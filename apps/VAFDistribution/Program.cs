using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using System.IO;
using System.Diagnostics;
using ASELib;

namespace VAFDistribution
{
    class Program
    {
        static Dictionary<string, ASETools.Case> cases;
        static ASETools.Configuration configuration;
        static Dictionary<string, ASETools.ASEMapPerGeneLine> perGeneASEMap;
        static ASETools.GeneMap geneMap;
        static ASETools.ASERepetitiveRegionMap repetitiveRegionMap;

        static Dictionary<string, ASETools.PreBucketedHistogram> global_one_mutation_histograms = new Dictionary<string, ASETools.PreBucketedHistogram>();
        static Dictionary<string, ASETools.PreBucketedHistogram> global_multiple_mutation_histograms = new Dictionary<string, ASETools.PreBucketedHistogram>();
        static Dictionary<string, ASETools.PreBucketedHistogram> global_normal_histograms = new Dictionary<string, ASETools.PreBucketedHistogram>();
        static Dictionary<string, ASETools.PreBucketedHistogram> global_silent_histograms = new Dictionary<string, ASETools.PreBucketedHistogram>();
        static Dictionary<string, ASETools.PreBucketedHistogram> global_idh1_one_mutation_histograms = new Dictionary<string, ASETools.PreBucketedHistogram>();  // Maps IDH1 mutation type->historgam
        static Dictionary<string, ASETools.PreBucketedHistogram> global_idh1_normal_histograms = new Dictionary<string, ASETools.PreBucketedHistogram>();  // Maps IDH1 mutation type->historgam
        static Dictionary<string, ASETools.PreBucketedHistogram> global_idh1_multiple_mutation_histograms = new Dictionary<string, ASETools.PreBucketedHistogram>();  // Maps IDH1 mutation type->historgam

        class HistogramsAndPValue : IComparable<HistogramsAndPValue>
        {
            public readonly double p;
            public readonly ASETools.PreBucketedHistogram oneMutationHistogram;
            public readonly ASETools.PreBucketedHistogram multipleMutationHistogram;
            public readonly ASETools.PreBucketedHistogram normalHistogram;
            public readonly ASETools.PreBucketedHistogram silentHistogram;
            public readonly string hugo_symbol;

            public bool isSignificant()
            {
                return p <= 0.01;
            }

            public int CompareTo(HistogramsAndPValue peer)
            {
                return p.CompareTo(peer.p);
            }

            public HistogramsAndPValue(double p_, ASETools.PreBucketedHistogram oneMutationHistogram_, ASETools.PreBucketedHistogram multipleMutationHistogram_, ASETools.PreBucketedHistogram normalHistogram_, ASETools.PreBucketedHistogram silentHistogram_, string hugo_symbol_)
            {
                p = p_;
                oneMutationHistogram = oneMutationHistogram_;
                multipleMutationHistogram = multipleMutationHistogram_;
                normalHistogram = normalHistogram_;
                silentHistogram = silentHistogram_;
                hugo_symbol = hugo_symbol_;
            }
        }

        static void Main(string[] args)
        {
            var timer = new Stopwatch();
            timer.Start();

            configuration = ASETools.Configuration.loadFromFile(args);

            if (null == configuration)
            {
                Console.WriteLine("Unable to load configuration");
                return;
            }

            cases = ASETools.Case.LoadCases(configuration.casesFilePathname);
            if (null == cases)
            {
                Console.WriteLine("Unable to load cases");
                return;
            }

            var geneLocationInformation = new ASETools.GeneLocationsByNameAndChromosome(ASETools.readKnownGeneFile(configuration.geneLocationInformationFilename));
            geneMap = new ASETools.GeneMap(geneLocationInformation.genesByName);
            perGeneASEMap = ASETools.ASEMapPerGeneLine.ReadFromFileToDictionary(configuration.finalResultsDirectory + ASETools.PerGeneASEMapFilename);

            repetitiveRegionMap = ASETools.ASERepetitiveRegionMap.loadFromFile(configuration.redundantChromosomeRegionFilename);
            if (null == repetitiveRegionMap)
            {
                Console.WriteLine("Unable to load repetitive region map from " + configuration.redundantChromosomeRegionFilename);
                return;
            }

            var casesToProcess = cases.Select(x => x.Value).Where(x => x.annotated_selected_variants_filename != "").ToList();

            Console.WriteLine("Processing " + casesToProcess.Count() + " cases, one dot/100 cases: ");
            ASETools.PrintNumberBar(casesToProcess.Count() / 100);

            var threading = new ASETools.ASVThreadingHelper<List<Dictionary<string, ASETools.PreBucketedHistogram>>>(casesToProcess, null, (x, y) => true, HandleOneCase, FinishUp, null, 100);
            threading.run();

            Console.WriteLine();

            const int minCount = 20;
            int bonferroniCorrection = global_one_mutation_histograms.Where(x => x.Value.count() >= minCount).Count();      // The set of genes with at least 20 tumors with single mutations
            bonferroniCorrection += global_idh1_one_mutation_histograms.Where(x => x.Value.count() >= minCount).Count();     // Or of specific IDH1 mutations
            Console.WriteLine("Bonferroni correction " + bonferroniCorrection);


            var results = new List<HistogramsAndPValue>();

            foreach (var geneEntry in global_one_mutation_histograms.Where(x => x.Value.count() >= minCount))
            {
                var hugo_symbol = geneEntry.Key;
                var histogram = geneEntry.Value;

                results.Add(new HistogramsAndPValue(ASETools.binomialTest((int)histogram.nLessThan(0.5), (int)histogram.count(), 0.5, ASETools.BinomalTestType.TwoSided) * bonferroniCorrection, histogram, 
                    global_multiple_mutation_histograms[hugo_symbol], global_normal_histograms[hugo_symbol], global_silent_histograms[hugo_symbol], hugo_symbol));

            }

            foreach (var particularMutation in global_idh1_one_mutation_histograms)
            {
                var mutation = particularMutation.Key;
                var histogram = particularMutation.Value;

                if (histogram.count() < 1)
                {
                    continue;
                }

                results.Add(new HistogramsAndPValue(ASETools.binomialTest((int)histogram.nLessThan(0.5), (int)histogram.count(), 0.5, ASETools.BinomalTestType.TwoSided) * bonferroniCorrection, histogram,
                    global_idh1_multiple_mutation_histograms[mutation], global_idh1_normal_histograms[mutation], new ASETools.PreBucketedHistogram(0, 1.01, 0.01), "IDH1:" + mutation));
            }

            results.Sort();

            var outputFile = ASETools.CreateStreamWriterWithRetry(configuration.finalResultsDirectory + ASETools.vaf_histogram_filename);

            foreach (var result in results)
            {
                if (result.p <= 0.01)
                {
                    Console.WriteLine(result.hugo_symbol + " has p=" + result.p + ", median 1 mutation VAF of " + result.oneMutationHistogram.ComputeHistogram()[50].cdfValue);
                }

                var oneMutationLines = result.oneMutationHistogram.ComputeHistogram();
                var multiMutationLines = result.multipleMutationHistogram.ComputeHistogram();
                var normalLines = result.normalHistogram.ComputeHistogram();
                var silentLines = result.silentHistogram.ComputeHistogram();

                outputFile.WriteLine("hugo_symbol\tp (one mutation not coin flip by binomial test, bonferroni corrected)\tn (one mutation)\tmean (one mutation)");
                outputFile.WriteLine(result.hugo_symbol + "\t" + result.p + "\t" + result.oneMutationHistogram.count() + "\t" + result.oneMutationHistogram.mean());
                outputFile.WriteLine("minValue\tone mutation (n = " + result.oneMutationHistogram.count() + ")\tmultiple mutations (n = " + result.multipleMutationHistogram.count() + ")\tnormal (n = " + result.normalHistogram.count() + 
                    ")\tsilent (n = " + result.silentHistogram.count() + ")");
                for (int i = 0; i < oneMutationLines.Count(); i++)
                {
                    outputFile.WriteLine(oneMutationLines[i].minValue + "\t" + oneMutationLines[i].cdfValue + "\t" + multiMutationLines[i].cdfValue + "\t" + normalLines[i].cdfValue + "\t" + silentLines[i].cdfValue);
                }
 
                outputFile.WriteLine();
            } // foreach result

            outputFile.Close();

            Console.WriteLine("Run time " + ASETools.ElapsedTimeInSeconds(timer));
        } // Main

        static string getHugoSymbolForVariant(ASETools.AnnotatedVariant variant)
        {
            if (variant.Hugo_symbol != "")
            {
                return variant.Hugo_symbol;
            }

            var mappedGenes = geneMap.getGenesMappedTo(variant.contig, variant.locus);

            if (mappedGenes.Count() == 1)
            {
                return mappedGenes[0].hugoSymbol;
            }

            return "";
        }

        const int nHistograms = 7;

        static void HandleOneCase(ASETools.Case case_, List<Dictionary<string, ASETools.PreBucketedHistogram>> histograms, List<ASETools.AnnotatedVariant> variantsForThisCase, ASETools.ASECorrection aseCorrection, Dictionary<bool, List<ASETools.CopyNumberVariation>> copyNumber)
        {
            var variantsByGene = variantsForThisCase.GroupByToDict(x => getHugoSymbolForVariant(x));

            if (histograms.Count() == 0)
            {
                for (int i = 0; i < nHistograms; i++)
                {
                    histograms.Add(new Dictionary<string, ASETools.PreBucketedHistogram>());
                }
            }

            foreach (var geneEntry in variantsByGene)
            {
                var variants = geneEntry.Value;
                var hugo_symbol = geneEntry.Key;

                if (hugo_symbol == "")
                {
                    //
                    // These are germline variants that aren't in exactly one gene.
                    //
                    continue;
                }

                if (!histograms[0].ContainsKey(hugo_symbol))
                {
                    for (int i = 0; i < nHistograms; i++)
                    {
                        histograms[i].Add(hugo_symbol, new ASETools.PreBucketedHistogram(0, 1.01, 0.01));
                    }
                }

                var mutations = variants.Where(x => x.somaticMutation && !x.isSilent()).ToList();
                int mutationCount = mutations.Count();

                if (hugo_symbol == "IDH1")
                {
                    foreach (var mutation in mutations)
                    {
                        var mutationType = mutation.IDH1MutationType();

                        if (!histograms[4].ContainsKey(mutationType))
                        {
                            histograms[4].Add(mutationType, new ASETools.PreBucketedHistogram(0, 1.01, .01));
                            histograms[5].Add(mutationType, new ASETools.PreBucketedHistogram(0, 1.01, .01));
                            histograms[6].Add(mutationType, new ASETools.PreBucketedHistogram(0, 1.01, .01));
                        }
                    }
                }

                if (mutationCount == 1 && mutations[0].IsASECandidate(true, copyNumber, configuration, perGeneASEMap, geneMap) && !mutations[0].CausesNonsenseMediatedDecay() && !mutations[0].isSilent() && 
                    !repetitiveRegionMap.isInRepetitiveRegion(mutations[0].contig, mutations[0].locus))
                {
#if false
                    if (mutations[0].Hugo_symbol == "COL1A2")
                    {
                        Console.WriteLine("COL1A2 case " + case_.case_id + " single mutation, tumor DNA " + mutations[0].tumorDNAReadCounts.nMatchingReference + " matching reference, " + mutations[0].tumorDNAReadCounts.nMatchingAlt +
                            " matching alt, tumor RNA " + mutations[0].tumorRNAReadCounts.nMatchingReference + " matching reference (" + mutations[0].reference_allele +"), " + mutations[0].tumorRNAReadCounts.nMatchingAlt + 
                            " matching alt (" + mutations[0].alt_allele + "), locus " + 
                            mutations[0].contig + ":" + ASETools.NumberWithCommas(mutations[0].locus));
                    }
#endif
                    histograms[0][hugo_symbol].addValue(mutations[0].GetTumorAltAlleleFraction());

                    if (hugo_symbol == "IDH1")
                    {
                        var mutationType = mutations[0].IDH1MutationType();

                        histograms[4][mutationType].addValue(mutations[0].GetTumorAltAlleleFraction());
                    }
                }

                if (mutationCount > 1 && mutations.Where(x => x.somaticMutation && x.IsASECandidate(true, copyNumber, configuration, perGeneASEMap, geneMap) && !x.isSilent() && !x.CausesNonsenseMediatedDecay()
                    && !repetitiveRegionMap.isInRepetitiveRegion(x.contig, x.locus)).Count() > 0)
                {
                    var vaf = mutations.Where(x => x.somaticMutation && x.IsASECandidate(true, copyNumber, configuration, perGeneASEMap, geneMap) && !x.isSilent() && !x.CausesNonsenseMediatedDecay()
                    && !repetitiveRegionMap.isInRepetitiveRegion(x.contig, x.locus)).Select(x => x.GetTumorAltAlleleFraction()).Average();

                    histograms[1][hugo_symbol].addValue(vaf);

                    if (hugo_symbol == "IDH1")
                    {
                        foreach (var mutation in mutations)
                        {
                            if (mutation.somaticMutation && mutation.IsASECandidate(true, copyNumber, configuration, perGeneASEMap, geneMap) && !mutation.isSilent() && !mutation.CausesNonsenseMediatedDecay() &&
                                !repetitiveRegionMap.isInRepetitiveRegion(mutation.contig, mutation.locus))
                            {
                                histograms[5][mutation.IDH1MutationType()].addValue(mutation.GetTumorAltAlleleFraction());
                            } // if it's qualified
                        } // each mutation
                    } // IDH1
                }

                if (variants.Where(x => !x.somaticMutation && x.IsASECandidate(false, copyNumber, configuration, perGeneASEMap, geneMap)).Count() > 0)  // Don't need to check repetitveRegionMap for germline variants, that happened in selectVariants
                {
                    var vaf = variants.Where(x => !x.somaticMutation && x.IsASECandidate(false, copyNumber, configuration, perGeneASEMap, geneMap)).Select(x => x.GetNormalAltAlleleFraction()).Average();
                    histograms[2][hugo_symbol].addValue(vaf);

                    if (hugo_symbol == "IDH1")
                    {
                        if (!histograms[4].ContainsKey("Other"))
                        {
                            histograms[4].Add("Other", new ASETools.PreBucketedHistogram(0, 1.01, .01));
                            histograms[5].Add("Other", new ASETools.PreBucketedHistogram(0, 1.01, .01));
                            histograms[6].Add("Other", new ASETools.PreBucketedHistogram(0, 1.01, .01));
                        }
                        histograms[6]["Other"].addValue(vaf);
                    }
                }

                if (mutationCount == 0 && variants.Where(x => x.isSilent() && x.somaticMutation && x.IsASECandidate(true, copyNumber, configuration, perGeneASEMap, geneMap)).Count() > 0) // how can this happen?  There are no silent/non-silent tags on germline variants
                {
                    var vaf = variants.Where(x => x.isSilent() && x.somaticMutation && x.IsASECandidate(true, copyNumber, configuration, perGeneASEMap, geneMap)).Select(x => x.GetTumorAltAlleleFraction()).Average();
                    histograms[3][hugo_symbol].addValue(vaf);
                }
            } //foreach gene
        } // HandleOneCase

        static void mergeHistograms(Dictionary<string, ASETools.PreBucketedHistogram> global, Dictionary<string, ASETools.PreBucketedHistogram> local)
        {
            foreach (var geneEntry in local)
            {
                var hugo_symbol_or_mutation = geneEntry.Key;
                var histogram = geneEntry.Value;

                if (!global.ContainsKey(hugo_symbol_or_mutation))
                {
                    global.Add(hugo_symbol_or_mutation, histogram);
                } else
                {
                    global[hugo_symbol_or_mutation].merge(histogram);
                }
            }
        }

        static void FinishUp(List<Dictionary<string, ASETools.PreBucketedHistogram>> histograms)
        {
            lock (global_one_mutation_histograms)
            {
                mergeHistograms(global_one_mutation_histograms, histograms[0]);
                mergeHistograms(global_multiple_mutation_histograms, histograms[1]);
                mergeHistograms(global_normal_histograms, histograms[2]);
                mergeHistograms(global_silent_histograms, histograms[3]);
                mergeHistograms(global_idh1_one_mutation_histograms, histograms[4]);
                mergeHistograms(global_idh1_multiple_mutation_histograms, histograms[5]);
                mergeHistograms(global_idh1_normal_histograms, histograms[6]);
            } // lock
        } // FinishUp
    }
}
