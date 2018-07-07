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

        static Dictionary<string, ASETools.PreBucketedHistogram> global_one_mutation_histograms = new Dictionary<string, ASETools.PreBucketedHistogram>();
        static Dictionary<string, ASETools.PreBucketedHistogram> global_multiple_mutation_histograms = new Dictionary<string, ASETools.PreBucketedHistogram>();
        static Dictionary<string, ASETools.PreBucketedHistogram> global_normal_histograms = new Dictionary<string, ASETools.PreBucketedHistogram>();
        static Dictionary<string, ASETools.PreBucketedHistogram> global_silent_histograms = new Dictionary<string, ASETools.PreBucketedHistogram>();

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

            var casesToProcess = cases.Select(x => x.Value).Where(x => x.annotated_selected_variants_filename != "").ToList();

            Console.Write("Processing " + casesToProcess.Count() + " cases, one dot/100 cases: ");

            var threading = new ASETools.ASVThreadingHelper<List<Dictionary<string, ASETools.PreBucketedHistogram>>>(casesToProcess, null, (x, y) => true, HandleOneCase, FinishUp, null, 100);
            threading.run(1);

            Console.WriteLine();

            const int minCount = 20;
            int bonferroniCorrection = global_one_mutation_histograms.Where(x => x.Value.count() >= minCount).Count(); // The set of genes with at least 20 tumors with single mutations
            Console.WriteLine("Bonferroni correction " + bonferroniCorrection);


            var results = new List<HistogramsAndPValue>();

            foreach (var geneEntry in global_one_mutation_histograms.Where(x => x.Value.count() >= minCount))
            {
                var hugo_symbol = geneEntry.Key;
                var histogram = geneEntry.Value;

                results.Add(new HistogramsAndPValue(ASETools.binomialTest((int)histogram.nLessThan(0.5), (int)histogram.count(), 0.5, ASETools.BinomalTestType.TwoSided) * bonferroniCorrection, histogram, 
                    global_multiple_mutation_histograms[hugo_symbol], global_normal_histograms[hugo_symbol], global_silent_histograms[hugo_symbol], hugo_symbol));

            }

            results.Sort();

            var outputFile = ASETools.CreateStreamWriterWithRetry(configuration.finalResultsDirectory + ASETools.vaf_histogram_filename);

            foreach (var result in results)
            {
                if (result.p <= 0.01)
                {
                    Console.WriteLine(result.hugo_symbol + " has p=" + result.p);
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

        const int nHistograms = 4;

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

                int mutationCount = variants.Where(x => x.somaticMutation && !x.isSilent()).Count();

                if (mutationCount == 1 && variants[0].IsASECandidate(true, copyNumber, configuration, perGeneASEMap, geneMap) && !variants[0].CausesNonsenseMediatedDecay() && !variants[0].isSilent())
                {
                    histograms[0][hugo_symbol].addValue(variants[0].GetTumorAltAlleleFraction());
                }

                if (mutationCount > 1 && variants.Where(x => x.somaticMutation && x.IsASECandidate(true, copyNumber, configuration, perGeneASEMap, geneMap) && !x.isSilent()).Count() > 0)
                {
                    var vaf = variants.Where(x => x.somaticMutation && x.IsASECandidate(true, copyNumber, configuration, perGeneASEMap, geneMap) && !x.isSilent()).Select(x => x.GetTumorAltAlleleFraction()).Average();

                    histograms[1][hugo_symbol].addValue(vaf);
                }

                if (variants.Where(x => !x.somaticMutation && x.IsASECandidate(false, copyNumber, configuration, perGeneASEMap, geneMap)).Count() > 0)
                {
                    var vaf = variants.Where(x => !x.somaticMutation && x.IsASECandidate(false, copyNumber, configuration, perGeneASEMap, geneMap)).Select(x => x.GetNormalAltAlleleFraction()).Average();
                    histograms[2][hugo_symbol].addValue(vaf);
                }

                if (mutationCount == 0 && variants.Where(x => x.isSilent() && x.somaticMutation && x.IsASECandidate(true, copyNumber, configuration, perGeneASEMap, geneMap)).Count() > 0)
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
                var hugo_symbol = geneEntry.Key;
                var histogram = geneEntry.Value;

                if (!global.ContainsKey(hugo_symbol))
                {
                    global.Add(hugo_symbol, histogram);
                } else
                {
                    global[hugo_symbol].merge(histogram);
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
            } // lock
        } // FinishUp
    }
}
