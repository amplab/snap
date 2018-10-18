using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using System.IO;
using System.Diagnostics;
using ASELib;

namespace SpliceosomeAllelicImbalance
{
    class Program
    {
        static ASETools.CommonData commonData;
        static ASETools.IsoformMap isoformMap;

        // tumor->disease->hugo_symbol->isoform->spliceosomemutant->frac of all reads mapped to this gene that are in this isoform.  Empty string disease is all diseases.
        static Dictionary<bool, Dictionary<string, Dictionary<string, Dictionary<string, Dictionary<bool, List<double>>>>>> perDiseaseIsoformBalanceBySpliceosome = 
            new Dictionary<bool, Dictionary<string, Dictionary<string, Dictionary<string, Dictionary<bool, List<double>>>>>>();

        static List<string> diseases;
        static List<string> diseasesAndAll = new List<string>();

        class DoubleAndBool : IComparer<DoubleAndBool>
        {
            public bool whichClass;
            public double value;

            public DoubleAndBool(bool whichClass_, double value_)
            {
                whichClass = whichClass_;
                value = value_;
            }

            public int Compare(DoubleAndBool x, DoubleAndBool y)
            {
                return x.value.CompareTo(y.value);
            }
        }

        class Result
        {
            public string name;
            public double pValue;
            public double meanNoMutation;
            public double meanWithMutation;
            public int nNoMutation;
            public int nWithMutation;

            public Result(string name_, double pValue_, double meanNoMutation_, double meanWithMutation_, int nNoMutation_, int nWithMutation_)
            {
                name = name_;
                pValue = pValue_;
                meanNoMutation = meanNoMutation_;
                meanWithMutation = meanWithMutation_;
                nNoMutation = nNoMutation_;
                nWithMutation = nWithMutation_;
            }
        }

        static void Main(string[] args)
        {
            commonData = ASETools.CommonData.LoadCommonData(args);

            if (null == commonData)
            {
                return;
            }
            int minSamplesPerClass = 50;

            if (commonData.listOfCases.Any(_ => _.isoform_read_counts_filename == ""))
            {
                Console.WriteLine("Not all cases have isoform reads counts.");
                return;
            }

            var cases = commonData.cases;

#if false
cases = new Dictionary<string, ASETools.Case>();
foreach (var case_ in commonData.listOfCases.Where(_ => _.disease() == "laml"))
{
    cases.Add(case_.case_id, case_);
}
minSamplesPerClass = 1;
#endif // false

            if (cases.Select(_ => _.Value).Any(_ => _.tumor_rna_allcount_filename == ""))
            {
                Console.WriteLine("Not all cases have tumor RNA allcount files.");
                return;
            }

            Console.Write("Computing isoform map...");
            isoformMap = new ASETools.IsoformMap(commonData.geneLocationInformation.genesByName);
            Console.WriteLine(ASETools.ElapsedTimeInSeconds(commonData.timer));

            // initialize the data structures
            diseases = ASETools.GetListOfDiseases(cases);
            diseasesAndAll.Add(""); // the "all" is the empty string.
            diseases.ForEach(_ => diseasesAndAll.Add(_));

            Console.Write("Allocating data structures...");

            foreach (var tumor in ASETools.BothBools)
            {
                perDiseaseIsoformBalanceBySpliceosome.Add(tumor, new Dictionary<string, Dictionary<string, Dictionary<string, Dictionary<bool, List<double>>>>>());
                foreach (var disease in diseasesAndAll)
                {
                    perDiseaseIsoformBalanceBySpliceosome[tumor].Add(disease, new Dictionary<string, Dictionary<string, Dictionary<bool, List<double>>>>());

                    foreach (var geneInfo in commonData.geneLocationInformation.genesByName.Select(x => x.Value).ToList())
                    {
                        perDiseaseIsoformBalanceBySpliceosome[tumor][disease].Add(geneInfo.hugoSymbol, new Dictionary<string, Dictionary<bool, List<double>>>());

                        foreach (var isoformName in geneInfo.isoforms.Select(x => x.ucscId))
                        {
                            perDiseaseIsoformBalanceBySpliceosome[tumor][disease][geneInfo.hugoSymbol].Add(isoformName, new Dictionary<bool, List<double>>());

                            foreach (var mutantSpliceosome in ASETools.BothBools)
                            {
                                perDiseaseIsoformBalanceBySpliceosome[tumor][disease][geneInfo.hugoSymbol][isoformName].Add(mutantSpliceosome, new List<double>());
                            } // mutant
                        } // isoforms
                    } // genes
                } // disease
            } // tumor
            Console.WriteLine(ASETools.ElapsedTimeInSeconds(commonData.timer) + " overall elapsed time.");


            int nPerDot;
            ASETools.PrintMessageAndNumberBar("Processing", "cases", cases.Count(), out nPerDot);
            var threading = new ASETools.WorkerThreadHelper<ASETools.Case, PerThreadState>(cases.Select(_ => _.Value).ToList(), HandleOneCase, FinishUp, null, nPerDot);
            threading.run();

            Console.WriteLine(ASETools.ElapsedTimeInSeconds(commonData.timer));

            // all of these are tumor->disease->value
            var nTooFewSamples = new Dictionary<bool, Dictionary<string, int>>();
            var nConsidered = new Dictionary<bool, Dictionary<string, int>>();
            var results = new Dictionary<bool, Dictionary<string, List<Result>>>();
            var nSignificant = new Dictionary<bool, Dictionary<string, int>>();

            ASETools.BothBools.ToList().ForEach(tumor => { nTooFewSamples.Add(tumor, new Dictionary<string, int>()); nConsidered.Add(tumor, new Dictionary<string, int>()); results.Add(tumor, new Dictionary<string, List<Result>>()); nSignificant.Add(tumor, new Dictionary<string, int>()); });
            ASETools.BothBools.ToList().ForEach(tumor => diseasesAndAll.ForEach(disease => { nTooFewSamples[tumor].Add(disease, 0); nConsidered[tumor].Add(disease, 0); results[tumor].Add(disease, new List<Result>()); nSignificant[tumor].Add(disease, 0); }));

            var pValueHistograms = new Dictionary<bool, Dictionary<string, ASETools.PreBucketedHistogram>>();
            foreach (var tumor in ASETools.BothBools)
            {
                pValueHistograms.Add(tumor, new Dictionary<string, ASETools.PreBucketedHistogram>());
                foreach (var disease in diseasesAndAll)
                {
                    pValueHistograms[tumor].Add(disease, new ASETools.PreBucketedHistogram(0, 1, 0.01));
                }
            }

            ASETools.PrintMessageAndNumberBar("Processing", "genes", commonData.geneLocationInformation.genesByName.Count(), out nPerDot);

            int nGenesProcessed = 0;

            foreach (var geneInfo in commonData.geneLocationInformation.genesByName.Select(x => x.Value).ToList())
            {
                nGenesProcessed++;

                if (nGenesProcessed % nPerDot == 0)
                {
                    Console.Write(".");
                }

                foreach (var tumor in ASETools.BothBools)
                {
                    foreach (var disease in diseasesAndAll)
                    {
                        foreach (var isoformName in geneInfo.isoforms.Select(x => x.ucscId))
                        {
                            if (perDiseaseIsoformBalanceBySpliceosome[tumor][disease][geneInfo.hugoSymbol][isoformName].Select(_ => _.Value.Count()).Min() < minSamplesPerClass)
                            {
                                nTooFewSamples[tumor][disease]++;
                                continue;
                            }   // This is also necessarily true for each disease, so skipping the isoform entirely is OK.

                            var items = new List<DoubleAndBool>();
                            var perDiseaseItems = new Dictionary<string, List<DoubleAndBool>>();
                            diseases.ForEach(_ => perDiseaseItems.Add(_, new List<DoubleAndBool>()));

                            foreach (var mutant in ASETools.BothBools)
                            {
                                perDiseaseIsoformBalanceBySpliceosome[tumor][disease][geneInfo.hugoSymbol][isoformName][mutant].ForEach(_ => items.Add(new DoubleAndBool(mutant, _)));
                            }

                            nConsidered[tumor][disease]++;

                            bool enoughData, reversed;
                            double nFirstGroup, nSecondGroup, U, z;

                            double p = ASETools.MannWhitney<DoubleAndBool>.ComputeMannWhitney(items, items[0], _ => _.whichClass, _ => _.value, out enoughData, out reversed, out nFirstGroup, out nSecondGroup, out U, out z);

                            pValueHistograms[tumor][disease].addValue(p);

                            var nonMutantItems = items.Where(_ => !_.whichClass).Select(_ => _.value);
                            var mutantItems = items.Where(_ => _.whichClass).Select(_ => _.value);

                            results[tumor][disease].Add(new Result(geneInfo.hugoSymbol + ":" + isoformName, p, nonMutantItems.Average(), mutantItems.Average(), nonMutantItems.Count(), mutantItems.Count()));

                        } // for each isoform
                    } // for each disease
                } // tumor
            } // for each gene

            Console.WriteLine();

            var pValueHistogramOutputFilename = commonData.configuration.finalResultsDirectory + ASETools.IsoformBalancePValueHistogramFilename;
            var pValueHistogramFile = ASETools.CreateStreamWriterWithRetry(pValueHistogramOutputFilename);
            if (null == pValueHistogramFile)
            {
                Console.WriteLine("Unable to open p value histogram file " + pValueHistogramOutputFilename);
                pValueHistogramFile = StreamWriter.Null;
                // Fall through
            }

            foreach (var tumor in ASETools.BothBools)
            {
                foreach (var disease in diseasesAndAll)
                {
                    nSignificant[tumor][disease] = results[tumor][disease].Where(_ => _.pValue * nConsidered[tumor][disease] <= .01).Count();

                    results[tumor][disease].Sort((x, y) => y.pValue.CompareTo(x.pValue));   // Backwards sort, so smallest is last

                    double maxSignificantUncorrectedPValue = 0.01 / nConsidered[tumor][disease];

                    string outputFilename;
                    string tumorString = tumor ? "tumor_" : "normal_";
                    if (disease == "")
                    {
                        outputFilename = commonData.configuration.finalResultsDirectory + tumorString + ASETools.IsoformBalanceFilename;
                    }
                    else
                    {
                        outputFilename = commonData.configuration.finalResultsDirectory + disease + "_" + tumorString + ASETools.IsoformBalanceFilename;
                    }

                    var outputFile = ASETools.CreateStreamWriterWithRetry(outputFilename);
                    if (null == outputFile)
                    {
                        Console.WriteLine("Unable to open output file " + outputFilename);
                        //
                        // Just let it go and print the last line anyway.
                        //
                    }
                    else
                    {
                        //
                        // Sort so significant results go first.  Within significant results, sort by ascending p value.  For insignificant results, sort by name.
                        //
                        results[tumor][disease].Sort((x, y) => (x.pValue < maxSignificantUncorrectedPValue && y.pValue < maxSignificantUncorrectedPValue) ? x.pValue.CompareTo(y.pValue) : (x.pValue < maxSignificantUncorrectedPValue ? -1 : (y.pValue < maxSignificantUncorrectedPValue ? 1 : x.name.CompareTo(y.name))));

                        outputFile.WriteLine("Name\tCorrected pValue\tn with no mutation\tmean with no mutation\tn with mutation\tmean with mutation");
                        results[tumor][disease].ForEach(_ => outputFile.WriteLine(_.name + "\t" + (_.pValue * nConsidered[tumor][disease]) + "\t" + _.nNoMutation + "\t" + _.meanNoMutation + "\t" + _.nWithMutation + "\t" + _.meanWithMutation));
                        outputFile.WriteLine("**done**");
                        outputFile.Close();
                    }

                    if (disease == "")
                    {
                        pValueHistogramFile.WriteLine("p value histogram for entire dataset for " + tumorString + ".");
                    }
                    else
                    {
                        pValueHistogramFile.WriteLine();
                        pValueHistogramFile.WriteLine("p value histogram for " + disease + " for " + tumorString);
                    }
                    pValueHistogramFile.WriteLine(ASETools.HistogramResultLine.Header());
                    pValueHistograms[tumor][disease].ComputeHistogram().ToList().ForEach(_ => pValueHistogramFile.WriteLine(_.ToString()));
                } // disease
            } // tumor

            pValueHistogramFile.WriteLine("**done**");
            pValueHistogramFile.Close();

            Console.WriteLine("Considered " + nConsidered[true][""] + " isoforms, of which " + nSignificant[true][""] + " are have significant differences in balance based on spliceosome mutations at the 0.01 level after Bonferroni correction in " + ASETools.ElapsedTimeInSeconds(commonData.timer));
        } // Main

        class PerThreadState
        {
            // tumor->hugo_symbol->isoform->spliceosomemutant->frac of all reads mapped to this gene that are in this isoform.  Empty string disease is all diseases.
            public Dictionary<bool,Dictionary<string, Dictionary<string, Dictionary<bool, List<double>>>>> overallIsoformBalanceBySpliceosome =
                new Dictionary<bool, Dictionary<string, Dictionary<string, Dictionary<bool, List<double>>>>>();

            public PerThreadState()
            {
                foreach (var tumor in ASETools.BothBools)
                {
                    overallIsoformBalanceBySpliceosome.Add(tumor, new Dictionary<string, Dictionary<string, Dictionary<bool, List<double>>>>());

                    foreach (var geneInfo in commonData.geneLocationInformation.genesByName.Select(x => x.Value).ToList())
                    {
                        overallIsoformBalanceBySpliceosome[tumor].Add(geneInfo.hugoSymbol, new Dictionary<string, Dictionary<bool, List<double>>>());

                        foreach (var isoformName in geneInfo.isoforms.Select(x => x.ucscId))
                        {
                            overallIsoformBalanceBySpliceosome[tumor][geneInfo.hugoSymbol].Add(isoformName, new Dictionary<bool, List<double>>());

                            foreach (var mutantSpliceosome in ASETools.BothBools)
                            {
                               overallIsoformBalanceBySpliceosome[tumor][geneInfo.hugoSymbol][isoformName].Add(mutantSpliceosome, new List<double>());
                            } // mutant
                        } // isoforms
                    } // genes
                } // tumor
            } // ctor
        } // PerThreadState

        static void HandleOneCase(ASETools.Case case_, PerThreadState state)
        {
            var mafLines = ASETools.MAFLine.ReadFile(case_.extracted_maf_lines_filename, "", false);
            bool anySpliceosomeMutations = mafLines.Any(_ => !_.IsSilent() && ASETools.spliceosome_genes.Contains(_.Hugo_Symbol));

            //var splicesomeMutations = mafLines.Where(_ => !_.IsSilent() && ASETools.spliceosome_genes.Contains(_.Hugo_Symbol)).ToList();

            var isoformReadCounts = ASETools.IsoformReadCounts.LoadFromFile(case_.isoform_read_counts_filename);
            if (isoformReadCounts == null)
            {
                Console.WriteLine("Unable to load isoform read counts for case " + case_.case_id + " from " + case_.isoform_read_counts_filename);
                return;
            }

            foreach (var tumor in ASETools.BothBools)
            {
                var totalIsoformReadsByGene = new Dictionary<string, long>();
                foreach (var geneInfo in commonData.geneLocationInformation.genesByName.Select(_ => _.Value).ToList())
                {
                    long totalIsoformReads = 0;
                    geneInfo.isoforms.Select(_ => _.ucscId).ToList().ForEach(_ => totalIsoformReads += isoformReadCounts[_].getCount(tumor));   // This is in place of the more natural .Sum() form to avoid integer overflow

                    totalIsoformReadsByGene.Add(geneInfo.hugoSymbol, totalIsoformReads);
                } // per gene

                lock (perDiseaseIsoformBalanceBySpliceosome[tumor][case_.disease()])
                {
                    foreach (var geneInfo in commonData.geneLocationInformation.genesByName.Select(_ => _.Value).ToList())
                    {
                        if (totalIsoformReadsByGene[geneInfo.hugoSymbol] < 100)
                        {
                            continue;
                        }

                        foreach (var ucsdName in geneInfo.isoforms.Select(_ => _.ucscId))
                        {
                            perDiseaseIsoformBalanceBySpliceosome[tumor][case_.disease()][geneInfo.hugoSymbol][ucsdName][anySpliceosomeMutations].Add((double)isoformReadCounts[ucsdName].getCount(tumor) / totalIsoformReadsByGene[geneInfo.hugoSymbol]);

                            state.overallIsoformBalanceBySpliceosome[tumor][geneInfo.hugoSymbol][ucsdName][anySpliceosomeMutations].Add((double)isoformReadCounts[ucsdName].getCount(tumor) / totalIsoformReadsByGene[geneInfo.hugoSymbol]);
                        }
                    } // per gene
                } // lock
            } // tumor

        } // HandleOneCase

        static void FinishUp(PerThreadState state)
        {
            foreach (var tumor in ASETools.BothBools)
            {
                lock (perDiseaseIsoformBalanceBySpliceosome[tumor][""])
                {
                    foreach (var geneInfo in commonData.geneLocationInformation.genesByName.Select(_ => _.Value).ToList())
                    {
                        foreach (var ucsdName in geneInfo.isoforms.Select(_ => _.ucscId))
                        {
                            foreach (var anySpliceosomeMutations in ASETools.BothBools)
                            {
                                if (state.overallIsoformBalanceBySpliceosome[tumor][geneInfo.hugoSymbol][ucsdName].ContainsKey(anySpliceosomeMutations))
                                {
                                    perDiseaseIsoformBalanceBySpliceosome[tumor][""][geneInfo.hugoSymbol][ucsdName][anySpliceosomeMutations].AddRange(
                                        state.overallIsoformBalanceBySpliceosome[tumor][geneInfo.hugoSymbol][ucsdName][anySpliceosomeMutations]);
                                } // anything to add

                            } // mutated or not
                        } // ucsdName
                    } // per gene
                } // lock
            } // tumor
        } // FinishUp
    }
}
