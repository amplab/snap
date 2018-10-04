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
        // hugo_symbol->isoform->spliceosomemutant->frac of all reads mapped to this gene that are in this isoform
        static Dictionary<string, Dictionary<string, Dictionary<bool, List<double>>>> isoformBalanceBySpliceosome = new Dictionary<string, Dictionary<string, Dictionary<bool, List<double>>>>();   


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

            if (commonData.cases.Select(_ => _.Value).Any(_ => _.tumor_rna_allcount_filename == ""))
            {
                Console.WriteLine("Not all cases have tumor RNA allcount files.");
                return;
            }

            isoformMap = new ASETools.IsoformMap(commonData.geneLocationInformation.genesByName);

            // initialize the data structure
            foreach (var geneInfo in commonData.geneLocationInformation.genesByName.Select(x => x.Value).ToList())
            {
                isoformBalanceBySpliceosome.Add(geneInfo.hugoSymbol, new Dictionary<string, Dictionary<bool, List<double>>>());
                foreach (var isoformName in geneInfo.isoforms.Select(x => x.ucscId))
                {
                    isoformBalanceBySpliceosome[geneInfo.hugoSymbol].Add(isoformName, new Dictionary<bool, List<double>>());
                    foreach (var mutantSpliceosome in ASETools.BothBools)
                    {
                        isoformBalanceBySpliceosome[geneInfo.hugoSymbol][isoformName].Add(mutantSpliceosome, new List<double>());
                    } // mutant
                } // isoforms
            } // genes

            var threading = new ASETools.WorkerThreadHelper<ASETools.Case, int>(commonData.cases.Select(_ => _.Value).ToList(), HandleOneCase, null, null, 100);
            Console.WriteLine("Processing allcount files, 1 dot/100 cases:");
            ASETools.PrintNumberBar(commonData.cases.Count() / 100);
            threading.run();

            int minSamplesPerClass = 10;

            int nTooFewSamples = 0;
            int nConsidered = 0;
            var results = new List<Result>();

            var pValueHistogram = new ASETools.PreBucketedHistogram(0, 1, 0.001);

            int nPerDot;
            ASETools.PrintMessageAndNumberBar("Processing", "genes", commonData.geneLocationInformation.genesByName.Count(), out nPerDot);

            int nGenesProcessed = 0;

            foreach (var geneInfo in commonData.geneLocationInformation.genesByName.Select(x => x.Value).ToList())
            {
                nGenesProcessed++;

                if (nGenesProcessed % nPerDot == 0)
                {
                    Console.Write(".");
                }
                foreach (var isoformName in geneInfo.isoforms.Select(x => x.ucscId))
                {
                    if (isoformBalanceBySpliceosome[geneInfo.hugoSymbol][isoformName].Select(_ => _.Value.Count()).Min() < minSamplesPerClass)
                    {
                        nTooFewSamples++;
                        continue;
                    }

                    var items = new List<DoubleAndBool>();
                    foreach (var mutant in ASETools.BothBools)
                    {
                        isoformBalanceBySpliceosome[geneInfo.hugoSymbol][isoformName][mutant].ForEach(_ => items.Add(new DoubleAndBool(mutant, _)));
                    }

                    nConsidered++;

                    bool enoughData, reversed;
                    double nFirstGroup, nSecondGroup, U, z;

                    double p = ASETools.MannWhitney<DoubleAndBool>.ComputeMannWhitney(items, items[0], _ => _.whichClass, _ => _.value, out enoughData, out reversed, out nFirstGroup, out nSecondGroup, out U, out z);

                    pValueHistogram.addValue(p);

                    var nonMutantItems = items.Where(_ => !_.whichClass).Select(_ => _.value);
                    var mutantItems = items.Where(_ => _.whichClass).Select(_ => _.value);

                    results.Add(new Result(geneInfo.hugoSymbol + ":" + isoformName, p, nonMutantItems.Average(), mutantItems.Average(), nonMutantItems.Count(), mutantItems.Count()));
                } // for each isoform
            } // for each gene

            Console.WriteLine();

            int nSignificant = results.Where(_ => _.pValue * nConsidered <= .01).Count();

            results.Sort((x, y) => y.pValue.CompareTo(x.pValue));   // Backwards sort, so smallest is last

            double maxSignificantUncorrectedPValue = 0.01 / nConsidered;
            results.Where(_ => _.pValue <= maxSignificantUncorrectedPValue).ToList().ForEach(_ => Console.WriteLine(_.name + " has pValue " + (_.pValue * nConsidered) + " after Bonferroni correction."));

            var outputFilename = commonData.configuration.finalResultsDirectory + ASETools.IsoformBalanceFilename;
            var outputFile = ASETools.CreateStreamWriterWithRetry(outputFilename);
            if (null == outputFile)
            {
                Console.WriteLine("Unable to open output file " + outputFilename);
                //
                // Just let it go and print the last line anyway.
                //
            } else
            {
                //
                // Sort so significant results go first.  Within significant results, sort by ascending p value.  For insignificant results, sort by name.
                //
                results.Sort((x, y) => (x.pValue < maxSignificantUncorrectedPValue && y.pValue < maxSignificantUncorrectedPValue) ? x.pValue.CompareTo(y.pValue) : (x.pValue < maxSignificantUncorrectedPValue ? -1 : (y.pValue < maxSignificantUncorrectedPValue ? 1 : x.name.CompareTo(y.name))));

                outputFile.WriteLine("Name\tCorrected pValuet\tn with no mutation\tmean with no mutation\tn with mutation\tmean with mutation");
                results.ForEach(_ => outputFile.WriteLine(_.name + "\t" + (_.pValue * nConsidered) + "\t" + _.nNoMutation + "\t" + _.meanNoMutation+ "\t" + _.nWithMutation  + "\t" + _.meanWithMutation));
                outputFile.WriteLine("**done**");
                outputFile.Close();
            }

            var pValueHistogramOutputFilename = commonData.configuration.finalResultsDirectory + ASETools.IsoformBalancePValueHistogramFilename;
            var pValueHistogramFile = ASETools.CreateStreamWriterWithRetry(pValueHistogramOutputFilename);
            if (null == pValueHistogramFile)
            {
                Console.WriteLine("Unable to open p value histogram file " + pValueHistogramOutputFilename);
                // Fall through
            } else
            {
                pValueHistogramFile.WriteLine(ASETools.HistogramResultLine.Header());
                pValueHistogram.ComputeHistogram().ToList().ForEach(_ => pValueHistogramFile.WriteLine(_.ToString()));
                pValueHistogramFile.WriteLine("**done**");
                pValueHistogramFile.Close();
            }

            Console.WriteLine("Considered " + nConsidered + " isoforms, of which " + nSignificant + " are have significant differences in balance based on spliceosome mutations at the 0.01 level after Bonferroni correction in " + ASETools.ElapsedTimeInSeconds(commonData.timer));
        } // Main

        static void HandleOneCase(ASETools.Case case_, int state)
        {
            var allcountReader = new ASETools.AllcountReader(case_.tumor_rna_allcount_filename);

            if (!allcountReader.openFile())
            {
                Console.WriteLine("Unable to open allcount file " + case_.tumor_rna_allcount_filename);
                return;
            }

            var valuesForThisCase = new Dictionary<string, Dictionary<string, int>>();  // gene->isoform->totalReadCount
            foreach (var geneInfo in commonData.geneLocationInformation.genesByName.Select(x => x.Value).ToList())
            {
                valuesForThisCase.Add(geneInfo.hugoSymbol, new Dictionary<string, int>());
                foreach (var isoformName in geneInfo.isoforms.Select(x => x.ucscId))
                {
                    valuesForThisCase[geneInfo.hugoSymbol].Add(isoformName, 0);
                } // isoform
            } // gene

            allcountReader.ReadAllcountFile((contigName, location, currentMappedReadCount) => processBase(valuesForThisCase, contigName, location, currentMappedReadCount));

            var mafLines = ASETools.MAFLine.ReadFile(case_.all_maf_lines_filename, "", false);
            bool anySpliceosomeMutations = mafLines.Any(_ => !_.IsSilent() && ASETools.spliceosome_genes.Contains(_.Hugo_Symbol));

            lock (isoformBalanceBySpliceosome)
            {
                foreach (var hugo_symbol in valuesForThisCase.Select(_ => _.Key).ToList())
                {
                    long totalIsoformReads = 0;  // This is more than the actual mapped reads, since some exons are in more than one isoform.  Hence "isoformReads"
                    valuesForThisCase[hugo_symbol].ToList().ForEach(_ => totalIsoformReads += _.Value); // This is instead of the more obvious .Sum() because it will overflow int
                        
                    foreach (var isoformName in valuesForThisCase[hugo_symbol].Select(_ => _.Key))
                    {
                        isoformBalanceBySpliceosome[hugo_symbol][isoformName][anySpliceosomeMutations].Add((double)valuesForThisCase[hugo_symbol][isoformName] / totalIsoformReads);
                    }
                }
            }

        }   // HandleOneCase

        static void processBase(Dictionary<string, Dictionary<string, int>> valuesForThisCase, string contigName, int location, int currentMappedReadCount)
        {
            foreach (var isoform in isoformMap.getIsoformsMappedTo(contigName, location))
            {
                valuesForThisCase[isoform.hugo_symbol][isoform.ucscId] += currentMappedReadCount;
            } // foreach isoform
        } // processBase
    }
}
