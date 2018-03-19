using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using System.IO;
using ASELib;
using System.Diagnostics;

namespace ApplyBonferroniCorrection
{
    class Program
    {
        static ASETools.Configuration configuration;
        static ASETools.GeneLocationsByNameAndChromosome geneLocationInformation;


        static int CountValidPValuesInFile(string filename)
        {
            int overallValidCount = 0;

			var results = ASETools.ExpressionResultsLine.readFile(filename);

			if (null == results)
            {
                Console.WriteLine("Unable to load results from file " + filename);
                return 0;
            }

                foreach (var result in results)
            {
                for (int i = 0; i < ASETools.nRegions; i++)
                {
                    if (result.nonExclusiveResultsByRange[i].oneVsMany != double.NegativeInfinity) overallValidCount++;
                    if (result.nonExclusiveResultsByRange[i].oneVsNotOne != double.NegativeInfinity) overallValidCount++;
                    if (result.nonExclusiveResultsByRange[i].oneSilentVsOneNotSilent != double.NegativeInfinity) overallValidCount++;
                    if (result.exclusiveResultsByRange[i].oneVsMany != double.NegativeInfinity) overallValidCount++;
                    if (result.exclusiveResultsByRange[i].oneVsNotOne != double.NegativeInfinity) overallValidCount++;
                    if (result.exclusiveResultsByRange[i].oneSilentVsOneNotSilent != double.NegativeInfinity) overallValidCount++;
                }
            }

            Console.WriteLine(ASETools.GetFileNameFromPathname(filename, true) + " has overall valid count (i.e. Bonferroni correction contribution) of " + overallValidCount);
            //
            // We're just computing overallValidCount, which we have done.
            //
            return overallValidCount;
        } // CountValidPValuesInFile
        
        class PValueStats
        {
            public bool foundAny = false;
            public double minP = 10000000;
            public double bestZeroVsOne = double.MinValue;
            public int bestZeroVsOneAt = -1;
            public bool bestZeroVsOneIsExclusive = false;
            public double bestLocalZeroVsOne = 1;
            public int bestLocalZeroVsOneAt = -1;
            public bool bestLocalZeroVsOneIsExclusive = false;
            public bool significant = false;
            public List<string> significantAt = new List<string>();

            public int minPAt = -1;
            public bool minPIsExclusive = false;
            public bool minPIsOneVsNotOne = false;
        }

        public static string histogramSetToName(int whichHistogramSet)
        {
            switch(whichHistogramSet)
            {
                case 0: return "1 vs. not 1";
                case 1: return "1 vs. many";
                case 2: return "1 silent vs. 1 not silent";
                default: throw new IndexOutOfRangeException();
            }
        }

        //
        // Apply the bonferroni correction and update our stats.
        //
        static int ProcessSinglePValue(string hugoSymbol, int nOne, int nOther, double ASE, double ASENotOne, string inputFilename, ref double p, int bonferroniCorrection, 
            bool exclusive, int whichHistogramSet /* 0 - one v. not 1, 1 - 1 v many, 2 - one silent v. one not silent*/, int regionIndex, PValueStats pValueStats, HistogramSet histogramSet, StreamWriter allSignificantResultsFile)
        {
            if (p == double.NegativeInfinity)   // Not a p-value because not enough data
            {
                return 0;
            }

            histogramSet.allResults.addValue(p);

            if (whichHistogramSet == 0)
            { 
                histogramSet.oneVsNotOne.addValue(p);
            } else if (whichHistogramSet == 1)
            {
                histogramSet.oneVsMany.addValue(p);
            } else
            {
                histogramSet.oneSilentVsOneNotSilent.addValue(p);
            }

            if (ASETools.GetFileNameFromPathname(inputFilename).Contains("_"))
            {
                histogramSet.singleDisease.addValue(p);
            } else
            {
                histogramSet.panCancer.addValue(p);
            }

            if (nOne > 20 && nOther > 20)
            {
                histogramSet.nGreaterThan20.addValue(p);
            } else
            {
                histogramSet.nLessThanOrEqualTo20.addValue(p);
            }

            if (exclusive)
            {
                histogramSet.exclusive.addValue(p);
            } else
            {
                histogramSet.inclusive.addValue(p);
            }

            if (regionIndex > 11)
            {
                histogramSet.atLeastAMegabase.addValue(p);
            } else
            {
                histogramSet.lessThanAMegabase.addValue(p);
            }

            p *= bonferroniCorrection;

            if (p < pValueStats.minP && whichHistogramSet < 2)
            {
                pValueStats.minP = p;
                pValueStats.minPAt = regionIndex;
                pValueStats.minPIsExclusive = exclusive;
                pValueStats.minPIsOneVsNotOne = whichHistogramSet == 0;

                pValueStats.significant = p <= configuration.significanceLevel;  // Recall, we're at the lowest p we've seen so far, so we can never set this to false if it's already true.
            }

            if (p <= configuration.significanceLevel)
            {
                pValueStats.significantAt.Add(Convert.ToString(regionIndex) + (exclusive ? "E" : "") + ((whichHistogramSet == 2) ? "S" : ""));

if (ASE == double.NegativeInfinity)
{
    Console.WriteLine("Here!");
}

                allSignificantResultsFile.WriteLine(hugoSymbol + "\t" + ASE + "\t" + (ASENotOne == double.NegativeInfinity ? "*" : Convert.ToString(ASENotOne)) + "\t" + 
                    ASETools.GetFileNameFromPathname(inputFilename) + "\t" + exclusive + "\t" + histogramSetToName(whichHistogramSet) + "\t" + regionIndex + "\t" + ASETools.regionIndexToString(regionIndex) + "\t" + p);
                return 1;
            } else
            {
                return 0;
            }            
        }

        static int ProcessSingleResult(string hugoSymbol, string inputFilename, ASETools.SingleExpressionResult singleResult, int bonferroniCorrection, bool exclusive, int regionIndex, 
            PValueStats pValueStats, HistogramSet histogramSet, StreamWriter allSignificantResultsFile)
        {
            int nSignificantResults = 0;

            double ASENotOne;
            if (singleResult.zeroMutationStats.n + singleResult.moreThanOneMutationStats.n == 0)
            {
                ASENotOne = double.NegativeInfinity;
            } else
            {
                ASENotOne = (singleResult.zeroMutationStats.n * singleResult.zeroMutationStats.mean + singleResult.moreThanOneMutationStats.n * singleResult.moreThanOneMutationStats.mean) / (singleResult.zeroMutationStats.n + singleResult.moreThanOneMutationStats.n);
            }

            nSignificantResults += ProcessSinglePValue(hugoSymbol, singleResult.oneMutationStats.n, singleResult.moreThanOneMutationStats.n, singleResult.oneMutationStats.mean, singleResult.moreThanOneMutationStats.mean, inputFilename, ref singleResult.oneVsMany, bonferroniCorrection, exclusive, 1, regionIndex, pValueStats, histogramSet, allSignificantResultsFile);
            nSignificantResults += ProcessSinglePValue(hugoSymbol, singleResult.oneMutationStats.n, singleResult.moreThanOneMutationStats.n + singleResult.zeroMutationStats.n, singleResult.oneMutationStats.mean, ASENotOne, inputFilename, ref singleResult.oneVsNotOne, bonferroniCorrection, exclusive, 0, regionIndex, pValueStats, histogramSet, allSignificantResultsFile);
            nSignificantResults += ProcessSinglePValue(hugoSymbol, singleResult.onlyOneSilentMutationStats.n, singleResult.oneMutationStats.n, singleResult.onlyOneSilentMutationStats.mean, singleResult.oneMutationStats.mean, inputFilename, ref singleResult.oneSilentVsOneNotSilent, bonferroniCorrection, exclusive, 2, regionIndex, pValueStats, histogramSet, allSignificantResultsFile);

            if (singleResult.oneMutationStats.n > 0 && singleResult.zeroMutationStats.n > 0 && singleResult.oneMutationStats.mean != 0 && 
                (singleResult.oneVsNotOne <= configuration.significanceLevel && singleResult.oneVsNotOne != double.NegativeInfinity|| singleResult.oneVsMany <= configuration.significanceLevel && singleResult.oneVsMany != double.NegativeInfinity))
            {
                //
                // We have a significant zero-vs-one.
                //

                double zeroVsOne = singleResult.oneMutationStats.mean / singleResult.zeroMutationStats.mean;
                if (zeroVsOne > pValueStats.bestZeroVsOne)
                {
                    pValueStats.bestZeroVsOne = zeroVsOne;
                    pValueStats.bestZeroVsOneAt = regionIndex;
                    pValueStats.bestZeroVsOneIsExclusive = exclusive;
                }
            }

            return nSignificantResults;
        }

        class SortKeyAndLine : IComparable
        {
            public SortKeyAndLine(bool significant_, double amount_, string hugo_symbol_, string outputLine_)
            {
                significant = significant_;
                amount = amount_;
                hugo_symbol = hugo_symbol_;
                outputLine = outputLine_;
            }
            public readonly bool significant;
            public readonly double amount;
            public readonly string hugo_symbol;
            public readonly string outputLine;

            public int CompareTo(object peerObject)
            {
                SortKeyAndLine peer = (SortKeyAndLine)peerObject;
                // Sort order is significant before not, then backward by amount (so larger amounts sort first).

                if (significant && !peer.significant) return -1;
                if (peer.significant && !significant) return 1;

                if (significant)
                {
                    if (amount != peer.amount) return peer.amount.CompareTo(amount);  // I cleverely reversed the order in the CompareTo to get it to sort higher first.
                }

                return hugo_symbol.CompareTo(peer.hugo_symbol);
            }
        }

        static int ProcessFile(string filename, StreamWriter allSignificantResultsFile, int bonferroniCorrection, HistogramSet histogramSet)
        {
            if (!filename.EndsWith(".txt"))
            {
                Console.WriteLine("ProcessFile: filename " + filename + " does not end in .txt");
                return 0;
            }

            string outputFilename = filename.Substring(0, filename.Count() - 4) + ASETools.bonferroniExtension;

			var results = ASETools.ExpressionResultsLine.readFile(filename);

			//
			// We rearrange the order of the columns somewhat.  The input is geneName followed by sets of columns with uncorrected p values, raw value means or standard 
			// deviations or counts, followed by four columns that have a count of tumors by mutation count for this gene.  We then add in min p and min p at columns.
			// The final output format is gene name, min p, min p at, the four tumor count columns from the end of the input and then the rest of the columns from
			// the input in order.
			//

            var sortableLines = new List<SortKeyAndLine>();

            int nSignificantReuslts = 0;
            foreach (var result in results)
            {
                var pValueStats = new PValueStats();
				
                for (int i = 0; i < ASETools.nRegions; i++)
                {
                    nSignificantReuslts += ProcessSingleResult(result.hugo_symbol, filename, result.nonExclusiveResultsByRange[i], bonferroniCorrection, false, i, pValueStats, histogramSet, allSignificantResultsFile);
                    nSignificantReuslts += ProcessSingleResult(result.hugo_symbol, filename, result.exclusiveResultsByRange[i], bonferroniCorrection, true, i, pValueStats, histogramSet, allSignificantResultsFile);
                }

                // hugo symbol
                var outputLine = ASETools.ConvertToExcelString(result.hugo_symbol);

                // alt fraction single
                if (result.nSingleContibutingToAltFraction < 10)
                {
                    outputLine += "\t*";
                } else
                {
                    outputLine += "\t" + ASETools.stringForDouble(result.singleAltFraction);
                }

                // alt fraction multiple
                if (result.nMultipleContributingToAltFraction < 10)
                {
                    outputLine += "\t*\t";
                } else
                {
                    outputLine += "\t" + ASETools.stringForDouble(result.multipleAltFraction) + "\t";
                }
 
                // min p, min P at, signficant, best 0 vs. 1, best 0 vs 1 ratio
                if (pValueStats.minPAt == -1)
                {
                    outputLine += "*\t*\tfalse\t*\t*\t";
                } else
                {
                    outputLine +=  pValueStats.minP + "\t" + ASETools.regionIndexToString(pValueStats.minPAt) + (pValueStats.minPIsExclusive ? " exclusive\t" : "\t") + (pValueStats.significant ? "true\t" : "false\t");
                    if (pValueStats.bestZeroVsOneAt != -1)
                    {
                        outputLine += pValueStats.bestZeroVsOne + "\t" + ASETools.regionIndexToString(pValueStats.bestZeroVsOneAt) + (pValueStats.bestZeroVsOneIsExclusive ? " exclusive\t" : "\t");
                    } else
                    {
                        outputLine += "*\t*\t";
                    }
                }

                // significant at
                pValueStats.significantAt.ForEach(x => outputLine += x + ",");

                // gene size
                if (geneLocationInformation.genesByName.ContainsKey(result.hugo_symbol))
                {
                    outputLine += "\t" + geneLocationInformation.genesByName[result.hugo_symbol].size() + "\t" + geneLocationInformation.genesByName[result.hugo_symbol].codingSize();
                } else
                {
                    outputLine += "\t*\t*";
                }

                outputLine += "\t" + ASETools.ConvertToExcelString(result.hugo_symbol);    // Yes, this is here twice.  We want it as the first column, and then having it again lets us use ASETools.ExpressionResultsLine.getHeaderString(), which is convenient

                bool anyValid = false;

                foreach (bool exclusive in ASETools.BothBools)
                {
                    for (int region = 0; region < ASETools.nRegions; region++)
                    {
                        outputLine += "\t";
                        result.resultsByRange[exclusive][region].appendToString(ref outputLine);
                        anyValid |= result.resultsByRange[exclusive][region].anyValid();
                    }
                }

                if (!anyValid)
                {
                    continue;   // Just cull these from the output, they're useless.
                }

                foreach (bool inclusive in ASETools.BothBools)
                {
                    for (int zeroOneManySilent = 0; zeroOneManySilent < 4; zeroOneManySilent++)
                    {
                       for (int region = 0; region < ASETools.nRegions; region++)
                        {
                            outputLine += "\t" + result.resultsByRange[inclusive][region].nMutationMean(zeroOneManySilent);
                        } // region
                    } // zeroOneMany
                } // exclusive

                outputLine += "\t" + result.nTumorsExcluded + "\t" + result.nZero + "\t" + result.nOne + "\t" + result.nMore + "\t" + result.nSingleContibutingToAltFraction + "\t" + result.nMultipleContributingToAltFraction + "\t" + 
                    ASETools.stringForDouble(result.singleAltFraction) + "\t" + ASETools.stringForDouble(result.multipleAltFraction);

                sortableLines.Add(new SortKeyAndLine(pValueStats.significant, pValueStats.bestZeroVsOne, result.hugo_symbol, outputLine));
            }
            var outputFile = ASETools.CreateStreamWriterWithRetry(outputFilename);
            outputFile.WriteLine("Hugo_Symbol\tAlt fraction for single mutations\tAlt fraction for multiple mutations\tMin p\tMin p at\tSignificant@.01\tBest 0 vs. 1 ratio for significant results\tBest 0 vs. 1 ratio at\tSignificant At\tGene Size\tCoding Size\t" + ASETools.ExpressionResultsLine.getHeaderString());

            sortableLines.Sort();

            sortableLines.ForEach(x => outputFile.WriteLine(x.outputLine));

            outputFile.WriteLine("**done**");
            outputFile.Close();
            Console.WriteLine("input file " + filename + " has " + nSignificantReuslts + " significant result" + ((nSignificantReuslts == 1) ? "." : "s."));

            return nSignificantReuslts;
        }

        static void ProcessFileGroup(List<string> inputFilenames, StreamWriter allSignificantResultsFile, out HistogramSet histogramSet)
        {
            int bonferonniCorrection = 0;
            foreach (var inputFilename in inputFilenames)
            {
                if (!inputFilename.ToLower().Contains("_by_range_graphs"))
                {
                    bonferonniCorrection += CountValidPValuesInFile(inputFilename);
                }
            }

            Console.WriteLine();
            Console.WriteLine("Total set Bonferonni correction is " + bonferonniCorrection);
            Console.WriteLine();

            int totalSignificantResults = 0;
            histogramSet = new HistogramSet(inputFilenames[0]);
            foreach (var inputFilename in inputFilenames)
            {
                totalSignificantResults += ProcessFile(inputFilename, allSignificantResultsFile, bonferonniCorrection, histogramSet);
            }

            Console.WriteLine("Set total of " + totalSignificantResults + " signficant results");
            Console.WriteLine();
        }

        class HistogramSet
        {
            public HistogramSet(string namePrefix)
            {
                allResults = new ASETools.Histogram(namePrefix + " all results");
                oneVsMany = new ASETools.Histogram(namePrefix + " one vs. many");
                oneVsNotOne = new ASETools.Histogram(namePrefix + " one vs. not one");
                oneSilentVsOneNotSilent = new ASETools.Histogram(namePrefix + " one silent vs. none at all");
                panCancer = new ASETools.Histogram(namePrefix + " pan cancer");
                singleDisease = new ASETools.Histogram(namePrefix + " single disease");
                nGreaterThan20 = new ASETools.Histogram(namePrefix + " n > 20");
                nLessThanOrEqualTo20 = new ASETools.Histogram(namePrefix + " n <= 20");
                lessThanAMegabase = new ASETools.Histogram(namePrefix + " less than a megabase");
                atLeastAMegabase = new ASETools.Histogram(namePrefix + " at least a megabase");
                inclusive = new ASETools.Histogram(namePrefix + " inclusive");
                exclusive = new ASETools.Histogram(namePrefix + " exclusive");

                allHistograms.Add(allResults);
                allHistograms.Add(oneVsMany);
                allHistograms.Add(oneVsNotOne);
                allHistograms.Add(oneSilentVsOneNotSilent);
                allHistograms.Add(panCancer);
                allHistograms.Add(singleDisease);
                allHistograms.Add(nGreaterThan20);
                allHistograms.Add(nLessThanOrEqualTo20);
                allHistograms.Add(lessThanAMegabase);
                allHistograms.Add(atLeastAMegabase);
                allHistograms.Add(inclusive);
                allHistograms.Add(exclusive);
            }
            public ASETools.Histogram allResults;
            public ASETools.Histogram oneVsMany;
            public ASETools.Histogram oneVsNotOne;
            public ASETools.Histogram oneSilentVsOneNotSilent;
            public ASETools.Histogram panCancer;
            public ASETools.Histogram singleDisease;
            public ASETools.Histogram nGreaterThan20;
            public ASETools.Histogram nLessThanOrEqualTo20;
            public ASETools.Histogram lessThanAMegabase;
            public ASETools.Histogram atLeastAMegabase;
            public ASETools.Histogram inclusive;
            public ASETools.Histogram exclusive;
            public List<ASETools.Histogram> allHistograms = new List<ASETools.Histogram>();
        }

        static void Main(string[] args)
        {
            var timer = new Stopwatch();
            timer.Start();

            configuration = ASETools.Configuration.loadFromFile(args);
            if (null == configuration)
            {
                Console.WriteLine("Couldn't load configuration.");
                return;
            }

            geneLocationInformation = new ASETools.GeneLocationsByNameAndChromosome(ASETools.readKnownGeneFile(configuration.geneLocationInformationFilename));

            double minASE = -1;
            double maxASE = 2;
            double minMaxChromASE = -1;
            bool dependOnTP53 = false;
            bool excludeTP53Mutants = false;

            for (int i = 0; i < configuration.commandLineArgs.Count(); i++)
            {
                if (configuration.commandLineArgs[i] == "-MaxASE" && configuration.commandLineArgs.Count() > i + 1)
                {
                    maxASE = Convert.ToDouble(configuration.commandLineArgs[i + 1]);
                    i++;
                } else if (configuration.commandLineArgs[i] == "-MinASE" && configuration.commandLineArgs.Count() > i + 1)
                {
                    minASE = Convert.ToDouble(configuration.commandLineArgs[i + 1]);
                    i++;
                } else if (configuration.commandLineArgs[i] == "-MinMaxChromASE" && configuration.commandLineArgs.Count() > i + 1)
                {
                    minMaxChromASE = Convert.ToDouble(configuration.commandLineArgs[i + 1]);
                    i++;
                }
                else if (configuration.commandLineArgs[i] == "-TP53Mutant")
                {
                    dependOnTP53 = true;
                    excludeTP53Mutants = false;
                }
                else if (configuration.commandLineArgs[i] == "-NoTP53Mutant")
                {
                    dependOnTP53 = true;
                    excludeTP53Mutants = true;
                } else
                {
                    Console.WriteLine("usage: ApplyBonferonniCorrection {-configuration configuration} {-MaxASE maxASE} {-MinASE minASE} {-TP53Mutant|-NoTP53Mutant}");
                    return;
                }
            }

            string filenameExtra = "";

            if (maxASE <= 1)
            {
                filenameExtra = "-maxASE-" + maxASE;
            }

            if (minASE >= 0)
            {
                filenameExtra += "-minASE-" + minASE;
            }

            if (minMaxChromASE >= 0)
            {
                filenameExtra += "-minMaxChromASE-";
            }
            if (dependOnTP53)
            {
                if (excludeTP53Mutants)
                {
                    filenameExtra += "-NoTP53Mutants-";
                }
                else
                {
                    filenameExtra += "-OnlyTP53Mutants-";
                }
            }

            //
            // Get the list of input files.  Note that these lists will also include the output files, so we filter them with a .Where().
            //
            List<List<string>> inputFileSets = new List<List<string>>();

            //inputFileSets.Add(Directory.GetFiles(configuration.finalResultsDirectory, "ExpressionDistributionByMutationCount*.txt").Where(x => !x.Contains("bonferroni")).ToList());
            inputFileSets.Add(Directory.GetFiles(configuration.finalResultsDirectory, ASETools.AlleleSpecificExpressionDistributionByMutationCountFilenameBase + filenameExtra + "*.txt").
                Where(x => !x.Contains("bonferroni") && (maxASE <= 1 || !x.Contains("-maxASE-")) && (minASE >= 0 || !x.Contains("-minASE-")) && (minMaxChromASE >= 0 || !x.Contains("-minMaxChromASE-")) && (dependOnTP53 || !x.Contains("TP53"))).ToList());

            var allSignificantResultsFile = ASETools.CreateStreamWriterWithRetry(configuration.finalResultsDirectory + ASETools.AllSignificantResultsFilenameBase + filenameExtra + ".txt");
            allSignificantResultsFile.WriteLine("Hugo Symbol\tASE (one mutation)\tASE (not one mutation)\tinput file\texclusive\tOneVsMany\trange index\trange\tp");

            var histogramSets = new List<HistogramSet>();

            foreach (var inputFileSet in inputFileSets)
            {
                if (inputFileSet.Count() == 0)
                {
                    continue;
                }

                HistogramSet histogramSet;

                ProcessFileGroup(inputFileSet, allSignificantResultsFile, out histogramSet);

                histogramSets.Add(histogramSet);
            }

            allSignificantResultsFile.WriteLine("**done**");
            allSignificantResultsFile.Close();

            var pValueHistogramFile = ASETools.CreateStreamWriterWithRetry(configuration.finalResultsDirectory + ASETools.pValueHistogramFilename);

            foreach (var histogramSet in histogramSets)
            {
                foreach (var histogram in histogramSet.allHistograms)
                {

                    pValueHistogramFile.WriteLine(histogram.name);
                    pValueHistogramFile.WriteLine(ASETools.HistogramResultLine.Header());
                    histogram.ComputeHistogram(0, 1, .01, "N2").ToList().ForEach(x => pValueHistogramFile.WriteLine(x));
                }
            }

            pValueHistogramFile.Close();
 
            Console.WriteLine("ApplyBonferroniCorrection took " + ASETools.ElapsedTimeInSeconds(timer));
        }
    }
}
