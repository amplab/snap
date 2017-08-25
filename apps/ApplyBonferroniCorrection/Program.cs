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
                    if (result.exclusiveResultsByRange[i].oneVsMany != double.NegativeInfinity) overallValidCount++;
                    if (result.exclusiveResultsByRange[i].oneVsNotOne != double.NegativeInfinity) overallValidCount++;
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

            public int minPAt = -1;
            public bool minPIsExclusive = false;
            public bool minPIsOneVsNotOne = false;
        }

        //
        // Apply the bonferroni correction and update our stats.
        //
        static int ProcessSinglePValue(string hugoSymbol, double ASE, double ASENotOne, string inputFilename, ref double p, int bonferroniCorrection, bool exclusive, bool oneVsNotOne, int regionIndex, PValueStats pValueStats, ASETools.Histogram histogramOfPValues, StreamWriter allSignificantResultsFile)
        {
            if (p == double.NegativeInfinity)   // Not a p-value because not enough data
            {
                return 0;
            }

            histogramOfPValues.addValue(p);

            p *= bonferroniCorrection;

            if (p < pValueStats.minP)
            {
                pValueStats.minP = p;
                pValueStats.minPAt = regionIndex;
                pValueStats.minPIsExclusive = exclusive;
                pValueStats.minPIsOneVsNotOne = oneVsNotOne;

                pValueStats.significant = p <= configuration.significanceLevel;  // Recall, we're at the lowest p we've seen so far, so we can never set this to false if it's already true.
            }

            if (p <= configuration.significanceLevel)
            {
                allSignificantResultsFile.WriteLine(hugoSymbol + "\t" + ASE + "\t" + (ASENotOne == double.NegativeInfinity ? "*" : Convert.ToString(ASENotOne)) + "\t" + 
                    ASETools.GetFileNameFromPathname(inputFilename) + "\t" + exclusive + "\t" + !oneVsNotOne + "\t" + regionIndex + "\t" + ASETools.regionIndexToString(regionIndex) + "\t" + p);
                return 1;
            } else
            {
                return 0;
            }            
        }

        static int ProcessSingleResult(string hugoSymbol, string inputFilename, ASETools.SingleExpressionResult singleResult, int bonferroniCorrection, bool exclusive, int regionIndex, PValueStats pValueStats, ASETools.Histogram histogramOfPValues, StreamWriter allSignificantResultsFile)
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

            nSignificantResults += ProcessSinglePValue(hugoSymbol, singleResult.oneMutationStats.mean, ASENotOne, inputFilename, ref singleResult.oneVsMany, bonferroniCorrection, exclusive, false, regionIndex, pValueStats, histogramOfPValues, allSignificantResultsFile);
            nSignificantResults += ProcessSinglePValue(hugoSymbol, singleResult.oneMutationStats.mean, ASENotOne, inputFilename, ref singleResult.oneVsNotOne, bonferroniCorrection, exclusive, true, regionIndex, pValueStats, histogramOfPValues, allSignificantResultsFile);

            if (singleResult.oneMutationStats.n > 0 && singleResult.zeroMutationStats.n > 0 && singleResult.oneMutationStats.mean != 0 && (singleResult.oneVsNotOne <= configuration.significanceLevel || singleResult.oneVsNotOne <= configuration.significanceLevel))
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


        static int ProcessFile(string filename, StreamWriter allSignificantResultsFile, int bonferroniCorrection, ASETools.Histogram histogramOfPValues)
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
			var outputFile = ASETools.CreateStreamWriterWithRetry(outputFilename);
            outputFile.WriteLine("Hugo Symbol\tMin p\tMin p at\tSignificant@.01\tBest 0 vs. 1 ratio for significant results\tBest 0 vs. 1 ratio at\t"+ ASETools.ExpressionResultsLine.getHeaderString());

            int nSignificantReuslts = 0;
            foreach (var result in results)
            {
                var pValueStats = new PValueStats();
				
                for (int i = 0; i < ASETools.nRegions; i++)
                {
                    nSignificantReuslts += ProcessSingleResult(result.hugo_symbol, filename, result.nonExclusiveResultsByRange[i], bonferroniCorrection, false, i, pValueStats, histogramOfPValues, allSignificantResultsFile);
                    nSignificantReuslts += ProcessSingleResult(result.hugo_symbol, filename, result.exclusiveResultsByRange[i], bonferroniCorrection, true, i, pValueStats, histogramOfPValues, allSignificantResultsFile);
                }

                outputFile.Write(ASETools.ConvertToExcelString(result.hugo_symbol) + "\t");

                if (pValueStats.minPAt == -1)
                {
                    outputFile.Write("*\t*\tfalse\t*\t*\t");
                } else
                {
                    outputFile.Write(pValueStats.minP + "\t" + ASETools.regionIndexToString(pValueStats.minPAt) + (pValueStats.minPIsExclusive ? " exclusive\t" : "\t") + (pValueStats.significant ? "true\t" : "false\t"));
                    if (pValueStats.bestZeroVsOneAt != -1)
                    {
                        outputFile.Write(pValueStats.bestZeroVsOne + "\t" + ASETools.regionIndexToString(pValueStats.bestZeroVsOneAt) + (pValueStats.bestZeroVsOneIsExclusive ? " exclusive\t" : "\t"));
                    } else
                    {
                        outputFile.Write("*\t*\t");
                    }

                }

                outputFile.Write(result.hugo_symbol);    // Yes, this is here twice.  We want it as the first column, and then having it again lets us use ASETools.ExpressionResultsLine.getHeaderString(), which is convenient

                foreach (bool exclusive in ASETools.BothBools)
                {
                     for (int region = 0; region < ASETools.nRegions; region++)
                    {
                        outputFile.Write("\t");
                        result.resultsByRange[exclusive][region].writeToFile(outputFile);
                    }
                }

                outputFile.WriteLine("\t" + result.nTumorsExcluded + "\t" + result.nZero + "\t" + result.nOne + "\t" + result.nMore);
            }

            outputFile.Close();
            Console.WriteLine("input file " + filename + " has " + nSignificantReuslts + " significant results.");

            return nSignificantReuslts;
        }

        static void ProcessFileGroup(List<string> inputFilenames, StreamWriter allSignificantResultsFile, out ASETools.Histogram histogramOfPValues)
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
            histogramOfPValues = new ASETools.Histogram(inputFilenames[0]);
            foreach (var inputFilename in inputFilenames)
            {
                totalSignificantResults += ProcessFile(inputFilename, allSignificantResultsFile, bonferonniCorrection, histogramOfPValues);
            }

            Console.WriteLine("Set total of " + totalSignificantResults + " signficant results");
            Console.WriteLine();
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

            if (configuration.commandLineArgs.Count() != 0)
            {
                Console.WriteLine("usage: ApplyBonferonniCorrection {-configuration configuration}");
                return;
            }

            //
            // Get the list of input files.  Note that these lists will also include the output files, so we filter them with a .Where().
            //
            List<List<string>> inputFileSets = new List<List<string>>();

            inputFileSets.Add(Directory.GetFiles(configuration.finalResultsDirectory, "ExpressionDistributionByMutationCount*.txt").Where(x => !x.Contains(ASETools.bonferroniExtension)).ToList());
            inputFileSets.Add(Directory.GetFiles(configuration.finalResultsDirectory, "AlleleSpecificExpressionDistributionByMutationCount*.txt").Where(x => !x.Contains(ASETools.bonferroniExtension)).ToList());

            var allSignificantResultsFile = ASETools.CreateStreamWriterWithRetry(configuration.finalResultsDirectory + ASETools.AllSignificantResultsFilename);
            allSignificantResultsFile.WriteLine("Hugo Symbol\tASE (one mutation)\tASE (not one mutation)\tinput file\texclusive\tOneVsMany\trange index\trange\tp");

            var histogramsOfPValues = new List<ASETools.Histogram>();

            foreach (var inputFileSet in inputFileSets)
            {
                if (inputFileSet.Count() == 0)
                {
                    continue;
                }
                ASETools.Histogram histogramOfPValues;
                ProcessFileGroup(inputFileSet, allSignificantResultsFile, out histogramOfPValues);
                histogramsOfPValues.Add(histogramOfPValues);
            }

            allSignificantResultsFile.WriteLine("**done**");
            allSignificantResultsFile.Close();

            var pValueHistogramFile = ASETools.CreateStreamWriterWithRetry(configuration.finalResultsDirectory + ASETools.pValueHistogramFilename);

            foreach (var histogramOfPValues in histogramsOfPValues)
            {
                pValueHistogramFile.WriteLine(histogramOfPValues.name);
                pValueHistogramFile.WriteLine("Min\tCount\tpdf\tcdf");
                var histogramResults = histogramOfPValues.ComputeHistogram(0, 1, .01, "N2");
                for (int i = 0; i < histogramResults.Length; i++)
                {
                    pValueHistogramFile.WriteLine(histogramResults[i].minValue + "\t" + histogramResults[i].count + "\t" + histogramResults[i].pdfValue + "\t" + histogramResults[i].cdfValue);
                }
            }

            pValueHistogramFile.Close();
 
            Console.WriteLine("ApplyBonferroniCorrection took " + ASETools.ElapsedTimeInSeconds(timer));
        }
    }
}
