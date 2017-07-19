using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using System.IO;
using ASELib;

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
            public double bestZeroVsOne = 1;
            public int bestZeroVsOneAt = -1;
            public bool bestZeroVsOneIsExclusive = false;
            public double bestLocalZeroVsOne = 1;
            public int bestLocalZeroVsOneAt = -1;
            public bool bestLocalZeroVsOneIsExclusive = false;
            public bool significant = false;
            public double zeroVsOne = 0;

            public int minPAt = -1;
            public bool minPIsExclusive = false;
            public bool minPIsOneVsNotOne = false;
        }

        //
        // Apply the bonferroni correction and update our stats.
        //
        static void ProcessSinglePValue(ref double p, int bonferroniCorrection, bool exclusive, bool oneVsNotOne, int regionIndex, PValueStats pValueStats, ASETools.Histogram histogramOfPValues)
        {
            if (p == double.NegativeInfinity)   // Not a p-value because not enough data
            {
                return;
            }

            histogramOfPValues.addValue(p);

            p *= bonferroniCorrection;

            if (p < pValueStats.minP)
            {
                pValueStats.minP = p;
                pValueStats.minPAt = regionIndex;
                pValueStats.minPIsExclusive = exclusive;
                pValueStats.minPIsOneVsNotOne = oneVsNotOne;

                pValueStats.significant = p < configuration.significanceLevel;  // Recall, we're at the lowest p we've seen so far, so we can never set this to false if it's already true.
            }
        }

        static void ProcessSingleResult(ASETools.SingleExpressionResult singleResult, int bonferroniCorrection, bool exclusive, int regionIndex, PValueStats pValueStats, ASETools.Histogram histogramOfPValues)
        {
            ProcessSinglePValue(ref singleResult.oneVsMany, bonferroniCorrection, exclusive, false, regionIndex, pValueStats, histogramOfPValues);
            ProcessSinglePValue(ref singleResult.oneVsNotOne, bonferroniCorrection, exclusive, true, regionIndex, pValueStats, histogramOfPValues);
        }


        static void ProcessFile(string filename, StreamWriter allSignificantResultsFile, int bonferroniCorrection, out int overallValidCount, ASETools.Histogram histogramOfPValues)
        {
            overallValidCount = 0;
            if (!filename.EndsWith(".txt"))
            {
                Console.WriteLine("ProcessFile: filename " + filename + " does not end in .txt");
                return;
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
            outputFile.WriteLine("Hugo Symbol\tMin p\tMin p at\t0 vs. 1 ratio at MinP\tSignificant@.01\tBest 0 vs. 1 ratio for significant results\tBest 0 vs. 1 ratio at\t"+ ASETools.ExpressionResultsLine.getHeaderString());

            int nSignificantReuslts = 0;
            foreach (var result in results)
            {
                var pValueStats = new PValueStats();

                for (int i = 0; i < ASETools.nRegions; i++)
                {
                    ProcessSingleResult(result.nonExclusiveResultsByRange[i], bonferroniCorrection, false, i, pValueStats, histogramOfPValues);
                    ProcessSingleResult(result.exclusiveResultsByRange[i], bonferroniCorrection, true, i, pValueStats, histogramOfPValues);
                }

                outputFile.Write(ASETools.ConvertToExcelString(result.hugo_symbol) + "\t");

                if (pValueStats.minPAt == -1)
                {
                    outputFile.Write("*\t*\t*\tfalse\t*\t*\t")
                } else
                {
                    outputFile.Write(pValueStats.minP + "\t" + ASETools.regionIndexToString(pValueStats.minPAt) + (pValueStats.minPIsExclusive ? " exclusive\t" : "\t");
                    if (result.resultsByRange[pValueStats.minPIsExclusive][pValueStats.minpAt].)
                }

                + pValueStats.minP + "\t" + ASETools.regionIndexToString(pValueStats.minPAt) + (pValueStats.minPIsExclusive ? " exclusive" : "") + "\t");

                if (result.result.resultsByRange[pValueStats.minPIsExclusive][])

                string restOfOutputLine = "";   // Since we want to put minP at the beginning but don't know its value until the end, build up the output in this string.

                for (int whichField = nFields - 4; whichField < nFields; whichField++)
                {
                    restOfOutputLine += "\t" + fields[whichField];
                }
         


                for (int whichField = 1; whichField < nFields - 4; whichField++)    // First field is gene name, last four are nTumors, nZero, nOne, nMore
                {
                    if (fields[whichField] == "*" || (!fieldsToConvert[whichField] && !zeroValueFields[whichField] && !oneValueFields[whichField]) )
                    {
                        restOfOutputLine += "\t" + fields[whichField];
                    }
                    else
                    {
                        try
                        {
                            double value = Convert.ToDouble(fields[whichField]);
                            if (fieldsToConvert[whichField])
                            {
                                histogramOfPValues.addValue(value); // NB: Before Bonferroni correction.

                                value *= bonferroniCorrection;
                            }

                            if (fieldsToConvert[whichField] && value < minP)
                            {
                                minP = value;
                                minPAt = whichField;
                                justSetMinP = true;
                                zeroValue = -1;
                                oneValue = -1;
                                zeroVsOne = -1;
                                foundAny = true;
                            }

                            if (fieldsToConvert[whichField] && value <= .01)
                            {
                                significant = true;
                                zeroValue = -1;
                                oneValue = -1;


                            }

                            if (zeroValueFields[whichField])
                            {
                                zeroValue = value;
                            }

                            if (oneValueFields[whichField])
                            {
                                oneValue = value;
                            }

                            if (-1 != zeroValue && -1 != oneValue)
                            {
                                if (justSetMinP)
                                {
                                    if (0 == oneValue)
                                    {
                                        zeroVsOne = 100000; // Something big.
                                    }
                                    else
                                    {
                                        zeroVsOne = zeroValue / oneValue;
                                    }
                                }

                                if (significant && 0 != oneValue)
                                {
                                    double candidate = zeroValue / oneValue;
                                    if (candidate > 1)
                                    {
                                        candidate = 1 / candidate;
                                    }

                                    if (candidate < bestZeroVsOne)
                                    {
                                        bestZeroVsOne = candidate;
                                        bestZeroVsOneAt = whichField;
                                    }

                                    if (zeroValueLocalFields[whichField] || oneValueLocalFields[whichField] && candidate < bestLocalZeroVsOne)
                                    {
                                        bestLocalZeroVsOne = candidate;
                                        bestLocalZeroVsOneAt = whichField;
                                    }

                                    allSignificantResultsFile.WriteLine(value + "\t" + fields[0] + "\t" + ASETools.GetFileNameFromPathname(filename, true) + "\t" + headerFields[whichField] + "\t" + candidate);
                                    nSignificantReuslts++;
                                }
                                else if (significant)
                                {
                                    allSignificantResultsFile.WriteLine(value + "\t" + fields[0] + "\t" + ASETools.GetFileNameFromPathname(filename, true) + "\t" + headerFields[whichField]);
                                    nSignificantReuslts++;
                                }

                                justSetMinP = false;
                                zeroValue = -1;
                                oneValue = -1;
                                significant = false;
                            }

                            restOfOutputLine += "\t" + Convert.ToString(value);
                        } catch (FormatException) {
                            Console.WriteLine("Unparsable field " + fields[whichField] + " in line " + whichLine + " of file " + filename);
                            outputFile.Close();
                            File.Delete(outputFilename);
                            return;
                        }
                    }
                } // for each value field.

                if (foundAny)
                {
                    outputFile.Write("\t" + minP + "\t" + headerFields[minPAt] + "\t" + zeroVsOne + "\t" + ((minP < .01) ? "yes" : "no"));
                }
                else
                {
                    outputFile.Write("\t*\t*\t*\t*");
                }

                if (bestZeroVsOneAt != -1)
                {
                    outputFile.Write("\t" + bestZeroVsOne + "\t" + headerFields[bestZeroVsOneAt]);
                }
                else
                {
                    outputFile.Write("\t*\t*");
                }

                if (bestLocalZeroVsOneAt != -1)
                {
                    outputFile.Write("\t" + bestLocalZeroVsOne + "\t" + headerFields[bestLocalZeroVsOneAt]);
                }
                else
                {
                    outputFile.Write("\t*\t*");
                }

                outputFile.WriteLine(restOfOutputLine);
            }


            Console.WriteLine(ASETools.GetFileNameFromPathname(filename, true) + " has " + nSignificantReuslts + " signficant results");
            overallValidCount = nSignificantReuslts;
            outputFile.Close();
        }

        static void ProcessFileGroup(List<string> inputFiles, StreamWriter allSignificantResultsFile, out ASETools.Histogram histogramOfPValues)
        {
            int bonferonniCorrection = 0;
            foreach (var inputFile in inputFiles)
            {
                if (!inputFile.ToLower().Contains("_by_range_graphs"))
                {
                    int overallValidCount;
                    ProcessFile(inputFile, allSignificantResultsFile, false, 0, out overallValidCount, null);

                    bonferonniCorrection += overallValidCount;
                }
            }

            Console.WriteLine();
            Console.WriteLine("Total set Bonferonni correction is " + bonferonniCorrection);
            Console.WriteLine();

            int totalSignificantResults = 0;
            histogramOfPValues = new ASETools.Histogram();
            foreach (var inputFile in inputFiles)
            {
                if (!inputFile.ToLower().Contains("_by_range_graphs"))
                {
                    int overallValidCount;
                    ProcessFile(inputFile, allSignificantResultsFile, true, bonferonniCorrection, out overallValidCount, histogramOfPValues);

                    totalSignificantResults += overallValidCount;
                }
            }

            Console.WriteLine("Set total of " + totalSignificantResults + " signficant results");
            Console.WriteLine();
        }

        static void Main(string[] args)
        {
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

            inputFileSets.Add(Directory.GetFiles(configuration.finalResultsDirectory, "ExpressionDistributionByMutationCount*.txt").Where(x => !x.EndsWith(ASETools.bonferroniExtension)).ToList());
            inputFileSets.Add(Directory.GetFiles(configuration.finalResultsDirectory, "AlleleSpecificExpressionDistributionByMutationCount*.txt").Where(x => !x.EndsWith(ASETools.bonferroniExtension)).ToList());

            var allSignificantResultsFile = ASETools.CreateStreamWriterWithRetry(configuration.finalResultsDirectory + "AllSignificantResults.txt");
            allSignificantResultsFile.WriteLine("Effect size\tGene\tinput file\tp");

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

            foreach (var histogramOfPValues in histogramsOfPValues)
            {
                allSignificantResultsFile.WriteLine();
                allSignificantResultsFile.WriteLine("Min\tCount\tpdf\tcdf");
                var histogramResults = histogramOfPValues.ComputeHistogram(0, 1, .01, "N2");
                for (int i = 0; i < histogramResults.Length; i++)
                {
                    allSignificantResultsFile.WriteLine(histogramResults[i].minValue + "\t" + histogramResults[i].count + "\t" + histogramResults[i].pdfValue + "\t" + histogramResults[i].cdfValue);
                }
            }

            allSignificantResultsFile.Close();
        }
    }
}
