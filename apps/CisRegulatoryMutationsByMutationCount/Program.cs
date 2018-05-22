using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using System.IO;
using System.Diagnostics;
using ASELib;

namespace CisRegulatoryMutationsByMutationCount
{
    class Program
    {
        static ASETools.Configuration configuration;
        static StreamWriter outputFile;

        class PerThreadState
        {
            public Dictionary<string, List<double>[,]> measurements = new Dictionary<string, List<double>[,]>();  // Hugo symbol -> array of [0,1,2] for low, middle, high count -> array by region index -> measurements
        }
        
        static Dictionary<string, List<double>[,]>  globalMeasurements = new Dictionary<string, List<double>[,]>();

        class DoubleAndBool : IComparer<DoubleAndBool>
        {
            public readonly double value;
            public readonly bool isHigh;

            public DoubleAndBool(double value_, bool isOne_)
            {
                value = value_;
                isHigh = isOne_;
            }

            public int Compare(DoubleAndBool a, DoubleAndBool b)
            {
                return a.value.CompareTo(b.value);
            }
        }

        static int countList(List<double> list)
        {
            if (list == null) return 0;
            return list.Count();
        }

        static string valueOrStar(double value)
        {
            if (value == double.NegativeInfinity)
            {
                return "*";
            }
            return value.ToString();
        }

        static string meanOrStar(List<double> values)
        {
            if (countList(values) == 0)
            {
                return "*";
            }

            return values.Average().ToString();
        }

        static string ratioOrStar(List<double> numerator, List<double> denominator)
        {
            if (denominator.Count() == 0 || numerator.Count() == 0)
            {
                return "*";
            }

            double denominatorValue = denominator.Average();
            if (denominatorValue == 0)
            {
                return "*";
            }

            return (numerator.Average() / denominatorValue).ToString();
        }

        const int nValuesNeeded = 10;
        static int totalValidPValues = 0;
        static double smallestPValue = 2;

        static int nRegionsToUse = 13;  // Only go this far, don't look absurdly far away.

        static void Main(string[] args)
        {

            var timer = new Stopwatch();
            timer.Start();

            configuration = ASETools.Configuration.loadFromFile(args);
            if (null == configuration)
            {
                Console.WriteLine("Unable to load configuration.");
                return;
            }

            var cases = ASETools.Case.LoadCases(configuration.casesFilePathname);
            if (null == cases)
            {
                Console.WriteLine("Unable to load cases.");
                return;
            }

            var casesToProcess = cases.Select(x => x.Value).Where(x => x.regulatory_mutations_near_mutations_filename != "").ToList();

            var threading = new ASETools.WorkerThreadHelper<ASETools.Case, PerThreadState>(casesToProcess, HandleOneCase, FinishUp, null, 100);
            Console.WriteLine("Processing " + casesToProcess.Count() + " cases, 1 dot/100 cases:");
            ASETools.PrintNumberBar(casesToProcess.Count() / 100);

            threading.run();

            outputFile = ASETools.CreateStreamWriterWithRetry(configuration.finalResultsDirectory + ASETools.cisRegulatoryMutationsByVAFFilename);
            if (null == outputFile)
            {
                Console.WriteLine("Unable to open output file " + configuration.finalResultsDirectory + ASETools.cisRegulatoryMutationsByVAFFilename);
                return;
            }

            outputFile.Write("Hugo symbol\tmin P");
            for (int region = 1; region < nRegionsToUse; region++)
            {
                var regionName = ASETools.regionIndexToString(region);
                outputFile.Write("\t" + regionName + " nLow\t" + regionName + " nModerate\t" + regionName + " nHigh\t" + regionName + " P high vs. not high\t" + regionName + " P high vs. moderate\t" + regionName + " low mean\t" +
                    regionName + " moderate mean\t" + regionName + " high mean\t" + regionName + " high mean/moderate mean");
            }
            outputFile.WriteLine();

            var genesToProcess = globalMeasurements.ToList();
            var geneThreading = new ASETools.WorkerThreadHelper<KeyValuePair<string, List<double>[,]>, int>(genesToProcess, HandleOneGene, null, null, 100);
            Console.WriteLine("Processing " + genesToProcess.Count()+  " genes, one dot/100 genes:");
            ASETools.PrintNumberBar(genesToProcess.Count() / 100);
            geneThreading.run();

            outputFile.WriteLine("**done**");
            outputFile.Close();

            Console.WriteLine(ASETools.ElapsedTimeInSeconds(timer));
            Console.WriteLine("Total of " + totalValidPValues + " P values, smallest " + smallestPValue);
        } // Main

        static void HandleOneCase(ASETools.Case case_, PerThreadState state)
        {
            var regulatoryMutationsByVAF = ASETools.RegulatoryMutationsNearGene.readFile(case_.regulatory_mutations_near_mutations_filename);

            if (regulatoryMutationsByVAF == null)
            {
                Console.WriteLine("Unable to read " + case_.regulatory_mutations_near_mutations_filename);
                return;
            }

            foreach (var geneMeasurement in regulatoryMutationsByVAF.Select(x => x.Value).ToList())
            {
                if (geneMeasurement.nNonSilentMutations != 1)
                {
                    continue;   // We only care about single mutation tumors
                }

                if (!state.measurements.ContainsKey(geneMeasurement.hugo_symbol))
                {
                    state.measurements.Add(geneMeasurement.hugo_symbol, new List<double>[3, nRegionsToUse]);
                }

                int vafIndex;

                if (geneMeasurement.nNonSilentMutationsWithHighVAF + geneMeasurement.nNonSilentMutationsWithLowVAF + geneMeasurement.nNonSilentMutationsWithModerateVAF != 1)
                {
                    throw new Exception("Values don't add up in for VAF stratification.  Case ID " + case_.case_id + " filename " + case_.regulatory_mutations_near_mutations_filename + " low " + geneMeasurement.nNonSilentMutationsWithLowVAF + 
                        " moderate " + geneMeasurement.nNonSilentMutationsWithModerateVAF + " high " + geneMeasurement.nNonSilentMutationsWithHighVAF + " total " + geneMeasurement.nNonSilentMutations + " silent " + geneMeasurement.nSilentMutations + 
                        " gene " + geneMeasurement.hugo_symbol);
                }

                if (geneMeasurement.nNonSilentMutationsWithLowVAF == 1)
                {
                    vafIndex = 0;
                }
                else if (geneMeasurement.nNonSilentMutationsWithModerateVAF == 1)
                {
                    vafIndex = 1;
                }
                else
                {
                    vafIndex = 2;
                }

                for (int range = 1; range < nRegionsToUse; range++)
                {
                    if (geneMeasurement.mutationsPerCoveredBaseByRange[range] != double.NegativeInfinity)
                    {
                        if (state.measurements[geneMeasurement.hugo_symbol][vafIndex, range] == null)
                        {
                            state.measurements[geneMeasurement.hugo_symbol][vafIndex, range] = new List<double>();
                        }
                        state.measurements[geneMeasurement.hugo_symbol][vafIndex, range].Add(geneMeasurement.mutationsPerCoveredBaseByRange[range]);
                    }
                }
            } // foreach gene
        } // HandleOneCase

        static void FinishUp(PerThreadState state)
        {
            lock (globalMeasurements)
            {
                foreach (var measurementEntry in state.measurements)
                {
                    var hugo_symbol = measurementEntry.Key;
                    var data = measurementEntry.Value;

                    if (globalMeasurements.ContainsKey(hugo_symbol))
                    {
                        for (int mutationCount = 0; mutationCount <= 2; mutationCount++)
                        {
                            for (int range = 0; range < nRegionsToUse; range++)
                            {
                                if (data[mutationCount, range] == null)
                                {
                                    continue;
                                }

                                if (globalMeasurements[hugo_symbol][mutationCount, range] == null)
                                {
                                    globalMeasurements[hugo_symbol][mutationCount, range] = data[mutationCount, range];
                                } else
                                {
                                    globalMeasurements[hugo_symbol][mutationCount, range].AddRange(data[mutationCount, range]);
                                } // If the mutationCount,range exists in the global state
                            } // for each range
                        } // for each mutation count
                    } else
                    {
                        globalMeasurements.Add(hugo_symbol, data);
                    }
                } // foreach measurement
            } // lock
        } // FinishUp

        static void HandleOneGene(KeyValuePair<string, List<double>[,]> geneDataEntry, int state)
        {
            var hugo_symbol = geneDataEntry.Key;
            var geneData = geneDataEntry.Value;

            var highVsNotHigh = new double[nRegionsToUse];
            var highVsModerate = new double[nRegionsToUse];

            double smallestPValueThisGene = 2;
            int nPValues = 0;

            for (int region = 1; region < nRegionsToUse; region++)
            {
                highVsModerate[region] = double.NegativeInfinity;
                highVsNotHigh[region] = double.NegativeInfinity;

                int[] counts = new int[3];
                for (int vafIndex = 0; vafIndex < 3; vafIndex++)
                {
                    counts[vafIndex] = countList(geneData[vafIndex, region]);
                }

                if (counts[2] < nValuesNeeded || counts[0] + counts[1] < nValuesNeeded)
                {
                    continue;   // Not enough data, skip it
                }

                var data = new List<DoubleAndBool>();

                geneData[2, region].ForEach(x => data.Add(new DoubleAndBool(x, true)));
                if (counts[1] > 0)
                {
                    geneData[1, region].ForEach(x => data.Add(new DoubleAndBool(x, false)));
                }

                if (counts[1] >= nValuesNeeded)
                {
                    bool enoughData;
                    bool reversed;
                    double nFirstGroup, nSecondGroup, U, z;

                    highVsModerate[region] = ASETools.MannWhitney<DoubleAndBool>.ComputeMannWhitney(data, data[0], x => x.isHigh, x => x.value, out enoughData, out reversed, out nFirstGroup, out nSecondGroup, out U, out z, true, nValuesNeeded);

                    if (enoughData)
                    {
                        nPValues++;
                        smallestPValueThisGene = Math.Min(smallestPValueThisGene, highVsModerate[region]);
                    }
                    else
                    {
                        highVsModerate[region] = double.NegativeInfinity;
                    }
                } // if enough for one vs. more than one

                if (counts[0] > 0)
                {
                    geneData[0, region].ForEach(x => data.Add(new DoubleAndBool(x, false)));
                }

                if (counts[0] + counts[1] >= nValuesNeeded)
                {
                    bool enoughData;
                    bool reversed;
                    double nFirstGroup, nSecondGroup, U, z;

                    highVsNotHigh[region] = ASETools.MannWhitney<DoubleAndBool>.ComputeMannWhitney(data, data[0], x => x.isHigh, x => x.value, out enoughData, out reversed, out nFirstGroup, out nSecondGroup, out U, out z, true, nValuesNeeded);

                    if (enoughData)
                    {
                        nPValues++;
                        smallestPValueThisGene = Math.Min(smallestPValueThisGene, highVsNotHigh[region]);
                    }
                    else
                    {
                        highVsNotHigh[region] = double.NegativeInfinity;
                    }
                } // if enough for one vs. not one

            } // for each region

            if (nPValues > 0)
            {
                lock (outputFile)
                {
                    outputFile.Write(hugo_symbol + "\t" + smallestPValueThisGene);
                    for (int region = 1; region < nRegionsToUse; region++)
                    {
                        outputFile.Write("\t" + countList(geneData[0, region]) + "\t" + countList(geneData[1, region]) + "\t" + countList(geneData[2, region]) + "\t" + valueOrStar(highVsNotHigh[region]) + "\t" +
                            valueOrStar(highVsModerate[region]) + "\t" + meanOrStar(geneData[0, region]) + "\t" + meanOrStar(geneData[1, region]) + "\t" + meanOrStar(geneData[2, region]) +
                            "\t" + ratioOrStar(geneData[2, region], geneData[1, region]));
                    }
                    outputFile.WriteLine();

                    smallestPValue = Math.Min(smallestPValue, smallestPValueThisGene);
                    totalValidPValues += nPValues;
                } // lock outputFile
            } // if we found any P values
        } // HandleOneGene

    }
}
