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

        class Measurements
        {
            public Dictionary<string, List<double>[,]> byRegion = new Dictionary<string, List<double>[,]>(); // Hugo symbol -> array of [0,1,2] for low, middle, high count -> array by region index -> measurements
            public Dictionary<string, Dictionary<string, List<double>[]>> byGeneHancer = new Dictionary<string, Dictionary<string, List<double>[]>>(); // Hugo symbol -> gene Hancer -> array of [0,1,2] for low, middle, high count->measurements 
            public Dictionary<string, Dictionary<string, List<double>[]>> byGeneHancerForAllTumors = new Dictionary<string, Dictionary<string, List<double>[]>>(); // Hugo symbol->gene hancer-> array of [0,1,2] for mutation count in gene
        }

        class PerThreadState
        {
            public Dictionary<string, List<double>[,]> measurements = new Dictionary<string, List<double>[,]>();  // Hugo symbol -> array of [0,1,2] for low, middle, high count -> array by region index -> measurements
            public Dictionary<string, Dictionary<string, List<double>[]>> byGeneHancerForSingleMutationTumors = new Dictionary<string, Dictionary<string, List<double>[]>>(); // Hugo symbol -> gene Hancer -> array of [0,1,2] for low, middle, high count->measurements 
            public Dictionary<string, Dictionary<string, List<double>[]>> byGeneHancerForAllTumors = new Dictionary<string, Dictionary<string, List<double>[]>>(); // Hugo symbol->gene hancer-> array of [0,1,2] for mutation count in gene

            void merge(PerThreadState peer)
            {
                foreach (var entry in peer.measurements)
                {
                    var hugo_symbol = entry.Key;
                    var lists = entry.Value;

                    if (!measurements.ContainsKey(hugo_symbol))
                    {
                        measurements.Add(hugo_symbol, new List<double>[3, nRegionsToUse]);
                    }

                    for (int i = 0; i < 3; i++)
                    {
                        for (int region = 0; region < nRegionsToUse; region++)
                        {
                            if (measurements[hugo_symbol][i, region] == null)
                            {
                                measurements[hugo_symbol][i, region] = peer.measurements[hugo_symbol][i, region];
                            } else
                            {
                                measurements[hugo_symbol][i, region].AddRange(peer.measurements[hugo_symbol][i, region]);
                            }
                        } // region
                    } // 0, 1, 2
                } // gene

                mergeGeneHancers(byGeneHancerForSingleMutationTumors, peer.byGeneHancerForSingleMutationTumors);
                mergeGeneHancers(byGeneHancerForAllTumors, peer.byGeneHancerForAllTumors);
            } // merge

            static void mergeGeneHancers(Dictionary<string, Dictionary<string, List<double>[]>> destination, Dictionary<string, Dictionary<string, List<double>[]>> peer)
            {
                foreach (var entry in peer)
                {
                    var hugo_symbol = entry.Key;
                    var geneHancers = entry.Value;

                    if (!destination.ContainsKey(hugo_symbol))
                    {
                        destination.Add(hugo_symbol, new Dictionary<string, List<double>[]>());
                    }

                    foreach (var geneHancerEntry in geneHancers)
                    {
                        var geneHancerId = geneHancerEntry.Key;
                        var lists = geneHancerEntry.Value;

                        if (!destination[hugo_symbol].ContainsKey(geneHancerId))
                        {
                            destination[hugo_symbol].Add(geneHancerId, new List<double>[3]);
                        }

                        for (int i = 0; i < 3; i++)
                        {
                            if (destination[hugo_symbol][geneHancerId][i] == null)
                            {
                                destination[hugo_symbol][geneHancerId][i] = peer[hugo_symbol][geneHancerId][i];
                            } else
                            {
                                destination[hugo_symbol][geneHancerId][i].AddRange(peer[hugo_symbol][geneHancerId][i]);
                            }
                        } // 0, 1, 2

                    } // geneHancer Entry
                } // gene
            } // mergeGeneHancers
        }


        static PerThreadState globalState = new PerThreadState();

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
            if (denominator == null || numerator == null || denominator.Count() == 0 || numerator.Count() == 0)
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

            outputFile.Write("Hugo symbol\tmin P from regions");
            for (int region = 1; region < nRegionsToUse; region++)
            {
                var regionName = ASETools.regionIndexToString(region);
                outputFile.Write("\t" + regionName + " nLow\t" + regionName + " nModerate\t" + regionName + " nHigh\t" + regionName + " P high vs. not high\t" + regionName + " P high vs. moderate\t" + regionName + " low mean\t" +
                    regionName + " moderate mean\t" + regionName + " high mean\t" + regionName + " high mean/moderate mean");
            }

            outputFile.WriteLine("\tGene Hancers by mutation count (id=n0, n1, n>1, p 1 != >1, p 1 != !1)\tSmallest P value for gene hancers by mutation count\tSmallest P Value for gene hancers by mutation count where 1 mutation > max(0 mutation, >1 mutation)\t" +
                "Gene Hancers for single mutant tumors (id = nLow, nMid, nHigh, p Mid != High, p Mid != !Mid)\tSmallest P value for gene hancers for single mutant tumors");

            var genesToProcess = globalState.measurements.Select(x => x.Key).ToList();
            var geneThreading = new ASETools.WorkerThreadHelper<string, int>(genesToProcess, HandleOneGene, null, null, 100);
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
                var hugo_symbol = geneMeasurement.hugo_symbol;

                updateGeneHancer(state.byGeneHancerForAllTumors, hugo_symbol, geneMeasurement.geneHancerMutationsPerBase, ASETools.ZeroOneMany(geneMeasurement.nNonSilentMutations));

                if (geneMeasurement.nNonSilentMutations != 1)
                {
                    continue;   // We only care about single mutation tumors
                }

                int vafIndex;

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

                updateGeneHancer(state.byGeneHancerForSingleMutationTumors, hugo_symbol, geneMeasurement.geneHancerMutationsPerBase, vafIndex);

                if (!state.measurements.ContainsKey(geneMeasurement.hugo_symbol))
                {
                    state.measurements.Add(geneMeasurement.hugo_symbol, new List<double>[3, nRegionsToUse]);
                }


                if (geneMeasurement.nNonSilentMutationsWithHighVAF + geneMeasurement.nNonSilentMutationsWithLowVAF + geneMeasurement.nNonSilentMutationsWithModerateVAF != 1)
                {
                    throw new Exception("Values don't add up in for VAF stratification.  Case ID " + case_.case_id + " filename " + case_.regulatory_mutations_near_mutations_filename + " low " + geneMeasurement.nNonSilentMutationsWithLowVAF + 
                        " moderate " + geneMeasurement.nNonSilentMutationsWithModerateVAF + " high " + geneMeasurement.nNonSilentMutationsWithHighVAF + " total " + geneMeasurement.nNonSilentMutations + " silent " + geneMeasurement.nSilentMutations + 
                        " gene " + geneMeasurement.hugo_symbol);
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

        static void updateGeneHancer(Dictionary<string, Dictionary<string, List<double>[]>>byGeneHancer, string hugo_symbol, Dictionary<string, double> geneHancerMutationsPerBase, int index)
        {
            if (!byGeneHancer.ContainsKey(hugo_symbol))
            {
                byGeneHancer.Add(hugo_symbol, new Dictionary<string, List<double>[]>());
            }

            foreach (var geneHancerEntry in geneHancerMutationsPerBase)
            {
                var geneHancerId = geneHancerEntry.Key;
                var mutationsPerCoveredBase = geneHancerEntry.Value;

                if (!byGeneHancer[hugo_symbol].ContainsKey(geneHancerId))
                {
                    byGeneHancer[hugo_symbol].Add(geneHancerId, new List<double>[3]);
                    for (int i = 0; i < 3; i++)
                    {
                        byGeneHancer[hugo_symbol][geneHancerId][i] = new List<double>();
                    }
                }

                byGeneHancer[hugo_symbol][geneHancerId][index].Add(mutationsPerCoveredBase);
            } // foreach geneHancer entry
        } // updateGeneHancer

        static void FinishUp(PerThreadState state)
        {
            lock (globalState)
            {
                foreach (var measurementEntry in state.measurements)
                {
                    var hugo_symbol = measurementEntry.Key;
                    var data = measurementEntry.Value;

                    if (globalState.measurements.ContainsKey(hugo_symbol))
                    {
                        for (int mutationCount = 0; mutationCount <= 2; mutationCount++)
                        {
                            for (int range = 0; range < nRegionsToUse; range++)
                            {
                                if (data[mutationCount, range] == null)
                                {
                                    continue;
                                }

                                if (globalState.measurements[hugo_symbol][mutationCount, range] == null)
                                {
                                    globalState.measurements[hugo_symbol][mutationCount, range] = data[mutationCount, range];
                                } else
                                {
                                    globalState.measurements[hugo_symbol][mutationCount, range].AddRange(data[mutationCount, range]);
                                } // If the mutationCount,range exists in the global state
                            } // for each range
                        } // for each mutation count
                    } else
                    {
                        globalState.measurements.Add(hugo_symbol, data);
                    }
                } // foreach measurement

                mergeGeneHancer(globalState.byGeneHancerForAllTumors, state.byGeneHancerForAllTumors);
                mergeGeneHancer(globalState.byGeneHancerForSingleMutationTumors, state.byGeneHancerForSingleMutationTumors);
            } // lock
        } // FinishUp

        static void mergeGeneHancer(Dictionary<string, Dictionary<string, List<double>[]>> mergeInto, Dictionary<string, Dictionary<string, List<double>[]>> mergeFrom)
        {
            foreach (var hugo_symbol in mergeFrom.Select(x => x.Key))
            {
                if (!mergeInto.ContainsKey(hugo_symbol))
                {
                    mergeInto.Add(hugo_symbol, mergeFrom[hugo_symbol]);
                }
                else
                {
                    foreach (var geneHancerId in mergeFrom[hugo_symbol].Select(x => x.Key))
                    {
                        if (!mergeInto[hugo_symbol].ContainsKey(geneHancerId))
                        {
                            mergeInto[hugo_symbol].Add(geneHancerId, mergeFrom[hugo_symbol][geneHancerId]);
                        }
                        else
                        {
                            for (int i = 0; i < 3; i++)
                            {
                                mergeInto[hugo_symbol][geneHancerId][i].AddRange(mergeFrom[hugo_symbol][geneHancerId][i]);
                            } // 0, 1, 2
                        }
                    } // geneHancerId
                }
            } // hugo_symbol
        } // mergeGeneHancer

        static void HandleOneGene(string hugo_symbol, int state)
        {
            var geneData = globalState.measurements[hugo_symbol];

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


            var outputLine = hugo_symbol + "\t" + smallestPValueThisGene;
            for (int region = 1; region < nRegionsToUse; region++)
            {
                outputLine += "\t" + countList(geneData[0, region]) + "\t" + countList(geneData[1, region]) + "\t" + countList(geneData[2, region]) + "\t" + valueOrStar(highVsNotHigh[region]) + "\t" +
                    valueOrStar(highVsModerate[region]) + "\t" + meanOrStar(geneData[0, region]) + "\t" + meanOrStar(geneData[1, region]) + "\t" + meanOrStar(geneData[2, region]) +
                    "\t" + ratioOrStar(geneData[2, region], geneData[1, region]);
            }

            outputLine += "\t";
            outputGeneHancers(ref outputLine, hugo_symbol, globalState.byGeneHancerForAllTumors, ref nPValues, ref smallestPValue, 1, 2, 0, true);
            outputLine += "\t";
            outputGeneHancers(ref outputLine, hugo_symbol, globalState.byGeneHancerForSingleMutationTumors, ref nPValues, ref smallestPValue, 2, 1, 0, false);

            if (nPValues > 0)
            {
                lock (outputFile)
                {
                    smallestPValue = Math.Min(smallestPValue, smallestPValueThisGene);
                    totalValidPValues += nPValues;

                    outputFile.WriteLine(outputLine);
                }
            } // if we found any P values
        } // HandleOneGene

        static void outputGeneHancers(ref string outputLine, string hugo_symbol, Dictionary<string, Dictionary<string, List<double>[]>> allHeneHancers, ref int nPValues, ref double smallestPValue, int falseGroup, int trueGroup1, int trueGroup2, bool doZeroOneTwoPValue)
        {
            if (!allHeneHancers.ContainsKey(hugo_symbol))
            {
                outputLine += "\t*\t*";
                return;
            }

            var geneHancersForThisGene = allHeneHancers[hugo_symbol];

            bool anyWritten = false;

            double minPThisGeneHancerSet = 2;
            double minPOneLessThanZeroAndTwo = 2;

            foreach (var geneHancerId in geneHancersForThisGene.Select(x => x.Key))
            {
                var thisGeneHancer = geneHancersForThisGene[geneHancerId];

                if (anyWritten)
                {
                    outputLine += ";";
                } else
                {
                    anyWritten = true;
                }

                outputLine += geneHancerId + "=";

                for (int i = 0; i < 3; i++)
                {
                    if (thisGeneHancer[i].Count() == 0)
                    {
                        outputLine += "0(*,*),";
                    } else
                    {
                        outputLine += thisGeneHancer[i].Count() + "(" + thisGeneHancer[i].Average() + "," + ((double)thisGeneHancer[i].Where(x => x != 0).Count() / thisGeneHancer[i].Count()) + "),";
                    }
                }

                var data = new List<DoubleAndBool>();

                thisGeneHancer[falseGroup].ForEach(x => data.Add(new DoubleAndBool(x, false)));
                thisGeneHancer[trueGroup1].ForEach(x => data.Add(new DoubleAndBool(x, true)));

                bool zeroOneTwoThisGeneHancer;

                if (doZeroOneTwoPValue && thisGeneHancer[0].Count() > 0 && thisGeneHancer[1].Count() > 0 && thisGeneHancer[2].Count() > 0)
                {
                    var zeroAverage = thisGeneHancer[0].Average();
                    var oneAverage = thisGeneHancer[1].Average();
                    var twoAverage = thisGeneHancer[2].Average();

                    zeroOneTwoThisGeneHancer = oneAverage > zeroAverage && oneAverage > twoAverage;
                } else
                {
                    zeroOneTwoThisGeneHancer = false;
                }

                var p = outputGeneHancerPValue(ref outputLine, data, ref nPValues, ref smallestPValue);
                minPThisGeneHancerSet = Math.Min(minPThisGeneHancerSet, p);

                if (zeroOneTwoThisGeneHancer)
                {
                    minPOneLessThanZeroAndTwo = Math.Min(minPOneLessThanZeroAndTwo, p);
                }

                outputLine += ",";

                thisGeneHancer[trueGroup2].ForEach(x => data.Add(new DoubleAndBool(x, true)));
                p = outputGeneHancerPValue(ref outputLine, data, ref nPValues, ref smallestPValue);

                minPThisGeneHancerSet = Math.Min(minPThisGeneHancerSet, p);

                if (zeroOneTwoThisGeneHancer)
                {
                    minPOneLessThanZeroAndTwo = Math.Min(minPOneLessThanZeroAndTwo, p);
                }

            }

            if (minPThisGeneHancerSet < 2)
            {
                outputLine += "\t" + minPThisGeneHancerSet;
            } else
            {
                outputLine += "\t*";
            }

            if (doZeroOneTwoPValue)
            {
                if (minPOneLessThanZeroAndTwo < 2)
                {
                    outputLine += "\t" + minPOneLessThanZeroAndTwo;
                } else
                {
                    outputLine += "\t*";
                }
            }
        } // outputGeneHancers

        static double outputGeneHancerPValue(ref string outputLine, List<DoubleAndBool> data, ref int nPValues, ref double smallestPValue)
        {
            if (!(data.Where(x => x.isHigh).Count() >= nValuesNeeded && data.Where(x => !x.isHigh).Count() >= nValuesNeeded))
            {
                outputLine += "*";
                return 2;
            }

            bool enoughData;
            bool reversed;
            double nFirstGroup, nSecondGroup, U, z;

            var p = ASETools.MannWhitney<DoubleAndBool>.ComputeMannWhitney(data, data[0], x => x.isHigh, x => x.value, out enoughData, out reversed, out nFirstGroup, out nSecondGroup, out U, out z, true, nValuesNeeded);

            if (enoughData)
            {
                outputLine += p;
                nPValues++;
                smallestPValue = Math.Min(smallestPValue, p);
                return p;
            }
            else
            {
                outputLine += "*";
                return 2;
            }
        } // outputGeneHancerPValue
    }
}
