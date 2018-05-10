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

        class PerThreadState
        {
            public Dictionary<string, List<double>[,]> measurements = new Dictionary<string, List<double>[,]>();  // Hugo symbol -> array of [0,1,2] for mutation count -> array by region index -> measurements
        }
        
        static Dictionary<string, List<double>[,]>  globalMeasurements = new Dictionary<string, List<double>[,]>();

        class DoubleAndBool : IComparer<DoubleAndBool>
        {
            public readonly double value;
            public readonly bool isOne;

            public DoubleAndBool(double value_, bool isOne_)
            {
                value = value_;
                isOne = isOne_;
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

        const int nValuesNeeded = 10;
        static int totalValidPValues = 0;
        static double smallestPValue = 2;

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

            var outputFile = ASETools.CreateStreamWriterWithRetry(configuration.finalResultsDirectory + ASETools.cisRegulatoryMutationsByMutationCountFilename);
            if (null == outputFile)
            {
                Console.WriteLine("Unable to open output file " + configuration.finalResultsDirectory + ASETools.cisRegulatoryMutationsByMutationCountFilename);
                return;
            }

            outputFile.Write("Hugo symbol\tmin P\t");
            for (int region = 1; region < ASETools.nRegions; region++)
            {
                var regionName = ASETools.regionIndexToString(region);
                outputFile.Write("\t" + regionName + " nZero\t" + regionName + " nOne\t" + regionName + " nMoreThanOne\t" + regionName + " P one vs. not one\t" + regionName + " P one vs. more than one\t" + regionName + " zero mean\t" +
                    regionName + " one mean\t" + regionName + " more than one mean");
            }
            outputFile.WriteLine();
            outputFile.WriteLine(ASETools.ElapsedTimeInSeconds(timer));

            var genesToProcess = globalMeasurements.ToList();
            var geneThreading = new ASETools.WorkerThreadHelper<KeyValuePair<string, List<double>[,]>, PerGeneThreadState>(genesToProcess, HandleOneGene, FinishGeneThread, null, 100);
            Console.WriteLine("Processing " + genesToProcess.Count()+  " genes, one dot/100 genes:");




            foreach (var geneDataEntry in globalMeasurements)
            {
                var hugo_symbol = geneDataEntry.Key;
                var geneData = geneDataEntry.Value;

                var oneVsNotOne = new double[ASETools.nRegions];
                var oneVsMany = new double[ASETools.nRegions];

                bool anyPValues = false;    // If we don't have enough data to compute even one P value, then we skip this gene.
                double smallestPValueThisGene = 2;

                for (int region = 1; region < ASETools.nRegions; region++)
                {
                    oneVsMany[region] = double.NegativeInfinity;
                    oneVsNotOne[region] = double.NegativeInfinity;

                    int []counts = new int[3];
                    for (int zeroOneMore = 0; zeroOneMore < 3; zeroOneMore++)
                    {
                        counts[zeroOneMore] = countList(geneData[zeroOneMore, region]);
                    }

                    if (counts[1] < nValuesNeeded || counts[0] + counts[2] < nValuesNeeded)
                    {
                        continue;   // Not enough data, skip it
                    }

                    var data = new List<DoubleAndBool>();

                    geneData[1, region].ForEach(x => data.Add(new DoubleAndBool(x, true)));
                    if (counts[2] > 0)
                    {
                        geneData[2, region].ForEach(x => data.Add(new DoubleAndBool(x, false)));
                    }

                    if (counts[2] >= nValuesNeeded)
                    {
                        bool enoughData;
                        bool reversed;
                        double nFirstGroup, nSecondGroup, U, z;

                        oneVsMany[region] = ASETools.MannWhitney<DoubleAndBool>.ComputeMannWhitney(data, data[0], x => x.isOne, x => x.value, out enoughData, out reversed, out nFirstGroup, out nSecondGroup, out U, out z, true, nValuesNeeded);

                        if (enoughData)
                        {
                            anyPValues = true;
                            smallestPValueThisGene = Math.Min(smallestPValueThisGene, oneVsMany[region]);
                        } else
                        {
                            oneVsMany[region] = double.NegativeInfinity;
                        }
                    } // if enough for one vs. more than one

                    if (counts[0] > 0)
                    {
                        geneData[0, region].ForEach(x => data.Add(new DoubleAndBool(x, false)));
                    }

                    if (counts[0] + counts[2] >= nValuesNeeded)
                    {
                        bool enoughData;
                        bool reversed;
                        double nFirstGroup, nSecondGroup, U, z;

                        oneVsNotOne[region] = ASETools.MannWhitney<DoubleAndBool>.ComputeMannWhitney(data, data[0], x => x.isOne, x => x.value, out enoughData, out reversed, out nFirstGroup, out nSecondGroup, out U, out z, true, nValuesNeeded);

                        if (enoughData)
                        {
                            anyPValues = true;
                            smallestPValueThisGene = Math.Min(smallestPValueThisGene, oneVsNotOne[region]);
                        }
                        else
                        {
                            oneVsNotOne[region] = double.NegativeInfinity;
                        }
                    } // if enough for one vs. not one
                    
                } // for each region

                if (anyPValues)
                {
                    outputFile.Write(hugo_symbol + "\t" + smallestPValue);
                    for (int region = 1; region < ASETools.nRegions; region++)
                    {
                        outputFile.Write("\t" + countList(geneData[0, region]) + "\t" + countList(geneData[1, region]) + "\t" + countList(geneData[2, region]) + "\t" + valueOrStar(oneVsNotOne[region]) + "\t" +
                            valueOrStar(oneVsMany[region]) + "\t" + meanOrStar(geneData[0, region]) + "\t" + meanOrStar(geneData[1, region]) + "\t" + meanOrStar(geneData[2, region]));
                    }
                    outputFile.WriteLine();

                    smallestPValue = Math.Min(smallestPValue, smallestPValueThisGene);
                }
            } // foreach gene

            outputFile.WriteLine("**done**");
            outputFile.Close();

            Console.WriteLine(ASETools.ElapsedTimeInSeconds(timer));
            Console.WriteLine("Total of " + totalValidPValues + " P values, smallest " + smallestPValue);
        } // Main

        static void HandleOneCase(ASETools.Case case_, PerThreadState state)
        {
            var regulatoryMutationsByMutationCount = ASETools.RegulatoryMutationsNearGene.readFile(case_.regulatory_mutations_near_mutations_filename);

            if (regulatoryMutationsByMutationCount == null)
            {
                Console.WriteLine("Unable to read " + case_.regulatory_mutations_near_mutations_filename);
                return;
            }

            foreach (var geneMeasurement in regulatoryMutationsByMutationCount.Select(x => x.Value).ToList())
            {
                if (!state.measurements.ContainsKey(geneMeasurement.hugo_symbol))
                {
                    state.measurements.Add(geneMeasurement.hugo_symbol, new List<double>[3, ASETools.nRegions]);
                }

                int mutationCount = ASETools.ZeroOneMany(geneMeasurement.nNonSilentMutations);  // Maps >2 -> 2

                for (int range = 1; range < ASETools.nRegions; range++)
                {
                    if (geneMeasurement.mutationsPerCoveredBaseByRange[range] != double.NegativeInfinity)
                    {
                        if (state.measurements[geneMeasurement.hugo_symbol][mutationCount, range] == null)
                        {
                            state.measurements[geneMeasurement.hugo_symbol][mutationCount, range] = new List<double>();
                        }
                        state.measurements[geneMeasurement.hugo_symbol][mutationCount, range].Add(geneMeasurement.mutationsPerCoveredBaseByRange[range]);
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
                            for (int range = 0; range < ASETools.nRegions; range++)
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

    }
}
