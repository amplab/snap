using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using System.Threading;
using System.IO;
using ASELib;
using System.Diagnostics;

namespace ASEHeatMap
{
    class Program
    {
        const int nASEs = 101;      // 0 - 1 inclusive by .01
        const int nRegions = 20;    // Powers of two * 1Kb starting with 1 KB (ie., no 0) and extending to 256MB (insclusive), plus one for rest of autosome.  This is selected because 256MB is bigger than chr1.
        class HeatMapEntry
        {
            public double totalASE = 0;
            public double totalASESquared = 0;
            public int n = 0;
            public int totalSamples = 0;
            public double maxASE = -1;
            public double minASE = 2;
        }

        static HeatMapEntry [,] globalHeatMap = new HeatMapEntry[nASEs, nRegions];
        static int[] globalASECount = new int[nASEs];

        static void initializeHeatMap(HeatMapEntry[,] heatMap, int[] aseCount)
        {
            for (int aseIndex = 0; aseIndex < nASEs; aseIndex++)
            {
                aseCount[aseIndex] = 0;

                for (int i = 0; i < nRegions; i++)
                {
                    heatMap[aseIndex, i] = new HeatMapEntry();
                }
            }
        }
 
        static void Main(string[] args)
        {
            var timer = new Stopwatch();
            timer.Start();

            var configuration = ASETools.ASEConfirguation.loadFromFile(args);

            if (null == configuration)
            {
                Console.WriteLine("Unable to load configuration.");
                return;
            }

            if (configuration.commandLineArgs.Count() != 0)
            {
                Console.WriteLine("usage: ASEHeatMap");
                return;
            }

            var cases = ASETools.Case.LoadCases(configuration.casesFilePathname);

            if (null == cases)
            {
                Console.WriteLine("Unable to load cases.");
                return;
            }

            List<string> caseIdsToProcess = new List<string>();

            foreach (var caseEntry in cases)
            {
                caseIdsToProcess.Add(caseEntry.Key);
            }

            initializeHeatMap(globalHeatMap, globalASECount);

            Console.Write("Progress (1 dot/100 cases): ");

            var threads = new List<Thread>();
            for (int i = 0; i < Environment.ProcessorCount * 2 /* *2 to account for hyperthreading*/; i++)
            {
                threads.Add(new Thread(() => WorkerThread(caseIdsToProcess, cases)));
            }

            threads.ForEach(t => t.Start());
            threads.ForEach(t => t.Join());

            Console.WriteLine();    // Newline after progress dots

            var outputFile = ASETools.CreateStreamWriterWithRetry(configuration.finalResultsDirectory + ASETools.heatMapFilename);

            if (null == outputFile)
            {
                Console.WriteLine("Error opening output file.");
                return;
            }

            outputFile.Write("Allele Specific Expression\tCount of samples at this ASE");
            for (int i = 0; i < nRegions; i++)
            {
                outputFile.Write("\tMean ASE Within " + distanceIndexToString(i));
            }

            for (int i = 0; i < nRegions; i++)
            {
                outputFile.Write("\tNumber of ASE samples that have anything within " + distanceIndexToString(i));
            }

            for (int i = 0; i < nRegions; i++)
            {
                outputFile.Write("\tMean number of ASE samples within " + distanceIndexToString(i));
            }

            for (int i = 0; i < nRegions; i++)
            {
                outputFile.Write("\t Standard Deviation of ASE within " + distanceIndexToString(i));
            }

            outputFile.WriteLine();

            for (int aseIndex = 0; aseIndex < nASEs; aseIndex++)
            {
                outputFile.Write(aseIndex * .01 + "\t" + globalASECount[aseIndex]);

                for (int i = 0; i < nRegions; i++)
                {
                    if (globalHeatMap[aseIndex, i].n == 0)
                    {
                        outputFile.Write("\t0");
                    } else
                    {
                        outputFile.Write("\t" + globalHeatMap[aseIndex, i].totalASE / globalHeatMap[aseIndex, i].n);
                    }
                }

                for (int i = 0; i < nRegions; i++)
                {
                    outputFile.Write("\t" + globalHeatMap[aseIndex, i].n);
                }

                for (int i = 0; i < nRegions; i++)
                {
                    if (globalHeatMap[aseIndex, i].n == 0)
                    {
                        outputFile.Write("\t0");
                    }
                    else
                    {
                        outputFile.Write("\t" + (double)globalHeatMap[aseIndex, i].totalSamples / globalHeatMap[aseIndex, i].n);
                    }
                }

                for (int i = 0; i < nRegions; i++)
                {
                    if (globalHeatMap[aseIndex, i].n == 0)
                    {
                        outputFile.Write("\t0");
                    }
                    else
                    {
                        outputFile.Write("\t" + Math.Sqrt(globalHeatMap[aseIndex, i].n  * globalHeatMap[aseIndex, i].totalASE - globalHeatMap[aseIndex,i].totalASESquared) / globalHeatMap[aseIndex,i].n);
                    }
                }

                outputFile.WriteLine();
            } // foreach ASE index

            outputFile.WriteLine("**done**");
            outputFile.Close();

            Console.WriteLine("Processed " + cases.Count() + " cases in " + ASETools.ElapsedTimeInSeconds(timer));

        } // Main

        static string distanceIndexToString(int distanceIndex)
        {
            switch (distanceIndex)
            {
                case 0: return "1KB";
                case 1: return "2KB";
                case 2: return "4KB";
                case 3: return "8KB";
                case 4: return "16KB";
                case 5: return "32KB";
                case 6: return "64KB";
                case 7: return "128KB";
                case 8: return "256KB";
                case 9: return "512KB";
                case 10: return "1MB";
                case 11: return "2MB";
                case 12: return "4MB";
                case 13: return "8MB";
                case 14: return "16MB";
                case 15: return "32MB";
                case 16: return "64MB";
                case 17: return "128MB";
                case 18: return "256MB";
                case 19: return "whole autosome";
            }

            return "Illegal index??";
        }

        static int nCasesProcessed = 0;

        static void WorkerThread(List<string> caseIdsToProcess, Dictionary<string, ASETools.Case> cases)
        {
            HeatMapEntry[,] localHeatMap = new HeatMapEntry[nASEs, nRegions];
            int[] localASECount = new int[nASEs];

            initializeHeatMap(localHeatMap, localASECount);

            while (true)
            {
                ASETools.Case case_;
                lock (caseIdsToProcess)
                {
                    if (caseIdsToProcess.Count() == 0)
                    {
                        //
                        // Merge our local heat map into the global one.
                        //

                        for (int aseIndex = 0; aseIndex < nASEs; aseIndex++)
                        {
                            globalASECount[aseIndex] += localASECount[aseIndex];

                            for (int i = 0; i < nRegions; i++)
                            {
                                globalHeatMap[aseIndex, i].n += localHeatMap[aseIndex, i].n;
                                globalHeatMap[aseIndex, i].totalASE += localHeatMap[aseIndex, i].totalASE;
                                globalHeatMap[aseIndex, i].totalASESquared += localHeatMap[aseIndex, i].totalASESquared;
                                globalHeatMap[aseIndex, i].totalSamples += localHeatMap[aseIndex, i].totalSamples;
                                globalHeatMap[aseIndex, i].minASE = Math.Min(globalHeatMap[aseIndex, i].minASE, localHeatMap[aseIndex, i].minASE);
                                globalHeatMap[aseIndex, i].maxASE = Math.Max(globalHeatMap[aseIndex, i].maxASE, localHeatMap[aseIndex, i].maxASE);
                            }
                        }
                        return;
                    }

                    case_ = cases[caseIdsToProcess[0]];
                    caseIdsToProcess.RemoveAt(0);
                }

                if (case_.annotated_selected_variants_filename == "")
                {
                    Console.WriteLine("Skipping case " + case_.case_id + " because it does not have annotated selected variants.");
                    continue;
                }

                var annotatedSelectedVariants = ASETools.AnnotatedVariant.readFile(case_.annotated_selected_variants_filename);

                foreach (var variant in annotatedSelectedVariants)
                {
                    if (variant.somaticMutation)
                    {
                        continue;   // We're only considering germline loci, at least for now, mostly because we don't want to have to filter subclones/homozygous loci, etc.
                    }

                    if (variant.tumorRNAReadCounts.nMatchingAlt + variant.tumorRNAReadCounts.nMatchingReference < 10 || !ASETools.isChromosomeAutosomal(variant.contig))
                    {
                        continue;
                    }

                    int aseValue = (int)(variant.tumorRNAReadCounts.AlleleSpecificValue() * 100 + .5); // +.5 is for rounding.

                    localASECount[aseValue]++;

                    for (int index = 0; index < nRegions; index++)
                    {
                        int distance = 1000 * (1 << index);

                        int n = 0;
                        double totalASE = 0;

                        foreach (var peerVariant in annotatedSelectedVariants)
                        {
                            if (peerVariant.somaticMutation || peerVariant.locus == variant.locus && peerVariant.contig == variant.contig || peerVariant.tumorRNAReadCounts.nMatchingReference + peerVariant.tumorRNAReadCounts.nMatchingAlt < 10)
                            {
                                continue;   // Either somatic variant or self or too few reads
                            }

                            if (peerVariant.contig == variant.contig && Math.Abs(peerVariant.locus - variant.locus) <= distance ||
                                index == nRegions - 1 && ASETools.isChromosomeAutosomal(peerVariant.contig))                            // Remember, the last region is "whole autosome"
                            {
                                n++;
                                totalASE += peerVariant.tumorRNAReadCounts.AlleleSpecificValue();
                            }
                        } // foreach peer variant

                        if (n != 0)
                        {
                            double meanASE = totalASE / n;

                            localHeatMap[aseValue, index].n++;
                            localHeatMap[aseValue, index].totalSamples += n;
                            localHeatMap[aseValue, index].totalASE += meanASE;
                            localHeatMap[aseValue, index].totalASESquared += meanASE * meanASE;
                            localHeatMap[aseValue, index].maxASE = Math.Max(globalHeatMap[aseValue, index].maxASE, meanASE);
                            localHeatMap[aseValue, index].minASE = Math.Min(globalHeatMap[aseValue, index].minASE, meanASE);
                        } // If we found anything
                    } // foreach region
                } // foreach variant

                if (Interlocked.Increment(ref nCasesProcessed) % 100 == 0)
                {
                    Console.Write(".");
                }
            } // while true
        }
    }
}
