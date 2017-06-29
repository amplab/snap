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

            public void merge(HeatMapEntry peer)
            {
                totalASE += peer.totalASE;
                totalASESquared += peer.totalASESquared;
                n += peer.n;
                totalSamples += peer.totalSamples;
                maxASE = Math.Max(maxASE, peer.maxASE);
                minASE = Math.Min(minASE, peer.minASE);
            }
        }

        class HeatMapAndCount
        {
            public HeatMapAndCount()
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

            public void merge(HeatMapAndCount peer)
            {
                for (int aseIndex = 0; aseIndex < nASEs; aseIndex++)
                {
                    aseCount[aseIndex] += peer.aseCount[aseIndex];

                    for (int i = 0; i < nRegions; i++)
                    {
                        heatMap[aseIndex, i].merge(peer.heatMap[aseIndex, i]);
                    }
                }
            }

            public void writeToFile(string filename)
            {
                var outputFile = ASETools.CreateStreamWriterWithRetry(filename);

                if (null == outputFile)
                {
                    Console.WriteLine("Error opening output file.");
                    return;
                }

                outputFile.Write("Allele Specific Expression\tCount of samples at this ASE\tpdf of ASE count");
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

                int totalSamples = aseCount.Sum();

                for (int aseIndex = 0; aseIndex < nASEs; aseIndex++)
                {
                    outputFile.Write(aseIndex * .01 + "\t" + aseCount[aseIndex] + "\t" + (double)aseCount[aseIndex] / totalSamples);

                    for (int i = 0; i < nRegions; i++)
                    {
                        if (heatMap[aseIndex, i].n == 0)
                        {
                            outputFile.Write("\t0");
                        }
                        else
                        {
                            outputFile.Write("\t" + heatMap[aseIndex, i].totalASE / heatMap[aseIndex, i].n);
                        }
                    }

                    for (int i = 0; i < nRegions; i++)
                    {
                        outputFile.Write("\t" + heatMap[aseIndex, i].n);
                    }

                    for (int i = 0; i < nRegions; i++)
                    {
                        if (heatMap[aseIndex, i].n == 0)
                        {
                            outputFile.Write("\t0");
                        }
                        else
                        {
                            outputFile.Write("\t" + (double)heatMap[aseIndex, i].totalSamples / heatMap[aseIndex, i].n);
                        }
                    }

                    for (int i = 0; i < nRegions; i++)
                    {
                        if (heatMap[aseIndex, i].n == 0)
                        {
                            outputFile.Write("\t0");
                        }
                        else
                        {
                            outputFile.Write("\t" + Math.Sqrt(heatMap[aseIndex, i].n * heatMap[aseIndex, i].totalASE - heatMap[aseIndex, i].totalASESquared) / heatMap[aseIndex, i].n);
                        }
                    }

                    outputFile.WriteLine();
                } // foreach ASE index

                outputFile.WriteLine("**done**");
                outputFile.Close();
            }


            public HeatMapEntry[,] heatMap = new HeatMapEntry[nASEs, nRegions];
            public int[] aseCount = new int[nASEs];
        }

        static HeatMapAndCount globalTumorHeatMapAndCount = new HeatMapAndCount();
        static HeatMapAndCount globalNormalHeatMapAndCount = new HeatMapAndCount();

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

            Console.Write("Progress (1 dot/100 cases): ");

            var threads = new List<Thread>();
            for (int i = 0; i < Environment.ProcessorCount ; i++)
            {
                threads.Add(new Thread(() => WorkerThread(caseIdsToProcess, cases)));
            }

            threads.ForEach(t => t.Start());
            threads.ForEach(t => t.Join());

            Console.WriteLine();    // Newline after progress dots

            globalTumorHeatMapAndCount.writeToFile(configuration.finalResultsDirectory + ASETools.tumorHeatMapFilename);
            globalNormalHeatMapAndCount.writeToFile(configuration.finalResultsDirectory + ASETools.normalHeatMapFilename);

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

        delegate ASETools.ReadCounts GetReadCounts(ASETools.AnnotatedVariant variant);

        static void ProcessAnnotatedSelectedVariants(List<ASETools.AnnotatedVariant> annotatedSelectedVariants, HeatMapAndCount heatMapAndCount, GetReadCounts getReadCounts)
        {
            if (annotatedSelectedVariants.Count() == 0 || getReadCounts(annotatedSelectedVariants[0]) == null)
            {
                return;
            }

            foreach (var variant in annotatedSelectedVariants)
            {
                int aseValue = (int)(getReadCounts(variant).AlleleSpecificValue() * 100 + .5); // +.5 is for rounding.

                heatMapAndCount.aseCount[aseValue]++;

                for (int index = 0; index < nRegions; index++)
                {
                    int distance = 1000 * (1 << index);

                    int n = 0;
                    double totalASE = 0;

                    foreach (var peerVariant in annotatedSelectedVariants)
                    {
                        if (peerVariant.contig == variant.contig && Math.Abs(peerVariant.locus - variant.locus) <= distance && peerVariant != variant ||
                            index == nRegions - 1)                            // The last region is "whole autosome" and we excluded non-autosomal variants above
                        {
                            n++;
                            totalASE += getReadCounts(peerVariant).AlleleSpecificValue();
                        }
                    } // foreach peer variant

                    if (n != 0)
                    {
                        double meanASE = totalASE / n;

                        heatMapAndCount.heatMap[aseValue, index].n++;
                        heatMapAndCount.heatMap[aseValue, index].totalSamples += n;
                        heatMapAndCount.heatMap[aseValue, index].totalASE += meanASE;
                        heatMapAndCount.heatMap[aseValue, index].totalASESquared += meanASE * meanASE;
                        heatMapAndCount.heatMap[aseValue, index].maxASE = Math.Max(heatMapAndCount.heatMap[aseValue, index].maxASE, meanASE);
                        heatMapAndCount.heatMap[aseValue, index].minASE = Math.Min(heatMapAndCount.heatMap[aseValue, index].minASE, meanASE);
                    } // If we found anything
                } // foreach region
            } // foreach variant
        }

        static void WorkerThread(List<string> caseIdsToProcess, Dictionary<string, ASETools.Case> cases)
        {
            var localTumorHeatMapAndASECount = new HeatMapAndCount();
            var localNormalHeatMapAndASECount = new HeatMapAndCount();

            HeatMapEntry[,] localHeatMap = new HeatMapEntry[nASEs, nRegions];
            int[] localASECount = new int[nASEs];

            while (true)
            {
                ASETools.Case case_;
                lock (caseIdsToProcess)
                {
                    if (caseIdsToProcess.Count() == 0)
                    {
                        //
                        // Merge our local heat maps into the global ones and finish up.
                        //
                        globalTumorHeatMapAndCount.merge(localTumorHeatMapAndASECount);
                        globalNormalHeatMapAndCount.merge(localNormalHeatMapAndASECount);
 
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

                //
                // Get all the selected variants that are germline, have high enough coverage in tumor RNA and DNA, and don't show loss-of-heterozygosity in tumor DNA.
                //
                var annotatedSelectedVariants = ASETools.AnnotatedVariant.readFile(case_.annotated_selected_variants_filename);
                var nUnfiltered = annotatedSelectedVariants.Count();

                annotatedSelectedVariants = annotatedSelectedVariants.Where(x => !x.somaticMutation).ToList();
                annotatedSelectedVariants = annotatedSelectedVariants.Where(x => x.tumorRNAReadCounts.nMatchingAlt + x.tumorRNAReadCounts.nMatchingReference >= 10).ToList();
                annotatedSelectedVariants = annotatedSelectedVariants.Where(x => ASETools.isChromosomeAutosomal(x.contig)).ToList();
                annotatedSelectedVariants = annotatedSelectedVariants.Where(x => x.tumorDNAReadCounts.nMatchingAlt + x.tumorDNAReadCounts.nMatchingReference >= 10).ToList();
                annotatedSelectedVariants = annotatedSelectedVariants.Where(x => x.tumorDNAReadCounts.nMatchingAlt * 3 > x.tumorDNAReadCounts.nMatchingReference * 2).ToList();
                annotatedSelectedVariants = annotatedSelectedVariants.Where(x => x.tumorDNAReadCounts.nMatchingReference * 3 > x.tumorDNAReadCounts.nMatchingAlt * 2).ToList();

                annotatedSelectedVariants = annotatedSelectedVariants.Where(x => !x.somaticMutation && x.tumorRNAReadCounts.nMatchingAlt + x.tumorRNAReadCounts.nMatchingReference >= 10 && ASETools.isChromosomeAutosomal(x.contig) &&
                x.tumorDNAReadCounts.nMatchingAlt + x.tumorDNAReadCounts.nMatchingReference >= 10 && x.tumorDNAReadCounts.nMatchingAlt * 3 > x.tumorDNAReadCounts.nMatchingReference * 2 && x.tumorDNAReadCounts.nMatchingReference * 3 > x.tumorDNAReadCounts.nMatchingAlt * 2).ToList();

                //Console.WriteLine("Left with " + annotatedSelectedVariants.Count() + " of " + nUnfiltered);

                ProcessAnnotatedSelectedVariants(annotatedSelectedVariants, localTumorHeatMapAndASECount, x => x.tumorRNAReadCounts);
                ProcessAnnotatedSelectedVariants(annotatedSelectedVariants, localNormalHeatMapAndASECount, x => x.normalRNAReadCounts);

                if (Interlocked.Increment(ref nCasesProcessed) % 100 == 0)
                {
                    Console.Write(".");
                }
            } // while true
        }
    }
}
