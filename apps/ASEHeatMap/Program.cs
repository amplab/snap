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
        const int nASEs = 21;      // 0 - 1 inclusive by .05
        const int nRegions = 20;    // Powers of two * 1Kb starting with 1 KB (ie., no 0) and extending to 256MB (insclusive), plus one for rest of autosome.  This is selected because 256MB is bigger than chr1.
        class HeatMapEntry
        {
            public double totalASE = 0;
            public double totalASESquared = 0;
            public int n = 0;
            public int totalSamples = 0;
            public double maxASE = -1;
            public double minASE = 2;
            public ASETools.Histogram histogram = new ASETools.Histogram();

            public void merge(HeatMapEntry peer)
            {
                totalASE += peer.totalASE;
                totalASESquared += peer.totalASESquared;
                n += peer.n;
                totalSamples += peer.totalSamples;
                maxASE = Math.Max(maxASE, peer.maxASE);
                minASE = Math.Min(minASE, peer.minASE);
                histogram.merge(peer.histogram);
            }
        }

        class HeatMapsAndCount
        {
            public HeatMapsAndCount()
            {
                for (int aseIndex = 0; aseIndex < nASEs; aseIndex++)
                {
                    aseCount[aseIndex] = 0;

                    for (int i = 0; i < nRegions; i++)
                    {
                        heatMap[aseIndex, i] = new HeatMapEntry();
                        exclusiveHeatMap[aseIndex, i] = new HeatMapEntry();
                        exclusiveAndExclusiveOfGeneHeatMap[aseIndex, i] = new HeatMapEntry();
                        selectedGenes[aseIndex, i] = new HeatMapEntry();
                    }
                }
            }

            public void merge(HeatMapsAndCount peer)
            {
                for (int aseIndex = 0; aseIndex < nASEs; aseIndex++)
                {
                    aseCount[aseIndex] += peer.aseCount[aseIndex];

                    for (int i = 0; i < nRegions; i++)
                    {
                        heatMap[aseIndex, i].merge(peer.heatMap[aseIndex, i]);
                        exclusiveHeatMap[aseIndex, i].merge(peer.exclusiveHeatMap[aseIndex, i]);
                        exclusiveAndExclusiveOfGeneHeatMap[aseIndex, i].merge(peer.exclusiveAndExclusiveOfGeneHeatMap[aseIndex, i]);
                        selectedGenes[aseIndex, i].merge(peer.selectedGenes[aseIndex, i]);
                    }
                }
            }

            void writeColumnHeaders(StreamWriter outputFile, string exclusivelyString, bool exclusive)
            {
                for (int i = 0; i < nRegions; i++)
                {
                    outputFile.Write("\tMean ASE " + exclusivelyString + "Within " + distanceIndexToString(i, exclusive));
                }

                for (int i = 0; i < nRegions; i++)
                {
                    outputFile.Write("\tNumber of ASE samples that have anything " + exclusivelyString + "within " + distanceIndexToString(i, exclusive));
                }

                for (int i = 0; i < nRegions; i++)
                {
                    outputFile.Write("\tMean number of ASE samples " + exclusivelyString + "within " + distanceIndexToString(i, exclusive));
                }

                for (int i = 0; i < nRegions; i++)
                {
                    outputFile.Write("\t Standard Deviation of ASE " + exclusivelyString + "within " + distanceIndexToString(i, exclusive));
                }
            }

            void writeColumns(StreamWriter outputFile, int aseIndex, HeatMapEntry[,] heatMapToUse)
            {
                for (int i = 0; i < nRegions; i++)
                {
                    if (heatMapToUse[aseIndex, i].n < configuration.minSampleCountForHeatMap)
                    {
                        outputFile.Write("\t"); // Just leave these blank
                    }
                    else
                    {
                        outputFile.Write("\t" + heatMapToUse[aseIndex, i].totalASE / heatMapToUse[aseIndex, i].n);
                    }
                }

                for (int i = 0; i < nRegions; i++)
                {
                    outputFile.Write("\t" + heatMapToUse[aseIndex, i].n);
                }

                for (int i = 0; i < nRegions; i++)
                {
                    if (heatMapToUse[aseIndex, i].n == 0)
                    {
                        outputFile.Write("\t0");
                    }
                    else
                    {
                        outputFile.Write("\t" + (double)heatMapToUse[aseIndex, i].totalSamples / heatMapToUse[aseIndex, i].n);
                    }
                }

                for (int i = 0; i < nRegions; i++)
                {
                    if (heatMapToUse[aseIndex, i].n < configuration.minSampleCountForHeatMap)
                    {
                        outputFile.Write("\t");
                    }
                    else
                    {
                        outputFile.Write("\t" + Math.Sqrt(heatMapToUse[aseIndex, i].n * heatMapToUse[aseIndex, i].totalASE - heatMapToUse[aseIndex, i].totalASESquared) / heatMapToUse[aseIndex, i].n);
                    }
                }
            }

            public void writeToFile(StreamWriter outputFile, bool writeHeader)
            {

                if (writeHeader)
                {
                    outputFile.Write("Allele Specific Expression\tCount of samples at this ASE\tpdf of ASE count");
                    writeColumnHeaders(outputFile, "Near selected genes ", true);
                    writeColumnHeaders(outputFile, "exclusive and exlcusive of gene ", true);
                    writeColumnHeaders(outputFile, "exclusive ", true);
                    writeColumnHeaders(outputFile, "", false);
                }
 
                outputFile.WriteLine();


                int totalSamples = aseCount.Sum();

                for (int aseIndex = 0; aseIndex < nASEs; aseIndex++)
                {
                    outputFile.Write((double)aseIndex * 1 / (nASEs - 1) + "\t" + aseCount[aseIndex] + "\t" + (double)aseCount[aseIndex] / totalSamples);

                    writeColumns(outputFile, aseIndex, selectedGenes);
                    writeColumns(outputFile, aseIndex, exclusiveAndExclusiveOfGeneHeatMap);
                    writeColumns(outputFile, aseIndex, exclusiveHeatMap);
                    writeColumns(outputFile, aseIndex, heatMap);

                    outputFile.WriteLine();
                } // foreach ASE index
            }

            public void writeHistogramsToFile(string filename)
            {
                var outputFile = ASETools.CreateStreamWriterWithRetry(filename);

                if (null == outputFile)
                {
                    Console.WriteLine("Error opening histogram output file.");
                    return;
                }

                outputFile.WriteLine("Class\tASE of variant site\tRange\tASE of histogram bucket\tcount\ttotal\tpdf\tcdf");

                dumpClassHistogramToFile(outputFile, selectedGenes, "Selected genes exclusive of gene", true);
                dumpClassHistogramToFile(outputFile, exclusiveAndExclusiveOfGeneHeatMap, "Exclusive and exclusive of gene", true);
                dumpClassHistogramToFile(outputFile, exclusiveAndExclusiveOfGeneHeatMap, "Exclusive", true);
                dumpClassHistogramToFile(outputFile, exclusiveAndExclusiveOfGeneHeatMap, "Inclusive", false);
            }

            void dumpClassHistogramToFile(StreamWriter outputFile, HeatMapEntry[,] heatMap, string className, bool exclusive) 
            {
                for (int aseIndex = 0; aseIndex < nASEs; aseIndex++)
                {
                    for (int regionIndex = 0; regionIndex < nRegions; regionIndex++)
                    {
                        var histogramLines = heatMap[aseIndex, regionIndex].histogram.ComputeHistogram(0, 1, (double)1 / (nASEs-1));
                        foreach (var line in histogramLines)
                        {
                            outputFile.WriteLine(className + "\t" + aseIndex * (double)1 / (nASEs - 1) + "\t" + distanceIndexToString(regionIndex, exclusive) + "\t" + line.minValue + "\t" + line.count + "\t" + line.total + "\t" + line.pdfValue + "\t" + line.cdfValue);
                        } // foreach histogram line
                    } // foreach region
                } // foreach ase
            }


            public HeatMapEntry[,] heatMap = new HeatMapEntry[nASEs, nRegions];
            public HeatMapEntry[,] exclusiveHeatMap = new HeatMapEntry[nASEs, nRegions];
            public HeatMapEntry[,] exclusiveAndExclusiveOfGeneHeatMap = new HeatMapEntry[nASEs, nRegions];
            public HeatMapEntry[,] selectedGenes = new HeatMapEntry[nASEs, nRegions];
            public int[] aseCount = new int[nASEs];
        }

        static HeatMapsAndCount globalTumorHeatMapAndCount = new HeatMapsAndCount();
        static HeatMapsAndCount globalNormalHeatMapAndCount = new HeatMapsAndCount();
        static ASETools.GeneMap geneMap = null;
        static ASETools.Configuration configuration;

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

            var geneLocationInformation = new ASETools.GeneLocationsByNameAndChromosome(ASETools.readKnownGeneFile(configuration.geneLocationInformationFilename));
            geneMap = new ASETools.GeneMap(geneLocationInformation.genesByName);

            var allSignificantResults = ASETools.AllSignificantResultsLine.ReadFile(configuration.finalResultsDirectory + ASETools.AllSignificantResultsFilename);
            if (null == allSignificantResults)
            {
                Console.WriteLine("You must first have an AllSignificantResults.txt file, which is generated by ApplyBonferroniCorrection");
                return;
            }

            var selectedGeneNames = allSignificantResults.Where(x => x.inputFile == ASETools.AllSignificantResultsFilename + ".txt" && x.range_index > 5).Select(x => x.hugo_symbol).Distinct().ToArray();

            var selectedGenes = new List<ASETools.GeneLocationInfo>();
            foreach (var geneName in selectedGeneNames)
            {
                selectedGenes.Add(geneLocationInformation.genesByName[geneName]);
            }

            Console.Write("Progress (1 dot/100 cases): ");

            var threads = new List<Thread>();
            for (int i = 0; i < Environment.ProcessorCount ; i++)
            {
                threads.Add(new Thread(() => WorkerThread(caseIdsToProcess, cases, selectedGenes)));
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

            globalTumorHeatMapAndCount.writeToFile(outputFile, true);

            outputFile.WriteLine();
            outputFile.WriteLine("Matched normal");

            globalNormalHeatMapAndCount.writeToFile(outputFile, false);

            outputFile.WriteLine("**done**");
            outputFile.Close();

            globalTumorHeatMapAndCount.writeHistogramsToFile(configuration.finalResultsDirectory + ASETools.tumorHeatMapHistogramFilename);
            globalNormalHeatMapAndCount.writeHistogramsToFile(configuration.finalResultsDirectory + ASETools.normalHeatMapHistogramFilename);

            Console.WriteLine("Processed " + cases.Count() + " cases in " + ASETools.ElapsedTimeInSeconds(timer));
        } // Main

        static string distanceIndexToString(int distanceIndex, bool exclusive)
        {
            if (exclusive)
            {
                switch (distanceIndex)
                {
                    case 0: return "<=1KB";
                    case 1: return "1KB-2KB";
                    case 2: return "2KB-4KB";
                    case 3: return "4KB-8KB";
                    case 4: return "8KB-16KB";
                    case 5: return "16KB-32KB";
                    case 6: return "32KB-64KB";
                    case 7: return "64KB-128KB";
                    case 8: return "128KB-256KB";
                    case 9: return "256KB-512KB";
                    case 10: return "512KB-1MB";
                    case 11: return "1MB-2MB";
                    case 12: return "2MB-4MB";
                    case 13: return "4MB-8MB";
                    case 14: return "8MB-16MB";
                    case 15: return "16MB-32MB";
                    case 16: return "32MB-64MB";
                    case 17: return "64MB-128MB";
                    case 18: return "chromosome";
                    case 19: return "autosome";
                }
            }

            switch (distanceIndex)
            {
                case 0: return "<=1KB";
                case 1: return "<=2KB";
                case 2: return "<=4KB";
                case 3: return "<=8KB";
                case 4: return "<=16KB";
                case 5: return "<=32KB";
                case 6: return "<=64KB";
                case 7: return "<=128KB";
                case 8: return "<=256KB";
                case 9: return "<=512KB";
                case 10: return "<=1MB";
                case 11: return "<=2MB";
                case 12: return "<=4MB";
                case 13: return "<=8MB";
                case 14: return "<=16MB";
                case 15: return "<=32MB";
                case 16: return "<=64MB";
                case 17: return "<=128MB";
                case 18: return "chromosome";
                case 19: return "autosome";
            }

            return "Illegal index??";
        }

        static int nCasesProcessed = 0;

        delegate ASETools.ReadCounts GetReadCounts(ASETools.AnnotatedVariant variant);

        static void ProcessAnnotatedSelectedVariants(List<ASETools.AnnotatedVariant> annotatedSelectedVariants, HeatMapsAndCount heatMapAndCount, GetReadCounts getReadCounts, List<ASETools.GeneLocationInfo> selectedGenes)
        {
            if (annotatedSelectedVariants.Count() == 0 || getReadCounts(annotatedSelectedVariants[0]) == null)
            {
                return;
            }

            foreach (var variant in annotatedSelectedVariants)
            {
                int aseIndex = (int)(getReadCounts(variant).AlleleSpecificValue() * (nASEs -1) + .5); // +.5 is for rounding.

                heatMapAndCount.aseCount[aseIndex]++;

                int minEncompassingLocus, maxEncompassingLocus;
                bool inAGene = geneMap.encompassingGeneRange(variant.contig, variant.locus, out minEncompassingLocus, out maxEncompassingLocus);   // min and max get filled in as -1 if the variant isn't in a gene

                bool inACoveredGene = false;
                selectedGenes.ForEach(x => inACoveredGene |= x.containsLocus(variant.contig, variant.locus));

                for (int index = 0; index < nRegions; index++)
                {
                    int distance = 1000 * (1 << index);
                    int exclusiveMinimumDistance = (0 == index) ? 0 : (distance / 2);

                    int n = 0;
                    double totalASE = 0;

                    int nExclusive = 0;
                    double totalASEExclusive = 0;

                    int nExclusiveAndExclusiveOfGene = 0;
                    double totalASEExclusiveAndExclusiveOfGene = 0;

                    foreach (var peerVariant in annotatedSelectedVariants)
                    {
                        if (peerVariant.contig == variant.contig && (Math.Abs(peerVariant.locus - variant.locus) <= distance && peerVariant != variant || index == nRegions - 2) || // nRegions-2 is whole chromosome
                            index == nRegions - 1)                            // The last two regions are "whole chromosome" and "whole autosome" and we excluded non-autosomal variants above
                        {
                            double ase = getReadCounts(peerVariant).AlleleSpecificValue();

                            n++;
                            totalASE += ase;

                            if (Math.Abs(peerVariant.locus - variant.locus) >= exclusiveMinimumDistance || index >= nRegions - 2)
                            {
                                nExclusive++;
                                totalASEExclusive += ase;

                                if (peerVariant.locus < minEncompassingLocus || peerVariant.locus > maxEncompassingLocus)
                                {
                                    nExclusiveAndExclusiveOfGene++;
                                    totalASEExclusiveAndExclusiveOfGene += ase;
                                }
                            }
                        }
                    } // foreach peer variant

                    if (n != 0)
                    {
                        double meanASE = totalASE / n;

                        heatMapAndCount.heatMap[aseIndex, index].n++;
                        heatMapAndCount.heatMap[aseIndex, index].totalSamples += n;
                        heatMapAndCount.heatMap[aseIndex, index].totalASE += meanASE;
                        heatMapAndCount.heatMap[aseIndex, index].totalASESquared += meanASE * meanASE;
                        heatMapAndCount.heatMap[aseIndex, index].maxASE = Math.Max(heatMapAndCount.heatMap[aseIndex, index].maxASE, meanASE);
                        heatMapAndCount.heatMap[aseIndex, index].minASE = Math.Min(heatMapAndCount.heatMap[aseIndex, index].minASE, meanASE);
                        heatMapAndCount.heatMap[aseIndex, index].histogram.addValue(meanASE);

                        if (nExclusive > 0)
                        {
                            double meanExclusiveASE = totalASEExclusive / nExclusive;

                            heatMapAndCount.exclusiveHeatMap[aseIndex, index].n++;
                            heatMapAndCount.exclusiveHeatMap[aseIndex, index].totalSamples += nExclusive;
                            heatMapAndCount.exclusiveHeatMap[aseIndex, index].totalASE += meanExclusiveASE;
                            heatMapAndCount.exclusiveHeatMap[aseIndex, index].totalASESquared += meanExclusiveASE * meanExclusiveASE;
                            heatMapAndCount.exclusiveHeatMap[aseIndex, index].maxASE = Math.Max(heatMapAndCount.exclusiveHeatMap[aseIndex, index].maxASE, meanExclusiveASE);
                            heatMapAndCount.exclusiveHeatMap[aseIndex, index].minASE = Math.Min(heatMapAndCount.exclusiveHeatMap[aseIndex, index].minASE, meanExclusiveASE);
                            heatMapAndCount.exclusiveHeatMap[aseIndex, index].histogram.addValue(meanExclusiveASE);

                            if (nExclusiveAndExclusiveOfGene > 0 && inAGene)
                            {
                                double meanExclusiveAndExclusiveOfGeneASE = totalASEExclusiveAndExclusiveOfGene / nExclusiveAndExclusiveOfGene;

                                heatMapAndCount.exclusiveAndExclusiveOfGeneHeatMap[aseIndex, index].n++;
                                heatMapAndCount.exclusiveAndExclusiveOfGeneHeatMap[aseIndex, index].totalSamples += nExclusiveAndExclusiveOfGene;
                                heatMapAndCount.exclusiveAndExclusiveOfGeneHeatMap[aseIndex, index].totalASE += meanExclusiveAndExclusiveOfGeneASE;
                                heatMapAndCount.exclusiveAndExclusiveOfGeneHeatMap[aseIndex, index].totalASESquared += meanExclusiveAndExclusiveOfGeneASE * meanExclusiveAndExclusiveOfGeneASE;
                                heatMapAndCount.exclusiveAndExclusiveOfGeneHeatMap[aseIndex, index].maxASE = Math.Max(heatMapAndCount.exclusiveHeatMap[aseIndex, index].maxASE, meanExclusiveAndExclusiveOfGeneASE);
                                heatMapAndCount.exclusiveAndExclusiveOfGeneHeatMap[aseIndex, index].minASE = Math.Min(heatMapAndCount.exclusiveHeatMap[aseIndex, index].minASE, meanExclusiveAndExclusiveOfGeneASE);
                                heatMapAndCount.exclusiveAndExclusiveOfGeneHeatMap[aseIndex, index].histogram.addValue(meanExclusiveAndExclusiveOfGeneASE);

                                if (inACoveredGene)
                                {
                                    heatMapAndCount.selectedGenes[aseIndex, index].n++;
                                    heatMapAndCount.selectedGenes[aseIndex, index].totalSamples += nExclusiveAndExclusiveOfGene;
                                    heatMapAndCount.selectedGenes[aseIndex, index].totalASE += meanExclusiveAndExclusiveOfGeneASE;
                                    heatMapAndCount.selectedGenes[aseIndex, index].totalASESquared += meanExclusiveAndExclusiveOfGeneASE * meanExclusiveAndExclusiveOfGeneASE;
                                    heatMapAndCount.selectedGenes[aseIndex, index].maxASE = Math.Max(heatMapAndCount.exclusiveHeatMap[aseIndex, index].maxASE, meanExclusiveAndExclusiveOfGeneASE);
                                    heatMapAndCount.selectedGenes[aseIndex, index].minASE = Math.Min(heatMapAndCount.exclusiveHeatMap[aseIndex, index].minASE, meanExclusiveAndExclusiveOfGeneASE);
                                    heatMapAndCount.selectedGenes[aseIndex, index].histogram.addValue(meanExclusiveAndExclusiveOfGeneASE);
                                }
                            }
                        }
                    } // If we found anything
                } // foreach region
            } // foreach variant
        }

        static void WorkerThread(List<string> caseIdsToProcess, Dictionary<string, ASETools.Case> cases, List<ASETools.GeneLocationInfo> selectedGenes)
        {
            var localTumorHeatMapAndASECount = new HeatMapsAndCount();
            var localNormalHeatMapAndASECount = new HeatMapsAndCount();

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

                annotatedSelectedVariants = annotatedSelectedVariants.Where(x => !x.somaticMutation && 
				ASETools.isChromosomeAutosomal(x.contig) && x.IsASECandidate()).ToList();

                //Console.WriteLine("Left with " + annotatedSelectedVariants.Count() + " of " + nUnfiltered);

                ProcessAnnotatedSelectedVariants(annotatedSelectedVariants, localTumorHeatMapAndASECount, x => x.tumorRNAReadCounts, selectedGenes);
                ProcessAnnotatedSelectedVariants(annotatedSelectedVariants, localNormalHeatMapAndASECount, x => x.normalRNAReadCounts, selectedGenes);

                if (Interlocked.Increment(ref nCasesProcessed) % 100 == 0)
                {
                    Console.Write(".");
                }
            } // while true
        }
    }
}
