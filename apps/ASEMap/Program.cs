using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using System.Threading;
using ASELib;
using System.Diagnostics;
using System.IO;

namespace ASEMap
{
    class Program
    {
        const int regionSize = 1000000;   // Break into megabase chunks.

        class MapEntry
        {
            public int nCases = 0;
            public double totalASE = 0;
            public double totalASESquared = 0;
        }

        class ASEMap
        {
            Dictionary<string, Dictionary<int, MapEntry>> map = new Dictionary<string, Dictionary<int, MapEntry>>();

            public void addASE(string chromosome, int locus, double ase)
            {
                if (!map.ContainsKey(chromosome))
                {
                    lock (map)
                    {
                        if (!map.ContainsKey(chromosome))
                        {

                            map.Add(chromosome, new Dictionary<int, MapEntry>());
                        }
                    }
                }

                lock (map[chromosome])
                {
                    int regionBase = locus - locus % regionSize;
                    if (!map[chromosome].ContainsKey(regionBase))
                    {
                        map[chromosome].Add(regionBase, new MapEntry());
                    }

                    var mapEntry = map[chromosome][regionBase];

                    mapEntry.nCases++;
                    mapEntry.totalASE += ase;
                    mapEntry.totalASESquared += ase * ase;
                }
            }

            public static string Header()
            {
                return "Chromsome\tlocus\ttumor\tn cases\tindex\tmean\tstandard deviation";
            }

            public void WriteToFile(StreamWriter outputFile, bool tumor)
            {
                foreach (var chromosomeEntry in map)
                {
                    foreach (var regionEntry in chromosomeEntry.Value)
                    {
                        var region = regionEntry.Value;

                        if (region.nCases < 100)
                        {
                            continue;   // Just skip the low-n ones to avoid noise.
                        }

                        outputFile.WriteLine(chromosomeEntry.Key + "\t" + regionEntry.Key + "\t" + tumor + "\t" + region.nCases + "\t" + ((long)(ASETools.ChromosomeNameToIndex(chromosomeEntry.Key) - 1) * 300000000 + regionEntry.Key) / regionSize + "\t" +
                            region.totalASE / region.nCases + "\t" + Math.Sqrt(region.nCases * region.totalASESquared - region.totalASE * region.totalASE) / region.nCases);
                    }
                }
            }
        }

        static ASEMap normalMap = new ASEMap();
        static ASEMap tumorMap = new ASEMap();
        static PerGeneASE globalNormalPerGeneASE = new PerGeneASE();
        static PerGeneASE globalTumorPerGeneASE = new PerGeneASE();

        static void WorkerThread(List<ASETools.Case> casesToProcess)
        {
            var thisThreadMeasurements = new List<ASEMeasurement>();
            var normalPerGeneASEThisThread = new PerGeneASE();
            var tumorPerGeneASEThisThread = new PerGeneASE();

            while (true)
            {
                ASETools.Case case_;

                lock (casesToProcess)
                {
                    if (casesToProcess.Count() == 0)
                    {
                        allMeasurements.AddRange(thisThreadMeasurements);
                        globalNormalPerGeneASE.mergeEqually(normalPerGeneASEThisThread);
                        globalTumorPerGeneASE.mergeEqually(tumorPerGeneASEThisThread);

                        return;
                    }

                    case_ = casesToProcess[0];
                    casesToProcess.RemoveAt(0);

                    nCasesProcessed++;
                    if (nCasesProcessed % 100 == 0)
                    {
                        Console.Write(".");
                    }
                }

                var normalPerGeneASEThisCase = new PerGeneASE();
                var tumorPerGeneASEThisCase = new PerGeneASE();

                var annotatedVariants = ASETools.AnnotatedVariant.readFile(case_.annotated_selected_variants_filename);
                var tumorCopyNumberVariation = ASETools.CopyNumberVariation.ReadFile(case_.tumor_copy_number_filename, case_.tumor_copy_number_file_id).Where(r => Math.Abs(r.Segment_Mean) > 1.0).ToList();
                List<ASETools.CopyNumberVariation> normalCopyNumberVariation = null;
                if (case_.normal_copy_number_filename != "")
                {
                    normalCopyNumberVariation = ASETools.CopyNumberVariation.ReadFile(case_.normal_copy_number_filename, case_.normal_copy_number_file_id).Where(r => Math.Abs(r.Segment_Mean) > 1.0).ToList();
                }
               

                if (null == annotatedVariants)
                {
                    Console.WriteLine("Unable to read annotated selected variants from " + case_.annotated_selected_variants_filename);
                    continue;
                }

                foreach (var annotatedVariant in annotatedVariants)
                {
                    if (annotatedVariant.somaticMutation)
                    {
                        continue;   // Only use germline variants for the map.
                    }

                    if (!annotatedVariant.somaticMutation && annotatedVariant.IsASECandidate(false, normalCopyNumberVariation, null, null, geneMap)) // null out the per-gene ASE, since we're where it comes from
                    {
 
                        normalMap.addASE(annotatedVariant.contig, annotatedVariant.locus, annotatedVariant.GetNormalAlleleSpecificExpression());
                        thisThreadMeasurements.Add(new ASEMeasurement(annotatedVariant.GetNormalAlleleSpecificExpression(), false));
                        normalPerGeneASEThisCase.recordSample(annotatedVariant.contig, annotatedVariant.locus, annotatedVariant.GetNormalAlleleSpecificExpression());

                    }

                    if (!annotatedVariant.somaticMutation && annotatedVariant.IsASECandidate(true, tumorCopyNumberVariation, null, null, geneMap)) // null out the per-gene ASE, since we're where it comes from
                    {
 
                        tumorMap.addASE(annotatedVariant.contig, annotatedVariant.locus, annotatedVariant.GetTumorAlleleSpecificExpression());
                        thisThreadMeasurements.Add(new ASEMeasurement(annotatedVariant.GetTumorAlleleSpecificExpression(), true));
                        tumorPerGeneASEThisCase.recordSample(annotatedVariant.contig, annotatedVariant.locus, annotatedVariant.GetTumorAlleleSpecificExpression());
                    }
                } // foreach variant

                normalPerGeneASEThisThread.mergeInOneSample(normalPerGeneASEThisCase);
                tumorPerGeneASEThisThread.mergeInOneSample(tumorPerGeneASEThisCase);
            } // while true (loop over cases)
        } // WorkerThread

        static int nCasesProcessed = 0;

        class ASEMeasurement : IComparer<ASEMeasurement>
        {
            public ASEMeasurement(double ase_, bool tumor_)
            {
                ase = ase_;
                tumor = tumor_;
            }

            public readonly double ase;
            public readonly bool tumor;

            public int Compare(ASEMeasurement a, ASEMeasurement b)
            {
                if (a.ase < b.ase) return -1;
                if (a.ase > b.ase) return 1;
                return 0;
            }
        }

        static List<ASEMeasurement> allMeasurements = new List<ASEMeasurement>();

        class PerGeneASE
        {
            class PerGeneData
            {
                public int n = 0;
                public double totalASE = 0;
                public double totalASESquared = 0;
            }

            Dictionary<string, PerGeneData> genes = new Dictionary<string, PerGeneData>();

            public void mergeEqually(PerGeneASE peer)
            {
                foreach (var geneEntry in peer.genes)
                {
                    var hugoSymbol = geneEntry.Key;
                    var perGeneData = geneEntry.Value;

                    if (!genes.ContainsKey(hugoSymbol))
                    {
                        genes.Add(hugoSymbol, new PerGeneData());
                    }

                    genes[hugoSymbol].n += perGeneData.n;
                    genes[hugoSymbol].totalASE += perGeneData.totalASE;
                    genes[hugoSymbol].totalASESquared += perGeneData.totalASESquared;
                }
            }

            //
            // Merges in a PerGeneASE that represents just one sample, so if it has multiple ASE value for a gene, just treat them as
            // a single value at the gene's mean.
            //
            public void mergeInOneSample(PerGeneASE sample)
            {
                foreach (var geneEntry in sample.genes)
                {
                    var hugoSymbol = geneEntry.Key;
                    var perGeneData = geneEntry.Value;

                    if (!genes.ContainsKey(hugoSymbol))
                    {
                        genes.Add(hugoSymbol, new PerGeneData());
                    }

                    double sampleASE = perGeneData.totalASE / perGeneData.n;
                    genes[hugoSymbol].n++;
                    genes[hugoSymbol].totalASE += sampleASE;
                    genes[hugoSymbol].totalASESquared += sampleASE * sampleASE;
                }
            }

            public void recordSample(string chromosome, int locus, double ase)
            {
                foreach (var geneInfo in geneMap.getGenesMappedTo(chromosome, locus))
                {
                    if (!genes.ContainsKey(geneInfo.hugoSymbol))
                    {
                        genes.Add(geneInfo.hugoSymbol, new PerGeneData());
                    }

                    genes[geneInfo.hugoSymbol].n++;
                    genes[geneInfo.hugoSymbol].totalASE += ase;
                    genes[geneInfo.hugoSymbol].totalASESquared += ase * ase;
                }
            }

            public void WriteToFile(StreamWriter outputFile, bool tumor, PerGeneASE peer)
            {
                int minSamplesToPrint = 10;

                foreach (var geneInfo in genes)
                {
                    if (geneInfo.Value.n < minSamplesToPrint)
                    {
                        continue;
                    }

                    var hugo_symbol = geneInfo.Key;
                    var geneLocation = geneLocationInformation.genesByName[hugo_symbol];

                    outputFile.WriteLine(ASETools.ConvertToExcelString(hugo_symbol) + "\t" + tumor + "\t" + geneInfo.Value.n + "\t" + geneInfo.Value.totalASE / geneInfo.Value.n + "\t" +
                        Math.Sqrt(geneInfo.Value.n * geneInfo.Value.totalASESquared - geneInfo.Value.totalASE * geneInfo.Value.totalASE) / geneInfo.Value.n + "\t" +
                        ((peer.genes.ContainsKey(hugo_symbol) && peer.genes[hugo_symbol].n >= minSamplesToPrint) ? Convert.ToString(geneInfo.Value.totalASE / geneInfo.Value.n - peer.genes[hugo_symbol].totalASE / peer.genes[hugo_symbol].n) : "*") + "\t" +
                        geneLocation.chromosome + "\t" + geneLocation.minLocus + "\t" + geneLocation.maxLocus
                        );
                }
            }
        }

        static ASETools.GeneMap geneMap;
        static ASETools.GeneLocationsByNameAndChromosome geneLocationInformation;

        static void Main(string[] args)
        {
            var timer = new Stopwatch();
            timer.Start();

            var configuration = ASETools.Configuration.loadFromFile(args);

            if (null == configuration)
            {
                Console.WriteLine("Giving up because we were unable to load configuration.");
                return;
            }

            geneLocationInformation = new ASETools.GeneLocationsByNameAndChromosome(ASETools.readKnownGeneFile(configuration.geneLocationInformationFilename));
            geneMap = new ASETools.GeneMap(geneLocationInformation.genesByName);

            var cases = ASETools.Case.LoadCases(configuration.casesFilePathname);

            if (null == cases)
            {
                Console.WriteLine("Unable to load cases.");
            }

            var casesToProcess = cases.Where(x => x.Value.annotated_selected_variants_filename != "" && x.Value.tumor_copy_number_filename != "").Select(x => x.Value).ToList();
            Console.Write("Processing " + casesToProcess.Count() + " cases, 1 dot/100 cases: ");

            var threads = new List<Thread>();
            for (int i = 0; i < Environment.ProcessorCount; i++)
            {
                threads.Add(new Thread(() => WorkerThread(casesToProcess)));
            }

            threads.ForEach(t => t.Start());
            threads.ForEach(t => t.Join());

            var outputFile = ASETools.CreateStreamWriterWithRetry(configuration.finalResultsDirectory + ASETools.ASEMapFilename);
            outputFile.WriteLine(ASEMap.Header());

            tumorMap.WriteToFile(outputFile, true);
            normalMap.WriteToFile(outputFile, false);

            outputFile.WriteLine("**done**");
            outputFile.Close();

            var perGeneOutputFile = ASETools.CreateStreamWriterWithRetry(configuration.finalResultsDirectory + ASETools.PerGeneASEMapFilename);
            perGeneOutputFile.WriteLine("Hugo Symbol\tTumor\tn Samples\tmean ASE\tstandard deviation of ASE\tASE difference with opposite tumor/normal\tchromosome\tmin locus\tmax locus");

            globalTumorPerGeneASE.WriteToFile(perGeneOutputFile, true, globalNormalPerGeneASE);
            globalNormalPerGeneASE.WriteToFile(perGeneOutputFile, false, globalTumorPerGeneASE);
            perGeneOutputFile.WriteLine("**done**");
            perGeneOutputFile.Close();

            Console.WriteLine();

            if (false) // For some reason, it takes more than a day to run Mann-Whitney on this, and then the p value is 0.  So, skip it.
            {
                bool enoughData;
                bool reversed;
                double nFirstGroup;
                double nSecondGroup;
                double U;
                double z;

                double p = ASETools.MannWhitney<ASEMeasurement>.ComputeMannWhitney(allMeasurements, allMeasurements[0], x => x.tumor, x => x.ase, out enoughData, out reversed, out nFirstGroup, out nSecondGroup, out U,
                    out z);

                Console.WriteLine((int)(nFirstGroup + nSecondGroup) + " total ASE measurements between tumor and normal.  Distributions differ with p = " + p + ".  Mean for tumor is " +
                    allMeasurements.Where(x => x.tumor).Select(x => x.ase).Average() + " mean for normal is " + allMeasurements.Where(x => !x.tumor).Select(x => x.ase).Average());
            }
            Console.Write("Took " + ASETools.ElapsedTimeInSeconds(timer));

        }
    }
}
