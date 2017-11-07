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
                double tumorMappedBaseCount = (double)(ASETools.MappedBaseCount.readFromFile(case_.tumor_rna_mapped_base_count_filename).mappedBaseCount);
                double normalMappedBaseCount = 0;

                var copyNumber = ASETools.CopyNumberVariation.ReadBothFiles(case_);

                if (case_.normal_rna_mapped_base_count_filename != "")
                {
                    normalMappedBaseCount = (double)(ASETools.MappedBaseCount.readFromFile(case_.normal_rna_mapped_base_count_filename).mappedBaseCount);
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

                    if (!annotatedVariant.somaticMutation && annotatedVariant.IsASECandidate(false, copyNumber, null, null, geneMap)) // null out the per-gene ASE, since we're where it comes from
                    {
                        normalMap.addASE(annotatedVariant.contig, annotatedVariant.locus, annotatedVariant.GetNormalAlleleSpecificExpression());
                        thisThreadMeasurements.Add(new ASEMeasurement(annotatedVariant.GetNormalAlleleSpecificExpression(), false));
                        normalPerGeneASEThisCase.recordSample(annotatedVariant.contig, annotatedVariant.locus, annotatedVariant.GetNormalAlleleSpecificExpression(), annotatedVariant.normalRNAReadCounts.totalReads() / normalMappedBaseCount);
                    }

                    if (!annotatedVariant.somaticMutation && annotatedVariant.IsASECandidate(true, copyNumber, null, null, geneMap)) // null out the per-gene ASE, since we're where it comes from
                    {
 
                        tumorMap.addASE(annotatedVariant.contig, annotatedVariant.locus, annotatedVariant.GetTumorAlleleSpecificExpression());
                        thisThreadMeasurements.Add(new ASEMeasurement(annotatedVariant.GetTumorAlleleSpecificExpression(), true));
                        tumorPerGeneASEThisCase.recordSample(annotatedVariant.contig, annotatedVariant.locus, annotatedVariant.GetTumorAlleleSpecificExpression(), annotatedVariant.tumorRNAReadCounts.totalReads() / tumorMappedBaseCount);
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
                public double totalFractionOfRNAReads = 0;
                public double totalFractionOfRNAReadsSquared = 0;

                public List<double> allASE = new List<double>();
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
                    genes[hugoSymbol].totalFractionOfRNAReads += perGeneData.totalFractionOfRNAReads;
                    genes[hugoSymbol].totalFractionOfRNAReadsSquared += perGeneData.totalFractionOfRNAReadsSquared;
                    genes[hugoSymbol].allASE.AddRange(perGeneData.allASE);
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
                    double totalFractionOfRNAReads = perGeneData.totalFractionOfRNAReads / perGeneData.n;
                    genes[hugoSymbol].totalFractionOfRNAReads += totalFractionOfRNAReads;
                    genes[hugoSymbol].totalFractionOfRNAReadsSquared += totalFractionOfRNAReads * totalFractionOfRNAReads;
                    genes[hugoSymbol].allASE.Add(sampleASE);
                }
            }

            public void recordSample(string chromosome, int locus, double ase, double RNAFraction)
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
                    genes[geneInfo.hugoSymbol].totalFractionOfRNAReads += RNAFraction;
                    genes[geneInfo.hugoSymbol].totalFractionOfRNAReadsSquared += RNAFraction * RNAFraction;
                    genes[geneInfo.hugoSymbol].allASE.Add(ase);
                }
            }


            class ASEValueAndTumor : IComparer<ASEValueAndTumor>
            {
                public ASEValueAndTumor(bool tumor_, double ase_)
                {
                    tumor = tumor_;
                    ase = ase_;
                }

                public readonly bool tumor;
                public readonly double ase;

                public int Compare(ASEValueAndTumor a, ASEValueAndTumor b)
                {
                    return a.ase.CompareTo(b.ase);
                }
            }
            public void WriteToFile(StreamWriter outputFile, PerGeneASE peer)
            {
                int minSamplesToPrint = 10;


                var mwPerGene = new Dictionary<string, double>();

                ASETools.MannWhitney<ASEValueAndTumor>.GetValue getValue = new ASETools.MannWhitney<ASEValueAndTumor>.GetValue(m => m.ase);
                ASETools.MannWhitney<ASEValueAndTumor>.WhichGroup whichGroup = new ASETools.MannWhitney<ASEValueAndTumor>.WhichGroup(m => m.tumor);

                foreach (var geneInfo in genes)
                {
                    var hugo_symbol = geneInfo.Key;
                    if (geneInfo.Value.n < minSamplesToPrint || !peer.genes.ContainsKey(hugo_symbol) || peer.genes[hugo_symbol].n <= minSamplesToPrint)
                    {
                        continue;
                    }

                    var values = new List<ASEValueAndTumor>();
                    foreach (var value in geneInfo.Value.allASE)
                    {
                        values.Add(new ASEValueAndTumor(true, value));
                    }

                    foreach (var value in peer.genes[hugo_symbol].allASE)
                    {
                        values.Add(new ASEValueAndTumor(false, value));
                    }

                    bool enoughData;
                    bool reversed;
                    double nFirstGroup;
                    double nSecondGroup;
                    double U;
                    double z;

                    var mw = ASETools.MannWhitney<ASEValueAndTumor>.ComputeMannWhitney(values, values[0], whichGroup, getValue,
                        out enoughData, out reversed, out nFirstGroup, out nSecondGroup, out U, out z);

                    if (enoughData)
                    {
                        mwPerGene.Add(hugo_symbol, mw);
                    }
                }

                int nMW = mwPerGene.Count();

                foreach (var mwEntry in mwPerGene)
                {

                    var hugo_symbol = mwEntry.Key;
                    var geneInfo = genes[hugo_symbol];

                    var geneLocation = geneLocationInformation.genesByName[hugo_symbol];

                    var peerGeneInfo = peer.genes[hugo_symbol];

                    outputFile.WriteLine(ASETools.ConvertToExcelString(hugo_symbol) + "\t" + geneInfo.n + "\t" + geneInfo.totalASE / geneInfo.n + "\t" + geneInfo.totalFractionOfRNAReads / geneInfo.n + "\t" +
                        Math.Sqrt(geneInfo.n * geneInfo.totalASESquared - geneInfo.totalASE * geneInfo.totalASE) / geneInfo.n + "\t" +
                        peerGeneInfo.n + "\t" + peerGeneInfo.totalASE / peerGeneInfo.n + "\t" + peerGeneInfo.totalFractionOfRNAReadsSquared / peerGeneInfo.n + "\t" +
                        Math.Sqrt(peerGeneInfo.n * peerGeneInfo.totalASESquared - peerGeneInfo.totalASE * peerGeneInfo.totalASE) / peerGeneInfo.n + "\t" +
                        (geneInfo.totalASE / geneInfo.n - peer.genes[hugo_symbol].totalASE / peer.genes[hugo_symbol].n) + "\t" + Math.Abs((geneInfo.totalASE / geneInfo.n - peer.genes[hugo_symbol].totalASE / peer.genes[hugo_symbol].n)) + "\t" +
                        mwEntry.Value * nMW + "\t" + 
                        geneLocation.chromosome + "\t" + geneLocation.minLocus + "\t" + geneLocation.maxLocus
                        );
                }
            }
        } // PerGeneASE

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
            perGeneOutputFile.WriteLine("Hugo Symbol\tn Tumor Samples\tmean Tumor ASE\tTumor Fraction of RNA at Variant Sites\tstandard deviation of Tumor ASE\tn Normal Samples\tmean Normal ASE\tNormal Fraction of RNA at Variant Sites\tstandard deviation of Normal ASE\t" +
                "Tumor ASE minus Normal ASE\tAbsolute value of Tumor ASE minus Normal ASE\tBonferroni Corrected MW p value for normal differing from tumor\tchromosome\tmin locus\tmax locus");

            globalTumorPerGeneASE.WriteToFile(perGeneOutputFile, globalNormalPerGeneASE);

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
