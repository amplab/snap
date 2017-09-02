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

                        outputFile.WriteLine(chromosomeEntry.Key + "\t" + regionEntry.Key +  "\t" + tumor + "\t" + region.nCases + "\t" + ((long)(ASETools.ChromosomeNameToIndex(chromosomeEntry.Key) - 1) * 300000000 + regionEntry.Key)/ regionSize + "\t" +
                            region.totalASE / region.nCases + "\t" + Math.Sqrt(region.nCases * region.totalASESquared - region.totalASE * region.totalASE) / region.nCases);
                    }
                }
            }
        }

        static ASEMap normalMap = new ASEMap();
        static ASEMap tumorMap = new ASEMap();

        static void WorkerThread(List<ASETools.Case> casesToProcess)
        {
            while (true)
            {
                ASETools.Case case_;

                lock (casesToProcess)
                {
                    if (casesToProcess.Count() == 0)
                    {
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

                var annotatedVariants = ASETools.AnnotatedVariant.readFile(case_.annotated_selected_variants_filename);
                var copyNumberVariation = ASETools.CopyNumberVariation.ReadFile(case_.tumor_copy_number_filename, case_.tumor_copy_number_file_id).Where(r => Math.Abs(r.Segment_Mean) > 1.0).ToList();
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

                    if (annotatedVariant.normalRNAReadCounts != null && annotatedVariant.IsASECandidate(false))
                    {
                        bool overlapsCNV = false;
                        if (normalCopyNumberVariation != null)
                        {
                            var overlappingCNV = copyNumberVariation.Where(cnv =>
                                    cnv.OverlapsLocus(ASETools.chromosomeNameToNonChrForm(annotatedVariant.contig),
                                    annotatedVariant.locus, annotatedVariant.locus + 1));
                            overlapsCNV = overlappingCNV.Count() > 0;
                        }

                        if (!overlapsCNV)
                        {
                            normalMap.addASE(annotatedVariant.contig, annotatedVariant.locus, annotatedVariant.GetNormalAlleleSpecificExpression());
                        }
                    }

                    if (annotatedVariant.IsASECandidate())
                    {
                        var overlappingCNV = copyNumberVariation.Where(cnv =>
                            cnv.OverlapsLocus(ASETools.chromosomeNameToNonChrForm(annotatedVariant.contig),
                            annotatedVariant.locus, annotatedVariant.locus + 1));
                        var overlapsCNV = overlappingCNV.Count() > 0;

                        if (!overlapsCNV)
                        {
                            tumorMap.addASE(annotatedVariant.contig, annotatedVariant.locus, annotatedVariant.GetTumorAlleleSpecificExpression());
                        }
                   }
                }
            }
        }

        static int nCasesProcessed = 0;

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

            Console.WriteLine();
            Console.Write("Took " + ASETools.ElapsedTimeInSeconds(timer));

        }
    }
}
