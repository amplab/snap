using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using System.Threading;
using ASELib;
using System.Diagnostics;

namespace FindASERuns
{
    class Program
    {
        static ASETools.Configuration configuration;
        static ASETools.GeneMap geneMap;
        static Dictionary<bool, Dictionary<string, ASETools.ASEMapPerGeneLine>> perGeneASEMap;

        static double minRangeASE = 0.3;
        static double ASERangeIncrement = 0.05;

        static int nRanges = (int)((1 - minRangeASE) / ASERangeIncrement);

        static int nCasesProcessed = 0;

        static void Main(string[] args)
        {
            var timer = new Stopwatch();
            timer.Start();

            configuration = ASETools.Configuration.loadFromFile(args);

            if (null == configuration)
            {
                Console.WriteLine("Giving up because we were unable to load configuration.");
                return;
            }

            var geneLocationInformation = new ASETools.GeneLocationsByNameAndChromosome(ASETools.readKnownGeneFile(configuration.geneLocationInformationFilename));
            geneMap = new ASETools.GeneMap(geneLocationInformation.genesByName);
            perGeneASEMap = ASETools.ASEMapPerGeneLine.ReadFromFileToDictionary(configuration.finalResultsDirectory + ASETools.PerGeneASEMapFilename);

            var cases = ASETools.Case.LoadCases(configuration.casesFilePathname);

            if (null == cases)
            {
                Console.WriteLine("Unable to load cases.");
            }

            var casesToProcess = cases.Select(x => x.Value).Where(x => x.annotated_selected_variants_filename != "" && x.tumor_copy_number_filename != "").ToList();

            Console.Write("Processing " + casesToProcess.Count() + " cases, 1 dot/1K cases: ");

            var threads = new List<Thread>();
            for (int i = 0; i < Environment.ProcessorCount; i++)
            {
                threads.Add(new Thread(() => WorkerThread(casesToProcess)));
            }

            threads.ForEach(t => t.Start());
            threads.ForEach(t => t.Join());
        } // Main

        class ASERun : IComparable
        {
            public ASERun(double minASEAllowed_, string chromosome_, int startLocus_)
            {
                minASEAllowed = minASEAllowed_;
                chromosome = chromosome_;
                startLocus = startLocus_;
                endLocus = startLocus;
            }

            public ASERun(double minASEAllowed_)
            {
                minASEAllowed = minASEAllowed_;
                chromosome = "Not a chromosome";
                endLocus = startLocus = 0;
            }

            public void recordObservation(int locus, double ase)
            {
                minASE = Math.Min(ase, minASE);
                maxASE = Math.Max(ase, maxASE);
                totalASE += ase;
                nLoci++;
                endLocus = locus;
            }

            public int CompareTo(object peerObject)
            {
                var peer = (ASERun)peerObject;

                if (nLoci < 10 || peer.nLoci < 10) return nLoci.CompareTo(peer.nLoci);

                return (endLocus - startLocus).CompareTo(peer.endLocus - peer.startLocus);
            }

            public readonly double minASEAllowed;
            public double minASE = 2;
            public double maxASE = -1;
            public double totalASE = 0;
            public readonly string chromosome;
            public readonly int startLocus;
            public int endLocus;
            public int nLoci = 0;
        }

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
                } // lock casesToProcess

                var tumorCopyNumberVariation = ASETools.CopyNumberVariation.ReadFile(case_.tumor_copy_number_filename, case_.tumor_copy_number_file_id).Where(r => Math.Abs(r.Segment_Mean) > 1.0).ToList();

                //
                // Get all of the ASE candidates 
                var variants = ASETools.AnnotatedVariant.readFile(case_.annotated_selected_variants_filename).Where(x => !x.somaticMutation && x.IsASECandidate(true, tumorCopyNumberVariation, configuration, perGeneASEMap, geneMap)).ToList();

                if (variants.Count() == 0)
                {
                    continue;
                }

                variants.Sort();

                var bestRunsByASE = new ASERun[nRanges];
                var currentRunsByASE = new ASERun[nRanges];
                for (int i = 0; i < nRanges; i++)
                {
                    bestRunsByASE[i] = new ASERun(minRangeASE + i * ASERangeIncrement);
                    currentRunsByASE[i] = new ASERun(minRangeASE + i * ASERangeIncrement);
                }

                int nVariants = variants.Count();
                for (int i = 0; i < nVariants; i++)
                {
                    var ase = variants[i].GetTumorAlleleSpecificExpression();

                    for (int range = 0; range < nRanges; range++)
                    {
                        if (ase < minRangeASE + range * ASERangeIncrement || currentRunsByASE[range].chromosome != variants[i].contig)
                        {
                            //
                            // End the current range because of too little ASE or end-of-chromosome.
                            //
                            if (bestRunsByASE[range].CompareTo(currentRunsByASE[range]) < 0) {
                                bestRunsByASE[range] = currentRunsByASE[range];
                            }
                        } else
                        {
                            bestRunsByASE[range].recordObservation(variants[i].locus, ase);
                        }
                    } // for range
                } // for variant



            } // while true
        }
    }
}
