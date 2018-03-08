using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using System.Diagnostics;
using ASELib;
using System.Threading;

namespace ComputePerCaseASE
{
    class Program
    {
        static ASETools.GeneMap geneMap;
        static ASETools.Configuration configuration;
        static Dictionary<string, ASETools.ASEMapPerGeneLine> perGeneASEMap;
        static ASETools.ASECorrection ASECorrection;


        static Dictionary<string, int> perCaseRawMAFCount = new Dictionary<string, int>();

        static List<string> outputLines = new List<string>();

        class MAFCountPerThreadState
        {
            public MAFCountPerThreadState()
            {

            }

            public static void HandleOneItem(List<ASETools.Case> cases, MAFCountPerThreadState perThreadState)
            {
                perThreadState.HandleOneItem(cases);
            }

            void HandleOneItem(List<ASETools.Case> cases)
            {
                var mafLines =  ASETools.MAFLine.ReadFile(cases[0].maf_filename, cases[0].maf_file_id, true);

                foreach (var case_ in cases)
                {
                    localPerCaseMAFCount.Add(case_.case_id, mafLines.Where(x => x.tumor_bam_uuid == case_.tumor_dna_file_id).Count());
                }
            }

            public static void FinishUp(MAFCountPerThreadState state)
            {
                foreach (var perCaseEntry in state.localPerCaseMAFCount)
                {
                    perCaseRawMAFCount.Add(perCaseEntry.Key, perCaseEntry.Value);
                }
            }

            Dictionary<string, int> localPerCaseMAFCount = new Dictionary<string, int>();
        }

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
            ASECorrection = ASETools.ASECorrection.LoadFromFile(configuration.finalResultsDirectory + ASETools.ASECorrectionFilename);
            if (null == ASECorrection)
            {
                Console.WriteLine("Unable to load ASE correction");
                return;
            }


            var cases = ASETools.Case.LoadCases(configuration.casesFilePathname);

            if (null == cases)
            {
                Console.WriteLine("Unable to load cases.");
            }

            var casesToProcess = cases.Where(x => x.Value.annotated_selected_variants_filename != "" && x.Value.tumor_copy_number_filename != "").Select(x => x.Value).ToList();

            var casesByMAFFile = casesToProcess.GroupByToDict(x => x.maf_filename).Select(x => x.Value).ToList();

            var mafCountThreading = new ASETools.WorkerThreadHelper<List<ASETools.Case>, MAFCountPerThreadState>(casesByMAFFile, MAFCountPerThreadState.HandleOneItem, MAFCountPerThreadState.FinishUp, null, 1);
            Console.Write("Loading raw maf lines from " + casesByMAFFile.Count() + " maf files, 1 dot/file: ");
            mafCountThreading.run();


            int nCases = casesToProcess.Count();
            Console.Write("Processing " + nCases + " cases, 1 dot/100 cases: ");

            var threads = new List<Thread>();
            for (int i = 0; i < Environment.ProcessorCount; i++)
            {
                threads.Add(new Thread(() => WorkerThread(casesToProcess)));
            }

            threads.ForEach(t => t.Start());
            threads.ForEach(t => t.Join());

            var outputFilename = configuration.finalResultsDirectory + ASETools.PerCaseASEFilename;
            var outputFile = ASETools.CreateStreamWriterWithRetry(outputFilename);
            if (outputFile == null)
            {
                Console.WriteLine("Unable to open output file " + outputFilename + ". Aborting.");
                return;
            }

            outputFile.Write("Case ID");

            foreach (var notTumor in ASETools.BothBools)
            {
                outputFile.Write("\t" + ASETools.tumorToString[!notTumor] + " ASE");
                for (int i = 1; i <= ASETools.nHumanAutosomes; i++)
                {
                    outputFile.Write("\t" + ASETools.tumorToString[!notTumor] + " chr" + i + " ASE");
                }
                outputFile.Write("\t" + ASETools.tumorToString[!notTumor] + " median chromosome ASE\t" + ASETools.tumorToString[!notTumor] + " min chromosome ASE\t" + ASETools.tumorToString[!notTumor] + " max chromosome ASE");
                outputFile.Write("\t" + ASETools.tumorToString[!notTumor] + " 10MB region count" + "\t" + ASETools.tumorToString[!notTumor] + " 10MB region median " + "\t" + ASETools.tumorToString[!notTumor] + " 10MB region min" + "\t" +
                    ASETools.tumorToString[!notTumor] + " 10MB region max" + "\t" + ASETools.tumorToString[!notTumor] + " 10MB region mean");
            }
            outputFile.WriteLine("\tSelected MAF Line Count\tRaw MAF Line Count");

            foreach (var outputLine in outputLines)
            {
                outputFile.WriteLine(outputLine);
            }

            outputFile.WriteLine("**done**");
            outputFile.Close();

            Console.WriteLine();
            Console.WriteLine("Processed " + nCases + " in " + ASETools.ElapsedTimeInSeconds(timer));
        }

        static int nCasesProcessed = 0;

        class Measurements
        {
            List<double> observations = new List<double>();

            public void recordObservation(double observation)
            {
                observations.Add(observation);
            }

            public double getASE()
            {
                return observations.Average();
            }

            public bool hasEnoughObservations()
            {
                return observations.Count() >= 30;
            }

            public double getMedian()
            {
                int n = observations.Count();

                if (n == 0)
                {
                    return 0;
                }

                observations.Sort();

                if (n % 2 == 1)
                {
                    return observations[n / 2];
                } else
                {
                    return (observations[n / 2] + observations[n / 2 + 1]) / 2;
                }
            }
        }

        static void WorkerThread(List<ASETools.Case> casesToProcess)
        {
            var localOutputLines = new List<string>();

            while (true)
            {
                ASETools.Case case_;

                lock (casesToProcess)
                {
                    if (casesToProcess.Count() == 0)
                    {
                        outputLines.AddRange(localOutputLines);
                        return;
                    }

                    case_ = casesToProcess[0];
                    casesToProcess.RemoveAt(0);

                    nCasesProcessed++;
                    if (nCasesProcessed % 100 == 0)
                    {
                        Console.Write('.');
                    }
                } // lock

                var annotatedSelectedVariants = ASETools.AnnotatedVariant.readFile(case_.annotated_selected_variants_filename).Where(x => !x.CausesNonsenseMediatedDecay()).ToList();

                var copyNumber = ASETools.CopyNumberVariation.ReadBothFiles(case_);
 
                if (annotatedSelectedVariants == null)
                {
                    Console.WriteLine("Unable to read annotated selected variants from " + case_.annotated_selected_variants_filename);
                    continue;
                }


                var overall = new Dictionary<bool, Measurements>();
                var perChromosome = new Dictionary<bool, Measurements[]>();
                var per10MB = new Dictionary<bool, List<Measurements>>();
 
                foreach (var tumor in ASETools.BothBools)
                {
                    overall.Add(tumor, new Measurements());
                    perChromosome.Add(tumor, new Measurements[ASETools.nHumanNuclearChromosomes]);
                    for (int i = 0; i < ASETools.nHumanNuclearChromosomes; i++)
                    {
                        perChromosome[tumor][i] = new Measurements();
                    }
                    per10MB.Add(tumor, new List<Measurements>());
                }

                string currentChromosome = "";
                int currentLocus = 0;
                annotatedSelectedVariants.Sort();

                for (int i = 0; i < annotatedSelectedVariants.Count(); i++)
                {
                    var variant = annotatedSelectedVariants[i];
                    int index = ASETools.ChromosomeNameToIndex(variant.contig);

                    if (variant.contig != currentChromosome || variant.locus / 10000000 != currentLocus / 10000000)
                    {
                        foreach (bool tumor in ASETools.BothBools)
                        {
                            per10MB[tumor].Insert(0, new Measurements());
                        }

                        currentChromosome = variant.contig;
                        currentLocus = variant.locus;
                    }

                    foreach (bool tumor in ASETools.BothBools)
                    {
                        if (!variant.somaticMutation && variant.IsASECandidate(tumor, copyNumber, configuration, perGeneASEMap, geneMap))
                        {
                            overall[tumor].recordObservation(variant.GetAlleleSpecificExpression(tumor, ASECorrection));
                            per10MB[tumor][0].recordObservation(variant.GetAlleleSpecificExpression(tumor, ASECorrection));
                            if (index != -1)
                            {
                                perChromosome[tumor][index].recordObservation(variant.GetAlleleSpecificExpression(tumor, ASECorrection));
                            }
                        }
                    } // foreach tumor/normal
                } // foreach variant

                string outputLine = case_.case_id;

                foreach (var notTumor in ASETools.BothBools) // We use notTumor here because we want to have normals first, since it makes the output file easier to work with
                {
                    if (overall[!notTumor].hasEnoughObservations())
                    {
                        outputLine += "\t" + overall[!notTumor].getASE();
                    } else
                    {
                        outputLine += "\t*";
                    }

                    var perChromosomeASEs = new List<double>();

                    for (int i = 1; i <= ASETools.nHumanAutosomes; i++)
                    {
                        if (perChromosome[!notTumor][i].hasEnoughObservations())
                        {
                            perChromosomeASEs.Add(perChromosome[!notTumor][i].getASE());
                            outputLine += "\t" + perChromosome[!notTumor][i].getASE();
                        } else
                        {
                            outputLine += "\t*";
                        }
                    }  // for chromosomes

 
                    int n = perChromosomeASEs.Count();
                    if (n == 0)
                    {
                        outputLine += "\t*\t*\t*";
                    }
                    else
                    {
                        perChromosomeASEs.Sort();
                        if (n % 2 == 1)
                        {
                            outputLine += "\t" + perChromosomeASEs[n / 2];
                        }
                        else
                        {
                            outputLine += "\t" + (perChromosomeASEs[n / 2] + perChromosomeASEs[n / 2 - 1]) / 2;
                        }
                        outputLine += "\t" + perChromosomeASEs[0] + "\t" + perChromosomeASEs[n - 1];        // Min & Max
                    }

                    var good10MBRegions = per10MB[!notTumor].Where(x => x.hasEnoughObservations()).Select(x => x.getASE()).ToList();
                    n = good10MBRegions.Count();
                    outputLine += "\t" + n;
                    if (n < 10)
                    {
                        outputLine += "\t*\t*\t*\t*";
                    } else
                    {
                        good10MBRegions.Sort();
                        // median
                        if (n % 2 == 1)
                        {
                            outputLine += "\t" + good10MBRegions[n / 2];
                        } else
                        {
                            outputLine += "\t" + (good10MBRegions[n / 2] + good10MBRegions[n / 2 - 1]) / 2;
                        }

                        // min, max, mean
                        outputLine += "\t" + good10MBRegions[0] + "\t" + good10MBRegions[n - 1] + "\t" + good10MBRegions.Average();
                    }
                } // foreach notTumor
                outputLine += "\t" + ASETools.MAFLine.ReadFile(case_.extracted_maf_lines_filename, case_.case_id, false).Count() + "\t" + perCaseRawMAFCount[case_.case_id];
                localOutputLines.Add(outputLine);
            } // while true
        } // WorkerThread
    }
}
