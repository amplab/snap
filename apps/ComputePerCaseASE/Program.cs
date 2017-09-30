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
        static Dictionary<bool, Dictionary<string, ASETools.ASEMapPerGeneLine>> perGeneASEMap;

        static List<string> outputLines = new List<string>();

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

            var casesToProcess = cases.Where(x => x.Value.per_case_ase_filename == "" && x.Value.annotated_selected_variants_filename != "" && x.Value.tumor_copy_number_filename != "").Select(x => x.Value).ToList();
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
            }
            outputFile.WriteLine();

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

            public bool hasAnyObservations()
            {
                return observations.Count() != 0;
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

                var annotatedSelectedVariants = ASETools.AnnotatedVariant.readFile(case_.annotated_selected_variants_filename);
                var tumorCopyNumberVariation = ASETools.CopyNumberVariation.ReadFile(case_.tumor_copy_number_filename, case_.tumor_copy_number_file_id).Where(r => Math.Abs(r.Segment_Mean) > 1.0).ToList();
                List<ASETools.CopyNumberVariation> normalCopyNumberVariation = null;
                if (case_.normal_copy_number_filename != "")
                {
                    normalCopyNumberVariation = ASETools.CopyNumberVariation.ReadFile(case_.normal_copy_number_filename, case_.normal_copy_number_file_id).Where(r => Math.Abs(r.Segment_Mean) > 1.0).ToList();
                }

                if (annotatedSelectedVariants == null)
                {
                    Console.WriteLine("Unable to read annotated selected variants from " + case_.annotated_selected_variants_filename);
                    continue;
                }


                var overall = new Dictionary<bool, Measurements>();
                var perChromosome = new Dictionary<bool, Measurements[]>();
                foreach (var tumor in ASETools.BothBools)
                {
                    overall.Add(tumor, new Measurements());
                    perChromosome.Add(tumor, new Measurements[ASETools.nHumanNuclearChromosomes]);
                    for (int i = 0; i < ASETools.nHumanNuclearChromosomes; i++)
                    {
                        perChromosome[tumor][i] = new Measurements();
                    }
                }

                foreach (var variant in annotatedSelectedVariants)
                {
                    int index = ASETools.ChromosomeNameToIndex(variant.contig);

                    foreach (bool tumor in ASETools.BothBools)
                    {
                        if (!variant.somaticMutation && variant.IsASECandidate(tumor, normalCopyNumberVariation, configuration, perGeneASEMap, geneMap))
                        {
                            overall[tumor].recordObservation(variant.GetAlleleSpecificExpression(tumor));
                            if (index != -1)
                            {
                                perChromosome[tumor][index].recordObservation(variant.GetAlleleSpecificExpression(tumor));
                            }
                        }
                    } // foreach tumor/normal
                } // foreach variant

                string outputLine = case_.case_id;

                foreach (var notTumor in ASETools.BothBools) // We use notTumor here because we want to have normals first, since it makes the output file easier to work with
                {
                    if (overall[!notTumor].hasAnyObservations())
                    {
                        outputLine += "\t" + overall[!notTumor].getASE();
                    } else
                    {
                        outputLine += "\t*";
                    }

                    var perChromosomeASEs = new List<double>();

                    for (int i = 1; i <= ASETools.nHumanAutosomes; i++)
                    {
                        if (perChromosome[!notTumor][i].hasAnyObservations())
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
                        outputLine += "\t" + perChromosomeASEs[0] + "\t" + perChromosomeASEs[perChromosomeASEs.Count() - 1];        // Min & Max
                    }
                } // foreach notTumor
                localOutputLines.Add(outputLine);
            } // while true
        } // WorkerThread
    }
}
