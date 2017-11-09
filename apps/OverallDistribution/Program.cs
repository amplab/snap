using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using ASELib;
using System.Diagnostics;
using System.IO;

namespace OverallDistribution
{
    class Program
    {

        static ASETools.Configuration configuration;
        static ASETools.GeneMap geneMap;
        static Dictionary<string, ASETools.ASEMapPerGeneLine> perGeneASEMap;
        static ASETools.GeneLocationsByNameAndChromosome geneLocationInformation;

        class PerThreadState
        {
            public PerThreadState()
            {
                foreach (var tumor in ASETools.BothBools)
                {
                    perThreadResult.Add(tumor, new Dictionary<bool, ASETools.Histogram>());
                    foreach (var somatic in ASETools.BothBools)
                    {
                        perThreadResult[tumor].Add(somatic, new ASETools.Histogram());
                    }
                }
            }

            public Dictionary<bool, Dictionary<bool, ASETools.Histogram>> perThreadResult = new Dictionary<bool, Dictionary<bool, ASETools.Histogram>>();

            public static void HandleOneCase(ASETools.Case case_, PerThreadState perThreadState)
            {
                perThreadState.HandleOneCase(case_);
            }

            void HandleOneCase(ASETools.Case case_)
            {
                var annotatedSelectedVariants = ASETools.AnnotatedVariant.readFile(case_.annotated_selected_variants_filename);
                var copyNumber = ASETools.CopyNumberVariation.ReadBothFiles(case_);

                foreach (var variant in annotatedSelectedVariants)
                {
                    foreach (var tumor in ASETools.BothBools)
                    {
                        foreach (var somatic in ASETools.BothBools)
                        {
                            if (!variant.somaticMutation == somatic)
                            {
                                continue;
                            }

                            if (variant.IsASECandidate(tumor, copyNumber, configuration, perGeneASEMap, geneMap))
                            {
                                perThreadResult[tumor][somatic].addValue(variant.GetAlleleSpecificExpression(tumor));
                            }
                        } // somatic
                    } // tumor
                } // variant
            } // HandleOneCase


            public static void FinishUp(PerThreadState perThreadState)
            {
                perThreadState.FinishUp();
            }

            void FinishUp()
            {
                foreach (var somatic in ASETools.BothBools)
                {
                    foreach (var tumor in ASETools.BothBools)
                    {
                        overallResult[tumor][somatic].merge(perThreadResult[tumor][somatic]);
                    }
                }
            } // FinishUp
        } // PerThreadState

        static Dictionary<bool, Dictionary<bool, ASETools.Histogram>> overallResult = new Dictionary<bool, Dictionary<bool, ASETools.Histogram>>();

        static void Main(string[] args)
        {
            var timer = new Stopwatch();
            timer.Start();

            foreach (var tumor in ASETools.BothBools)
            {
                overallResult.Add(tumor, new Dictionary<bool, ASETools.Histogram>());
                foreach (var somatic in ASETools.BothBools)
                {
                    overallResult[tumor].Add(somatic, new ASETools.Histogram());
                }
            }

            configuration = ASETools.Configuration.loadFromFile(args);

            if (null == configuration)
            {
                Console.WriteLine("Giving up because we were unable to load configuration.");
                return;
            }

            geneLocationInformation = new ASETools.GeneLocationsByNameAndChromosome(ASETools.readKnownGeneFile(configuration.geneLocationInformationFilename));
            geneMap = new ASETools.GeneMap(geneLocationInformation.genesByName);

            perGeneASEMap = ASETools.ASEMapPerGeneLine.ReadFromFileToDictionary(configuration.finalResultsDirectory + ASETools.PerGeneASEMapFilename);

            if (null == perGeneASEMap)
            {
                Console.WriteLine("You must first create the per-gene ASE map in " + configuration.finalResultsDirectory + ASETools.PerGeneASEMapFilename);
                return;
            }

            var cases = ASETools.Case.LoadCases(configuration.casesFilePathname);

            if (null == cases)
            {
                Console.WriteLine("Unable to load cases.");
            }

            var casesToProcess = cases.Select(x => x.Value).Where(x => x.annotated_selected_variants_filename != "" && x.tumor_copy_number_filename != "").ToList();
            int n = casesToProcess.Count();


            var threading = new ASETools.WorkerThreadHelper<ASETools.Case, PerThreadState>(casesToProcess, PerThreadState.HandleOneCase, PerThreadState.FinishUp, null, 100);

            Console.Write("Processing " + casesToProcess.Count() + " cases, 1 dot/100 cases: ");
            threading.run();

            var outputFile = ASETools.CreateStreamWriterWithRetry(configuration.finalResultsDirectory + ASETools.OverallASEFilename);

            foreach (var germline in ASETools.BothBools)    // Doing germline rather than somatic here is in order to get the germlines ones to come out first in the output file.
            {
                foreach (var tumor in ASETools.BothBools)
                {
                    outputFile.WriteLine("Overall ASE distribution somatic: " + !germline + " tumor: " + tumor);
                    outputFile.WriteLine(ASETools.HistogramResultLine.Header());
                    overallResult[tumor][!germline].ComputeHistogram(0, 1, 0.01).ToList().ForEach(x => outputFile.WriteLine(x.ToString()));
                    outputFile.WriteLine();
                }
            }

            outputFile.WriteLine("**done**");
            outputFile.Close();

            Console.WriteLine("Completed in " + ASETools.ElapsedTimeInSeconds(timer));
        } // Main
    }
}
