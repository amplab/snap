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

        static int nProcessed = 0;
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

            public static PerThreadState create() { return new PerThreadState(); }

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

            public static void ItemDequeued()
            {
                if (nProcessed++ % 100 == 0)
                {
                    Console.Write(".");
                }
            }

            public static void FinishUp(PerThreadState perThreadState)
            {
                perThreadState.FinishUp();
            }

            void FinishUp()
            {
                ParallelMergeOptions into global
            }
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

            Console.Write("Processing " + casesToProcess.Count() + " cases, 1 dot/100 cases: ");

            var threading = new ASETools.WorkerThreadHelper<ASETools.Case, PerThreadState>(casesToProcess, PerThreadState.HandleOneCase, PerThreadState.FinishUp, PerThreadState.create,
                PerThreadState.ItemDqueued);
        }
    }
}
