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
        static Dictionary<bool, Dictionary<bool, Dictionary<bool, ASETools.PreBucketedHistogram>>> readDepthDistribution = new Dictionary<bool, Dictionary<bool, Dictionary<bool, ASETools.PreBucketedHistogram>>>();  // tumor, somatic, DNA in that order
        static Dictionary<bool,Dictionary<bool, Dictionary<int, ASETools.PreBucketedHistogram>>> ASEByReadDepthDistribution = new Dictionary<bool, Dictionary<bool, Dictionary<int, ASETools.PreBucketedHistogram>>>();  // Maps tumor->somatic->read depth->ASE distribution.

        class PerThreadState
        {
            public PerThreadState()
            {
                foreach (var tumor in ASETools.BothBools)
                {
                    perThreadResult.Add(tumor, new Dictionary<bool, ASETools.PreBucketedHistogram>());
                    perThreadReadDepthDistribution.Add(tumor, new Dictionary<bool, Dictionary<bool, ASETools.PreBucketedHistogram>>());
                    perThreadASEByReadDepthDistribution.Add(tumor, new Dictionary<bool, Dictionary<int, ASETools.PreBucketedHistogram>>());
                    foreach (var somatic in ASETools.BothBools)
                    {
                        perThreadResult[tumor].Add(somatic, new ASETools.PreBucketedHistogram(0, 1, 0.01));
                        perThreadReadDepthDistribution[tumor].Add(somatic, new Dictionary<bool, ASETools.PreBucketedHistogram>());
                        perThreadASEByReadDepthDistribution[tumor].Add(somatic, new Dictionary<int, ASETools.PreBucketedHistogram>());

                        foreach (var dna in ASETools.BothBools)
                        {
                            perThreadReadDepthDistribution[tumor][somatic].Add(dna, new ASETools.PreBucketedHistogram(0, 1000, 1));
                        }
                    }
                }
            }

            Dictionary<bool, Dictionary<bool, ASETools.PreBucketedHistogram>> perThreadResult = new Dictionary<bool, Dictionary<bool, ASETools.PreBucketedHistogram>>();
            Dictionary<bool, Dictionary<bool, Dictionary<bool, ASETools.PreBucketedHistogram>>> perThreadReadDepthDistribution = new Dictionary<bool, Dictionary<bool, Dictionary<bool, ASETools.PreBucketedHistogram>>>();    // tumor, somatic, DNA in that order
            Dictionary<bool, Dictionary<bool, Dictionary<int, ASETools.PreBucketedHistogram>>> perThreadASEByReadDepthDistribution = new Dictionary<bool, Dictionary<bool, Dictionary<int, ASETools.PreBucketedHistogram>>>();  // Maps tumor->somatic->read depth->ASE distribution.

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

                            if (variant.IsASECandidate(tumor, copyNumber, configuration, perGeneASEMap, geneMap, 1))   // Force read depth to 1, since we're measuring the read depth distribution (though < 10 was filtered upstream)
                            {
                                foreach (var dna in ASETools.BothBools)
                                {
                                    perThreadReadDepthDistribution[tumor][somatic][dna].addValue(variant.getReadCount(tumor, dna).usefulReads());
                                }
                            }

                            if (variant.IsASECandidate(tumor, copyNumber, configuration, perGeneASEMap, geneMap))
                            {
                                var ASE = variant.GetAlleleSpecificExpression(tumor);
                                perThreadResult[tumor][somatic].addValue(ASE);

                                int readDepth = variant.getReadCount(tumor, false).usefulReads();
                                if (!perThreadASEByReadDepthDistribution[tumor][somatic].ContainsKey(readDepth))
                                {
                                    perThreadASEByReadDepthDistribution[tumor][somatic].Add(readDepth, new ASETools.PreBucketedHistogram(0, 1, 0.01));
                                }

                                perThreadASEByReadDepthDistribution[tumor][somatic][readDepth].addValue(ASE);
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

                        var readDepths = perThreadASEByReadDepthDistribution[tumor][somatic].Select(x => x.Key).ToList();    // We need to enumerate this way so we can modify perThreadASEByReadDepthDistribution inside the loop
                        foreach (var readDepth in readDepths)
                        {
                            var histogram = perThreadASEByReadDepthDistribution[tumor][somatic][readDepth];
                            if (!ASEByReadDepthDistribution[tumor][somatic].ContainsKey(readDepth))
                            {
                                ASEByReadDepthDistribution[tumor][somatic].Add(readDepth, histogram);
                            }
                            else
                            {
                                ASEByReadDepthDistribution[tumor][somatic][readDepth].merge(histogram);
                            }
                        } // for each read depth

                        foreach (var dna in ASETools.BothBools)
                        {
                            readDepthDistribution[tumor][somatic][dna].merge(perThreadReadDepthDistribution[tumor][somatic][dna]);
                        }
                    }
                }
            } // FinishUp
        } // PerThreadState

        static Dictionary<bool, Dictionary<bool, ASETools.PreBucketedHistogram>> overallResult = new Dictionary<bool, Dictionary<bool, ASETools.PreBucketedHistogram>>();

        static void Main(string[] args)
        {
            var timer = new Stopwatch();
            timer.Start();

            foreach (var tumor in ASETools.BothBools)
            {
                overallResult.Add(tumor, new Dictionary<bool, ASETools.PreBucketedHistogram>());
                readDepthDistribution.Add(tumor, new Dictionary<bool, Dictionary<bool, ASETools.PreBucketedHistogram>>());
                ASEByReadDepthDistribution.Add(tumor, new Dictionary<bool, Dictionary<int, ASETools.PreBucketedHistogram>>());
                foreach (var somatic in ASETools.BothBools)
                {
                    overallResult[tumor].Add(somatic, new ASETools.PreBucketedHistogram(0, 1, 0.01));
                    readDepthDistribution[tumor].Add(somatic, new Dictionary<bool, ASETools.PreBucketedHistogram>());
                    ASEByReadDepthDistribution[tumor].Add(somatic, new Dictionary<int, ASETools.PreBucketedHistogram>());

                    foreach (var dna in ASETools.BothBools)
                    {
                        readDepthDistribution[tumor][somatic].Add(dna, new ASETools.PreBucketedHistogram(0, 1000, 1));
                    }
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
                    overallResult[tumor][!germline].ComputeHistogram().ToList().ForEach(x => outputFile.WriteLine(x.ToString()));
                    outputFile.WriteLine();
                }
            }

            //
            // Write the tumor germline ASE distribution into its own file, beause it's needed programmatically elsewhere.
            //
            var tumorGermlineASEDistributionFile = ASETools.CreateStreamWriterWithRetry(configuration.finalResultsDirectory + ASETools.TumorGermlineASEDistributionFilename);
            tumorGermlineASEDistributionFile.WriteLine(ASETools.HistogramResultLine.Header());
            overallResult[true][false].ComputeHistogram().ToList().ForEach(x => tumorGermlineASEDistributionFile.WriteLine(x.ToString()));
            tumorGermlineASEDistributionFile.WriteLine("**done**");
            tumorGermlineASEDistributionFile.Close();

            foreach (var tumor in ASETools.BothBools)
            {
                foreach (var somatic in ASETools.BothBools)
                {
                    foreach (var dna in ASETools.BothBools)
                    {
                        outputFile.WriteLine("Read depth distribution tumor: " + tumor + ", somatic: " + somatic + ", dna: " + dna);
                        outputFile.WriteLine(ASETools.HistogramResultLine.Header());
                        readDepthDistribution[tumor][somatic][dna].ComputeHistogram().ToList().ForEach(x => outputFile.WriteLine(x.ToString()));
                        outputFile.WriteLine();
                    }
                }
            }

            //
            // Write the tumor RNA read depth into its own file, because it's needed programmatically elsewhere.
            //
            var tumorRNAReadDepthFile = ASETools.CreateStreamWriterWithRetry(configuration.finalResultsDirectory + ASETools.TumorRNAReadDepthDistributionFilename);
            tumorRNAReadDepthFile.WriteLine(ASETools.HistogramResultLine.Header());
            readDepthDistribution[true][false][false].ComputeHistogram().ToList().ForEach(x => tumorRNAReadDepthFile.WriteLine(x.ToString()));
            tumorRNAReadDepthFile.WriteLine("**done**");
            tumorRNAReadDepthFile.Close();


            foreach (var tumor in ASETools.BothBools) {
                foreach (var somatic in ASETools.BothBools)
                {
                    if (!tumor && somatic)
                    {
                        //
                        // Germline somatic is meaningless/just pipeline errors. Skip it.
                        //
                        continue;
                    }
                    for (int readDepth = 10; readDepth <= 500; readDepth++)
                    {
                        if (ASEByReadDepthDistribution[tumor][somatic].ContainsKey(readDepth))
                        {
                            outputFile.WriteLine("ASE distributon for tumor: " + tumor + " somatic: " + somatic + " sites with read depth exactly " + readDepth);
                            outputFile.WriteLine(ASETools.HistogramResultLine.Header());
                            ASEByReadDepthDistribution[tumor][somatic][readDepth].ComputeHistogram().ToList().ForEach(x => outputFile.WriteLine(x.ToString()));
                            outputFile.WriteLine();

                        }
                    }
                }
            }

            outputFile.WriteLine("**done**");
            outputFile.Close();

            Console.WriteLine("Completed in " + ASETools.ElapsedTimeInSeconds(timer));
        } // Main
    }
}
