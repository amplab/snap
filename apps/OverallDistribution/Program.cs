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

        const int maxReadDepthInHistograms = 1000;

        static ASETools.Configuration configuration;
        static ASETools.GeneMap geneMap;
        static Dictionary<string, ASETools.ASEMapPerGeneLine> perGeneASEMap;
        static ASETools.GeneLocationsByNameAndChromosome geneLocationInformation;
        static Dictionary<bool, Dictionary<bool, Dictionary<bool, ASETools.PreBucketedHistogram>>> readDepthDistribution = new Dictionary<bool, Dictionary<bool, Dictionary<bool, ASETools.PreBucketedHistogram>>>();  // tumor, somatic, DNA in that order
        static Dictionary<bool,Dictionary<bool, Dictionary<int, ASETools.PreBucketedHistogram>>> ASEByReadDepthDistribution = new Dictionary<bool, Dictionary<bool, Dictionary<int, ASETools.PreBucketedHistogram>>>();  // Maps tumor->somatic->read depth->ASE distribution.
        static Dictionary<string, ASETools.PreBucketedHistogram> germlineTumorDNAReadDepthByDiseaseType = new Dictionary<string, ASETools.PreBucketedHistogram>();
        static Dictionary<string, ASETools.PreBucketedHistogram> validGermlineASESitesByDisease = new Dictionary<string, ASETools.PreBucketedHistogram>();
        static ASETools.PreBucketedHistogram overallValidGermlineASESites = new ASETools.PreBucketedHistogram(0, 2000, 1);

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
                        perThreadResult[tumor].Add(somatic, new ASETools.PreBucketedHistogram(0, 1, 0.01, "per thread ASE"));
                        perThreadReadDepthDistribution[tumor].Add(somatic, new Dictionary<bool, ASETools.PreBucketedHistogram>());
                        perThreadASEByReadDepthDistribution[tumor].Add(somatic, new Dictionary<int, ASETools.PreBucketedHistogram>());

                        foreach (var dna in ASETools.BothBools)
                        {
                            perThreadReadDepthDistribution[tumor][somatic].Add(dna, new ASETools.PreBucketedHistogram(0, maxReadDepthInHistograms, 1, "per thread read depth"));
                        }
                    }
                }
            }

            Dictionary<bool, Dictionary<bool, ASETools.PreBucketedHistogram>> perThreadResult = new Dictionary<bool, Dictionary<bool, ASETools.PreBucketedHistogram>>();
            Dictionary<bool, Dictionary<bool, Dictionary<bool, ASETools.PreBucketedHistogram>>> perThreadReadDepthDistribution = new Dictionary<bool, Dictionary<bool, Dictionary<bool, ASETools.PreBucketedHistogram>>>();    // tumor, somatic, DNA in that order
            Dictionary<bool, Dictionary<bool, Dictionary<int, ASETools.PreBucketedHistogram>>> perThreadASEByReadDepthDistribution = new Dictionary<bool, Dictionary<bool, Dictionary<int, ASETools.PreBucketedHistogram>>>();  // Maps tumor->somatic->read depth->ASE distribution.
            Dictionary<string, ASETools.PreBucketedHistogram> localGermlineTumorDNAReadDepthByDiseaseType = new Dictionary<string, ASETools.PreBucketedHistogram>();
            Dictionary<string, ASETools.PreBucketedHistogram> localValidGermlineASESitesByDisease = new Dictionary<string, ASETools.PreBucketedHistogram>();
            ASETools.PreBucketedHistogram localOverallValidGermlineASESites = new ASETools.PreBucketedHistogram(0, 2000, 1, "local overall valid germline ASE sites");

            public static void HandleOneCase(ASETools.Case case_, PerThreadState perThreadState)
            {
                perThreadState.HandleOneCase(case_);
            }

            void HandleOneCase(ASETools.Case case_)
            {
                var annotatedSelectedVariants = ASETools.AnnotatedVariant.readFile(case_.annotated_selected_variants_filename);
                var copyNumber = ASETools.CopyNumberVariation.ReadBothFiles(case_);
                var disease = case_.disease();
                int nValidASESites = 0;

                if (!localGermlineTumorDNAReadDepthByDiseaseType.ContainsKey(disease))
                {
                    localGermlineTumorDNAReadDepthByDiseaseType.Add(disease, new ASETools.PreBucketedHistogram(0, maxReadDepthInHistograms, 1, "per thread read depth by disease"));
                }

                if (!localValidGermlineASESitesByDisease.ContainsKey(disease))
                {
                    localValidGermlineASESitesByDisease.Add(disease, new ASETools.PreBucketedHistogram(0, 2000, 1, "per thread germline ASE site count by disease"));
                }

                foreach (var variant in annotatedSelectedVariants)
                {
                    foreach (var tumor in ASETools.BothBools)
                    {
                        foreach (var somatic in ASETools.BothBools)
                        {
                            if (!variant.somaticMutation == somatic || variant.isSilent())
                            {
                                continue;
                            }

                            if (variant.IsASECandidate(tumor, copyNumber, configuration, perGeneASEMap, geneMap, 1))   // Force read depth to 1, since we're measuring the read depth distribution (though < 10 was filtered upstream)
                            {
                                foreach (var dna in ASETools.BothBools)
                                {
                                    perThreadReadDepthDistribution[tumor][somatic][dna].addValue(variant.getReadCount(tumor, dna).usefulReads());
                                }

                                if (!somatic && tumor)
                                {
                                    localGermlineTumorDNAReadDepthByDiseaseType[disease].addValue(variant.getReadCount(true, true).usefulReads());
                                }
                            }

                            if (variant.IsASECandidate(tumor, copyNumber, configuration, perGeneASEMap, geneMap))
                            {
                                var ASE = variant.GetAlleleSpecificExpression(tumor, aseCorrection);
                                perThreadResult[tumor][somatic].addValue(ASE);

                                int readDepth = variant.getReadCount(tumor, false).usefulReads();
                                if (!perThreadASEByReadDepthDistribution[tumor][somatic].ContainsKey(readDepth))
                                {
                                    perThreadASEByReadDepthDistribution[tumor][somatic].Add(readDepth, new ASETools.PreBucketedHistogram(0, 1, 0.01, "per thread ASE by read depth"));
                                }

                                perThreadASEByReadDepthDistribution[tumor][somatic][readDepth].addValue(ASE);

                                if (!somatic)
                                {
                                    nValidASESites++;
                                }
                            }
                        } // somatic
                    } // tumor
                } // variant

                localValidGermlineASESitesByDisease[disease].addValue(nValidASESites);
                localOverallValidGermlineASESites.addValue(nValidASESites);
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

                foreach (var perDiseaseEntry in localGermlineTumorDNAReadDepthByDiseaseType)
                {
                    var disease = perDiseaseEntry.Key;
                    var histogram = perDiseaseEntry.Value;
                    if (!germlineTumorDNAReadDepthByDiseaseType.ContainsKey(disease))
                    {
                        germlineTumorDNAReadDepthByDiseaseType.Add(disease, histogram);
                    } else
                    {
                        germlineTumorDNAReadDepthByDiseaseType[disease].merge(histogram);
                    }
                }

                foreach (var perDiseaseEntry in localValidGermlineASESitesByDisease)
                {
                    var disease = perDiseaseEntry.Key;
                    var histogram = perDiseaseEntry.Value;
                    if (!validGermlineASESitesByDisease.ContainsKey(disease))
                    {
                        validGermlineASESitesByDisease.Add(disease, histogram);
                    }
                    else
                    {
                        validGermlineASESitesByDisease[disease].merge(histogram);
                    }
                }

                overallValidGermlineASESites.merge(localOverallValidGermlineASESites);
            } // FinishUp
        } // PerThreadState

        static Dictionary<bool, Dictionary<bool, ASETools.PreBucketedHistogram>> overallResult = new Dictionary<bool, Dictionary<bool, ASETools.PreBucketedHistogram>>();
        static ASETools.ASECorrection aseCorrection = null;

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
                    overallResult[tumor].Add(somatic, new ASETools.PreBucketedHistogram(0, 1, 0.01, "overall ASE"));
                    readDepthDistribution[tumor].Add(somatic, new Dictionary<bool, ASETools.PreBucketedHistogram>());
                    ASEByReadDepthDistribution[tumor].Add(somatic, new Dictionary<int, ASETools.PreBucketedHistogram>());

                    foreach (var dna in ASETools.BothBools)
                    {
                        readDepthDistribution[tumor][somatic].Add(dna, new ASETools.PreBucketedHistogram(0, maxReadDepthInHistograms, 1, "overall read depth"));
                    }
                }
            }

            configuration = ASETools.Configuration.loadFromFile(args);

            if (null == configuration)
            {
                Console.WriteLine("Giving up because we were unable to load configuration.");
                return;
            }

            if (configuration.commandLineArgs.Count() != 0 && configuration.commandLineArgs.Count() != 1 || configuration.commandLineArgs.Count() == 1 && configuration.commandLineArgs[0] != "-c")
            {
                Console.WriteLine("usage: OverallDistribution {-c}");
                Console.WriteLine("-c means to generate corrected ASE values.");
                return;
            }

            bool useCorrection = configuration.commandLineArgs.Count() == 1 && configuration.commandLineArgs[0] == "-c";

            if (useCorrection)
            {
                var filename = configuration.finalResultsDirectory + ASETools.ASECorrectionFilename;
                aseCorrection = ASETools.ASECorrection.LoadFromFile(filename);
                if (null == aseCorrection)
                {
                    Console.WriteLine("Unable to load ASECorrection from " + filename);
                    return;
                }
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

            var outputFile = ASETools.CreateStreamWriterWithRetry(configuration.finalResultsDirectory + (useCorrection ? ASETools.CorrectedOverallASEFilename :  ASETools.UncorrectedOverallASEFilename));

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


            if (!useCorrection)
            {
                //
                // Write the tumor germline ASE distribution into its own file, beause it's needed programmatically elsewhere.
                //
                var tumorGermlineASEDistributionFile = ASETools.CreateStreamWriterWithRetry(configuration.finalResultsDirectory + ASETools.TumorGermlineASEDistributionFilename);
                tumorGermlineASEDistributionFile.WriteLine(ASETools.HistogramResultLine.Header());
                overallResult[true][false].ComputeHistogram().ToList().ForEach(x => tumorGermlineASEDistributionFile.WriteLine(x.ToString()));
                tumorGermlineASEDistributionFile.WriteLine("**done**");
                tumorGermlineASEDistributionFile.Close();

                //
                // Write the tumor RNA read depth into its own file, because it's needed programmatically elsewhere.
                //
                var tumorRNAReadDepthFile = ASETools.CreateStreamWriterWithRetry(configuration.finalResultsDirectory + ASETools.TumorRNAReadDepthDistributionFilename);
                tumorRNAReadDepthFile.WriteLine(ASETools.HistogramResultLine.Header());
                readDepthDistribution[true][false][false].ComputeHistogram().ToList().ForEach(x => tumorRNAReadDepthFile.WriteLine(x.ToString()));
                tumorRNAReadDepthFile.WriteLine("**done**");
                tumorRNAReadDepthFile.Close();
            }

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

            var allDNAReadDepthsByDisease = germlineTumorDNAReadDepthByDiseaseType.Select(x => new DiseaseAndHistogramLines(x.Key, x.Value)).ToList();
            allDNAReadDepthsByDisease.Sort();
            allDNAReadDepthsByDisease.Insert(0, new DiseaseAndHistogramLines("overall", readDepthDistribution[true][false][true]));

            outputFile.WriteLine("Tumor DNA read depth distribution at germline variant sites by disease type");
            outputFile.Write("minValue\tpdf");
            allDNAReadDepthsByDisease.ForEach(x => outputFile.Write("\t" + x.disease));
            outputFile.Write("\tcdf");
            allDNAReadDepthsByDisease.ForEach(x => outputFile.Write("\t" + x.disease));
            outputFile.Write("\tcount");
            allDNAReadDepthsByDisease.ForEach(x => outputFile.Write("\t" + x.disease));
            outputFile.WriteLine();

            for (int i = 0; i < allDNAReadDepthsByDisease[0].histogramLines.Count(); i++)
            {
                outputFile.Write(allDNAReadDepthsByDisease[0].histogramLines[i].minValue +"\t");    // Skip the "pdf" column, which is meant to be blank
                allDNAReadDepthsByDisease.ForEach(x => outputFile.Write("\t" + x.histogramLines[i].pdfValue));
                outputFile.Write("\t"); // Skip the "cdf" column, which is meant to be blank
                allDNAReadDepthsByDisease.ForEach(x => outputFile.Write("\t" + x.histogramLines[i].cdfValue));
                outputFile.Write("\t");// Skip the "count" column, which is meant to be blank
                allDNAReadDepthsByDisease.ForEach(x => outputFile.Write("\t" + x.histogramLines[i].count));
                outputFile.WriteLine();
            }

            var allCountOfASESitesByDisease = validGermlineASESitesByDisease.Select(x => new DiseaseAndHistogramLines(x.Key, x.Value)).ToList();
            allCountOfASESitesByDisease.Sort();
            allCountOfASESitesByDisease.Insert(0, new DiseaseAndHistogramLines("overall", overallValidGermlineASESites));

            outputFile.WriteLine("Distribution of valid gemline ASE sites at germline variant sites by disease type");
            outputFile.Write("minValue\tpdf");
            allCountOfASESitesByDisease.ForEach(x => outputFile.Write("\t" + x.disease));
            outputFile.Write("\tcdf");
            allCountOfASESitesByDisease.ForEach(x => outputFile.Write("\t" + x.disease));
            outputFile.Write("\tcount");
            allCountOfASESitesByDisease.ForEach(x => outputFile.Write("\t" + x.disease));
            outputFile.WriteLine();

            for (int i = 0; i < allCountOfASESitesByDisease[0].histogramLines.Count(); i++)
            {
                outputFile.Write(allCountOfASESitesByDisease[0].histogramLines[i].minValue + "\t");    // Skip the "pdf" column, which is meant to be blank
                allCountOfASESitesByDisease.ForEach(x => outputFile.Write("\t" + x.histogramLines[i].pdfValue));
                outputFile.Write("\t"); // Skip the "cdf" column, which is meant to be blank
                allCountOfASESitesByDisease.ForEach(x => outputFile.Write("\t" + x.histogramLines[i].cdfValue));
                outputFile.Write("\t");// Skip the "count" column, which is meant to be blank
                allCountOfASESitesByDisease.ForEach(x => outputFile.Write("\t" + x.histogramLines[i].count));
                outputFile.WriteLine();
            }

            outputFile.WriteLine("**done**");
            outputFile.Close();

            Console.WriteLine("Completed in " + ASETools.ElapsedTimeInSeconds(timer));
        } // Main

        class DiseaseAndHistogramLines : IComparable
        {
            public readonly string disease;
            public readonly ASETools.HistogramResultLine[] histogramLines;

            public DiseaseAndHistogramLines(string disease_, ASETools.PreBucketedHistogram histogram)
            {
                disease = disease_;
                histogramLines = histogram.ComputeHistogram();
            }

            public int CompareTo(object peerObject)
            {
                DiseaseAndHistogramLines peer = (DiseaseAndHistogramLines)peerObject;

                return disease.CompareTo(peer.disease);
            }
        }
    }
}
