using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using ASELib;
using System.Diagnostics;
using System.Threading;
using System.IO;

namespace ASEConsistency
{
    class Program
    {
        static ASETools.Configuration configuration;
        static ASETools.GeneMap geneMap;
        static Dictionary<string, ASETools.ASEMapPerGeneLine> perGeneASEMap;
        static ASETools.GeneLocationsByNameAndChromosome geneLocationInformation;
        static Dictionary<bool, Dictionary<string, ASETools.Histogram>> perGeneResults = new Dictionary<bool, Dictionary<string, ASETools.Histogram>>();
        static Dictionary<bool, ASETools.Histogram> overallResults = new Dictionary<bool, ASETools.Histogram>();
        static Dictionary<bool, ASETools.Histogram> referenceFraction = new Dictionary<bool, ASETools.Histogram>();

        static StreamWriter badInstancesFile;
        static StreamWriter questionableInstancesFile;

        static int nProcessed = 0;

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

            badInstancesFile = ASETools.CreateStreamWriterWithRetry("ReallyInconsistentASELocations.txt");
            questionableInstancesFile = ASETools.CreateStreamWriterWithRetry("KindaInconsistentASELocations.txt");
            string header = "Case ID\tgene\tchrom\tmin locus\tmax locus\ttumor";
            badInstancesFile.WriteLine(header);
            questionableInstancesFile.WriteLine(header);

            foreach (var tumor in ASETools.BothBools)
            {
                overallResults.Add(tumor, new ASETools.Histogram());
                perGeneResults.Add(tumor, new Dictionary<string, ASETools.Histogram>());
                referenceFraction.Add(tumor, new ASETools.Histogram("reference fraction: tumor " + tumor));
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

            var threads = new List<Thread>();
            for (int i = 0; i < Environment.ProcessorCount; i++)
            {
                threads.Add(new Thread(() => WorkerThread(casesToProcess)));
            }

            threads.ForEach(t => t.Start());
            threads.ForEach(t => t.Join());

            Console.WriteLine();

            var outputFile = ASETools.CreateStreamWriterWithRetry(configuration.finalResultsDirectory +  ASETools.ASEConsistencyFilename);

            foreach (var tumor in ASETools.BothBools)
            {
                outputFile.WriteLine("Overall results: tumor " + tumor);
                outputFile.WriteLine(ASETools.HistogramResultLine.Header());
                overallResults[tumor].ComputeHistogram(0, 1, 0.01).ToList().ForEach(x => outputFile.WriteLine(x.ToString()));
            }

            foreach (var tumor in ASETools.BothBools)
            {
                outputFile.WriteLine();
                outputFile.WriteLine("Reference fraction: tumor " + tumor);
                outputFile.WriteLine(ASETools.HistogramResultLine.Header());
                referenceFraction[tumor].ComputeHistogram(0, 1, 0.01).ToList().ForEach(x => outputFile.WriteLine(x.ToString()));
            }

            foreach (var tumor in ASETools.BothBools)
            {
                var interestingPerGene = perGeneResults[tumor].Select(x => x.Value).Where(x => x.count() >= 100).ToList();
                interestingPerGene.Sort((x, y) => y.mean().CompareTo(x.mean()));    // NB: Backward comparison to put the highest mean ones first.

                foreach (var perGeneHistogram in interestingPerGene)
                {
                    outputFile.WriteLine();
                    outputFile.WriteLine(perGeneHistogram.name);
                    perGeneHistogram.ComputeHistogram(0, 1, 0.01).ToList().ForEach(x => outputFile.WriteLine(x.ToString()));
                }
            }

            outputFile.WriteLine("**done**");
            outputFile.Close();

            badInstancesFile.WriteLine("**done**");
            badInstancesFile.Close();

            questionableInstancesFile.WriteLine("**done**");
            questionableInstancesFile.Close();

            Console.WriteLine("Took " + ASETools.ElapsedTimeInSeconds(timer));

        } // Main

        static void WorkerThread(List<ASETools.Case> casesToProcess)
        {
            var localPerGeneResults = new Dictionary<bool, Dictionary<string, ASETools.Histogram>>();
            var localOverallResults = new Dictionary<bool, ASETools.Histogram>();
            var localReferenceFraction = new Dictionary<bool, ASETools.Histogram>();

            foreach (var tumor in ASETools.BothBools)
            {
                localPerGeneResults.Add(tumor, new Dictionary<string, ASETools.Histogram>());
                localOverallResults.Add(tumor, new ASETools.Histogram());
                localReferenceFraction.Add(tumor, new ASETools.Histogram());
            }

            while (true)
            {
                ASETools.Case case_;

                lock (casesToProcess)
                {
                    if (casesToProcess.Count() == 0)
                    {
                        foreach (var tumor in ASETools.BothBools)
                        {
                            foreach (var localResultEntry in localPerGeneResults[tumor])
                            {
                                var hugo_symbol = localResultEntry.Key;
                                var localResult = localResultEntry.Value;

                                if (!perGeneResults[tumor].ContainsKey(hugo_symbol))
                                {
                                    perGeneResults[tumor].Add(hugo_symbol, localResult);
                                }
                                else
                                {
                                    perGeneResults[tumor][hugo_symbol].merge(localResult);
                                }
                            }

                            overallResults[tumor].merge(localOverallResults[tumor]);
                            referenceFraction[tumor].merge(localReferenceFraction[tumor]);
                        } // tumor/normal
                        return;
                    }

                    case_ = casesToProcess[0];
                    casesToProcess.RemoveAt(0);
                    nProcessed++;
                    if (nProcessed % 100 == 0)
                    {
                        Console.Write(".");
                    }
                } // lock (casesToProcess)

                var annotatedSelectedVariants = ASETools.AnnotatedVariant.readFile(case_.annotated_selected_variants_filename);
                var tumorCopyNumberVariation = ASETools.CopyNumberVariation.ReadFile(case_.tumor_copy_number_filename, case_.tumor_copy_number_file_id).Where(r => Math.Abs(r.Segment_Mean) > 1.0).ToList();
                List<ASETools.CopyNumberVariation> normalCopyNumberVariation = null;
                if (case_.normal_copy_number_filename != "")
                {
                    normalCopyNumberVariation = ASETools.CopyNumberVariation.ReadFile(case_.normal_copy_number_filename, case_.normal_copy_number_file_id).Where(r => Math.Abs(r.Segment_Mean) > 1.0).ToList();
                }

                var map = new Dictionary<bool, Dictionary<string, List<double>>>(); // tumor/normal->gene->ASEs
                map.Add(false, new Dictionary<string, List<double>>());
                map.Add(true, new Dictionary<string, List<double>>());

                foreach (var variant in annotatedSelectedVariants)
                {
                    if (variant.somaticMutation)
                    {
                        continue;
                    }

                    foreach (var tumor in ASETools.BothBools)
                    {
                        if (variant.IsASECandidate(tumor, tumorCopyNumberVariation, configuration, perGeneASEMap, geneMap))
                        {
                            localReferenceFraction[tumor].addValue(1 - variant.GetAltAlleleFraction(tumor));

                            if (geneMap.getGenesMappedTo(variant.contig, variant.locus).Count() != 1)   // Skip loci in more than one gene (or in none)
                            {
                                continue;
                            }

                            foreach (var gene in geneMap.getGenesMappedTo(variant.contig, variant.locus)) {
                                if (!map[tumor].ContainsKey(gene.hugoSymbol))
                                {
                                    map[tumor].Add(gene.hugoSymbol, new List<double>());
                                }

                                map[tumor][gene.hugoSymbol].Add(variant.GetAlleleSpecificExpression(tumor));
                            } // foreach mapped gene


                        } // if it's an ASE candidate
                    } // tumor in BothBools
                } // foreach variant

                foreach (var tumor in ASETools.BothBools)
                {
                    foreach (var gene in map[tumor])
                    {
                        if (gene.Value.Count() < 2)
                        {
                            continue;
                        }

                        var hugo_symbol = gene.Key;
                        if (!localPerGeneResults[tumor].ContainsKey(hugo_symbol))
                        {
                            localPerGeneResults[tumor].Add(hugo_symbol, new ASETools.Histogram(hugo_symbol + " tumor: " + tumor));
                        }

                        var spread = gene.Value.Max() - gene.Value.Min();
                        localPerGeneResults[tumor][hugo_symbol].addValue(spread);

                        localOverallResults[tumor].addValue(spread);

                        StreamWriter outputFile = null;
                        if (spread >= 0.98)
                        {
                            outputFile = badInstancesFile;
                        } else if (spread >= 0.330 && spread <= 0.332)
                        {
                            outputFile = questionableInstancesFile;
                        }

                        if (null != outputFile) { 

                            lock (outputFile)
                            {
                                // "Case ID\tgene\tchrom\tmin locus\tmax locus\ttumor"
                                outputFile.WriteLine(case_.case_id + "\t" + hugo_symbol + "\t" + geneLocationInformation.genesByName[hugo_symbol].chromosome + "\t" + geneLocationInformation.genesByName[hugo_symbol].minLocus + "\t" +
                                    geneLocationInformation.genesByName[hugo_symbol].maxLocus + "\t" + tumor);
                            }
                        } 
                    } // foreach gene
                } // foreach tumor in bothbools
            } // while (true)
        }
    } // Program
} // ASEConsistency
