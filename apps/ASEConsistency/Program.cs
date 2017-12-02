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
        static Dictionary<bool, Dictionary<string, ASETools.RunningMeanAndStdDev>> perGeneRefFraction = new Dictionary<bool, Dictionary<string, ASETools.RunningMeanAndStdDev>>();
        static Dictionary<bool, Dictionary<int, ASETools.Histogram>> overallResults = new Dictionary<bool, Dictionary<int, ASETools.Histogram>>(); // tumor, min read depth
        static Dictionary<bool, Dictionary<int, ASETools.Histogram>> overallResultsWithCorrection = new Dictionary<bool, Dictionary<int, ASETools.Histogram>>(); // tumor, min read depth
        static Dictionary<bool, ASETools.Histogram> overallSameExonResults = new Dictionary<bool, ASETools.Histogram>();
        static Dictionary<bool, ASETools.Histogram> overallReadDepthDifference = new Dictionary<bool, ASETools.Histogram>();
        static Dictionary<bool, ASETools.Histogram> overallReadDepthDifferenceSameExon = new Dictionary<bool, ASETools.Histogram>();
        static Dictionary<bool, ASETools.Histogram> referenceFraction = new Dictionary<bool, ASETools.Histogram>();
        static ASETools.ASECorrection ASECorrection;

        static int[] readDepths = { 10, 25, 45, 110, 135, 250 };    // The values must be in ascending order for the code to work properly.

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

            foreach (var tumor in ASETools.BothBools)
            {
                overallResults.Add(tumor, new Dictionary<int, ASETools.Histogram>());
                overallResultsWithCorrection.Add(tumor, new Dictionary<int, ASETools.Histogram>());
                foreach (var readDepth in readDepths)
                {
                    overallResults[tumor].Add(readDepth, new ASETools.Histogram());
                    overallResultsWithCorrection[tumor].Add(readDepth, new ASETools.Histogram());
                }

                perGeneResults.Add(tumor, new Dictionary<string, ASETools.Histogram>());
                perGeneRefFraction.Add(tumor, new Dictionary<string, ASETools.RunningMeanAndStdDev>());
                referenceFraction.Add(tumor, new ASETools.Histogram("reference fraction: tumor " + tumor));
                overallSameExonResults.Add(tumor, new ASETools.Histogram());
                overallReadDepthDifference.Add(tumor, new ASETools.Histogram());
                overallReadDepthDifferenceSameExon.Add(tumor, new ASETools.Histogram());
            }

            geneLocationInformation = new ASETools.GeneLocationsByNameAndChromosome(ASETools.readKnownGeneFile(configuration.geneLocationInformationFilename));
            geneMap = new ASETools.GeneMap(geneLocationInformation.genesByName);

            perGeneASEMap = ASETools.ASEMapPerGeneLine.ReadFromFileToDictionary(configuration.finalResultsDirectory + ASETools.PerGeneASEMapFilename);

            if (null == perGeneASEMap)
            {
                Console.WriteLine("You must first create the per-gene ASE map in " + configuration.finalResultsDirectory + ASETools.PerGeneASEMapFilename);
                return;
            }

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
                outputFile.WriteLine("Read depth difference distribution (same exon), tumor: " + tumor);
                outputFile.WriteLine(ASETools.HistogramResultLine.Header());
                overallReadDepthDifferenceSameExon[tumor].ComputeHistogram(0, 200, 1).ToList().ForEach(x => outputFile.WriteLine(x.ToString()));
            }

            foreach (var tumor in ASETools.BothBools)
            {
                outputFile.WriteLine("Read depth difference distribution, tumor: " + tumor);
                outputFile.WriteLine(ASETools.HistogramResultLine.Header());
                overallReadDepthDifference[tumor].ComputeHistogram(0, 200, 1).ToList().ForEach(x => outputFile.WriteLine(x.ToString()));
            }

            foreach (var tumor in ASETools.BothBools)
            {
                foreach (var readDepth in readDepths)
                {
                    outputFile.WriteLine("Overall results: tumor " + tumor + ", min read depth: " + readDepth);
                    outputFile.WriteLine(ASETools.HistogramResultLine.Header());
                    overallResults[tumor][readDepth].ComputeHistogram(0, 1, 0.01).ToList().ForEach(x => outputFile.WriteLine(x.ToString()));

                    outputFile.WriteLine("Overall results with correction: tumor " + tumor + ", min read depth: " + readDepth);
                    outputFile.WriteLine(ASETools.HistogramResultLine.Header());
                    overallResultsWithCorrection[tumor][readDepth].ComputeHistogram(0, 1, 0.01).ToList().ForEach(x => outputFile.WriteLine(x.ToString()));
                }
            }

            foreach (var tumor in ASETools.BothBools)
            {
                outputFile.WriteLine("Overall result (same exon): tumor " + tumor);
                outputFile.WriteLine(ASETools.HistogramResultLine.Header());
                overallSameExonResults[tumor].ComputeHistogram(0, 1, 0.01).ToList().ForEach(x => outputFile.WriteLine(x.ToString()));
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
                var extremeRefFraction = perGeneRefFraction[tumor].Where(x => x.Value.getCount() >= 20 && (x.Value.getMeanAndStdDev().mean > 0.7 || x.Value.getMeanAndStdDev().mean < 0.3)).ToList();
                extremeRefFraction.Sort((x, y) => y.Value.getMeanAndStdDev().mean.CompareTo(x.Value.getMeanAndStdDev().mean));// NB: Backward comparison to put the highest mean ones first.

                outputFile.WriteLine("Genes with min 20 variants and ref fraction > 0.7, tumor: " + tumor);
                outputFile.WriteLine("Gene\tn\tReference Fraction");
                extremeRefFraction.ForEach(x => outputFile.WriteLine(x.Key + "\t" + x.Value.getCount() + "\t" + x.Value.getMeanAndStdDev().mean));
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

            Console.WriteLine("Took " + ASETools.ElapsedTimeInSeconds(timer));

        } // Main

        struct ASEAndLocus
        {
            public readonly int locus;
            public readonly double ASE;

            public ASEAndLocus(int locus_, double ASE_)
            {
                locus = locus_;
                ASE = ASE_;
            }
        }

        static void WorkerThread(List<ASETools.Case> casesToProcess)
        {
            var localPerGeneResults = new Dictionary<bool, Dictionary<string, ASETools.Histogram>>();
            var localPerGeneRefFraction = new Dictionary<bool, Dictionary<string, ASETools.RunningMeanAndStdDev>>();
            var localOverallResults = new Dictionary<bool, Dictionary<int, ASETools.Histogram>>();
            var localOverallResultsWithCorrection = new Dictionary<bool, Dictionary<int, ASETools.Histogram>>();
            var localReferenceFraction = new Dictionary<bool, ASETools.Histogram>();
            var localOverallSameExonResults = new Dictionary<bool, ASETools.Histogram>();
            var localOverallReadDepthDifference = new Dictionary<bool, ASETools.Histogram>();
            var localOverallReadDepthDifferenceSameExon = new Dictionary<bool, ASETools.Histogram>();

            foreach (var tumor in ASETools.BothBools)
            {
                localPerGeneResults.Add(tumor, new Dictionary<string, ASETools.Histogram>());
                localPerGeneRefFraction.Add(tumor, new Dictionary<string, ASETools.RunningMeanAndStdDev>());
                localOverallResults.Add(tumor, new Dictionary<int, ASETools.Histogram>());
                localOverallResultsWithCorrection.Add(tumor, new Dictionary<int, ASETools.Histogram>());
                foreach (var readDepth in readDepths)
                {
                    localOverallResults[tumor].Add(readDepth, new ASETools.Histogram());
                    localOverallResultsWithCorrection[tumor].Add(readDepth, new ASETools.Histogram());
                }
                localReferenceFraction.Add(tumor, new ASETools.Histogram());
                localOverallSameExonResults.Add(tumor, new ASETools.Histogram());
                localOverallReadDepthDifference.Add(tumor, new ASETools.Histogram());
                localOverallReadDepthDifferenceSameExon.Add(tumor, new ASETools.Histogram());
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

                            foreach (var localPerGeneRefFractionEntry in localPerGeneRefFraction[tumor])
                            {
                                var hugo_symbol = localPerGeneRefFractionEntry.Key;
                                var localResult = localPerGeneRefFractionEntry.Value;

                                if (!perGeneRefFraction[tumor].ContainsKey(hugo_symbol))
                                {
                                    perGeneRefFraction[tumor].Add(hugo_symbol, localResult);
                                }
                                else
                                {
                                    perGeneRefFraction[tumor][hugo_symbol].merge(localResult);
                                }
                            }

                            foreach (var readDepth in readDepths)
                            {
                                overallResults[tumor][readDepth].merge(localOverallResults[tumor][readDepth]);
                                overallResultsWithCorrection[tumor][readDepth].merge(localOverallResultsWithCorrection[tumor][readDepth]);
                            }
                            referenceFraction[tumor].merge(localReferenceFraction[tumor]);
                            overallSameExonResults[tumor].merge(localOverallSameExonResults[tumor]);
                            overallReadDepthDifference[tumor].merge(localOverallReadDepthDifference[tumor]);
                            overallReadDepthDifferenceSameExon[tumor].merge(localOverallReadDepthDifference[tumor]);
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


                var copyNumber = ASETools.CopyNumberVariation.ReadBothFiles(case_);
                var annotatedSelectedVariants = ASETools.AnnotatedVariant.readFile(case_.annotated_selected_variants_filename);

                var map = new Dictionary<bool, Dictionary<int, Dictionary<string, List<ASEAndLocus>>>>(); // tumor/normal->min read depth->gene->ASEs
                var mapWithCorrection = new Dictionary<bool, Dictionary<int, Dictionary<string, List<ASEAndLocus>>>>(); // tumor/normal->min read depth->gene->ASEs
                var readDepthMap = new Dictionary<bool, Dictionary<string, List<int>>>(); // tumor->gene name->read depth
                foreach (var tumor in ASETools.BothBools)
                {
                    map.Add(tumor, new Dictionary<int, Dictionary<string, List<ASEAndLocus>>>());
                    mapWithCorrection.Add(tumor, new Dictionary<int, Dictionary<string, List<ASEAndLocus>>>());
                    foreach (var readDepth in readDepths)
                    {
                        map[tumor].Add(readDepth, new Dictionary<string, List<ASEAndLocus>>());
                        mapWithCorrection[tumor].Add(readDepth, new Dictionary<string, List<ASEAndLocus>>());
                    }
                    readDepthMap.Add(tumor, new Dictionary<string, List<int>>());
                }

                foreach (var variant in annotatedSelectedVariants)
                {
                    if (variant.somaticMutation)
                    {
                        continue;
                    }

                    foreach (var tumor in ASETools.BothBools)
                    {
                        foreach (var readDepth in readDepths)
                        {
                            if (variant.IsASECandidate(tumor, copyNumber, configuration, perGeneASEMap, geneMap, readDepth))
                            {
                                if (readDepth == readDepths[0])
                                {
                                    localReferenceFraction[tumor].addValue(1 - variant.GetAltAlleleFraction(tumor));
                                }

                                if (geneMap.getGenesMappedTo(variant.contig, variant.locus).Count() != 1)   // Skip loci in more than one gene (or in none)
                                {
                                    continue;
                                }

                                foreach (var gene in geneMap.getGenesMappedTo(variant.contig, variant.locus))
                                {
                                    if (!map[tumor][readDepth].ContainsKey(gene.hugoSymbol))
                                    {
                                        map[tumor][readDepth].Add(gene.hugoSymbol, new List<ASEAndLocus>());
                                        mapWithCorrection[tumor][readDepth].Add(gene.hugoSymbol, new List<ASEAndLocus>());
                                    }

                                    map[tumor][readDepth][gene.hugoSymbol].Add(new ASEAndLocus(variant.locus, variant.GetAlleleSpecificExpression(tumor)));
                                    mapWithCorrection[tumor][readDepth][gene.hugoSymbol].Add(new ASEAndLocus(variant.locus, variant.GetAlleleSpecificExpression(tumor, ASECorrection)));

                                    if (readDepth == readDepths[0])
                                    {
                                        if (!readDepthMap[tumor].ContainsKey(gene.hugoSymbol))
                                        {
                                            readDepthMap[tumor].Add(gene.hugoSymbol, new List<int>());
                                        }

                                        readDepthMap[tumor][gene.hugoSymbol].Add(variant.getReadCount(tumor, false).usefulReads());

                                        if (!localPerGeneRefFraction[tumor].ContainsKey(gene.hugoSymbol))
                                        {
                                            localPerGeneRefFraction[tumor].Add(gene.hugoSymbol, new ASETools.RunningMeanAndStdDev());
                                        }

                                        localPerGeneRefFraction[tumor][gene.hugoSymbol].addValue(1 - variant.GetAltAlleleFraction(tumor));
                                    }
                                } // foreach mapped gene
                            } // if it's an ASE candidate
                        }// read depth
                    } // tumor in BothBools
                } // foreach variant

                foreach (var tumor in ASETools.BothBools)
                {
                    foreach (var gene in map[tumor][readDepths[0]])
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

                        var spread = gene.Value.Select(x => x.ASE).Max() - gene.Value.Select(x => x.ASE).Min();
                        localPerGeneResults[tumor][hugo_symbol].addValue(spread);

                        localOverallReadDepthDifference[tumor].addValue(readDepthMap[tumor][hugo_symbol].Max() - readDepthMap[tumor][hugo_symbol].Min());

                        //
                        // See if we have any results in the same exon.
                        //
                        double maxPerExonValue = -1;
                        foreach (var isoform in geneLocationInformation.genesByName[hugo_symbol].isoforms)
                        {
                            foreach (ASETools.Exon exon in isoform.exons)
                            {
                                var matchingResults = gene.Value.Where(x => x.locus >= exon.start && x.locus <= exon.end).ToList();

                                if (matchingResults.Count() < 2)
                                {
                                    continue;
                                }

                                maxPerExonValue = Math.Max(maxPerExonValue, matchingResults.Select(x => x.ASE).Max() - matchingResults.Select(x => x.ASE).Min());
                            }
                        }

                        if (maxPerExonValue >= 0)
                        {
                            localOverallSameExonResults[tumor].addValue(maxPerExonValue);
                        }

                        foreach (var readDepth in readDepths)
                        {
                            if (!map[tumor][readDepth].ContainsKey(hugo_symbol) || map[tumor][readDepth][hugo_symbol].Count() < 2)
                            {
                                break;
                            }

                            var spreadAtThisReadDepth = map[tumor][readDepth][hugo_symbol].Select(x => x.ASE).Max() - map[tumor][readDepth][hugo_symbol].Select(x => x.ASE).Min();
                            var spreadAtThisReadDepthWithCorrection = mapWithCorrection[tumor][readDepth][hugo_symbol].Select(x => x.ASE).Max() - mapWithCorrection[tumor][readDepth][hugo_symbol].Select(x => x.ASE).Min();

                            localOverallResults[tumor][readDepth].addValue(spreadAtThisReadDepth);
                            localOverallResultsWithCorrection[tumor][readDepth].addValue(spreadAtThisReadDepthWithCorrection);
                        }
                    } // foreach gene
                } // foreach tumor in bothbools
            } // while (true)
        }
    } // Program
} // ASEConsistency
