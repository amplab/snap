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
        static Dictionary<string, int[,]> scatterChartByGene = new Dictionary<string, int[,]>();
        static Dictionary<string, int[,]> scatterChartByGeneSameExon = new Dictionary<string, int[,]>();
        static ASETools.PreBucketedHistogram overallSomaticGermlineDifference = new ASETools.PreBucketedHistogram(-1, 1, 0.01);
        static ASETools.PreBucketedHistogram overallSomaticGermlineDifferenceSameExon = new ASETools.PreBucketedHistogram(-1, 1, 0.01);
        static Dictionary<bool, ASETools.PreBucketedHistogram> overallSomaticGermlineDifferenceByWhetherGeneIsInterestingForLowVAF = new Dictionary<bool, ASETools.PreBucketedHistogram>();
        static ASETools.ASECorrection ASECorrection;
        static ASETools.CommonData commonData;
        static Dictionary<string, int> genesWithHighSomaticAndLowGermlineASE = new Dictionary<string, int> ();
        static Dictionary<string, List<string>> interestingGeneSomaticGermlinePairs = new Dictionary<string, List<string>>();

        static int[] readDepths = { 10, 25, 45, 110, 135, 250 };    // The values must be in ascending order for the code to work properly.

        const int nScatterGraphChunks = 101; // 0 -> 1 in units of 0.01

        static int[,] SummarizeScatterCharts(Dictionary<string, int[,]>scatterCharts)
        {
            return SummarizeScatterCharts(scatterCharts.ToList());
        } // SummarizeScatterCharts

        static int[,] SummarizeScatterCharts(List<KeyValuePair<string, int[,]>>scatterChartEntries)
        {
            var output = new int[nScatterGraphChunks, nScatterGraphChunks];

            foreach (var scatterChart in scatterChartEntries.Select(_ => _.Value))
            {
                for (int i = 0; i < nScatterGraphChunks; i++)
                {
                    for (int j = 0; j < nScatterGraphChunks; j++)
                    {
                        output[i, j] += scatterChart[i, j];
                    }
                }
            }

            return output;
        }

        static void WriteScatterChart(StreamWriter output, int [,]scatterChart)
        {
            output.Write("SomaticASE/GermlineASE");
            for (int j = 0; j < nScatterGraphChunks; j++)
            {
                output.Write("\t" + (double)j / (nScatterGraphChunks - 1));
            }
            output.WriteLine();

            for (int i = 0; i < nScatterGraphChunks; i++)
            {
                output.Write((double)i / (nScatterGraphChunks - 1));
                for (int j = 0; j < nScatterGraphChunks; j++)
                {
                    output.Write("\t" + scatterChart[i, j]);
                } // j
                output.WriteLine();
            }// i
        } // WriteScatterChart

        static List<string> interestingGenesForLowVAF = new List<string>();

        static void Main(string[] args)
        {
            commonData = ASETools.CommonData.LoadCommonData(args);
            if (null == commonData)
            {
                return;
            }

            interestingGenesForLowVAF.Add("PRKCSH");
            interestingGenesForLowVAF.Add("COL3A1");
            interestingGenesForLowVAF.Add("COL1A2");
            interestingGenesForLowVAF.Add("COL6A3");
            interestingGenesForLowVAF.Add("COL6A2");
            interestingGenesForLowVAF.Add("COL5A2");
            interestingGenesForLowVAF.Add("COL1A1");
            interestingGenesForLowVAF.Add("VWF");
            interestingGenesForLowVAF.Add("COL4A1");
            interestingGenesForLowVAF.Add("CSTA");
            interestingGenesForLowVAF.Add("VPS41");

            foreach (var interesting in ASETools.BothBools)
            {
                overallSomaticGermlineDifferenceByWhetherGeneIsInterestingForLowVAF.Add(interesting, new ASETools.PreBucketedHistogram(-1, 1, 0.01));
            }

            configuration = commonData.configuration;
            geneLocationInformation = commonData.geneLocationInformation;
            geneMap = commonData.geneMap;
            perGeneASEMap = commonData.perGeneASEMap;
            if (null == perGeneASEMap)
            {
                Console.WriteLine("You must first create the per-gene ASE map in " + configuration.finalResultsDirectory + ASETools.PerGeneASEMapFilename);
                return;
            }

            ASECorrection = commonData.aseCorrection;
            if (null == ASECorrection)
            {
                Console.WriteLine("Unable to load ASE correction");
                return;
            }
            var cases = commonData.cases;

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

            var casesToProcess = cases.Select(x => x.Value).Where(x => x.annotated_selected_variants_filename != "" && x.tumor_copy_number_filename != "").ToList();

            int nCasesPerDot;
            ASETools.PrintMessageAndNumberBar("Processing", "cases", casesToProcess.Count(), out nCasesPerDot);
            var threading = new ASETools.WorkerThreadHelper<ASETools.Case, WorkerThreadState>(casesToProcess, HandleOneCase, FinishUp, null, nCasesPerDot);
            threading.run();

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
                    outputFile.WriteLine(ASETools.HistogramResultLine.Header());
                    perGeneHistogram.ComputeHistogram(0, 1, 0.01).ToList().ForEach(x => outputFile.WriteLine(x.ToString()));
                }
            }

            outputFile.WriteLine();
            outputFile.WriteLine("Genes with at least 0.9 somatic ASE and no more than 0.1 germline ASE in the same tumor");
            outputFile.WriteLine("Hugo Symbol\tCount");
            foreach (var entry in genesWithHighSomaticAndLowGermlineASE)
            {
                outputFile.WriteLine(entry.Key + "\t" + entry.Value);
            }

            outputFile.WriteLine();
            outputFile.WriteLine("Distribution of somaticASE - germlineASE for variants in the same gene");
            overallSomaticGermlineDifference.WriteHistogram(outputFile);

            outputFile.WriteLine();
            outputFile.WriteLine("Distribution of somaticASE - germlineASE for variants in the same exon");
            overallSomaticGermlineDifferenceSameExon.WriteHistogram(outputFile);

            var summaryScatterGraph = SummarizeScatterCharts(scatterChartByGene);
            outputFile.WriteLine("Summary scatter chart for all genes.");
            WriteScatterChart(outputFile, summaryScatterGraph);

            foreach (var interesting in ASETools.BothBools)
            {
                var scatterGraph = SummarizeScatterCharts(scatterChartByGene.Where(_ => interestingGenesForLowVAF.Contains(_.Key) == interesting).ToList());
                outputFile.WriteLine();
                outputFile.WriteLine("Summary scatter chart for all genes that are " + (interesting ? "" : "not ") + " on the list of genes with aberrantly low somatic VAF");
                WriteScatterChart(outputFile, scatterGraph);

                outputFile.WriteLine();
                outputFile.WriteLine("Distribution of somatic ASE - germlineASE for variants in the same gene that " + (interesting ? "" : "not ") + " on the list of genes with aberrantly low somatic VAF");
                overallSomaticGermlineDifferenceByWhetherGeneIsInterestingForLowVAF[interesting].WriteHistogram(outputFile);
            }

            var sameExonSummaryGraph = SummarizeScatterCharts(scatterChartByGeneSameExon);
            outputFile.WriteLine();
            outputFile.WriteLine("Summary scatter chart for all genes with somatic and germline variants in the same exon");
            WriteScatterChart(outputFile, sameExonSummaryGraph);

            int nPointsToEmitScatterChart = 20;

 


            foreach (var hugo_symbol in scatterChartByGene.Select(_ => _.Key))
            {
                if (scatterChartByGene[hugo_symbol].Cast<int>().Sum() >= nPointsToEmitScatterChart || interestingGenesForLowVAF.Contains(hugo_symbol))
                {
                    outputFile.WriteLine();
                    outputFile.WriteLine("Scatter chart for " + hugo_symbol + " regardless of exon");
                    WriteScatterChart(outputFile, scatterChartByGene[hugo_symbol]);

                    if (scatterChartByGeneSameExon.ContainsKey(hugo_symbol) && (scatterChartByGeneSameExon[hugo_symbol].Cast<int>().Sum() >= nPointsToEmitScatterChart || interestingGenesForLowVAF.Contains(hugo_symbol)))
                    {
                        outputFile.WriteLine();
                        outputFile.WriteLine("Scatter chart for " + hugo_symbol + " all points in same exon");
                        WriteScatterChart(outputFile, scatterChartByGeneSameExon[hugo_symbol]);
                    }
                }
            }

            foreach (var hugo_symbol in interestingGeneSomaticGermlinePairs.Select(_ => _.Key))
            {
                outputFile.WriteLine();
                outputFile.WriteLine(hugo_symbol + " cases with both somatic and germline ASE measurements in the same tumor");

                interestingGeneSomaticGermlinePairs[hugo_symbol].ForEach(_ => outputFile.WriteLine(_));
            }

            outputFile.WriteLine("**done**");
            outputFile.Close();

            Console.WriteLine("Took " + ASETools.ElapsedTimeInSeconds(commonData.timer));

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

        class WorkerThreadState
        {
            public Dictionary<bool, Dictionary<string, ASETools.Histogram>> perGeneResults = new Dictionary<bool, Dictionary<string, ASETools.Histogram>>();
            public Dictionary<bool, Dictionary<string, ASETools.RunningMeanAndStdDev>> perGeneRefFraction = new Dictionary<bool, Dictionary<string, ASETools.RunningMeanAndStdDev>>();
            public Dictionary<bool, Dictionary<int, ASETools.Histogram>> overallResults = new Dictionary<bool, Dictionary<int, ASETools.Histogram>>();
            public Dictionary<bool, Dictionary<int, ASETools.Histogram>> overallResultsWithCorrection = new Dictionary<bool, Dictionary<int, ASETools.Histogram>>();
            public Dictionary<bool, ASETools.Histogram> referenceFraction = new Dictionary<bool, ASETools.Histogram>();
            public Dictionary<bool, ASETools.Histogram> overallSameExonResults = new Dictionary<bool, ASETools.Histogram>();
            public Dictionary<bool, ASETools.Histogram> overallReadDepthDifference = new Dictionary<bool, ASETools.Histogram>();
            public Dictionary<bool, ASETools.Histogram> overallReadDepthDifferenceSameExon = new Dictionary<bool, ASETools.Histogram>();
            public Dictionary<string, int[,]> scatterChartByGene = new Dictionary<string, int[,]>();
            public Dictionary<string, int[,]> scatterChartByGeneSameExon = new Dictionary<string, int[,]>();
            public ASETools.PreBucketedHistogram overallSomaticGermlineDifference = new ASETools.PreBucketedHistogram(-1,1,0.01);
            public ASETools.PreBucketedHistogram overallSomaticGermlineDifferenceSameExon = new ASETools.PreBucketedHistogram(-1, 1, 0.01);
            public Dictionary<bool, ASETools.PreBucketedHistogram> overallSomaticGermlineDifferenceByWhetherGeneIsInterestingForLowVAF = new Dictionary<bool, ASETools.PreBucketedHistogram>();
            public Dictionary<string, int> genesWithHighSomaticAndLowGermlineASE = new Dictionary<string, int>();
            public Dictionary<string, List<string>> interestingGeneSomaticGermlinePairs = new Dictionary<string, List<string>>();

            public WorkerThreadState()
            {
                foreach (var tumor in ASETools.BothBools)
                {
                    perGeneResults.Add(tumor, new Dictionary<string, ASETools.Histogram>());
                    perGeneRefFraction.Add(tumor, new Dictionary<string, ASETools.RunningMeanAndStdDev>());
                    overallResults.Add(tumor, new Dictionary<int, ASETools.Histogram>());
                    overallResultsWithCorrection.Add(tumor, new Dictionary<int, ASETools.Histogram>());
                    foreach (var readDepth in readDepths)
                    {
                        overallResults[tumor].Add(readDepth, new ASETools.Histogram());
                        overallResultsWithCorrection[tumor].Add(readDepth, new ASETools.Histogram());
                    }
                    referenceFraction.Add(tumor, new ASETools.Histogram());
                    overallSameExonResults.Add(tumor, new ASETools.Histogram());
                    overallReadDepthDifference.Add(tumor, new ASETools.Histogram());
                    overallReadDepthDifferenceSameExon.Add(tumor, new ASETools.Histogram());
                } // tumor/normal

                foreach (var interesting in ASETools.BothBools)
                {
                    overallSomaticGermlineDifferenceByWhetherGeneIsInterestingForLowVAF.Add(interesting, new ASETools.PreBucketedHistogram(-1, 1, 0.01));
                }
            } // ctor
        } // WorkerThreadState

        static void mergeScatterGraphs(Dictionary<string,int[,]> into, Dictionary<string,int [,]> from)
        {
            foreach (var hugo_symbol in from.Select(_ => _.Key)) {
                if (!into.ContainsKey(hugo_symbol))
                {
                    into.Add(hugo_symbol, new int[nScatterGraphChunks, nScatterGraphChunks]);
                }

                for (int i = 0; i < nScatterGraphChunks; i++)
                {
                    for (int j = 0; j < nScatterGraphChunks; j++)
                    {
                        into[hugo_symbol][i, j] += from[hugo_symbol][i, j];
                    } // j
                } // i
            } // hugo symbol
        } // mergeScatterGraphs

        static void FinishUp(WorkerThreadState state)
        {
            lock (perGeneResults)
            {
                mergeScatterGraphs(scatterChartByGene, state.scatterChartByGene);
                mergeScatterGraphs(scatterChartByGeneSameExon, state.scatterChartByGeneSameExon);
                overallSomaticGermlineDifference.merge(state.overallSomaticGermlineDifference);
                overallSomaticGermlineDifferenceSameExon.merge(overallSomaticGermlineDifferenceSameExon);

                foreach (var tumor in ASETools.BothBools)
                {
                    foreach (var resultEntry in state.perGeneResults[tumor])
                    {
                        var hugo_symbol = resultEntry.Key;
                        var result = resultEntry.Value;

                        if (!perGeneResults[tumor].ContainsKey(hugo_symbol))
                        {
                            perGeneResults[tumor].Add(hugo_symbol, result);
                        }
                        else
                        {
                            perGeneResults[tumor][hugo_symbol].merge(result);
                        }
                    }

                    foreach (var perGeneRefFractionEntry in state.perGeneRefFraction[tumor])
                    {
                        var hugo_symbol = perGeneRefFractionEntry.Key;
                        var result = perGeneRefFractionEntry.Value;

                        if (!perGeneRefFraction[tumor].ContainsKey(hugo_symbol))
                        {
                            perGeneRefFraction[tumor].Add(hugo_symbol, result);
                        }
                        else
                        {
                            perGeneRefFraction[tumor][hugo_symbol].merge(result);
                        }
                    }

                    foreach (var readDepth in readDepths)
                    {
                        overallResults[tumor][readDepth].merge(state.overallResults[tumor][readDepth]);
                        overallResultsWithCorrection[tumor][readDepth].merge(state.overallResultsWithCorrection[tumor][readDepth]);
                    }
                    referenceFraction[tumor].merge(state.referenceFraction[tumor]);
                    overallSameExonResults[tumor].merge(state.overallSameExonResults[tumor]);
                    overallReadDepthDifference[tumor].merge(state.overallReadDepthDifference[tumor]);
                    overallReadDepthDifferenceSameExon[tumor].merge(state.overallReadDepthDifference[tumor]);
                } // tumor/normal

                foreach (var interesting in ASETools.BothBools)
                {
                    overallSomaticGermlineDifferenceByWhetherGeneIsInterestingForLowVAF[interesting].merge(state.overallSomaticGermlineDifferenceByWhetherGeneIsInterestingForLowVAF[interesting]);
                }

                foreach (var entry in state.genesWithHighSomaticAndLowGermlineASE)
                {
                    if (!genesWithHighSomaticAndLowGermlineASE.ContainsKey(entry.Key))
                    {
                        genesWithHighSomaticAndLowGermlineASE.Add(entry.Key, entry.Value);
                    } else
                    {
                        genesWithHighSomaticAndLowGermlineASE[entry.Key] += entry.Value;
                    }
                }

                foreach (var hugo_symbol in state.interestingGeneSomaticGermlinePairs.Select(_ => _.Key))
                {
                    if (!interestingGeneSomaticGermlinePairs.ContainsKey(hugo_symbol))
                    {
                        interestingGeneSomaticGermlinePairs.Add(hugo_symbol, state.interestingGeneSomaticGermlinePairs[hugo_symbol]);
                    } else
                    {
                        interestingGeneSomaticGermlinePairs[hugo_symbol].AddRange(state.interestingGeneSomaticGermlinePairs[hugo_symbol]);
                    }
                }


            } // lock
        } // FinishUp

        static void HandleOneCase(ASETools.Case case_, WorkerThreadState state)
        {

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

            var genesAlreadyInScatterChart = new Dictionary<string, int>(); // The value isn't meaningful, it's just to keep track of which genes we've seen

            var genesForGermineVariants = new Dictionary<ASETools.AnnotatedVariant, List<string>>();
            foreach (var variant in annotatedSelectedVariants.Where(_ => !_.somaticMutation && _.IsASECandidate(true, copyNumber, commonData)).ToList())  // germline variants don't have 
            {
                genesForGermineVariants.Add(variant, commonData.geneMap.getGenesMappedTo(variant.contig, variant.locus).Select(_ => _.hugoSymbol).ToList());
            }

            foreach (var variant in annotatedSelectedVariants)
            {
                if (variant.somaticMutation && variant.IsASECandidate(true, copyNumber, commonData) && 
                    !genesAlreadyInScatterChart.ContainsKey(variant.Hugo_symbol) && 
                    annotatedSelectedVariants.Any(_ => !_.somaticMutation && _.IsASECandidate(true, copyNumber, commonData) && genesForGermineVariants[_].Contains(variant.Hugo_symbol)))
                {
                    //
                    // This one has (at least) two variants in the same gene covering both somatic and germline variants.  It's a candidate for the scatter chart.
                    //
                    var hugo_symbol = variant.Hugo_symbol;
                    genesAlreadyInScatterChart.Add(hugo_symbol, 42);    // Value doesn't matter, this is just to keep track of the genes we've already done
                    var somaticVariants = annotatedSelectedVariants.Where(_ => _.somaticMutation && _.Hugo_symbol == hugo_symbol && _.IsASECandidate(true, copyNumber, commonData)).ToList();
                    var germlineVariants = annotatedSelectedVariants.Where(_ => !_.somaticMutation && _.IsASECandidate(true, copyNumber, commonData) && genesForGermineVariants[_].Contains(hugo_symbol)).ToList();

                    if (!state.scatterChartByGene.ContainsKey(hugo_symbol))
                    {
                        state.scatterChartByGene.Add(hugo_symbol, new int[nScatterGraphChunks, nScatterGraphChunks]);
                    }

                    var somaticASE = somaticVariants.Select(_ => _.GetAlleleSpecificExpression(true, commonData.aseCorrection)).Average();
                    var germlineASE = germlineVariants.Select(_ => _.GetAlleleSpecificExpression(true, commonData.aseCorrection)).Average();

                    if (somaticASE >= 0.9 && germlineASE <= 0.1)
                    {  if (!state.genesWithHighSomaticAndLowGermlineASE.ContainsKey(hugo_symbol))
                        {
                            state.genesWithHighSomaticAndLowGermlineASE.Add(hugo_symbol, 0);
                        }

                        state.genesWithHighSomaticAndLowGermlineASE[hugo_symbol]++;
                    }

                    state.overallSomaticGermlineDifference.addValue(somaticASE - germlineASE);
                    state.overallSomaticGermlineDifferenceByWhetherGeneIsInterestingForLowVAF[interestingGenesForLowVAF.Contains(hugo_symbol)].addValue(somaticASE - germlineASE);

                    state.scatterChartByGene[hugo_symbol][(int)(somaticASE * (nScatterGraphChunks - 1)), (int)(germlineASE * (nScatterGraphChunks - 1))]++;

                    int distanceInSameExon = -1;

                    foreach (var isoform in commonData.geneLocationInformation.genesByName[hugo_symbol].isoforms)
                    {
                        foreach (var exon in isoform.exons)
                        {
                            var somaticInThisExon = somaticVariants.Where(_ => exon.start <= _.locus && exon.end >= _.locus).ToList();
                            var germlineInThisExon = germlineVariants.Where(_ => exon.start <= _.locus && exon.end >= _.locus).ToList();

                            if (somaticInThisExon.Count() > 0 && germlineInThisExon.Count() > 0)
                            {
                                distanceInSameExon = Math.Abs(somaticInThisExon[0].locus - germlineInThisExon[0].locus);    // Technically, they could be closer, but in practice there's rarely more than one of each, so this is close enough.

                                var somaticASEInThisExon = somaticInThisExon.Select(_ => _.GetAlleleSpecificExpression(true, commonData.aseCorrection)).Average();
                                var germlineASEInThisExon = germlineInThisExon.Select(_ => _.GetAlleleSpecificExpression(true, commonData.aseCorrection)).Average();

                                if (!state.scatterChartByGeneSameExon.ContainsKey(hugo_symbol))
                                {
                                    state.scatterChartByGeneSameExon.Add(hugo_symbol, new int[nScatterGraphChunks, nScatterGraphChunks]);
                                }

                                state.scatterChartByGeneSameExon[hugo_symbol][(int)(somaticASEInThisExon * (nScatterGraphChunks - 1)), (int)(germlineASEInThisExon * (nScatterGraphChunks - 1))]++;
                                state.overallSomaticGermlineDifferenceSameExon.addValue(somaticASEInThisExon - germlineASEInThisExon);

                                //
                                // Now remove them from the list, so that if there are overlapping exons (or isoforms) that contain them that we won't double count.
                                //
                                foreach (var variantToRemove in somaticInThisExon)
                                {
                                    somaticVariants.Remove(variantToRemove);
                                }

                                foreach (var variantToRemove in germlineInThisExon)
                                {
                                    germlineVariants.Remove(variantToRemove);
                                }
                            }  // If we have both somatic and germline variants in this exon
                        } // for each exon in this isoform
                    } // for each isoform 

                    if (interestingGenesForLowVAF.Contains(hugo_symbol))
                    {
                        if (!state.interestingGeneSomaticGermlinePairs.ContainsKey(hugo_symbol))
                        {
                            state.interestingGeneSomaticGermlinePairs.Add(hugo_symbol, new List<string>());
                        }

                        state.interestingGeneSomaticGermlinePairs[hugo_symbol].Add(case_.case_id + " somatic ASE " + somaticASE + " germline ASE " + germlineASE + ((distanceInSameExon != -1) ? " in same exon at distance " + distanceInSameExon : ""));
                    }
                } // If we have both germline and somatic variants in this gene

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
                                state.referenceFraction[tumor].addValue(1 - variant.GetAltAlleleFraction(tumor));
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

                                    if (!state.perGeneRefFraction[tumor].ContainsKey(gene.hugoSymbol))
                                    {
                                        state.perGeneRefFraction[tumor].Add(gene.hugoSymbol, new ASETools.RunningMeanAndStdDev());
                                    }

                                    state.perGeneRefFraction[tumor][gene.hugoSymbol].addValue(1 - variant.GetAltAlleleFraction(tumor));
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
                    if (!state.perGeneResults[tumor].ContainsKey(hugo_symbol))
                    {
                        state.perGeneResults[tumor].Add(hugo_symbol, new ASETools.Histogram(hugo_symbol + " tumor: " + tumor));
                    }

                    var spread = gene.Value.Select(x => x.ASE).Max() - gene.Value.Select(x => x.ASE).Min();
                    state.perGeneResults[tumor][hugo_symbol].addValue(spread);

                    state.overallReadDepthDifference[tumor].addValue(readDepthMap[tumor][hugo_symbol].Max() - readDepthMap[tumor][hugo_symbol].Min());

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
                        state.overallSameExonResults[tumor].addValue(maxPerExonValue);
                    }

                    foreach (var readDepth in readDepths)
                    {
                        if (!map[tumor][readDepth].ContainsKey(hugo_symbol) || map[tumor][readDepth][hugo_symbol].Count() < 2)
                        {
                            break;
                        }

                        var spreadAtThisReadDepth = map[tumor][readDepth][hugo_symbol].Select(x => x.ASE).Max() - map[tumor][readDepth][hugo_symbol].Select(x => x.ASE).Min();
                        var spreadAtThisReadDepthWithCorrection = mapWithCorrection[tumor][readDepth][hugo_symbol].Select(x => x.ASE).Max() - mapWithCorrection[tumor][readDepth][hugo_symbol].Select(x => x.ASE).Min();

                        state.overallResults[tumor][readDepth].addValue(spreadAtThisReadDepth);
                        state.overallResultsWithCorrection[tumor][readDepth].addValue(spreadAtThisReadDepthWithCorrection);
                    }
                } // foreach gene
            } // foreach tumor in bothbools 
        } // HandleOneCase

    } // Program
} // ASEConsistency
