using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using System.IO;
using System.Diagnostics;
using ASELib;

namespace OverallGeneExpressionByMutationCount
{
    class Program
    {
        static ASETools.CommonData commonData;

        class GeneData
        {
            public Dictionary<string, ASETools.PreBucketedHistogram[]> geneDataByGene = new Dictionary<string, ASETools.PreBucketedHistogram[]>();
 
            public void merge(GeneData peer)
            {
                foreach (var geneHistogramsEntry in peer.geneDataByGene)
                {
                    var hugo_symbol = geneHistogramsEntry.Key;
                    var geneHistograms = geneHistogramsEntry.Value;

                    if (!geneDataByGene.ContainsKey(hugo_symbol))
                    {
                        geneDataByGene.Add(hugo_symbol, geneHistograms);
                    }
                    else
                    {
                        foreach (var i in histogramTypes)
                        {
                            geneDataByGene[hugo_symbol][i].merge(geneHistograms[i]);
                        }
                    }
                }
            } // merge
        } // GeneData

        static GeneData overallGeneData = new GeneData();
        const int silentMutationHistogramType = 3;
        static int[] histogramTypes = { 0, 1, 2, 3 }; // 0 mutations, 1 mutation, >1 mutation, 1 silent mutation
        static int nHistogramTypes = histogramTypes.Count();

        static void Main(string[] args)
        {
            commonData = ASETools.CommonData.LoadCommonData(args);

            if (null == commonData)
            {
                return;
            }

            if (commonData.configuration.commandLineArgs.Count() != 0)
            {
                Console.WriteLine("usage: OverallGeneExpressionByMutationCount");
                return;
            }

            if (commonData.cases.Any(x => x.Value.expression_by_gene_filename == "" || x.Value.annotated_selected_variants_filename == "")) {
                Console.WriteLine("At least one case doesn't have an expression by gene or annotated selected variants file.");
                return;
            }

            var casesToProcess = commonData.cases.Select(x => x.Value).Where(x => x.expression_by_gene_filename != "" && x.annotated_selected_variants_filename != "").ToList();

            Console.WriteLine("Processing " + casesToProcess.Count() + " cases, 1 dot/100:");
            ASETools.PrintNumberBar(casesToProcess.Count() / 100);

            var threading = new ASETools.ASVThreadingHelper<GeneData>(casesToProcess, commonData.aseCorrection, (x, y) => true, HandleOneCase, MergeIntoGlobal, null, 100);
            threading.run();

            //
            // Now generate the output file.
            //
            var outputFile = ASETools.CreateStreamWriterWithRetry(commonData.configuration.finalResultsDirectory + ASETools.Configuration.PerGeneExpressionHistogramsFilename);
            if (null == outputFile)
            {
                Console.WriteLine("Unable to open output file " + commonData.configuration.finalResultsDirectory + ASETools.Configuration.PerGeneExpressionHistogramsFilename);
                return;
            }

            foreach (var entry in overallGeneData.geneDataByGene)
            {
                var hugo_symbol = entry.Key;
                var histograms = entry.Value;

                outputFile.WriteLine(hugo_symbol);
                outputFile.WriteLine("Bucket\t0 mutations (n = " + histograms[0].count() + ")\t1 mutation (excluding NMD, n = " + histograms[1].count() + ")\t>1 mutation (excluding NMD, n = " + histograms[2].count() 
                    + ")\t1 silent mutation and no others (n = " + histograms[silentMutationHistogramType].count() + ")");

                var histogramLines = new List<ASETools.HistogramResultLine[]>();
                for (int i = 0; i < nHistogramTypes; i++)
                {
                    histogramLines.Add(histograms[i].ComputeHistogram());
                }

                for (int histogramBucket = 0; histogramBucket < histogramLines[0].Count(); histogramBucket++)
                {
                    outputFile.Write(histogramLines[0][histogramBucket].minValue);

                    foreach (int whichHistogram in histogramTypes)
                    {
                        outputFile.Write("\t" + histogramLines[whichHistogram][histogramBucket].cdfValue);
                    }

                    outputFile.WriteLine();
                } // foreach line in the histogram.

                outputFile.WriteLine(); // Blank line between genes.
                 
            } // for each gene

            outputFile.WriteLine("**done**");
            outputFile.Close();

            Console.WriteLine(ASETools.ElapsedTimeInSeconds(commonData.timer));

        } // Main

        static void HandleOneCase(ASETools.Case case_, GeneData state, List<ASETools.AnnotatedVariant> annotatedSelectedVariants, ASETools.ASECorrection aseCorrection, Dictionary<bool, List<ASETools.CopyNumberVariation>> copyNumber)
        {
            var expressionByGene = ASETools.ExpressionByGeneLine.ReadFromFile(case_.expression_by_gene_filename);

            if (null == expressionByGene)
            {
                Console.WriteLine("Unable to load expression by gene data for " + case_.expression_by_gene_filename);
                throw new Exception();
            }

            var asvByGene = annotatedSelectedVariants.Where(x => x.somaticMutation).GroupByToDict(x => x.Hugo_symbol);

            foreach (var entry in expressionByGene)
            {
                var hugo_symbol = entry.Key;
                var fractionOfMeanExpression = entry.Value.fractionOfMeanExpression;

                if (double.IsNaN(fractionOfMeanExpression))  // == double.NaN doens't work.
                {
                    continue;   // ignore this gene.
                }

                int nonSilentMutationCount;
                int silentMutationCount;
                if (asvByGene.ContainsKey(hugo_symbol))
                {
                    if (asvByGene[hugo_symbol].Any(x => x.CausesNonsenseMediatedDecay()))
                    {
                        // Just skip these for obvious reasons.
                        continue;
                    }

                    nonSilentMutationCount = asvByGene[hugo_symbol].Where(x => !x.isSilent()).Count();
                    silentMutationCount = asvByGene[hugo_symbol].Where(x => x.isSilent()).Count();
                } else
                {
                    nonSilentMutationCount = 0;
                    silentMutationCount = 0;
                }
                
                if (!state.geneDataByGene.ContainsKey(hugo_symbol))
                {
                    state.geneDataByGene.Add(hugo_symbol, new ASETools.PreBucketedHistogram[nHistogramTypes]);

                    foreach (var i in histogramTypes)
                    {
                        state.geneDataByGene[hugo_symbol][i] = new ASETools.PreBucketedHistogram(0, 3, 0.025);
                    }
                } // If this is the first time for this gene on this thread.

                state.geneDataByGene[hugo_symbol][ASETools.ZeroOneMany(nonSilentMutationCount)].addValue(fractionOfMeanExpression);
                if (silentMutationCount == 1 && nonSilentMutationCount == 0)
                {
                    state.geneDataByGene[hugo_symbol][silentMutationHistogramType].addValue(fractionOfMeanExpression);
                }
            } // foreach gene
        } // HandleOneCase

        static void MergeIntoGlobal(GeneData state) // Thread ended, merge its state into the global state
        {
            lock (overallGeneData)
            {
                overallGeneData.merge(state);
            }
        }
    }
}
