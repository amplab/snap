using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using ASELib;
using System.Diagnostics;
using System.IO;

namespace MethylationAndExpressionCorrelation
{
    class Program
    {

        static Dictionary<string, double[]> geneExpressionByCase = new Dictionary<string, double[]>();
        static Dictionary<string, ASETools.MeanAndStdDev> geneExpressionMeanAndStandardDeviation = new Dictionary<string, ASETools.MeanAndStdDev>();
        static string[] selectedGenes;
        static int nSelectedGenes;
        static ASETools.MethylationDistribution methylationDistribution;
        static ASETools.CommonData commonData;
        static string chromosomeToProcess = null;

        static Dictionary<string, Dictionary<int, double[]>> correlationState = new Dictionary<string, Dictionary<int, double[]>>(); // contig->locus->gene index->running total of product of differences from the mean.  Locking is per locus.

        static void Main(string[] args)
        {
            commonData = ASETools.CommonData.LoadCommonData(args);
            if (null == commonData)
            {
                return;
            }

            if (commonData.configuration.commandLineArgs.Count() != 0)
            {
                if (commonData.configuration.commandLineArgs.Count() > 1)
                {
                    Console.WriteLine("usage: MethylationAndExpressionCorrelation <chromosome>");
                    return;
                }

                chromosomeToProcess = commonData.configuration.commandLineArgs[0];
            }

            var casesToProcess = commonData.cases.Select(_ => _.Value).Where(_ => _.tumor_methylation_filename != "" && _.expression_by_gene_filename != "").ToList()./*BJB*/GetRange(0, 100).ToList();
            if (casesToProcess.Count() != commonData.cases.Where(_ => _.Value.tumor_methylation_file_id != "").Count()) // Don't complain about ones that don't have methylation data in TCGA.
            {
                Console.WriteLine("**warning: only " + casesToProcess.Count() + " cases of " + commonData.cases.Count() + " have the requisite data to process.");
                if (casesToProcess.Count() == 0)
                {
                    return;
                }
            }

            selectedGenes = ASETools.ExpressionByGene.LoadFromFile(casesToProcess[0].expression_by_gene_filename).GetAllHugoSymbols().ToArray();
            nSelectedGenes = selectedGenes.Length;

            for (int geneIndex = 0; geneIndex < nSelectedGenes; geneIndex++)
            {
                var runningValue = new ASETools.RunningMeanAndStdDev();
                geneExpressionByCase.Select(_ => _.Value[geneIndex]).ToList().ForEach(_ => runningValue.addValue(_));

                geneExpressionMeanAndStandardDeviation.Add(selectedGenes[geneIndex], runningValue.getMeanAndStdDev());
            }

            //
            // Load one methylation so we can initialize the map.
            //
            bool foundBig = false;
            bool foundSmall = false;
            ulong nMethylationSites = 0;

            foreach (var case_ in casesToProcess)
            {
                var methylations = ASETools.MethylationAnnotationLine.ReadFile(case_.tumor_methylation_filename);
                var big = methylations.Count() > 400000;
 
                if (big && foundBig || !big && foundSmall)
                {
                    continue;
                }


                foreach (var methylation in methylations.Where(_ =>  (chromosomeToProcess == null || chromosomeToProcess ==  _.compositeREF.Chromosome)))
                {
                    var contig = methylation.compositeREF.Chromosome;
                    int locus = methylation.compositeREF.Start;

                    if (!correlationState.ContainsKey(contig))
                    {
                        correlationState.Add(contig, new Dictionary<int, double[]>());
                    }

                    if (!correlationState[contig].ContainsKey(locus))
                    {
                        correlationState[contig].Add(locus, new double[nSelectedGenes]);
                        var thisLocus = correlationState[contig][locus];
                        nMethylationSites++;
                    } // locus
                } // methylation

                if (big)
                {
                    foundBig = true;
                } else
                {
                    foundSmall = true;
                }

                if (foundSmall && foundBig)
                {
                    break;
                }
                Console.WriteLine();
            } // cases

            if (!(foundBig && foundSmall))
            {
                Console.WriteLine("Didn't find both a big and a small methylation array.");
                return;
            }

            Console.WriteLine();
            Console.WriteLine("Total of " + nMethylationSites + " methylation sites, " + selectedGenes.Count() + " genes and " + casesToProcess.Count() +" cases, yields " + ASETools.SizeToUnits(3 * nMethylationSites * (ulong)selectedGenes.Count() * (ulong)casesToProcess.Count()) + " floating point ops (more-or-less; not all cases have all methylation sites). Elapsed time " + ASETools.ElapsedTimeInSeconds(commonData.timer)); // 3x is because there's a subtraction, multiplication and addition for each one.

            methylationDistribution = ASETools.MethylationDistribution.ReadFromFile(commonData.configuration.baseDirectory + ASETools.MethylationDistributionsFilename);
            if (null == methylationDistribution)
            {
                Console.WriteLine("Unable to read methylation distributions.");
                return;
            }

            Console.WriteLine("Loading expression by gene for " + casesToProcess.Count() + " cases, 1 dot/100 cases:");
            ASETools.PrintNumberBar(casesToProcess.Count() / 100);

            var threading = new ASETools.WorkerThreadHelper<ASETools.Case, int>(casesToProcess, LoadOneGeneExpression, null, null, 100);
            threading.run();

            var standardDeviationsOfGeneExpressions = new Dictionary<string, ASETools.RunningMeanAndStdDev>();
            foreach (var hugo_symbol in selectedGenes)
            {
                standardDeviationsOfGeneExpressions.Add(hugo_symbol, new ASETools.RunningMeanAndStdDev());
            }

            foreach (var caseId in casesToProcess.Select(_ => _.case_id))
            {
                for (int geneIndex = 0; geneIndex < nSelectedGenes; geneIndex++)
                {
                    standardDeviationsOfGeneExpressions[selectedGenes[geneIndex]].addValue(geneExpressionByCase[caseId][geneIndex]);
                }
            }

            Console.WriteLine("Elapsed time " + ASETools.ElapsedTimeInSeconds(commonData.timer) + ". Processing methylation data for " +
                casesToProcess.Count() + " cases, 1 dot/100 cases.");
            ASETools.PrintNumberBar(casesToProcess.Count() / 100);

            threading = new ASETools.WorkerThreadHelper<ASETools.Case, int>(casesToProcess, LoadOneMethylation, null, null, 100);
            threading.run();

            Console.WriteLine("Elapsed time " + ASETools.ElapsedTimeInSeconds(commonData.timer) + ".  Writing output.");

            var outputFilename = commonData.configuration.baseDirectory + ASETools.MethylationCorrelationsFilename + (chromosomeToProcess == null ? "" : "." + chromosomeToProcess);
            var outputFile = ASETools.CreateStreamWriterWithRetry(outputFilename);
            outputFile.WriteLine("Contig\tlocus\thugo symbol\tcorrelation coefficient");

            var histogram = new ASETools.PreBucketedHistogram(-1, 1, 0.01);

            foreach (var contigEntry in correlationState)
            {
                var contig = contigEntry.Key;

                foreach (var locusEntry in contigEntry.Value)
                {
                    var locus = locusEntry.Key;
                    var sumsOfCovariants = locusEntry.Value;    // Array indexed by geneIndex

                    for (int geneIndex = 0; geneIndex < nSelectedGenes; geneIndex++)
                    {
                        var hugoSymbol = selectedGenes[geneIndex];
                        var sumOfCovariants = sumsOfCovariants[geneIndex];

                        var meanOfCovariants = sumOfCovariants / methylationDistribution.forLocus(contig, locus).n;

                        var correlationCoefficient = meanOfCovariants / (methylationDistribution.forLocus(contig, locus).stdDev * standardDeviationsOfGeneExpressions[hugoSymbol].getMeanAndStdDev().stddev);

                        outputFile.WriteLine(contig + "\t" + locus + "\t" + hugoSymbol + "\t" + correlationCoefficient);

                        if (!double.IsNaN(correlationCoefficient))
                        {
                            histogram.addValue(correlationCoefficient);
                        }
                    } // gene
                } // locus
            } // contig

            outputFile.WriteLine("**done**");
            outputFile.Close();

            var histogramFilename = commonData.configuration.baseDirectory + ASETools.MethylagtionCorrelationsHistogramFilename + (chromosomeToProcess == null ? "" : "." + chromosomeToProcess);
            var histogramFile = ASETools.CreateStreamWriterWithRetry(histogramFilename);
            histogramFile.WriteLine(ASETools.HistogramResultLine.Header());
            histogram.ComputeHistogram().ToList().ForEach(_ => histogramFile.WriteLine(_));
            histogramFile.WriteLine("**done**");
            histogramFile.Close();

            Console.WriteLine("Total elapsed time " + ASETools.ElapsedTimeInSeconds(commonData.timer));

        } // Main

        static void LoadOneGeneExpression(ASETools.Case case_, int state)
        {
            var geneExpression = ASETools.ExpressionByGene.LoadFromFile(case_.expression_by_gene_filename);

            lock (geneExpressionByCase)
            {
                geneExpressionByCase.Add(case_.case_id, new double[nSelectedGenes]);

                var thisCase = geneExpressionByCase[case_.case_id];

                for (int geneIndex = 0; geneIndex < nSelectedGenes; geneIndex++)
                {
                    thisCase[geneIndex] = geneExpression.getMeanExpressionForGene(selectedGenes[geneIndex]);
                }
            }
        } // LoadOneGeneExpression

        static void LoadOneMethylation(ASETools.Case case_, int state)
        {
            var methylations = ASETools.MethylationAnnotationLine.ReadFile(case_.tumor_methylation_filename).Where(_ => (chromosomeToProcess == null || chromosomeToProcess == _.compositeREF.Chromosome)).ToList().OrderRandomly().ToList();
            var caseData = geneExpressionByCase[case_.case_id].ToList();

            foreach (var methylation in methylations)   // Random order avoids convoy effects on the lock.
            {
                var methylationDistributionLine = methylationDistribution.forLocus(methylation.compositeREF.Chromosome, methylation.compositeREF.Start);
                var methylationDifferenceFromMean = methylation.Beta_Value - methylationDistributionLine.mean;

                var locusMap = correlationState[methylation.compositeREF.Chromosome][methylation.compositeREF.Start];

                lock (locusMap)
                {
                    for (int geneIndex = 0; geneIndex < nSelectedGenes; geneIndex++)
                    {
                        var hugo_symbol = selectedGenes[geneIndex];

                        var meanExpressionForThisCase = caseData[geneIndex];
                        locusMap[geneIndex] += (meanExpressionForThisCase - 1) * methylationDifferenceFromMean; /* The -1 is because the expression is in units of the mean expression for that gene, so the mean of all of them is 1*/
                    } // gene index
                } // lock
            } // methylation
        } // LoadOneMethylation
    } // Program
} // namespace
