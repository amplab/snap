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

        static Dictionary<string, Dictionary<string, ASETools.ExpressionByGeneLine>> geneExpressionByCase = new Dictionary<string, Dictionary<string, ASETools.ExpressionByGeneLine>>();
        static Dictionary<string, ASETools.MeanAndStdDev> geneExpressionMeanAndStandardDeviation = new Dictionary<string, ASETools.MeanAndStdDev>();
        static Dictionary<string, ASETools.SelectedGene> selectedGenes;
        static ASETools.CommonData commonData;
        struct NAndSum
        {
            public int n;
            public double sum;
        }

        static Dictionary<string, Dictionary<int, Dictionary<string, NAndSum>>> correlationState = new Dictionary<string, Dictionary<int, Dictionary<string, NAndSum>>>(); // contig->locus->gene->state.  Locking is per locus.

        static void Main(string[] args)
        {
            commonData = ASETools.CommonData.LoadCommonData(args);
            if (null == commonData)
            {
                return;
            }

            selectedGenes = ASETools.SelectedGene.LoadFromFile(commonData.configuration.selectedGenesFilename).GroupByToDictUnique(_ => _.Hugo_Symbol);

            var casesToProcess = commonData.cases.Select(_ => _.Value).Where(_ => _.tumor_methylation_filename != "").ToList();

           foreach (var hugo_symbol in selectedGenes.Select(_ => _.Key))
            {
                var runningValue = new ASETools.RunningMeanAndStdDev();
                geneExpressionByCase.Select(_ => _.Value[hugo_symbol].fractionOfMeanExpression).ToList().ForEach(_ => runningValue.addValue(_));

                geneExpressionMeanAndStandardDeviation.Add(hugo_symbol, runningValue.getMeanAndStdDev());
            }

            //
            // Load one methylation so we can initialize the map.
            //
            foreach (var case_ in casesToProcess)
            {
                var methylations = ASETools.MethylationAnnotationLine.ReadFile(case_.tumor_methylation_filename);
                if (methylations.Count() < 9000)
                {
                    //
                    // There are two methylation array sizes.  We want the big one.  This isn't it.
                    //
                    continue;
                }

                foreach (var methylation in methylations)
                {
                    var contig = methylation.compositeREF.Chromosome;
                    int locus = methylation.compositeREF.Start;

                    if (!correlationState.ContainsKey(contig))
                    {
                        correlationState.Add(contig, new Dictionary<int, Dictionary<string, NAndSum>>());
                    }

                    if (!correlationState[contig].ContainsKey(locus))
                    {
                        correlationState[contig].Add(locus, new Dictionary<string, NAndSum>());
                        var thisLocus = correlationState[contig][locus];
                        foreach (var gene in selectedGenes)
                        {
                            var state = new NAndSum();
                            state.n = 0;
                            state.sum = 0;
                            thisLocus.Add(gene.Key, state);
                        } // gene
                    } // locus
                } // contig
                break;
            } // cases

            if (correlationState.Count() == 0)
            {
                Console.WriteLine("No cases had a big methylation array.");
                return;
            }
            Console.WriteLine("Loading expression by gene for " + casesToProcess.Count() + " cases, 1 dot/100 cases:");
            ASETools.PrintNumberBar(casesToProcess.Count() / 100);

            var threading = new ASETools.WorkerThreadHelper<ASETools.Case, int>(casesToProcess, LoadOneGeneExpression, null, null, 100);
            threading.run();
            Console.WriteLine("Elapsed time " + ASETools.ElapsedTimeInSeconds(commonData.timer) + ".  Computing distributions for " + selectedGenes.Count() + " genes.");



            Console.WriteLine("Elapsed time " + ASETools.ElapsedTimeInSeconds(commonData.timer) + ". Processing methylation data for " +
                casesToProcess.Count() + " cases, 1 dot/100 cases.");
            ASETools.PrintNumberBar(casesToProcess.Count() / 100);

            threading = new ASETools.WorkerThreadHelper<ASETools.Case, int>(casesToProcess, LoadOneMethylation, null, null, 100);
            threading.run();

            Console.WriteLine("Elapsed time " + ASETools.ElapsedTimeInSeconds(commonData.timer) + ".  Writing output.");

        } // Main

        static void LoadOneGeneExpression(ASETools.Case case_, int state)
        {
            var geneExpression = ASETools.ExpressionByGeneLine.ReadFromFile(case_.gene_expression_filename);

            lock (geneExpressionByCase)
            {
                geneExpressionByCase.Add(case_.case_id, geneExpression);
            }
        } // LoadOneGeneExpression

        static void LoadOneMethylation(ASETools.Case case_, int state)
        {
            var methylation = ASETools.MethylationAnnotationLine.ReadFile(case_.tumor_methylation_filename, case_.case_id, false);

            foreach 
        }
    }
}
