using System;
using System.CodeDom.Compiler;
using System.Collections.Generic;
using System.Linq;
using System.Net;
using System.Runtime.InteropServices;
using System.Security.Cryptography;
using System.Text;
using System.Threading.Tasks;
using System.Xml.Linq;
using ASELib;

namespace AnalyzeMiRNAExpression
{
    class Program
    {
        static ASETools.CommonData commonData;
        static List<string> miRNANames;
        static List<string> geneNames;

        class miRNAAndMutationCounts
        {
            public miRNAAndMutationCounts(Dictionary<string, ASETools.miRNAExpressionQuantification> miRNA_, List<ASETools.AnnotatedVariant> annotatedVariants)
            {
                miRNA = miRNA_;
                foreach (var geneName in geneNames)
                {
                    var mutationIndex = ASETools.ZeroOneMany(annotatedVariants.Where(_ => _.Hugo_symbol == geneName && _.somaticMutation && !_.isSilent()).Count());
                    mutationCounts.Add(geneName, mutationIndex);
                } // gene
            }  // ctor

            public readonly Dictionary<string, int> mutationCounts = new Dictionary<string, int>();
            public readonly Dictionary<string, ASETools.miRNAExpressionQuantification> miRNA;
        }

        static Dictionary<string, miRNAAndMutationCounts> miRNAAndMutationCountsByCase = new Dictionary<string, miRNAAndMutationCounts>();

        static void LoadOneCase(ASETools.Case case_, int state)
        {
            var miRNAExpressionQuantification = ASETools.miRNAExpressionQuantification.LoadFromFile(case_.tumor_miRNA_expression_quantification_filename);
            var asv = ASETools.AnnotatedVariant.readFile(case_.annotated_selected_variants_filename);

            var data = new miRNAAndMutationCounts(miRNAExpressionQuantification, asv);

            lock (miRNAAndMutationCountsByCase)
            {
                miRNAAndMutationCountsByCase.Add(case_.case_id, data);
            }
        }

        static void HandleOneGene(string geneName, List<ResultAndPValue> state)
        {
            var tumorsByMutationCount = miRNAAndMutationCountsByCase.Keys.GroupByToDict(caseId => miRNAAndMutationCountsByCase[caseId].mutationCounts[geneName]);
            foreach (var mutationCount in ASETools.ZeroOneTwo)
            {
                if (!tumorsByMutationCount.ContainsKey(mutationCount))
                {
                    tumorsByMutationCount.Add(mutationCount, new List<string>());
                }
            }

            int[] countOfTumorsByNMutations = { tumorsByMutationCount[0].Count(), tumorsByMutationCount[1].Count(), tumorsByMutationCount[2].Count() };

            if (countOfTumorsByNMutations[1] < ASETools.Configuration.minInstancesOfEachClassToComputeMannWhitney ||
                countOfTumorsByNMutations[0] + countOfTumorsByNMutations[2] < ASETools.Configuration.minInstancesOfEachClassToComputeMannWhitney)
            {
                //
                // Hopeless for this gene.  Stop now.
                //
                return;
            }

            foreach (var miRNAName in miRNANames)
            {
                //
                // Make a list of the exactly one mutation and > 1 mutation cases
                //
                var expressionLevels = tumorsByMutationCount[1].Select(caseId => new ExpressionLevelAndClass(miRNAAndMutationCountsByCase[caseId].miRNA[miRNAName].reads_mer_milion_miRNA_mapped, true)).ToList();
                expressionLevels.AddRange(tumorsByMutationCount[2].Select(caseId => new ExpressionLevelAndClass(miRNAAndMutationCountsByCase[caseId].miRNA[miRNAName].reads_mer_milion_miRNA_mapped, false)).ToList());

                bool enoughData;
                bool reversed;
                double nInFirstGroup, nInSecondGroup, U, z, p;

                if (countOfTumorsByNMutations[2] >= ASETools.Configuration.minInstancesOfEachClassToComputeMannWhitney)
                {
                    //
                    // We have enough for the 1 vs. >1 comparison.
                    //
                    p = ASETools.MannWhitney<ExpressionLevelAndClass>.ComputeMannWhitney(expressionLevels, expressionLevels[0], a => a.oneMutation, a => a.expressionLevel,
                                                                                         out enoughData, out reversed, out nInFirstGroup, out nInSecondGroup, out U, out z, 
                                                                                         true, ASETools.Configuration.minInstancesOfEachClassToComputeMannWhitney);

                    if (!enoughData)
                    {
                        throw new Exception("HandleOneGene: MannWhitney said not enough data when we verified it going in.");
                    }

                    state.Add(new ResultAndPValue(geneName, miRNAName , false, p, expressionLevels.Where(_ => _.oneMutation).Select(_ => _.expressionLevel).Average(), countOfTumorsByNMutations[1],
                                                  expressionLevels.Where(_ => !_.oneMutation).Select(_ => _.expressionLevel).Average(), countOfTumorsByNMutations[2]));
                }
                continue; // BJB - just do 1 vs 2 for now.

                expressionLevels.AddRange(tumorsByMutationCount[0].Select(caseId => new ExpressionLevelAndClass(miRNAAndMutationCountsByCase[caseId].miRNA[miRNAName].reads_mer_milion_miRNA_mapped, false)).ToList());
                p = ASETools.MannWhitney<ExpressionLevelAndClass>.ComputeMannWhitney(expressionLevels, expressionLevels[0], a => a.oneMutation, a => a.expressionLevel,
                        out enoughData, out reversed, out nInFirstGroup, out nInSecondGroup, out U, out z, true, ASETools.Configuration.minInstancesOfEachClassToComputeMannWhitney);

                if (!enoughData)
                {
                    throw new Exception("HandleOneGene: MannWhitney said not enough data when we verified it going in.");
                }

                state.Add(new ResultAndPValue(geneName, miRNAName, true, p, expressionLevels.Where(_ => _.oneMutation).Select(_ => _.expressionLevel).Average(), countOfTumorsByNMutations[1], 
                                              expressionLevels.Where(_ => !_.oneMutation).Select(_ => _.expressionLevel).Average(), countOfTumorsByNMutations[0] + countOfTumorsByNMutations[2]));
            } // miRNA
        } // HandleOneGene

        static void FinishUp(List<ResultAndPValue> state)
        {
            lock (results)
            {
                results.AddRange(state);
            }
        } // FinishUp

        static List<ResultAndPValue> results = new List<ResultAndPValue>();

        static void Main(string[] args)
        {
            commonData = ASETools.CommonData.LoadCommonData(args);

            if (null == commonData)
            {
                return;
            }

            if (commonData.listOfCases.Any(_ => _.tumor_miRNA_expression_quantification_filename == "" && _.tumor_miRNA_expression_quantification_file_id != "" || _.annotated_selected_variants_filename == ""))
            {
                Console.WriteLine("At least one case that has tumor miRNA quantification is missing an input.  Download it and try again.");
                return;
            }

            var casesToRun = commonData.listOfCases.Where(_ => _.tumor_miRNA_expression_quantification_filename != "" && _.annotated_selected_variants_filename != "").ToList();
            if (casesToRun.Count() == 0)
            {
                Console.WriteLine("No cases have data (?)");
                return;
            }

            //
            // We need to initialize the list of miRNANames.  Read in a single case's miRNA quantification.
            //
            var miRNAExpressionQuantification = ASETools.miRNAExpressionQuantification.LoadFromFile(casesToRun[0].tumor_miRNA_expression_quantification_filename);
            miRNANames = miRNAExpressionQuantification.Select(_ => _.Key).ToList();
            geneNames = ASETools.SelectedGene.LoadFromFile(commonData.configuration.selectedGenesFilename).Select(_ => _.Hugo_Symbol).ToList();

            int nPerDot;
            ASETools.PrintMessageAndNumberBar("Loading data for", "cases", casesToRun.Count(), out nPerDot);

            var threading = new ASETools.WorkerThreadHelper<ASETools.Case, int>(casesToRun, LoadOneCase, null, null, nPerDot);

            threading.run();
            Console.WriteLine(ASETools.ElapsedTimeInSeconds(commonData.timer));

            ASETools.PrintMessageAndNumberBar("Processing", "genes", geneNames.Count(), out nPerDot);
            var threading2 = new ASETools.WorkerThreadHelper<string, List<ResultAndPValue>>(geneNames, HandleOneGene, FinishUp, null, nPerDot);

            threading2.run();
            Console.WriteLine("total time to now " + ASETools.ElapsedTimeInSeconds(commonData.timer));


            results.Sort((a, b) => a.pValue.CompareTo(b.pValue));

            int bonferroniCorrection = results.Count();

            Console.WriteLine("Found " + bonferroniCorrection + " miRNA, gene and comparison combinations with enough instances to compute a p Value.  " + results.Where(_ => _.pValue * bonferroniCorrection <= 0.01).Count() + " are significant at 0.01.");


            var outputFile = ASETools.CreateStreamWriterWithRetry(commonData.configuration.finalResultsDirectory + ASETools.miRNAExpressionSummaryFilename);
            if (null == outputFile)
            {
                Console.WriteLine("Unable to open " + commonData.configuration.finalResultsDirectory + ASETools.miRNAExpressionSummaryFilename);
                return;
            }

            outputFile.WriteLine("Hugo Symbol\tmiRNA Name\tOne vs Not One\tCorrected p\tn With One Mutation\tMean with One Mutation\tn In Other Class\tMean for other class\tRatio Of Larger to Smaller\tSignificant");
            foreach (var resultAndPValue in results)
            {
                outputFile.WriteLine(resultAndPValue.HugoSymbol + "\t" + resultAndPValue.miRNAName + "\t" + resultAndPValue.oneVsNotOne + "\t" + 
                    (resultAndPValue.pValue * bonferroniCorrection) + "\t" + resultAndPValue.nWithOneMutation + "\t" + resultAndPValue.meanWithOneMutation + "\t" + 
                    resultAndPValue.nInOtherClass + "\t" + resultAndPValue.meanForOtherClass + "\t" + 
                    Math.Max(resultAndPValue.meanWithOneMutation, resultAndPValue.meanForOtherClass) / Math.Min(resultAndPValue.meanWithOneMutation, resultAndPValue.meanForOtherClass) + "\t" +
                    (((resultAndPValue.pValue * bonferroniCorrection) <= 0.01) ? "Y" : "N"));
            }

            outputFile.WriteLine("**done**");
            outputFile.Close();

            ASETools.PreBucketedHistogram pValueHistogram = new ASETools.PreBucketedHistogram(0, 1, 0.01);
            results.ForEach(_ => pValueHistogram.addValue(_.pValue));

            var pValueHisogramFile = ASETools.CreateStreamWriterWithRetry(commonData.configuration.finalResultsDirectory + ASETools.miRNAExpressionPValueHistogramFilename);
            if (null != pValueHisogramFile)
            {
                pValueHistogram.WriteHistogram(pValueHisogramFile);
                pValueHisogramFile.WriteLine("**done**");
                pValueHisogramFile.Close();
            }


            Console.WriteLine(ASETools.ElapsedTimeInSeconds(commonData.timer));
        } // Main

        class ResultAndPValue 
        {
            public ResultAndPValue(string HugoSymbol_, string miRNAName_, bool oneVsNotOne_, double pValue_, double meanWithOneMutation_, int nWithOneMutation_, double meanForOtherClass_, int nInOtherClass_)
            {
                HugoSymbol = HugoSymbol_;
                miRNAName = miRNAName_;
                oneVsNotOne = oneVsNotOne_;
                pValue = pValue_;
                meanWithOneMutation = meanWithOneMutation_;
                nWithOneMutation = nWithOneMutation_;
                meanForOtherClass = meanForOtherClass_;
                nInOtherClass = nInOtherClass_;
            }

            public readonly string HugoSymbol;
            public readonly string miRNAName;
            public readonly bool oneVsNotOne;
            public readonly double pValue;
            public readonly double meanWithOneMutation;
            public readonly int nWithOneMutation;
            public readonly double meanForOtherClass;
            public readonly int nInOtherClass;
        } // ResultAndPValue

        class ExpressionLevelAndClass : IComparer<ExpressionLevelAndClass>
        {
            public ExpressionLevelAndClass(double expressionLevel_, bool oneMutation_)
            {
                expressionLevel = expressionLevel_;
                oneMutation = oneMutation_;
            }

            public int Compare(ExpressionLevelAndClass a, ExpressionLevelAndClass b)
            {
                return a.expressionLevel.CompareTo(b.expressionLevel);
            }

            public double expressionLevel;
            public bool oneMutation;
        } // ExpressionLevelAndClass

    }
}
