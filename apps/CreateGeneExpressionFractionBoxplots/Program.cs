using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using ASELib;
using System.IO;
using System.CodeDom.Compiler;

namespace CreateGeneExpressionFractionBoxplots
{
    internal class Program
    {
        static ASETools.CommonData commonData;
        static Dictionary<string, List<double>> result = new Dictionary<string, List<double>>();

        static void Main(string[] args)
        {
            commonData = ASETools.CommonData.LoadCommonData(args);
            if (null == commonData)
            {
                return;
            }

            if ((commonData.configuration.commandLineArgs.Count() != 1 || !commonData.diseases.Contains(commonData.configuration.commandLineArgs[0].ToLower())) ||
                (commonData.configuration.commandLineArgs.Count() < 2 || commonData.configuration.commandLineArgs[0] != "-f"))
            {
                Console.WriteLine("usage: CreateGeneExpressionFractionBoxplots disease");
                Console.WriteLine("-or-   CreateGeneExpressionFractionBoxplots -f listOfInputAllcountFiles");
                return;
            }

            var disease = commonData.configuration.commandLineArgs[0].ToLower();
            List<string> allcountFilesToProcess;

            if (commonData.configuration.commandLineArgs.Count() == 1)
            {
                var casesToProcess = commonData.listOfCases.Where(_ => _.disease() == disease).ToList();

                if (casesToProcess.Any(_ => _.gene_expression_fraction_filename == ""))
                {
                    Console.Write("At least one case is missing a gene expression fraction file.");
                    return;
                }

                allcountFilesToProcess = commonData.listOfCases.Where(_ => _.disease() == commonData.configuration.commandLineArgs[1]).Select(_ => _.tumor_rna_allcount_filename).ToList();
            } else
            {
                allcountFilesToProcess = commonData.configuration.commandLineArgs.ToList().GetRange(1, commonData.configuration.commandLineArgs.Count() - 1);
            }



            int nPerDot;
            ASETools.PrintMessageAndNumberBar("Processing", "cases", casesToProcess.Count(), out nPerDot);
            var threading = new ASETools.WorkerThreadHelper<string, Dictionary<string, List<double>>>(allcountFilesToProcess, ProcessOneCase, FinishUp, null, nPerDot);
            threading.run();

            var outputFile = ASETools.CreateStreamWriterWithRetry(@"\temp\GeneExpressionFractionBoxplot_" + disease + ".txt");
            if (null == outputFile)
            {
                Console.WriteLine("Unable to open output file");
                return;
            }

            outputFile.Write("Gene/Case");
            for (int i = 0; i < allcountFilesToProcess.Count(); i++)
            {
                outputFile.Write("\t" + i);
            }
            outputFile.WriteLine();


            foreach (var hugo_symbol in result.Select(_ => _.Key))
            {
                outputFile.Write(hugo_symbol);
                result[hugo_symbol].ForEach(_ => outputFile.Write("\t" + _));
                outputFile.WriteLine();
            } // gene

            outputFile.WriteLine("**done**");
            outputFile.Close();

            outputFile = ASETools.CreateStreamWriterWithRetry(@"\temp\GeneExpressionFractionRanged_" + disease + ".txt");
            if (null == outputFile)
            {
                Console.WriteLine("Unable to open output file");
                return;
            }

            outputFile.WriteLine("gene\tmin\t10%ile\t25%ile\tmedian\t75%ile\t90%ile\tmax\tinter-quartile range\tinter-quartile range multiple");
            foreach (var hugo_symbol in result.Select(_ => _.Key))
            {
                result[hugo_symbol].Sort();
                var n = result[hugo_symbol].Count();

                outputFile.WriteLine(hugo_symbol + "\t" + result[hugo_symbol][0] + "\t" + result[hugo_symbol][n / 10] + "\t" + result[hugo_symbol][n / 4] + "\t" + result[hugo_symbol][n / 2] + "\t" + result[hugo_symbol][3 * n / 4] + "\t" +
                    result[hugo_symbol][9 * n / 10] + "\t" + result[hugo_symbol][n - 1] + "\t" + (result[hugo_symbol][3 * n / 4] - result[hugo_symbol][n / 4]) + "\t" + 
                    (result[hugo_symbol][3 * n / 4] - result[hugo_symbol][n / 4]) / result[hugo_symbol][n / 4]);

            } // gene

            outputFile.WriteLine("**done**");
            outputFile.Close();


            Console.WriteLine("Processed " + casesToProcess.Count() + " cases in " + ASETools.ElapsedTime(commonData.timer));
        } // Main

        static void ProcessOneCase(string allcountFilename, Dictionary<string, List<double>> state)
        {
            var geneExpressionFractions = ASETools.GeneExpressionFraction.LoadFromFile(allcountFilename);

            foreach (var expressionFraction in geneExpressionFractions.Select(_ => _.Value))
            {
                if (!state.ContainsKey(expressionFraction.Hugo_symbol))
                {
                    state.Add(expressionFraction.Hugo_symbol, new List<double>());
                }
                state[expressionFraction.Hugo_symbol].Add(expressionFraction.expressionByTumor[true]);
            } // foreach gene
        } // ProcessOneCase

        static void FinishUp(Dictionary<string, List<double>> state)
        {
            lock (result)
            {
                foreach (var geneEntry in state)
                {
                    if (!result.ContainsKey(geneEntry.Key))
                    {
                        result.Add(geneEntry.Key, geneEntry.Value);
                    } else
                    {
                        result[geneEntry.Key].AddRange(geneEntry.Value);
                    }
                } // each gene
            } // lock
        } // FinishUp

    } // Program
} // namespace
