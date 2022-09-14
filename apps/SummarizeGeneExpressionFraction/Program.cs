using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using ASELib;

namespace SummarizeGeneExpressionFraction
{
    internal class Program
    {
        static ASETools.CommonData commonData;

        static Dictionary<string, Dictionary<string, ASETools.RunningMeanAndStdDev>> meanExpressionFractionByDisease = new Dictionary<string, Dictionary<string, ASETools.RunningMeanAndStdDev>>();
        static Dictionary<string, Dictionary<string, ASETools.RunningMeanAndStdDev>> meanUnnormalizedExpressionFractionByDisease = new Dictionary<string, Dictionary<string, ASETools.RunningMeanAndStdDev>>();

        static void Main(string[] args)
        {
            commonData = ASETools.CommonData.LoadCommonData(args);

            if (commonData == null)
            {
                return;
            }

            if (commonData.configuration.commandLineArgs.Any(_ => !_.EndsWith(ASETools.geneExpressionFractionExtension)))
            {
                Console.WriteLine("usage: SummarizeGeneExpressionFraction <list of expression files>");
                Console.WriteLine("If you don't specify the list of files, it will run on all cases by disease");
            }


            List<List<string>> expressionFileOrDiseaseSets;
            if (commonData.configuration.commandLineArgs.Count() == 0)
            {
                if (commonData.listOfCases.Any(_ => _.gene_expression_fraction_filename == ""))
                {
                    Console.WriteLine("Not all cases have gene expression fraction data.");
                }

                expressionFileOrDiseaseSets = new List<List<string>>();
                foreach (var disease in commonData.diseases)
                {
if (disease != "laml") continue; // BJB only do AML for now
                    var listOfDisease = new List<string>();
                    listOfDisease.Add(disease);

                    expressionFileOrDiseaseSets.Add(listOfDisease);
                }
            } else
            {
                expressionFileOrDiseaseSets = new List<List<string>>();
                expressionFileOrDiseaseSets.Add(commonData.configuration.commandLineArgs.ToList());
            }

            int nPerDot;
            ASETools.PrintMessageAndNumberBar("Processing", "diseases/input sets", expressionFileOrDiseaseSets.Count(), out nPerDot);
            var threading = new ASETools.WorkerThreadHelper<List<string>, int>(expressionFileOrDiseaseSets, ProcessOneDisease, null, null, nPerDot);
            threading.run();

            string outfilename;
            if (commonData.configuration.commandLineArgs.Count() == 0)
            {
                outfilename = @"\temp\GeneExpressionFractionSummary.txt";
            } else
            {
                outfilename = @"\temp\GeneExpressionFractionSummarySpecial.txt";
            }

            var outputFile = ASETools.CreateStreamWriterWithRetry(outfilename);  // XXX - switch this to the global output directory
            outputFile.Write("Gene");
            var diseasesProcessed = meanExpressionFractionByDisease.Select(_ => _.Key).ToList();
            foreach (var disease in diseasesProcessed)
            {
                outputFile.Write("\t" + disease + " mean\t" + disease + " std dev\t" + disease + " z\t" + disease + " unnormalized\t" + disease + " unnormalized std dev\t" + disease + " unnormalized z"); ;
            }
            outputFile.WriteLine();

            foreach (var hugo_symbol in meanExpressionFractionByDisease[diseasesProcessed[0]].Select(_ => _.Key))
            {
                outputFile.Write(ASETools.ConvertToExcelString(hugo_symbol));

                foreach (var disease in diseasesProcessed)
                {
                    outputFile.Write("\t" + meanExpressionFractionByDisease[disease][hugo_symbol].getMeanAndStdDev().mean + "\t" + meanExpressionFractionByDisease[disease][hugo_symbol].getMeanAndStdDev().stddev + "\t" +
                        meanExpressionFractionByDisease[disease][hugo_symbol].getMeanAndStdDev().stddev / meanExpressionFractionByDisease[disease][hugo_symbol].getMeanAndStdDev().mean + "\t" +
                        meanUnnormalizedExpressionFractionByDisease[disease][hugo_symbol].getMeanAndStdDev().mean + "\t" + meanUnnormalizedExpressionFractionByDisease[disease][hugo_symbol].getMeanAndStdDev().stddev + "\t" +
                        meanUnnormalizedExpressionFractionByDisease[disease][hugo_symbol].getMeanAndStdDev().stddev / meanUnnormalizedExpressionFractionByDisease[disease][hugo_symbol].getMeanAndStdDev().mean);
                }

                outputFile.WriteLine();
            } // foreach gene

            outputFile.WriteLine("**done**");
            outputFile.Close();

            Console.WriteLine("Processed " + commonData.listOfCases.Where(_ => _.gene_expression_fraction_filename != "").Count() + " cases in " + ASETools.ElapsedTime(commonData.timer));

        } // Main

        static void ProcessOneDisease(List<string> expressionFractionFilenamesOrDisease, int state)
        {
            List<Dictionary<string, ASETools.GeneExpressionFraction>> expressionFractions;
            string diseaseName;
            bool tumor; // We use tumor for TCGA and normal for special
            if (expressionFractionFilenamesOrDisease.Count() == 1)
            {
                tumor = true;
                diseaseName = expressionFractionFilenamesOrDisease[0];
                expressionFractions = commonData.listOfCases.Where(_ => _.disease() == diseaseName).Select(_ => ASETools.GeneExpressionFraction.LoadFromFile(_.gene_expression_fraction_filename)).ToList();
            }
            else
            {
                tumor = false;
                diseaseName = "";
                expressionFractions = expressionFractionFilenamesOrDisease.Select(_ => ASETools.GeneExpressionFraction.LoadFromFile(_)).ToList();
            }

            if (expressionFractions.Count() == 0)
            {
                Console.WriteLine("Didn't find any cases");
                return;
            }

            var meanExpressionFractions = new Dictionary<string, ASETools.RunningMeanAndStdDev>();
            var meanUnnormalizedExpressionFractions = new Dictionary<string, ASETools.RunningMeanAndStdDev>();

            foreach (var hugo_symbol in expressionFractions[0].Keys)
            {         
                var runningMeanAndStdDev = new ASETools.RunningMeanAndStdDev();

                expressionFractions.Where(_ => _.ContainsKey(hugo_symbol)).ToList().ForEach(_ => runningMeanAndStdDev.addValue(_[hugo_symbol].normalizedExpressionByTumor[tumor]));
                meanExpressionFractions.Add(hugo_symbol, runningMeanAndStdDev);

                runningMeanAndStdDev = new ASETools.RunningMeanAndStdDev();
                expressionFractions.Where(_ => _.ContainsKey(hugo_symbol)).ToList().ForEach(_ => runningMeanAndStdDev.addValue(_[hugo_symbol].expressionByTumor[tumor]));
                meanUnnormalizedExpressionFractions.Add(hugo_symbol, runningMeanAndStdDev);
            }

            lock (meanExpressionFractionByDisease)
            {
                meanExpressionFractionByDisease.Add(diseaseName, meanExpressionFractions);
                meanUnnormalizedExpressionFractionByDisease.Add(diseaseName, meanUnnormalizedExpressionFractions);
            } // lock

        } // ProcessOneDisease

    } // Program
} // namespace
