using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using ASELib;
using System.Diagnostics;
using System.Threading;

namespace GenerateHistogramsForSignificantResults
{
    class Program
    {
        class SignificantResultAndHistograms
        {
            public SignificantResultAndHistograms(ASETools.AllSignificantResultsLine result_)
            {
                result = result_;
 
            }

            public readonly ASETools.AllSignificantResultsLine result;
            public ASETools.Histogram oneHistogram = new ASETools.Histogram(); // Histogram for tumors with exactly one mutation
            public ASETools.Histogram otherHistogram = new ASETools.Histogram();  // The histogram for the other category, which may be not-1 or may be > 1 depending on the result
        }

        static void WorkerThread(List<ASETools.Case> casesToProcess, List<SignificantResultAndHistograms> resultsAndHistograms, List<string> genesToProcess, ref int nProcessed, ref int nMissing)
        {
            while (true)
            {
                ASETools.Case case_;
                lock (casesToProcess)
                {
                    if (casesToProcess.Count() == 0)
                    {
                        return;
                    }

                    case_ = casesToProcess[0];
                    casesToProcess.RemoveAt(0);

                    nProcessed++;
                    if (nProcessed % 100 == 0)
                    {
                        Console.Write(".");
                    }

                    if (case_.tumor_allele_specific_gene_expression_filename == "")
                    {
                        nMissing++;
                        continue;
                    }
                } // lock

                var expressionForThisCase = ASETools.ASESignalLine.ReadFile(case_.tumor_allele_specific_gene_expression_filename, genesToProcess);
                var aseLines = ASETools.ASESignalLine.ReadFile(case_.tumor_allele_specific_gene_expression_filename);

                foreach (var resultAndHistograms in resultsAndHistograms)
                {
                    var hugo_symbol = resultAndHistograms.result.hugo_symbol;

                    if (resultAndHistograms.result.inputFile != ASETools.AlleleSpecificExpressionDistributionByMutationCountFilenameBase + ".txt")
                    {
                        if (!resultAndHistograms.result.inputFile.Contains("_") || !resultAndHistograms.result.inputFile.EndsWith(".txt") ||
                            !ASETools.GetFileNameFromPathname(resultAndHistograms.result.inputFile).StartsWith(ASETools.AlleleSpecificExpressionDistributionByMutationCountFilenameBase))
                        {
                            Console.WriteLine("Can't parse input filename " + resultAndHistograms.result.inputFile);
                            continue;
                        }

                        var diseaseDotText = resultAndHistograms.result.inputFile.Substring(resultAndHistograms.result.inputFile.LastIndexOf("_") + 1).ToLower();
                        if (case_.disease() + ".txt" != diseaseDotText)
                        {
                            continue;
                        }
                    } // if this is a per-disease result

                    if (!aseLines.ContainsKey(hugo_symbol))
                    {
                        continue;   // This case doens't have data for this gene.
                    }

                    if (aseLines[hugo_symbol].nMutations == 0 && resultAndHistograms.result.one_vs_many)
                    {
                        continue;   // This result is for one mutation vs. many, but the tumor has zero.  Ignore it.
                    }

                    double[] ase = resultAndHistograms.result.exclusive ? aseLines[hugo_symbol].exclusiveASE : aseLines[hugo_symbol].inclusiveASE;

                    if (ase[resultAndHistograms.result.range_index] == double.NegativeInfinity)
                    {
                        continue;   // There wasn't a result for this range because of no germline variant site in the range.
                    }

                    lock (resultsAndHistograms)
                    {
                        if (aseLines[hugo_symbol].nMutations == 1)
                        {
                            resultAndHistograms.oneHistogram.addValue(ase[resultAndHistograms.result.range_index]);
                        } else
                        {
                            resultAndHistograms.otherHistogram.addValue(ase[resultAndHistograms.result.range_index]);
                        }
                    }
                }
            } // while true
        }


        static void Main(string[] args)
        {
            var congfiguration = ASETools.Configuration.loadFromFile(args);

            if (null == congfiguration)
            {
                Console.WriteLine("Unable to load configuration");
                return;
            }

            if (congfiguration.commandLineArgs.Count() != 0)
            {
                Console.WriteLine("usage: GenerateHistogramsForSignificantResults");
                return;
            }

            var allSignificantResults = ASETools.AllSignificantResultsLine.ReadFile(congfiguration.finalResultsDirectory + ASETools.AllSignificantResultsFilename);

            if (null == allSignificantResults)
            {
                Console.WriteLine("Unable to load AllSignificantResults.txt file from " + congfiguration.finalResultsDirectory + ASETools.AllSignificantResultsFilename);
                return;
            }

            var resultsAndHistograms = allSignificantResults.Select(x => new SignificantResultAndHistograms(x)).ToList();
            resultsAndHistograms.Sort((x, y) => String.Compare(x.result.hugo_symbol, y.result.hugo_symbol));

            var cases = ASETools.Case.LoadCases(congfiguration.casesFilePathname);
            if (null == cases)
            {
                Console.WriteLine("Unable to load cases file.");
                return;
            }

            var genesToProcess = new List<string>();
            foreach (var result in allSignificantResults)
            {
                if (!genesToProcess.Contains(result.hugo_symbol))
                {
                    genesToProcess.Add(result.hugo_symbol);
                }
            }

            int nMissing = 0;
            int nProcessed = 0;

            var timer = new Stopwatch();
            timer.Start();

            Console.Write("Processing " + cases.Count() + " cases, 1 dot/100 cases: ");

            var casesToProcess = cases.Select(x => x.Value).ToList();

            var threads = new List<Thread>();
            for (int i = 0; i < Environment.ProcessorCount; i++)
            {
                threads.Add(new Thread(() => WorkerThread(casesToProcess, resultsAndHistograms, genesToProcess, ref nProcessed, ref nMissing)));
            }

            threads.ForEach(t => t.Start());
            threads.ForEach(t => t.Join());

            Console.WriteLine();
            Console.WriteLine("Took " + ASETools.ElapsedTimeInSeconds(timer));

            var outputFile = ASETools.CreateStreamWriterWithRetry(congfiguration.finalResultsDirectory + ASETools.HistogramsForSignficantResultsFilename);

            string[] mutationCountNames = { "0", "1", ">1" };

            foreach (var resultAndHistograms in resultsAndHistograms)
            {
                outputFile.WriteLine(resultAndHistograms.result.hugo_symbol + " " + resultAndHistograms.result.inputFile + " exclusive: " + resultAndHistograms.result.exclusive +
                    " one vs. many: " + resultAndHistograms.result.one_vs_many + " 1 mutation, " + resultAndHistograms.result.range);
                outputFile.WriteLine(ASETools.HistogramResultLine.Header());
                resultAndHistograms.oneHistogram.ComputeHistogram(0, 1, .01).ToList().ForEach(x => outputFile.WriteLine(x));
                outputFile.WriteLine();

                outputFile.WriteLine(resultAndHistograms.result.hugo_symbol + " " + resultAndHistograms.result.inputFile + " exclusive: " + resultAndHistograms.result.exclusive +
                    " one vs. many: " + resultAndHistograms.result.one_vs_many + (resultAndHistograms.result.one_vs_many ? " > 1 mutations, " : " not one mutation, ") + resultAndHistograms.result.range);
                outputFile.WriteLine(ASETools.HistogramResultLine.Header());
                resultAndHistograms.otherHistogram.ComputeHistogram(0, 1, .01).ToList().ForEach(x => outputFile.WriteLine(x));
                outputFile.WriteLine();


            }

            outputFile.WriteLine("**done**");
            outputFile.Close();

        }
    }
}
