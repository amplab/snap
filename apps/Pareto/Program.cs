using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using ASELib;
using System.IO;
using System.Globalization;

namespace Pareto
{
    class Program
    {

        class Result
        {
            public Result(string resultName_, ASETools.SNAPRunTiming timing_, ASETools.ConcordanceResults concordance_)
            {
                resultName = resultName_;
                timing = timing_;
                concordance = concordance_;

                runTimeInDays = ((double)timing.alignTime + timing.loadingTime) / 24 / 3600;    // Dividing converts from seconds to days
            }

            public Result(string resultName_, double runTimeInDays_, ASETools.ConcordanceResults concordance_)
            {
                resultName = resultName_;
                runTimeInDays = runTimeInDays_;
                concordance = concordance_;
            }

            readonly public ASETools.SNAPRunTiming timing;
            readonly public ASETools.ConcordanceResults concordance;
            readonly public double runTimeInDays;
            readonly public string resultName;
        }
        static void Main(string[] args)
        {
            int[] hValues = { 250, 500, 1000, 2000, 4000, 8000, 16000, 32000 };
            int[] seedSizes = { 16, 20, 22, 25, 27, 30, 32 };

            var timingDirectory = @"d:\temp\timings\";
            var concordanceDiretory = @"d:\temp\concordance\";
            var runName = "hg003";

            var results = new Dictionary<int, Dictionary<int, Result>>();

            int nRunsFound = 0;
            int nRunsWithConcordance = 0;

            foreach (var H in hValues)
            {
                results.Add(H, new Dictionary<int, Result>());
                foreach (var seedSize in seedSizes)
                {
                    var timingFilename = timingDirectory + runName + ".snap.H" + H + ".s" + seedSize + "_timings.txt";
                    var concordanceFilename = concordanceDiretory + runName + ".snap.H" + H + ".s" + seedSize + ".concordance.tar";

                    if (File.Exists(timingFilename))
                    {
                        nRunsFound++;

                        var runTiming = ASETools.SNAPRunTiming.LoadFromFile(timingFilename);

                        ASETools.ConcordanceResults concordance = null;

                        if (File.Exists(concordanceFilename))
                        {
                            concordance = new ASETools.ConcordanceResults(concordanceFilename);
                            nRunsWithConcordance++;
                        }

                        results[H].Add(seedSize, new Result("" + H + "\t" + seedSize,runTiming, concordance));
                    } else
                    {
                        results[H].Add(seedSize, new Result("" + H + "\t" + seedSize, 0, null));
                    }

                } // seed size
            } // H

            Console.WriteLine(nRunsFound + " total runs, of which " + nRunsWithConcordance + " have concordance.");


            var bwa = new Result("BWA\t", 0.832175926, new ASETools.ConcordanceResults(new ASETools.SingleConcordanceResult("SNP", 0.999568, 0.982369), new ASETools.SingleConcordanceResult("INDEL", 0.961477, 0.955432)));
            var bowtie = new Result("Bowtie\t", 1.045185185, new ASETools.ConcordanceResults(new ASETools.SingleConcordanceResult("SNP", 0.978907, 0.993347), new ASETools.SingleConcordanceResult("INDEL", 0.95743, 0.965235)));
            var novoalign = new Result("Novoalign\t", 8.247222222, new ASETools.ConcordanceResults(new ASETools.SingleConcordanceResult("SNP", 0.999556, 0.988173), new ASETools.SingleConcordanceResult("INDEL", 0.968255, 0.963358)));
            var gem = new Result("GEM\t", 0.222291667, new ASETools.ConcordanceResults(new ASETools.SingleConcordanceResult("SNP", 0.983147, 0.992375), new ASETools.SingleConcordanceResult("INDEL", 0.950528, 0.963512)));

            var allResults = new List<Result>();
            allResults.Add(bwa);
            allResults.Add(bowtie);
            allResults.Add(novoalign);
            allResults.Add(gem);


            foreach (var H in hValues)
            {
                foreach (var seedSize in seedSizes)
                {
                    allResults.Add(results[H][seedSize]);
                } // seedSize
            } // H

            var outputFilename = @"d:\temp\pareto_output.txt";
            var outputFile = ASETools.CreateStreamWriterWithRetry(outputFilename);

            if (null == outputFile) return;

            outputFile.WriteLine("H\tseed\tRun time\tSNV F1\tIndel F1\tMean F1\tSNV recall\tSNV precision\tIndel recall\tIndel precision");
            foreach (var result in allResults)
            {
                outputFile.Write(result.resultName + "\t" + result.runTimeInDays);
                if (result.concordance != null)
                {
                    outputFile.Write("\t" + result.concordance.results[ASETools.VariantType.SNV].F1_score + "\t" + result.concordance.results[ASETools.VariantType.Indel].F1_score + "\t" +
                        (result.concordance.results[ASETools.VariantType.SNV].F1_score + result.concordance.results[ASETools.VariantType.Indel].F1_score)/2 + "\t" +
                        result.concordance.results[ASETools.VariantType.SNV].recall + "\t" + result.concordance.results[ASETools.VariantType.SNV].precision + "\t" +
                        result.concordance.results[ASETools.VariantType.Indel].recall + "\t" + result.concordance.results[ASETools.VariantType.Indel].precision);
                }
                outputFile.WriteLine();
            } // all results

            outputFile.Close();
        } // Main


        } // Program
} // namespace
