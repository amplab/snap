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
            public Result(string resultName_, ASETools.RunTiming timing_, ASETools.ConcordanceResults concordance_)
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

            public bool Dominates2D(Result peer)
            {
                return timing.alignTime <= peer.timing.alignTime && meanF1() >= peer.meanF1();
            }

            public bool Dominates3D(Result peer)
            {
                return timing.alignTime <= peer.timing.alignTime && F1(ASETools.VariantType.Indel) >= peer.F1(ASETools.VariantType.Indel) && F1(ASETools.VariantType.SNV) >= peer.F1(ASETools.VariantType.SNV);
            }

            public bool Dominates5D(Result peer)
            {
                return timing.alignTime <= peer.timing.alignTime && recall(ASETools.VariantType.Indel) >= peer.recall(ASETools.VariantType.Indel) && precision(ASETools.VariantType.Indel) >= peer.precision(ASETools.VariantType.Indel) &&
                       recall(ASETools.VariantType.SNV) >= peer.recall(ASETools.VariantType.SNV) && precision(ASETools.VariantType.SNV) >= peer.precision(ASETools.VariantType.SNV);
            }

            public double meanF1()
            {
                return (F1(ASETools.VariantType.Indel) + F1(ASETools.VariantType.SNV)) / 2;
            }

            public double F1(ASETools.VariantType variantType)
            {
                return concordance.results[variantType].F1_score;
            }

            public double recall(ASETools.VariantType variantType)
            {
                return concordance.results[variantType].recall;
            }

            public double precision(ASETools.VariantType variantType)
            {
                return concordance.results[variantType].precision;
            }

            readonly public ASETools.RunTiming timing;
            readonly public ASETools.ConcordanceResults concordance;
            readonly public double runTimeInDays;
            readonly public string resultName;
        }

        static void Pareto(string runName)
        {
            int[] hValues = { 250, 500, 1000, 2000, 4000, 8000, 16000, 32000 };
            int[] seedSizes = { 16, 20, 22, 25, 27, 30, 32 };

            var timingDirectory = @"d:\temp\timings\";
            var concordanceDirectory = @"d:\temp\concordance\";

            var results = new Dictionary<int, Dictionary<int, Result>>();

            int nRunsFound = 0;
            int nRunsWithConcordance = 0;

            foreach (var H in hValues)
            {
                results.Add(H, new Dictionary<int, Result>());
                foreach (var seedSize in seedSizes)
                {
                    var timingFilename = timingDirectory + runName + ".snap.H" + H + ".s" + seedSize + "_timings.txt";
                    var concordanceFilename = concordanceDirectory + runName + ".snap.H" + H + ".s" + seedSize + ".concordance.tar";

                    //
                    // Special case H4000 s25, which is the default.
                    //
                    if (H == 4000 && seedSize == 25)
                    {
                        timingFilename = timingDirectory + runName + ".snap_timings.txt";
                        concordanceFilename = concordanceDirectory + runName + ".snap.concordance.tar";
                    }

                    if (File.Exists(timingFilename))
                    {
                        nRunsFound++;

                        var runTiming = ASETools.RunTiming.LoadFromSNAPFile(timingFilename);

                        ASETools.ConcordanceResults concordance = null;

                        if (File.Exists(concordanceFilename))
                        {
                            concordance = new ASETools.ConcordanceResults(concordanceFilename);
                            nRunsWithConcordance++;
                        }

                        results[H].Add(seedSize, new Result("" + H + "\t" + seedSize, runTiming, concordance));
                    }
                    else
                    {
                        results[H].Add(seedSize, new Result("" + H + "\t" + seedSize, 0, null));
                    }

                } // seed size
            } // H

            Console.WriteLine(runName + ": " + nRunsFound + " total SNAP runs, of which " + nRunsWithConcordance + " have concordance.");

            var allResults = new List<Result>();

            foreach (var aligner in ASETools.EnumUtil.GetValues<ASETools.Aligner>())
            {
                if (aligner == ASETools.Aligner.SNAP) continue;

                var timingFilename = timingDirectory + runName + "." + ASETools.alignerName[aligner].ToLower() + "_timings.txt";
                if (File.Exists(timingFilename))
                {
                    ASETools.ConcordanceResults concordance = null;

                    var concordanceFilename = concordanceDirectory + runName + "." + ASETools.alignerName[aligner].ToLower() + ".concordance.tar";
                    if (File.Exists(concordanceFilename))
                    {
                        concordance = new ASETools.ConcordanceResults(concordanceFilename);
                    }

                    allResults.Add(new Result(ASETools.alignerName[aligner] + "\t", ASETools.RunTiming.LoadFromLinuxFile(timingFilename), concordance));
                }
            }

#if false

            var bwa = new Result("BWA\t", 0.832175926, new ASETools.ConcordanceResults(new ASETools.SingleConcordanceResult("SNP", 0.999568, 0.982369), new ASETools.SingleConcordanceResult("INDEL", 0.961477, 0.955432)));
            var bowtie = new Result("Bowtie\t", 1.045185185, new ASETools.ConcordanceResults(new ASETools.SingleConcordanceResult("SNP", 0.978907, 0.993347), new ASETools.SingleConcordanceResult("INDEL", 0.95743, 0.965235)));
            var novoalign = new Result("Novoalign\t", 8.247222222, new ASETools.ConcordanceResults(new ASETools.SingleConcordanceResult("SNP", 0.999556, 0.988173), new ASETools.SingleConcordanceResult("INDEL", 0.968255, 0.963358)));
            var gem = new Result("GEM\t", 0.222291667, new ASETools.ConcordanceResults(new ASETools.SingleConcordanceResult("SNP", 0.983147, 0.992375), new ASETools.SingleConcordanceResult("INDEL", 0.950528, 0.963512)));

            allResults.Add(bwa);
            allResults.Add(bowtie);
            allResults.Add(novoalign);
            allResults.Add(gem);
#endif // false

            var snapResultsWithConcordance = new List<Result>();

            foreach (var H in hValues)
            {
                foreach (var seedSize in seedSizes)
                {
                    allResults.Add(results[H][seedSize]);
                    if (results[H][seedSize].concordance != null)
                    {
                        snapResultsWithConcordance.Add(results[H][seedSize]);
                    }
                } // seedSize
            } // H

            var outputFilename = @"d:\temp\pareto_" + runName + ".txt";
            var outputFile = ASETools.CreateStreamWriterWithRetry(outputFilename);

            if (null == outputFile) return;

            outputFile.WriteLine("H\tseed\tRun time\tSNV F1\tIndel F1\tMean F1\tSNV recall\tSNV precision\tIndel recall\tIndel precision\t2D dominated\t3D dominated\t5D dominated\t2D dominating\t3D dominating\t5D dominating");
            foreach (var result in allResults)
            {
                outputFile.Write(result.resultName + "\t" + result.runTimeInDays);
                if (result.concordance != null)
                {
                    var otherSNAPResults = snapResultsWithConcordance.Where(_ => _ != result).ToList();
                    outputFile.Write("\t" + result.F1(ASETools.VariantType.SNV) + "\t" + result.F1(ASETools.VariantType.Indel) + "\t" +
                        result.meanF1() + "\t" +
                        result.recall(ASETools.VariantType.SNV) + "\t" + result.precision(ASETools.VariantType.SNV) + "\t" +
                        result.recall(ASETools.VariantType.Indel) + "\t" + result.precision(ASETools.VariantType.Indel) + "\t" +
                        otherSNAPResults.Any(_ => _.Dominates2D(result)) + "\t" + otherSNAPResults.Any(_ => _.Dominates3D(result)) + "\t" + otherSNAPResults.Any(_ => _.Dominates5D(result)) + "\t" +
                        otherSNAPResults.All(_ => result.Dominates2D(_)) + "\t" + otherSNAPResults.All(_ => result.Dominates3D(_)) + "\t" + otherSNAPResults.All(_ => result.Dominates5D(_)));
                }
                outputFile.WriteLine();
            } // all results

            outputFile.Close();
        }
        static void Main(string[] args)
        {
            Pareto("hg003");
            Pareto("ERR194147");
        } // Main


        } // Program
} // namespace
