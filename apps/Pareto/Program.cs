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

        enum ParetoRunType { HAndS, D, DandI};
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

        static void Pareto(string runName, ParetoRunType runType)
        {

            var timingDirectory = @"d:\temp\timings\";
            var concordanceDirectory = @"d:\temp\concordance\";

            int nRunsFound = 0;
            int nRunsWithConcordance = 0;

            var allResults = new List<Result>();

            switch (runType) 
            {
                case ParetoRunType.HAndS:

                    int[] hValues = { 250, 500, 1000, 2000, 4000, 8000, 16000, 32000 };
                    int[] seedSizes = { 16, 20, 22, 25, 27, 30, 32 };

                    foreach (var H in hValues)
                    {
                        foreach (var seedSize in seedSizes)
                        {
                            var timingFilename = timingDirectory + runName + ".snap.H" + H + ".s" + seedSize + "_timings.txt";
                            var concordanceFilename = concordanceDirectory + runName + ".snap.H" + H + ".s" + seedSize + ".concordance.tar";

                            if (File.Exists(timingFilename))
                            {
                                nRunsFound++;

                                var runTiming = ASETools.RunTiming.LoadFromSNAPFile(timingFilename);

                                ASETools.ConcordanceResults concordance = null;

                                if (File.Exists(concordanceFilename) && new FileInfo(concordanceFilename).Length > 10240)
                                {
                                    string sampleNameOverride = null;
                                    if (H == 4000 && seedSize == 25) // This used to be the default, and I renamed the tar files without rerunning them.
                                    {
                                        sampleNameOverride = runName + ".snap.concordance";
                                    }
                                    concordance = new ASETools.ConcordanceResults(concordanceFilename, sampleNameOverride);
                                    nRunsWithConcordance++;
                                }

                                allResults.Add(new Result("" + H + "\t" + seedSize, runTiming, concordance));
                            }
                            else
                            {
                                allResults.Add(new Result("" + H + "\t" + seedSize, 0, null));
                            }

                        } // seed size
                    } // H
                    break;

                case ParetoRunType.D:
                    {
                        int[] dValues = { 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25, 26, 27, 28, 29, 30, 31, 32 };

                        foreach (var d in dValues)
                        {
                            var timingFilename = timingDirectory + runName + ".snap.d" + d + ".s25_timings.txt";
                            var concordanceFilename = concordanceDirectory + runName + ".snap.d" + d + ".s25.concordance.tar";

                            if (File.Exists(timingFilename))
                            {

                                var runTiming = ASETools.RunTiming.LoadFromSNAPFile(timingFilename);

                                if (runTiming != null)
                                {
                                    nRunsFound++;

                                    ASETools.ConcordanceResults condordance = null;

                                    if (File.Exists(concordanceFilename) && new FileInfo(concordanceFilename).Length > 10240)
                                    {
                                        condordance = new ASETools.ConcordanceResults(concordanceFilename);
                                        nRunsWithConcordance++;
                                    }

                                    allResults.Add(new Result("" + d, runTiming, condordance));
                                }
                            }
                        }
                    }

                    break;

                case ParetoRunType.DandI:
                    {
                        int[] dValues = { 15, 17, 19, 21, 23, 25, 27, 29 };
                        int[] iValues = { 27, 30, 35, 40, 45, 50, 55, 60 };

                        foreach (var d in dValues)
                        {
                            foreach (var i in iValues)
                            {
                                var timingFilename = timingDirectory + runName + ".snap.d" + d + ".i" + i + ".s22_timings.txt";
                                var concordanceFilename = concordanceDirectory + runName + ".snap.d" + d + ".i" + i + ".s22.concordance.tar";

                                if (File.Exists(timingFilename))
                                {
                                    var runTiming = ASETools.RunTiming.LoadFromSNAPFile(timingFilename);

                                    if (runTiming != null)
                                    {
                                        nRunsFound++;

                                        ASETools.ConcordanceResults concordance = null;

                                        if (File.Exists(concordanceFilename) && new FileInfo(concordanceFilename).Length > 10240)
                                        {
                                            concordance = new ASETools.ConcordanceResults(concordanceFilename);
                                            nRunsWithConcordance++;
                                        } // if the concordance file exists

                                        allResults.Add(new Result("" + d + "\t" + i, runTiming, concordance));
                                    } // if we loaded a good run timing
                                } // if the timing file exists
                                else { Console.WriteLine("Didn't find d and I file " + timingFilename); }
                            } // for i
                        } // for d
                    }
                    break;

            } // switch (runType)

            Console.WriteLine(runName + ": " + nRunsFound + " total SNAP runs, of which " + nRunsWithConcordance + " have concordance.");

            var outputFilename = @"d:\temp\pareto_" + runName + (runType == ParetoRunType.D ? "-d" : (runType == ParetoRunType.DandI ? "-dAndI" : "")) + ".txt";
            var outputFile = ASETools.CreateStreamWriterWithRetry(outputFilename);

            if (null == outputFile) return;

            switch (runType)
            {
                case ParetoRunType.HAndS:
                    outputFile.Write("H\tseed\t");
                    break;
                case ParetoRunType.D:
                    outputFile.Write("D\t");
                    break;
                case ParetoRunType.DandI:
                    outputFile.Write("D\tI\t");
                    break;

            }
            outputFile.WriteLine("Run time\tSNV F1\tIndel F1\tMean F1\tSNV recall\tSNV precision\tIndel recall\tIndel precision\t2D dominated\t3D dominated\t5D dominated\t2D dominating\t3D dominating\t5D dominating");
            foreach (var result in allResults)
            {
                outputFile.Write(result.resultName + "\t" + result.runTimeInDays);
                if (result.concordance != null)
                {
                    var otherSNAPResults = allResults.Where(_ => _ != result && _.concordance != null).ToList();
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
            Pareto("hg003", ParetoRunType.HAndS);
            Pareto("ERR194147", ParetoRunType.HAndS);
            Pareto("hg007", ParetoRunType.HAndS);
            Pareto("hg003", ParetoRunType.D);
            Pareto("ERR194147", ParetoRunType.D);
            Pareto("hg003", ParetoRunType.DandI);
        } // Main


    } // Program
} // namespace
