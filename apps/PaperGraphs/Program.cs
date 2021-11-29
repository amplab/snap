using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using ASELib;
using System.IO;
using System.Runtime.CompilerServices;

namespace PaperGraphs
{
    class Program
    {

        class Result
        {
            public readonly ASETools.RunTiming runTiming;
            public readonly ASETools.ConcordanceResults concordance;
            public readonly List<ASETools.RunTiming> replicas;

            public Result(ASETools.RunTiming runTiming_, ASETools.ConcordanceResults concordance_, List<ASETools.RunTiming> replicas_)
            {
                runTiming = runTiming_;
                concordance = concordance_;
                replicas = replicas_;
            }
        }

        delegate double ExtractConcordanceValue(ASETools.ConcordanceResults concordance);

        class ScatterGraphInfo
        {
            public readonly string name;
            public readonly ExtractConcordanceValue extractor;

            public ScatterGraphInfo(string name_, ExtractConcordanceValue extractor_)
            {
                name = name_;
                extractor = extractor_;
            } // ctor

        } // ScatterGraphInfo

        static string timingsDirectory = @"d:\temp\timings\";
        static string concordanceDirectory = @"d:\temp\concordance-dragen3.8VC\";

        static void RunSeedSizeSweep(StreamWriter outputFile, string sample)
        {
            var seedSizesToRun = new List<int>();

            for (int seedSize = 13; seedSize <= 32; seedSize++)
            {
                if (seedSize == 24)
                {
                    if (File.Exists(timingsDirectory + sample + ".snap_timings.txt") && File.Exists(concordanceDirectory + sample + ".snap.dragen3.8VC.concordance.tar"))
                    {
                        seedSizesToRun.Add(seedSize);
                    }
                } 
                else if (File.Exists(timingsDirectory + sample + ".snap_s" + seedSize + "_timings.txt") && File.Exists(concordanceDirectory + sample + ".snap_s" + seedSize + ".dragen3.8VC.concordance.tar"))
                {
                    seedSizesToRun.Add(seedSize);
                }
            } // for each possible seed size

            if (seedSizesToRun.Count() == 0)
            {
                return;
            }

            outputFile.WriteLine();
            outputFile.WriteLine(sample + " seed size sweep");
            outputFile.WriteLine("seed size\ttime\tmean F1\tSNV F1\tindel F1");

            foreach (var seedSize in seedSizesToRun)
            {
                string timingFilename;
                string concordanceFilename;

                if (seedSize == 24)
                {
                    timingFilename = timingsDirectory + sample + ".snap_timings.txt";
                    concordanceFilename = concordanceDirectory + sample + ".snap.dragen3.8VC.concordance.tar";
                } else
                {
                    timingFilename = timingsDirectory + sample + ".snap_s" + seedSize + "_timings.txt";
                    concordanceFilename = concordanceDirectory + sample + ".snap_s" + seedSize + ".dragen3.8VC.concordance.tar";
                }

                var timing = ASETools.RunTiming.LoadFromSNAPFile(timingFilename);
                var concordance = new ASETools.ConcordanceResults(concordanceFilename);

                outputFile.WriteLine(seedSize + "\t" + ((double)timing.alignTime + timing.loadingTime) / 3600 / 24 + "\t" + concordance.Mean_F1_Score() + "\t" + concordance.results[ASETools.VariantType.SNV].F1_score + "\t" + concordance.results[ASETools.VariantType.Indel].F1_score);
            }

        } // RunSeedSizeSweep


        static int nDigits = 4;

        static void Main(string[] args)
        {


            string[] dataSets = { "ERR194146", "ERR194147" , "hg001", "hg002", "hg003", "hg004", "hg005", "hg006", "hg007",  };

            Dictionary<string, Dictionary<ASETools.Aligner, Result>> resultsByDataSetAndAligner = new Dictionary<string, Dictionary<ASETools.Aligner, Result>>();
            Dictionary<string, Dictionary<ASETools.Aligner, List<Result>>> replicasByDataSetAndAligner = new Dictionary<string, Dictionary<ASETools.Aligner, List<Result>>>();

            var aligners = new List<ASETools.Aligner>();    // The same as the aligners from ASETools, but in a good order for the paper
            aligners.Add(ASETools.Aligner.SNAP);
            aligners.Add(ASETools.Aligner.Dragen);
            aligners.Add(ASETools.Aligner.Novoalign);
            aligners.Add(ASETools.Aligner.BWA);
            aligners.Add(ASETools.Aligner.GEM);
            aligners.Add(ASETools.Aligner.Bowtie);


            int nTotal = 0;
            int nWithTimings = 0;
            int nWithConcordance = 0;

            foreach (var dataSet in dataSets) {
                resultsByDataSetAndAligner.Add(dataSet, new Dictionary<ASETools.Aligner, Result>());

                foreach (var aligner in aligners)
                {
                    nTotal++;

                    var timingFilename = timingsDirectory + dataSet + "." + ASETools.alignerName[aligner] + (aligner == ASETools.Aligner.Bowtie ? "_s1000" : (aligner == ASETools.Aligner.Dragen ? "3.8" : "")) + "_timings.txt";
                    var splitTimingsFilenames = new List<string>();
                    for (int i = 0; i < 9; i++)
                    {
                        splitTimingsFilenames.Add(timingsDirectory + dataSet + "." + i + "." + ASETools.alignerName[aligner]  + "_timings.txt");
                    }

                    var secondStageFilename = timingsDirectory + dataSet + "." + ASETools.alignerName[aligner] + "_second_stage_timings.txt";

                    ASETools.RunTiming timings = null;
                    ASETools.ConcordanceResults concordance = null;
                    var replicas = new List<ASETools.RunTiming>();

                    if (File.Exists(timingFilename) || splitTimingsFilenames.All(_ => File.Exists(_) && File.Exists(secondStageFilename))) {
                        nWithTimings++;


                        if (aligner == ASETools.Aligner.SNAP)
                        {
                            timings = ASETools.RunTiming.LoadFromSNAPFile(timingFilename);
                        } else
                        {
                            if (File.Exists(timingFilename))
                            {
                                timings = ASETools.RunTiming.LoadFromLinuxFile(timingFilename);
                            } else
                            {
                                timings = ASETools.RunTiming.LoadFromSplitFiles(splitTimingsFilenames, secondStageFilename);
                            }
                        }

                        replicas.Add(timings);

                        foreach (var replicaFilename in Directory.GetFiles(timingsDirectory, dataSet + "." + ASETools.alignerName[aligner] +  (aligner == ASETools.Aligner.Bowtie ? "_s1000" : "") + "_timings_*.txt"))
                        {
                            if (replicaFilename.Contains("_timings_s"))
                            {
                                continue;   //SNAP runs with different seed sizes.
                            }

                            var replica = ASETools.RunTiming.LoadFromReplicaFile(replicaFilename);
                            if (replica != null)
                            {
                                replicas.Add(replica);
                            }
                        }

                    }

                    var concordanceFilename = concordanceDirectory + dataSet + "." + ASETools.alignerName[aligner].ToLower() + (aligner == ASETools.Aligner.Bowtie ? "_s1000" : (aligner == ASETools.Aligner.Dragen ? "3.8" : "")) + ".dragen3.8VC.concordance.tar";
                    //Console.WriteLine(concordanceFilename);

                    if (File.Exists(concordanceFilename))
                    {
                        nWithConcordance++;
                        concordance = new ASETools.ConcordanceResults(concordanceFilename);
                    }

                    resultsByDataSetAndAligner[dataSet].Add(aligner, new Result(timings, concordance, replicas));

                } // aligner
            } // dataSet

            Console.WriteLine("Of " + nTotal + " possible runs, " + nWithTimings + " have timings and " + nWithConcordance + " have concordance");             

            var outputFile = ASETools.CreateStreamWriterWithRetry(@"d:\temp\Aligner Concordance and Timings.txt");
            if (null == outputFile)
            {
                return;
            }

            //
            // Write the first header line
            //
            string oneTabPerAligner = "";
            foreach (var aligner in aligners) oneTabPerAligner += "\t";

            outputFile.WriteLine("\tMean F1" + oneTabPerAligner + "SNV F1" + oneTabPerAligner + "Indel F1" + oneTabPerAligner + "SNV Precision" + oneTabPerAligner + "SNV Recall" + oneTabPerAligner + "Indel Precision" + oneTabPerAligner + "Indel Recall");

            for (int i = 0; i < 7; i++)
            {
                foreach (var aligner in aligners)
                {
                    outputFile.Write("\t" + ASETools.alignerName[aligner]);
                }
            }

            outputFile.WriteLine();

            foreach (var dataSet in dataSets)
            {
                outputFile.Write(dataSet);
                foreach (var aligner in aligners)
                {
                    outputFile.Write("\t");
                    var concordance = resultsByDataSetAndAligner[dataSet][aligner].concordance;
                    if (concordance != null)
                    {
                        outputFile.Write(ASETools.DoubleRounded(resultsByDataSetAndAligner[dataSet][aligner].concordance.Mean_F1_Score(), nDigits));
                    }
                } // aligner

                foreach (var aligner in aligners)
                {
                    outputFile.Write("\t");
                    var concordance = resultsByDataSetAndAligner[dataSet][aligner].concordance;
                    if (concordance != null)
                    {
                        outputFile.Write(ASETools.DoubleRounded(resultsByDataSetAndAligner[dataSet][aligner].concordance.results[ASETools.VariantType.SNV].F1_score, nDigits));
                    }
                } // aligner

                foreach (var aligner in aligners)
                {
                    outputFile.Write("\t");
                    var concordance = resultsByDataSetAndAligner[dataSet][aligner].concordance;
                    if (concordance != null)
                    {
                        outputFile.Write(ASETools.DoubleRounded(resultsByDataSetAndAligner[dataSet][aligner].concordance.results[ASETools.VariantType.Indel].F1_score, nDigits));
                    }
                } // aligner
                
                foreach (var aligner in aligners)
                {
                    outputFile.Write("\t");
                    var concordance = resultsByDataSetAndAligner[dataSet][aligner].concordance;
                    if (concordance != null)
                    {
                        outputFile.Write(ASETools.DoubleRounded(resultsByDataSetAndAligner[dataSet][aligner].concordance.results[ASETools.VariantType.SNV].precision, nDigits));
                    }
                } // aligner

                foreach (var aligner in aligners)
                {
                    outputFile.Write("\t");
                    var concordance = resultsByDataSetAndAligner[dataSet][aligner].concordance;
                    if (concordance != null)
                    {
                        outputFile.Write(ASETools.DoubleRounded(resultsByDataSetAndAligner[dataSet][aligner].concordance.results[ASETools.VariantType.SNV].recall, nDigits));
                    }
                } // aligner

                foreach (var aligner in aligners)
                {
                    outputFile.Write("\t");
                    var concordance = resultsByDataSetAndAligner[dataSet][aligner].concordance;
                    if (concordance != null)
                    {
                        outputFile.Write(ASETools.DoubleRounded(resultsByDataSetAndAligner[dataSet][aligner].concordance.results[ASETools.VariantType.Indel].precision, nDigits));
                    }
                } // aligner

                foreach (var aligner in aligners)
                {
                    outputFile.Write("\t");
                    var concordance = resultsByDataSetAndAligner[dataSet][aligner].concordance;
                    if (concordance != null)
                    {
                        outputFile.Write(ASETools.DoubleRounded(resultsByDataSetAndAligner[dataSet][aligner].concordance.results[ASETools.VariantType.Indel].recall, nDigits));
                    }
                } // aligner

                outputFile.WriteLine();
            } // dataSet

            outputFile.WriteLine();

            outputFile.WriteLine("\tSNV" + oneTabPerAligner + "Indel" + oneTabPerAligner + "Mean");
            outputFile.Write("9s");

            // 9s (= -1Log1(1-F1 score))
            for (int i = 0; i < 3; i++)
            {
                foreach (var aligner in aligners)
                {
                    var alignerName = ASETools.alignerName[aligner];
                    outputFile.Write("\t" + alignerName);
                } // aligner
            }
            outputFile.WriteLine();

            foreach (var dataSet in dataSets)
            {
                outputFile.Write(dataSet);
                foreach (var aligner in aligners)
                {
                    var concordance = resultsByDataSetAndAligner[dataSet][aligner].concordance;
                    if (concordance != null)
                    {
                        outputFile.Write("\t" + ASETools.DoubleRounded(-Math.Log10(1 - concordance.results[ASETools.VariantType.SNV].F1_score), nDigits));
                    }
                    else
                    {
                        outputFile.Write("\t");
                    }// If we have concordance
                } // aligner

                foreach (var aligner in aligners)
                {
                    var concordance = resultsByDataSetAndAligner[dataSet][aligner].concordance;
                    if (concordance != null)
                    {
                        outputFile.Write("\t" + ASETools.DoubleRounded(-Math.Log10(1 - concordance.results[ASETools.VariantType.Indel].F1_score), nDigits));
                    }
                    else
                    {
                        outputFile.Write("\t");
                    }// If we have concordance
                } // aligner

                foreach (var aligner in aligners)
                {
                    var concordance = resultsByDataSetAndAligner[dataSet][aligner].concordance;
                    if (concordance != null)
                    {
                        outputFile.Write("\t" + ASETools.DoubleRounded(-Math.Log10(1 - concordance.Mean_F1_Score()), nDigits));
                    }
                    else
                    {
                        outputFile.Write("\t");
                    }// If we have concordance
                } // aligner
                outputFile.WriteLine();
            } // dataSet

            outputFile.WriteLine();

            outputFile.WriteLine("\tLoad & Align" + oneTabPerAligner + "Load & Align error bar" + oneTabPerAligner + "Whole Pipeline");
            outputFile.Write("Run Time");
            for (int i = 0; i < 3; i++)
            {
                foreach (var aligner in aligners)
                {
                    var alignerName = ASETools.alignerName[aligner];
                    outputFile.Write("\t" + alignerName);
                } // aligner
            }
            outputFile.WriteLine();

            foreach (var dataSet in dataSets)
            {
                outputFile.Write(dataSet);

                // Load & Align
                foreach (var aligner in aligners)
                {
                    var timing = resultsByDataSetAndAligner[dataSet][aligner].runTiming;
                    if (timing != null)
                    {
                        int runTime;

                        if (resultsByDataSetAndAligner[dataSet][aligner].replicas != null && resultsByDataSetAndAligner[dataSet][aligner].replicas.Count() >= 5)
                        {
                            runTime = (resultsByDataSetAndAligner[dataSet][aligner].replicas.Select(_ => _.alignTime + _.loadingTime).Sum() + timing.alignTime + timing.loadingTime) / (1 + resultsByDataSetAndAligner[dataSet][aligner].replicas.Count());
                        }
                        else
                        {
                            runTime = timing.alignTime + timing.loadingTime;
                        }
                        outputFile.Write("\t" + ASETools.DoubleRounded(((double)runTime / 3600 / 24), nDigits)); // Math converts seconds to days.
                    }
                    else
                    {
                        outputFile.Write("\t");
                    }// If we have timing

                } // aligner

                // Standard error for load & align
                foreach (var aligner in aligners)
                {
                    var timing = resultsByDataSetAndAligner[dataSet][aligner].runTiming;
                    var replicas = resultsByDataSetAndAligner[dataSet][aligner].replicas;

                    if (timing != null && replicas.Count() >= 4)
                    {
                        var times = replicas.Select(_ => ((double)_.alignTime + _.loadingTime) / 3600 /24).ToList();
 
                        var sigma = ASETools.StandardDeviationOfList(times);

                        outputFile.Write("\t" + ASETools.DoubleRounded(sigma / Math.Sqrt(times.Count()), nDigits));
                    } else
                    {
                        outputFile.Write("\t");
                    }
                }

                foreach (var aligner in aligners)
                {
                    var timing = resultsByDataSetAndAligner[dataSet][aligner].runTiming;
                    if (timing != null)
                    {
                        outputFile.Write("\t" + ASETools.DoubleRounded((double)timing.overallRuntime / 3600 / 24, nDigits)); // Math converts seconds to days.
                    }
                    else
                    {
                        outputFile.Write("\t");
                    }// If we have concordance
                } // aligner
                outputFile.WriteLine();
            } // dataSet

            outputFile.WriteLine();

            outputFile.WriteLine("\tLoad & Align" + oneTabPerAligner + "Whole Pipeline");
            outputFile.Write("Normalized run time");
            for (int i = 0; i < 2; i++)
            {
                foreach (var aligner in aligners)
                {
                    var alignerName = ASETools.alignerName[aligner];
                    outputFile.Write("\t" + alignerName);
                } // aligner
            }
            outputFile.WriteLine();

            foreach (var dataSet in dataSets)
            {
                outputFile.Write(dataSet);
                var snapTiming = resultsByDataSetAndAligner[dataSet][ASETools.Aligner.SNAP].runTiming;

                foreach (var aligner in aligners)
                {
                    var timing = resultsByDataSetAndAligner[dataSet][aligner].runTiming;
                    if (timing != null && snapTiming != null)
                    {
                        outputFile.Write("\t" + ASETools.DoubleRounded((double)(timing.alignTime + timing.loadingTime) / (snapTiming.alignTime + snapTiming.loadingTime), 2)); // Math converts seconds to days.
                    }
                    else
                    {
                        outputFile.Write("\t");
                    }// If we have concordance
                } // aligner

                foreach (var aligner in aligners)
                {
                    var timing = resultsByDataSetAndAligner[dataSet][aligner].runTiming;
                    if (timing != null && snapTiming != null)
                    {
                        outputFile.Write("\t" + ASETools.DoubleRounded((double)timing.overallRuntime / snapTiming.overallRuntime, 2)); // Math converts seconds to days.
                    }
                    else
                    {
                        outputFile.Write("\t");
                    }// If we have concordance
                } // aligner
                outputFile.WriteLine();
            } // dataSet

            //
            // Scatter graphs
            //

            var scatterGraphs = new List<ScatterGraphInfo>();
            scatterGraphs.Add(new ScatterGraphInfo("Mean F1", _ => _.Mean_F1_Score()));
            scatterGraphs.Add(new ScatterGraphInfo("SNV F1", _ => _.results[ASETools.VariantType.SNV].F1_score));
            scatterGraphs.Add(new ScatterGraphInfo("Indel F1", _ => _.results[ASETools.VariantType.Indel].F1_score));
            scatterGraphs.Add(new ScatterGraphInfo("SNV Recall", _ => _.results[ASETools.VariantType.SNV].recall));
            scatterGraphs.Add(new ScatterGraphInfo("SNV Precision", _ => _.results[ASETools.VariantType.SNV].precision));
            scatterGraphs.Add(new ScatterGraphInfo("Indel Recall", _ => _.results[ASETools.VariantType.Indel].recall));
            scatterGraphs.Add(new ScatterGraphInfo("Indel Precision", _ => _.results[ASETools.VariantType.Indel].precision));

            foreach (var dataSet in dataSets)
            {
                //
                // Write the graph headers
                //
                outputFile.WriteLine();
                outputFile.WriteLine("Scatter graphs for " + dataSet);
                foreach (var scatterGraph in scatterGraphs)
                {
                    outputFile.Write("\t" + scatterGraph.name + "\t");
                }
                outputFile.WriteLine();

                outputFile.Write("Aligner");
                foreach (var scatterGraph in scatterGraphs)
                {
                    outputFile.Write("\tRun Time\tConcordance"); ;
                }
                outputFile.WriteLine();

                foreach (var aligner in aligners)
                {
                    outputFile.Write(ASETools.alignerName[aligner]);
                    foreach (var scatterGraph in scatterGraphs)
                    {
                        if (resultsByDataSetAndAligner[dataSet][aligner].runTiming != null && resultsByDataSetAndAligner[dataSet][aligner].concordance != null)
                        {
                            outputFile.Write("\t" + (double)(resultsByDataSetAndAligner[dataSet][aligner].runTiming.alignTime + resultsByDataSetAndAligner[dataSet][aligner].runTiming.loadingTime) / 3600 / 24 + "\t" +
                                ASETools.DoubleRounded(scatterGraph.extractor(resultsByDataSetAndAligner[dataSet][aligner].concordance), nDigits));
                        } else
                        {
                            outputFile.Write("\t\t");
                        }
                    } // scatter graph
                    outputFile.WriteLine();
                } // aligner
            }

            RunSeedSizeSweep(outputFile, "hg003");
            RunSeedSizeSweep(outputFile, "hg007");
            RunSeedSizeSweep(outputFile, "ERR194147");

            outputFile.WriteLine("**done**");
            outputFile.Close();
        } // Main
    }
}
