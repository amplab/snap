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


        static int nDigits = 4;

        static void Main(string[] args)
        {
            var timingsDirectory = @"d:\temp\timings\";
            var concordanceDirectory = @"d:\temp\concordance\";

            string[] dataSets = { "hg001", "hg002", "hg003", "hg004", "hg005", "hg006", "hg007", "ERR194146", "ERR194147", "mp002", "mp003", "mp004", "mp005", "mp006", "mp007",  };

            Dictionary<string, Dictionary<ASETools.Aligner, Result>> resultsByDataSetAndAligner = new Dictionary<string, Dictionary<ASETools.Aligner, Result>>();
            Dictionary<string, Dictionary<ASETools.Aligner, List<Result>>> replicasByDataSetAndAligner = new Dictionary<string, Dictionary<ASETools.Aligner, List<Result>>>();


            int nTotal = 0;
            int nWithTimings = 0;
            int nWithConcordance = 0;

            foreach (var dataSet in dataSets) {
                resultsByDataSetAndAligner.Add(dataSet, new Dictionary<ASETools.Aligner, Result>());

                foreach (var aligner in ASETools.allAligners)
                {
                    nTotal++;

                    var timingFilename = timingsDirectory + dataSet + "." + ASETools.alignerName[aligner] + "_timings.txt";
                    var splitTimingsFilenames = new List<string>();
                    for (int i = 0; i < 9; i++)
                    {
                        splitTimingsFilenames.Add(timingsDirectory + dataSet + "." + i + "." + ASETools.alignerName[aligner] + "_timings.txt");
                    }

                    if (File.Exists(timingFilename) || splitTimingsFilenames.All(_ => File.Exists(_))) {
                        nWithTimings++;

                        var replicas = new List<ASETools.RunTiming>();

                        ASETools.RunTiming timings;
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
                                timings = ASETools.RunTiming.LoadFromSplitFiles(splitTimingsFilenames);
                            }
                        }

                        var concordanceFilename = concordanceDirectory + dataSet + "." + ASETools.alignerName[aligner].ToLower() + ".concordance.tar";
                        //Console.WriteLine(concordanceFilename);
                        ASETools.ConcordanceResults concordance = null;

                        if (File.Exists(concordanceFilename))
                        {
                            nWithConcordance++;
                            concordance = new ASETools.ConcordanceResults(concordanceFilename);
                        }

                        foreach (var replicaFilename in Directory.GetFiles(timingsDirectory, dataSet + "." + ASETools.alignerName[aligner].ToLower() + "_timings_*.txt"))
                        {
                            var replica = ASETools.RunTiming.LoadFromReplicaFile(replicaFilename);
                            if (replica != null)
                            {
                                replicas.Add(replica);
                            }
                        }

                        resultsByDataSetAndAligner[dataSet].Add(aligner, new Result(timings, concordance, replicas));
                    } else
                    {
                        resultsByDataSetAndAligner[dataSet].Add(aligner, new Result(null, null, null));
                    }


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
            outputFile.WriteLine("\tSNV Recall\t\t\t\t\tSNV Precision\t\t\t\t\tIndel Recall\t\t\t\t\tIndel Precision\t\t\t\t\tSNV F1\t\t\t\t\tIndel F1\t\t\t\t\tMean F1");

            for (int i = 0; i < 7; i++)
            {
                foreach (var aligner in ASETools.allAligners)
                {
                    outputFile.Write("\t" + ASETools.alignerName[aligner]);
                }
            }

            outputFile.WriteLine();

            foreach (var dataSet in dataSets)
            {
                outputFile.Write(dataSet);
                foreach (var aligner in ASETools.allAligners)
                {
                    outputFile.Write("\t");
                    var concordance = resultsByDataSetAndAligner[dataSet][aligner].concordance;
                    if (concordance != null)
                    {
                        outputFile.Write(ASETools.DoubleRounded(resultsByDataSetAndAligner[dataSet][aligner].concordance.results[ASETools.VariantType.SNV].recall, nDigits));
                    }
                } // aligner

                foreach (var aligner in ASETools.allAligners)
                {
                    outputFile.Write("\t");
                    var concordance = resultsByDataSetAndAligner[dataSet][aligner].concordance;
                    if (concordance != null)
                    {
                        outputFile.Write(ASETools.DoubleRounded(resultsByDataSetAndAligner[dataSet][aligner].concordance.results[ASETools.VariantType.SNV].precision, nDigits));
                    }
                } // aligner

                foreach (var aligner in ASETools.allAligners)
                {
                    outputFile.Write("\t");
                    var concordance = resultsByDataSetAndAligner[dataSet][aligner].concordance;
                    if (concordance != null)
                    {
                        outputFile.Write(ASETools.DoubleRounded(resultsByDataSetAndAligner[dataSet][aligner].concordance.results[ASETools.VariantType.Indel].recall, nDigits));
                    }
                } // aligner

                foreach (var aligner in ASETools.allAligners)
                {
                    outputFile.Write("\t");
                    var concordance = resultsByDataSetAndAligner[dataSet][aligner].concordance;
                    if (concordance != null)
                    {
                        outputFile.Write(ASETools.DoubleRounded(resultsByDataSetAndAligner[dataSet][aligner].concordance.results[ASETools.VariantType.Indel].precision, nDigits));
                    }
                } // aligner

                foreach (var aligner in ASETools.allAligners)
                {
                    outputFile.Write("\t");
                    var concordance = resultsByDataSetAndAligner[dataSet][aligner].concordance;
                    if (concordance != null)
                    {
                        outputFile.Write(ASETools.DoubleRounded(resultsByDataSetAndAligner[dataSet][aligner].concordance.results[ASETools.VariantType.SNV].F1_score, nDigits));
                    }
                } // aligner

                foreach (var aligner in ASETools.allAligners)
                {
                    outputFile.Write("\t");
                    var concordance = resultsByDataSetAndAligner[dataSet][aligner].concordance;
                    if (concordance != null)
                    {
                        outputFile.Write(ASETools.DoubleRounded(resultsByDataSetAndAligner[dataSet][aligner].concordance.results[ASETools.VariantType.Indel].F1_score, nDigits));
                    }
                } // aligner

                foreach (var aligner in ASETools.allAligners)
                {
                    outputFile.Write("\t");
                    var concordance = resultsByDataSetAndAligner[dataSet][aligner].concordance;
                    if (concordance != null)
                    {
                        outputFile.Write(ASETools.DoubleRounded(resultsByDataSetAndAligner[dataSet][aligner].concordance.Mean_F1_Score(), nDigits));
                    }
                } // aligner

                outputFile.WriteLine();
            } // dataSet

            outputFile.WriteLine();

            outputFile.WriteLine("\tSNV\t\t\t\t\tIndel\t\t\t\t\tMean");
            outputFile.Write("9s");

            // 9s (= -1Log1(1-F1 score))
            for (int i = 0; i < 3; i++)
            {
                foreach (var aligner in ASETools.allAligners)
                {
                    var alignerName = ASETools.alignerName[aligner];
                    outputFile.Write("\t" + alignerName);
                } // aligner
            }
            outputFile.WriteLine();

            foreach (var dataSet in dataSets)
            {
                outputFile.Write(dataSet);
                foreach (var aligner in ASETools.allAligners)
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

                foreach (var aligner in ASETools.allAligners)
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

                foreach (var aligner in ASETools.allAligners)
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

            outputFile.WriteLine("\tLoad & Align\t\t\t\t\tLoad & Align error bar\t\t\t\t\tWhole Pipeline");
            outputFile.Write("Run Time");
            for (int i = 0; i < 3; i++)
            {
                foreach (var aligner in ASETools.allAligners)
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
                foreach (var aligner in ASETools.allAligners)
                {
                    var timing = resultsByDataSetAndAligner[dataSet][aligner].runTiming;
                    if (timing != null)
                    {
                        int runTime;

                        if (resultsByDataSetAndAligner[dataSet][aligner].replicas != null)
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
                foreach (var aligner in ASETools.allAligners)
                {
                    var timing = resultsByDataSetAndAligner[dataSet][aligner].runTiming;
                    var replicas = resultsByDataSetAndAligner[dataSet][aligner].replicas;

                    if (timing != null && replicas.Count() >= 4)
                    {
                        var times = replicas.Select(_ => ((double)_.alignTime + _.loadingTime) / 3600 /24).ToList();
                        times.Add(((double)timing.alignTime + timing.loadingTime) / 3600 / 24);

                        var sigma = ASETools.StandardDeviationOfList(times);

                        outputFile.Write("\t" + ASETools.DoubleRounded(sigma / Math.Sqrt(times.Count()), nDigits));
                    } else
                    {
                        outputFile.Write("\t");
                    }
                }

                foreach (var aligner in ASETools.allAligners)
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

            outputFile.WriteLine("\tLoad & Align\t\t\t\t\tWhole Pipeline");
            outputFile.Write("Normalized run time");
            for (int i = 0; i < 2; i++)
            {
                foreach (var aligner in ASETools.allAligners)
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

                foreach (var aligner in ASETools.allAligners)
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

                foreach (var aligner in ASETools.allAligners)
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
            outputFile.WriteLine();
            outputFile.WriteLine("Scatter Graphs");
            string[] scatterGraphTypes = { "SNV Recall", "SNV Precision", "SNV F1", "Indel Recall", "Indel Precision", "Indel F1", "Mean F1" };
            foreach (var scatterGraphType in scatterGraphTypes)
            {
                outputFile.Write("\t" + scatterGraphType + "\t\t\t\t");                
            } // scatter graph type
            outputFile.WriteLine();

            outputFile.Write("Sample");
            foreach (var scatterGraphType in scatterGraphTypes)
            {
                foreach (var aligner in ASETools.allAligners)
                {

                } // aligners
            } // scatter graph type
            outputFile.WriteLine();

            outputFile.WriteLine("**done**");
            outputFile.Close();
        } // Main
    }
}
