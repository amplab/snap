using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using ASELib;
using System.IO;
using System.Diagnostics;

namespace MethylationDistributions
{
    class Program
    {
        class BetaValues
        {
            Dictionary<string, Dictionary<int, ASETools.RunningMeanAndStdDev>> betas = new Dictionary<string, Dictionary<int, ASETools.RunningMeanAndStdDev>>();

            public void merge (BetaValues peer)
            {
                foreach (var contigEntry in peer.betas)
                {
                    var contig = contigEntry.Key;
                    var loci = contigEntry.Value;

                    if (!betas.ContainsKey(contig))
                    {
                        betas.Add(contig, new Dictionary<int, ASETools.RunningMeanAndStdDev>());
                    }

                    foreach (var locusEntry in loci)
                    {
                        var locus = locusEntry.Key;
                        var runningMeanAndStdDev = locusEntry.Value;

                        if (!betas[contig].ContainsKey(locus))
                        {
                            betas[contig].Add(locus, runningMeanAndStdDev);
                        } else
                        {
                            betas[contig][locus].merge(runningMeanAndStdDev);
                        }
                    } // locus
                } // contig
            } // merge

            public void addValue(string contig, int locus, double value)
            {
                if (!betas.ContainsKey(contig))
                {
                    betas.Add(contig, new Dictionary<int, ASETools.RunningMeanAndStdDev>());
                }

                if (!betas[contig].ContainsKey(locus))
                {
                    betas[contig].Add(locus, new ASETools.RunningMeanAndStdDev());
                }

                betas[contig][locus].addValue(value);
            }

            public void WriteToFile(StreamWriter outputFile)
            {
                foreach (var contigEntry in betas)
                {
                    var contig = contigEntry.Key;
                    foreach (var locusEntry in contigEntry.Value)
                    {
                        var locus = locusEntry.Key;
                        var runningMeanAndStdDev = locusEntry.Value;

                        var meanAndStdDev = runningMeanAndStdDev.getMeanAndStdDev();

                        outputFile.WriteLine(contig + "\t" + locus + "\t" + runningMeanAndStdDev.getCount() + "\t" + runningMeanAndStdDev.getMin() + "\t" + runningMeanAndStdDev.getMax() + "\t" + meanAndStdDev.mean + "\t" + meanAndStdDev.stddev);
                    }
                }
            }
        } // BetaValues

        static BetaValues overallBetaValues = new BetaValues();

        static ASETools.CommonData commonData;
        static void Main(string[] args)
        {
            commonData = ASETools.CommonData.LoadCommonData(args);

            if (null == commonData)
            {
                return;
            }

            if (commonData.cases.Any(x => x.Value.tumor_methylation_filename == "" && x.Value.tumor_methylation_file_id != ""))
            {
                Console.WriteLine("At least one case doesn't have its methylation data downloaded.");
                return;
            }

            var casesToProcess = commonData.cases.Select(_ => _.Value).Where(_ => _.tumor_methylation_filename != "").ToList();

            Console.WriteLine("Processing " + casesToProcess.Count() + " cases, one dot/100.");
            ASETools.PrintNumberBar(casesToProcess.Count() / 100);

            var threading = new ASETools.WorkerThreadHelper<ASETools.Case, BetaValues>(casesToProcess, HandleOneCase, FinishThread, null, 100);
            threading.run();

            Console.WriteLine("Computing and writing distributions");

            var outputFile = ASETools.CreateStreamWriterWithRetry(commonData.configuration.baseDirectory + ASETools.MethylationDistributionsFilename);
            if (outputFile == null)
            {
                Console.WriteLine("Unable to open output file " + commonData.configuration.baseDirectory + ASETools.MethylationDistributionsFilename);
                return;
            }

            outputFile.WriteLine("contig\tlocus\tn\tmin beta\tmax beta\tmean beta\tbeta standard deviation");
            overallBetaValues.WriteToFile(outputFile); 

            outputFile.WriteLine("**done**");
            outputFile.Close();

            Console.WriteLine(ASETools.ElapsedTimeInSeconds(commonData.timer));
        }

        static void HandleOneCase(ASETools.Case case_, BetaValues state)
        {
            var methylationLines = ASETools.MethylationAnnotationLine.ReadFile(case_.tumor_methylation_filename, case_.tumor_methylation_file_id, false);
            if (null == methylationLines)
            {
                throw new Exception("Failed to read methylation data for case " + case_.case_id);
            }

            foreach (var methylationLine in methylationLines)
            {
                state.addValue(methylationLine.compositeREF.Chromosome, methylationLine.compositeREF.Start, methylationLine.Beta_Value);
            }
        }

        static void FinishThread(BetaValues state)
        {
            lock (overallBetaValues)
            {
                overallBetaValues.merge(state);
            }
        }
    }
}
