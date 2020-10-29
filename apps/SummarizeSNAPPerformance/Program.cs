using System;
using System.Collections.Generic;
using System.Diagnostics;
using System.Linq;
using System.Runtime.InteropServices;
using System.Runtime.Remoting;
using System.Text;
using System.Threading;
using System.Threading.Tasks;
using ASELib;

namespace SummarizeSNAPPerformance
{
    class Program
    {

        static Dictionary<string, ASETools.CaseMetadata> metadataByCase;

        static int[] speedByReadCountMins = { 0, 60000000, 70000000, 80000000, 90000000, 100000000, 110000000, 120000000, 130000000, 140000000, 150000000, 160000000, 170000000, 180000000, 190000000, 200000000 };

        class Histograms
        {
            public ASETools.PreBucketedHistogram alignTimeHistogram = new ASETools.PreBucketedHistogram(0, ASETools.maxTimeForRealignHistograms, ASETools.timeStepForRealignHistograms);
            public ASETools.PreBucketedHistogram indexLoadTimeHistogram = new ASETools.PreBucketedHistogram(0, ASETools.maxTimeForRealignHistograms, ASETools.timeStepForRealignHistograms);
            public ASETools.PreBucketedHistogram loadAndAlignTimeHistogram = new ASETools.PreBucketedHistogram(0, ASETools.maxTimeForRealignHistograms, ASETools.timeStepForRealignHistograms);
            public ASETools.PreBucketedHistogram sortTimeHistogram = new ASETools.PreBucketedHistogram(0, ASETools.maxTimeForRealignHistograms, ASETools.timeStepForRealignHistograms);
            public ASETools.PreBucketedHistogram totalSNAPTimeHistogram = new ASETools.PreBucketedHistogram(0, ASETools.maxTimeForRealignHistograms, ASETools.timeStepForRealignHistograms);
            public ASETools.PreBucketedHistogram copyInTimeHistogram = new ASETools.PreBucketedHistogram(0, ASETools.maxTimeForRealignHistograms, ASETools.timeStepForRealignHistograms);
            public ASETools.PreBucketedHistogram copyOutTimeHistogram = new ASETools.PreBucketedHistogram(0, ASETools.maxTimeForRealignHistograms, ASETools.timeStepForRealignHistograms);
            public ASETools.PreBucketedHistogram endToEndTimeHistogram = new ASETools.PreBucketedHistogram(0, ASETools.maxTimeForRealignHistograms, ASETools.timeStepForRealignHistograms);
            public ASETools.PreBucketedHistogram readsPerSecondReportedHistogram = new ASETools.PreBucketedHistogram(0, ASETools.maxSpeedForRealignmentHistograms, ASETools.stepForSpeedRealignmentHistograms);
            public ASETools.PreBucketedHistogram readsPerLoadAndAlignSecondHistogram = new ASETools.PreBucketedHistogram(0, ASETools.maxSpeedForRealignmentHistograms, ASETools.stepForSpeedRealignmentHistograms);
            public ASETools.PreBucketedHistogram readsPerSecondAllSnapHistogram = new ASETools.PreBucketedHistogram(0, ASETools.maxSpeedForRealignmentHistograms, ASETools.stepForSpeedRealignmentHistograms);
            public ASETools.PreBucketedHistogram totalReadsHistogram = new ASETools.PreBucketedHistogram(0, 1000000000, 1000000);
            public ASETools.PreBucketedHistogram[] speedByReadCount = new ASETools.PreBucketedHistogram[speedByReadCountMins.Count()];

            public Histograms()
            {
                for (int i = 0; i < speedByReadCountMins.Count(); i++)
                {
                    speedByReadCount[i] = new ASETools.PreBucketedHistogram(0, ASETools.maxSpeedForRealignmentHistograms, ASETools.stepForSpeedRealignmentHistograms);
                }
            }

            public int timeForSlowestRun = 0;
            public string caseIdForSlowestRun = "";

            public void merge(Histograms peer)
            {
                alignTimeHistogram.merge(peer.alignTimeHistogram);
                indexLoadTimeHistogram.merge(peer.indexLoadTimeHistogram);
                loadAndAlignTimeHistogram.merge(peer.loadAndAlignTimeHistogram);
                sortTimeHistogram.merge(peer.sortTimeHistogram);
                totalSNAPTimeHistogram.merge(peer.totalSNAPTimeHistogram);
                copyInTimeHistogram.merge(peer.copyInTimeHistogram);
                copyOutTimeHistogram.merge(peer.copyOutTimeHistogram);
                endToEndTimeHistogram.merge(peer.endToEndTimeHistogram);
                readsPerSecondReportedHistogram.merge(peer.readsPerSecondReportedHistogram);
                readsPerLoadAndAlignSecondHistogram.merge(peer.readsPerLoadAndAlignSecondHistogram);
                readsPerSecondAllSnapHistogram.merge(peer.readsPerSecondAllSnapHistogram);
                totalReadsHistogram.merge(peer.totalReadsHistogram);

                if (timeForSlowestRun < peer.timeForSlowestRun)
                {
                    timeForSlowestRun = peer.timeForSlowestRun;
                    caseIdForSlowestRun = peer.caseIdForSlowestRun;
                }

                for (int i = 0; i < speedByReadCountMins.Count(); i++)
                {
                    speedByReadCount[i].merge(peer.speedByReadCount[i]);
                }

            }
        } // Histograms

        class PerThreadState
        {
            public Dictionary<bool, Dictionary<bool, Histograms>> histograms = new Dictionary<bool, Dictionary<bool, Histograms>>();    // tumor -> paired -> histograms

            public PerThreadState()
            {
                foreach (var tumor in ASETools.BothBools)
                {
                    histograms.Add(tumor, new Dictionary<bool, Histograms>());
                    foreach (var paired in ASETools.BothBools)
                    {
                        histograms[tumor].Add(paired, new Histograms());
                    }
                }
            } // ctor

            public void merge(PerThreadState peer)
            {
                foreach (var tumor in ASETools.BothBools)
                {
                    foreach (var paired in ASETools.BothBools)
                    {
                        histograms[tumor][paired].merge(peer.histograms[tumor][paired]);
                    }
                }
            }
        } // PerThreadState

        static PerThreadState globalState = new PerThreadState();
        static void ProcessOneCase(ASETools.Case case_, PerThreadState state)
        {
            foreach (var tumor in ASETools.BothBools)
            {
                var inputFilename = case_.realignments[ASETools.Aligner.SNAP][tumor].dna_statistics_filename;
 
                if (inputFilename == "")
                {
                    continue;
                }

                var paired = metadataByCase[case_.case_id].getBAMMetadata(tumor, true).isPaired;

                var stats = ASETools.RunTiming.LoadFromSNAPFile(inputFilename);

                state.histograms[tumor][paired].alignTimeHistogram.addValue(stats.alignTime);
                state.histograms[tumor][paired].indexLoadTimeHistogram.addValue(stats.loadingTime);
                state.histograms[tumor][paired].loadAndAlignTimeHistogram.addValue(stats.alignTime + stats.loadingTime);
                state.histograms[tumor][paired].sortTimeHistogram.addValue(stats.sortTime);
                state.histograms[tumor][paired].totalSNAPTimeHistogram.addValue(stats.overallRuntime);
                state.histograms[tumor][paired].copyInTimeHistogram.addValue(stats.copyInTime);
                state.histograms[tumor][paired].copyOutTimeHistogram.addValue(stats.copyOutTime);
                state.histograms[tumor][paired].endToEndTimeHistogram.addValue(stats.overallRuntime + stats.copyInTime + stats.copyOutTime);
                state.histograms[tumor][paired].readsPerSecondReportedHistogram.addValue(stats.readsPerSecond);
                state.histograms[tumor][paired].readsPerLoadAndAlignSecondHistogram.addValue(stats.totalReads / (stats.alignTime + stats.loadingTime));
                state.histograms[tumor][paired].readsPerSecondAllSnapHistogram.addValue(stats.totalReads / stats.overallRuntime);
                state.histograms[tumor][paired].totalReadsHistogram.addValue(stats.totalReads);

                if (stats.overallRuntime > state.histograms[tumor][paired].timeForSlowestRun)
                {
                    state.histograms[tumor][paired].timeForSlowestRun = stats.overallRuntime;
                    state.histograms[tumor][paired].caseIdForSlowestRun = case_.case_id;
                }

                int whichSpeedBucket = 0;
                for (int i = 0; i < speedByReadCountMins.Count(); i++)
                {
                    if (stats.totalReads > speedByReadCountMins[i])
                    {
                        whichSpeedBucket = i;
                    }
                }

                state.histograms[tumor][paired].speedByReadCount[whichSpeedBucket].addValue(stats.readsPerSecond);
            } // tumor/normal
        } // ProcessOneCase

        static void FinishUp(PerThreadState state)
        {
            lock (globalState)
            {
                globalState.merge(state);
            } // lock(globalState)
        } // FinishUp

        static void Main(string[] args)
        {
            var timer = new Stopwatch();
            timer.Start();

            var configuration = ASETools.Configuration.loadFromFile(args);
            if (configuration == null)
            {
                return;
            }

            var cases = ASETools.Case.LoadCases(configuration.casesFilePathname);
            if (cases == null)
            {
                Console.WriteLine("Unable to load cases from " + configuration.casesFilePathname);
                return;
            }

            var listOfCases = cases.Select(_ => _.Value).ToList();

            if (listOfCases.Any(case_ => case_.case_metadata_filename == "" || case_.realignments[ASETools.Aligner.SNAP].Any(_ => _.Value.dna_statistics_filename == "")))
            {
                Console.WriteLine("Some cases are missing data.");
                //BJB return;
            }

            var casesToRun = listOfCases.Where(case_ => case_.case_metadata_filename != "" && case_.realignments[ASETools.Aligner.SNAP].Any(_ => _.Value.dna_statistics_filename != "")).ToList();

            metadataByCase = ASETools.CaseMetadata.ReadConsolodatedCaseMetadata(configuration.finalResultsDirectory + ASETools.ConsolodatedCaseMetadataFilename);

            int nPerDot;
            ASETools.PrintMessageAndNumberBar("Processing", "cases", casesToRun.Count(), out nPerDot);

            var threading = new ASETools.WorkerThreadHelper<ASETools.Case, PerThreadState>(casesToRun, ProcessOneCase, FinishUp, null, nPerDot);
            threading.run();

            var outputFilename = configuration.finalResultsDirectory + ASETools.SNAPSummaryFilename;

            var outputFile = ASETools.CreateStreamWriterWithRetry(outputFilename);
            if (outputFile == null)
            {
                return;
            }

            foreach (var paired in ASETools.BothBools)
            {
                foreach (var tumor in ASETools.BothBools)
                {
                    var timeHistograms = new List<KeyValuePair<string, ASETools.PreBucketedHistogram>>();
                    timeHistograms.Add(new KeyValuePair<string, ASETools.PreBucketedHistogram>("align time", globalState.histograms[tumor][paired].alignTimeHistogram));
                    timeHistograms.Add(new KeyValuePair<string, ASETools.PreBucketedHistogram>("index load", globalState.histograms[tumor][paired].indexLoadTimeHistogram));
                    timeHistograms.Add(new KeyValuePair<string, ASETools.PreBucketedHistogram>("load and align", globalState.histograms[tumor][paired].loadAndAlignTimeHistogram));
                    timeHistograms.Add(new KeyValuePair<string, ASETools.PreBucketedHistogram>("sort time", globalState.histograms[tumor][paired].sortTimeHistogram));
                    timeHistograms.Add(new KeyValuePair<string, ASETools.PreBucketedHistogram>("total SNAP time", globalState.histograms[tumor][paired].totalSNAPTimeHistogram));
                    timeHistograms.Add(new KeyValuePair<string, ASETools.PreBucketedHistogram>("copy in time", globalState.histograms[tumor][paired].copyInTimeHistogram));
                    timeHistograms.Add(new KeyValuePair<string, ASETools.PreBucketedHistogram>("copy out time", globalState.histograms[tumor][paired].copyOutTimeHistogram));
                    timeHistograms.Add(new KeyValuePair<string, ASETools.PreBucketedHistogram>("end-to-end time", globalState.histograms[tumor][paired].endToEndTimeHistogram));

                    var rateHistograms = new List<KeyValuePair<string, ASETools.PreBucketedHistogram>>();
                    rateHistograms.Add(new KeyValuePair<string, ASETools.PreBucketedHistogram>("reads per second reported", globalState.histograms[tumor][paired].readsPerSecondReportedHistogram));
                    rateHistograms.Add(new KeyValuePair<string, ASETools.PreBucketedHistogram>("reads per second of load and align time", globalState.histograms[tumor][paired].readsPerLoadAndAlignSecondHistogram));
                    rateHistograms.Add(new KeyValuePair<string, ASETools.PreBucketedHistogram>("reads per second in SNAP", globalState.histograms[tumor][paired].readsPerSecondAllSnapHistogram));

                    outputFile.WriteLine("Data for " + (tumor ? "tumor" : "normal") + " " + (paired ? "paired-end" : "single-end"));
                    outputFile.WriteLine("Slowest run " + globalState.histograms[tumor][paired].timeForSlowestRun + " for " + globalState.histograms[tumor][paired].caseIdForSlowestRun);
                    timeHistograms.ForEach(_ => outputFile.WriteLine(_.Key + " mean " + _.Value.mean()));
                    rateHistograms.ForEach(_ => outputFile.WriteLine(_.Key + " mean " + _.Value.mean()));

                    outputFile.WriteLine();
                    timeHistograms.ForEach(_ => outputFile.WriteLine(_.Key + " max " + _.Value.max()));
                    rateHistograms.ForEach(_ => outputFile.WriteLine(_.Key + " max " + _.Value.max()));

                    outputFile.WriteLine();
                    ASETools.PreBucketedHistogram.WriteBatchOfHistogramCDFs(outputFile, timeHistograms);

                    outputFile.WriteLine();
                    ASETools.PreBucketedHistogram.WriteBatchOfHistogramCDFs(outputFile, rateHistograms);

                    outputFile.WriteLine();
                    outputFile.WriteLine("Total reads per sample");
                    globalState.histograms[tumor][paired].totalReadsHistogram.WriteHistogram(outputFile);

                    var speedByReadCountHistograms = new List<KeyValuePair<string, ASETools.PreBucketedHistogram>>();
                    for (int i = 0; i < speedByReadCountMins.Count(); i++)
                    {
                        speedByReadCountHistograms.Add(new KeyValuePair<string, ASETools.PreBucketedHistogram>(speedByReadCountMins[i].ToString() + ((i == speedByReadCountMins.Count() - 1) ? "+" : "-" + speedByReadCountMins[i + 1]), globalState.histograms[tumor][paired].speedByReadCount[i]));
                    }

                    outputFile.WriteLine();
                    outputFile.WriteLine("Reads/s reported in SNAP by size of input file.");
                    ASETools.PreBucketedHistogram.WriteBatchOfHistogramCDFs(outputFile, speedByReadCountHistograms);

                    outputFile.WriteLine();
                } // tumor/normal
            } // paired

            outputFile.WriteLine("**done**");
            outputFile.Close();

            Console.WriteLine("Processed " + casesToRun.Count() + " cases in " + ASETools.ElapsedTimeInSeconds(timer));

        } // Main
    } // Program
} // SummarizeSNAPPerformance
