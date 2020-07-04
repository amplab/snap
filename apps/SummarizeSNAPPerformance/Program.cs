using System;
using System.Collections.Generic;
using System.Diagnostics;
using System.Linq;
using System.Runtime.InteropServices;
using System.Text;
using System.Threading;
using System.Threading.Tasks;
using ASELib;

namespace SummarizeSNAPPerformance
{
    class Program
    {

        class PerThreadState
        {
            public ASETools.PreBucketedHistogram alignTimeHistogram = new ASETools.PreBucketedHistogram(0, ASETools.maxTimeForRealignHistograms, 60);
            public ASETools.PreBucketedHistogram indexLoadTimeHistogram = new ASETools.PreBucketedHistogram(0, ASETools.maxTimeForRealignHistograms, 60);
            public ASETools.PreBucketedHistogram loadAndAlignTimeHistogram = new ASETools.PreBucketedHistogram(0, ASETools.maxTimeForRealignHistograms, 60);
            public ASETools.PreBucketedHistogram sortTimeHistogram = new ASETools.PreBucketedHistogram(0, ASETools.maxTimeForRealignHistograms, 60);
            public ASETools.PreBucketedHistogram totalSNAPTimeHistogram = new ASETools.PreBucketedHistogram(0, ASETools.maxTimeForRealignHistograms, 60);
            public ASETools.PreBucketedHistogram copyInTimeHistogram = new ASETools.PreBucketedHistogram(0, ASETools.maxTimeForRealignHistograms, 60);
            public ASETools.PreBucketedHistogram copyOutTimeHistogram = new ASETools.PreBucketedHistogram(0, ASETools.maxTimeForRealignHistograms, 60);
            public ASETools.PreBucketedHistogram endToEndTimeHistogram = new ASETools.PreBucketedHistogram(0, ASETools.maxTimeForRealignHistograms, 60);
            public ASETools.PreBucketedHistogram readsPerSecondReportedHistogram = new ASETools.PreBucketedHistogram(0, ASETools.maxSpeedForRealignmentHistograms, ASETools.stepForSpeedRealignmentHistograms);
            public ASETools.PreBucketedHistogram readsPerSecondAllSnapHistogram = new ASETools.PreBucketedHistogram(0, ASETools.maxSpeedForRealignmentHistograms, ASETools.stepForSpeedRealignmentHistograms);

            public int timeForSlowestRun = 0;
            public string caseIdForSlowestRun = "";

            public void merge(PerThreadState peer)
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
                readsPerSecondAllSnapHistogram.merge(peer.readsPerSecondAllSnapHistogram);

                if (timeForSlowestRun < peer.timeForSlowestRun)
                {
                    timeForSlowestRun = peer.timeForSlowestRun;
                    caseIdForSlowestRun = peer.caseIdForSlowestRun;
                }
            }
        }

        static Dictionary<bool, PerThreadState> globalState = new Dictionary<bool, PerThreadState>();   // tumor->state

        static void ProcessOneCase(ASETools.Case case_, Dictionary<bool, PerThreadState> state)
        {
            foreach (var tumor in ASETools.BothBools)
            {
                if (!state.ContainsKey(tumor))
                {
                    state.Add(tumor, new PerThreadState());
                }

                var inputFilename = tumor ? case_.snap_realigned_tumor_dna_statictics_filename : case_.snap_realigned_normal_dna_statictics_filename;

                if (inputFilename == "")
                {
                    continue;
                }

                var stats = ASETools.SNAPRunTiming.LoadFromFile(inputFilename);

                state[tumor].alignTimeHistogram.addValue(stats.alignTime);
                state[tumor].indexLoadTimeHistogram.addValue(stats.loadingTime);
                state[tumor].loadAndAlignTimeHistogram.addValue(stats.alignTime + stats.loadingTime);
                state[tumor].sortTimeHistogram.addValue(stats.sortTime);
                state[tumor].totalSNAPTimeHistogram.addValue(stats.overallRuntime);
                state[tumor].copyInTimeHistogram.addValue(stats.copyInTime);
                state[tumor].copyOutTimeHistogram.addValue(stats.copyOutTime);
                state[tumor].endToEndTimeHistogram.addValue(stats.overallRuntime + stats.copyInTime + stats.copyOutTime);
                state[tumor].readsPerSecondReportedHistogram.addValue(stats.readsPerSecond);
                state[tumor].readsPerSecondAllSnapHistogram.addValue(stats.totalReads / stats.overallRuntime);

                if (stats.overallRuntime > state[tumor].timeForSlowestRun)
                {
                    state[tumor].timeForSlowestRun = stats.overallRuntime;
                    state[tumor].caseIdForSlowestRun = case_.case_id;
                }
            } // tumor/normal
        } // ProcessOneCase

        static void FinishUp(Dictionary<bool, PerThreadState> state)
        {
            lock (globalState)
            {
                foreach (var tumor in ASETools.BothBools)
                {
                    globalState[tumor].merge(state[tumor]);
                }
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

            foreach (var tumor in ASETools.BothBools)
            {
                globalState.Add(tumor, new PerThreadState());
            }

            if (listOfCases.Any(_ => _.snap_realigned_normal_dna_statictics_filename == "" || _.snap_realigned_tumor_dna_statictics_filename == ""))
            {
                Console.WriteLine("Some cases are missing data.");
                //BJB return;
            }

            var casesToRun = listOfCases.Where(_ => _.snap_realigned_normal_dna_statictics_filename != "" || _.snap_realigned_tumor_dna_statictics_filename != "").ToList();

            int nPerDot;
            ASETools.PrintMessageAndNumberBar("Processing", "cases", casesToRun.Count(), out nPerDot);

            var threading = new ASETools.WorkerThreadHelper<ASETools.Case, Dictionary<bool, PerThreadState>>(casesToRun, ProcessOneCase, FinishUp, null, nPerDot);
            threading.run();

            var outputFilename = configuration.finalResultsDirectory + ASETools.SNAPSummaryFilename;

            var outputFile = ASETools.CreateStreamWriterWithRetry(outputFilename);
            if (outputFile == null)
            {
                return;
            }

            foreach (var tumor in ASETools.BothBools)
            {
                var timeHistograms = new List<KeyValuePair<string, ASETools.PreBucketedHistogram>>();
                timeHistograms.Add(new KeyValuePair<string, ASETools.PreBucketedHistogram>("align time", globalState[tumor].alignTimeHistogram));
                timeHistograms.Add(new KeyValuePair<string, ASETools.PreBucketedHistogram>("index load", globalState[tumor].indexLoadTimeHistogram));
                timeHistograms.Add(new KeyValuePair<string, ASETools.PreBucketedHistogram>("load and align", globalState[tumor].loadAndAlignTimeHistogram));
                timeHistograms.Add(new KeyValuePair<string, ASETools.PreBucketedHistogram>("sort time", globalState[tumor].sortTimeHistogram));
                timeHistograms.Add(new KeyValuePair<string, ASETools.PreBucketedHistogram>("total SNAP time", globalState[tumor].totalSNAPTimeHistogram));
                timeHistograms.Add(new KeyValuePair<string, ASETools.PreBucketedHistogram>("copy in time", globalState[tumor].copyInTimeHistogram));
                timeHistograms.Add(new KeyValuePair<string, ASETools.PreBucketedHistogram>("copy out time", globalState[tumor].copyOutTimeHistogram));
                timeHistograms.Add(new KeyValuePair<string, ASETools.PreBucketedHistogram>("end-to-end time", globalState[tumor].endToEndTimeHistogram));

                var rateHistograms = new List<KeyValuePair<string, ASETools.PreBucketedHistogram>>();
                rateHistograms.Add(new KeyValuePair<string, ASETools.PreBucketedHistogram>("reads per second reported", globalState[tumor].readsPerSecondReportedHistogram));
                rateHistograms.Add(new KeyValuePair<string, ASETools.PreBucketedHistogram>("reads per second in SNAP", globalState[tumor].readsPerSecondAllSnapHistogram));

                outputFile.WriteLine("Data for " + (tumor ? "tumor" : "normal"));
                outputFile.WriteLine("Slowest run " + globalState[tumor].timeForSlowestRun + " for " + globalState[tumor].caseIdForSlowestRun);
                timeHistograms.ForEach(_ => outputFile.WriteLine(_.Key + " mean " + _.Value.mean()));
                rateHistograms.ForEach(_ => outputFile.WriteLine(_.Key + " mean " + _.Value.mean()));

                outputFile.WriteLine();
                timeHistograms.ForEach(_ => outputFile.WriteLine(_.Key + " max " + _.Value.max()));
                rateHistograms.ForEach(_ => outputFile.WriteLine(_.Key + " max " + _.Value.max()));

                outputFile.WriteLine();
                ASETools.PreBucketedHistogram.WriteBatchOfHistogramCDFs(outputFile, timeHistograms);
                ASETools.PreBucketedHistogram.WriteBatchOfHistogramCDFs(outputFile, rateHistograms);

                outputFile.WriteLine();
            } // tumor/normal

            outputFile.WriteLine("**done**");
            outputFile.Close();

            Console.WriteLine("Processed " + casesToRun.Count() + " cases in " + ASETools.ElapsedTimeInSeconds(timer));

        } // Main
    } // Program
} // SummarizeSNAPPerformance
