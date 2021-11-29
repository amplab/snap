using System;
using System.Collections.Generic;
using System.Linq;
using System.Runtime.Remoting.Messaging;
using System.Text;
using System.Threading.Tasks;
using System.Xml.Serialization;
using ASELib;
using System.Diagnostics;
using System.Data;
using System.Runtime.Remoting;

namespace SummarizeBowtiePerformance // This is a misnomer: it all Linux aligners.  They have identical formats.
{
    class Program
    {
        static Dictionary<string, ASETools.CasePairedness> casePairedness = null;
        class HistogramState
        { 
            public ASETools.PreBucketedHistogram copyTimeHistogram = new ASETools.PreBucketedHistogram(0, ASETools.maxTimeForRealignHistograms, ASETools.timeStepForRealignHistograms);
            public ASETools.PreBucketedHistogram alignerTimeHistogram = new ASETools.PreBucketedHistogram(0, ASETools.maxTimeForRealignHistograms, ASETools.timeStepForRealignHistograms);
            public ASETools.PreBucketedHistogram bowtieIndexTimeHistogram = new ASETools.PreBucketedHistogram(0, ASETools.maxTimeForRealignHistograms, ASETools.timeStepForRealignHistograms);
            public ASETools.PreBucketedHistogram samToBamTimeHistogram = new ASETools.PreBucketedHistogram(0, ASETools.maxTimeForRealignHistograms, ASETools.timeStepForRealignHistograms);
            public ASETools.PreBucketedHistogram sortTimeHistogram = new ASETools.PreBucketedHistogram(0, ASETools.maxTimeForRealignHistograms, ASETools.timeStepForRealignHistograms);
            public ASETools.PreBucketedHistogram duplicateMarkHistogram = new ASETools.PreBucketedHistogram(0, ASETools.maxTimeForRealignHistograms, ASETools.timeStepForRealignHistograms);
            public ASETools.PreBucketedHistogram indexHistogram = new ASETools.PreBucketedHistogram(0, ASETools.maxTimeForRealignHistograms, ASETools.timeStepForRealignHistograms);
            public ASETools.PreBucketedHistogram copyOutTimeHistogram = new ASETools.PreBucketedHistogram(0, ASETools.maxTimeForRealignHistograms, ASETools.timeStepForRealignHistograms);
            public ASETools.PreBucketedHistogram totalTimeHistogram = new ASETools.PreBucketedHistogram(0, ASETools.maxTimeForRealignHistograms, ASETools.timeStepForRealignHistograms);
            public ASETools.PreBucketedHistogram copiedToFinalBamTimeHistogram = new ASETools.PreBucketedHistogram(0, ASETools.maxTimeForRealignHistograms, ASETools.timeStepForRealignHistograms); // The time we charge them with--doesn't count copy in/out

            public ASETools.PreBucketedHistogram readsPerSecondReportedHistogram = new ASETools.PreBucketedHistogram(0, ASETools.maxSpeedForRealignmentHistograms, ASETools.stepForSpeedRealignmentHistograms);
            public ASETools.PreBucketedHistogram readsPerSecondAllTimeHistogram = new ASETools.PreBucketedHistogram(0, ASETools.maxSpeedForRealignmentHistograms, ASETools.stepForSpeedRealignmentHistograms);

            public void mergeWith(HistogramState peer)
            {
                copyTimeHistogram.merge(peer.copyTimeHistogram);
                alignerTimeHistogram.merge(peer.alignerTimeHistogram);
                bowtieIndexTimeHistogram.merge(peer.bowtieIndexTimeHistogram);
                samToBamTimeHistogram.merge(peer.samToBamTimeHistogram);
                sortTimeHistogram.merge(peer.sortTimeHistogram);
                duplicateMarkHistogram.merge(peer.duplicateMarkHistogram);
                indexHistogram.merge(peer.indexHistogram);
                copyOutTimeHistogram.merge(peer.copyOutTimeHistogram);
                totalTimeHistogram.merge(peer.totalTimeHistogram);
                copiedToFinalBamTimeHistogram.merge(peer.copiedToFinalBamTimeHistogram);
                readsPerSecondReportedHistogram.merge(peer.readsPerSecondReportedHistogram);
                readsPerSecondAllTimeHistogram.merge(peer.readsPerSecondAllTimeHistogram);
            }
        } // HistogramState

        class PerThreadState
        {
            public Dictionary<ASETools.EndCount, HistogramState> histograms = new Dictionary<ASETools.EndCount, HistogramState>();

            public PerThreadState()
            {
                foreach (var endCount in ASETools.EnumUtil.GetValues<ASETools.EndCount>())
                {
                    histograms.Add(endCount, new HistogramState());
                }
            }

            public void mergeWith(PerThreadState peer)
            {
                foreach (var endCount in ASETools.EnumUtil.GetValues<ASETools.EndCount>())
                {
                    histograms[endCount].mergeWith(peer.histograms[endCount]);
                }
            } // mergeWith
        } // PerThreadState

        static Dictionary<bool, PerThreadState> globalState = new Dictionary<bool, PerThreadState>(); // tumor -> state
        static ASETools.Aligner aligner;


        static void HandleOneCase(ASETools.Case case_, Dictionary<bool, PerThreadState> state)
        {
            var readStats = ASETools.ReadStatistics.ReadAllFromFile(case_.read_statictics_filename);

            foreach (var tumor in ASETools.BothBools)
            {
                if (!state.ContainsKey(tumor))
                {
                    state.Add(tumor, new PerThreadState());
                }

                var inputFilename = case_.realignments[aligner][tumor].dna_statistics_filename;
                if (inputFilename == "")
                {
                    continue;
                }

                var inputFile = ASETools.CreateStreamReaderWithRetry(inputFilename);
                if (inputFile == null)
                {
                    throw new Exception("Unable to open " + inputFilename);
                }

                var inputLines = new List<string>();
                string inputLine;
                while (null != (inputLine = inputFile.ReadLine()))
                {
                    inputLines.Add(inputLine);
                }

                var linesWithDates = Enumerable.Range(0, inputLines.Count()).Where(_ => ASETools.daysOfTheWeek.Any(day => inputLines[_].StartsWith(day))).ToList();
                if (linesWithDates.Count() != 7)
                {
                    Console.WriteLine("Incorrect number of lines with dates for " + inputFilename);
                    return;
                }

                List<DateTime> dates;
                try
                {
                    dates = linesWithDates.Select(_ => ASETools.LinuxDateStringToDateTime(inputLines[_])).ToList();
                }
                catch (Exception e)
                {
                    Console.WriteLine("Exception.  Lines with dates: ");
                    linesWithDates.ForEach(_ => Console.WriteLine((inputLines[_])));
                    Console.WriteLine("input filename " + inputFilename);
                    throw e;
                }

                if (dates[2].Subtract(dates[1]).TotalSeconds < 30)
                {
                    Console.WriteLine("Extremely fast aligner run for case " + case_.case_id + ", stats file " + inputFilename);
                }

                bool paired = casePairedness[case_.case_id].pairedness[tumor][true];
                ASETools.EndCount endCount;
                if (paired)
                {
                    endCount = ASETools.EndCount.Paired;
                } else
                {
                    endCount = ASETools.EndCount.Single;
                }

                state[tumor].histograms[endCount].copyTimeHistogram.addValue(dates[1].Subtract(dates[0]).TotalSeconds);
                state[tumor].histograms[endCount].alignerTimeHistogram.addValue(dates[2].Subtract(dates[1]).TotalSeconds);
                state[tumor].histograms[endCount].samToBamTimeHistogram.addValue(dates[3].Subtract(dates[2]).TotalSeconds);
                state[tumor].histograms[endCount].sortTimeHistogram.addValue(dates[4].Subtract(dates[3]).TotalSeconds);
                state[tumor].histograms[endCount].duplicateMarkHistogram.addValue(dates[5].Subtract(dates[4]).TotalSeconds);
                state[tumor].histograms[endCount].indexHistogram.addValue(dates[6].Subtract(dates[5]).TotalSeconds);
                state[tumor].histograms[endCount].totalTimeHistogram.addValue(dates[6].Subtract(dates[0]).TotalSeconds);
                state[tumor].histograms[endCount].copiedToFinalBamTimeHistogram.addValue(dates[6].Subtract(dates[1]).TotalSeconds);

                state[tumor].histograms[endCount].readsPerSecondReportedHistogram.addValue(readStats[tumor][true].totalReads / dates[2].Subtract(dates[1]).TotalSeconds);    // true means dna
                state[tumor].histograms[endCount].readsPerSecondAllTimeHistogram.addValue(readStats[tumor][true].totalReads / dates[6].Subtract(dates[1]).TotalSeconds);// true means dna
            } // tumor or normal
        } // HandleOneCase

        static void FinishUp(Dictionary<bool, PerThreadState> state)
        {
            lock (globalState)
            {
                foreach (var tumor in ASETools.BothBools)
                {
                    if (state.ContainsKey(tumor) && state[tumor] != null && globalState.ContainsKey(tumor) && globalState[tumor] != null)
                    {
                        globalState[tumor].mergeWith(state[tumor]);
                    }
                }
            }
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

            if (configuration.commandLineArgs.Count() != 1 || !ASETools.isStringAnAlignerName(configuration.commandLineArgs[0]))
            { 
                Console.WriteLine("usage: SummarizeBowtiePerformance linuxAligner");
                return;
            }

            casePairedness = ASETools.CasePairedness.ReadConsolodatedCasePairedness(configuration.finalResultsDirectory + ASETools.ConsolodatedCasePairednessFilename);
            if (null == casePairedness)
            {
                Console.WriteLine("Unable to read case pairedness from " + configuration.finalResultsDirectory + ASETools.ConsolodatedCasePairednessFilename);
                return;
            }

            aligner = ASETools.AlignerNameToAligner(configuration.commandLineArgs[0]);

            if (aligner == ASETools.Aligner.SNAP)
            {
                Console.WriteLine("Sorry, SNAP is not a Linux Aligner (at least as we run it).");
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

            if (listOfCases.Any(_ => _.read_statictics_filename == "")) {
                Console.WriteLine("Some read statistics files are missing.  Try again once you've run all the realignments.");
                //return;
                listOfCases = listOfCases.Where(_ => _.read_statictics_filename != "").ToList();
            }

            int nPerDot;
            ASETools.PrintMessageAndNumberBar("Processing", "cases", listOfCases.Count(), out nPerDot);
            var threading = new ASETools.WorkerThreadHelper<ASETools.Case, Dictionary<bool, PerThreadState>>(listOfCases, HandleOneCase, FinishUp, null, nPerDot);
            threading.run();

            var outputFile = ASETools.CreateStreamWriterWithRetry(configuration.finalResultsDirectory + ASETools.alignerName[aligner] + "Histograms.txt");
            if (outputFile == null)
            {
                return;
            }

            foreach (var tumor in ASETools.BothBools) {
                foreach (var ends in ASETools.EnumUtil.GetValues<ASETools.EndCount>())
                {
                    if (globalState[tumor].histograms[ends].copyTimeHistogram.count() == 0)
                    {
                        //
                        // No data, just skip it.
                        //
                        continue;
                    }


                    if (!tumor)
                    {
                        outputFile.WriteLine(); // Space after the tumor results
                    }

                    outputFile.WriteLine(ASETools.tumorToString[tumor] + " " + ASETools.endCountToString[ends] + " statistics");

                    var histogramsToWrite = new List<KeyValuePair<string, ASETools.PreBucketedHistogram>>();
                    histogramsToWrite.Add(new KeyValuePair<string, ASETools.PreBucketedHistogram>("Copy In", globalState[tumor].histograms[ends].copyTimeHistogram));
                    histogramsToWrite.Add(new KeyValuePair<string, ASETools.PreBucketedHistogram>(ASETools.alignerName[aligner], globalState[tumor].histograms[ends].alignerTimeHistogram));
                    histogramsToWrite.Add(new KeyValuePair<string, ASETools.PreBucketedHistogram>("sam->bam", globalState[tumor].histograms[ends].samToBamTimeHistogram));
                    histogramsToWrite.Add(new KeyValuePair<string, ASETools.PreBucketedHistogram>("sort", globalState[tumor].histograms[ends].sortTimeHistogram));
                    histogramsToWrite.Add(new KeyValuePair<string, ASETools.PreBucketedHistogram>("Mark Duplicates", globalState[tumor].histograms[ends].duplicateMarkHistogram));
                    histogramsToWrite.Add(new KeyValuePair<string, ASETools.PreBucketedHistogram>("index", globalState[tumor].histograms[ends].indexHistogram));
                    histogramsToWrite.Add(new KeyValuePair<string, ASETools.PreBucketedHistogram>("total time", globalState[tumor].histograms[ends].totalTimeHistogram));
                    histogramsToWrite.Add(new KeyValuePair<string, ASETools.PreBucketedHistogram>(ASETools.alignerName[aligner] + "->finished", globalState[tumor].histograms[ends].copiedToFinalBamTimeHistogram));

                    histogramsToWrite.ForEach(_ => outputFile.WriteLine(_.Key + " mean " + _.Value.mean()));
                    outputFile.WriteLine();

                    histogramsToWrite.ForEach(_ => outputFile.WriteLine(_.Key + " max " + _.Value.max()));
                    outputFile.WriteLine();

                    ASETools.PreBucketedHistogram.WriteBatchOfHistogramCDFs(outputFile, histogramsToWrite);

                    histogramsToWrite = new List<KeyValuePair<string, ASETools.PreBucketedHistogram>>();
                    histogramsToWrite.Add(new KeyValuePair<string, ASETools.PreBucketedHistogram>("Reads per second in " + ASETools.alignerName[aligner], globalState[tumor].histograms[ends].readsPerSecondReportedHistogram));
                    histogramsToWrite.Add(new KeyValuePair<string, ASETools.PreBucketedHistogram>("Reads per second to final BAM", globalState[tumor].histograms[ends].readsPerSecondAllTimeHistogram));

                    outputFile.WriteLine();

                    histogramsToWrite.ForEach(_ => outputFile.WriteLine(_.Key + " mean " + _.Value.mean()));
                    outputFile.WriteLine();

                    histogramsToWrite.ForEach(_ => outputFile.WriteLine(_.Key + " max " + _.Value.max()));
                    outputFile.WriteLine();
                    ASETools.PreBucketedHistogram.WriteBatchOfHistogramCDFs(outputFile, histogramsToWrite);

                } // single/paired ended
            } // tumor/normal

            outputFile.WriteLine("**done**");
            outputFile.Close();

            Console.WriteLine();
            Console.WriteLine("Processed " + listOfCases.Count() + " cases in " + ASETools.ElapsedTimeInSeconds(timer) + ", " + listOfCases.Count() / timer.Elapsed.TotalSeconds + " cases/s");
            
        } // Main
    }
}
