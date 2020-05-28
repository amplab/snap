using System;
using System.Collections.Generic;
using System.Linq;
using System.Runtime.Remoting.Messaging;
using System.Text;
using System.Threading.Tasks;
using System.Xml.Serialization;
using ASELib;

namespace SummarizeBowtiePerformance
{
    class Program
    {
        class PerThreadState
        {
            public const int maxTime = 7200;

            public ASETools.PreBucketedHistogram copyTimeHistogram = new ASETools.PreBucketedHistogram(0, maxTime, 60);
            public ASETools.PreBucketedHistogram bowtieTimeHistogram = new ASETools.PreBucketedHistogram(0, maxTime, 60);
            public ASETools.PreBucketedHistogram bowtieIndexTimeHistogram = new ASETools.PreBucketedHistogram(0, maxTime, 60);
            public ASETools.PreBucketedHistogram samToBamTimeHistogram = new ASETools.PreBucketedHistogram(0, maxTime, 60);
            public ASETools.PreBucketedHistogram sortTimeHistogram = new ASETools.PreBucketedHistogram(0, maxTime, 60);
            public ASETools.PreBucketedHistogram duplicateMarkHistogram = new ASETools.PreBucketedHistogram(0, maxTime, 60);
            public ASETools.PreBucketedHistogram indexHistogram = new ASETools.PreBucketedHistogram(0, maxTime, 60);
            public ASETools.PreBucketedHistogram copyOutTimeHistogram = new ASETools.PreBucketedHistogram(0, maxTime, 60);
            public ASETools.PreBucketedHistogram totalTimeHistogram = new ASETools.PreBucketedHistogram(0, maxTime, 60);
            public ASETools.PreBucketedHistogram copiedToFinalBamTimeHistogram = new ASETools.PreBucketedHistogram(0, maxTime, 60); // The time we charge them with--doesn't count copy in/out

            public void mergeWith(PerThreadState peer)
            {
                copyTimeHistogram.merge(peer.copyTimeHistogram);
                bowtieTimeHistogram.merge(peer.copyTimeHistogram);
                bowtieIndexTimeHistogram.merge(peer.bowtieIndexTimeHistogram);
                samToBamTimeHistogram.merge(peer.samToBamTimeHistogram);
                sortTimeHistogram.merge(peer.sortTimeHistogram);
                duplicateMarkHistogram.merge(peer.duplicateMarkHistogram);
                indexHistogram.merge(peer.indexHistogram);
                copyOutTimeHistogram.merge(peer.copyOutTimeHistogram);
                totalTimeHistogram.merge(peer.totalTimeHistogram);
                copiedToFinalBamTimeHistogram.merge(peer.copiedToFinalBamTimeHistogram);
            }
        } // PerThreadState

        static Dictionary<bool, PerThreadState> globalState = new Dictionary<bool, PerThreadState>(); // tumor -> state

        static void HandleOneCase(ASETools.Case case_, PerThreadState state)
        {

        }

        static void Main(string[] args)
        {
            var commonData = ASETools.CommonData.LoadCommonData(args);
            if  (commonData == null)
            {
                return;
            }

            if (commonData.listOfCases.Any(_ => _.bowtie_realigned_normal_dna_statictics_filename == "" || _.bowtie_realigned_tumor_dna_statictics_filename == "")) {
                Console.WriteLine("Some bowtie realignment statistics files are missing.  Try again once you've run all the realignments.");
                return;
            }

            string[] daysOfTheWeek = { "Mon", "Tue", "Wed", "Thu", "Fri", "Sat", "Sun" };   // This is the beginning of the date output.

            var outputFile = ASETools.CreateStreamWriterWithRetry(commonData.configuration.finalResultsDirectory + ASETools.BowtieHistogramsFilename);
            if (outputFile == null)
            {
                return;
            }

            int nPerDot;
            ASETools.PrintMessageAndNumberBar("Processing", "BAM files", commonData.listOfCases.Count() * 2, out nPerDot);

            int nProcessed = 0;

            foreach (var tumor in ASETools.BothBools) {

                ASETools.Case.ColumnGetter getInputFilename = c => tumor ? c.bowtie_realigned_tumor_dna_statictics_filename : c.bowtie_realigned_normal_dna_statictics_filename;

                var copyTimeHistogram = new ASETools.PreBucketedHistogram(0, 7200, 60);
                var bowtieTimeHistogram = new ASETools.PreBucketedHistogram(0, 7200, 60);
                var bowtieIndexTimeHistogram = new ASETools.PreBucketedHistogram(0, 7200, 60);
                var samToBamTimeHistogram = new ASETools.PreBucketedHistogram(0, 7200, 60);
                var sortTimeHistogram = new ASETools.PreBucketedHistogram(0, 7200, 60);
                var duplicateMarkHistogram = new ASETools.PreBucketedHistogram(0, 7200, 60);
                var indexHistogram = new ASETools.PreBucketedHistogram(0, 7200, 60);
                var copyOutTimeHistogram = new ASETools.PreBucketedHistogram(0, 7200, 60);
                var totalTimeHistogram = new ASETools.PreBucketedHistogram(0, 7200, 60);
                var copiedToFinalBamTimeHistogram = new ASETools.PreBucketedHistogram(0, 7200, 60); // The time we charge them with--doesn't count copy in/out

                var inputFilenames = commonData.listOfCases.Where(_ => getInputFilename(_) != "").Select(_ => getInputFilename(_)).ToList();

                foreach (var inputFilename in inputFilenames)
                {
                    var inputFile = ASETools.CreateStreamReaderWithRetry(inputFilename);
                    if (inputFile == null)
                    {
                        Console.WriteLine("Unable to open " + inputFilename);
                        return;
                    }

                    var inputLines = new List<string>();
                    string inputLine;
                    while (null != (inputLine = inputFile.ReadLine()))
                    {
                        inputLines.Add(inputLine);
                    }

                    var linesWithDates = Enumerable.Range(0, inputLines.Count()).Where(_ => daysOfTheWeek.Any(day => inputLines[_].StartsWith(day))).ToList();
                    if (linesWithDates.Count() != 7)
                    {
                        Console.WriteLine("Incorrect number of lines with dates for " + inputFilename);
                        return;
                    }

                    List<DateTime> dates;
                    try
                    {
                        dates = linesWithDates.Select(_ => ASETools.LinuxDateStringToDateTime(inputLines[_])).ToList();
                    } catch (Exception e)
                    {
                        Console.WriteLine("Exception.  Lines with dates: ");
                        linesWithDates.ForEach(_ => Console.WriteLine((inputLines[_])));
                        Console.WriteLine("input filename " + inputFilename);
                        throw e;
                    }

                    copyTimeHistogram.addValue(dates[1].Subtract(dates[0]).TotalSeconds);
                    bowtieTimeHistogram.addValue(dates[2].Subtract(dates[1]).TotalSeconds);
                    samToBamTimeHistogram.addValue(dates[3].Subtract(dates[2]).TotalSeconds);
                    sortTimeHistogram.addValue(dates[4].Subtract(dates[3]).TotalSeconds);
                    duplicateMarkHistogram.addValue(dates[5].Subtract(dates[4]).TotalSeconds);
                    indexHistogram.addValue(dates[6].Subtract(dates[5]).TotalSeconds);
                    totalTimeHistogram.addValue(dates[6].Subtract(dates[0]).TotalSeconds);
                    copiedToFinalBamTimeHistogram.addValue(dates[6].Subtract(dates[1]).TotalSeconds);

                    nProcessed++;
                    if (nProcessed % nPerDot == 0)
                    {
                        Console.Write(".");
                    }

                } // cases

                var histogramsToWrite = new List<KeyValuePair<string, ASETools.PreBucketedHistogram>>();
                histogramsToWrite.Add(new KeyValuePair<string, ASETools.PreBucketedHistogram>("Copy In", copyTimeHistogram));
                histogramsToWrite.Add(new KeyValuePair<string, ASETools.PreBucketedHistogram>("Bowtie", bowtieTimeHistogram));
                histogramsToWrite.Add(new KeyValuePair<string, ASETools.PreBucketedHistogram>("sam->bam", samToBamTimeHistogram));
                histogramsToWrite.Add(new KeyValuePair<string, ASETools.PreBucketedHistogram>("sort", sortTimeHistogram));
                histogramsToWrite.Add(new KeyValuePair<string, ASETools.PreBucketedHistogram>("Mark Duplicates", duplicateMarkHistogram));
                histogramsToWrite.Add(new KeyValuePair<string, ASETools.PreBucketedHistogram>("index", indexHistogram));
                histogramsToWrite.Add(new KeyValuePair<string, ASETools.PreBucketedHistogram>("total time", totalTimeHistogram));
                histogramsToWrite.Add(new KeyValuePair<string, ASETools.PreBucketedHistogram>("Bowtie->finished", copiedToFinalBamTimeHistogram));

                ASETools.PreBucketedHistogram.WriteBatchOfHistogramCDFs(outputFile, histogramsToWrite);

            } // tumor/normal

            outputFile.WriteLine("**done**");
            outputFile.Close();

            Console.WriteLine();
        } // Main
    }
}
