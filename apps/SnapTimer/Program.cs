using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using ASELib;
using System.IO;
using System.Diagnostics;
using System.Net;

namespace SnapTimer
{
    class Program
    {
        static void Main(string[] args)
        {
            var configuration = ASETools.Configuration.loadFromFile(args);
            if (configuration.commandLineArgs.Count() < 7)
            {
                Console.WriteLine("Usage: SnapTimer outputFile tempDirectory <snap args>");
                return;
            }

            var tempDirectory = configuration.commandLineArgs[1];
            if (!tempDirectory.EndsWith(@"\"))
            {
                Console.WriteLine("Temp directory must end with a backslash.");
                return;
            }

            int outputArgIndex = -1;

            string initialOutputFile = "";
            string tempOutputFile = "";
            string tempOutputBAIFile = "";

            var initialInputFiles = new List<string>();
            var tempInputFiles = new List<string>();
            var tempFilesToDelete = new List<string>();

            var snapArgs = "";
            for (int i = 2; i < configuration.commandLineArgs.Count(); i++)
            {
                var thisArg = configuration.commandLineArgs[i];
                if (thisArg.EndsWith(".bam") || thisArg.EndsWith(".fastq"))
                {
                    if (outputArgIndex != i) 
                    {
                        initialInputFiles.Add(thisArg);
                        string tempFilename;
                        if (thisArg.StartsWith(@"\\"))
                        {
                            tempFilename = tempDirectory + ASETools.GetFileNameFromPathname(thisArg);
                            tempFilesToDelete.Add(tempFilename);
                        } else
                        {
                            tempFilename = thisArg;
                        }
                        snapArgs += tempFilename + " ";
                        tempInputFiles.Add(tempFilename);                        
                    }
                    else
                    {
                        if (thisArg.StartsWith(@"\\"))
                        {
                            tempOutputFile = tempDirectory + ASETools.GetFileNameFromPathname(thisArg);
                            tempOutputBAIFile = tempOutputFile + ".bai";

                            tempFilesToDelete.Add(tempOutputFile);
                            tempFilesToDelete.Add(tempOutputBAIFile);
                        } else
                        {
                            tempOutputFile = thisArg;
                            tempOutputBAIFile = thisArg + ".bai";
                        }
                        initialOutputFile = thisArg;
                        snapArgs += tempOutputFile + " ";
                    }
                }
                else
                {
                    snapArgs += thisArg + " ";
                }

                if (thisArg == "-o")
                {
                    outputArgIndex = i + 1;
                }
            }

            if (outputArgIndex == -1 || outputArgIndex >= configuration.commandLineArgs.Count())
            {
                Console.WriteLine("Couldn't properly process the output arg for SNAP.");
                return;
            }

            if (initialInputFiles.Count() == 0)
            {
                Console.WriteLine("Couldn't find snap input file.");
                return;
            }

            //
            // SNAP's output directory may not exist.  Create it if needed.
            //
            Directory.CreateDirectory(ASETools.GetDirectoryFromPathname(configuration.commandLineArgs[outputArgIndex]));

            //
            // Copy in the input files.
            //
            var copyInTimer = new Stopwatch();
            copyInTimer.Start();
            Console.Write("Copying input files...");
            foreach (var inputFile in initialInputFiles.Where(_ => _.StartsWith(@"\\")))
            {
                var tempFilename = tempDirectory + ASETools.GetFileNameFromPathname(inputFile);
                File.Copy(inputFile, tempFilename, true); // true allows overwrite
            }
            copyInTimer.Stop();
            Console.WriteLine(ASETools.ElapsedTimeInSeconds(copyInTimer));

            //
            // Run SNAP.
            //
            Console.WriteLine("Running snap " + snapArgs);
            var startInfo = new ProcessStartInfo(configuration.binariesDirectory + "snap.exe", snapArgs);
            startInfo.RedirectStandardOutput = true;
            startInfo.UseShellExecute = false;

            Stopwatch runTimer = new Stopwatch();
            runTimer.Start();

            var startTime = DateTime.Now;

            Process process;
            try
            {
                process = Process.Start(startInfo);
            }
            catch (Exception e)
            {
                Console.WriteLine("Error trying to start SNAP process ");
                Console.WriteLine("Exception message: " + e.Message);
                tempInputFiles.ForEach(_ => File.Delete(_));

                throw e;
            }

            var snapOutput = new List<string>();

            string outputLine;
            while (null != (outputLine = process.StandardOutput.ReadLine()))
            {
                snapOutput.Add(outputLine);
                Console.WriteLine(outputLine);
            }

            process.WaitForExit();
            runTimer.Stop();

            var stopTime = DateTime.Now;

            if (snapOutput.Count() < 5)
            {
                Console.WriteLine("Too few lines of output from snap.");
                Console.WriteLine("Output: ");
                for (int i = 0; i < snapOutput.Count(); i++)
                {
                    Console.WriteLine(i + ": " + snapOutput[i]);
                }
                tempInputFiles.ForEach(_ => File.Delete(_));
                File.Delete(tempOutputFile);
                File.Delete(tempOutputBAIFile);
                return;
            }

            Console.WriteLine("Run time: " + ASETools.ElapsedTimeInSeconds(runTimer));

            string loadingString = "Loading index from directory... ";
            if (!snapOutput[0].StartsWith(loadingString) || snapOutput[0].IndexOf('s') <= 0)
            {
                Console.WriteLine("Can't parse line 0");
                snapOutput.ForEach(_ => Console.WriteLine(_));
                tempInputFiles.ForEach(_ => File.Delete(_));
                File.Delete(tempOutputFile);
                File.Delete(tempOutputBAIFile);
                return;
            }

            var loadingTimeString = snapOutput[0].Substring(loadingString.Count());
            loadingTimeString = loadingTimeString.Substring(0, loadingTimeString.IndexOf('s'));

            var loadingTime = Convert.ToInt32(loadingTimeString);

            if (snapOutput[2].IndexOf(',') == -1)
            {
                Console.WriteLine("Can't parse sort line");
                snapOutput.ForEach(_ => Console.WriteLine(_));
                tempInputFiles.ForEach(_ => File.Delete(_));
                File.Delete(tempOutputFile);
                File.Delete(tempOutputBAIFile);
                return;
            }
            var sortTimeString = snapOutput[2].Substring(snapOutput[2].IndexOf(',') + 2);
            sortTimeString = sortTimeString.Substring(0, sortTimeString.IndexOf(' '));

            var sortTime = Convert.ToInt32(sortTimeString);

            var headerLine = snapOutput[3];
            if (!headerLine.Contains("Time in Aligner (s)") || !headerLine.Contains("Total Reads") || !headerLine.Contains("Reads/s") || 
                !headerLine.Contains("Aligned, MAPQ >= 10") || !headerLine.Contains("Aligned, MAPQ < 10") || !headerLine.Contains("Unaligned") ||
                !headerLine.Contains("Too Short/Too Many Ns"))
            {
                Console.WriteLine("Can't parse header line");
                snapOutput.ForEach(_ => Console.WriteLine(_));
                tempInputFiles.ForEach(_ => File.Delete(_));
                File.Delete(tempOutputFile);
                File.Delete(tempOutputBAIFile);
                return;
            }

            var statsLine = snapOutput[4];
            long totalReads = ASETools.GetLongFromString(statsLine);
            long alignedHQ = ASETools.GetLongFromString(statsLine.Substring(headerLine.IndexOf("Aligned, MAPQ >= 10")));
            long alignedLQ = ASETools.GetLongFromString(statsLine.Substring(headerLine.IndexOf("Aligned, MAPQ < 10")));
            long unaligned = ASETools.GetLongFromString(statsLine.Substring(headerLine.IndexOf("Unaligned")));
            long tooShort = ASETools.GetLongFromString(statsLine.Substring(headerLine.IndexOf("Too Short/Too Many Ns")));
            long speed = ASETools.GetLongFromString(statsLine.Substring(headerLine.IndexOf("Reads/s")));
            long timeInAligner = ASETools.GetLongFromString(statsLine.Substring(headerLine.IndexOf("Time in Aligner (s)")));
            double agCalledSingle = 0;
            double agUsedSingle = 0;

            if (headerLine.IndexOf("%AgSingle") >= 0)
            {
                agCalledSingle = ASETools.GetDoubleFromString(statsLine.Substring(headerLine.IndexOf("%AgSingle"))) / 100.0;
                agUsedSingle = ASETools.GetDoubleFromString(statsLine.Substring(headerLine.IndexOf("%AgUsedSingle"))) / 100.0;
            }

            double pairs = 0;
            if (headerLine.Contains("%Pairs"))
            {
                pairs = ASETools.GetDoubleFromString(statsLine.Substring(headerLine.IndexOf("%Pairs")));
            }

            Console.Write("Copying back output...");
            var copyOutTimer = new Stopwatch();
            copyOutTimer.Start();
            if (tempOutputFile != initialOutputFile)
            {
                File.Copy(tempOutputFile, initialOutputFile, true);
                File.Copy(tempOutputBAIFile, initialOutputFile + ".bai", true);
            }

            copyOutTimer.Stop();
            Console.Write(ASETools.ElapsedTimeInSeconds(copyOutTimer));

            tempFilesToDelete.ForEach(_ => File.Delete(_));

            var outputFile = ASETools.CreateStreamWriterWithRetry(configuration.commandLineArgs[0]);
            if (null == outputFile)
            {
                return; // We already have an error message printed.
            }

            outputFile.WriteLine("Copy in Time (s)\tOverall Runtime (s)\tLoading Time (s)\tTime in Aligner (s)\tSort Time (s)\tTotal Reads\tAligned, MAPQ >= 10\tAligned, MAPQ < 10\tUnaligned\tToo Short/Too Many Ns\t%Pairs\tReads/s\tCopy out Time (s)\tStart Time\tStop Time\tAG called Single\tAG used single result");
            outputFile.WriteLine(
                ASETools.ElapsedTimeInSeconds(copyInTimer) + "\t" +
                ASETools.GetIntFromString(ASETools.ElapsedTimeInSeconds(runTimer)) + "\t" + 
                loadingTime + "\t" + 
                timeInAligner + "\t" + 
                sortTime + "\t" + 
                totalReads + "\t" + 
                alignedHQ + "\t" + 
                alignedLQ + "\t" + 
                unaligned + "\t" +
                tooShort + "\t" + 
                pairs + "\t" + 
                speed + "\t" + 
                ASETools.ElapsedTimeInSeconds(copyOutTimer) + "\t" +
                startTime.ToString() + "\t" + 
                stopTime.ToString() + "\t" +
                agCalledSingle + "\t" +
                agUsedSingle);

            snapOutput.ForEach(_ => outputFile.WriteLine(_));
            outputFile.WriteLine("**done**");
            outputFile.Close();
        } // Main
    } // Program
} // namespace SnapTimer
