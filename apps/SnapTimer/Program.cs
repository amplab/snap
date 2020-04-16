using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using ASELib;
using System.IO;
using System.Diagnostics;

namespace SnapTimer
{
    class Program
    {
        static void Main(string[] args)
        {
            var commonData = ASETools.CommonData.LoadCommonData(args);
            if (commonData.configuration.commandLineArgs.Count() < 6)
            {
                Console.WriteLine("Usage: SnapTimer outputFile <snap args>");
                return;
            }

            int outputArgIndex = -1;

            var snapArgs = "";
            for (int i = 1; i < commonData.configuration.commandLineArgs.Count(); i++)
            {
                snapArgs += commonData.configuration.commandLineArgs[i] + " ";

                if (commonData.configuration.commandLineArgs[i] == "-o")
                {
                    outputArgIndex = i + 1;
                }
            }

            if (outputArgIndex == -1 || outputArgIndex >= commonData.configuration.commandLineArgs.Count())
            {
                Console.WriteLine("Couldn't properly process the output arg for SNAP.");
                return;
            }

            //
            // SNAP's output directory may not exist.  Create it if needed.
            //
            Directory.CreateDirectory(ASETools.GetDirectoryFromPathname(commonData.configuration.commandLineArgs[outputArgIndex]));

            var startInfo = new ProcessStartInfo(commonData.configuration.binariesDirectory + "snap.exe", snapArgs);
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

            if (snapOutput.Count() < 7)
            {
                Console.WriteLine("Too few lines of output from snap.");
                Console.WriteLine("Output: ");
                for (int i = 0; i < snapOutput.Count(); i++)
                {
                    Console.WriteLine(i + ": " + snapOutput[i]);
                }
                return;
            }

            Console.WriteLine("Run time: " + ASETools.ElapsedTimeInSeconds(runTimer));

            string loadingString = "Loading index from directory... ";
            if (!snapOutput[0].StartsWith(loadingString) || snapOutput[0].IndexOf('s') <= 0)
            {
                Console.WriteLine("Can't parse line 0");
                snapOutput.ForEach(_ => Console.WriteLine(_));
                return;
            }

            var loadingTimeString = snapOutput[0].Substring(loadingString.Count());
            loadingTimeString = loadingTimeString.Substring(0, loadingTimeString.IndexOf('s'));

            var loadingTime = Convert.ToInt32(loadingTimeString);

            if (snapOutput[2].IndexOf(',') == -1)
            {
                Console.WriteLine("Can't parse sort line");
                snapOutput.ForEach(_ => Console.WriteLine(_));
                return;
            }
            var sortTimeString = snapOutput[2].Substring(snapOutput[2].IndexOf(',') + 2);
            sortTimeString = sortTimeString.Substring(0, sortTimeString.IndexOf(' '));

            var sortTime = Convert.ToInt32(sortTimeString);

            var headerLine = snapOutput[5];
            if (!headerLine.Contains("Time in Aligner (s)") || !headerLine.Contains("Total Reads") || !headerLine.Contains("Reads/s") || 
                !headerLine.Contains("Aligned, MAPQ >= 10") || !headerLine.Contains("Aligned, MAPQ < 10") || !headerLine.Contains("Unaligned") ||
                !headerLine.Contains("Too Short/Too Many Ns"))
            {
                Console.WriteLine("Can't parse header line");
                snapOutput.ForEach(_ => Console.WriteLine(_));
                return;
            }

            var statsLine = snapOutput[6];
            int totalReads = ASETools.GetIntFromString(statsLine);
            int alignedHQ = ASETools.GetIntFromString(statsLine.Substring(headerLine.IndexOf("Aligned, MAPQ >= 10")));
            int alignedLQ = ASETools.GetIntFromString(statsLine.Substring(headerLine.IndexOf("Aligned, MAPQ < 10")));
            int unaligned = ASETools.GetIntFromString(statsLine.Substring(headerLine.IndexOf("Unaligned")));
            int tooShort = ASETools.GetIntFromString(statsLine.Substring(headerLine.IndexOf("Too Short/Too Many Ns")));
            int speed = ASETools.GetIntFromString(statsLine.Substring(headerLine.IndexOf("Reads/s")));
            int timeInAligner = ASETools.GetIntFromString(statsLine.Substring(headerLine.IndexOf("Time in Aligner (s)")));

            double pairs = 0;
            if (headerLine.Contains("%Pairs"))
            {
                pairs = ASETools.GetDoubleFromString(statsLine.Substring(headerLine.IndexOf("%Pairs")));
            }

            var outputFile = ASETools.CreateStreamWriterWithRetry(commonData.configuration.commandLineArgs[0]);
            if (null == outputFile)
            {
                return; // We already have an error message printed.
            }

            outputFile.WriteLine("Overall Runtime (s)\tLoading Time (s)\tTime in Aligner (s)\tSort Time (s)\tTotal Reads\tAligned, MAPQ >= 10\tAligned, MAPQ < 10\tUnaligned\tToo Short/Too Many Ns\t%Pairs\tReads/s\tStart Time\tStop Time");
            outputFile.WriteLine(ASETools.GetIntFromString(ASETools.ElapsedTimeInSeconds(runTimer)) + "\t" + loadingTime + "\t" + timeInAligner + "\t" + sortTime + "\t" + totalReads + "\t" + alignedHQ + "\t" + alignedLQ + "\t" + unaligned + "\t" +
                                 tooShort + "\t" + pairs + "\t" + speed + "\t" + startTime.ToString() + "\t" + stopTime.ToString());

            snapOutput.ForEach(_ => outputFile.WriteLine(_));
            outputFile.WriteLine("**done**");
            outputFile.Close();
        } // Main
    } // Program
} // namespace SnapTimer
