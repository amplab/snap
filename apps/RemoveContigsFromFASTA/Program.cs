using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using ASELib;
using System.IO;
using System.Diagnostics;

namespace RemoveContigsFromFASTA
{
    class Program
    {
        static void Main(string[] args)
        {
            var timer = new Stopwatch();
            timer.Start();

            long totalLines = 0;
            long emittedLines = 0;
            long inputContigs = 0;
            long outputContigs = 0;

            if (args.Count() != 2)
            {
                Console.WriteLine("Usage: RemoveContigsFromFASTA inputFile outputFile");
                return;
            }

            var inputFile = ASETools.CreateStreamReaderWithRetry(args[0]);
            var outputFile = ASETools.CreateStreamWriterWithRetry(args[1]);

            if (inputFile == null || outputFile == null)
            {
                return; // The error message was already printed.
            }

            string currentContig = "";

            string inputLine;
            bool skipping = false;

            while (null != (inputLine = inputFile.ReadLine()))
            {
                totalLines++;
                if (inputLine.StartsWith(">"))
                {
                    inputContigs++;

                    currentContig = inputLine.Substring(1);
                    skipping = currentContig.Contains("_"); // Could generalize this if you cared

                    if (!skipping)
                    {
                        outputContigs++;
                    }
                }

                if (!skipping)
                {
                    emittedLines++;
                    outputFile.WriteLine(inputLine);
                }
            } // for each line

            inputFile.Close();
            outputFile.Close();

            Console.WriteLine("Processed " + ASETools.NumberWithCommas(totalLines) + " lines in " + inputContigs + " contigs, of which we kept " + ASETools.NumberWithCommas(emittedLines) + " in " + outputContigs + " contigs, in " + ASETools.ElapsedTimeInSeconds(timer));

        }
    }
}
