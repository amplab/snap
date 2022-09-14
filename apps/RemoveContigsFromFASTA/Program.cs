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
            bool removeALTs = false;

            if (args.Count() != 2 && (args.Count() != 3 || args[2] != "-a"))
            {
                Console.WriteLine("Usage: RemoveContigsFromFASTA inputFile outputFile {-a}");
                Console.WriteLine("-a means to remove ALTs (contigs whose name starts with HLA- or contains _alt, regardless of case)");
                return;
            }

            removeALTs = args.Count() > 2 && args[2] == "-a";

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
                    if (removeALTs)
                    {
                        skipping = currentContig.ToLower().StartsWith("hla-") || currentContig.ToLower().Contains("_alt");  // Contains rather than EndsWith because there's extra junk on the contig line
                    }
                    else
                    {
                        skipping = currentContig.Contains("_"); // Could generalize this if you cared
                    }

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
