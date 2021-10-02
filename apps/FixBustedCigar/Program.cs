using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using System.IO;
using ASELib;
using System.Diagnostics;

namespace FixBustedCigar
{
    class Program
    {
        static void Main(string[] args)
        {
            var timer = new Stopwatch();
            timer.Start();

            if (args.Count() != 2)
            {
                Console.WriteLine("usage: FixBustedCigar inputFile outputFile");
                Console.WriteLine("Fixes SAM lines that have a cigar string that soft clips the entire read by replacing them with an unmapped read.");
                return;
            }

            var inputFile = ASETools.CreateStreamReaderWithRetry(args[0]);
            if (null == inputFile)
            {
                return; // It already printed an error message
            }

            var outputFile = ASETools.CreateStreamWriterWithRetry(args[1]);
            if (null == outputFile)
            {
                return; // It already printed an error message
            }


            long nReads = 0;
            long nFixed = 0;

            string line = "";

            while (null != (line = inputFile.ReadLine()))
            {
                if (line.StartsWith("@"))
                {
                    outputFile.WriteLine(line);
                }

                nReads++;

                var SAMLine = new ASETools.SAMLine(line);
                if (!SAMLine.cigarAllClipping)
                {
                    outputFile.WriteLine(line);
                } else
                {
                    var fields = line.Split('\t');

                    var outputLine = fields[0];

                    for (int i = 1; i < fields.Count(); i++)
                    {
                        // We change flags to add in segment unmapped, MAPQ to 0 and cigar to "*", leaving the rest as is
                        if (i == 1)
                        {
                            outputLine += "\t" + (SAMLine.flag | ASETools.SAMLine.Unmapped);
                        } else if (i == 4)
                        {
                            outputLine += "\t0";    // MAPQ
                        } else if (i == 5)
                        {
                            outputLine += "\t*";    // CIGAR
                        } else
                        {
                            outputLine += "\t" + fields[i];
                        }
                    }

                    outputFile.WriteLine(outputLine);

                    nFixed++;
                    Console.WriteLine("Changed SAM line " + line);
                    Console.WriteLine("to: " + outputLine);
                }
            } // while we have an input line

            outputFile.Close();
            inputFile.Close();

            Console.WriteLine("Processed " + nReads + " reads of which we fixed " + nFixed + " in " + ASETools.ElapsedTimeInSeconds(timer));
        } // Main
    } // Program
} // namespace
