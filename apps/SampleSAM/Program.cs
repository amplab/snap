using System;
using System.CodeDom;
using System.Collections.Generic;
using System.Diagnostics;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using ASELib;

namespace SampleSAM
{
    class Program
    {
        static void Main(string[] args)
        {
            var timer = new Stopwatch();
            timer.Start();

            if (args.Count() != 3)
            {
                Console.WriteLine("usage: SampleSAM inputFile outputFile percentage (expressed as an int)");
                return;
            }

            var inputFile = ASETools.CreateStreamReaderWithRetry(args[0]);
            if (inputFile == null)
            {
                Console.WriteLine("Unable to open input file " + args[0]);
                return;
            }

            int percentage = Convert.ToInt32(args[2]);
            if (percentage < 1 || percentage > 99)
            {
                Console.WriteLine("Percentage must be between 1 and 99");
                return;
            }

            var outputFile = ASETools.CreateStreamWriterWithRetry(args[1]);
            if (outputFile == null)
            {
                Console.WriteLine("Unable to open output file " + args[1]);
                return;
            }

            var halfSeenReads = new Dictionary<string, string>();  // Maps name -> read

            int nTotalReads = 0;
            int nSingleEndReadsWritten = 0;
            int nPairedEndReadsWritten = 0;
            var random = new Random();

            string line;
            while (null != (line = inputFile.ReadLine()))
            {
                if (line == "")
                {
                    Console.WriteLine("blank input line");
                    return;
                }

                if (line[0] == '@')
                {
                    // Header line
                    outputFile.WriteLine(line);
                    continue;
                }

                nTotalReads++;

                var fields = line.Split('\t');
                if (fields.Count() < 11)
                {
                    Console.WriteLine("Too short SAM line: " + line);
                    return;
                }

                var flags = Convert.ToInt32(fields[1]);
                if ((flags & 0x900) != 0)
                {
                    // Secondary and/or supplentary alignment.  Skip it.
                    continue;
                }

                if ((flags & 0x1) == 0)
                {
                    //
                    // A single-ended read
                    //
                    if (random.Next(100) < percentage)
                    {
                        outputFile.WriteLine(line);
                        nSingleEndReadsWritten++;
                    }

                    continue;
                } // single end

                if (halfSeenReads.ContainsKey(fields[0]))
                {
                    if (random.Next(99) < percentage)
                    {
                        outputFile.WriteLine(halfSeenReads[fields[0]]);
                        outputFile.WriteLine(line);
                        nPairedEndReadsWritten += 2;
                    }

                    halfSeenReads.Remove(fields[0]);
                } else
                {
                    halfSeenReads.Add(fields[0], line);
                }
            } // while we have input lines

            outputFile.Close();
            inputFile.Close();

            Console.WriteLine(nTotalReads + " total reads, " + nPairedEndReadsWritten + " paired-end reads written, " + nSingleEndReadsWritten + " single end reads written and " + halfSeenReads.Count() + " with unfound mates in " +
                ASETools.ElapsedTimeInSeconds(timer));
        } // Main
    } // Program
} // namespace
