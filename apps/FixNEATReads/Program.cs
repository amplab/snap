using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using System.Diagnostics;
using ASELib;

// NEAT seems to produce reads with all the same read name.  Fix that.

//
// This is deterministic (just adds .# to the end of the read name), so it works fine running independently on each file of a paired-end FASTQ
//

namespace FixNEATReads
{
    internal class Program
    {
        static void Main(string[] args)
        {
            var timer = new Stopwatch();
            timer.Start();

            if (args.Length != 2)
            {
                Console.WriteLine("usage: FixNEATReads input.fq (or .fq.gz) output.fq");
                return;
            }

            var inputFilename = args[0];
            var outputFilename = args[1];

            var inputFile = ASETools.CreateStreamReaderWithRetryCompressedBasedOnFilename(inputFilename);
            if (null == inputFile)
            {
                Console.WriteLine("unable to open input file " + inputFilename);
                return;
            }

            var outputFile = ASETools.CreateStreamWriterWithRetry(outputFilename);
            if (null == outputFile)
            {
                Console.WriteLine("unable to open output file " + outputFilename);
                return;
            }

            long readsSeen = 0;

            string line;

            Console.Write("Progress (1 dot/10M reads): ");

            while (null != (line = inputFile.ReadLine()))
            {
                if (!line.StartsWith("@"))
                {
                    Console.WriteLine("input read name line doesn't start with @ (after " + readsSeen + " reads): " + line);
                    return;
                }

                var firstSpace = line.IndexOf(' ');
                if (firstSpace == -1)
                {
                    outputFile.WriteLine(line + "." + readsSeen);
                } else
                {
                    outputFile.WriteLine(line.Substring(0, firstSpace) + "." + readsSeen + line.Substring(firstSpace));
                }
                
                // Just copy the remaining three lines of the read

                for (int i = 0; i < 3; i++)
                {
                    line = inputFile.ReadLine();
                    if (null == line)
                    {
                        Console.WriteLine("EOF in the middle of a read after " + readsSeen + " reads");
                        return;
                    }

                    outputFile.WriteLine(line);
                }

                readsSeen++;
                if (readsSeen % 10000000 == 0)
                {
                    Console.Write(".");
                }
            } // while we read the first line of a read

            outputFile.Close();
            inputFile.Close();

            Console.WriteLine("Processed " + readsSeen + " reads in " + ASETools.ElapsedTime(timer));
        } // Main
    } // Program
} // namespace
