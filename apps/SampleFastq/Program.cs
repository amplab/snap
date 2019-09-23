using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using System.IO;
using ASELib;
using System.Diagnostics;

namespace SampleFastq
{
    class Program
    {
        static void Main(string[] args)
        {
            if (args.Count() != 3 && args.Count() != 5)
            {
                Console.WriteLine("usage: SampleFastq nReadsToWritePerInput input1 {input2} output1 {output2}");
                return;
            }

            Stopwatch timer = new Stopwatch();
            timer.Start();

            string input1Name = args[1];
            string input2Name = (args.Count() == 5) ? args[2] : null;
            string output1Name = (args.Count() == 5) ? args[3] : args[2];
            string output2Name = (args.Count() == 5) ? args[4] : null;
 
            var nReadsToWrite = Convert.ToInt64(args[0]);
            var nReadsWritten = 0;

            var input1 = ASETools.CreateStreamReaderWithRetryCompressedBasedOnFilename(input1Name);
            if (null == input1)
            {
                Console.WriteLine("Unable to open input file " + input1Name);
                return;
            }

            StreamReader input2;
            if (input2Name != null)
            {
                input2 = ASETools.CreateStreamReaderWithRetryCompressedBasedOnFilename(input2Name);
                if (input2 == null)
                {
                    Console.WriteLine("Unable to open file " + input2Name);
                    return;
                }
            }
            else
            {
                input2 = null;
            }

            var output1 = ASETools.CreateStreamWriterWithRetry(output1Name);
            if (null == output1)
            {
                Console.WriteLine("Unable to open output file " + output1Name);
                return;
            }

            StreamWriter output2;
            if (output2Name != null)
            {
                output2 = ASETools.CreateStreamWriterWithRetry(output2Name);
                if (null == output2)
                {
                    Console.WriteLine("Unable to open output file " + output2Name);
                    return;
                }
            } else
            {
                output2 = null;
            }

            var input1Size = new FileInfo(input1Name).Length;

            long input1SizeConsumed = 0;
            long input1LineGroupsConsumed = 0;

            const int nLinesPerFastqRead = 4;
            string[] lines1 = new string[nLinesPerFastqRead];
            string[] lines2 = new string[nLinesPerFastqRead];
            for (int i = 0; i < nLinesPerFastqRead; i++)
            {
                lines2[i] = null;
            }

            var random = new Random();

            int nPerDot;
            ASETools.PrintMessageAndNumberBar("Writing", "reads", nReadsToWrite, out nPerDot);



            while (true)
            {
                for (int i = 0; i < nLinesPerFastqRead; i++)
                {
                    lines1[i] = input1.ReadLine();
                    if (input2 != null)
                    {
                        lines2[i] = input2.ReadLine();
                    }
                }

                if (lines1[0] == null)
                {
                    if (lines1.Any(_ => _!= null) || lines2.Any(_ => _ != null))
                    {
                        Console.WriteLine("Input file wasn't a multiple of " + nLinesPerFastqRead + " lines.");
                        return;
                    }

                    break;
                }

                input1LineGroupsConsumed++;
                input1SizeConsumed += lines1.Select(_ => _.Length).Sum();

                var estimatedLinesRemaining = (input1Size - input1SizeConsumed) / ((double)input1SizeConsumed/input1LineGroupsConsumed);
                if (random.NextDouble()  < (nReadsToWrite - nReadsWritten) / estimatedLinesRemaining)
                {
                    lines1.ToList().ForEach(_ => output1.WriteLine(_));
                    if (output2 != null)
                    {
                        lines2.ToList().ForEach(_ => output2.WriteLine(_));
                    }
                    nReadsWritten++;

                    if (nReadsWritten % nPerDot == 0)
                    {
                        Console.Write(".");
                    }

                    if (nReadsWritten == nReadsToWrite)
                    {
                        break;
                    }
                }


            } // while (true)

            output1.Close();

            if (null != output2)
            {
                output2.Close();
            }

            Console.WriteLine();

            Console.WriteLine("Wrote " + nReadsWritten + " FASTQ reads in " + ASETools.ElapsedTimeInSeconds(timer));
        } // Main
    }
}
