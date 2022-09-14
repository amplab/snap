using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using ASELib;
using System.IO;
using System.Diagnostics;
using System.Security.Cryptography;

namespace SampleFastqWithKnownSize
{
    class Program
    {
        static Random random = new Random();
        static long LongRandom(long nPossibleValues) // returns a long in the range 0..nPossibleValues - 1.
        {
            //
            // This could be done better to avoid a round-off bias, but I don't care that much for this application.  This ain't cryptography.
            //
            // Adapted from the first answer in https://stackoverflow.com/questions/677373/generate-random-values-in-c-sharp
            //

            var buffer = new byte[sizeof(Int64)];
            random.NextBytes(buffer);
            var rawValue =  BitConverter.ToInt64(buffer, 0);
            if (rawValue < 0)
            {
                rawValue = -rawValue;
            }

            return rawValue % nPossibleValues;
        }

        static void Main(string[] args)
        {
            var timer = new Stopwatch();
            timer.Start();

            if (args.Count() != 4 && args.Count() != 6)
            {
                Console.WriteLine("usage: SampleFastqWithKnownSize nWanted nInInput input {input2} output {output2}");
                return;
            }

            long nToEmit = Convert.ToInt64(args[0]);
            long nInTotal = Convert.ToInt64(args[1]);

            if (nToEmit > nInTotal)
            {
                Console.WriteLine("Can't emit more reads than are in the input.");
                return;
            }

            var nInputs = args.Count() == 4 ? 1 : 2;

            if (nInputs == 2 && (nToEmit % 2 != 0 || nInTotal % 2 != 0))
            {
                Console.WriteLine("You have paired-end reads and an odd number of total reads or reads to emit, which doesn't make sense.");
            }

            var inputs = new StreamReader[nInputs];
            var outputs = new StreamWriter[nInputs];

            inputs[0] = ASETools.CreateStreamReaderWithRetryCompressedBasedOnFilename(args[2]);
            outputs[0] = ASETools.CreateStreamWriterWithRetry((nInputs == 1) ? args[3] : args[4]);

            if (nInputs == 2)
            {
                inputs[1] = ASETools.CreateStreamReaderWithRetryCompressedBasedOnFilename(args[3]);
                outputs[1] = ASETools.CreateStreamWriterWithRetry(args[5]);
            }

            if (inputs.Any(_ => _ == null) || outputs.Any(_ => _ == null))
            {
                Console.WriteLine("Didn't open all files.");
                return;
            }

            //
            // These are expressed in pairs (if we're doing paired-end).
            //
            var nRemaining = nInTotal / nInputs;
            var nRemainingToEmit = nToEmit / nInputs;

            int nPerDot;
            ASETools.PrintMessageAndNumberBar("Processing", "reads", nInTotal, out nPerDot);

            long nProcessed = 0;

            while (nRemainingToEmit > 0)
            {
                bool emitThisOne = LongRandom(nRemaining) < nRemainingToEmit;

                for (int whichInput = 0; whichInput < nInputs; whichInput++)
                {
                    for (int whichLine = 0; whichLine < 4; whichLine++) // 4 lines/read in FASTQ
                    {
                        var line = inputs[whichInput].ReadLine();
                        if (line == null)
                        {
                            Console.WriteLine("Premature EOF.  Are you sure you gave the correct number of reads in the input?");
                            return;
                        }

                        if (emitThisOne)
                        {
                            outputs[whichInput].WriteLine(line);
                        }
                    } // lines/read
                } // input files

                if (emitThisOne)
                {
                    nRemainingToEmit--;
                }

                nRemaining--;
                nProcessed += nInputs;

                if (nProcessed % nPerDot == 0)  // We know nPerDot will always be even (unless it's 1, in which case we'll only print half the needed dots, but really, who wants to sample a FASTQ with only a handful of reads?)
                {
                    Console.Write(".");
                }
            }

            Console.WriteLine();

            inputs.ToList().ForEach(_ => _.Close());
            outputs.ToList().ForEach(_ => _.Close());

            Console.WriteLine("Total run time " + ASETools.ElapsedTimeInSeconds(timer));

        } // Main
    } // Program
} // Namespace
