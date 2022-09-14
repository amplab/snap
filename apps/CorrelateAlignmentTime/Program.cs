using System;
using System.Collections.Generic;
using System.IO;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using ASELib;
using System.Diagnostics;

// Compute a scatter graph of alignment times for the same reads aligned by different aligners (or different versions)

namespace CorrelateAlignmentTime
{
    internal class Program
    {
        static int cheezyLog2(int value)
        {
            if (value < 0)
            {
                throw new ArgumentException("cheezyLog2: value may not be negavive: " + value);
            }

            int retVal = 0;
            value = value >> 1; // Because log2(1) == 0 (this also makes 0 go into 0)

            while (value > 0)
            {
                value = value >> 1;
                retVal++;
            }

            return retVal;
        } // cheezyLog2


        static void Main(string[] args)
        {
            var timer = new Stopwatch();
            timer.Start();

            if (args.Length != 3 && args.Length != 4 && args.Length != 5)
            {
                Console.WriteLine("usage: CorrelateAlignmentTime input1.SAM input2.SAM output.txt {pairs_output.txt {downsample factor}}");
                return;
            }

            const int nInputFiles = 2;  // Can't really change this without totally redoing the code, it's just here for clarity

            var inputFile = new StreamReader[nInputFiles];
            for (int i = 0; i < nInputFiles; i++)
            {
                inputFile[i] = ASETools.CreateStreamReaderWithRetry(args[i]);
                if (inputFile[i] == null)
                {
                    Console.WriteLine("Unable to open " + args[i]);
                    return;
                }
            } // for each input file.

            var outputFile = ASETools.CreateStreamWriterWithRetry(args[2]);
            if (outputFile == null)
            {
                return;
            }

            StreamWriter pairsOutputFile = null;
            if (args.Length > 3)
            {
                pairsOutputFile = ASETools.CreateStreamWriterWithRetry(args[3]);
                if (pairsOutputFile == null)
                {
                    Console.WriteLine("Unable to open " + args[3]);
                    return;
                }
            }

            int downsample = 1;
            if (args.Length > 4)
            {
                downsample = Convert.ToInt32(args[4]);
                if (downsample < 1)
                {
                    Console.WriteLine("downsample must be >= 1 (and it's meaningless if it's 1)");
                    return;
                }
            }

            var rand = new Random();

            //
            // We use decoratedQname() to add /1 or /2 to paired-end reads to make sure they match correctly.
            //

            const int maxTimeInUs = 1024 * 1024;    // Anything over this we'll just treat as 1s
            int nBuckets = cheezyLog2(maxTimeInUs) + 1;
            long[,] counts = new long[nBuckets, nBuckets];
            var unmatchedReads = new Dictionary<string, ASETools.SAMLine>[nInputFiles];

            bool[] notEOF = new bool[nInputFiles];
            for (int i = 0; i < nInputFiles; i++)
            {
                notEOF[i] = true;
                unmatchedReads[i] = new Dictionary<string, ASETools.SAMLine>();
            }

            long processedPairs = 0;
            var pearsonR = new ASETools.PearsonR();
            var pearsonROfLogTime = new ASETools.PearsonR();

            Console.Write("Processing input lines (one dot/million): ");
            long nLinesRead = 0;

            while (notEOF.Any(_ => _))
            {
                for (int i = 0; i < nInputFiles; i++)
                {
                    if (notEOF[i])
                    {
                        var line = inputFile[i].ReadLine();
                        if (line == null)
                        {
                            notEOF[i] = false;
                            continue;
                        }

                        nLinesRead++;
                        if (nLinesRead % 1000000 == 0)
                        {
                            Console.Write(".");
                        }

                        if (!line.StartsWith("@"))
                        {
                            var samLine = new ASETools.SAMLine(line);

                            if (samLine.isSecondaryAlignment() || samLine.isSupplementaryAlignment())
                            {
                                continue;
                            }

                            if (unmatchedReads[1 - i].ContainsKey(samLine.decoratedQname()))
                            {
                                var samLines = new ASETools.SAMLine[2];
                                samLines[i] = samLine;
                                samLines[1-i] = unmatchedReads[1 - i][samLine.decoratedQname()];

                                var times = new int[nInputFiles];

                                for (int j = 0; j < nInputFiles; j++)
                                {
                                    times[j] = samLines[j].AT();
                                } // input files

                                counts[Math.Min(cheezyLog2(times[0]), nBuckets - 1), Math.Min(cheezyLog2(times[1]), nBuckets - 1)]++;
                                pearsonR.addSamplePair(times[0], times[1]);
                                pearsonROfLogTime.addSamplePair(Math.Min(cheezyLog2(times[0]), nBuckets - 1), Math.Min(cheezyLog2(times[1]), nBuckets - 1));

                                unmatchedReads[1 - i].Remove(samLine.decoratedQname());

                                if (pairsOutputFile != null && rand.Next(downsample) == 0)
                                {
                                    pairsOutputFile.WriteLine(times[0] + ", " + times[1]);
                                }

                                processedPairs++;
                            }  else
                            {
                                unmatchedReads[i].Add(samLine.decoratedQname(), samLine);
                            }
                        } // if it's not a header line
                    } // if not EOF on this input file
                } // for each input file
            } // while we have data

            Console.WriteLine();
            Console.WriteLine("Processed " + nLinesRead + " lines in " +ASETools.ElapsedTimeInSeconds(timer));

            outputFile.WriteLine("Correctly matched " + processedPairs + " reads from each input");
            Console.WriteLine("Correctly matched " + processedPairs + " reads from each input");

            outputFile.WriteLine("Pearson's correlation coefficient " + pearsonR.r() +".  Pearson's correlation coefficient of log times is " + pearsonROfLogTime.r() + ".");
            Console.WriteLine("Pearson's correlation coefficient " + pearsonR.r() + ".  Pearson's correlation coefficient of log times is " + pearsonROfLogTime.r() + ".");

            if (unmatchedReads[0].Count() > 0 || unmatchedReads[1].Count() > 0)
            {
                outputFile.WriteLine("The first input file has " + unmatchedReads[0].Count() + " unmatched reads and the second has " + unmatchedReads[1].Count() + ".");
                Console.WriteLine("The first input file has " + unmatchedReads[0].Count() + " unmatched reads and the second has " + unmatchedReads[1].Count() + ".");
            }

            outputFile.Write("File1/File2");
            for (int i = 0; i < nBuckets; i++)
            {
                outputFile.Write("\t" +(1 << i) + "us");
            }
            outputFile.WriteLine();

            for (int i = 0; i < nBuckets; i++)
            {
                outputFile.Write((1 << i) + "us");
                for (int j = 0; j < nBuckets; j++)
                {
                    outputFile.Write("\t" + counts[i, j]);
                }
                outputFile.WriteLine();
            }

            outputFile.Close();
            if (pairsOutputFile != null)
            {
                pairsOutputFile.Close();
            }
        } // Main
    } // Program
} // namespace
