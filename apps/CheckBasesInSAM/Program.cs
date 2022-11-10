using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using ASELib;
using System.Diagnostics;

namespace CheckBasesInSAM
{
    internal class Program
    {
        static long readsProcessed = 0;

        static void validateBases(ASETools.SAMLine samLine)
        {
            bool printed = false;
            foreach (var base_ in samLine.seq) {
                if (printed)
                {
                    break;
                }

                switch (base_)
                {
                    case 'A':
                    case 'a':
                    case 'T':
                    case 't':
                    case 'C':
                    case 'c':
                    case 'G':
                    case 'g':
                    case 'N':
                    case 'n':
                    case 'U':
                    case 'u':
                    case 'R':
                    case 'r':
                    case 'Y':
                    case 'y':
                    case 'K':
                    case 'k':
                    case 'M':
                    case 'm':
                    case 'S':
                    case 's':
                    case 'W':
                    case 'w':
                    case 'B':
                    case 'b':
                    case 'D':
                    case 'd':
                    case 'H':
                    case 'h':
                    case 'V':
                    case 'v':
                    case '*':
                        break;
                    default:
                        Console.WriteLine("Read with strange base '" + base_ + "; in seq: " + samLine.line);
                        printed = true;
                        break;
                } // switch
            } // foreach base

            readsProcessed++;
            if (readsProcessed % 10000000 == 0)
            {
                Console.Write(".");
            }
        } // validateBases

        static void Main(string[] args)
        {
            var timer = new Stopwatch();
            timer.Start();

            if (args.Length != 1)
            {
                Console.WriteLine("usage: CheckBasesInSAM input.sam");
                return;
            }

            var inputFilename = args[0];
            var inputFile = ASETools.CreateStreamReaderWithRetry(inputFilename);
            if (null == inputFile)
            {
                Console.WriteLine("unable to open " + inputFilename);
                return;
            }

            Console.Write("Processing reads, one dot/10M: ");

            ASETools.SAMLine.ReadSAMLinesFromFile(inputFile, validateBases);

            Console.WriteLine();
            Console.WriteLine("Processed " + readsProcessed + " reads in " + ASETools.ElapsedTime(timer) + ", " + readsProcessed * 1000 / timer.ElapsedMilliseconds + " reads/s");



        } // Main
    } // Program
} // namespace
