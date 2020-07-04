using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using ASELib;
using System.Diagnostics;
using System.Threading;

namespace FixFASTA
{
    class Program
    {
        static void Main(string[] args)
        {
            var timer = new Stopwatch();
            timer.Start();

            if (args.Count() != 2)
            {
                Console.WriteLine("usage: FixFASTA inputFile outputFile");
                return;
            }

            var inputFile = ASETools.CreateStreamReaderWithRetry(args[0]);
            if (inputFile == null)
            {
                Console.WriteLine("Unable to open " + args[0] + " for read.");
                return;
            }

            var outputFile = ASETools.CreateStreamWriterWithRetry(args[1]);
            if (outputFile == null)
            {
                Console.WriteLine("Unable to open " + args[0] + " for write.");
                return;
            }

            string inputLine;
            int nContigs = 0;
            int nLines = 0;
            long nBases = 0;
            var baseCounts = new Dictionary<char, long>();

            while (null != (inputLine = inputFile.ReadLine()))
            {
                nLines++;

                if (inputLine.Length == 0)
                {
                    Console.WriteLine("Blank input line");
                    return;
                }

                if (inputLine[0] == '>')
                {
                    nContigs++;
                    outputFile.WriteLine(inputLine);
                    continue;
                }

                string outputLine = "";
                for (int i = 0; i < inputLine.Length; i++)
                {
                    if (!baseCounts.ContainsKey(inputLine[i]))
                    {
                        baseCounts.Add(inputLine[i], 0);
                    }

                    baseCounts[inputLine[i]]++;
                    nBases++;
                    switch (inputLine[i])
                    {
                        case 'A':
                        case 'C':
                        case 'T':
                        case 'G':
                        case 'N':
                        case 'a':
                        case 't':
                        case 'c':
                        case 'g':
                        case 'n':
                            outputLine += inputLine[i];
                            break;
                        default:
                            outputLine += 'N';
                            break;
                    }
                } // switch

                outputFile.WriteLine(outputLine);
            } // each line of the file

            inputFile.Close();
            outputFile.Close();

            Console.WriteLine("Base\tCount");
            foreach (var base_ in baseCounts.Select(_ => _.Key)) 
            {
                if (baseCounts[base_] != 0)
                {
                    Console.WriteLine(base_ + "\t" + baseCounts[base_]);
                }
            }

            Console.WriteLine();
            Console.WriteLine("Processed " + nContigs + " contigs in " + nLines + " lines with " + nBases + " bases in " + ASETools.ElapsedTimeInSeconds(timer));

        } // Main
    }
}
