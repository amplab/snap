using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using System.IO;

namespace SplitIntoPieces
{
    class Program
    {
        static void PrintUsage()
        {
            Console.WriteLine("usage: SplitIntoPieces inputFile nPieces lineGranularity outputNameBase outputNameExtension {-u} {-h headerLine}");
            Console.WriteLine(" -u means to use unix-style line terminators.");
            Console.WriteLine(" -h provides a header line to use in each output file.");
        }
        static void Main(string[] args)
        {
            if (args.Count() < 5 || args.Count() > 8)
            {
                Console.WriteLine("Incorrect arg count " + args.Count() + " must be between 5 and 8");
                PrintUsage();
                return;
            }

            string inputFilename = args[0];
            string outputNameBase = args[3];
            string outputNameExtension = args[4];
            string headerLine = null;
            bool useUnixLineEndings = false;

            int i = 5;
            while (i < args.Count())
            {
                if ("-u" == args[i])
                {
                    useUnixLineEndings = true;
                }
                else if ("-h" == args[i] && i < args.Count() - 1)
                {
                    i++;
                    headerLine = args[i];
                }
                else
                {
                    Console.WriteLine("Couldn't parse arg " + args[i]);
                    PrintUsage();
                    return;
                }

                i++;
            }

            int nPieces = Convert.ToInt32(args[1]);
            if (nPieces < 1)
            {
                Console.WriteLine("nPieces must be >= 1");
                PrintUsage();
                return;
            }

            int lineGranularity = Convert.ToInt32(args[2]);
            if (lineGranularity < 1)
            {
                Console.WriteLine("lineGranularity must be >= 1");
                PrintUsage();
                return;
            }

            string[] lines = File.ReadAllLines(args[0]);
            int nLines = lines.Count();

            int nLinesPerPiece = nLines / nPieces;
            if (nLinesPerPiece % lineGranularity != 0)
            {
                nLinesPerPiece += lineGranularity - nLinesPerPiece % lineGranularity;
            }

            for (i = 0; i < nPieces; i++)
            {
                StreamWriter outputFile = new StreamWriter(outputNameBase + i + outputNameExtension);
                if (headerLine != null)
                {
                    if (useUnixLineEndings)
                    {
                        outputFile.Write(headerLine + "\n");
                    }
                    else
                    {
                        outputFile.WriteLine(headerLine);
                    }
                }

                for (int j = i * nLinesPerPiece; j < (i + 1) * nLinesPerPiece && j < nLines; j++)
                {
                    if (useUnixLineEndings)
                    {
                        outputFile.Write(lines[j] + "\n");
                    }
                    else
                    {
                        outputFile.WriteLine(lines[j]);
                    }
                }

                outputFile.Close();
            }
        }
    }
}
