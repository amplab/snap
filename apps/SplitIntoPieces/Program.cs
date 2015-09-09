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
            Console.WriteLine("usage: SplitIntoPieces inputFile nPieces lineGranularity outputNameBase outputNameExtension");
        }
        static void Main(string[] args)
        {
            if (args.Count() != 5)
            {
                PrintUsage();
                return;
            }

            string inputFilename = args[0];
            string outputNameBase = args[3];
            string outputNameExtension = args[4];

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

            for (int i = 0; i < nPieces; i++)
            {
                StreamWriter outputFile = new StreamWriter(outputNameBase + i + outputNameExtension);

                for (int j = i * nLinesPerPiece; j < (i + 1) * nLinesPerPiece && j < nLines; j++)
                {
                    outputFile.WriteLine(lines[j]);
                }

                outputFile.Close();
            }
        }
    }
}
