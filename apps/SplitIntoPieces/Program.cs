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
            Console.WriteLine("usage: SplitIntoPieces inputFile nPieces lineGranularity outputNameBase outputNameExtension {-u} {-h headerLine} {-r}");
            Console.WriteLine(" -u means to use unix-style line terminators.");
            Console.WriteLine(" -h provides a header line to use in each output file.");
            Console.WriteLine(" -r means to randomize the lines both within and between groups");
        }
        static void Main(string[] args)
        {
            if (args.Count() < 5 || args.Count() > 9)
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
            bool randomize = false;

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
                else if ("-r" == args[i])
                {
                    randomize = true;
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


            if (randomize)
            {
                Random random = new Random();
                int nLineGroups = (lines.Count() + lineGranularity - 1) / lineGranularity;

                string[,] lineGroups = new string[nLineGroups, lineGranularity];
                for (i = 0; i < nLineGroups; i++)
                {
                    for (int j = 0; j < lineGranularity; j++)
                    {
                        if (i * lineGranularity + j >= nLines)
                        {
                            lineGroups[i, j] = "";
                        }
                        else
                        {
                            lineGroups[i, j] = lines[i * lineGranularity + j];
                        }
                    }
                }
                int nLineGroupsRemaining = nLineGroups;

                StreamWriter[] outputFiles = new StreamWriter[nPieces];
                for (i = 0; i < nPieces; i++)
                {
                    outputFiles[i] = new StreamWriter(outputNameBase + i + outputNameExtension);
                    if (headerLine != null)
                    {
                        if (useUnixLineEndings)
                        {
                            outputFiles[i].Write(headerLine + "\n");
                        }
                        else
                        {
                            outputFiles[i].WriteLine(headerLine);
                        }
                    }
                }

                int whichOutputFile = 0;
                while (nLineGroupsRemaining > 0)
                {
                    int whichLineGroup = random.Next() % nLineGroupsRemaining;
                    for (int j = 0; j < lineGranularity; j++)
                    {
                        if (useUnixLineEndings)
                        {
                            outputFiles[whichOutputFile].Write(lineGroups[whichLineGroup, j] + "\n"); ;
                        }
                        else
                        {
                            outputFiles[whichOutputFile].WriteLine(lineGroups[whichLineGroup, j]);
                        }
                        lineGroups[whichLineGroup, j] = lineGroups[nLineGroupsRemaining - 1, j];
                    }
                    whichOutputFile = (whichOutputFile + 1) % nPieces;
                    nLineGroupsRemaining--;
                }

                for (i = 0; i < nPieces; i++)
                {
                    outputFiles[i].Close();
                }
            }
            else
            {
                // Not randomized

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
}
