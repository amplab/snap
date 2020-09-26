using System;
using System.Collections.Generic;
using System.ComponentModel.Design;
using System.IO;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using ASELib;

namespace ExtractReadsFromFASTQ
{
    class Program
    {
        static void Main(string[] args)
        {
            if (args.Count() != 3 && args.Count() != 5)
            {
                Console.WriteLine("usage: ExtractReadsFromFASTQ readIDsToKeepFile  {input output | input1 input2 output1 output2}");
                return;
            }

            var IDsToKeep = File.ReadAllLines(args[0]);
            if (IDsToKeep == null || IDsToKeep.Count() == 0)
            {
                Console.WriteLine("Couldn't read IDs to keep from " + args[0]);
                return;
            }

            int nInputs = (args.Count() == 3) ? 1 : 2;

            var inputFiles = new StreamReader[nInputs];
            var outputFiles = new StreamWriter[nInputs];

            for (int i = 0; i < nInputs; i++)
            {
                inputFiles[i] = ASETools.CreateStreamReaderWithRetry(args[i + 1]);

                if (inputFiles[i] == null)
                {
                    Console.WriteLine("Unable to open input file " + args[i + 1]);
                    return;
                }
             }

            for (int i = 0; i < nInputs; i++)
            {
                outputFiles[i] = ASETools.CreateStreamWriterWithRetry(args[i + 1 + nInputs]);
                if (outputFiles[i] == null)
                {
                    Console.WriteLine("Can't open output file " + args[i + 1 + nInputs]);
                }
            }

            var inputLines = new string[nInputs, 4];
            while (true)
            {
                for (int i = 0; i < nInputs; i++)
                {
                    for (int j = 0; j < 4; j++)
                    {
                        inputLines[i, j] = inputFiles[i].ReadLine();
                    }
                }

                if (inputLines[0,0] == null)
                {
                    break;
                }

                for (int i = 0; i < nInputs; i++)
                {
                    if (!inputLines[i,0].StartsWith("@"))
                    {
                        Console.WriteLine("Malformed FASTQ: line 1 of a read doesn't start with @");
                        return;
                    }

                    if (inputLines[i,2] != "+")
                    {
                        Console.WriteLine("Malformed FASTQ: line 3 of a read isn't +");
                        return;
                    }
                }

                if (IDsToKeep.Any(_ => inputLines[0,0].Contains("@" + _ + " ") || inputLines[0,0] == "@" + _))
                {
                    for (int i = 0; i < nInputs; i++)
                    {
                        for (int j = 0; j < 4; j++)
                        {
                            outputFiles[i].WriteLine(inputLines[i, j]);
                        }
                    }
                }
            }

            for (int i = 0; i < nInputs; i++)
            {
                outputFiles[i].Close();
                inputFiles[i].Close();
            }


        } // Main
    } // Program
} // Namespace
