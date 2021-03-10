using System;
using System.Collections.Generic;
using System.IO;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using ASELib;

namespace SplitSAM
{
    class Program
    {
        static void Main(string[] args)
        {
            if (args.Count() != 3)
            {
                Console.WriteLine("usage: SplitSAM input.sam outputBaseName nPieces");
                return;
            }

            int nPieces = Convert.ToInt32(args[2]);
            if (nPieces < 1)
            {
                Console.WriteLine("Invalid n Pieces");
                return;
            }

            var inputFile = ASETools.CreateStreamReaderWithRetry(args[0]);
            if (null == inputFile)
            {
                return;
            }

            var outputFiles = new List<StreamWriter>();

            for (int i = 0; i < nPieces; i++)
            {
                var outputFile = ASETools.CreateStreamWriterWithRetry(args[1] + i + ".sam");
                if (null == outputFile)
                {
                    return;
                }

                outputFiles.Add(outputFile);
            }

            var random = new Random();
            var halfSeenReads = new Dictionary<string, string>();  // Maps name -> read

            string line;
            while (null != (line = inputFile.ReadLine()))
            {
                if (line == "")
                {
                    Console.WriteLine("Empty input line");
                    return;
                }

                if (line[0] == '@')
                {
                    //
                    // Copy header lines into each output file
                    //
                    foreach (var outputFile in outputFiles)
                    {
                        outputFile.WriteLine(line);
                    }

                    continue;
                }


                var fields = line.Split('\t');
                if (fields.Count() < 11)
                {
                    Console.WriteLine("Too few fields: " + line);
                    return;
                }

                var flags = Convert.ToInt32(fields[1]);
                if ((flags & 0x900) != 0)
                {
                    //
                    // Supplementary and/or secondary alignments.  Skip.
                    //
                    continue;
                }

                if ((flags & 0x1) != 0)
                {
                    if (halfSeenReads.ContainsKey(fields[0]))
                    {
                        var outputFile = outputFiles[random.Next(nPieces)];

                        outputFile.WriteLine(line);
                        outputFile.WriteLine(halfSeenReads[fields[0]]);

                        halfSeenReads.Remove(fields[0]);
                    } else
                    {
                        halfSeenReads.Add(fields[0], line);
                    }
                }
                else
                {
                    //
                    // Single-end read, just write to a random file.
                    //
                    outputFiles[random.Next(nPieces)].WriteLine(line);
                }

            } // while we have input lines

            foreach (var outputFile in outputFiles)
            {
                outputFile.Close();
            }

            inputFile.Close();
        } // Main
    } // Program
} // namespace
