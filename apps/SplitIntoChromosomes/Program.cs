using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using System.IO;

namespace SplitIntoChromosomes
{
    class Program
    {
        static void Main(string[] args)
        {
            if (args.Count() != 1 && args.Count() != 2)
            {
                Console.WriteLine("usage: SplitIntoChromosomes input.fa {linesPerChunk}");
                return;
            }

            int linesPerChunk = 0;
            if (args.Count() == 2)
            {
                linesPerChunk = Convert.ToInt32(args[1]);
            }

            StreamWriter outputFile = null;
            StreamReader inputFile = new StreamReader(args[0]);

            string inputLine;
            int linesThisChunk = 0;
            int chunkNumber = 0;
            string chromosomeName = null;
            string lineRemainder = "";
            while ((inputLine = inputFile.ReadLine()) != null)
            {
                if (inputLine.Count() > 0 && inputLine[0] == '>' || (linesPerChunk > 0 && linesThisChunk >= linesPerChunk))
                {
                    if (inputLine[0] == '>')
                    {
                        if (lineRemainder != "")
                        {
                            outputFile.WriteLine(lineRemainder.ToUpper());
                            lineRemainder = "";
                        }

                        chunkNumber = 0;
                        chromosomeName = inputLine.Substring(1);
                        if (chromosomeName.Contains(' '))
                        {
                            chromosomeName = chromosomeName.Substring(0, chromosomeName.IndexOf(' '));
                        }
                        if (chromosomeName.Count() < 3 || chromosomeName.Substring(0, 3).ToLower() != "chr")
                        {
                            chromosomeName = "chr" + chromosomeName;
                        }
                    }
                    else
                    {
                        chunkNumber++;
                    }

                    if (outputFile != null)
                    {
                        outputFile.Close();
                    }

                    if (linesPerChunk == 0)
                    {
                        outputFile = new StreamWriter(chromosomeName);
                    }
                    else
                    {
                        outputFile = new StreamWriter(chromosomeName + "." + chunkNumber);
                    }
                    linesThisChunk = 0;
                }

                if (outputFile == null)
                {
                    Console.WriteLine("No output file.  Either input FASTA doesn't start with a >, or we failed to open the output file.");
                    return;
                }

                if (inputLine[0] == '>')
                {
                    outputFile.WriteLine(inputLine);
                    lineRemainder = "";
                    linesThisChunk++;
                }
                else
                {
                    lineRemainder = lineRemainder + inputLine;
                    if (lineRemainder.Count() >= 100)
                    {
                        outputFile.WriteLine(lineRemainder.Substring(0, 100).ToUpper());
                        lineRemainder = lineRemainder.Substring(100);
                        linesThisChunk++;
                    }
                }

            }
            if (outputFile != null)
            {
                outputFile.Close();
            }
        }
    }
}
