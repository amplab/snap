using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using System.IO;
using ASELib;
using System.Diagnostics;

namespace AddFASTQComments
{
    internal class Program
    {
        static void Main(string[] args)
        {
            var timer = new Stopwatch();
            timer.Start();

            if (args.Count() < 3)
            {
                Console.WriteLine("usage: addFASTQComments input.FASTQ output.FASTQ comment");
                Console.WriteLine("The comment can be multiple args, which are then tab separated in the output.");
                Console.WriteLine("Any existing FASTQ comments are not preserved.");
                Console.WriteLine("Input files can be compressed if they end in an appropriate extension, but output will be");
                Console.WriteLine("uncompressed regardless of filename.");
                return;
            }

            var inputFile = ASETools.CreateStreamReaderWithRetryCompressedBasedOnFilename(args[0]);
            if (null == inputFile)
            {
                Console.WriteLine("Can't open input file " + args[0]);
                return;
            }

            var outputFile = ASETools.CreateStreamWriterWithRetry(args[1]);
            if (null == outputFile)
            {
                Console.WriteLine("Can't open output file " + args[1]);
                return;
            }

            string comment = args[2];
            for (int i = 3; i < args.Count(); i++)
            {
                comment += "\t" + args[i];
            }


            Console.Write("Processing, 1 dot/100K reads: ");

            long lineCount = 0;
            string line;
            while (null != (line = inputFile.ReadLine()))
            {
                if (lineCount % 4 == 0)
                {
                    string leadingPart;
                    if (line.Contains(" ") || line.Contains("\t"))
                    {
                        leadingPart = line.Substring(0, line.IndexOf(' '));
                    } else
                    {
                        leadingPart = line;
                    }

                    if (leadingPart.Contains("\t"))
                    {
                        leadingPart = line.Substring(0, line.IndexOf('\t'));
                    }

                    outputFile.WriteLine(leadingPart + " " + comment);
                } 
                else
                {
                    outputFile.WriteLine(line);
                }

                lineCount++;

                if (lineCount % 400000 == 0)    // 4 input lines/read
                {
                    Console.Write(".");
                }
            } // while we have a line

            outputFile.Close();
            inputFile.Close();

            Console.WriteLine();
            if (lineCount % 4 != 0)
            {
                Console.WriteLine("Malformed input FASTQ doesn't have a multiple of 4 lines");
            }

            Console.WriteLine("Processed " + lineCount / 4 + " reads in " + ASETools.ElapsedTime(timer));

        } // Main
    } // Program
} // namespace
