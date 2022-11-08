using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using System.Diagnostics;
using System.IO;
using ASELib;

namespace SimulateShortReads
{
    internal class Program
    {
        static void WriteFASTQLine(StreamWriter outputFile, ASETools.FASTA fasta, int readLength, bool RC, string contigName, int offsetInContig, long readNumber)
        {
            outputFile.WriteLine("@simulated." + readNumber);
            string seq = "";
            string qual = "";
            var contigBytes = fasta.contigs[contigName];
            for (int offset = offsetInContig; offset < offsetInContig + readLength; offset++)
            {
                seq += contigBytes[offset];
                qual += "M";    // High quality (since, after all, we know it's perfect).
            }

            if (RC)
            {
                outputFile.WriteLine(ASETools.ReverseCompliment(seq));
            } else
            {
                outputFile.WriteLine(seq);
            }

            outputFile.WriteLine("+");
            outputFile.WriteLine(qual);
        } // WriteFASTQLine

        static void Main(string[] args)
        {
            var timer = new Stopwatch();
            timer.Start();

            if (args.Length != 4 && args.Length != 5)
            {
                Console.WriteLine("usage: SimulateShortReads input.fasta output_base nReads length {average inter-read gap}");
                Console.WriteLine("If you want paired reads, specify the average inter-read gap and the generated gap will");
                Console.WriteLine("be uniformly distributed between .5x-1.5x");
                return;
            }

            var inputFilename = args[0];
            var outputBaseFilename = args[1];
            long nReads = Convert.ToInt64(args[2]);
            int readLength = Convert.ToInt32(args[3]);
            bool paired = args.Length == 5;
            int interReadGap = 0;
            if (paired)
            {
                interReadGap = Convert.ToInt32(args[4]);
            }

            var fasta = ASETools.FASTA.loadFromFile(inputFilename);
            if (fasta == null)
            {
                Console.WriteLine("Unable to load " + inputFilename);
                return;
            }

            StreamWriter outputFile;
            StreamWriter outputFile2 = null;

            if (paired)
            {
                var outputFilename = outputBaseFilename + "_1.fastq";
                outputFile = ASETools.CreateStreamWriterWithRetry(outputFilename);
                if (outputFile == null)
                {
                    Console.WriteLine("Unable to open output file " + outputFilename);
                    return;
                }

                outputFilename = outputBaseFilename + "_2.fastq";
                outputFile2 = ASETools.CreateStreamWriterWithRetry(outputFilename);
                if (outputFile2 == null)
                {
                    Console.WriteLine("Unable to open output file " + outputFilename);
                    return;
                }
            }
            else
            {
                var outputFilename = outputBaseFilename + ".fastq";
                outputFile = ASETools.CreateStreamWriterWithRetry(outputFilename);
                if (outputFile == null)
                {
                    Console.WriteLine("Unable to open output file " + outputFilename);
                    return;
                }
            }

            var random = new Random();

            int nPerDot;
            if (nReads > 1000)
            {
                ASETools.PrintMessageAndNumberBar("Generating", "reads", nReads, out nPerDot);
            } else
            {
                nPerDot = 1000;
            }

            long nGenerated = 0;
            while(nGenerated < nReads)
            {
                string contigName;
                int offsetInContig;

                int gap;
                int offsetOfEndOfReads;
                if (paired)
                {
                    gap = random.Next(interReadGap / 2, interReadGap * 3 / 2);
                    offsetOfEndOfReads = gap + readLength * 2;
                } else
                {
                    gap = 0;
                    offsetOfEndOfReads = readLength;
                }

                do
                {
                    fasta.getRandomLocation(out contigName, out offsetInContig);
                } while (fasta.contigs[contigName].Length < offsetOfEndOfReads + offsetInContig);   // Keep going until we get one that's not off the end of the contig

                var oldNDots = nGenerated / nPerDot;

                WriteFASTQLine(outputFile, fasta, readLength, false, contigName, offsetInContig, nGenerated);

                if (paired)
                {
                    WriteFASTQLine(outputFile2, fasta, readLength, true, contigName, offsetInContig + gap + readLength, nGenerated);
                    nGenerated++;
                    if (nGenerated % nPerDot == 0) Console.Write(".");
                }

                nGenerated++;
                if (nGenerated % nPerDot == 0) Console.Write(".");
            } // for each read (pair) to generate

            outputFile.Close();
            if (paired)
            {
                outputFile2.Close();
            }

            if (nReads > 1000)
            {
                Console.WriteLine();
                Console.WriteLine("Generated " + nReads + " in " + ASETools.ElapsedTime(timer) + ", " + nReads * 1000 / timer.ElapsedMilliseconds + " reads/s");
            }

        } // Main
    } // Program
} // namespace
