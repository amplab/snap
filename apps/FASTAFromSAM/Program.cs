using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using ASELib;
using System.Diagnostics;

//
// Take a SAM file and generate a FASTA for each read.  The contig name is the QNAME of the read.
// This is intended for long reads.
//


namespace FASTAFromSAM
{
    internal class Program
    {
        static int nSAMLinesProcessed = 0;
        static long nBasesProcessed = 0;

        const int nBasesPerFASTALine = 100;

        static void ProcessSAMLine(string outputBaseName, ASETools.SAMLine samLine)
        {
            if (samLine.seq == "*" || (samLine.flag & (ASETools.SAMLine.SecondaryAligment | ASETools.SAMLine.SupplementaryAlignment)) != 0)
            {
                return;
            }

            var outputFilename = outputBaseName + samLine.qname + ".fasta";
            var outputFile = ASETools.CreateStreamWriterWithRetry(outputFilename);
            if (null == outputFile)
            {
                Console.WriteLine("Unable to open output file " + outputFilename);
                Environment.Exit(1);
            }

            outputFile.WriteLine(">" + samLine.qname);

            int offsetInSEQ = 0;
            while (offsetInSEQ < samLine.seq.Length)
            {
                int amountToWrite = Math.Min(nBasesPerFASTALine, samLine.seq.Length - offsetInSEQ);
                outputFile.WriteLine(samLine.seq.Substring(offsetInSEQ, amountToWrite));
                offsetInSEQ += amountToWrite;
            }

            outputFile.Close();

            nSAMLinesProcessed++;
            nBasesProcessed += samLine.seq.Length;
        } // ProcessSAMLine


        static void Main(string[] args)
        {
            var timer = new Stopwatch();
            timer.Start();

            if (args.Length != 2)
            {
                Console.WriteLine("usage: FASTAFromSAM inputSAM outputBaseName");
                return;
            }

            var inputFilename = args[0];
            var outputBaseName = args[1];

            var inputFile = ASETools.CreateStreamReaderWithRetryCompressedBasedOnFilename(inputFilename);
            if (null == inputFile)
            {
                Console.WriteLine("Unable to to open " + inputFilename);
                return;
            }

            ASETools.SAMLine.ReadSAMLinesFromFile(inputFile, _ => ProcessSAMLine(outputBaseName, _));

            inputFile.Close();

            Console.WriteLine("Processed " + nBasesProcessed + " in " + nSAMLinesProcessed + " SAM lines in " + ASETools.ElapsedTime(timer));
        } // Main
    } // Program
} // namespace
