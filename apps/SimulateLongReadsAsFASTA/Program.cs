using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using ASELib;

namespace SimulateLongReadsAsFASTA
{
    internal class Program
    {
        static void Main(string[] args)
        {
            if (args.Count() != 4 && args.Count() != 5 || args.Count() == 5 && args[4] != "-n")
            {
                Console.WriteLine("usage: SimulateLongReadsAsFASTA inputFASTA outputFASTABase readLength coverage {-n}");
                Console.WriteLine("-n means not to generate FASTAs, just to make the single FASTQ output.");
                return;
            }

            bool nullFASTA = args.Count() == 5; // Yick.

            var coverage = Convert.ToDouble(args[3]);
            if (coverage <= 0)
            {
                Console.WriteLine("Coverage must be strictly positive");
                return;
            }

            var readLength = Convert.ToInt32(args[2]);
            if (readLength <= 0)
            {
                Console.WriteLine("read length must be strictly positive");
            }

            var inputFasta = ASETools.FASTA.loadFromFile(args[0]);
            if (null == inputFasta)
            {
                Console.WriteLine("unable to load FASTA from " + args[0]);
                return;
            }

            long totalBasesEmitted = 0;
            int outputFASTANumber = 0;

            var outputFastq = ASETools.CreateStreamWriterWithRetry(args[1] + ".fastq");
            if (outputFastq == null)
            {
                return;
            }

            var qualityString = ""; // for the FASTQ
            for (int i = 0; i < readLength; i++)
            {
                qualityString += "M";   // A good quality (since we know it's actually perfect)
            }

            while ((double)totalBasesEmitted < inputFasta.totalSize * coverage)
            {
                string contigName;
                int contigOffset;

                //
                // Find a random location that doesn't hang off the end of the contig
                //
                do
                {
                    inputFasta.getRandomLocation(out contigName, out contigOffset);
                } while (inputFasta.contigs[contigName].Length <= contigOffset + readLength);

                var outputFilename = args[1] + "_" + outputFASTANumber + ".fasta";
                var outputFile = nullFASTA ? System.IO.StreamWriter.Null : ASETools.CreateStreamWriterWithRetry(outputFilename);
                if (null == outputFile)
                {
                    Console.WriteLine("Unable to create output file " + outputFilename);
                    return;
                }

                outputFile.WriteLine(">" + contigName + "_" + contigOffset);

                outputFastq.WriteLine("@" + contigName + "_" + contigOffset);

                string outputLine = "";
                for (int i = contigOffset - 1; i <= contigOffset + readLength - 1; i++) // -1 is because contigs are 1 based, but the array is 0 based (feh)
                {
                    outputLine += inputFasta.contigs[contigName][i];
                    if (outputLine.Length >= 100)
                    {
                        outputFile.WriteLine(outputLine);
                        outputFastq.Write(outputLine);
                        outputLine = "";
                    }
                }

                if (outputLine != "")
                {
                    outputFile.WriteLine(outputLine);
                }

                outputFastq.WriteLine();
                outputFastq.WriteLine("+");
                outputFastq.WriteLine(qualityString);

                outputFile.Close();

                outputFASTANumber++;

                totalBasesEmitted += readLength;
            } // while we're emitting reads

            outputFastq.Close();
        } // Main
    } // Program
} // namespace
