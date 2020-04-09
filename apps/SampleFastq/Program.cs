using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using System.IO;
using ASELib;
using System.Diagnostics;
using System.Security.Cryptography;

namespace SampleFastq
{
    class Program
    {
        class FASTQRead
        {
            byte[] data;
            int idLength;
            int queryLength; // equals quality length

            static public FASTQRead ReadFromStream(StreamReader inputFile)
            {
                var idLine = inputFile.ReadLine();
                var queryLine = inputFile.ReadLine();
                var plusLine = inputFile.ReadLine();
                var qualityLine = inputFile.ReadLine();

                if (qualityLine == null)
                {
                    if (idLine != null)
                    {
                        throw new Exception("Ill-formatted FASTQ file or read error.");
                    }

                    return null;
                }

                if (!idLine.StartsWith("@"))
                {
                    throw new Exception("FASTQ ID line doens't start with an @");
                }

                if (plusLine != "+")
                {
                    throw new Exception("FASTQ plus line isn't a plus.");
                }

                if (queryLine.Length != qualityLine.Length)
                {
                    throw new Exception("FASTQ query and quality lines don't match in length");
                }

                var read = new FASTQRead();
                read.data = new byte[idLine.Length - 1 + 2 * queryLine.Length]; // -1 is because we don't bother to store the @
                read.idLength = idLine.Length - 1;
                read.queryLength = queryLine.Length;

                int offsetInData = 0;
                for (int i = 1; i < idLine.Length; i++)
                {
                    read.data[offsetInData++] = Convert.ToByte(idLine[i]);
                }

                for (int i = 0; i < queryLine.Length; i++)
                {
                    read.data[offsetInData++] = Convert.ToByte(queryLine[i]);
                }

                for (int i = 0; i < qualityLine.Length; i++)
                {
                    read.data[offsetInData++] = Convert.ToByte(qualityLine[i]);
                }
                return read;
            } // factory

            public void WriteToStream(StreamWriter outputFile)
            {
                string idLine = "@";
                for (int i = 0; i < idLength; i++)
                {
                    idLine += Convert.ToChar(data[i]);
                }
                outputFile.WriteLine(idLine);

                string queryString = "";
                for (int i = idLength; i < idLength + queryLength; i++)
                {
                    queryString += Convert.ToChar(data[i]);
                }

                outputFile.WriteLine(queryString);
                outputFile.WriteLine("+");

                string quality = "";
                for (int i = idLength + queryLength; i < idLength + 2 * queryLength; i++)
                {
                    quality += Convert.ToChar(data[i]);
                }

                outputFile.WriteLine(quality);
            }
        }
        static void Main(string[] args)
        {
            if (args.Count() != 3 && args.Count() != 5)
            {
                Console.WriteLine("usage: SampleFastq nReadsToWritePerInput input1 {input2} output1 {output2}");
                return;
            }

            Stopwatch timer = new Stopwatch();
            timer.Start();

            int nInputFiles = args.Count() == 5 ? 2 : 1;

            string input0Name = args[1];
            string input1Name = (nInputFiles == 2) ? args[2] : null;
            string output0Name = (nInputFiles == 2) ? args[3] : args[2];
            string output1Name = (nInputFiles == 2) ? args[4] : null;
 
            var nReadsToWrite = Convert.ToInt64(args[0]);

            var input0 = ASETools.CreateStreamReaderWithRetryCompressedBasedOnFilename(input0Name);
            if (null == input0)
            {
                Console.WriteLine("Unable to open input file " + input0Name);
                return;
            }

            StreamReader input1;
            if (nInputFiles == 2)
            {
                input1 = ASETools.CreateStreamReaderWithRetryCompressedBasedOnFilename(input1Name);
                if (input1 == null)
                {
                    Console.WriteLine("Unable to open file " + input1Name);
                    return;
                }
            }
            else
            {
                input1 = null;
            }

            var output0 = ASETools.CreateStreamWriterWithRetry(output0Name);
            if (null == output0)
            {
                Console.WriteLine("Unable to open output file " + output0Name);
                return;
            }

            StreamWriter output1;
            if (nInputFiles == 2)
            {
                output1 = ASETools.CreateStreamWriterWithRetry(output1Name);
                if (null == output1)
                {
                    Console.WriteLine("Unable to open output file " + output1Name);
                    return;
                }
            } else
            {
                output1 = null;
            }

            StreamReader[] inputFiles = { input0, input1};
            StreamWriter[] outputFiles = { output0, output1 };

            var reads = new List<FASTQRead>[nInputFiles];
            for (int whichFile = 0; whichFile < nInputFiles; whichFile++)
            {
                reads[whichFile] = new List<FASTQRead>();
            }

            var random = new Random();

            Console.Write("Reading input files, 1 dot/1,000,000 reads (in each file): ");
            int nRead = 0;

            while (true)
            {
                int nFailedReads = 0;
                for (int whichFile = 0; whichFile < nInputFiles; whichFile++)
                {
                    var read = FASTQRead.ReadFromStream(inputFiles[whichFile]);
                    if (read == null)
                    {
                        nFailedReads++;
                    } else
                    {
                        reads[whichFile].Add(read);
                    }
                }

                if (nFailedReads != 0)
                {
                    if (nFailedReads != nInputFiles)
                    {
                        Console.WriteLine("Read error or input files didn't have the same number of reads");
                        return;
                    }
                    break;
                }

                nRead++;
                if (nRead % 1000000 == 0)
                {
                    Console.Write(".");
                }
            } // while we have something to read
            Console.WriteLine();

            int []indexArray = new int[nRead];
            for (int i = 0; i < nRead; i++)
            {
                indexArray[i] = i;
            }

            int nPerDot;
            ASETools.PrintMessageAndNumberBar("Writing", "reads", nReadsToWrite, out nPerDot);

            int nRemaining = nRead;
            for (int nWritten = 0; nWritten < nReadsToWrite; nWritten++)
            {
                int whichIndex = random.Next(0, nRemaining - 1);
                int whichRead = indexArray[whichIndex];

                for (int whichFile = 0; whichFile < nInputFiles; whichFile++)
                {
                    reads[whichFile][whichRead].WriteToStream(outputFiles[whichFile]);
                }

                indexArray[whichIndex] = indexArray[nRemaining - 1];
                nRemaining--;

                if ((nWritten + 1) % nPerDot == 0)
                {
                    Console.Write(".");
                }
            }

            foreach (var outputFile in outputFiles)
            {
                outputFile.Close();
            }

            Console.WriteLine();
            Console.WriteLine("Wrote " + nReadsToWrite + " FASTQ reads in " + ASETools.ElapsedTimeInSeconds(timer));
        } // Main
    }
}
