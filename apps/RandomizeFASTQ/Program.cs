using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using System.IO;
using System.Diagnostics;
using System.Runtime.Remoting.Metadata;

namespace RandomizeFASTQ
{
    class Program
    {
        class FASTQRead
        {
            public string[] lines = new string[3];  // Skip the third line, which is just a "+"

            public static FASTQRead readFromFile(StreamReader inputFile)
            {
                var read = new FASTQRead();
                read.lines[0] = inputFile.ReadLine();
                if (read.lines[0] == null)
                {
                    return null;    // EOF
                }
                if (read.lines[0][0] != '@')
                {
                    throw new Exception("FASTQRead.readFromFile: first line (read name) doesn't start with an @: " + read.lines[0]);
                }
                read.lines[1] = inputFile.ReadLine();
                var plusLine = inputFile.ReadLine();
                if (plusLine[0] != '+')
                {
                    throw new Exception("FASTQRead.readFromFile: third line doesn't start with a  +:" + plusLine[0]);
                }

                read.lines[2] = inputFile.ReadLine();
                if (read.lines[1].Length != read.lines[2].Length)
                {
                    throw new Exception("FASTQRead.readFromFile: bases and base call quality are of different length: " + read.lines[1] + " and " + read.lines[2]);
                }

                return read;
            } // readFromFile

            public void writeToFile(StreamWriter outputFile)
            {
                outputFile.WriteLine(lines[0]);
                outputFile.WriteLine(lines[1]);
                outputFile.WriteLine("+");
                outputFile.WriteLine(lines[2]);
            }
        } // class FASTQRead


        static void Main(string[] args)
        {
            if (args.Count() != 2 && args.Count() != 4)
            {
                Console.WriteLine("usage: randomizeFastq {input output} {input1 input2 output1 output2}");
                Console.WriteLine("Use the second form for paired reads.");
                return;
            }

            StreamReader[] inputFiles;
            StreamWriter[] outputFiles;

            if (args.Count() == 2)
            {
                inputFiles = new StreamReader[1];
                outputFiles = new StreamWriter[1];

                inputFiles[0] = new StreamReader(args[0]);
                outputFiles[1] = new StreamWriter(args[1]);
            } else
            {
                inputFiles = new StreamReader[2];
                outputFiles = new StreamWriter[2];

                inputFiles[0] = new StreamReader(args[0]);
                inputFiles[1] = new StreamReader(args[1]);

                outputFiles[0] = new StreamWriter(args[2]);
                outputFiles[1] = new StreamWriter(args[3]);
            }

            int nInputFiles = inputFiles.Count();

            int nReads = 0;
            FASTQRead newRead;
            List<FASTQRead>[] reads = new List<FASTQRead>[nInputFiles];
            for (int i = 0; i < nInputFiles; i++)
            {
                reads[i] = new List<FASTQRead>();
            }

            Console.Write("Reading inputs, one dot/100000 read (pairs): ");

            while (null != (newRead = FASTQRead.readFromFile(inputFiles[0])))
            {
                reads[0].Add(newRead);
                for (int i = 1; i < nInputFiles; i++)
                {
                    newRead = FASTQRead.readFromFile(inputFiles[i]);
                    if (newRead == null)
                    {
                        throw new Exception("Unequal number of reads in input files.");
                    }

                    reads[i].Add(newRead);
                }

                nReads++;
                if (nReads % 100000 == 0)
                {
                    Console.Write(".");
                }
            }

            for (int i = 1; i < nInputFiles; i++)
            {
                if (null != FASTQRead.readFromFile(inputFiles[i]))
                {
                    throw new Exception("Unequal number of reads in input files (at end).");
                }
            }

            Console.WriteLine("");
            Console.Write("Writing " + nReads + " read (pairs), one dot/100000:");
            int nWritten = 0;
            var random = new Random();
            for (int maxRead = nReads - 1; maxRead >= 0; maxRead--)
            {
                var indexToEmit = random.Next(maxRead);

                for (int i = 0; i < nInputFiles; i++)
                {
                    reads[i][indexToEmit].writeToFile(outputFiles[i]);
                    reads[i][indexToEmit] = reads[i][maxRead];  // Move the max read into the slot we just wrote.
                }

                nWritten++;
                if (nWritten % 100000 == 0)
                {
                    Console.Write(".");
                }
            }
            Console.WriteLine();

            for (int i = 0; i < nInputFiles; i++)
            {
                inputFiles[i].Close();
                outputFiles[i].Close();
            }
        } // Main
    }
}
