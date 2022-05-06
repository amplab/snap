using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using System.IO;
using ASELib;
using System.Threading;
using System.Diagnostics;


//
// This fixes some of the Novoalign runs where the read group line has PLILLUMINA instead of PL:ILLUMINA
//

namespace FixBrokenPL
{
    internal class Program
    {



        static bool fixedPL = false;

        static void GotLine(string line, Process writerProcess)
        {
            if (!fixedPL && line.StartsWith("@RG"))
            {
                fixedPL = true;
                writerProcess.StandardInput.WriteLine(line.Replace("PLILLUMINA", "PL:ILLUMINA"));
            } else
            {
                writerProcess.StandardInput.WriteLine(line);
            }
        } // GotLine

        static void OutputThread(StreamReader stdout, StreamWriter bamFile)
        {
            var binaryReader = new BinaryReader(stdout.BaseStream);
            var binaryWriter = new BinaryWriter(bamFile.BaseStream);

            var bufferSize = 20 * 1024 * 1024;

            long totalRead = 0;

            while (true)
            {
                var buffer = binaryReader.ReadBytes(bufferSize);
                if (buffer == null || buffer.Length == 0)
                {
                    break;  // EOF
                }

                binaryWriter.Write(buffer);

                if (totalRead / (1024 * 1024 * 1024) != (totalRead + buffer.Length) / (1024 * 1024 * 1024) && totalRead != 0)
                {
                    Console.Write(".");
                }

                totalRead += buffer.Length;
            }
        }

        static void Main(string[] args)
        {
            var timer = new Stopwatch();
            timer.Start();


            if (args.Length != 2)
            {
                Console.WriteLine("usage: FixBrokenPL input.bam output.bam");
                return;
            }

            Console.Write("Processing (one dot/GB):");

            StreamWriter bamFile = ASETools.CreateStreamWriterWithRetry(args[1]);
            if (bamFile == null)
            {
                return;
            }

            var samtoolsBinary = @"c:\bolosky\bin\samtools.exe";

            var startInfo = new ProcessStartInfo(samtoolsBinary, "view -h -b -S -");

            startInfo.UseShellExecute = false;
            startInfo.RedirectStandardOutput = true;
            startInfo.RedirectStandardInput = true;

            Process writerProcess;
            try
            {
                writerProcess = Process.Start(startInfo);
            }
            catch (Exception e)
            {
                Console.WriteLine("Error trying to start process " + samtoolsBinary);
                Console.WriteLine("Exception message: " + e.Message);

                return;
            }

            var outputThread = new Thread(() => OutputThread(writerProcess.StandardOutput, bamFile));   // This takes the output from the writer and writes it to the file system.
            outputThread.Start();

            ASETools.RunProcess(samtoolsBinary, "view -h " + args[0], _ => GotLine(_, writerProcess));  // This is the reader.  GotLine will fix the line if necessary and pass it to the writer

            outputThread.Join();

            bamFile.Close();

            Console.WriteLine();
            Console.WriteLine("Indexing...");

            ASETools.RunAndWaitForProcess(samtoolsBinary, "index " + args[1] + " " + args[1] + ".bai");

            Console.WriteLine();
            Console.WriteLine("Total elapsed time: " + ASETools.ElapsedTime(timer));

        } // Main
    } // program
} // namespace
