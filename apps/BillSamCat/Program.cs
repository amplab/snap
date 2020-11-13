using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using System.IO;
using ASELib;
using System.Diagnostics;

namespace BillSamCat
{
    class Program
    {
        static void Main(string[] args)
        {
            if (args.Count() < 3)
            {
                Console.WriteLine("usage: BillSamCat outputFile.sam inputFile0.bam inputFile1.bam ...");
                return;
            }


            var outputFilename = args[0];
            if (File.Exists(outputFilename))
            {
                Console.WriteLine("Cannot write to an existing file (so you don't accidentally overwrite something)");
                return;
            }

            var outputFile = new StreamWriter(outputFilename);

            bool fixedRG = false;

            for (int i = 1; i < args.Count(); i++)
            {
                var startInfo = new ProcessStartInfo("samtools.exe", "view " + ((i == 1) ? "-h " : "") + args[i]);

                startInfo.UseShellExecute = false;
                startInfo.RedirectStandardOutput = true;

                Process process;
                try
                {
                    process = Process.Start(startInfo);
                }
                catch (Exception e)
                {
                    Console.WriteLine("Error trying to start process samtools.exe");
                    Console.WriteLine("Exception message: " + e.Message);

                    throw e;
                }

                string outputLine;
                while (null != (outputLine = process.StandardOutput.ReadLine()))
                {
                    if (!fixedRG && outputLine.StartsWith("@RG"))
                    {
                        outputLine = outputLine.Replace("PLILLUMINA", "PL:ILLUMINA");
                        fixedRG = true;
                    }

                    outputFile.WriteLine(outputLine);
                }

                process.WaitForExit();
            } // for each input file

            outputFile.Close();

        } // Main
    }// Program
} // namespace
