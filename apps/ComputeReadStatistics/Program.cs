using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using System.IO;
using System.Diagnostics;
using ASELib;
using System.Data;

namespace ComputeReadStatistics
{
    class Program
    {
        static ASETools.CommonData commonData;

        static void Main(string[] args)
        {
            commonData = ASETools.CommonData.LoadCommonData(args);

            if (null == commonData)
            {
                return;
            }

            if (commonData.configuration.commandLineArgs.Count() == 0 ||  commonData.configuration.commandLineArgs.Any( _ => !commonData.cases.ContainsKey(_)))
            {
                Console.WriteLine("usage: ComputeReadStatistics {caseID}");
                return;
            }

            var casesToProcess = commonData.listOfCases.Where(_ => commonData.configuration.commandLineArgs.Contains(_.case_id)).ToList();

            if (casesToProcess.Any(_ => _.normal_dna_filename == "" || _.tumor_dna_filename == "" || _.tumor_rna_filename == "" || _.normal_rna_file_id != "" && _.normal_rna_filename == ""))
            {
                Console.WriteLine("At least one case is missing a BAM.");
                return;
            }

            int nPerDot;
            ASETools.PrintMessageAndNumberBar("Processing", "cases", casesToProcess.Count(), out nPerDot);

            var threading = new ASETools.WorkerThreadHelper<ASETools.Case, int>(casesToProcess, HandleOneCase, null, null, nPerDot);
            threading.run();
            Console.WriteLine(ASETools.ElapsedTimeInSeconds(commonData.timer));
        } // Main


        static ASETools.ReadStatistics ComputeReadStaticticsForBAM(string inputFilename)
        {
            var readStatistics = new ASETools.ReadStatistics();

            var startInfo = new ProcessStartInfo(/*BJB commonData.configuration.binariesDirectory*/ @"c:\bolosky\bin\" + "samtools.exe", "view " + inputFilename);
            startInfo.RedirectStandardOutput = true;
            startInfo.UseShellExecute = false;

            Process process = Process.Start(startInfo);

            var inputFile = process.StandardOutput;

            string inputLine = "";
            while (null != (inputLine = inputFile.ReadLine()))
            {
                readStatistics.addSAMLine(inputLine, inputFilename);
            }

            return readStatistics;
        } // HandleOneBAM

        static void HandleOneCase(ASETools.Case case_, int state)
        {
            //try
            {
                var outputDirectory = ASETools.GetDirectoryFromPathname(case_.normal_dna_filename) + @"\..\..\derived_files\" + case_.case_id;
                Directory.CreateDirectory(outputDirectory);
                var outputFilename = outputDirectory + @"\" + case_.case_id + ASETools.readStatisticsExtension;

                var outputFile = ASETools.CreateStreamWriterWithRetry(outputFilename);
                if (null == outputFile)
                {
                    Console.WriteLine("Unable to open output file " + outputFilename);
                    return;
                }

                outputFile.WriteLine("Normal DNA");
                ComputeReadStaticticsForBAM(case_.normal_dna_filename).WriteToFile(outputFile);

                outputFile.WriteLine();
                outputFile.WriteLine("Tumor DNA");
                ComputeReadStaticticsForBAM(case_.tumor_dna_filename).WriteToFile(outputFile);

                outputFile.WriteLine();
                outputFile.WriteLine("Tumor RNA");
                ComputeReadStaticticsForBAM(case_.tumor_rna_filename).WriteToFile(outputFile);

                if (case_.normal_rna_filename != "")
                {
                    outputFile.WriteLine();
                    ComputeReadStaticticsForBAM(case_.normal_rna_filename).WriteToFile(outputFile);
                }

                outputFile.WriteLine();
                outputFile.WriteLine("**done**");
                outputFile.Close();
            } /*catch (Exception e)
            {
                Console.WriteLine("Case " + case_.case_id + " aborted because of exception: " + e.Message + ".");
                Console.WriteLine("Stack trace: ");
                Console.WriteLine(e.StackTrace);
            }*/
        } // HandleOneCase


    }
}
