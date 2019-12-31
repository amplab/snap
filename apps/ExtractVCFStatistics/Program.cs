using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using ASELib;
using System.Diagnostics;

namespace ExtractVCFStatistics
{
    class Program
    {
        static ASETools.Configuration configuration;


        static void Main(string[] args)
        {
            var timer = new Stopwatch();
            timer.Start();

            configuration = ASETools.Configuration.loadFromFile(args);
            if (null == configuration)
            {
                Console.WriteLine("Unable to load configuration.");
                return;
            }


            var cases = ASETools.Case.LoadCases(configuration.casesFilePathname).Select(_ => _.Value).ToList();
            if (null == cases)
            {
                Console.WriteLine("Unable to load cases.");
                return;
            }

            if (configuration.commandLineArgs.Count() == 0 || configuration.commandLineArgs.Any(arg => !cases.Select(_ => _.case_id).Contains(arg)))
            {
                Console.WriteLine("usage: ExtractVCFStatistics {caseId}");
                return;
            }

            var casesToProcess = cases.Where(_ => configuration.commandLineArgs.Contains(_.case_id)).ToList();

            int nPerDot;
            ASETools.PrintMessageAndNumberBar("Processing", "cases", casesToProcess.Count(), out nPerDot);

            var threading = new ASETools.WorkerThreadHelper<ASETools.Case, int>(casesToProcess, HandleOneCase, null, null, nPerDot);
            threading.run();
            Console.WriteLine(ASETools.ElapsedTimeInSeconds(timer));
        } // Main

        static void HandleOneCase(ASETools.Case case_, int state)
        {
            var resultsThisCase = new Dictionary<string, int>();
            foreach (var chr in ASETools.chromosomes)
            {
                resultsThisCase.Add(chr, 0);
            }

            resultsThisCase.Add("other", 0);

            var inputFile = ASETools.CreateStreamReaderWithRetry(case_.vcf_filename);
            if (null == inputFile)
            {
                Console.WriteLine("Unable to open " + case_.vcf_filename);
                return;
            }

            string inputLine;
            while (null != (inputLine = inputFile.ReadLine()))
            {
                if (inputLine.StartsWith("#"))
                {
                    continue;   // Ignore the headers
                }

                var fields = inputLine.Split('\t');
                if (fields.Length < 9)
                {
                    Console.WriteLine(case_.vcf_filename + " has a line with too few fields: " + inputLine);
                    return;
                }

                var chr = fields[0].ToLower();
                if (!ASETools.chromosomes.Contains(chr))
                {
                    chr = "other";
                }

                resultsThisCase[chr]++;
            } // while we have an input line

            inputFile.Close();

            var outputFilename = ASETools.GetDirectoryFromPathname(case_.vcf_filename) + @"\" + case_.normal_dna_file_id + ASETools.vcfStatisticsExtension;
            var outputFile = ASETools.CreateStreamWriterWithRetry(outputFilename);

            outputFile.WriteLine("Chromosome\tCount of variants");
            foreach (var chr in resultsThisCase.Select(_ => _.Key))
            {
                outputFile.WriteLine(chr + "\t" + resultsThisCase[chr]);
            }

            outputFile.WriteLine("**done**");

            outputFile.Close();
        } // HandleOneCase

    } // Program
} // namespace ExtractVCFStatistics
