using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using System.IO;
using ASELib;
using System.Diagnostics;

namespace MappedBaseCountDistribution
{
    class Program
    {
        static void HandleOne(string filename, ASETools.PreBucketedHistogram histogram)
        {
            if (filename == "")
            {
                return;
            }

            var mappedBaseCount = ASETools.MappedBaseCount.readFromFile(filename);
            if (mappedBaseCount == null)
            {
                Console.WriteLine("Unable to read mapped base count from " + filename);
                return;
            }

            histogram.addValue(mappedBaseCount.basesCovered);
        }

        static void WriteOne(string name, ASETools.HistogramResultLine[] lines, StreamWriter outputFile)
        {
            outputFile.WriteLine();
            outputFile.WriteLine(name);
            outputFile.WriteLine(ASETools.HistogramResultLine.Header());
            for (int i = 0; i < lines.Count(); i++)
            {
                outputFile.WriteLine(lines[i]);
            }
        }

        static void Main(string[] args)
        {

            var timer = new Stopwatch();
            timer.Start();

            var configuration = ASETools.Configuration.loadFromFile(args);
            if (null == configuration)
            {
                Console.WriteLine("Unable to load configuration");
                return;
            }

            var cases = ASETools.Case.LoadCases(configuration.casesFilePathname).Select(x => x.Value).ToList(); ;
            if (null == cases)
            {
                Console.WriteLine("Unable to load cases");
                return;
            }

            const double maxbases = 4000000000;
            const double increment = 50000000;
            var tumorDNAHistogram = new ASETools.PreBucketedHistogram(0, maxbases, increment);
            var normalDNAHistogram = new ASETools.PreBucketedHistogram(0, maxbases, increment);
            var tumorRNAHistogram = new ASETools.PreBucketedHistogram(0, maxbases, increment);
            var normalRNAHistogram = new ASETools.PreBucketedHistogram(0, maxbases, increment);

            Console.Write("Processing " + cases.Count() + " cases, 1 dot/100: ");
            int nCasesProcessed = 0;

            foreach (var case_ in cases)
            {
                HandleOne(case_.tumor_dna_mapped_base_count_filename, tumorDNAHistogram);
                HandleOne(case_.normal_dna_mapped_base_count_filename, normalDNAHistogram);
                HandleOne(case_.tumor_rna_mapped_base_count_filename, tumorRNAHistogram);
                HandleOne(case_.normal_rna_mapped_base_count_filename, normalRNAHistogram);

                nCasesProcessed++;
                if (nCasesProcessed % 100 == 0)
                {
                    Console.Write(".");
                }
            }

            var outputFile = ASETools.CreateStreamWriterWithRetry(configuration.finalResultsDirectory + ASETools.mappedBaseCountDistributionFilename);
            if (outputFile == null)
            {
                Console.WriteLine("Unable to open output file " + configuration.finalResultsDirectory + ASETools.mappedBaseCountDistributionFilename);
                return;
            }

            outputFile.WriteLine("Bucket\tnormal DNA\ttumor DNA\tnormal RNA\ttumorRNA");
            var normalDNALines = normalDNAHistogram.ComputeHistogram();
            var tumorDNALines = tumorDNAHistogram.ComputeHistogram();
            var normalRNALines = normalRNAHistogram.ComputeHistogram();
            var tumorRNALines = tumorRNAHistogram.ComputeHistogram();

            for (int i = 0; i < normalDNALines.Count(); i++)
            {
                outputFile.WriteLine(normalDNALines[i].minValue + "\t" + normalDNALines[i].cdfValue + "\t" + tumorDNALines[i].cdfValue + "\t" + normalRNALines[i].cdfValue + "\t" + tumorRNALines[i].cdfValue);
            }

            WriteOne("normal DNA", normalDNALines, outputFile);
            WriteOne("tumor DNA", tumorDNALines, outputFile);
            WriteOne("normal RNA", normalRNALines, outputFile);
            WriteOne("tumor RNA", tumorRNALines, outputFile);

            outputFile.WriteLine("**done**");
            outputFile.Close();
            Console.WriteLine(ASETools.ElapsedTimeInSeconds(timer));
        }
    }
}
