using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using ASELib;
using System.IO;


namespace SeedSweep
{
    class Program
    {
 
        static void Main(string[] args)
        {
            var timingDirectory = @"d:\temp\timings\";
            var concordanceDirectory = @"d:\temp\concordance-dragen3.8VC\";

            string[] samplesToRun = { "hg003", "hg007", "ERR194147" };

            var outputFilename = @"\temp\SNAPSeedSweep.txt";
            var outputFile = ASETools.CreateStreamWriterWithRetry(outputFilename);

            if (null == outputFile)
            {
                Console.WriteLine("Unable to open output file " + outputFilename);
                return;
            }

            foreach (var sample in samplesToRun)
            {
                outputFile.WriteLine(sample);
                outputFile.WriteLine("Seed size\tRun time (days)\tmean F1\tSNV F1\tIndel F1\tSNV Recall\tSNV precision\tIndel Recall\tIndel precision");

                for (int seedSize = 16; seedSize <= 32; seedSize++)
                {
                    var timingFilename = timingDirectory + sample + ".snap_s" + seedSize + "_timings.txt";
                    var concordanceFilename = concordanceDirectory + sample + ".snap_s" + seedSize + ".dragen3.8VC.concordance.tar";
//Console.WriteLine(timingFilename + " " + File.Exists(timingFilename) + " " + concordanceFilename + File.Exists(concordanceFilename));

                    if (File.Exists(timingFilename) && File.Exists(concordanceFilename))
                    {
                        var timing = ASETools.RunTiming.LoadFromSNAPFile(timingFilename);
                        var concordance = new ASETools.ConcordanceResults(concordanceFilename);

                        outputFile.WriteLine(seedSize + "\t" + ((double)timing.alignTime + timing.loadingTime) / 3600 / 24 + "\t" + concordance.Mean_F1_Score() + "\t" +
                            concordance.results[ASETools.VariantType.SNV].F1_score + "\t" + concordance.results[ASETools.VariantType.Indel].F1_score + "\t" +
                            concordance.results[ASETools.VariantType.SNV].recall + "\t" + concordance.results[ASETools.VariantType.SNV].precision + "\t" +
                            concordance.results[ASETools.VariantType.Indel].recall + "\t" + concordance.results[ASETools.VariantType.Indel].precision + "\t");
                    }
                } // seed size
                outputFile.WriteLine();
            } // samples

            outputFile.Close();

        } // Main
    } // Program
} // namespace
