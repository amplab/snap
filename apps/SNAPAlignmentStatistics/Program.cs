using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using ASELib;
using System.IO;

namespace SNAPAlignmentStatistics
{
    class Program
    {
        static void Main(string[] args)
        {
            string[] samples = { "ERR194146", "ERR194147", "hg001", "hg002", "hg003", "hg004", "hg005", "hg006", "hg007", };

            var outputFilename = @"d:\gdc\final_results\snap_alignment_statictics.txt";
            var outputFile = ASETools.CreateStreamWriterWithRetry(outputFilename);

            if (outputFile == null)
            {
                Console.WriteLine("Unable to open output file " + outputFilename);
                return;
            }

            outputFile.WriteLine("Sample\ttotal reads\tAligned MAPQ >= 10\tAligned, MAPQ < 10\tUnaligned\tToo Short or too many Ns\t%Pairs\tTime In Aligner (s)\t% Aligned Reads aligned as single");

            foreach (var sample in samples)
            {
                var timingFilename = @"d:\temp\timings\" + sample + ".snap_timings.txt";

                if (!File.Exists(timingFilename))
                {
                    Console.WriteLine("Timing file " + timingFilename + " doesn't exist");
                    continue;
                }

                var timing = ASETools.RunTiming.LoadFromSNAPFile(timingFilename);
                outputFile.WriteLine(sample + "\t" + ASETools.NumberWithCommas(timing.totalReads) + "\t" + ASETools.NumberWithCommas(timing.alignedMapQAtLeast10) + "\t" +
                    ASETools.NumberWithCommas(timing.alignedMapQLessThan10) + "\t" + ASETools.NumberWithCommas(timing.unaligned) + "\t" +
                    ASETools.NumberWithCommas(timing.tooShortTooManyNs) + "\t" + timing.fractionPairs * 100.0 + "%\t" + ASETools.NumberWithCommas(timing.overallRuntime) + "\t" +
                    (1.0 - (timing.fractionPairs * timing.totalReads / (timing.alignedMapQAtLeast10 + timing.alignedMapQLessThan10))) * 100.0 + "%");
            } // sample

            outputFile.WriteLine("**done**");
            outputFile.Close();
        } // Main
    } // Program
} // namespace
