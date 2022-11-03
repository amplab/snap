using System;
using System.CodeDom;
using System.Collections.Generic;
using System.IO;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using ASELib;

namespace CollectResults
{
    internal class Program
    {
        static string nonInfogainFilenamePrefix = "mason_hg3_";
        static string nonInfogainFilenamePostfix = "x.concordance.tar";

        static int getCoverageFromPathname(string pathname)
        {
            var filename = ASETools.GetFileNameFromPathname(pathname);
            if (!filename.StartsWith(nonInfogainFilenamePrefix))
            {
                throw new Exception("getCoverageFromPathname: filename doesn't start correctly " + pathname);
            }

            return ASETools.GetNextNumberFromString(filename.Substring(nonInfogainFilenamePrefix.Length));
        } // getCoverageFromPathname

        static void Main(string[] args)
        {
            //
            // First do the non-Infogain baseline.
            //
            var nonInfogainResults = new Dictionary<int, ASETools.ConcordanceResults>(); // Maps coverage to concordance

            foreach (var giabMachine in ASETools.GIABMachines)
            {
                foreach (var pathname in Directory.EnumerateFiles(@"\\" + giabMachine + @"\d$\temp\", nonInfogainFilenamePrefix + "*" + nonInfogainFilenamePostfix))
                {
                    var coverage = getCoverageFromPathname(pathname);

                    var concordance = new ASETools.ConcordanceResults(pathname);

                    if (nonInfogainResults.ContainsKey(coverage))
                    {
                        Console.WriteLine("Coverage " + coverage + " appears twice for non-Infogain.  Second time at " + pathname);
                    } else
                    {
                        nonInfogainResults.Add(coverage, concordance);
                    }                    
                } // foreach concordance file on this machine
            } // foreach machine

            Console.WriteLine("Found " + nonInfogainResults.Count + " non-Infogain results");

            var outputFile = ASETools.CreateStreamWriterWithRetry(@"\temp\InfoGainResults.txt");

            outputFile.WriteLine("Coverage\tMean F1\tSNV F1\tIndel F1\tSNV precision\tSNV recall\tIndel precision\tIndel recall");
            for (int i = 1; i < nonInfogainResults.Select(_ => _.Key).Max(); i++)
            {
                if (!nonInfogainResults.ContainsKey(i))
                {
                    continue;
                }

                outputFile.WriteLine(i + "\t" + nonInfogainResults[i].Mean_F1_Score() + "\t" + nonInfogainResults[i].results[ASETools.VariantType.SNV].F1_score + "\t" + nonInfogainResults[i].results[ASETools.VariantType.Indel].F1_score + "\t" +
                                     nonInfogainResults[i].results[ASETools.VariantType.SNV].precision + "\t" + nonInfogainResults[i].results[ASETools.VariantType.SNV].recall + "\t" +
                                     nonInfogainResults[i].results[ASETools.VariantType.Indel].precision + "\t" + nonInfogainResults[i].results[ASETools.VariantType.Indel].recall);
            } // foreach potential coverage

            outputFile.WriteLine("**done**");
            outputFile.Close();

        } // Main
    } // Program
} // namespace
