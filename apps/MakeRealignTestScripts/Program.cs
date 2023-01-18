using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using ASELib;
using System.IO;

namespace MakeRealignTestScripts
{
    internal class Program
    {
        static ASETools.CommonData commonData;
        static void Main(string[] args)
        {
            //
            // Command args are file(s) containing normal DNA file IDs of already completed runs
            //

            commonData = ASETools.CommonData.LoadCommonData(args);
            if (null == commonData)
            {
                return;
            }

            var alreadyCompleted = new List<string>();

            if (commonData.configuration.commandLineArgs.Count() != 0)
            {
                foreach (var inputFilename in commonData.configuration.commandLineArgs)
                {
                    alreadyCompleted.AddRange(File.ReadAllLines(inputFilename));
                }
            }

            var outputFile = ASETools.CreateStreamWriterWithRetry("RealignTest.cmd");

            int nGenerated = 0;

            foreach (var case_ in commonData.listOfCases.Where(_ => _.normal_dna_filename != "" && _.case_pairedness_filename != ""))
            {
                if (alreadyCompleted.Contains(case_.normal_dna_file_id))
                {
                    continue;
                }

                nGenerated++;

                var pairedness = ASETools.CasePairedness.ReadFromFile(case_.case_pairedness_filename);
                
                outputFile.WriteLine(@"d:\gdc\bin\snap.exe " + (pairedness.pairedness[false][true] ? "paired" : "single")   // maps tumor->dna->pairedness
                   + @" \sequence\indices\Homo_sapiens_assembly38-22-large-liftover " + case_.normal_dna_filename + " -hc- -o delete-me." + case_.normal_dna_file_id + ".bam -d 20 -i 40 -so -sm 20" + (pairedness.pairedness[false][true] ? " -eh" : ""));              
            }

            outputFile.Close();

            Console.WriteLine("Generated script with " + nGenerated + " in " + ASETools.ElapsedTime(commonData.timer));
        } // Main
    } // Program
} // namespace
