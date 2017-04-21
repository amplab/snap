using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using System.IO;
using ASELib;
using System.Diagnostics;

namespace ExtractMAFLines
{
    class Program
    {
        static void Main(string[] args)
        {
            var stopwatch = new Stopwatch();
            stopwatch.Start();

            var configuration = ASETools.ASEConfirguation.loadFromFile(args);
            
            if (configuration.commandLineArgs.Count() != 0) {
                Console.WriteLine("usage: ExtractMAFLines {-configuration configurationFileName}");
                return;
            }

                
            var cases = ASETools.Case.LoadCases(configuration.casesFilePathname);


            var casesByMAF = new Dictionary<string, List<ASETools.Case>>();

            int nCasesToProcess = 0;

            foreach (var caseEntry in cases)
            {
                var case_ = caseEntry.Value;

                if (case_.extracted_maf_lines != null && case_.extracted_maf_lines != "" || case_.maf_filename == null || case_.maf_filename == "")
                {
                    continue;
                }

                if (!casesByMAF.ContainsKey(case_.maf_filename))
                {
                    casesByMAF.Add(case_.maf_filename, new List<ASETools.Case>());
                }

                casesByMAF[case_.maf_filename].Add(case_);
                nCasesToProcess++;
            }

            Console.WriteLine("Processing " + nCasesToProcess + " cases with " + casesByMAF.Count() + " MAF files.");



            foreach (var casesForThisMAFEntry in casesByMAF)
            {
                var casesForThisMAF = casesForThisMAFEntry.Value;

                var mafFilename = casesForThisMAF[0].maf_filename;

                Console.Write("Loading MAF " + mafFilename + "...");

                var loadMAFTimer = new Stopwatch();
                loadMAFTimer.Start();

                var mafLines = ASETools.MAFLine.ReadFile(casesForThisMAF[0].maf_filename, casesForThisMAF[0].maf_file_id);

                Console.WriteLine(ASETools.ElapsedTimeInSeconds(loadMAFTimer));

                var writeMAFsTimer = new Stopwatch();
                writeMAFsTimer.Start();
                Console.Write("Writing " + casesForThisMAF.Count() + "...");

                foreach (var case_ in casesForThisMAF)
                {
                    var caseDirectory = ASETools.GetDirectoryFromPathname(mafFilename) + @"\..\..\" + configuration.derivedFilesDirectory + @"\" + case_.case_id + @"\";

                    Directory.CreateDirectory(caseDirectory);  // This is a no-op if the directory already exists.

                    var selectedLines = mafLines.Where(x => case_.tumor_dna_file_id == x.tumor_bam_uuid).ToList();

                    if (selectedLines.Count() == 0)
                    {
                        Console.WriteLine("Found no MAF lines for case " + case_.case_id);
                        continue;
                    }

                    ASETools.MAFLine.WriteToFile(caseDirectory + case_.case_id + ASETools.extractedMAFLinesExtension, selectedLines);
                }

                Console.WriteLine(ASETools.ElapsedTimeInSeconds(writeMAFsTimer));
            } // foereach MAF

            Console.WriteLine("Processed " + nCasesToProcess + " cases in " + casesByMAF.Count() + " MAFs in " + ASETools.ElapsedTimeInSeconds(stopwatch));
        }
    }
}
