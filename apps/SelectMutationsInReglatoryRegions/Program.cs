using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using System.IO;
using System.Diagnostics;
using ASELib;

namespace SelectMutationsInReglatoryRegions
{
    class Program
    {
        static ASETools.Configuration configuration;
        static ASETools.BEDFile bedFile;
        static Dictionary<string, ASETools.Case> cases;

        static void Main(string[] args)
        {
            Stopwatch timer = new Stopwatch();
            timer.Start();

            configuration = ASETools.Configuration.loadFromFile(args);

            if (configuration == null)
            {
                Console.WriteLine("Unable to load configuration");
                return;
            }

            cases = ASETools.Case.LoadCases(configuration.casesFilePathname);
            if (null == cases)
            {
                Console.WriteLine("Unable to load cases");
                return;
            }

            if (configuration.commandLineArgs.Count() == 0 || configuration.commandLineArgs.Any(x => !cases.ContainsKey(x)))
            {
                Console.WriteLine("usage: SelectMutationsInRegulatoryRegions <cases>");
                return;
            }

            bedFile = ASETools.BEDFile.ReadFromFile(configuration.encodeBEDFile, !configuration.usesChr);

            var casesToProcess = configuration.commandLineArgs.Select(x => cases[x]).ToList();

            var casesWithoutData = casesToProcess.Where(x => x.all_maf_lines_filename == "").ToList();
            if (casesWithoutData.Count() != 0)
            {
                Console.Write("Cases");
                casesWithoutData.ForEach(x => Console.Write(" " + x.case_id));
                Console.WriteLine(" do not have all maf line files.  Ignoring.");
            }

            var threading = new ASETools.WorkerThreadHelper<ASETools.Case, int>(casesToProcess, HandleOneCase, null, null, 1);

            Console.Write("Processing " + casesToProcess.Count() + " cases, one dot/case: ");

            threading.run();

            Console.WriteLine(ASETools.ElapsedTimeInSeconds(timer));
        }

        static void HandleOneCase(ASETools.Case case_, int state)
        {
            if (case_.all_maf_lines_filename == "")
            {
                return;
            }

            var mafLines = ASETools.MAFLine.ReadFile(case_.all_maf_lines_filename, case_.maf_file_id, false);

            if (mafLines == null)
            {
                Console.WriteLine("Unable to read MAF lines from " + case_.all_maf_lines_filename);
                return;
            }

            //
            // Find the MAF lines that are in regulatory regions.
            //
            var selectedMAFLines = mafLines.Where(x => (x.n_alt_count + x.n_ref_count >= configuration.minDNAReadCoverage) && (bedFile.isInRegulatoryRegion(x.Chromosome, x.Start_Position) || bedFile.isInRegulatoryRegion(x.Chromosome, x.End_Positon))).ToList();

            if (!case_.all_maf_lines_filename.EndsWith(ASETools.allMAFLinesExtension))
            {
                Console.WriteLine("all maf lines file " + case_.all_maf_lines_filename + " doesn't end with the right extention.");
                return;
            }

            var outputFilename = case_.all_maf_lines_filename.Substring(0, case_.all_maf_lines_filename.Length - ASETools.allMAFLinesExtension.Length) + ASETools.selectedRegulatoryMAFLinesExtension;

            ASETools.MAFLine.WriteToFile(outputFilename, selectedMAFLines);
        }
    }
}
