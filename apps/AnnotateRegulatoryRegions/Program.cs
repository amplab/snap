using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using System.Diagnostics;
using System.IO;
using ASELib;

namespace AnnotateRegulatoryRegions
{
    class Program
    {

        static ASETools.Configuration configuration;
        static Dictionary<string, ASETools.Case> cases;
        static ASETools.BEDFile regulatoryRegions;

        static void Main(string[] args)
        {
            var timer = new Stopwatch();
            timer.Start();

            configuration = ASETools.Configuration.loadFromFile(args);
            if (null == configuration)
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

            regulatoryRegions = ASETools.BEDFile.ReadFromFile(configuration.encodeBEDFile);
            if (null == regulatoryRegions)
            {
                Console.WriteLine("Unable to read regulatory regions BED file.");
                return;
            }

            if (configuration.commandLineArgs.Count() == 0 || configuration.commandLineArgs.Any(x => !cases.ContainsKey(x)))
            {
                Console.WriteLine("usage: AnnotateRegulatoryRegions <caseId>");
                return;
            }

            var casesToProcess = configuration.commandLineArgs.Select(x => cases[x]).ToList();
            Console.Write("Processing " + casesToProcess.Count() + " cases, 1 dot/case:");

            var threading = new ASETools.WorkerThreadHelper<ASETools.Case, int>(casesToProcess, HandleOneCase, null, null, 1);
            threading.run();

            Console.WriteLine(ASETools.ElapsedTimeInSeconds(timer));
        } // Main

        static void HandleOneCase(ASETools.Case case_, int state)
        {
            if (case_.tumor_dna_allcount_filename == "")
            {
                Console.WriteLine("Case " + case_.case_id + " does not have a tumor DNA allcount file.  Skipping.");
                return;
            }

            var allcountReader = new ASETools.AllcountReader(case_.tumor_dna_allcount_filename);
            if (!allcountReader.openFile())
            {
                Console.WriteLine("Failed to open allcount file for case " + case_.case_id + ".  SKipping.");
                return;
            }

            if (case_.selected_regulatory_maf_filename == "")
            {
                Console.WriteLine("Case " + case_.case_id + " doesn't have a selected regulatory maf lines file.  Skipping.");
                return;
            }

            var mafLines = ASETools.MAFLine.ReadFile(case_.selected_regulatory_maf_filename, case_.case_id, false);
            if (null == mafLines)
            {
                Console.WriteLine("Unable to read MAF lines for case " + case_.case_id);
                return;
            }

            var mafTree = new ASETools.AVLTree<ASETools.MAFLine>();
            foreach (var mafLine in mafLines)
            {
                if (!mafTree.ContainsKey(mafLine))
                {
                    mafTree.Insert(mafLine);
                }
            }

            var annotatedBEDLines = new ASETools.AVLTree<ASETools.AnnotatedBEDLine>();
            regulatoryRegions.ForEachLine(x => annotatedBEDLines.Insert(new ASETools.AnnotatedBEDLine(x)));

            allcountReader.ReadAllcountFile((string contigName, int location, int currentMappedReadCount) => ProcessOneBase(contigName, location, currentMappedReadCount, annotatedBEDLines));

            var outputFilename = ASETools.GetDirectoryFromPathname(case_.tumor_dna_allcount_filename) + @"\" + case_.case_id + ASETools.annotatedBEDLinesExtension;
            var outputFile = ASETools.CreateStreamWriterWithRetry(outputFilename);
            if (null == outputFile)
            {
                Console.WriteLine("Unable to open output file " + outputFilename);
                return;
            }

            outputFile.WriteLine("Chrom\tChromStart\tChromEnd\tname\tscore\tstrand\tthickStart\tthickEnd\titemRGB\tnBasesWithCoverage\tnMutations\tnMutationsBelow40Percent\tnMutationsAbove60Percent");
            for (var line = annotatedBEDLines.min(); null != line; line = annotatedBEDLines.FindFirstGreaterThan(line)) // I just can't get IEnumerable to work for AVLTree, so no foreach()
            {
                if (line.nBasesWithCoverage == 0)
                {
                    continue;   // Skip these, they're no coverage so no use
                }

                var mutations = mafLines.Where(x => x.Chromosome == line.chrom && x.Start_Position <= line.chromEnd && x.End_Positon >= line.chromStart).ToList();
                line.nMutations = mutations.Count();
                line.nMutationsAbove60Percent = mutations.Where(x => (double)x.n_alt_count / (x.n_alt_count + x.n_ref_count) > 0.6).Count();
                line.nMutationsBelow40Percent = mutations.Where(x => (double)x.n_alt_count / (x.n_alt_count + x.n_ref_count) < 0.4).Count();

                outputFile.WriteLine(line.chrom + "\t" + line.chromStart + "\t" + line.chromEnd + "\t" + line.name + "\t" + line.score + "\t" + line.strand + "\t" + line.thickStart + "\t" + line.thickEnd + "\t" +
                    line.itemRGB + "\t" + line.nBasesWithCoverage + "\t" + line.nMutations + "\t" + line.nMutationsBelow40Percent + "\t" + line.nMutationsAbove60Percent);
            }

            outputFile.WriteLine("**done**");
            outputFile.Close();

        } // HandleOneCase

        static void ProcessOneBase(string contigName, int location, int currentMappedReadCount, ASETools.AVLTree<ASETools.AnnotatedBEDLine> annotatedBEDLines)
        {
            if (currentMappedReadCount < configuration.minDNAReadCoverage)
            {
                return;
            }

            ASETools.AnnotatedBEDLine line;
            ASETools.AnnotatedBEDLine key = new ASETools.AnnotatedBEDLine(new ASETools.BEDLine(contigName, location));

            if (!annotatedBEDLines.FindFirstLessThanOrEqualTo(key, out line) || line.chrom != contigName || line.chromEnd < location) {
                return;
            }

            line.nBasesWithCoverage++;
        }
    }
}
