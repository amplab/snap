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
        static List<ASETools.GeneHancerLine> geneHancers;


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

            regulatoryRegions = ASETools.BEDFile.ReadFromFile(configuration.encodeBEDFile, !configuration.usesChr);
            if (null == regulatoryRegions)
            {
                Console.WriteLine("Unable to read regulatory regions BED file.");
                return;
            }

            geneHancers = ASETools.GeneHancerLine.ReadFromFile(configuration.geneHancerFilename);
            if (geneHancers == null)
            {
                Console.WriteLine("Unable to read gene hancers.");
                return;
            }


            if (configuration.commandLineArgs.Count() == 0 || configuration.commandLineArgs.Any(x => !cases.ContainsKey(x)))
            {
                Console.WriteLine("usage: AnnotateRegulatoryRegions <caseId>");
                return;
            }

            var casesToProcess = configuration.commandLineArgs.Select(x => cases[x]).ToList();
            Console.WriteLine("Processing " + casesToProcess.Count() + " cases, 1 dot/case:");
            ASETools.PrintNumberBar(casesToProcess.Count());

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

            var annotatedGeneHancerLines = new ASETools.AVLTree<ASETools.AnnotatedGeneHancerLine>();
            geneHancers.ForEach(x => annotatedGeneHancerLines.Insert(new ASETools.AnnotatedGeneHancerLine(x)));

            allcountReader.ReadAllcountFile((string contigName, int location, int currentMappedReadCount) => ProcessOneBase(contigName, location, currentMappedReadCount, annotatedBEDLines, annotatedGeneHancerLines));

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

                HandleOneLine(mafLines, line.chrom, line.chromStart, line.chromEnd, out line.nMutations, out line.nMutationsAbove60Percent, out line.nMutationsBelow40Percent);

                outputFile.WriteLine(line.chrom + "\t" + line.chromStart + "\t" + line.chromEnd + "\t" + line.name + "\t" + line.score + "\t" + line.strand + "\t" + line.thickStart + "\t" + line.thickEnd + "\t" +
                    line.itemRGB + "\t" + line.nBasesWithCoverage + "\t" + line.nMutations + "\t" + line.nMutationsBelow40Percent + "\t" + line.nMutationsAbove60Percent);
            }

            outputFile.WriteLine("**done**");
            outputFile.Close();

            outputFilename = ASETools.GetDirectoryFromPathname(case_.tumor_dna_allcount_filename) + @"\" + case_.case_id + ASETools.annotatedGeneHancerLinesExtension;
            outputFile = ASETools.CreateStreamWriterWithRetry(outputFilename);
            if (null == outputFile)
            {
                Console.WriteLine("Unable to open output file " + outputFilename);
                return;
            }

            outputFile.WriteLine("Chromosome\tStart\tEnd\tFeature\tScore\tAttributes\tn Bases With Coverage\tn Mutations\tn Mutations Below 40 Percent\tn Mutations Above 60 Percent");
            for (var line = annotatedGeneHancerLines.min(); null != line; line = annotatedGeneHancerLines.FindFirstGreaterThan(line))
            {
                if (line.nBasesWithCoverage == 0)
                {
                    continue;
                }

                HandleOneLine(mafLines, line.chromosome, line.start, line.end, out line.nMutations, out line.nMutationsAbove60Percent, out line.nMutationsBelow40Percent);

                outputFile.WriteLine(line.chromosome + "\t" + line.start + "\t" + line.end + "\t" + line.feature + "\t" + line.score
                    + "\t" + line.attributes + "\t" + line.nBasesWithCoverage + "\t" + line.nMutations + "\t" +
                    line.nMutationsBelow40Percent + "\t" + line.nMutationsAbove60Percent);
            }

            outputFile.WriteLine("**done**");
            outputFile.Close();
        } // HandleOneCase

        static void HandleOneLine(List<ASETools.MAFLine> mafLines, string chrom, int start, int end, out int nMutations,
            out int nMutationsAbove60Percent, out int nMutationsBelow40Percent)
        {
            var mutations = mafLines.Where(x => x.Chromosome == chrom && x.Start_Position <= end && x.End_Positon >= start).ToList();
            nMutations = mutations.Count();
            nMutationsAbove60Percent = mutations.Where(x => (double)x.n_alt_count / (x.n_alt_count + x.n_ref_count) > 0.6).Count();
            nMutationsBelow40Percent = mutations.Where(x => (double)x.n_alt_count / (x.n_alt_count + x.n_ref_count) < 0.4).Count();
        }

        static void ProcessOneBase(string contigName, int location, int currentMappedReadCount, ASETools.AVLTree<ASETools.AnnotatedBEDLine> annotatedBEDLines, 
            ASETools.AVLTree<ASETools.AnnotatedGeneHancerLine> annotatedGeneHancerLines)
        {
            if (currentMappedReadCount < configuration.minDNAReadCoverage)
            {
                return;
            }

            ASETools.AnnotatedBEDLine bedLine;
            ASETools.AnnotatedBEDLine bedKey = new ASETools.AnnotatedBEDLine(new ASETools.BEDLine(contigName, location));

            if (annotatedBEDLines.FindFirstLessThanOrEqualTo(bedKey, out bedLine) && bedLine.chrom == contigName && bedLine.chromEnd < location) {
                bedLine.nBasesWithCoverage++;
            }

            ASETools.AnnotatedGeneHancerLine geneHancerLine;
            ASETools.AnnotatedGeneHancerLine geneHancerKey = new ASETools.AnnotatedGeneHancerLine(new ASETools.GeneHancerLine(contigName, location));

            if (annotatedGeneHancerLines.FindFirstLessThanOrEqualTo(geneHancerKey, out geneHancerLine) && geneHancerLine.chromosome == contigName &&
                geneHancerLine.start <= location && geneHancerLine.end >= location)
            {
                geneHancerLine.nBasesWithCoverage++;
            }
        } // ProcessOneBase
    }
}
