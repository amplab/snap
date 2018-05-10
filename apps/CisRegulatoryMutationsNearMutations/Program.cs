using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using System.Diagnostics;
using System.IO;
using ASELib;

namespace CisRegulatoryMutationsNearMutations
{
    class Program
    {

        static ASETools.Configuration configuration;
        static List<ASETools.SelectedGene> selectedGenes;
        static ASETools.GeneLocationsByNameAndChromosome geneLocationInformation;
        static ASETools.GeneMap geneMap;

        static void Main(string[] args)
        {
            var timer = new Stopwatch();
            timer.Start();

            configuration = ASETools.Configuration.loadFromFile(args);
            if (null == configuration)
            {
                Console.WriteLine("Unable to read configuration");
                return;
            }

            var cases = ASETools.Case.LoadCases(configuration.casesFilePathname);
            if (null == cases)
            {
                Console.WriteLine("Unable to load cases");
                return;
            }

            if (configuration.commandLineArgs.Count() == 0 || configuration.commandLineArgs.Any(x => !cases.ContainsKey(x)))
            {
                Console.WriteLine("usage: AnaylzeCisRegulatoryRegions {caseId}");
                return;
            }

            var casesToProcess = cases.Select(x => x.Value).Where(x => configuration.commandLineArgs.Contains(x.case_id)).ToList();

            if (casesToProcess.Any(x => x.annotated_regulatory_regions_filename == "" || x.extracted_maf_lines_filename == ""))
            {
                Console.Write("Skipping cases because of missing input data:");
                casesToProcess.Where(x => x.annotated_regulatory_regions_filename == "" || x.extracted_maf_lines_filename == "").ToList().ForEach(x => Console.Write(" " + x.case_id));
                Console.WriteLine();
                casesToProcess = casesToProcess.Where(x => x.annotated_regulatory_regions_filename != "" && x.extracted_maf_lines_filename != "").ToList();

                if (casesToProcess.Count() == 0)
                {
                    return;
                }
            }

            geneLocationInformation = new ASETools.GeneLocationsByNameAndChromosome(ASETools.readKnownGeneFile(configuration.geneLocationInformationFilename));
            geneMap = new ASETools.GeneMap(geneLocationInformation.genesByName);

            selectedGenes = ASETools.SelectedGene.LoadFromFile(configuration.selectedGenesFilename).Where(x => geneLocationInformation.genesByName.ContainsKey(x.Hugo_Symbol)).ToList();
            if (null == selectedGenes)
            {
                Console.WriteLine("Unable to load selected genes.");
                return;
            }

            var threading = new ASETools.WorkerThreadHelper<ASETools.Case, int>(casesToProcess, HandleOneCase, null, null, 1);
            Console.Write("Processing " + casesToProcess.Count() + " cases, one dot/case: ");
            threading.run();

            Console.WriteLine(ASETools.ElapsedTimeInSeconds(timer));
        } // Main

        static void HandleOneCase(ASETools.Case case_, int state)
        {
            var mafLines = ASETools.MAFLine.ReadFile(case_.extracted_maf_lines_filename, case_.case_id, false);

            if (null == mafLines)
            {
                Console.WriteLine("Unable to read maf lines from " + case_.extracted_maf_lines_filename);
                return;
            }

            var annotatedBEDLines = ASETools.AnnotatedBEDLine.ReadFromFile(case_.annotated_regulatory_regions_filename);
            if (null == annotatedBEDLines)
            {
                Console.WriteLine("Unable to read BED lines from " + case_.annotated_regulatory_regions_filename);
                return;
            }


            string outputFilename = ASETools.GetDirectoryFromPathname(case_.annotated_regulatory_regions_filename) + @"\" + case_.case_id + ASETools.regulatoryMutationsNearMutationsExtension;
            var outputFile = ASETools.CreateStreamWriterWithRetry(outputFilename);
            if (null == outputFile)
            {
                Console.WriteLine("Unable to open output file " + outputFile);
                return;
            }

            outputFile.Write("Hugo_Symbol\tNon-Silent Mutation Count");
            for (int regionIndex = 1; regionIndex < ASETools.nRegions; regionIndex++)
            {
                outputFile.Write("\t" + ASETools.regionIndexToString(regionIndex));
            }
            outputFile.WriteLine();

            foreach (var selectedGene in selectedGenes)
            {
                outputFile.Write(selectedGene.Hugo_Symbol + "\t" + mafLines.Where(x => x.Hugo_Symbol == selectedGene.Hugo_Symbol && x.Variant_Classification != "Silent").Count());
                var geneLocation = geneLocationInformation.genesByName[selectedGene.Hugo_Symbol];

                var bedLinesThisChromosome = annotatedBEDLines.Where(x => geneLocation.chromosome == x.chrom).ToList();

                for (int regionIndex = 1; regionIndex < ASETools.nRegions; regionIndex++) {
                    List<ASETools.AnnotatedBEDLine> bedLinesInRange;

                    if (regionIndex == ASETools.nRegions - 1)
                    {
                        bedLinesInRange = annotatedBEDLines;
                    } else
                    {
                        bedLinesInRange = bedLinesThisChromosome.Where(x => geneLocation.minDistanceFromTSS(x.chrom, x.chromStart, x.chromEnd) <= ASETools.regionIndexToSizeInBases(regionIndex)).ToList();
                    }
    
                    int nBasesWithCoverage = bedLinesInRange.Select(x => x.nBasesWithCoverage).Sum();

                    if (nBasesWithCoverage < 10)
                    {
                        outputFile.Write("\t*");
                    } else {
                        int nCisRegulatoryMutations = bedLinesInRange.Select(x => x.nMutations).Sum();
                        outputFile.Write("\t" + (double)nCisRegulatoryMutations / nBasesWithCoverage);
                    }
                    
                } // for each region
                outputFile.WriteLine();
            }

            outputFile.WriteLine("**done**");
            outputFile.Close();
        } // HandleOneCase


    } // Program
} // namespace CisRegulatoryMutationsNearMutations
