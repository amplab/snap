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
        static Dictionary<string, ASETools.ASEMapPerGeneLine> perGeneASEMap;

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

            perGeneASEMap = ASETools.ASEMapPerGeneLine.ReadFromFileToDictionary(configuration.finalResultsDirectory + ASETools.PerGeneASEMapFilename);
            if (perGeneASEMap == null)
            {
                Console.WriteLine("Unable to load per-Gene ASE map from " + configuration.finalResultsDirectory + ASETools.PerGeneASEMapFilename);
                return;
            }

            var casesToProcess = cases.Select(x => x.Value).Where(x => configuration.commandLineArgs.Contains(x.case_id)).ToList();

            if (casesToProcess.Any(x => x.annotated_regulatory_regions_filename == "" || x.annotated_selected_variants_filename == ""))
            {
                Console.Write("Skipping cases because of missing input data:");
                casesToProcess.Where(x => x.annotated_regulatory_regions_filename == "" || x.annotated_selected_variants_filename == "").ToList().ForEach(x => Console.Write(" " + x.case_id));
                Console.WriteLine();
                casesToProcess = casesToProcess.Where(x => x.annotated_regulatory_regions_filename != "" && x.annotated_selected_variants_filename != "").ToList();

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
            Console.WriteLine("Processing " + casesToProcess.Count() + " cases, one dot/case");
            ASETools.PrintNumberBar(casesToProcess.Count());
            threading.run();

            Console.WriteLine(ASETools.ElapsedTimeInSeconds(timer));
        } // Main

        static void HandleOneCase(ASETools.Case case_, int state)
        {

            var annotatedSelectedVariants = ASETools.AnnotatedVariant.readFile(case_.annotated_selected_variants_filename);

            if (null == annotatedSelectedVariants)
            {
                Console.WriteLine("Unable to read annotated selected variants from " + case_.annotated_selected_variants_filename);
                return;
            }

            var copyNumber = ASETools.CopyNumberVariation.ReadBothFiles(case_);
            // It's OK for this to be null, a few didn't have them in TCGA.

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

            outputFile.Write("Hugo_Symbol\tNon-Silent Mutation Count\tNon - Silent Mutations with Low VAF\tNon - Silent Mutations with Moderate VAF\tNon - Silent Mutations with High VAF\tSilent Mutations\tNot ASE Candidates");
            for (int regionIndex = 1; regionIndex < ASETools.nRegions; regionIndex++)
            {
                outputFile.Write("\t" + ASETools.regionIndexToString(regionIndex));
            }
            outputFile.WriteLine();

            var variantsByGene = annotatedSelectedVariants.GroupByToDict(x => x.Hugo_symbol);
            var bedLinesByChromosome = annotatedBEDLines.GroupByToDict(x => x.chrom);

            foreach (var selectedGene in selectedGenes)
            {
                List<ASETools.AnnotatedVariant> variantsThisGene;
                if (variantsByGene.ContainsKey(selectedGene.Hugo_Symbol))
                {
                    variantsThisGene = variantsByGene[selectedGene.Hugo_Symbol];
                } else
                {
                    variantsThisGene = new List<ASETools.AnnotatedVariant>();
                }
 
                var nonSilentVariantsThisGene = variantsThisGene.Where(x => !x.isSilent() && x.IsASECandidate(true, copyNumber, configuration, perGeneASEMap, geneMap)).ToList();
                outputFile.Write(
                    selectedGene.Hugo_Symbol + "\t" +
                    nonSilentVariantsThisGene.Count() + "\t" +
                    nonSilentVariantsThisGene.Where(x => x.GetTumorAltAlleleFraction() < 0.4).Count() + "\t" +
                    nonSilentVariantsThisGene.Where(x => x.GetTumorAltAlleleFraction() >= 0.4 && x.GetTumorAltAlleleFraction() <= 0.6).Count() + "\t" +
                    nonSilentVariantsThisGene.Where(x => x.GetTumorAltAlleleFraction() > 0.6).Count() + "\t" +
                    variantsThisGene.Where(x => x.isSilent()).Count() + "\t" +
                    variantsThisGene.Where(x => !x.isSilent() && !x.IsASECandidate(true, copyNumber, configuration, perGeneASEMap, geneMap)).Count());

                var geneLocation = geneLocationInformation.genesByName[selectedGene.Hugo_Symbol];

                List<ASETools.AnnotatedBEDLine> bedLinesThisChromosome;
                if (bedLinesByChromosome.ContainsKey(geneLocation.chromosome))
                {
                    bedLinesThisChromosome = bedLinesByChromosome[geneLocation.chromosome];
                } else
                {
                    bedLinesThisChromosome = new List<ASETools.AnnotatedBEDLine>();
                }

                int nBedLinesThisChromosome = bedLinesThisChromosome.Count();

                for (int regionIndex = 1; regionIndex < ASETools.nRegions; regionIndex++) {
                    List<ASETools.AnnotatedBEDLine> bedLinesInRange;

                    if (regionIndex == ASETools.nRegions - 1)
                    {
                        bedLinesInRange = annotatedBEDLines;
                    } else
                    {
                        bedLinesInRange = new List<ASETools.AnnotatedBEDLine>();

                        for (int i = 0; i < nBedLinesThisChromosome; i++)
                        {
                            var bedLine = bedLinesThisChromosome[i];
                            if (geneLocation.minDistanceFromTSS(bedLine.chrom, bedLine.chromStart, bedLine.chromEnd) < ASETools.regionIndexToSizeInBases(regionIndex))
                            {
                                bedLinesInRange.Add(bedLine);
                            }
                        }
                        // too slow bedLinesInRange = bedLinesThisChromosome.Where(x => geneLocation.minDistanceFromTSS(x.chrom, x.chromStart, x.chromEnd) <= ASETools.regionIndexToSizeInBases(regionIndex)).ToList();
                    }

                    int nBasesWithCoverage = 0;
                    int nBedLinesInRange = bedLinesInRange.Count();
                    int nCisRegulatoryMutations = 0;

                    for (int i = 0; i < nBedLinesInRange; i++)
                    {
                        nBasesWithCoverage += bedLinesInRange[i].nBasesWithCoverage;
                        nCisRegulatoryMutations += bedLinesInRange[i].nMutationsWithModerateVAF();
                    }

                    // unrolled into for loop above for speed    bedLinesInRange.Select(x => x.nBasesWithCoverage).Sum();

                    if (nBasesWithCoverage < 10)
                    {
                        outputFile.Write("\t*");
                    } else {
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
