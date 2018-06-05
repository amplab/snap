using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using ASELib;
using System.Diagnostics;

namespace CisRegulatoryMutationsInKnownRegions
{
    class Program
    {
        static ASETools.Configuration configuration;
        static Dictionary<string, List<ASETools.GeneHancerLine>> geneHancersByGene;
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
                Console.WriteLine("Unable to load configuration.");
                return;
            }

            if (configuration.commandLineArgs.Count() == 0)
            {
                Console.WriteLine("Usage: CisRegulatoryMutationsInKnownRegions <caseId>");
            }

            var cases = ASETools.Case.LoadCases(configuration.casesFilePathname);
            if (null == cases)
            {
                Console.WriteLine("Unable to load cases.");
                return;
            }

            geneHancersByGene = ASETools.GeneHancerLine.ReadFromFileToDict(configuration.geneHancerFilename);
            if (geneHancersByGene == null)
            {
                Console.WriteLine("Unable to read gene hancers.");
                return;
            }

            perGeneASEMap = ASETools.ASEMapPerGeneLine.ReadFromFileToDictionary(configuration.finalResultsDirectory + ASETools.PerGeneASEMapFilename);
            if (perGeneASEMap == null)
            {
                Console.WriteLine("Unable to load per-Gene ASE map from " + configuration.finalResultsDirectory + ASETools.PerGeneASEMapFilename);
                return;
            }

            geneLocationInformation = new ASETools.GeneLocationsByNameAndChromosome(ASETools.readKnownGeneFile(configuration.geneLocationInformationFilename));
            geneMap = new ASETools.GeneMap(geneLocationInformation.genesByName);

            var casesToProcess = cases.Select(x => x.Value).Where(x => configuration.commandLineArgs.Contains(x.case_id)).ToList();

            if (casesToProcess.Count() != configuration.commandLineArgs.Count())
            {
                Console.Write("The following args don't appear to be case IDs, ignoring:");
                configuration.commandLineArgs.Where(x => !cases.Select(y => y.Value.case_id).Contains(x)).ToList().ForEach(x => Console.Write(x + " "));
            }

            Console.WriteLine("Processing " + casesToProcess.Count() + " cases, one dot/case: ");
            ASETools.PrintNumberBar(casesToProcess.Count());

            var threading = new ASETools.WorkerThreadHelper<ASETools.Case, int>(casesToProcess, HandleOneCase, null, null, 1);
            threading.run();

            Console.WriteLine();
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

            var mafLines = ASETools.MAFLine.ReadFile(case_.all_maf_lines_filename, "", false);
            if (null == mafLines)
            {
                Console.WriteLine("Unable to load MAF lines for case " + case_.case_id + " from file " + case_.all_maf_lines_filename);
                return;
            }

            var copyNumber = ASETools.CopyNumberVariation.ReadBothFiles(case_);
            // It's OK for this to be null, a few didn't have them in TCGA.

            var outputFilename = ASETools.GetDirectoryFromPathname(case_.all_maf_lines_filename) + @"\" + case_.case_id + ASETools.cisRegulatoryMutationsInKnownRegionsExtension;
            var outputFile = ASETools.CreateStreamWriterWithRetry(outputFilename);
            if (null == outputFile)
            {
                Console.WriteLine("Unable to open output file " + outputFilename);
                return;
            }

            outputFile.WriteLine("Hugo Symbol\tMutation Count\tn Silent\tn Not ASE Candidates\tn < 0.4 VAF\tn 0.4 <= VAF <= 0.6\tn > 0.6 VAF\tList of GeneHancers and Mutation Counts for this gene");
            foreach (var geneEntry in geneHancersByGene)
            {
                var hugo_symbol = geneEntry.Key;
                var geneHancers = geneEntry.Value;

                var annotatedSelectedVariantsThisGene = annotatedSelectedVariants.Where(x => x.somaticMutation && x.Hugo_symbol == hugo_symbol).ToList();
                var nonSlientASVThisGene = annotatedSelectedVariantsThisGene.Where(x => !x.isSilent());
                var nonSilentASECandidates = nonSlientASVThisGene.Where(x => x.IsASECandidate(true, copyNumber, configuration, perGeneASEMap, geneMap));

                outputFile.Write(hugo_symbol + "\t"  + annotatedSelectedVariantsThisGene.Count() + "\t" + (annotatedSelectedVariantsThisGene.Count() - nonSlientASVThisGene.Count()) + "\t" + 
                    nonSlientASVThisGene.Where(x => x.V)
            }


            outputFile.Close();
        } // HandleOneCase

    }
}
