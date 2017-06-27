using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using System.Threading;
using System.IO;
using ASELib;

namespace CountVariantsInGenes
{
    class Program
    {

        static void Main(string[] args)
        {
            var configuration = ASETools.ASEConfirguation.loadFromFile(args);

            if (null == configuration)
            {
                Console.WriteLine("Unable to load configuration.");
                return;
            }

            if (configuration.commandLineArgs.Count() == 0)
            {
                Console.WriteLine("usage: CountVariantsInGenes <one or more case IDs>");
                return;
            }

            var cases = ASETools.Case.LoadCases(configuration.casesFilePathname);

            if (null == cases)
            {
                Console.WriteLine("Unable to load cases.");
                return;
            }

            List<string> caseIdsToProcess = new List<string>();

            foreach (var caseId in configuration.commandLineArgs)
            {
                if (!cases.ContainsKey(caseId))
                {
                    Console.WriteLine(caseId + " does not appear to be a valid case ID");
                } else
                {
                    caseIdsToProcess.Add(caseId);
                }
            }

            var selectedGenes = ASETools.SelectedGene.LoadFromFile(configuration.selectedGenesFilename);
            if (selectedGenes == null)
            {
                Console.WriteLine("Unable to load selected genes from " + configuration.selectedGenesFilename);
                return;
            }

            var geneLocationInformation = new ASETools.GeneLocationsByNameAndChromosome(ASETools.readKnownGeneFile(ASETools.ASEConfirguation.defaultGeneLocationInformation));
            var selectedGenesWithLocationInformation = selectedGenes.Where(x => geneLocationInformation.genesByName.ContainsKey(x.Hugo_Symbol)).ToList();

            var threads = new List<Thread>();
            for (int i = 0; i < Environment.ProcessorCount; i++)
            {
                threads.Add(new Thread(() => WorkerThread(caseIdsToProcess, cases, selectedGenes, geneLocationInformation)));
            }

            threads.ForEach(t => t.Start());
            threads.ForEach(t => t.Join());

        } // Main

        static void WorkerThread(List<string> caseIdsToProcess, Dictionary<string, ASETools.Case> cases, List<ASETools.SelectedGene> selectedGenesWithLocationInformation, ASETools.GeneLocationsByNameAndChromosome geneLocationInformation)
        {
            while (true)
            {
                ASETools.Case case_;

                lock (caseIdsToProcess)
                {
                    if (caseIdsToProcess.Count() == 0)
                    {
                        return;
                    }

                    case_ = cases[caseIdsToProcess[0]];
                    caseIdsToProcess.RemoveAt(0);
                }

                if (case_.annotated_selected_variants_filename == "")
                {
                    Console.WriteLine("Case " + case_.case_id + " does not have annotated selected variants.  Skipping.");
                    continue;
                }

                var annotatedSelectedVariants = ASETools.AnnotatedVariant.readFile(case_.annotated_selected_variants_filename);
                if (null == annotatedSelectedVariants)
                {
                    Console.WriteLine("Unable to load annotated selected variants for case " + case_.case_id + " from " + case_.annotated_selected_variants_filename);
                    continue;
                }

                var countsByGene = new Dictionary<string, int>();
                selectedGenesWithLocationInformation.ForEach(x => countsByGene.Add(x.Hugo_Symbol, 0));

                int totalSelectedVariantsInGenes = 0;

                foreach (var selectedVariant in annotatedSelectedVariants)
                {
                    if (selectedVariant.somaticMutation)
                    {
                        continue;
                    }

                    foreach (var gene in selectedGenesWithLocationInformation)
                    {
                        if (geneLocationInformation.genesByName.ContainsKey(gene.Hugo_Symbol) && 
                            !geneLocationInformation.genesByName[gene.Hugo_Symbol].inconsistent &&
                            geneLocationInformation.genesByName[gene.Hugo_Symbol].chromosome == selectedVariant.contig && 
                            geneLocationInformation.genesByName[gene.Hugo_Symbol].maxLocus >= selectedVariant.locus && 
                            geneLocationInformation.genesByName[gene.Hugo_Symbol].minLocus <= selectedVariant.locus)
                        {
                            var foo = geneLocationInformation.genesByName[gene.Hugo_Symbol];
                            countsByGene[gene.Hugo_Symbol]++;
                            totalSelectedVariantsInGenes++;
                        }
                    } // for every gene
                } // for every selected variant.

                var outputFilename = ASETools.GetDirectoryFromPathname(case_.annotated_selected_variants_filename) + @"\" + case_.case_id + ASETools.selectedVariantCountByGeneExtension;
                var outputFile = ASETools.CreateStreamWriterWithRetry(outputFilename);
                if (null == outputFile)
                {
                    Console.WriteLine("Unable to open output file " + outputFilename);
                    continue;
                }

                outputFile.WriteLine("Hugo_Symbol\tnVariantsInGene\tnSomaticMutations");

                foreach (var geneEntry in countsByGene)
                {
                    outputFile.WriteLine(geneEntry.Key + "\t" + geneEntry.Value + "\t" + annotatedSelectedVariants.Where(x => x.Hugo_symbol == geneEntry.Key && x.somaticMutation).Count());
                }

                outputFile.WriteLine("**done**");
                outputFile.Close();

                Console.WriteLine("Case " + case_.case_id + " has " + totalSelectedVariantsInGenes + " of its " + annotatedSelectedVariants.Where(x => !x.somaticMutation).Count() + " germline variants in selected genes.");

            } // while true
        } // WorkerThread
    } // Program
}
