using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using ASELib;
using System.Diagnostics;

namespace SelectGenes
{
    class Program
    {

        class Gene
        {
            public string Hugo_Symbol = "";
            public Dictionary<int, int> tumorsByMutationCount = new Dictionary<int, int>();
            public int nFlankingMutations = 0;
            public int nRNAMutations = 0;
            public int nTumorsWithFlankingMutations = 0;
            public int nTumorsWithNonFlankingMutations = 0;
        }

        static Dictionary<string, Gene> genesByHugoSymbol = new Dictionary<string, Gene>(); // Cheezy to make this global, but I'm lazy

        static void ProcessSample(string Hugo_Symbol, int nNonFlankingMutations, int nFlankingMutations, int nRNAMutations)
        {
            if (!genesByHugoSymbol.ContainsKey(Hugo_Symbol))
            {
                genesByHugoSymbol.Add(Hugo_Symbol, new Gene());
                genesByHugoSymbol[Hugo_Symbol].Hugo_Symbol = Hugo_Symbol;
            }

            var gene = genesByHugoSymbol[Hugo_Symbol];
            if (!gene.tumorsByMutationCount.ContainsKey(nNonFlankingMutations))
            {
                gene.tumorsByMutationCount.Add(nNonFlankingMutations, 0);
            }

            gene.tumorsByMutationCount[nNonFlankingMutations]++;
            gene.nFlankingMutations += nFlankingMutations;
            if (nFlankingMutations > 0)
            {
                gene.nTumorsWithFlankingMutations++;
            }

            if (nNonFlankingMutations > 0)
            {
                gene.nTumorsWithNonFlankingMutations++;
            }

            gene.nRNAMutations += nRNAMutations;
        }

        static int Main(string[] args)
        {
            var timer = new Stopwatch();
            timer.Start();

            var configuration = ASETools.ASEConfirguation.loadFromFile(args);

            if (configuration.commandLineArgs.Count() != 0)
            {
                Console.WriteLine("usage: SelectGenes {-configuration configFileName}");
                return 1;
            }

            var cases = ASETools.Case.LoadCases(configuration.casesFilePathname);

            if (null == cases)
            {
                Console.WriteLine("You must have selected cases before you can select genes.");
                return 1;
            }

            int nMafsProcessed = 0;

            foreach (var caseEntry in cases)
            {
                var case_ = caseEntry.Value;
                if (case_.extracted_maf_lines_filename == "")
                {
                    Console.WriteLine("Case " + case_.case_id + " does not have an extracted MAF lines file.  Aborting.");
                    return 1;
                }

                var mafLines = ASETools.MAFLine.ReadFile(case_.extracted_maf_lines_filename, case_.maf_file_id, false);

                mafLines.Sort((a, b) => String.Compare(a.Hugo_Symbol, b.Hugo_Symbol));

                string currentHugoSymbol = "";
                int nNonFlankingMutations = 0;
                int nFlankingMutations = 0;
                int nRNAMutations = 0;

                for (int i = 0; i < mafLines.Count(); i++)
                {
                    var mafLine = mafLines[i];
                    if (mafLine.Hugo_Symbol != currentHugoSymbol)
                    {
                        if (currentHugoSymbol != "")
                        {
                            ProcessSample(currentHugoSymbol, nNonFlankingMutations, nFlankingMutations, nRNAMutations);
                        }
                        nFlankingMutations = 0;
                        nNonFlankingMutations = 0;
                        nRNAMutations = 0;
                        currentHugoSymbol = mafLine.Hugo_Symbol;
                    }

                    if (mafLine.Variant_Classification == "5'Flank" || mafLine.Variant_Classification == "3'Flank")
                    {
                        nFlankingMutations++;
                    }
                    else
                    {
                        nNonFlankingMutations++;

                        if (mafLine.Variant_Classification == "RNA")
                        {
                            nRNAMutations++;
                        }
                    }
                } // foreach maf line

                if (currentHugoSymbol != "")
                {
                    ProcessSample(currentHugoSymbol, nNonFlankingMutations, nFlankingMutations, nRNAMutations);
                }

                if (++nMafsProcessed % 100 == 0)
                {
                    Console.Write(".");
                }
            } // foreach case
            Console.WriteLine();

            var selectedGenes = new List<ASETools.SelectedGene>();
            
            foreach (var geneEntry in genesByHugoSymbol)
            {
                var gene = geneEntry.Value;

                if (gene.nTumorsWithNonFlankingMutations < configuration.nTumorsToIncludeGene)
                {
                    continue;
                }

                selectedGenes.Add(new ASETools.SelectedGene(gene.Hugo_Symbol, gene.nFlankingMutations, gene.nRNAMutations, gene.tumorsByMutationCount));
            } // foreach gene

            ASETools.SelectedGene.SaveAllToFile(selectedGenes, configuration.selectedGenesFilename);

            Console.WriteLine("Selected " + selectedGenes.Count() + " genes from a pool of " + genesByHugoSymbol.Count() + " possible genes in " + ASETools.ElapsedTimeInSeconds(timer));

            return 0;
        }
    }
}
