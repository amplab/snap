using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using System.Threading;
using System.IO;
using ASELib;
using System.Diagnostics;

namespace SummarizeSelectedVariantCountsByGene
{
    class Program
    {

        class PerGeneData
        {
            public PerGeneData(string Hugo_Symbol_)
            {
                HugoSymbol = Hugo_Symbol_;
                nSelectedVariantsInGeneByMutationCount[0] = 0;
                nSelectedVariantsInGeneByMutationCount[1] = 0;
                nSelectedVariantsInGeneByMutationCount[2] = 0;
            }

            public string HugoSymbol;
            public int nSelectedVariantsInGene = 0;
            public int[] nSelectedVariantsInGeneByMutationCount = new int[3];  // Really, by 0, 1, or 2 where 2 means > 1.
        }

        static void Main(string[] args)
        {
            var timer = new Stopwatch();
            timer.Start();

            var configuration = ASETools.ASEConfirguation.loadFromFile(args);

            if (null == configuration)
            {
                Console.WriteLine("Unable to load configuration.");
                return;
            }

            if (configuration.commandLineArgs.Count() != 0)
            {
                Console.WriteLine("usage: SummarizeSelectedVariantCountsByGene");
                return;
            }

            var cases = ASETools.Case.LoadCases(configuration.casesFilePathname);

            if (null == cases)
            {
                Console.WriteLine("Unable to load cases.");
                return;
            }

            int nMissingInput = 0;

            var perGene = new Dictionary<string, PerGeneData>();

            foreach (var caseEntry in cases)
            {
                var case_ = caseEntry.Value;

                if (case_.selected_variant_counts_by_gene_filename == "")
                {
                    nMissingInput++;
                    continue;
                }

                var selectedVariantCounts = ASETools.SelectedVariantCountByGene.LoadFromFile(case_.selected_variant_counts_by_gene_filename);

                if (null == selectedVariantCounts)
                {
                    Console.WriteLine("Couldn't load selected variant counts from " + case_.selected_variant_counts_by_gene_filename + ".  Ignoring.");
                    nMissingInput++;
                    continue;
                }

                foreach (var selectedVariantCount in selectedVariantCounts)
                {
                    if (!perGene.ContainsKey(selectedVariantCount.Hugo_Symbol))
                    {
                        perGene.Add(selectedVariantCount.Hugo_Symbol, new PerGeneData(selectedVariantCount.Hugo_Symbol));
                    }

                    perGene[selectedVariantCount.Hugo_Symbol].nSelectedVariantsInGene += selectedVariantCount.nVariantsInGene;
                    perGene[selectedVariantCount.Hugo_Symbol].nSelectedVariantsInGeneByMutationCount[ASETools.ZeroOneMany(selectedVariantCount.nSomaticMutations)] += selectedVariantCount.nVariantsInGene;
                }
            } // foreach case

            var output = ASETools.CreateStreamWriterWithRetry(configuration.finalResultsDirectory + ASETools.genesWithSelectedVariantsFilename);

            if (null == output)
            {
                Console.WriteLine("Unable to open output file: " + configuration.finalResultsDirectory + ASETools.genesWithSelectedVariantsFilename);
                return;
            }

            var geneLocationInformation = new ASETools.GeneLocationsByNameAndChromosome(ASETools.readKnownGeneFile(ASETools.ASEConfirguation.defaultGeneLocationInformation));

            output.WriteLine("Hugo Symbol\tSize of gene in bases (max tx end - min tx start)\tn Selected Variants in the gene\tn in tumors with no somatic mutations in this gene\tn in tumors with exactly 1 somatic mutation in this gene\tn in tumors with more than one somatic ");

            var perGeneList = new List<PerGeneData>();
            foreach (var perGeneEntry in perGene)
            {
                perGeneList.Add(perGeneEntry.Value);
            }

            perGeneList.Sort((a, b) => String.Compare(a.HugoSymbol, b.HugoSymbol));

            for (int i =0; i < perGeneList.Count(); i++)
            {
                var Hugo_Symbol = perGeneList[i].HugoSymbol;
                int size = -1;
                if (geneLocationInformation.genesByName.ContainsKey(Hugo_Symbol))
                {
                    size = (geneLocationInformation.genesByName[Hugo_Symbol].maxLocus - geneLocationInformation.genesByName[Hugo_Symbol].minLocus);
                }
                output.Write(Hugo_Symbol + "\t" + size + "\t" + perGeneList[i].nSelectedVariantsInGene);
                for (int j = 0; j < 3; j++)
                {
                    output.Write("\t" + perGeneList[i].nSelectedVariantsInGeneByMutationCount[j]);
                }
                output.WriteLine();
            }

            output.WriteLine("**done**");
            output.Close();

            Console.WriteLine("Processed " + cases.Count() + " cases (of which " + nMissingInput + " were damaged or missing input) in " + ASETools.ElapsedTimeInSeconds(timer));
        } // Main 
    }
}
