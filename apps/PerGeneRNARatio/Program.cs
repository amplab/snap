using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using ASELib;

namespace PerGeneRNARatio
{
    class Program
    {
        class GeneRNAExpression
        {
            public int n = 0;
            public double total = 0;
            public readonly string hugo_symbol;
            public readonly bool significant;
            public int nOverPoint5 = 0;

            public GeneRNAExpression(string hugo_symbol_, bool significant_)
            {
                hugo_symbol = hugo_symbol_;
                significant = significant_;
            }

            public void addSample(double sample)
            {
                n++;
                total += sample;
                if (sample > 0.5) nOverPoint5++;
            }

            public double mean()
            {
                return total / n;
            }

            public double fractionOverPoint5()
            {
                return (double)nOverPoint5 / n;
            }
        }

        static void Main(string[] args)
        {
            var configuration = ASETools.Configuration.loadFromFile(args);

            if (null == configuration)
            {
                Console.WriteLine("Giving up because we were unable to load configuration.");
                return;
            }

            var significantGenes = ASETools.AllSignificantResultsLine.ReadFile(configuration.finalResultsDirectory + ASETools.AllSignificantResultsFilename).Select(x => x.hugo_symbol).ToList().Distinct().ToList();

            var scatterGraphLines = ASETools.GeneScatterGraphLine.LoadAllGeneScatterGraphLines(configuration.geneScatterGraphsDirectory, false, "*");

            var perGeneExpression = new Dictionary<string, GeneRNAExpression>();
            foreach (var scatterGraphLine in scatterGraphLines)
            {
                if (scatterGraphLine.nMutationsThisGene != 1 || scatterGraphLine.tumorDNAFraction < .4 || scatterGraphLine.tumorDNAFraction > .6 || 
                    scatterGraphLine.tumorDNAReadCounts.nMatchingAlt + scatterGraphLine.tumorDNAReadCounts.nMatchingReference < 10 || 
                    scatterGraphLine.tumorRNAReadCounts.nMatchingAlt + scatterGraphLine.tumorRNAReadCounts.nMatchingReference < 10)
                {
                    continue;
                }

                if (!perGeneExpression.ContainsKey(scatterGraphLine.Hugo_Symbol))
                {
                    perGeneExpression.Add(scatterGraphLine.Hugo_Symbol, new GeneRNAExpression(scatterGraphLine.Hugo_Symbol, significantGenes.Contains(scatterGraphLine.Hugo_Symbol)));
                }

                perGeneExpression[scatterGraphLine.Hugo_Symbol].addSample(scatterGraphLine.tumorRNAFraction);
            }

            var outputLines = perGeneExpression.Select(x => x.Value).ToList();
            outputLines.Sort((x, y) => string.Compare(x.hugo_symbol, y.hugo_symbol));

            var outputFile = ASETools.CreateStreamWriterWithRetry(configuration.finalResultsDirectory + ASETools.PerGeneRNARatioFilename);
            outputFile.WriteLine("hugo symbol\tsignificant\tn\ttotal\tfraction over 0.5\tmean");
            foreach (var outputLine in outputLines)
            {
                outputFile.WriteLine(outputLine.hugo_symbol + "\t" + outputLine.significant + "\t" + outputLine.n + "\t" + outputLine.total + "\t" + outputLine.fractionOverPoint5() + "\t" + outputLine.mean());
            }

            outputFile.WriteLine("**done**");
            outputFile.Close();


        }
    }
}
