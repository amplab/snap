using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using System.IO;
using MathNet.Numerics;
using ASELib;
using System.Diagnostics;


namespace MannWhitney
{
    class Program
    {
        class Mutation: IComparer<Mutation>
        {
            public bool isSingle;
            public double RatioOfRatios;    // (RNAMut / RNANormal) / (DNAMut / DNANormal), or...how far above the main diagonal it is
            public int Compare(Mutation a, Mutation b)
            {
                return xCompare(a, b);
            }

            static public int xCompare(Mutation a, Mutation b)
            {
                if (a.RatioOfRatios > b.RatioOfRatios) return 1;
                if (a.RatioOfRatios < b.RatioOfRatios) return -1;
                return 0;
            }

        }

        class OutputLine
        {
            public string line;
            public double p;
        }
        static void Main(string[] args)
        {
            var endToEndTimer = new Stopwatch();
            endToEndTimer.Start();

            var configuration = ASETools.Configuration.loadFromFile(args);
            if (null == configuration)
            {
                Console.WriteLine("Unable to load configuration.");
            }

            //
            // Get a list of sex chromosome genes by reading the _summary.txt file.
            //

            List<string> sexGenes = new List<string>();
            var summaryLines = ASETools.SummaryLine.ReadFromFile(configuration.geneScatterGraphsDirectory + ASETools.scatterGraphsSummaryFilename);

            if (null == summaryLines)
            {
                Console.WriteLine("Unable to read summary file " + configuration.geneScatterGraphsDirectory + ASETools.scatterGraphsSummaryFilename);
                return;
            }

            foreach (var summaryLine in summaryLines)
            {
                if (summaryLine.sex)
                {
                    sexGenes.Add(summaryLine.gene.ToLower());
                }
            }

            var outputLines = new List<OutputLine>();

            ASETools.MannWhitney<Mutation>.GetValue getValue = new ASETools.MannWhitney<Mutation>.GetValue(m => m.RatioOfRatios);
            ASETools.MannWhitney<Mutation>.WhichGroup whichGroup = new ASETools.MannWhitney<Mutation>.WhichGroup(m => m.isSingle);
            var scatterGraphLines = ASETools.GeneScatterGraphLine.LoadAllGeneScatterGraphEntries(configuration.geneScatterGraphsDirectory, false, "*");

            var allHugoSymbols = new HashSet<string>();
            scatterGraphLines.ForEach(s => allHugoSymbols.Add(s.Hugo_Symbol));

            int nGenesProcessed = 0;

            foreach (var hugoSymbol in allHugoSymbols)
            {
                nGenesProcessed++;
                if (nGenesProcessed % 100 == 0)
                {
                    Console.Write(".");
                }

                List<Mutation> mutations = new List<Mutation>();
                List<double> allRatioOfRatios = new List<double>();

                foreach (var scatterGraphLine in scatterGraphLines.Where(g => g.Hugo_Symbol == hugoSymbol))
                {

                    if (scatterGraphLine.tumorDNAFraction < 1.0 / 3.0 || scatterGraphLine.tumorDNAFraction > 2.0 / 3.0)
                    {
                        continue;   // Either too small of a subclone (fraction < 1/3), or probable loss of heterozygosity (fraction > 2/3)
                    }

                    var mutation = new Mutation();
                    mutation.isSingle = !scatterGraphLine.MultipleMutationsInThisGene;
                    mutation.RatioOfRatios = scatterGraphLine.ratioOfRatios;
                    allRatioOfRatios.Add(mutation.RatioOfRatios);

                    mutations.Add(mutation);
                } // input lines


                double nSingle = 0;
                double nMultiple = 0;
                double p;
                bool reversed;
                double U;
                double z;

                bool enoughData;
                if (mutations.Count() == 0)
                {
                    continue;
                }
                p = ASETools.MannWhitney<Mutation>.ComputeMannWhitney(mutations, mutations[0], whichGroup, getValue, out enoughData, out reversed, out nSingle, out nMultiple, out U, out z, true, 10);
                if (!enoughData)
                {
                    continue;
                }

                allRatioOfRatios.Sort();

                var outputLine = new OutputLine();
                outputLine.line = ASETools.ConvertToExcelString(hugoSymbol) + "\t" + nSingle + "\t" + nMultiple + "\t" + U + "\t" + z + "\t" + reversed + "\t" + p + "\t" + allRatioOfRatios[allRatioOfRatios.Count() / 2] + "\t" + sexGenes.Contains(hugoSymbol.ToLower());
                outputLine.p = p;
                outputLines.Add(outputLine);
                // Probably should tweak median for even-sized distributions

            } // foreach hugo symhol

            var output = ASETools.CreateStreamWriterWithRetry(configuration.geneScatterGraphsDirectory + ASETools.mannWhitneyFilename);
            output.WriteLine("HugoSymbol\tnSingle\tnMultiple\tU\tz\treversed\tp (Pre-Bonferroni)\tmedian ratio-of-ratios\tsex\tp (post-Bonferroni)");

            foreach (var outputLine in outputLines)
            {
                output.WriteLine(outputLine.line + "\t" + outputLine.p * outputLines.Count());
            }
            output.Close();

            Console.WriteLine("Processed " + allHugoSymbols.Count() + " genes in " + ASETools.ElapsedTimeInSeconds(endToEndTimer));
        } // Main
    }
}
