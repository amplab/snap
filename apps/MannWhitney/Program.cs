using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using System.IO;
using MathNet.Numerics;
using ExpressionLib;


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
            //
            // Get a list of sex chromosome genes by reading the _summary.txt file.
            //

            List<string> sexGenes = new List<string>();
            string[] summaryLines = File.ReadAllLines(ExpressionTools.geneScatterGraphsDirectory + "_summary.txt");
            foreach (var summaryLine in summaryLines)
            {
                string[] fields = summaryLine.Split('\t');
                if (fields[5].ToUpper() == "TRUE")
                {
                    sexGenes.Add(fields[0].ToLower());
                }
            }

            var outputLines = new List<OutputLine>();

            ExpressionTools.MannWhitney<Mutation>.GetValue getValue = new ExpressionTools.MannWhitney<Mutation>.GetValue(m => m.RatioOfRatios);
            ExpressionTools.MannWhitney<Mutation>.WhichGroup whichGroup = new ExpressionTools.MannWhitney<Mutation>.WhichGroup(m => m.isSingle);

            foreach (var file in Directory.EnumerateFiles(ExpressionTools.geneScatterGraphsDirectory, "*.txt"))
            {
                string []pathComponents = file.Split('\\');
                string filename = pathComponents[pathComponents.Count() - 1];
                if (filename == "_summary.txt" || filename == "_MannWhitney.txt")
                {
                    // Not a gene
                    continue;
                }

                StreamReader reader = new StreamReader(file);
                reader.ReadLine();  // Skip header

                List<Mutation> mutations = new List<Mutation>();

                string hugoSymbol = "";
                string line;
                List<double> allRatioOfRatios = new List<double>();
                while (null != (line = reader.ReadLine()))
                {
                    string[] fields = line.Split('\t');

                    double tumorDNAFraction = Convert.ToDouble(fields[47]); // The tumorDNAFraction column

                    if (tumorDNAFraction < 1.0/3.0 || tumorDNAFraction > 2.0 / 3.0)
                    {
                        continue;   // Either too small of a subclone (fraction < 1/3), or probable loss of heterozygosity (fraction > 2/3)
                    }

                    var mutation = new Mutation();
                    mutation.isSingle = fields[56].ToLower() == "true";
                    mutation.RatioOfRatios = Convert.ToDouble(fields[52]) / Convert.ToDouble(fields[51]);
                    allRatioOfRatios.Add(mutation.RatioOfRatios);

                    mutations.Add(mutation);

                    hugoSymbol = ExpressionTools.ConvertToNonExcelString(fields[2]); 
                } // input lines
                reader.Close();

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
                p = ExpressionTools.MannWhitney<Mutation>.ComputeMannWhitney(mutations, mutations[0], whichGroup, getValue, out enoughData, out reversed, out nSingle, out nMultiple, out U, out z);
                if (!enoughData)
                {
                    continue;
                }                

                allRatioOfRatios.Sort();

                var outputLine = new OutputLine();
                outputLine.line = ExpressionTools.ConvertToExcelString(hugoSymbol) + "\t" + nSingle + "\t" + nMultiple + "\t" + U + "\t" + z + "\t" + reversed + "\t" + p + "\t" + allRatioOfRatios[allRatioOfRatios.Count() / 2] + "\t" + sexGenes.Contains(hugoSymbol.ToLower());
                outputLine.p = p;
                outputLines.Add(outputLine);
                // Probably should tweak median for even-sized distributions

            }
            var output = new StreamWriter(ExpressionTools.geneScatterGraphsDirectory + "_MannWhitney.txt");
            output.WriteLine("HugoSymbol\tnSingle\tnMultiple\tU\tz\treversed\tp (Pre-Bonferroni)\tmedian ratio-of-ratios\tsex\tp (post-Bonferroni)");

            foreach (var outputLine in outputLines)
            {
                output.WriteLine(outputLine.line + "\t" + outputLine.p * outputLines.Count());
            }
            output.Close();
        }
    }
}
