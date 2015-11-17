using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using System.IO;
using MathNet.Numerics;


namespace MannWhitney
{
    class Program
    {
        class Mutation
        {
            public bool isSingle;
            public double RatioOfRatios;    // (RNAMut / RNANormal) / (DNAMut / DNANormal), or...how far above the main diagonal it is
            static public int Compare(Mutation a, Mutation b)
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
            string[] summaryLines = File.ReadAllLines(@"f:\temp\gene_scatter_graphs\_summary.txt");
            foreach (var summaryLine in summaryLines)
            {
                string[] fields = summaryLine.Split('\t');
                if (fields[5].ToUpper() == "TRUE")
                {
                    sexGenes.Add(fields[0].ToLower());
                }
            }

            var outputLines = new List<OutputLine>();

            foreach (var file in Directory.EnumerateFiles(@"f:\temp\gene_scatter_graphs", "*.txt"))
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

                    hugoSymbol = fields[2]; 
                } // input lines
                reader.Close();

                mutations.Sort(Mutation.Compare);

                double Rsingle = 0;
                double nSingle = 0;
                double nMultiple = 0;
                for (int i = 0; i < mutations.Count(); i++)
                {
                    if (mutations[i].isSingle)
                    {
                        int cumulativeR = i + 1;// +1 because Mann-Whitney is one-based, while C# is 0-based.  BUGBUG: Handle ties properly (which would be more interesting to do if this turned up more than TP53)
                        int n = 1;
                        //
                        // Now add in adjascent indices if there's a tie.  For ties, we use the mean of all the indices in the tied region for each of them (regardless of whether they're single or multiple).
                        //
                        for (int j = i - 1; j > 0 && mutations[j].RatioOfRatios == mutations[i].RatioOfRatios; j--)
                        {
                            cumulativeR += j + 1;
                            n++;
                        }

                        for (int j = i + 1; j < mutations.Count() && mutations[j].RatioOfRatios == mutations[i].RatioOfRatios; j++) {
                            cumulativeR += j+1;
                            n++;
                        }


                        Rsingle += cumulativeR / n;
                        nSingle++;
                    }
                    else
                    {
                        nMultiple++;
                    }
                }

                if (nSingle == 0 || nMultiple == 0) {
                    //
                    // Not enough data, reject this gene
                    //
                    continue;
                }

                double U = (double)Rsingle - (double)nSingle * (nSingle + 1.0) / 2.0;

                double z = Math.Abs(U - nSingle * nMultiple / 2) / Math.Sqrt(nMultiple * nSingle * (nMultiple + nSingle + 1) / 12);

                double p = MathNet.Numerics.Distributions.Normal.CDF(1.0, 1.0, z);
                //
                // We're doing two-tailed, so we need to see if this is on the other end]
                bool reversed = false;
                if (p > 0.5)
                {
                    p = 1.0 - p;
                    reversed = true;
                }

                //
                // And then multiply by two.
                //
                p *= 2.0;

                allRatioOfRatios.Sort();

                var outputLine = new OutputLine();
                outputLine.line = hugoSymbol + "\t" + nSingle + "\t" + nMultiple + "\t" + U + "\t" + z + "\t" + reversed + "\t" + p + "\t" + allRatioOfRatios[allRatioOfRatios.Count() / 2] + "\t" + sexGenes.Contains(hugoSymbol.ToLower());
                outputLine.p = p;
                outputLines.Add(outputLine);
                // Probably should tweak median for even-sized distributions

            }
            var output = new StreamWriter(@"f:\temp\gene_scatter_graphs\_MannWhitney.txt");
            output.WriteLine("HugoSymbol\tnSingle\tnMultiple\tU\tz\treversed\tp (Pre-Bonferroni)\tmedian ratio-of-ratios\tsex\tp (post-Bonferroni)");

            foreach (var outputLine in outputLines)
            {
                output.WriteLine(outputLine.line + "\t" + outputLine.p * outputLines.Count());
            }
            output.Close();



        }
    }
}
