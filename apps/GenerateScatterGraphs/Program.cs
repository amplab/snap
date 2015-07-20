using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using System.IO;

namespace GenerateScatterGraphs
{
    class Program
    {
        class GeneState
        {
            public GeneState()
            {
                for (int i = 0; i < 2; i++)
                {
                    indices[i] = new Dictionary<string, bool>();
                }
            }

            Dictionary<string, bool>[] indices = new Dictionary<string, bool>[2];

        }
        static void Main(string[] args)
        {
            string[] tumorDNA = File.ReadAllLines(@"c:\temp\tumor_dna.txt");
            string[] tumorRNA = File.ReadAllLines(@"c:\temp\tumor_rna.txt");

            string[][] tumorDNAFields = new string[tumorDNA.Count()][];
            string[][] tumorRNAFields = new string[tumorRNA.Count()][];

            for (int i = 0; i < tumorDNA.Count(); i++)
            {
                tumorDNAFields[i] = tumorDNA[i].Split('\t');
            }

            for (int i = 0; i < tumorRNA.Count(); i++)
            {
                tumorRNAFields[i] = tumorRNA[i].Split('\t');
            }

            string[][] measurements = new string[2][];
            measurements[0] = tumorDNA;
            measurements[1] = tumorRNA;

            var genes = new List<string>();

            if (false)
            {
                genes.Add("tp53");
                genes.Add("cdkn2a");
                genes.Add("egfr");
                genes.Add("idh1");
                genes.Add("idh2");
                genes.Add("brca1");
                genes.Add("brca2");
                genes.Add("cebpa");
                genes.Add("npm1");
            }
            else
            {
                foreach (var line in tumorDNA)
                {
                    string[] fields = line.Split('\t');
                    if (!genes.Contains(fields[2].ToLower()))
                    {
                        genes.Add(fields[2].ToLower());
                    }
                }
            }

            foreach (var gene in genes)
            {
                //
                // We need to find all mutations for this gene that are present in both the DNA and RNA sets.  Make indices of them.
                //
                var indices = new Dictionary<string, bool>[2];
                for (int i = 0; i < 2; i++)
                {
                    indices[i] = new Dictionary<string, bool>();
                    foreach (var line in measurements[i])
                    {
                        string[] fields = line.Split('\t');
                        if (fields[2].ToLower() != gene)
                        {
                            continue;   // Wrong gene; ignore it
                        }
                        if (fields[10].ToLower() == "silent")
                        {
                            //
                            // Ignore silent mutations.
                            //
                            continue;
                        }

                        if (!indices[i].ContainsKey(fields[34])) // fields[34] is normal sample ID
                        {
                            indices[i].Add(fields[34], true);
                        }
                    }
                }

                //
                // Now run through each one and keep only the ones that match.
                //
                List<string>[] matchedMeasurements = new List<string>[2];
                for (int i = 0; i < 2; i++)
                {
                    matchedMeasurements[i] = new List<string>();
                    foreach (var line in measurements[i])
                    {
                        string[] fields = line.Split('\t');

                        if (fields[2].ToLower() != gene)
                        {
                            continue;   // Wrong gene; ignore it
                        }

                        if (indices[1 - i].ContainsKey(fields[34]))
                        {
                            matchedMeasurements[i].Add(line);
                        }
                    }
                }


                StreamWriter outFile = new StreamWriter(@"c:\temp\gene_expression_graphs\" + gene + ".txt");
                int nMultiple = 0;
                int nEmitted = 0;
                int nInteresting = 0;
                for (int i = 0; i < matchedMeasurements[0].Count(); i++)
                {
                    //
                    // It it single?
                    //
                    string [] currentFields = matchedMeasurements[0][i].Split('\t');

                    double nNormalDNA = Convert.ToDouble(currentFields[39]);
                    double nTumorDNA = Convert.ToDouble(currentFields[40]);
                    double nNeitherDNA = Convert.ToDouble(currentFields[41]);

                    string[] currentRNA = matchedMeasurements[1][i].Split('\t');
                    double nNormalRNA = Convert.ToDouble(currentRNA[39]);
                    double nTumorRNA = Convert.ToDouble(currentRNA[40]);
                    double nNeitherRNA = Convert.ToDouble(currentRNA[41]);

                    //
                    // Give up on ones that don't have at least 10 DNA and RNA reads.
                    //
                    if (nNormalDNA + nTumorDNA + nNeitherDNA < 10) continue;
                    if (nNormalRNA + nTumorRNA + nNeitherRNA < 10) continue;

                    nEmitted++;
                    outFile.Write(matchedMeasurements[0][i] + "\t" + matchedMeasurements[1][i]);
                    bool single = true;
                    if (i != 0)
                    {
                        string[] priorFields = matchedMeasurements[0][i - 1].Split('\t');
                        if (priorFields[34] == currentFields[34])
                        {
                            single = false;
                            nMultiple++;
                        }
                    }

                    if (single && i != matchedMeasurements[0].Count()-1)
                    {
                        string[] nextFields = matchedMeasurements[0][i + 1].Split('\t');
                        if (nextFields[34] == currentFields[34])
                        {
                            single = false;
                            nMultiple++;
                        }
                    }

                    outFile.Write("\t" + single + "\t" + (nTumorDNA / (nTumorDNA + nNormalDNA + nNeitherDNA)) + "\t" + (nTumorRNA / (nTumorRNA + nNormalRNA + nNeitherRNA)));
                    if ((nTumorDNA / (nTumorDNA + nNormalDNA + nNeitherDNA)) > 0)
                    {
                        outFile.WriteLine("\t" + (nTumorRNA / (nTumorRNA + nNormalRNA + nNeitherRNA) / nTumorDNA / (nTumorDNA + nNormalDNA + nNeitherDNA)));
                    }
                    else
                    {
                        outFile.WriteLine();
                    }

                    if ((nTumorDNA / (nTumorDNA + nNormalDNA + nNeitherDNA)) <= .55 && (nTumorRNA / (nTumorRNA + nNormalRNA + nNeitherRNA)) >= .67)
                    {
                        nInteresting++;
                    }

                }
                outFile.Close();
                Console.WriteLine(gene + " has " + matchedMeasurements[0].Count() + " total matched samples, of which we emitted " + nEmitted + " of which " + nMultiple + " were multiple mutations, " +
                    nInteresting + " were interesting.");
            }
        }
    }
}
