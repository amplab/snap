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
            public string gene;
            public List<string> single = new List<string>();
            public List<string> multiple = new List<string>();

            public Dictionary<string, List<string[]>> dnaMutationsByTumorAnalysis = new Dictionary<string, List<string []>>();

            public int nSingle = 0;
            public int nMultiple = 0;
            public int nInteresting = 0;

            public bool sex;
            public string chrom;
        }
 
        static void Main(string[] args)
        {
            const int minReads = 10;    // Minimum reads (in each RNA and DNA) to keep a sample.

            string[] tumorDNA = File.ReadAllLines(@"f:\temp\tumor_dna.txt");

            var geneStates = new Dictionary<string, GeneState>();

            foreach (var line in tumorDNA)
            {
                string[] fields = line.Split('\t');

                string variantClassification = fields[10];
                if (variantClassification.ToLower() == "silent")
                {
                    //
                    // Just ignore silent mutations.
                    //
                    continue;
                }

                string gene = fields[2].ToLower();

                if (gene == "unknown")
                {
                    continue;
                }

                string tumorSampleID = fields[34];
                string chrom = fields[6].ToLower();
                bool sex = chrom == "x" || chrom == "y" || chrom == "chrx" || chrom == "chry";

                if (!geneStates.ContainsKey(gene))
                {
                    geneStates.Add(gene, new GeneState());
                    geneStates[gene].gene = gene;
                    geneStates[gene].sex = sex;
                    geneStates[gene].chrom = chrom;
                }

                if (geneStates[gene].sex != sex)
                {
                    Console.WriteLine("Disagreement on sex for gene " + gene + " this record is " + chrom + " and first one was " + geneStates[gene].chrom);
                }

                if (!geneStates[gene].dnaMutationsByTumorAnalysis.ContainsKey(tumorSampleID))
                {
                    geneStates[gene].dnaMutationsByTumorAnalysis.Add(tumorSampleID, new List<string[]>());
                }

                geneStates[gene].dnaMutationsByTumorAnalysis[tumorSampleID].Add(fields);
            }

            tumorDNA = null;    // Let it GC this.

            string[] tumorRNA = File.ReadAllLines(@"f:\temp\tumor_rna.txt");

            //
            // Now run through the RNA and look for matches.
            //
            foreach (var line in tumorRNA)
            {
                string [] fields = line.Split('\t');
                string gene = fields[2].ToLower();

                if (gene == "unknown")
                {
                    continue;
                }

                string tumorSampleID = fields[34];

                double nRNAReference = Convert.ToInt32(fields[39]);
                double nRNATumor = Convert.ToInt32(fields[40]);
                double nRNANeither = Convert.ToInt32(fields[41]);

                if (nRNAReference + nRNATumor + nRNANeither < minReads)
                {
                    continue;
                }

 
                if (!geneStates.ContainsKey(gene) || !geneStates[gene].dnaMutationsByTumorAnalysis.ContainsKey(tumorSampleID))
                {
                    continue;
                }

                var geneState = geneStates[gene];

                //
                // Look for this particular mutation (they should always be there if we get this far).  Fields[2..37] should match (2 is a filename, after 37
                // are counts; all differ between RNA and DNA.
                //
                foreach (var dnaMutation in geneState.dnaMutationsByTumorAnalysis[tumorSampleID])
                {
                    bool foundIt = true;
                    for (int i = 2; i <= 37; i++)
                    {
                        if (fields[i] != dnaMutation[i])
                        {
                            foundIt = false;
                            break;
                        }
                    }


                    if (foundIt)
                    {
                        double nDNAReference = Convert.ToInt32(dnaMutation[39]);
                        double nDNATumor = Convert.ToInt32(dnaMutation[40]);
                        double nDNANeither = Convert.ToInt32(dnaMutation[41]);

                        if (nDNAReference + nDNATumor + nDNANeither < minReads)
                        {
                            continue;
                        }

                        string thisLine = fields[1] + "\t" + dnaMutation[1];

                        for (int i = 2; i <= 42; i++)
                        {
                            thisLine += "\t" + dnaMutation[i];
                        }

                        //
                        // Write the RNA counts after the DNA counts
                        //
                        for (int i = 39; i <= 42; i++)
                        {
                            thisLine += "\t" + fields[i];
                        }

                        bool isSingle = geneState.dnaMutationsByTumorAnalysis[tumorSampleID].Count() == 1;
                        double tumorDNAFraction = nDNATumor / (nDNATumor + nDNAReference + nDNANeither);
                        double tumorRNAFraction = nRNATumor / (nRNATumor + nRNAReference + nRNANeither);
                        thisLine += "\t" + tumorDNAFraction + "\t" + tumorRNAFraction + "\t" + isSingle;

                        if (isSingle)
                        {
                            geneState.nSingle++;
                            geneState.single.Add(thisLine);
                        }
                        else
                        {
                            geneState.nMultiple++;
                            geneState.multiple.Add(thisLine);
                        }

                        if (tumorDNAFraction < .6 && tumorDNAFraction > .1 && tumorRNAFraction >= .67 && tumorRNAFraction < .99)
                        {
                            geneState.nInteresting++;
                        }
                        break;
                    }
                } // for each possible matching DNA mutation
            } // foreach line in the RNA input file.

            StreamWriter summaryFile = new StreamWriter(@"f:\temp\gene_scatter_graphs\_summary.txt");
            summaryFile.WriteLine("Gene\tnSingle\tnMultiple\tnInteresting\t%Interesting\tsex");

            foreach (var entry in geneStates)
            {
                var geneState = entry.Value;
                if (geneState.single.Count() + geneState.multiple.Count()< 10)
                {
                    continue;   // Not enough patients to be interesting
                }
                StreamWriter file = new StreamWriter(@"f:\temp\gene_scatter_graphs\" + geneState.gene + ".txt");
                file.WriteLine("RNAFile\tDNAFile\t" +
                    "Hugo_Symbol\tEntrez_Gene_Id\tCenter\tNCBI_Build\tChromosome\tStart_Position\tEnd_Position\tStrand\tVariant_Classification\tVariant_Type\tReference_Allele\tTumor_Seq_Allele_1\tTumor_Seq_Allele_2\tdbSNP_RS\tdbSNP_Val_Status\t" +
                    "Tumor_Sample_Barcode\tMatched_Norm_Sample_Barcode\tMatch_Norm_Seq_Allele1\tMatch_Norm_Seq_Allele2\tTumor_Validation_Allele1\tTumor_Validation_Allele2\tMatch_Norm_Validation_Allele1\tMatch_Norm_Validation_Allele2\tVerification_Status\t" +
                    "Validation_Status\tMutation_Status\tSequencing_Phase\tSequence_Source\tValidation_Method\tScore\tBAM_File\tSequencer\tTumor_Sample_UUID\tMatched_Norm_Sample_UUID\tFile_Name\tArchive_Name\tLine_Number\t" +
                    "n_DNA_Matching_Reference\tn_DNA_Matching_Tumor\tn_DNA_Matching_Neither\tn_DNA_Matching_Both\t" +
                    "n_RNA_Matching_Reference\tn_RNA_Matching_Tumor\tn_RNA_Matching_Neither\tn_RNA_Matching_Both\t" +
                    "tumorDNAFraction\ttumorRNAFraction\tIsSingle");
                foreach (var line in geneState.multiple)
                {
                    file.WriteLine(line);
                }
                foreach (var line in geneState.single)
                {
                    file.WriteLine(line);
                }
                file.Close();
                summaryFile.WriteLine(geneState.gene + "\t" + geneState.nSingle + "\t" + geneState.nMultiple + "\t" + geneState.nInteresting + "\t" + ((geneState.nInteresting * 100) / (geneState.nSingle + geneState.nMultiple)) + "\t" + geneState.sex);
            }
            summaryFile.Close();

        } // Main
    }
}
