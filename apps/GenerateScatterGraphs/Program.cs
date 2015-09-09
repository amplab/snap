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
            public List<string> output = new List<string>();

            public Dictionary<string, int> countByCancerType = new Dictionary<string, int>();

            public Dictionary<string, List<string[]>> dnaMutationsByTumorAnalysis = new Dictionary<string, List<string []>>();

            public int nSingle = 0;
            public int nMultiple = 0;
            public int nInteresting = 0;
            public int nWithDNA = 0;

            public bool sex;
            public string chrom;
        }

        const double maxMultiple = 30.0;

        static double ComputeMultiple(double nTumor, double nNormal)
        {
            if (nTumor == 0) return - maxMultiple;

            if (nNormal == 0) return maxMultiple;

            if (nTumor >= nNormal) return nTumor / nNormal;
            return -nNormal / nTumor;
        }

        static double ComputeRatio(double nTumor, double nNormal)
        {
            if (nTumor == 0) return .0001;

            if (nNormal == 0) return 1000;

            return nTumor / nNormal;
        }

        public struct NormalizedExpressionAndZScore
        {
            public double normalizedExpression;
            public double z;
        }

        public static Dictionary<string, Dictionary<string, NormalizedExpressionAndZScore>> BuildExpressionMultipleMapping() // maps gene->(tumor RNA Analysis ID -> multiple of mean expression of this gene in this disease for this sample)
        {
            //
            // Makes a mapping from tumor type->gene->mean expression for that gene in that tumor type's tumor RNA.
            //

            string[] experiments = File.ReadAllLines(@"f:\temp\expression\experiments.txt");

            var diseaseByTumorRNA = new Dictionary<string, string>();

            for (int i = 1; i < experiments.Count(); i++) // start at 1 to skip the header line
            {
                string[] fields = experiments[i].Split('\t');
                diseaseByTumorRNA.Add(fields[3], fields[0]);
            }

            experiments = null;

            var isoformToGene = new Dictionary<string, string>();
            string[] lines = File.ReadAllLines(@"f:\temp\expression\genes_and_isoforms.txt");
            foreach (var line in lines)
            {
                string[] fields = line.Split('\t');
                isoformToGene.Add(fields[1].ToLower(), fields[0].ToLower());
            }

            var expressionsByTumorTypeAndGene = new Dictionary<string, Dictionary<string, List<double>>>();
            var expressionByGeneAndAnalysisID = new Dictionary<string, Dictionary<string, double>>();  // The same shape as the final mapping, but not normalized
            string[] isoformCounts = File.ReadAllLines(@"f:\sequence\reads\tcga\isoform_counts.txt");

            foreach (var line in isoformCounts)
            {
                string[] fields = line.Split('\t');
                //
                // The first field is actually filename:isoform (it came that way from findstr).
                //
                string[] pathnameAndIsoform = fields[0].Split(':');
                string pathname = pathnameAndIsoform[0];
                string extendedIsoform = pathnameAndIsoform[1];

                string[] isoformPieces = extendedIsoform.Split('.');  // isoform.version
                string isoform = isoformPieces[0];

                if (!isoformToGene.ContainsKey(isoform)) {
                    Console.WriteLine("Saw unexpected isoform " + isoform);
                    continue;
                }
                string gene = isoformToGene[isoform];

                //
                // Furthermore, the filename is of the form (prefix)\<analysis-id>-isoforms.txt.
                //
                string[] pathnameComponents = pathname.Split('\\');
                string filename = pathnameComponents[pathnameComponents.Count() - 1];

                //
                // And the filename is <analysisID>-isoforms.txt.
                //
                if (filename.Count() != 49 || filename.Substring(36) != "-isoforms.txt")
                {
                    Console.WriteLine("Found incorrectly-formatted isoform file name in isoform_counts.txt" + pathname + " (" + filename + "?");
                    continue;
                }

                string tumorRNAAnalysisId = filename.Substring(0,36);

                if (!diseaseByTumorRNA.ContainsKey(tumorRNAAnalysisId))
                {
                    Console.WriteLine("Can't find tumor RNA analysis ID in experiments list for " + tumorRNAAnalysisId);
                    continue;
                }

                string disease = diseaseByTumorRNA[tumorRNAAnalysisId];
                if (!expressionsByTumorTypeAndGene.ContainsKey(disease))
                {
                    expressionsByTumorTypeAndGene.Add(disease, new Dictionary<string, List<double>>());
                }

                if (!expressionsByTumorTypeAndGene[disease].ContainsKey(gene))
                {
                    expressionsByTumorTypeAndGene[disease].Add(gene, new List<double>());
                }

                double rawExpression = Convert.ToDouble(fields[14]);
                expressionsByTumorTypeAndGene[disease][gene].Add(rawExpression);

                if (!expressionByGeneAndAnalysisID.ContainsKey(gene))
                {
                    expressionByGeneAndAnalysisID.Add(gene, new Dictionary<string, double>());
                }

                expressionByGeneAndAnalysisID[gene].Add(tumorRNAAnalysisId, rawExpression);
            }

            var meanExpression = new Dictionary<string, Dictionary<string, double>>();
            var stdDeviation = new Dictionary<string, Dictionary<string, double>>();

            foreach (var entry in expressionsByTumorTypeAndGene)
            {
                string disease = entry.Key;
                meanExpression.Add(disease, new Dictionary<string, double>());
                stdDeviation.Add(disease, new Dictionary<string, double>());
                foreach (var geneEntry in entry.Value)
                {
                    string gene = geneEntry.Key;
                    double mean = expressionsByTumorTypeAndGene[disease][gene].Average();
                    meanExpression[disease].Add(gene, mean);
                    double sumOfSquareDeviations = 0;
                    foreach (var value in expressionsByTumorTypeAndGene[disease][gene]) {
                        double deviation = value - mean;
                        sumOfSquareDeviations += deviation * deviation;
                    }
                    stdDeviation[disease].Add(gene, Math.Sqrt(sumOfSquareDeviations / expressionsByTumorTypeAndGene[disease][gene].Count()));
                }
            }

            //
            // Finally, go back through the isoforms again and build the final mapping now that we have the means.
            //
            var normalizedExpressionByGeneAndAnalysisID = new Dictionary<string, Dictionary<string, NormalizedExpressionAndZScore>>();  // The same shape as the final mapping, but not normalized
            foreach (var geneEntry in expressionByGeneAndAnalysisID)
            {
                string gene = geneEntry.Key;
                normalizedExpressionByGeneAndAnalysisID.Add(gene, new Dictionary<string, NormalizedExpressionAndZScore>());
                foreach (var analysisEntry in geneEntry.Value) 
                { 
                    string analysisID = analysisEntry.Key;
                    var value = new NormalizedExpressionAndZScore();
                    value.normalizedExpression = expressionByGeneAndAnalysisID[gene][analysisID] / meanExpression[diseaseByTumorRNA[analysisID]][gene];
                    value.z = (expressionByGeneAndAnalysisID[gene][analysisID] - meanExpression[diseaseByTumorRNA[analysisID]][gene]) / stdDeviation[diseaseByTumorRNA[analysisID]][gene];

                    normalizedExpressionByGeneAndAnalysisID[gene].Add(analysisID, value);
                }
            }

            return normalizedExpressionByGeneAndAnalysisID;
        }
        static void Main(string[] args)
        {
            const int minReads = 10;    // Minimum reads (in each RNA and DNA) to keep a sample.

            var normalizedExpression = BuildExpressionMultipleMapping();

            string[] tumorDNA = File.ReadAllLines(@"f:\sequence\reads\tcga\tumor_dna.txt");

            var geneStates = new Dictionary<string, GeneState>();

            foreach (var line in tumorDNA)
            {
                string[] fields = line.Split('\t');

                if (fields[0] == "Cancer_Type")
                {
                    //
                    // This is a header line.  Ignore it.
                    //
                    continue;
                }

                string cancerType = fields[0].ToLower();

                string variantClassification = fields[10];
                if (variantClassification.ToLower() == "silent")
                {
                    //
                    // Just ignore silent mutations.
                    //
                    continue;
                }

                string gene = fields[2].ToLower();

                if (gene == "unknown" || gene == ".")
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
                    geneStates[gene].sex = true;    // Err on the side of being sex linked
                }

                if (!geneStates[gene].dnaMutationsByTumorAnalysis.ContainsKey(tumorSampleID))
                {
                    geneStates[gene].dnaMutationsByTumorAnalysis.Add(tumorSampleID, new List<string[]>());
                }

                geneStates[gene].dnaMutationsByTumorAnalysis[tumorSampleID].Add(fields);

                if (!geneStates[gene].countByCancerType.ContainsKey(cancerType))
                {
                    geneStates[gene].countByCancerType.Add(cancerType, 0);
                }
                geneStates[gene].countByCancerType[cancerType]++;
                geneStates[gene].nWithDNA++;
            }

            tumorDNA = null;    // Let it GC this.

            string[] tumorRNA = File.ReadAllLines(@"f:\sequence\reads\tcga\tumor_rna.txt");

            //
            // Now run through the RNA and look for matches.
            //
            foreach (var line in tumorRNA)
            {
                string [] fields = line.Split('\t');

                if (fields[0] == "Cancer_Type")
                {
                    //
                    // This is a header line.  Ignore it.
                    //
                    continue;
                }

                string cancerType = fields[0].ToLower();
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

                        string ouptutLine = fields[1] + "\t" + dnaMutation[1];

                        for (int i = 2; i <= 42; i++)
                        {
                            ouptutLine += "\t" + dnaMutation[i];
                        }

                        //
                        // Write the RNA counts after the DNA counts
                        //
                        for (int i = 39; i <= 42; i++)
                        {
                            ouptutLine += "\t" + fields[i];
                        }

                        bool isSingle = geneState.dnaMutationsByTumorAnalysis[tumorSampleID].Count() == 1;
                        if (isSingle)
                        {
                            geneState.nSingle++;
                        }
                        else
                        {
                            geneState.nMultiple++;
                        }
                        double tumorDNAFraction = nDNATumor / (nDNATumor + nDNAReference + nDNANeither);
                        double tumorRNAFraction = nRNATumor / (nRNATumor + nRNAReference + nRNANeither);
                        double tumorDNAMultiple = ComputeMultiple(nDNATumor, nDNAReference + nDNANeither);
                        double tumorRNAMultiple = ComputeMultiple(nRNATumor, nRNAReference + nRNANeither);
                        double tumorDNALog = ComputeRatio(nDNATumor, nDNAReference + nDNANeither);
                        double tumorRNALog = ComputeRatio(nRNATumor, nRNAReference + nRNANeither);

                        string [] filenameParts = fields[1].Split('\\');
                        if (filenameParts.Count() != 2 || filenameParts[1].Count() < 36) {
                            Console.WriteLine("Can't parse filename in rna input " + fields[1]);
                        }
                        string tumorRNAAnalysisID = filenameParts[1].Substring(0,36);

                        string normalizedExpressionString;
                        if (normalizedExpression.ContainsKey(gene) && normalizedExpression[gene].ContainsKey(tumorRNAAnalysisID))
                        {
                            normalizedExpressionString = Convert.ToString(normalizedExpression[gene][tumorRNAAnalysisID].normalizedExpression) + "\t" + Convert.ToString(normalizedExpression[gene][tumorRNAAnalysisID].z);
                        }
                        else
                        {
                            normalizedExpressionString = "\t";
                        }

                        ouptutLine += "\t" + tumorDNAFraction + "\t" + tumorRNAFraction + "\t" + tumorDNAMultiple + "\t" + tumorRNAMultiple + "\t" + tumorDNALog + "\t" + tumorRNALog + "\t" + normalizedExpressionString + "\t" + isSingle +"\t";
                        if (geneState.countByCancerType[cancerType] * 20 >= geneState.nWithDNA) // 5% or more
                        {
                            ouptutLine += cancerType;
                        }
                        else
                        {
                            ouptutLine += "other";
                        }

                        geneState.output.Add(ouptutLine);

                        if (tumorDNAFraction < .6 && tumorDNAFraction > .1 && tumorRNAFraction >= .67 && tumorRNAFraction < .999)
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
                if (geneState.output.Count() < 10)
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
                    "tumorDNAFraction\ttumorRNAFraction\ttumorDNAMultiple\ttumorRNAMultiple\ttumorDNARatio\ttumorRNARatio\tFractionOfMeanExpression\tzOfmeanExpression\tIsSingle\tCancerType");
                foreach (var line in geneState.output)
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
