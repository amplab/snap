using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using System.IO;
using System.IO.Compression;
using System.Diagnostics;
using ExpressionLib;

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

            //
            // The ratio of the RNA mutant expression ratio to the DNA mutant expression ratio for each mutation that makes the cut.
            //
            public List<double> singleRatioRatios = new List<double>();
            public List<double> multipleRatioRatios = new List<double>();
            public List<double> ratioRatios = new List<double>();
            public List<double> ratioRatiosMidDNA = new List<double>();

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

        public static Dictionary<string, Dictionary<string, NormalizedExpressionAndZScore>> BuildExpressionMultipleMapping(string[] experiments) // maps gene->(tumor RNA Analysis ID -> multiple of mean expression of this gene in this disease for this sample)
        {
            //
            // Makes a mapping from tumor type->gene->mean expression for that gene in that tumor type's tumor RNA.
            //


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

        public static Dictionary<string, string> MakeTumorSampleIdToGenderMap(string[] experiments)
        {
            var participantToGenderMap = new Dictionary<string, string>();

            foreach (var line in experiments)
            {
                string[] fields = line.Split('\t');
                if (fields[15] != "" && fields[0] != "disease_abbr")
                {
                    participantToGenderMap.Add(fields[15], fields[12]);
                }
            }

            return participantToGenderMap;
        }

        public static Dictionary<string, string> MakeTumorSampleToAllcountFileMap(string[] experiments)
        {
            var tumorSampleToAllcountMap = new Dictionary<string, string>();

            foreach (var line in experiments)
            {
                string[] fields = line.Split('\t');

                if (fields[16] != "" && fields[0] != "disease_abbr")
                {
                    tumorSampleToAllcountMap.Add(fields[3], fields[16]);
                }
            }

            return tumorSampleToAllcountMap;
        }

        static double Median(List<double> input, List<double> input2 = null)
        {
            int size = input.Count();
            if (null != input2)
            {
                size += input2.Count();
            }

            if (size == 0) {
                return 0;
            }
            var newList = new List<double>();
            foreach (var element in input)
            {
                newList.Add(element);
            }

            if (null != input2)
            {
                foreach (var element in input2)
                {
                    newList.Add(element);
                }
            }

            newList.Sort();

            if (size % 2 == 1)
            {
                return newList[size / 2];
            }

            return (newList[size / 2] + newList[size / 2 - 1])/2;
        }

        static string[] BuildHistogram(List<double> values, List<double> bucketMaxima)
        {
            int nBuckets = bucketMaxima.Count() + 1;
            int[] buckets = new int[nBuckets];  // +1 is for "More"

            for (int i = 0; i < nBuckets; i++)
            {
                buckets[i] = 0;
            }

            foreach (var value in values)
            {
                for (int i = 0; i < nBuckets; i++)
                {
                    if (i == nBuckets - 1 || value < bucketMaxima[i])   // Yes, we could do a binary search, but it's not worth the effort
                    {
                        buckets[i]++;
                        break;
                    }
                }
            }

            string[] output = new string[nBuckets + 1]; // +1 is for the header
            output[0] = "BucketMax\tCount\tpdf\tcdf";
            int runningTotal = 0;
            for (int i = 0; i < nBuckets; i++)
            {
                runningTotal += buckets[i];
                output[i + 1] = ((i + 1 == nBuckets) ? "More" : Convert.ToString(bucketMaxima[i])) + "\t" + buckets[i] + "\t" + ((double)buckets[i]) / values.Count() + "\t" + ((double)runningTotal) / values.Count();
            }

            return output;
        }



        static public string RNALineToDesignator(string rnaLine)
        {
            var fields = rnaLine.Split('\t');
            string cancer_type = fields[0].ToLower();

            if (cancer_type != "ov" || fields[5] != "36") return cancer_type;

            return "ov_hg18";

        }


        static void Main(string[] args)
        {
            bool useExpression = true;
            const int minReads = 30;    // Minimum reads (in each RNA and DNA) to keep a sample.  There's a similar constant in ProcessVCFs.

            if (args.Count() > 0)
            {
                if (args.Count() > 1 || args[0] != "--e")
                {
                    Console.WriteLine("usage: GenerateScatterGraphs {--e}");
                    return;
                }
                useExpression = false;
            }

            string[] experiments = File.ReadAllLines(@"f:\temp\expression\experiments.txt");
            var normalizedExpression = /*BuildExpressionMultipleMapping(experiments);*/ new Dictionary<string, Dictionary<string, NormalizedExpressionAndZScore>>();    // Don't bother with this anymore
            var tumorSampleIDToGenderMap = MakeTumorSampleIdToGenderMap(experiments);
            var tumorSampleToAllcountMap = MakeTumorSampleToAllcountFileMap(experiments);
            var tumorSampleToHighQualityNuclearRNAReads = new Dictionary<string, long>();

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
            Array.Sort(tumorRNA, (val1, val2) => string.Compare(RNALineToDesignator(val1), RNALineToDesignator(val2))); // Ordered by tumor type (with special ov_hg18)

            string loadedCancerType = "";
            Dictionary<string, Dictionary<int,ExpressionTools.MeanAndStdDev>> expression = null;
            Stopwatch stopwatch = null;

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
                string designator = RNALineToDesignator(line);
                string gene = fields[2].ToLower();

                if (gene == "unknown")
                {
                    continue;
                }

                if (designator != loadedCancerType && useExpression)
                {
                    loadedCancerType = designator;
                    if (null != stopwatch)
                    {
                        Console.WriteLine("" + (stopwatch.ElapsedMilliseconds + 30000) / 60000 + " minutes");
                    }
                    Console.Write("Processing " + loadedCancerType + "...");
                    stopwatch = new Stopwatch();
                    stopwatch.Start();

                    //
                    // Load in the expression file.  We do this the mininum number of times, because we sorted the RNA lines by tumor type already.
                    //
                    expression = ExpressionTools.LoadExpressionFile(@"f:\sequence\reads\tcga\expression\expression_" + designator);
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

                        string outputLine = fields[1] + "\t" + dnaMutation[1];

                        for (int i = 2; i <= 42; i++)
                        {
                            outputLine += "\t" + dnaMutation[i];
                        }

                        //
                        // Write the RNA counts after the DNA counts
                        //
                        for (int i = 39; i <= 42; i++)
                        {
                            outputLine += "\t" + fields[i];
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

                        string ratioOfRatiosString;
                        if (tumorDNALog == 0)
                        {
                            ratioOfRatiosString = "1000";
                        }
                        else
                        {
                            ratioOfRatiosString = Convert.ToString(tumorRNALog / tumorDNALog);
                        }

                        outputLine += "\t" + tumorDNAFraction + "\t" + tumorRNAFraction + "\t" + tumorDNAMultiple + "\t" + tumorRNAMultiple + "\t" + tumorDNALog + "\t" + tumorRNALog + "\t" + normalizedExpressionString + "\t" + ratioOfRatiosString + "\t" + isSingle + "\t";
                        if (geneState.countByCancerType[cancerType] * 20 >= geneState.nWithDNA) // 5% or more
                        {
                            outputLine += cancerType + "\t";
                        }
                        else
                        {
                            outputLine += "other\t";
                        }

                        if (tumorSampleIDToGenderMap.ContainsKey(tumorSampleID))
                        {
                            outputLine += tumorSampleIDToGenderMap[tumorSampleID] = "\t";
                        }
                        else
                        {
                            outputLine += "unknown\t";
                        }

                        //
                        // See if we have the mean & std deviation for this locus, and if so write the z values for normal and mutant expression.  Recall that the mean and std. dev differ by
                        // tumor type, but that we've loaded the values for this tumor type.
                        //
                        int startPosition = Convert.ToInt32(fields[7]);
                        if (expression != null && expression.ContainsKey(fields[6]) && expression[fields[6]].ContainsKey(startPosition))
                        {
                            if (!tumorSampleToHighQualityNuclearRNAReads.ContainsKey(tumorRNAAnalysisID) && tumorSampleToAllcountMap.ContainsKey(tumorRNAAnalysisID))
                            {
                                var reader = new StreamReader(new GZipStream(new StreamReader(tumorSampleToAllcountMap[tumorRNAAnalysisID]).BaseStream, CompressionMode.Decompress));
                                reader.ReadLine();
                                string countLine = reader.ReadLine();
                                string[] countLineFields = countLine.Split(' ');
                                long count = Convert.ToInt64(countLineFields[0]);
                                tumorSampleToHighQualityNuclearRNAReads.Add(tumorRNAAnalysisID, count);

                                reader.Close();
                            }

                            if (!tumorSampleToHighQualityNuclearRNAReads.ContainsKey(tumorRNAAnalysisID)) 
                            {
                                outputLine += "\t\t\t\t\t";
                            } 
                            else 
                            {
                                var stats = expression[fields[6]][startPosition];
                                outputLine += ((((double)nRNATumor / (double)tumorSampleToHighQualityNuclearRNAReads[tumorRNAAnalysisID]) - stats.mean) / stats.stddev) + "\t" +
                                    ((((double)(nRNAReference + nRNANeither) / (double)tumorSampleToHighQualityNuclearRNAReads[tumorRNAAnalysisID]) - stats.mean) / stats.stddev) + "\t" +
                                    ((((double)nRNATumor * 2.0 / (double)tumorSampleToHighQualityNuclearRNAReads[tumorRNAAnalysisID]) - stats.mean) / stats.stddev) + "\t" +
                                    ((((double)(nRNAReference + nRNANeither) * 2.0 / (double)tumorSampleToHighQualityNuclearRNAReads[tumorRNAAnalysisID]) - stats.mean) / stats.stddev) + "\t";

                                if (stats.mean != 0)
                                {
                                    outputLine += ((double)nRNATumor / (double)tumorSampleToHighQualityNuclearRNAReads[tumorRNAAnalysisID]) / stats.mean + "\t" +
                                                  ((double)(nRNAReference + nRNANeither) / (double)tumorSampleToHighQualityNuclearRNAReads[tumorRNAAnalysisID]) / stats.mean;
                                }
                                else
                                {
                                    outputLine += "\t";
                                }
                            }
                            
                        }
                        else
                        {
                            outputLine += "\t\t\t\t\t";
                        }


                        geneState.output.Add(outputLine);

                        if (tumorDNAFraction < .6 && tumorDNAFraction > .1 && tumorRNAFraction >= .67 && tumorRNAFraction < .999)
                        {
                            geneState.nInteresting++;
                        }

                        if (isSingle)
                        {
                            geneState.singleRatioRatios.Add(tumorRNALog / tumorDNALog);
                        }
                        else
                        {
                            geneState.multipleRatioRatios.Add(tumorRNALog / tumorDNALog);
                        }
                        geneState.ratioRatios.Add(tumorRNALog / tumorDNALog);
                        if (tumorDNALog >= 0.5 && tumorDNALog <= 2)
                        {
                            geneState.ratioRatiosMidDNA.Add(tumorRNALog / tumorDNALog);
                        }
                        break;
                    }
                } // for each possible matching DNA mutation
            } // foreach line in the RNA input file.

            StreamWriter summaryFile = new StreamWriter(@"f:\temp\gene_scatter_graphs\_summary.txt");
            summaryFile.WriteLine("Gene\tnSingle\tnMultiple\tnInteresting\t%Interesting\tsex\tmedian single\tmedian multi\tmedian combined\tmedian heterozygous");

            var histogramBucketMaxima = new List<double>();
            for (double value = .001; value <= 1000; value *= Math.Pow(10, .1))
            {
                histogramBucketMaxima.Add(value);
            }

            foreach (var entry in geneStates)
            {
                var geneState = entry.Value;
                if (geneState.output.Count() < 30)
                {
                    continue;   // Not enough patients to be interesting
                }

                summaryFile.WriteLine(geneState.gene + "\t" + geneState.nSingle + "\t" + geneState.nMultiple + "\t" + geneState.nInteresting + "\t" + ((geneState.nInteresting * 100) / (geneState.nSingle + geneState.nMultiple)) + "\t" + geneState.sex + "\t" +
                    Median(geneState.singleRatioRatios) + "\t" + Median(geneState.multipleRatioRatios) + "\t" + Median(geneState.singleRatioRatios, geneState.multipleRatioRatios) + "\t" +  Median(geneState.ratioRatiosMidDNA));

                string[] histogramLines = BuildHistogram(geneState.ratioRatios, histogramBucketMaxima);
                string[] histogram2Lines = BuildHistogram(geneState.ratioRatiosMidDNA, histogramBucketMaxima);

               StreamWriter file = new StreamWriter(@"f:\temp\gene_scatter_graphs\" + geneState.gene + ".txt");
                file.WriteLine("RNAFile\tDNAFile\t" +
                    "Hugo_Symbol\tEntrez_Gene_Id\tCenter\tNCBI_Build\tChromosome\tStart_Position\tEnd_Position\tStrand\tVariant_Classification\tVariant_Type\tReference_Allele\tTumor_Seq_Allele_1\tTumor_Seq_Allele_2\tdbSNP_RS\tdbSNP_Val_Status\t" +
                    "Tumor_Sample_Barcode\tMatched_Norm_Sample_Barcode\tMatch_Norm_Seq_Allele1\tMatch_Norm_Seq_Allele2\tTumor_Validation_Allele1\tTumor_Validation_Allele2\tMatch_Norm_Validation_Allele1\tMatch_Norm_Validation_Allele2\tVerification_Status\t" +
                    "Validation_Status\tMutation_Status\tSequencing_Phase\tSequence_Source\tValidation_Method\tScore\tBAM_File\tSequencer\tTumor_Sample_UUID\tMatched_Norm_Sample_UUID\tFile_Name\tArchive_Name\tLine_Number\t" +
                    "n_DNA_Matching_Reference\tn_DNA_Matching_Tumor\tn_DNA_Matching_Neither\tn_DNA_Matching_Both\t" +
                    "n_RNA_Matching_Reference\tn_RNA_Matching_Tumor\tn_RNA_Matching_Neither\tn_RNA_Matching_Both\t" +
                    "tumorDNAFraction\ttumorRNAFraction\ttumorDNAMultiple\ttumorRNAMultiple\ttumorDNARatio\ttumorRNARatio\tFractionOfMeanExpression\tzOfmeanExpression\tratioOfRatios\tIsSingle\tCancerType\tgender\tzTumor\tzNormal\tz2Tumor\tz2Normal\t%MeanTumor\t%MeanNormal\t" + 
                    histogramLines[0] + "\t" + histogram2Lines[0]);
                for (int i = 0; i < geneState.output.Count(); i++)
                {
                    if (i < histogramLines.Count()-1)
                    {
                        file.WriteLine(geneState.output[i] + "\t" + histogramLines[i+1] + "\t" + histogram2Lines[i+1]);
                    }
                    else
                    {
                        file.WriteLine(geneState.output[i]);
                    }
                }

                file.Close();
            }
            summaryFile.Close();

        } // Main
    }
}
