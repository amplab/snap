using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using System.IO;
using System.IO.Compression;
using System.Diagnostics;
using ASELib;

namespace GenerateScatterGraphs
{
    class Program
    {


        class GeneState
        {
            public string gene;
            public List<string> output = new List<string>();
            public List<string> outputUnfiltered = new List<string>();

            public Dictionary<string, int> countByDisease = new Dictionary<string, int>();

            public int nSingle = 0;
            public int nMultiple = 0;
            public int nInteresting = 0;

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
            var configuration = ASETools.ASEConfirguation.loadFromFile(args);
            if (null == configuration)
            {
                Console.WriteLine("Unable to load configuration.");
            }

            if (configuration.commandLineArgs.Count() > 0)
            {
                if (args.Count() > 1 || args[0] != "--e")
                {
                    Console.WriteLine("usage: GenerateScatterGraphs {--e}");
                    return;
                }
            }

            var useExpression = !configuration.commandLineArgs.Contains("--e");

            var cases = ASETools.Case.LoadCases(configuration.casesFilePathname);
            if (null == cases)
            {
                Console.WriteLine("You must generate cases first.");
                return;
            }

            var selectedGenes = ASETools.SeletedGene.LoadFromFile(configuration.selectedGenesFilename);
            if (null == selectedGenes)
            {
                Console.WriteLine("Unable to load Selected Genes file.");
                return;
            }

            var selectedGeneNames = new List<string>();
            foreach (var selectedGene in selectedGenes)
            {
                selectedGeneNames.Add(selectedGene.Hugo_Symbol.ToLower());
            }

            int nMissingInputs = 0;
            var geneStates = new Dictionary<string, GeneState>();

            foreach (var disease in ASETools.GetListOfDiseases(cases))
            {
                var expressionFilename = configuration.expressionFilesDirectory + "expression_" + disease;
                if (!File.Exists(expressionFilename))
                {
                    Console.WriteLine("Can't find expression file " + expressionFilename + ", skipping " + disease);
                    continue;
                }

                ASETools.ExpressionFile expression = null;

                if (useExpression)
                {
                    expression = new ASETools.ExpressionFile();
                    var timer = new Stopwatch();
                    timer.Start();
                    Console.Write("Loading expression file for " + disease);
                    expression.LoadFromFile(expressionFilename);
                    Console.WriteLine(ASETools.ElapsedTimeInSeconds(timer));
                }

                foreach (var caseEntry in cases.Where(c => c.Value.disease() == disease))
                {
                    var case_ = caseEntry.Value;

                    if (case_.annotated_selected_variants_filename == "" || case_.tumor_rna_mapped_base_count_filename == "")
                    {
                        nMissingInputs++;
                        continue;
                    }

                    var annotatedVariants = ASETools.AnnotatedVariant.readFile(case_.annotated_selected_variants_filename);
                    var mappedBaseCount = ASETools.MappedBaseCount.readFromFile(case_.tumor_rna_mapped_base_count_filename);
                    if (mappedBaseCount == null || annotatedVariants == null)
                    {
                        nMissingInputs++;
                        continue;
                    }

                    foreach (var annotatedVariant in annotatedVariants)
                    {
                        if (!annotatedVariant.somaticMutation)
                        {
                            continue;   // We only care about somatic mutations for this stage, so skip the germline ones.
                        }

                        string Hugo_symbol = annotatedVariant.Hugo_symbol;

                        if (Hugo_symbol == "")
                        {
                            Console.WriteLine("Found empty hugo symbol in annotated variant for somatic mutation on case " + case_.case_id);
                            break;
                        }

                        if (!selectedGeneNames.Contains(Hugo_symbol.ToLower()))
                        {
                            //
                            // It's not a gene we care about, skip it.
                            //
                            continue;
                        }

                        var sex = ASETools.isChromosomeSex(annotatedVariant.contig);

                        GeneState geneState;

                        if (!geneStates.ContainsKey(Hugo_symbol))
                        {
                            geneStates.Add(Hugo_symbol, new GeneState());
                            geneState = geneStates[Hugo_symbol];
                            geneState.gene = Hugo_symbol;
                            geneState.sex = sex;
                            geneState.chrom = annotatedVariant.contig;
                        } else
                        {
                            geneState = geneStates[Hugo_symbol];
                        }

                        if (geneState.sex != sex)
                        {
                            Console.WriteLine("Disagreement on sex for gene " + Hugo_symbol + " this record is " + annotatedVariant.contig + " and first one was " + geneState.chrom);
                            geneState.sex = true;    // Err on the side of being sex linked
                        }

                        if (!geneState.countByDisease.ContainsKey(disease))
                        {
                            geneState.countByDisease.Add(disease, 0);
                        }
                        geneState.countByDisease[disease]++;

                        var nDNAReference = annotatedVariant.tumorDNAReadCounts.nMatchingReference;
                        var nDNATumor = annotatedVariant.tumorDNAReadCounts.nMatchingAlt;
                        var nDNANeither = annotatedVariant.tumorDNAReadCounts.nMatchingNeither;

                        var nRNAReference = annotatedVariant.tumorRNAReadCounts.nMatchingReference;
                        var nRNATumor = annotatedVariant.tumorRNAReadCounts.nMatchingAlt;
                        var nRNANeither = annotatedVariant.tumorRNAReadCounts.nMatchingNeither;

                        int nMutationsThisGene = annotatedVariants.Where(e => e.Hugo_symbol == annotatedVariant.Hugo_symbol).Count(); // Yes, we could me more efficient here and make only one pass, but it really shouldn't matter much

                        if (nMutationsThisGene > 1)
                        {
                            geneState.nMultiple++;
                        }
                        else
                        {
                            geneState.nSingle++;
                        }

                        string outputLine = ASETools.ConvertToExcelString(Hugo_symbol) + "\t" + annotatedVariant.contig + "\t" + annotatedVariant.locus + "\t" + annotatedVariant.variantClassification + "\t" + annotatedVariant.variantType + "\t" +
                            annotatedVariant.reference_allele + "\t" + annotatedVariant.alt_allele + "\t" + disease + "\t" + case_.case_id + "\t" + case_.tumor_dna_file_id + "\t" + case_.tumor_rna_file_id + "\t" + case_.normal_dna_file_id + "\t" + case_.normal_rna_file_id + "\t" +
                            annotatedVariant.normalDNAReadCounts.nMatchingReference + "\t" + annotatedVariant.normalDNAReadCounts.nMatchingAlt + "\t" + annotatedVariant.normalDNAReadCounts.nMatchingNeither + "\t" + annotatedVariant.normalDNAReadCounts.nMatchingBoth + "\t" +
                            annotatedVariant.tumorDNAReadCounts.nMatchingReference + "\t" + annotatedVariant.tumorDNAReadCounts.nMatchingAlt + "\t" + annotatedVariant.tumorDNAReadCounts.nMatchingNeither + "\t" + annotatedVariant.tumorDNAReadCounts.nMatchingBoth + "\t";


                        if (annotatedVariant.normalRNAReadCounts != null)
                        {
                            outputLine += 
                            annotatedVariant.normalRNAReadCounts.nMatchingReference + "\t" + annotatedVariant.normalRNAReadCounts.nMatchingAlt + "\t" + annotatedVariant.normalRNAReadCounts.nMatchingNeither + "\t" + annotatedVariant.normalRNAReadCounts.nMatchingBoth;
                        }
                        else
                        {
                            outputLine += "\t\t\t";
                        }

                        outputLine += "\t" +
                            annotatedVariant.tumorRNAReadCounts.nMatchingReference + "\t" + annotatedVariant.tumorRNAReadCounts.nMatchingAlt + "\t" + annotatedVariant.tumorRNAReadCounts.nMatchingNeither + "\t" + annotatedVariant.tumorRNAReadCounts.nMatchingBoth + "\t" +
                            (nMutationsThisGene > 1) + "\t" + nMutationsThisGene + "\t";

                        geneState.outputUnfiltered.Add(outputLine);

                        if (nRNAReference + nRNATumor + nRNANeither < configuration.nReadsToIncludeVariant || nDNAReference + nDNATumor + nDNANeither < configuration.nReadsToIncludeVariant)
                        {
                            continue;
                        }

                        double tumorDNAFraction = nDNATumor / (nDNATumor + nDNAReference + nDNANeither);
                        double tumorRNAFraction = nRNATumor / (nRNATumor + nRNAReference + nRNANeither);
                        double tumorDNAMultiple = ComputeMultiple(nDNATumor, nDNAReference + nDNANeither);
                        double tumorRNAMultiple = ComputeMultiple(nRNATumor, nRNAReference + nRNANeither);
                        double tumorDNALog = ComputeRatio(nDNATumor, nDNAReference + nDNANeither);
                        double tumorRNALog = ComputeRatio(nRNATumor, nRNAReference + nRNANeither);

                        string ratioOfRatiosString;
                        if (tumorDNALog == 0)
                        {
                            ratioOfRatiosString = "1000";
                        }
                        else
                        {
                            ratioOfRatiosString = Convert.ToString(tumorRNALog / tumorDNALog);
                        }

                        outputLine += "\t" + tumorDNAFraction + "\t" + tumorRNAFraction + "\t" + tumorDNAMultiple + "\t" + tumorRNAMultiple + "\t" + tumorDNALog + "\t" + tumorRNALog + "\t" + ratioOfRatiosString + "\t";

                        //
                        // See if we have the mean & std deviation for this locus, and if so write the z values for normal and mutant expression.  Recall that the mean and std. dev differ by
                        // tumor type, but that we've loaded the values for this tumor type.
                        //
                        ASETools.MeanAndStdDev meanAndStdDev;

                        if (expression != null && mappedBaseCount.mappedBaseCount != 0 && expression.getValue(annotatedVariant.contig, annotatedVariant.locus, out meanAndStdDev))
                        { 
                            outputLine += ((((double)nRNATumor / (double)mappedBaseCount.mappedBaseCount) - meanAndStdDev.mean) / meanAndStdDev.stddev) + "\t" +
                                ((((double)(nRNAReference + nRNANeither) / (double)mappedBaseCount.mappedBaseCount) - meanAndStdDev.mean) / meanAndStdDev.stddev) + "\t" +
                                ((((double)nRNATumor * 2.0 / (double)mappedBaseCount.mappedBaseCount) - meanAndStdDev.mean) / meanAndStdDev.stddev) + "\t" +
                                ((((double)(nRNAReference + nRNANeither) * 2.0 / (double)mappedBaseCount.mappedBaseCount) - meanAndStdDev.mean) / meanAndStdDev.stddev) + "\t";

                            if (meanAndStdDev.mean != 0)
                            {
                                outputLine += ((double)nRNATumor / (double)mappedBaseCount.mappedBaseCount) / meanAndStdDev.mean + "\t" +
                                                ((double)(nRNAReference + nRNANeither) / (double)mappedBaseCount.mappedBaseCount) / meanAndStdDev.mean;
                            }
                            else
                            {
                                outputLine += "\t";
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

                        if (nMutationsThisGene == 1)
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
                    }
                } // foreach case in a given disease
            } // foreach disease


            StreamWriter summaryFile = new StreamWriter(configuration.geneScatterGraphsDirectory +  ASETools.scatterGraphsSummaryFilename);
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

                StreamWriter file = new StreamWriter(configuration.geneScatterGraphsDirectory + geneState.gene + ".txt");
                string unfilteredHeaderLine =
                    "Hugo_Symbol\tChromosome\tStart_Position\tVariant_Classification\tVariant_Type\tReference_Allele\tAlt_Allele\tdisease" +
                    "Case Id\tTumor DNA File ID\tTumor RNA File ID\tNormal DNA File ID\tNormal RNA File ID\t" +
                    "n_normal_DNA_Matching_Reference\tn_normal_DNA_Matching_Alt\tn_normal_DNA_Matching_Neither\tn_normal_DNA_Matching_Both\t" +
                    "n_tumor_DNA_Matching_Reference\tn_tumor_DNA_Matching_Alt\tn_tumor_DNA_Matching_Neither\tn_tumor_DNA_Matching_Both\t" +
                    "n_normal_RNA_Matching_Reference\tn_normal_RNA_Matching_Alt\tn_normal_RNA_Matching_Neither\tn_normal_RNA_Matching_Both\t" +
                    "n_tumor_RNA_Matching_Reference\tn_tumor_RNA_Matching_Alt\tn_tumor_RNA_Matching_Neither\tn_tumor_RNA_Matching_Both\t" +
                    "Multiple Mutations in this Gene\tn Mutations in this gene\t";
                string headerLine = unfilteredHeaderLine + 
                    "\ttumorDNAFraction\ttumorRNAFraction\ttumorDNAMultiple\ttumorRNAMultiple\ttumorDNARatio\ttumorRNARatio\tFractionOfMeanExpression\tzOfmeanExpression\tratioOfRatios\tIsSingle\tCancerType\tzTumor\tzNormal\tz2Tumor\tz2Normal\t%MeanTumor\t%MeanNormal\t";

                file.WriteLine(headerLine + histogramLines[0] + "\t" + histogram2Lines[0]);
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

                file.WriteLine("*done**");
                file.Close();

                var unfilteredFile = new StreamWriter(configuration.geneScatterGraphsDirectory + geneState.gene + ASETools.ASEConfirguation.unfilteredCountsExtention);
                unfilteredFile.WriteLine(unfilteredHeaderLine);
                foreach (var line in geneState.outputUnfiltered)
                {
                    unfilteredFile.WriteLine(line);
                }
                unfilteredFile.WriteLine("**done**");
                unfilteredFile.Close();
            }

            summaryFile.WriteLine("**done**");
            summaryFile.Close();

            Console.WriteLine("" + nMissingInputs + " cases are missing input files.");

        } // Main
    }
}
