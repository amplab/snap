using MathNet.Numerics;
using MathNet.Numerics.Distributions;
using System;
using System.Collections.Generic;
using System.Data.SqlTypes;
using System.Diagnostics;
using System.Diagnostics.Contracts;
using System.Globalization;
using System.IO;
using System.IO.Compression;
using System.Linq;
using System.Net;
using System.Net.Sockets;
using System.Runtime.CompilerServices;
using System.Runtime.InteropServices;
using System.Runtime.Remoting.Metadata.W3cXsd2001;
using System.Runtime.Serialization;
using System.Runtime.Serialization.Json;
using System.Text;
using System.Threading;
using System.Web.UI;

namespace ASELib
{
    public class ASETools
    {
        static ASETools()   // This is the static initializer, which gets called at program start
        {
            tumorToString.Add(true, "Tumor");
            tumorToString.Add(false, "Normal");

            dnaToString.Add(true, "DNA");
            dnaToString.Add(false, "RNA");

            chromosomeSizes.ToList().ForEach(x => chromosomeSizesByName.Add(x.name, x));
            chromosomeSizes.ToList().ForEach(x => chromosomeSizesByName.Add(chromosomeNameToNonChrForm(x.name), x));
        }

        public const string urlPrefix = @"https://api.gdc.cancer.gov/";

        public const int GuidStringLength = 36;

        public const int nHumanNuclearChromosomes = 24;   // 1-22, X and Y.
        public const int nHumanAutosomes = 22;

        public const int Million = 1000000;
        public const int Billion = 1000000000;

        public static readonly string[] chromosomes = { "chr1", "chr2", "chr3", "chr4", "chr5", "chr6", "chr7", "chr8", "chr9", "chr10", "chr11", "chr12", "chr13", "chr14", "chr15", "chr16", "chr17", "chr18", "chr19", "chr20", "chr21", "chr22", "chrx", "chry" };// Leaves off the mitochondrial genes
        public static readonly string[] chromosomesWithMitochondria = { "chr1", "chr2", "chr3", "chr4", "chr5", "chr6", "chr7", "chr8", "chr9", "chr10", "chr11", "chr12", "chr13", "chr14", "chr15", "chr16", "chr17", "chr18", "chr19", "chr20", "chr21", "chr22", "chrx", "chry", "chrm" };
        public static readonly string[] autosomes = { "chr1", "chr2", "chr3", "chr4", "chr5", "chr6", "chr7", "chr8", "chr9", "chr10", "chr11", "chr12", "chr13", "chr14", "chr15", "chr16", "chr17", "chr18", "chr19", "chr20", "chr21", "chr22" };

        // Used for writing out Bonferroni corrected p values in Mann Whitney tests
        public class OutputLine
        {
            public string line;
            public double p;
        }

        static public bool isChromosomeMitochondrial(string chromosome)
        {
            string lowerChromosome = chromosome.ToLower();
            return lowerChromosome == "m" || lowerChromosome == "mt" || lowerChromosome == "chrm" || lowerChromosome == "chrmt";
        }

        static public bool isChromosomeSex(string chromosome)
        {
            string lowerChromosome = chromosome.ToLower();
            return lowerChromosome == "x" || lowerChromosome == "chrx" || lowerChromosome == "y" || lowerChromosome == "chry";
        }

        static public bool isChromosomeAutosomal(string chromosome)
        {
            string lowerChromosome = chromosome.ToLower();
            return autosomes.Contains(lowerChromosome) || autosomes.Contains("chr" + lowerChromosome);
        }

        public static int ChromosomeNameToIndex(string chromosomeName)
        {
            var name = chromosomeName.ToLower();

            for (int i = 1; i <= 22; i++)
            {
                var comp = Convert.ToString(i);
                if (name == comp || name == "chr" + comp) return i;
            }

            //
            // Use 0 and 23 for the sex chromosomes, so the array indices correspond to the chromosome names.
            //
            if (name == "chrx" || name == "x") return 0;
            if (name == "chry" || name == "y") return 23;

            return -1;
        }


        public static string ChromosomeIndexToName(int index, bool useChr = false)
        {
            if (index >= 1 && index <= 22)
            {
                return (useChr ? "chr" : "") + index;
            }

            if (index == 0) return (useChr ? "chr" : "") + "X";

            if (index == 23) return (useChr ? "chr" : "") + "Y";

            return "Invalid chromosome index: " + index;
        }

        public static string chromosomeNameToNonChrForm(string rawName)
        {
            if (rawName.Count() < 4 || rawName.Substring(0, 3).ToLower() != "chr")
            {
                return rawName.ToLower();
            }

            return rawName.Substring(3).ToLower();
        }


        public class MannWhitney<T>
        {
            public delegate bool WhichGroup(T element);
            public delegate double GetValue(T element);
            public static double ComputeMannWhitney(List<T> elements, IComparer<T> comparer, WhichGroup whichGroup, GetValue getValue, out bool enoughData, out bool reversed,
                out double nFirstGroup, out double nSecondGroup, out double U, out double z, bool twoTailed = true, int minGroupSize = 1)
            {
                elements.Sort(comparer);

                reversed = false;

                double RfirstGroup = 0; // The rank sum for the first group
                nFirstGroup = 0;
                nSecondGroup = 0;
                U = 0;
                z = 0;
                int n = elements.Count();

                for (int i = 0; i < n; i++)
                {
                    if (whichGroup(elements[i]))
                    {
                        int cumulativeR = n - i;
                        int nTied = 1;
                        //
                        // Now add in adjascent indices if there's a tie.  For ties, we use the mean of all the indices in the tied region for each of them (regardless of whether they're single or multiple).
                        //
                        for (int j = i - 1; j >= 0 && getValue(elements[j]) == getValue(elements[i]); j--)
                        {
                            cumulativeR += n - j;
                            nTied++;
                        }

                        for (int j = i + 1; j < elements.Count() && getValue(elements[j]) == getValue(elements[i]); j++)
                        {
                            cumulativeR += n - j;
                            nTied++;
                        }


                        RfirstGroup += cumulativeR / nTied;
                        nFirstGroup++;
                    }
                    else
                    {
                        nSecondGroup++;
                    }
                }

                if (nFirstGroup < minGroupSize || nSecondGroup < minGroupSize)
                {
                    //
                    // Not enough data, reject this gene
                    //
                    enoughData = false;
                    return -1;
                }

                U = (double)RfirstGroup - (double)nFirstGroup * (nFirstGroup + 1.0) / 2.0;

                z = (U - nFirstGroup * nSecondGroup / 2) / Math.Sqrt(nSecondGroup * nFirstGroup * (nSecondGroup + nFirstGroup + 1) / 12);

                double p = MathNet.Numerics.Distributions.Normal.CDF(0, 1.0, z);


                if (twoTailed)
                {
                    //
                    // We're doing two-tailed, so we need to see if this is on the other end
                    if (p > 0.5)
                    {
                        p = 1.0 - p;
                        reversed = true;
                    }

                    //
                    // And then multiply by two (because this is a two-tailed test).
                    //
                    p *= 2.0;
                }

                enoughData = true;
                return p;
            }
        } // MannWhitney

        public class WelchsTTest
        {
            public static double T(IEnumerable<double> X1, IEnumerable<double> X2)   // Using 1 and 2 to correspond to the notation in Wikipedia: https://en.wikipedia.org/wiki/Welch%27s_t-test
            {
                var X1bar = MathNet.Numerics.Statistics.Statistics.Mean(X1);
                var X2bar = MathNet.Numerics.Statistics.Statistics.Mean(X2);

                var s1Squared = MathNet.Numerics.Statistics.Statistics.Variance(X1);    // s is standard deviation, so s^2 is variance, which we get directly from the library
                var s2Squared = MathNet.Numerics.Statistics.Statistics.Variance(X2);    // s is standard deviation, so s^2 is variance, which we get directly from the library

                var N1 = X1.Count();
                var N2 = X2.Count();

                return T(X1bar, X2bar, s1Squared, s2Squared, N1, N2);
            }

            static double T(double X1bar, double X2bar, double s1Squared, double s2Squared, int N1, int N2)
            {
                return (X1bar - X2bar) / Math.Sqrt(s1Squared / N1 + s2Squared / N2);
            }

            public static double OneSidedTTest(IEnumerable<double> X1, IEnumerable<double> X2)
            {

                var X1bar = MathNet.Numerics.Statistics.Statistics.Mean(X1);
                var X2bar = MathNet.Numerics.Statistics.Statistics.Mean(X2);

                var s1Squared = MathNet.Numerics.Statistics.Statistics.Variance(X1);    // s is standard deviation, so s^2 is variance, which we get directly from the library
                var s2Squared = MathNet.Numerics.Statistics.Statistics.Variance(X2);    // s is standard deviation, so s^2 is variance, which we get directly from the library

                var N1 = X1.Count();
                var N2 = X2.Count();

                var t = T(X1bar, X2bar, s1Squared, s2Squared, N1, N2);
                var nu = Math.Pow((s1Squared / N1 + s2Squared / N2), 2) / (s1Squared * s1Squared / ((double)N1 * N1 * (N1 - 1) * (N1 - 1)) + s2Squared * s2Squared / ((double)N2 * N2 * (N2 - 1) * (N2 - 1))); // Casting N1 and N2 to double avoids integer overflow

                double p = MathNet.Numerics.Distributions.StudentT.CDF(0, 1, nu, t);

                return p;
            }

            public static double TwoSidedTTest(IEnumerable<double> X1, IEnumerable<double> X2)
            {
                var p = OneSidedTTest(X1, X2);

                if (p > 0.5)
                {
                    p = 1.0 - p;
                }

                p *= 2; // Because 2 tailed

                return p;
            }
        } // WelchsTTest

        public class GeneScatterGraphLine
        {
            public readonly string Hugo_Symbol;
            public readonly string Chromosome;
            public readonly int Start_Position;
            public readonly string Variant_Classification;
            public readonly string Variant_Type;
            public readonly string Reference_Allele;
            public readonly string Alt_Allele;
            public readonly string disease;
            public readonly string case_id;
            public readonly string tumor_dna_file_id;
            public readonly string normal_dna_file_id;
            public readonly string tumor_rna_file_id;
            public readonly string normal_rna_file_id;
            public readonly ReadCounts normalDNAReadCounts;
            public readonly ReadCounts tumorDNAReadCounts;
            public readonly ReadCounts normalRNAReadCounts;
            public readonly ReadCounts tumorRNAReadCounts;
            public readonly int nMutationsThisGene; // In this tumor, of course
            public readonly double tumorDNAFraction;
            public readonly double tumorRNAFraction;
            public readonly double tumorDNAMultiple;
            public readonly double tumorRNAMultiple;
            public readonly double tumorDNARatio;
            public readonly double tumorRNARatio;
            public readonly double ratioOfRatios;
            public readonly bool MultipleMutationsInThisGene;

            public readonly bool zKnown; // The next 6 fields are valid iff zKnown.

            public readonly double zTumor;
            public readonly double zNormal;
            public readonly double z2Tumor;
            public readonly double z2Normal;
            public readonly double percentMeanTumor; // Expressed as a fraction (i.e., 100% -> 1.0)
            public readonly double percentMeanNormal;// Expressed as a fraction (i.e., 100% -> 1.0)

            public readonly bool fromUnfilteredFile; // If this is true, then the fields after the read counts are uninitialized and invalid.

            public readonly string rawInputLine;

            //
            // This is only filled in in the annotated and semi-annotated version.
            //
            public double[] tumorRNAFracAltPercentile = null;  // min, 10th%ile, 20th%ile...90th%ile, max
            public double[] tumorRNAFracRefPercentile = null;
            public double[] tumorRNAFracAllPercentile = null;

            public const string unfilteredHeaderLine = "Hugo_Symbol\tChromosome\tStart_Position\tVariant_Classification\tVariant_Type\tReference_Allele\tAlt_Allele\tdisease\t" +
                    "Case Id\tTumor DNA File ID\tTumor RNA File ID\tNormal DNA File ID\tNormal RNA File ID\t" +
                    "n_normal_DNA_Matching_Reference\tn_normal_DNA_Matching_Alt\tn_normal_DNA_Matching_Neither\tn_normal_DNA_Matching_Both\t" +
                    "n_tumor_DNA_Matching_Reference\tn_tumor_DNA_Matching_Alt\tn_tumor_DNA_Matching_Neither\tn_tumor_DNA_Matching_Both\t" +
                    "n_normal_RNA_Matching_Reference\tn_normal_RNA_Matching_Alt\tn_normal_RNA_Matching_Neither\tn_normal_RNA_Matching_Both\t" +
                    "n_tumor_RNA_Matching_Reference\tn_tumor_RNA_Matching_Alt\tn_tumor_RNA_Matching_Neither\tn_tumor_RNA_Matching_Both\t" +
                    "Multiple Mutations in this Gene\tn Mutations in this gene";

            public const string filteredHeaderLine = unfilteredHeaderLine +
                    "\ttumorDNAFraction\ttumorRNAFraction\ttumorDNAMultiple\ttumorRNAMultiple\ttumorDNARatio\ttumorRNARatio\tRatioOfRatios\tzOftotalExpression\tzTumor\tzNormal\tz2Tumor\tz2Normal\t%MeanTumor\t%MeanNormal\t";

            public const string percentileHeaderLine = filteredHeaderLine +
                "frac min 2ref\tfrac 10th %ile 2ref\tfrac 20th %ile 2ref\tfrac 30th %ile 2ref\tfrac 40th %ile 2ref\tfrac 50th %ile 2ref\tfrac 60th %ile 2ref\tfrac 70th %ile 2ref\tfrac 80th %ile 2ref\tfrac 90th %ile 2ref\tfrac max 2ref\t" +
                "frac min 2alt\tfrac 10th %ile 2alt\tfrac 20th %ile 2alt\tfrac 30th %ile 2alt\tfrac 40th %ile 2alt\tfrac 50th %ile 2alt\tfrac 60th %ile 2alt\tfrac 70th %ile 2alt\tfrac 80th %ile 2alt\tfrac 90th %ile 2alt\tfrac max 2alt\t" +
                "frac min all\tfrac 10th %ile all\tfrac 20th %ile all\tfrac 30th %ile all\tfrac 40th %ile all\tfrac 50th %ile all\tfrac 60th %ile all\tfrac 70th %ile all\tfrac 80th %ile all\tfrac 90th %ile all\tfrac max all";

            public const string annotatedHeaderLine = percentileHeaderLine + "\tase candidate?\twhy not?\twhy not (other than nonsense mediated decay)?\t";



            public string toUnfilteredLine()
            {
                var line = ASETools.ConvertToExcelString(Hugo_Symbol) + "\t" + Chromosome + "\t" + Start_Position + "\t" + Variant_Classification + "\t" + Variant_Type + "\t" +
                    Reference_Allele + "\t" + Alt_Allele + "\t" + disease + "\t" + case_id + "\t" + tumor_dna_file_id + "\t" + tumor_rna_file_id + "\t" + normal_dna_file_id + "\t" + normal_rna_file_id + "\t" +
                    normalDNAReadCounts.nMatchingReference + "\t" + normalDNAReadCounts.nMatchingAlt + "\t" + normalDNAReadCounts.nMatchingNeither + "\t" + normalDNAReadCounts.nMatchingBoth + "\t" +
                    tumorDNAReadCounts.nMatchingReference + "\t" + tumorDNAReadCounts.nMatchingAlt + "\t" + tumorDNAReadCounts.nMatchingNeither + "\t" + tumorDNAReadCounts.nMatchingBoth + "\t";


                if (normalRNAReadCounts != null)
                {
                    line += normalRNAReadCounts.nMatchingReference + "\t" + normalRNAReadCounts.nMatchingAlt + "\t" + normalRNAReadCounts.nMatchingNeither + "\t" + normalRNAReadCounts.nMatchingBoth;
                }
                else
                {
                    line += "\t\t\t";
                }

                line += "\t" +
                    tumorRNAReadCounts.nMatchingReference + "\t" + tumorRNAReadCounts.nMatchingAlt + "\t" + tumorRNAReadCounts.nMatchingNeither + "\t" + tumorRNAReadCounts.nMatchingBoth + "\t" +
                    (nMutationsThisGene > 1) + "\t" + nMutationsThisGene;

                return line;
            }

            public string toFilteredLine()
            {
                var line = toUnfilteredLine() + "\t" + tumorDNAFraction + "\t" + tumorRNAFraction + "\t" + tumorDNAMultiple + "\t" + tumorRNAMultiple + "\t" + tumorDNARatio + "\t" + tumorRNARatio + "\t" + ratioOfRatios + "\t";

                if (zKnown)
                {
                    line += zTumor + "\t" + zNormal + "\t" + z2Tumor + "\t" + z2Normal + "\t" + percentMeanTumor + "\t" + percentMeanNormal + "\t";

                }
                else
                {
                    line += "\t\t\t\t\t\t";
                }

                return line;
            }

            static string lineForPercentiles(double[] percentiles)
            {
                string line = "";

                for (int i = 0; i < 11; i++)
                {
                    if (percentiles[i] == double.NegativeInfinity)
                    {
                        line += "\t*";
                    }
                    else
                    {
                        line += "\t" + percentiles[i];
                    }
                }
                return line;
            }

            public string toPercentileLine()
            {
                return toFilteredLine() + lineForPercentiles(tumorRNAFracRefPercentile) + lineForPercentiles(tumorRNAFracAltPercentile) + lineForPercentiles(tumorRNAFracAllPercentile);
            }




            public int CompareByDiseaseAndChromosome(GeneScatterGraphLine peer)
            {
                if (disease != peer.disease)
                {
                    return disease.CompareTo(peer.disease);
                }

                return (Chromosome.CompareTo(peer.Chromosome));
            }

            GeneScatterGraphLine(string Hugo_Symbol_, string Chromosome_, int Start_Position_, string Variant_Classification_, string Variant_Type_, string Reference_Allele_, string Alt_Allele_, string disease_,
                string case_id_, string normal_dna_file_id_, string tumor_dna_file_id_, string normal_rna_file_id_, string tumor_rna_file_id_, ReadCounts normalDNAReadCounts_, ReadCounts tumorDNAReadCounts_,
                ReadCounts normalRNAReadCounts_, ReadCounts tumorRNAReadCounts_, bool MultipleMutationsInThisGene_, int nMutationsThisGene_, double tumorDNAFraction_, double tumorRNAFraction_, double tumorDNAMultiple_, double tumorRNAMultiple_, double tumorDNARatio_,
                double tumorRNARatio_, double ratioOfRatios_, bool zKnown_, double zTumor_, double zNormal_, double z2Tumor_, double z2Normal_, double percentMeanTumor_,
                double percentMeanNormal_, bool fromUnfilteredFile_, string rawInputLine_, double[] refPercentiles_, double[] altPercentiles_, double[] allPercentiles_)
            {
                Hugo_Symbol = Hugo_Symbol_;
                Chromosome = Chromosome_;
                Start_Position = Start_Position_;
                Variant_Classification = Variant_Classification_;
                Variant_Type = Variant_Type_;
                Reference_Allele = Reference_Allele_;
                Alt_Allele = Alt_Allele_;
                disease = disease_;
                case_id = case_id_;
                normal_dna_file_id = normal_dna_file_id_;
                tumor_dna_file_id = tumor_dna_file_id_;
                normal_rna_file_id = normal_rna_file_id_;
                tumor_rna_file_id = tumor_rna_file_id_;
                normalDNAReadCounts = normalDNAReadCounts_;
                tumorDNAReadCounts = tumorDNAReadCounts_;
                normalRNAReadCounts = normalRNAReadCounts_;
                tumorRNAReadCounts = tumorRNAReadCounts_;
                MultipleMutationsInThisGene = MultipleMutationsInThisGene_;
                nMutationsThisGene = nMutationsThisGene_;
                tumorDNAFraction = tumorDNAFraction_;
                tumorRNAFraction = tumorRNAFraction_;
                tumorDNAMultiple = tumorDNAMultiple_;
                tumorRNAMultiple = tumorRNAMultiple_;
                tumorDNARatio = tumorDNARatio_;
                tumorRNARatio = tumorRNARatio_;
                ratioOfRatios = ratioOfRatios_;
                zKnown = zKnown_;
                zTumor = zTumor_;
                zNormal = zNormal_;
                z2Tumor = z2Tumor_;
                z2Normal = z2Normal_;
                percentMeanTumor = percentMeanTumor_;
                percentMeanNormal = percentMeanNormal_;
                fromUnfilteredFile = fromUnfilteredFile_;
                rawInputLine = rawInputLine_;
                tumorRNAFracRefPercentile = refPercentiles_;
                tumorRNAFracAltPercentile = altPercentiles_;
                tumorRNAFracAllPercentile = allPercentiles_;
            }

            static string[] wantedFieldsArrayUnfilteredVersion =
 {
                    "Hugo_Symbol",
                    "Chromosome",
                    "Start_Position",
                    "Variant_Classification",
                    "Variant_Type",
                    "Reference_Allele",
                    "Alt_Allele",
                    "disease",
                    "Case Id",
                    "Tumor DNA File ID",
                    "Tumor RNA File ID",
                    "Normal DNA File ID",
                    "Normal RNA File ID",
                    "n_normal_DNA_Matching_Reference",
                    "n_normal_DNA_Matching_Alt",
                    "n_normal_DNA_Matching_Neither",
                    "n_normal_DNA_Matching_Both",
                    "n_tumor_DNA_Matching_Reference",
                    "n_tumor_DNA_Matching_Alt",
                    "n_tumor_DNA_Matching_Neither",
                    "n_tumor_DNA_Matching_Both",
                    "n_normal_RNA_Matching_Reference",
                    "n_normal_RNA_Matching_Alt",
                    "n_normal_RNA_Matching_Neither",
                    "n_normal_RNA_Matching_Both",
                    "n_tumor_RNA_Matching_Reference",
                    "n_tumor_RNA_Matching_Alt",
                    "n_tumor_RNA_Matching_Neither",
                    "n_tumor_RNA_Matching_Both",
                    "Multiple Mutations in this Gene",
                    "n Mutations in this gene",
            };

            static string[] wantedFieldsArrayFieldsAdditionalForFiltered =
            {
                    "tumorDNAFraction",
                    "tumorRNAFraction",
                    "tumorDNAMultiple",
                    "tumorRNAMultiple",
                    "tumorDNARatio",
                    "tumorRNARatio",
                    "RatioOfRatios",
                    "zTumor",
                    "zNormal",
                    "z2Tumor",
                    "z2Normal",
                    "%MeanTumor",
                    "%MeanNormal",
            };

            static string[] wantedFieldsArrayFieldsAdditionalForPercentiles =
            {
                @"frac min 2ref",
                @"frac 10th %ile 2ref",
                @"frac 20th %ile 2ref",
                @"frac 30th %ile 2ref",
                @"frac 40th %ile 2ref",
                @"frac 50th %ile 2ref",
                @"frac 60th %ile 2ref",
                @"frac 70th %ile 2ref",
                @"frac 80th %ile 2ref",
                @"frac 90th %ile 2ref",
                @"frac max 2ref",
                @"frac min 2alt",
                @"frac 10th %ile 2alt",
                @"frac 20th %ile 2alt",
                @"frac 30th %ile 2alt",
                @"frac 40th %ile 2alt",
                @"frac 50th %ile 2alt",
                @"frac 60th %ile 2alt",
                @"frac 70th %ile 2alt",
                @"frac 80th %ile 2alt",
                @"frac 90th %ile 2alt",
                @"frac max 2alt",
                @"frac min all",
                @"frac 10th %ile all",
                @"frac 20th %ile all",
                @"frac 30th %ile all",
                @"frac 40th %ile all",
                @"frac 50th %ile all",
                @"frac 60th %ile all",
                @"frac 70th %ile all",
                @"frac 80th %ile all",
                @"frac 90th %ile all",
                @"frac max all"
            };

            public static List<GeneScatterGraphLine> LoadAllGeneScatterGraphLines(string directoryName, bool fromUnfiltered, string hugoSymbol /* this may be * to load all*/, bool includePercentiles = false)
            {
                var geneScatterGraphEntries = new List<GeneScatterGraphLine>();

                var wantedFieldsFilteredVersion = wantedFieldsArrayUnfilteredVersion.ToList();
                wantedFieldsFilteredVersion.AddRange(wantedFieldsArrayFieldsAdditionalForFiltered.ToList());

                var wantedFields = fromUnfiltered ? wantedFieldsArrayUnfilteredVersion.ToList() : wantedFieldsFilteredVersion;

                foreach (var filename in Directory.EnumerateFiles(directoryName, hugoSymbol + (fromUnfiltered ? Configuration.unfilteredCountsExtention : ".txt")))
                {
                    if (filename.Contains(ASETools.annotated_scatter_graph_filename_extension) || filename.Contains(annotated_scatter_graphs_histogram_filename) || filename.Contains(raw_median_data_extension))
                    {
                        continue;    // Skip the annotated ones for now.
                    }

                    if (filename.Count() == 0 || GetFileNameFromPathname(filename)[0] == '_')
                    {
                        continue;   // Summary file like _MannWhitney rather than a gene file
                    }

                    if (!fromUnfiltered && filename.Contains(Configuration.unfilteredCountsExtention))
                    {
                        continue;   // This happens when hugo_symbol is *
                    }

                    var inputFile = CreateStreamReaderWithRetry(filename);
                    if (null == inputFile)
                    {
                        Console.WriteLine("Unable to open " + filename);
                        throw new FileNotFoundException();
                    }
                    var headerizedFile = new HeaderizedFile<GeneScatterGraphLine>(inputFile, false, true, "", wantedFields, inputFilename_: filename);

                    List<GeneScatterGraphLine> linesFromThisFile;

                    headerizedFile.ParseFile(x => ParseLine(x, fromUnfiltered, includePercentiles), out linesFromThisFile);

                    geneScatterGraphEntries.AddRange(linesFromThisFile);
                }

                return geneScatterGraphEntries;
            }

            public static List<GeneScatterGraphLine> LoadAllLinesFromPercentilesDirectory(string directoryName)
            {
                var geneScatterGraphEntries = new List<GeneScatterGraphLine>();

                var wantedFields = wantedFieldsArrayUnfilteredVersion.ToList();
                wantedFields.AddRange(wantedFieldsArrayFieldsAdditionalForFiltered.ToList());
                wantedFields.AddRange(wantedFieldsArrayFieldsAdditionalForPercentiles.ToList());

                foreach (var filename in Directory.EnumerateFiles(directoryName, GeneScatterGraphLinesWithPercentilesPrefix + "*"))
                {
                    if (filename == "")
                    {
                        continue;
                    }

                    var inputFile = CreateStreamReaderWithRetry(filename);
                    if (null == inputFile)
                    {
                        throw new Exception("Unable to open input file " + filename);
                    }

                    var headerizedFile = new HeaderizedFile<GeneScatterGraphLine>(inputFile, false, true, "", wantedFields, inputFilename_: filename);
                    List<GeneScatterGraphLine> linesFromThisFile;

                    headerizedFile.ParseFile(x => ParseLine(x, false, true), out linesFromThisFile);
                    geneScatterGraphEntries.AddRange(linesFromThisFile);
                } // each file

                return geneScatterGraphEntries;
            }

            static GeneScatterGraphLine ParseLine(ASETools.HeaderizedFile<GeneScatterGraphLine>.FieldGrabber fieldGrabber, bool fromUnfilteredFile, bool includePercentiles)
            {
                double[] refPercentiles = null;
                double[] altPercentiles = null;
                double[] allPercentiles = null;

                ReadCounts normalRNAReadCounts;

                if (fieldGrabber.AsString("n_normal_RNA_Matching_Reference") != "")
                {
                    normalRNAReadCounts = new ReadCounts(fieldGrabber.AsInt("n_normal_RNA_Matching_Reference"), fieldGrabber.AsInt("n_normal_RNA_Matching_Alt"), fieldGrabber.AsInt("n_normal_RNA_Matching_Neither"),
                        fieldGrabber.AsInt("n_normal_RNA_Matching_Both"));
                }
                else
                {
                    normalRNAReadCounts = null;
                }

                double tumorDNAFraction = -1, tumorRNAFraction = -1, tumorDNAMultiple = -1, tumorRNAMultiple = -1, tumorDNARatio = -1, tumorRNARatio = -1, ratioOfRatios = -1, zTumor = -1, zNormal = -1, z2Tumor = -1, z2Normal = -1, percentMeanTumor = -1, percentMeanNormal = -1;
                bool zKnown = false;
                if (!fromUnfilteredFile) {
                    tumorDNAFraction = fieldGrabber.AsDouble("tumorDNAFraction");
                    tumorRNAFraction = fieldGrabber.AsDouble("tumorRNAFraction");
                    tumorDNAMultiple = fieldGrabber.AsDouble("tumorDNAMultiple");
                    tumorRNAMultiple = fieldGrabber.AsDouble("tumorRNAMultiple");
                    tumorDNARatio = fieldGrabber.AsDouble("tumorDNARatio");
                    tumorRNARatio = fieldGrabber.AsDouble("tumorRNARatio");
                    ratioOfRatios = fieldGrabber.AsDouble("RatioOfRatios");
                    if (fieldGrabber.AsString("zTumor") != "")
                    {
                        zKnown = true;
                        zTumor = fieldGrabber.AsDouble("zTumor");
                        zNormal = fieldGrabber.AsDouble("zNormal");
                        z2Tumor = fieldGrabber.AsDouble("z2Tumor");
                        z2Normal = fieldGrabber.AsDouble("z2Normal");
                        percentMeanTumor = fieldGrabber.AsDouble("%MeanTumor");
                        percentMeanNormal = fieldGrabber.AsDoubleNegativeInfinityIfStarOrEmptyString("%MeanNormal");
                    }
                }

                if (includePercentiles)
                {
                    refPercentiles = LoadPercentileSet(fieldGrabber, "2ref");
                    altPercentiles = LoadPercentileSet(fieldGrabber, "2alt");
                    allPercentiles = LoadPercentileSet(fieldGrabber, "all");
                }

                return new GeneScatterGraphLine(
                  fieldGrabber.AsString("Hugo_Symbol"), fieldGrabber.AsString("Chromosome"), fieldGrabber.AsInt("Start_Position"), fieldGrabber.AsString("Variant_Classification"), fieldGrabber.AsString("Variant_Type"),
                    fieldGrabber.AsString("Reference_Allele"), fieldGrabber.AsString("Alt_Allele"), fieldGrabber.AsString("disease"), fieldGrabber.AsString("Case Id"), fieldGrabber.AsString("Normal DNA File ID"), fieldGrabber.AsString("Tumor DNA File ID"),
                    fieldGrabber.AsString("Normal RNA File ID"), fieldGrabber.AsString("Tumor RNA File ID"),
                    new ReadCounts(fieldGrabber.AsInt("n_normal_DNA_Matching_Reference"), fieldGrabber.AsInt("n_normal_DNA_Matching_Alt"), fieldGrabber.AsInt("n_normal_DNA_Matching_Neither"),
                        fieldGrabber.AsInt("n_normal_DNA_Matching_Both")),
                    new ReadCounts(fieldGrabber.AsInt("n_tumor_DNA_Matching_Reference"), fieldGrabber.AsInt("n_tumor_DNA_Matching_Alt"), fieldGrabber.AsInt("n_tumor_DNA_Matching_Neither"),
                        fieldGrabber.AsInt("n_tumor_DNA_Matching_Both")),
                    normalRNAReadCounts,
                    new ReadCounts(fieldGrabber.AsInt("n_tumor_RNA_Matching_Reference"), fieldGrabber.AsInt("n_tumor_RNA_Matching_Alt"), fieldGrabber.AsInt("n_tumor_RNA_Matching_Neither"),
                        fieldGrabber.AsInt("n_tumor_RNA_Matching_Both")),
                    fieldGrabber.AsBool("Multiple Mutations in this Gene"), fieldGrabber.AsInt("n Mutations in this gene"), tumorDNAFraction, tumorRNAFraction, tumorDNAMultiple, tumorRNAMultiple, tumorDNARatio, tumorRNARatio,
                    ratioOfRatios, zKnown, zTumor, zNormal, z2Tumor, z2Normal, percentMeanTumor, percentMeanNormal, fromUnfilteredFile, fieldGrabber.rawLine(), refPercentiles, altPercentiles, allPercentiles);
            } // ParseLine

            static double[] LoadPercentileSet(ASETools.HeaderizedFile<GeneScatterGraphLine>.FieldGrabber fieldGrabber, string fieldSuffix)
            {
                double[] percentiles = new double[11];

                percentiles[0] = fieldGrabber.AsDoubleNegativeInfinityIfStarOrEmptyString("frac min " + fieldSuffix);

                for (int i = 1; i < 10; i++)
                {
                    percentiles[i] = fieldGrabber.AsDoubleNegativeInfinityIfStarOrEmptyString("frac " + i + "0th %ile " + fieldSuffix);
                }

                percentiles[10] = fieldGrabber.AsDoubleNegativeInfinityIfStarOrEmptyString("frac max " + fieldSuffix);

                return percentiles;
            }

            public bool isASECandidate(List<CopyNumberVariation> copyNumber, Configuration configuration, Dictionary<string, ASEMapPerGeneLine> perGeneASEMap, GeneMap geneMap, ASERepetitiveRegionMap repetitiveRegionMap)
            {
                string whyNot;
                return isASECandidate(out whyNot, copyNumber, configuration, perGeneASEMap, geneMap, repetitiveRegionMap);
            }

            public bool isASECandidate(out string whyNot, List<CopyNumberVariation> copyNumber, Configuration configuration, Dictionary<string, ASEMapPerGeneLine> perGeneASEMap, GeneMap geneMap, ASERepetitiveRegionMap repetitiveRegionMap,
                bool ignoreNonseMediatedDecay = false)

            {
                if (NonsenseMediatedDecayCausingVariantClassifications.Contains(Variant_Classification) && !ignoreNonseMediatedDecay)
                {
                    nNonsenseMediatedDecay++;
                    whyNot = "Nonsense mediated decay";
                    return false;
                }

                if (tumorDNAReadCounts == null || tumorRNAReadCounts == null)
                {
                    nNoReadCounts++;
                    whyNot = "Missing read counts";
                    return false;
                }

                if (!checkReadCountsForASECandidacy(out whyNot, tumorDNAReadCounts, tumorRNAReadCounts, configuration.minRNAReadCoverage, configuration.minDNAReadCoverage))
                {
                    nBadReadCounts++;
                    return false;
                }

                //
                // If this is in a gene that has too much ASE in the matched normals, then it's not an ASE candidate.
                //
                if (perGeneASEMap != null)
                {
                    //
                    // If any genes that contain this locus have ASE > 0.5, then reject this place.
                    //
                    foreach (var gene in geneMap.getGenesMappedTo(Chromosome, Start_Position))
                    {
                        if (perGeneASEMap.ContainsKey(gene.hugoSymbol) && perGeneASEMap[gene.hugoSymbol].sampleData[false].meanASE >= configuration.ASEInNormalAtWhichToExcludeGenes)
                        {
                            whyNot = "Gene has too high overall ASE";
                            nBadGene++;
                            return false;
                        }
                    }
                }

                if (null != repetitiveRegionMap && repetitiveRegionMap.isCloseToRepetitiveRegion(ChromosomeNameToIndex(Chromosome), Start_Position, 150))
                {
                    whyNot = "In a repetitive region";
                    nRepetitive++;
                    return false;
                }

                if (copyNumber == null)
                {
                    whyNot = "ASE candidate";
                    nCandidate++;
                    return true;
                }

                var overlappingCNV = copyNumber.Where(cnv =>
                    cnv.OverlapsLocus(ASETools.chromosomeNameToNonChrForm(Chromosome),
                    Start_Position, Start_Position + 1)).Where(r => Math.Abs(r.Segment_Mean) > 1.0).ToList();

                if (overlappingCNV.Count() == 0)
                {
                    whyNot = "ASE candidate";
                    nCandidate++;
                    return true;
                } else
                {
                    whyNot = "Copy number variant";
                    nBadCopyNumber++;
                    return false;
                }
            } // isASECandidate
        } // GeneScatterGraphLine

        public static int nNonsenseMediatedDecay = 0, nNoReadCounts = 0, nBadReadCounts = 0, nBadGene = 0, nNoCopyNumber = 0, nBadCopyNumber = 0, nRepetitive = 0, nCandidate = 0; // BJB


        static bool checkReadCountsForASECandidacy(ReadCounts DNAReadCounts, ReadCounts RNAReadCounts, int minDNACoverage, int minRNACoverage)
        {
            string whyNot;
            return checkReadCountsForASECandidacy(out whyNot, DNAReadCounts, RNAReadCounts, minDNACoverage, minRNACoverage);
        }

        static bool checkReadCountsForASECandidacy(out string whyNot, ReadCounts DNAReadCounts, ReadCounts RNAReadCounts, int minDNACoverage, int minRNACoverage)
        {
            // check there is sufficient coverage for DNA and RNA
            // check annotated variant is not a possible minor subclone
            if (!(DNAReadCounts.nMatchingReference + DNAReadCounts.nMatchingAlt >= minDNACoverage))
            {
                whyNot = "Low DNA Read Count";
                return false;
            }

            if (!(RNAReadCounts.nMatchingReference + RNAReadCounts.nMatchingAlt >= minRNACoverage))
            {
                whyNot = "Low RNA read count";
                return false;
            }

            if (!(DNAReadCounts.nMatchingReference * 3 >= DNAReadCounts.nMatchingAlt * 2))
            {
                whyNot = "Low reference DNA";
                return false;
            }

            if (!(DNAReadCounts.nMatchingAlt * 3 >= DNAReadCounts.nMatchingReference * 2))
            {
                whyNot = "High reference DNA";
                return false;
            }

            if (!((DNAReadCounts.nMatchingNeither + DNAReadCounts.nMatchingBoth) * 20 < DNAReadCounts.totalReads()))
            {
                whyNot = "Too many funny DNA reads";
                return false;
            }

            if (!((RNAReadCounts.nMatchingNeither + RNAReadCounts.nMatchingBoth) * 20 < RNAReadCounts.totalReads())) {
                whyNot = "Too many funny RNA reads";
                return false;
            }

            whyNot = "ASE candidate";
            return true;
        }

        public static double MeanOfList(List<double> values)
        {
            return values.Sum() / values.Count();
        }

        public static double StandardDeviationOfList(List<double> values)
        {
            var mean = MeanOfList(values);

            double variance = 0;
            foreach (var value in values)
            {
                var difference = value - mean;
                variance = difference * difference;
            }

            return Math.Sqrt(variance);
        }


        public class Exon : IEquatable<Exon>    // IEquatable allows this to be a dictionary key.
        {
            public Exon(string startString, string endString)
            {
                start = Convert.ToInt32(startString);
                end = Convert.ToInt32(endString);
            }

            public int start;
            public int end;

            public bool Equals(Exon other)
            {
                return (start == other.start && end == other.end);
            }

            public int length()
            {
                return end - start;
            }

            public override int GetHashCode()
            {
                return start ^ (end << 4);  // 4 is the max shift that won't lose any bits for the largest chromosome offset, which is the end of chr1 at ~249Mb = 0x0e4e1c00
            }
        }

        public class Isoform    // These come from the knownGene file
        {
            public string hugo_symbol;
            public string ucscId;
            public string chromosome;   // The non-chr version
            public string strand;
            public int txStart;
            public int txEnd;
            public int cdsStart;
            public int cdsEnd;
            public string proteinId;
            public string alignId;
            public Exon[] exons;

            public static Isoform fromFileLine(string fileLine)
            {
                var fields = fileLine.Split('\t');
                if (fields.Count() < 12)
                {
                    Console.WriteLine("Isoform.fromFileLine: line had wrong number of fields, " + fields.Count() + " != 12");
                    return null;
                }

                var isoform = new Isoform();

                isoform.hugo_symbol = ConvertToNonExcelString(fields[16]);
                isoform.ucscId = ConvertToNonExcelString(fields[0]);
                isoform.chromosome = ConvertToNonExcelString(fields[1]);
                isoform.strand = ConvertToNonExcelString(fields[2]);

                int exonCount;
                try
                {
                    isoform.txStart = Convert.ToInt32(ConvertToNonExcelString(fields[3]));
                    isoform.txEnd = Convert.ToInt32(ConvertToNonExcelString(fields[4]));

                    if (isoform.txEnd <= isoform.txStart)
                    {
                        Console.WriteLine("Isoform.fromFileLine: warning: isoform " + isoform.ucscId + " has empty or negative transcription region " + fileLine);
                    }

                    isoform.cdsStart = Convert.ToInt32(ConvertToNonExcelString(fields[5]));
                    isoform.cdsEnd = Convert.ToInt32(ConvertToNonExcelString(fields[6]));
                    exonCount = Convert.ToInt32(ConvertToNonExcelString(fields[7]));


                    var exonStartStrings = ConvertToNonExcelString(fields[8]).Split(',');
                    var exonEndStrings = ConvertToNonExcelString(fields[9]).Split(',');

                    //
                    // They have trailing commas, so they should have one more field than there are exons, and also should have an empty string
                    // as their last element.
                    //
                    if (exonStartStrings.Count() != exonCount + 1 || exonEndStrings.Count() != exonCount + 1 || exonStartStrings[exonCount] != "" || exonEndStrings[exonCount] != "")
                    {
                        Console.WriteLine("Isoform.fromFileLine: Bad exon start/end: " + fileLine);
                    }

                    isoform.exons = new Exon[exonCount];
                    for (int i = 0; i < exonCount; i++)
                    {
                        isoform.exons[i] = new Exon(exonStartStrings[i], exonEndStrings[i]);
                    }


                }
                catch (FormatException)
                {
                    Console.WriteLine("Isoform.fromFileLine: Format exception parsing a numeric field in line: " + fileLine);
                    return null;
                }

                isoform.proteinId = fields[10];
                isoform.alignId = fields[11];

                return isoform;
            }

            public int codingSize()
            {
                return exons.Select(x => x.end - x.start).Sum();
            }
        } // Isoform


        // Holds Gene location and isoform information
        public class GeneLocationInfo
        {
            public string hugoSymbol;
            public string chromosome;   // in non-chr form
            public int minLocus;
            public int maxLocus;

            public bool inconsistent = false;

            public bool overlapsRange(string locusChromosome, int start, int end)
            {
                if (start < maxLocus && end > minLocus)
                    return true;
                else
                    return false;
            }

            public bool containsLocus(string locusChromosome, int locus)
            {
                return chromosome == locusChromosome && minLocus <= locus && maxLocus >= locus;
            }

            public bool containsLocusInCodingAndKnownExpressionRegion(string locusChromosome, int locus, ExpressionFile expressionFile)
            {
                return expressionFile.valueKnown(locusChromosome, locus) && isoforms.Any(isoform => isoform.chromosome == locusChromosome && isoform.exons.Any(exon => exon.start <= locus && exon.end >= locus));
            }

            public int size()
            {
                return maxLocus - minLocus;
            }

            public int codingSizeOfLargestIsoform() // The sum of the exons of the largest coding isoform
            {
                if (isoforms.Count() == 0)
                {
                    return 0;
                }

                return isoforms.Select(x => x.codingSize()).ToList().Max();
            }

            public int totalBasesInCodingAndKnownExpressionRegion(ExpressionFile expressionFile) // Just what it says.
            {
                int basesInCodingRegion = 0;
                for (int locus = minLocus; locus <= maxLocus; locus++)
                {
                    if (containsLocusInCodingAndKnownExpressionRegion(chromosome, locus, expressionFile))
                    {
                        basesInCodingRegion++;
                    }
                }

                return basesInCodingRegion;
            }

            public int minDistanceFromTSS(string locusChromosome, int locus)    // How far away is this locus from any transcription start site in this gene (different isoforms may have different TSSs, take the closest one).
            {
                if (chromosome != locusChromosome)
                    return int.MaxValue;

                var nIsoforms = isoforms.Count();

                int minValue = int.MaxValue;
                for (int i = 0; i < nIsoforms; i++)
                {
                    minValue = Math.Min(minValue, Math.Abs(isoforms[i].txStart - locus));
                }

                return minValue;
                // This seemed absurdly slow, so it's now just an explicit for loop.  I think maybe the runtime is taking an unnecessary lock here or something.   return isoforms.Select(x => Math.Abs(x.txStart - locus)).ToList().Min();
            }

            public int minDistanceFromTSS(string locusChromosome, int start, int end)
            {
                return Math.Min(minDistanceFromTSS(locusChromosome, start), minDistanceFromTSS(locusChromosome, end));
            }

            public List<Isoform> isoforms = new List<Isoform>();
        } // GeneLocationInfo

        // Stores a dictionary of genes for each chromosome
        public class GeneLocationsByNameAndChromosome
        {
            public GeneLocationsByNameAndChromosome(Dictionary<string, GeneLocationInfo> genesByName_)
            {
                genesByName = genesByName_;

                foreach (var geneEntry in genesByName)
                {
                    var gene = geneEntry.Value;
                    if (!genesByChromosome.ContainsKey(gene.chromosome))
                    {
                        genesByChromosome.Add(gene.chromosome, new List<GeneLocationInfo>());
                    }

                    genesByChromosome[gene.chromosome].Add(gene);
                }
            }

            public Dictionary<string, GeneLocationInfo> genesByName;
            public Dictionary<string, List<GeneLocationInfo>> genesByChromosome = new Dictionary<string, List<GeneLocationInfo>>();    // chromosome is in non-chr form.
        } // GeneLocationsByNameAndChromosome

        public class GeneMap
        {
            public GeneMap(Dictionary<string, GeneLocationInfo> genesByName)
            {
                var key = new Range();

                foreach (var geneEntry in genesByName)
                {
                    var gene = geneEntry.Value;

                    if (gene.inconsistent)
                    {
                        continue;
                    }

                    key.chromosome = gene.chromosome;
                    key.minLocus = gene.minLocus;

                    //
                    // We need to split any regions that cross the beginning and end of the gene, and add the gene to any regions entirely
                    // contained within it.
                    //

                    int handledUpTo = gene.minLocus;

                    Range overlapRange;
                    if (map.FindFirstLessThanOrEqualTo(key, out overlapRange) && overlapRange.chromosome == gene.chromosome && overlapRange.minLocus < gene.minLocus && overlapRange.maxLocus >= gene.minLocus)
                    {
                        //
                        // We overlap at the beginning.  Delete the existing range and split it in two, with the first half not containing us and the second half containing us.
                        //
                        map.Delete(overlapRange);
                        Range firstHalf = new Range();
                        firstHalf.chromosome = gene.chromosome;
                        firstHalf.minLocus = overlapRange.minLocus;
                        firstHalf.maxLocus = gene.minLocus - 1;
                        firstHalf.genes = overlapRange.genes;

                        map.Insert(firstHalf);

                        Range secondHalf = new Range();
                        secondHalf.chromosome = gene.chromosome;
                        secondHalf.minLocus = gene.minLocus;
                        secondHalf.maxLocus = Math.Min(overlapRange.maxLocus, gene.maxLocus);
                        overlapRange.genes.ForEach(x => secondHalf.genes.Add(x));
                        secondHalf.genes.Add(gene);

                        map.Insert(secondHalf);
                        handledUpTo = secondHalf.maxLocus + 1;

                        if (handledUpTo < overlapRange.maxLocus)
                        {
                            //
                            // The new gene fit entirely within an old one.  We need to add the end of the overlap range.
                            //
                            Range thirdHalf = new Range();
                            thirdHalf.chromosome = gene.chromosome;
                            thirdHalf.minLocus = handledUpTo;
                            thirdHalf.maxLocus = overlapRange.maxLocus;
                            overlapRange.genes.ForEach(x => thirdHalf.genes.Add(x));

                            map.Insert(thirdHalf);
                            handledUpTo = thirdHalf.maxLocus + 1;
                        }
                    }

                    //
                    // Now cruise along looking at ranges that are entirely within the gene.
                    //

                    while (handledUpTo < gene.maxLocus)
                    {
                        key.minLocus = handledUpTo;

                        Range containedRange;
                        bool foundContainedRange;
                        if ((foundContainedRange = map.FindFirstGreaterThanOrEqualTo(key, out containedRange)) && containedRange.minLocus == handledUpTo && containedRange.chromosome == gene.chromosome)
                        {
                            //
                            // Two cases: the contained range is entirely within gene, and the contained range hangs off the end of gene.
                            // In the first case, we just add ourself to the list.  In the second, we need to split.
                            //
                            if (containedRange.maxLocus <= gene.maxLocus)
                            {
                                containedRange.genes.Add(gene);
                            } else
                            {
                                //
                                // Split.  This is the symmetric case for the split at the beginning of the gene above this loop.
                                //
                                var newRange = new Range();
                                newRange.minLocus = containedRange.minLocus;
                                newRange.chromosome = containedRange.chromosome;
                                newRange.maxLocus = gene.maxLocus;
                                containedRange.genes.ForEach(x => newRange.genes.Add(x));
                                newRange.genes.Add(gene);

                                map.Delete(containedRange);
                                map.Insert(newRange);

                                containedRange.minLocus = gene.maxLocus + 1;
                                map.Insert(containedRange);
                            }
                            handledUpTo = containedRange.maxLocus + 1;
                        } else
                        {
                            //
                            // This is an empty range.  Create a new one that goes up to either the end of gene or the beginnig of the range we found.
                            //
                            var newRange = new Range();
                            newRange.chromosome = gene.chromosome;
                            newRange.minLocus = handledUpTo;
                            if (foundContainedRange && containedRange.minLocus <= gene.maxLocus && containedRange.chromosome == gene.chromosome)
                            {
                                newRange.maxLocus = containedRange.minLocus - 1;
                            } else
                            {
                                newRange.maxLocus = gene.maxLocus;
                            }
                            newRange.genes.Add(gene);

                            handledUpTo = newRange.maxLocus + 1;

                            map.Insert(newRange);
                        }
                    } // While we have genome to cover.
                } // foreach gene
            } // GeneMap constructor

            class Range : IComparable
            {
                public string chromosome;
                public int minLocus;
                public int maxLocus;

                public int CompareTo(object peerObject)
                {
                    var peer = (Range)peerObject;

                    var chromosomesComparison = chromosome.CompareTo(peer.chromosome);
                    if (chromosomesComparison != 0) return chromosomesComparison;

                    return minLocus.CompareTo(peer.minLocus);
                }

                public List<GeneLocationInfo> genes = new List<GeneLocationInfo>();
            } // Range

            public List<GeneLocationInfo> getGenesMappedTo(string chromosome, int locus)
            {
                var key = new Range();
                key.chromosome = chromosome;
                key.minLocus = locus;

                Range mapRange;
                if (!map.FindFirstLessThanOrEqualTo(key, out mapRange) || mapRange.chromosome != chromosome || mapRange.maxLocus < locus)
                {
                    return new List<GeneLocationInfo>();
                }

                return mapRange.genes;
            } // getGenesMappedTo

            //
            // Figure out the maximum range around locus that is in any gene that includes locus.  Return false iff locus isn't
            // in any gene.
            //
            public bool encompassingGeneRange(string chromosome, int locus, out int minEncompassingLocus, out int maxEncompassingLocus)
            {
                var mappedGenes = getGenesMappedTo(chromosome, locus);

                if (mappedGenes.Count() == 0)
                {
                    minEncompassingLocus = maxEncompassingLocus = -1;
                    return false;
                }

                minEncompassingLocus = mappedGenes.Select(x => x.minLocus).Min();
                maxEncompassingLocus = mappedGenes.Select(x => x.maxLocus).Max();

                return true;
            }

            AVLTree<Range> map = new AVLTree<Range>();
        } // GeneMap

        public class IsoformMap
        {
            public IsoformMap(Dictionary<string, GeneLocationInfo> genesByName)
            {
                foreach (var geneInfo in genesByName.Select(x => x.Value).ToList())
                {
                    foreach (var isoform in geneInfo.isoforms)
                    {
                        isoformsByName.Add(isoform.ucscId, isoform);
                        if (!map.ContainsKey(isoform.chromosome))
                        {
                            map.Add(isoform.chromosome, new Dictionary<int, List<Isoform>>());
                        }

                        foreach (var exon in isoform.exons)
                        {
                            for (int locus = exon.start; locus <= exon.end; locus++)
                            {
                                if (!map[isoform.chromosome].ContainsKey(locus))
                                {
                                    map[isoform.chromosome].Add(locus, new List<Isoform>());
                                }

                                map[isoform.chromosome][locus].Add(isoform);
                            } // for each locus in the exon
                        } // for each exon
                    }  // for each isoform
                } // for each gene
            } // ctor

            public List<Isoform> getIsoformsMappedTo(string contig, int locus)
            {
                if (!map.ContainsKey(contig) || !map[contig].ContainsKey(locus))
                {
                    return new List<Isoform>();
                }

                return map[contig][locus];
            }

            public Isoform getIsoform(string ucsdId)
            {
                return isoformsByName[ucsdId];
            }

            public List<Isoform> getAllIsoforms()
            {
                return isoformsByName.Select(_ => _.Value).ToList();
            }

            Dictionary<string, Dictionary<int, List<Isoform>>> map = new Dictionary<string, Dictionary<int, List<Isoform>>>();  // chr->locus->list of isoforms
            Dictionary<string, Isoform> isoformsByName = new Dictionary<string, Isoform>();
        }

        //
        // A Case is a person with TCGA data.
        //
        public class Case
        {
            //
            // Mandatory metadata that's created by GenerateCases.
            //
            public string case_id;
            public string normal_dna_file_id;
            public string tumor_dna_file_id;
            public string normal_rna_file_id = "";    // This is optional
            public string tumor_rna_file_id;
            public string maf_file_id;
            public string tumor_methylation_file_id = "";
            public string normal_methylation_file_id = "";
            public string tumor_copy_number_file_id = "";
            public string normal_copy_number_file_id = "";
            public string tumor_fpkm_file_id = "";
            public string normal_fpkm_file_id = "";
            public string project_id;   // This is TCGA_<DiseaseType> for TCGA.
            public string clinical_supplement_file_id = "";
            public List<string> sample_ids = new List<string>();
            public string normal_miRNA_file_id = "";
            public string tumor_miRNA_file_id = "";
            public string tumor_miRNA_expression_quantification_file_id = "";
            public string normal_miRNA_expression_quantification_file_id = "";
            public string tumor_isoform_expression_quantification_file_id = "";
            public string normal_isoform_expression_quantification_file_id = "";

            //
            // Pathnames for downloaded files.
            //
            public string normal_dna_filename = "";
            public string tumor_dna_filename = "";
            public string normal_rna_filename = "";
            public string tumor_rna_filename = "";
            public string maf_filename = "";
            public string decompressed_maf_filename = "";
            public string tumor_methylation_filename = "";
            public string normal_methylation_filename = "";
            public string tumor_copy_number_filename = "";
            public string normal_copy_number_filename = "";
            public string tumor_fpkm_filename = "";
            public string normal_fpkm_filename = "";
            public string clinical_supplement_filename = "";
            public string normal_miRNA_filename = "";
            public string tumor_miRNA_filename = "";
            public string tumor_miRNA_expression_quantification_filename = "";
            public string normal_miRNA_expression_quantification_filename = "";
            public string tumor_isoform_expression_quantification_filename = "";
            public string normal_isoform_expression_quantification_filename = "";

            //
            // Sizes for downloaded files.
            //
            public long normal_dna_size = 0;
            public long tumor_dna_size = 0;
            public long normal_rna_size = 0;
            public long tumor_rna_size = 0;
            public long tumor_methylation_size = 0;
            public long normal_methylation_size = 0;
            public long tumor_copy_number_size = 0;
            public long normal_copy_number_size = 0;
            public long tumor_fpkm_size = 0;
            public long normal_fpkm_size = 0;
            public long clinical_supplement_size = 0;
            public long normal_miRNA_size = 0;
            public long tumor_miRNA_size = 0;
            public long tumor_miRNA_expression_quantification_size = 0;
            public long normal_miRNA_expression_quantification_size = 0;
            public long tumor_isoform_expression_quantification_size = 0;
            public long normal_isoform_expression_quantification_size = 0;

            //
            // Pathnames for derived files.
            //
            public string normal_dna_allcount_filename = "";
            public string tumor_dna_allcount_filename = "";
            public string normal_rna_allcount_filename = "";
            public string tumor_rna_allcount_filename = "";
            public string regional_expression_filename = "";
            public string gene_expression_filename = "";
            public string tentative_selected_variants_filename = "";
            //public string selected_variants_filename = "";
            public string tumor_dna_reads_at_tentative_selected_variants_filename = "";
            public string tumor_dna_reads_at_tentative_selected_variants_index_filename = "";
            public string tumor_rna_reads_at_tentative_selected_variants_filename = "";
            public string tumor_rna_reads_at_tentative_selected_variants_index_filename = "";
            public string normal_rna_reads_at_tentative_selected_variants_filename = "";
            public string normal_rna_reads_at_tentative_selected_variants_index_filename = "";
            public string normal_dna_reads_at_tentative_selected_variants_filename = "";
            public string normal_dna_reads_at_tentative_selected_variants_index_filename = "";
            public string tentative_annotated_selected_variants_filename = "";
            public string annotated_selected_variants_filename = "";
            public string tumor_allele_specific_gene_expression_filename = "";
            public string normal_allele_specific_gene_expression_filename = "";
            public string tumor_dna_gene_coverage_filname = "";
            public string raw_vcf_filename = "";
            public string vcf_filename = "";
            public string extracted_maf_lines_filename = "";
            public string all_maf_lines_filename = "";
            public string normal_dna_mapped_base_count_filename = "";
            public string tumor_dna_mapped_base_count_filename = "";
            public string normal_rna_mapped_base_count_filename = "";
            public string tumor_rna_mapped_base_count_filename = "";
            public string selected_variant_counts_by_gene_filename = "";
            public string selected_regulatory_maf_filename = "";
            public string annotated_regulatory_regions_filename = "";
            public string annotated_geneHancer_lines_filename = "";
            public string regulatory_mutations_near_mutations_filename = "";
            public string expression_by_gene_filename = "";
            public string isoform_read_counts_filename = "";
            public string compressed_vcf_filename = "";
            public string case_metadata_filename = "";
            public string tentative_asv_without_cnvs_filename = "";
            public string variant_phasing_filename = "";
            public string vcf_statistics_filename = "";
            public string read_statictics_filename = "";
            public string snap_realigned_normal_dna_filename = "";
            public string snap_realigned_normal_dna_bai_filename = "";
            public string snap_realigned_tumor_dna_filename = "";
            public string snap_realigned_tumor_dna_bai_filename = "";
            public string snap_realigned_normal_dna_statictics_filename = "";
            public string snap_realigned_tumor_dna_statictics_filename = "";
            public string bowtie_realigned_normal_dna_filename = "";
            public string bowtie_realigned_normal_dna_bai_filename = "";
            public string bowtie_realigned_tumor_dna_filename = "";
            public string bowtie_realigned_tumor_dna_bai_filename = "";
            public string bowtie_realigned_normal_dna_statictics_filename = "";
            public string bowtie_realigned_tumor_dna_statictics_filename = "";
            public string normal_dna_fastq_filename = "";
            public string normal_dna_fastq_second_end = "";
            public string tumor_dna_fastq = "";
            public string tumor_dna_fastq_second_end = "";

            //
            // If you add another drived file type and it has a **done** terminator, please add it to the CheckDone tool.     
            //

            //
            // Checksums for downloaded files. (NB: I'm not sure how to get the BAI md5s from the metadata, so they're left blank and the
            // BAIs aren't checked.  You have to notice if they're corrupt by hand when samtools can't read them and delete & re-download them.)
            //
            public string normal_rna_file_bam_md5 = "";
            public string normal_rna_file_bai_md5 = "";
            public string tumor_rna_file_bam_md5 = "";
            public string tumor_rna_file_bai_md5 = "";
            public string normal_dna_file_bam_md5 = "";
            public string normal_dna_file_bai_md5 = "";
            public string tumor_dna_file_bam_md5 = "";
            public string tumor_dna_file_bai_md5 = "";
            public string tumor_methylation_file_md5 = "";
            public string normal_methylation_file_md5 = "";
            public string tumor_copy_number_file_md5 = "";
            public string normal_copy_number_file_md5 = "";
            public string tumor_fpkm_file_md5 = "";
            public string normal_fpkm_file_md5 = "";
            public string clinical_supplement_md5 = "";
            public string normal_miRNA_md5 = "";
            public string tumor_miRNA_md5 = "";
            public string tumor_miRNA_expression_quantification_md5 = "";
            public string normal_miRNA_expression_quantification_md5 = "";
            public string tumor_isoform_expression_quantification_md5 = "";
            public string normal_isoform_expression_quantification_md5 = "";

            //
            // Sizes for derived files.
            //
            public long normal_dna_allcount_size = 0;
            public long tumor_dna_allcount_size = 0;
            public long normal_rna_allcount_size = 0;
            public long tumor_rna_allcount_size = 0;
            public long regional_expression_size = 0;
            public long gene_expression_size = 0;
            public long tentative_selected_variants_size = 0;
            public long selected_variants_size = 0;
            public long tumor_dna_reads_at_selected_variants_size = 0;
            public long tumor_dna_reads_at_selected_variants_index_size = 0;
            public long normal_dna_reads_at_selected_variants_size = 0;
            public long normal_dna_reads_at_selected_variants_index_size = 0;
            public long tumor_rna_reads_at_selected_variants_size = 0;
            public long tumor_rna_reads_at_selected_variants_index_size = 0;
            public long normal_rna_reads_at_selected_variants_size = 0;
            public long normal_rna_reads_at_selected_variants_index_size = 0;
            public long tentative_annotated_selected_variants_size = 0;
            public long annotated_selected_variants_size = 0;
            public long normal_allele_specific_gene_expression_size = 0;
            public long tumor_allele_specific_gene_expression_size = 0;
            public long tumor_dna_gene_coverage_size = 0;
            public long raw_vcf_size = 0;
            public long vcf_size = 0;
            public long extracted_maf_lines_size = 0;
            public long all_maf_lines_size = 0;
            public long normal_dna_mapped_base_count_size = 0;
            public long tumor_dna_mapped_base_count_size = 0;
            public long normal_rna_mapped_base_count_size = 0;
            public long tumor_rna_mapped_base_count_size = 0;
            public long selected_variant_counts_by_gene_size = 0;
            public long tumor_regional_methylation_size = 0;
            public long normal_regional_methylation_size = 0;
            public long selected_regulatory_maf_lines_size = 0;
            public long annotated_regulatory_regions_size = 0;
            public long annotated_geneHancer_lines_size = 0;
            public long regulatory_mutations_near_mutations_size = 0;
            public long expression_by_gene_size = 0;
            public long isoform_read_counts_file_size = 0;
            public long compressed_vcf_file_size = 0;
            public long case_metadata_file_size = 0;
            public long tentative_asv_without_cnvs_size = 0;
            public long variant_phasing_size = 0;
            public long vcf_statistics_size = 0;
            public long read_statictics_size = 0;
            public long snap_realigned_normal_dna_size = 0;
            public long snap_realigned_normal_dna_bai_size = 0;
            public long snap_realigned_tumor_dna_size = 0;
            public long snap_realigned_tumor_dna_bai_size = 0;
            public long snap_realigned_normal_dna_statictics_size = 0;
            public long snap_realigned_tumor_dna_statictics_size = 0;
            public long bowtie_realigned_normal_dna_size = 0;
            public long bowtie_realigned_normal_dna_bai_size = 0;
            public long bowtie_realigned_tumor_dna_size = 0;
            public long bowtie_realigned_tumor_dna_bai_size = 0;
            public long bowtie_realigned_normal_dna_statictics_size = 0;
            public long bowtie_realigned_tumor_dna_statictics_size = 0;
            public long normal_dna_fastq_size = 0;
            public long normal_dna_fastq_second_end_size = 0;
            public long tumor_dna_fastq_size = 0;
            public long tumor_dna_fastq_second_end_size = 0;
            //
            // The column numbers from the cases file for these fields.  They're used by C++ programs, which don't have access to the HeaderizedFile class,
            // so so just take the column numbers on the command line.  (That seemed better than compiling file format knowledge into them.)
            //
            static public int ProjectColumn = -1;
            static public int TumorRNAAllcountFilenameColumn = -1;
            static public int TumorRNAMappedBaseCountColumn = -1;

            public string disease()
            {
                if (project_id.Contains('-'))
                {
                    return project_id.Substring(project_id.LastIndexOf('-') + 1).ToLower();
                }

                if (project_id.Contains('_'))
                {
                    return project_id.Substring(project_id.LastIndexOf('_') + 1).ToLower();
                }

                return project_id.ToLower();
            }

            public static Dictionary<string, Case> LoadCases(string inputFilename)
            {
                if (!File.Exists(inputFilename))
                {
                    Console.WriteLine("Cases file " + inputFilename + " doesn't exist.");

                    return null;   // Nothing to load because we haven't generated a cases file yet.
                }

                var wantedFields = new List<string>();

                foreach (var info in AllFields)
                {
                    wantedFields.Add(info.columnName);
                }

                var inputFile = CreateStreamReaderWithRetry(inputFilename);
                if (null == inputFile)
                {
                    return null;
                }

                var headerizedFile = new HeaderizedFile<Case>(inputFile, false, true, "", wantedFields);

                List<Case> listOfCases;
                Dictionary<string, int> fieldMappings;
                if (!headerizedFile.ParseFile(fromSaveFileLine, out listOfCases, out fieldMappings)) {
                    inputFile.Close();
                    return null;
                }

                ProjectColumn = fieldMappings["Project ID"];
                TumorRNAAllcountFilenameColumn = fieldMappings["Tumor RNA Allcount Filename"];
                TumorRNAMappedBaseCountColumn = fieldMappings["Tumor RNA Mapped Base Count Filename"];

                inputFile.Close();

                var cases = new Dictionary<string, Case>();
                foreach (var case_ in listOfCases) {
                    cases.Add(case_.case_id, case_);
                }

                return cases;
            } // LoadCases

            public delegate string ColumnGetter(Case case_);
            public delegate void ColumnSetter(Case case_, string value);
            public delegate string ExpectedIdGetter(Case case_);
            public delegate long SizeColumnGetter(Case case_);
            public delegate void SizeColumnSetter(Case case_, long value);


            public class DownloadableFileType
            {
                public DownloadableFileType(string name_, ColumnGetter fileIdGetter_, ColumnSetter fileIdSetter_, ColumnGetter filenameGetter_, ColumnSetter filenameSetter_,
                                            SizeColumnGetter sizeGetter_, SizeColumnSetter sizeSetter_, ColumnGetter md5Getter_, ColumnSetter md5Setter_)
                {
                    name = name_;
                    fileIdGetter = fileIdGetter_;
                    fileIdSetter = fileIdSetter_;
                    filenameGetter = filenameGetter_;
                    filenameSetter = filenameSetter_;
                    sizeGetter = sizeGetter_;
                    sizeSetter = sizeSetter_;
                    md5Getter = md5Getter_;
                    md5Setter = md5Setter_;
                }

                public DownloadableFileType(string name_, ColumnGetter fileIdGetter_, ColumnSetter fileIdSetter_, ColumnGetter filenameGetter_, ColumnSetter filenameSetter_,
                            SizeColumnGetter sizeGetter_, SizeColumnSetter sizeSetter_, ColumnGetter bamMD5Getter_, ColumnSetter bamMD5Setter_, ColumnGetter baiMD5Getter_, ColumnSetter baiMD5Setter_)
                {
                    name = name_;
                    fileIdGetter = fileIdGetter_;
                    fileIdSetter = fileIdSetter_;
                    filenameGetter = filenameGetter_;
                    filenameSetter = filenameSetter_;
                    sizeGetter = sizeGetter_;
                    sizeSetter = sizeSetter_;
                    bamMD5Getter = bamMD5Getter_;
                    bamMD5Setter = bamMD5Setter_;
                    baiMD5Getter = baiMD5Getter_;
                    baiMD5Setter = baiMD5Setter_;
                }

                public readonly string name;
                public readonly ColumnGetter fileIdGetter;
                public readonly ColumnSetter fileIdSetter;
                public readonly ColumnGetter filenameGetter;
                public readonly ColumnSetter filenameSetter;
                public readonly SizeColumnGetter sizeGetter;
                public readonly SizeColumnSetter sizeSetter;
                public readonly ColumnGetter md5Getter = null;
                public readonly ColumnSetter md5Setter = null;
                public readonly ColumnGetter bamMD5Getter = null;
                public readonly ColumnSetter bamMD5Setter = null;
                public readonly ColumnGetter baiMD5Getter = null;
                public readonly ColumnSetter baiMD5Setter = null;

            }

            public class DerivedFileType
            {
                public DerivedFileType(string name_, ColumnGetter filenameGetter_, ColumnSetter filenameSetter_, DerivedFile.Type derivedFileType_, string extension_, ExpectedIdGetter fileIDGetter_, SizeColumnGetter sizeColumnGetter_, SizeColumnSetter sizeColumnSetter_)
                {
                    name = name_;
                    filenameGetter = filenameGetter_;
                    filenameSetter = filenameSetter_;
                    derivedFileType = derivedFileType_;
                    extension = extension_;
                    fileIDGetter = fileIDGetter_;
                    sizeColumnGetter = sizeColumnGetter_;
                    sizeColumnSetter = sizeColumnSetter_;
                }

                public readonly string name;
                public readonly ColumnGetter filenameGetter;
                public readonly ColumnSetter filenameSetter;
                public readonly DerivedFile.Type derivedFileType;
                public readonly string extension;
                public readonly ExpectedIdGetter fileIDGetter;
                public readonly SizeColumnGetter sizeColumnGetter;
                public readonly SizeColumnSetter sizeColumnSetter;
            } // DerivedFileType

            public static List<DownloadableFileType> downloadableFileTypes;
            public static List<DerivedFileType> derivedFileTypes;

            public static Dictionary<string, DownloadableFileType> downloadableFileTypesByName;
            public static Dictionary<string, DerivedFileType> derivedFileTypesByName;

            class TypeAndFieldInformationInitializer
            {
                public TypeAndFieldInformationInitializer()
                {
                    downloadableFileTypes = new List<DownloadableFileType>();
                    derivedFileTypes = new List<DerivedFileType>();

                    downloadableFileTypes.Add(new DownloadableFileType("Normal DNA",
                                                                        c => c.normal_dna_file_id, (c, v) => c.normal_dna_file_id = v,
                                                                        c => c.normal_dna_filename, (c, v) => c.normal_dna_filename = v,
                                                                        c => c.normal_dna_size, (c, v) => c.normal_dna_size = v,
                                                                        c => c.normal_dna_file_bam_md5, (c, v) => c.normal_dna_file_bam_md5 = v,
                                                                        c => c.normal_dna_file_bai_md5, (c, v) => c.normal_dna_file_bai_md5 = v));

                    downloadableFileTypes.Add(new DownloadableFileType("Tumor DNA",
                                                                        c => c.tumor_dna_file_id, (c, v) => c.tumor_dna_file_id = v,
                                                                        c => c.tumor_dna_filename, (c, v) => c.tumor_dna_filename = v,
                                                                        c => c.tumor_dna_size, (c, v) => c.tumor_dna_size = v,
                                                                        c => c.tumor_dna_file_bam_md5, (c, v) => c.tumor_dna_file_bam_md5 = v,
                                                                        c => c.tumor_dna_file_bai_md5, (c, v) => c.tumor_dna_file_bai_md5 = v));

                    downloadableFileTypes.Add(new DownloadableFileType("Normal RNA",
                                                                        c => c.normal_rna_file_id, (c, v) => c.normal_rna_file_id = v,
                                                                        c => c.normal_rna_filename, (c, v) => c.normal_rna_filename = v,
                                                                        c => c.normal_rna_size, (c, v) => c.normal_rna_size = v,
                                                                        c => c.normal_rna_file_bam_md5, (c, v) => c.normal_rna_file_bam_md5 = v,
                                                                        c => c.normal_rna_file_bai_md5, (c, v) => c.normal_rna_file_bai_md5 = v));

                    downloadableFileTypes.Add(new DownloadableFileType("Tumor RNA",
                                                                        c => c.tumor_rna_file_id, (c, v) => c.tumor_rna_file_id = v,
                                                                        c => c.tumor_rna_filename, (c, v) => c.tumor_rna_filename = v,
                                                                        c => c.tumor_rna_size, (c, v) => c.tumor_rna_size = v,
                                                                        c => c.tumor_rna_file_bam_md5, (c, v) => c.tumor_rna_file_bam_md5 = v,
                                                                        c => c.tumor_rna_file_bai_md5, (c, v) => c.tumor_rna_file_bai_md5 = v));

                    downloadableFileTypes.Add(new DownloadableFileType("MAF",
                                                                        c => c.maf_file_id, (c, v) => c.maf_file_id = v,
                                                                        c => c.maf_filename, (c, v) => c.maf_filename = v,
                                                                        c => 0, (c, v) => { },
                                                                        c => "", (c, v) => { })); // MAF doesn't have size or MD5

                    downloadableFileTypes.Add(new DownloadableFileType("Tumor Methylation",
                                                                        c => c.tumor_methylation_file_id, (c, v) => c.tumor_methylation_file_id = v,
                                                                        c => c.tumor_methylation_filename, (c, v) => c.tumor_methylation_filename = v,
                                                                        c => c.tumor_methylation_size, (c, v) => c.tumor_methylation_size = v,
                                                                        c => c.tumor_methylation_file_md5, (c, v) => c.tumor_methylation_file_md5 = v));

                    downloadableFileTypes.Add(new DownloadableFileType("Normal Methylation",
                                                                        c => c.normal_methylation_file_id, (c, v) => c.normal_methylation_file_id = v,
                                                                        c => c.normal_methylation_filename, (c, v) => c.normal_methylation_filename = v,
                                                                        c => c.normal_methylation_size, (c, v) => c.normal_methylation_size = v,
                                                                        c => c.normal_methylation_file_md5, (c, v) => c.normal_methylation_file_md5 = v));

                    downloadableFileTypes.Add(new DownloadableFileType("Tumor Copy Number",
                                                                        c => c.tumor_copy_number_file_id, (c, v) => c.tumor_copy_number_file_id = v,
                                                                        c => c.tumor_copy_number_filename, (c, v) => c.tumor_copy_number_filename = v,
                                                                        c => c.tumor_copy_number_size, (c, v) => c.tumor_copy_number_size = v,
                                                                        c => c.tumor_copy_number_file_md5, (c, v) => c.tumor_copy_number_file_md5 = v));

                    downloadableFileTypes.Add(new DownloadableFileType("Normal Copy Number",
                                                                        c => c.normal_copy_number_file_id, (c, v) => c.normal_copy_number_file_id = v,
                                                                        c => c.normal_copy_number_filename, (c, v) => c.normal_copy_number_filename = v,
                                                                        c => c.normal_copy_number_size, (c, v) => c.normal_copy_number_size = v,
                                                                        c => c.normal_copy_number_file_md5, (c, v) => c.normal_copy_number_file_md5 = v));

                    downloadableFileTypes.Add(new DownloadableFileType("Tumor FPKM",
                                                                        c => c.tumor_fpkm_file_id, (c, v) => c.tumor_fpkm_file_id = v,
                                                                        c => c.tumor_fpkm_filename, (c, v) => c.tumor_fpkm_filename = v,
                                                                        c => c.tumor_fpkm_size, (c, v) => c.tumor_fpkm_size = v,
                                                                        c => c.tumor_fpkm_file_md5, (c, v) => c.tumor_fpkm_file_md5 = v));

                    downloadableFileTypes.Add(new DownloadableFileType("Normal FPKM",
                                                                        c => c.normal_fpkm_file_id, (c, v) => c.normal_fpkm_file_id = v,
                                                                        c => c.normal_fpkm_filename, (c, v) => c.normal_fpkm_filename = v,
                                                                        c => c.normal_fpkm_size, (c, v) => c.normal_fpkm_size = v,
                                                                        c => c.normal_fpkm_file_md5, (c, v) => c.normal_fpkm_file_md5 = v));

                    downloadableFileTypes.Add(new DownloadableFileType("Clinical Supplement",
                                                                        c => c.clinical_supplement_file_id, (c, v) => c.clinical_supplement_file_id = v,
                                                                        c => c.clinical_supplement_filename, (c, v) => c.clinical_supplement_filename = v,
                                                                        c => c.clinical_supplement_size, (c, v) => c.clinical_supplement_size = v,
                                                                        c => c.clinical_supplement_md5, (c, v) => c.clinical_supplement_md5 = v));

                    downloadableFileTypes.Add(new DownloadableFileType("Tumor miRNA-seq",
                                                                        c => c.tumor_miRNA_file_id, (c, v) => c.tumor_miRNA_file_id = v,
                                                                        c => c.tumor_miRNA_filename, (c, v) => c.tumor_miRNA_filename = v,
                                                                        c => c.tumor_miRNA_size, (c, v) => c.tumor_miRNA_size = v,
                                                                        c => c.tumor_miRNA_md5, (c, v) => c.tumor_miRNA_md5 = v));

                    downloadableFileTypes.Add(new DownloadableFileType("Normal miRNA-seq",
                                                                        c => c.normal_miRNA_file_id, (c, v) => c.normal_miRNA_file_id = v,
                                                                        c => c.normal_miRNA_filename, (c, v) => c.normal_miRNA_filename = v,
                                                                        c => c.normal_miRNA_size, (c, v) => c.normal_miRNA_size = v,
                                                                        c => c.normal_miRNA_md5, (c, v) => c.normal_miRNA_md5 = v));

                    downloadableFileTypes.Add(new DownloadableFileType("Tumor miRNA Expression Quantification",
                                                                        c => c.tumor_miRNA_expression_quantification_file_id, (c, v) => c.tumor_miRNA_expression_quantification_file_id = v,
                                                                        c => c.tumor_miRNA_expression_quantification_filename, (c, v) => c.tumor_miRNA_expression_quantification_filename = v,
                                                                        c => c.tumor_miRNA_expression_quantification_size, (c, v) => c.tumor_miRNA_expression_quantification_size = v,
                                                                        c => c.tumor_miRNA_expression_quantification_md5, (c, v) => c.tumor_miRNA_expression_quantification_md5 = v));

                    downloadableFileTypes.Add(new DownloadableFileType("Normal miRNA Expression Quantification",
                                                                        c => c.normal_miRNA_expression_quantification_file_id, (c, v) => c.normal_miRNA_expression_quantification_file_id = v,
                                                                        c => c.normal_miRNA_expression_quantification_filename, (c, v) => c.normal_miRNA_expression_quantification_filename = v,
                                                                        c => c.normal_miRNA_expression_quantification_size, (c, v) => c.normal_miRNA_expression_quantification_size = v,
                                                                        c => c.normal_miRNA_expression_quantification_md5, (c, v) => c.normal_miRNA_expression_quantification_md5 = v));

                    downloadableFileTypes.Add(new DownloadableFileType("Tumor Isoform Expression Quantification",
                                                                        c => c.tumor_isoform_expression_quantification_file_id, (c, v) => c.tumor_isoform_expression_quantification_file_id = v,
                                                                        c => c.tumor_isoform_expression_quantification_filename, (c, v) => c.tumor_isoform_expression_quantification_filename = v,
                                                                        c => c.tumor_isoform_expression_quantification_size, (c, v) => c.tumor_isoform_expression_quantification_size = v,
                                                                        c => c.tumor_isoform_expression_quantification_md5, (c, v) => c.tumor_isoform_expression_quantification_md5 = v));

                    downloadableFileTypes.Add(new DownloadableFileType("Normal Isoform Expression Quantification",
                                                                        c => c.normal_isoform_expression_quantification_file_id, (c, v) => c.normal_isoform_expression_quantification_file_id = v,
                                                                        c => c.normal_isoform_expression_quantification_filename, (c, v) => c.normal_isoform_expression_quantification_filename = v,
                                                                        c => c.normal_isoform_expression_quantification_size, (c, v) => c.normal_isoform_expression_quantification_size = v,
                                                                        c => c.normal_isoform_expression_quantification_md5, (c, v) => c.normal_isoform_expression_quantification_md5 = v));


                    derivedFileTypes.Add(new DerivedFileType("Normal DNA Allcount", c => c.normal_dna_allcount_filename, (c, v) => c.normal_dna_allcount_filename = v, DerivedFile.Type.NormalDNAAllcount, normalDNAAllcountExtension, c => c.normal_dna_file_id, c => c.normal_dna_allcount_size, (c, v) => c.normal_dna_allcount_size = v));
                    derivedFileTypes.Add(new DerivedFileType("Tumor DNA Allcount", c => c.tumor_dna_allcount_filename, (c, v) => c.tumor_dna_allcount_filename = v, DerivedFile.Type.TumorDNAAllcount, tumorDNAAllcountExtension, c => c.tumor_dna_file_id, c => c.tumor_dna_allcount_size, (c, v) => c.tumor_dna_allcount_size = v));
                    derivedFileTypes.Add(new DerivedFileType("Normal RNA Allcount", c => c.normal_rna_allcount_filename, (c, v) => c.normal_rna_allcount_filename = v, DerivedFile.Type.NormalRNAAllcount, normalRNAAllcountExtension, c => c.normal_rna_file_id, c => c.normal_rna_allcount_size, (c, v) => c.normal_rna_allcount_size = v));
                    derivedFileTypes.Add(new DerivedFileType("Tumor RNA Allcount", c => c.tumor_rna_allcount_filename, (c, v) => c.tumor_rna_allcount_filename = v, DerivedFile.Type.TumorRNAAllcount, tumorRNAAllcountExtension, c => c.tumor_rna_file_id, c => c.tumor_rna_allcount_size, (c, v) => c.tumor_rna_allcount_size = v));
                    derivedFileTypes.Add(new DerivedFileType("Regional Expression", c => c.regional_expression_filename, (c, v) => c.regional_expression_filename = v, DerivedFile.Type.RegionalExpression, regionalExpressionExtension, c => c.tumor_rna_file_id, c => c.regional_expression_size, (c, v) => c.regional_expression_size = v));
                    derivedFileTypes.Add(new DerivedFileType("Gene Expression", c => c.gene_expression_filename, (c, v) => c.gene_expression_filename = v, DerivedFile.Type.GeneExpression, geneExpressionExtension, c => c.case_id, c => c.gene_expression_size, (c, v) => c.gene_expression_size = v));
                    derivedFileTypes.Add(new DerivedFileType("Tentative Selected Variants", c => c.tentative_selected_variants_filename, (c, v) => c.tentative_selected_variants_filename = v, DerivedFile.Type.TentativeSelectedVariants, tentativeSelectedVariantsExtension, c => c.normal_dna_file_id, c => c.tentative_selected_variants_size, (c, v) => c.tentative_selected_variants_size = v));
                    derivedFileTypes.Add(new DerivedFileType("Normal DNA Reads At Selected Variants", c => c.normal_dna_reads_at_tentative_selected_variants_filename, (c, v) => c.normal_dna_reads_at_tentative_selected_variants_filename = v, DerivedFile.Type.NormalDNAReadsAtSelectedVariants, normalDNAReadsAtTentativeSelectedVariantsExtension, c => c.normal_dna_file_id, c => c.normal_dna_reads_at_selected_variants_size, (c, v) => c.normal_dna_reads_at_selected_variants_size = v));
                    derivedFileTypes.Add(new DerivedFileType("Normal DNA Reads At Selected Variants Index", c => c.normal_dna_reads_at_tentative_selected_variants_index_filename, (c, v) => c.normal_dna_reads_at_tentative_selected_variants_index_filename = v, DerivedFile.Type.NormalDNAReadsAtSelectedVariantsIndex, normalDNAReadsAtTentativeSelectedVariantsIndexExtension, c => c.normal_dna_file_id, c => c.normal_dna_reads_at_selected_variants_index_size, (c, v) => c.normal_dna_reads_at_selected_variants_index_size = v));
                    derivedFileTypes.Add(new DerivedFileType("Tumor DNA Reads At Selected Variants", c => c.tumor_dna_reads_at_tentative_selected_variants_filename, (c, v) => c.tumor_dna_reads_at_tentative_selected_variants_filename = v, DerivedFile.Type.TumorDNAReadsAtSelectedVariants, tumorDNAReadsAtTentativeSelectedVariantsExtension, c => c.tumor_dna_file_id, c => c.tumor_dna_reads_at_selected_variants_size, (c, v) => c.tumor_dna_reads_at_selected_variants_size = v));
                    derivedFileTypes.Add(new DerivedFileType("Tumor DNA Reads At Selected Variants Index", c => c.tumor_dna_reads_at_tentative_selected_variants_index_filename, (c, v) => c.tumor_dna_reads_at_tentative_selected_variants_index_filename = v, DerivedFile.Type.TumorDNAReadsAtSelectedVariantsIndex, tumorDNAReadsAtTentativeSelectedVariantsIndexExtension, c => c.tumor_dna_file_id, c => c.tumor_dna_reads_at_selected_variants_index_size, (c, v) => c.tumor_dna_reads_at_selected_variants_index_size = v));
                    derivedFileTypes.Add(new DerivedFileType("Normal RNA Reads At Selected Variants", c => c.normal_rna_reads_at_tentative_selected_variants_filename, (c, v) => c.normal_rna_reads_at_tentative_selected_variants_filename = v, DerivedFile.Type.NormalRNAReadsAtSelectedVariants, normalRNAReadsAtTentativeSelectedVariantsExtension, c => c.normal_rna_file_id, c => c.normal_rna_reads_at_selected_variants_size, (c, v) => c.normal_rna_reads_at_selected_variants_size = v));
                    derivedFileTypes.Add(new DerivedFileType("Normal RNA Reads At Selected Variants Index", c => c.normal_rna_reads_at_tentative_selected_variants_index_filename, (c, v) => c.normal_rna_reads_at_tentative_selected_variants_index_filename = v, DerivedFile.Type.NormalRNAReadsAtSelectedVariantsIndex, normalRNAReadsAtTentativeSelectedVariantsIndexExtension, c => c.normal_rna_file_id, c => c.normal_rna_reads_at_selected_variants_index_size, (c, v) => c.normal_rna_reads_at_selected_variants_index_size = v));
                    derivedFileTypes.Add(new DerivedFileType("Tumor RNA Reads At Selected Variants", c => c.tumor_rna_reads_at_tentative_selected_variants_filename, (c, v) => c.tumor_rna_reads_at_tentative_selected_variants_filename = v, DerivedFile.Type.TumorRNAReadsAtSelectedVariants, tumorRNAReadsAtTentativeSelectedVariantsExtension, c => c.tumor_rna_file_id, c => c.tumor_rna_reads_at_selected_variants_size, (c, v) => c.tumor_rna_reads_at_selected_variants_size = v));
                    derivedFileTypes.Add(new DerivedFileType("Tumor RNA Reads At Selected Variants Index", c => c.tumor_rna_reads_at_tentative_selected_variants_index_filename, (c, v) => c.tumor_rna_reads_at_tentative_selected_variants_index_filename = v, DerivedFile.Type.TumorRNAReadsAtSelectedVariantsIndex, tumorRNAReadsAtTentativeSelectedVariantsIndexExtension, c => c.tumor_rna_file_id, c => c.tumor_rna_reads_at_selected_variants_index_size, (c, v) => c.tumor_rna_reads_at_selected_variants_index_size = v));
                    derivedFileTypes.Add(new DerivedFileType("Tentative Annotated Selected Variants", c => c.tentative_annotated_selected_variants_filename, (c, v) => c.tentative_annotated_selected_variants_filename = v, DerivedFile.Type.TentativeAnnotatedSelectedVariants, tentativeAnnotatedSelectedVariantsExtension, c => c.case_id, c => c.tentative_annotated_selected_variants_size, (c, v) => c.tentative_annotated_selected_variants_size = v));
                    derivedFileTypes.Add(new DerivedFileType("Annotated Selected Variants", c => c.annotated_selected_variants_filename, (c, v) => c.annotated_selected_variants_filename = v, DerivedFile.Type.AnnotatedSelectedVariants, annotatedSelectedVariantsExtension, c => c.case_id, c => c.annotated_selected_variants_size, (c, v) => c.annotated_selected_variants_size = v));
                    derivedFileTypes.Add(new DerivedFileType("Normal Allele Specific Gene Expression", c => c.normal_allele_specific_gene_expression_filename, (c, v) => c.normal_allele_specific_gene_expression_filename = v, DerivedFile.Type.NormalAlleleSpecificGeneExpression, normalAlleleSpecificGeneExpressionExtension, c => c.case_id, c => c.normal_allele_specific_gene_expression_size, (c, v) => c.normal_allele_specific_gene_expression_size = v));
                    derivedFileTypes.Add(new DerivedFileType("Tumor Allele Specific Gene Expression", c => c.tumor_allele_specific_gene_expression_filename, (c, v) => c.tumor_allele_specific_gene_expression_filename = v, DerivedFile.Type.TumorAlleleSpecificGeneExpression, tumorAlleleSpecificGeneExpressionExtension, c => c.case_id, c => c.tumor_allele_specific_gene_expression_size, (c, v) => c.tumor_allele_specific_gene_expression_size = v));
                    derivedFileTypes.Add(new DerivedFileType("Tumor DNA Gene Coverage", c => c.tumor_dna_gene_coverage_filname, (c, v) => c.tumor_dna_gene_coverage_filname = v, DerivedFile.Type.TumorDNAGeneCoverage, tumorDNAGeneCoverageExtension, c => c.tumor_dna_file_id, c => c.tumor_dna_gene_coverage_size, (c, v) => c.tumor_dna_gene_coverage_size = v));
                    derivedFileTypes.Add(new DerivedFileType("VCF", c => c.vcf_filename, (c, v) => c.vcf_filename = v, DerivedFile.Type.VCF, vcfExtension, c => c.normal_dna_file_id, c => c.vcf_size, (c, v) => c.vcf_size = v));
                    derivedFileTypes.Add(new DerivedFileType("Extracted MAF Lines", c => c.extracted_maf_lines_filename, (c, v) => c.extracted_maf_lines_filename = v, DerivedFile.Type.ExtractedMAFLines, extractedMAFLinesExtension, c => c.case_id, c => c.extracted_maf_lines_size, (c, v) => c.extracted_maf_lines_size = v));
                    derivedFileTypes.Add(new DerivedFileType("All MAF Lines", c => c.all_maf_lines_filename, (c, v) => c.all_maf_lines_filename = v, DerivedFile.Type.AllMAFLines, allMAFLinesExtension, c => c.case_id, c => c.all_maf_lines_size, (c, v) => c.all_maf_lines_size = v));
                    derivedFileTypes.Add(new DerivedFileType("Normal DNA Mapped Base Count", c => c.normal_dna_mapped_base_count_filename, (c, v) => c.normal_dna_mapped_base_count_filename = v, DerivedFile.Type.NormalDNAMappedBaseCount, normalDNAMappedBaseCountExtension, c => c.normal_dna_file_id, c => c.normal_dna_mapped_base_count_size, (c, v) => c.normal_dna_mapped_base_count_size = v));
                    derivedFileTypes.Add(new DerivedFileType("Tumor DNA Mapped Base Count", c => c.tumor_dna_mapped_base_count_filename, (c, v) => c.tumor_dna_mapped_base_count_filename = v, DerivedFile.Type.TumorDNAMappedBaseCount, tumorDNAMappedBaseCountExtension, c => c.tumor_dna_file_id, c => c.tumor_dna_mapped_base_count_size, (c, v) => c.tumor_dna_mapped_base_count_size = v));
                    derivedFileTypes.Add(new DerivedFileType("Normal RNA Mapped Base Count", c => c.normal_rna_mapped_base_count_filename, (c, v) => c.normal_rna_mapped_base_count_filename = v, DerivedFile.Type.NormalRNAMappedBaseCount, normalRNAMappedBaseCountExtension, c => c.normal_rna_file_id, c => c.normal_rna_mapped_base_count_size, (c, v) => c.normal_rna_mapped_base_count_size = v));
                    derivedFileTypes.Add(new DerivedFileType("Tumor RNA Mapped Base Count", c => c.tumor_rna_mapped_base_count_filename, (c, v) => c.tumor_rna_mapped_base_count_filename = v, DerivedFile.Type.TumorRNAMappedBaseCount, tumorRNAMappedBaseCountExtension, c => c.tumor_rna_file_id, c => c.tumor_rna_mapped_base_count_size, (c, v) => c.tumor_rna_mapped_base_count_size = v));
                    derivedFileTypes.Add(new DerivedFileType("Selected Variant Counts By Gene", c => c.selected_variant_counts_by_gene_filename, (c, v) => c.selected_variant_counts_by_gene_filename = v, DerivedFile.Type.SelectedVariantCountByGene, selectedVariantCountByGeneExtension, c => c.case_id, c => c.selected_variant_counts_by_gene_size, (c, v) => c.selected_variant_counts_by_gene_size = v));
                    derivedFileTypes.Add(new DerivedFileType("Selected Regulatory MAF Lines", c => c.selected_regulatory_maf_filename, (c, v) => c.selected_regulatory_maf_filename = v, DerivedFile.Type.SelectedRegulatoryMAFLines, selectedRegulatoryMAFLinesExtension, c => c.case_id, c => c.selected_regulatory_maf_lines_size, (c, v) => c.selected_regulatory_maf_lines_size = v));
                    derivedFileTypes.Add(new DerivedFileType("Annotated Regulatory Regions", c => c.annotated_regulatory_regions_filename, (c, v) => c.annotated_regulatory_regions_filename = v, DerivedFile.Type.AnnotatedRegulatoryRegions, annotatedBEDLinesExtension, c => c.case_id, c => c.annotated_regulatory_regions_size, (c, v) => c.annotated_regulatory_regions_size = v));
                    derivedFileTypes.Add(new DerivedFileType("Annotated GeneHancer Lines", c => c.annotated_geneHancer_lines_filename, (c, v) => c.annotated_geneHancer_lines_filename = v, DerivedFile.Type.AnnotatedGeneHancer, annotatedGeneHancerLinesExtension, c => c.case_id, c => c.annotated_geneHancer_lines_size, (c, v) => c.annotated_geneHancer_lines_size = v));
                    derivedFileTypes.Add(new DerivedFileType("Regulatory Mutations Near Mutations", c => c.regulatory_mutations_near_mutations_filename, (c, v) => c.regulatory_mutations_near_mutations_filename = v, DerivedFile.Type.RegulatoryMutationsNearMutations, regulatoryMutationsNearMutationsExtension, c => c.case_id, c => c.regulatory_mutations_near_mutations_size, (c, v) => c.regulatory_mutations_near_mutations_size = v));
                    derivedFileTypes.Add(new DerivedFileType("Expression By Gene", c => c.expression_by_gene_filename, (c, v) => c.expression_by_gene_filename = v, DerivedFile.Type.ExpressionByGene, expressionByGeneExtension, c => c.case_id, c => c.expression_by_gene_size, (c, v) => c.expression_by_gene_size = v));
                    derivedFileTypes.Add(new DerivedFileType("Isoform Read Counts", c => c.isoform_read_counts_filename, (c, v) => c.isoform_read_counts_filename = v, DerivedFile.Type.IsoformReadCounts, isoformReadCountsExtension, c => c.case_id, c => c.isoform_read_counts_file_size, (c, v) => c.isoform_read_counts_file_size = v));
                    derivedFileTypes.Add(new DerivedFileType("Compressed VCF", c => c.compressed_vcf_filename, (c, v) => c.compressed_vcf_filename = v, DerivedFile.Type.CompressedVCF, compressedVCFExtension, c => c.normal_dna_file_id, c => c.compressed_vcf_file_size, (c, v) => c.compressed_vcf_file_size = v));
                    derivedFileTypes.Add(new DerivedFileType("Case Metadata", c => c.case_metadata_filename, (c, v) => c.case_metadata_filename = v, DerivedFile.Type.CaseMetadata, caseMetadataExtension, c => c.case_id, c => c.case_metadata_file_size, (c, v) => c.case_metadata_file_size = v));
                    derivedFileTypes.Add(new DerivedFileType("Tentative ASVs without CNVs", c => c.tentative_asv_without_cnvs_filename, (c, v) => c.tentative_asv_without_cnvs_filename = v, DerivedFile.Type.TentativeASVsWithoutCNVs, tentativeASVsWithoutCNVsExtension, c => c.case_id, c => c.tentative_asv_without_cnvs_size, (c, v) => c.tentative_asv_without_cnvs_size = v));
                    derivedFileTypes.Add(new DerivedFileType("Variant Phasing", c => c.variant_phasing_filename, (c, v) => c.variant_phasing_filename = v, DerivedFile.Type.VariantPhasing, variantPhasingExtension, c => c.case_id, c => c.variant_phasing_size, (c, v) => c.variant_phasing_size = v));
                    derivedFileTypes.Add(new DerivedFileType("VCF Statistics", c => c.vcf_statistics_filename, (c, v) => c.vcf_statistics_filename = v, DerivedFile.Type.VCFStatistics, vcfStatisticsExtension, c => c.normal_dna_file_id, c => c.vcf_statistics_size, (c, v) => c.vcf_statistics_size = v));
                    derivedFileTypes.Add(new DerivedFileType("Read Statistics", c => c.read_statictics_filename, (c, v) => c.read_statictics_filename = v, DerivedFile.Type.ReadStatictics, readStatisticsExtension, c => c.case_id, c => c.read_statictics_size, (c, v) => c.read_statictics_size = v));
                    derivedFileTypes.Add(new DerivedFileType("SNAP Realigned Normal DNA", c => c.snap_realigned_normal_dna_filename, (c, v) => c.snap_realigned_normal_dna_filename = v, DerivedFile.Type.SnapRealignedNormalDNA, snapRealignedNormalDNAExtension, c => c.normal_dna_file_id, c => c.snap_realigned_normal_dna_size, (c, v) => c.snap_realigned_normal_dna_size = v));
                    derivedFileTypes.Add(new DerivedFileType("SNAP Realigned Normal DNA BAI", c => c.snap_realigned_normal_dna_bai_filename, (c, v) => c.snap_realigned_normal_dna_bai_filename = v, DerivedFile.Type.SnapRealignedNormalDNABai, snapRealignedNormalDNABaiExtension, c => c.normal_dna_file_id, c => c.snap_realigned_normal_dna_bai_size, (c, v) => c.snap_realigned_normal_dna_bai_size = v));
                    derivedFileTypes.Add(new DerivedFileType("SNAP Realigned Tumor DNA", c => c.snap_realigned_tumor_dna_filename, (c, v) => c.snap_realigned_tumor_dna_filename = v, DerivedFile.Type.SnapRealignedTumorDNA, snapRealignedTumorDNAExtension, c => c.tumor_dna_file_id, c => c.snap_realigned_tumor_dna_size, (c, v) => c.snap_realigned_tumor_dna_size = v));
                    derivedFileTypes.Add(new DerivedFileType("SNAP Realigned Tumor DNA BAI", c => c.snap_realigned_tumor_dna_bai_filename, (c, v) => c.snap_realigned_tumor_dna_bai_filename = v, DerivedFile.Type.SnapRealignedTumorDNABai, snapRealignedTumorDNABaiExtension, c => c.tumor_dna_file_id, c => c.snap_realigned_tumor_dna_bai_size, (c, v) => c.snap_realigned_tumor_dna_bai_size = v));
                    derivedFileTypes.Add(new DerivedFileType("SNAP Realigned Normal DNA Statictics", c => c.snap_realigned_normal_dna_statictics_filename, (c, v) => c.snap_realigned_normal_dna_statictics_filename = v, DerivedFile.Type.SnapNormalDNAStatictics, snapRealignedNormalDNAStaticticsExtension, c => c.normal_dna_file_id, c => c.snap_realigned_normal_dna_statictics_size, (c, v) => c.snap_realigned_normal_dna_statictics_size = v));
                    derivedFileTypes.Add(new DerivedFileType("SNAP Realigned Tumor DNA Statictics", c => c.snap_realigned_tumor_dna_statictics_filename, (c, v) => c.snap_realigned_tumor_dna_statictics_filename = v, DerivedFile.Type.SnapTumorDNAStatictics, snapRealignedTumorDNAStatisticsExtension, c => c.tumor_dna_file_id, c => c.snap_realigned_tumor_dna_statictics_size, (c, v) => c.snap_realigned_tumor_dna_statictics_size = v));
                    derivedFileTypes.Add(new DerivedFileType("Bowtie Realigned Normal DNA", c => c.bowtie_realigned_normal_dna_filename, (c, v) => c.bowtie_realigned_normal_dna_filename = v, DerivedFile.Type.BowtieRealignedNormalDNA, bowtieRealignedNormalDNAExtension, c => c.normal_dna_file_id, c => c.bowtie_realigned_normal_dna_size, (c, v) => c.bowtie_realigned_normal_dna_size = v));
                    derivedFileTypes.Add(new DerivedFileType("Bowtie Realigned Normal DNA BAI", c => c.bowtie_realigned_normal_dna_bai_filename, (c, v) => c.bowtie_realigned_normal_dna_bai_filename = v, DerivedFile.Type.BowtieRealignedNormalDNABai, bowtieRealignedNormalDNABaiExtension, c => c.normal_dna_file_id, c => c.bowtie_realigned_normal_dna_bai_size, (c, v) => c.bowtie_realigned_normal_dna_bai_size = v));
                    derivedFileTypes.Add(new DerivedFileType("Bowtie Realigned Tumor DNA", c => c.bowtie_realigned_tumor_dna_filename, (c, v) => c.bowtie_realigned_tumor_dna_filename = v, DerivedFile.Type.BowtieRealignedTumorDNA, bowtieRealignedTumorDNAExtension, c => c.tumor_dna_file_id, c => c.bowtie_realigned_tumor_dna_size, (c, v) => c.bowtie_realigned_tumor_dna_size = v));
                    derivedFileTypes.Add(new DerivedFileType("Bowtie Realigned Tumor DNA BAI", c => c.bowtie_realigned_tumor_dna_bai_filename, (c, v) => c.bowtie_realigned_tumor_dna_bai_filename = v, DerivedFile.Type.BowtieRealignedTumorDNABai, bowtieRealignedTumorDNABaiExtension, c => c.tumor_dna_file_id, c => c.bowtie_realigned_tumor_dna_bai_size, (c, v) => c.bowtie_realigned_tumor_dna_bai_size = v));
                    derivedFileTypes.Add(new DerivedFileType("Bowtie Realigned Normal DNA Statictics", c => c.bowtie_realigned_normal_dna_statictics_filename, (c, v) => c.bowtie_realigned_normal_dna_statictics_filename = v, DerivedFile.Type.BowtieNormalDNAStatictics, bowtieRealignedNormalDNAStaticticsExtension, c => c.normal_dna_file_id, c => c.bowtie_realigned_normal_dna_statictics_size, (c, v) => c.bowtie_realigned_normal_dna_statictics_size = v));
                    derivedFileTypes.Add(new DerivedFileType("Bowtie Realigned Tumor DNA Statictics", c => c.bowtie_realigned_tumor_dna_statictics_filename, (c, v) => c.bowtie_realigned_tumor_dna_statictics_filename = v, DerivedFile.Type.BowtieTumorDNAStatictics, bowtieRealignedTumorDNAStatisticsExtension, c => c.tumor_dna_file_id, c => c.bowtie_realigned_tumor_dna_statictics_size, (c, v) => c.bowtie_realigned_tumor_dna_statictics_size = v));
                    derivedFileTypes.Add(new DerivedFileType("Normal DNA FASTQ", c => c.normal_dna_fastq_filename, (c, v) => c.normal_dna_fastq_filename = v, DerivedFile.Type.NormalDNAFastq, normalFastqExtension, c => c.normal_dna_file_id, c => c.normal_dna_fastq_size, (c, v) => c.normal_dna_fastq_size = v));
                    derivedFileTypes.Add(new DerivedFileType("Normal DNA FASTQ Second End", c => c.normal_dna_fastq_filename, (c, v) => c.normal_dna_fastq_second_end = v, DerivedFile.Type.NormalDNAFastqSecondEnd, normalFastqExtension, c => c.normal_dna_file_id, c => c.normal_dna_fastq_size, (c, v) => c.normal_dna_fastq_size = v));

                    downloadableFileTypesByName = downloadableFileTypes.GroupByToDictUnique(_ => _.name);
                    derivedFileTypesByName = derivedFileTypes.GroupByToDictUnique(_ => _.name);

                    //
                    // Now make the FieldInformations.
                    //

                    AllFields = new List<FieldInformation>();
                    AllFields.Add(new FieldInformation("Case ID", c => c.case_id.ToLower(), (c, v) => c.case_id = v));

                    foreach (var downloadableFileType in downloadableFileTypes)
                    {
                        AllFields.Add(new FieldInformation(downloadableFileType.name + " File ID", downloadableFileType.fileIdGetter, downloadableFileType.fileIdSetter));
                    }

                    //
                    // Special metadata that doesn't correspond to a file.
                    //
                    AllFields.Add(new FieldInformation("Project ID", c => c.project_id, (c, v) => c.project_id = v));
                    AllFields.Add(new FieldInformation("Sample IDs", c => c.sampleIdsInCommaSeparatedList(), (c, v) => c.sample_ids = v.Split(',').ToList()));

                    foreach (var downloadableFileType in downloadableFileTypes)
                    {
                        AllFields.Add(new FieldInformation(downloadableFileType.name + " Filename", downloadableFileType.filenameGetter, downloadableFileType.filenameSetter));
                    }

                    foreach (var derivedFileType in derivedFileTypes)
                    {
                        AllFields.Add(new FieldInformation(derivedFileType.name + " Filename", derivedFileType.filenameGetter, derivedFileType.filenameSetter, derivedFileType.derivedFileType, derivedFileType.extension, derivedFileType.fileIDGetter, derivedFileType.name + " File Size", derivedFileType.sizeColumnGetter, derivedFileType.sizeColumnSetter));
                    }


                    foreach (var downloadableFileType in downloadableFileTypes)
                    {
                        AllFields.Add(new FieldInformation(downloadableFileType.name + " Size", c => Convert.ToString(downloadableFileType.sizeGetter(c)), (c, v) => downloadableFileType.sizeSetter(c, LongFromString(v))));
                    } //ctor

                    foreach (var downloadableFileType in downloadableFileTypes)
                    {
                        if (downloadableFileType.md5Getter != null)
                        {
                            AllFields.Add(new FieldInformation(downloadableFileType.name + " MD5", downloadableFileType.md5Getter, downloadableFileType.md5Setter));
                        } else
                        {
                            AllFields.Add(new FieldInformation(downloadableFileType.name + " BAM MD5", downloadableFileType.bamMD5Getter, downloadableFileType.bamMD5Setter));
                            AllFields.Add(new FieldInformation(downloadableFileType.name + " BAI MD5", downloadableFileType.baiMD5Getter, downloadableFileType.baiMD5Setter));
                        }
                    }

                } // ctor
            } // TypeAndFieldInformationInitializer

            static TypeAndFieldInformationInitializer initializer = new TypeAndFieldInformationInitializer();

            public static List<FieldInformation> AllFields;


            public class FieldInformation
            {
                public FieldInformation(string columnName_, ColumnGetter getter_, ColumnSetter setter_, DerivedFile.Type type_, string extension_, ExpectedIdGetter idGetter_, string sizeColumnName_, SizeColumnGetter sizeColumnGetter_, SizeColumnSetter sizeColumnSetter_)
                {
                    columnName = columnName_;
                    getter = getter_;
                    setter = setter_;

                    //
                    // These fields apply only to derived files.  For other fields, use the other constructor.
                    //
                    type = type_;
                    extension = extension_;
                    idGetter = idGetter_;
                    sizeColumnName = sizeColumnName_;
                    sizeColumnGetter = sizeColumnGetter_;
                    sizeColumnSetter = sizeColumnSetter_;
                }

                public FieldInformation(string columnName_, ColumnGetter getter_, ColumnSetter setter_)
                {
                    columnName = columnName_;
                    getter = getter_;
                    setter = setter_;
                }

                public string getValue(Case case_)
                {
                    return getter(case_);
                }

                public void setValue(Case case_, string value)
                {
                    setter(case_, value);
                }

                public void setSize(Case case_, long size)
                {
                    sizeColumnSetter(case_, size);
                }

                public long getSize(Case case_)
                {
                    return sizeColumnGetter(case_);
                }

                public string getExpectedId(Case case_)
                {
                    return idGetter(case_);
                }

                public readonly string columnName;
                ColumnGetter getter;
                ColumnSetter setter;
                ExpectedIdGetter idGetter = c => "";
                public readonly DerivedFile.Type type = DerivedFile.Type.Unknown;
                public readonly string extension = "";
                public readonly string sizeColumnName = ""; // For the size of derived file
                SizeColumnGetter sizeColumnGetter;
                SizeColumnSetter sizeColumnSetter;
            } // FieldInformation 

#if false
            public static FieldInformation[] AllFields =
            {
                new FieldInformation("Case ID",                                             c => c.case_id.ToLower(), (c, v) => c.case_id = v),
                new FieldInformation("Normal DNA File ID",                                  c => c.normal_dna_file_id, (c, v) => c.normal_dna_file_id = v),
                new FieldInformation("Tumor DNA File ID",                                   c => c.tumor_dna_file_id, (c,v) => c.tumor_dna_file_id = v),
                new FieldInformation("Normal RNA File ID",                                  c => c.normal_rna_file_id, (c,v) => c.normal_rna_file_id = v),
                new FieldInformation("Tumor RNA File ID",                                   c => c.tumor_rna_file_id, (c,v) => c.tumor_rna_file_id = v),
                new FieldInformation("MAF File ID",                                         c => c.maf_file_id, (c,v) => c.maf_file_id = v),
                new FieldInformation("Tumor Methylation File ID",                           c => c.tumor_methylation_file_id, (c,v) => c.tumor_methylation_file_id = v),
                new FieldInformation("Normal Methylation File ID",                          c => c.normal_methylation_file_id, (c,v) => c.normal_methylation_file_id = v),
                new FieldInformation("Tumor Copy Number File ID",                           c => c.tumor_copy_number_file_id, (c,v) => c.tumor_copy_number_file_id = v),
                new FieldInformation("Normal Copy Number File ID",                          c => c.normal_copy_number_file_id, (c,v) => c.normal_copy_number_file_id = v),
                new FieldInformation("Tumor FPKM File ID",                                  c => c.tumor_fpkm_file_id, (c,v) => c.tumor_fpkm_file_id = v),
                new FieldInformation("Normal FPKM File ID",                                 c => c.normal_fpkm_file_id, (c,v) => c.normal_fpkm_file_id = v),
                new FieldInformation("Clinical Supplement File ID",                         c => c.clinical_supplement_file_id, (c, v) => c.clinical_supplement_file_id = v),
                new FieldInformation("Project ID",                                          c => c.project_id, (c,v) => c.project_id = v),
                new FieldInformation("Sample IDs",                                          c => c.sampleIdsInCommaSeparatedList(), (c,v) => c.sample_ids = v.Split(',').ToList()),
                new FieldInformation("Normal miRNA-seq File ID",                            c => c.normal_miRNA_file_id, (c, v) => c.normal_miRNA_file_id = v),
                new FieldInformation("Tumor miRNA-seq File ID",                             c => c.tumor_miRNA_file_id, (c, v) => c.tumor_miRNA_file_id = v),
                new FieldInformation("Tumor miRNA Expression Quantification File ID",       c => c.tumor_miRNA_expression_quantification_file_id, (c, v) => c.tumor_miRNA_expression_quantification_file_id = v),
                new FieldInformation("Normal miRNA Expression Quantification File ID",      c => c.normal_miRNA_expression_quantification_file_id, (c, v) => c.normal_miRNA_expression_quantification_file_id = v),
                new FieldInformation("Tumor Isoform Expression Quantification File ID",     c => c.tumor_isoform_expression_quantification_file_id, (c, v) => c.tumor_isoform_expression_quantification_file_id = v),
                new FieldInformation("Normal Isoform Expression Quantification File ID",    c => c.normal_isoform_expression_quantification_file_id, (c, v) => c.normal_isoform_expression_quantification_file_id = v),

                new FieldInformation("Normal DNA Filename",                                 c => c.normal_dna_filename, (c,v) => c.normal_dna_filename = v),
                new FieldInformation("Tumor DNA Filename",                                  c => c.tumor_dna_filename, (c,v) => c.tumor_dna_filename = v),
                new FieldInformation("Normal RNA Filename",                                 c => c.normal_rna_filename, (c,v) => c.normal_rna_filename = v),
                new FieldInformation("Tumor RNA Filename",                                  c => c.tumor_rna_filename, (c,v) => c.tumor_rna_filename = v),
                new FieldInformation("Tumor Methylation Filename",                          c => c.tumor_methylation_filename, (c,v) => c.tumor_methylation_filename = v),
                new FieldInformation("Normal Methylation Filename",                         c => c.normal_methylation_filename, (c,v) => c.normal_methylation_filename = v),
                new FieldInformation("Tumor Copy Number Filename",                          c => c.tumor_copy_number_filename, (c,v) => c.tumor_copy_number_filename = v),
                new FieldInformation("Normal Copy Number Filename",                         c => c.normal_copy_number_filename, (c,v) => c.normal_copy_number_filename = v),
                new FieldInformation("Tumor FPKM Filename",                                 c => c.tumor_fpkm_filename, (c,v) => c.tumor_fpkm_filename = v),
                new FieldInformation("Normal FPKM Filename",                                c => c.normal_fpkm_filename, (c,v) => c.normal_fpkm_filename = v),
                new FieldInformation("MAF Filename",                                        c => c.maf_filename, (c,v) => c.maf_filename = v),
                new FieldInformation("Decompressed MAF Filename",                           c => c.decompressed_maf_filename, (c,v) => c.decompressed_maf_filename = v),
                new FieldInformation("Clinical Supplement Filename",                        c => c.clinical_supplement_filename, (c, v) => c.clinical_supplement_filename = v),
                new FieldInformation("Normal miRNA-seq Filename",                           c => c.normal_miRNA_filename, (c,v) => c.normal_miRNA_filename = v),
                new FieldInformation("Tumor miRNA-seq Filename",                            c => c.tumor_miRNA_filename, (c,v) => c.tumor_miRNA_filename = v),
                new FieldInformation("Tumor miRNA Expression Quantification Filename",      c => c.tumor_miRNA_expression_quantification_filename, (c, v) => c.tumor_miRNA_expression_quantification_filename = v),
                new FieldInformation("Normal miRNA Expression Quantification Filename",     c => c.normal_miRNA_expression_quantification_filename, (c, v) => c.normal_miRNA_expression_quantification_filename = v),
                new FieldInformation("Tumor Isoform Expression Quantification Filename",    c => c.tumor_isoform_expression_quantification_filename, (c, v) => c.tumor_isoform_expression_quantification_filename = v),
                new FieldInformation("Normal Isoform Expression Quantification Filename",   c => c.normal_isoform_expression_quantification_filename, (c, v) => c.normal_isoform_expression_quantification_filename = v),

                new FieldInformation("Normal DNA Size",                                     c => Convert.ToString(c.normal_dna_size), (c,v) => c.normal_dna_size = LongFromString(v)),
                new FieldInformation("Tumor DNA Size",                                      c => Convert.ToString(c.tumor_dna_size), (c,v) => c.tumor_dna_size = LongFromString(v)),
                new FieldInformation("Normal RNA Size",                                     c => Convert.ToString(c.normal_rna_size), (c,v) => c.normal_rna_size = LongFromString(v)),
                new FieldInformation("Tumor RNA Size",                                      c => Convert.ToString(c.tumor_rna_size), (c,v) => c.tumor_rna_size = LongFromString(v)),
                new FieldInformation("Tumor Methylation Size",                              c => Convert.ToString(c.tumor_methylation_size), (c,v) => c.tumor_methylation_size = LongFromString(v)),
                new FieldInformation("Normal Methylation Size",                             c => Convert.ToString(c.normal_methylation_size), (c,v) => c.normal_methylation_size = LongFromString(v)),
                new FieldInformation("Tumor Copy Number Size",                              c => Convert.ToString(c.tumor_copy_number_size), (c,v) => c.tumor_copy_number_size = LongFromString(v)),
                new FieldInformation("Normal Copy Number Size",                             c => Convert.ToString(c.normal_copy_number_size), (c,v) => c.normal_copy_number_size = LongFromString(v)),
                new FieldInformation("Tumor FPKM Size",                                     c => Convert.ToString(c.tumor_fpkm_size), (c,v) => c.tumor_fpkm_size = LongFromString(v)),
                new FieldInformation("Normal FPKM Size",                                    c => Convert.ToString(c.normal_fpkm_size), (c,v) => c.normal_fpkm_size = LongFromString(v)),
                new FieldInformation("Clinical Supplement Size",                            c => Convert.ToString(c.clinical_supplement_size), (c, v) => c.clinical_supplement_size = LongFromString(v)),
                new FieldInformation("Normal miRNA-seq Size",                               c => Convert.ToString(c.normal_miRNA_size), (c, v) => c.normal_miRNA_size = LongFromString(v)),
                new FieldInformation("Tumor miRNA-seq Size",                                c => Convert.ToString(c.tumor_miRNA_size), (c, v) => c.tumor_miRNA_size = LongFromString(v)),
                new FieldInformation("Tumor miRNA Expression Quantification Size",          c => Convert.ToString(c.tumor_miRNA_expression_quantification_size), (c, v) => c.tumor_miRNA_expression_quantification_size = LongFromString(v)),
                new FieldInformation("Normal miRNA Expression Quantification Size",         c => Convert.ToString(c.normal_miRNA_expression_quantification_size), (c, v) => c.normal_miRNA_expression_quantification_size = LongFromString(v)),
                new FieldInformation("Tumor Isoform Expression Quantification Size",        c => Convert.ToString(c.tumor_isoform_expression_quantification_size), (c, v) => c.tumor_isoform_expression_quantification_size = LongFromString(v)),
                new FieldInformation("Normal Isoform Expression Quantification Size",       c => Convert.ToString(c.normal_isoform_expression_quantification_size), (c, v) => c.normal_isoform_expression_quantification_size = LongFromString(v)),

                new FieldInformation("Normal DNA Allcount Filename",                        c => c.normal_dna_allcount_filename, (c,v) => c.normal_dna_allcount_filename = v, DerivedFile.Type.NormalDNAAllcount, normalDNAAllcountExtension, c => c.normal_dna_file_id, "Normal DNA Allcount File Size", c => c.normal_dna_allcount_size, (c,v) => c.normal_dna_allcount_size = v),
                new FieldInformation("Tumor DNA Allcount Filename",                         c => c.tumor_dna_allcount_filename, (c,v) => c.tumor_dna_allcount_filename = v, DerivedFile.Type.TumorDNAAllcount, tumorDNAAllcountExtension, c => c.tumor_dna_file_id, "Tumor DNA Allcount File Size", c => c.tumor_dna_allcount_size, (c, v) => c.tumor_dna_allcount_size = v),
                new FieldInformation("Normal RNA Allcount Filename",                        c => c.normal_rna_allcount_filename, (c,v) => c.normal_rna_allcount_filename = v, DerivedFile.Type.NormalRNAAllcount, normalRNAAllcountExtension, c => c.normal_rna_file_id, "Normal RNA Allcount File Size", c => c.normal_rna_allcount_size, (c, v) => c.normal_rna_allcount_size = v),
                new FieldInformation("Tumor RNA Allcount Filename",                         c => c.tumor_rna_allcount_filename, (c,v) => c.tumor_rna_allcount_filename = v, DerivedFile.Type.TumorRNAAllcount, tumorRNAAllcountExtension, c => c.tumor_rna_file_id, "Tumor RNA Allcount File Size", c => c.tumor_rna_allcount_size, (c, v) => c.tumor_rna_allcount_size = v),
                new FieldInformation("Regional Expression Filename",                        c => c.regional_expression_filename, (c,v) => c.regional_expression_filename = v, DerivedFile.Type.RegionalExpression, regionalExpressionExtension, c => c.tumor_rna_file_id, "Regional Expression File Size", c => c.regional_expression_size, (c, v) => c.regional_expression_size = v),
                new FieldInformation("Gene Expression Filename",                            c => c.gene_expression_filename, (c,v) => c.gene_expression_filename = v, DerivedFile.Type.GeneExpression, geneExpressionExtension, c => c.case_id, "Gene Expression File Size", c => c.gene_expression_size, (c, v) => c.gene_expression_size = v),
                new FieldInformation("Tentative Selected Variants",                         c => c.tentative_selected_variants_filename, (c, v) => c.tentative_selected_variants_filename = v, DerivedFile.Type.TentativeSelectedVariants, tentativeSelectedVariantsExtension, c => c.normal_dna_file_id, "Tentative Selected Variants Size", c => c.tentative_selected_variants_size, (c, v) => c.tentative_selected_variants_size = v),
                //new FieldInformation("Selected Variants Filename",                          c => c.selected_variants_filename, (c,v) => c.selected_variants_filename = v, DerivedFile.Type.SelectedVariants, selectedVariantsExtension, c => c.normal_dna_file_id, "Selected Variants File Size", c=> c.selected_variants_size, (c, v) => c.selected_variants_size = v),
                new FieldInformation("Normal DNA Reads At Selected Variants Filename",      c => c.normal_dna_reads_at_tentative_selected_variants_filename, (c,v) => c.normal_dna_reads_at_tentative_selected_variants_filename = v, DerivedFile.Type.NormalDNAReadsAtSelectedVariants, normalDNAReadsAtTentativeSelectedVariantsExtension, c => c.normal_dna_file_id, "Normal DNA Reads At Selected Variants File Size", c => c.normal_dna_reads_at_selected_variants_size, (c, v) => c.normal_dna_reads_at_selected_variants_size = v),
                new FieldInformation("Normal DNA Reads At Selected Variants Index Filename",c => c.normal_dna_reads_at_tentative_selected_variants_index_filename, (c,v) => c.normal_dna_reads_at_tentative_selected_variants_index_filename = v, DerivedFile.Type.NormalDNAReadsAtSelectedVariantsIndex, normalDNAReadsAtTentativeSelectedVariantsIndexExtension, c => c.normal_dna_file_id, "Normal DNA Reads At Selected Variants Index File Size", c => c.normal_dna_reads_at_selected_variants_index_size, (c, v) => c.normal_dna_reads_at_selected_variants_index_size = v),
                new FieldInformation("Tumor DNA Reads At Selected Variants Filename",       c => c.tumor_dna_reads_at_tentative_selected_variants_filename, (c,v) => c.tumor_dna_reads_at_tentative_selected_variants_filename = v, DerivedFile.Type.TumorDNAReadsAtSelectedVariants, tumorDNAReadsAtTentativeSelectedVariantsExtension, c => c.tumor_dna_file_id, "Tumor DNA Reads At Selected Variants File Size", c => c.tumor_dna_reads_at_selected_variants_size, (c, v) => c.tumor_dna_reads_at_selected_variants_size = v),
                new FieldInformation("Tumor DNA Reads At Selected Variants Index Filename", c => c.tumor_dna_reads_at_tentative_selected_variants_index_filename, (c,v) => c.tumor_dna_reads_at_tentative_selected_variants_index_filename = v, DerivedFile.Type.TumorDNAReadsAtSelectedVariantsIndex, tumorDNAReadsAtTentativeSelectedVariantsIndexExtension, c => c.tumor_dna_file_id, "Tumor DNA Reads At Selected Variants Index File Size", c => c.tumor_dna_reads_at_selected_variants_index_size, (c, v) => c.tumor_dna_reads_at_selected_variants_index_size = v),
                new FieldInformation("Normal RNA Reads At Selected Variants Filename",      c => c.normal_rna_reads_at_tentative_selected_variants_filename, (c,v) => c.normal_rna_reads_at_tentative_selected_variants_filename = v, DerivedFile.Type.NormalRNAReadsAtSelectedVariants, normalRNAReadsAtTentativeSelectedVariantsExtension, c => c.normal_rna_file_id, "Normal RNA Reads At Selected Variants File Size", c => c.normal_rna_reads_at_selected_variants_size, (c, v) => c.normal_rna_reads_at_selected_variants_size = v),
                new FieldInformation("Normal RNA Reads At Selected Variants Index Filename",c => c.normal_rna_reads_at_tentative_selected_variants_index_filename, (c,v) => c.normal_rna_reads_at_tentative_selected_variants_index_filename = v, DerivedFile.Type.NormalRNAReadsAtSelectedVariantsIndex, normalRNAReadsAtTentativeSelectedVariantsIndexExtension, c => c.normal_rna_file_id, "Normal RNA Reads At Selected Variants Index File Size", c => c.normal_rna_reads_at_selected_variants_index_size, (c, v) => c.normal_rna_reads_at_selected_variants_index_size = v),
                new FieldInformation("Tumor RNA Reads At Selected Variants Filename",       c => c.tumor_rna_reads_at_tentative_selected_variants_filename, (c,v) => c.tumor_rna_reads_at_tentative_selected_variants_filename = v, DerivedFile.Type.TumorRNAReadsAtSelectedVariants, tumorRNAReadsAtTentativeSelectedVariantsExtension, c => c.tumor_rna_file_id, "Tumor RNA Reads At Selected Variants File Size", c => c.tumor_rna_reads_at_selected_variants_size, (c, v) => c.tumor_rna_reads_at_selected_variants_size = v),
                new FieldInformation("Tumor RNA Reads At Selected Variants Index Filename", c => c.tumor_rna_reads_at_tentative_selected_variants_index_filename, (c,v) => c.tumor_rna_reads_at_tentative_selected_variants_index_filename = v, DerivedFile.Type.TumorRNAReadsAtSelectedVariantsIndex, tumorRNAReadsAtTentativeSelectedVariantsIndexExtension, c => c.tumor_rna_file_id, "Tumor RNA Reads At Selected Variants Index File Size", c => c.tumor_rna_reads_at_selected_variants_index_size, (c, v) => c.tumor_rna_reads_at_selected_variants_index_size = v),
                new FieldInformation("Tentative Annotated Selected Variants Filename",      c => c.tentative_annotated_selected_variants_filename, (c, v) => c.tentative_annotated_selected_variants_filename = v, DerivedFile.Type.TentativeAnnotatedSelectedVariants, tentativeAnnotatedSelectedVariantsExtension, c => c.case_id, "Tentative Annotated Selected Variants", c => c.tentative_annotated_selected_variants_size, (c, v) => c.tentative_annotated_selected_variants_size = v),
                new FieldInformation("Annotated Selected Variants Filename",                c => c.annotated_selected_variants_filename, (c,v) => c.annotated_selected_variants_filename = v, DerivedFile.Type.AnnotatedSelectedVariants, annotatedSelectedVariantsExtension, c => c.case_id, "Annotated Selected Variants File Size", c => c.annotated_selected_variants_size, (c, v) => c.annotated_selected_variants_size = v),
                new FieldInformation("Normal Allele Specific Gene Expression Filename",     c => c.normal_allele_specific_gene_expression_filename, (c,v) => c.normal_allele_specific_gene_expression_filename = v, DerivedFile.Type.NormalAlleleSpecificGeneExpression, normalAlleleSpecificGeneExpressionExtension, c => c.case_id, "Normal Allele Specific Gene Expression File Size", c => c.normal_allele_specific_gene_expression_size, (c, v) => c.normal_allele_specific_gene_expression_size = v),
                new FieldInformation("Tumor Allele Specific Gene Expression Filename",      c => c.tumor_allele_specific_gene_expression_filename, (c,v) => c.tumor_allele_specific_gene_expression_filename = v, DerivedFile.Type.TumorAlleleSpecificGeneExpression, tumorAlleleSpecificGeneExpressionExtension, c => c.case_id, "Tumor Allele Specific Gene Expression File Size", c => c.tumor_allele_specific_gene_expression_size, (c, v) => c.tumor_allele_specific_gene_expression_size = v),
                new FieldInformation("Tumor DNA Gene Coverage Filename",                    c => c.tumor_dna_gene_coverage_filname, (c,v) => c.tumor_dna_gene_coverage_filname = v, DerivedFile.Type.TumorDNAGeneCoverage, tumorDNAGeneCoverageExtension, c => c.tumor_dna_file_id, "Tumor DNA Gene Coverage File Size", c => c.tumor_dna_gene_coverage_size, (c, v) => c.tumor_dna_gene_coverage_size = v),
                new FieldInformation("VCF Filename",                                        c => c.vcf_filename, (c,v) => c.vcf_filename = v, DerivedFile.Type.VCF, vcfExtension, c => c.normal_dna_file_id, "VCF File Size", c => c.vcf_size, (c, v) => c.vcf_size = v),
                new FieldInformation("Extracted MAF Lines Filename",                        c => c.extracted_maf_lines_filename, (c,v) => c.extracted_maf_lines_filename = v, DerivedFile.Type.ExtractedMAFLines, extractedMAFLinesExtension, c => c.case_id, "Extracted MAF Lines File Size", c => c.extracted_maf_lines_size, (c, v) => c.extracted_maf_lines_size = v),
                new FieldInformation("All MAF Lines Filename",                              c => c.all_maf_lines_filename, (c,v) => c.all_maf_lines_filename = v, DerivedFile.Type.AllMAFLines, allMAFLinesExtension, c => c.case_id, "All MAF Lines File Size", c => c.all_maf_lines_size, (c, v) => c.all_maf_lines_size = v),
                new FieldInformation("Normal DNA Mapped Base Count Filename",               c => c.normal_dna_mapped_base_count_filename, (c, v) => c.normal_dna_mapped_base_count_filename = v, DerivedFile.Type.NormalDNAMappedBaseCount, normalDNAMappedBaseCountExtension, c => c.normal_dna_file_id, "Normal DNA Mapped Base Count File Size", c => c.normal_dna_mapped_base_count_size, (c, v) => c.normal_dna_mapped_base_count_size = v),
                new FieldInformation("Tumor DNA Mapped Base Count Filename",                c => c.tumor_dna_mapped_base_count_filename, (c, v) => c.tumor_dna_mapped_base_count_filename = v, DerivedFile.Type.TumorDNAMappedBaseCount, tumorDNAMappedBaseCountExtension, c => c.tumor_dna_file_id, "Tumor DNA Mapped Base Count File Size", c => c.tumor_dna_mapped_base_count_size, (c, v) => c.tumor_dna_mapped_base_count_size = v),
                new FieldInformation("Normal RNA Mapped Base Count Filename",               c => c.normal_rna_mapped_base_count_filename, (c, v) => c.normal_rna_mapped_base_count_filename = v, DerivedFile.Type.NormalRNAMappedBaseCount, normalRNAMappedBaseCountExtension, c => c.normal_rna_file_id, "Normal RNA Mapped Base Count File Size", c => c.normal_rna_mapped_base_count_size, (c, v) => c.normal_rna_mapped_base_count_size = v),
                new FieldInformation("Tumor RNA Mapped Base Count Filename",                c => c.tumor_rna_mapped_base_count_filename, (c, v) => c.tumor_rna_mapped_base_count_filename = v, DerivedFile.Type.TumorRNAMappedBaseCount, tumorRNAMappedBaseCountExtension, c => c.tumor_rna_file_id, "Tumor RNA Mapped Base Count File Size", c => c.tumor_rna_mapped_base_count_size, (c, v) => c.tumor_rna_mapped_base_count_size = v),
                new FieldInformation("Selected Variant Counts By Gene Filename",            c => c.selected_variant_counts_by_gene_filename, (c, v) => c.selected_variant_counts_by_gene_filename = v, DerivedFile.Type.SelectedVariantCountByGene, selectedVariantCountByGeneExtension, c => c.case_id, "Selected Variant Counts By Gene File Size", c => c.selected_variant_counts_by_gene_size, (c, v) => c.selected_variant_counts_by_gene_size = v),
                new FieldInformation("Selected Regulatory MAF Lines Filename",              c => c.selected_regulatory_maf_filename, (c, v) => c.selected_regulatory_maf_filename = v, DerivedFile.Type.SelectedRegulatoryMAFLines, selectedRegulatoryMAFLinesExtension, c => c.case_id, "Selected Regulatory MAF Lines Size", c => c.selected_regulatory_maf_lines_size, (c, v) => c.selected_regulatory_maf_lines_size = v),
                new FieldInformation("Annotated Regulatory Regions Filename",               c => c.annotated_regulatory_regions_filename, (c, v) => c.annotated_regulatory_regions_filename = v, DerivedFile.Type.AnnotatedRegulatoryRegions, annotatedBEDLinesExtension, c => c.case_id, "Annotated Regulatory Regions Size", c => c.annotated_regulatory_regions_size, (c, v) => c.annotated_regulatory_regions_size = v),
                new FieldInformation("Annotated GeneHancer Lines Filename",                 c => c.annotated_geneHancer_lines_filename, (c, v) => c.annotated_geneHancer_lines_filename = v, DerivedFile.Type.AnnotatedGeneHancer, annotatedGeneHancerLinesExtension, c => c.case_id, "Annotated GeneHancer Lines Size", c => c.annotated_geneHancer_lines_size, (c, v) => c.annotated_geneHancer_lines_size = v),            
                new FieldInformation("Regulatory Mutations Near Mutations Filename",        c => c.regulatory_mutations_near_mutations_filename, (c, v) => c.regulatory_mutations_near_mutations_filename = v, DerivedFile.Type.RegulatoryMutationsNearMutations, regulatoryMutationsNearMutationsExtension, c => c.case_id, "Regulatory Mutations Near Mutations Size", c => c.regulatory_mutations_near_mutations_size, (c, v) => c.regulatory_mutations_near_mutations_size = v),
                new FieldInformation("Expression By Gene",                                  c => c.expression_by_gene_filename, (c, v) => c.expression_by_gene_filename = v, DerivedFile.Type.ExpressionByGene, expressionByGeneExtension, c => c.case_id, "Expression By Gene Size", c => c.expression_by_gene_size, (c, v) => c.expression_by_gene_size = v),
                new FieldInformation("Isoform Read Counts",                                 c => c.isoform_read_counts_filename, (c, v) => c.isoform_read_counts_filename = v, DerivedFile.Type.IsoformReadCounts, isoformReadCountsExtension, c => c.case_id, "Isoform Gene Counts Size", c => c.isoform_read_counts_file_size, (c, v) => c.isoform_read_counts_file_size = v),
                new FieldInformation("Compressed VCF",                                      c => c.compressed_vcf_filename, (c, v) => c.compressed_vcf_filename = v, DerivedFile.Type.CompressedVCF, compressedVCFExtension, c => c.normal_dna_file_id, "Compressed VCF Size", c => c.compressed_vcf_file_size, (c, v) => c.compressed_vcf_file_size = v),
                new FieldInformation("Case Metadata",                                       c => c.case_metadata_filename, (c, v) => c.case_metadata_filename = v, DerivedFile.Type.CaseMetadata, caseMetadataExtension, c => c.case_id, "Case Metadata Size", c => c.case_metadata_file_size, (c, v) => c.case_metadata_file_size = v),
                new FieldInformation("Tentative ASVs without CNVs Filename",                c => c.tentative_asv_without_cnvs_filename, (c, v) => c.tentative_asv_without_cnvs_filename = v, DerivedFile.Type.TentativeASVsWithoutCNVs, tentativeASVsWithoutCNVsExtension, c => c.case_id, "Tentative ASVs Without CNVs Size", c => c.tentative_asv_without_cnvs_size, (c, v) => c.tentative_asv_without_cnvs_size = v),
                new FieldInformation("Variant Phasing Filename",                            c => c.variant_phasing_filename, (c, v) => c.variant_phasing_filename = v, DerivedFile.Type.VariantPhasing, variantPhasingExtension, c => c.case_id, "Variant Phasing Size", c => c.variant_phasing_size, (c, v) => c.variant_phasing_size = v),
                new FieldInformation("VCF Statistics Filename",                             c => c.vcf_statistics_filename, (c,v) => c.vcf_statistics_filename = v, DerivedFile.Type.VCFStatistics, vcfStatisticsExtension, c => c.normal_dna_file_id, "VCF Statistics Size", c=> c.vcf_statistics_size, (c, v) => c.vcf_statistics_size = v),
                new FieldInformation("Read Statistics Filename",                            c => c.read_statictics_filename, (c, v) => c.read_statictics_filename = v, DerivedFile.Type.ReadStatictics, readStatisticsExtension, c => c.case_id, "Read Statictics Size", c => c.read_statictics_size, (c, v) => c.read_statictics_size = v)                        ,   

                new FieldInformation("Normal RNA BAM MD5",                                  c => c.normal_rna_file_bam_md5, (c,v) => c.normal_rna_file_bam_md5 = v),
                new FieldInformation("Normal RNA BAI MD5",                                  c => c.normal_rna_file_bai_md5, (c,v) => c.normal_rna_file_bai_md5 = v),
                new FieldInformation("Tumor RNA BAM MD5",                                   c => c.tumor_rna_file_bam_md5, (c,v) => c.tumor_rna_file_bam_md5 = v),
                new FieldInformation("Tumor RNA BAI MD5",                                   c => c.tumor_rna_file_bai_md5, (c,v) => c.tumor_rna_file_bai_md5 = v),
                new FieldInformation("Normal DNA BAM MD5",                                  c => c.normal_dna_file_bam_md5, (c,v) => c.normal_dna_file_bam_md5 = v),
                new FieldInformation("Normal DNA BAI MD5",                                  c => c.normal_dna_file_bai_md5, (c,v) => c.normal_dna_file_bai_md5 = v),
                new FieldInformation("Tumor DNA BAM MD5",                                   c => c.tumor_dna_file_bam_md5, (c,v) => c.tumor_dna_file_bam_md5 = v),
                new FieldInformation("Tumor DNA BAI MD5",                                   c => c.tumor_dna_file_bai_md5, (c,v) => c.tumor_dna_file_bai_md5 = v),
                new FieldInformation("Tumor Methylation MD5",                               c => c.tumor_methylation_file_md5, (c,v) => c.tumor_methylation_file_md5 = v),
                new FieldInformation("Normal Methylation MD5",                              c => c.normal_methylation_file_md5, (c,v) => c.normal_methylation_file_md5 = v),
                new FieldInformation("Tumor Copy Number MD5",                               c => c.tumor_copy_number_file_md5, (c,v) => c.tumor_copy_number_file_md5 = v),
                new FieldInformation("Normal Copy Number MD5",                              c => c.normal_copy_number_file_md5, (c,v) => c.normal_copy_number_file_md5 = v),
                new FieldInformation("Tumor FPKM MD5",                                      c => c.tumor_fpkm_file_md5, (c,v) => c.tumor_fpkm_file_md5 = v),
                new FieldInformation("Normal FPKM MD5",                                     c => c.normal_fpkm_file_md5, (c,v) => c.normal_fpkm_file_md5 = v),
                new FieldInformation("Clinical Supplement MD5",                             c => c.clinical_supplement_md5, (c,v) => c.clinical_supplement_md5 = v),
                new FieldInformation("Normal miRNA MD5",                                    c => c.normal_miRNA_md5, (c,v) => c.normal_miRNA_md5 = v),
                new FieldInformation("Tumor miRNA MD5",                                     c => c.tumor_miRNA_md5, (c,v) => c.tumor_miRNA_md5 = v),
                new FieldInformation("Tumor miRNA Expression Quantification MD5",           c => c.tumor_miRNA_expression_quantification_md5, (c, v) => c.tumor_miRNA_expression_quantification_md5 = v),
                new FieldInformation("Normal miRNA Expression Quantification MD5",          c => c.normal_miRNA_expression_quantification_md5, (c, v) => c.normal_miRNA_expression_quantification_md5 = v),
                new FieldInformation("Tumor Isoform Expression Quantification MD5",         c => c.tumor_isoform_expression_quantification_md5, (c, v) => c.tumor_isoform_expression_quantification_md5 = v),
                new FieldInformation("Normal Isoform Expression Quantification MD5",        c => c.normal_isoform_expression_quantification_md5, (c, v) => c.normal_isoform_expression_quantification_md5 = v),

            }; // fieldInformation
#endif // false


            static public Case fromSaveFileLine(Dictionary<string, int> fieldMappings, string[] fields)
            {
                var case_ = new Case();

                foreach (var info in AllFields)
                {
                    info.setValue(case_, fields[fieldMappings[info.columnName]]);
                }

                return case_;
            } // fromSaveFile

            string sampleIdsInCommaSeparatedList()
            {
                if (sample_ids.Count() == 0) {
                    return "";
                }

                var retVal = sample_ids[0];

                for (int i = 1; i < sample_ids.Count(); i++) {
                    retVal += "," + sample_ids[i];
                }

                return retVal;

            }   // sampleIdsInCommaSeparatedList    

            public string GenerateLine()
            {

                string value = AllFields[0].getValue(this);

                for (int i = 1; i < AllFields.Count(); i++)
                {
                    value += "\t" + AllFields[i].getValue(this);
                }

                for (int i = 0; i < AllFields.Count(); i++)
                {
                    if (AllFields[i].type != DerivedFile.Type.Unknown)
                    {
                        value += "\t" + AllFields[i].getSize(this);
                    }
                }

                return value;

            } // GenerateLine

            public static void SaveToFile(Dictionary<string, Case> cases, string filename)
            {
                var file = CreateStreamWriterWithRetry(filename);

                //
                // Header
                //
                for (int i = 0; i < AllFields.Count() - 1; i++) // We go to the next-to-last to get the tabs right
                {
                    file.Write(AllFields[i].columnName + "\t");
                }

                file.Write(AllFields[AllFields.Count() - 1].columnName);

                for (int i = 0; i < AllFields.Count(); i++)
                {
                    if (AllFields[i].type != DerivedFile.Type.Unknown)
                    {
                        file.Write("\t" + AllFields[i].sizeColumnName);
                    }
                }

                file.WriteLine();

                //
                // Cases
                //
                foreach (var case_ in cases)
                {
                    file.WriteLine(case_.Value.GenerateLine());
                }
                file.WriteLine("**done**");
                file.Close();
            } // SaveToFile

            public static void loadAllFileLocations(Configuration configuration, Dictionary<string, Case> cases, Dictionary<string, DownloadedFile> downloadedFiles, Dictionary<string, List<DerivedFile>> derivedFiles)
            {
                foreach (var caseEntry in cases)
                {
                    caseEntry.Value.loadFileLocations(configuration, downloadedFiles, derivedFiles);
                }
            }

            delegate void Setter(string value);
            static void addDownloadedFileLocation(Dictionary<string, DownloadedFile> downloadedFiles, string fileID, Setter setter)
            {
                if (fileID != "" && downloadedFiles.ContainsKey(fileID))
                {
                    setter(downloadedFiles[fileID].fileInfo.FullName);
                } else
                {
                    setter("");
                }
            }
            //
            // Fill in all the derived file locations if they exist.
            //
            public void loadFileLocations(Configuration configuration, Dictionary<string, DownloadedFile> downloadedFiles, Dictionary<string, List<DerivedFile>> derivedFiles)
            {

                addDownloadedFileLocation(downloadedFiles, normal_dna_file_id, x => normal_dna_filename = x);

                addDownloadedFileLocation(downloadedFiles, tumor_dna_file_id, x => tumor_dna_filename = x);
                if (normal_rna_file_id != "" && downloadedFiles.ContainsKey(normal_rna_file_id))
                {
                    normal_rna_filename = downloadedFiles[normal_rna_file_id].fileInfo.FullName;
                }
                else
                {
                    normal_rna_filename = "";
                }

                addDownloadedFileLocation(downloadedFiles, tumor_rna_file_id, x => tumor_rna_filename = x);
                if (!configuration.isBeatAML)
                {
                    if (downloadedFiles.ContainsKey(maf_file_id))
                    {
                        maf_filename = downloadedFiles[maf_file_id].fileInfo.FullName;
                    }
                    else
                    {
                        maf_filename = "";
                    }

                    if (downloadedFiles.ContainsKey(disease() + ".maf"))
                    {
                        decompressed_maf_filename = downloadedFiles[disease() + ".maf"].fileInfo.FullName;
                    }
                }

                addDownloadedFileLocation(downloadedFiles, tumor_methylation_file_id, x => tumor_methylation_filename = x);
                addDownloadedFileLocation(downloadedFiles, normal_methylation_file_id, x => normal_methylation_filename = x);
                addDownloadedFileLocation(downloadedFiles, tumor_copy_number_file_id, x => tumor_copy_number_filename = x);
                addDownloadedFileLocation(downloadedFiles, normal_copy_number_file_id, x => normal_copy_number_filename = x);
                addDownloadedFileLocation(downloadedFiles, tumor_fpkm_file_id, x => tumor_fpkm_filename = x);
                addDownloadedFileLocation(downloadedFiles, normal_fpkm_file_id, x => normal_fpkm_filename = x);
                addDownloadedFileLocation(downloadedFiles, tumor_isoform_expression_quantification_file_id, x => tumor_isoform_expression_quantification_filename = x);
                addDownloadedFileLocation(downloadedFiles, normal_isoform_expression_quantification_file_id, x => normal_isoform_expression_quantification_filename = x);
                addDownloadedFileLocation(downloadedFiles, tumor_miRNA_expression_quantification_file_id, x => tumor_miRNA_expression_quantification_filename = x);
                addDownloadedFileLocation(downloadedFiles, tumor_miRNA_file_id, x => tumor_miRNA_filename = x);
                addDownloadedFileLocation(downloadedFiles, normal_miRNA_file_id, x => normal_miRNA_filename = x);


                if (!derivedFiles.ContainsKey(case_id))
                {
                    tumor_rna_allcount_filename = "";
                    tumor_allele_specific_gene_expression_filename = "";
                    normal_allele_specific_gene_expression_filename = "";
                    annotated_selected_variants_filename = "";
                    tumor_dna_reads_at_tentative_selected_variants_filename = "";
                    tumor_dna_reads_at_tentative_selected_variants_index_filename = "";
                    gene_expression_filename = "";
                    tumor_dna_reads_at_tentative_selected_variants_filename = "";
                    regional_expression_filename = "";
                    tumor_rna_reads_at_tentative_selected_variants_filename = "";
                    tumor_rna_reads_at_tentative_selected_variants_index_filename = "";
                    //selected_variants_filename = "";
                    tentative_selected_variants_filename = "";
                    tentative_annotated_selected_variants_filename = "";
                    vcf_filename = "";
                    return;
                }

                Dictionary<DerivedFile.Type, List<DerivedFile>> derivedFilesByType = new Dictionary<DerivedFile.Type, List<DerivedFile>>();
                var derivedFilesForThisCase = derivedFiles[case_id];
                foreach (var type in (DerivedFile.Type[])Enum.GetValues(typeof(DerivedFile.Type))) {
                    derivedFilesByType.Add(type, derivedFilesForThisCase.Where(x => x.type == type).ToList());
                    if (type != DerivedFile.Type.Unknown && derivedFilesByType[type].Count() > 1) {
                        Console.Write("Found more than one file of type " + type + " for case " + case_id + ":");
                        foreach (var file in derivedFilesByType[type]) {
                            Console.Write(" " + file.fileinfo.FullName);
                        }
                        Console.WriteLine();
                    }
                }

                foreach (var field in AllFields)
                {
                    if (field.type == DerivedFile.Type.Unknown)
                    {
                        continue;
                    }

                    if (derivedFilesByType[field.type].Count() > 0)
                    {
                        var filename = derivedFilesByType[field.type][0].fileinfo.FullName;
                        field.setValue(this, filename);
                        field.setSize(this, derivedFilesByType[field.type][0].fileinfo.Length);
                        if (GetFileIdFromPathname(filename) != field.getExpectedId(this))
                        {
                            Console.WriteLine("Found derived file with unexpected file ID: " + filename + ", expected file id: " + field.getExpectedId(this));
                            derivedFilesByType[field.type] = new List<DerivedFile>();   // Just blow it away and don't use it.
                        }
                    }
                    else
                    {
                        field.setValue(this, "");
                    }
                }
            } // loadFileLocations

            public string getAllcountFilename(bool tumor, bool dna)
            {
                if (tumor)
                {
                    if (dna) return tumor_dna_allcount_filename;
                    return tumor_rna_allcount_filename;
                }
                if (dna) return normal_dna_allcount_filename;
                return normal_rna_allcount_filename;
            } // getAllcountFilename

            public string getReadsAtTentativeSelectedVariantsFilename(bool tumor, bool dna)
            {
                if (tumor)
                {
                    if (dna) return tumor_dna_reads_at_tentative_selected_variants_filename;
                    else return tumor_rna_reads_at_tentative_selected_variants_filename;
                } else
                {
                    if (dna) return normal_dna_reads_at_tentative_selected_variants_filename;
                    return normal_rna_reads_at_tentative_selected_variants_filename;
                }
            } // getReadsAtSelectedVariantsFilename


            public string getReadsAtTentativeSelectedVariantsIndexFilename(bool tumor, bool dna)
            {
                if (tumor)
                {
                    if (dna) return tumor_dna_reads_at_tentative_selected_variants_index_filename;
                    else return tumor_rna_reads_at_tentative_selected_variants_index_filename;
                }
                else
                {
                    if (dna) return normal_dna_reads_at_tentative_selected_variants_index_filename;
                    else return normal_rna_reads_at_tentative_selected_variants_index_filename;
                }
            } // getReadsAtSelectedVariantsFilename

            public string getDownloadedReadsFileId(bool tumor, bool dna)
            {
                if (tumor)
                {
                    if (dna) return tumor_dna_file_id;
                    return tumor_rna_file_id;
                } else
                {
                    if (dna) return normal_dna_file_id;
                    return normal_rna_file_id;
                }
            } // getDownloadedReadsFileId

            public string getDownloadedReadsFilename(bool tumor, bool dna)
            {
                if (tumor)
                {
                    if (dna) return tumor_dna_filename;
                    return tumor_rna_filename;
                } else
                {
                    if (dna) return normal_dna_filename;
                    return normal_rna_filename;
                }
            } // getDownloadedReadsFilename

            public string getDownloadedReadsBamMD5(bool tumor, bool dna)
            {
                if (tumor)
                {
                    if (dna) return tumor_dna_file_bam_md5;
                    return tumor_rna_file_bam_md5;
                } else
                {
                    if (dna) return normal_dna_file_bam_md5;
                    return normal_rna_file_bam_md5;
                }
            }

        } // Case

        public static long LongFromString(string String)
        {
            if (String == null || String == "") return 0;

            try
            {
                return Convert.ToInt64(String);
            }
            catch (FormatException)
            {
                Console.WriteLine("LongFromString: string '" + String + "' didn't parse.");
                return 0;
            }
        }

        static bool isExceptionFatalForRetryOpen(Exception e)
        {
            return e is FileNotFoundException || e is DirectoryNotFoundException || e.Message.Contains("The network path was not found");
        }

        public static StreamWriter CreateStreamWriterWithRetry(string filename, int bufferSize = 0)
        {
            while (true)
            {
                try
                {
                    if (bufferSize == 0)
                    {
                        return new StreamWriter(filename);
                    } else
                    {
                        return new StreamWriter(filename, false, new UTF8Encoding(), bufferSize);  // false means don't append
                    }
                }
                catch (IOException e)
                {
                    if (isExceptionFatalForRetryOpen(e))
                    {
                        return null;
                    }
                    Console.WriteLine("IOException opening " + filename + " for write.  Sleeping and retrying.");
                    Thread.Sleep(10 * 1000);
                }
            }
        }

        public static StreamWriter CreateAppendingStreamWriterWithRetry(string filename)
        {
            while (true)
            {
                try
                {
                    var writer = File.AppendText(filename);
                    return writer;
                }
                catch (IOException e)
                {
                    if (isExceptionFatalForRetryOpen(e))
                    {
                        return null;
                    }
                    Console.WriteLine("IOException opening " + filename + " for append.  Sleeping and retrying.");
                    Thread.Sleep(10 * 1000);
                }
            }
        }

        public static StreamReader CreateStreamReaderWithRetry(string filename)
        {
            while (true)
            {
                try
                {
                    var reader = new StreamReader(filename);
                    return reader;
                }
                catch (IOException e)
                {
                    if (isExceptionFatalForRetryOpen(e))
                    {
                        return null;
                    }
                    Console.WriteLine("IOException opening " + filename + " for read.  Sleeping and retrying.");
                    Thread.Sleep(10 * 1000);
                }
            }
        }

        public static StreamReader CreateCompressedStreamReaderWithRetry(string filename)
        {
            var innerReader = CreateStreamReaderWithRetry(filename);
            if (innerReader == null) {
                return null;
            }
            return new StreamReader(new GZipStream(innerReader.BaseStream, CompressionMode.Decompress));
        }

        public static string[] ReadAllLinesWithRetry(string filename)
        {
            while (true)
            {
                try
                {
                    var lines = File.ReadAllLines(filename);
                    return lines;
                }
                catch (IOException e)
                {
                    if (isExceptionFatalForRetryOpen(e))
                    {
                        return null;
                    }

                    Console.WriteLine("IOException reading " + filename + ".  Sleeping and retrying.");
                    Thread.Sleep(10 * 1000);
                }
            }
        }

        public static StreamReader CreateStreamReaderWithRetryCompressedBasedOnFilename(string filename)
        {
            if (filename.ToLower().EndsWith(".gz"))
            {
                return CreateCompressedStreamReaderWithRetry(filename);
            } else
            {
                return CreateStreamReaderWithRetry(filename);
            }
        }

        public class Configuration
        {
            // If you add anything that references defaultBaseDirectory, be sure to add it to
            // reinitializeWithBaseDirectory()

            public const string defaultBaseDirectory = @"\\air-k26-01\d$\gdc\";
            public const string defaultConfigurationFilePathame = defaultBaseDirectory + "configuration.txt";

            public string baseDirectory = defaultBaseDirectory;

            public string accessTokenPathname = defaultBaseDirectory + @"access_token.txt";

            public const string defaultGenomeBuild = "hg38";
            public string ensemblToGeneFilename = defaultBaseDirectory + @"ensemblGeneNames.txt";

            public string hg19GeneLocationInformation = defaultBaseDirectory + "knownGene-hg19.txt";

            public List<string> dataDirectories = new List<string>();
            public string mafManifestPathname = defaultBaseDirectory + "mafManifest.txt";
            public string mutationCaller = "mutect";
            public List<string> programNames = new List<string>();
            public string binariesDirectory = defaultBaseDirectory + @"bin\";
            public string configuationFilePathname = defaultConfigurationFilePathame;
            public string casesFilePathname = defaultBaseDirectory + "cases.txt";

            public string indexDirectory = defaultBaseDirectory + @"indices\hg38-16";
            public string localIndexDirectory
                = @"d:\sequence\indices\hg38-20-large";
            public string derivedFilesDirectory = "derived_files";    // This is relative to each download directory
            public string hpcScriptFilename = "";    // The empty string says to black hole this script
            public string hpcScheduler = "gcr";
            public string hpcBinariesDirectory = @"\\gcr\scratch\b99\bolosky\";
            public string hpcIndexDirectory = defaultBaseDirectory + @"indices\hg38-16";
            public string azureScriptFilename = ""; // The empty string says to black hole this script
            public string expressionFilesDirectory = defaultBaseDirectory + @"expression\";
            public string basesInKnownCodingRegionsDirectory = defaultBaseDirectory + @"bases_in_known_coding_regions\";
            public string completedVCFsDirectory = "";  // Where completed VCFs from Azure are dropped
            public List<string> excludedFileIDs = new List<string>();
            public int regionalExpressionRegionSize = 1000;
            public int nTumorsToIncludeGene = 30;   // How many tumors must have at least one mutation in a gene in order to include it.
            public int nReadsToIncludeVariant = 30; // How much coverage do we need of a mutation or somatic variant to consider it?
            public string selectedGenesFilename = defaultBaseDirectory + @"seleted_genes.txt";
            public string scriptOutputDirectory = @"\temp\";
            public string finalResultsDirectory = defaultBaseDirectory + @"final_results\"; // includes the trailing backslash
            public double significanceLevel = 0.01;
            public string geneLocationInformationFilename = defaultBaseDirectory + "knownGene-" + defaultGenomeBuild + ".txt";
            public int minSampleCountForHeatMap = 100;
            public double ASEInNormalAtWhichToExcludeGenes = 0.5;  // If there's at least this much ASE overall in the matched normal for a gene, then don't use selected variants for it.
            public bool downloadedFilesHaveMD5Sums = true;  // The BeatAML data set doesn't have them (though the synapse client does internally).
            public string samplesSummaryPathname = @"\\fds-k24-09\d$\BeatAML-20180411\Samples Summary.txt";
            public string synapseDirectory = @"\\fds-k24-09\d$\BeatAML-20180411\";
            public string patient_metadata_directory = defaultBaseDirectory + @"patient_metadata\";

            public string geneScatterGraphsDirectory = defaultBaseDirectory + @"gene_scatter_graphs\";
            public string geneScatterGraphsLinesWithPercentilesDirectory = defaultBaseDirectory + @"gene_scatter_graphs_with_percentiles\";
            public string encodeBEDFile = defaultBaseDirectory + @"encode\ENCFF621SFC_hg38.bed";
            // items used in bisulfite analysis
            public string bisulfiteDirectory = defaultBaseDirectory + @"bisulfate\";
            public string bisulfiteCasesFilePathname = defaultBaseDirectory + @"bisulfite\cases_bisulfite.txt";

            public int readLengthForRepetitiveRegionDetection = 48;

            public string chromosomeMapsDirectory = defaultBaseDirectory + @"chromosome_maps\";
            public string redundantChromosomeRegionFilename = defaultBaseDirectory + "redundantRegions.txt";

            // chain files
            public string hg38Tohg19ChainFile = defaultBaseDirectory + @"chain\hg38ToHg19.over.chain";
            public string hg19Tohg38ChainFile = defaultBaseDirectory + @"chain\hg19ToHg38.over.chain";

            public string regionalExpressionDirectory = defaultBaseDirectory + @"regional_expression\";

            public string unfilteredCountsDirectory = defaultBaseDirectory + @"gene_mutations_with_counts\";
            public const string unfilteredCountsExtention = @"_unfiltered_counts.txt";
            public string methylationREFsFilename = defaultBaseDirectory + "compositeREFs450.txt";

            public const string geneCategorizationFilename = "GeneCategorization.txt";

            public string zero_one_two_directory = defaultBaseDirectory + @"012graphs\";
            //public string expression_distribution_directory = defaultBaseDirectory + @"expression_distribution\";
            // string expression_distribution_by_chromosome_directory = defaultBaseDirectory + @"expression_distribution_by_chromosome\";
            public string expression_distribution_by_chromosome_map_filename = defaultBaseDirectory + @"expression_distribution_by_chromosomes_map.txt";

            public string geneHancerFilename = defaultBaseDirectory + @"genehancer.txt";

            public const string GlobalScatterGraphFilename = "ASEByCase.txt";

            public const string PerGeneExpressionHistogramsFilename = "PerGeneExpressionHistograms.txt";

            public bool usesChr = true; // Do the data files use "chr" in the chromosome names?

            public bool configurationFileLocationExplicitlySpecified = false;

            public int minRNAReadCoverage = 10;
            public int minDNAReadCoverage = 10;

            public int maxProximityForReferenceBiasCheck = 25;
            public int minDepthForReferenceBiasCheck = 30;      // This is an interpretation of the phrase "a few dozen" in the wikipedia article on confidence intervals.
            public int minDistanceFromRepetitiveRegion = 50;    // How far from a region with multiply mapped reads must we be to include a germline variant (if it's not decidable by actual reference bias measurements)
            public double repetitiveRegionConfidence = 0.95;   // 
            public int minDistanceBetweenGermlineVariants = 1000;

            public const int minInstancesOfEachClassToComputeMannWhitney = 10;  // This should be a parameter rather than a compiled-in constant.

            public int nWorkerMachines = 120;                 // Over how many machines would we like to evenly split our work?

            public bool isBeatAML = false;

            void reinitialzieWithBaseDirectory()
            {
                accessTokenPathname = baseDirectory + @"access_token.txt";
                ensemblToGeneFilename = baseDirectory + @"ensemblGeneNames.txt";
                binariesDirectory = baseDirectory + @"bin\";
                hg19GeneLocationInformation = baseDirectory + "knownGene-hg19.txt";
                mafManifestPathname = baseDirectory + "mafManifest.txt";
                binariesDirectory = baseDirectory + @"bin\";
                casesFilePathname = baseDirectory + "cases.txt";
                indexDirectory = baseDirectory + @"indices\hg38-16";
                hpcIndexDirectory = baseDirectory + @"indices\hg38-16";
                expressionFilesDirectory = baseDirectory + @"expression\";
                basesInKnownCodingRegionsDirectory = baseDirectory + @"bases_in_known_coding_regions\";
                selectedGenesFilename = baseDirectory + @"seleted_genes.txt";
                finalResultsDirectory = baseDirectory + @"final_results\";
                geneLocationInformationFilename = baseDirectory + "knownGene-" + defaultGenomeBuild + ".txt";
                geneScatterGraphsDirectory = baseDirectory + @"gene_scatter_graphs\";
                geneScatterGraphsLinesWithPercentilesDirectory = baseDirectory + @"gene_scatter_graphs_with_percentiles\";
                encodeBEDFile = baseDirectory + @"encode\ENCFF621SFC_hg38.bed";
                bisulfiteDirectory = baseDirectory + @"bisulfate\";
                bisulfiteCasesFilePathname = baseDirectory + @"bisulfite\cases_bisulfite.txt";
                chromosomeMapsDirectory = baseDirectory + @"chromosome_maps\";
                redundantChromosomeRegionFilename = baseDirectory + "redundantRegions.txt";
                hg38Tohg19ChainFile = baseDirectory + @"chain\hg38ToHg19.over.chain";
                hg19Tohg38ChainFile = baseDirectory + @"chain\hg19ToHg38.over.chain";
                regionalExpressionDirectory = baseDirectory + @"regional_expression\";
                unfilteredCountsDirectory = baseDirectory + @"gene_mutations_with_counts\";
                methylationREFsFilename = baseDirectory + "compositeREFs450.txt";
                zero_one_two_directory = baseDirectory + @"012graphs\";
                geneHancerFilename = baseDirectory + @"genehancer.txt";
                patient_metadata_directory = baseDirectory + @"patient_metadata\";


                dataDirectories = new List<string>();
                dataDirectories.Add(baseDirectory + @"downloaded_files\");
            }

            public List<string> neededPerReplicaDirectories()
            {
                var neededDirectories = new List<string>();

                foreach (var dataDirectory in dataDirectories)
                {
                    neededDirectories.Add(dataDirectory);
                    neededDirectories.Add(dataDirectory + @"..\" + derivedFilesDirectory);
                }

                return neededDirectories;
            }

            public List<string> neededGlobalDirectories()
            {
                var neededDirectories = new List<string>();

                neededDirectories.Add(zero_one_two_directory);
                neededDirectories.Add(basesInKnownCodingRegionsDirectory);
                neededDirectories.Add(chromosomeMapsDirectory);
                neededDirectories.Add(completedVCFsDirectory);
                neededDirectories.Add(GetDirectoryFromPathname(encodeBEDFile));
                neededDirectories.Add(expressionFilesDirectory);
                neededDirectories.Add(finalResultsDirectory);
                neededDirectories.Add(unfilteredCountsDirectory);
                neededDirectories.Add(geneScatterGraphsDirectory);
                neededDirectories.Add(geneScatterGraphsLinesWithPercentilesDirectory);
                neededDirectories.Add(indexDirectory + @"..\");  // One level up, because indexDirectory includes hg38-20, which we don't want here.
                neededDirectories.Add(regionalExpressionDirectory);

                return neededDirectories;
            }

            public int getMinReadCoverage(bool dna)
            {
                if (dna) return minDNAReadCoverage;
                return minRNAReadCoverage;
            }


            public string[] commandLineArgs = null;    // The args excluding -configuration <filename>

            Configuration()
            {
                programNames.Add("TCGA");   // The default value
                dataDirectories.Add(defaultBaseDirectory + @"downloaded_files\");
            }

            //
            // Parse the args to find the configuration file pathname, and then load from that path (or the default if it's not present).
            //
            public static Configuration loadFromFile(string[] args)
            {
                string pathname = @"\\air-k26-01\d$\gdc\configuration.txt";

                bool fromCommandLine = false;
                var nonConsumedArgs = new List<string>();
                for (int i = 0; i < args.Count(); i++) {
                    if (args[i] == "-configuration") {
                        if (i >= args.Count() - 1) {
                            Console.WriteLine("-configuation can't be the last parameter, it needs to be followed by a file path.  Ignoring.");
                        } else {
                            pathname = args[i + 1];
                            i++;
                            fromCommandLine = true;
                        }
                    }
                    else
                    {
                        nonConsumedArgs.Add(args[i]);
                    }
                }

                var retVal = new Configuration();
                retVal.commandLineArgs = nonConsumedArgs.ToArray();
                retVal.configurationFileLocationExplicitlySpecified = fromCommandLine;

                if (!File.Exists(pathname))
                {
                    //
                    // No config file means to use the default.
                    //
                    if (fromCommandLine)
                    {
                        Console.WriteLine("Unable to open configuration file " + pathname);
                        return null;
                    }
                    return retVal;
                }

                retVal.configuationFilePathname = pathname; // Don't load this from the file, just keep track of where we got it from.

                var lines = ReadAllLinesWithRetry(pathname);
                bool seenProgramNames = false;
                bool seenDataDirectories = false;
                bool firstLine = true; // This is used because we can only reset the base directory in the first line, since resetting it changes lots of other fields.

                foreach (var line in lines)
                {
                    if (line.Count() == 0 || line[0] == '#')
                    {
                        continue;   // Ignore blank lines or comments.
                    }

                    var fields = line.Split('\t');
                    if (fields.Count() != 2) {
                        Console.WriteLine("ASEConfiguration.loadFromFile: configuration file " + pathname + " contains a line that doesn't have exactly two tab separated fields: " + line + ".  Ignoring.");
                        continue;
                    }

                    var type = fields[0].ToLower();
                    if (type == "access token") {
                        retVal.accessTokenPathname = fields[1];
                    } else if (type == "data directory") {
                        if (!fields[1].EndsWith(@"\"))
                        {
                            Console.WriteLine("The configuration contains a data directory that doesn't end with a backslash: " + fields[1] + ".  Ignoring.");
                            continue;
                        }
                        if (!seenDataDirectories)   // Don't keep the default
                        {
                            retVal.dataDirectories = new List<string>();
                            seenDataDirectories = true;
                        }
                        retVal.dataDirectories.Add(fields[1]);
                    } else if (type == "maf manifest") {
                        retVal.mafManifestPathname = fields[1];
                    } else if (type == "mutation caller") {
                        retVal.mutationCaller = fields[1];
                    } else if (type == "program name") {
                        if (!seenProgramNames)  // If we've got the default value, override it.
                        {
                            retVal.programNames = new List<string>();
                            seenProgramNames = true;
                        }
                        retVal.programNames.Add(fields[1]);
                    } else if (type == "binary directory") {
                        retVal.binariesDirectory = fields[1];
                    } else if (type == "cases") {
                        retVal.casesFilePathname = fields[1];
                    } else if (type == "index directory") {
                        retVal.indexDirectory = fields[1];
                    } else if (type == "derived files") {
                        retVal.derivedFilesDirectory = fields[1];
                    } else if (type == "hpc script name") {
                        retVal.hpcScriptFilename = fields[1];
                    } else if (type == "hpc scheduler") {
                        retVal.hpcScheduler = fields[1];
                    } else if (type == "hpc binaries directory") {
                        retVal.hpcBinariesDirectory = fields[1];
                    } else if (type == "hpc index directory") {
                        retVal.hpcIndexDirectory = fields[1];
                    } else if (type == "expression files directory") {
                        retVal.expressionFilesDirectory = fields[1];
                    } else if (type == "azure script name") {
                        retVal.azureScriptFilename = fields[1];
                    } else if (type == "completed vcfs directory") {
                        retVal.completedVCFsDirectory = fields[1];
                    } else if (type == "regional expression region size") {
                        retVal.regionalExpressionRegionSize = Convert.ToInt32(fields[1]);   // should probably wrap this in a try
                    } else if (type == "excluded file id") {
                        retVal.excludedFileIDs.Add(fields[1]);
                    } else if (type == "selected genes filename") {
                        retVal.selectedGenesFilename = fields[1];
                    } else if (type == "count of tumors to include gene") {
                        retVal.nTumorsToIncludeGene = Convert.ToInt32(fields[1]);
                    } else if (type == "script output directory") {
                        retVal.scriptOutputDirectory = fields[1];
                    } else if (type == "final results directory") {
                        retVal.finalResultsDirectory = fields[1];
                    } else if (type == "gene scatter graphs directory") {
                        retVal.geneScatterGraphsDirectory = fields[1];
                    } else if (type == "gene scatter graphs with percentiles directory") {
                        retVal.geneScatterGraphsLinesWithPercentilesDirectory = fields[1];
                    } else if (type == "read count to include variant") {
                        retVal.nReadsToIncludeVariant = Convert.ToInt32(fields[1]);
                    } else if (type == "regional expression directory") {
                        retVal.regionalExpressionDirectory = fields[1];
                    } else if (type == "significance level") {
                        retVal.significanceLevel = Convert.ToDouble(fields[1]);
                    } else if (type == "gene location information filename") {
                        retVal.geneLocationInformationFilename = fields[1];
                    } else if (type == "min sample count for heatmap") {
                        retVal.minSampleCountForHeatMap = Convert.ToInt32(fields[1]);
                    } else if (type == "max normal ase for using variants") {
                        retVal.ASEInNormalAtWhichToExcludeGenes = Convert.ToDouble(fields[1]);
                    } else if (type == "min rna read coverage") {
                        retVal.minRNAReadCoverage = Convert.ToInt32(fields[1]);
                    } else if (type == "min dna read coverage") {
                        retVal.minDNAReadCoverage = Convert.ToInt32(fields[1]);
                    } else if (type == "encode bed file") {
                        retVal.encodeBEDFile = fields[1];
                    } else if (type == "read length for repetitive region detection") {
                        retVal.readLengthForRepetitiveRegionDetection = Convert.ToInt32(fields[1]);
                    } else if (type == "downloaded files have md5 sums") {
                        if (fields[1].ToLower() == "false")
                        {
                            retVal.downloadedFilesHaveMD5Sums = false;
                        } else if (fields[1].ToLower() == "true")
                        {
                            retVal.downloadedFilesHaveMD5Sums = true;
                        } else
                        {
                            Console.WriteLine("'downloaded files have md5 sums' configuration parameter must have value true or false.");
                            return null;
                        }
                    } else if (type == "beat aml") {
                        if (fields[1].ToLower() == "true")
                        {
                            retVal.isBeatAML = true;
                        } else
                        {
                            retVal.isBeatAML = false;
                        }
                    } else if (type == "base directory") {
                        if (!firstLine)
                        {
                            Console.WriteLine("You can only change the base directory as the first line of the configuration file.");
                            return null;
                        }
                        retVal.baseDirectory = fields[1];
                        retVal.reinitialzieWithBaseDirectory();
                    } else if (type == "synapse directory") {
                        retVal.synapseDirectory = fields[1];
                    } else if (type == "genehancer filename") {
                        retVal.geneHancerFilename = fields[1];
                    } else if (type == "local index directory") {
                        retVal.localIndexDirectory = fields[1];
                    } else if (type == "n worker machines") {
                        retVal.nWorkerMachines = Convert.ToInt32(fields[1]);
                    } else if (type == "uses chr") {
                        if (fields[1] == "yes")
                        {
                            retVal.usesChr = true;
                        } else if (fields[1] == "no")
                        {
                            retVal.usesChr = false;
                        } else
                        {
                            Console.WriteLine("value in configuration file after 'uses chr' must be 'yes' or 'no'.");
                            continue;
                        }
                    } else if (type == "max distance for reference bias proximity check") {
                        retVal.maxProximityForReferenceBiasCheck = Convert.ToInt32(fields[1]);
                    } else if (type == "min depth for reference bias check") {
                        retVal.minDepthForReferenceBiasCheck = Convert.ToInt32(fields[1]);
                    } else {
                        Console.WriteLine("ASEConfiguration.loadFromFile: configuration file " + pathname + " contains a line with an unknown configuration parameter type: " + line + ".  Ignoring.");
                        continue;
                    }

                    firstLine = false;
                } // foreach (var line in lines) 

                return retVal;
            }
        } // Configuration

        class ExpressionDistributionByChromosomeMapLine
        {
            public readonly string disease;
            public readonly string chromosome;
            public readonly string filename;

            public ExpressionDistributionByChromosomeMapLine(string disease_, string chromosome_, string filename_)
            {
                disease = disease_;
                chromosome = chromosome_;
                filename = filename_;
            }
        }
        public class ExpressionDistributionByChromosomeMap
        {
            public static ExpressionDistributionByChromosomeMap LoadFromFile(string filename)
            {
                var retVal = new ExpressionDistributionByChromosomeMap();

                if (!File.Exists(filename))
                {
                    return retVal;  // No file -> empty map
                }

                var inputFile = CreateStreamReaderWithRetry(filename);
                if (null == inputFile)
                {
                    return null;
                }

                var headerizedFile = new HeaderizedFile<ExpressionDistributionByChromosomeMapLine>(inputFile, false, true, "", wantedFields.ToList());

                List<ExpressionDistributionByChromosomeMapLine> lines;
                headerizedFile.ParseFile(ParseOneLine, out lines);

                foreach (var line in lines)
                {
                    if (!retVal.map.ContainsKey(line.disease))
                    {
                        retVal.map.Add(line.disease, new Dictionary<string, string>());
                    }

                    retVal.map[line.disease].Add(line.chromosome, line.filename);
                }

                return retVal;
            }

            static ExpressionDistributionByChromosomeMap LoadFromFile(Configuration configuration)
            {
                return LoadFromFile(configuration.expression_distribution_by_chromosome_map_filename);
            }

            static ExpressionDistributionByChromosomeMapLine ParseOneLine(HeaderizedFile<ExpressionDistributionByChromosomeMapLine>.FieldGrabber fieldGrabber)
            {
                return new ExpressionDistributionByChromosomeMapLine(fieldGrabber.AsString("disease"), fieldGrabber.AsString("chromosome"), fieldGrabber.AsString("filename"));
            }

            public void WriteToFile(string outputFilename)
            {
                var outputFile = CreateStreamWriterWithRetry(outputFilename);

                if (null == outputFile)
                {
                    throw new Exception("ExpressionDistributionByChromosomeMap.WriteToFile: unable to open output file " + outputFilename);
                }

                outputFile.WriteLine("disease\tchromosome\tfilename");

                foreach (var diseaseEntry in map)
                {
                    var disease = diseaseEntry.Key;
                    foreach (var chromosomeEntry in diseaseEntry.Value)
                    {
                        var chromosome = chromosomeEntry.Key;
                        var filename = chromosomeEntry.Value;
                        outputFile.WriteLine(disease + "\t" + chromosome + "\t" + filename);
                    } // chromosome
                } // disease

                outputFile.WriteLine("**done**");
                outputFile.Close();
            } // WriteToFile

            public List<string> allFiles()
            {
                var retVal = new List<string>();

                foreach (var diseaseEntry in map)
                {
                    foreach (var chromosomeEntry in diseaseEntry.Value)
                    {
                        retVal.Add(chromosomeEntry.Value);
                    }
                }

                return retVal;
            } // allFiles()

            public Dictionary<string, Dictionary<string, string>> map = new Dictionary<string, Dictionary<string, string>>(); // disease->chromosome->filename

            static string[] wantedFields = { "disease", "chromosome", "filename" };

        } // ExpressionDistributionByChromosomeMap

        public class SelectedGene
        {
            public SelectedGene(string Hugo_Symbol_, int nFlankingMutations_, int nRNAMutations_, Dictionary<int, int> tumorsByMutationCount)
            {
                Hugo_Symbol = Hugo_Symbol_;
                nFlankingMutations = nFlankingMutations_;
                nRNAMutations = nRNAMutations_;

                int maxMutations = 0;
                foreach (var countEntry in tumorsByMutationCount)
                {
                    maxMutations = Math.Max(maxMutations, countEntry.Key);
                }

                tumorsByCountOfNonSilentMutations = new int[maxMutations + 1];
                nMutations = 0;

                for (int i = 0; i <= maxMutations; i++)
                {
                    if (tumorsByMutationCount.ContainsKey(i))
                    {
                        tumorsByCountOfNonSilentMutations[i] = tumorsByMutationCount[i];
                        nMutations += tumorsByMutationCount[i] * i;
                    } else
                    {
                        tumorsByCountOfNonSilentMutations[i] = 0;
                    }
                }

                nTumors = tumorsByCountOfNonSilentMutations.Sum();
            }

            public int countOfTumorsByMutation(int nMutations)
            {
                if (nMutations >= tumorsByCountOfNonSilentMutations.Count())
                {
                    return 0;
                }

                return tumorsByCountOfNonSilentMutations[nMutations];
            }

            public static void SaveAllToFile(List<SelectedGene> selectedGenes, string filename)
            {
                var writer = CreateStreamWriterWithRetry(filename);

                writer.WriteLine("Hugo_Symbol\tnTumors\tnMutations\tnFlankingMutations\tnRNAMutations\tnTumorsWithOneMutation\tnTumorsWithMoreThanOneMutation\tTumorsByMutationCount");


                foreach (var selectedGene in selectedGenes)
                {
                    writer.Write(ConvertToExcelString(selectedGene.Hugo_Symbol) + "\t" + selectedGene.nTumors + "\t" + selectedGene.nMutations + "\t" + selectedGene.nFlankingMutations + "\t" + selectedGene.nRNAMutations + "\t" +
                        selectedGene.countOfTumorsByMutation(1) + "\t");
                    int nMoreThanOne = 0;
                    for (int i = 2; i < selectedGene.tumorsByCountOfNonSilentMutations.Count(); i++)
                    {
                        nMoreThanOne += selectedGene.countOfTumorsByMutation(i);
                    }
                    writer.Write(nMoreThanOne + "\t");
                    for (int i = 0; i < selectedGene.tumorsByCountOfNonSilentMutations.Count(); i++)
                    {
                        writer.Write(selectedGene.tumorsByCountOfNonSilentMutations[i] + ",");
                    }
                    writer.WriteLine();
                }

                writer.WriteLine("**done**");
                writer.Close();
            }

            static SelectedGene parser(Dictionary<string, int> fieldMappings, string[] fields)
            {
                var tumorsByMutationCountList = fields[fieldMappings["TumorsByMutationCount"]].Split(',');

                var tumorsByMutationCount = new Dictionary<int, int>();
                for (int i = 0; i < tumorsByMutationCountList.Count() - 1; i++) // -1 is because the list ends with a comma, so there's an empty field at the end that we ignore
                {
                    tumorsByMutationCount.Add(i, Convert.ToInt32(tumorsByMutationCountList[i]));
                }

                return new SelectedGene(ConvertToNonExcelString(fields[fieldMappings["Hugo_Symbol"]]), Convert.ToInt32(fields[fieldMappings["nFlankingMutations"]]), Convert.ToInt32(fields[fieldMappings["nRNAMutations"]]), tumorsByMutationCount);
            }

            public static List<SelectedGene> LoadFromFile(string filename)
            {
                var wantedFields = new List<string>();
                wantedFields.Add("Hugo_Symbol");
                wantedFields.Add("nTumors");
                wantedFields.Add("nMutations");
                wantedFields.Add("nRNAMutations");
                wantedFields.Add("nFlankingMutations");
                wantedFields.Add("TumorsByMutationCount");

                var inputFile = CreateStreamReaderWithRetry(filename);

                var file = new HeaderizedFile<SelectedGene>(inputFile, false, true, "", wantedFields);

                List<SelectedGene> result;
                file.ParseFile(parser, out result);

                inputFile.Close();

                return result;
            }

            public readonly string Hugo_Symbol;
            public readonly int nMutations;
            public readonly int nTumors;
            public readonly int nFlankingMutations;
            public readonly int nRNAMutations;
            public readonly int[] tumorsByCountOfNonSilentMutations;
        } // SelectedGene

        public static System.Net.WebClient getWebClient()
        {
            ServicePointManager.SecurityProtocol = SecurityProtocolType.Tls12;  // NIH uses TLS v 1.2, and for some reason the auto negotiate doesn't pick that up, so we'll just set it explicitly

            var webClient = new WebClient();

            //
            // Check the status of gdc to make sure that we're in the right version.
            //

            var statusSerializer = new DataContractJsonSerializer(typeof(ASETools.Status));
            ASETools.Status status = (ASETools.Status)statusSerializer.ReadObject(new MemoryStream(webClient.DownloadData(ASETools.urlPrefix + "status")));

            if (status.status != "OK")
            {
                Console.WriteLine("GDC returned not OK status of " + status.status + ", aborting");
                return null;
            }

            if (status.version != "1")
            {
                Console.WriteLine("GDC returned a version number that we don't understand (" + status.version + ", you probably need to update this program.");
                return null;
            }

            return webClient;
        }

        [DataContract]
        public class Status
        {
            [DataMember]
            public string commit = "";  // Default values are to void having Visual Studio generate a warning, since it can't see that they're set by the JSON serializer

            [DataMember]
            public string status = "";

            [DataMember]
            public string tag = "";

            [DataMember]
            public string version = "";
        }

        [DataContract]
        public class GDCPagination
        {
            [DataMember]
            public int count = 0;

            [DataMember]
            public string sort = "";

            [DataMember]
            public int from = 0;

            [DataMember]
            public int page = 0;

            [DataMember]
            public int total = 0;

            [DataMember]
            public int pages = 0;

            [DataMember]
            public string size = "";

            static public GDCPagination extractFromString(string inputString)
            {
                int paginationIndex = inputString.IndexOf("\"pagination\":");

                if (-1 == paginationIndex)
                {
                    return null;
                }

                int openCurlyBraceIndex = paginationIndex;

                int size = inputString.Count();

                while (openCurlyBraceIndex < size && inputString[openCurlyBraceIndex] != '{')
                {
                    openCurlyBraceIndex++;
                }

                if (openCurlyBraceIndex >= size)
                {
                    return null;
                }

                int currentIndex = openCurlyBraceIndex + 1;
                int openCurlyBraceCount = 1;

                while (currentIndex < size && openCurlyBraceCount > 0)
                {
                    if (inputString[currentIndex] == '}')
                    {
                        openCurlyBraceCount--;
                    }
                    else if (inputString[currentIndex] == '{')
                    {
                        openCurlyBraceCount++;
                    }
                    currentIndex++;
                }

                if (0 != openCurlyBraceCount)
                {
                    throw new FormatException();
                }


                var paginationString = inputString.Substring(openCurlyBraceIndex, currentIndex - openCurlyBraceIndex);
                var serializer = new DataContractJsonSerializer(typeof(GDCPagination));
                return (GDCPagination)serializer.ReadObject(new MemoryStream(Encoding.ASCII.GetBytes(paginationString)));
            }
        }

        [DataContract]
        public class GDCCaseProject
        {
            [DataMember]
            public string project_id;
        }

        [DataContract]
        public class GDCCase
        {
            [DataMember]
            public string[] sample_ids = { };

            [DataMember]
            public string[] portion_ids = { };

            [DataMember]
            public string updated_datetime = "";

            [DataMember]
            public string created_datetime = "";

            //            [DataMember]
            //            public string[] submitter_aliquot_ids = { };

            //            [DataMember]
            //            public string[] submitter_portion_ids = { };

            //            [DataMember]
            //            public string[] submitter_analyte_ids = { };

            //            [DataMember]
            //            public string[] analyte_ids = { };

            [DataMember]
            public string submitter_id = "";

            [DataMember]
            public string case_id = "";

            [DataMember]
            public string state = "";

            //            [DataMember]
            //            public string[] aliquot_ids = { };

            //            [DataMember]
            //            public string[] slide_ids = { };

            //            [DataMember]
            //            public string[] submitter_sample_ids = { };

            //            [DataMember]
            //            public string[] submitter_slide_ids = { };

            [DataMember]
            public GDCCaseProject project = null;

            public const string fields = //"sample_ids,portion_ids,updated_datetime,created_datetime,submitter_aliquot_ids,submitter_portion_ids,submitter_analyte_ids,analyte_ids,submitter_id,case_id,state,aliquot_ids,slide_ids,submitter_sample_ids,submitter_slide_ids,project.project_id";
                "sample_ids,portion_ids,updated_datetime,created_datetime,submitter_id,case_id,state,project.project_id";

            public string debugString()
            {
                var retVal = "sample IDs (" + sample_ids.Count() + "): {";
                foreach (var sampleId in sample_ids)
                {
                    retVal += sampleId + " ";
                }
                retVal += "} updated datetime: " + updated_datetime + ", created datetime:  " + created_datetime + ", submitter_id: " + submitter_id + ", case_id: " + case_id + ", state: " + state + ", project: " + project.project_id;

                return retVal;
            }
        }

        [DataContract]
        public class GDCHits<containedClass>
        {
            [DataMember]
            public containedClass[] hits = { };
        }

        [DataContract]
        public class GDCData<containedClass>
        {
            [DataMember]
            public GDCHits<containedClass> data = null;
        }

        [DataContract]
        public class GDCSamples
        {
            [DataMember]
            public string sample_type_id = "";

            [DataMember]
            public string sample_type = "";

            [DataMember]
            public string sample_id = "";
        }

        [DataContract]
        public class GDCFileCases
        {
            [DataMember]
            public GDCSamples[] samples = { };
        }

        [DataContract]
        public class GDCFileAnalysis
        {
            [DataMember]
            public string workflow_link = "";
        }

        [DataContract]
        public class GDCFile
        {
            [DataMember]
            public string data_type = "";

            [DataMember]
            public string updated_datetime = "";

            [DataMember]
            public string created_datetime = "";

            [DataMember]
            public string file_name = "";

            [DataMember]
            public string md5sum = "";

            [DataMember]
            public string data_format = "";

            [DataMember]
            public string[] acl = { };

            [DataMember]
            public string access = "";

            [DataMember]
            public string platform = "";

            [DataMember]
            public string state = "";

            [DataMember]
            public string file_id = "";

            [DataMember]
            public string data_category = "";

            [DataMember]
            public long file_size = 0;

            [DataMember]
            public string submitter_id = "";

            [DataMember]
            public string type = "";

            [DataMember]
            public string file_state = "";

            [DataMember]
            public string experimental_strategy = "";

            [DataMember]
            public GDCFileCases[] cases = { };

            [DataMember]
            public GDCFileAnalysis analysis = null;

            public static GDCFile selectNewestUpdated(GDCFile one, GDCFile two)
            {
                if (null == one) return two;
                if (null == two) return one;

                if (Convert.ToDateTime(one.updated_datetime) < Convert.ToDateTime(two.updated_datetime))
                {
                    return two;
                }
                else
                {
                    return one;
                }
            }
        } // GDCFile


        public static string generateFilterList(List<string> itemsToAdd)
        {
            var outputString = "[";

            bool addedOneAlready = false;
            foreach (var item in itemsToAdd)
            {
                if (addedOneAlready)
                {
                    outputString += ",";
                }
                else
                {
                    addedOneAlready = true;
                }
                outputString += "\"" + item + "\"";
            }

            outputString += "]";

            return outputString;
        }

        public static string ElapsedTimeInSeconds(Stopwatch stopwatch)
        {
            return "" + (stopwatch.ElapsedMilliseconds + 500) / 1000 + "s";
        }

        public static string ShareFromPathname(string pathname)
        {
            //
            // A pathname is of the form \\computer\share\...
            // Get up to share, discard the rest.
            //
            string[] fields = pathname.Split('\\');
            if (fields.Count() < 4 || fields[0] != "" || fields[1] != "")
            {
                throw new Exception("ShareFromPathname: can't parse pathname " + pathname);
            }

            string share = "";
            for (int i = 1; i < 4; i++) // Skip 0 because we know it's the empty string, and it makes the \ part easier.
            {
                share += @"\" + fields[i];
            }

            return share;
        }

        public static string ComputerFromPathname(string pathname)
        {
            //
            // A pathname is of the form \\computer\share\...
            // Return the computer name
            //

            var fields = pathname.Split('\\');
            if (fields.Count() < 4 || fields[0] != "" || fields[1] != "")
            {
                throw new Exception("ComputerFromPathname: can't parse pathname " + pathname);
            }

            return fields[2];
        }

        public static string PathnameToLinuxPathname(string pathname)
        {
            return pathname.Replace('\\', '/'); // Just switch the directory separators.  You're responsible for the UNC path and whatnot.
        }

        public static string PathnameWithoutUNC(string pathname)
        {
            var fields = pathname.Split('\\');
            if (fields.Count() < 4 || fields[0] != "" || fields[1] != "")
            {
                throw new Exception("PathnameWithoutUNC: can't parse pathname " + pathname);
            }

            string retVal = "";
            for (int i = 4; i < fields.Count(); i++)
            {
                retVal += @"\" + fields[i];
            }

            return retVal;
        } // PathnameWithoutUNC

        public static string GetFileNameFromPathname(string pathname, bool excludeExtension = false)
        {
            string lastComponent;
            if (pathname.LastIndexOf('\\') == -1)
            {
                lastComponent = pathname;
            }
            else
            {
                lastComponent = pathname.Substring(pathname.LastIndexOf('\\') + 1);
            }

            if (!excludeExtension || lastComponent.LastIndexOf('.') == -1)
            {
                return lastComponent;
            }

            return lastComponent.Substring(0, lastComponent.LastIndexOf('.'));
        }

        public static string GetDirectoryFromPathname(string pathname)
        {
            if (!pathname.Contains('\\'))
            {
                throw new FormatException();
            }

            return pathname.Substring(0, pathname.LastIndexOf('\\'));
        }

        public static string GetFileIdFromPathname(string pathname)
        {
            string filename = GetFileNameFromPathname(pathname);

            if (filename.Count() < FileIdLength)
            {
                throw new FormatException();
            }

            return filename.Substring(0, FileIdLength);
        }

        static readonly int FileIdLength = "01c5f902-ec3c-4eb7-9e38-1d29ae6ab959".Count();
        //
        // For filenames of the from ...\analysis-id\filename
        //
        static public string GetAnalysisIdFromPathname(string pathname)
        {
            if (!pathname.Contains('\\'))
            {
                Console.WriteLine("ASETools.GetAnalysisIdFromPathname: invalid pathname input " + pathname);
                return "";
            }

            string directoryName;
            if (pathname.LastIndexOf('\\') != pathname.Count())
            {
                //
                // Does not end in a trailing \.
                //
                directoryName = pathname.Substring(0, pathname.LastIndexOf('\\'));
            }
            else
            {
                //
                // Ends in a trailing backslash, just trim it off.
                //
                directoryName = pathname.Substring(0, pathname.Count() - 1);
            }

            if (directoryName.Contains('\\'))
            {
                directoryName = directoryName.Substring(directoryName.LastIndexOf('\\') + 1);
            }

            if (directoryName.Count() != 36)
            {
                Console.WriteLine("ASETools.GetAnalysisIdFromPathname: invalid pathname input does not include analysis id as last directory " + pathname);
                return "";
            }


            return directoryName;
        }

        public class DownloadedFile
        {
            public readonly string file_id;
            public readonly FileInfo fileInfo;
            public readonly string storedMD5;    // This is the null string for files for which we don't know the sum (i.e., they haven't yet been computed or we don't know the right answer).
            public readonly FileInfo md5FileInfo;

            public DownloadedFile(string file_id_, string pathname, string storedMd5Sum_, string md5FilePathname)
            {
                file_id = file_id_;
                storedMD5 = storedMd5Sum_;

                fileInfo = new FileInfo(pathname);

                if (md5FilePathname == null || md5FilePathname == "")
                {
                    md5FileInfo = null;
                }
                else
                {
                    md5FileInfo = new FileInfo(md5FilePathname);
                }
            }
        } // DownloadedFile

        public static string getReadsAtTentativeSelectedVariantsExtension(bool tumor, bool dna)
        {
            if (tumor)
            {
                if (dna) return tumorDNAReadsAtTentativeSelectedVariantsExtension;
                return tumorRNAReadsAtTentativeSelectedVariantsExtension;
            }
            if (dna) return normalDNAReadsAtTentativeSelectedVariantsExtension;
            return normalRNAReadsAtTentativeSelectedVariantsExtension;
        } // getReadsAtTentativeSelectedVariantsExtension


        public const string tumorRNAAllcountExtension = ".allcount.gz";
        public const string normalRNAAllcountExtension = ".normal_rna_allcount.gz";
        public const string normalDNAAllcountExtension = ".normal_dna_allcount.gz";
        public const string tumorDNAAllcountExtension = ".tumor_dna_allcount.gz";
        public const string tentativeSelectedVariantsExtension = ".tentativeSelectedVariants";
        public const string selectedVariantsExtension = ".selectedVariants";
        public const string tentativeAnnotatedSelectedVariantsExtension = ".tentativeAnnotatedSelectedVariants";
        public const string annotatedSelectedVariantsExtension = ".annotatedSeletedVariants";
        public const string regionalExpressionExtension = ".regional_expression.txt";
        public const string geneExpressionExtension = ".gene_expression.txt";
        public const string tumorAlleleSpecificGeneExpressionExtension = ".tumor-allele-specific_gene_expression.txt";
        public const string normalAlleleSpecificGeneExpressionExtension = ".normal-allele-specific_gene_expression.txt";
        public const string tumorRNAReadsAtTentativeSelectedVariantsExtension = ".tumor-rna-reads-at-tentative-selected-variants.txt";
        public const string tumorRNAReadsAtTentativeSelectedVariantsIndexExtension = ".tumor-rna-reads-at-tentative-selected-variants.txt.index";
        public const string normalRNAReadsAtTentativeSelectedVariantsExtension = ".normal-rna-reads-at-tentative-selected-variants.txt";
        public const string normalRNAReadsAtTentativeSelectedVariantsIndexExtension = ".normal-rna-reads-at-tentative-selected-variants.txt.index";
        public const string normalDNAReadsAtTentativeSelectedVariantsExtension = ".normal-dna-reads-at-tentative-selected-variants.txt";
        public const string normalDNAReadsAtTentativeSelectedVariantsIndexExtension = ".normal-dna-reads-at-tentative-selected-variants.txt.index";
        public const string tumorDNAReadsAtTentativeSelectedVariantsExtension = ".tumor-dna-reads-at-tentative-selected-variants.txt";
        public const string tumorDNAReadsAtTentativeSelectedVariantsIndexExtension = ".tumor-dna-reads-at-tentative-selected-variants.txt.index";
        public const string vcfExtension = ".vcf";
        public const string extractedMAFLinesExtension = ".extracted_maf_lines.txt";
        public const string allMAFLinesExtension = ".all_maf_lines";
        public const string tumorDNAGeneCoverageExtension = ".tumor_dna_gene_coverage.txt";
        public const string normalDNAMappedBaseCountExtension = ".normal_dna_mapped_base_count.txt";
        public const string tumorDNAMappedBaseCountExtension = ".tumor_dna_mapped_base_count.txt";
        public const string normalRNAMappedBaseCountExtension = ".normal_rna_mapped_base_count.txt";
        public const string tumorRNAMappedBaseCountExtension = ".tumor_rna_mapped_base_count.txt";
        public const string selectedVariantCountByGeneExtension = ".selected_variant_count_by_gene.txt";
        public const string bonferroniExtension = "_bonferroni.txt";
        public const string regulatoryMutationsNearMutationsExtension = ".regulatory_mutations_near_mutations.txt";
        public const string expressionByGeneExtension = ".expressionByGene";
        public const string isoformReadCountsExtension = ".isoformReadCounts";
        public const string compressedVCFExtension = ".vcf.gz";
        public const string caseMetadataExtension = ".case_metadata.txt";
        public const string tentativeASVsWithoutCNVsExtension = ".tentative_asvs_without_cnvs.txt";

        public const string scatterGraphsSummaryFilename = "_summary.txt";
        public const string mannWhitneyFilename = "_MannWhitney.txt";
        public const string genesWithSelectedVariantsFilename = "GenesWithSelectedVariantCounts.txt";
        public const string heatMapFilename = "AlleleSpecificExpressionHeatMap.txt";
        public const string regionalMethylationExtension = ".regional_methylation.txt";

        public const string selectedRegulatoryMAFLinesExtension = ".regulatory_maf_lines";
        public const string variantPhasingExtension = ".variant_phasing.txt";
        public const string vcfStatisticsExtension = ".vcf_statistics";

        public const string tumorHeatMapFilename = "TumorAlleleSpecificExpressionHeatMap.txt";
        public const string normalHeatMapFilename = "NormalAlleleSpecificExpressionHeatMap.txt";
        public const string tumorHeatMapHistogramFilename = "TumorAlleleSpecificExpressionHeatMapHistograms.txt";
        public const string normalHeatMapHistogramFilename = "NormalAlleleSpecificExpressionHeatMapHistograms.txt";
        public const string ASEConsistencyFilename = "ASEConsistency.txt";
        public const string GenesByFunnyASECountFilename = "ASEInconsistencyByGene.txt";
        public const string AllSignificantResultsFilenameBase = "AllSignificantResults";
        public const string AllSignificantResultsFilename = AllSignificantResultsFilenameBase + ".txt";
        public const string pValueHistogramFilename = "pValueHistogram.txt";
        public const string AlleleSpecificExpressionDistributionByMutationCountFilenameBase = "AlleleSpecificExpressionDistributionByMutationCount";
        public const string HistogramsForSignficantResultsFilename = "HistogramsForSignificantResults.txt";
        public const string PerGeneRNARatioFilename = "PerGeneRNARatio.txt";
        public const string ASEMapFilename = "ASEMap.txt";
        public const string ASEDifferenceMapFilename = "ASEDifferenceMap.txt";
        public const string PerGeneASEMapFilename = "ASEMap-PerGene.txt";
        public const string PerCaseASEFilename = "PerCaseASE.txt";
        public const string UncorrectedOverallASEFilename = "UncorrectedOverallASE.txt";
        public const string CorrectedOverallASEFilename = "CorrectedOverallASE.txt";
        public const string ASECorrectionFilename = "ase_correction.txt";
        public const string TumorRNAReadDepthDistributionFilename = "TumorRNAReadDepth.txt";
        public const string TumorGermlineASEDistributionFilename = "TumorGermlineASEDistribution.txt";
        public const string allLociExtension = "-all-loci.sam";
        public const string allLociAlignedExtension = "-all-loci-aligned.sam";
        public const string transcriptomeFastaFilename = "transcriptome.fasta";
        public const string generatedTranscriptomeIndexName = "TranscriptomeSNAPIndex";
        public const string AllSitesReadDepthFilename = "AllSitesReadDepth.txt";
        public const string Expression_distribution_filename_base = "expression_distribution_";
        public const string Expression_filename_base = "expression_";
        public const string annotated_scatter_graph_filename_extension = "_annotated_scatter_lines.txt";
        public const string raw_median_data_extension = "_raw_median_data.txt";
        public const string annotated_scatter_graphs_histogram_filename = "_annotated_scatter_graphs_histograms.txt";
        public const string vaf_histogram_filename = "VAFHistograms.txt";
        public const string simulatedResultsFilename = "SimulatedASEError.txt";
        public const string mappedBaseCountDistributionFilename = "MappedBaseCountDistribution.txt";
        public const string annotatedBEDLinesExtension = ".annotated_bed_lines.txt";
        public const string annotatedGeneHancerLinesExtension = ".annotated_gene_hancers.txt";
        public const string cisRegulatoryMutationsByVAFFilename = "CisRegulatoryMutationsByVAF.txt";
        public const string cisRegulatoryMutationsInKnownRegionsExtension = ".cis_regulatory_mutations_in_known_regions.txt";
        public const string readStatisticsExtension = ".read_statistics";
        public const string snapRealignedNormalDNAExtension = ".snap-normal-dna.bam";
        public const string snapRealignedNormalDNABaiExtension = ".snap-normal-dna.bam.bai";
        public const string snapRealignedTumorDNAExtension = ".snap-tumor-dna.bam";
        public const string snapRealignedTumorDNABaiExtension = ".snap-tumor-dna.bam.bai";
        public const string snapRealignedNormalDNAStaticticsExtension = ".snap-normal-dna-statictics.txt";
        public const string snapRealignedTumorDNAStatisticsExtension = ".snap-tumor-dna-statictics.txt";
        public const string bowtieRealignedNormalDNAExtension = ".bowtie-normal-dna.bam";
        public const string bowtieRealignedNormalDNABaiExtension = ".bowtie-normal-dna.bam.bai";
        public const string bowtieRealignedTumorDNAExtension = ".bowtie-tumor-dna.bam";
        public const string bowtieRealignedTumorDNABaiExtension = ".bowtie-tumor-dna.bam.bai";
        public const string bowtieRealignedNormalDNAStaticticsExtension = ".bowtie-normal-dna-statictics.txt";
        public const string bowtieRealignedTumorDNAStatisticsExtension = ".bowtie-tumor-dna-statictics.txt";
        public const string normalFastqExtension = ".normal_1.fastq";
        public const string normalSecondEndFastqExtension = ".normal_2.fastq";
        public const string tumorFastqExtension = ".tumor_1.fastq";
        public const string tumorSecondEndFastqExtension = ".tumor_2.fastq";

        public const string MethylationDistributionsFilename = "MethylationDistributions.txt";
        public const string MethylationCorrelationsFilename = "MethylationCorrelations.txt";
        public const string MethylagtionCorrelationsHistogramFilename = "MethylationCorrelationsHistogram.txt";
        public const string IsoformBalanceFilenameBase = "IsoformBalance";
        public const string IsoformBalancePValueHistogramFilename = "IsoformBalance-pValueHistogram.txt";
        public const string ReadLengthHistogramFilename = "ReadLengthHistograms.txt";
        public const string variantSelectionSummaryFilename = "VariantSelectionSummary.txt";
        public const string DistanceBetweenMutationsFilename = "DistanceBetweenMutations.txt";
        public const string SingleReadPhasingFilename = "SingleReadPhasing.txt";
        public const string GeneScatterGraphLinesWithPercentilesPrefix = "GeneScatterGraphsWithPercentiles_";   // Followed by disease _ chromosome
        public const string ExpressionDistrbutionByChromosomeDirectory = @"expression_distribution_by_chromosome\";
        public const string UniparentalDisomyFilename = "UniparentalDisomy.txt";
        public const string UniparentalDisomyHistogramsFilename = "UniparentalDisomyHistograms.txt";
        public const string PhasingForNearbyVariantsFilename = "PhasingForNearbyVariants.txt";
        public const string VCFStatisticsFilename = "VCFStatistics.txt";
        public const string miRNAExpressionSummaryFilename = "miRNASummary.txt";
        public const string miRNAExpressionPValueHistogramFilename = "miRNAPValueHistogram.txt";
        public const string CaseMetadataSummaryFilename = "CaseMetadataSummary.txt";

        public const string basesInKnownCodingRegionsPrefix = "BasesInKnownCodingRegions_";



        public class DerivedFile
        {
            public readonly string derived_from_file_id;
            public readonly Type type;
            public readonly FileInfo fileinfo;
            public readonly string case_id;

            public DerivedFile(string pathname, string case_id_)
            {
                fileinfo = new FileInfo(pathname);
                type = Type.Unknown;

                var filename = GetFileNameFromPathname(pathname).ToLower();
                case_id = case_id_;

                if (filename.Count() <= FileIdLength + 1) {
                    Console.WriteLine("Derived file has too short of a name: " + pathname);
                    derived_from_file_id = "";
                    return;
                }

                derived_from_file_id = filename.Substring(0, FileIdLength);

                var extension = filename.Substring(FileIdLength).ToLower();

                foreach (var field in Case.AllFields)
                {
                    if (field.type != Type.Unknown && extension == field.extension.ToLower())
                    {
                        type = field.type;
                    }
                }

                if (Type.Unknown == type)
                {
                    Console.WriteLine("Derived file with unknown extension: " + pathname);
                }
            } // DerviedFile.DerivedFile


            public enum Type { Unknown, NormalRNAAllcount, TumorRNAAllcount, NormalDNAAllcount, TumorDNAAllcount, RegionalExpression, GeneExpression, TumorDNAGeneCoverage, TentativeSelectedVariants,
                SelectedVariants, NormalDNAReadsAtSelectedVariants, NormalDNAReadsAtSelectedVariantsIndex, TumorDNAReadsAtSelectedVariants, TumorDNAReadsAtSelectedVariantsIndex, TumorRNAReadsAtSelectedVariants,
                TumorRNAReadsAtSelectedVariantsIndex, NormalRNAReadsAtSelectedVariants, NormalRNAReadsAtSelectedVariantsIndex, AnnotatedSelectedVariants, NormalAlleleSpecificGeneExpression, TumorAlleleSpecificGeneExpression, VCF, ExtractedMAFLines, AllMAFLines,
                NormalDNAMappedBaseCount, TumorDNAMappedBaseCount, NormalRNAMappedBaseCount, TumorRNAMappedBaseCount, SelectedVariantCountByGene, SelectedRegulatoryMAFLines, AnnotatedRegulatoryRegions, RegulatoryMutationsNearMutations,
                AnnotatedGeneHancer, ExpressionByGene, TentativeAnnotatedSelectedVariants, IsoformReadCounts, CompressedVCF, CaseMetadata, TentativeASVsWithoutCNVs, VariantPhasing, VCFStatistics, ReadStatictics,
                SnapRealignedNormalDNA, SnapRealignedNormalDNABai, SnapRealignedTumorDNA, SnapRealignedTumorDNABai, SnapNormalDNAStatictics, SnapTumorDNAStatictics,
                BowtieRealignedNormalDNA, BowtieRealignedNormalDNABai, BowtieRealignedTumorDNA, BowtieRealignedTumorDNABai, BowtieNormalDNAStatictics, BowtieTumorDNAStatictics, NormalDNAFastq, NormalDNAFastqSecondEnd, TumorDNAFastq, TumorDNAFastqSecondEnd
            };
        } // DerivedFile

        public class ScannedFilesystem
        {
            public ulong totalFreeBytes = 0;
            public ulong totalStorageBytes = 0;
            public ulong totalBytesInDownloadedFiles = 0;
            public ulong totalBytesInDerivedFiles = 0;
            public int nDownloadedFiles = 0;
            public int nDerivedFiles = 0;

            public List<DownloadedFile> downloadedFiles = new List<DownloadedFile>();
            public List<DerivedFile> derivedFiles = new List<DerivedFile>();

            public void mergeWith(ScannedFilesystem peer)
            {
                totalFreeBytes += peer.totalFreeBytes;
                totalStorageBytes += peer.totalStorageBytes;
                totalBytesInDownloadedFiles += peer.totalBytesInDownloadedFiles;
                totalBytesInDerivedFiles += peer.totalBytesInDerivedFiles;
                nDownloadedFiles += peer.nDownloadedFiles;
                nDerivedFiles += peer.nDerivedFiles;

                downloadedFiles.AddRange(peer.downloadedFiles);
                derivedFiles.AddRange(peer.derivedFiles);
            }
        }

        static void ScanOneFilesystem(Configuration configuration, string downloadedFilesDirectory, ScannedFilesystem state, Stopwatch stopwatch, int directoryFieldLength)
        {
            if (!Directory.Exists(downloadedFilesDirectory))
            {
                Console.WriteLine("Warning: can't find download directory " + downloadedFilesDirectory + ".  Ignoring.");
                return;
            }

            int nDownloadedFiles = 0;
            int nDerivedFiles = 0;
            ulong totalBytesInDownloadedFiles = 0;
            ulong totalBytesInDerivedFiles = 0;

            ulong freeBytesAvailable, totalBytes, totalNumberOfFreeBytes;

            var downloadedFiles = new List<DownloadedFile>();
            var derivedFiles = new List<DerivedFile>();

            try
            {
                foreach (var subdir in Directory.EnumerateDirectories(downloadedFilesDirectory))
                {
                    var file_id = GetFileNameFromPathname(subdir).ToLower();    // The directory is the same as the file id.

                    if (file_id.Count() != GuidStringLength)
                    {
                        lock (state)
                        {
                            Console.WriteLine("Found subdirectory of a downloaded files dierectory whose name doesn't appear to be a guid, ignoring: " + subdir);
                        }
                        continue;
                    }

                    //
                    // Look through the subdirectory to find the downloaded file and also any .md5 versions of it.
                    // Also look for decompressed MAF files
                    //
                    string candidatePathname = null;
                    string md5Pathname = null;
                    bool sawBAI = false;
                    bool sawBAM = false;
                    string unsortedBAMPathmame = null;
                    string extraMAF = null;
                    foreach (var pathname in Directory.EnumerateFiles(subdir))
                    {
                        var filename = GetFileNameFromPathname(pathname).ToLower();
                        if (filename == "annotations.txt")
                        {
                            continue;
                        }

                        if (filename.EndsWith(".bai"))
                        {
                            sawBAI = true;
                            continue;
                        }


                        if (filename.EndsWith(".md5"))
                        {
                            if (null != md5Pathname)
                            {
                                Console.WriteLine("Directory " + subdir + " contains more than one file with extension .md5");
                                continue;
                            }
                            md5Pathname = pathname;
                            continue;
                        }

                        if (filename.EndsWith(".unsorted-bam") && configuration.isBeatAML)
                        {
                            unsortedBAMPathmame = pathname;
                            continue;
                        }

                        if (filename.EndsWith(".maf"))
                        {
                            extraMAF = pathname;
                            continue;
                        }

                        if (null != candidatePathname)
                        {
                            Console.WriteLine("Found more than one file candidate in " + subdir);
                            continue;
                        }

                        candidatePathname = pathname;

                        if (filename.EndsWith(".bam"))
                        {
                            sawBAM = true;
                        }
                    }

                    if (candidatePathname == null && unsortedBAMPathmame != null)
                    {
                        candidatePathname = unsortedBAMPathmame;    // Until we sort it, we'll just use that.
                    }

                    if (candidatePathname == null && md5Pathname != null)
                    {
                        Console.WriteLine("Found md5 file in directory without accompanying downloaded file (?): " + md5Pathname);
                        continue;
                    }

                    if (sawBAM != sawBAI && (!configuration.isBeatAML || !sawBAM)) // Some of the BeatAML files aren't indexed, so we run a stage for that.
                    {
                        Console.WriteLine("Downloaded file directory " + subdir + " has exactly one BAM and BAI file.");
                        continue;
                    }

                    string md5Value = "";
                    if (md5Pathname != null)
                    {
                        if (md5Pathname.ToLower() != candidatePathname.ToLower() + ".md5")
                        {
                            Console.WriteLine("md5 file has wrong basename: " + md5Pathname);
                        }
                        else
                        {
                            try
                            {
                                var lines = File.ReadAllLines(md5Pathname); // Don't use the with retry version, because we don't want to get stuck behind a worker hashing a file.
                                if (lines.Count() != 1 || lines[0].Count() != 32)
                                {
                                    Console.WriteLine("md5 file " + md5Pathname + " has a non-md5 content.  Ignoring.");
                                }
                                else
                                {
                                    md5Value = lines[0];
                                }
                            }
                            catch (IOException)
                            {
                                Console.WriteLine("Skipping md5 file " + md5Pathname + " because of an IO error.  It's most likely being created now.");
                            }
                        }
                    }


                    if (null == candidatePathname)
                    {
                        Console.WriteLine("Couldn't find downloaded file in directory " + subdir);
                    }
                    else
                    {
                        nDownloadedFiles++;
                        var downloadedFile = new DownloadedFile(file_id, candidatePathname, md5Value, md5Pathname);
                        downloadedFiles.Add(downloadedFile);

                        totalBytesInDownloadedFiles += (ulong)downloadedFile.fileInfo.Length;
                    }

                    if (extraMAF != null)
                    {
                        nDownloadedFiles++;
                        var filename = GetFileNameFromPathname(extraMAF);
                        var disease = filename.Substring(0, filename.Length - 4);   // -1 is for .maf.  These files are named disease.maf
                        var downloadedFile = new DownloadedFile(disease + ".maf", extraMAF, "", "");

                        downloadedFiles.Add(downloadedFile);

                        totalBytesInDownloadedFiles += (ulong)downloadedFile.fileInfo.Length;
                    }
                } // foreach subdir of downloaded

                string derivedFilesDirectory = downloadedFilesDirectory + @"..\" + configuration.derivedFilesDirectory;
                foreach (var derivedCase in Directory.EnumerateDirectories(derivedFilesDirectory))
                {
                    var caseId = ASETools.GetFileNameFromPathname(derivedCase).ToLower();

                    if (caseId.Count() != GuidStringLength)
                    {
                        lock (state)
                        {
                            Console.WriteLine("Found subdirectory of a derived files dierectory whose name doesn't appear to be a guid, ignoring: " + derivedCase);
                        }
                        continue;
                    }

                    foreach (var derivedFilePathname in Directory.EnumerateFiles(derivedCase).Select(x => x.ToLower()))
                    {
                        var derivedFile = new DerivedFile(derivedFilePathname, caseId);
                        nDerivedFiles++;
                        totalBytesInDerivedFiles += (ulong)derivedFile.fileinfo.Length;
                        derivedFiles.Add(derivedFile);
                    }
                }

                GetDiskFreeSpaceEx(downloadedFilesDirectory, out freeBytesAvailable, out totalBytes, out totalNumberOfFreeBytes);
            } catch (Exception e)
            {
                Console.WriteLine("Exception in ScanOneFilesystem for " + downloadedFilesDirectory);
                throw e;
            }

            lock (state)
            {
                Console.WriteLine(String.Format("{0," + directoryFieldLength + "}", downloadedFilesDirectory) + " " + String.Format("{0,16}", "" + nDownloadedFiles + " (" + SizeToUnits(totalBytesInDownloadedFiles) + "B)") + " " +
                    String.Format("{0,13}", "" + nDerivedFiles + " (" + SizeToUnits(totalBytesInDerivedFiles) + "B)") + " " + String.Format("{0,10}", SizeToUnits(freeBytesAvailable) + "B") + " " +
                    String.Format("{0,9}", ElapsedTimeInSeconds(stopwatch)));

                state.nDownloadedFiles += nDownloadedFiles;
                state.nDerivedFiles += nDerivedFiles;
                state.totalBytesInDownloadedFiles += totalBytesInDownloadedFiles;
                state.totalBytesInDerivedFiles += totalBytesInDerivedFiles;
                state.totalFreeBytes += freeBytesAvailable;
                state.totalStorageBytes += totalBytes;
                state.downloadedFiles.AddRange(downloadedFiles);
                state.derivedFiles.AddRange(derivedFiles);
            }
        }

        //
        // Grovel the file system(s) to look for downloaded and derived files.  Returns a dictionary that maps from file_id -> DownloadedFile object.
        //
        public static void ScanFilesystems(Configuration configuration, out Dictionary<string, DownloadedFile> downloadedFiles, out Dictionary<string, List<DerivedFile>> derivedFiles, out Dictionary<string, ScannedFilesystem> fileSystems)
        {
            var stopwatch = new Stopwatch();
            stopwatch.Start();

            int longestDirectoryNameLength = 0;

            foreach (var directory in configuration.dataDirectories)
            {
                longestDirectoryNameLength = Math.Max(longestDirectoryNameLength, directory.Count());
            }

            const string directoryHeader = "Directory";
            int paddingLength = Math.Max(0, longestDirectoryNameLength - directoryHeader.Count());

            Console.Write("Directory");
            for (int i = 0; i < paddingLength; i++)
            {
                Console.Write(" ");
            }

            Console.WriteLine(" Downloaded Files Derived Files Free Size  Scan Time");
            for (int i = 0; i < directoryHeader.Count() + paddingLength; i++)
            {
                Console.Write("-");
            }
            Console.WriteLine(" ---------------- ------------- ---------- ---------");

            var perDirectoryState = new Dictionary<string, ScannedFilesystem>();
            configuration.dataDirectories.ForEach(_ => perDirectoryState.Add(_, new ScannedFilesystem()));

            var threads = new List<Thread>();
            
            foreach (var directory in configuration.dataDirectories)
            {
                threads.Add(new Thread(() => ScanOneFilesystem(configuration, directory, perDirectoryState[directory], stopwatch, paddingLength + directoryHeader.Count())));
            } // foreach data directory

            threads.ForEach(t => t.Start());
            threads.ForEach(t => t.Join());

            fileSystems = perDirectoryState;
            var state = new ScannedFilesystem();

            downloadedFiles = new Dictionary<string, DownloadedFile>();
            derivedFiles = new Dictionary<string, List<DerivedFile>>();

            foreach (var downloadedFile in state.downloadedFiles)
            {
                if (downloadedFiles.ContainsKey(downloadedFile.file_id))
                {
                    Console.WriteLine("Found multiple directories with the downloaded file id " + downloadedFile.file_id + ": " + downloadedFile.fileInfo.FullName + " and " + downloadedFiles[downloadedFile.file_id].fileInfo.FullName);
                    continue;
                }
                downloadedFiles.Add(downloadedFile.file_id, downloadedFile);
            }

            foreach (var derivedFile in state.derivedFiles)
            {
                if (!derivedFiles.ContainsKey(derivedFile.case_id))
                {
                    derivedFiles.Add(derivedFile.case_id, new List<DerivedFile>());
                }

                derivedFiles[derivedFile.case_id].Add(derivedFile);
            }

            Console.WriteLine("Scanned " + configuration.dataDirectories.Count() + " data directories in " + ElapsedTimeInSeconds(stopwatch) + ", containing " + state.nDownloadedFiles + " (" + SizeToUnits(state.totalBytesInDownloadedFiles) +
                "B) downloaded and " + state.nDerivedFiles + " (" + SizeToUnits(state.totalBytesInDerivedFiles) + "B) derived files.  " + SizeToUnits(state.totalFreeBytes) + "B remain free of " + SizeToUnits(state.totalStorageBytes) + "B total.");
        } // ScanFilesystems

        public class MAFInfo
        {
            public MAFInfo(string file_id_, string md5Sum_, string filename_)
            {
                file_id = file_id_;
                md5Sum = md5Sum_;
                filename = filename_;
            }
            public static MAFInfo fromSaveFileLine(string saveFileLine)
            {
                var fields = saveFileLine.Split('\t');

                if (fields.Count() != 3)
                {
                    Console.WriteLine("Incorrect number of fields in maf configuration file line: " + saveFileLine);
                    return null;
                }

                return new MAFInfo(fields[0], fields[1], fields[2]);
            }

            public static Dictionary<string, MAFInfo> LoadMAFManifest(string filename)
            {
                if (!File.Exists(filename))
                {
                    return null;
                }

                var lines = ReadAllLinesWithRetry(filename);

                if (lines.Count() < 2)
                {
                    Console.WriteLine("MAF manifest file " + filename + " has too few lines.  Ignoring.");
                    return null;
                }

                var headerPrefix = "MAF Manifest v1.0 generated at ";
                if (lines[0].Count() < headerPrefix.Count() || lines[0].Substring(0, headerPrefix.Count()) != headerPrefix) {
                    Console.WriteLine("Corrupt or unrecognized version in maf manifest header, ignoring: " + lines[0]);
                    return null;
                }

                var retVal = new Dictionary<string, MAFInfo>();

                if (lines[lines.Count() - 1] != "**done**") {
                    Console.WriteLine("maf manifest file " + filename + " does not end with **done**, and so is probably truncated.  Ignoring.");
                    return null;
                }

                for (int i = 1; i < lines.Count() - 1; i++) {
                    var mafInfo = MAFInfo.fromSaveFileLine(lines[i]);

                    if (null == mafInfo) {
                        Console.WriteLine("failed to parse maf manifest line " + lines[i] + ".  Ignoring configuration file.");
                        return null;
                    }

                    if (retVal.ContainsKey(mafInfo.file_id)) {
                        Console.WriteLine("Duplicate file id " + mafInfo.file_id + " in maf manifest file " + filename);
                        return null;
                    }

                    retVal.Add(mafInfo.file_id, mafInfo);
                }

                return retVal;
            }

            public string file_id;
            public string md5Sum;
            public string filename;
        } // MAFInfo

        //
        // Code for GetVolumeFreeSpace adapted from http://stackoverflow.com/questions/14465187/get-available-disk-free-space-for-a-given-path-on-windows
        //
        [DllImport("kernel32.dll", SetLastError = true, CharSet = CharSet.Auto)]
        [return: MarshalAs(UnmanagedType.Bool)]
        static extern bool GetDiskFreeSpaceEx(string lpDirectoryName,
           out ulong lpFreeBytesAvailable,
           out ulong lpTotalNumberOfBytes,
           out ulong lpTotalNumberOfFreeBytes);
        public static ulong GetVolumeFreeSpace(string pathname)
        {
            ulong FreeBytesAvailable;
            ulong TotalNumberOfBytes;
            ulong TotalNumberOfFreeBytes;

            if (!GetDiskFreeSpaceEx(pathname,
                                        out FreeBytesAvailable,
                                        out TotalNumberOfBytes,
                                        out TotalNumberOfFreeBytes))
            {
                throw new System.ComponentModel.Win32Exception();
            }

            return FreeBytesAvailable;
        } // GetVolumeFreeSpace

        //
        // Code to get computer memory size adapted from http://stackoverflow.com/questions/105031/how-do-you-get-total-amount-of-ram-the-computer-has
        //
        [StructLayout(LayoutKind.Sequential, CharSet = CharSet.Auto)]
        private class MEMORYSTATUSEX
        {
            public uint dwLength;
            public uint dwMemoryLoad;
            public ulong ullTotalPhys;
            public ulong ullAvailPhys;
            public ulong ullTotalPageFile;
            public ulong ullAvailPageFile;
            public ulong ullTotalVirtual;
            public ulong ullAvailVirtual;
            public ulong ullAvailExtendedVirtual;
            public MEMORYSTATUSEX()
            {
                this.dwLength = (uint)Marshal.SizeOf(typeof(MEMORYSTATUSEX));
            }
        }

        [DllImport("kernel32.dll", SetLastError = true, CharSet = CharSet.Auto)]
        [return: MarshalAs(UnmanagedType.Bool)]
        static extern bool GlobalMemoryStatusEx([In, Out] MEMORYSTATUSEX lpBuffer);

        static public ulong GetTotalComputerMemory()
        {
            MEMORYSTATUSEX status = new MEMORYSTATUSEX();
            if (!GlobalMemoryStatusEx(status))
            {
                return 0;
            }

            return status.ullTotalPhys;
        }

        public class HeaderizedFile<outputType>
        {
            public delegate outputType Parse(Dictionary<string, int> fieldMappings, string[] fields);
            public delegate outputType FieldGrabberParser(FieldGrabber fieldGrabber);
            public delegate void ProcessNewItem(outputType newItem);

            public HeaderizedFile(StreamReader inputFile_, bool hasVersion_, bool hasDone_, string expectedVersion_, List<string> wantedFields_,
                bool skipHash_ = false, bool allowMissingColumnsInData_ = false, int skippedRows_ = 0, bool allowMissingRowsInData_ = false,
                char separator_ = '\t', bool stopAtBlankLine_ = false, string inputFilename_ = "", int rowsToSkipAfterHeader_ = 0, List<string> additionalValuesToTreatAsStar_ = null)
            {
                inputFile = inputFile_;
                hasVersion = hasVersion_;
                hasDone = hasDone_;
                expectedVersion = expectedVersion_;
                wantedFields = wantedFields_;
                skipHash = skipHash_;
                allowMissingColumnsInData = allowMissingColumnsInData_;
                skippedRows = skippedRows_;
                allowMissingRowsInData = allowMissingRowsInData_;
                separator = separator_;
                stopAtBlankLine = stopAtBlankLine_;
                inputFilename = inputFilename_;
                rowsToSkipAfterHeader = rowsToSkipAfterHeader_;
                if (additionalValuesToTreatAsStar_ != null)
                {
                    additionalValuesToTreatAsStar = additionalValuesToTreatAsStar_;
                }
            }

            public bool ParseFile(FieldGrabberParser parser, ProcessNewItem processNewItem)
            {
                Dictionary<string, int> fieldMappings;

                return ParseFile(null, parser, out fieldMappings, processNewItem);
            }

            //
            // Just eat the fieldMappings output the clunky way, since out parameters can't have default values.
            //
            public bool ParseFile(Parse parser, out List<outputType> result)
            {
                Dictionary<string, int> fieldMappings;

                return ParseFile(parser, out result, out fieldMappings);
            }

            public bool ParseFile(FieldGrabberParser parser, out List<outputType> result)
            {
                Dictionary<string, int> fieldMappings;

                return ParseFile(parser, out result, out fieldMappings);
            }

            public bool ParseFile(Parse parser, out List<outputType> result, out Dictionary<string, int> fieldMappings_out)
            {
                return ParseFile(parser, null, out result, out fieldMappings_out);
            }

            public bool ParseFile(FieldGrabberParser parser, out List<outputType> result, out Dictionary<string, int> fieldMappings_out)
            {
                return ParseFile(null, parser, out result, out fieldMappings_out);
            }

            bool ParseFile(Parse parser, FieldGrabberParser fieldGrabbingParser, out List<outputType> result, out Dictionary<string, int> fieldMappings_out)
            {
                result = new List<outputType>();
                var resultCopy = result;    // This is necessary because you can't use an out parameter in a lambda

                bool retVal = ParseFile(parser, fieldGrabbingParser, out fieldMappings_out, x => resultCopy.Add(x));

                if (!retVal)
                {
                    result = null;
                }

                return retVal;
            }

            bool ParseFile(Parse parser, FieldGrabberParser fieldGrabbingParser, out Dictionary<string, int> fieldMappings_out, ProcessNewItem processNewItem)
            {
                fieldMappings_out = null;

                if (wantedFields.Any(x => wantedFields.Where(y => x == y).Count() != 1))
                {
                    Console.Write("ASETools.HeaderizedFile.ParseFile: wantedFields contains the following duplicate columns (which will appear more than once in this list):");

                    foreach (var duplicateField in wantedFields.Where(x => wantedFields.Where(y => x == y).Count() != 1))
                    {
                        Console.Write(" " + duplicateField);
                    }
                    return false;
                }

                if (hasVersion)
                {
                    var versionString = inputFile.ReadLine();

                    if (versionString == null || expectedVersion != null && versionString != expectedVersion)
                    {
                        return false;
                    }
                }

                for (var i = 0; i < skippedRows; i++)
                {
                    inputFile.ReadLine();
                }

                var header = inputFile.ReadLine();

                if (skipHash)
                {
                    while (header != null && header.StartsWith("#"))
                    {
                        header = inputFile.ReadLine();
                    }
                }

                if (null == header)
                {
                    return false;
                }

                for (int i = 0; i < rowsToSkipAfterHeader; i++)
                {
                    inputFile.ReadLine();
                }

                var columns = header.Split(separator);
                var fieldMappings = new Dictionary<string, int>();
                int maxNeededField = -1;

                for (int i = 0; i < columns.Count(); i++)
                {
                    if (wantedFields.Contains(columns[i]))
                    {
                        if (fieldMappings.ContainsKey(columns[i]))
                        {
                            Console.WriteLine("Duplicate needed column in headerized file " + inputFilename + " (or code bug or something): " + columns[i]);
                            return false;
                        }

                        fieldMappings.Add(columns[i], i);
                        maxNeededField = i;
                    }
                }

                if (fieldMappings.Count() != wantedFields.Count())
                {
                    var missingColumns = new List<string>();
                    foreach (var wantedField in wantedFields)
                    {
                        if (!fieldMappings.ContainsKey(wantedField))
                        {
                            missingColumns.Add(wantedField);
                        }
                    }

                    if (fieldMappings.Count() + missingColumns.Count() != wantedFields.Count())
                    {
                        Console.WriteLine("Got the wrong number of missing fields.  Code bug.");
                        return false;

                    }
                    Console.Write("Headerized file (" + inputFilename + "): missing columns:");
                    foreach (var missingColumn in missingColumns)
                    {
                        Console.Write(" " + missingColumn);
                        fieldMappings.Add(missingColumn, maxNeededField + 1);
                    }
                    Console.WriteLine(".  Filling in with the empty string.");
                    hasMissingFields = true;
                }

                string inputLine;
                bool sawDone = false;
                while (null != (inputLine = inputFile.ReadLine())) {
                    if (sawDone) {
                        Console.WriteLine("HeaderizedFile (" + inputFilename + "): Saw data after **done**");
                        return false;
                    }

                    if ("**done**" == inputLine.Trim()) {
                        sawDone = true;
                        continue;
                    }

                    if (inputLine == "" && stopAtBlankLine)
                    {
                        break;
                    }

                    var fields = inputLine.Split('\t');
                    if (fields.Count() <= maxNeededField && !allowMissingColumnsInData)
                    {
                        Console.WriteLine("HeaderizedFile.Parse (" + inputFilename + "): input line didn't include a needed field " + inputLine);
                        return false;
                    }
                    else if ((hasMissingFields || allowMissingColumnsInData) && fields.Count() <= maxNeededField + 1)
                    {
                        var extendedFields = new string[maxNeededField + 2];
                        for (int i = 0; i <= maxNeededField; i++)
                        {
                            if (i < fields.Count())
                            {
                                extendedFields[i] = fields[i];
                            } else
                            {
                                extendedFields[i] = "";
                            }
                        }
                        fields = extendedFields;
                    }

                    if (hasMissingFields || allowMissingColumnsInData && fields.Count() <= maxNeededField + 1)
                    {
                        fields[maxNeededField + 1] = "";    // This is for all missing fields
                    }

                    if (null != parser)
                    {
                        processNewItem(parser(fieldMappings, fields));
                    } else
                    {
                        processNewItem(fieldGrabbingParser(new FieldGrabber(fieldMappings, fields, inputLine, additionalValuesToTreatAsStar)));
                    }
                }

                if (hasDone && !sawDone)
                {
                    Console.WriteLine("HeaderizedFile.Parse (" + inputFilename + "): missing **done**");
                    return false;
                }
                else if (!hasDone && sawDone)
                {
                    Console.WriteLine(inputFilename + "Saw unepected **done**.  Ignoring.");
                    return false;
                }

                fieldMappings_out = fieldMappings;
                return true;
            } // ParseFile

            public class FieldGrabber
            {
                public FieldGrabber(Dictionary<string, int> fieldMappings_, string[] fields_, string rawInputLine_, List<string> additionalValuesToTreatAsStar_)
                {
                    fieldMappings = fieldMappings_;
                    fields = fields_;
                    rawInputLine = rawInputLine_;
                    additionalValuesToTreatAsStar = additionalValuesToTreatAsStar_;
                }

                public string AsString(string fieldName)
                {
                    var x = ConvertToNonExcelString(fields[fieldMappings[fieldName]]);
                    return x;
                }

                public string AsStringEmptyForNA(string fieldName)
                {
                    var x = ConvertToNonExcelString(fields[fieldMappings[fieldName]]);

                    if (x == "NA")
                    {
                        return "";
                    }

                    return x;
                }
                public int AsInt(string fieldName)
                {
                    return Convert.ToInt32(AsString(fieldName));
                }

                public string AsStringFromSquareBracketString(string fieldName)
                {
                    var rawValue = AsString(fieldName);
                    if (rawValue == "")
                    {
                        return "";
                    }

                    if (rawValue.First() != '[' || rawValue.Last() != ']')
                    {
                        throw new FormatException("FieldGrabber.AsStringFromSquareBracketString: doesn't have square brackets: " + rawValue);
                    }

                    return rawValue.Substring(1, rawValue.Length - 2);
                }

                public string AsStringFromSquareBracketSingleQuoteString(string fieldName)
                {
                    var deBracketed = AsStringFromSquareBracketString(fieldName);

                    if (deBracketed == "")
                    {
                        return "";
                    }

                    if (deBracketed.First() != '\'' || deBracketed.Last() != '\'')
                    {
                        throw new FormatException("FieldGrabber.AsStringFromSquareBracketSingleQuoteString: missing single quotes: " + deBracketed);
                    }

                    return deBracketed.Substring(1, deBracketed.Length - 2);
                }

                public int AsIntFromSquareBracketInt(string fieldName)
                {
                    if (AsString(fieldName) == "") return 0;

                    return Convert.ToInt32(AsStringFromSquareBracketString(fieldName));
                }

                public long AsLongFromSquareBracketLong(string fieldName)
                {
                    if (AsString(fieldName) == "") return 0;

                    return Convert.ToInt64(AsStringFromSquareBracketString(fieldName));
                }

                public long AsLong(string fieldName)
                {
                    return Convert.ToInt64(AsString(fieldName));
                }

                public int AsIntMinusOneIfStarOrEmptyString(string fieldName)
                {
                    var value = fields[fieldMappings[fieldName]];
                    if (value == "*" || value == "" || additionalValuesToTreatAsStar.Contains(value))
                    {
                        return -1;
                    }

                    return AsInt(fieldName);
                }

                public double AsDouble(string fieldName)
                {
                    return Convert.ToDouble(AsString(fieldName));
                }

                public double AsDoubleNegativeInfinityIfStarOrEmptyString(string fieldName)
                {
                    var value = fields[fieldMappings[fieldName]];
                    if (value == "*" || value == "" || additionalValuesToTreatAsStar.Contains(value))
                    {
                        return double.NegativeInfinity;
                    }

                    return AsDouble(fieldName);
                }

                public bool AsBool(string fieldName)
                {
                    return Convert.ToBoolean(AsString(fieldName));
                }

                public string rawLine()
                {
                    var result = fields[0];
                    for (int i = 1; i < fields.Count(); i++)
                    {
                        result += "\t" + fields[i];
                    }

                    return result;
                }

                public bool YOrNToBool(string fieldName)
                {
                    if (fields[fieldMappings[fieldName]] == "y" || fields[fieldMappings[fieldName]] == "Y")
                    {
                        return true;
                    }

                    if (fields[fieldMappings[fieldName]] == "n" || fields[fieldMappings[fieldName]] == "N")
                    {
                        return false;
                    }

                    throw new FormatException();
                }

                Dictionary<string, int> fieldMappings;
                string[] fields;
                readonly string rawInputLine;
                List<string> additionalValuesToTreatAsStar;
            } // FieldGrabber

            StreamReader inputFile;
            bool hasVersion;
            bool hasDone;
            string expectedVersion;
            List<string> wantedFields;
            bool hasMissingFields = false;
            bool skipHash;
            bool allowMissingColumnsInData;
            bool allowMissingRowsInData;
            int skippedRows;
            char separator;
            bool stopAtBlankLine;
            string inputFilename;
            int rowsToSkipAfterHeader;
            List<string> additionalValuesToTreatAsStar = new List<string>();
        } // HeaderizedFile

        public static int ConvertToInt32TreatingNullStringAsZero(string value)
        {
            if (value == "") return 0;
            return Convert.ToInt32(value);
        }

        // Gives information on CompositeREFs from methylation and corresponding gene information
        public class CompositeREFInfoLine
        {
            public string Composite_Element_REF;
            public string Chromosome;
            public int Position;
            public string Hugo_Symbol;

            static KeyValuePair<string, GeneLocationInfo> ParseLine(Dictionary<string, int> fieldMappings, string[] fields)
            {
                var composite = fields[fieldMappings["compositeREF"]];
                var geneInfo = new GeneLocationInfo();
                geneInfo.chromosome = fields[fieldMappings["chr"]];
                geneInfo.minLocus = Convert.ToInt32(fields[fieldMappings["start"]]);
                geneInfo.maxLocus = Convert.ToInt32(fields[fieldMappings["end"]]);

                return new KeyValuePair<string, GeneLocationInfo>(composite, geneInfo);

            }


            static public List<KeyValuePair<string, GeneLocationInfo>> ReadFile(string filename)
            {
                StreamReader inputFile = CreateStreamReaderWithRetry(filename);

                var neededFields = new List<string>();
                neededFields.Add("compositeREF");
                neededFields.Add("chr");
                neededFields.Add("start");
                neededFields.Add("end");

                bool hasDone = false;
                var headerizedFile = new HeaderizedFile<KeyValuePair<string, GeneLocationInfo>>(inputFile, false, hasDone, "#version gdc-1.0.0", neededFields);

                List<KeyValuePair<string, GeneLocationInfo>> result;

                if (!headerizedFile.ParseFile((a, b) => ParseLine(a, b), out result))
                {
                    Console.WriteLine("Error reading Composite REF File " + filename);
                    return null;
                }

                inputFile.Close();
                return result;
            } // ReadFile

        } // CompositeREFInfoFile


        public class CompositeREF
        {
            public string Composite_Element_REF;
            public string Chromosome;
            public int Start;
            public int End;
            public string[] Gene_Symbol;
            public string[] Gene_Type;
            public string[] Transcript_ID;
            public int[] Position_to_TSS;
            public string CGI_Coordinate;
            public FeatureType Feature_Type;

            // Methylation types
            public enum FeatureType
            {
                // Feature type not specified
                Other,
                // methylation found in 2kb region upstream from island
                N_Shore,
                // methylation found in 2kb region downstream from island
                S_Shore,
                // CpG island
                Island,
            }

            public CompositeREF(string Composite_Element_REF_,
                string Chromosome_,
                int Start_,
                int End_,
                string[] Gene_Symbol_,
                string[] Gene_Type_,
                string[] Transcript_ID_,
                int[] Position_to_TSS_,
                string CGI_Coordinate_,
                FeatureType Feature_Type_)
            {

                var geneCount = Gene_Symbol_.Length;

                if (Gene_Type_.Length != geneCount || Transcript_ID_.Length != geneCount || Position_to_TSS_.Length != geneCount)
                {
                    Console.WriteLine("Error: Gene Count does not match supporting gene data");
                }

                Composite_Element_REF = Composite_Element_REF_;
                Chromosome = Chromosome_;
                Start = Start_;
                End = End_;
                Gene_Symbol = Gene_Symbol_;
                Gene_Type = Gene_Type_;
                Transcript_ID = Transcript_ID_;
                Position_to_TSS = Position_to_TSS_;
                CGI_Coordinate = CGI_Coordinate_;
                Feature_Type = Feature_Type_;

            }
        }

        // contains annotation data from one line in a methylation file
        public class MethylationAnnotationLine
        {
            public static List<MethylationAnnotationLine> filterByTSS(List<MethylationAnnotationLine> sites, int maxDistance = 2000)
            {
                var filteredSites = sites.Where(r =>
                {
                    var geneAbsAverages =
                    r.compositeREF.Gene_Symbol.Zip(r.compositeREF.Position_to_TSS, (first, second) => new Tuple<string, int>(first, second))
                        .GroupBy(i => i.Item1)  // group by Gene symbol
                        .Select(i => Math.Abs(i.Select(k => k.Item2).Average()));   // get average distance from TSS for each gene

                    return geneAbsAverages.Where(i => i <= maxDistance).Count() > 0;
                });

                return filteredSites.ToList();
            }

            public static double hemimethylation(double value)
            {
                return Math.Abs(2 * value - 1);
            }

            public static double M2Beta(double M_Value)
            {
                var expM = Math.Pow(2, M_Value);
                return expM / (expM + 1);
            }

            public static double betaToM(double beta) {

                // correct from 15 digit precision NaNs
                var epsilon = 0.000000000000001;
                if (beta == 1)
                    beta -= epsilon;
                if (beta == 0)
                    beta += epsilon;

                // M_value calculated as shown in Du et. al. 2010
                return Math.Log(beta / (1 - beta));
            }

            public double Beta_Value;
            public CompositeREF compositeREF;


            // this is not in the original file, but calculated from beta_value
            public readonly double M_Value;

            MethylationAnnotationLine(CompositeREF compositeREF_, double Beta_Value_)
            {
                //var geneCount = Gene_Symbol_.Length;

                //if (Gene_Type_.Length != geneCount || Transcript_ID_.Length != geneCount || Position_to_TSS_.Length != geneCount)
                //{
                //	Console.WriteLine("Error: Gene Count does not match supporting gene data");
                //}

                compositeREF = compositeREF_;
                Beta_Value = Beta_Value_;

                // M_value calculated as shown in Du et. al. 2010
                M_Value = betaToM(Beta_Value);
            }

            // Converts strings to doubles, returning 1 if empty string.
            // We return 1 here because this function parses beta values, which must be between 0 and 1
            static Double ConvertToDoubleTreatingNullStringAsOne(string value)
            {
                Double n;
                if (!Double.TryParse(value, out n))
                    n = 1;
                return n;
            }

            // splits string to array by specified delimiter. This is used mainly when gene symbols are left blank,
            // which is specified by a '.'
            static string[] splitPotentialString(string value, char delim)
            {
                string[] parsed = new string[] { };
                if (value != ".")
                    parsed = value.Split(delim);
                return parsed;
            }

            static MethylationAnnotationLine ParseLine(Dictionary<string, int> fieldMappings, string[] fields)
            {
                // delimiter separating columns
                var delim = ';';

                Enum.TryParse(fields[fieldMappings["Feature_Type"]], out CompositeREF.FeatureType Feature_Type);

                var compositeREF = new CompositeREF(
                fields[fieldMappings["Composite Element REF"]],
                fields[fieldMappings["Chromosome"]],
                Convert.ToInt32(fields[fieldMappings["Start"]]),
                Convert.ToInt32(fields[fieldMappings["End"]]),
                splitPotentialString(fields[fieldMappings["Gene_Symbol"]], delim),
                splitPotentialString(fields[fieldMappings["Gene_Type"]], delim),
                splitPotentialString(fields[fieldMappings["Transcript_ID"]], delim),
                splitPotentialString(fields[fieldMappings["Position_to_TSS"]], delim).Select(s => Convert.ToInt32(s)).ToArray(),
                fields[fieldMappings["CGI_Coordinate"]],
                Feature_Type);

                return new MethylationAnnotationLine(compositeREF,
                ConvertToDoubleTreatingNullStringAsOne(fields[fieldMappings["Beta_value"]]));
            }

            static public List<MethylationAnnotationLine> ReadFile(string filename)
            {
                StreamReader inputFile;

                inputFile = CreateStreamReaderWithRetry(filename);
                if (null == inputFile)
                {
                    Console.WriteLine("Unable to open " + filename);
                    return null;
                }

                var neededFields = new List<string>();
                neededFields.Add("Composite Element REF");
                neededFields.Add("Beta_value");
                neededFields.Add("Chromosome");
                neededFields.Add("Start");
                neededFields.Add("End");
                neededFields.Add("Gene_Symbol");
                neededFields.Add("Gene_Type");
                neededFields.Add("Transcript_ID");
                neededFields.Add("Position_to_TSS");
                neededFields.Add("CGI_Coordinate");
                neededFields.Add("Feature_Type");

                bool hasDone = false;
                var headerizedFile = new HeaderizedFile<MethylationAnnotationLine>(inputFile, false, hasDone, "#version gdc-1.0.0", neededFields);

                List<MethylationAnnotationLine> result;

                if (!headerizedFile.ParseFile((a, b) => ParseLine(a, b), out result))
                {
                    Console.WriteLine("Error reading Annotation File " + filename);
                    return null;
                }

                inputFile.Close();

                // filter out invalid results. Invalid results include records without valid M values and records without gene symbols.
                result = result.Where(c => (c.M_Value != double.PositiveInfinity) && (c.compositeREF.Gene_Symbol.Length > 0)).ToList();

                return result;
            } // ReadFile
        } // MethylationAnnotationLine

        public class MethylationDistributionLine
        {
            public string contig;
            public int locus;
            public int n;
            public double min;
            public double max;
            public double mean;
            public double stdDev;

            MethylationDistributionLine(string contig_, int locus_, int n_, double min_, double max_, double mean_, double stdDev_)
            {
                contig = contig_;
                locus = locus_;
                n = n_;
                min = min_;
                max = max_;
                mean = mean_;
                stdDev = stdDev_;
            }

            static public List<MethylationDistributionLine> ReadFromFile(string filename)
            {

                var inputFile = CreateStreamReaderWithRetry(filename);
                if (null == inputFile)
                {
                    Console.WriteLine("MethylationDistributionLine.ReadFromFile: unable to open input file " + filename);
                    return null;
                }

                string[] wantedFields = { "contig", "locus", "n", "min beta", "max beta", "mean beta", "beta standard deviation" };

                var headerizedFile = new HeaderizedFile<MethylationDistributionLine>(inputFile, false, true, "", wantedFields.ToList());
                List<MethylationDistributionLine> retVal;

                if (!headerizedFile.ParseFile(Parse, out retVal))
                {
                    Console.WriteLine("Unable to parse " + filename);
                    inputFile.Close();
                    return null;
                }

                inputFile.Close();

                return retVal;
            }

            static MethylationDistributionLine Parse(HeaderizedFile<MethylationDistributionLine>.FieldGrabber fieldGrbber)
            {
                return new MethylationDistributionLine(fieldGrbber.AsString("contig"), fieldGrbber.AsInt("locus"), fieldGrbber.AsInt("n"), fieldGrbber.AsDouble("min beta"), fieldGrbber.AsDouble("max beta"), fieldGrbber.AsDouble("mean beta"),
                    fieldGrbber.AsDouble("beta standard deviation"));
            }
        }
        public class MethylationDistribution
        {

            public static MethylationDistribution ReadFromFile(string filename)
            {
                var lines = MethylationDistributionLine.ReadFromFile(filename);
                if (null == lines)
                {
                    Console.WriteLine("Failed to read methylation lines from " + filename);
                    return null;
                }

                var retVal = new MethylationDistribution();

                foreach (var line in lines)
                {
                    if (!retVal.lineMap.ContainsKey(line.contig))
                    {
                        retVal.lineMap.Add(line.contig, new Dictionary<int, MethylationDistributionLine>());
                    }

                    retVal.lineMap[line.contig].Add(line.locus, line);
                }

                return retVal;
            }

            public MethylationDistributionLine forLocus(string contig, int locus)
            {
                if (!lineMap.ContainsKey(contig) || !lineMap[contig].ContainsKey(locus))
                {
                    return null;
                }

                return lineMap[contig][locus];
            }

            Dictionary<string, Dictionary<int, MethylationDistributionLine>> lineMap = new Dictionary<string, Dictionary<int, MethylationDistributionLine>>();
        } // MethylationDistribution



        public class BedLine
        {
            public readonly string Chromosome;
            public readonly int Start_Position;
            public readonly int End_Position;
            public readonly string name;
            public readonly int score; // value 0-1000. for methylation score/1000 gives you methylation percentage
            public readonly char strand; // can be '.', '+' or '-'
            public readonly double thickStart;
            public readonly double thickEnd;

            public BedLine(string Chromosome_,
                int Start_Position_,
                int End_Position_,
                string name_,
                int score_,
                char strand_,
                double thickStart_ = 0,
                double thickEnd_ = 1)
            {
                Chromosome = Chromosome_;
                Start_Position = Start_Position_;
                End_Position = End_Position_;
                name = name_;
                score = score_;
                strand = strand_;
                thickStart = thickStart_;
                thickEnd = thickEnd_;
            }

            public bool overlaps(string otherChromosome, int otherStart, int otherEnd)
            {
                if (ASETools.chromosomeNameToNonChrForm(otherChromosome) == ASETools.chromosomeNameToNonChrForm(Chromosome))
                {
                    if (otherStart < End_Position && otherEnd > Start_Position)
                    {
                        return true;
                    }
                }
                return false;
            }

            public string toString() {
                return this.Chromosome + "\t" + this.Start_Position + "\t" + this.End_Position + "\t" + this.name + "\t" + this.score +
                    "\t" + this.strand + "\t" + this.thickStart + "\t" + this.thickEnd;
            }

            public static void writeFile(string filename, List<BedLine> bedLines, string header = "")
            {
                var file = CreateStreamWriterWithRetry(filename);

                if (null == file)
                {
                    Console.WriteLine("BedLine.writeFile: unable to create file " + filename);
                    return;
                }

                // Write header
                file.WriteLine(header);

                foreach (var line in bedLines)
                {
                    file.WriteLine(line.toString());
                }
                file.Close();
            }

            static public List<BedLine> ReadFile(string filename, bool hasHeader = true)
            {
                StreamReader inputFile = CreateStreamReaderWithRetry(filename);

                if (hasHeader)
                {
                    var header = inputFile.ReadLine();
                }

                List<BedLine> result = new List<BedLine>();
                string line;
                while ((line = inputFile.ReadLine()) != null)
                {
                    var split = line.Split('\t');
                    result.Add(new BedLine(split[0], Convert.ToInt32(split[1]), Convert.ToInt32(split[2]), split[3], Convert.ToInt32(split[4]), split[5].ToCharArray()[0], Convert.ToDouble(split[6]), Convert.ToDouble(split[7])));
                }

                inputFile.Close();

                return result;
            }
        }

        public class MAFLine : IComparable
        {
            public readonly string Hugo_Symbol;
            public readonly string NCBI_Build;
            public readonly string Chromosome;
            public int Start_Position;	// these are sometimes reset when switching builds
            public int End_Positon;
            public readonly string Variant_Classification;
            public readonly string Variant_Type;
            public readonly string Reference_Allele;
            public readonly string Tumor_Seq_Allele1;
            public readonly string Tumor_Seq_Allele2;
            public readonly string Match_Norm_Seq_Allele1;
            public readonly string Match_Norm_Seq_Allele2;
            public readonly string Tumor_Sample_UUID;
            public readonly string Matched_Norm_Sample_UUID;
            public readonly string tumor_bam_uuid;
            public readonly string normal_bam_uuid;
            public readonly int t_depth;
            public readonly int t_ref_count;
            public readonly int t_alt_count;
            public readonly int n_depth;
            public readonly int n_ref_count;
            public readonly int n_alt_count;

            public readonly string maf_file_id; // Of the MAF file

            public MAFLine(string Chromosome_, int Start_Position_) // A constructor for making keys for comparison in AVL Trees
            {
                Chromosome = Chromosome_;
                Start_Position = Start_Position_;
            }

            public int CompareTo(object peerObject)
            {
                MAFLine peer = (MAFLine)peerObject;

                if (Chromosome != peer.Chromosome)
                {
                    return Chromosome.CompareTo(peer.Chromosome);
                }

                return Start_Position.CompareTo(peer.Start_Position);
            }

            public bool IsNonsenseMediatedDecayCausing()
            {
                return NonsenseMediatedDecayCausingVariantClassifications.Contains(Variant_Classification);
            }

            public bool IsSilent()
            {
                return Variant_Classification == "Silent";
            }

            public string IDH1MutationType()
            {
                if (Hugo_Symbol.ToUpper() != "IDH1")
                {
                    throw new Exception("MAFLine.IHD1MutationType: this is not an IDH1 mutation, hugo symbol = " + Hugo_Symbol);
                }

                if (Variant_Type != "SNP")
                {
                    return "Other";
                }

                return IDH1MutantDescription(Start_Position, Tumor_Seq_Allele2[0]);
            }

            MAFLine(string Hugo_Symbol_,
             string NCBI_Build_,
             string Chromosome_,
             int Start_Position_,
             int End_Positon_,
             string Variant_Classification_,
             string Variant_Type_,
             string Reference_Allele_,
             string Tumor_Seq_Allele1_,
             string Tumor_Seq_Allele2_,
             string Match_Norm_Seq_Allele1_,
             string Match_Norm_Seq_Allele2_,
             string Tumor_Sample_UUID_,
             string Matched_Norm_Sample_UUID_,
             string tumor_bam_uuid_,
             string normal_bam_uuid_,
             int t_depth_,
             int t_ref_count_,
             int t_alt_count_,
             int n_depth_,
             int n_ref_count_,
             int n_alt_count_,
             string maf_file_id_)
            {
                Hugo_Symbol = ConvertToNonExcelString(Hugo_Symbol_);
                NCBI_Build = NCBI_Build_;
                Chromosome = Chromosome_;
                Start_Position = Start_Position_;
                End_Positon = End_Positon_;
                Variant_Classification = Variant_Classification_;
                Variant_Type = Variant_Type_;
                Reference_Allele = Reference_Allele_;
                Tumor_Seq_Allele1 = Tumor_Seq_Allele1_;
                Tumor_Seq_Allele2 = Tumor_Seq_Allele2_;
                Match_Norm_Seq_Allele1 = Match_Norm_Seq_Allele1_;
                Match_Norm_Seq_Allele2 = Match_Norm_Seq_Allele2_;
                Tumor_Sample_UUID = Tumor_Sample_UUID_;
                Matched_Norm_Sample_UUID = Matched_Norm_Sample_UUID_;
                tumor_bam_uuid = tumor_bam_uuid_;
                normal_bam_uuid = normal_bam_uuid_;
                t_depth = t_depth_;
                t_ref_count = t_ref_count_;
                t_alt_count = t_alt_count_;
                n_depth = n_depth_;
                n_ref_count = n_ref_count_;
                n_alt_count = n_alt_count_;
                maf_file_id = maf_file_id_;
            }

            static MAFLine ParseLine(Dictionary<string, int> fieldMappings, string[] fields, string maf_file_id, bool isBeatAML)
            {
                return new MAFLine(
                    fields[fieldMappings["Hugo_Symbol"]],
                    fields[fieldMappings["NCBI_Build"]],
                    fields[fieldMappings["Chromosome"]],
                    Convert.ToInt32(fields[fieldMappings["Start_Position"]]),
                    Convert.ToInt32(fields[fieldMappings["End_Position"]]),
                    fields[fieldMappings["Variant_Classification"]],
                    fields[fieldMappings["Variant_Type"]],
                    fields[fieldMappings["Reference_Allele"]],
                    fields[fieldMappings["Tumor_Seq_Allele1"]],
                    fields[fieldMappings["Tumor_Seq_Allele2"]],
                    fields[fieldMappings["Match_Norm_Seq_Allele1"]],
                    fields[fieldMappings["Match_Norm_Seq_Allele2"]],
                    fields[fieldMappings["Tumor_Sample_UUID"]],
                    fields[fieldMappings["Matched_Norm_Sample_UUID"]],
                    (isBeatAML ? "" : fields[fieldMappings["tumor_bam_uuid"]]),
                    (isBeatAML ? "" : fields[fieldMappings["normal_bam_uuid"]]),
                    ConvertToInt32TreatingNullStringAsZero(fields[fieldMappings["t_depth"]]),
                    ConvertToInt32TreatingNullStringAsZero(fields[fieldMappings["t_ref_count"]]),
                    ConvertToInt32TreatingNullStringAsZero(fields[fieldMappings["t_alt_count"]]),
                    ConvertToInt32TreatingNullStringAsZero(fields[fieldMappings["n_depth"]]),
                    ConvertToInt32TreatingNullStringAsZero(fields[fieldMappings["n_ref_count"]]),
                    ConvertToInt32TreatingNullStringAsZero(fields[fieldMappings["n_alt_count"]]),
                    maf_file_id
                    );
            }

            static public Dictionary<string, int> GetMutationCounts(List<MAFLine> mafLines)
            {
                // group MAFLines by hugo symbol to count mutations over symbol
                var mutationCounts = mafLines.GroupBy(r => r.Hugo_Symbol)
                    .Select(group => new
                    {
                        Hugo_Symbol = group.Key,
                        Count = group.Count()
                    });

                return mutationCounts.ToDictionary(x => x.Hugo_Symbol, x => x.Count);
            }

            static public List<MAFLine> ReadFile(string filename, string file_id, bool fileHasVersion, bool isBeatAML = false)
            {
                StreamReader inputFile;

                if (filename.Count() > 2 && filename.Substring(filename.Count() - 3, 3) == ".gz")
                {
                    inputFile = CreateCompressedStreamReaderWithRetry(filename);
                }
                else
                {
                    inputFile = CreateStreamReaderWithRetry(filename);
                }

                var neededFields = new List<string>();
                neededFields.Add("Hugo_Symbol");
                neededFields.Add("NCBI_Build");
                neededFields.Add("Chromosome");
                neededFields.Add("Start_Position");
                neededFields.Add("End_Position");
                neededFields.Add("Variant_Classification");
                neededFields.Add("Variant_Type");
                neededFields.Add("Reference_Allele");
                neededFields.Add("Tumor_Seq_Allele1");
                neededFields.Add("Tumor_Seq_Allele2");
                neededFields.Add("Match_Norm_Seq_Allele1");
                neededFields.Add("Match_Norm_Seq_Allele2");
                neededFields.Add("Tumor_Sample_UUID");
                neededFields.Add("Matched_Norm_Sample_UUID");
                if (!isBeatAML)
                {
                    neededFields.Add("tumor_bam_uuid");
                    neededFields.Add("normal_bam_uuid");
                }
                neededFields.Add("t_depth");
                neededFields.Add("t_ref_count");
                neededFields.Add("t_alt_count");
                neededFields.Add("n_depth");
                neededFields.Add("n_ref_count");
                neededFields.Add("n_alt_count");


                var headerizedFile = new HeaderizedFile<MAFLine>(inputFile, fileHasVersion, !fileHasVersion, (isBeatAML ? "#version 2.4" : "#version gdc-1.0.0"), neededFields, fileHasVersion);

                List<MAFLine> result;

                if (!headerizedFile.ParseFile((a, b) => ParseLine(a, b, file_id, isBeatAML), out result))
                {
                    Console.WriteLine("Error reading MAF File " + filename);
                    return null;
                }

                inputFile.Close();

                return result;
            } // ReadFile

            static public List<MAFLine> ReadFileOnlyOnePerLocus(string filename, string file_id, bool fileHasVersion)
            {
                var rawList = ReadFile(filename, file_id, fileHasVersion);
                rawList.Sort();

                int n = rawList.Count();

                var cookedList = new List<MAFLine>();

                bool inRun = false;

                for (int i = 0; i < n; i++)
                {
                    if (inRun && rawList[i - 1].Chromosome == rawList[i].Chromosome && rawList[i - 1].Start_Position == rawList[i].Start_Position)
                    {
                        continue;
                    }

                    inRun = false;

                    int endPosition = rawList[i].End_Positon;
                    int biggestInRun = i;

                    for (int j = i + 1; j < n && rawList[j].Chromosome == rawList[i].Chromosome && rawList[j].Start_Position == rawList[i].Start_Position; j++)
                    {
                        inRun = true;
                        if (rawList[j].End_Positon > endPosition)
                        {
                            endPosition = rawList[j].End_Positon;
                            biggestInRun = j;
                        }
                    }

                    cookedList.Add(rawList[biggestInRun]);
                } // for all in raw list

                return cookedList;
            }

            public static void WriteHeaderLine(StreamWriter outputStream)
            {
                outputStream.WriteLine("Hugo_Symbol\tNCBI_Build\tChromosome\tStart_Position\tEnd_Position\tVariant_Classification\tVariant_Type\tReference_Allele\tTumor_Seq_Allele1\tTumor_Seq_Allele2\tMatch_Norm_Seq_Allele1\tMatch_Norm_Seq_Allele2\tTumor_Sample_UUID\tMatched_Norm_Sample_UUID\ttumor_bam_uuid\tnormal_bam_uuid\tt_depth\tt_ref_count\tt_alt_count\tn_depth\tn_ref_count\tn_alt_count");
            }

            public void WriteToStream(StreamWriter output)
            {
                output.WriteLine(
                    ConvertToExcelString(Hugo_Symbol) + "\t" +
                    NCBI_Build + "\t" +
                    Chromosome + "\t" +
                    Start_Position + "\t" +
                    End_Positon + "\t" +
                    Variant_Classification + "\t" +
                    Variant_Type + "\t" +
                    Reference_Allele + "\t" +
                    Tumor_Seq_Allele1 + "\t" +
                    Tumor_Seq_Allele2 + "\t" +
                    Match_Norm_Seq_Allele1 + "\t" +
                    Match_Norm_Seq_Allele2 + "\t" +
                    Tumor_Sample_UUID + "\t" +
                    Matched_Norm_Sample_UUID + "\t" +
                    tumor_bam_uuid + "\t" +
                    normal_bam_uuid + "\t" +
                    t_depth + "\t" +
                    t_ref_count + "\t" +
                    t_alt_count + "\t" +
                    n_depth + "\t" +
                    n_ref_count + "\t" +
                    n_alt_count
                    );
            }

            public static void WriteToFile(string filename, List<MAFLine> linesToWrite)
            {
                var output = CreateStreamWriterWithRetry(filename);

                WriteHeaderLine(output);

                foreach (var mafLine in linesToWrite)
                {
                    mafLine.WriteToStream(output);
                }

                output.WriteLine("**done**");
                output.Close();
            }

            public string getExtractedReadsExtension()
            {
                return "-" + Chromosome + "-" + Math.Max(1, Start_Position - 200) + "-" + (End_Positon + 10);
            }

        } // MAFLine


        public static string SizeToUnits(ulong size)
        {
            if (size < 1024) return "" + size;

            ulong divisor;

            string suffix;

            if (size < 1024 * 1024)
            {
                divisor = 1024;
                suffix = "K";
            }
            else if (size < 1024 * 1024 * 1024)
            {
                divisor = 1024 * 1024;
                suffix = "M";
            }
            else if (size < (ulong)1024 * 1024 * 1024 * 1024)
            {
                divisor = 1024 * 1024 * 1024;
                suffix = "G";
            }
            else if (size < (ulong)1024 * 1024 * 1024 * 1024 * 1024)
            {
                divisor = (ulong)1024 * 1024 * 1024 * 1024;
                suffix = "T";
            } else
            {
                divisor = (ulong)1024 * 1024 * 1024 * 1024 * 1024;
                suffix = "P";
            }

            if (size / divisor > 9)
            {
                return "" + (size + divisor / 2) / divisor + suffix;
            }

            return "" + size / divisor + "." + ((size * 10 + divisor / 20) / divisor) % 10 + suffix;
        }

        public class AllcountReader
        {
            public AllcountReader(string filename_)
            {
                filename = filename_;
            }

            public bool openFile()
            {
                long mappedHQNuclearReads;
                int numContigs;

                return openFile(out mappedHQNuclearReads, out numContigs);
            }

            public bool openFile(out long mappedHQNuclearReads, out int out_numContigs)
            {
                if (null != allcountReader)
                {
                    Console.WriteLine("ASETools.AllcountReader.openFile(): file is already open.  You can only do this once per allcount file.");
                }

                allcountReader = null;
                out_numContigs = numContigs = 0;
                contigs = null;
                mappedHQNuclearReads = 0;

                try
                {
                    compressedStreamReader = CreateStreamReaderWithRetry(filename);

                    if (compressedStreamReader == null)
                    {
                        Console.WriteLine("AllcountReader.openFile: Failed to open compressed stream reader.");
                        return false;
                    }
                    allcountReader = new StreamReader(new GZipStream(compressedStreamReader.BaseStream, CompressionMode.Decompress));


                    //
                    // The format is two header lines like:
                    //      CountReadsCovering v1.1 \\msr-genomics-0\d$\sequence\indices\grch37-24 -a \\msr-genomics-0\e$\tcga\hnsc\tumor\rna\0415c9cc-48e6-46e7-9076-06b6b56bb4be\PCAWG.e7697e88-96df-4c96-ad8d-6c54df95d29b.STAR.v1.bam -
                    //      83303682 mapped high quality nuclear reads, 9891899 low quality reads, 0 mitochondrial reads, 100740272 total reads
                    // followed by a blank line followed by
                    //      NumContigs: 93
                    //      ContigName\tLength
                    //      <ContigName\tLength> x numContigs
                    // Then for each contig:
                    //      >ContigName
                    // And within each contig one of three line types:
                    //      offset\tcount
                    // Which is an offset in the contig (in HEX) followed by the count of high quality reads that map there (also in hex) or:
                    //      xnumber
                    // which is the number of consecutive loci with the same number of HQ reads mapped as the previous locus (x is literal, count is in hex) or:
                    //      count
                    // which is the new count of reads mapped for the next locus (in hex)
                    //
                    // And the final line in the file is
                    //      **done**
                    // 
                    //

                    var line = allcountReader.ReadLine();

                    const string headerBeginning = "CountReadsCovering v";

                    if (null == line || line.Count() < headerBeginning.Count() + 1 || line.Substring(0, headerBeginning.Count()) != headerBeginning)
                    {
                        Console.WriteLine("Empty or corrupt allcount file " + filename);
                        return false;
                    }

                    if (line[headerBeginning.Count()] != '1')
                    {
                        Console.WriteLine("Unsupported major version of allcount file " + filename + ".  Header line: " + line);
                        return false;
                    }

                    line = allcountReader.ReadLine();
                    if (null == line)
                    {
                        Console.WriteLine("Corrupt or tuncated allcount file " + filename);
                        return false;
                    }

                    var fields = line.Split(' ');
                    if (fields.Count() != 16)
                    {
                        Console.WriteLine("Corrupt or tuncated allcount file " + filename + ".  Second line has " + fields.Count() + " fields: " + line);
                        return false;
                    }

                    try
                    {
                        mappedHQNuclearReads = Convert.ToInt64(fields[0]);
                    }
                    catch (FormatException)
                    {
                        Console.WriteLine("Format exception parsing mapped HQ read count for file " + filename + " from line: " + line);
                        return false;
                    }

                    line = allcountReader.ReadLine();   // The blank line

                    line = allcountReader.ReadLine();
                    if (null == line)
                    {
                        Console.WriteLine("Allcount file truncated before contig count: " + filename);
                        return false;
                    }

                    const string numContigsLineBeginning = "NumContigs: ";
                    if (line.Count() < numContigsLineBeginning.Count() + 1 || line.Substring(0, numContigsLineBeginning.Count()) != numContigsLineBeginning)
                    {
                        Console.WriteLine("Malformed NumContigs line in " + filename + ": " + line);
                        return false;
                    }

                    try
                    {
                        out_numContigs = numContigs = Convert.ToInt32(line.Substring(numContigsLineBeginning.Count()));
                    }
                    catch (FormatException)
                    {
                        Console.WriteLine("Couldn't parse NumContigs line in file " + filename + ": " + line);
                        return false;
                    }

                    if (numContigs < 1)
                    {
                        Console.WriteLine("Invalid numContigs in " + filename + ": " + line);
                        return false;
                    }

                    line = allcountReader.ReadLine();   // The header line for the contigs.

                    contigs = new Contig[numContigs];

                    int whichContig;

                    for (whichContig = 0; whichContig < numContigs; whichContig++)
                    {
                        contigs[whichContig] = new Contig();

                        line = allcountReader.ReadLine();
                        if (null == line)
                        {
                            Console.WriteLine("File truncated in contig list " + filename);
                            return false;
                        }

                        fields = line.Split('\t');
                        if (fields.Count() != 2)
                        {
                            Console.WriteLine("Incorrect contig line format in file " + filename + ": " + line);
                            return false;
                        }

                        contigs[whichContig].name = fields[0].ToLower();
                        try
                        {
                            contigs[whichContig].length = Convert.ToInt64(fields[1]);
                        }
                        catch (FormatException)
                        {
                            Console.WriteLine("Incorrect contig line format in file " + filename + ": " + line);
                            return false;
                        }
                    } // for all expected contigs
                } // try
                catch (Exception e)
                {
                    if (e is IOException || e is InvalidDataException)
                    {
                        Console.WriteLine("IOException or InvalidDataException opening allcount file " + filename);
                        return false;
                    }
                    else
                    {
                        throw e;
                    }
                }

                return true;
            } // openFile

            public delegate void ProcessBase(string contigName, int location, int currentMappedReadCount);

            public bool ReadAllcountFile(ProcessBase processBase, string onlyThisContig = null)
            {
                int currentOffset = -1;
                int currentMappedReadCount = -1;
                int whichContig = -1;

                bool sawDone = false;
                string contigName = "";

                bool processThisContig = onlyThisContig == null;

                string line;

                while (null != (line = allcountReader.ReadLine()))
                {
                    if (sawDone)
                    {
                        Console.WriteLine("File " + filename + " continues after **done** line: " + line);
                        return false;
                    }

                    if ("**done**" == line)
                    {
                        sawDone = true;
                        continue;
                    }

                    if (line.Count() == 0)
                    {
                        Console.WriteLine("Unexpected blank line in " + filename);
                        return false;
                    }

                    if (line[0] == '>')
                    {
                        whichContig++;
                        if (whichContig >= numContigs)
                        {
                            Console.WriteLine("Saw too many contigs in " + filename + ": " + line);
                            return false;
                        }

                        if (line.Substring(1).ToLower() != contigs[whichContig].name)
                        {
                            Console.WriteLine("Unexpected contig in " + filename + ".  Expected " + contigs[whichContig].name + ", got ", line.Substring(1));
                            return false;
                        }

                        if (onlyThisContig != null && chromosomeNameToNonChrForm(contigName) == chromosomeNameToNonChrForm(onlyThisContig))
                        {
                            //
                            // We've finished the one contig we're supposed to read.  Just quit now.
                            //
                            return true;
                        }

                        contigName = line.Substring(1).ToLower();

                        processThisContig = onlyThisContig == null || chromosomeNameToNonChrForm(onlyThisContig) == chromosomeNameToNonChrForm(contigName);

                        currentOffset = -1;
                        currentMappedReadCount = -1;
                        continue;
                    }

                    if (-1 == whichContig)
                    {
                        Console.WriteLine("Expected contig line after list of contigs, got " + line);
                        return false;
                    }

                    if (!processThisContig)
                    {
                        continue;
                    }

                    var fields = line.Split('\t');
                    if (fields.Count() == 1)
                    {
                        //
                        // Either xRepeatCount or newCount
                        //
                        if (line[0] == 'x')
                        {
                            if (currentMappedReadCount < 0 || currentOffset < 0)
                            {
                                Console.WriteLine("Got unexpected x line " + line);
                                return false;
                            }

                            int repeatCount;
                            try
                            {
                                repeatCount = Convert.ToInt32(line.Substring(1), 16); // 16 means the string is in hex
                                if (repeatCount <= 1)
                                {
                                    Console.WriteLine("Bogus repeat count " + line);
                                    return false;
                                }
                            }
                            catch (SystemException ex)
                            {
                                if (ex is FormatException || ex is ArgumentOutOfRangeException)
                                {
                                    Console.WriteLine("Format exception processing x line " + line);
                                    return false;
                                }
                                else
                                {
                                    throw ex;
                                }
                            }
                            for (; repeatCount > 1; repeatCount--)  // > 1 because this count includes the locus that specified the mapped read count (the previous non-x line) and we already emitted that one.
                            {
                                currentOffset++;
                                processBase(contigName, currentOffset, currentMappedReadCount);
                            }
                        }
                        else
                        {
                            //
                            // A new mapped read count line
                            //
                            if (currentOffset <= 0)
                            {
                                Console.WriteLine("Got unexpected mapped read count line " + line);
                                return false;
                            }
                            try
                            {
                                currentMappedReadCount = Convert.ToInt32(line, 16); // 16 means the string is in hex
                                if (currentMappedReadCount <= 0)
                                {
                                    Console.WriteLine("Bogus current count " + line);
                                    return false;
                                }
                            }
                            catch (SystemException ex)
                            {
                                if (ex is FormatException || ex is ArgumentOutOfRangeException)
                                {
                                    Console.WriteLine("Format exception processing count line " + line);
                                    return false;
                                }
                                else
                                {
                                    throw ex;
                                }
                            }

                            currentOffset++;
                            processBase(contigName, currentOffset, currentMappedReadCount);
                        }

                        continue;
                    }

                    if (fields.Count() != 2)
                    {
                        Console.WriteLine("Saw too many fields in line " + line);
                        return false;
                    }

                    //
                    // An offset + count line
                    //
                    try
                    {
                        currentOffset = Convert.ToInt32(fields[0], 16); // 16 means the string is in hex
                        currentMappedReadCount = Convert.ToInt32(fields[1], 16); // 16 means the string is in hex

                        if (currentOffset <= 0 || currentMappedReadCount <= 0)
                        {
                            Console.WriteLine("Bogus offset + count line " + line);
                            return false;
                        }
                    }
                    catch (SystemException ex)
                    {
                        if (ex is FormatException || ex is ArgumentOutOfRangeException)
                        {
                            Console.WriteLine("Unable to parse offset + count line " + line);
                            return false;
                        }
                        else
                        {
                            throw ex;
                        }
                    }

                    processBase(contigName, currentOffset, currentMappedReadCount);
                }

                if (!sawDone)
                {
                    Console.WriteLine("Truncated allcount file " + filename);
                    return false;
                }

                return true;
            }

            public void Close()
            {
                if (null != allcountReader)
                {
                    allcountReader.Close();
                    compressedStreamReader.Close();
                }
            }

            StreamReader compressedStreamReader = null;
            StreamReader allcountReader = null;
            int numContigs = 0;
            Contig[] contigs = null;

            public readonly string filename;


            class Contig
            {
                public string name = "";
                public long length = -1;
            }

        } // AllcountReader

        public static string DownloadStringWithRetry(WebClient webclient, string address)
        {
            while (true)
            {
                try
                {
                    var reply = webclient.DownloadString(address);
                    if (reply == "")
                    {
                        Console.WriteLine("Server returned the empty string.  Pausing 10 seconds and retrying.");
                        Thread.Sleep(10000);
                    } else
                    {
                        return reply;
                    }
                }
                catch (WebException e)
                {
                    Console.WriteLine("DownloadStringWithRetry: caught WebException: " + e.Message + ".  Pausing 10 seconds and retrying.");
                    Thread.Sleep(10000);
                }
            }
        } // DownloadStringWithRetry

        // Loads gene information from a reference file from UCSC Genome Browser
        public static Dictionary<string, GeneLocationInfo> readKnownGeneFile(string knownGenesFilename)
        {
            var results = new Dictionary<string, GeneLocationInfo>();

            var refFile = CreateStreamReaderWithRetry(knownGenesFilename);

            refFile.ReadLine();    // Skip the header

            string line;
            while (null != (line = refFile.ReadLine()))
            {
                // Parse required fields
                var fields = line.Split('\t');

                if (fields.Count() != 22)
                {
                    Console.WriteLine("LoadGeneLocationInfo: wrong number of fields in file " + fields.Count() + " != 22 " + ": " + line);
                    continue;
                }

                var hugoSymbol = ConvertToNonExcelString(fields[16]);
                var isoform = Isoform.fromFileLine(line);

                if (!results.ContainsKey(hugoSymbol))
                {
                    results.Add(hugoSymbol, new GeneLocationInfo());
                    results[hugoSymbol].hugoSymbol = hugoSymbol;
                    results[hugoSymbol].minLocus = 2000000000;
                    results[hugoSymbol].maxLocus = 0;
                    results[hugoSymbol].chromosome = "";
                }
                results[hugoSymbol].isoforms.Add(isoform);
            }

            //
            // We've collected all of the isoforms with their genes, now compute the min and max loci for each.
            //

            foreach (var resultEntry in results)
            {
                var result = resultEntry.Value;

                foreach (var isoform in result.isoforms)
                {
                    if (!isoform.chromosome.ToLower().Contains("alt") && !isoform.chromosome.ToLower().Contains("random") && !isoform.chromosome.ToLower().Contains("chrun"))
                    {
                        if (result.chromosome != "" && result.chromosome != isoform.chromosome)
                        {
                            //Console.WriteLine("Inconsistent chromosomes for " + result.hugoSymbol + ": " + result.chromosome + " and " + isoform.chromosome);
                            result.inconsistent = true;
                            break;
                        } else
                        {
                            result.chromosome = isoform.chromosome;
                        }

                        result.maxLocus = Math.Max(result.maxLocus, isoform.txEnd);
                        result.minLocus = Math.Min(result.minLocus, isoform.txStart);
                    }
                } // foreach isoform

                if (result.chromosome == "")
                {
                    //Console.WriteLine("No non-alt chromosome found for gene " + result.hugoSymbol);
                    result.inconsistent = true;
                }
            }

            //
            // Mark off "genes" that we know to be inconsistent for whatever reason.
            //
            foreach (var resultEntry in results.Where(x => knownInconsistentGenes.Contains(x.Value.hugoSymbol)))
            {
                resultEntry.Value.inconsistent = true;
            }

            refFile.Close();
            return results;
        } // readKnownGeneFile

        public static string[] knownInconsistentGenes = { "SNORA44", "U6atac", "SNORD78", "FAM95B1", "PGM5-AS1", "SNORA35", "LINC00864", "H19_3", "SNORA79", "RP11-420N3.3", };

        public static string WindowsToLinuxPathname(string windowsPathname)
        {
            string[] components = windowsPathname.Split('\\');


            if (components.Count() < 5)
            {
                Console.WriteLine("WindowsToLinuxPathname: unable to parse " + windowsPathname);
                return @"/dev/null";
            }

            string linuxPathname = @"/mnt/" + components[2] + @"/" + components[3][0];  // components[3][0] is the drive letter, and strips off the $

            for (int i = 4; i < components.Count(); i++)
            {
                linuxPathname += @"/" + components[i];
            }

            return linuxPathname;
        } // WindowsToLinuxPathname

        static public string GetDirectoryPathFromFullyQualifiedFilename(string filename)
        {
            if (filename.Count() < 1 || filename[0] != '\\')
            {
                Console.WriteLine("GetDirectoryPathFromFullyQualifiedFilename -- invalid input " + filename);
                return "";
            }

            string[] pathnameComponents = filename.Split('\\');
            if (pathnameComponents.Count() < 4 || pathnameComponents[0] != "" || pathnameComponents[1] != "" || pathnameComponents[3].Count() != 2 || pathnameComponents[3].Substring(1) != "$")
            {
                Console.WriteLine("GetDirectoryPathFromFullyQualifiedFilename: expected \\ pathname, got " + filename);
                return "";
            }

            string dirName = @"";

            for (int i = 0; i < pathnameComponents.Count() - 1; i++)
            {
                dirName += pathnameComponents[i] + @"\";
            }

            return dirName;
        } // GetDirectoryPathFromFullyQualifiedFilename

        class MAFLoadStatus
        {
            public int nToLoad;
            public int nLoaded = 0;
            public int nFailed = 0;
        } // MAFLoadStatus

        static void ReadMafFileAndAppendToList(string filename, string file_id, List<ASETools.MAFLine> allMAFLines, MAFLoadStatus loadStatus, Configuration configuration)
        {
            var mafLinesForThisFile = ASETools.MAFLine.ReadFile(filename, file_id, true, configuration.isBeatAML);

            lock (loadStatus)
            {
                loadStatus.nLoaded++;

                Console.Write(loadStatus.nLoaded + "/" + loadStatus.nToLoad + " ");

                if (null == mafLinesForThisFile)
                {
                    loadStatus.nFailed++;
                    Console.WriteLine("MAF load of " + filename + " failed.");
                }

                allMAFLines.AddRange(mafLinesForThisFile);
            }
        } // ReadMafFileAndAppendToList

        public static List<MAFLine> LoadMAFs(Configuration configuration, Dictionary<string, DownloadedFile> downloadedFiles, out Dictionary<string, List<MAFLine>> byTumorSampleId)
        {
            var stopwatch = new Stopwatch();
            stopwatch.Start();

            byTumorSampleId = null;

            if (configuration.mafManifestPathname == null || configuration.mafManifestPathname == "")
            {
                Console.WriteLine("The MAF manifest must exist in order to run GenerateCases.  Go back and run the process manager to generate the script to make the maf manifest.");
                return null;
            }

            var mafManifest = ASETools.MAFInfo.LoadMAFManifest(configuration.mafManifestPathname);
            if (null == mafManifest)
            {
                Console.WriteLine("Unable to load maf manifest " + configuration.mafManifestPathname);
                return null;
            }


            foreach (var mafEntry in mafManifest)
            {
                var mafInfo = mafEntry.Value;

                if (!downloadedFiles.ContainsKey(mafInfo.file_id))
                {
                    Console.WriteLine("Missing MAF file " + mafInfo.file_id + ".  Download it (which should be in the script generated by ASEProcessManager) and rerun.");
                    return null;
                }
            }

            Console.Write("Loading MAF entries..");
            var allMAFLines = new List<MAFLine>();

            var threads = new List<Thread>();

            var loadStatus = new MAFLoadStatus();
            loadStatus.nToLoad = mafManifest.Count();


            foreach (var mafEntry in mafManifest)
            {
                var mafInfo = mafEntry.Value;

                threads.Add(new Thread(() => ReadMafFileAndAppendToList(downloadedFiles[mafInfo.file_id].fileInfo.FullName, mafInfo.file_id, allMAFLines, loadStatus, configuration)));
            }


            threads.ForEach(t => t.Start());
            threads.ForEach(t => t.Join());

            if (loadStatus.nFailed > 0)
            {
                Console.WriteLine("Giving up due to failed MAF load(s).");
                return null;
            }

            byTumorSampleId = new Dictionary<string, List<MAFLine>>();

            foreach (MAFLine mafLine in allMAFLines) {
                if (!byTumorSampleId.ContainsKey(mafLine.Tumor_Sample_UUID))
                {
                    byTumorSampleId.Add(mafLine.Tumor_Sample_UUID, new List<MAFLine>());
                }

                byTumorSampleId[mafLine.Tumor_Sample_UUID].Add(mafLine);
            }

            Console.WriteLine(ElapsedTimeInSeconds(stopwatch) + " to load " + allMAFLines.Count() + " mutations.");

            return allMAFLines;
        } // LoadMAFs

        public struct MeanAndStdDev
        {
            public MeanAndStdDev(double mean_, double stdDev_)
            {
                mean = mean_;
                stddev = stdDev_;
            }

            public readonly double mean;
            public readonly double stddev;
        }

        public class RunningMeanAndStdDev
        {
            public RunningMeanAndStdDev() { }

            public void addValue(double value)
            {
                if (0 == n)
                {
                    min = max = value;
                } else
                {
                    min = Math.Min(min, value);
                    max = Math.Max(max, value);
                }
                n++;
                total += value;
                totalOfSquares += value * value;
            }

            public void merge(RunningMeanAndStdDev peer)
            {
                n += peer.n;
                total += peer.total;
                totalOfSquares += peer.totalOfSquares;
            }

            public MeanAndStdDev getMeanAndStdDev()
            {
                if (n == 0) return new MeanAndStdDev(0, 0);

                return new MeanAndStdDev(total / n, Math.Sqrt(n * totalOfSquares - total * total) / n);
            }

            public int getCount()
            {
                return n;
            }

            public double getMin()
            {
                return min;
            }

            public double getMax()
            {
                return max;
            }

            int n = 0;
            double total = 0;
            double totalOfSquares = 0;
            double min;
            double max;
        }

        //
        // Maps chromosomeName -> (offset -> MeanAndStdDev)
        //
        public struct ChromosomeSize : IComparable<ChromosomeSize>
        {
            public ChromosomeSize(string name_, int size_, int centromere_)
            {
                name = name_;
                size = size_;
                centromere = centromere_;
            }

            public int CompareTo(ChromosomeSize other)
            {
                //
                // This is backwards from what you'd ordinarily think because we want them ordered largest to smallest in hopes that that helps out the
                // memory allocator.
                //
                if (size > other.size) return -1;
                if (size == other.size) return 0;
                return 1;
            }

            public readonly string name;
            public readonly int size;
            public readonly int centromere;
        } // ChromosomeSize

        // GRCh38 chromosome sizes.  Centromere locations from Wikipedia.
        static public readonly ChromosomeSize[] chromosomeSizes = {new ChromosomeSize("chr1",  248956422, 123400000), new ChromosomeSize("chr2",  242193529, 93900000), new ChromosomeSize("chr3",  198295559, 90900000),
                                                                   new ChromosomeSize("chr4",  190214555,  50000000), new ChromosomeSize("chr5",  181538259, 48800000), new ChromosomeSize("chr6",  170805979, 59800000),
                                                                   new ChromosomeSize("chr7",  159345973,  60100000), new ChromosomeSize("chr8",  145138636, 45200000), new ChromosomeSize("chr9",  138394717, 43000000),
                                                                   new ChromosomeSize("chr10", 133797422,  39800000), new ChromosomeSize("chr11", 135086622, 53400000), new ChromosomeSize("chr12", 133275309, 35500000),
                                                                   new ChromosomeSize("chr13", 114364328,  17700000), new ChromosomeSize("chr14", 107043718, 17200000), new ChromosomeSize("chr15", 101991189, 19000000),
                                                                   new ChromosomeSize("chr16",  90338345,  36800000), new ChromosomeSize("chr17",  83257441, 25100000), new ChromosomeSize("chr18",  80373285, 18500000),
                                                                   new ChromosomeSize("chr19",  58617616,  26200000), new ChromosomeSize("chr20",  64444167, 28100000), new ChromosomeSize("chr21",  46709983, 12000000),
                                                                   new ChromosomeSize("chr22",  50818468,  15000000), new ChromosomeSize("chrx",  156040895, 61000000), new ChromosomeSize("chry",   57227415, 10400000) };

        static public readonly Dictionary<string, ChromosomeSize> chromosomeSizesByName = new Dictionary<string, ChromosomeSize>();

        public class GeneExpression
        {
            static GeneExpression() // This is a static initializer that runs once at program start time, it's not a constructor.
            {
                regionSizeByRegionSizeIndex[0] = 0;
                regionSizeByRegionSizeIndex[1] = 1000;
                for (int i = 2; i < nRegionSizes; i++)
                {
                    regionSizeByRegionSizeIndex[i] = regionSizeByRegionSizeIndex[i - 1] * 2;
                }

                comparer = StringComparer.OrdinalIgnoreCase;
            }

            public GeneExpression(ASETools.GeneLocationInfo gene_)
            {
                geneLocationInfo = gene_;

                for (int sizeIndex = 0; sizeIndex < nRegionSizes; sizeIndex++)
                {
                    regionalExpressionState[sizeIndex] = new ASETools.RegionalExpressionState();
                    exclusiveRegionalExpressionState[sizeIndex] = new ASETools.RegionalExpressionState();
                }
            }

            public void AddRegionalExpression(int locus, double z, double mu, bool isTumor)
            {
                int distance;
                if (locus >= geneLocationInfo.minLocus && locus <= geneLocationInfo.maxLocus)
                {
                    distance = 0;
                }
                else if (locus < geneLocationInfo.minLocus)
                {
                    distance = geneLocationInfo.minLocus - locus;
                }
                else
                {
                    distance = locus - geneLocationInfo.maxLocus;
                }

                for (int sizeIndex = nRegionSizes - 1; sizeIndex >= 0; sizeIndex--)
                {
                    if (regionSizeByRegionSizeIndex[sizeIndex] < distance)
                    {
                        if (sizeIndex != nRegionSizes - 1)
                        {
                            if (isTumor)
                            {
                                exclusiveRegionalExpressionState[sizeIndex + 1].AddTumorExpression(z, mu);

                            }
                            else
                            {
                                exclusiveRegionalExpressionState[sizeIndex + 1].AddNormalExpression(z, mu);
                            }
                            break;
                        }
                    }
                    if (isTumor)
                    {
                        regionalExpressionState[sizeIndex].AddTumorExpression(z, mu);

                    }
                    else
                    {
                        regionalExpressionState[sizeIndex].AddNormalExpression(z, mu);
                    }
                }

                if (0 == distance)  // Have to special case this, since exclusive gets added when we're one smaller, and there is nothing smaller than sizeIndex 0.
                {
                    if (isTumor)
                    {
                        exclusiveRegionalExpressionState[0].AddTumorExpression(z, mu);

                    }
                    else
                    {
                        exclusiveRegionalExpressionState[0].AddNormalExpression(z, mu);
                    }
                }
            }


            public static int CompareByGeneName(GeneExpression a, GeneExpression b)
            {
                return comparer.Compare(a.geneLocationInfo.hugoSymbol, b.geneLocationInfo.hugoSymbol);
            }

            public const int nRegionSizes = 20;    // Because we have 0 (in the gene), this range is 2^(20 - 2) * 1000 = 262 Mbases on either side, i.e., the entire chromosome
            public static readonly int[] regionSizeByRegionSizeIndex = new int[nRegionSizes];

            public ASETools.RegionalExpressionState[] regionalExpressionState = new ASETools.RegionalExpressionState[nRegionSizes]; // Dimension is log2(regionSize) - 1
            public ASETools.RegionalExpressionState[] exclusiveRegionalExpressionState = new ASETools.RegionalExpressionState[nRegionSizes];  // Expression in this region but not closer, so from log2(regionSize - 1) to log2(regionSize) - 1.  The zero element is the same as regionalExpressionState

            public ASETools.GeneLocationInfo geneLocationInfo;
            public int nonSilentMutationCount = 0;
            public int silentMutationCount = 0;
            public static StringComparer comparer;
        } // GeneExpression

        public class RegionalExpressionState
        {
            // Tumor expression data
            public int nRegionsIncludedTumor = 0;
            public double minTumorExpression = 100000;
            public double maxTumorExpression = -100000;
            public double totalTumorExpression = 0;

            public double minMeanTumorExpression = 100000;
            public double maxMeanTumorExpression = -1;
            public double totalMeanTumorExpression = 0;

            // Normal expression data
            public int nRegionsIncludedNormal = 0;
            public double minNormalExpression = 100000;
            public double maxNormalExpression = -100000;
            public double totalNormalExpression = 0;

            public double minMeanNormalExpression = 100000;
            public double maxMeanNormalExpression = -1;
            public double totalMeanNormalExpression = 0;

            public void AddTumorExpression(double z, double mu)
            {
                nRegionsIncludedTumor++;
                totalTumorExpression += z;

                minTumorExpression = Math.Min(minTumorExpression, z);
                maxTumorExpression = Math.Max(maxTumorExpression, z);

                totalMeanTumorExpression += mu;
                minMeanTumorExpression = Math.Min(minMeanTumorExpression, mu);
                maxMeanTumorExpression = Math.Max(maxMeanTumorExpression, mu);
            }

            public void AddNormalExpression(double z, double mu)
            {
                nRegionsIncludedNormal++;
                totalNormalExpression += z;
                minNormalExpression = Math.Min(minNormalExpression, z);
                maxNormalExpression = Math.Max(maxNormalExpression, z);

                totalMeanNormalExpression += mu;
                minMeanNormalExpression = Math.Min(minMeanNormalExpression, mu);
                maxMeanNormalExpression = Math.Max(maxMeanNormalExpression, mu);
            }
        } // RegionalExpressionState

        // file that contains distance over gene names. Files created from ExpressionNearMutations
        public class RegionalSignalFile
        {
            string fileHeader;
            List<ASETools.GeneExpression> allExpressions;
            ASETools.RegionalExpressionState wholeAutosomeRegionalExpression;
            Dictionary<string, ASETools.RegionalExpressionState> allButThisChromosomeAutosomalRegionalExpressionState;
            ASETools.RegionalExpressionState[] perChromosomeRegionalExpressionState;

            // fileHeader is label for first line in file. For many uses, we print the pipeline name as well as the case id.
            public RegionalSignalFile(string fileHeader_,
                List<ASETools.GeneExpression> allExpressions_,                                                               // regional expression
                ASETools.RegionalExpressionState wholeAutosomeRegionalExpression_,                                           // whole autosome
                Dictionary<string, ASETools.RegionalExpressionState> allButThisChromosomeAutosomalRegionalExpressionState_,  // dictionary of chromosome, expression information
                ASETools.RegionalExpressionState[] perChromosomeRegionalExpressionState_)
            {
                fileHeader = fileHeader_;
                allExpressions = allExpressions_;
                wholeAutosomeRegionalExpression = wholeAutosomeRegionalExpression_;
                allButThisChromosomeAutosomalRegionalExpressionState = allButThisChromosomeAutosomalRegionalExpressionState_;
                perChromosomeRegionalExpressionState = perChromosomeRegionalExpressionState_;
            }

            public void WriteFile(string outputFilename, bool printMu, int minExamplesPerRegion, string columnSuffix, bool isTumor)
            {
                var outputFile = ASETools.CreateStreamWriterWithRetry(outputFilename);

                outputFile.WriteLine(fileHeader);
                outputFile.Write("Gene name\tnon-silent mutation count\tsilent mutation count");

                writeRow(outputFile, allExpressions[0], printMu, minExamplesPerRegion, true, columnSuffix, isTumor);

                outputFile.WriteLine();

                for (int i = 0; i < allExpressions.Count(); i++)
                {
                    outputFile.Write(ASETools.ConvertToExcelString(allExpressions[i].geneLocationInfo.hugoSymbol) + "\t" + allExpressions[i].nonSilentMutationCount + "\t" + allExpressions[i].silentMutationCount);

                    writeRow(outputFile, allExpressions[i], printMu, minExamplesPerRegion, false, "", isTumor);

                    outputFile.WriteLine();
                } // for each gene

                outputFile.WriteLine("**done**");
                outputFile.Close();
            }

            // Write order:
            //	1. Regional values (20)
            //  2. Autosome value (1)
            //  3. If not ASE, regional mean values (20)
            //  4. If not ASE, autosome value (1)
            //  5. Exclusive regional values (20)
            //  6. Per chromosome values (
            void writeRow(
                StreamWriter outputFile,
                ASETools.GeneExpression allExpression,
                bool printMu,
                int minExamplesPerRegion,
                Boolean header,
                string columnSuffix,
                bool isTumor)
            {
                // format suffix for header names
                string columnSuffix_mu = "(" + columnSuffix + " mu)";
                columnSuffix = "(" + columnSuffix + ")";

                // 1. Chromosome specific regional values
                for (int sizeIndex = 0; sizeIndex < ASETools.GeneExpression.nRegionSizes; sizeIndex++)
                {
                    if (header)
                    {
                        outputFile.Write("\t" + ASETools.GeneExpression.regionSizeByRegionSizeIndex[sizeIndex] + columnSuffix);
                        continue;
                    }

                    if (isTumor)
                    {
                        if (allExpression.regionalExpressionState[sizeIndex].nRegionsIncludedTumor >= minExamplesPerRegion)
                        {
                            outputFile.Write("\t" + allExpression.regionalExpressionState[sizeIndex].totalTumorExpression / allExpression.regionalExpressionState[sizeIndex].nRegionsIncludedTumor);
                        }
                        else
                        {
                            outputFile.Write("\t*");
                        }
                    }
                    else if (!isTumor)
                    {
                        if (allExpression.regionalExpressionState[sizeIndex].nRegionsIncludedNormal >= minExamplesPerRegion)
                        {
                            outputFile.Write("\t" + allExpression.regionalExpressionState[sizeIndex].totalNormalExpression / allExpression.regionalExpressionState[sizeIndex].nRegionsIncludedNormal);
                        }
                        else
                        {
                            outputFile.Write("\t*");
                        }
                    }


                }

                // 2. Whole autosome
                if (header)
                {
                    outputFile.Write("\tWhole Autosome " + columnSuffix);
                }
                else if (isTumor)
                {
                    if (wholeAutosomeRegionalExpression.nRegionsIncludedTumor >= minExamplesPerRegion)
                    {
                        outputFile.Write("\t" + wholeAutosomeRegionalExpression.totalTumorExpression / wholeAutosomeRegionalExpression.nRegionsIncludedTumor);
                    }
                    else
                    {
                        outputFile.Write("\t*");
                    }
                }
                else if (!isTumor)
                {
                    if (wholeAutosomeRegionalExpression.nRegionsIncludedNormal >= minExamplesPerRegion)
                    {
                        outputFile.Write("\t" + wholeAutosomeRegionalExpression.totalNormalExpression / wholeAutosomeRegionalExpression.nRegionsIncludedNormal);
                    }
                    else
                    {
                        outputFile.Write("\t*");
                    }
                }
                // 3. Write regional mean values
                if (printMu)
                {

                    for (int sizeIndex = 0; sizeIndex < ASETools.GeneExpression.nRegionSizes; sizeIndex++)
                    {
                        if (header)
                        {
                            outputFile.Write("\t" + ASETools.GeneExpression.regionSizeByRegionSizeIndex[sizeIndex] + columnSuffix_mu);
                            continue;
                        }

                        if (isTumor)
                        {
                            if (allExpression.regionalExpressionState[sizeIndex].nRegionsIncludedTumor >= minExamplesPerRegion)
                            {
                                outputFile.Write("\t" + allExpression.regionalExpressionState[sizeIndex].totalMeanTumorExpression / allExpression.regionalExpressionState[sizeIndex].nRegionsIncludedTumor);
                            }
                            else
                            {
                                outputFile.Write("\t*");
                            }

                        }
                        else
                        {
                            if (allExpression.regionalExpressionState[sizeIndex].nRegionsIncludedNormal >= minExamplesPerRegion)
                            {
                                outputFile.Write("\t" + allExpression.regionalExpressionState[sizeIndex].totalMeanNormalExpression / allExpression.regionalExpressionState[sizeIndex].nRegionsIncludedNormal);
                            }
                            else
                            {
                                outputFile.Write("\t*");
                            }
                        }
                    }

                    // 4. Write mean autosome value
                    if (header)
                    {
                        outputFile.Write("\tWhole Autosome " + columnSuffix_mu);
                    }
                    else if (isTumor)
                    {
                        if (wholeAutosomeRegionalExpression.nRegionsIncludedTumor >= minExamplesPerRegion)
                        {
                            outputFile.Write("\t" + wholeAutosomeRegionalExpression.totalMeanTumorExpression / wholeAutosomeRegionalExpression.nRegionsIncludedTumor);
                        }
                        else
                        {
                            outputFile.Write("\t*");
                        }
                    }
                    else if (!isTumor)
                    {
                        if (wholeAutosomeRegionalExpression.nRegionsIncludedNormal >= minExamplesPerRegion)
                        {
                            outputFile.Write("\t" + wholeAutosomeRegionalExpression.totalNormalExpression / wholeAutosomeRegionalExpression.nRegionsIncludedNormal);
                        }
                        else
                        {
                            outputFile.Write("\t*");
                        }
                    }
                }

                // 5. Write exclusive regional values
                for (int sizeIndex = 0; sizeIndex < ASETools.GeneExpression.nRegionSizes; sizeIndex++)
                {
                    if (header)
                    {
                        outputFile.Write("\t" + ASETools.GeneExpression.regionSizeByRegionSizeIndex[sizeIndex] + " exclusive " + columnSuffix);
                        continue;
                    }
                    if (isTumor)
                    {
                        if (allExpression.exclusiveRegionalExpressionState[sizeIndex].nRegionsIncludedTumor >= minExamplesPerRegion)
                        {
                            outputFile.Write("\t" + allExpression.exclusiveRegionalExpressionState[sizeIndex].totalTumorExpression / allExpression.exclusiveRegionalExpressionState[sizeIndex].nRegionsIncludedTumor);

                        }
                        else
                        {
                            outputFile.Write("\t*");
                        }
                    }
                    else if (!isTumor)
                    {
                        if (allExpression.exclusiveRegionalExpressionState[sizeIndex].nRegionsIncludedNormal >= minExamplesPerRegion)
                        {
                            outputFile.Write("\t" + allExpression.exclusiveRegionalExpressionState[sizeIndex].totalNormalExpression / allExpression.exclusiveRegionalExpressionState[sizeIndex].nRegionsIncludedNormal);
                        }
                        else
                        {
                            outputFile.Write("\t*");
                        }
                    }
                }

                // 6. Whole autosome exclusive values
                if (header)
                {
                    outputFile.Write("\tWhole Autosome exclusive " + columnSuffix);
                }
                else if (isTumor)
                {
                    if (allButThisChromosomeAutosomalRegionalExpressionState[allExpression.geneLocationInfo.chromosome].nRegionsIncludedTumor >= minExamplesPerRegion)
                    {
                        outputFile.Write("\t" + allButThisChromosomeAutosomalRegionalExpressionState[allExpression.geneLocationInfo.chromosome].totalTumorExpression / allButThisChromosomeAutosomalRegionalExpressionState[allExpression.geneLocationInfo.chromosome].nRegionsIncludedTumor);
                    }
                    else
                    {
                        outputFile.Write("\t*");
                    }
                }
                else if (!isTumor)
                {
                    if (allButThisChromosomeAutosomalRegionalExpressionState[allExpression.geneLocationInfo.chromosome].nRegionsIncludedNormal >= minExamplesPerRegion)
                    {
                        outputFile.Write("\t" + allButThisChromosomeAutosomalRegionalExpressionState[allExpression.geneLocationInfo.chromosome].totalNormalExpression / allButThisChromosomeAutosomalRegionalExpressionState[allExpression.geneLocationInfo.chromosome].nRegionsIncludedNormal);
                    }
                    else
                    {
                        outputFile.Write("\t*");
                    }
                }

                if (printMu)
                {
                    // 7. Exlusive regional means
                    for (int sizeIndex = 0; sizeIndex < ASETools.GeneExpression.nRegionSizes; sizeIndex++)
                    {
                        if (header)
                        {
                            outputFile.Write("\t" + ASETools.GeneExpression.regionSizeByRegionSizeIndex[sizeIndex] + " exclusive " + columnSuffix_mu);
                            continue;
                        }
                        if (isTumor)
                        {
                            if (allExpression.exclusiveRegionalExpressionState[sizeIndex].nRegionsIncludedTumor >= minExamplesPerRegion)
                            {
                                outputFile.Write("\t" + allExpression.exclusiveRegionalExpressionState[sizeIndex].totalMeanTumorExpression / allExpression.exclusiveRegionalExpressionState[sizeIndex].nRegionsIncludedTumor);
                            }
                            else
                            {
                                outputFile.Write("\t*");
                            }
                        }
                        else if (!isTumor)
                        {
                            if (allExpression.exclusiveRegionalExpressionState[sizeIndex].nRegionsIncludedNormal >= minExamplesPerRegion)
                            {
                                outputFile.Write("\t" + allExpression.exclusiveRegionalExpressionState[sizeIndex].totalMeanNormalExpression / allExpression.exclusiveRegionalExpressionState[sizeIndex].nRegionsIncludedNormal);
                            }
                            else
                            {
                                outputFile.Write("\t*");
                            }
                        }
                    }

                    // 8. Write out chromosome exclusive autosome means
                    if (header)
                    {
                        outputFile.Write("\tWhole Autosome exclusive " + columnSuffix_mu);
                    }
                    else if (isTumor)
                    {
                        if (allButThisChromosomeAutosomalRegionalExpressionState[allExpression.geneLocationInfo.chromosome].nRegionsIncludedTumor >= minExamplesPerRegion)
                        {
                            outputFile.Write("\t" + allButThisChromosomeAutosomalRegionalExpressionState[allExpression.geneLocationInfo.chromosome].totalMeanTumorExpression / allButThisChromosomeAutosomalRegionalExpressionState[allExpression.geneLocationInfo.chromosome].nRegionsIncludedTumor);
                        }
                        else
                        {
                            outputFile.Write("\t*");
                        }
                    }
                    else if (!isTumor)
                    {
                        if (allButThisChromosomeAutosomalRegionalExpressionState[allExpression.geneLocationInfo.chromosome].nRegionsIncludedNormal >= minExamplesPerRegion)
                        {
                            outputFile.Write("\t" + allButThisChromosomeAutosomalRegionalExpressionState[allExpression.geneLocationInfo.chromosome].totalMeanNormalExpression / allButThisChromosomeAutosomalRegionalExpressionState[allExpression.geneLocationInfo.chromosome].nRegionsIncludedNormal);
                        }
                        else
                        {
                            outputFile.Write("\t*");
                        }
                    }
                }

                // 9. Write out per chromosome
                for (int whichChromosome = 0; whichChromosome < ASETools.nHumanNuclearChromosomes; whichChromosome++)
                {
                    if (header)
                    {
                        outputFile.Write("\t" + ASETools.ChromosomeIndexToName(whichChromosome, true) + columnSuffix);
                        continue;
                    }
                    if (isTumor)
                    {
                        if (perChromosomeRegionalExpressionState[whichChromosome].nRegionsIncludedTumor >= minExamplesPerRegion)
                        {
                            outputFile.Write("\t" + perChromosomeRegionalExpressionState[whichChromosome].totalTumorExpression / perChromosomeRegionalExpressionState[whichChromosome].nRegionsIncludedTumor);
                        }
                        else
                        {
                            outputFile.Write("\t*");
                        }
                    }
                    else if (!isTumor)
                    {
                        if (perChromosomeRegionalExpressionState[whichChromosome].nRegionsIncludedNormal >= minExamplesPerRegion)
                        {
                            outputFile.Write("\t" + perChromosomeRegionalExpressionState[whichChromosome].totalNormalExpression / perChromosomeRegionalExpressionState[whichChromosome].nRegionsIncludedNormal);
                        }
                        else
                        {
                            outputFile.Write("\t*");
                        }
                    }
                }

                // 10. Write out per chromosome mean
                if (printMu)
                {
                    for (int whichChromosome = 0; whichChromosome < ASETools.nHumanNuclearChromosomes; whichChromosome++)
                    {
                        if (header)
                        {
                            outputFile.Write("\t" + ASETools.ChromosomeIndexToName(whichChromosome, true) + columnSuffix_mu);
                            continue;
                        }
                        if (isTumor)
                        {
                            if (perChromosomeRegionalExpressionState[whichChromosome].nRegionsIncludedTumor >= minExamplesPerRegion)
                            {
                                outputFile.Write("\t" + perChromosomeRegionalExpressionState[whichChromosome].totalMeanTumorExpression / perChromosomeRegionalExpressionState[whichChromosome].nRegionsIncludedTumor);
                            }
                            else
                            {
                                outputFile.Write("\t*");
                            }
                        }
                        else if (!isTumor)
                        {
                            if (perChromosomeRegionalExpressionState[whichChromosome].nRegionsIncludedNormal >= minExamplesPerRegion)
                            {
                                outputFile.Write("\t" + perChromosomeRegionalExpressionState[whichChromosome].totalMeanNormalExpression / perChromosomeRegionalExpressionState[whichChromosome].nRegionsIncludedNormal);
                            }
                            else
                            {
                                outputFile.Write("\t*");
                            }
                        }
                    }
                }
            }

            // Returns a tuple of ExpressionMap, and list of labels corresponding to the expression values in the dictionary
            public static Tuple<Dictionary<string, double[]>, List<string>> ReadFile(string filename, bool skipFirstLine = true, bool includeMutationCount = false, List<string> genesToInclude = null)
            {
                // keeps track of index of distances
                List<string> index = new List<string>();

                // for each gene, maintain a list of distance, expression tuples
                Dictionary<string, double[]> expressionMap = new Dictionary<string, double[]>();

                var reader = CreateStreamReaderWithRetry(filename);

                string line;

                List<string> excelizedGenesToInclude;
                if (genesToInclude == null)
                {
                    excelizedGenesToInclude = null;
                } else
                {
                    excelizedGenesToInclude = genesToInclude.Select(x => ConvertToExcelString(x)).ToList();
                }

                // skip first line if version or other info resides there
                if (skipFirstLine)
                    reader.ReadLine();

                // read in header
                var header = reader.ReadLine().Split('\t');

                // skip gene name and mutation count. Get all distance labels and keep them in the index.
                for (var i = 2; i < header.Count(); i++)
                {
                    index.Add(header[i]);
                }

                while (null != (line = reader.ReadLine()))
                {
                    if (line.Trim() == "**done**")
                    {
                        break;
                    }

                    if (excelizedGenesToInclude != null)
                    {
                        bool includeThisGene = false;

                        foreach (var gene in excelizedGenesToInclude)
                        {
                            if (line.StartsWith(gene))
                            {
                                includeThisGene = true;
                                break;
                            }
                        }

                        if (!includeThisGene)
                        {
                            continue;
                        }
                    }

                    string[] fields = line.Split('\t');

                    var hugoSymbol = ConvertToNonExcelString(fields[0]);

                    // get rid of gene name and optionally mutation count
                    fields = fields.Skip(includeMutationCount ? 1 : 3).ToArray();

                    // the rest of the fields match the index fields

                    if (fields.Length != index.Count() + (includeMutationCount ? 1 : 0))
                    {
                        throw new Exception("header does not match field length for file " + filename);
                    }

                    double[] numericFields = fields.Select(r => {
                        // Set default value
                        double value = double.NegativeInfinity;

                        if (r != "*" && r != "-Infinity")
                        {
                            value = Convert.ToDouble(r);
                        }
                        return value;
                    }).ToArray();

                    // add gene to dictionary
                    expressionMap.Add(hugoSymbol, numericFields);
                }
                return new Tuple<Dictionary<string, double[]>, List<string>>(expressionMap, index);
            }

            //
            // The other direction of our conversion of * to -Infinity
            //
            public static string stringForDouble(double value)
            {
                if (value == double.NegativeInfinity) return "*";
                return value.ToString();
            }
        } // RegionalSignalFile

        public static string stringForDouble(double value)
        {
            if (value == double.NegativeInfinity)
            {
                return "*";
            }
            return Convert.ToString(value);
        }

        public class ASESignalLine // Should probably add the per-chromosome stuff here, too.
        {
            public static Dictionary<string, ASESignalLine> ReadFile(string filename, List<string> genesToInclude = null)
            {
                var regionalSignalFile = RegionalSignalFile.ReadFile(filename, true, true, genesToInclude);

                var retVal = new Dictionary<string, ASESignalLine>();

                foreach (var regionalSignalEntry in regionalSignalFile.Item1)
                {
                    var regionalSignal = regionalSignalEntry.Value;

                    var inclusiveASE = new double[ASETools.nRegions];
                    Array.Copy(regionalSignal, 1, inclusiveASE, 0, ASETools.nRegions);

                    var exclusiveASE = new double[ASETools.nRegions];
                    Array.Copy(regionalSignal, 1 + ASETools.nRegions, exclusiveASE, 0, ASETools.nRegions);

                    retVal.Add(regionalSignalEntry.Key, new ASESignalLine((int)regionalSignal[0], regionalSignalEntry.Key, inclusiveASE, exclusiveASE));
                }

                return retVal;
            }

            ASESignalLine(int nMutations_, string hugo_symbol_, double[] inclusiveASE_, double[] exclusiveASE_)
            {
                nMutations = nMutations_;
                nMutationsIndex = (nMutations < 2) ? nMutations : 2;
                hugo_symbol = hugo_symbol_;
                inclusiveASE = inclusiveASE_;
                exclusiveASE = exclusiveASE_;
            }

            public readonly int nMutations;
            public readonly int nMutationsIndex;    // 0, 1 or 2 depending on 0, 1, or > 1 mutations
            public readonly string hugo_symbol;
            public readonly double[] inclusiveASE;
            public readonly double[] exclusiveASE;
        } // ASESignalLine


        //
        // This breaks out the array of doubles from the entires in a RegionalSignalFile into a class that names them for the allele-specific expression case
        //
        public class AlleleSpecificSignal
        {
            public readonly string Hugo_Symbol;
            public readonly int nNonSilentMutations;
            public readonly int nSilentMutations;
            public readonly double[] nonExclusive = new double[nRegions];
            public readonly double[] exclusive = new double[nRegions];
            public readonly double[] chromosomes = new double[nHumanNuclearChromosomes];

            public AlleleSpecificSignal(string Hugo_Symbol_, double[] numericValues)
            {
                Hugo_Symbol = Hugo_Symbol_;
                if (numericValues.Count() != 2 * nRegions + nHumanNuclearChromosomes + 2)
                {
                    throw new FormatException();
                }

                int nextValueIndex = 0;
                nNonSilentMutations = (int)numericValues[nextValueIndex];
                nextValueIndex++;

                nSilentMutations = (int)numericValues[nextValueIndex];
                nextValueIndex++;

                for (int i = 0; i < nRegions; i++)
                {
                    nonExclusive[i] = numericValues[nextValueIndex];
                    nextValueIndex++;
                }

                for (int i = 0; i < nRegions; i++)
                {
                    exclusive[i] = numericValues[nextValueIndex];
                    nextValueIndex++;
                }

                for (int i = 0; i < nHumanNuclearChromosomes; i++)
                {
                    chromosomes[i] = numericValues[nextValueIndex];
                    nextValueIndex++;
                }
            }

            public static Dictionary<string, AlleleSpecificSignal> fromRegionalSignalFile(Dictionary<string, double[]> regionalSignalFile)
            {
                var retVal = new Dictionary<string, AlleleSpecificSignal>();
                foreach (var entry in regionalSignalFile)
                {
                    retVal.Add(entry.Key, new AlleleSpecificSignal(entry.Key, entry.Value));
                }

                return retVal;
            }

            public string OutputString()
            {
                string retVal = Hugo_Symbol + "\t" + nNonSilentMutations + "\t" + nSilentMutations;
                for (int i = 0; i < nRegions; i++)
                {
                    retVal += "\t" + RegionalSignalFile.stringForDouble(nonExclusive[i]);
                }

                for (int i = 0; i < nRegions; i++)
                {
                    retVal += "\t" + RegionalSignalFile.stringForDouble(exclusive[i]);
                }

                for (int i = 0; i < nHumanNuclearChromosomes; i++)
                {
                    retVal += "\t" + RegionalSignalFile.stringForDouble(chromosomes[i]);
                }

                return retVal;
            }
        }


        //
        // The content of an expression file.  We used to just return the more obvious Dictionary<string, Dictionary<int, MeanAndStdDev>>, but the inside dictionaries got big enough to make
        // something bad happen in the CLR runtime, even on 64 bit with gcAllowVeryLargeObjects set.  So, instead we add one more level of indirection to get around dictionaries with
        // tens of millions of entries.
        //
        public class ExpressionFile
        {
            public ExpressionFile() { }

            public void LoadFromFile(string filename)
            {
                string currentContig = null;

                var reader = CreateStreamReaderWithRetry(filename);

                string statLine;

                while (null != (statLine = reader.ReadLine()))
                {
                    string[] statLineFields = statLine.Split('\t');
                    if (statLineFields.Count() == 1)
                    {
                        //
                        // It's a chromosome header.
                        //

                        currentContig = statLineFields[0].ToLower();
                        if (expression.ContainsKey(currentContig))
                        {
                            throw new Exception("Duplicate contig name in expression file: " + currentContig);
                        }
                        expression.Add(currentContig, new Dictionary<int, Dictionary<int, MeanAndStdDev>>());
                        higestOffsetForEachContig.Add(currentContig, 0);
                    } else
                    {
                        int chromosomeOffset = Convert.ToInt32(statLineFields[0]);
                        if (!expression[currentContig].ContainsKey(chromosomeOffset / subchunkSize))
                        {
                            expression[currentContig].Add(chromosomeOffset / subchunkSize, new Dictionary<int, MeanAndStdDev>());
                        }

                        expression[currentContig][chromosomeOffset / subchunkSize].Add(chromosomeOffset, new MeanAndStdDev(Convert.ToDouble(statLineFields[2]), Convert.ToDouble(statLineFields[3])));
                        higestOffsetForEachContig[currentContig] = Math.Max(higestOffsetForEachContig[currentContig], chromosomeOffset);
                    }
                }
            }
            public bool valueKnown(string contigName, int chromosomeOffset)
            {
                contigName = contigName.ToLower();
                return !(!expression.ContainsKey(contigName) || !expression[contigName].ContainsKey(chromosomeOffset / subchunkSize) || !expression[contigName][chromosomeOffset / subchunkSize].ContainsKey(chromosomeOffset));
            }

            public bool getValue(string contigName, int chromosomeOffset, out MeanAndStdDev meanAndStdDev)
            {
                contigName = contigName.ToLower();
                if (!expression.ContainsKey(contigName) || !expression[contigName].ContainsKey(chromosomeOffset / subchunkSize) || !expression[contigName][chromosomeOffset / subchunkSize].ContainsKey(chromosomeOffset))
                {
                    meanAndStdDev = new MeanAndStdDev(-1, -1);
                    return false;
                }

                meanAndStdDev = expression[contigName][chromosomeOffset / subchunkSize][chromosomeOffset];
                return true;
            }

            public int getHigestOffsetForContig(string contigName)
            {
                return higestOffsetForEachContig[contigName];
            }

            Dictionary<string, Dictionary<int, Dictionary<int, MeanAndStdDev>>> expression = new Dictionary<string, Dictionary<int, Dictionary<int, MeanAndStdDev>>>();
            public Dictionary<string, int> higestOffsetForEachContig = new Dictionary<string, int>();
            const int subchunkSize = 100000;
        }

        //
        // Excel has the bad habit of taking strings that looks like dates and converting them into actual dates.  This is a problem for some gene names ("MARCH1"), which is
        // bad enough that it's actually made into the literature (not to mention the mafs that I imported).  This takes a string and converts it into a format that Excel
        // will not munge.
        //
        public static string ConvertToExcelString(string input)
        {
            if (input.Count() < 3 || input.Substring(0, 2) != "=\"" || input[input.Count() - 1] != '\"')
            {
                return "=\"" + input + "\"";
            }

            //
            // It's already in excel format, just keep it.
            //
            return input;
        } // ConvertToExcelString

        public static string ConvertToNonExcelString(string input)
        {
            if (input.Count() < 3 || input.Substring(0, 2) != "=\"" || input[input.Count() - 1] != '\"')
            {
                //
                // It's not an excel string, just return it.
                //
                return input;
            }

            return input.Substring(2, input.Count() - 3);
        } // ConvertToNonExcelString

        public class ConsolodatedFileReader
        {
            public ConsolodatedFileReader()
            {
            }

            public bool open(string filename_) {
                filename = filename_;
                filestream = File.OpenRead(filename);

                var indexReader = CreateStreamReaderWithRetry(filename + ".index");

                string line;
                bool sawDone = false;
                while (null != (line = indexReader.ReadLine())) {
                    if (sawDone)
                    {
                        Console.WriteLine("ASETools.ConsolodatedFileReader: index file continues after **done**: " + filename);
                        return false;
                    }
                    if (line.Trim() == "**done**")
                    {
                        sawDone = true;
                        continue;
                    }
                    var fields = line.Split('\t');
                    if (fields.Count() != 3)
                    {
                        Console.WriteLine("ASETools.ConsolodatedFileReader: saw index line with an incorrect number of fields: " + line + ": " + filename);
                        return false;
                    }

                    if (subfiles.ContainsKey(fields[0]))
                    {
                        //
                        // We sometimes get more than one instance of the same subfile if, for instance, there are two MAF entries for the
                        // same mutation locus.  It's OK as long as they have the same subfile size.
                        //
                        if (subfiles[fields[0]].size == Convert.ToInt64(fields[2]))
                        {
                            continue;
                        }
                        Console.WriteLine("ASETools.ConsolodatedFileReader: index file has more than one instance of subfile " + fields[0] + " with different sizes " + subfiles[fields[0]].size + " and " + Convert.ToInt32(fields[2]));
                        return false;
                    }

                    try
                    {
                        subfiles.Add(fields[0], new SubFile(fields[0], Convert.ToInt64(fields[1]), Convert.ToInt64(fields[2])));
                    } catch (FormatException)
                    {
                        Console.WriteLine("ASETools.ConsolodatedFileReader:Unparsable line in index file: " + line + ": " + filename);
                        return false;
                    }
                }

                indexReader.Close();

                if (!sawDone)
                {
                    Console.WriteLine("ASETools.ConsolodatedFile: truncated index file: " + filename);
                    return false;
                }

                return true;
            }

            public bool isSubfileTooBigToRead(string subfileName)
            {
                if (!subfiles.ContainsKey(subfileName))
                {
                    return false;   // Doesn't exist, so it's not too big to read
                }

                return subfiles[subfileName].size > int.MaxValue;
            }

            public StreamReader getSubfile(string subfileName)
            {
                if (!subfiles.ContainsKey(subfileName))
                {
                    return null;
                }

                filestream.Seek(subfiles[subfileName].offset, SeekOrigin.Begin);

                long totalToRead = subfiles[subfileName].size;
                if (totalToRead > int.MaxValue)
                {
                    Console.WriteLine("ConsolodatedFileReader.GetSubfile: Unable to create subfile because the stream is too big: " + totalToRead + " > " + int.MaxValue + ".  File " + filename + " subfile " + subfileName + ".");
                    return null;
                }
                var buffer = new byte[totalToRead];
                long amountRead = 0;
                int maxReadSize = 2146435072;   // 2 GiB - 1MiB

                if (totalToRead < maxReadSize)
                {
                    int amountActuallyRead = filestream.Read(buffer, 0, (int)totalToRead);
                    if (amountActuallyRead != totalToRead)
                    {
                        Console.WriteLine("getSubfile: Filestream.Read didn't read the whole requested size, " + amountActuallyRead + " != " + amountRead); // If this happens, you can always fix the code to loop.
                        return null;
                    }
                }
                else
                {
                    //
                    // Because Filestream.Read will only read an int's worth at a time, we have to break it into
                    // chunks and read it into a temp buffer and then copy.
                    //
                    while (amountRead < totalToRead)
                    {
                        int amountToRead;
                        if (totalToRead - amountRead < maxReadSize)
                        {
                            amountToRead = (int)(totalToRead - amountRead);
                        }
                        else
                        {
                            amountToRead = maxReadSize;
                        }

                        var tempBuffer = new byte[amountToRead];
                        int amountActuallyRead = filestream.Read(tempBuffer, 0, amountToRead);
                        if (0 == amountActuallyRead)
                        {
                            Console.WriteLine("GetSubfile: error reading from subfile, got 0 back from Filestream.Read()");
                            return null;
                        }
                        Array.Copy(tempBuffer, 0, buffer, amountRead, amountActuallyRead);
                        amountRead += amountActuallyRead;
                    } // While we have more to read
                } // if we can do it in one read or not

                return new StreamReader(new MemoryStream(buffer));
            }

            class SubFile
            {
                public SubFile(string name_, long offset_, long size_)
                {
                    name = name_;
                    offset = offset_;
                    size = size_;
                }

                public string name;
                public long offset;
                public long size;
            }

            public List<string> allSubfiles()
            {
                return subfiles.Select(_ => _.Key).ToList();
            }

            FileStream filestream;
            Dictionary<string, SubFile> subfiles = new Dictionary<string, SubFile>();
            string filename;
        } // ConsolodatedFileReader

        public class Genome
        {
            public Genome()
            {
            }

            private string formatContig(string contig)
            {
                //
                // Try adding "chr" to the contig.
                //
                if (!(contig.Count() > 3 && contig.Substring(0, 3) == "chr"))
                {
                    contig = "chr" + contig;
                }
                // use lower case to avoid issues with _random contigs
                return contig.ToLower();
            }

            public bool load(string snapIndexDirectory)
            {
                var metadataFilename = snapIndexDirectory + @"\GenomeIndex";

                if (!File.Exists(metadataFilename))
                {
                    Console.WriteLine("Genome.load: can't find SNAP index metadata in expected file " + metadataFilename);
                    return false;
                }

                var metadata = File.ReadAllLines(snapIndexDirectory + @"\GenomeIndex");

                if (metadata.Count() != 1)
                {
                    Console.WriteLine("Genome.load: SNAP index metadata has other than one line (" + metadataFilename + ")");
                    return false;
                }

                var fields = metadata[0].Split(' ');
                if (fields.Count() < 2 || fields[0] != "5")
                {
                    Console.WriteLine("Genome.load: Wrong version or corrupt snap index in directory " + snapIndexDirectory);
                    return false;
                }

                if (fields.Count() != 10)
                {
                    Console.WriteLine("Genome.load: wrong field count in snap index metadata: " + metadataFilename);
                    return false;
                }

                try
                {
                    chromosomePadding = Convert.ToInt32(fields[5]);
                } catch (FormatException)
                {
                    Console.WriteLine("Genome.load: error parsing snap index metadata line in " + metadataFilename);
                    return false;
                }

                var genomeFilename = snapIndexDirectory + @"\Genome";

                var genomeReader = CreateStreamReaderWithRetry(genomeFilename);
                var line = genomeReader.ReadLine();

                if (null == line)
                {
                    Console.WriteLine("Genome.load: empty genome file " + genomeFilename);
                    return false;
                }

                fields = line.Split(' ');

                if (fields.Count() != 2)
                {
                    Console.WriteLine("Genome.load: wrong field count on genome header line in file " + genomeFilename);
                    return false;
                }

                try
                {
                    genomeLength = Convert.ToInt64(fields[0]);
                    nContigs = Convert.ToInt32(fields[1]);

                    for (int i = 0; i < nContigs; i++)
                    {
                        line = genomeReader.ReadLine();

                        if (null == line)
                        {
                            Console.WriteLine("Genome.load: truncated genome file " + genomeFilename);
                            return false;
                        }

                        fields = line.Split(' ');
                        if (fields.Count() != 2)
                        {
                            Console.WriteLine("Genome.load: corrupt contig line in " + genomeFilename);
                            return false;
                        }

                        var contig = formatContig(fields[1]);

                        if (contigsByName.ContainsKey(contig))
                        {
                            Console.WriteLine("Genome.load: duplicate contig (" + contig + " in header of " + genomeFilename);
                            return false;
                        }

                        contigsByName.Add(contig, new Contig(contig, Convert.ToInt64(fields[0])));
                        contigsInOrder.Add(contigsByName[contig]);
                    }
                } catch (FormatException)
                {
                    Console.WriteLine("Genome.load: Format error parsing header of genome file " + genomeFilename);
                    return false;
                }

                for (int i = 0; i < contigsInOrder.Count() - 1; i++)
                {
                    contigsInOrder[i].size = (int)(contigsInOrder[i + 1].offset - contigsInOrder[i].offset - chromosomePadding);
                }
                contigsInOrder[contigsInOrder.Count() - 1].size = (int)(genomeLength - contigsInOrder[contigsInOrder.Count() - 1].offset);    // Yes, this throws an exception for a single contig genome, but we know we don't have that.

                long currentOffset = 0;
                for (int i = 0; i < contigsInOrder.Count(); i++)
                {
                    int charsRead;
                    if (currentOffset < contigsInOrder[i].offset)
                    {
                        int sizeToSkip = (int)(contigsInOrder[i].offset - currentOffset);
                        var garbageBuffer = new char[sizeToSkip];
                        charsRead = genomeReader.ReadBlock(garbageBuffer, 0, sizeToSkip);

                        if (charsRead != sizeToSkip)
                        {
                            Console.WriteLine("Genome.load: failed to read in as much padding as expected, " + charsRead + " != " + sizeToSkip);
                            return false;
                        }

                        currentOffset += charsRead;
                    }

                    contigsInOrder[i].data = new char[contigsInOrder[i].size];

                    charsRead = genomeReader.ReadBlock(contigsInOrder[i].data, 0, contigsInOrder[i].size);

                    if (charsRead != contigsInOrder[i].size)
                    {
                        Console.WriteLine("Genome.load: read unexpected number of bytes " + charsRead + " != " + contigsInOrder[i].size);
                        return false;
                    }

                    currentOffset += charsRead;
                }

                return true;
            }

            public char getBase(string contigName, int offset)
            {
                contigName = formatContig(contigName);
                return contigsByName[contigName].data[offset - 1];   // Offset - 1 because genome coordinates are 1 based, while C# is 0 based.
            }

            public bool isAContig(string contigName)
            {
                return contigsByName.ContainsKey(formatContig(contigName));
            }

            public int getContigLength(string contigName)
            {
                return contigsByName[formatContig(contigName)].size;
            }

            int chromosomePadding;
            long genomeLength;
            int nContigs;
            List<Contig> contigsInOrder = new List<Contig>();
            Dictionary<string, Contig> contigsByName = new Dictionary<string, Contig>();

            class Contig
            {
                public Contig(string name_, long offset_)
                {
                    name = name_;
                    offset = offset_;
                }

                public string name;
                public long offset;
                public int size;

                public char[] data;
            }
        }// Genome

        public class RandomizingStreamWriter
        {
            public RandomizingStreamWriter(StreamWriter underlyingWriter_)
            {
                underlyingWriter = underlyingWriter_;
            }

            public void Write(string newData)
            {
                unfinishedLine += newData;
            }

            public void WriteLine(string newData)
            {
                unfinishedLine += newData;
                queuedLines.Add(unfinishedLine);
                unfinishedLine = "";
            }

            public void WriteLine()
            {
                WriteLine("");
            }

            public void Close()
            {
                if (unfinishedLine != "")
                {
                    WriteLine("");
                }

                var random = new Random();

                while (queuedLines.Count() > 0)
                {
                    int lineToEmit = random.Next(queuedLines.Count());
                    underlyingWriter.WriteLine(queuedLines[lineToEmit]);
                    queuedLines.RemoveAt(lineToEmit);
                }

                underlyingWriter.Close();
            }

            StreamWriter underlyingWriter;
            string unfinishedLine = "";
            List<string> queuedLines = new List<string>();
        }

        public class SelectedVariant
        {
            SelectedVariant(string contig_, int locus_, char referenceBase_, char altBase_)
            {
                contig = contig_;
                locus = locus_;
                referenceBase = referenceBase_;
                altBase = altBase_;
            }

            public static List<SelectedVariant> LoadFromFile(string filename)
            {
                //
                // Alas, this file type doesn't have a header, so we can't use HeaderizedFile.
                //

                var retVal = new List<SelectedVariant>();

                bool seenDone = false;

                var reader = CreateStreamReaderWithRetry(filename);

                var line = reader.ReadLine();

                if (null == line)
                {
                    Console.WriteLine("Empty selected variants file " + filename);
                    return null;
                }

                while (null != (line = reader.ReadLine()))
                {
                    if (seenDone)
                    {
                        Console.WriteLine("Saw data after **done** in " + filename);
                        return null;
                    }

                    if (line == "**done**")
                    {
                        seenDone = true;
                        continue;
                    }

                    var fields = line.Split('\t');
                    if (fields.Count() < 7)
                    {
                        Console.WriteLine("Not enough fields in a line of " + filename + ".  It's probably truncated.");
                        return null;
                    }

                    try
                    {
                        retVal.Add(new SelectedVariant(fields[0], Convert.ToInt32(fields[1]), fields[5][0], fields[6][0]));
                    } catch (FormatException) {
                        Console.WriteLine("Format exception in " + filename + ".  It's probably truncated.");
                        return null;
                    }
                }

                if (!seenDone)
                {
                    Console.WriteLine("Truncated selected variants file " + filename);
                    return null;
                }

                return retVal;
            }

            public string getExtractedReadsExtension()
            {
                return "-" + contig + "-" + locus;
            }

            public readonly string contig;
            public int locus;
            public readonly char referenceBase;
            public readonly char altBase;
        } // SelectedVariant

        public class MethylationCounts
        {
            public MethylationCounts(int nMethylated_, int nUnmethylated_, int nNeither_)
            {
                nMethylated = nMethylated_;
                nUnmethylated = nUnmethylated_;
                nNeither = nNeither_;
            }

            public readonly int nMethylated;
            public readonly int nUnmethylated;
            public readonly int nNeither;

            public override string ToString()
            {
                return (nMethylated.ToString() + '\t' + nUnmethylated.ToString() + '\t' + nNeither.ToString());
            }


            // Returns ReadCounts for methylation counts of each allele. First allele is ref, second allele is alt
            public static Tuple<MethylationCounts, MethylationCounts> ComputeMethylationReadCounts(string selectedReadsFilename,
                string contig,
                int variant_start_position,
                string variant_reference_allele,
                string variant_alt_allele,
                string variantType,
                int cpg_start_position,
                string subfileName,
                Genome genome)
            {

                // holds all the sequences matching ref
                List<Dictionary<int, char>> matchingRef = new List<Dictionary<int, char>>();

                // holds all the sequences matching alt
                List<Dictionary<int, char>> matchingAlt = new List<Dictionary<int, char>>();


                var consolodatedFile = new ASETools.ConsolodatedFileReader();

                if (!consolodatedFile.open(selectedReadsFilename))
                {
                    Console.WriteLine("Unable to open reads at selected variants file " + selectedReadsFilename);
                    return null;
                }

                var subfileReader = consolodatedFile.getSubfile(subfileName);
                if (null == subfileReader)
                {
                    Console.WriteLine("ComputeReadCounts: no subfile named " + subfileName);
                    return null;
                }

                // Divide reads up into alt and ref
                var padding = 60;
                int[] posArray = new int[padding];
                for (int i = 0; i < padding; i++)
                {
                    posArray[i] = variant_start_position - padding / 2 + i;
                }

                // get reference and alt sequences
                var reference_bases = posArray
                    .Select(r => new Tuple<int, char>(r, genome.getBase(contig, r))).ToDictionary(x => x.Item1, x => x.Item2);
                var alt_bases = posArray
                    .Select(r => new Tuple<int, char>(r, genome.getBase(contig, r))).ToDictionary(x => x.Item1, x => x.Item2);

                // If the variant we are looking at is an INS, we have to shift the variant position 
                // by one to match the SamLine insertion location
                int start = variantType == "INS" ? variant_start_position + 1 : variant_start_position;

                // modify alt_bases based on VariantType
                if (variantType == "DEL")
                {
                    // replace missing bases with "N" so we can do a direct comparison
                    for (var i = 0; i < variant_reference_allele.Length; i++)
                    {
                        alt_bases[start + i] = 'N';
                    }
                }
                else if (variantType != "INS")
                {
                    // for NP: replace bases for alt
                    for (var i = 0; i < variant_alt_allele.Length; i++)
                    {
                        alt_bases[start + i] = variant_alt_allele[i];
                    }
                }
                else
                {
                    // insert INS into dictionary
                    KeyValuePair<int, char>[] ins = new KeyValuePair<int, char>[variant_alt_allele.Length];

                    var chArr = variant_alt_allele.ToCharArray();
                    for (int i = 0; i < chArr.Length; i++)
                    {
                        ins[i] = new KeyValuePair<int, char>(i + start, chArr[i]);
                    }

                    alt_bases =
                        alt_bases.Where(r => r.Key < start).ToArray().Concat(ins.ToArray())
                        .Concat(alt_bases.Where(r => r.Key >= start).Select(r => new KeyValuePair<int, char>(r.Key + variant_alt_allele.Length, r.Value)).ToArray())
                        .ToDictionary(x => x.Key, x => x.Value);
                }

                string line;
                while (null != (line = subfileReader.ReadLine()))
                {
                    ASETools.SAMLine samLine;

                    int nMatchingRef = 0;
                    int nMatchingAlt = 0;
                    int nMatchingNeither = 0;
                    int nMatchingBoth = 0;

                    try
                    {
                        samLine = new ASETools.SAMLine(line);
                    }
                    catch (FormatException)
                    {
                        Console.WriteLine("Unable to parse sam line in extracted reads main file " + selectedReadsFilename + " subfile " + subfileName + ": " + line);
                        subfileReader.Close();
                        return null;
                    }

                    if (samLine.isUnmapped() || samLine.mapq < 10)
                    {
                        //
                        // Probably half of a paired-end read with the other end mapped.  Ignore it.
                        //
                        continue;
                    }

                    if (contig != samLine.rname)
                    {
                        // Wrong contig. Ignore it.
                        continue;
                    }

                    // make sure variant and cpg can be found on this SamLine
                    if (!samLine.mappedBases.ContainsKey(start) || !samLine.mappedBases.ContainsKey(cpg_start_position))
                    {
                        //
                        // This read does not overlap the variant.  Ignore it.
                        // The 'before' case differs between variant types, and must be dealt with on a casewise basis
                        // below in the switch statement.
                        //
                        continue;
                    }
                    if (!(start > 1 + samLine.mappedBases.Keys.Min()) || !(start < samLine.mappedBases.Keys.Max() - 1))
                    {
                        //
                        // There is no padding on the variant. Without it, its hard to make a match call. Ignore.
                        //
                        continue;
                    }

                    bool matchesRef = false;
                    bool matchesAlt = false;

                    switch (variantType)
                    {
                        case "INS":
                        case "DEL":
                        case "SNP":
                        case "DNP":
                        case "TNP":

                            // match to ref on both sides of variant
                            var refMatchesLeft = ReadCounts.matchSequence(reference_bases, samLine.mappedBases, start, false);
                            var refMatchesRight = ReadCounts.matchSequence(reference_bases, samLine.mappedBases, start, true);

                            // match to alt on both sides of variant
                            var altMatchesLeft = ReadCounts.matchSequence(alt_bases, samLine.mappedBases, start, false);
                            var altMatchesRight = ReadCounts.matchSequence(alt_bases, samLine.mappedBases, start, true);

                            // match criteria: there are some matching bases on both sides, with a total >= 10
                            matchesRef = refMatchesLeft > 0 && refMatchesRight > 0 && refMatchesLeft + refMatchesRight >= 10;
                            matchesAlt = altMatchesLeft > 0 && altMatchesRight > 0 && altMatchesLeft + altMatchesRight >= 10;

                            break;
                        default:
                            Console.WriteLine("Unknown variant type: " + variantType);
                            subfileReader.Close();
                            return null;
                    }

                    // increment matches
                    if (matchesRef)
                    {
                        if (matchesAlt)
                        {
                            nMatchingBoth++;
                        }
                        else
                        {
                            matchingRef.Add(samLine.mappedBases);
                            nMatchingRef++;
                        }
                    }
                    else if (matchesAlt)
                    {
                        matchingAlt.Add(samLine.mappedBases);
                        nMatchingAlt++;
                    }
                    else
                    {
                        nMatchingNeither++;
                    }
                } // foreach SAMLine

                subfileReader.Close();

                // Part 2: take the sequences matching ref and alt and calculate ASM at the cpg
                var nMethylatedRef = 0;
                var nMethylatedAlt = 0;

                var nUnmethylatedRef = 0;
                var nUnmethylatedAlt = 0;

                var nNeitherRef = 0;
                var nNeitherAlt = 0;

                foreach (var read in matchingRef)
                {
                    if (read[cpg_start_position] == 'C')
                    {
                        nMethylatedRef++;
                    }
                    else if (read[cpg_start_position] == 'T')
                    {
                        nUnmethylatedRef++;
                    }
                    else
                    {
                        nNeitherRef++;
                    }
                }

                foreach (var read in matchingAlt)
                {
                    if (read[cpg_start_position] == 'C')
                    {
                        nMethylatedAlt++;
                    }
                    else if (read[cpg_start_position] == 'T')
                    {
                        nUnmethylatedAlt++;
                    }
                    else
                    {
                        nNeitherAlt++;
                    }
                }

                // return methylated (C) and  unmethylated (T)
                return new Tuple<MethylationCounts, MethylationCounts>(new MethylationCounts(nMethylatedRef, nUnmethylatedRef, nNeitherRef),
                    new MethylationCounts(nMethylatedAlt, nUnmethylatedAlt, nNeitherAlt));

            }
        } // MethylationCounts


        public class ReadCounts
        {

            enum MappedReadResult { MatchesRef, MatchesAlt, MatchesNeither, MatchesBoth };

            public double AlleleSpecificValue() // I'm not calling this "allele specific expression" because DNA doesn't have expression.
            {
                if (nMatchingAlt + nMatchingReference == 0)
                {
                    return 0;
                }

                return ((double)Math.Abs(nMatchingReference - nMatchingAlt)) / (nMatchingReference + nMatchingAlt);
            }

            public double AltFraction()
            {
                if (nMatchingAlt + nMatchingReference == 0)
                {
                    return 0.5; // Whatever
                }

                return (double)nMatchingAlt / (nMatchingReference + nMatchingAlt);
            }

            public static ReadCounts ComputeReadCounts(string selectedReadsFilename, string contig,
                int start_position,
                string reference_allele,
                string alt_allele,
                string variantType,
                string subfileName,
                Genome genome)
            {
                int nMatchingRef = 0;
                int nMatchingAlt = 0;
                int nMatchingNeither = 0;
                int nMatchingBoth = 0;

                var consolodatedFile = new ASETools.ConsolodatedFileReader();

                if (!consolodatedFile.open(selectedReadsFilename))
                {
                    Console.WriteLine("Unable to open reads at selected variants file " + selectedReadsFilename);
                    return null;
                }

                var subfileReader = consolodatedFile.getSubfile(subfileName);
                if (null == subfileReader)
                {
                    Console.WriteLine("ComputeReadCounts: no subfile named " + subfileName);
                    return null;
                }

                var padding = 20;
                int[] posArray = new int[padding];
                for (int i = 0; i < padding; i++)
                {
                    posArray[i] = start_position - padding / 2 + i;
                }

                // get reference and alt sequences
                var reference_bases = posArray
                    .Select(r => new Tuple<int, char>(r, genome.getBase(contig, r))).ToDictionary(x => x.Item1, x => x.Item2);
                var alt_bases = posArray
                    .Select(r => new Tuple<int, char>(r, genome.getBase(contig, r))).ToDictionary(x => x.Item1, x => x.Item2);

                // If the variant we are looking at is an INS, we have to shift the variant position 
                // by one to match the SamLine insertion location
                int start = variantType == "INS" ? start_position + 1 : start_position;

                // modify alt_bases based on VariantType
                if (variantType == "DEL")
                {
                    // replace missing bases with "N" so we can do a direct comparison
                    for (var i = 0; i < reference_allele.Length; i++)
                    {
                        alt_bases[start + i] = 'N';
                    }
                }
                else if (variantType != "INS")
                {
                    // for NP: replace bases for alt
                    for (var i = 0; i < alt_allele.Length; i++)
                    {
                        alt_bases[start + i] = alt_allele[i];
                    }
                }
                else
                {
                    // insert INS into dictionary
                    KeyValuePair<int, char>[] ins = new KeyValuePair<int, char>[alt_allele.Length];

                    var chArr = alt_allele.ToCharArray();
                    for (int i = 0; i < chArr.Length; i++)
                    {
                        ins[i] = new KeyValuePair<int, char>(i + start, chArr[i]);
                    }

                    alt_bases =
                        alt_bases.Where(r => r.Key < start).ToArray().Concat(ins.ToArray())
                        .Concat(alt_bases.Where(r => r.Key >= start).Select(r => new KeyValuePair<int, char>(r.Key + alt_allele.Length, r.Value)).ToArray())
                        .ToDictionary(x => x.Key, x => x.Value);
                }

                //
                // Keep track of the reads we've mapped, so if we wind up with both ends of a single read covering the variant
                // that we only count it once (or not at all if it doesn't map the same way).
                //
                var readResults = new Dictionary<string, MappedReadResult>();

                string line;
                while (null != (line = subfileReader.ReadLine()))
                {
                    ASETools.SAMLine samLine;

                    try
                    {
                        samLine = new ASETools.SAMLine(line);
                    }
                    catch (Exception e)
                    {
                        Console.WriteLine("Unable to parse sam line in extracted reads main file + " + selectedReadsFilename + " subfile " + subfileName + ": " + line);
                        throw e;
                    }

                    if (samLine.isUnmapped() || samLine.mapq < 10 /*|| samLine.isDuplicate()*/)
                    {
                        //
                        // Probably half of a paired-end read with the other end mapped.  Ignore it.
                        //
                        continue;
                    }

                    if (contig.ToLower() != samLine.rname.ToLower())
                    {
                        // Wrong contig. Ignore it.
                        continue;
                    }

                    // make sure variant can be found on this SamLine
                    if (!samLine.mappedBases.ContainsKey(start))
                    {
                        //
                        // This read does not overlap the variant.  Ignore it.
                        // The 'before' case differs between variant types, and must be dealt with on a casewise basis
                        // below in the switch statement.
                        //
                        continue;
                    }

                    if (!(start > 1 + samLine.mappedBases.Keys.Min()) || !(start < samLine.mappedBases.Keys.Max() - 1))
                    {
                        //
                        // There is no padding on the variant. Without it, its hard to make a match call. Ignore.
                        //
                        continue;
                    }

                    if (samLine.mappedQual[start] < 10)
                    {
                        //
                        // Too low of a base call quality at the variant site.  Ignore this read.
                        //
                        continue;
                    }

                    bool matchesRef = false;
                    bool matchesAlt = false;

                    switch (variantType)
                    {
                        case "INS":
                        case "DEL":
                        case "SNP":
                        case "DNP":
                        case "TNP":

                            // match to ref on both sides of variant
                            var refMatchesLeft = matchSequence(reference_bases, samLine.mappedBases, start, false);
                            var refMatchesRight = matchSequence(reference_bases, samLine.mappedBases, start, true);

                            // match to alt on both sides of variant
                            var altMatchesLeft = matchSequence(alt_bases, samLine.mappedBases, start, false);
                            var altMatchesRight = matchSequence(alt_bases, samLine.mappedBases, start, true);

                            // match criteria: there are some matching bases on both sides, with a total >= 10
                            matchesRef = refMatchesLeft > 0 && refMatchesRight > 0 && refMatchesLeft + refMatchesRight >= 10;
                            matchesAlt = altMatchesLeft > 0 && altMatchesRight > 0 && altMatchesLeft + altMatchesRight >= 10;

                            break;
                        default:
                            Console.WriteLine("Unknown variant type: " + variantType);
                            subfileReader.Close();
                            return null;
                    }

                    MappedReadResult result;

                    // increment matches
                    if (matchesRef)
                    {
                        if (matchesAlt)
                        {
                            result = MappedReadResult.MatchesBoth;
                        }
                        else
                        {
                            result = MappedReadResult.MatchesRef;
                        }
                    }
                    else if (matchesAlt)
                    {
                        result = MappedReadResult.MatchesAlt;
                    }
                    else
                    {
                        result = MappedReadResult.MatchesNeither;
                    }

                    if (samLine.qname == "*" || !readResults.ContainsKey(samLine.qname))
                    {
                        switch (result)
                        {
                            case MappedReadResult.MatchesRef: nMatchingRef++; break;
                            case MappedReadResult.MatchesAlt: nMatchingAlt++; break;
                            case MappedReadResult.MatchesBoth: nMatchingBoth++; break;
                            case MappedReadResult.MatchesNeither: nMatchingNeither++; break;
                        }

                        if (samLine.qname != "*")
                        {
                            readResults.Add(samLine.qname, result);
                        }
                    } else if (readResults[samLine.qname] != result)
                    {
                        //
                        // The other end of this read also mapped to the same location, and it has a different result.  Just ignore the read.
                        //
                        switch (readResults[samLine.qname])
                        {
                            case MappedReadResult.MatchesRef: nMatchingRef--; break;
                            case MappedReadResult.MatchesAlt: nMatchingAlt--; break;
                            case MappedReadResult.MatchesBoth: nMatchingBoth--; break;
                            case MappedReadResult.MatchesNeither: nMatchingNeither--; break;
                        }
                    }
                } // foreach SAMLine

                subfileReader.Close();
                return new ASETools.ReadCounts(nMatchingRef, nMatchingAlt, nMatchingNeither, nMatchingBoth);
            }

            // recursively matches a sequence and counts number of matches until the sequence ends or a mismatch is found
            public static int matchSequence(Dictionary<int, char> reference, Dictionary<int, char> sequence, int position, bool goRight, int matchCount = 0, int threshold = 10)
            {
                if (matchCount >= threshold)
                {
                    // We have already seen enough bases that match. Return.
                    return matchCount;
                }
                // base case 1: strings have been consumed
                if (!reference.ContainsKey(position) || !sequence.ContainsKey(position))
                {
                    // Reached end of sequence. Return. 
                    return matchCount;
                }

                bool match = reference[position] == sequence[position];
                if (match)
                {
                    if (goRight)
                        position += 1;
                    else
                        position -= 1;

                    matchCount += 1;
                    return matchSequence(reference, sequence, position, goRight, matchCount);
                }
                else
                {
                    // base case 2: hit a mismatch
                    return matchCount;
                }
            }


            public ReadCounts(int nMatchingReference_, int nMatchingAlt_, int nMatchingNeither_, int nMatchingBoth_)
            {
                nMatchingReference = nMatchingReference_;
                nMatchingAlt = nMatchingAlt_;
                nMatchingNeither = nMatchingNeither_;
                nMatchingBoth = nMatchingBoth_;
            }

            public readonly int nMatchingReference;
            public readonly int nMatchingAlt;
            public readonly int nMatchingNeither;
            public readonly int nMatchingBoth;

            public override string ToString()
            {
                return (nMatchingReference.ToString() + '\t' + nMatchingAlt.ToString() + '\t' + nMatchingNeither.ToString() + '\t' + nMatchingBoth.ToString());
            }

            public int totalReads()
            {
                return nMatchingAlt + nMatchingBoth + nMatchingNeither + nMatchingReference;
            }

            public int usefulReads()
            {
                return nMatchingReference + nMatchingAlt;
            }

        } // ReadCounts


        public class RegionalExpressionLine
        {
            public string contig;
            public int contig_offset;
            public int nbases_expressed;                                    // n Bases Expressed
            public int nbases_expressed_baseline;                           // n Bases Expressed With Baseline Expression
            public int nbases_expressed_no_baseline;                        // n Bases Expressed Without Baseline Expression
            public int nreads_baseline;                                     // Total Reads Mapped To Bases With Baseline Expression
            public int nreads_no_baseline;                                  // Total Reads Mapped To Bases Without Baseline Expression
            public int nbases_baseline_other_samples;                       // Count of bases with baseline expression but not in this sample
            public double z;                                                // Total z For Bases With Baseline Expression
            public double min_z;                                            // Min z For Bases With Baseline Expression
            public double max_z;                                            // Max z For Bases With Baseline Expression
            public double mean_z;                                           // Mean z for Bases With Baseline Expression
            public double mean_mu;                                          // Mean mu for Bases with Baseline Expression

            public RegionalExpressionLine(string contig_, int contig_offset_,
                int nbases_expressed_,
                int nbases_expressed_baseline_,
                int nbases_expressed_no_baseline_,
                int nreads_baseline_,
                int nreads_no_baseline_,
                int nbases_baseline_other_samples_,
                double z_,
                double min_z_,
                double max_z_,
                double mean_z_,
                double mean_mu_)
            {
                contig = contig_;
                contig_offset = contig_offset_;
                nbases_expressed = nbases_expressed_;
                nbases_expressed_baseline = nbases_expressed_baseline_;
                nbases_expressed_no_baseline = nbases_expressed_no_baseline_;
                nreads_baseline = nreads_baseline_;
                nreads_no_baseline = nreads_no_baseline_;
                nbases_baseline_other_samples = nbases_baseline_other_samples_;
                z = z_;
                min_z = min_z_;
                max_z = max_z_;
                mean_z = mean_z_;
                mean_mu = mean_mu_;
            }
        } // RegionalExpressionLine

        public class Region
        {
            static List<string> getHeaders()
            {
                List<string> headers = new List<string>();
                // Add headers
                headers.Add("Contig");
                headers.Add("Contig Offset");
                headers.Add("n Bases Expressed");
                headers.Add("n Bases Expressed With Baseline Expression");
                headers.Add("n Bases Expressed Without Baseline Expression");
                headers.Add("Total Reads Mapped To Bases With Baseline Expression");
                headers.Add("Total Reads Mapped To Bases Without Baseline Expression");
                headers.Add("Count of bases with baseline expression but not in this sample");
                headers.Add("Total z For Bases With Baseline Expression");
                headers.Add("Min z For Bases With Baseline Expression");
                headers.Add("Max z For Bases With BaselineExpression");
                headers.Add("Mean z for Bases With Baseline Expression");
                headers.Add("Mean mu for Bases with Baseline Expression");

                return headers;
            }

            public Region(int regionSize_, ASETools.ExpressionFile expression_, Dictionary<string, int> highestOffsetForEachContig_, long nHighQualityMappedNuclearReads_, StreamWriter outputFile_)
            {
                regionSize = regionSize_;
                expression = expression_;
                highestOffsetForEachContig = highestOffsetForEachContig_;
                nHighQualityMappedNuclearReads = nHighQualityMappedNuclearReads_;
                outputFile = outputFile_;
            }

            ASETools.ExpressionFile expression;
            Dictionary<string, int> highestOffsetForEachContig;
            int regionSize;
            long nHighQualityMappedNuclearReads;
            StreamWriter outputFile;
            string currentContig = "";

            int baseOffset = 0;
            int lastBaseSeen = 0;
            long nBasesExpressed = 0;
            long nBasesExpressedWithBaselineExpression = 0;
            long nBasesExpressedWithoutBaselineExpression = 0;
            long totalReadsMappedToBasesWithoutBaselineExpression = 0;
            long totalReadsMappedToBasesWithBaselineExpression = 0;
            long nBasesWithBaselineButNoLocalExpression = 0;
            double totalZForBasesWithBaselineExpression = 0;
            double totalMuForBasesWithBaselineExpression = 0;   // Instead of standard deviations above/below the mean, just means.  This is 0-based (0 expression => 0 mu)
            double minZForBasesWithBaselineExpression = 1000000000;
            double maxZForBasesWithBaselineExpression = -10000000000;

            public static List<RegionalExpressionLine> readFile(string filename)
            {
                var file = ASETools.CreateStreamReaderWithRetry(filename);

                if (null == file)
                {
                    Console.WriteLine("Region.readFile: unable to open file " + filename);
                    return null;
                }

                // Skip first 2 rows
                var headerizedFile = new ASETools.HeaderizedFile<RegionalExpressionLine>(file, false, true, "", getHeaders(), false, false, 2);

                List<RegionalExpressionLine> retVal;
                if (!headerizedFile.ParseFile(parse, out retVal))
                {
                    return null;
                }

                return retVal;
            }

            static RegionalExpressionLine parse(Dictionary<string, int> fieldMappings, string[] fields)
            {
                var contig = fields[fieldMappings["Contig"]];
                var contig_offset = Convert.ToInt32(fields[fieldMappings["Contig Offset"]]);
                var nbases_expressed = Convert.ToInt32(fields[fieldMappings["n Bases Expressed"]]);
                var nbases_expressed_baseline = Convert.ToInt32(fields[fieldMappings["n Bases Expressed With Baseline Expression"]]);
                var nbases_expressed_no_baseline = Convert.ToInt32(fields[fieldMappings["n Bases Expressed Without Baseline Expression"]]);
                var nreads_baseline = Convert.ToInt32(fields[fieldMappings["Total Reads Mapped To Bases With Baseline Expression"]]);
                var nreads_no_baseline = Convert.ToInt32(fields[fieldMappings["Total Reads Mapped To Bases Without Baseline Expression"]]);
                var nbases_baseline_other_samples = Convert.ToInt32(fields[fieldMappings["Count of bases with baseline expression but not in this sample"]]);
                var z = Convert.ToDouble(fields[fieldMappings["Total z For Bases With Baseline Expression"]]);
                var min_z = Convert.ToDouble(fields[fieldMappings["Min z For Bases With Baseline Expression"]]);
                var max_z = Convert.ToDouble(fields[fieldMappings["Max z For Bases With BaselineExpression"]]);
                var mean_z = Convert.ToDouble(fields[fieldMappings["Mean z for Bases With Baseline Expression"]]);
                var mean_mu = Convert.ToDouble(fields[fieldMappings["Mean mu for Bases with Baseline Expression"]]);

                return new RegionalExpressionLine(contig, contig_offset, nbases_expressed, nbases_expressed_baseline, nbases_expressed_no_baseline, nreads_baseline,
                    nreads_no_baseline, nbases_baseline_other_samples, z, min_z, max_z, mean_z, mean_mu);
            }

            void ZeroRegion(int regionBaseOffset)
            {
                baseOffset = regionBaseOffset;
                lastBaseSeen = baseOffset - 1;

                closeRegion();
            }

            public void processBase(string contig, int offset, long mappedReadCount)
            {
                if (contig != currentContig)
                {
                    //
                    // Zero out the rest of the contig.
                    //
                    if (currentContig != "")
                    {
                        closeRegion();

                        for (int i = baseOffset + regionSize; i <= highestOffsetForEachContig[currentContig]; i += regionSize)
                        {
                            ZeroRegion(i);
                        }
                    }

                    baseOffset = offset - offset % regionSize;
                    lastBaseSeen = baseOffset - 1;
                    currentContig = contig;

                }
                else if (baseOffset + regionSize < offset)

                {
                    //
                    // Finish this region, and zero any with no expression.
                    //
                    closeRegion();

                    int baseOffsetForNextRegionWithExpression = offset - offset % regionSize;
                    for (int i = baseOffset + regionSize; i < baseOffsetForNextRegionWithExpression; i += regionSize)
                    {
                        ZeroRegion(i);
                    }

                    baseOffset = offset - offset % regionSize;
                    lastBaseSeen = baseOffset - 1;
                }


                for (int i = lastBaseSeen + 1; i < offset; i++)
                {
                    processSkippedBase(i);
                }

                nBasesExpressed++;
                ASETools.MeanAndStdDev meanAndStdDev;
                if (expression.getValue(contig, offset, out meanAndStdDev))
                {
                    nBasesExpressedWithBaselineExpression++;
                    double z = (((double)mappedReadCount / (double)nHighQualityMappedNuclearReads) - meanAndStdDev.mean) / meanAndStdDev.stddev;

                    totalZForBasesWithBaselineExpression += z;
                    minZForBasesWithBaselineExpression = Math.Min(z, minZForBasesWithBaselineExpression);
                    maxZForBasesWithBaselineExpression = Math.Max(z, maxZForBasesWithBaselineExpression);
                    totalReadsMappedToBasesWithBaselineExpression += mappedReadCount;

                    totalMuForBasesWithBaselineExpression += ((double)mappedReadCount / (double)nHighQualityMappedNuclearReads) / meanAndStdDev.mean;
                }
                else
                {
                    nBasesExpressedWithoutBaselineExpression++;
                    totalReadsMappedToBasesWithoutBaselineExpression += mappedReadCount;
                }

                lastBaseSeen = offset;
            }

            public void closeRegion()
            {
                for (int i = lastBaseSeen + 1; i < baseOffset + regionSize; i++)
                {
                    processSkippedBase(i);
                }

                outputFile.WriteLine(currentContig + "\t" + baseOffset + "\t" + nBasesExpressed + "\t" + nBasesExpressedWithBaselineExpression + "\t" + nBasesExpressedWithoutBaselineExpression + "\t" + totalReadsMappedToBasesWithBaselineExpression + "\t" +
                    totalReadsMappedToBasesWithoutBaselineExpression + "\t" + nBasesWithBaselineButNoLocalExpression + "\t" + totalZForBasesWithBaselineExpression + "\t" + minZForBasesWithBaselineExpression + "\t" + maxZForBasesWithBaselineExpression + "\t" +
                    totalZForBasesWithBaselineExpression / (double)regionSize + "\t" + totalMuForBasesWithBaselineExpression / regionSize);

                nBasesExpressed = 0;
                nBasesExpressedWithBaselineExpression = 0;
                nBasesExpressedWithoutBaselineExpression = 0;
                totalReadsMappedToBasesWithoutBaselineExpression = 0;
                totalReadsMappedToBasesWithBaselineExpression = 0;
                totalZForBasesWithBaselineExpression = 0;
                totalMuForBasesWithBaselineExpression = 0;
                nBasesWithBaselineButNoLocalExpression = 0;
                minZForBasesWithBaselineExpression = 1000000000;
                maxZForBasesWithBaselineExpression = -10000000000;
            }

            void processSkippedBase(int offset)
            {
                ASETools.MeanAndStdDev meanAndStdDev;
                if (expression.getValue(currentContig, offset, out meanAndStdDev))
                {
                    nBasesWithBaselineButNoLocalExpression++;
                    totalZForBasesWithBaselineExpression -= meanAndStdDev.mean / meanAndStdDev.stddev; // It's one mean below the mean: ie. 0 expression
                                                                                                       // No need to update totalMu, since 0 expression adds 0 there.
                }
            }

            static public void printHeader(StreamWriter outputFile)
            {
                outputFile.WriteLine(String.Join("\t", getHeaders()));
            }
        } // Region

        public readonly static string[] NonsenseMediatedDecayCausingVariantClassifications = { "Splice_Site", "Nonsense_Mutation", "Frame_Shift_Del", "Frame_Shift_Ins", "Nonstop_Mutation", "Splice_Region" };

        public class AnnotatedVariant : IComparable
        {
            public AnnotatedVariant() { }

            public string toString()
            {
                var str = ConvertToExcelString(Hugo_symbol) + "\t" + contig + '\t' + locus + '\t' + reference_allele + '\t' + alt_allele + '\t' + variantType + '\t' + variantClassification + '\t' + somaticMutation + '\t' + tumorDNAReadCounts.ToString() + '\t' + normalDNAReadCounts.ToString() + '\t' + tumorRNAReadCounts.ToString();

                if (normalRNAReadCounts != null)
                {
                    str += '\t' + normalRNAReadCounts.ToString();
                }
                else
                {
                    str += "\t\t\t\t";
                }

                return str;
            }

            public string getExtractedReadsExtension()
            {
                if (somaticMutation)
                {
                    return "-" + contig + "-" + Math.Max(1, locus - 200) + "-" + (locus + 10);    // Is this right for indels??
                }
                else
                {
                    return "-" + contig + "-" + locus;
                }
            }

            public bool isSilent()
            {
                return variantClassification == "Silent";
            }


            public bool CausesNonsenseMediatedDecay()
            {
                return NonsenseMediatedDecayCausingVariantClassifications.Contains(variantClassification);
            }

            public double GetTumorAlleleSpecificExpression(ASECorrection aseCorrection = null)
            {
                return GetAlleleSpecificExpression(tumorRNAReadCounts, aseCorrection);
            }

            public double GetNormalAlleleSpecificExpression(ASECorrection aseCorrection = null)
            {
                return GetAlleleSpecificExpression(normalRNAReadCounts, aseCorrection);
            }

            public double GetAlleleSpecificExpression(bool tumor, ASECorrection aseCorrection = null)
            {
                if (tumor) return GetTumorAlleleSpecificExpression(aseCorrection);
                return GetNormalAlleleSpecificExpression(aseCorrection);
            }

            public double GetTumorAltAlleleFraction()
            {
                return GetAltAlleleFraction(tumorRNAReadCounts);
            }

            public double GetNormalAltAlleleFraction()
            {
                return GetAltAlleleFraction(normalRNAReadCounts);
            }

            public double GetAltAlleleFraction(bool tumor)
            {
                if (tumor) return GetTumorAltAlleleFraction();
                return GetNormalAltAlleleFraction();
            }

            private static double GetAlleleSpecificExpression(ReadCounts readCounts, ASECorrection aseCorrection)
            {
                //
                // This used to work by doing:
                // var rnaFractionTumor = GetAltAlleleFraction(readCounts);
                // double rawASE = Math.Abs(rnaFractionTumor * 2.0 - 1.0);
                //
                // But due to the vaugaries of floating point, it would produce slightly different values
                // for the same difference in read counts (i.e., |a - b|) depending on whether
                // a or b was bigger.
                //
                int a = Math.Max(readCounts.nMatchingAlt, readCounts.nMatchingReference);
                int b = Math.Min(readCounts.nMatchingAlt, readCounts.nMatchingReference);

                double rawASE = ((double)(a - b)) / (a + b);

                if (aseCorrection == null)
                {
                    return rawASE;
                } else
                {
                    return aseCorrection.getCorrectedASE(rawASE, readCounts.usefulReads());
                }
            }

            private static double GetAltAlleleFraction(ReadCounts readCounts)
            {
                return (double)readCounts.nMatchingAlt / (readCounts.nMatchingReference + readCounts.nMatchingAlt);
            }

            public AnnotatedVariant(string Hugo_symbol_, bool somaticMutation_, string contig_, int locus_, string reference_allele_, string alt_allele_, string variantType_, string variantClassification_, ReadCounts tumorDNAReadCounts_, ReadCounts tumorRNAReadCounts_, ReadCounts normalDNAReadCounts_, ReadCounts normalRNAReadCounts_)
            {
                Hugo_symbol = Hugo_symbol_;
                somaticMutation = somaticMutation_;
                contig = contig_;
                locus = locus_;
                reference_allele = reference_allele_;
                alt_allele = alt_allele_;
                variantType = variantType_;
                variantClassification = variantClassification_;

                tumorDNAReadCounts = tumorDNAReadCounts_;
                tumorRNAReadCounts = tumorRNAReadCounts_;
                normalDNAReadCounts = normalDNAReadCounts_;
                normalRNAReadCounts = normalRNAReadCounts_;
            }

            static AnnotatedVariant parse(Dictionary<string, int> fieldMappings, string[] fields)
            {
                var Hugo_symbol = ConvertToNonExcelString(fields[fieldMappings["Hugo_symbol"]]);
                var contig = fields[fieldMappings["Chromosome"]];
                var loc = Convert.ToInt32(fields[fieldMappings["Position"]]);
                var Ref = fields[fieldMappings["Ref_allele"]];
                var alt = fields[fieldMappings["Alt_allele"]];
                var variantType = fields[fieldMappings["Variant_type"]];
                var variantClassification = fields[fieldMappings["Variant_classification"]];
                var isSomatic = Convert.ToBoolean(fields[fieldMappings["Somatic"]]);

                var nMatchingTumorReferenceDNA = Convert.ToInt32(fields[fieldMappings["tumorDNAmatchingRef"]]);
                var nMatchingTumorAltDNA = Convert.ToInt32(fields[fieldMappings["tumorDNAmatchingAlt"]]);
                var nMatchingTumorNeitherDNA = Convert.ToInt32(fields[fieldMappings["tumorDNAmatchingNeither"]]);
                var nMatchingTumorBothDNA = Convert.ToInt32(fields[fieldMappings["tumorDNAmatchingBoth"]]);

                var nMatchingNormalReferenceDNA = Convert.ToInt32(fields[fieldMappings["normalDNAmatchingRef"]]);
                var nMatchingNormalAltDNA = Convert.ToInt32(fields[fieldMappings["normalDNAmatchingAlt"]]);
                var nMatchingNormalNeitherDNA = Convert.ToInt32(fields[fieldMappings["normalDNAmatchingNeither"]]);
                var nMatchingNormalBothDNA = Convert.ToInt32(fields[fieldMappings["normalDNAmatchingBoth"]]);

                var nMatchingTumorReferenceRNA = Convert.ToInt32(fields[fieldMappings["tumorRNAmatchingRef"]]);
                var nMatchingTumorAltRNA = Convert.ToInt32(fields[fieldMappings["tumorRNAmatchingAlt"]]);
                var nMatchingTumorNeitherRNA = Convert.ToInt32(fields[fieldMappings["tumorRNAmatchingNeither"]]);
                var nMatchingTumorBothRNA = Convert.ToInt32(fields[fieldMappings["tumorRNAmatchingBoth"]]);


                var normalDNAReadCounts = new ReadCounts(nMatchingNormalReferenceDNA, nMatchingNormalAltDNA, nMatchingNormalNeitherDNA, nMatchingNormalBothDNA);
                var tumorDNAReadCounts = new ReadCounts(nMatchingTumorReferenceDNA, nMatchingTumorAltDNA, nMatchingTumorNeitherDNA, nMatchingTumorBothDNA);
                var tumorRNAReadCounts = new ReadCounts(nMatchingTumorReferenceRNA, nMatchingTumorAltRNA, nMatchingTumorNeitherRNA, nMatchingTumorBothRNA);

                ReadCounts normalRNAReadCounts = null;
                // check if normal RNA is present
                if (fields[fieldMappings["normalRNAmatchingRef"]] != "") {
                    var nMatchingNormalReferenceRNA = Convert.ToInt32(fields[fieldMappings["normalRNAmatchingRef"]]);
                    var nMatchingNormalAltRNA = Convert.ToInt32(fields[fieldMappings["normalRNAmatchingAlt"]]);
                    var nMatchingNormalNeitherRNA = Convert.ToInt32(fields[fieldMappings["normalRNAmatchingNeither"]]);
                    var nMatchingNormalBothRNA = Convert.ToInt32(fields[fieldMappings["normalRNAmatchingBoth"]]);

                    normalRNAReadCounts = new ReadCounts(nMatchingNormalReferenceRNA, nMatchingNormalAltRNA, nMatchingNormalNeitherRNA, nMatchingNormalBothRNA);
                } else
                {
                    normalRNAReadCounts = null;
                }

                return new AnnotatedVariant(Hugo_symbol, isSomatic, contig, loc, Ref, alt, variantType, variantClassification, tumorDNAReadCounts, tumorRNAReadCounts, normalDNAReadCounts, normalRNAReadCounts);
            }

            public static List<AnnotatedVariant> readFile(string filename)
            {
                var file = CreateStreamReaderWithRetry(filename);

                if (null == file)
                {
                    Console.WriteLine("AnnotatedVariant.readFile: unable to open file " + filename);
                    return null;
                }

                var wantedFields = new List<string>();
                wantedFields.Add("Hugo_symbol");
                wantedFields.Add("Chromosome");
                wantedFields.Add("Position");
                wantedFields.Add("Ref_allele");
                wantedFields.Add("Alt_allele");
                wantedFields.Add("Variant_type");
                wantedFields.Add("Variant_classification");
                wantedFields.Add("Somatic");
                wantedFields.Add("tumorDNAmatchingRef");
                wantedFields.Add("tumorDNAmatchingAlt");
                wantedFields.Add("tumorDNAmatchingNeither");
                wantedFields.Add("tumorDNAmatchingBoth");
                wantedFields.Add("normalDNAmatchingRef");
                wantedFields.Add("normalDNAmatchingAlt");
                wantedFields.Add("normalDNAmatchingNeither");
                wantedFields.Add("normalDNAmatchingBoth");
                wantedFields.Add("tumorRNAmatchingRef");
                wantedFields.Add("tumorRNAmatchingAlt");
                wantedFields.Add("tumorRNAmatchingNeither");
                wantedFields.Add("tumorRNAmatchingBoth");
                wantedFields.Add("normalRNAmatchingRef");
                wantedFields.Add("normalRNAmatchingAlt");
                wantedFields.Add("normalRNAmatchingNeither");
                wantedFields.Add("normalRNAmatchingBoth");


                var headerizedFile = new HeaderizedFile<AnnotatedVariant>(file, false, true, "", wantedFields);

                List<AnnotatedVariant> retVal;
                if (!headerizedFile.ParseFile(parse, out retVal))
                {
                    return null;
                }

                return retVal;
            }

            public ReadCounts getReadCount(bool tumor, bool dna)
            {
                if (tumor)
                {
                    if (dna)
                    {
                        return tumorDNAReadCounts;
                    } else
                    {
                        return tumorRNAReadCounts;
                    }
                } else
                {
                    if (dna)
                    {
                        return normalDNAReadCounts;
                    } else
                    {
                        return normalRNAReadCounts;
                    }
                }
            }

            public static void writeFile(string filename, List<AnnotatedVariant> annotatedVariants)
            {
                var file = CreateStreamWriterWithRetry(filename);

                if (null == file)
                {
                    Console.WriteLine("AnnotatedVariant.writeFile: unable to create file " + filename);
                    return;
                }

                // Write header
                file.WriteLine(headerString);

                foreach (var annotatedVariant in annotatedVariants)
                {
                    file.WriteLine(annotatedVariant.toString());
                }
                file.WriteLine("**done**");
                file.Close();

            }

            public static readonly string headerString = "Hugo_symbol\tChromosome\tPosition\tRef_allele\tAlt_allele\tVariant_type\tVariant_classification\tSomatic\t" +
                    "tumorDNAmatchingRef\ttumorDNAmatchingAlt\ttumorDNAmatchingNeither\ttumorDNAmatchingBoth\t" +
                    "normalDNAmatchingRef\tnormalDNAmatchingAlt\tnormalDNAmatchingNeither\tnormalDNAmatchingBoth\t" +
                    "tumorRNAmatchingRef\ttumorRNAmatchingAlt\ttumorRNAmatchingNeither\ttumorRNAmatchingBoth\t" +
                    "normalRNAmatchingRef\tnormalRNAmatchingAlt\tnormalRNAmatchingNeither\tnormalRNAmatchingBoth";

            public readonly string Hugo_symbol;                 // Only present for somatic mutations, otherwise ""
            bool IsASECandidate(bool isTumor, int minRNACoverage, int minDNACoverage) // Has at least 10 tumor DNA & RNA reads, and has variant allele frequency between .4 and .6 in tumor DNA
            {
                string whyNot;
                return IsASECandidate(out whyNot, isTumor, minRNACoverage, minDNACoverage);
            }

            bool IsASECandidate(out string whyNot, bool isTumor, int minRNACoverage, int minDNACoverage) // Has at least 10 tumor DNA & RNA reads, and has variant allele frequency between .4 and .6 in tumor DNA

            {
                // get read counts for either tumor or normal
                var DNAReadCounts = isTumor ? this.tumorDNAReadCounts : this.normalDNAReadCounts;
                var RNAReadCounts = isTumor ? this.tumorRNAReadCounts : this.normalRNAReadCounts;

                // normal readcounts are often not available. Throw exception if this is the case.
                if (DNAReadCounts == null || RNAReadCounts == null)
                {
                    throw new Exception(isTumor ? "tumor " : "normal " + "read counts are null");
                }

                return checkReadCountsForASECandidacy(out whyNot, DNAReadCounts, RNAReadCounts, minDNACoverage, minRNACoverage);
            }

            //
            // A version that rolls in the copy number file check.
            //
            public bool IsASECandidate(bool isTumor, Dictionary<bool, List<CopyNumberVariation>> copyNumber, Configuration configuration, Dictionary<string, ASEMapPerGeneLine> perGeneASEMap, GeneMap geneMap, int inputMinRNACoverage = -1, int inputMinDNACoverage = -1)
            {
                string whyNot;
                return IsASECandidate(out whyNot, isTumor, copyNumber, configuration, perGeneASEMap, geneMap, inputMinRNACoverage, inputMinDNACoverage);
            }

            public bool IsASECandidate(bool isTumor, Dictionary<bool, List<CopyNumberVariation>> copyNumber, CommonData commonData)
            {
                return IsASECandidate(isTumor, copyNumber, commonData.configuration, commonData.perGeneASEMap, commonData.geneMap);
            }

            public bool IsASECandidate(out string whyNot, bool isTumor, Dictionary<bool, List<CopyNumberVariation>> copyNumber, Configuration configuration, Dictionary<string, ASEMapPerGeneLine> perGeneASEMap, GeneMap geneMap, int inputMinRNACoverage = -1, int inputMinDNACoverage = -1)
            {
                if (!isTumor && normalRNAReadCounts == null)
                {
                    whyNot = "No Normal RNA";
                    return false;
                }

                //
                // If provided, the input min coverage numbers override the ones in the configuration.  If not, stick with the configuration.
                //
                int minRNACoverage = (inputMinRNACoverage == -1) ? configuration.minRNAReadCoverage : inputMinRNACoverage;
                int minDNACoverage = (inputMinDNACoverage == -1) ? configuration.minDNAReadCoverage : inputMinDNACoverage;

                if (!IsASECandidate(out whyNot, isTumor, minRNACoverage, minDNACoverage))
                {
                    return false;
                }

                //
                // If this is in a gene that has too much ASE in the matched normals, then it's not an ASE candidate.
                //
                if (perGeneASEMap != null)
                {
                    //
                    // If any genes that contain this locus have ASE > 0.5, then reject this place.
                    //
                    foreach (var gene in geneMap.getGenesMappedTo(contig, locus))
                    {
                        if (perGeneASEMap.ContainsKey(gene.hugoSymbol) && perGeneASEMap[gene.hugoSymbol].sampleData[false].meanASE >= configuration.ASEInNormalAtWhichToExcludeGenes)
                        {
                            whyNot = "Gene with too high mean normal ASE";
                            return false;
                        }
                    }
                }

                if (copyNumber == null || copyNumber[isTumor] == null)
                {
                    whyNot = "ASE Candidate";
                    return true;
                }

                var overlappingCNV = copyNumber[isTumor].Where(cnv =>
                    cnv.OverlapsLocus(ASETools.chromosomeNameToNonChrForm(contig),
                    locus, locus + 1)).Where(r => Math.Abs(r.Segment_Mean) > 1.0).ToList();

                if (overlappingCNV.Count() == 0)
                {
                    whyNot = "ASE Candidate";
                    return true;
                }

                whyNot = "Copy number variant";
                return false;
            }

            public static int CompareByLocus(AnnotatedVariant one, AnnotatedVariant two)
            {
                if (one.contig != two.contig)
                {
                    return string.Compare(one.contig, two.contig);
                }

                if (one.locus < two.locus)
                {
                    return -1;
                }

                if (one.locus == two.locus)
                {
                    return 0;
                }

                return 1;
            }

            public int CompareTo(object peer)
            {
                return CompareByLocus(this, (AnnotatedVariant)peer);
            }

            public string IDH1MutationType()
            {
                if (Hugo_symbol.ToUpper() != "IDH1")
                {
                    throw new Exception("Annotatedariant.IHD1MutationType: this is not an IDH1 mutation, hugo symbol = " + Hugo_symbol);
                }

                if (variantType != "SNP")
                {
                    return "Other";
                }

                return IDH1MutantDescription(locus, alt_allele[0]);
            }

            public readonly bool somaticMutation;
            public readonly string contig;
            public readonly int locus;
            public readonly string reference_allele;
            public readonly string alt_allele;
            public readonly string variantType;
            public readonly string variantClassification;       // Only present for somatic mutations, otherwise ""

            public readonly ReadCounts tumorDNAReadCounts;
            public readonly ReadCounts tumorRNAReadCounts;
            public readonly ReadCounts normalDNAReadCounts;
            public readonly ReadCounts normalRNAReadCounts; // This will be null if there is no normal RNA for this case.
        } // AnnotatedVariant

        public static int PhredToInt(char phred)
        {
            return (int)phred - 33;
        }

        public class SAMLine
        {
            static public List<SAMLine> ReadFromFile(StreamReader inputFile)
            {
                var retVal = new List<SAMLine>();

                string rawLine;
                while (null != (rawLine = inputFile.ReadLine()))
                {
                    retVal.Add(new SAMLine(rawLine));
                }

                return retVal;
            }

            public static SAMLine attemptSAMLine(string rawline)
            {
                try
                {
                    return new SAMLine(rawline);
                } catch (Exception e)
                {
                    Console.WriteLine("Unable to parse ill-formatted SAM line: '" + rawline + ", got exception " + e.Message);
                    return null;
                }
            }

            public SAMLine(string rawline)
            {
                var fields = rawline.Split('\t');

                if (fields.Count() < 11)
                {
                    throw new FormatException();
                }

                qname = fields[0];
                flag = Convert.ToInt32(fields[1]);
                rname = fields[2];
                pos = Convert.ToInt32(fields[3]);
                mapq = Convert.ToInt32(fields[4]);
                cigar = fields[5];
                rnext = fields[6];
                pnext = Convert.ToInt32(fields[7]);
                tlen = Convert.ToInt32(fields[8]);
                seq = fields[9];
                qual = fields[10];
                optionalFields = new string[fields.Count() - 11];
                for (int i = 11; i < fields.Count(); i++)
                {
                    optionalFields[i - 11] = fields[i];
                }

                nonclippedBases = 0;

                //
                // Break down the bases in seq based on their mapped location.  Orindarily, this is one
                // base per location, but in the cases of insertions or deletions (or skipped bases for RNA) it could be zero or more than one.
                //
                int currentPos = pos;
                int offsetInCigarString = 0;
                int offsetInSeq = 0;

                // keep track of cigar elements
                int cigarElementStart = 0;
                int cigarElementLength = 0;

                //
                // A small number of reads are incorrectly formatted with |qual| < |seq|.  Just pad them out with 2's.
                //
                while (qual.Length < seq.Length)
                {
                    qual += "2";
                }

                //try
                {

                    while (offsetInCigarString < cigar.Count())
                    {
                        switch (cigar[offsetInCigarString])
                        {
                            case 'M':
                            case 'I':
                            case '=':
                            case 'X':
                                offsetInCigarString++;  // Consume the M, I, = or X
                                int count = GetNextNumberFromString(cigar.Substring(cigarElementStart, cigarElementLength));

                                if (0 == count)
                                {
                                    throw new FormatException();
                                }

                                for (int i = 0; i < count; i++)
                                {
                                    mappedBases.Add(currentPos, seq[offsetInSeq]);
                                    mappedQual.Add(currentPos, PhredToInt(qual[offsetInSeq]));

                                    currentPos++;
                                    offsetInSeq++;
                                }

                                // reset cigar element information
                                cigarElementLength = 0;
                                cigarElementStart = offsetInCigarString;

                                nonclippedBases += count;

                                break;

                            case 'D':
                                offsetInCigarString++;  // Consume the D
                                count = GetNextNumberFromString(cigar.Substring(cigarElementStart, cigarElementLength));

                                if (0 == count)
                                {
                                    throw new FormatException();
                                }

                                for (int i = 0; i < count; i++)
                                {
                                    mappedBases.Add(currentPos, 'N'); // Add placeholder for DEL
                                    mappedQual.Add(currentPos, 70);     // Just use a high quality for a missing base
                                    currentPos++;
                                }

                                // reset cigar element information
                                cigarElementLength = 0;
                                cigarElementStart = offsetInCigarString;

                                break;

                            case 'N':
                            case 'H':
                                offsetInCigarString++;  // Consume the N or H
                                count = GetNextNumberFromString(cigar.Substring(cigarElementStart, cigarElementLength));

                                // skip region. Reset
                                currentPos += count;

                                cigarElementLength = 0;
                                cigarElementStart = offsetInCigarString;
                                break;

                            case 'S':
                                offsetInCigarString++;  // Consume the S
                                count = GetNextNumberFromString(cigar.Substring(cigarElementStart, cigarElementLength));

                                // skip region. Reset
                                offsetInSeq += count;

                                cigarElementLength = 0;
                                cigarElementStart = offsetInCigarString;
                                break;

                            default:
                                offsetInCigarString++;
                                cigarElementLength++;
                                break;
                        } // switch over character in the CIGAR string
                    } // while over the CIGAR string
                }  /*catch (Exception e) 
                {
                    Console.WriteLine("SAMLine.SAMLine: error parsing cigar string " + cigar + " for read " + rawline);
                    throw e;
                } // catch */


            } // SAMLine

            public bool isUnmapped()
            {
                return (flag & Unmapped) == Unmapped;
            }

            public bool isRC()
            {
                return (flag & RC) == RC;
            }

            public bool isSecondaryAlignment()
            {
                return (flag & SecondaryAligment) == SecondaryAligment;
            }

            public bool isSupplementaryAlignment()
            {
                return (flag & SupplementaryAlignment) == SupplementaryAlignment;
            }

            public bool isDuplicate()
            {
                return (flag & Duplicate) == Duplicate;
            }

            public bool isPaired()
            {
                return (flag & MultipleSegments) == MultipleSegments;
            }

            public bool bothHalvesOfPairMapped()
            {
                return (flag & AllSegmentsProperlyAligned) == AllSegmentsProperlyAligned;
            }

            public bool isNextUnmapped()
            {
                return (flag & NextUnmapped) == NextUnmapped;
            }

#if false
            static public bool operator==(SAMLine one, SAMLine two)
            {
                return one.qname == two.qname && one.flag == two.flag;  // This is sufficient as long as there are only pairs and not more than that with the same name.
            }

            static public bool operator!=(SAMLine one, SAMLine two)
            {
                return !(one == two);   // Why you have to define this and it's not automatic I'll never know.
            }
#endif

            public bool NMKnown()
            {
                return (optionalFields.Any(_ => _.StartsWith("NM:i:")));
            }

            public int NM()
            {
                var NMField = optionalFields.Where(_ => _.StartsWith("NM:i:")).ToList();
                if (NMField.Count != 1)
                {
                    throw new Exception("SAMLine.NM: line does not have exactly one NM:i: field");
                }

                return Convert.ToInt32(NMField[0].Substring(5));
            } // NM

            public readonly string qname;
            public readonly int flag;
            public readonly string rname;
            public readonly int pos;
            public readonly int mapq;
            public readonly string cigar;
            public readonly string rnext;
            public readonly int pnext;
            public readonly int tlen;
            public readonly string seq;
            public readonly string qual;
            public readonly int nonclippedBases;
            public readonly string[] optionalFields;

            public const int MultipleSegments = 0x1;    // i.e., paired
            public const int AllSegmentsProperlyAligned = 0x2;
            public const int Unmapped = 0x4;
            public const int NextUnmapped = 0x8;
            public const int RC = 0x10;
            public const int MateRC = 0x20;
            public const int FirstSegment = 0x40;
            public const int LastSegment = 0x80;
            public const int SecondaryAligment = 0x100;
            public const int Duplicate = 0x400;
            public const int SupplementaryAlignment = 0x800;

            // Dictionary of position and bases at position
            public Dictionary<int, char> mappedBases = new Dictionary<int, char>();
            public Dictionary<int, int> mappedQual = new Dictionary<int, int>();            // base call quality, converted from phred to int

        } // SAMLine

        //
        // Take a string of form ###<otherstuff> and return ### as an int.
        //
        public static int GetNextNumberFromString(string input)
        {
            int nDigits = 0;
            while (nDigits < input.Count() && input[nDigits] >= '0' && input[nDigits] <= '9')
            {
                nDigits++;
            }

            if (0 == nDigits)
            {
                throw new FormatException();
            }

            return Convert.ToInt32(input.Substring(0, nDigits));
        }

        //
        // Generate a pathname that's the equivalent of pathname\..\..\..\.. (where there levelsToGoUp copies of \..) by truncating the top levels
        // rather than adding the \..'s.  Returns a pathname with a trailing backslash.
        //
        public static string GoUpFilesystemLevels(string pathname, int levelsToGoUp)
        {
            string outputPathname = pathname;

            if (outputPathname.EndsWith(@"\"))
            {
                outputPathname = outputPathname.Substring(0, outputPathname.Count() - 1);
            }

            for (int i = 0; i < levelsToGoUp; i++)
            {
                int lastBackslash = outputPathname.LastIndexOf('\\');

                if (lastBackslash <= 0)
                {
                    Console.WriteLine("GoUpFilesystemLevels: not " + levelsToGoUp + " levels to go up in " + pathname);
                    return @"\";
                }

                outputPathname = outputPathname.Substring(0, lastBackslash);
            }

            return outputPathname + @"\";
        }

        public static string GetDataDirectoryFromFilename(string filename, Configuration configuration)
        {
            foreach (var dataDirectory in configuration.dataDirectories)
            {
                if (filename.ToLower().StartsWith(dataDirectory.ToLower()))
                {
                    return dataDirectory;
                }
            }

            return "none";
        }

        public static string GetDerivedFiledDirectoryFromFilename(string filename, Configuration configuration)
        {
            var dataDirectory = GetDataDirectoryFromFilename(filename, configuration);
            if (dataDirectory == "none")
            {
                return dataDirectory;
            }

            return dataDirectory.Substring(0, dataDirectory.LastIndexOf(@"downloaded_files\")) + @"derived_files\";
        }

        public static List<string> GetListOfDiseases(Dictionary<string, Case> cases)
        {
            var diseases = new List<string>();

            foreach (var caseEntry in cases)
            {
                var case_ = caseEntry.Value;

                if (!diseases.Contains(case_.disease()))
                {
                    diseases.Add(case_.disease());
                }
            }

            return diseases;
        }

        public class MappedBaseCount
        {
            public static MappedBaseCount readFromFile(string filename)
            {
                var reader = CreateStreamReaderWithRetry(filename);
                if (null == reader)
                {
                    Console.WriteLine("Unable to read mapped base count file " + filename);
                    return null;
                }

                var line = reader.ReadLine();
                if (null == line)
                {
                    Console.WriteLine("Mapped base count file " + filename + " is empty.");
                    reader.Close();
                    return null;
                }

                long mappedBaseCount, basesCovered;

                var fields = line.Split('\t');
                if (fields.Count() != 3)
                {
                    Console.WriteLine("Format error in mapped base count file " + filename);
                    reader.Close();
                    return null;
                }

                try
                {
                    mappedBaseCount = Convert.ToInt64(fields[0]);
                    basesCovered = Convert.ToInt64(fields[2]);

                } catch (FormatException)
                {
                    Console.WriteLine("Mapped base count file " + filename + " contains an unparsable count: " + fields[0] + " or " + fields[2]);
                    reader.Close();
                    return null;
                }

                line = reader.ReadLine();
                if (line != "**done**")
                {
                    Console.WriteLine("Mapped base count file " + filename + " doens't have **done** where it should.");
                    reader.Close();
                    return null;
                }

                if (null != reader.ReadLine())
                {
                    Console.WriteLine("Mapped base count file " + filename + " continues after **done**");
                    reader.Close();
                    return null;
                }

                reader.Close();
                return new MappedBaseCount(mappedBaseCount, basesCovered);
            }


            MappedBaseCount(long mappedBaseCount_, long basesCovered_)
            {
                mappedBaseCount = mappedBaseCount_;
                basesCovered = basesCovered_;
            }

            public readonly long mappedBaseCount;       // The total length of all of the reads
            public readonly long basesCovered;      // The total number of bases with at least one read mapped to it.
        }

        public class SummaryLine
        {
            static public List<SummaryLine> ReadFromFile(string filename)
            {
                List<SummaryLine> retVal;

                string[] wantedFields =
                {
                    "Gene",
                    "nSingle",
                    "nMultiple",
                    "nInteresting",
                    "%Interesting",
                    "sex",
                    "median single",
                    "median multi",
                    "median combined",
                    "median heterozygous",
                };

                var reader = CreateStreamReaderWithRetry(filename);
                if (null == reader)
                {
                    Console.WriteLine("Unable to open file " + filename);
                    return null;
                }

                var headerizedFile = new HeaderizedFile<SummaryLine>(reader, false, true, "", wantedFields.ToList());

                if (!headerizedFile.ParseFile(parser, out retVal))
                {
                    Console.WriteLine("Failed to parse " + filename);
                    return null;
                }

                return retVal;
            }

            public SummaryLine(string gene_, int nSingle_, int nMultiple_, int nInteresting_, double percentInteresting_, bool sex_, double medianSingle_, double medianMulti_, double medianCombined_, double medianHeterozygous_)
            {
                gene = gene_;
                nSingle = nSingle_;
                nMultiple = nMultiple_;
                nInteresting = nInteresting_;
                percentInteresting = percentInteresting_;
                sex = sex_;
                medianSingle = medianSingle_;
                medianMulti = medianMulti_;
                medianCombined = medianCombined_;
                medianHeterozygous = medianHeterozygous_;
            }

            static SummaryLine parser(ASETools.HeaderizedFile<SummaryLine>.FieldGrabber fields)
            {
                return new SummaryLine(
                    fields.AsString("Gene"),
                    fields.AsInt("nSingle"),
                    fields.AsInt("nMultiple"),
                    fields.AsInt("nInteresting"),
                    fields.AsDouble("%Interesting"),
                    fields.AsBool("sex"),
                    fields.AsDouble("median single"),
                    fields.AsDouble("median multi"),
                    fields.AsDouble("median combined"),
                    fields.AsDouble("median heterozygous")
                    );
            }

            public readonly string gene;
            public readonly int nSingle;
            public readonly int nMultiple;
            public readonly int nInteresting;
            public readonly double percentInteresting;
            public readonly bool sex;
            public readonly double medianSingle;
            public readonly double medianMulti;
            public readonly double medianCombined;
            public readonly double medianHeterozygous;
        } // SummaryLine

        public class SelectedVariantCountByGene
        {
            SelectedVariantCountByGene(string Hugo_Symbol_, int nVariantsInGene_, int nSomaticMutations_)
            {
                Hugo_Symbol = Hugo_Symbol_;
                nVariantsInGene = nVariantsInGene_;
                nSomaticMutations = nSomaticMutations_;
            }

            static public List<SelectedVariantCountByGene> LoadFromFile(string filename)
            {
                List<SelectedVariantCountByGene> retVal;
                string[] wantedFields =
                {
                    "Hugo_Symbol",
                    "nVariantsInGene",
                    "nSomaticMutations"
                };

                var reader = CreateStreamReaderWithRetry(filename);
                if (reader == null)
                {
                    Console.WriteLine("Unable to open " + filename);
                    return null;
                }

                var headerizedFile = new HeaderizedFile<SelectedVariantCountByGene>(reader, false, true, "", wantedFields.ToList());

                if (!headerizedFile.ParseFile(parser, out retVal))
                {
                    Console.WriteLine("Unable to parse " + filename);
                    return null;
                }

                return retVal;
            }

            static SelectedVariantCountByGene parser(ASETools.HeaderizedFile<SelectedVariantCountByGene>.FieldGrabber fields)
            {
                return new SelectedVariantCountByGene(fields.AsString("Hugo_Symbol"), fields.AsInt("nVariantsInGene"), fields.AsInt("nSomaticMutations"));
            }

            public readonly string Hugo_Symbol;
            public readonly int nVariantsInGene;
            public readonly int nSomaticMutations;
        } // SelectedVariantCountByGene

        //
        // Just a cute thing for classifying tumors that you can use as an array index.
        //
        public static int ZeroOneMany(int value)
        {
            if (value == 0) return value;
            if (value == 1) return value;
            return 2;
        }

        public static string ZeroOneManyString(int value)
        {
            switch (ZeroOneMany(value))
            {
                case 0: return "0";
                case 1: return "1";
                case 2: return ">1";
            }

            return "DANGER, WILL ROBINSON!  ASETools.ZeroOneMany() returned an unexpected value.";    // Just to make it compile
        }


        public class BisulfateCase
        {
            public string case_id;


            public string bam_tumor_filename;
            public string bed_tumor_filename;
            public string bed_normal_filename;
            public string bam_normal_filename;

            public string bam_tumor_file_id;
            public string bed_tumor_file_id;
            public string bed_normal_file_id;
            public string bam_normal_file_id;

            public string realigned_maf_filename;
            public string realigned_selected_variants_filename;

            public BisulfateCase()
            {
            }

            public static Dictionary<string, BisulfateCase> loadCases(string filename)
            {

                StreamReader inputFile = ASETools.CreateStreamReaderWithRetry(filename);

                // format filename, location, file id, data format, access, data category, file size, project id, case id
                var header = inputFile.ReadLine();

                string line;

                Dictionary<string, BisulfateCase> cases = new Dictionary<string, BisulfateCase>();

                while ((line = inputFile.ReadLine()) != null)
                {
                    var split = line.Split('\t');

                    // build bisulfate case
                    var case_id = split[0];


                    if (!cases.ContainsKey(case_id))
                    {
                        cases.Add(case_id, new BisulfateCase());
                    }
                    else
                    {
                        throw new Exception("Case " + case_id + " already in file");
                    }

                    cases[case_id].case_id = case_id;
                    cases[case_id].realigned_maf_filename = split[1];
                    cases[case_id].realigned_selected_variants_filename = split[2];
                    cases[case_id].bam_tumor_filename = split[3];
                    cases[case_id].bam_tumor_file_id = split[4];
                    cases[case_id].bam_normal_filename = split[5];
                    cases[case_id].bam_normal_file_id = split[6];

                    cases[case_id].bed_tumor_filename = split[7];
                    cases[case_id].bed_tumor_file_id = split[8];
                    cases[case_id].bed_normal_file_id = split[9];
                    cases[case_id].bed_normal_file_id = split[10];
                }

                return cases;
            }
        }


        public class MethylationPoint : IComparer<MethylationPoint>
        {
            public string identifier;
            public double adjustedBValue;
            public double mValue;
            public bool hasOneMutation;

            public MethylationPoint(string identifier_, double mValue_, bool hasOneMutation_)
            {
                mValue = mValue_;
                identifier = identifier_;
                hasOneMutation = hasOneMutation_;

                var bValue = MethylationAnnotationLine.M2Beta(mValue);
                // Values of 1 have no/full methylation. values of partial methylation have 0 score
                adjustedBValue = Math.Abs(2 * bValue - 1);
            }

            public int Compare(MethylationPoint a, MethylationPoint b)
            {
                return xCompare(a, b);
            }

            static public int xCompare(MethylationPoint a, MethylationPoint b)
            {
                if (a.adjustedBValue > b.adjustedBValue) return 1;
                if (a.adjustedBValue < b.adjustedBValue) return -1;
                return 0;
            }
        }

        public static bool[] BothBools = { true, false };   // There's probably something like this in the runtime, but whatever.
        public static Dictionary<bool, string> tumorToString = new Dictionary<bool, string>();  // true -> Tumor, false -> Normal
        public static Dictionary<bool, string> dnaToString = new Dictionary<bool, string>();
        public static string[] BothGenders = { "male", "female" };

        public class HistogramResultLine
        {
            public string minValue;
            public double minValueAsDouble;
            public long count = 0;
            public double total = 0;    // The sum of all of the values in this line
            public double pdfValue = 0;
            public double cdfValue = 0;

            public override string ToString()
            {
                return minValue + "\t" + count + "\t" + total + "\t" + pdfValue + "\t" + cdfValue;
            }

            public static string Header()
            {
                return "minValue\tcount\ttotal\tpdf\tcdf";
            }

            static public HistogramResultLine[] ReadFromFile(string filename)
            {
                var file = CreateStreamReaderWithRetry(filename);
                if (null == file)
                {
                    Console.WriteLine("ASETools.Histogram.ReadFromFile: unable to open " + filename);
                    return null;
                }

                var retVal = ReadFromStream(file, true);
                file.Close();

                return retVal;
            }

            static public HistogramResultLine[] ReadFromStream(StreamReader file, bool wholeFile)
            {

                string[] wantedFields = { "minValue", "count", "total", "pdf", "cdf" };

                var headerizedFile = new HeaderizedFile<HistogramResultLine>(file, false, wholeFile, "", wantedFields.ToList(), false, false, 0, false, '\t', !wholeFile); // If not whole file, then don't expect **done** and stop at the first blank line

                List<HistogramResultLine> lines;

                if (!headerizedFile.ParseFile(parse, out lines))
                {
                    Console.WriteLine("ASETools.Histogram.ReadFromFile: unable to parse input file.");
                    file.Close();
                    return null;
                }


                return lines.ToArray();
            }

            static HistogramResultLine parse(HeaderizedFile<HistogramResultLine>.FieldGrabber fieldGrabber)
            {
                double minValue;
                if (fieldGrabber.AsString("minValue") == "More")
                {
                    minValue = 1;
                } else
                {
                    minValue = fieldGrabber.AsDouble("minValue");
                }
                return new HistogramResultLine(fieldGrabber.AsString("minValue"), minValue, fieldGrabber.AsInt("count"), fieldGrabber.AsDouble("total"), fieldGrabber.AsDouble("pdf"), fieldGrabber.AsDouble("cdf"));
            }

            public HistogramResultLine(string minValue_, double minValueAsDouble_, int count_, double total_, double pdfValue_, double cdfValue_)
            {
                minValue = minValue_;
                minValueAsDouble = minValueAsDouble_;
                count = count_;
                total = total_;
                pdfValue = pdfValue_;
                cdfValue = cdfValue_;
            }

            public HistogramResultLine() { }
        } // HistogramResultLine

        public class Histogram
        {
            public Histogram(string name_ = "")
            {
                name = name_;
            }

            public void addValue(double value)
            {
                values.Add(value);
            }

            public double min() // Throws ArgumentNullExcpetion if no data jas been added
            {
                return values.Min();
            }

            public double max()
            {
                return values.Max();
            }

            public double mean()
            {
                return values.Average();
            }

            public double median()
            {
                values.Sort();
                if (values.Count() % 2 == 0)
                {
                    return (values[values.Count() / 2] + values[values.Count() / 2 - 1]) / 2;   // Take the mean of the middle values for an even-sized distribution.
                } else
                {
                    return values[values.Count() / 2];
                }
            }

            public void merge(Histogram peer)
            {
                values.AddRange(peer.values);
            }
            public int count()
            {
                return values.Count();
            }
            public HistogramResultLine[] ComputeHistogram(double min, double max, double increment, string format = "G")
            {
                int nBuckets = (int)((max - min) / increment) + 1;  // +1 is for "more"

                var result = new HistogramResultLine[nBuckets];

                double x = min;
                for (int i = 0; i < nBuckets - 1; i++)
                {
                    result[i] = new HistogramResultLine();
                    result[i].minValue = x.ToString(format);
                    x += increment;
                }
                result[nBuckets - 1] = new HistogramResultLine();
                result[nBuckets - 1].minValue = "More";

                foreach (var value in values)
                {
                    int whichBucket;
                    if (value > max)
                    {
                        whichBucket = nBuckets - 1;
                    }
                    else
                    {
                        whichBucket = (int)((value - min) / increment);
                    }

                    result[whichBucket].count++;
                    result[whichBucket].total += value;
                }

                int overallCount = values.Count();
                long runningCount = 0;
                for (int whichBucket = 0; whichBucket < nBuckets; whichBucket++)
                {
                    runningCount += result[whichBucket].count;

                    result[whichBucket].pdfValue = ((double)result[whichBucket].count) / overallCount;
                    result[whichBucket].cdfValue = ((double)runningCount) / overallCount;
                }

                return result;
            } // ComputeHistogram

            public static double[] ReadStreamToPDFValues(StreamReader inputStream, bool wholeFile)
            {
                var lines = HistogramResultLine.ReadFromStream(inputStream, wholeFile);
                if (null == lines) return null;

                return lines.Select(x => x.pdfValue).ToArray();
            }

            public static double[] ReadFileToPDFValues(string filename)
            {
                var lines = HistogramResultLine.ReadFromFile(filename);
                if (null == lines) return null;

                return lines.Select(x => x.pdfValue).ToArray();
            }

            public static double[] ReadStreamToCDFValues(StreamReader inputStream, bool wholeFile)
            {
                var lines = HistogramResultLine.ReadFromStream(inputStream, wholeFile);
                if (null == lines) return null;

                return lines.Select(x => x.cdfValue).ToArray();
            }

            public static double[] ReadFileToCDFValues(string filename)
            {
                var lines = HistogramResultLine.ReadFromFile(filename);
                if (null == lines) return null;

                return lines.Select(x => x.cdfValue).ToArray();
            }

            public List<double> getValues()
            {
                return values;
            }

            List<double> values = new List<double>();
            public readonly string name;
        } // Histogram

        public class PreBucketedHistogram
        {
            double minBucket;
            double maxBucket;
            double increment;
            int nBuckets;
            HistogramResultLine[] buckets;
            public readonly string name;

            double minSeenValue = double.MaxValue;
            double maxSeenValue = double.MinValue;
            double totalSeenValue = 0;
            long nValues = 0;

            public PreBucketedHistogram(double minBucket_, double maxBucket_, double increment_, string name_ = "")
            {
                minBucket = minBucket_;
                maxBucket = maxBucket_;
                increment = increment_;
                nBuckets = (int)((maxBucket - minBucket) / increment) + 1;  // +1 is for "more"
                buckets = new HistogramResultLine[nBuckets];
                for (int i = 0; i < nBuckets; i++)
                {
                    buckets[i] = new HistogramResultLine();
                }

                name = name_;
            }

            public void addValue(double value)
            {
                if (value < minBucket || double.IsNaN(value))
                {
                    throw new FormatException("PreBucketedHistogram " + name + " : value (" + value + ") smaller than minBucket (" + minBucket + ") or NaN");
                }

                int whichBucket = bucketForValue(value);

                /*BJB*/
                if (whichBucket >= nBuckets || whichBucket < 0)
                {
                    Console.WriteLine(whichBucket + " >= " + nBuckets + " value = " + value);
                }

                buckets[whichBucket].count++;
                buckets[whichBucket].total += value;

                minSeenValue = Math.Min(minSeenValue, value);
                maxSeenValue = Math.Max(maxSeenValue, value);
                totalSeenValue += value;
                nValues++;
            }

            public double min()
            {
                return minSeenValue;
            }

            public double max()
            {
                return maxSeenValue;
            }

            public double mean()
            {
                return totalSeenValue / nValues;
            }

            public double percentile(int desiredPercentile)
            {
                if (desiredPercentile < 1 || desiredPercentile > 99)
                {
                    throw new Exception("PreBuckedHistogram: percentiles must be between 1 and 99 inclusive");
                }

                var lines = ComputeHistogram();

                double desiredPercentileAsDouble = (double)desiredPercentile / 100;

                for (int i = 0; i < lines.Count() - 1; i++)
                {
                    if (lines[i].cdfValue >= desiredPercentileAsDouble)
                    {
                        if (i == 0)
                        {
                            return lines[i].minValueAsDouble;
                        }
                    }
                }

                /*BJB*/
                return 0;   // writeme
            }

            public void merge(PreBucketedHistogram peer)
            {
                if (peer.minBucket != minBucket || peer.maxBucket != maxBucket || peer.increment != increment || peer.nBuckets != nBuckets)
                {
                    throw new FormatException("Can't merge pre-bucketed histograms with different parameters");
                }

                for (int i = 0; i < nBuckets; i++)
                {
                    buckets[i].count += peer.buckets[i].count;
                    buckets[i].total += peer.buckets[i].total;
                }

                minSeenValue = Math.Min(minSeenValue, peer.minSeenValue);
                maxSeenValue = Math.Max(maxSeenValue, peer.maxSeenValue);
                totalSeenValue += peer.totalSeenValue;
                nValues += peer.nValues;
            }

            static public PreBucketedHistogram ReadFromSteam(StreamReader inputFile, string name, bool wholeFile)
            {
                var resultLines = HistogramResultLine.ReadFromStream(inputFile, wholeFile);

                if (resultLines.Count() < 2 || resultLines[resultLines.Count() - 1].minValue != "More")
                {
                    throw new FormatException("ASELib.PreBucketedHistogram.ReadFromStream: found too few results lines, or the last one isn't 'More'");
                }

                int nBuckets = resultLines.Count() - 1;
                double minValue = resultLines[0].minValueAsDouble;
                double maxValue = resultLines[nBuckets - 1].minValueAsDouble;
                double increment = (maxValue - minValue) / (nBuckets - 1);

                var retVal = new PreBucketedHistogram(minValue, maxValue, increment, name);

                //
                // Now fill in the actual values in the buckets.
                //
                for (int i = 0; i < nBuckets; i++)
                {
                    retVal.buckets[i].count = resultLines[i].count;
                    retVal.buckets[i].total = resultLines[i].total;
                    retVal.totalSeenValue += resultLines[i].total;
                    retVal.nValues += resultLines[i].count;
                }

                return retVal;
            }

            public long count()
            {
                return nValues;

            }

            int bucketForValue(double value)
            {
                if (value > maxBucket)
                {
                    return nBuckets - 1;
                }
                else
                {
                    return (int)((value - minBucket) / increment);
                }
            }

            public long nLessThan(double value)
            {
                long retVal = 0;

                for (int i = 0; i < bucketForValue(value); i++)
                {
                    retVal += buckets[i].count;
                }

                return retVal;
            }

            public HistogramResultLine[] ComputeHistogram(string format = "G")
            {
                var result = new HistogramResultLine[nBuckets];

                double x = minBucket;
                long runningCount = 0;
                for (int i = 0; i < nBuckets; i++)
                {
                    result[i] = new HistogramResultLine();
                    result[i].minValue = x.ToString(format);
                    result[i].count = buckets[i].count;
                    result[i].total = buckets[i].total;
                    result[i].pdfValue = ((double)result[i].count) / nValues;
                    runningCount += result[i].count;
                    result[i].cdfValue = ((double)runningCount) / nValues;
                    x += increment;
                }

                result[nBuckets - 1].minValue = "More";

                return result;
            }// ComputeHistogram

            public void WriteHistogram(StreamWriter outputFile, string format = "G")
            {
                outputFile.WriteLine(HistogramResultLine.Header());
                ComputeHistogram(format).ToList().ForEach(_ => outputFile.WriteLine(_));
            }

            public delegate double ValueGetter(HistogramResultLine h);
            static public void WriteBatchOfHistogramValues(StreamWriter outputFile, List<KeyValuePair<string, PreBucketedHistogram>> headersAndHistograms, ValueGetter valueGetter)
            {
                if (headersAndHistograms.Count() == 0)
                {
                    return; // No histograms, no output.
                }

                var first = headersAndHistograms[0].Value;
                if (headersAndHistograms.Select(_ => _.Value).Any(_ => _.minBucket != first.minBucket || _.maxBucket != first.maxBucket || _.increment != first.increment))
                {
                    throw new Exception("ASETools.PreBucketedHistogram.WriteBatchOfHistogramCDFs: not all input histograms have the same shape.");
                }

                outputFile.Write("minValue");
                foreach (var headerAndHistogram in headersAndHistograms)
                {
                    outputFile.Write("\t" + headerAndHistogram.Key + " (n = " + headerAndHistogram.Value.count() + ")");
                }
                outputFile.WriteLine();

                var lines = new HistogramResultLine[headersAndHistograms.Count][];
                for (int i = 0; i < headersAndHistograms.Count(); i++)
                {
                    lines[i] = headersAndHistograms[i].Value.ComputeHistogram();
                }

                for (int whichLine = 0; whichLine < lines[0].Count(); whichLine++)
                {
                    outputFile.Write(lines[0][whichLine].minValue);

                    for (int whichHistogram = 0; whichHistogram < headersAndHistograms.Count(); whichHistogram++)
                    {
                        outputFile.Write("\t" + valueGetter(lines[whichHistogram][whichLine]));
                    } // which histogram
                    outputFile.WriteLine();
                } // Which line

            } // WriteBatchOfHistogramValues

            static public void WriteBatchOfHistogramCDFs(StreamWriter outputFile, List<KeyValuePair<string, PreBucketedHistogram>> headersAndHistograms)
            {
                WriteBatchOfHistogramValues(outputFile, headersAndHistograms, line => line.cdfValue);
            } // WriteBatchOfHistogramCDFs

            static public void WriteBatchOfHistogramPDFs(StreamWriter outputFile, List<KeyValuePair<string, PreBucketedHistogram>> headersAndHistograms)
            {
                WriteBatchOfHistogramValues(outputFile, headersAndHistograms, line => line.pdfValue);
            } // WriteBat

        } // PreBucketedHistogram

        public class NMeandAndStdDev
        {
            public NMeandAndStdDev(int n_, double mean_, double stdDev_)
            {
                n = n_;
                mean = mean_;
                stdDev = stdDev_;
            }

            public readonly int n;
            public readonly double mean;
            public readonly double stdDev;

            public static List<string> getHeaders(string prefix)
            {
                string[] headers = { prefix + " N", prefix + " mean", prefix + " stDdev" };
                return headers.ToList();
            }

            public static string getHeaderString(string prefix)
            {
                return prefix + "N\t" + prefix + " mean\t" + prefix + " stdDev";
            }

            public string outputString()
            {
                return (n == -1 ? "*" : Convert.ToString(n)) + "\t" + (mean == double.NegativeInfinity ? "*" : Convert.ToString(mean)) + "\t" + (stdDev == double.NegativeInfinity ? "*" : Convert.ToString(stdDev));
            }
            public bool valid()
            {
                return n != -1 && mean != double.NegativeInfinity && stdDev != double.NegativeInfinity;
            }
        }

        public class SingleExpressionResult
        {
            public readonly int rangeIndex;     // i.e., 0 -> nRanges - 1
            public readonly int rangeInBases;   // int.MaxValue for whole autosome

            //
            // p values.  double.NegativeInfinity means no value.
            //
            public double oneVsMany;
            public double oneVsNotOne;
            public double oneSilentVsOneNotSilent;

            public readonly NMeandAndStdDev zeroMutationStats;
            public readonly NMeandAndStdDev oneMutationStats;
            public readonly NMeandAndStdDev moreThanOneMutationStats;
            public readonly NMeandAndStdDev onlyOneSilentMutationStats;

            public static List<string> getHeaders(bool mu, bool exclusive, string range)
            {
                //
                // This is pretty much a match for the code that generates the headers in the first place, which results in strange double spaces in some places.
                //
                var result = new List<string>();
                string muString = ((!mu) ? "" : " mu") + (!exclusive ? "" : " exclusive");
                result.Add(range + " 1 vs. many" + muString);
                result.Add(range + " 1 vs. not 1" + muString);
                result.Add(range + " 1 silent vs. 1 not silent" + muString);

                string[] mutationSets = { "0", "1", ">1", "1 Silent" };

                for (int i = 0; i < mutationSets.Count(); i++)
                {
                    result.Add(range + " " + mutationSets[i] + " mutation" + muString + " N");
                    result.Add(range + " " + mutationSets[i] + " mutation " + muString + " mean");
                    result.Add(range + " " + mutationSets[i] + " mutation " + muString + " stdDev");
                }

                return result;
            }

            public static string getHeaderString(bool exclusive, int rangeIndex)
            {
                var regionString = regionIndexToString(rangeIndex) + (exclusive ? " exclusive " : " ");
                return regionString + "1 vs. many\t" +
                       regionString + "1 vs. not 1\t" +
                       regionString + "1 silent vs. 1 not silent\t" +
                       NMeandAndStdDev.getHeaderString(regionString + "0 mutation ") + "\t" +
                       NMeandAndStdDev.getHeaderString(regionString + "1 mutation ") + "\t" +
                       NMeandAndStdDev.getHeaderString(regionString + ">1 mutation ") + "\t" +
                       NMeandAndStdDev.getHeaderString(regionString + "1 Silent");
            }

            public SingleExpressionResult(int rangeIndex_, double oneVsMany_, double oneVsNotOne_, double oneSilentVsOneNotSilent_, NMeandAndStdDev zeroMutationStats_, NMeandAndStdDev oneMutationStats_, NMeandAndStdDev moreThanOneMutationStats_, NMeandAndStdDev onlyOneSilentMutationStats_)
            {
                rangeIndex = rangeIndex_;
                if (rangeIndex == 0)
                {
                    rangeInBases = 0;
                } else if (rangeIndex == nRegions - 1)
                {
                    rangeInBases = int.MaxValue;
                } else
                {
                    rangeInBases = (1 << rangeIndex - 1) * 1000;
                }

                oneVsMany = oneVsMany_;
                oneVsNotOne = oneVsNotOne_;
                oneSilentVsOneNotSilent = oneSilentVsOneNotSilent_;
                zeroMutationStats = zeroMutationStats_;
                oneMutationStats = oneMutationStats_;
                moreThanOneMutationStats = moreThanOneMutationStats_;
                onlyOneSilentMutationStats = onlyOneSilentMutationStats_;
            }

            public static SingleExpressionResult Parse(int rangeIndex, HeaderizedFile<ExpressionResultsLine>.FieldGrabber fieldGrabber, bool mu, bool exclusive)
            {
                //
                // This is pretty much a match for the code that generates the headers in the first place, which results in strange double spaces in some places.
                //
                string rangeString = regionIndexToString(rangeIndex);
                string muString = ((!mu) ? "" : " mu") + (!exclusive ? "" : " exclusive");

                return new SingleExpressionResult(
                    rangeIndex,
                    fieldGrabber.AsDoubleNegativeInfinityIfStarOrEmptyString(rangeString + " 1 vs. many" + muString),
                    fieldGrabber.AsDoubleNegativeInfinityIfStarOrEmptyString(rangeString + " 1 vs. not 1" + muString),
                    fieldGrabber.AsDoubleNegativeInfinityIfStarOrEmptyString(rangeString + " 1 silent vs. 1 not silent" + muString),
                    new NMeandAndStdDev(fieldGrabber.AsIntMinusOneIfStarOrEmptyString(rangeString + " 0 mutation" + muString + " N"), fieldGrabber.AsDoubleNegativeInfinityIfStarOrEmptyString(rangeString + " 0 mutation " + muString + " mean"), fieldGrabber.AsDoubleNegativeInfinityIfStarOrEmptyString(rangeString + " 0 mutation " + muString + " stdDev")),
                    new NMeandAndStdDev(fieldGrabber.AsIntMinusOneIfStarOrEmptyString(rangeString + " 1 mutation" + muString + " N"), fieldGrabber.AsDoubleNegativeInfinityIfStarOrEmptyString(rangeString + " 1 mutation " + muString + " mean"), fieldGrabber.AsDoubleNegativeInfinityIfStarOrEmptyString(rangeString + " 1 mutation " + muString + " stdDev")),
                    new NMeandAndStdDev(fieldGrabber.AsIntMinusOneIfStarOrEmptyString(rangeString + " >1 mutation" + muString + " N"), fieldGrabber.AsDoubleNegativeInfinityIfStarOrEmptyString(rangeString + " >1 mutation " + muString + " mean"), fieldGrabber.AsDoubleNegativeInfinityIfStarOrEmptyString(rangeString + " >1 mutation " + muString + " stdDev")),
                    new NMeandAndStdDev(fieldGrabber.AsIntMinusOneIfStarOrEmptyString(rangeString + " 1 Silent mutation" + muString + " N"), fieldGrabber.AsDoubleNegativeInfinityIfStarOrEmptyString(rangeString + " 1 Silent mutation " + muString + " mean"), fieldGrabber.AsDoubleNegativeInfinityIfStarOrEmptyString(rangeString + " 1 Silent mutation " + muString + " stdDev"))
                    );
            }

            public void appendToString(ref string outputLine)
            {
                outputLine += ASETools.stringForDouble(oneVsMany) + "\t" + ASETools.stringForDouble(oneVsNotOne) + "\t" + ASETools.stringForDouble(oneSilentVsOneNotSilent) + "\t" +
                    zeroMutationStats.outputString() + "\t" + oneMutationStats.outputString() + "\t" + moreThanOneMutationStats.outputString() + "\t" + onlyOneSilentMutationStats.outputString();
            }

            public bool anyValid()
            {
                return oneVsMany != double.NegativeInfinity || oneVsNotOne != double.NegativeInfinity || oneSilentVsOneNotSilent != double.NegativeInfinity || zeroMutationStats.valid() || oneMutationStats.valid() || moreThanOneMutationStats.valid() || onlyOneSilentMutationStats.valid();
            }

            public string nMutationMean(int n)  // n is 0, 1, 2 (which means >=2) or 3 (which means 1 silent)
            {
                NMeandAndStdDev stats = (n == 0) ? zeroMutationStats : ((n == 1) ? oneMutationStats : ((n == 2) ? moreThanOneMutationStats : onlyOneSilentMutationStats));

                if (stats.n < 10 || stats.mean == double.NegativeInfinity)
                {
                    return "";      // Empty string rather than the usual *, because the * just confuses Excel when it's trying to make graphs
                }

                return Convert.ToString(stats.mean);
            }
        } // SingleExpressionResult

        public class ExpressionResultsLine
        {
            static public List<ExpressionResultsLine> readFile(string filename)
            {
                string[] wantedFieldsOtherThanExpression =
                {
                    "Hugo Symbol",
                    "nTumorsExcluded",
                    "nZero",
                    "nOne",
                    "nMore",
                    "nSingleContibutingToAltFraction",
                    "nMultipleContributingToAltFraction",
                    "single alt fraction",
                    "multiple alt fraction"
                };

                var wantedFields = wantedFieldsOtherThanExpression.ToList();
                for (var i = 0; i < nRegions; i++)
                {
                    foreach (bool exclusive in BothBools) {
                        wantedFields.AddRange(SingleExpressionResult.getHeaders(false, exclusive, "" + regionIndexToString(i)));
                    }
                }

                var inputFile = CreateStreamReaderWithRetry(filename);
                if (null == inputFile)
                {
                    return null;
                }
                var headerizedFile = new HeaderizedFile<ExpressionResultsLine>(inputFile, false, false, "", wantedFields, false, true);  // We allow missing columns because genes with no data are truncated

                List<ExpressionResultsLine> retVal;
                headerizedFile.ParseFile(Parse, out retVal);

                return retVal;
            }

            public ExpressionResultsLine(string hugo_symbol_, SingleExpressionResult[] nonExclusiveResultsByRange_, SingleExpressionResult[] exclusiveResultsByRange_, int nTumorsExcluded_, int nZero_, int nOne_, int nMore_,
                int nSingleContribuingToAltFraction_, int nMultipleContributingToAltFraction_, double singleAltFraction_, double multipleAltFraction_, string rawLine_)
            {
                hugo_symbol = hugo_symbol_;
                nonExclusiveResultsByRange = nonExclusiveResultsByRange_;
                exclusiveResultsByRange = exclusiveResultsByRange_;
                nTumorsExcluded = nTumorsExcluded_;
                nZero = nZero_;
                nOne = nOne_;
                nMore = nMore_;
                nSingleContibutingToAltFraction = nSingleContribuingToAltFraction_;
                nMultipleContributingToAltFraction = nMultipleContributingToAltFraction_;
                singleAltFraction = singleAltFraction_;
                multipleAltFraction = multipleAltFraction_;
                rawLine = rawLine_;

                resultsByRange.Add(true, exclusiveResultsByRange);
                resultsByRange.Add(false, nonExclusiveResultsByRange);
            }

            static ExpressionResultsLine Parse(HeaderizedFile<ExpressionResultsLine>.FieldGrabber fieldGrabber)
            {
                var nonExclusiveResults = new SingleExpressionResult[nRegions];
                var exclusiveResults = new SingleExpressionResult[nRegions];

                for (int i = 0; i < nRegions; i++)
                {
                    nonExclusiveResults[i] = SingleExpressionResult.Parse(i, fieldGrabber, false, false);
                    exclusiveResults[i] = SingleExpressionResult.Parse(i, fieldGrabber, false, true);
                }

                return new ExpressionResultsLine(fieldGrabber.AsString("Hugo Symbol"), nonExclusiveResults, exclusiveResults, fieldGrabber.AsIntMinusOneIfStarOrEmptyString("nTumorsExcluded"),
                    fieldGrabber.AsIntMinusOneIfStarOrEmptyString("nZero"), fieldGrabber.AsIntMinusOneIfStarOrEmptyString("nOne"), fieldGrabber.AsIntMinusOneIfStarOrEmptyString("nMore"),
                    fieldGrabber.AsInt("nSingleContibutingToAltFraction"), fieldGrabber.AsInt("nMultipleContributingToAltFraction"),
                    fieldGrabber.AsDoubleNegativeInfinityIfStarOrEmptyString("single alt fraction"), fieldGrabber.AsDoubleNegativeInfinityIfStarOrEmptyString("multiple alt fraction"),
                    fieldGrabber.rawLine());
            }

            public static string getHeaderString()
            {
                var result = "Hugo Symbol";
                for (int inclusive = 0; inclusive < 2; inclusive++)
                {
                    for (int i = 0; i < nRegions; i++)
                    {
                        result += "\t" + SingleExpressionResult.getHeaderString(inclusive == 0, i);
                    }
                }

                foreach (bool inclusive in BothBools)
                {
                    for (int zeroOneManySilent = 0; zeroOneManySilent < 4; zeroOneManySilent++)
                    {
                        for (int i = 0; i < nRegions; i++)
                        {
                            result += "\t" + ((zeroOneManySilent == 3) ? "one silent" : ZeroOneManyString(zeroOneManySilent)) + " mutation " + regionIndexToString(i) + (inclusive ? "" : " exclusive");
                        }
                    }
                }

                return result + "\tnTumorsExcluded\tnZero\tnOne\tnMore\tnSingleContibutingToAltFraction\tnMultipleContributingToAltFraction\tsingle alt fraction\tmultiple alt fraction";
            }

            public readonly string hugo_symbol;
            public readonly SingleExpressionResult[] nonExclusiveResultsByRange = new SingleExpressionResult[nRegions];
            public readonly SingleExpressionResult[] exclusiveResultsByRange = new SingleExpressionResult[nRegions];
            public readonly Dictionary<bool, SingleExpressionResult[]> resultsByRange = new Dictionary<bool, SingleExpressionResult[]>();

            //
            // These are -1 if not filled in.
            //
            public readonly int nTumorsExcluded;
            public readonly int nZero;
            public readonly int nOne;
            public readonly int nMore;
            public readonly int nSingleContibutingToAltFraction;
            public readonly int nMultipleContributingToAltFraction;
            public readonly double singleAltFraction;
            public readonly double multipleAltFraction;
            public readonly string rawLine;
        } // ExpressionResultsLine

        //
        // Handling power-of-two * 1Kb regions.
        //
        public const int nRegions = 21;    // 0, 1 Kb, 2Kb, 4Kb .. 256Kb, Whole Autosome

        public static string regionIndexToString(int regionIndex)
        {
            switch (regionIndex)
            {
                case 0: return "0Kbp";
                case 1: return "1Kbp";
                case 2: return "2Kbp";
                case 3: return "4Kbp";
                case 4: return "8Kbp";
                case 5: return "16Kbp";
                case 6: return "32Kbp";
                case 7: return "64Kbp";
                case 8: return "128Kbp";
                case 9: return "256Kbp";
                case 10: return "512Kbp";
                case 11: return "1Mbp";
                case 12: return "2Mbp";
                case 13: return "4Mbp";
                case 14: return "8Mbp";
                case 15: return "16Mbp";
                case 16: return "32Mbp";
                case 17: return "64Mbp";
                case 18: return "128Mbp";
                case 19: return "256Mbp";
                case 20: return "whole autosome";
                default: return "Incorerect region index " + regionIndex;
            }
        }

        public static int regionIndexToSizeInBases(int regionIndex)
        {
            if (0 == regionIndex) return 0;
            if (nRegions - 1 == regionIndex) return int.MaxValue;
            return (1 << regionIndex) * 1000;
        }


        public class AVLTree<TValue> where TValue : IComparable
        {
            public static AVLTree<TValue> CreateFromList(List<TValue> list)
            {
                var tree = new AVLTree<TValue>();
                list.ForEach(x => tree.Insert(x));
                return tree;
            }

            class Node
            {
                public TValue value;

                public Node left = null, right = null, parent = null;

                public enum Balance { AVLLeft, AVLBalanced, AVLRight, AVLNew };
                public Balance balance = Balance.AVLNew;

                public int checkAndReturnDepth()
                {
                    //
                    // Links match up.
                    //
                    if (left != null && left.parent != this ||
                        right != null && right.parent != this)
                    {
                        throw new InvalidDataException();   // More-or-less arbitrary exception.
                    }

                    //
                    // The binary search tree ordering property applies.
                    //
                    if (left != null && value.CompareTo(left.value) <= 0 ||
                        right != null && value.CompareTo(right.value) >= 0)
                    {
                        throw new InvalidDataException();   // More-or-less arbitrary exception.
                    }

                    int leftDepth;
                    if (left != null)
                    {
                        leftDepth = left.checkAndReturnDepth();
                    } else
                    {
                        leftDepth = 0;
                    }

                    int rightDepth;
                    if (right != null) {
                        rightDepth = right.checkAndReturnDepth();
                    } else
                    {
                        rightDepth = 0;
                    }

                    if (leftDepth == rightDepth)
                    {
                        if (balance != Balance.AVLBalanced)
                        {
                            throw new InvalidDataException();
                        }
                        return leftDepth + 1;
                    }

                    if (leftDepth == rightDepth + 1)
                    {
                        if (balance != Balance.AVLLeft)
                        {
                            throw new InvalidDataException();
                        }
                        return leftDepth + 1;
                    }

                    if (leftDepth + 1 == rightDepth)
                    {
                        if (balance != Balance.AVLRight)
                        {
                            throw new InvalidDataException();
                        }
                        return rightDepth + 1;
                    }

                    //
                    // We're out of balance.
                    //
                    throw new InvalidDataException();
                }
            } //Node


            public void check()
            {
                if (root == null)
                {
                    return;
                }

                if (root.parent != null)
                {
                    throw new InvalidDataException();
                }

                root.checkAndReturnDepth();
            }

            //
            // This is the way to allow tree[element] and tree[element] = newElement;  Note that the second one is the only way to update a value in the tree without first deleting it and re-inserting it.
            // However, if it doesn't exist then it will insert it new.  So, tree[element] = element is the same as "update if existing, insert otherwise."
            //
            public TValue this[TValue key]
            {
                get
                {
                    TValue retVal;
                    if (!Lookup(key, out retVal))
                    {
                        throw new ArgumentException("Indexed value not found in AVL tree");
                    }

                    return retVal;
                }

                set
                {
                    Node currentCandiate = root;
                    while (currentCandiate != null)
                    {
                        int comparison = key.CompareTo(currentCandiate.value);
                        if (comparison == 0)
                        {
                            currentCandiate.value = value;
                            return;
                        }

                        if (comparison < 0)
                        {
                            currentCandiate = currentCandiate.left;
                        } else
                        {
                            currentCandiate = currentCandiate.right;
                        }
                    }

                    Insert(value);
                }
            }

            public bool findMin(out TValue returnValue)
            {
                if (root == null)
                {
                    returnValue = default(TValue);
                    return false;
                }

                Node currentNode = root;
                while (currentNode.left != null)
                {
                    currentNode = currentNode.left;
                }

                returnValue = currentNode.value;
                return true;
            }

            public bool findMax(out TValue returnValue)
            {
                if (root == null)
                {
                    returnValue = default(TValue);
                    return false;
                }

                Node currentNode = root;
                while (currentNode.right != null)
                {
                    currentNode = currentNode.right;
                }

                returnValue = currentNode.value;
                return true;
            }

            //
            // min and max return the default value if the tree is empty.
            //
            public TValue min()
            {
                TValue retVal;
                findMin(out retVal);
                return retVal;
            }

            public TValue max()
            {
                TValue retVal;
                findMax(out retVal);
                return retVal;
            }

            public bool FindFirstGreaterThanOrEqualTo(TValue key, out TValue returnValue)
            {
                Node currentNode = root;
                Node lastLeft = null;

                while (currentNode != null)
                {
                    int comparison = key.CompareTo(currentNode.value);

                    if (comparison == 0)
                    {
                        returnValue = currentNode.value;
                        return true;
                    }

                    if (comparison > 0)
                    {
                        currentNode = currentNode.right;
                    }
                    else
                    {
                        lastLeft = currentNode;
                        currentNode = currentNode.left;
                    }
                } // while (currentNode != null)

                if (lastLeft == null)
                {
                    returnValue = default(TValue);
                    return false;
                }

                returnValue = lastLeft.value;
                return true;
            } // FindFirstGreaterThanOrEqualTo

            public bool FindFirstLessThanOrEqualTo(TValue key, out TValue returnValue)
            {
                Node currentNode = root;
                Node lastRight = null;

                while (currentNode != null)
                {
                    int comparison = key.CompareTo(currentNode.value);

                    if (comparison == 0)
                    {
                        returnValue = currentNode.value;
                        return true;
                    }

                    if (comparison < 0)
                    {
                        currentNode = currentNode.left;
                    }
                    else
                    {
                        lastRight = currentNode;
                        currentNode = currentNode.right;
                    }
                } // while (currentNode != null)

                if (lastRight == null)
                {
                    returnValue = default(TValue);
                    return false;
                }

                returnValue = lastRight.value;
                return true;
            } // FindFirstLessThanOrEqualTo

            public TValue FindFirstLessThan(TValue key)
            {
                TValue retVal;
                if (!FindFirstLessThan(key, out retVal))
                {
                    return default(TValue);
                }
                return retVal;
            } // FindFirstLessThan

            public TValue FindFirstLessThanOrEqualTo(TValue key)
            {
                TValue retVal;
                if (!FindFirstLessThanOrEqualTo(key, out retVal))
                {
                    return default(TValue);
                }
                return retVal;
            } // FindFirstLessThanOrEqualTo

            public TValue FindFirstGreaterThan(TValue key)
            {
                TValue retVal;
                if (!FindFirstGreaterThan(key, out retVal))
                {
                    return default(TValue);
                }
                return retVal;
            } // FindFirstGreaterThan

            public TValue FindNextGreaterThanOrEqualTo(TValue key)
            {
                TValue retVal;
                if (!FindFirstGreaterThanOrEqualTo(key, out retVal))
                {
                    return default(TValue);
                }

                return retVal;
            } // FindNextGreaterThanOrEqualTo


            public bool FindFirstLessThan(TValue key, out TValue retVal)
            {
                Node currentNode = root;
                Node lastLessThan = null;

                while (currentNode != null)
                {
                    int comparison = key.CompareTo(currentNode.value);
                    if (comparison <= 0)
                    {
                        currentNode = currentNode.left;
                    } else
                    {
                        lastLessThan = currentNode;
                        currentNode = currentNode.right;
                    }
                }

                if (lastLessThan != null)
                {
                    retVal = lastLessThan.value;
                    return true;
                }

                retVal = default(TValue);
                return false;
            }

            public bool FindFirstGreaterThan(TValue key, out TValue retVal)
            {
                Node currentNode = root;
                Node lastGreaterThan = null;

                while (currentNode != null)
                {
                    int comparison = key.CompareTo(currentNode.value);
                    if (comparison >= 0)
                    {
                        currentNode = currentNode.right;
                    }
                    else
                    {
                        lastGreaterThan = currentNode;
                        currentNode = currentNode.left;
                    }
                }

                if (lastGreaterThan != null)
                {
                    retVal = lastGreaterThan.value;
                    return true;
                }

                retVal = default(TValue);
                return false;
            }
            public bool Lookup(TValue key, out TValue returnValue)
            {
                if (!FindFirstLessThanOrEqualTo(key, out returnValue) || returnValue.CompareTo(key) != 0)
                {
                    returnValue = default(TValue);
                    return false;
                }

                return true;
            } // Lookup

            public bool ContainsKey(TValue key)
            {
                TValue value;
                return Lookup(key, out value);
            }

            public void Insert(TValue newValue)
            {
                inserts++;

                var newNode = new Node();
                newNode.value = newValue;
                newNode.balance = Node.Balance.AVLBalanced;

                if (root == null)
                {
                    root = newNode;
                    return;
                }

                Node currentNode = root;
                Node previousNode = null;

                while (currentNode != null)
                {
                    previousNode = currentNode;
                    int comparison = currentNode.value.CompareTo(newValue);
                    if (comparison < 0)
                    {
                        currentNode = currentNode.right;
                    } else if (comparison > 0)
                    {
                        currentNode = currentNode.left;
                    } else
                    {
                        throw new ArgumentException("Trying to insert a duplicate element into an AVL tree.");
                    }
                } // while (currentNode != null)

                newNode.parent = previousNode;
                if (previousNode.value.CompareTo(newValue) < 0)
                {
                    previousNode.right = newNode;
                    rightAdded(previousNode);
                } else
                {
                    previousNode.left = newNode;
                    leftAdded(previousNode);
                }
            }

            void rightAdded(Node node)
            {
                if (node.balance == Node.Balance.AVLLeft)
                {
                    node.balance = Node.Balance.AVLBalanced;
                    //
                    // The depth of the subtree rooted here hasn't changed, we're done.
                    //
                    return;
                }

                if (node.balance == Node.Balance.AVLBalanced)
                {
                    //
                    // We've just gotten one deeper, but are still balanced.  Update and recurse up the tree.
                    //
                    node.balance = Node.Balance.AVLRight;
                    if (node == root)
                    {
                        return;
                    }

                    if (node.parent.right == node)
                    {
                        rightAdded(node.parent);
                    }
                    else
                    {
                        leftAdded(node.parent);
                    }
                    return;

                } // If we went to right heavy

                //
                // We went to double right.  Rotate.
                //
                if (node.right.balance == Node.Balance.AVLRight)
                {
                    singleRotate(node, node.right, Node.Balance.AVLRight);
                } else
                {
                    doubleRotate(node, node.right, node.right.left, Node.Balance.AVLRight);
                }
            } // rightAdded

            void leftAdded(Node node)
            {
                if (node.balance == Node.Balance.AVLRight)
                {
                    node.balance = Node.Balance.AVLBalanced;
                    //
                    // The depth of the subtree rooted here hasn't changed, we're done.
                    //
                    return;
                }

                if (node.balance == Node.Balance.AVLBalanced)
                {
                    //
                    // We've just gotten one deeper, but are still balanced.  Update and recurse up the tree.
                    //
                    node.balance = Node.Balance.AVLLeft;
                    if (node == root)
                    {
                        return;
                    }

                    if (node.parent.right == node)
                    {
                        rightAdded(node.parent);
                    }
                    else
                    {
                        leftAdded(node.parent);
                    }
                    return;
                } // if node was balanced

                //
                // We went to double left.  Rotate.
                //
                if (node.left.balance == Node.Balance.AVLLeft)
                {
                    singleRotate(node, node.left, Node.Balance.AVLLeft);
                }
                else
                {
                    doubleRotate(node, node.left, node.left.right, Node.Balance.AVLLeft);
                }
            }
            void singleRotate(Node node, Node child, Node.Balance whichSide)
            {
                //
                // Promote the child to our position in the tree.
                //

                if (node.parent != null)
                {
                    if (node.parent.left == node)
                    {
                        node.parent.left = child;
                    } else
                    {
                        node.parent.right = child;
                    }
                    child.parent = node.parent;
                }
                else
                {
                    //
                    // We're the root of the tree.
                    //
                    root = child;
                    child.parent = null;
                }

                //
                // Attach the child's light subtree to our heavy side (ie., where the child is attched now)
                // Then, attach us to the child's light subtree.
                //
                if (whichSide == Node.Balance.AVLRight)
                {
                    node.right = child.left;
                    if (node.right != null)
                    {
                        node.right.parent = node;
                    }
                    child.left = node;
                }
                else
                {
                    node.left = child.right;
                    if (node.left != null)
                    {
                        node.left.parent = node;
                    }
                    child.right = node;
                }

                node.parent = child;

                //
                // Finally, now both our and our (former) child are balanced.
                //
                node.balance = Node.Balance.AVLBalanced;
                child.balance = Node.Balance.AVLBalanced;
                //
                // NB: One of the cases in delete will result in the above balance settings being incorrect.
                // That case fixes up the settings after we return.
                //
            }

            void doubleRotate(Node node, Node child, Node grandchild, Node.Balance whichSide)
            {
                //
                // Write down a copy of all of the subtrees; see Knuth v3 p454 for the picture.
                // NOTE: The alpha and delta trees are never moved, so we don't store them.
                //

                Node beta, gamma;
                if (whichSide == Node.Balance.AVLRight)
                {
                    beta = grandchild.left;
                    gamma = grandchild.right;
                } else
                {
                    beta = grandchild.right;
                    gamma = grandchild.left;
                }

                //
                // Promote grandchild to our position.
                //
                if (root != node)
                {
                    if (node.parent.left == node)
                    {
                        node.parent.left = grandchild;
                    }
                    else
                    {
                        node.parent.right = grandchild;
                    }
                } else
                {
                    root = grandchild;
                }
                grandchild.parent = node.parent;

                //
                // Attach appropriate children to grandchild.
                //
                if (whichSide == Node.Balance.AVLRight)
                {
                    grandchild.right = child;
                    grandchild.left = node;
                }
                else
                {
                    grandchild.right = node;
                    grandchild.left = child;
                }
                node.parent = grandchild;
                child.parent = grandchild;

                //
                // Attach beta and gamma to node and child.
                //
                if (whichSide == Node.Balance.AVLRight)
                {
                    node.right = beta;
                    if (beta != null)
                    {
                        beta.parent = node;
                    }
                    child.left = gamma;
                    if (gamma != null)
                    {
                        gamma.parent = child;
                    }
                }
                else
                {
                    node.left = beta;
                    if (beta != null)
                    {
                        beta.parent = node;
                    }
                    child.right = gamma;
                    if (gamma != null)
                    {
                        gamma.parent = child;
                    }
                }

                //
                // Now upate the balance fields.
                //
                switch (grandchild.balance)
                {
                    case Node.Balance.AVLLeft:
                        if (whichSide == Node.Balance.AVLRight)
                        {
                            node.balance = Node.Balance.AVLBalanced;
                            child.balance = Node.Balance.AVLRight;
                        } else
                        {
                            node.balance = Node.Balance.AVLRight;
                            child.balance = Node.Balance.AVLBalanced;
                        }
                        break;

                    case Node.Balance.AVLBalanced:
                        node.balance = Node.Balance.AVLBalanced;
                        child.balance = Node.Balance.AVLBalanced;
                        break;

                    case Node.Balance.AVLRight:
                        if (whichSide == Node.Balance.AVLRight)
                        {
                            node.balance = Node.Balance.AVLLeft;
                            child.balance = Node.Balance.AVLBalanced;
                        } else
                        {
                            node.balance = Node.Balance.AVLBalanced;
                            child.balance = Node.Balance.AVLLeft;
                        }
                        break;
                } // switch

                grandchild.balance = Node.Balance.AVLBalanced;
            } // doubleRotate

            public void Delete(TValue key)
            {
                deletes++;

                Node nodeToRemove = root;
                while (nodeToRemove != null)
                {
                    int comparison = key.CompareTo(nodeToRemove.value);
                    if (comparison == 0)
                    {
                        break;
                    } else if (comparison > 0)
                    {
                        nodeToRemove = nodeToRemove.right;
                    } else
                    {
                        nodeToRemove = nodeToRemove.left;
                    }
                } // while nodeToRemove != null

                if (nodeToRemove == null)
                {
                    throw new ArgumentException("Tried to remove a value from an AVL tree that's not in the tree.");
                }

                if (nodeToRemove.left == null)
                {
                    //
                    // The right child either doesn't exist or is a leaf (because of the AVL balance property).
                    //
                    if (nodeToRemove.right != null)
                    {
                        nodeToRemove.right.parent = nodeToRemove.parent;
                    }

                    if (nodeToRemove.parent != null)
                    {
                        if (nodeToRemove.parent.left == nodeToRemove)
                        {
                            nodeToRemove.parent.left = nodeToRemove.right;
                            gotOneShorter(nodeToRemove.parent, Node.Balance.AVLLeft);
                        } else
                        {
                            nodeToRemove.parent.right = nodeToRemove.right;
                            gotOneShorter(nodeToRemove.parent, Node.Balance.AVLRight);
                        }
                    } else
                    {
                        root = nodeToRemove.right;
                    }
                } else if (nodeToRemove.right == null)
                {
                    //
                    // The left child must be a leaf because of the AVL balance property
                    //
                    nodeToRemove.left.parent = nodeToRemove.parent;
                    if (nodeToRemove.parent != null)
                    {
                        if (nodeToRemove.parent.left == nodeToRemove)
                        {
                            nodeToRemove.parent.left = nodeToRemove.left;
                            gotOneShorter(nodeToRemove.parent, Node.Balance.AVLLeft);
                        } else
                        {
                            nodeToRemove.parent.right = nodeToRemove.left;
                            gotOneShorter(nodeToRemove.parent, Node.Balance.AVLRight);
                        }
                    } else
                    {
                        root = nodeToRemove.left;
                    }
                } else // Node has both children
                {
                    //
                    // Find the symmetric successor and promote it.  The symmetric successor is the smallest element in the right
                    // subtree; it's found by following all left links in the right subtree until we find a node with no left link.
                    // That node may be promoted to the place of this without corrupting the binary tree ordering properties. (We could
                    // just as easily use the symmetric predecessor by finding the largest element in the left subtree, but there's
                    // no point.)
                    //

                    Node successorCandidate = nodeToRemove.right;
                    while (successorCandidate.left != null)
                    {
                        successorCandidate = successorCandidate.left;
                    }

                    Node shorterRoot;
                    Node.Balance shorterSide;

                    if (successorCandidate.parent.left == successorCandidate)
                    {
                        shorterRoot = successorCandidate.parent;
                        shorterSide = Node.Balance.AVLLeft;
                        successorCandidate.parent.left = successorCandidate.right;
                        if (successorCandidate.right != null)
                        {
                            successorCandidate.right.parent = successorCandidate.parent;
                        }
                        successorCandidate.right = nodeToRemove.right;
                        successorCandidate.left = nodeToRemove.left;
                        successorCandidate.balance = nodeToRemove.balance;
                        successorCandidate.right.parent = successorCandidate;
                        successorCandidate.left.parent = successorCandidate;
                        if (nodeToRemove.parent != null)
                        {
                            if (nodeToRemove.parent.left == nodeToRemove)
                            {
                                nodeToRemove.parent.left = successorCandidate;
                            } else
                            {
                                nodeToRemove.parent.right = successorCandidate;
                            }
                        } else
                        {
                            root = successorCandidate;
                        }
                        successorCandidate.parent = nodeToRemove.parent;
                    } else
                    {
                        //
                        // The successor was our child.  Just directly promote it.
                        //
                        if (nodeToRemove.parent != null)
                        {
                            if (nodeToRemove.parent.right == nodeToRemove)
                            {
                                nodeToRemove.parent.right = successorCandidate;
                            } else
                            {
                                nodeToRemove.parent.left = successorCandidate;
                            }
                        } else
                        {
                            root = successorCandidate;
                        }
                        successorCandidate.parent = nodeToRemove.parent;
                        successorCandidate.left = nodeToRemove.left;
                        if (nodeToRemove.left != null)
                        {
                            nodeToRemove.left.parent = successorCandidate;
                        }
                        //
                        // We just made our right subtree shorter.
                        //
                        successorCandidate.balance = nodeToRemove.balance;
                        shorterRoot = successorCandidate;
                        shorterSide = Node.Balance.AVLRight;
                    }
                    if (shorterRoot != null)
                    {
                        gotOneShorter(shorterRoot, shorterSide);
                    }
                } // Whether nodeToRemove has right only, left only or both children
            } // Delete

            void gotOneShorter(Node node, Node.Balance whichSide)
            {
                if (node.balance == Node.Balance.AVLBalanced)
                {
                    //
                    // We've just shrunk one subttree, but our depth has stayed the same.
                    // Reset our balance indicator and punt.
                    //
                    if (whichSide == Node.Balance.AVLRight)
                    {
                        node.balance = Node.Balance.AVLLeft;
                    } else
                    {
                        node.balance = Node.Balance.AVLRight;
                    }
                    return;
                } else if (node.balance == whichSide)
                {
                    //
                    // We just shrunk our heavy side; set our balance to neutral and recurse up the tree
                    //
                    node.balance = Node.Balance.AVLBalanced;
                    if (node.parent != null)
                    {
                        if (node.parent.right == node)
                        {
                            gotOneShorter(node.parent, Node.Balance.AVLRight);
                        } else
                        {
                            gotOneShorter(node.parent, Node.Balance.AVLLeft);
                        }
                    } // else we were the root and we're done.
                    return;
                } else
                {
                    //
                    // We've just gone out of balance.  Figure out a rotation to do.  This is almost like having added a
                    // node to the opposide side, except that the opposite side might be balanced.
                    //
                    Node.Balance heavySide;
                    Node heavyChild;
                    Node replacement;

                    if (whichSide == Node.Balance.AVLRight)
                    {
                        heavySide = Node.Balance.AVLLeft;
                        heavyChild = node.left;
                    } else
                    {
                        heavySide = Node.Balance.AVLRight;
                        heavyChild = node.right;
                    }

                    if (heavyChild.balance == heavySide)
                    {
                        //
                        // Typical single rotation case
                        //
                        singleRotate(node, heavyChild, heavySide);
                        replacement = heavyChild;
                    } else if (heavyChild.balance == whichSide)
                    {
                        //
                        // Typical double rotation case.
                        //
                        Node grandchild;
                        if (heavySide == Node.Balance.AVLRight)
                        {
                            grandchild = heavyChild.left;
                        } else
                        {
                            grandchild = heavyChild.right;
                        }

                        doubleRotate(node, heavyChild, grandchild, heavySide);
                        replacement = grandchild;
                    } else
                    {
                        singleRotate(node, heavyChild, heavySide);
                        //
                        // singleRotate has incorrectly set the balances; reset them.
                        //
                        node.balance = heavySide;
                        heavyChild.balance = whichSide;
                        // Overall depth hasn't changed; we're done.
                        return;
                    }

                    // NB: we have now changed position in the tree, so parent, right & left have changed!
                    if (replacement.parent == null)
                    {
                        return;
                    }

                    if (replacement.parent.right == replacement)
                    {
                        gotOneShorter(replacement.parent, Node.Balance.AVLRight);
                    } else
                    {
                        gotOneShorter(replacement.parent, Node.Balance.AVLLeft);
                    }
                }
            } // gotOneShorter

            public void ForEach(Action<TValue> action)
            {
                for (TValue x = min(); x != null; x = FindFirstGreaterThan(x))
                {
                    action(x);
                }
            }

            Node root = null;
            int inserts = 0;
            int deletes = 0;

        } // AVLTree

        public static void Shuffle<T>(IList<T> list)
        {
            Random rng = new Random();

            int n = list.Count;
            while (n > 1)
            {
                n--;
                int k = rng.Next(n + 1);
                T value = list[k];
                list[k] = list[n];
                list[n] = value;
            }
        }

        public class EnsembleGeneFile
        {
            static Tuple<string, string> ParseLine(Dictionary<string, int> fieldMappings, string[] fields)
            {
                var ensemblId = fields[fieldMappings["Gene stable ID"]];
                var geneName = fields[fieldMappings["Gene name"]];

                return new Tuple<string, string>(ensemblId, geneName);
            }

            public static Dictionary<string, string> ReadFile(string filename)
            {
                StreamReader inputFile = CreateStreamReaderWithRetry(filename);

                var wantedFields = new List<string>();
                wantedFields.Add("Gene stable ID");
                wantedFields.Add("Gene name");

                var headerizedFile = new HeaderizedFile<Tuple<string, string>>(inputFile, false, false, "", wantedFields);

                List<Tuple<string, string>> linesFromThisFile;

                headerizedFile.ParseFile(ParseLine, out linesFromThisFile);

                return linesFromThisFile.ToDictionary(r => r.Item1, r => r.Item2);
            }
        }

        public class CopyNumberVariation
        {
            //public string Sample; This changed to GDC_aliquot in the new data, but we're not using it, so I just deleted it.
            public string Chromosome;
            public int Start;
            public int End;
            public int Num_Probes;
            public double Segment_Mean;

            CopyNumberVariation(
                string Chromosome_,
                int Start_,
                int End_,
                int Num_Probes_,
                double Segment_Mean_)
            {
                Chromosome = Chromosome_;
                Start = Start_;
                End = End_;
                Num_Probes = Num_Probes_;
                Segment_Mean = Segment_Mean_;
            }
            static CopyNumberVariation ParseLine(Dictionary<string, int> fieldMappings, string[] fields)
            {
                // Start and End can be formatted in scientific notation, so we convert them accordingly
                return new CopyNumberVariation(
                    fields[fieldMappings["Chromosome"]],
                    Convert.ToInt32(Double.Parse(fields[fieldMappings["Start"]], System.Globalization.NumberStyles.Float)),
                    Convert.ToInt32(Double.Parse(fields[fieldMappings["End"]], System.Globalization.NumberStyles.Float)),
                    Convert.ToInt32(fields[fieldMappings["Num_Probes"]]),
                    Convert.ToDouble(fields[fieldMappings["Segment_Mean"]])
                    );
            }

            public bool OverlapsLocus(string Chromosome_, int Start_, int End_)
            {
                if (this.Chromosome != Chromosome_)
                {
                    return false;
                }

                if (this.Start < End_ && this.End > Start_)
                {
                    return true;
                }
                else
                {
                    return false;
                }
            }

            static public int getBasesinCNV(List<CopyNumberVariation> CNVs, double threshold = 0.2)
            {
                return CNVs.Select(r => r.End - r.Start).Sum();
            }

            static public List<CopyNumberVariation> ReadFile(string filename)
            {
                if (filename == "")
                {
                    return null;
                }

                StreamReader inputFile;

                inputFile = CreateStreamReaderWithRetry(filename);

                var neededFields = new List<string>();
                neededFields.Add("Chromosome");
                neededFields.Add("Start");
                neededFields.Add("End");
                neededFields.Add("Num_Probes");
                neededFields.Add("Segment_Mean");

                var headerizedFile = new HeaderizedFile<CopyNumberVariation>(inputFile, false, false, "#version gdc-1.0.0", neededFields);

                List<CopyNumberVariation> result;

                if (!headerizedFile.ParseFile((a, b) => ParseLine(a, b), out result))
                {
                    Console.WriteLine("Error reading Annotation File " + filename);
                    return null;
                }

                inputFile.Close();

                return result;
            } // ReadFile

            static public Dictionary<bool, List<CopyNumberVariation>> ReadBothFiles(Case case_)
            {
                var retVal = new Dictionary<bool, List<CopyNumberVariation>>();
                if (case_.tumor_copy_number_filename == "")
                {
                    retVal.Add(true, null);
                }
                else
                {
                    retVal.Add(true, ReadFile(case_.tumor_copy_number_filename).Where(r => Math.Abs(r.Segment_Mean) > 1.0).ToList());
                }
                if (case_.normal_copy_number_filename == "") {
                    retVal.Add(false, null);
                } else {
                    retVal.Add(false, ReadFile(case_.normal_copy_number_filename).Where(r => Math.Abs(r.Segment_Mean) > 1.0).ToList());
                }
                return retVal;
            }
        } // CopyNumberVariation

        public class IGV
        {
            // writes batch file that visualizes CNV, normal and tumor DNA, as well as tumor RNA
            // specify case, output file and optional position, in the form of "Gene_name"
            // or "chrX:start-end"
            public static void writeBatchFile(ASETools.Case case_, string filename, string position = "")
            {
                var outFile = ASETools.CreateStreamWriterWithRetry(filename);

                // start new session
                outFile.WriteLine("new");

                // Write normal CNV data
                if (case_.normal_copy_number_filename != "")
                {
                    outFile.WriteLine("load " + case_.normal_copy_number_filename);
                }

                // Write normal CNV data
                if (case_.tumor_copy_number_filename != "")
                {
                    outFile.WriteLine("load " + case_.tumor_copy_number_filename);
                }

                // load normal DNA, tumor DNA and tumor RNA
                outFile.WriteLine("load " + case_.normal_dna_filename);
                outFile.WriteLine("load " + case_.tumor_dna_filename);
                outFile.WriteLine("load " + case_.tumor_rna_filename);

                if (position != "")
                {
                    outFile.WriteLine("goto " + position);
                }

                outFile.Close();

            }
        }


        public class FPKMFile
        {
            public static Dictionary<string, double> ReadFile(string filename, Configuration configuration)
            {
                // first, load in default ensembl id => gene names
                var idMap = EnsembleGeneFile.ReadFile(configuration.ensemblToGeneFilename);

                StreamReader inputFile;

                if (filename.Count() > 2 && filename.Substring(filename.Count() - 3, 3) == ".gz")
                {
                    inputFile = CreateCompressedStreamReaderWithRetry(filename);
                }
                else
                {
                    inputFile = CreateStreamReaderWithRetry(filename);
                }

                var result = new Dictionary<string, double>();

                string line;
                while ((line = inputFile.ReadLine()) != null)
                {
                    var split = line.Split('\t');
                    var ens_id = split[0].Split('.').First();
                    var fpkm = Convert.ToDouble(split[1]);

                    string geneName;
                    if (idMap.TryGetValue(ens_id, out geneName))
                    {
                        // there are some inconsistent genes appearing more than once
                        if (!result.ContainsKey(geneName))
                        {
                            result.Add(geneName, fpkm);
                        }
                    }
                }

                inputFile.Close();
                return result;
            }
        } // FPKMFile

        public class AllSignificantResultsLine
        {
            public AllSignificantResultsLine(string hugo_symbol_, double ase_one_mutation_, double ase_not_one_mutation_, string inputFile_, bool exclusive_, bool one_vs_many_, int range_index_, string range_, double p_)
            {
                hugo_symbol = hugo_symbol_;
                ase_one_mutation = ase_one_mutation_;
                ase_not_one_mutation = ase_not_one_mutation_;
                inputFile = inputFile_;
                exclusive = exclusive_;
                one_vs_many = one_vs_many_;
                range_index = range_index_;
                range = range_;
                p = p_;
            }
            public static List<AllSignificantResultsLine> ReadFile(string filename)
            {
                var inputFile = ASETools.CreateStreamReaderWithRetry(filename);

                if (null == inputFile)
                {
                    return null;
                }

                var wantedFields = new List<string>();
                wantedFields.Add("Hugo Symbol");
                wantedFields.Add("ASE (one mutation)");
                wantedFields.Add("ASE (not one mutation)");
                wantedFields.Add("input file");
                wantedFields.Add("exclusive");
                wantedFields.Add("OneVsMany");
                wantedFields.Add("range index");
                wantedFields.Add("range");
                wantedFields.Add("p");

                var headerizedFile = new HeaderizedFile<AllSignificantResultsLine>(inputFile, false, true, "", wantedFields);
                List<AllSignificantResultsLine> retVal;

                headerizedFile.ParseFile(parser, out retVal);

                inputFile.Close();

                return retVal;
            }

            static AllSignificantResultsLine parser(HeaderizedFile<AllSignificantResultsLine>.FieldGrabber fieldGrabber)
            {
                return new AllSignificantResultsLine(fieldGrabber.AsString("Hugo Symbol"), fieldGrabber.AsDouble("ASE (one mutation)"), fieldGrabber.AsDouble("ASE (not one mutation)"),
                    fieldGrabber.AsString("input file"), fieldGrabber.AsBool("exclusive"), fieldGrabber.AsBool("OneVsMany"), fieldGrabber.AsInt("range index"), fieldGrabber.AsString("range"), fieldGrabber.AsDouble("p"));
            }

            public readonly string hugo_symbol;
            public readonly double ase_one_mutation;
            public readonly double ase_not_one_mutation;
            public readonly string inputFile;
            public readonly bool exclusive;
            public readonly bool one_vs_many;
            public readonly int range_index;
            public readonly string range;
            public readonly double p;

        } // AllSignificantResultsLine

        public class ASEMapLine
        {
            public static List<ASEMapLine> ReadFile(string pathnmame)
            {
                string[] wantedFields = { "Chromosome", "locus", "tumor", "n cases", "index", "mean", "standard deviation" };

                var inputFile = CreateStreamReaderWithRetry(pathnmame);
                if (inputFile == null)
                {
                    Console.WriteLine("Unable to open " + pathnmame);
                    return null;
                }

                var headerizedFile = new HeaderizedFile<ASEMapLine>(inputFile, false, true, "", wantedFields.ToList());

                List<ASEMapLine> retVal;
                headerizedFile.ParseFile(Parse, out retVal);

                inputFile.Close();

                return retVal;
            }

            static ASEMapLine Parse(HeaderizedFile<ASEMapLine>.FieldGrabber fieldGrabber)
            {
                return new ASEMapLine(fieldGrabber.AsString("Chromosome"), fieldGrabber.AsInt("locus"), fieldGrabber.AsBool("tumor"), fieldGrabber.AsInt("n cases"),
                    fieldGrabber.AsInt("index"), fieldGrabber.AsDouble("mean"), fieldGrabber.AsDouble("standard deviation"));
            }

            ASEMapLine(string chromosome_, int locus_, bool tumor_, int nCases_, int index_, double mean_, double stdDev_)
            {
                chromosome = chromosome_;
                locus = locus_;
                tumor = tumor_;
                nCases = nCases_;
                index = index_;
                mean = mean_;
                stdDev = stdDev_;
            }

            public readonly string chromosome;
            public readonly int locus;
            public readonly bool tumor;
            public readonly int nCases;
            public readonly int index;              // This is the x coordinate for the scatter graph.
            public readonly double mean;
            public readonly double stdDev;
        } // ASEMapLine

        public class ASEMapPerGeneLine
        {
            public static Dictionary<string, ASEMapPerGeneLine> ReadFromFileToDictionary(string filename)
            {
                var allLines = ReadFromFile(filename);

                if (null == allLines)
                {
                    return null;
                }

                var result = new Dictionary<string, ASEMapPerGeneLine>();

                allLines.ForEach(x => result.Add(x.hugo_symbol, x));

                return result;
            }

            public static List<ASEMapPerGeneLine> ReadFromFile(string filename)
            {
                var inputFile = CreateStreamReaderWithRetry(filename);
                if (inputFile == null)
                {
                    Console.WriteLine("ASEMapPerGeneLine.ReadFromFile: unable to open file " + filename);
                    return null;
                }

                string[] wantedFields = { "Hugo Symbol", "n Tumor Samples", "mean Tumor ASE",  "Tumor Fraction of RNA at Variant Sites", "standard deviation of Tumor ASE",
                    "n Normal Samples", "mean Normal ASE", "Normal Fraction of RNA at Variant Sites", "standard deviation of Normal ASE", "Tumor ASE minus Normal ASE",
                    "Absolute value of Tumor ASE minus Normal ASE", "Bonferroni Corrected MW p value for normal differing from tumor",
                    "chromosome", "min locus", "max locus" };

                var headerizedFile = new HeaderizedFile<ASEMapPerGeneLine>(inputFile, false, true, "", wantedFields.ToList());

                List<ASEMapPerGeneLine> result;
                if (!headerizedFile.ParseFile(parseLine, out result))
                {
                    return null;
                }

                inputFile.Close();

                return result;
            }

            static ASEMapPerGeneLine parseLine(ASETools.HeaderizedFile<ASEMapPerGeneLine>.FieldGrabber fieldGrabber)
            {
                return new ASEMapPerGeneLine(fieldGrabber.AsString("Hugo Symbol"), fieldGrabber.AsInt("n Tumor Samples"), fieldGrabber.AsDouble("mean Tumor ASE"), fieldGrabber.AsDouble("standard deviation of Tumor ASE"), fieldGrabber.AsDouble("Tumor Fraction of RNA at Variant Sites"),
                    fieldGrabber.AsInt("n Normal Samples"), fieldGrabber.AsDouble("mean Normal ASE"), fieldGrabber.AsDouble("standard deviation of Normal ASE"), fieldGrabber.AsDouble("Normal Fraction of RNA at Variant Sites"),
                    fieldGrabber.AsDouble("Tumor ASE minus Normal ASE"),
                    fieldGrabber.AsString("chromosome"), fieldGrabber.AsInt("min locus"),
                    fieldGrabber.AsInt("max locus"));
            }

            ASEMapPerGeneLine(string hugo_symbol_, int nSamplesTumor, double meanASETumor, double stdDeviationTumor, double meanFractionRNAReadsTumor,
                int nSamplesNormal, double meanASENormal, double stdDeviationNormal, double meanFractionRNAReadsNormal, double TumorExcessASEOverNormal_, string chromosome_, int minLocus_, int maxLocus_)
            {
                hugo_symbol = ConvertToNonExcelString(hugo_symbol_);
                sampleData.Add(true, new SampleData(nSamplesTumor, meanASETumor, stdDeviationTumor, meanFractionRNAReadsTumor));
                sampleData.Add(false, new SampleData(nSamplesNormal, meanASENormal, stdDeviationNormal, meanFractionRNAReadsNormal));


                TumorExcessASEOverNormal = TumorExcessASEOverNormal_;
                chromosome = chromosome_;
                minLocus = minLocus_;
                maxLocus = maxLocus_;
            }

            public class SampleData
            {
                public SampleData(int nSamples_, double meanASE_, double stdDeviation_, double meanFractionRNAReads_)
                {
                    nSamples = nSamples_;
                    meanASE = meanASE_;
                    stdDeviation = stdDeviation_;
                    meanFractionRNAReads = meanFractionRNAReads_;
                }
                public readonly int nSamples;
                public readonly double meanASE;
                public readonly double stdDeviation;
                public readonly double meanFractionRNAReads;
            }

            public readonly string hugo_symbol;
            public readonly Dictionary<bool, SampleData> sampleData = new Dictionary<bool, SampleData>();

            public readonly double TumorExcessASEOverNormal;   // Signed difference between tumor meanASE and normal meanASE
            public readonly string chromosome;
            public readonly int minLocus;
            public readonly int maxLocus;
        } // ASEMapPerGeneLine


        public interface ClinicalSummaryLine
        {
            string getVitalStatus();
            int getOverallSurvivalInDays();
            string getPatientId();
            string getGender();
        }
        //
        // This for the BeatAML data, not TCGA
        //
        public class BeatAMLClinicalSummaryLine : ClinicalSummaryLine
        {
            static BeatAMLClinicalSummaryLine Parse(HeaderizedFile<BeatAMLClinicalSummaryLine>.FieldGrabber fieldGrabber)
            {
                return new BeatAMLClinicalSummaryLine(
                    fieldGrabber.AsString("patientId"),
                    fieldGrabber.AsString("labId"),
                    fieldGrabber.AsString("gender"),
                    fieldGrabber.AsIntMinusOneIfStarOrEmptyString("ageAtDiagnosis"),
                    fieldGrabber.YOrNToBool("priorMalignancyNonMyeloid"),
                    fieldGrabber.AsString("dxAtInclusion"),
                    fieldGrabber.AsString("specificDxAtInclusion"),
                    fieldGrabber.AsString("riskGroup"),
                    fieldGrabber.AsString("karyotype"),
                    fieldGrabber.AsString("dxAtSpecimenAcquisition"),
                    fieldGrabber.AsIntMinusOneIfStarOrEmptyString("ageAtSpecimenAcquisition"),
                    fieldGrabber.AsIntMinusOneIfStarOrEmptyString("timeOfSampleCollectionRelativeToInclusion"),
                    fieldGrabber.YOrNToBool("rnaSeq"),
                    fieldGrabber.YOrNToBool("exomeSeq"),
                    fieldGrabber.AsString("vitalStatus"),
                    fieldGrabber.AsIntMinusOneIfStarOrEmptyString("overallSurvival"),
                    fieldGrabber.AsString("causeOfDeath"),
                    fieldGrabber.AsString("TP53"),
                    fieldGrabber.AsString("TP53_GeneTrails"),
                    fieldGrabber.AsString("priorTreatmentRegimens"),
                    fieldGrabber.AsString("currentRegimen")
                    );
            }

            public int getOverallSurvivalInDays()
            {
                return overallSurvival;
            }

            public string getVitalStatus()
            {
                return vitalStatus;
            }

            public string getPatientId()
            {
                return patientId;
            }

            public string getGender()
            {
                return gender;
            }

            public static List<BeatAMLClinicalSummaryLine> readFile(string filename)
            {
                var inputFile = CreateStreamReaderWithRetry(filename);

                if (null == inputFile)
                {
                    return null;
                }

                string[] wantedFields =
                {
                    "patientId",
                    "labId",
                    "gender",
                    "ageAtDiagnosis",
                    "priorMalignancyNonMyeloid",
                    "priorMalignancyType",
                    "priorMalignancyChemo",
                    "priorMalignancyRadiationTx",
                    "priorMDS",
                    "priorMDSMoreThanTwoMths",
                    "priorMDSMPN",
                    "priorMDSMPNMoreThanTwoMths",
                    "priorMPN",
                    "priorMPNMoreThanTwoMths",
                    "performanceStatus",
                    "dxAtInclusion",
                    "specificDxAtInclusion",
                    "fabBlastMorphology",
                    "riskGroup",
                    "karyotype",
                    "otherCytogenetics",
                    "dxAtSpecimenAcquisition",
                    "ageAtSpecimenAcquisition",
                    "timeOfSampleCollectionRelativeToInclusion",
                    "specimenGroups",
                    "specimenType",
                    "rnaSeq",
                    "exomeSeq",
                    "ic50",
                    "priorTreatmentTypeCount",
                    "priorTreatmentTypes",
                    "priorTreatmentRegimenCount",
                    "priorTreatmentRegimens",
                    "priorTreatmentStageCount",
                    "priorTreatmentStages",
                    "responseToInductionTx",
                    "typeInductionTx",
                    "responseDurationToInductionTx",
                    "currentTreatmentType",
                    "currentRegimen",
                    "currentStage",
                    "currentTreatmentDuration",
                    "vitalStatus",
                    "overallSurvival",
                    "causeOfDeath",
                    "percentBlastsBM",
                    "percentBlastsPB",
                    "percentAbnormalPlasmaBM",
                    "percentBandsPB",
                    "percentBasophilsPB",
                    "percentEosinophilsPB",
                    "percentImmatureGranulocytesPB",
                    "percentLymphocytesPB",
                    "percentMetamyelocytesPB",
                    "percentMonocytesPB",
                    "percentMyelocytesPB",
                    "percentNeutrophilsPB",
                    "percentNucleatedRBCsPB",
                    "percentPromonocytes",
                    "percentPromyelocytes",
                    "percentPromyelocytesPB",
                    "percentReactiveLymphocytesPB",
                    "percentWBC",
                    "wbcCount",
                    "plateletCount",
                    "albumin",
                    "bCellGeneRearrangement",
                    "bilirubin",
                    "creatinine",
                    "hemoglobin",
                    "surfaceAntigensImmunohistochemicalStains",
                    "tCellReceptorGene",
                    "ALT",
                    "AST",
                    "ASXL1",
                    "ASXL1_GeneTrails",
                    "ATM",
                    "BCL_1_t11_14",
                    "BCL_2_14_18",
                    "BCOR",
                    "BCOR_GeneTrails",
                    "BCORL1",
                    "BCR_ABL",
                    "BCR_ABL_International_Scale",
                    "BCR_ABL_FISH_",
                    "BCR_ABL_LOG_DROP",
                    "BCR_ABL_RNA_Quant",
                    "BRAF",
                    "BRAF_Sequenome",
                    "CAD",
                    "Calreticulin",
                    "CBL",
                    "CBL_GeneTrails",
                    "CCND2",
                    "CD36",
                    "CEBPA",
                    "CEBPA_GeneTrails",
                    "CEBPA_Sequenome",
                    "CHEK2",
                    "CIITA",
                    "CREBBP",
                    "CREBBP_GeneTrails",
                    "CSF3R",
                    "CSF3R_GeneTrails",
                    "CUX1",
                    "DNMT3A",
                    "DNMT3A_GeneTrails",
                    "ETV6",
                    "ETV6_GeneTrails",
                    "EZH2",
                    "EZH2_GeneTrails",
                    "FBXW7_GeneTrails",
                    "FLT3",
                    "FLT3_GeneTrails",
                    "FLT3_Sequenome",
                    "FLT3_D835",
                    "FLT3_D835_GeneTrails",
                    "FLT3_D835_Sequenome",
                    "FLT3_ITD",
                    "FLT3_ITD_GeneTrails",
                    "FLT3_ITD_Internal",
                    "FLT3_N676K_GeneTrails",
                    "FLT3_V491L",
                    "GATA1",
                    "GATA1_GeneTrails",
                    "GATA2",
                    "GATA2_Genetrails",
                    "GNAS",
                    "IDH1",
                    "IDH1_GeneTrails",
                    "IDH1_Sequenome",
                    "IDH2",
                    "IDH2_GeneTrails",
                    "IDH2_Sequenome",
                    "IKZF1",
                    "IKZF1_GeneTrails",
                    "JAK1_GeneTrails",
                    "JAK2",
                    "JAK2_GeneTrails",
                    "JAK2_Sequenome",
                    "JAK3",
                    "JAK3_GeneTrails",
                    "KDM6A",
                    "KDM6A_GeneTrails",
                    "KIT",
                    "KIT_GeneTrails",
                    "KIT_Sequenome",
                    "KMT2A",
                    "KRAS",
                    "KRAS_GeneTrails",
                    "KRAS_Sequenome",
                    "LDH",
                    "MCV",
                    "MEN1",
                    "MLL",
                    "MLL_GeneTrails",
                    "MLL2",
                    "MLL2_GeneTrails",
                    "MPL",
                    "MPL_GeneTrails",
                    "MUTYH",
                    "FGFR1_ZNF198",
                    "HNF1A",
                    "NF1",
                    "NOTCH1",
                    "NOTCH1_GeneTrails",
                    "NPM",
                    "NPM_GeneTrails",
                    "NPM_Sequenome",
                    "NPM_Internal",
                    "NRAS",
                    "NRAS_GeneTrails",
                    "NRAS_Sequenome",
                    "NRBC",
                    "NTRK2",
                    "PAX5",
                    "PAX5_GeneTrails",
                    "PHF6",
                    "POT1",
                    "PRDM1_GeneTrails",
                    "PTPN11",
                    "PTPN11_GeneTrails",
                    "PTPN11_Sequenome",
                    "ROS1",
                    "RUNX1",
                    "RUNX1_GeneTrails",
                    "SETBP1",
                    "CSF1R",
                    "SF1",
                    "SF3B1",
                    "SF3B1_GeneTrails",
                    "SRSF2",
                    "SRSF2_GeneTrails",
                    "STAG2",
                    "STAG2_GeneTrails",
                    "STAT3",
                    "STAT3_Genetrails",
                    "STK11",
                    "SUZ12",
                    "SUZ12_GeneTrails",
                    "TET2",
                    "TET2_GeneTrails",
                    "TP53",
                    "TP53_GeneTrails",
                    "U2AF1",
                    "U2AF1_GeneTrails",
                    "WT1",
                    "WT1_GeneTrails",
                    "ZRSR2",
                    "ZRSR2_GeneTrails",
                    "PML_RAR_Alpha_t_15_17",
                    "PML_RARA_Quantitivate"
                };

                var headerizedFile = new HeaderizedFile<BeatAMLClinicalSummaryLine>(inputFile, false, false, "", wantedFields.ToList());

                List<BeatAMLClinicalSummaryLine> result;
                headerizedFile.ParseFile(Parse, out result);

                inputFile.Close();

                return result;
            }



            BeatAMLClinicalSummaryLine(
                string patientId_,
                string labId_,
                string gender_,
                int ageAtDiag_,
                bool priorMalignencyNonMyeloid_,
                string dxAtInclusion_,
                string specificDxAtInclusion_,
                string riskGroup_,
                string karyotype_,
                string dXAtSpeciminAcquisition_,
                int ageAtSpeciminAcquisition_,
                int timeOfSampleCollectionRelativeToInclusion_,
                bool rnaSeq_,
                bool exomeSeq_,
                string vitalStatus_,
                int overallSurvival_,
                string causeOfDeath_,
                string TP53_,
                string TP53_Genetrails_,
                string priorTreatmentRegimens_,
                string currentRegimen_
            )
            {
                patientId = patientId_;
                labId = labId_;
                gender = gender_;
                ageAtDiag = ageAtDiag_;
                priorMalignencyNonMyeloid = priorMalignencyNonMyeloid_;
                dxAtInclusion = dxAtInclusion_;
                specificDxAtInclusion = specificDxAtInclusion_;
                riskGroup = riskGroup_;
                karyotype = karyotype_;
                dXAtSpeciminAcquisition = dXAtSpeciminAcquisition_;
                ageAtSpeciminAcquisition = ageAtSpeciminAcquisition_;
                timeOfSampleCollectionRelativeToInclusion = timeOfSampleCollectionRelativeToInclusion_;
                rnaSeq = rnaSeq_;
                exomeSeq = exomeSeq_;
                vitalStatus = vitalStatus_;
                overallSurvival = overallSurvival_;
                causeOfDeath = causeOfDeath_;
                TP53 = TP53_;
                TP53_Genetrails = TP53_Genetrails_;
                priorTreatmentRegimens = priorTreatmentRegimens_;
                currentRegimen = currentRegimen_;
            }

            public bool p53Mutant()
            {
                return TP53.ToLower().Contains("positive") || TP53_Genetrails.ToLower().Contains("positive");
            }

            public readonly string patientId;
            public readonly string labId;
            public readonly string gender;
            public readonly int ageAtDiag;
            public readonly bool priorMalignencyNonMyeloid;
            public readonly string dxAtInclusion;
            public readonly string specificDxAtInclusion;
            public readonly string riskGroup;
            public readonly string karyotype;
            public readonly string dXAtSpeciminAcquisition;
            public readonly int ageAtSpeciminAcquisition;
            public readonly int timeOfSampleCollectionRelativeToInclusion;
            public readonly bool rnaSeq;
            public readonly bool exomeSeq;
            public /*readonly*/ string vitalStatus; // Turn off readonly for a hack in AnalyzeBeatAMLCases
            public readonly int overallSurvival;
            public readonly string causeOfDeath;
            public readonly string TP53;
            public readonly string TP53_Genetrails;
            public readonly string priorTreatmentRegimens;  // This is actually a list separated by "|", but for now we just get the raw string
            public readonly string currentRegimen;
        } // ClinicalSummaryLine

        public class TCGAClinicalSummaryLine : ClinicalSummaryLine
        {
            public static List<TCGAClinicalSummaryLine> readFromFile(string filename)
            {
                var inputFile = CreateStreamReaderWithRetry(filename);

                if (null == inputFile)
                {
                    throw new Exception("TCGAClinicalSummaryLine.readFromFile(" + filename + "): unable to open file.");
                }

                //
                // These things don't have consistent headers, of course.  So, read the first line, see what it has, then choose the appropriate option.
                //

                string daysToDeathFieldName;
                string ageAtDiagnosisFieldName;
                var header = inputFile.ReadLine();
                if (header == null)
                {
                    return null;
                }

                if (header.Contains("days_to_death"))
                {
                    daysToDeathFieldName = "days_to_death";
                } else if (header.Contains("death_days_to"))
                {
                    daysToDeathFieldName = "death_days_to";
                } else
                {
                    Console.Write("No recognized days to death field name in header.  These are the fields that contain the word 'death':");
                    header.Split('\t').Where(x => x.Contains("death")).ToList().ForEach(x => Console.WriteLine(" " + x));
                    return null;
                }

                if (header.Contains("age_at_diagnosis"))
                {
                    ageAtDiagnosisFieldName = "age_at_diagnosis";
                } else if (header.Contains("age_at_initial_pathologic_diagnosis"))
                {
                    ageAtDiagnosisFieldName = "age_at_initial_pathologic_diagnosis";
                } else
                {
                    Console.Write("No recognized age at diagnosis field name in header.  These are the fields that contain the word 'age':");
                    header.Split('\t').Where(x => x.Contains("age")).ToList().ForEach(x => Console.WriteLine(" " + x));
                    return null;
                }

                //
                // Now reopen the input file.
                //
                inputFile.Close();
                inputFile = CreateStreamReaderWithRetry(filename);

                if (null == inputFile)
                {
                    return null;
                }

                string[] wantedFields =
                {
                    "bcr_patient_uuid",
                    "vital_status",
                    "gender",
                };  // Tons more stuff we could include.

                var wantedFieldsList = wantedFields.ToList();
                wantedFieldsList.Add(daysToDeathFieldName);
                wantedFieldsList.Add(ageAtDiagnosisFieldName);

                // 
                // These files are weird: they have two copies of the header line and then a second, slightly different header, and then a third really weird line with CDE_ID on it.
                // So we skip the two weird ones.
                //
                string[] additionalValuesToTreatAsStar = { "[Not Applicable]", "[Discrepancy]", "[Not Available]", "[Completed]" };
                var headerizedFile = new HeaderizedFile<TCGAClinicalSummaryLine>(inputFile, false, false, "", wantedFieldsList, rowsToSkipAfterHeader_: 2, additionalValuesToTreatAsStar_: additionalValuesToTreatAsStar.ToList());
                List<TCGAClinicalSummaryLine> result;
                headerizedFile.ParseFile(x => parse(x, daysToDeathFieldName, ageAtDiagnosisFieldName), out result);

                return result;
            }

            static TCGAClinicalSummaryLine parse(HeaderizedFile<TCGAClinicalSummaryLine>.FieldGrabber fieldGrabber, string daysToDeathFieldName, string ageAtDiagnosisFieldName)
            {
                return new TCGAClinicalSummaryLine(fieldGrabber.AsString("bcr_patient_uuid").ToLower(), fieldGrabber.AsIntMinusOneIfStarOrEmptyString(daysToDeathFieldName), fieldGrabber.AsString("vital_status"), fieldGrabber.AsString("gender").ToLower(),
                    fieldGrabber.AsIntMinusOneIfStarOrEmptyString(ageAtDiagnosisFieldName));
            }

            TCGAClinicalSummaryLine(string patientId_, int days_to_death_, string vitalStatus_, string gender_, int age_at_diagnosis_)
            {
                patientId = patientId_;
                days_to_death = days_to_death_;
                vitalStatus = vitalStatus_;
                gender = gender_;
                age_at_diagnosis = age_at_diagnosis_;
            }

            public string getPatientId()
            {
                return patientId;
            }

            public int getOverallSurvivalInDays()
            {
                return days_to_death;
            }

            public string getVitalStatus()
            {
                return vitalStatus;
            }

            public string getGender()
            {
                return gender;
            }

            public static Dictionary<string, TCGAClinicalSummaryLine> readAllToDictionary(Configuration configuration)
            {
                var retVal = new Dictionary<string, TCGAClinicalSummaryLine>();

                foreach (var file in Directory.EnumerateFiles(configuration.patient_metadata_directory))
                {
                    if (!file.EndsWith(".txt"))
                    {
                        continue;   // ., .., etc.
                    }

                    readFromFile(file).ForEach(_ => retVal.Add(_.getPatientId(), _));
                }

                return retVal;
            }

            public static List<string> getInputFilenames(Configuration configuration)
            {
                return Directory.EnumerateFiles(configuration.patient_metadata_directory).Where(_ => _.EndsWith(".txt")).ToList();
            }

            public readonly string patientId;
            public readonly int days_to_death;
            public readonly string vitalStatus;
            public readonly string gender;
            public readonly int age_at_diagnosis;
        } // TCGAClinicalSummaryLine

        public class KaplanMeierPoint
        {
            public KaplanMeierPoint(int daysSinceInclusion_, double fractionAlive_)
            {
                daysSinceInclusion = daysSinceInclusion_;
                fractionAlive = fractionAlive_;

            }
            public int daysSinceInclusion;
            public double fractionAlive;
        }

        static public List<KaplanMeierPoint> KaplanMeier(List<ClinicalSummaryLine> summaryLines, out int n, out double slope, out double rSquared)
        {
            var oneLinePerPatient = summaryLines.GroupByToDict<ClinicalSummaryLine, string>(x => x.getPatientId()).Select(x => x.Value[0]).Where(x => x.getOverallSurvivalInDays() > 0).ToList();
            n = oneLinePerPatient.Count();
            oneLinePerPatient.Sort((x, y) => x.getOverallSurvivalInDays().CompareTo(y.getOverallSurvivalInDays()));   // This isn't necessary, but it's nice for debugging 

            var longestSurvivor = summaryLines.Select(x => x.getOverallSurvivalInDays()).Max();

            var fractionAliveByDay = new double[longestSurvivor + 1];
            fractionAliveByDay[0] = 1;

            for (int day = 1; day <= longestSurvivor; day++)
            {
                fractionAliveByDay[day] = (double)oneLinePerPatient.Where(x => x.getOverallSurvivalInDays() >= day).Count() / oneLinePerPatient.Where(x => x.getOverallSurvivalInDays() >= day || x.getVitalStatus() == "Dead").Count();
            }

            var retVal = new List<KaplanMeierPoint>();
            for (int day = 0; day <= longestSurvivor; day++)
            {
                if (day == 0 || fractionAliveByDay[day] != fractionAliveByDay[day - 1]) {
                    retVal.Add(new KaplanMeierPoint(day, fractionAliveByDay[day]));
                }
            }

            var nonZeroLogPoints = new List<Tuple<double, double>>();   // Maps days since inclusion -> log2(fractionAlive), excluding fractionAlive = 0 point if any (since you can't take a log of 0).
            for (int i = 0; i < retVal.Count(); i++)
            {
                if (retVal[i].fractionAlive > 0)
                {
                    nonZeroLogPoints.Add(new Tuple<double, double>(retVal[i].daysSinceInclusion, Math.Log(retVal[i].fractionAlive, 2)));
                }
            }

            double slopeValue = MathNet.Numerics.LinearRegression.SimpleRegression.FitThroughOrigin(nonZeroLogPoints); // slopeValue is here because you can't pass an out parameter to a method as its value (in this case CoefficientToDetermination)
            slope = slopeValue;
            rSquared = MathNet.Numerics.GoodnessOfFit.CoefficientOfDetermination(nonZeroLogPoints.Select(x => x.Item1 * slopeValue), nonZeroLogPoints.Select(x => x.Item2));

            return retVal;
        } // KaplanMeier

        public class BeatAMLSampleMapping
        {
            public static Dictionary<string, List<BeatAMLSampleMapping>> readFromFile(string filename)
            {
                var inputFile = CreateStreamReaderWithRetry(filename);
                if (null == inputFile)
                {
                    return null;
                }

                string[] wantedFields =
                {
                    "PatientID",
                    "AML_Original_LabID",
                    "Normal_Original_LabID",
                    "AML_SeqID",
                    "Normal_SeqID",
                    "AML_BAM",
                    "Normal_BAM",
                    "AML_Outliers",
                    "Normal_Outliers",
                    "RNA_SeqID",
                    "RNA_BAM",
                    "RNA_Outliers"
                };

                var headerizedFile = new HeaderizedFile<BeatAMLSampleMapping>(inputFile, false, false, "", wantedFields.ToList());

                List<BeatAMLSampleMapping> allMappings;
                headerizedFile.ParseFile(parse, out allMappings);

                return allMappings.GroupByToDict<BeatAMLSampleMapping, string>(x => x.patientId);
            }

            static BeatAMLSampleMapping parse(HeaderizedFile<BeatAMLSampleMapping>.FieldGrabber fieldGrabber)
            {
                return new BeatAMLSampleMapping(fieldGrabber.AsString("PatientID"), fieldGrabber.AsStringEmptyForNA("AML_Original_LabID"), fieldGrabber.AsStringEmptyForNA("Normal_Original_LabID"),
                    fieldGrabber.AsStringEmptyForNA("AML_SeqID"), fieldGrabber.AsStringEmptyForNA("Normal_SeqID"), fieldGrabber.AsStringEmptyForNA("AML_BAM"), fieldGrabber.AsStringEmptyForNA("Normal_BAM"),
                    fieldGrabber.AsStringEmptyForNA("AML_Outliers"), fieldGrabber.AsStringEmptyForNA("Normal_Outliers"), fieldGrabber.AsStringEmptyForNA("RNA_SeqID"), fieldGrabber.AsStringEmptyForNA("RNA_BAM"),
                    fieldGrabber.AsStringEmptyForNA("RNA_Outliers"));
            }

            BeatAMLSampleMapping(string patientId_, string AML_Original_LabID_, string Normal_Original_LabID_, string AML_SeqID_, string Normal_SeqID_, string AML_BAM_, string Normal_BAM_, string AML_Outliers_, string Normal_Outliers_,
                string RNA_SeqID_, string RNA_BAM_, string RNA_Outliers_)
            {
                patientId = patientId_;
                AML_Original_LabID = AML_Original_LabID_;
                Normal_Original_LabID = Normal_Original_LabID_;
                AML_SeqID = AML_SeqID_;
                Normal_SeqID = Normal_SeqID_;
                AML_BAM = AML_BAM_;
                Normal_BAM = Normal_BAM_;
                AML_Outliers = AML_Outliers_;
                Normal_Outliers = Normal_Outliers_;
                RNA_SeqID = RNA_SeqID_;
                RNA_BAM = RNA_BAM_;
                RNA_Outliers = RNA_Outliers_;
            }

            public readonly string patientId;
            public readonly string AML_Original_LabID;
            public readonly string Normal_Original_LabID;
            public readonly string AML_SeqID;
            public readonly string Normal_SeqID;
            public readonly string AML_BAM;
            public readonly string Normal_BAM;
            public readonly string AML_Outliers;
            public readonly string Normal_Outliers;
            public readonly string RNA_SeqID;
            public readonly string RNA_BAM;
            public readonly string RNA_Outliers;
        } // BeatAMLSampleMapping

        public class PerCaseASE
        {
            static public Dictionary<string, PerCaseASE> loadAll(Configuration configuration)
            {
                var retVal = new Dictionary<string, PerCaseASE>();

                string[] wantedFields =
                {
                    "Case ID",
                    "Normal ASE",
                    "Tumor ASE",
                    "Normal median chromosome ASE",
                    "Normal min chromosome ASE",
                    "Normal max chromosome ASE",
                    "Tumor median chromosome ASE",
                    "Tumor min chromosome ASE",
                    "Tumor max chromosome ASE"
                };

                var wantedFieldsList = wantedFields.ToList();

                for (int i = 1; i <= nHumanAutosomes; i++)
                {
                    wantedFieldsList.Add("Normal chr" + i + " ASE");
                    wantedFieldsList.Add("Tumor chr" + i + " ASE");
                }


                var inputFile = CreateStreamReaderWithRetry(configuration.finalResultsDirectory + PerCaseASEFilename);
                var headerizedFile = new HeaderizedFile<PerCaseASE>(inputFile, false, true, "", wantedFieldsList);

                List<PerCaseASE> allPerCaseASE;
                headerizedFile.ParseFile(parse, out allPerCaseASE);

                allPerCaseASE.ForEach(x => retVal.Add(x.caseId, x));

                return retVal;
            } // loadAll

            static PerCaseASE parse(HeaderizedFile<PerCaseASE>.FieldGrabber fieldGrabber)
            {
                var perChromosomeNormalASE = new double[nHumanAutosomes + 1];   // +1 because indexed by chromosome number.  [0] is unused
                var perChromosomeTumorASE = new double[nHumanAutosomes + 1];   // +1 because indexed by chromosome number.  [0] is unused

                perChromosomeNormalASE[0] = perChromosomeTumorASE[0] = double.NegativeInfinity;

                for (int i = 1; i <= nHumanAutosomes; i++)
                {
                    perChromosomeNormalASE[i] = fieldGrabber.AsDoubleNegativeInfinityIfStarOrEmptyString("Normal chr" + i + " ASE");
                    perChromosomeTumorASE[i] = fieldGrabber.AsDoubleNegativeInfinityIfStarOrEmptyString("Tumor chr" + i + " ASE");
                }
                return new PerCaseASE(
                    fieldGrabber.AsString("Case ID"), fieldGrabber.AsDoubleNegativeInfinityIfStarOrEmptyString("Normal ASE"),
                    fieldGrabber.AsDoubleNegativeInfinityIfStarOrEmptyString("Tumor ASE"), perChromosomeNormalASE, perChromosomeTumorASE,
                    fieldGrabber.AsDoubleNegativeInfinityIfStarOrEmptyString("Tumor min chromosome ASE"),
                    fieldGrabber.AsDoubleNegativeInfinityIfStarOrEmptyString("Tumor median chromosome ASE"),
                    fieldGrabber.AsDoubleNegativeInfinityIfStarOrEmptyString("Tumor max chromosome ASE"),
                    fieldGrabber.AsDoubleNegativeInfinityIfStarOrEmptyString("Normal min chromosome ASE"),
                    fieldGrabber.AsDoubleNegativeInfinityIfStarOrEmptyString("Normal median chromosome ASE"),
                    fieldGrabber.AsDoubleNegativeInfinityIfStarOrEmptyString("Normal max chromosome ASE"));
            }

            PerCaseASE(string caseId_, double normalASE_, double tumorASE_, double[] perChromosomeNormalASE_, double[] perChromosomeTumorASE_,
                double tumorMinChromosome_, double tumorMedianChromosome_, double tumorMaxChromosome_, double normalMinChromosome_, double normalMedianChromosome_,
                double normalMaxChromosome_)
            {
                caseId = caseId_;
                normalASE = normalASE_;
                tumorASE = tumorASE_;
                perChromosomeNormalASE = perChromosomeNormalASE_;
                perChromosomeTumorASE = perChromosomeTumorASE_;
                tumorMinChromosome = tumorMinChromosome_;
                tumorMedianChromosome = tumorMedianChromosome_;
                tumorMaxChromosome = tumorMaxChromosome_;
                normalMinChromosome = normalMinChromosome_;
                normalMedianChromosome = normalMedianChromosome_;
                normalMaxChromosome = normalMaxChromosome_;
            }

            public readonly string caseId;
            public readonly double normalASE;
            public readonly double tumorASE;

            public readonly double[] perChromosomeNormalASE;
            public readonly double[] perChromosomeTumorASE;

            public readonly double tumorMinChromosome, tumorMedianChromosome, tumorMaxChromosome;
            public readonly double normalMinChromosome, normalMedianChromosome, normalMaxChromosome;
        }  // PerCaseASE

        public class WorkerThreadHelper<TQueueItem, TPerThreadState> where TPerThreadState : new()
        {
            public delegate void HandleOneItem(TQueueItem item, TPerThreadState perThreadState);
            public delegate void FinishUp(TPerThreadState perThreadState);
            public delegate void ItemDequeued();
            public WorkerThreadHelper(List<TQueueItem> queue_, HandleOneItem handleOneItem_, FinishUp finishUp_, ItemDequeued itemDequeued_, int nItemsPerDot_ = 0)
            {
                queue = new List<TQueueItem>();
                queue_.ForEach(_ => queue.Add(_));  // Copy the queue, so that we don't modify the one the caller proovided.
                handleOneItem = handleOneItem_;
                finishUp = finishUp_;
                itemDequeued = itemDequeued_;
                nItemsPerDot = nItemsPerDot_;
            }

            public void run()
            {
                run(Environment.ProcessorCount);
            }

            public void run(int nThreads)
            {
                if (nThreads <= 0)
                {
                    throw new Exception("ASETools.WorkerThreadHelper.run: thread count must be strictly positive.");
                }
                start(nThreads);
                join();
            }

            public void start(int nThreads)
            {
                nItemsProcessed = 0;
                threads = new List<Thread>();
                for (int i = 0; i < nThreads; i++)
                {
                    threads.Add(new Thread(() => WorkerThread(this)));
                }

                threads.ForEach(t => t.Start());
            }

            public void join()
            {
                threads.ForEach(t => t.Join());

                if (nItemsPerDot != 0)
                {
                    Console.WriteLine();
                }
            }

            static void WorkerThread(WorkerThreadHelper<TQueueItem, TPerThreadState> helper)
            {
                helper.WorkerThread();
            }

            void WorkerThread()
            {
                var perThreadState = new TPerThreadState();

                while (true)
                {
                    TQueueItem queueItem;
                    lock (queue)
                    {
                        if (queue.Count() == 0)
                        {
                            if (finishUp != null)
                            {
                                finishUp(perThreadState);
                            }
                            return;
                        }

                        queueItem = queue[0];
                        queue.RemoveAt(0);

                        if (itemDequeued != null)
                        {
                            itemDequeued();
                        }
                    }

                    if (null != handleOneItem)  // It's not clear what use it would be not to have one, but just in case
                    {
                        handleOneItem(queueItem, perThreadState);
                    }

                    lock (queue)
                    {
                        nItemsProcessed++;

                        if (nItemsPerDot != 0 && nItemsProcessed % nItemsPerDot == 0)
                        {
                            Console.Write(".");
                        }
                    }
                } // while (true)
            }

            List<TQueueItem> queue;
            HandleOneItem handleOneItem;
            FinishUp finishUp;
            ItemDequeued itemDequeued;
            int nItemsPerDot = 0;
            int nItemsProcessed = 0;
            List<Thread> threads;
        } // WorkerThreadHelper

        public class ASVThreadingHelper<TPerThreadState> where TPerThreadState : new()
        {

            public delegate void HandleOneCase(Case case_, TPerThreadState state, List<AnnotatedVariant> variantsForThisCase, ASECorrection aseCorrection, Dictionary<bool, List<CopyNumberVariation>> copyNumber);
            public delegate bool ASVDecider(AnnotatedVariant asv, Dictionary<bool, List<CopyNumberVariation>> copyNumber);
            public ASVThreadingHelper(List<Case> casesToProcess, ASECorrection aseCorrection_, ASVDecider asvDecider_, HandleOneCase handleOneCase_, WorkerThreadHelper<Case, TPerThreadState>.FinishUp finishUp_, WorkerThreadHelper<Case, TPerThreadState>.ItemDequeued itemDequeued_,
                int nItemsPerDot = 0)
            {
                aseCorrection = aseCorrection_;
                asvDecider = asvDecider_;
                handleOneCase = handleOneCase_;
                finishUp = finishUp_;
                itemDequeued = itemDequeued_;

                workerThreadHelper = new WorkerThreadHelper<Case, TPerThreadState>(casesToProcess,
                            (x, y) => this.HandleOneItem(x, y),
                            x => this.FinishUp(x),
                            () => this.ItemDequeued(),
                            nItemsPerDot);
            }
            public static ASVThreadingHelper<TPerThreadState> create(List<Case> casesToProcess, bool useSomatic, bool useGermline, bool tumor, Configuration configuration, Dictionary<string, ASEMapPerGeneLine> perGeneASEMap,
                GeneMap geneMap, ASECorrection aseCorrection,
                HandleOneCase handleOneCase, WorkerThreadHelper<Case, TPerThreadState>.FinishUp finishUp, WorkerThreadHelper<Case, TPerThreadState>.ItemDequeued itemDequeued,
                int inputMinRNACoverage = -1, int inputMinDNACoverage = -1,
                int nItemsPerDot = 0)
            {
                return new ASVThreadingHelper<TPerThreadState>(casesToProcess, aseCorrection,
                    (v, c) => ((v.somaticMutation && useSomatic && !v.isSilent()) || (!v.somaticMutation && useGermline)) &&
                        v.IsASECandidate(tumor, c, configuration, perGeneASEMap, geneMap, inputMinRNACoverage, inputMinDNACoverage),
                    handleOneCase, finishUp, itemDequeued, nItemsPerDot);
            }

            static bool loadState(Configuration configuration, out Dictionary<string, ASEMapPerGeneLine> perGeneASEMap, out GeneMap geneMap, out ASECorrection aseCorrection, out GeneLocationsByNameAndChromosome geneLocationInformation)
            {
                aseCorrection = null;
                geneLocationInformation = null;
                geneMap = null;

                perGeneASEMap = ASEMapPerGeneLine.ReadFromFileToDictionary(configuration.finalResultsDirectory + ASETools.PerGeneASEMapFilename);

                if (null == perGeneASEMap)
                {
                    Console.WriteLine("You must first create the per-gene ASE map in " + configuration.finalResultsDirectory + ASETools.PerGeneASEMapFilename);
                    return false;
                }

                aseCorrection = ASECorrection.LoadFromFile(configuration.finalResultsDirectory + ASETools.ASECorrectionFilename);
                if (null == aseCorrection)
                {
                    Console.WriteLine("Unable to load ASE correction");
                    return false;
                }

                geneLocationInformation = new GeneLocationsByNameAndChromosome(ASETools.readKnownGeneFile(configuration.geneLocationInformationFilename));
                geneMap = new GeneMap(geneLocationInformation.genesByName);

                return true;
            }

            public static ASVThreadingHelper<TPerThreadState> create(List<Case> casesToProcess, bool useSomatic, bool useGermline, bool tumor, Configuration configuration,
                 HandleOneCase handleOneCase, WorkerThreadHelper<Case, TPerThreadState>.FinishUp finishUp, WorkerThreadHelper<Case, TPerThreadState>.ItemDequeued itemDequeued,
                int inputMinRNACoverage = -1, int inputMinDNACoverage = -1,
                int nItemsPerDot = 0)
            {
                Dictionary<string, ASEMapPerGeneLine> perGeneASEMap;
                GeneMap geneMap;
                ASECorrection aseCorrection;
                GeneLocationsByNameAndChromosome geneLocationInformation;

                if (!loadState(configuration, out perGeneASEMap, out geneMap, out aseCorrection, out geneLocationInformation))
                {
                    return null;
                }

                return ASVThreadingHelper<TPerThreadState>.create(casesToProcess, useSomatic, useGermline, tumor, configuration, perGeneASEMap, geneMap, aseCorrection, handleOneCase, finishUp, itemDequeued, inputMinRNACoverage, inputMinDNACoverage, nItemsPerDot);
            }

            public static ASVThreadingHelper<TPerThreadState> create(List<Case> casesToProcess, Configuration configuration, ASVDecider asvDecider,
                HandleOneCase handleOneCase, WorkerThreadHelper<Case, TPerThreadState>.FinishUp finishUp, WorkerThreadHelper<Case, TPerThreadState>.ItemDequeued itemDequeued, int nItemsPerDot = 0)
            {
                Dictionary<string, ASEMapPerGeneLine> perGeneASEMap;
                GeneMap geneMap;
                ASECorrection aseCorrection;
                GeneLocationsByNameAndChromosome geneLocationInformation;

                if (!loadState(configuration, out perGeneASEMap, out geneMap, out aseCorrection, out geneLocationInformation))
                {
                    return null;
                }

                return new ASVThreadingHelper<TPerThreadState>(casesToProcess, aseCorrection, asvDecider, handleOneCase, finishUp, itemDequeued, nItemsPerDot);
            }

            public void run()
            {
                workerThreadHelper.run(Environment.ProcessorCount);
            }

            public void run(int nThreads)
            {
                workerThreadHelper.run(nThreads);
            }

            void HandleOneItem(Case case_, TPerThreadState perThreadState)
            {
                if (case_.annotated_selected_variants_filename == "" || case_.tumor_copy_number_filename == "")
                {
                    Console.WriteLine("ASV Threading Helper: case " + case_.case_id + " is missing annotated selected variants and/or copy numbers.  Skipping it.");
                    return;
                }

                var asv = AnnotatedVariant.readFile(case_.annotated_selected_variants_filename);
                if (null == asv)
                {
                    Console.WriteLine("ASV Threading Helper: case unable to read ASV from " + case_.annotated_selected_variants_filename + " Skipping case.");
                    return;
                }

                var copyNumber = CopyNumberVariation.ReadBothFiles(case_);
                if (null == copyNumber)
                {
                    Console.WriteLine("ASV Threading Helper: case unable to read ASV from " + case_.annotated_selected_variants_filename + " Skipping case.");
                    return;
                }

                var variantsToProcess = asv.Where(x => asvDecider(x, copyNumber)).ToList();

                handleOneCase(case_, perThreadState, variantsToProcess, aseCorrection, copyNumber);
            }

            void FinishUp(TPerThreadState perThreadState)
            {
                if (null != finishUp)
                {
                    finishUp(perThreadState);
                }
            }

            void ItemDequeued()
            {
                if (null != itemDequeued)
                {
                    itemDequeued();
                }
            }

            WorkerThreadHelper<Case, TPerThreadState> workerThreadHelper;
            HandleOneCase handleOneCase;
            WorkerThreadHelper<Case, TPerThreadState>.FinishUp finishUp;
            WorkerThreadHelper<Case, TPerThreadState>.ItemDequeued itemDequeued;
            ASVDecider asvDecider;


            ASECorrection aseCorrection;
        }

        public class ASECorrection
        {
            public static ASECorrection LoadFromFile(string filename)
            {
                var inputFile = CreateStreamReaderWithRetry(filename);

                if (null == inputFile)
                {
                    Console.WriteLine("Unable to open ASECorrection file " + filename);
                    return null;
                }

                var lines = new List<string>();

                inputFile.ReadLine();   // Skip the header
                {
                    string line;
                    while (null != (line = inputFile.ReadLine())) {
                        lines.Add(line);
                    }
                }
                inputFile.Close();

                if (lines.Count() == 0)
                {
                    Console.WriteLine("ASECorrection file " + filename + " doesn't appear to have any data.");
                    return null;
                }

                int maxReadDepth = lines.Select(x => Convert.ToInt32(x.Split('\t')[0])).ToList().Max();

                var retVal = new ASECorrection();
                retVal.maxReadDepth = maxReadDepth;
                retVal.correction = new double[retVal.maxReadDepth + 1, 101];

                for (int i = 0; i <= maxReadDepth; i++)
                {
                    for (int j = 0; j < 101; j++)
                    {
                        retVal.correction[i, j] = 0;
                    }
                }

                foreach (var line in lines)
                {
                    var fields = line.Split('\t');
                    if (fields.Count() != 102)
                    {
                        Console.WriteLine("ASECorrection file " + filename + " has a line with the wrong number of fields " + fields.Count() + " != 102: " + line);
                        return null;
                    }

                    int readDepth = Convert.ToInt32(fields[0]);
                    for (int i = 0; i < 101; i++)
                    {
                        if (fields[i + 1] == "*")
                        {
                            retVal.correction[readDepth, i] = 0;
                        }
                        else
                        {
                            retVal.correction[readDepth, i] = Convert.ToDouble(fields[i + 1]);
                        }
                    } // for each field

                } // for each line (read depth)
                return retVal;
            } // LoadFromFile

            public double getCorrectedASE(double rawASE, int readDepth)
            {
                if (readDepth > maxReadDepth)
                {
                    return rawASE;
                }

                return rawASE + correction[readDepth, (int)(Math.Round(rawASE * 100))];
            }

            int maxReadDepth = 0;

            double[,] correction = null;
        } // ASECorrection

        public class ASERepetitiveRegionMap
        {
            class SingleRepetitiveRegion : IComparable
            {
                public int chromosomeNumber;
                public int start;
                public int end;

                public SingleRepetitiveRegion(int chromosomeNumber_, int start_, int end_)
                {
                    chromosomeNumber = chromosomeNumber_;
                    start = start_;
                    end = end_;
                }

                public int CompareTo(object peerObject)
                {
                    SingleRepetitiveRegion peer = (SingleRepetitiveRegion)peerObject;

                    if (chromosomeNumber != peer.chromosomeNumber)
                    {
                        return chromosomeNumber.CompareTo(peer.chromosomeNumber);
                    }

                    return start.CompareTo(peer.start);
                }
            }
            static public ASERepetitiveRegionMap loadFromFile(string filename)
            {
                var inputFile = CreateStreamReaderWithRetry(filename);
                if (null == inputFile)
                {
                    Console.WriteLine("ASERepetitiveRegionMap: unable to open input file " + filename);
                    return null;
                }

                string[] wantedFields =
                {
                    "Chromosome",
                    "Begin",
                    "End"
                };

                var headerizedFile = new HeaderizedFile<SingleRepetitiveRegion>(inputFile, false, true, "", wantedFields.ToList());

                List<SingleRepetitiveRegion> lines;
                headerizedFile.ParseFile(parser, out lines);

                var tree = new AVLTree<SingleRepetitiveRegion>();
                lines.ForEach(x => tree.Insert(x));

                return new ASERepetitiveRegionMap(tree);
            }

            AVLTree<SingleRepetitiveRegion> tree;

            ASERepetitiveRegionMap(AVLTree<SingleRepetitiveRegion> tree_)
            {
                tree = tree_;
            }

            public bool isInRepetitiveRegion(string chromosome, int locus)
            {
                return isInRepetitiveRegion(chromosomeNameToNonChrForm(chromosome)[0], locus);
            }

            public bool isInRepetitiveRegion(char chromosomeNumber, int locus)
            {
                var key = new SingleRepetitiveRegion(chromosomeNumber, locus, locus);
                SingleRepetitiveRegion result;
                if (!tree.FindFirstLessThanOrEqualTo(key, out result))
                {
                    return false;
                }

                return result.chromosomeNumber == chromosomeNumber && result.end >= locus;
            }

            public bool isCloseToRepetitiveRegion(int chromosomeNumber, int locus, int range)
            {
                var key = new SingleRepetitiveRegion(chromosomeNumber, locus, locus);

                SingleRepetitiveRegion result;
                if (tree.FindFirstLessThanOrEqualTo(key, out result) && chromosomeNumber == result.chromosomeNumber && result.end + range >= locus)
                {
                    return true;
                }

                if (tree.FindFirstGreaterThanOrEqualTo(key, out result) && chromosomeNumber == result.chromosomeNumber && locus + range >= result.start)
                {
                    return true;
                }

                return false;
            }

            static SingleRepetitiveRegion parser(HeaderizedFile<SingleRepetitiveRegion>.FieldGrabber fieldGrabber)
            {
                return new SingleRepetitiveRegion(fieldGrabber.AsInt("Chromosome"), fieldGrabber.AsInt("Begin"), fieldGrabber.AsInt("End"));
            }
        } // ASERepetitiveRegionMap

        //
        // PDFs of ASE at each read depth
        //
        public class MeasuredASEMatrix
        {
            public MeasuredASEMatrix(int nDiscreteASEValues_, int maxReadDepth_)
            {
                nDiscreteASEValues = nDiscreteASEValues_;
                maxReadDepth = maxReadDepth_;
                pdfMatrix = new double[nDiscreteASEValues, maxReadDepth + 1];
            }

            double[,] pdfMatrix;
            int nDiscreteASEValues;
            int maxReadDepth;

            public double getValue(double ase, int readDepth)
            {
                int ASEIndex = (int)Math.Round(ase * (nDiscreteASEValues - 1));

                int column = Math.Min(maxReadDepth, readDepth);

                return pdfMatrix[ASEIndex, column];
            }

            static public MeasuredASEMatrix readFromFile(string filename, bool tumor, bool somatic)
            {
                var inputFile = CreateStreamReaderWithRetry(filename);

                if (null == inputFile)
                {
                    Console.WriteLine("MeasuredASEMatrix: unable to open " + filename);
                    return null;
                }

                var rows = new Dictionary<int, double[]>();

                int maxReadDepth = -1;

                string inputLine;
                while (null != (inputLine = inputFile.ReadLine()))
                {
                    if (inputLine.StartsWith("ASE distributon for tumor: "))
                    {
                        var words = inputLine.Split(' ');
                        if (words.Count() != 13)
                        {
                            Console.WriteLine("MeasuredASEMatrix: failed to parse header line " + inputLine);
                            continue;
                        }

                        if (tumor && words[4] != "True" || !tumor && words[4] != "False")
                        {
                            continue;
                        }

                        if (somatic && words[6] != "True" || !somatic && words[6] != "False")
                        {
                            continue;
                        }

                        int readDepth = Convert.ToInt32(words[12]);

                        rows.Add(readDepth, Histogram.ReadStreamToPDFValues(inputFile, false));

                        maxReadDepth = Math.Max(maxReadDepth, readDepth);
                    }
                }

                if (rows.Count() == 0)
                {
                    Console.WriteLine("MeasuredASEMatrix.readFromFile: Unable to find any data in " + filename);
                    return null;
                }

                int nDiscreteASEValues = rows[maxReadDepth].Count();
                if (rows.Any(x => x.Value.Count() != nDiscreteASEValues))
                {
                    Console.WriteLine("MeasuredASEMatrix.readFromFile: not all rows (pdf's) have the same number of values (discrete ASE values)");
                    return null;
                }

                var retVal = new MeasuredASEMatrix(nDiscreteASEValues, maxReadDepth);
                foreach (var rowEntry in rows)
                {
                    for (int i = 0; i < nDiscreteASEValues; i++)
                    {
                        retVal.pdfMatrix[i, rowEntry.Key] = rowEntry.Value[i];
                    }
                }

                return retVal;

            }
        } // MeasuredASEMatrix

        public class BonferroniCorrectedASEDistributionLine // This class doesn't contain all of the actual measurements, just the metadata.  Feel free to add the measurements if you need them.
        {
            public static List<BonferroniCorrectedASEDistributionLine> readFromFile(string filename)
            {
                string[] wantedFields =
                {
                    "Hugo Symbol",
                    "Alt fraction for single mutations",
                    "Alt fraction for multiple mutations",
                    "Min p",
                    "Significant@.01",
                    "Best 0 vs. 1 ratio for significant results",
                    "Best 0 vs. 1 ratio at",
                    "Significant At",
                    "Gene Size",
                    "nTumorsExcluded",
                    "nZero",
                    "nOne",
                    "nMore",
                    "nSingleContibutingToAltFraction",
                    "nMultipleContributingToAltFraction",
                };

                var inputFile = CreateStreamReaderWithRetry(filename);
                if (null == inputFile)
                {
                    Console.WriteLine("Unable to open " + filename);
                    return null;
                }

                var headerizedFile = new HeaderizedFile<BonferroniCorrectedASEDistributionLine>(inputFile, false, true, "", wantedFields.ToList());
                List<BonferroniCorrectedASEDistributionLine> retval;

                if (!headerizedFile.ParseFile(parser, out retval))
                {
                    inputFile.Close();
                    Console.WriteLine("BonferroniCorrectedASEDistributionLine: failed to parse file " + filename);
                    return null;
                }

                inputFile.Close();

                return retval;
            } // readFromFile

            static BonferroniCorrectedASEDistributionLine parser(HeaderizedFile<BonferroniCorrectedASEDistributionLine>.FieldGrabber fieldGrabber)
            {
                return new BonferroniCorrectedASEDistributionLine(
                    fieldGrabber.AsString("Hugo Symbol"), fieldGrabber.AsDoubleNegativeInfinityIfStarOrEmptyString("Alt fraction for single mutations"), fieldGrabber.AsDoubleNegativeInfinityIfStarOrEmptyString("Alt fraction for multiple mutations"),
                    fieldGrabber.AsDoubleNegativeInfinityIfStarOrEmptyString("Min p"), fieldGrabber.AsBool("Significant@.01"), fieldGrabber.AsDoubleNegativeInfinityIfStarOrEmptyString("Best 0 vs. 1 ratio for significant results"),
                    fieldGrabber.AsString("Best 0 vs. 1 ratio at"), fieldGrabber.AsString("Significant At"), fieldGrabber.AsIntMinusOneIfStarOrEmptyString("Gene Size"), fieldGrabber.AsInt("nTumorsExcluded"), fieldGrabber.AsInt("nZero"),
                    fieldGrabber.AsInt("nOne"), fieldGrabber.AsInt("nMore"), fieldGrabber.AsInt("nSingleContibutingToAltFraction"), fieldGrabber.AsInt("nMultipleContributingToAltFraction"));
            }

            BonferroniCorrectedASEDistributionLine(string hugo_symbol_, double altFractionForSingleMutations_, double altFractionForMultipleMutations_, double minP_, bool significant01_, double best0vs1Ratio_, string best0vs1RatioAt_,
                string significantAt_, int geneSize_, int nTumorsExcluded_, int nZero_, int nOne_, int nMore_, int nSingleContributingToAltFraction_, int nMultipleContributingToAltFraction_)
            {
                hugo_symbol = hugo_symbol_;
                altFractionForSingleMutations = altFractionForSingleMutations_;
                altFractionForMultipleMutations = altFractionForMultipleMutations_;
                minP = minP_;
                significant01 = significant01_;
                best0vs1Ratio = best0vs1Ratio_;
                best0vs1RatioAt = best0vs1RatioAt_;
                significantAt = significantAt_;
                significantAtArray = significantAt.Split(',');
                geneSize = geneSize_;
                nTumorsExcluded = nTumorsExcluded_;
                nZero = nZero_;
                nOne = nOne_;
                nMore = nMore_;
                nSingleContributingToAltFraction = nSingleContributingToAltFraction_;
                nMultipleContributingToAltFraction = nMultipleContributingToAltFraction_;
            }

            public readonly string hugo_symbol;
            public readonly double altFractionForSingleMutations;
            public readonly double altFractionForMultipleMutations;
            public readonly double minP;
            public readonly bool significant01;
            public readonly double best0vs1Ratio;
            public readonly string best0vs1RatioAt;
            public readonly string significantAt;
            public readonly string[] significantAtArray;
            public readonly int geneSize;
            public readonly int nTumorsExcluded;
            public readonly int nZero;
            public readonly int nOne;
            public readonly int nMore;
            public readonly int nSingleContributingToAltFraction;
            public readonly int nMultipleContributingToAltFraction;
        } // BonferroniCorrectedASEDistributionLine

        static public string RangeToDescriptiveString(string rangeString)
        {
            if ("0" == rangeString)
            {
                return "In gene";
            }

            bool exclusiveRange = rangeString.EndsWith("E");
            int numericRange;
            if (exclusiveRange)
            {
                numericRange = Convert.ToInt32(rangeString.Substring(0, rangeString.Length - 1));
            } else
            {
                numericRange = Convert.ToInt32(rangeString);
            }

            int rangeInBases = 1000 * (1 << numericRange);

            if (exclusiveRange)
            {
                return "" + SizeToUnits(((ulong)rangeInBases / 2)) + "bp-" + SizeToUnits((ulong)rangeInBases) + "bp";
            }

            return "<= " + SizeToUnits((ulong)rangeInBases) + "bp";
        } // RangeToDescriptiveString

        public static int RangeValueFromRange(string range)
        {
            if (range.EndsWith("E"))
            {
                return Convert.ToInt32(range.Substring(0, range.Count() - 1));
            }

            return Convert.ToInt32(range);
        }

        public class ExpressionDistribution
        {
            static List<string> wantedFields()
            {
                string[] wantedFieldsArray = { "Chromosome", "Locus", "min", "max" };
                var wantedFields = wantedFieldsArray.ToList();
                for (int i = 10; i < 100; i += 10)
                {
                    wantedFields.Add(i + "th %ile");
                }
                return wantedFields;
            }
            public static List<ExpressionDistribution> ReadFromFile(string filename)
            {
                var inputFile = CreateStreamReaderWithRetry(filename);

                if (null == inputFile)
                {
                    Console.WriteLine("ExpressionDistribution: unable to open input file " + filename);
                    return null;
                }

                var headerizedFile = new HeaderizedFile<ExpressionDistribution>(inputFile, false, true, "", wantedFields());

                List<ExpressionDistribution> result;
                headerizedFile.ParseFile(Parse, out result);

                return result;
            }

            public static bool ReadFromFile(string filename, HeaderizedFile<ExpressionDistribution>.ProcessNewItem processNewItem)
            {
                var inputFile = CreateStreamReaderWithRetry(filename);

                if (null == inputFile)
                {
                    Console.WriteLine("ExpressionDistribution: unable to open input file " + filename);
                    return false;
                }

                var headerizedFile = new HeaderizedFile<ExpressionDistribution>(inputFile, false, true, "", wantedFields());

                return headerizedFile.ParseFile(Parse, processNewItem);
            }

            static ExpressionDistribution Parse(HeaderizedFile<ExpressionDistribution>.FieldGrabber fieldGrabber)
            {
                var retVal = new ExpressionDistribution(fieldGrabber.AsString("Chromosome"), fieldGrabber.AsInt("Locus"), fieldGrabber.AsDouble("min"), fieldGrabber.AsDouble("max"));

                for (int i = 10; i < 100; i += 10)
                {
                    retVal.percentiles[i / 10 - 1] = fieldGrabber.AsDouble(i + "th %ile");
                }

                return retVal;
            }

            ExpressionDistribution(string chr_, int locus_, double min_, double max_)
            {
                chr = chr_;
                locus = locus_;
                min = min_;
                max = max_;
            }

            public readonly string chr;
            public readonly int locus;
            public readonly double min;
            public readonly double max;
            double[] percentiles = new double[9];

            public double getPercentile(int whichPercentile)
            {
                if (whichPercentile == 0)
                {
                    return min;
                }

                if (whichPercentile == 100)
                {
                    return max;
                }

                if (whichPercentile % 10 != 0 || whichPercentile < 10 || whichPercentile > 90)
                {
                    throw new InvalidParameterException();
                }

                return percentiles[whichPercentile / 10 - 1];
            }

        } // ExpressionDistribution

        public class ExpressionDistributionMap
        {
            void addToMap(ExpressionDistribution line)
            {
                if (!map.ContainsKey(line.chr))
                {
                    map.Add(line.chr, new Dictionary<int, ExpressionDistribution>());
                }

                map[line.chr].Add(line.locus, line);    // Assumes that each locus occurs at most once in the input.
            } // addToMap

            public ExpressionDistributionMap(string filename_)
            {
                filename = filename_;

                var worked = ExpressionDistribution.ReadFromFile(filename, x => this.addToMap(x));
                if (!worked)
                {
                    throw new Exception("Failed to load ExpressionDistributionMap.");
                }
            }

            public readonly string filename;
            public readonly Dictionary<string, Dictionary<int, ExpressionDistribution>> map = new Dictionary<string, Dictionary<int, ExpressionDistribution>>();
        } // ExpressionDistributionMap

#if false
        public class ExpressionDistributionMaps
        {
            public ExpressionDistributionMaps(Dictionary<string, Case> cases)
            {
                foreach (var disease in cases.Select(x => x.Value.disease()))
                {
                    if (!maps.ContainsKey(disease))
                    {
                        maps.Add(disease, new ExpressionDistributionMap(configuration.expression_distribution_directory + Expression_distribution_filename_base + disease));
                    }
                }
            }

            public readonly Dictionary<string, ExpressionDistributionMap> maps = new Dictionary<string, ExpressionDistributionMap>();   // Maps disease -> ExpressionDistributionMap
        } // ExpressionDistributionMaps
#endif  // If we're actually using this, we need to pass in the configuration

        public enum BinomalTestType
        {
            Less,
            Greater,
            TwoSided
        }


        //
        // binomialTest is adapted from the R binom.test code at https://github.com/SurajGupta/r-source/blob/master/src/library/stats/R/binom.test.R
        // which is under GPL2 or later. Copyright notice from that code follows:
        //

        /*
# File src/library/stats/R/binom.test.R
# Part of the R package, https://www.R-project.org
#
# Copyright (C) 1995-2012 The R Core Team
#
# This program is free software; you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation; either version 2 of the License, or
#  (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# A copy of the GNU General Public License is available at
# https://www.R-project.org/Licenses/
         */

        public static double binomialTest(int x, int n, double p, BinomalTestType alternative)
        {
            if (x < 0 || n < 1 || x > n || p < 0 || p > 1)
            {
                throw new InvalidParameterException();
            }

            try
            {
                switch (alternative)
                {
                    case BinomalTestType.Less:
                        return MathNet.Numerics.Distributions.Binomial.CDF(p, n, x);

                    case BinomalTestType.Greater:
                        return 1 - MathNet.Numerics.Distributions.Binomial.CDF(p, n, x - 1);

                    case BinomalTestType.TwoSided:
                        if (p == 0)
                        {
                            if (x == 0)
                            {
                                return 1;
                            } else
                            {
                                return 0;
                            }
                        }

                        if (p == 1)
                        {
                            if (x == n)
                            {
                                return 1;
                            } else
                            {
                                return 0;
                            }
                        }

                        /*
## Do
##   d <- dbinom(0 : n, n, p)
##   sum(d[d <= dbinom(x, n, p)])
## a bit more efficiently ...

## Note that we need a little fuzz.
                         */


                        double d = MathNet.Numerics.Distributions.Binomial.PMF(p, n, x);
                        double m = n * p;
                        if (x == m) // As if floating point multiplication will ever do this
                        {
                            return 1;
                        }

                        double relErr = 1 + 1e-7;

                        if (x < m)
                        {
                            int y = 0;
                            for (int i = (int)Math.Ceiling(m); i <= n; i++)
                            {
                                if (MathNet.Numerics.Distributions.Binomial.PMF(p, n, i) <= d * relErr)
                                {
                                    y++;
                                }
                            }

                            return MathNet.Numerics.Distributions.Binomial.CDF(p, n, x) + (1 - MathNet.Numerics.Distributions.Binomial.CDF(p, n, n - y));
                        }
                        else
                        {
                            int y = 0;
                            for (int i = 0; i <= Math.Floor(m); i++)
                            {
                                if (MathNet.Numerics.Distributions.Binomial.PMF(p, n, i) <= d * relErr)
                                {
                                    y++;
                                }
                            }

                            return MathNet.Numerics.Distributions.Binomial.CDF(p, n, y - 1) + (1 - MathNet.Numerics.Distributions.Binomial.CDF(p, n, x - 1));
                        }


                    default:
                        throw new InvalidParameterException();
                }  // switch
            } catch (Exception e)
            {
                Console.WriteLine("Binomial test: exception (probably in MathDotNet.Numerics), p = " + p + ", n = " + n + ", x = " + x + ".");
                throw e;
            }
        } // binomialTest

        // BED file format at http://uswest.ensembl.org/info/website/upload/bed.html
        public class BEDLine : IComparable
        {
            public readonly string chrom;
            public readonly int chromStart;
            public readonly int chromEnd;
            public readonly string name;
            public readonly int score;
            public readonly string strand;         // probably could be a char
            public readonly int thickStart;
            public readonly int thickEnd;
            public readonly string itemRGB;

            // The rest of the fields aren't filled in for the files we're using, so ignore them.

            public BEDLine(string chrom_, int chromStart_, int chromEnd_, string name_, int score_, string strand_, int thickStart_, int thickEnd_, string itemRGB_)
            {
                chrom = chrom_;
                chromStart = chromStart_;
                chromEnd = chromEnd_;
                name = name_;
                score = score_;
                strand = strand_;
                thickStart = thickStart_;
                thickEnd = thickEnd_;
                itemRGB = itemRGB_;
            }

            public BEDLine(string chrom_, int chromStart_)  // this is to make lookup keys for the AVL tree
            {
                chrom = chrom_;
                chromStart = chromStart_;
            }

            public static BEDLine fromLine(string line, bool stripChr)
            {
                var fields = line.Split('\t');

                string chrom;
                if (stripChr && fields[0].ToLower().StartsWith("chr"))
                {
                    chrom = fields[0].Substring(3);
                } else
                {
                    chrom = fields[0];
                }

                return new BEDLine(chrom, Convert.ToInt32(fields[1]), Convert.ToInt32(fields[2]), fields[3], Convert.ToInt32(fields[4]), fields[5], Convert.ToInt32(fields[6]), Convert.ToInt32(fields[7]), fields[8]);
            }

            public int CompareTo(object o_peer)
            {
                var peer = (BEDLine)o_peer;

                var stringComp = chrom.CompareTo(peer.chrom);
                if (stringComp != 0)
                {
                    return stringComp;
                }

                return chromStart.CompareTo(peer.chromStart);
            }
        }

        public class AnnotatedBEDLine : BEDLine
        {
            public AnnotatedBEDLine(BEDLine non_annotated) : base(non_annotated.chrom, non_annotated.chromStart, non_annotated.chromEnd, non_annotated.name, non_annotated.score, non_annotated.strand, non_annotated.thickStart, non_annotated.thickEnd, non_annotated.itemRGB)
            {
            }

            public int nBasesWithCoverage = 0;
            public int nMutations = 0;
            public int nMutationsBelow40Percent = 0;
            public int nMutationsAbove60Percent = 0;

            public int nMutationsWithModerateVAF()
            {
                return nMutations - nMutationsAbove60Percent - nMutationsBelow40Percent;
            }

            public double fractionCovered()
            {
                return (double)nBasesWithCoverage / (chromEnd - chromStart);
            }

            public double mutationsPerCoveredBase()
            {
                return (double)nMutations / (chromEnd - chromStart);

            }

            public static List<AnnotatedBEDLine> ReadFromFile(string filename)
            {
                string[] wantedFields = { "Chrom", "ChromStart", "ChromEnd", "name", "score", "strand", "thickStart", "thickEnd", "itemRGB", "nBasesWithCoverage", "nMutations", "nMutationsBelow40Percent", "nMutationsAbove60Percent" };

                var inputFile = CreateStreamReaderWithRetry(filename);
                if (null == inputFile)
                {
                    Console.WriteLine("AnnotatedBEDLine.ReadFromFile(): unable to open input file " + filename);
                    return null;
                }

                var headerizedFile = new HeaderizedFile<AnnotatedBEDLine>(inputFile, false, true, "", wantedFields.ToList());

                List<AnnotatedBEDLine> retVal;
                var worked = headerizedFile.ParseFile(ParseOneLine, out retVal);

                inputFile.Close();
                if (!worked)
                {
                    return null;
                }

                return retVal;
            }

            AnnotatedBEDLine(string chrom_, int chromStart_, int chromEnd_, string name_, int score_, string strand_, int thickStart_, int thickEnd_, string itemRGB_, int nBasesWithCoverage_, int nMutations_, int nMutationsBelow40Percent_, int nMutationsAbove60Percent_) :
                base(chrom_, chromStart_, chromEnd_, name_, score_, strand_, thickStart_, thickEnd_, itemRGB_)
            {
                nBasesWithCoverage = nBasesWithCoverage_;
                nMutations = nMutations_;
                nMutationsBelow40Percent = nMutationsBelow40Percent_;
                nMutationsAbove60Percent = nMutationsAbove60Percent_;
            }

            static AnnotatedBEDLine ParseOneLine(HeaderizedFile<AnnotatedBEDLine>.FieldGrabber fieldGrabber)
            {
                return new AnnotatedBEDLine(fieldGrabber.AsString("Chrom"), fieldGrabber.AsInt("ChromStart"), fieldGrabber.AsInt("ChromEnd"), fieldGrabber.AsString("name"), fieldGrabber.AsInt("score"), fieldGrabber.AsString("strand"),
                                            fieldGrabber.AsInt("thickStart"), fieldGrabber.AsInt("thickEnd"), fieldGrabber.AsString("itemRGB"), fieldGrabber.AsInt("nBasesWithCoverage"), fieldGrabber.AsInt("nMutations"),
                                            fieldGrabber.AsInt("nMutationsBelow40Percent"), fieldGrabber.AsInt("nMutationsAbove60Percent"));
            }
        } // AnnotatedBEDLine

        public class BEDFile
        {
            public static BEDFile ReadFromFile(string filename, bool stripChr)
            {
                //
                // This is not a headerized file, so we just break it into columns and parse it directly.
                //

                StreamReader inputFile;

                if (filename.ToLower().EndsWith(".gz"))
                {
                    inputFile = CreateCompressedStreamReaderWithRetry(filename);
                } else
                {
                    inputFile = CreateStreamReaderWithRetry(filename);
                }

                if (null == inputFile)
                {
                    Console.WriteLine("Unable to open BED file " + filename);
                    return null;
                }

                var bedFile = new BEDFile();

                string line;
                while (null != (line = inputFile.ReadLine()))
                {
                    var bedLine = BEDLine.fromLine(line, stripChr);
                    BEDLine containedLine;
                    if (bedFile.lines.Lookup(bedLine, out containedLine))
                    {
                        if (containedLine.chromEnd < bedLine.chromEnd)
                        {
                            bedFile.lines.Delete(containedLine);    // Keep the bigger one (there is one of these in the file).
                        }
                        else
                        {
                            continue;   // Skip the duplicate
                        }
                    }

                    bedFile.lines.Insert(bedLine);
                }

                inputFile.Close();
                return bedFile;
            }

            public bool isInRegulatoryRegion(string chrom, int locus)
            {
                var key = new BEDLine(chrom, locus);

                var treeEntry = lines.FindFirstLessThanOrEqualTo(key);

                if (treeEntry == null || treeEntry.chrom != chrom || treeEntry.chromEnd < locus)
                {
                    return false;
                }

                return true;
            }

            public void checkOverlap()
            {
                //
                // Make sure that no regions overlap with one another.
                //

                BEDLine prevLine = null;
                lines.findMin(out prevLine);
                BEDLine line;
                while (null != (line = lines.FindFirstGreaterThan(prevLine)))
                {
                    if (prevLine.chrom == line.chrom && prevLine.chromEnd > line.chromStart)
                    {
                        Console.WriteLine("Overlapping BED regions at " + prevLine.chrom + " " + prevLine.chromStart + "-" + prevLine.chromEnd + " and " + line.chromStart + "-" + line.chromEnd);
                    }
                    prevLine = line;
                }
            }

            public void ForEachLine(Action<BEDLine> action)
            {
                lines.ForEach(action);
            }

            AVLTree<BEDLine> lines = new AVLTree<BEDLine>();
        } // BEDFile


        public class BeatAMLSamplesSummaryLine
        {
            static public List<BeatAMLSamplesSummaryLine> ReadFile(string filename)
            {
                var inputFile = CreateStreamReaderWithRetry(filename);
                if (inputFile == null)
                {
                    return null;
                }

                string[] wantedFields = {"SeqID", "Original_LabID", "PatientID",  "Design_PatientID", "Diagnosis", "SpecificDiagnosis",  "Secondary_Specific_Diagnosis", "Tertiary_Specific_Diagnosis", "DiagnosisClass",
                                          "SpecimenType", "SampleGroup", "SpecimenSite", "CaptureGroup", "Sequential_CaptureGroup", "FlowCell", "Lane", "CoreProjectID", "CoreRunID",  "CoreFlowCellID", "CoreSampleID", "Outliers",
                                           "Specimen_HasRNASeq", "Patient_HasRNASeq"};

                var headerizedFile = new HeaderizedFile<BeatAMLSamplesSummaryLine>(inputFile, false, false, "", wantedFields.ToList());
                List<BeatAMLSamplesSummaryLine> retVal;
                if (!headerizedFile.ParseFile(ParseOneLine, out retVal))
                {
                    inputFile.Close();
                    return null;
                }

                inputFile.Close();
                return retVal;
            }

            static BeatAMLSamplesSummaryLine ParseOneLine(HeaderizedFile<BeatAMLSamplesSummaryLine>.FieldGrabber fieldGrabber)
            {
                return new BeatAMLSamplesSummaryLine(
                    fieldGrabber.AsString("SeqID"), fieldGrabber.AsString("Original_LabID"), fieldGrabber.AsInt("PatientID"), fieldGrabber.AsInt("Design_PatientID"), fieldGrabber.AsString("Diagnosis"),
                    fieldGrabber.AsString("SpecificDiagnosis"), fieldGrabber.AsString("Secondary_Specific_Diagnosis"), fieldGrabber.AsString("Tertiary_Specific_Diagnosis"), fieldGrabber.AsString("DiagnosisClass"),
                    fieldGrabber.AsString("SpecimenType"), fieldGrabber.AsString("SpecimenSite"), fieldGrabber.AsInt("SampleGroup"), fieldGrabber.AsString("CaptureGroup"), fieldGrabber.AsString("Sequential_CaptureGroup"),
                    fieldGrabber.AsString("FlowCell"), fieldGrabber.AsString("Lane"), fieldGrabber.AsString("CoreProjectID"), fieldGrabber.AsString("CoreRunID"), fieldGrabber.AsString("CoreFlowCellID"),
                    fieldGrabber.AsString("CoreSampleID"), fieldGrabber.AsString("Outliers"), fieldGrabber.AsBool("Specimen_HasRNASeq"), fieldGrabber.AsBool("Patient_HasRNASeq"));
            }

            BeatAMLSamplesSummaryLine(string seq_id_, string original_lab_id_, int patient_id_, int design_patient_id_, string diagnosis_, string specific_diagnosis_, string secondary_specific_diagnosis_,
                                      string tertiary_specific_diagnosis_, string diagnosis_class_, string specimen_type_, string specimen_site_, int sample_group_, string capture_group_, string sequential_capture_group_,
                                      string flow_cell_, string lane_, string core_projectId_, string core_runId_, string core_flow_cell_id_, string core_sample_id_, string outliers_, bool specimen_has_rna_seq_,
                                      bool patient_has_rna_seq_)
            {
                seq_id = seq_id_;
                original_lab_id = original_lab_id_;
                patient_id = patient_id_;
                design_patient_id = design_patient_id_;
                diagnosis = diagnosis_;
                specific_diagnosis = specific_diagnosis_;
                secondary_specific_diagnosis = secondary_specific_diagnosis_;
                tertiary_specific_diagnosis = tertiary_specific_diagnosis_;
                diagnosis_class = diagnosis_class_;
                specimen_type = specimen_type_;
                specimen_site = specimen_site_;
                sample_group = sample_group_;
                capture_group = capture_group_;
                sequential_capture_group = sequential_capture_group_;
                flow_cell = flow_cell_;
                lane = lane_;
                core_projectId = core_projectId_;
                core_runId = core_runId_;
                core_flow_cell_id = core_flow_cell_id_;
                core_sample_id = core_sample_id_;
                outliers = outliers_;
                specimen_has_rna_seq = specimen_has_rna_seq_;
                patient_has_rna_seq = patient_has_rna_seq_;
            }

            public readonly string seq_id;
            public readonly string original_lab_id;
            public readonly int patient_id;
            public readonly int design_patient_id;
            public readonly string diagnosis;
            public readonly string specific_diagnosis;
            public readonly string secondary_specific_diagnosis;
            public readonly string tertiary_specific_diagnosis;
            public readonly string diagnosis_class;
            public readonly string specimen_type;
            public readonly string specimen_site;
            public readonly int sample_group;
            public readonly string capture_group;
            public readonly string sequential_capture_group;
            public readonly string flow_cell;
            public readonly string lane;
            public readonly string core_projectId;
            public readonly string core_runId;
            public readonly string core_flow_cell_id;
            public readonly string core_sample_id;
            public readonly string outliers;
            public readonly bool specimen_has_rna_seq;
            public readonly bool patient_has_rna_seq;

        } // BeatAMLSamplesSummaryLine

        public class RegulatoryMutationsNearGene
        {
            public readonly string hugo_symbol;
            public readonly int nNonSilentMutations;
            public readonly int nNonSilentMutationsWithLowVAF;  // < 0.4
            public readonly int nNonSilentMutationsWithModerateVAF; // 0.4 <= VAF <= 0.6
            public readonly int nNonSilentMutationsWithHighVAF; // > 0.6
            public readonly int nSilentMutations;
            public readonly int nNotASECandidates;
            public readonly double[] mutationsPerCoveredBaseByRange;
            public readonly Dictionary<string, double> geneHancerMutationsPerBase;  // maps geneHancer ID -> mutations per base

            public static Dictionary<string, RegulatoryMutationsNearGene> readFile(string filename)
            {
                var inputFile = CreateStreamReaderWithRetry(filename);
                if (inputFile == null)
                {
                    Console.WriteLine("Unable to open input file " + filename);
                    return null;
                }

                var retVal = new Dictionary<string, RegulatoryMutationsNearGene>();

                string[] wantedStaticFields = { "Hugo_Symbol", "Non-Silent Mutation Count", "Non - Silent Mutations with Low VAF", "Non - Silent Mutations with Moderate VAF",
                    "Non - Silent Mutations with High VAF", "Silent Mutations", "Not ASE Candidates", "GeneHancers and frac mutated" };

                var wantedFields = wantedStaticFields.ToList();
                for (int i = 1; i < nRegions; i++)
                {
                    wantedFields.Add(regionIndexToString(i));
                }

                var headerizedFile = new HeaderizedFile<RegulatoryMutationsNearGene>(inputFile, false, true, "", wantedFields);
                List<RegulatoryMutationsNearGene> list;
                headerizedFile.ParseFile(x => Parse(filename, x), out list);

                if (null == list)
                {
                    Console.WriteLine("RegulatoryMutationsNearGene.readFile: failed to parse " + filename);
                    return null;
                }

                inputFile.Close();

                list.ForEach(x => retVal.Add(x.hugo_symbol, x));
                return retVal;
            } // readFile

            static RegulatoryMutationsNearGene Parse(string filename, HeaderizedFile<RegulatoryMutationsNearGene>.FieldGrabber fieldGrabber)
            {
                var mutationsPerCoveredBaseByRange = new double[nRegions];  // 0 isn't used, but it's still there because that way the address space matches everywhere else

                mutationsPerCoveredBaseByRange[0] = double.NegativeInfinity;

                for (int i = 1; i < nRegions; i++)
                {
                    mutationsPerCoveredBaseByRange[i] = fieldGrabber.AsDoubleNegativeInfinityIfStarOrEmptyString(regionIndexToString(i));
                }

                var geneHancerMutationsPerBase = new Dictionary<string, double>();
                var allGeneHancers = fieldGrabber.AsString("GeneHancers and frac mutated");
                if (allGeneHancers != "")
                {
                    foreach (var geneHancerString in allGeneHancers.Split(';'))
                    {
                        var fields = geneHancerString.Split('=');
                        if (fields.Count() != 2)
                        {
                            throw new Exception("GeneHancer field has the wrong number of things separated by an equal sign: " + geneHancerString + " in " + filename);
                        }

                        geneHancerMutationsPerBase.Add(fields[0], Convert.ToDouble(fields[1]));
                    }
                }


                return new RegulatoryMutationsNearGene(fieldGrabber.AsString("Hugo_Symbol"), fieldGrabber.AsInt("Non-Silent Mutation Count"), fieldGrabber.AsInt("Non - Silent Mutations with Low VAF"),
                    fieldGrabber.AsInt("Non - Silent Mutations with Moderate VAF"), fieldGrabber.AsInt("Non - Silent Mutations with High VAF"), fieldGrabber.AsInt("Silent Mutations"),
                    fieldGrabber.AsInt("Not ASE Candidates"), mutationsPerCoveredBaseByRange, geneHancerMutationsPerBase);
            } // Parse

            RegulatoryMutationsNearGene(string hugo_symbol_, int nNonSilentMutations_, int nNonSilentMutationsWithLowVAF_, int nNonSilentMutationsWithModerateVAF_,
                int nNonSilentMutationsWithHighVAF_, int nSilentMutations_, int nNotASECandiates_, double[] mutationsPerCoveredBaseByRange_, Dictionary<string, double> geneHancerMutationsPerBase_)
            {
                if (mutationsPerCoveredBaseByRange_.Count() != nRegions)
                {
                    throw new InvalidParameterException();
                }

                hugo_symbol = hugo_symbol_;
                nNonSilentMutations = nNonSilentMutations_;
                nNonSilentMutationsWithLowVAF = nNonSilentMutationsWithLowVAF_;
                nNonSilentMutationsWithModerateVAF = nNonSilentMutationsWithModerateVAF_;
                nNonSilentMutationsWithHighVAF = nNonSilentMutationsWithHighVAF_;
                nSilentMutations = nSilentMutations_;
                nNotASECandidates = nNotASECandiates_;
                mutationsPerCoveredBaseByRange = mutationsPerCoveredBaseByRange_;
                geneHancerMutationsPerBase = geneHancerMutationsPerBase_;
            } // RegulatoryMutationsNearGene ctor
        } // RegulatoryMutationsNearGene

        public static void PrintNumberBar(int maxValue)
        {
            if (maxValue >= 1000) return;

            if (maxValue >= 100)
            {
                for (int i = 1; i <= maxValue; i++)
                {
                    if (i % 100 == 0)
                    {
                        Console.Write((i / 100) % 10);
                    } else
                    {
                        Console.Write(" ");
                    }
                }
                Console.WriteLine();
            } // If we need a hundreds bar

            if (maxValue >= 10)
            {
                for (int i = 1; i <= maxValue; i++)
                {
                    if (i % 10 == 0)
                    {
                        Console.Write((i / 10) % 10);
                    } else
                    {
                        Console.Write(" ");
                    }
                }
                Console.WriteLine();
            } // If we need a tens bar

            for (int i = 1; i <= maxValue; i++)
            {
                Console.Write(i % 10);
            }
            Console.WriteLine();
        } // PrintNumberBar

        public static void PrintMessageAndNumberBar(string messageBeforeCount, string messageAfterCount, long count, out int nPerDot)
        {
            nPerDot = 1;
            int maxDots = 200;

            while (count / nPerDot > maxDots)
            {
                nPerDot *= 10;
            }

            Console.WriteLine(messageBeforeCount + " " + NumberWithCommas(count) + " " + messageAfterCount + " (one dot/" + NumberWithCommas(nPerDot) + ")");
            PrintNumberBar((int)(count / nPerDot));
        }

        public static string PadStringToGuidLength(string inputString)
        {
            if (inputString.Count() > GuidStringLength)
            {
                throw new Exception("PadStringToGuidLength: input string too long: " + inputString);
            }

            if (inputString.Count() == GuidStringLength)
            {
                return inputString;
            }

            var retVal = inputString;
            retVal += "-";
            while (retVal.Count() < GuidStringLength)
            {
                retVal += "x";
            }
            return retVal;
        } // PadStringToGuidLength

        public static string ExtractUnderlyingStringFromGuidPaddedString(string paddedString)
        {
            if (paddedString.Count() != GuidStringLength)
            {
                throw new Exception("ExtractUnderlyingStringFromGuidPaddedString : input string is the wrong length: " + paddedString);
            }

            int nToCut = 0;
            while (paddedString[GuidStringLength - nToCut - 1] == 'x')
            {
                nToCut++;
            }

            if (paddedString[GuidStringLength - nToCut - 1] == '-')
            {
                nToCut++;
            } else
            {
                if (nToCut != 0)
                {
                    throw new Exception("ExtractUnderlyingStringFromGuidPaddedString : input string is not a padded string: " + paddedString);
                }

                //
                // Input string was the right length without any padding.
                //
                return paddedString;
            }

            return paddedString.Substring(0, GuidStringLength - nToCut);
        } // ExtractUnderlyingStringFromGuidPaddedString

        public class GeneHancerLine : IComparable
        {
            public readonly string chromosome;
            public readonly int start;
            public readonly int end;
            public readonly string feature;
            public readonly double score;
            public readonly string attributes;
            public readonly string geneHancerId;
            public readonly Dictionary<string, double> geneScores = new Dictionary<string, double>();

            public GeneHancerLine(string chromosome_, int start_) // ctor for AVLTree keys
            {
                chromosome = chromosome_;
                start = start_;
            }

            public GeneHancerLine(string chromosome_, int start_, int end_, string feature_, double score_, string attributes_)
            {
                chromosome = chromosome_;
                start = start_;
                end = end_;
                feature = feature_;
                score = score_;
                attributes = attributes_;

                if (attributes != null)
                {
                    var attributeFields = attributes.Split(';');
                    int nAttributeFields = attributeFields.Count();
                    //
                    // The format is that each field is of the form name=value.  In practice, the first one is always
                    // genehancer_id, and then all the rest are connected_gene followed by score, but we'll parse
                    // it in a slightly more general way.
                    //
                    int whichField = 0;
                    while (whichField < nAttributeFields)
                    {
                        var subfields = attributeFields[whichField].Split('=');
                        if (subfields.Count() != 2)
                        {
                            throw new Exception("Unable to properly parse GeneHancer attribute field, subfield has the wrong count: " + attributes);
                        }

                        if (subfields[0] == "genehancer_id")
                        {
                            geneHancerId = subfields[1];
                            whichField++;
                        }
                        else if (subfields[0] == "connected_gene")
                        {
                            if (whichField + 1 == nAttributeFields)
                            {
                                throw new Exception("Unable to properly parse GeneHancer attribute field, connected_gene is the last field: " + attributes);
                            }

                            var score_subfields = attributeFields[whichField + 1].Split('=');
                            if (score_subfields.Count() != 2)
                            {
                                throw new Exception("Unable to properly parse GeneHancer attribute field, score has wrong subfield count: " + attributes);
                            }

                            if (score_subfields[0] != "score")
                            {
                                throw new Exception("Unable to properly parse GeneHancer attribute field, expected score field missing: " + attributes);
                            }

                            if (geneScores.ContainsKey(subfields[1]))
                            {
                                throw new Exception("Unable to properly parse GeneHancer attribute field, same gene present twice for a single entry: " + attributes);
                            }

                            geneScores.Add(subfields[1], Convert.ToDouble(score_subfields[1]));
                            whichField += 2;
                        }
                        else
                        {
                            throw new Exception("Unable to properly parse GeneHancer attribute field, subfield has the unknown type: " + attributes);
                        }
                    }
                }
            } // GeneHancerEntry ctor

            public int CompareTo(object peerObject)
            {
                GeneHancerLine peer = (GeneHancerLine)peerObject;

                if (peer.chromosome != chromosome)
                {
                    return chromosome.CompareTo(peer.chromosome);
                }

                return start.CompareTo(peer.start);
            }

            public static List<GeneHancerLine> ReadFromFile(string filename)
            {
                var inputFile = CreateStreamReaderWithRetry(filename);
                if (null == inputFile)
                {
                    return null;
                }

                string[] wantedFields = { "chrom", "feature name", "start", "end", "score", "attributes" };

                var headerizedFile = new HeaderizedFile<GeneHancerLine>(inputFile, false, false, "", wantedFields.ToList());

                List<GeneHancerLine> retVal;
                headerizedFile.ParseFile(Parse, out retVal);

                return retVal;
            } // ReadFromFile

            static GeneHancerLine Parse(HeaderizedFile<GeneHancerLine>.FieldGrabber fieldGrabber)
            {
                return new GeneHancerLine(fieldGrabber.AsString("chrom"), fieldGrabber.AsInt("start"), fieldGrabber.AsInt("end"), fieldGrabber.AsString("feature name"), fieldGrabber.AsDouble("score"), fieldGrabber.AsString("attributes"));
            }

            public static Dictionary<string, List<GeneHancerLine>> ReadFromFileToDict(string filename)
            {
                var geneHancers = ReadFromFile(filename);

                var geneHancersByGene = new Dictionary<string, List<GeneHancerLine>>();

                foreach (var geneHancer in geneHancers)
                {
                    foreach (var hugo_symbol in geneHancer.geneScores.Select(x => x.Key))
                    {
                        if (!geneHancersByGene.ContainsKey(hugo_symbol))
                        {
                            geneHancersByGene.Add(hugo_symbol, new List<GeneHancerLine>());
                        }
                        geneHancersByGene[hugo_symbol].Add(geneHancer);
                    }
                }

                return geneHancersByGene;
            }
        } // GeneHancerLine

        public class AnnotatedGeneHancerLine : GeneHancerLine
        {
            public int nBasesWithCoverage = 0;
            public int nMutations = 0;
            public int nMutationsBelow40Percent = 0;
            public int nMutationsAbove60Percent = 0;

            public AnnotatedGeneHancerLine(GeneHancerLine geneHancerLine) : base(geneHancerLine.chromosome, geneHancerLine.start, geneHancerLine.end, geneHancerLine.feature, geneHancerLine.score, geneHancerLine.attributes)
            {
            }

            public int nMutationsWithModerateVAF()
            {
                return nMutations - nMutationsAbove60Percent - nMutationsBelow40Percent;
            }

            public double fractionCovered()
            {
                return (double)nBasesWithCoverage / (end - start);
            }

            public double mutationsPerCoveredBase()
            {
                return (double)nMutations / (end - start);
            }

            public static new List<AnnotatedGeneHancerLine> ReadFromFile(string filename)
            {
                var inputFile = CreateStreamReaderWithRetry(filename);

                if (inputFile == null)
                {
                    return null;
                }

                string[] wantedFields = { "Chromosome", "Start", "End", "Feature", "Score", "Attributes", "n Bases With Coverage", "n Mutations", "n Mutations Below 40 Percent", "n Mutations Above 60 Percent" };

                var headerizedFile = new HeaderizedFile<AnnotatedGeneHancerLine>(inputFile, false, true, "", wantedFields.ToList());

                List<AnnotatedGeneHancerLine> retVal;

                headerizedFile.ParseFile(parse, out retVal);
                return retVal;
            }

            public static new Dictionary<string, List<AnnotatedGeneHancerLine>> ReadFromFileToDict(string filename) // This is an exact copy of the same function in the base class, but with different assumed types (from ReadFromFile).  There's got to be a better way to do this.
            {
                var geneHancers = ReadFromFile(filename);

                if (null == geneHancers)
                {
                    Console.WriteLine("AnnotatedGeneHancerLine.ReadFromFileToDict: failed to read/parse input file " + filename);
                }

                var geneHancersByGene = new Dictionary<string, List<AnnotatedGeneHancerLine>>();

                foreach (var geneHancer in geneHancers)
                {
                    foreach (var hugo_symbol in geneHancer.geneScores.Select(x => x.Key))
                    {
                        if (!geneHancersByGene.ContainsKey(hugo_symbol))
                        {
                            geneHancersByGene.Add(hugo_symbol, new List<AnnotatedGeneHancerLine>());
                        }
                        geneHancersByGene[hugo_symbol].Add(geneHancer);
                    }
                }

                return geneHancersByGene;
            }

            static AnnotatedGeneHancerLine parse(HeaderizedFile<AnnotatedGeneHancerLine>.FieldGrabber fieldGrabber)
            {
                var line = new AnnotatedGeneHancerLine(new GeneHancerLine(fieldGrabber.AsString("Chromosome"), fieldGrabber.AsInt("Start"), fieldGrabber.AsInt("End"), fieldGrabber.AsString("Feature"), fieldGrabber.AsDouble("Score"), fieldGrabber.AsString("Attributes")));

                line.nBasesWithCoverage = fieldGrabber.AsInt("n Bases With Coverage");
                line.nMutations = fieldGrabber.AsInt("n Mutations");
                line.nMutationsBelow40Percent = fieldGrabber.AsInt("n Mutations Below 40 Percent");
                line.nMutationsAbove60Percent = fieldGrabber.AsInt("n Mutations Above 60 Percent");

                return line;
            }
        } // AnnotatedGeneHancerLine

        public class CommonData // Stuff that's typically used by most of the apps, all rolled up into one place so I don't have to keep copying the initialization code into each new app
        {
            public readonly Dictionary<string, Case> cases;
            public readonly List<Case> listOfCases;
            public readonly Configuration configuration;
            public readonly GeneLocationsByNameAndChromosome geneLocationInformation;
            public readonly GeneMap geneMap;
            public readonly Dictionary<string, ASETools.ASEMapPerGeneLine> perGeneASEMap;
            public readonly ASECorrection aseCorrection;
            public ExpressionDistributionByChromosomeMap expressionDistributionByChromosomeMap;
            public List<string> diseases;
            public Stopwatch timer; // Starts running as soon as this is created
            public Dictionary<string, TCGAClinicalSummaryLine> clinicalSummariesByPatientId;
            // If you add new things here, be sure to add them to DetermineStateOfTheWorld() in ASEProcessManager.

            public List<string> GetInputFiles()
            {
                var retVal = new List<string>();

                retVal.AddRange(ASETools.TCGAClinicalSummaryLine.getInputFilenames(configuration));
                retVal.Add(configuration.casesFilePathname);
                retVal.Add(configuration.finalResultsDirectory + ASETools.ASECorrectionFilename);
                retVal.Add(configuration.geneLocationInformationFilename);
                retVal.Add(configuration.finalResultsDirectory + ASETools.PerGeneASEMapFilename);

                return retVal;
            }

            public static CommonData LoadCommonData(string[] args)
            {
                var timer = new Stopwatch();
                timer.Start();

                var configuration = ASETools.Configuration.loadFromFile(args);
                if (null == configuration)
                {
                    Console.WriteLine("Unable to load configuration.");
                    return null;
                }

                Dictionary<string, TCGAClinicalSummaryLine> clinicalSummariesByPatientId;
                if (configuration.isBeatAML)
                {
                    clinicalSummariesByPatientId = null;
                } else
                {
                    clinicalSummariesByPatientId = TCGAClinicalSummaryLine.readAllToDictionary(configuration);
                }

                var cases = ASETools.Case.LoadCases(configuration.casesFilePathname);

                if (null == cases)
                {
                    Console.WriteLine("Unable to load cases.");
                    return null;
                }

                var aseCorrection = ASETools.ASECorrection.LoadFromFile(configuration.finalResultsDirectory + ASETools.ASECorrectionFilename);
                if (null == aseCorrection)
                {
                    Console.WriteLine("Unable to load ASE correction");
                    return null;
                }

                var geneLocationInformation = new ASETools.GeneLocationsByNameAndChromosome(ASETools.readKnownGeneFile(configuration.geneLocationInformationFilename));
                var geneMap = new ASETools.GeneMap(geneLocationInformation.genesByName);

                var perGeneASEMap = ASETools.ASEMapPerGeneLine.ReadFromFileToDictionary(configuration.finalResultsDirectory + ASETools.PerGeneASEMapFilename);

                if (null == perGeneASEMap)
                {
                    Console.WriteLine("You must first create the per-gene ASE map in " + configuration.finalResultsDirectory + ASETools.PerGeneASEMapFilename);
                    return null;
                }

                var diseases = new List<string>();

                foreach (var caseEntry in cases)
                {
                    var case_ = caseEntry.Value;

                    if (!diseases.Contains(case_.disease()))
                    {
                        diseases.Add(case_.disease());
                    }
                } // case

                ExpressionDistributionByChromosomeMap expressionDistributionByChromosomeMap;

                if (configuration.expression_distribution_by_chromosome_map_filename == "")
                {
                    expressionDistributionByChromosomeMap = null;
                } else
                {
                    expressionDistributionByChromosomeMap = ASETools.ExpressionDistributionByChromosomeMap.LoadFromFile(configuration.expression_distribution_by_chromosome_map_filename);
                }

                return new CommonData(cases, configuration, geneLocationInformation, geneMap, perGeneASEMap, aseCorrection, diseases, timer, expressionDistributionByChromosomeMap, clinicalSummariesByPatientId);
            } // LoadCommonData

            CommonData(Dictionary<string, Case> cases_, Configuration configuration_, GeneLocationsByNameAndChromosome geneLocationInformation_, GeneMap geneMap_,
                Dictionary<string, ASETools.ASEMapPerGeneLine> perGeneASEMap_, ASECorrection aseCorrection_, List<string> diseases_, Stopwatch timer_, ExpressionDistributionByChromosomeMap expressionDistributionByChromosomeMap_,
                Dictionary<string, TCGAClinicalSummaryLine> clinicalSummariesByPatientId_)
            {
                cases = cases_;
                listOfCases = cases.Select(x => x.Value).ToList();
                configuration = configuration_;
                geneLocationInformation = geneLocationInformation_;
                geneMap = geneMap_;
                perGeneASEMap = perGeneASEMap_;
                aseCorrection = aseCorrection_;
                diseases = diseases_;
                timer = timer_;
                expressionDistributionByChromosomeMap = expressionDistributionByChromosomeMap_;
                clinicalSummariesByPatientId = clinicalSummariesByPatientId_;
            }
        } // CommonData

        public class BasesInCodingAndKnownExpressionRegions
        {
            class GeneAndBaseCount
            {
                public readonly string hugo_symbol;
                public readonly int baseCount;

                public GeneAndBaseCount(string hugo_symbol_, int baseCount_)
                {
                    hugo_symbol = hugo_symbol_;
                    baseCount = baseCount_;
                }
            } // GeneAndBaseCount

            public bool isGeneKnown(string hugo_symbol)
            {
                return basesPerGene.ContainsKey(hugo_symbol);
            }

            public int basesForGene(string hugo_symbol)
            {
                return basesPerGene[hugo_symbol];
            }

            public static BasesInCodingAndKnownExpressionRegions LoadFromFile(string filename)
            {
                var inputFile = CreateStreamReaderWithRetry(filename);
                if (null == inputFile)
                {
                    Console.WriteLine("Unable to open file " + filename);
                    return null;
                }

                string[] wantedFields = { "Hugo Symbol", "Bases In Coding And Known Expression Region" };

                var parser = new HeaderizedFile<GeneAndBaseCount>(inputFile, false, true, "", wantedFields.ToList());

                List<GeneAndBaseCount> genesAndBases;
                if (!parser.ParseFile(ParseOne, out genesAndBases))
                {
                    Console.WriteLine("Unable to parse input file " + filename);
                    inputFile.Close();
                }

                var retVal = new BasesInCodingAndKnownExpressionRegions();

                genesAndBases.ForEach(x => retVal.basesPerGene.Add(x.hugo_symbol, x.baseCount));

                return retVal;
            }

            static GeneAndBaseCount ParseOne(HeaderizedFile<GeneAndBaseCount>.FieldGrabber fieldGrabber)
            {
                return new GeneAndBaseCount(fieldGrabber.AsString("Hugo Symbol"), fieldGrabber.AsInt("Bases In Coding And Known Expression Region"));
            }

            Dictionary<string, int> basesPerGene = new Dictionary<string, int>();
        } // BasesInCodingAndKnownExpressionRegions

        public readonly static int[] ZeroOneTwo = { 0, 1, 2 };

        public class ExpressionByGeneLine
        {
            public static Dictionary<string, ExpressionByGeneLine> ReadFromFile(string filename)
            {
                var inputFile = ASETools.CreateStreamReaderWithRetry(filename);
                if (null == inputFile)
                {
                    return null;
                }

                string[] wantedFields = { "Hugo Symbol", "Fraction of mean expression" };
                var headerizedFile = new HeaderizedFile<ExpressionByGeneLine>(inputFile, false, true, "", wantedFields.ToList());

                List<ExpressionByGeneLine> listOfResults;
                headerizedFile.ParseFile(parser, out listOfResults);

                return listOfResults.GroupByToDictUnique(x => x.hugo_symbol);
            }

            static ExpressionByGeneLine parser(ASETools.HeaderizedFile<ExpressionByGeneLine>.FieldGrabber fieldGrabber)
            {
                return new ExpressionByGeneLine(fieldGrabber.AsString("Hugo Symbol"), fieldGrabber.AsDoubleNegativeInfinityIfStarOrEmptyString("Fraction of mean expression"));
            }

            ExpressionByGeneLine(string hugo_symbol_, double fractionOfMeanExpression_)
            {
                hugo_symbol = hugo_symbol_;
                fractionOfMeanExpression = fractionOfMeanExpression_;
            }

            public readonly string hugo_symbol;
            public readonly double fractionOfMeanExpression;
        } // ExpressionByGeneLine

        public class ThrottledParallelQueue<TElement>
        {
            public ThrottledParallelQueue(int maxSize_, int nInitialWriters)
            {
                maxSize = maxSize_;
                nRemainingWriters = nInitialWriters;
            }

            //
            // Returns false iff the queue is terminated.
            //
            public bool Dequeue(out TElement element)
            {
                while (true)
                {
                    lock (this)
                    {
                        var size = queue.Count();

                        if (0 == size)
                        {
                            if (nRemainingWriters == 0)
                            {
                                element = default(TElement);
                                return false;
                            }

                            queueNotEmptyEvent.Reset();
                        }
                        else
                        {
                            element = queue[0];
                            queue.RemoveAt(0);
                            queueNotFullEvent.Set();

                            if (1 == size)
                            {
                                queueNotEmptyEvent.Reset();
                            }

                            return true;
                        }
                    } // lock

                    queueNotEmptyEvent.WaitOne();
                } // while(true) 
            } // Dequeue

            public void Enqueue(TElement element)
            {
                while (true)
                {
                    lock (this)
                    {
                        int size = queue.Count();

                        if (size >= maxSize)
                        {
                            queueNotFullEvent.Reset();
                        }
                        else
                        {
                            queue.Add(element);
                            if (size + 1 >= maxSize)
                            {
                                queueNotFullEvent.Reset();
                            }

                            queueNotEmptyEvent.Set();
                            return;
                        }
                    } // lock

                    queueNotFullEvent.WaitOne();
                } // while (true)
            } // Enqueue

            public void TerminateWriter()
            {
                lock (this)
                {
                    if (nRemainingWriters < 1)
                    {
                        throw new Exception("ThrottledParallelQueue: terminating writer when there are no remaining writers to terminate.");
                    } else
                    {
                        nRemainingWriters--;
                        if (nRemainingWriters == 0)
                        {
                            queueNotEmptyEvent.Set();   // Wake the waiters so they can exit.
                        }
                    }
                }
            } // TerminateWriter

            public void AddWriter()
            {
                lock (this)
                {
                    nRemainingWriters++;
                }
            } // AddWriter



            int maxSize;
            int nRemainingWriters;
            List<TElement> queue = new List<TElement>();
            ManualResetEvent queueNotFullEvent = new ManualResetEvent(true);
            ManualResetEvent queueNotEmptyEvent = new ManualResetEvent(false);

        } // ThrottledParallelQueue

        public class ExpressionByGene
        {
            class HugoSymbolAndExpression
            {
                public HugoSymbolAndExpression(string hugo_symbol_, double expression_)
                {
                    hugo_symbol = hugo_symbol_;
                    expression = expression_;
                }

                public string hugo_symbol;
                public double expression;
            }

            static public ExpressionByGene LoadFromFile(string filename)
            {
                string[] wantedFields = { "Hugo Symbol", "Fraction of mean expression" };

                var inputFile = CreateStreamReaderWithRetry(filename);
                if (null == inputFile)
                {
                    Console.WriteLine("Unable to open " + filename);
                    return null;
                }

                var headerizedFile = new HeaderizedFile<HugoSymbolAndExpression>(inputFile, false, true, "", wantedFields.ToList());

                List<HugoSymbolAndExpression> values;
                if (!headerizedFile.ParseFile(parse, out values))
                {
                    Console.WriteLine("Unable to parse " + filename);
                    return null;
                }

                inputFile.Close();

                var retVal = new ExpressionByGene();
                values.ForEach(_ => retVal.fractionOfAllReads.Add(_.hugo_symbol, _.expression));

                return retVal;
            }

            static HugoSymbolAndExpression parse(HeaderizedFile<HugoSymbolAndExpression>.FieldGrabber fieldGrabber)
            {
                return new HugoSymbolAndExpression(fieldGrabber.AsString("Hugo Symbol"), fieldGrabber.AsDouble("Fraction of mean expression"));
            }

            public double getMeanExpressionForGene(string hugo_symbol)
            {
                return fractionOfAllReads[hugo_symbol];
            }

            public List<string> GetAllHugoSymbols()
            {
                return fractionOfAllReads.Select(_ => _.Key).ToList();
            }

            Dictionary<string, double> fractionOfAllReads = new Dictionary<string, double>();
        } // ExpressionByGene

        // These are the genes from figure S11(a) of "Genomic and Epigenomic Landscapes of Adult De Novo Acute Myeloid Leukemia," N Engl J Med 2013; 368:2059-2074;  https://www.nejm.org/doi/suppl/10.1056/NEJMoa1301689/suppl_file/nejmoa1301689_appendix.pdf
        public static string[] spliceosome_genes = { "CSTF2T", "DDX1", "DDX23", "DHX32", "HNRNPK", "METTL3", "PLRG1", "POLR2A", "PRPF3", "PRPF8", "RBMX", "SF3B1", "SNRNP200", "SRRM2", "SRSF2", "SRSF6", "SUPT5H", "TRA2B", "U2AF1", "U2AF1L4", "U2AF2" };

        public static string NumberWithCommas(long value)
        {
            return String.Format("{0:n0}", value);
        }

        public static void ComputeConfidenceInterval(IEnumerable<double> enumerableValues, double desiredConfidence, out double mean, out double range)
        {
            //
            // Implemented from the Wikipedia page, https://en.wikipedia.org/wiki/Confidence_interval, "Theoretical Example."   This assumes a normal distribution of the underlying distribution, but
            // for n > "a few dozen" it's "quite good," unless the distribution's "CDF does not have any discontinuities and its skewness is moderate."
            //

            var values = enumerableValues.ToList();
            int n = values.Count();

            var meanAndStdDev = new RunningMeanAndStdDev();
            values.ForEach(_ => meanAndStdDev.addValue(_));

            var c = MathNet.Numerics.Distributions.StudentT.InvCDF(0, 1, n - 1, desiredConfidence / 2);

            mean = meanAndStdDev.getMeanAndStdDev().mean;
            range = c * meanAndStdDev.getMeanAndStdDev().stddev / Math.Sqrt(n);
        } // ComputeConfidenceInterval

        public class IsoformReadCounts
        {
            public IsoformReadCounts(string ucsdId_, int tumor_, int normal_)
            {
                ucsdId = ucsdId_;
                tumor = tumor_;
                normal = normal_;
            }

            public int getCount(bool forTumor)
            {
                if (forTumor)
                {
                    return tumor;
                }

                return normal;
            }

            public static Dictionary<string, IsoformReadCounts> LoadFromFile(string filename)
            {
                var inputFile = CreateStreamReaderWithRetry(filename);
                if (null == inputFile)
                {
                    Console.WriteLine("Unable to open input file " + filename);
                    return null;
                }

                string[] wantedFields = { "Uscd ID", "Normal Read Count", "Tumor Read Count" };

                var headerizedFile = new HeaderizedFile<IsoformReadCounts>(inputFile, false, true, "", wantedFields.ToList());
                List<IsoformReadCounts> listOfResults;
                if (!headerizedFile.ParseFile(parse, out listOfResults))
                {
                    Console.WriteLine("Unable to parse " + filename);
                    return null;
                }

                var retVal = new Dictionary<string, IsoformReadCounts>();
                listOfResults.ForEach(_ => retVal.Add(_.ucsdId, _));

                return retVal;
            } // LoadFromFile

            static IsoformReadCounts parse(HeaderizedFile<IsoformReadCounts>.FieldGrabber fieldGrabber)
            {
                return new IsoformReadCounts(fieldGrabber.AsString("Uscd ID"), fieldGrabber.AsInt("Tumor Read Count"), fieldGrabber.AsInt("Normal Read Count"));
            }

            public static void WriteToFile(string filename, List<IsoformReadCounts> readCounts)
            {
                var outputFile = CreateStreamWriterWithRetry(filename);
                if (null == outputFile)
                {
                    Console.WriteLine("Unable to open output file " + filename);
                    return;
                }

                outputFile.WriteLine("Uscd ID\tNormal Read Count\tTumor Read Count");

                foreach (var readCount in readCounts)
                {
                    outputFile.WriteLine(readCount.ucsdId + "\t" + readCount.normal + "\t" + readCount.tumor);
                }

                outputFile.WriteLine("**done**");
                outputFile.Close();
            } // WriteToFile


            public readonly string ucsdId;
            public readonly int tumor, normal;
        } // IsoformReadCounts

        public static char[] GeneticBases = { 'A', 'T', 'C', 'G' };
        static public char GeneticCode(string codon)
        {
            codon = codon.ToUpper();

            if (codon.Length != 3 || codon.ToArray().Any(_ => !GeneticBases.Contains(_)))
            {
                throw new Exception("GeneticCode: input isn't a codon: " + codon);
            }

            switch (codon)
            {
                case "AAA":
                case "AAG":
                    return 'K';
                case "AAC":
                case "AAT":
                    return 'N';
                case "ACA":
                case "ACT":
                case "ACG":
                case "ACC":
                    return 'T';
                case "AGA":
                case "AGG":
                    return 'R';
                case "AGC":
                case "AGT":
                    return 'S';
                case "ATA":
                case "ATC":
                case "ATT":
                    return 'I';
                case "ATG":
                    return 'M';
                case "CAA":
                case "CAG":
                    return 'Q';
                case "CAT":
                case "CAC":
                    return 'H';
                case "CCA":
                case "CCC":
                case "CCG":
                case "CCT":
                    return 'P';
                case "CGA":
                case "CGC":
                case "CGG":
                case "CGT":
                    return 'R';
                case "CTA":
                case "CTC":
                case "CTG":
                case "CTT":
                    return 'L';
                case "GAA":
                case "GAG":
                    return 'E';
                case "GAC":
                case "GAT":
                    return 'D';
                case "GCA":
                case "GCC":
                case "GCG":
                case "GCT":
                    return 'A';
                case "GGA":
                case "GGC":
                case "GGG":
                case "GGT":
                    return 'G';
                case "GTA":
                case "GTC":
                case "GTG":
                case "GTT":
                    return 'V';
                case "TAA":
                case "TAG":
                    return '*';
                case "TAC":
                case "TAT":
                    return 'Y';
                case "TCA":
                case "TCC":
                case "TCG":
                case "TCT":
                    return 'S';
                case "TGA":
                    return '*';
                case "TGC":
                case "TGT":
                    return 'C';
                case "TGG":
                    return 'W';
                case "TTA":
                case "TTG":
                    return 'L';
                case "TTC":
                case "TTT":
                    return 'F';

            }

            throw new Exception("GeneticCode: fell out bottom of switch, which shouldn't be possible.");
        } // GeneticCode

        public static string ReverseCompliment(string input)
        {
            string output = "";

            for (int i = 0; i < input.Length; i++)
            {
                switch (input[i])
                {
                    case 'A':
                    case 'a':
                        output = "T" + output;
                        break;

                    case 'T':
                    case 't':
                        output = "A" + output;
                        break;

                    case 'C':
                    case 'c':
                        output = "G" + output;
                        break;

                    case 'G':
                    case 'g':
                        output = "C" + output;
                        break;

                    case 'N':
                    case 'n':
                        output = "N" + output;
                        break;

                    default:
                        throw new Exception("ReverseCompliment: input contains non-base (or N): " + input);
                } // switch
            } // for

            return output;
        }  // ReverseCompliment

        //
        // Generate the friendly description of an IDH1 R132 mutation
        //
        static public string IDH1MutantDescription(int locus, char newBase)
        {
            if (locus < 208248387 || locus > 208248389)
            {
                return "Other";
            }

            StringBuilder codon = new StringBuilder("AGC");
            codon[locus - 208248387] = newBase;
            return "R132" + GeneticCode(ReverseCompliment(codon.ToString()));
        } // IDH1MutantDescription

        static public int ratioToPercentage(int numerator, int denominator)
        {
            return (int)(100 * (double)numerator / denominator);
        }

        public class BAMMetadata
        {
            public bool isPaired;
            public int minReadLength;
            public int maxReadLength;
            public double meanReadLength;
            public double medianReadLength;
            public int minInsert;
            public int maxInsert;
            public double meanInsert;
            public double medianInsert;

            public BAMMetadata(bool isPaired_, int minReadLength_, int maxReadLength_, double meanReadLength_, double medianReadLength_, int minInsert_, int maxInsert_, double meanInsert_, double medianInsert_)
            {
                isPaired = isPaired_;
                minReadLength = minReadLength_;
                maxReadLength = maxReadLength_;
                meanReadLength = meanReadLength_;
                medianReadLength = medianReadLength_;
                minInsert = minInsert_;
                maxInsert = maxInsert_;
                meanInsert = meanInsert_;
                medianInsert = medianInsert_;
            }
        } // BAMMetadata

        public class CaseMetadata
        {
            public readonly string caseID;
            Dictionary<bool, Dictionary<bool, BAMMetadata>> bamMetadata = null; // tumor->dna=>metadata
            Dictionary<bool, Dictionary<bool, string>> sample = null; // tumor->dna->sample ID

            public CaseMetadata(string caseID_, Dictionary<bool, Dictionary<bool, string>> sample_, Dictionary<bool, Dictionary<bool, BAMMetadata>> bamMetadata_)
            {
                caseID = caseID_;
                sample = sample_;
                bamMetadata = bamMetadata_;
            }

            public BAMMetadata getBAMMetadata(bool tumor, bool dna)
            {
                return bamMetadata[tumor][dna];
            }

            public string getSample(bool tumor, bool dna)
            {
                return sample[tumor][dna];
            }

            public static CaseMetadata ReadFromFile(string filename)
            {
                var inputFile = CreateStreamReaderWithRetry(filename);

                if (null == inputFile)
                {
                    throw new Exception("CaseMetadata.ReadFromFile: unable to open " + filename);
                }

                List<string> neededFields = new List<string>();
                neededFields.Add("Case ID");

                foreach (var tumor in BothBools)
                {
                    foreach (var dna in BothBools)
                    {
                        string specifier = getSpecifier(tumor, dna);

                        neededFields.Add(specifier + "is paired");
                        neededFields.Add(specifier + "min read length");
                        neededFields.Add(specifier + "max read length");
                        neededFields.Add(specifier + "mean read length");
                        neededFields.Add(specifier + "min insert");
                        neededFields.Add(specifier + "max insert");
                        neededFields.Add(specifier + "median insert");
                        neededFields.Add(specifier + "mean insert");
                        neededFields.Add(specifier + "median read length");

                        neededFields.Add(specifier + "Sample");
                    } // dna

                } // tumor

                var headerizedFile = new HeaderizedFile<CaseMetadata>(inputFile, false, true, "", neededFields);
                List<CaseMetadata> parsedMetadata;
                headerizedFile.ParseFile(parser, out parsedMetadata);

                inputFile.Close();

                if (parsedMetadata.Count() != 1)
                {
                    throw new Exception("CaseMetadata.ReadFromFile: file " + filename + " has wrong line count: " + parsedMetadata.Count());
                }

                return parsedMetadata[0];
            } // ReadFromFile

            public static string getSpecifier(bool tumor, bool dna)
            {
                string specifier;

                if (tumor)
                {
                    specifier = "Tumor ";
                }
                else
                {
                    specifier = "Normal ";
                }

                if (dna)
                {
                    specifier += "DNA ";
                }
                else
                {
                    specifier += "RNA ";
                }

                return specifier;
            }

            static CaseMetadata parser(HeaderizedFile<CaseMetadata>.FieldGrabber fieldGrabber)
            {
                Dictionary<bool, Dictionary<bool, BAMMetadata>> bamMetadata = new Dictionary<bool, Dictionary<bool, BAMMetadata>>();
                Dictionary<bool, Dictionary<bool, string>> sample = new Dictionary<bool, Dictionary<bool, string>>();

                foreach (var tumor in BothBools)
                {
                    bamMetadata.Add(tumor, new Dictionary<bool, BAMMetadata>());
                    sample.Add(tumor, new Dictionary<bool, string>());

                    foreach (var dna in BothBools)
                    {
                        string specifier = getSpecifier(tumor, dna);

                        bamMetadata[tumor].Add(dna, new BAMMetadata(
                                                        fieldGrabber.AsBool(specifier + "is paired"),
                                                        fieldGrabber.AsInt(specifier + "min read length"),
                                                        fieldGrabber.AsInt(specifier + "max read length"),
                                                        fieldGrabber.AsDoubleNegativeInfinityIfStarOrEmptyString(specifier + "mean read length"),
                                                        fieldGrabber.AsDoubleNegativeInfinityIfStarOrEmptyString(specifier + "median read length"),
                                                        fieldGrabber.AsInt(specifier + "min insert"),
                                                        fieldGrabber.AsInt(specifier + "max insert"),
                                                        fieldGrabber.AsDoubleNegativeInfinityIfStarOrEmptyString(specifier + "mean insert"),
                                                        fieldGrabber.AsDoubleNegativeInfinityIfStarOrEmptyString(specifier + "median insert")));

                        if (bamMetadata[tumor][dna].minReadLength == -1)
                        {
                            //
                            // This is a missing normal RNA sample.
                            //
                            if (tumor || dna)
                            {
                                throw new Exception("CaseMetadata.parser: other than normal RNA has missing read length");
                            }
                            bamMetadata[tumor][dna] = null;
                        }

                        sample[tumor].Add(dna, fieldGrabber.AsString(specifier + "Sample"));
                    } // dna
                } // tumor

                return new CaseMetadata(fieldGrabber.AsString("Case ID"), sample, bamMetadata);
            } // parser
        }// CaseMetadata

        public static void RunAndWaitForProcess(string binaryName, string commandLine)
        {
            var startInfo = new ProcessStartInfo(binaryName, commandLine);

            startInfo.UseShellExecute = false;

            Process process;
            try
            {
                process = Process.Start(startInfo);
            }
            catch (Exception e)
            {
                Console.WriteLine("Error trying to start process cmd.exe");
                Console.WriteLine("Exception message: " + e.Message);

                throw e;
            }

            process.WaitForExit();
        } // RunAndWaitForProcess
        public static List<string> RunProcessAndGetOutput(string binaryName, string commandLine)
        {
            var output = new List<string>();

            var startInfo = new ProcessStartInfo(binaryName, commandLine);

            startInfo.UseShellExecute = false;
            startInfo.RedirectStandardOutput = true;

            Process process;
            try
            {
                process = Process.Start(startInfo);
            }
            catch (Exception e)
            {
                Console.WriteLine("Error trying to start process cmd.exe");
                Console.WriteLine("Exception message: " + e.Message);

                throw e;
            }

            string outputLine;
            while (null != (outputLine = process.StandardOutput.ReadLine()))
            {
                output.Add(outputLine);
            }

            process.WaitForExit();

            return output;
        } // RunProcessAndGetOutput

#if false
        public static string ParsePossiblyQuotedString(string inputString, string terminatingCharacters, out int nCharacters)
        {
            //
            // Trim leading spaces.
            //
            nCharacters = 0;
            while (nCharacters < inputString.Length && inputString[nCharacters] == ' ')
            {
                nCharacters++;
            }

            if (inputString.Substring(nCharacters) == "")
            {
                return "";
            }

            var retVal = "";
            bool isQuoted = inputString[nCharacters] == '\"';
            if (isQuoted)
            {
                nCharacters++;
                int beginNonSpace = nCharacters;
                while (nCharacters < inputString.Length)
                {
                    if (inputString[nCharacters] == '\"')
                    {
                        nCharacters++;
                        return inputString.Substring()
                    }
                }
            }

        }

        public static string ParsePossiblyQuotedString(string inputString, string terminatingCharacters = "")
        {
            int nCharacters;
            return ParsePossiblyQuotedString(inputString, terminatingCharacters, out nCharacters);
        }

        public class VCFFile
        {
            public static VCFFile ReadFromFile(string filename)
            {

            }

            public class Info
            {
                public readonly string ID;
                public readonly string number;

                public enum Type { Integer, Float, Flag, Character, String};
                public readonly Type type;

                public readonly string description;

                class MetaInformationField
                {
                    public readonly string key;
                    public readonly string value;
                    public readonly int charsConsumed;  // How many characters of the input line were used in this field?

                    public MetaInformationField(string inputLine)
                    {
                        int firstEqual = inputLine.IndexOf('=');
                        if (firstEqual == -1)
                        {
                            throw new Exception("MetaInformationField: input line didn't contain an =: " + inputLine);
                        }

                        key = inputLine.Substring(0, firstEqual);


                    }
                }


                class MetaInformationLine
                {
                    public MetaInformationLine(string line)
                    {
                        if (!line.StartsWith("##") || !line.Contains("=") || line.Length == line.IndexOf('=') || line[line.IndexOf('=') + 1] == '<' && !line.EndsWith(">"))
                        {
                            throw new Exception("ASETools.VCFile.MetaInformationLine constructor: line is not properly formatted: " + line);
                        }

                        var indexOfEqualSign = line.IndexOf('=');
                        key = line.Substring(2, indexOfEqualSign - 2);
                        if (line[indexOfEqualSign+1] != '<')
                        {
                            //
                            // Single value line.
                            //
                            values.Add("", line.Substring(indexOfEqualSign + 1));
                            return;
                        }

                        var lineGuts = line.Substring(indexOfEqualSign + 2, line.Length - indexOfEqualSign - 3);

                        while (lineGuts.Length != 0)
                        {
                            var fields = lineGuts.Split('=');
                            
                        }
                    }

                    public readonly Dictionary<string, string> values = new Dictionary<string, string>();   // If there is only one value, then it's at the empty string key.
                }

                public Info(string infoLine)
                {
                    if (!infoLine.ToLower().StartsWith("##info=<") || !infoLine.EndsWith(">"))
                    {
                        throw new Exception("ASETools.VCFFile.Info constructor called with a non-info line");
                    }

                    string guts = infoLine.Substring(8, infoLine.Length - 9);
                    var fields = guts.Split(',');

                    var idFields = fields.Where(_ => _.ToLower().StartsWith("id=")).ToList();
                    var numberFields = fields.Where(_ => _.ToLower().StartsWith("number=")).ToList();
                    var typeFields = fields.Where(_ => _.ToLower().StartsWith("type=")).ToList();
                    var descriptionFields = fields.Where(_ => _.ToLower().StartsWith("description=\"") && _.EndsWith("\"")).ToList();

                    if (idFields.Count() != 1 || numberFields.Count() != 1 || typeFields.Count() != 1 || descriptionFields.Count() != 1)
                    {
                        throw new Exception("ASETools.VCFFIle.Info constructor: no (or more than one) of ID, Type, Number or Description fields.");
                    }

                    ID = idFields[0].Substring(3);
                    number = numberFields[0].Substring(7);
                    var typeString = typeFields[0].Substring(5);
                    switch (typeString.ToLower())
                    {
                        case "integer": type = Type.Integer; break;
                        case "float": type = Type.Float; break;
                        case "flag": type = Type.Flag; break;
                        case "characher": type = Type.Character; break;
                        case "string": type = Type.String; break;
                        default:
                            throw new Exception("ASETools.VCFFile.Info constructor: type field is of unknown value: " + typeString);
                    }

                    description = descriptionFields[0].Substring(13, descriptionFields[0].Length - 14); // Probably should check for internal quotes, escape sequences and properly handle quoted commas.  But not until it comes up.
                } // Info ctor
            } // Info class

            public class Filter
            {
                public readonly string ID;
                public readonly string description;

                public Filter(string filterLine)
                {
                    if (!filterLine.ToLower().StartsWith("##filter=<") || !filterLine.EndsWith("<"))
                    {
                        throw new Exception()
                    }
                }
            }

            Dictionary<string, Info> info = new Dictionary<string, Info>(); // Indexed by ID

            public readonly string fileFormat;
            public readonly string 
        }

        public class VCFLine
        {
            public static List<VCFLine> ReadFromFile(StreamReader inputFile, VCFFile file)
            {

            }
        }
#endif // false

        public static void WriteMatrixWithPercentages(StreamWriter outputFile, string title, int[,] matrix, bool normalize)
        {
            outputFile.WriteLine(title);

            var total = matrix.Cast<int>().Sum();  // https://stackoverflow.com/questions/19034970/sum-multidimensional-array-c-sharp

            for (int y = matrix.GetLength(1) - 1; y >= 0; y--)
            {
                outputFile.Write((double)y / (matrix.GetLength(1) - 1));
                for (int x = 0; x < matrix.GetLength(1); x++)
                {
                    if (normalize)
                    {
                        outputFile.Write("\t" + (double)matrix[x, y] / total);
                    }
                    else
                    {
                        outputFile.Write("\t" + matrix[x, y]);
                    }
                }
                outputFile.WriteLine();
            }

            //
            // Now write the column footers.
            //
            for (int x = 0; x < matrix.GetLength(0); x++)
            {
                outputFile.Write("\t" + (double)x / (matrix.GetLength(0) - 1));
            }
            outputFile.WriteLine();
        } // WriteMatrixWithPercentages

        public class FASTA
        {
            public FASTA(StreamReader inputFile)
            {
                string nextLine;

                List<string> currentContig = null;
                string currentContigName = null;

                while (null != (nextLine = inputFile.ReadLine()))
                {
                    if (nextLine.StartsWith(">"))
                    {
                        addContig(currentContigName, currentContig);
                        currentContig = new List<string>();
                        currentContigName = nextLine.Substring(1);  // Cut off the >
                        if (currentContigName.Contains(' '))
                        {
                            currentContigName = currentContigName.Substring(0, currentContigName.IndexOf(' '));
                        }
                    } else
                    {
                        if (currentContig == null)
                        {
                            Console.WriteLine("FASTA: input file does not start with a contig name");
                            throw new Exception("FASTA: input file does not start with a contig name");
                        }

                        currentContig.Add(nextLine.ToUpper());
                    }
                }

                addContig(currentContigName, currentContig);
            }

            void addContig(string contigName, List<string> contigStrings)
            {
                if (contigName == null)
                {
                    return;
                }

                if (contigs.ContainsKey(contigName))
                {
                    Console.WriteLine("FASTA: contig " + contigName + " appears at least twice");
                    throw new Exception("FASTA: contig " + contigName + " appears at least twice");
                }

                var contig = new char[contigStrings.Select(_ => _.Length).Sum()];
                int howFar = 0;
                foreach (var contigString in contigStrings)
                {
                    for (int i = 0; i < contigString.Length; i++)
                    {
                        contig[howFar + i] = contigString[i];
                    }
                    howFar += contigString.Length;
                }

                contigs.Add(contigName, contig);
            }

            public static FASTA loadFromFile(string inputFilename)
            {
                var inputFile = CreateStreamReaderWithRetry(inputFilename);
                if (null == inputFile)
                {
                    return null;
                }

                var retVal = new FASTA(inputFile);
                inputFile.Close();

                return retVal;
            }

            Dictionary<string, char[]> contigs = new Dictionary<string, char[]>();

            static public bool compare(FASTA one, FASTA two)
            {
                var fastas = new FASTA[2];
                fastas[0] = one;
                fastas[1] = two;

                bool equal = true;

                for (int i = 0; i < 2; i++) {
                    foreach (var contig in fastas[i].contigs)
                    {
                        if (!fastas[1 - i].contigs.ContainsKey(contig.Key))
                        {
                            Console.WriteLine("One FASTA has contig " + contig.Key + ", which isn't in the other");
                            equal = false;
                        }
                    }
                } // for each FASTA

                foreach (var contig in one.contigs)
                {
                    var contigName = contig.Key;
                    if (!two.contigs.ContainsKey(contigName))
                    {
                        continue;
                    }

                    if (contig.Value.Length != two.contigs[contigName].Length)
                    {
                        Console.WriteLine("Contig " + contigName + " differs in length: " + contig.Value.Length + " and " + two.contigs[contigName].Length);
                        equal = false;
                        continue;
                    }

                    for (int i = 0; i < contig.Value.Length; i++)
                    {
                        if (contig.Value[i] != two.contigs[contigName][i] && contig.Value[i] != 'N' && two.contigs[contigName][i] != 'N')
                        {
                            Console.WriteLine("Contig " + contigName + " first differs at " + i + " (zero based) with values: " + contig.Value[i] + " and " + two.contigs[contigName][i]);
                            equal = false;
                            break;
                        }
                    }
                }

                return equal;


            }
        } // FASTA

        public delegate V MergeDictionaryValues<V>(V one, V two);


        public static void MergeDictionaries<K, V>(Dictionary<K, V> into, Dictionary<K, V> from, MergeDictionaryValues<V> mergeFunction)
        {
            foreach (var k in from.Select(_ => _.Key))
            {
                if (!into.ContainsKey(k))
                {
                    into.Add(k, from[k]);
                }
                else
                {
                    into[k] = mergeFunction(into[k], from[k]);
                }
            } // for each key in from
        }

        public static void MergeDictionaries<K, V>(Dictionary<K, List<V>> into, Dictionary<K, List<V>> from)
        {
            foreach (var k in from.Select(_ => _.Key))
            {
                if (!into.ContainsKey(k))
                {
                    into.Add(k, from[k]);
                }
                else
                {
                    into[k].AddRange(from[k]);
                }
            }
        } // MergeDictionaries

        public static void MergeDictionaries<K>(Dictionary<K, int> into, Dictionary<K, int> from)
        {
            MergeDictionaries(into, from, (a, b) => a + b);
        } // MergeDictionaries

        public static void MergeDictionaries<K>(Dictionary<K, double> into, Dictionary<K, double> from)
        {
            MergeDictionaries(into, from, (a, b) => a + b);
        } // MergeDictionaries

        public static void MergeDictionaries<K>(Dictionary<K, PreBucketedHistogram> into, Dictionary<K, PreBucketedHistogram> from)
        {
            foreach (var k in from.Select(_ => _.Key))
            {
                if (!into.ContainsKey(k))
                {
                    into.Add(k, from[k]);
                } else
                {
                    into[k].merge(from[k]);
                }
            }
        } // MergeDictionaries

        public class VariantPairPhasing
        {
            public readonly string contig;
            public readonly int locus0, locus1;
            public readonly int[] loci = new int[2];
            public readonly bool somatic0, somatic1;
            public readonly bool[] somatic = new bool[2];
            public readonly string ref0, alt0, ref1, alt1;
            public readonly string[] ref_ = new string[2], alt = new string[2];
            public readonly Dictionary<bool, Dictionary<bool, Dictionary<bool, int>>> counts = new Dictionary<bool, Dictionary<bool, Dictionary<bool, int>>>();    // Doesn't include weird.  Maps dna->ref locus 0->ref locus 1->count
            public readonly Dictionary<bool, int> weirdCounts = new Dictionary<bool, int>();  // DNA->weird count
            public readonly bool dnaConsistent, rnaConsistent;
            public readonly Dictionary<bool, bool> consistent = new Dictionary<bool, bool>();  // DNA->consistent
            public readonly bool mutuallyConsistent;

            public static List<VariantPairPhasing> readFromFile(string filename)
            {
                var inputFile = CreateStreamReaderWithRetry(filename);
                if (null == inputFile)
                {
                    return null;
                }

                var headerizedFile = new HeaderizedFile<VariantPairPhasing>(inputFile, false, true, "", wantedFields.ToList());
                List<VariantPairPhasing> retVal;
                if (!headerizedFile.ParseFile(ParseOneLine, out retVal))
                {
                    return null;
                }

                return retVal;
            }

            public int nSomaticVariants()
            {
                int retVal = 0;
                if (somatic0) retVal++;
                if (somatic1) retVal++;

                return retVal;
            }

            public int countOfReadsForOneReadType(bool dna)
            {
                return counts[dna][true][true] + counts[dna][true][false] + counts[dna][false][true] + counts[dna][false][false];
            }

            public bool hasEnoughReadsForOneReadType(bool dna, int minReads)
            {
                return countOfReadsForOneReadType(dna) >= minReads;
            }

            public bool hasEnoughReadsOfEitherReadType(int minReads)
            {
                return ASETools.BothBools.Any(dna => hasEnoughReadsForOneReadType(dna, minReads));
            }



            VariantPairPhasing(string contig_, int locus0_, int locus1_, bool somatic0_, bool somatic1_, string ref0_, string ref1_, string alt0_, string alt1_, int dnaRR, int dnaAA, int dnaRA, int dnaAR, int dnaWeird,
                               int rnaRR, int rnaAA, int rnaRA, int rnaAR, int rnaWeird, bool dnaConsistent_, bool rnaConsistent_, bool mutuallyConsistent_)
            {
                contig = contig_;
                locus0 = locus0_;
                locus1 = locus1_;
                somatic0 = somatic0_;
                somatic1 = somatic1_;
                ref0 = ref0_;
                ref1 = ref1_;
                alt0 = alt0_;
                alt1 = alt1_;
                dnaConsistent = dnaConsistent_;
                rnaConsistent = rnaConsistent_;
                mutuallyConsistent = mutuallyConsistent_;

                loci[0] = locus0;
                loci[1] = locus1;

                somatic[0] = somatic0;
                somatic[1] = somatic1;

                ref_[0] = ref0;
                ref_[1] = ref1;

                alt[0] = alt0;
                alt[1] = alt1;

                foreach (var dna in ASETools.BothBools)
                {
                    counts.Add(dna, new Dictionary<bool, Dictionary<bool, int>>());

                    foreach (var refLocus0 in ASETools.BothBools)
                    {
                        counts[dna].Add(refLocus0, new Dictionary<bool, int>());

                        // Don't bother with the last Add, we do it below when we populate.
                    } // refLocus0
                } // dna

                counts[true][true].Add(true, dnaRR);
                counts[true][true].Add(false, dnaRA);
                counts[true][false].Add(true, dnaAR);
                counts[true][false].Add(false, dnaAA);
                counts[false][true].Add(true, rnaRR);
                counts[false][true].Add(false, rnaRA);
                counts[false][false].Add(true, rnaAR);
                counts[false][false].Add(false, rnaAA);

                weirdCounts.Add(true, dnaWeird);
                weirdCounts.Add(false, rnaWeird);

                consistent.Add(true, dnaConsistent);
                consistent.Add(false, rnaConsistent);

            }

            static VariantPairPhasing ParseOneLine(HeaderizedFile<VariantPairPhasing>.FieldGrabber fieldGrabber)
            {
                return new VariantPairPhasing(
                                contig_: fieldGrabber.AsString("Contig"),
                                locus0_: fieldGrabber.AsInt("locus 0"),
                                locus1_: fieldGrabber.AsInt("locus 1"),
                                somatic0_: fieldGrabber.AsBool("Somatic 0"),
                                somatic1_: fieldGrabber.AsBool("Somatic 1"),
                                ref0_: fieldGrabber.AsString("Ref 0"),
                                ref1_: fieldGrabber.AsString("Ref 1"),
                                alt0_: fieldGrabber.AsString("Alt 0"),
                                alt1_: fieldGrabber.AsString("Alt 1"),
                                dnaRR: fieldGrabber.AsInt("DNA ref ref"),
                                dnaRA: fieldGrabber.AsInt("DNA ref alt"),
                                dnaAR: fieldGrabber.AsInt("DNA alt ref"),
                                dnaAA: fieldGrabber.AsInt("DNA alt alt"),
                                rnaRR: fieldGrabber.AsInt("RNA ref ref"),
                                rnaRA: fieldGrabber.AsInt("RNA ref alt"),
                                rnaAR: fieldGrabber.AsInt("RNA alt ref"),
                                rnaAA: fieldGrabber.AsInt("RNA alt alt"),
                                dnaWeird: fieldGrabber.AsInt("DNA weird"),
                                rnaWeird: fieldGrabber.AsInt("RNA weird"),
                                dnaConsistent_: fieldGrabber.AsBool("DNA consistent"),
                                rnaConsistent_: fieldGrabber.AsBool("RNA consistent"),
                                mutuallyConsistent_: fieldGrabber.AsBool("DNA consistent with RNA")
                                );

            } // ParseOneLine

            static readonly string[] wantedFields =
            {
                "Contig",
                "locus 0",
                "locus 1",
                "Somatic 0",
                "Somatic 1",
                "Ref 0",
                "Alt 0",
                "Ref 1",
                "Alt 1",
                "DNA ref ref",
                "DNA alt alt",
                "DNA ref alt",
                "DNA alt ref",
                "DNA weird",
                "RNA ref ref",
                "RNA alt alt",
                "RNA ref alt",
                "RNA alt ref",
                "RNA weird",
                "DNA consistent",
                "RNA consistent",
                "DNA consistent with RNA"
            };
        } // VariantPairPhasing

        public class VCFStatistics
        {
            public static VCFStatistics LoadFromFile(string filename)
            {
                var inputFile = CreateStreamReaderWithRetry(filename);
                if (null == inputFile)
                    return null;

                var headerizedFile = new HeaderizedFile<KeyValuePair<string, int>>(inputFile, false, true, "", wantedFields.ToList());
                List<KeyValuePair<string, int>> listOfChromosomeStats;
                if (!headerizedFile.ParseFile(parser, out listOfChromosomeStats))
                {
                    return null;
                }

                var vcfCounts = new Dictionary<string, int>();
                listOfChromosomeStats.ForEach(_ => vcfCounts.Add(_.Key, _.Value));
                return new VCFStatistics(vcfCounts);

            } // LoadFromFile

            VCFStatistics(Dictionary<string, int> vcfCounts_)
            {
                vcfCounts = vcfCounts_;
            }

            static KeyValuePair<string, int> parser(ASETools.HeaderizedFile<KeyValuePair<string, int>>.FieldGrabber fieldGrabber)
            {
                return new KeyValuePair<string, int>(fieldGrabber.AsString("Chromosome"), fieldGrabber.AsInt("Count of variants"));
            }

            public int totalVariants()
            {
                return vcfCounts.Select(_ => _.Value).Sum();
            }

            static readonly string[] wantedFields =
            {
                "Chromosome",
                "Count of variants"
            }; // wantedFields

            public readonly Dictionary<string, int> vcfCounts;
        } // VCFStatistics

        const int MaxNMToTrack = 100;
        public class ReadStatistics
        {
            public PreBucketedHistogram distributionOfMappedReadsByNM = new PreBucketedHistogram(0, MaxNMToTrack, 1);
            public PreBucketedHistogram distributionOfMAPQ = new PreBucketedHistogram(0, 61, 1);
            public PreBucketedHistogram distributionOfPairedReadGap = new PreBucketedHistogram(0, 2000, 5);
            public PreBucketedHistogram distributionOfReadLengths = new PreBucketedHistogram(0, 200, 1);
            public long unmappedReads = 0;
            public long totalReads = 0;
            public long pairedReads = 0;
            public long bothHalvesOfPairMapped = 0;
            public long forwardReads = 0;
            public long RCReads = 0;
            public long crossChromosomeMappedPairs = 0;
            public long secondaryAlignments = 0;
            public long supplementaryAlignments = 0;

            public Dictionary<string, int> countOfMappedReadsByChromosome = new Dictionary<string, int>();
            public long readsInALTContigs = 0;
            public long readsInNonNormalNonALTContigs = 0; // unplaced parts of the primary assembly

            public ReadStatistics()
            {
                foreach (var chromosome in chromosomesWithMitochondria)
                {
                    countOfMappedReadsByChromosome.Add(chromosome, 0);
                }
            } // ctor

            public void addSAMLine(string samText, string filename)
            {
                var samLine = SAMLine.attemptSAMLine(samText);
                if (samLine == null)
                {
                    Console.WriteLine("Ignoring bad SAMLine in " + filename);
                    return;
                }
                addSAMLine(samLine);
            }

            public void addSAMLine(SAMLine samLine)
            {

                totalReads++;

                if (samLine.isPaired())
                {
                    pairedReads++;  // We count unmapped paired reads as paired

                    if (samLine.pnext != 0 && samLine.bothHalvesOfPairMapped())
                    {
                        bothHalvesOfPairMapped++;

                        if (samLine.rnext != "=" && samLine.rnext != "*" && samLine.rnext != samLine.rname)
                        {
                            crossChromosomeMappedPairs++;
                        } else
                        {
                            if (0 != samLine.tlen)
                            {
                                distributionOfPairedReadGap.addValue(Math.Abs(samLine.tlen));
                            }
                        }
                    } // both halves of pair mapped
                } // paired

                distributionOfReadLengths.addValue(samLine.seq.Length);

                if (samLine.isSecondaryAlignment())
                {
                    secondaryAlignments++;
                }

                if (samLine.isSupplementaryAlignment())
                {
                    supplementaryAlignments++;
                }

                if (samLine.isUnmapped())
                {
                    unmappedReads++;
                    return;
                }

                if (samLine.isRC())
                {
                    RCReads++;
                } else
                {
                    forwardReads++;
                }

                distributionOfMAPQ.addValue(samLine.mapq);

                var contig = samLine.rname.ToLower();

                if (chromosomes.Contains(contig))
                {
                    countOfMappedReadsByChromosome[contig]++;
                } else if (contig.EndsWith("_alt"))
                {
                    readsInALTContigs++;
                } else
                {
                    readsInNonNormalNonALTContigs++;
                }

                if (samLine.NMKnown())
                {
                    distributionOfMappedReadsByNM.addValue(samLine.NM());
                }

            } // addSAMLine

            public long mappedReads() { return totalReads - unmappedReads; }

            public void WriteToFile(StreamWriter outputFile)
            {
                outputFile.WriteLine("Total Reads\tUnmapped Reads\tPaired Reads\tReads in ALT contigs\tReads in Non-Normal, Non-ALT contigs\tBoth Halves of Pair Mapped\tForward\tRC\tCross Chromosome Pairs\tSecondary Alignments\tSupplementary Alignments");
                outputFile.WriteLine(totalReads + "\t" + unmappedReads + "\t" + pairedReads + "\t" + readsInALTContigs + "\t" + readsInNonNormalNonALTContigs + "\t" + bothHalvesOfPairMapped + "\t" + forwardReads
                                     + "\t" + RCReads + "\t" + crossChromosomeMappedPairs + " \t" + secondaryAlignments + "\t" + supplementaryAlignments);

                outputFile.WriteLine();
                outputFile.WriteLine("Chromosome\tMapped Reads");
                countOfMappedReadsByChromosome.ToList().ForEach(_ => outputFile.WriteLine(_.Key + "\t" + _.Value));

                outputFile.WriteLine();
                outputFile.WriteLine("Distribution of NM");
                distributionOfMappedReadsByNM.WriteHistogram(outputFile);

                outputFile.WriteLine();
                outputFile.WriteLine("Distribution of MAPQ");
                distributionOfMAPQ.WriteHistogram(outputFile);

                outputFile.WriteLine();
                outputFile.WriteLine("Distribution of Paired Read Gap");
                distributionOfPairedReadGap.WriteHistogram(outputFile);

                outputFile.WriteLine();
                outputFile.WriteLine("Distribution of Read Lengths");
                distributionOfReadLengths.WriteHistogram(outputFile);
            }


            public static ReadStatistics ReadFromFile(string inputFilename)
            {
                var inputFile = CreateStreamReaderWithRetry(inputFilename);
                if (null == inputFile)
                {
                    Console.WriteLine("ASELib.ReadStatictics.ReadFromFile: unable to open input file " + inputFilename);
                    return null;
                }

                var result = ReadFromFile(inputFile);
                inputFile.Close();

                return result;
            }


            public static ReadStatistics ReadFromFile(StreamReader inputFile)
            {
                //
                // The first two lines are a header and the counts (possibly preceeded by blank lines), followed by a blank line.
                //
                var headerLines = new List<string>();

                string inputLine;
                while ((inputLine = inputFile.ReadLine()) == "") { }

                headerLines.Add(inputFile.ReadLine());
                headerLines.Add(inputFile.ReadLine());

                inputFile.ReadLine();

                var headerizedFile = new HeaderizedFile<ReadStatistics>(StreamReaderFromStrings(headerLines), false, false, "", wantedFields.ToList());
                List<ReadStatistics> result;
                if (!headerizedFile.ParseFile(Parser, out result))
                {
                    return null;
                }

                if (result.Count() != 1)
                {
                    throw new Exception("ASELib.ReadStatistics.ReadFromFile: got more than one result from HeaderizedLine (??)");
                }

                inputLine = inputFile.ReadLine();
                if (inputLine != "Chromosome\tMapped Reads")
                {
                    throw new Exception("ASETools.ReadStatictics.ReadFromFile: didn't find Chromosome Mapped Reads line where expected. Got: " + inputLine + " instead.");
                }

                while ((inputLine = inputFile.ReadLine()) != "")
                {
                    var fields = inputLine.Split('\t');
                    if (fields.Count() != 2)
                    {
                        throw new Exception("ASETools.ReadStatictics.ReadFromFile: Chromosme Mapped Reads line didnt' have two fields: " + inputLine);
                    }

                    result[0].countOfMappedReadsByChromosome[fields[0]] = Convert.ToInt32(fields[1]);
                }

                inputLine = inputFile.ReadLine();
                if (inputLine != "Distribution of NM")
                {
                    throw new Exception("ASETools.ReadStatictics.ReadFromFile: didn't find Distribution of NM line where expected. Got: " + inputLine + " instead.");
                }
                result[0].distributionOfMappedReadsByNM = PreBucketedHistogram.ReadFromSteam(inputFile, "", false);

                return result[0];
            } // ReadFromFile

            static ReadStatistics Parser(HeaderizedFile<ReadStatistics>.FieldGrabber fieldGrabber)
            {
                var readStatistics = new ReadStatistics();
                readStatistics.totalReads = fieldGrabber.AsLong("Total Reads");
                readStatistics.unmappedReads = fieldGrabber.AsLong("Unmapped Reads");
                readStatistics.pairedReads = fieldGrabber.AsLong("Paired Reads");
                readStatistics.readsInALTContigs = fieldGrabber.AsLong("Reads in ALT contigs");
                readStatistics.readsInNonNormalNonALTContigs = fieldGrabber.AsLong("Reads in Non-Normal, Non-ALT contigs");
                readStatistics.bothHalvesOfPairMapped = fieldGrabber.AsLong("Both Halves of Pair Mapped");
                readStatistics.forwardReads = fieldGrabber.AsLong("Forward");
                readStatistics.RCReads = fieldGrabber.AsLong("RC");
                readStatistics.crossChromosomeMappedPairs = fieldGrabber.AsLong("Cross Chromosome Pairs");
                readStatistics.secondaryAlignments = fieldGrabber.AsLong("Secondary Alignments");
                readStatistics.supplementaryAlignments = fieldGrabber.AsLong("Supplementary Alignments");

                return readStatistics;
            }

            static string[] wantedFields =
            {
                "Total Reads", "Unmapped Reads", "Paired Reads", "Reads in ALT contigs", "Reads in Non-Normal, Non-ALT contigs", "Both Halves of Pair Mapped", "Forward", "RC", "Cross Chromosome Pairs", "Secondary Alignments", "Supplementary Alignments"
            };
        } // ReadStatistics

        public static StreamReader StreamReaderFromStrings(List<string> strings)
        {
            string concatenatedString = "";

            strings.ForEach(_ => concatenatedString += _ + "\n");

            return new StreamReader(new MemoryStream(System.Text.Encoding.UTF8.GetBytes(concatenatedString)));
        }

        public class miRNAExpressionQuantification
        {
            public readonly string miRNA_ID;
            public readonly int read_count;
            public readonly double reads_mer_milion_miRNA_mapped;
            public readonly bool cross_mapped;

            public static Dictionary<string, miRNAExpressionQuantification> LoadFromFile(string inputFilename) // maps name->miRNAExpressionQuantification
            {
                var inputFile = CreateStreamReaderWithRetry(inputFilename);
                if (inputFile == null)
                {
                    return null;
                }

                string[] wantedFields = { "miRNA_ID", "read_count", "reads_per_million_miRNA_mapped", "cross-mapped" };

                var headerizedFile = new HeaderizedFile<miRNAExpressionQuantification>(inputFile, false, false, "", wantedFields.ToList());
                List<miRNAExpressionQuantification> listOfmiRNAs;


                Dictionary<string, miRNAExpressionQuantification> retVal;

                if (!headerizedFile.ParseFile(parser, out listOfmiRNAs))
                {
                    retVal = null;
                } else
                {
                    retVal = listOfmiRNAs.GroupByToDictUnique(_ => _.miRNA_ID);
                }

                inputFile.Close();
                return retVal;
            }

            static miRNAExpressionQuantification parser(HeaderizedFile<miRNAExpressionQuantification>.FieldGrabber fieldGrabber)
            {
                bool crossMapped;
                if (fieldGrabber.AsString("cross-mapped") == "Y")
                {
                    crossMapped = true;
                }
                else if (fieldGrabber.AsString("cross-mapped") == "N")
                {
                    crossMapped = false;
                } else
                {
                    throw new Exception("ASETools.miRNAExpressionQuantification.parser: cross-mapped column contains neither Y nor N for line " + fieldGrabber.rawLine());
                }

                return new miRNAExpressionQuantification(fieldGrabber.AsString("miRNA_ID"), fieldGrabber.AsInt("read_count"), fieldGrabber.AsDouble("reads_per_million_miRNA_mapped"), crossMapped);
            }

            miRNAExpressionQuantification(string miRNA_ID_, int read_count_, double reads_mer_milion_miRNA_mapped_, bool cross_mapped_)
            {
                miRNA_ID = miRNA_ID_;
                read_count = read_count_;
                reads_mer_milion_miRNA_mapped = reads_mer_milion_miRNA_mapped_;
                cross_mapped = cross_mapped_;
            }
        } // miRNAExpressionQuantification

        //
        // This is like Convert.ToInt32, except it removes any commas and if the string extends beyond the number it will truncate it.  Strings that don't start
        // with a number will still throw a format exception.
        //
        public static int GetIntFromString(string value)
        {
            var withoutCommas = value.Replace(",", "");

            for (int i = 0; i < withoutCommas.Length; i++)
            {
                if (withoutCommas[i] < '0' || withoutCommas[i] > '9')
                {
                    withoutCommas = withoutCommas.Substring(0, i);
                    break;
                }
            }

            return Convert.ToInt32(withoutCommas);
        } // GetIntFromString

        //
        // This is like Convert.ToInt32, except it allows stuff to be after the value in the string.
        //
        public static double GetDoubleFromString(string value)
        {
            for (int i = 0; i < value.Length; i++)
            {
                if ((value[i] < '0' || value[i] > '9') && value[i] != '-' && value[i] != 'e' && value[i] != '+' && value[i] != '.')
                {
                    value = value.Substring(0, i);
                    break;
                }
            }

            return Convert.ToDouble(value);
        } // GetDoubleFromString

        public class BiasedRandom<KeyType>
        {
            public BiasedRandom(List<KeyValuePair<KeyType, ulong>> biasValues)
            {
                biasValues.ForEach(_ => sortedBiasValues.Add(_));

                sortedBiasValues.Sort(comparer);
                biasTotal = 0;
                for (int i = 0; i < sortedBiasValues.Count(); i++)  // What genius decided that ulong shoudln't have .Sum() implemented??
                {
                    biasTotal += sortedBiasValues[i].Value;
                }
            }

            public KeyType select()
            {
                //
                // Code for long-valued rand in C# cribbed from https://gist.github.com/subena22jf/c7bb027ea99127944981 and modified slightly
                //
 
                ulong ulongRand;
                do
                {
                    byte[] buf = new byte[8];
                    rand.NextBytes(buf);
                    ulongRand = (ulong)BitConverter.ToInt64(buf, 0);
                } while (ulongRand > ulong.MaxValue - ((ulong.MaxValue % biasTotal) + 1) % biasTotal);

                ulong runningTotal = 0;
                for (int i = 0; i < sortedBiasValues.Count(); i++)
                {
                    runningTotal += sortedBiasValues[i].Value;  // It's not at all clear that we need to sort for this to work...

                    if (ulongRand <= runningTotal)
                    {
                        return sortedBiasValues[i].Key;
                    }
                }

                throw new Exception("BaisedRandom.select(): shouldn't get here.  biasTotal = " + biasTotal + " ulongRand = " + ulongRand + " runningTotal = " + runningTotal);
            }


            ulong biasTotal;
            List<KeyValuePair<KeyType, ulong>> sortedBiasValues = new List<KeyValuePair<KeyType, ulong>>();
            Random rand = new Random();

            static int comparer(KeyValuePair<KeyType, ulong> first, KeyValuePair<KeyType, ulong> second)
            {
                return first.Value.CompareTo(second.Value);
            }
        }

    } // ASETools

    //
    // I love that C# lets you do stuff like this.  This is a simple groupBy method that takes an IEnumerable and a key extractor function, and returns a dictionary indexed on keys of lists 
    // of elements in the enumerable.  This is a slight tweak on the existing GroupBy methods that fits my needs a little better.
    //
    public static class GroupByExtension
    {
        public static Dictionary<TKey, List<TSource>> GroupByToDict<TSource, TKey>(this IEnumerable<TSource> source, Func<TSource, TKey> keyExtractor)
        {
            var retVal = new Dictionary<TKey, List<TSource>>();

            foreach (var element in source)
            {
                var key = keyExtractor(element);

                if (!retVal.ContainsKey(key))
                {
                    retVal.Add(key, new List<TSource>());
                }

                retVal[key].Add(element);
            }

            return retVal;
        } // GroupByToDict

        public static Dictionary<TKey, TSource> GroupByToDictUnique<TSource, TKey>(this IEnumerable<TSource> source, Func<TSource, TKey> keyExtractor)
        {
            var retVal = new Dictionary<TKey, TSource>();

            foreach (var element in source)
            {
                retVal.Add(keyExtractor(element), element); // If it's not unique, this throws an exception.
            }

            return retVal;
        }
    } // GroupByExtension

    public static class ListExtension
    {
        //
        // Make a string like elt0, elt1, elt2 
        //
        // If you want it to have a length other than its actual length, then set length.  If length is longer than
        // source.Count(), then it's padded with padding.
        //
        public static string EnumerateWithCommas<TValue>(this List<TValue> source, int length = -1, string padding = "")  // Probably could do this for IEnumerable...
        {
            var count = source.Count();
            var lengthToUse = (length < 0) ? count : length;
            var retVal = "";
            for (int i = 0; i < lengthToUse; i++)
            {
                if (i >= count)
                {
                    retVal += padding;
                }
                else
                {
                    retVal += source[i];
                }

                if (i != lengthToUse - 1)
                {
                    retVal += ", ";
                }
            }

            return retVal;
        } // EnumerateWithCommas
    } // ListExtension

    public static class IEnumerableExtension
    {
        // The code for OrderRandomly comes from https://stackoverflow.com/questions/254573/can-you-enumerate-a-collection-in-c-sharp-out-of-order, but I modified it to remove the RemoveAt() call.
        public static IEnumerable<T> OrderRandomly<T>(this IEnumerable<T> sequence)
        {
            Random random = new Random();
            var copy = sequence.ToArray();

            for (int i = copy.Count(); i > 0; i--)
            {
                int index = random.Next(i);
                yield return copy[index];
                copy[index] = copy[i-1];  // This removes the returned item from the active part of the array.  There's no need to swap, however.
            }
        } // OrderRandomly
    } // IEnumerableExtension

    public static class EnumerableMedian    // There's certainly a way to do this for generic numeric types rather than just double, but I haven't bothered to figure out how.
    {
        public static double Median(this IEnumerable<double> input) 
        {
            var list = input.ToList();
            var n = list.Count();
            if (n == 0)
            {
                return 0;
            }

            list.Sort();

            if (n % 2 == 1)
            {
                return list[n / 2];
            } else
            {
                return (list[n / 2] + list[n / 2 - 1]) / 2;
            }
        } // Median

        public static double Median(this IEnumerable<int> input)
        {
            var list = input.ToList();
            var n = list.Count();
            if (n == 0)
            {
                return 0;
            }

            list.Sort();

            if (n % 2 == 1)
            {
                return list[n / 2];
            }
            else
            {
                return ((double)list[n / 2] + list[n / 2 - 1]) / 2;
            }
        } // Median 
    } // EnumerableMedian
}
