using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using System.IO;
using System.Threading;
using System.IO.Compression;
using System.Net;
using System.Web;
using System.Runtime.Serialization;
using System.Runtime.Serialization.Json;
using System.Diagnostics;
using System.Runtime.InteropServices;
using MathNet.Numerics;

namespace ASELib
{
    public class ASETools
    {
        public const string urlPrefix = @"https://gdc-api.nci.nih.gov/";

        public const int GuidStringLength = 36;

		public const int nHumanNuclearChromosomes = 24;   // 1-22, X and Y.

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
			return !isChromosomeSex(chromosome) && !isChromosomeMitochondrial(chromosome);
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
		}

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

            GeneScatterGraphLine(string Hugo_Symbol_, string Chromosome_, int Start_Position_, string Variant_Classification_, string Variant_Type_, string Reference_Allele_, string Alt_Allele_, string disease_,
                string case_id_, string normal_dna_file_id_, string tumor_dna_file_id_, string normal_rna_file_id_, string tumor_rna_file_id_, ReadCounts normalDNAReadCounts_, ReadCounts tumorDNAReadCounts_,
                ReadCounts normalRNAReadCounts_, ReadCounts tumorRNAReadCounts_, bool MultipleMutationsInThisGene_, int nMutationsThisGene_, double tumorDNAFraction_, double tumorRNAFraction_, double tumorDNAMultiple_, double tumorRNAMultiple_, double tumorDNARatio_,
                double tumorRNARatio_, double ratioOfRatios_,  bool zKnown_, double zTumor_, double zNormal_, double z2Tumor_, double z2Normal_, double percentMeanTumor_,
                double percentMeanNormal_, bool fromUnfilteredFile_)
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
            }



			public static List<GeneScatterGraphLine> LoadAllGeneScatterGraphEntries(string directoryName, bool fromUnfiltered, string hugoSymbol /* this may be * to load all*/)
			{
				var geneScatterGraphEntries = new List<GeneScatterGraphLine>();



                string[] wantedFieldsArrayUnfilteredVersion =
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

                string[] wantedFieldsArrayAdditionalFields =
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

                var wantedFieldsFullVersion = wantedFieldsArrayUnfilteredVersion.ToList();
                wantedFieldsFullVersion.AddRange(wantedFieldsArrayAdditionalFields.ToList());


                var wantedFields = fromUnfiltered ? wantedFieldsArrayUnfilteredVersion.ToList() : wantedFieldsFullVersion;

				foreach (var filename in Directory.EnumerateFiles(directoryName, hugoSymbol + (fromUnfiltered ? Configuration.unfilteredCountsExtention : ".txt")))
				{
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
                    var headerizedFile = new HeaderizedFile<GeneScatterGraphLine>(inputFile, false, true, "", wantedFields);

                    List<GeneScatterGraphLine> linesFromThisFile;

                    headerizedFile.ParseFile(ParseLine, out linesFromThisFile);

                    geneScatterGraphEntries.AddRange(linesFromThisFile);
				}

				return geneScatterGraphEntries;
			}

            static GeneScatterGraphLine ParseLine(Dictionary<string, int> fieldMappings, string[] fields)
            {
                ReadCounts normalRNAReadCounts;

                if (fields[fieldMappings["n_normal_RNA_Matching_Reference"]] != "")
                {
                    normalRNAReadCounts = new ReadCounts(Convert.ToInt32(fields[fieldMappings["n_normal_RNA_Matching_Reference"]]), Convert.ToInt32(fields[fieldMappings["n_normal_RNA_Matching_Alt"]]), Convert.ToInt32(fields[fieldMappings["n_normal_RNA_Matching_Neither"]]),
                        Convert.ToInt32(fields[fieldMappings["n_normal_RNA_Matching_Both"]]));
                }
                else
                {
                    normalRNAReadCounts = null;
                }

                double tumorDNAFraction = -1, tumorRNAFraction = -1, tumorDNAMultiple = -1, tumorRNAMultiple = -1, tumorDNARatio = -1, tumorRNARatio = -1, ratioOfRatios = -1, zTumor = -1, zNormal = -1, z2Tumor = -1, z2Normal = -1, percentMeanTumor = -1, percentMeanNormal = -1;
                bool fromUnfilteredFile = !fieldMappings.ContainsKey("tumorDNAFraction");
                bool zKnown = false;
                if (!fromUnfilteredFile) {
                    tumorDNAFraction = Convert.ToDouble(fields[fieldMappings["tumorDNAFraction"]]);
                    tumorRNAFraction = Convert.ToDouble(fields[fieldMappings["tumorRNAFraction"]]);
                    tumorDNAMultiple = Convert.ToDouble(fields[fieldMappings["tumorDNAMultiple"]]);
                    tumorRNAMultiple = Convert.ToDouble(fields[fieldMappings["tumorRNAMultiple"]]);
                    tumorDNARatio = Convert.ToDouble(fields[fieldMappings["tumorDNARatio"]]);
                    tumorRNARatio = Convert.ToDouble(fields[fieldMappings["tumorRNARatio"]]);
                    ratioOfRatios = Convert.ToDouble(fields[fieldMappings["RatioOfRatios"]]);
                    if (fields[fieldMappings["zTumor"]] != "")
                    {
                        zKnown = true;
                        zTumor = Convert.ToDouble(fields[fieldMappings["zTumor"]]);
                        zNormal = Convert.ToDouble(fields[fieldMappings["zNormal"]]);
                        z2Tumor = Convert.ToDouble(fields[fieldMappings["z2Tumor"]]);
                        z2Normal = Convert.ToDouble(fields[fieldMappings["z2Normal"]]);
                        percentMeanTumor = Convert.ToDouble(fields[fieldMappings["%MeanTumor"]]);
                        percentMeanNormal = Convert.ToDouble(fields[fieldMappings["%MeanNormal"]]);
                    }
                }

                return new GeneScatterGraphLine(
                  ConvertToNonExcelString(fields[fieldMappings["Hugo_Symbol"]]), fields[fieldMappings["Chromosome"]], Convert.ToInt32(fields[fieldMappings["Start_Position"]]), fields[fieldMappings["Variant_Classification"]], fields[fieldMappings["Variant_Type"]],
                    fields[fieldMappings["Reference_Allele"]], fields[fieldMappings["Reference_Allele"]], fields[fieldMappings["disease"]], fields[fieldMappings["Case Id"]], fields[fieldMappings["Normal DNA File ID"]], fields[fieldMappings["Tumor DNA File ID"]],
                    fields[fieldMappings["Normal RNA File ID"]], fields[fieldMappings["Tumor RNA File ID"]],
                    new ReadCounts(Convert.ToInt32(fields[fieldMappings["n_normal_DNA_Matching_Reference"]]), Convert.ToInt32(fields[fieldMappings["n_normal_DNA_Matching_Alt"]]), Convert.ToInt32(fields[fieldMappings["n_normal_DNA_Matching_Neither"]]),
                        Convert.ToInt32(fields[fieldMappings["n_normal_DNA_Matching_Both"]])),
                    new ReadCounts(Convert.ToInt32(fields[fieldMappings["n_tumor_DNA_Matching_Reference"]]), Convert.ToInt32(fields[fieldMappings["n_tumor_DNA_Matching_Alt"]]), Convert.ToInt32(fields[fieldMappings["n_tumor_DNA_Matching_Neither"]]),
                        Convert.ToInt32(fields[fieldMappings["n_tumor_DNA_Matching_Both"]])),
                    normalRNAReadCounts,
                    new ReadCounts(Convert.ToInt32(fields[fieldMappings["n_tumor_RNA_Matching_Reference"]]), Convert.ToInt32(fields[fieldMappings["n_tumor_RNA_Matching_Alt"]]), Convert.ToInt32(fields[fieldMappings["n_tumor_RNA_Matching_Neither"]]),
                        Convert.ToInt32(fields[fieldMappings["n_tumor_RNA_Matching_Both"]])),
                    Convert.ToBoolean(fields[fieldMappings["Multiple Mutations in this Gene"]]), Convert.ToInt32(fields[fieldMappings["n Mutations in this gene"]]), tumorDNAFraction, tumorRNAFraction, tumorDNAMultiple, tumorRNAMultiple, tumorDNARatio, tumorRNARatio,
                    ratioOfRatios, zKnown, zTumor, zNormal, z2Tumor, z2Normal, percentMeanTumor, percentMeanNormal, fromUnfilteredFile);
            }


        } // GeneScatterGraphLine


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


		public class Exon
		{
			public Exon(string startString, string endString)
			{
				start = Convert.ToInt32(startString);
				end = Convert.ToInt32(endString);
			}

			public int start;
			public int end;
		}

		public class Isoform    // These come from the knownGene file
		{
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
		}


		// Holds Gene location and isoform information
		public class GeneLocationInfo
		{
			public string hugoSymbol;
			public string chromosome;   // in non-chr form
			public int minLocus;
			public int maxLocus;

			public bool inconsistent = false;

            public bool containsLocus(string locusChromosome, int locus)
            {
                return chromosome == locusChromosome && minLocus <= locus && maxLocus >= locus;
            }

			public List<Isoform> isoforms = new List<Isoform>();
		}

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
		}

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
                    // We need to split any regions that cross the beginnig and end of the gene, and add the gene to any regions entirely
                    // contained within it.
                    //

                    int handledUpTo = gene.minLocus;

                    Range overlapRange;
                    if (map.FindFirstLessThanOrEqualTo(key, out overlapRange) &&  overlapRange.chromosome == gene.chromosome && overlapRange.minLocus < gene.minLocus && overlapRange.maxLocus >= gene.minLocus)
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

            }

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
			public string project_id;   // This is TCGA_<DiseaseType> for TCGA.
            public List<string> sample_ids = new List<string>();

            //
            // Pathnames for downloaded files.
            //
            public string normal_dna_filename = "";
            public string tumor_dna_filename = "";
            public string normal_rna_filename = "";
            public string tumor_rna_filename = "";
	        public string maf_filename = "";
            public string tumor_methylation_filename = "";
			public string normal_methylation_filename = "";
			public string tumor_copy_number_filename = "";
			public string normal_copy_number_filename = "";

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


			//
			// Pathnames for derived files.
			//
			public string normal_dna_allcount_filename = "";
            public string tumor_dna_allcount_filename = "";
            public string normal_rna_allcount_filename = "";
            public string tumor_rna_allcount_filename = "";
            public string regional_expression_filename = "";
            public string gene_expression_filename = "";
            public string selected_variants_filename = "";
            public string tumor_dna_reads_at_selected_variants_filename = "";
            public string tumor_dna_reads_at_selected_variants_index_filename = "";
            public string tumor_rna_reads_at_selected_variants_filename = "";
            public string tumor_rna_reads_at_selected_variants_index_filename = "";
            public string normal_rna_reads_at_selected_variants_filename = "";
            public string normal_rna_reads_at_selected_variants_index_filename = "";
            public string normal_dna_reads_at_selected_variants_filename = "";
            public string normal_dna_reads_at_selected_variants_index_filename = "";
            public string annotated_selected_variants_filename = "";
            public string normal_allele_specific_gene_expression_filename = "";
			public string tumor_allele_specific_gene_expression_filename = "";
			public string tumor_dna_gene_coverage_filname = "";
            public string vcf_filename = "";
            public string extracted_maf_lines_filename = "";
            public string normal_dna_mapped_base_count_filename = "";
            public string tumor_dna_mapped_base_count_filename = "";
            public string normal_rna_mapped_base_count_filename = "";
            public string tumor_rna_mapped_base_count_filename = "";
            public string selected_variant_counts_by_gene_filename = "";
            public string tumor_regional_methylation_filename = "";
            public string normal_regional_methylation_filename = "";
            // If you add another drived file type and it has a **done** terminator, please add it to the CheckDone tool.     

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

            //
            // Sizes for derived files.
            //
            public long normal_dna_allcount_size = 0;
            public long tumor_dna_allcount_size = 0;
            public long normal_rna_allcount_size = 0;
            public long tumor_rna_allcount_size = 0;
            public long regional_expression_size = 0;
            public long gene_expression_size = 0;
            public long selected_variants_size = 0;
            public long tumor_dna_reads_at_selected_variants_size = 0;
            public long tumor_dna_reads_at_selected_variants_index_size = 0;
            public long normal_dna_reads_at_selected_variants_size = 0;
            public long normal_dna_reads_at_selected_variants_index_size = 0;
            public long tumor_rna_reads_at_selected_variants_size = 0;
            public long tumor_rna_reads_at_selected_variants_index_size = 0;
            public long normal_rna_reads_at_selected_variants_size = 0;
            public long normal_rna_reads_at_selected_variants_index_size = 0;
            public long annotated_selected_variants_size = 0;
            public long normal_allele_specific_gene_expression_size = 0;
            public long tumor_allele_specific_gene_expression_size = 0;
            public long tumor_dna_gene_coverage_size = 0;
            public long vcf_size = 0;
            public long extracted_maf_lines_size = 0;
            public long normal_dna_mapped_base_count_size = 0;
            public long tumor_dna_mapped_base_count_size = 0;
            public long normal_rna_mapped_base_count_size = 0;
            public long tumor_rna_mapped_base_count_size = 0;
            public long selected_variant_counts_by_gene_size = 0;
            public long tumor_regional_methylation_size = 0;
            public long normal_regional_methylation_size = 0;

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

            public static FieldInformation[] AllFields = 
            {
                new FieldInformation("Case ID",                                             c => c.case_id, (c, v) => c.case_id = v),
                new FieldInformation("Normal DNA File ID",                                  c => c.normal_dna_file_id, (c, v) => c.normal_dna_file_id = v),
                new FieldInformation("Tumor DNA File ID",                                   c => c.tumor_dna_file_id, (c,v) => c.tumor_dna_file_id = v),
                new FieldInformation("Normal RNA File ID",                                  c => c.normal_rna_file_id, (c,v) => c.normal_rna_file_id = v),
                new FieldInformation("Tumor RNA File ID",                                   c => c.tumor_rna_file_id, (c,v) => c.tumor_rna_file_id = v),
                new FieldInformation("MAF File ID",                                         c => c.maf_file_id, (c,v) => c.maf_file_id = v),
                new FieldInformation("Tumor Methylation File ID",                           c => c.tumor_methylation_file_id, (c,v) => c.tumor_methylation_file_id = v),
				new FieldInformation("Normal Methylation File ID",                          c => c.normal_methylation_file_id, (c,v) => c.normal_methylation_file_id = v),
				new FieldInformation("Tumor Copy Number File ID",                           c => c.tumor_copy_number_file_id, (c,v) => c.tumor_copy_number_file_id = v),
				new FieldInformation("Normal Copy Number File ID",                          c => c.normal_copy_number_file_id, (c,v) => c.normal_copy_number_file_id = v),
				new FieldInformation("Project ID",                                          c => c.project_id, (c,v) => c.project_id = v),
                new FieldInformation("Sample IDs",                                          c => c.sampleIdsInCommaSeparatedList(), (c,v) => c.sample_ids = v.Split(',').ToList()),

                new FieldInformation("Normal DNA Filename",                                 c => c.normal_dna_filename, (c,v) => c.normal_dna_filename = v),
                new FieldInformation("Tumor DNA Filename",                                  c => c.tumor_dna_filename, (c,v) => c.tumor_dna_filename = v),
                new FieldInformation("Normal RNA Filename",                                 c => c.normal_rna_filename, (c,v) => c.normal_rna_filename = v),
                new FieldInformation("Tumor RNA Filename",                                  c => c.tumor_rna_filename, (c,v) => c.tumor_rna_filename = v),
                new FieldInformation("Tumor Methylation Filename",                          c => c.tumor_methylation_filename, (c,v) => c.tumor_methylation_filename = v),
				new FieldInformation("Normal Methylation Filename",                         c => c.normal_methylation_filename, (c,v) => c.normal_methylation_filename = v),
				new FieldInformation("Tumor Copy Number Filename",                          c => c.tumor_copy_number_filename, (c,v) => c.tumor_copy_number_filename = v),
				new FieldInformation("Normal Copy Number Filename",                         c => c.normal_copy_number_filename, (c,v) => c.normal_copy_number_filename = v),
				new FieldInformation("MAF Filename",                                        c => c.maf_filename, (c,v) => c.maf_filename = v),

                new FieldInformation("Normal DNA Size",                                     c => Convert.ToString(c.normal_dna_size), (c,v) => c.normal_dna_size = LongFromString(v)),
                new FieldInformation("Tumor DNA Size",                                      c => Convert.ToString(c.tumor_dna_size), (c,v) => c.tumor_dna_size = LongFromString(v)),
                new FieldInformation("Normal RNA Size",                                     c => Convert.ToString(c.normal_rna_size), (c,v) => c.normal_rna_size = LongFromString(v)),
                new FieldInformation("Tumor RNA Size",                                      c => Convert.ToString(c.tumor_rna_size), (c,v) => c.tumor_rna_size = LongFromString(v)),
                new FieldInformation("Tumor Methylation Size",                              c => Convert.ToString(c.tumor_methylation_size), (c,v) => c.tumor_methylation_size = LongFromString(v)),
				new FieldInformation("Normal Methylation Size",                             c => Convert.ToString(c.normal_methylation_size), (c,v) => c.normal_methylation_size = LongFromString(v)),
				new FieldInformation("Tumor Copy Number Size",                              c => Convert.ToString(c.tumor_copy_number_size), (c,v) => c.tumor_copy_number_size = LongFromString(v)),
				new FieldInformation("Normal Copy Number Size",                             c => Convert.ToString(c.normal_copy_number_size), (c,v) => c.normal_copy_number_size = LongFromString(v)),

				new FieldInformation("Normal DNA Allcount Filename",                        c => c.normal_dna_allcount_filename, (c,v) => c.normal_dna_allcount_filename = v, DerivedFile.Type.NormalDNAAllcount, normalDNAAllcountExtension, c => c.normal_dna_file_id, "Normal DNA Allcount File Size", c => c.normal_dna_allcount_size, (c,v) => c.normal_dna_allcount_size = v),
                new FieldInformation("Tumor DNA Allcount Filename",                         c => c.tumor_dna_allcount_filename, (c,v) => c.tumor_dna_allcount_filename = v, DerivedFile.Type.TumorDNAAllcount, tumorDNAAllcountExtension, c => c.tumor_dna_file_id, "Tumor DNA Allcount File Size", c => c.tumor_dna_allcount_size, (c, v) => c.tumor_dna_allcount_size = v),
                new FieldInformation("Normal RNA Allcount Filename",                        c => c.normal_rna_allcount_filename, (c,v) => c.normal_rna_allcount_filename = v, DerivedFile.Type.NormalRNAAllcount, normalRNAAllcountExtension, c => c.normal_rna_file_id, "Normal RNA Allcount File Size", c => c.normal_rna_allcount_size, (c, v) => c.normal_rna_allcount_size = v),
                new FieldInformation("Tumor RNA Allcount Filename",                         c => c.tumor_rna_allcount_filename, (c,v) => c.tumor_rna_allcount_filename = v, DerivedFile.Type.TumorRNAAllcount, tumorRNAAllcountExtension, c => c.tumor_rna_file_id, "Tumor RNA Allcount File Size", c => c.tumor_rna_allcount_size, (c, v) => c.tumor_rna_allcount_size = v),
                new FieldInformation("Regional Expression Filename",                        c => c.regional_expression_filename, (c,v) => c.regional_expression_filename = v, DerivedFile.Type.RegionalExpression, regionalExpressionExtension, c => c.tumor_rna_file_id, "Regional Expression File Size", c => c.regional_expression_size, (c, v) => c.regional_expression_size = v),
                new FieldInformation("Gene Expression Filename",                            c => c.gene_expression_filename, (c,v) => c.gene_expression_filename = v, DerivedFile.Type.GeneExpression, geneExpressionExtension, c => c.case_id, "Gene Expression File Size", c => c.gene_expression_size, (c, v) => c.gene_expression_size = v),
                new FieldInformation("Selected Variants Filename",                          c => c.selected_variants_filename, (c,v) => c.selected_variants_filename = v, DerivedFile.Type.SelectedVariants, selectedVariantsExtension, c => c.normal_dna_file_id, "Selected Variants File Size", c=> c.selected_variants_size, (c, v) => c.selected_variants_size = v),
                new FieldInformation("Normal DNA Reads At Selected Variants Filename",      c => c.normal_dna_reads_at_selected_variants_filename, (c,v) => c.normal_dna_reads_at_selected_variants_filename = v, DerivedFile.Type.NormalDNAReadsAtSelectedVariants, normalDNAReadsAtSelectedVariantsExtension, c => c.normal_dna_file_id, "Normal DNA Reads At Selected Variants File Size", c => c.normal_dna_reads_at_selected_variants_size, (c, v) => c.normal_dna_reads_at_selected_variants_size = v),
                new FieldInformation("Normal DNA Reads At Selected Variants Index Filename",c => c.normal_dna_reads_at_selected_variants_index_filename, (c,v) => c.normal_dna_reads_at_selected_variants_index_filename = v, DerivedFile.Type.NormalDNAReadsAtSelectedVariantsIndex, normalDNAReadsAtSelectedVariantsIndexExtension, c => c.normal_dna_file_id, "Normal DNA Reads At Selected Variants Index File Size", c => c.normal_dna_reads_at_selected_variants_index_size, (c, v) => c.normal_dna_reads_at_selected_variants_index_size = v),
                new FieldInformation("Tumor DNA Reads At Selected Variants Filename",       c => c.tumor_dna_reads_at_selected_variants_filename, (c,v) => c.tumor_dna_reads_at_selected_variants_filename = v, DerivedFile.Type.TumorDNAReadsAtSelectedVariants, tumorDNAReadsAtSelectedVariantsExtension, c => c.tumor_dna_file_id, "Tumor DNA Reads At Selected Variants File Size", c => c.tumor_dna_reads_at_selected_variants_size, (c, v) => c.tumor_dna_reads_at_selected_variants_size = v),
                new FieldInformation("Tumor DNA Reads At Selected Variants Index Filename", c => c.tumor_dna_reads_at_selected_variants_index_filename, (c,v) => c.tumor_dna_reads_at_selected_variants_index_filename = v, DerivedFile.Type.TumorDNAReadsAtSelectedVariantsIndex, tumorDNAReadsAtSelectedVariantsIndexExtension, c => c.tumor_dna_file_id, "Tumor DNA Reads At Selected Variants Index File Size", c => c.tumor_dna_reads_at_selected_variants_index_size, (c, v) => c.tumor_dna_reads_at_selected_variants_index_size = v),
                new FieldInformation("Normal RNA Reads At Selected Variants Filename",      c => c.normal_rna_reads_at_selected_variants_filename, (c,v) => c.normal_rna_reads_at_selected_variants_filename = v, DerivedFile.Type.NormalRNAReadsAtSelectedVariants, normalRNAReadsAtSelectedVariantsExtension, c => c.normal_rna_file_id, "Normal RNA Reads At Selected Variants File Size", c => c.normal_rna_reads_at_selected_variants_size, (c, v) => c.normal_rna_reads_at_selected_variants_size = v),
                new FieldInformation("Normal RNA Reads At Selected Variants Index Filename",c => c.normal_rna_reads_at_selected_variants_index_filename, (c,v) => c.normal_rna_reads_at_selected_variants_index_filename = v, DerivedFile.Type.NormalRNAReadsAtSelectedVariantsIndex, normalRNAReadsAtSelectedVariantsIndexExtension, c => c.normal_rna_file_id, "Normal RNA Reads At Selected Variants Index File Size", c => c.normal_rna_reads_at_selected_variants_index_size, (c, v) => c.normal_rna_reads_at_selected_variants_index_size = v),
                new FieldInformation("Tumor RNA Reads At Selected Variants Filename",       c => c.tumor_rna_reads_at_selected_variants_filename, (c,v) => c.tumor_rna_reads_at_selected_variants_filename = v, DerivedFile.Type.TumorRNAReadsAtSelectedVariants, tumorRNAReadsAtSelectedVariantsExtension, c => c.tumor_rna_file_id, "Tumor RNA Reads At Selected Variants File Size", c => c.tumor_rna_reads_at_selected_variants_size, (c, v) => c.tumor_rna_reads_at_selected_variants_size = v),
                new FieldInformation("Tumor RNA Reads At Selected Variants Index Filename", c => c.tumor_rna_reads_at_selected_variants_index_filename, (c,v) => c.tumor_rna_reads_at_selected_variants_index_filename = v, DerivedFile.Type.TumorRNAReadsAtSelectedVariantsIndex, tumorRNAReadsAtSelectedVariantsIndexExtension, c => c.tumor_rna_file_id, "Tumor RNA Reads At Selected Variants Index File Size", c => c.tumor_rna_reads_at_selected_variants_index_size, (c, v) => c.tumor_rna_reads_at_selected_variants_index_size = v),
                new FieldInformation("Annotated Selected Variants Filename",                c => c.annotated_selected_variants_filename, (c,v) => c.annotated_selected_variants_filename = v, DerivedFile.Type.AnnotatedSelectedVariants, annotatedSelectedVariantsExtension, c => c.case_id, "Annotated Selected Variants File Size", c => c.annotated_selected_variants_size, (c, v) => c.annotated_selected_variants_size = v),
                new FieldInformation("Normal Allele Specific Gene Expression Filename",     c => c.normal_allele_specific_gene_expression_filename, (c,v) => c.normal_allele_specific_gene_expression_filename = v, DerivedFile.Type.NormalAlleleSpecificGeneExpression, normalAlleleSpecificGeneExpressionExtension, c => c.case_id, "Normal Allele Specific Gene Expression File Size", c => c.normal_allele_specific_gene_expression_size, (c, v) => c.normal_allele_specific_gene_expression_size = v),
                new FieldInformation("Tumor Allele Specific Gene Expression Filename",      c => c.tumor_allele_specific_gene_expression_filename, (c,v) => c.tumor_allele_specific_gene_expression_filename = v, DerivedFile.Type.TumorAlleleSpecificGeneExpression, tumorAlleleSpecificGeneExpressionExtension, c => c.case_id, "Tumor Allele Specific Gene Expression File Size", c => c.tumor_allele_specific_gene_expression_size, (c, v) => c.tumor_allele_specific_gene_expression_size = v),
                new FieldInformation("Tumor DNA Gene Coverage Filename",                    c => c.tumor_dna_gene_coverage_filname, (c,v) => c.tumor_dna_gene_coverage_filname = v, DerivedFile.Type.TumorDNAGeneCoverage, tumorDNAGeneCoverageExtension, c => c.tumor_dna_file_id, "Tumor DNA Gene Coverage File Size", c => c.tumor_dna_gene_coverage_size, (c, v) => c.tumor_dna_gene_coverage_size = v),
                new FieldInformation("VCF Filename",                                        c => c.vcf_filename, (c,v) => c.vcf_filename = v, DerivedFile.Type.VCF, vcfExtension, c => c.normal_dna_file_id, "VCF File Size", c => c.vcf_size, (c, v) => c.vcf_size = v),
                new FieldInformation("Extracted MAF Lines Filename",                        c => c.extracted_maf_lines_filename, (c,v) => c.extracted_maf_lines_filename = v, DerivedFile.Type.ExtractedMAFLines, extractedMAFLinesExtension, c => c.case_id, "Extracted MAF Lines File Size", c => c.extracted_maf_lines_size, (c, v) => c.extracted_maf_lines_size = v),
                new FieldInformation("Normal DNA Mapped Base Count Filename",               c => c.normal_dna_mapped_base_count_filename, (c, v) => c.normal_dna_mapped_base_count_filename = v, DerivedFile.Type.NormalDNAMappedBaseCount, normalDNAMappedBaseCountExtension, c => c.normal_dna_file_id, "Normal DNA Mapped Base Count File Size", c => c.normal_dna_mapped_base_count_size, (c, v) => c.normal_dna_mapped_base_count_size = v),
                new FieldInformation("Tumor DNA Mapped Base Count Filename",                c => c.tumor_dna_mapped_base_count_filename, (c, v) => c.tumor_dna_mapped_base_count_filename = v, DerivedFile.Type.TumorDNAMappedBaseCount, tumorDNAMappedBaseCountExtension, c => c.tumor_dna_file_id, "Tumor DNA Mapped Base Count File Size", c => c.tumor_dna_mapped_base_count_size, (c, v) => c.tumor_dna_mapped_base_count_size = v),
                new FieldInformation("Normal RNA Mapped Base Count Filename",               c => c.normal_rna_mapped_base_count_filename, (c, v) => c.normal_rna_mapped_base_count_filename = v, DerivedFile.Type.NormalRNAMappedBaseCount, normalRNAMappedBaseCountExtension, c => c.normal_rna_file_id, "Normal RNA Mapped Base Count File Size", c => c.normal_rna_mapped_base_count_size, (c, v) => c.normal_rna_mapped_base_count_size = v),
                new FieldInformation("Tumor RNA Mapped Base Count Filename",                c => c.tumor_rna_mapped_base_count_filename, (c, v) => c.tumor_rna_mapped_base_count_filename = v, DerivedFile.Type.TumorRNAMappedBaseCount, tumorRNAMappedBaseCountExtension, c => c.tumor_rna_file_id, "Tumor RNA Mapped Base Count File Size", c => c.tumor_rna_mapped_base_count_size, (c, v) => c.tumor_rna_mapped_base_count_size = v),
                new FieldInformation("Selected Variant Counts By Gene Filename",            c => c.selected_variant_counts_by_gene_filename, (c, v) => c.selected_variant_counts_by_gene_filename = v, DerivedFile.Type.SelectedVariantCountByGene, selectedVariantCountByGeneExtension, c => c.case_id, "Selected Variant Counts By Gene File Size", c => c.selected_variant_counts_by_gene_size, (c, v) => c.selected_variant_counts_by_gene_size = v),
                new FieldInformation("Tumor Regional Methylation Filename",                 c => c.tumor_regional_methylation_filename, (c, v) => c.tumor_regional_methylation_filename = v, DerivedFile.Type.TumorRegionalMethylation, tumorRegionalMethylationExtension, c => c.case_id, "Tumor Regional Methylation File Size", c => c.tumor_regional_methylation_size, (c, v) => c.tumor_regional_methylation_size = v),
                new FieldInformation("Normal Regional Methylation Filename",                c => c.normal_regional_methylation_filename, (c, v) => c.normal_regional_methylation_filename = v, DerivedFile.Type.NormalRegionalMethylation, normalRegionalMethylationExtension, c => c.case_id, "Normal Regional Methylation File Size", c => c.normal_regional_methylation_size, (c, v) => c.normal_regional_methylation_size = v),

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

			}; // fieldInformation

  
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

            public static void loadAllFileLocations(Dictionary<string, Case> cases, Dictionary<string, DownloadedFile> downloadedFiles, Dictionary<string, List<DerivedFile>> derivedFiles)
            {
                foreach (var caseEntry in cases)
                {
                    caseEntry.Value.loadFileLocations(downloadedFiles, derivedFiles);
                }
            }

            //
            // Fill in all the derived file locations if they exist.
            //
            public void loadFileLocations(Dictionary<string, DownloadedFile> downloadedFiles, Dictionary<string, List<DerivedFile>> derivedFiles)
            {

                if (downloadedFiles.ContainsKey(normal_dna_file_id))
                {
                    normal_dna_filename = downloadedFiles[normal_dna_file_id].fileInfo.FullName;
                }
                else
                {
                    normal_dna_filename = "";
                }

                if (downloadedFiles.ContainsKey(tumor_dna_file_id))
                {
                  tumor_dna_filename = downloadedFiles[tumor_dna_file_id].fileInfo.FullName;
                }
                else
                {
                    tumor_dna_filename = "";
                }

                if (normal_rna_file_id != "" && downloadedFiles.ContainsKey(normal_rna_file_id))
                {
                  normal_rna_filename = downloadedFiles[normal_rna_file_id].fileInfo.FullName;
                }
                else
                {
                    normal_rna_filename = "";
                }

                if (downloadedFiles.ContainsKey(tumor_rna_file_id))
                {
                  tumor_rna_filename = downloadedFiles[tumor_rna_file_id].fileInfo.FullName;
                }
                else
                {
                    tumor_rna_filename = "";
                }

                if (downloadedFiles.ContainsKey(maf_file_id))
                {
                  maf_filename  = downloadedFiles[maf_file_id].fileInfo.FullName;
                }
                else
                {
                    maf_filename = "";
                }

                if (tumor_methylation_file_id != "" && downloadedFiles.ContainsKey(tumor_methylation_file_id))
                {
                  tumor_methylation_filename  = downloadedFiles[tumor_methylation_file_id].fileInfo.FullName;
                }
                else
                {
                    tumor_methylation_filename = "";
                }

				if (normal_methylation_file_id != "" && downloadedFiles.ContainsKey(normal_methylation_file_id))
				{
					normal_methylation_filename = downloadedFiles[normal_methylation_file_id].fileInfo.FullName;
				}
				else
				{
					normal_methylation_filename = "";
				}

				if (tumor_copy_number_file_id != "" && downloadedFiles.ContainsKey(tumor_copy_number_file_id))
                {
                    tumor_copy_number_filename = downloadedFiles[tumor_copy_number_file_id].fileInfo.FullName;
                }
                else
                {
                    tumor_copy_number_filename = "";
                }

				if (normal_copy_number_file_id != "" && downloadedFiles.ContainsKey(normal_copy_number_file_id))
				{
					normal_copy_number_filename = downloadedFiles[normal_copy_number_file_id].fileInfo.FullName;
				}
				else
				{
					normal_copy_number_filename = "";
				}

				if (!derivedFiles.ContainsKey(case_id))
                {
                    tumor_rna_allcount_filename = "";
                    normal_allele_specific_gene_expression_filename = "";
					tumor_allele_specific_gene_expression_filename = "";
					annotated_selected_variants_filename = "";
                    tumor_dna_reads_at_selected_variants_filename = "";
                    tumor_dna_reads_at_selected_variants_index_filename = "";
                    gene_expression_filename = "";
                    tumor_dna_reads_at_selected_variants_filename = "";
                    regional_expression_filename = "";
                    tumor_rna_reads_at_selected_variants_filename = "";
                    tumor_rna_reads_at_selected_variants_index_filename = "";
                    selected_variants_filename = "";
                    vcf_filename = "";
                    return;
                }

                Dictionary<DerivedFile.Type, List<DerivedFile>> derivedFilesByType = new Dictionary<DerivedFile.Type,List<DerivedFile>>();
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

        public static StreamWriter CreateStreamWriterWithRetry(string filename)
        {
            while (true)
            {
                try
                {
                    var writer = new StreamWriter(filename);
                    return writer;
                }
                catch (IOException e)
                {
                    if (e is FileNotFoundException || e is DirectoryNotFoundException)
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
                    if (e is FileNotFoundException || e is DirectoryNotFoundException)
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
                    if (e is FileNotFoundException || e is DirectoryNotFoundException)
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
                catch (IOException)
                {
                    Console.WriteLine("IOException reading " + filename + ".  Sleeping and retrying.");
                    Thread.Sleep(10 * 1000);
                }
            }
        }

		public class Configuration
        {
            public const string defaultBaseDirectory = @"\\msr-genomics-0\d$\gdc\";
            public const string defaultConfigurationFilePathame = defaultBaseDirectory + "configuration.txt";

            public string accessTokenPathname = defaultBaseDirectory +  @"access_token.txt";

			public const string defaultGenomeBuild = "hg38";
			public const string defaultGeneLocationInformationFilename = defaultBaseDirectory + "knownGene-" + defaultGenomeBuild + ".txt";

            public List<string> dataDirectories = new List<string>();
            public string mafManifestPathname = defaultBaseDirectory + "mafManifest.txt";
            public string mutationCaller = "mutect";
            public List<string> programNames = new List<string>();
            public string binariesDirectory = defaultBaseDirectory + @"bin\";
            public string configuationFilePathname = defaultConfigurationFilePathame;
            public string casesFilePathname = defaultBaseDirectory + "cases.txt";
            public string indexDirectory = defaultBaseDirectory + @"indices\hg38-20";
            public string derivedFilesDirectory = "derived_files";    // This is relative to each download directory
            public string hpcScriptFilename = "";    // The empty string says to black hole this script
            public string hpcScheduler = "gcr";
            public string hpcBinariesDirectory = @"\\gcr\scratch\b99\bolosky\";
            public string hpcIndexDirectory = @"\\msr-genomics-0\d$\gdc\indices\hg38-20";
            public string azureScriptFilename = ""; // The empty string says to black hole this script
            public string expressionFilesDirectory = @"\\msr-genomics-0\d$\gdc\expression\";
            public string completedVCFsDirectory = "";  // Where completed VCFs from Azure are dropped
            public List<string> excludedFileIDs = new List<string>();
            public int regionalExpressionRegionSize = 1000;
            public int nTumorsToIncludeGene = 30;   // How many tumors must have at least one mutation in a gene in order to include it.
            public int nReadsToIncludeVariant = 30; // How much coverage do we need of a mutation or somatic variant to consider it?
            public string selectedGenesFilename = @"\\msr-genomics-0\d$\gdc\seleted_genes.txt";
            public string scriptOutputDirectory = @"\temp\";
            public string finalResultsDirectory = @"\\msr-genomics-0\d$\gdc\final_results\";
            public double significanceLevel = 0.01;
            public string geneLocationInformationFilename = defaultGeneLocationInformationFilename;
            public int minSampleCountForHeatMap = 100;

            public string geneScatterGraphsDirectory = defaultBaseDirectory + @"gene_scatter_graphs\";
            public string regionalExpressionDirectory = defaultBaseDirectory + @"regional_expression\";

			public const string unfilteredCountsDirectory = defaultBaseDirectory + @"gene_mutations_with_counts\";
			public const string unfilteredCountsExtention = @"_unfiltered_counts.txt";
			public const string methylationREFsFilename = defaultBaseDirectory + "compositeREFs.txt";


			public string[] commandLineArgs = null;    // The args excluding -configuration <filename>

            Configuration()
            {
                programNames.Add("TCGA");   // The default value
                dataDirectories.Add(defaultBaseDirectory + @"downloaded_files\");
            }

            //
            // Parse the args to find the configuration file pathname, and then load from that path (or the default if it's not present).
            //
            public static Configuration loadFromFile(string [] args) 
            {
                string pathname = @"\\msr-genomics-0\d$\gdc\configuration.txt";

                bool fromCommandLine = false;
                var nonConsumedArgs = new List<string>();
                for (int i = 0; i < args.Count(); i++) {
                    if (args[i] == "-configuration") {
                        if (i >= args.Count() - 1) {
                            Console.WriteLine("-configuation can't be the last parameter, it needs to be followed by a file path.  Ignoring.");
                        } else {
                            pathname = args[i+1];
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

                foreach (var line in lines) 
                {
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
                    } else {
                        Console.WriteLine("ASEConfiguration.loadFromFile: configuration file " + pathname + " contains a line with an unknown configuration parameter type: " + line + ".  Ignoring.");
                        continue;
                    }
                }

                return retVal;
            }

        }

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
            public readonly int []tumorsByCountOfNonSilentMutations;
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
                Console.WriteLine("ShareFromPathname: can't parse pathname " + pathname);
                return "";
            }

            string share = "";
            for (int i = 1; i < 4; i++) // Skip 0 because we know it's the empty string, and it makes the \ part easier.
            {
                share += @"\" + fields[i];
            }

            return share;
        }

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


        public const string tumorRNAAllcountExtension = ".allcount.gz";
        public const string normalRNAAllcountExtension = ".normal_rna_allcount.gz";
        public const string normalDNAAllcountExtension = ".normal_dna_allcount.gz";
        public const string tumorDNAAllcountExtension = ".tumor_dna_allcount.gz";
        public const string selectedVariantsExtension = ".selectedVariants";
        public const string annotatedSelectedVariantsExtension = ".annotatedSeletedVariants";
        public const string regionalExpressionExtension = ".regional_expression.txt";
        public const string geneExpressionExtension = ".gene_expression.txt";
		public const string normalAlleleSpecificGeneExpressionExtension = ".normal-allele-specific_gene_expression.txt";
		public const string tumorAlleleSpecificGeneExpressionExtension = ".tumor-allele-specific_gene_expression.txt";
        public const string tumorRNAReadsAtSelectedVariantsExtension = ".tumor-rna-reads-at-selected-variants.txt";
        public const string tumorRNAReadsAtSelectedVariantsIndexExtension = ".tumor-rna-reads-at-selected-variants.txt.index";
        public const string normalRNAReadsAtSelectedVariantsExtension = ".normal-rna-reads-at-selected-variants.txt";
        public const string normalRNAReadsAtSelectedVariantsIndexExtension = ".normal-rna-reads-at-selected-variants.txt.index";
        public const string normalDNAReadsAtSelectedVariantsExtension = ".normal-dna-reads-at-selected-variants.txt";
        public const string normalDNAReadsAtSelectedVariantsIndexExtension = ".normal-dna-reads-at-selected-variants.txt.index";
        public const string tumorDNAReadsAtSelectedVariantsExtension = ".tumor-dna-reads-at-selected-variants.txt";
        public const string tumorDNAReadsAtSelectedVariantsIndexExtension = ".tumor-dna-reads-at-selected-variants.txt.index";
        public const string vcfExtension = ".vcf";
        public const string extractedMAFLinesExtension = ".extracted_maf_lines.txt";
        public const string tumorDNAGeneCoverageExtension = ".tumor_dna_gene_coverage.txt";
        public const string normalDNAMappedBaseCountExtension = ".normal_dna_mapped_base_count.txt";
        public const string tumorDNAMappedBaseCountExtension = ".tumor_dna_mapped_base_count.txt";
        public const string normalRNAMappedBaseCountExtension = ".normal_rna_mapped_base_count.txt";
        public const string tumorRNAMappedBaseCountExtension = ".tumor_rna_mapped_base_count.txt";
        public const string selectedVariantCountByGeneExtension = ".selected_variant_count_by_gene.txt";
        public const string tumorRegionalMethylationExtension = ".tumor_regional_methylation.txt";
        public const string normalRegionalMethylationExtension = ".normal_regional_methylation.txt";
        public const string bonferroniExtension = "_bonferroni.txt";

        public const string scatterGraphsSummaryFilename = "_summary.txt";
        public const string mannWhitneyFilename = "_MannWhitney.txt";
        public const string genesWithSelectedVariantsFilename = "GenesWithSelectedVariantCounts.txt";
        public const string heatMapFilename = "AlleleSpecificExpressionHeatMap.txt";
        public const string tumorHeatMapHistogramFilename = "TumorAlleleSpecificExpressionHeatMapHistograms.txt";
        public const string normalHeatMapHistogramFilename = "NormalAlleleSpecificExpressionHeatMapHistograms.txt";
        public const string ASEConsistencyFilename = "ASEConsistency.txt";
        public const string GenesByFunnyASECountFilename = "ASEInconsistencyByGene.txt";

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

                derived_from_file_id = filename.Substring(0,FileIdLength);

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


            public enum Type { Unknown, NormalRNAAllcount, TumorRNAAllcount, NormalDNAAllcount, TumorDNAAllcount, RegionalExpression, GeneExpression, TumorDNAGeneCoverage,
                SelectedVariants, NormalDNAReadsAtSelectedVariants, NormalDNAReadsAtSelectedVariantsIndex, TumorDNAReadsAtSelectedVariants, TumorDNAReadsAtSelectedVariantsIndex, TumorRNAReadsAtSelectedVariants,
                TumorRNAReadsAtSelectedVariantsIndex, NormalRNAReadsAtSelectedVariants, NormalRNAReadsAtSelectedVariantsIndex, AnnotatedSelectedVariants, NormalAlleleSpecificGeneExpression, TumorAlleleSpecificGeneExpression, VCF, ExtractedMAFLines,
                NormalDNAMappedBaseCount, TumorDNAMappedBaseCount, NormalRNAMappedBaseCount, TumorRNAMappedBaseCount, SelectedVariantCountByGene, TumorRegionalMethylation, NormalRegionalMethylation,
            };
        } // DerivedFile

        class ScanFilesystemState
        {
            public ulong totalFreeBytes = 0;
            public ulong totalStorageBytes = 0;
            public ulong totalBytesInDownloadedFiles = 0;
            public ulong totalBytesInDerivedFiles = 0;
            public int nDownloadedFiles = 0;
            public int nDerivedFiles = 0;

            public List<DownloadedFile> downloadedFiles = new List<DownloadedFile>();
            public List<DerivedFile> derivedFiles = new List<DerivedFile>();
        }

        static void ScanOneFilesystem(Configuration configuration, string downloadedFilesDirectory, ScanFilesystemState state, Stopwatch stopwatch, int directoryFieldLength) 
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

            var downloadedFiles = new List<DownloadedFile>();
            var derivedFiles = new List<DerivedFile>();

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
                //
                string candidatePathname = null;
                string md5Pathname = null;
                bool sawBAI = false;
                bool sawBAM = false;
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

                if (candidatePathname == null && md5Pathname != null)
                {
                    Console.WriteLine("Found md5 file in directory without accompanying downloaded file (?): " + md5Pathname);
                    continue;
                }

                if (sawBAM != sawBAI)
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
            } // foreach subdir

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
 
                foreach (var derivedFilePathname in Directory.EnumerateFiles(derivedCase))
                {
                    var derivedFile = new DerivedFile(derivedFilePathname, caseId);
                    nDerivedFiles++;
                    totalBytesInDerivedFiles += (ulong)derivedFile.fileinfo.Length;
                    derivedFiles.Add(derivedFile);
                }
            }

            ulong freeBytesAvailable, totalBytes, totalNumberOfFreeBytes;
            GetDiskFreeSpaceEx(downloadedFilesDirectory, out freeBytesAvailable, out totalBytes, out totalNumberOfFreeBytes);

            lock (state)
            {
                Console.WriteLine(String.Format("{0," + directoryFieldLength + "}", downloadedFilesDirectory) + " " + String.Format("{0,16}", "" + nDownloadedFiles + " (" + SizeToUnits(totalBytesInDownloadedFiles) + "B)") + " " +
                    String.Format("{0,13}", "" + nDerivedFiles + " (" + SizeToUnits(totalBytesInDerivedFiles) + "B)") + " " + String.Format("{0,10}", SizeToUnits(freeBytesAvailable) +"B") + " " +
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
        public static void ScanFilesystems(Configuration configuration, out Dictionary<string, DownloadedFile> downloadedFiles, out Dictionary<string, List<DerivedFile>> derivedFiles)
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

            Console.WriteLine(" Downloaded Files Derived Files Free Size Scan Time");
            for (int i = 0; i < directoryHeader.Count() + paddingLength; i++ )
            {
                Console.Write("-");
            }
            Console.WriteLine(" ---------------- ------------- ---------- ---------");

            var state = new ScanFilesystemState();
            var threads = new List<Thread>();

            foreach (var directory in configuration.dataDirectories)
            {
                threads.Add(new Thread(() => ScanOneFilesystem(configuration, directory, state, stopwatch, paddingLength + directoryHeader.Count())));
            } // foreach data directory

            threads.ForEach(t => t.Start());
            threads.ForEach(t => t.Join());

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
                if (lines[0].Count() < headerPrefix.Count() || lines[0].Substring(0,headerPrefix.Count()) != headerPrefix) {
                    Console.WriteLine("Corrupt or unrecognized version in maf manifest header, ignoring: " + lines[0]);
                    return null;
                }

                var retVal = new Dictionary<string, MAFInfo>();

                if (lines[lines.Count() -1] != "**done**") {
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
        }

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
        }

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

            public HeaderizedFile(StreamReader inputFile_, bool hasVersion_, bool hasDone_, string expectedVersion_, List<string> wantedFields_, bool skipHash_ = false, bool allowMissingColumnsInData_ = false, int skippedRows_ = 0)
            {
                inputFile = inputFile_;
                hasVersion = hasVersion_;
                hasDone = hasDone_;
                expectedVersion = expectedVersion_;
                wantedFields = wantedFields_;
                skipHash = skipHash_;
                allowMissingColumnsInData = allowMissingColumnsInData_;
				skippedRows = skippedRows_;
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
                fieldMappings_out = null;

                if (hasVersion)
                {
                    var versionString = inputFile.ReadLine();

                    if (versionString == null || expectedVersion != null && versionString != expectedVersion)
                    {
                        result = null;
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
                    result = null;
                    return false;
                }

                var columns = header.Split('\t');
                var fieldMappings = new Dictionary<string, int>();
                int maxNeededField = -1;

                for (int i = 0; i < columns.Count(); i++)
                {
                    if (wantedFields.Contains(columns[i]))
                    {
                        if (fieldMappings.ContainsKey(columns[i]))
                        {
                            Console.WriteLine("Duplicate needed column in headerized file (or code bug or something): " + columns[i]);
                            result = null;
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
                        result = null;
                        return false;

                    }
                    Console.Write("Headerized file: missing columns:");
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
                result = new List<outputType>();
                while (null != (inputLine = inputFile.ReadLine())) {
                    if (sawDone) {
                        Console.WriteLine("HeaderizedFile: Saw data after **done**");
                        result = null;
                        return false;
                    }

                    if ("**done**" == inputLine) {
                        sawDone = true;
                        continue;
                    }

                    var fields = inputLine.Split('\t');
                    if (fields.Count() <= maxNeededField && !allowMissingColumnsInData)
                    {
                        Console.WriteLine("HeaderizedFile.Parse: input line didn't include a needed field " + inputLine);
                        result = null;
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
                        result.Add(parser(fieldMappings, fields));
                    } else
                    {
                        result.Add(fieldGrabbingParser(new FieldGrabber(fieldMappings, fields)));
                    }
                }

                if (hasDone && !sawDone)
                {
                    Console.WriteLine("HeaderizedFile.Parse: missing **done**");
                    result = null;
                    return false;
                }
                else if (!hasDone && sawDone)
                {
                    Console.WriteLine("Saw unepected **done**.  Ignoring.");
                    result = null;
                    return false;
                }


                fieldMappings_out = fieldMappings;
                return true;
            } // ParseFile

            public class FieldGrabber
            {
                public FieldGrabber(Dictionary<string, int> fieldMappings_, string[] fields_)
                {
                    fieldMappings = fieldMappings_;
                    fields = fields_;
                }

                public string AsString(string fieldName)
                {
                    return ConvertToNonExcelString(fields[fieldMappings[fieldName]]);
                }
                public int AsInt(string fieldName)
                {
                    return Convert.ToInt32(AsString(fieldName));
                }

                public int AsIntMinusOneIfStarOrEmptyString(string fieldName)
                {
                    if (fields[fieldMappings[fieldName]] == "*" || fields[fieldMappings[fieldName]] == "")
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
                    if (fields[fieldMappings[fieldName]] == "*" || fields[fieldMappings[fieldName]] == "")
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

                Dictionary<string, int> fieldMappings;
                string[] fields;

            } // FieldGrabber

            StreamReader inputFile;
            bool hasVersion;
            bool hasDone;
            string expectedVersion;
            List<string> wantedFields;
            bool hasMissingFields = false;
            bool skipHash;
            bool allowMissingColumnsInData;
			int skippedRows;
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
			public string Chromsome;
			public int Position;
			public string Hugo_Symbol;

			static KeyValuePair<string, GeneLocationInfo> ParseLine(Dictionary<string, int> fieldMappings, string[] fields)
			{
				var composite = fields[fieldMappings["CompositeREF"]];
				var geneInfo = new GeneLocationInfo();
				geneInfo.chromosome = fields[fieldMappings["Chromosome"]];
				geneInfo.minLocus = Convert.ToInt32(fields[fieldMappings["Position"]]);
				geneInfo.hugoSymbol = fields[fieldMappings["Hugo Symbol"]];

				return new KeyValuePair<string, GeneLocationInfo>(composite, geneInfo);

			}


			static public List<KeyValuePair<string, GeneLocationInfo>> ReadFile(string filename)
			{
				StreamReader inputFile;

				inputFile = CreateStreamReaderWithRetry(filename);

				var neededFields = new List<string>();
				neededFields.Add("CompositeREF");
				neededFields.Add("Chromosome");
				neededFields.Add("Position");
				neededFields.Add("Hugo Symbol");

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

		// contains annotation data from one line in a methylation file
		public class AnnotationLine
		{
			public static double betaToM(double beta) {
				// M_value calculated as shown in Du et. al. 2010
				return Math.Log(beta / (1 - beta));
			}

			// Methylation types
			public enum FeatureType {
				// Feature type not specified
				Other,
				// methylation found in 2kb region upstream from island
				N_Shore,
				// methylation found in 2kb region downstream from island
				S_Shore,
				// CpG island
				Island,
			}

			public string Composite_Element_REF;
			public double Beta_Value;
			public string Chromsome;
			public int Start;
			public int End;
			public string[] Gene_Symbol;
			public string[] Gene_Type;
			public string[] Transcript_ID;
			public int[] Position_to_TSS;
			public string CGI_Coordinate;
			public FeatureType Feature_Type;

			// this is not in the original file, but calculated from beta_value
			public readonly double M_Value;

			AnnotationLine(string Composite_Element_REF_,
			double Beta_Value_,
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
				Beta_Value = Beta_Value_;
				Chromsome = Chromosome_;
				Start = Start_;
				End = End_;
				Gene_Symbol = Gene_Symbol_;
				Gene_Type = Gene_Type_;
				Transcript_ID = Transcript_ID_;
				Position_to_TSS = Position_to_TSS_;
				CGI_Coordinate = CGI_Coordinate_;
				Feature_Type = Feature_Type_;

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

			static AnnotationLine ParseLine(Dictionary<string, int> fieldMappings, string[] fields)
			{
				// delimiter separating columns
				var delim = ';';

				Enum.TryParse(fields[fieldMappings["Feature_Type"]], out FeatureType Feature_Type);

				 return new AnnotationLine(
				fields[fieldMappings["Composite Element REF"]],
				ConvertToDoubleTreatingNullStringAsOne(fields[fieldMappings["Beta_value"]]),
				fields[fieldMappings["Chromosome"]],
				Convert.ToInt32(fields[fieldMappings["Start"]]),
				Convert.ToInt32(fields[fieldMappings["End"]]),
				splitPotentialString(fields[fieldMappings["Gene_Symbol"]], delim),
				splitPotentialString(fields[fieldMappings["Gene_Type"]], delim),
				splitPotentialString(fields[fieldMappings["Transcript_ID"]], delim),
				splitPotentialString(fields[fieldMappings["Position_to_TSS"]], delim).Select(s => Convert.ToInt32(s)).ToArray(),
				fields[fieldMappings["CGI_Coordinate"]],
				Feature_Type
				);

			}

			static public List<AnnotationLine> ReadFile(string filename, string file_id, bool fileHasVersion)
			{
				StreamReader inputFile;

				inputFile = CreateStreamReaderWithRetry(filename);

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
				var headerizedFile = new HeaderizedFile<AnnotationLine>(inputFile, fileHasVersion, hasDone, "#version gdc-1.0.0", neededFields);

				List<AnnotationLine> result;

				if (!headerizedFile.ParseFile((a, b) => ParseLine(a, b), out result))
				{
					Console.WriteLine("Error reading Annotation File " + filename);
					return null;
				}

				inputFile.Close();

				// filter out invalid results. Invalid results include records without valid M values and records without gene symbols.
				result = result.Where(c => (c.M_Value != double.PositiveInfinity) && (c.Gene_Symbol.Length > 0)).ToList();

				return result;
			} // ReadFile


		}

		public class MAFLine
        {
            public readonly string Hugo_Symbol;
            public readonly string NCBI_Build;
            public readonly string Chromosome;
            public readonly int Start_Position;
            public readonly int End_Positon;
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

            static MAFLine ParseLine(Dictionary<string, int> fieldMappings, string[] fields, string maf_file_id)
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
                    fields[fieldMappings["tumor_bam_uuid"]],
                    fields[fieldMappings["normal_bam_uuid"]],
                    ConvertToInt32TreatingNullStringAsZero(fields[fieldMappings["t_depth"]]),
                    ConvertToInt32TreatingNullStringAsZero(fields[fieldMappings["t_ref_count"]]),
                    ConvertToInt32TreatingNullStringAsZero(fields[fieldMappings["t_alt_count"]]),
                    ConvertToInt32TreatingNullStringAsZero(fields[fieldMappings["n_depth"]]),
                    ConvertToInt32TreatingNullStringAsZero(fields[fieldMappings["n_ref_count"]]),
                    ConvertToInt32TreatingNullStringAsZero(fields[fieldMappings["n_alt_count"]]),
                    maf_file_id
                    );
            }

            static public List<MAFLine> ReadFile(string filename, string file_id, bool fileHasVersion)
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
                neededFields.Add("tumor_bam_uuid");
                neededFields.Add("normal_bam_uuid");
                neededFields.Add("t_depth");
                neededFields.Add("t_ref_count");
                neededFields.Add("t_alt_count");
                neededFields.Add("n_depth");
                neededFields.Add("n_ref_count");
                neededFields.Add("n_alt_count");


                var headerizedFile = new HeaderizedFile<MAFLine>(inputFile, fileHasVersion, !fileHasVersion, "#version gdc-1.0.0", neededFields, fileHasVersion);

                List<MAFLine> result;

                if (!headerizedFile.ParseFile((a, b) => ParseLine(a, b, file_id), out result))
                {
                    Console.WriteLine("Error reading MAF File " + filename);
                    return null;
                }

                inputFile.Close();

                return result;
            } // ReadFile

            public static void WriteHeaderLine(StreamWriter outputStream)
            {
                outputStream.WriteLine("Hugo_Symbol\tNCBI_Build\tChromosome\tStart_Position\tEnd_Position\tVariant_Classification\tVariant_Type\tReference_Allele\tTumor_Seq_Allele1\tTumor_Seq_Allele2\tMatch_Norm_Seq_Allele1\tMatch_Norm_Seq_Allele2\tTumor_Sample_UUID\tMatched_Norm_Sample_UUID\ttumor_bam_uuid\tnormal_bam_uuid\tt_depth\tt_ref_count\tt_alt_count\tn_depth\tn_ref_count\tn_alt_count");
            }

            public void WriteToStream(StreamWriter output)
            {
                output.WriteLine(
                    ConvertToExcelString(Hugo_Symbol)   + "\t" +
                    NCBI_Build                          + "\t" +
                    Chromosome                          + "\t" +
                    Start_Position                      + "\t" +
                    End_Positon                         + "\t" +
                    Variant_Classification              + "\t" +
                    Variant_Type                        + "\t" +
                    Reference_Allele                    + "\t" +
                    Tumor_Seq_Allele1                   + "\t" +
                    Tumor_Seq_Allele2                   + "\t" +
                    Match_Norm_Seq_Allele1              + "\t" +
                    Match_Norm_Seq_Allele2              + "\t" +
                    Tumor_Sample_UUID                   + "\t" +
                    Matched_Norm_Sample_UUID            + "\t" +
                    tumor_bam_uuid                      + "\t" +
                    normal_bam_uuid                     + "\t" +
                    t_depth                             + "\t" +
                    t_ref_count                         + "\t" +
                    t_alt_count                         + "\t" +
                    n_depth                             + "\t" +
                    n_ref_count                         + "\t" +
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

            if (size < 1024 * 1024) return "" + ((size + 512) / 1024) + "K";

            if (size < 1024 * 1024 * 1024) return "" + ((size + 512 * 1024) / (1024 * 1024)) + "M";

            if (size < (ulong)1024 * 1024 * 1024 * 1024) return "" + ((size + 512 * 1024 * 1024) / (1024 * 1024 * 1024)) + "G";

            if (size < (ulong)1024 * 1024 * 1024 * 1024 * 1024) return "" + ((size + (ulong)512 * 1024 * 1204 * 1024) / ((ulong)1024 * 1024 * 1024 * 1024)) + "T";

            return "" + ((size + (ulong)512 * 1024 * 1204 * 1024 * 1024) / ((ulong)1024 * 1024 * 1024 * 1024 * 1024)) + "P";
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

            public bool ReadAllcountFile(ProcessBase processBase)
            {
                int currentOffset = -1;
                int currentMappedReadCount = -1;
                int whichContig = -1;

                bool sawDone = false;
                string contigName = "";

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


                        contigName = line.Substring(1).ToLower();

                        currentOffset = -1;
                        currentMappedReadCount = -1;
                        continue;
                    }

                    if (-1 == whichContig)
                    {
                        Console.WriteLine("Expected contig line after list of contigs, got " + line);
                        return false;
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

            foreach (var resultEntry in results.Where(x => !x.Value.inconsistent && x.Value.maxLocus - x.Value.minLocus > 1000000))
            {
                var result = resultEntry.Value;
                Console.WriteLine("Gene bigger than a megabase: " + result.hugoSymbol + " " + result.chromosome + ": " + result.minLocus + "-" + result.maxLocus + "(" + (result.maxLocus - result.minLocus) + ")");
            }

            refFile.Close();
			return results;
		}

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

        static void ReadMafFileAndAppendToList(string filename, string file_id, List<ASETools.MAFLine> allMAFLines, MAFLoadStatus loadStatus)
        {
            var mafLinesForThisFile = ASETools.MAFLine.ReadFile(filename, file_id, true);

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

                threads.Add(new Thread(() => ReadMafFileAndAppendToList(downloadedFiles[mafInfo.file_id].fileInfo.FullName, mafInfo.file_id, allMAFLines, loadStatus)));
            }
 

            threads.ForEach(t => t.Start());
            threads.ForEach(t => t.Join());

            if (loadStatus.nFailed > 0)
            {
                Console.WriteLine("Giving up due to failed MAF load(s).");
                return null;
            }

            byTumorSampleId = new Dictionary<string,List<MAFLine>>();

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

        //
        // Maps chromosomeName -> (offset -> MeanAndStdDev)
        //
        struct ChromosomeSize : IComparable<ChromosomeSize>
        {
            public ChromosomeSize(string name_, int size_)
            {
                name = name_;
                size = size_;
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

            public string name;
            public int size;
        }

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

		}
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
					regionalExpressionState[sizeIndex] = new RegionalExpressionState();
					exclusiveRegionalExpressionState[sizeIndex] = new RegionalExpressionState();
				}
			}

			public void AddRegionalExpression(int chromosomeOffset, double z, double mu, bool isTumor)
			{
				int distance;
				if (chromosomeOffset >= geneLocationInfo.minLocus && chromosomeOffset <= geneLocationInfo.maxLocus)
				{
					distance = 0;
				}
				else if (chromosomeOffset < geneLocationInfo.minLocus)
				{
					distance = geneLocationInfo.minLocus - chromosomeOffset;
				}
				else
				{
					distance = chromosomeOffset - geneLocationInfo.maxLocus;
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

			public RegionalExpressionState[] regionalExpressionState = new RegionalExpressionState[nRegionSizes]; // Dimension is log2(regionSize) - 1
			public RegionalExpressionState[] exclusiveRegionalExpressionState = new RegionalExpressionState[nRegionSizes];  // Expression in this region but not closer, so from log2(regionSize - 1) to log2(regionSize) - 1.  The zero element is the same as regionalExpressionState

			public ASETools.GeneLocationInfo geneLocationInfo;
			public int mutationCount = 0;
			public static StringComparer comparer;
		}


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
				outputFile.Write("Gene name\tnon-silent mutation count");

				writeRow(outputFile, allExpressions[0], printMu, minExamplesPerRegion, true, columnSuffix, isTumor);

				outputFile.WriteLine();

				for (int i = 0; i < allExpressions.Count(); i++)
				{
					outputFile.Write(ASETools.ConvertToExcelString(allExpressions[i].geneLocationInfo.hugoSymbol) + "\t" + allExpressions[i].mutationCount);

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

			public static Dictionary<string, double[]> ReadFile(string filename, bool skipFirstLine = true, bool includeMutationCount = false)
			{

				// keeps track of index of distances
				List<string> index = new List<string>();

				// for each gene, maintain a list of distance, expression tuples
				Dictionary<string, double[]> expressionMap = new Dictionary<string, double[]>();

				var reader = CreateStreamReaderWithRetry(filename);

				string line;

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

					string[] fields = line.Split('\t');

					var hugoSymbol = ConvertToNonExcelString(fields[0]);

					var mutationCount = Convert.ToInt32(fields[1]);

					// get rid of gene name and optionally mutation count
					fields = fields.Skip(includeMutationCount ? 1 : 2).ToArray();

					// the rest of the fields match the index fields

					if (fields.Length != index.Count() + (includeMutationCount ? 1 : 0))
					{
						throw new Exception("header does not match field length for file " + filename);
					}

					double[] numericFields = fields.Select(r => {
						// Set default value
						double value = double.NegativeInfinity;

						if (r != "*")
						{
							value = Convert.ToDouble(r);
						}
						return value;
					}).ToArray();

					// add gene to dictionary
					expressionMap.Add(hugoSymbol, numericFields);
				}
				return expressionMap;
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

        //
        // This breaks out the array of doubles from the entires in a RegionalSignalFile into a class that names them for the allele-specific expression case
        //
        public class AlleleSpecificSignal
        {
            public readonly string Hugo_Symbol;
            public readonly int mutationCount;
            public readonly double[] nonExclusive = new double[nRegions];
            public readonly double[] exclusive = new double[nRegions];
            public readonly double[] chromosomes = new double[nHumanNuclearChromosomes];

            public AlleleSpecificSignal(string Hugo_Symbol_, double[] numericValues)
            {
                Hugo_Symbol = Hugo_Symbol_;
                if (numericValues.Count() != 2 * nRegions + nHumanNuclearChromosomes + 1)
                {
                    throw new FormatException();
                }

                int nextValueIndex = 0;
                mutationCount = (int)numericValues[nextValueIndex];
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
                string retVal =  Hugo_Symbol + "\t" + mutationCount;
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

            public bool getValue(string contigName, int chromosomeOffset, out MeanAndStdDev meanAndStdDev)
            {
                contigName = contigName.ToLower();
                if (!expression.ContainsKey(contigName) || !expression[contigName].ContainsKey(chromosomeOffset / subchunkSize) || !expression[contigName][chromosomeOffset/subchunkSize].ContainsKey(chromosomeOffset))
                {
                    meanAndStdDev = new MeanAndStdDev(-1, -1);
                    return false;
                }

                meanAndStdDev =  expression[contigName][chromosomeOffset / subchunkSize][chromosomeOffset];
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

            public bool open(string filename) {
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

                if (!sawDone)
                {
                    Console.WriteLine("ASETools.ConsolodatedFile: truncated index file: " + filename);
                    return false;
                }

                return true;
            }

            public StreamReader getSubfile(string subfileName)
            {
                if (!subfiles.ContainsKey(subfileName))
                {
                    return null;
                }

                filestream.Seek(subfiles[subfileName].offset, SeekOrigin.Begin);

                long totalToRead = subfiles[subfileName].size;
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

            FileStream filestream;
            Dictionary<string, SubFile> subfiles = new Dictionary<string, SubFile>();
        }

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
            public readonly int locus;
            public readonly char referenceBase;
            public readonly char altBase;
        } // SelectedVariant

        public class ReadCounts
        {
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

            public double AlleleSpecificValue() // I'm not calling this "allele specific expression" because DNA doesn't have expression.
            {
                if (nMatchingAlt + nMatchingReference == 0)
                {
                    return 0;
                }

                return ((double)Math.Abs(nMatchingReference - nMatchingAlt)) / (nMatchingReference + nMatchingAlt);
            }

		}

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
			public double max_z;											// Max z For Bases With Baseline Expression
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
		}

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
		}

		public class AnnotatedVariant
        {

			public AnnotatedVariant(){}

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
                try
                {
                    var nMatchingNormalReferenceRNA = Convert.ToInt32(fields[fieldMappings["normalRNAmatchingRef"]]);
                    var nMatchingNormalAltRNA = Convert.ToInt32(fields[fieldMappings["normalRNAmatchingAlt"]]);
                    var nMatchingNormalNeitherRNA = Convert.ToInt32(fields[fieldMappings["normalRNAmatchingNeither"]]);
                    var nMatchingNormalBothRNA = Convert.ToInt32(fields[fieldMappings["normalRNAmatchingBoth"]]);

                    normalRNAReadCounts = new ReadCounts(nMatchingNormalReferenceRNA, nMatchingNormalAltRNA, nMatchingNormalNeitherRNA, nMatchingNormalBothRNA);
                }
                catch (Exception)
                {
                    // no op. No normal RNA
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

			
			public static void writeFile(string filename, List<AnnotatedVariant> annotatedVariants)
			{
				var file = CreateStreamWriterWithRetry(filename);

				if (null == file)
				{
					Console.WriteLine("AnnotatedVariant.writeFile: unable to create file " + filename);
					return;
				}

				// Write header
				file.WriteLine("Hugo_symbol\tChromosome\tPosition\tRef_allele\tAlt_allele\tVariant_type\tVariant_classification\tSomatic\t" +
					"tumorDNAmatchingRef\ttumorDNAmatchingAlt\ttumorDNAmatchingNeither\ttumorDNAmatchingBoth\t" +
					"normalDNAmatchingRef\tnormalDNAmatchingAlt\tnormalDNAmatchingNeither\tnormalDNAmatchingBoth\t" +
					"tumorRNAmatchingRef\ttumorRNAmatchingAlt\ttumorRNAmatchingNeither\ttumorRNAmatchingBoth\t" +
					"normalRNAmatchingRef\tnormalRNAmatchingAlt\tnormalRNAmatchingNeither\tnormalRNAmatchingBoth");

				foreach (var annotatedVariant in annotatedVariants)
				{
					file.WriteLine(annotatedVariant.toString());
				}
				file.WriteLine("**done**");
				file.Close();

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

            public bool IsASECandidate()    // Has at least 10 tumor DNA & RNA reads, and has variant allele frequency between .4 and .6 in tumor DNA
            {
                return tumorDNAReadCounts.nMatchingReference + tumorDNAReadCounts.nMatchingAlt >= 10 &&
                            tumorRNAReadCounts.nMatchingReference + tumorRNAReadCounts.nMatchingAlt >= 10 &&
                            tumorDNAReadCounts.nMatchingReference * 3 >= tumorDNAReadCounts.nMatchingAlt * 2 &&
                            tumorDNAReadCounts.nMatchingAlt * 3 >= tumorDNAReadCounts.nMatchingReference * 2;
            }

            public readonly string Hugo_symbol;                 // Only present for somatic mutations, otherwise ""
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

        public class SAMLine
        {
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

                while (offsetInCigarString < cigar.Count())
                {
                    switch (cigar[offsetInCigarString])
                    {
                        case 'M':
						case 'I':
						case '=':
                        case 'X':
                            offsetInCigarString++;  // Consume the M, = or X
                            int count = GetNextNumberFromString(cigar.Substring(cigarElementStart, cigarElementLength));

                            if (0 == count)
                            {
                                throw new FormatException();
                            }

                            for (int i = 0; i < count; i++)
                            {
                                mappedBases.Add(currentPos, seq[offsetInSeq]);

                                currentPos++;
                                offsetInSeq++;
                            }

							// reset cigar element information
							cigarElementLength = 0;
							cigarElementStart = offsetInCigarString;

                            break;
						case 'D':
							offsetInCigarString++;  // Consume the M, = or X
							count = GetNextNumberFromString(cigar.Substring(cigarElementStart, cigarElementLength));

							if (0 == count)
							{
								throw new FormatException();
							}
                            
							for (int i = 0; i < count; i++)
							{
								mappedBases.Add(currentPos, 'N'); // Add placeholder for DEL
								currentPos++;
							}

							// reset cigar element information
							cigarElementLength = 0;
							cigarElementStart = offsetInCigarString;

							break;
						case 'N':
						case 'H':
							offsetInCigarString++;  // Consume the M, = or X
							count = GetNextNumberFromString(cigar.Substring(cigarElementStart, cigarElementLength));

							// skip region. Reset
							currentPos += count;

							cigarElementLength = 0;
							cigarElementStart = offsetInCigarString;
							break;
						case 'S':
							offsetInCigarString++;  // Consume the M, = or X
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
            } // SAMLine

            public bool isUnmapped()
            {
                return (flag & Unmapped) == Unmapped;
            }

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

            public const int Unmapped = 0x4;

			// Dictionary of position and bases at position
            public Dictionary<int, char> mappedBases = new Dictionary<int, char>();

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
            bool isDerivedFilesDirectory = filename.Contains(configuration.derivedFilesDirectory);  // Derived files are stored in dataDirectory\..\derivedFiles
            foreach (var dataDirectory in configuration.dataDirectories)
            {
                if (filename.ToLower().StartsWith(dataDirectory.ToLower()) || isDerivedFilesDirectory && filename.ToLower().StartsWith((GoUpFilesystemLevels(dataDirectory, 1) + configuration.derivedFilesDirectory).ToLower()))
                {
                    return dataDirectory;
                }
            }

            return "none";
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

                long mappedBaseCount;

                var fields = line.Split('\t');
                if (fields.Count() != 2)
                {
                    Console.WriteLine("Format error in mapped base count file " + filename);
                    reader.Close();
                    return null;
                }

                try
                {
                    mappedBaseCount = Convert.ToInt64(fields[0]);
                } catch (FormatException)
                {
                    Console.WriteLine("Mapped base count file " + filename + " contains an unparsable count: " + fields[0]);
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
                return new MappedBaseCount(mappedBaseCount);
            }


            MappedBaseCount(long mappedBaseCount_)
            {
                mappedBaseCount = mappedBaseCount_;
            }

            public readonly long mappedBaseCount;
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

        public static bool[] BothBools = { true, false };   // There's probably something like this in the runtime, but whatever.

        public class HistogramResultLine
        {
            public string minValue;
            public int count = 0;
            public double total = 0;    // The sum of all of the values in this line
            public double pdfValue = 0;
            public double cdfValue = 0;
        } // HistogramResultLine

        public class Histogram
        {
            public Histogram() { }

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
                int runningCount = 0;
                for (int whichBucket = 0; whichBucket < nBuckets; whichBucket++)
                {
                    runningCount += result[whichBucket].count;

                    result[whichBucket].pdfValue = ((double)result[whichBucket].count) / overallCount;
                    result[whichBucket].cdfValue = ((double)runningCount) / overallCount;
                }

                return result;
            } // ComputeHistogram

            List<double> values = new List<double>();
        } // Histogram

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

            public readonly NMeandAndStdDev zeroMutationStats;
            public readonly NMeandAndStdDev oneMutationStats;
            public readonly NMeandAndStdDev moreThanOneMutationStats;

        public static List<string> getHeaders(bool mu, bool exclusive, string range)
            {
                //
                // This is pretty much a match for the code that generates the headers in the first place, which results in strange double spaces in some places.
                //
                var result = new List<string>();
                string muString = ((!mu) ? "" : " mu") + (!exclusive? "" : " exclusive");
                result.Add(range + " 1 vs. many" + muString);
                result.Add(range + " 1 vs. not 1" + muString);

                string[] mutationSets = { "0", "1", ">1" };

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
                       NMeandAndStdDev.getHeaderString(regionString + "0 mutation ") + "\t" +
                       NMeandAndStdDev.getHeaderString(regionString + "1 mutation ") + "\t" +
                       NMeandAndStdDev.getHeaderString(regionString + ">1 mutation ");
            }

            public SingleExpressionResult(int rangeIndex_, double oneVsMany_, double oneVsNotOne_, NMeandAndStdDev zeroMutationStats_, NMeandAndStdDev oneMutationStats_, NMeandAndStdDev moreThanOneMutationStats_)
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
                    rangeInBases = (1 << rangeIndex-1) * 1000;
                }

                oneVsMany = oneVsMany_;
                oneVsNotOne = oneVsNotOne_;
                zeroMutationStats = zeroMutationStats_;
                oneMutationStats = oneMutationStats_;
                moreThanOneMutationStats = moreThanOneMutationStats_;
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
                    new NMeandAndStdDev(fieldGrabber.AsIntMinusOneIfStarOrEmptyString(rangeString + " 0 mutation" + muString + " N"), fieldGrabber.AsDoubleNegativeInfinityIfStarOrEmptyString(rangeString + " 0 mutation " + muString + " mean"), fieldGrabber.AsDoubleNegativeInfinityIfStarOrEmptyString(rangeString + " 0 mutation " + muString + " stdDev")),
                    new NMeandAndStdDev(fieldGrabber.AsIntMinusOneIfStarOrEmptyString(rangeString + " 1 mutation" + muString + " N"), fieldGrabber.AsDoubleNegativeInfinityIfStarOrEmptyString(rangeString + " 1 mutation " + muString + " mean"), fieldGrabber.AsDoubleNegativeInfinityIfStarOrEmptyString(rangeString + " 1 mutation " + muString + " stdDev")),
                    new NMeandAndStdDev(fieldGrabber.AsIntMinusOneIfStarOrEmptyString(rangeString + " >1 mutation" + muString + " N"), fieldGrabber.AsDoubleNegativeInfinityIfStarOrEmptyString(rangeString + " >1 mutation " + muString + " mean"), fieldGrabber.AsDoubleNegativeInfinityIfStarOrEmptyString(rangeString + " >1 mutation " + muString + " stdDev"))
                    );
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
                    "nMore"
                };

                var wantedFields = wantedFieldsOtherThanExpression.ToList();
                for (var i = 0; i < nRegions; i++)
                {
                    foreach (bool exclusive in BothBools) {
                        wantedFields.AddRange(SingleExpressionResult.getHeaders(false, exclusive, regionIndexToString(i)));
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

            public ExpressionResultsLine(string hugo_symbol_, SingleExpressionResult[] nonExclusiveResultsByRange_, SingleExpressionResult[] exclusiveResultsByRange_, int nTumorsExcluded_, int nZero_, int nOne_, int nMore_, string rawLine_)
            {
                hugo_symbol = hugo_symbol_;
                nonExclusiveResultsByRange = nonExclusiveResultsByRange_;
                exclusiveResultsByRange = exclusiveResultsByRange_;
                nTumorsExcluded = nTumorsExcluded_;
                nZero = nZero_;
                nOne = nOne_;
                nMore = nMore_;
                rawLine = rawLine_;

                resultsByRange.Add(true, nonExclusiveResultsByRange);
                resultsByRange.Add(false, exclusiveResultsByRange);
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
                    fieldGrabber.rawLine());
            }

            public static string getHeaderString()
            {
                var result = "Hugo Symbol";
                foreach (var inclusive in BothBools)
                {
                    bool exclusive = !inclusive;
                    for (int i = 0; i < nRegions; i++)
                    {
                        result += "\t" + SingleExpressionResult.getHeaderString(exclusive, i);
                    }
                }

                return result + "\tnTumorsExcluded\tnZero\tnOne\tnMany";
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
            public readonly string rawLine;
        } // ExpressionResultsLine

        //
        // Handling power-of-two * 1Kb regions.
        //
        public const int nRegions = 21;    // 0, 1 Kb, 2Kb, 4Kb .. 256Kb, Whole Autosome
        
        public static string regionIndexToString(int regionIndex)
        {
            if (0 == regionIndex) return "0Kbp";
            if (nRegions - 1 == regionIndex) return "whole autosome";
            return "" + (1 << (regionIndex-1)) + "Kbp";
        }

        public static int regionIndexToSizeInBases(int regionIndex)
        {
            if (0 == regionIndex) return 0;
            if (nRegions - 1 == regionIndex) return int.MaxValue;
            return (1 << regionIndex) * 1000;
        }

        public class AVLTree<TValue>  where TValue : IComparable
        {
            class Node
            {
                public TValue value;

                public Node left = null, right = null, parent = null;

                public enum Balance { AVLLeft, AVLBalanced, AVLRight, AVLNew};
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

            public bool findFirstLessThan(TValue key, out TValue retVal)
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

            public bool findFirstGreaterThan(TValue key, out TValue retVal)
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

            Node root = null;
            int inserts = 0;
            int deletes = 0;

        } // AVLTree


    } // ASETools
}
