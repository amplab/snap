using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.IO;
using System.Threading.Tasks;
using ASELib;
using ExpressionNearMutations;

namespace BisulfiteMethylation
{
	class AnnotatedCpG
	{
		public AnnotatedCpG(string contig_, int locus_, ASETools.ReadCounts tumorDNAReadCounts_, ASETools.ReadCounts normalDNAReadCounts_)
		{
			contig = contig_;
			locus = locus_;

			tumorDNAReadCounts = tumorDNAReadCounts_;
			normalDNAReadCounts = normalDNAReadCounts_;
		}

		public readonly string contig;
		public readonly int locus;

		public readonly ASETools.ReadCounts tumorDNAReadCounts;
		public readonly ASETools.ReadCounts normalDNAReadCounts;

	}

	public class BisulfateCase {
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

		public BisulfateCase() {
		}

		public static Dictionary<string, BisulfateCase> loadCases(string filename) {

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
				else {
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


	class Program
	{
		static ASETools.Genome genome = new ASETools.Genome();

		static void ProcessCpGs(string filename, 
			string file_id, 
			string contig, 
			int locus, 
			string reference_allele, 
			string alt_allele, 
			string variantType,
			IEnumerable<ASETools.BedLine> cpgs)
		{
			// make sure the variant is not C. If it is modified by busulfite, we can't determine which strands are from which allele
			if (alt_allele.Contains('C') || reference_allele.Contains('C'))
			{
				return;
			}

			// calculate the ASM for each CPG (it will be the same for all CpGs overlapping this site)\
			foreach (var cpg in cpgs)
			{
				
				// count reads at CpG
				var methylationCounts = ASETools.MethylationCounts.ComputeMethylationReadCounts(filename, cpg.Chromosome,
					locus,
					reference_allele,
					alt_allele,
					variantType,
					cpg.Start_Position,
					file_id,
					genome);

				// TODO: possibly normalize by skew of normal dna in the first place?
				var nMatchingReferenceTumorDNA = methylationCounts.Item1.nMethylated + methylationCounts.Item1.nUnmethylated;
				var nMatchingVariantTumorDNA = methylationCounts.Item2.nMethylated + methylationCounts.Item2.nUnmethylated;

				// Make sure there are enough total reads, and that there is a valid mix of variant and reference DNA
				if (methylationCounts.Item1.nMethylated + methylationCounts.Item1.nUnmethylated < 3 ||
					methylationCounts.Item2.nMethylated + methylationCounts.Item2.nUnmethylated < 3 ||
					nMatchingReferenceTumorDNA * 3 < nMatchingVariantTumorDNA * 2 ||                // It's not more than 2/3 variant DNA
					nMatchingVariantTumorDNA * 3 < nMatchingReferenceTumorDNA * 2 ||
					methylationCounts.Item1.nNeither > 0 ||
					methylationCounts.Item2.nNeither > 0)
				{
					continue;
				}

				// ase1 (allelel1) = Cs with mutated read/cs with mutated + Ts with mutated
				// ase2 (allelel2) = Cs with ref read/cs with ref + Ts with ref
				// ase = ase1 / ase1 + ase2 (0 and 1 is high ASM, middle is equal amounts of methylation
				double methylationFraction_allele1 = methylationCounts.Item1.nMethylated / Convert.ToDouble(methylationCounts.Item1.nMethylated + methylationCounts.Item1.nUnmethylated);
				double methylationFraction_allele2 = methylationCounts.Item2.nMethylated / Convert.ToDouble(methylationCounts.Item2.nMethylated + methylationCounts.Item2.nUnmethylated);

				double asm_tumor = 0.5; // by default, no ASM
				if (methylationFraction_allele1 + methylationFraction_allele2 > 0) // avoid div by 0 error
				{
					asm_tumor = methylationFraction_allele2 / Convert.ToDouble(methylationFraction_allele1 + methylationFraction_allele2);
				}
				asm_tumor = Math.Abs(asm_tumor * 2.0 - 1.0);

				var asm_tumor_mu = 0;

				int chromosomeId = ASETools.ChromosomeNameToIndex(cpg.Chromosome);

				// add ASM value to autosome values
				if (ASETools.isChromosomeAutosomal(cpg.Chromosome))
				{
					wholeAutosomeRegionalExpression.AddTumorExpression(asm_tumor, asm_tumor_mu);

					foreach (var entry in allButThisChromosomeAutosomalRegionalExpressionState)
					{
						if (entry.Key != chromosomeId.ToString())
						{
							entry.Value.AddTumorExpression(asm_tumor, asm_tumor_mu);
						}
					}

					// TODO possibly add in normal expression
				}

				if (chromosomeId != -1)
				{
					if (perChromosomeRegionalExpressionState[chromosomeId] == null)
					{
						perChromosomeRegionalExpressionState[chromosomeId] = new ASETools.RegionalExpressionState();
					}

					perChromosomeRegionalExpressionState[chromosomeId].AddTumorExpression(asm_tumor, asm_tumor_mu);
					// TODO possibly add in normal expression
				}

				// add ASM value for genes
				foreach (var geneLocation in geneLocationInformation.genesByChromosome[ASETools.chromosomeNameToNonChrForm(cpg.Chromosome)])
				{
					if (!geneExpressions.ContainsKey(geneLocation.hugoSymbol))
					{
						geneExpressions.Add(geneLocation.hugoSymbol, new ASETools.GeneExpression(geneLocation));
					}

					geneExpressions[geneLocation.hugoSymbol].AddRegionalExpression(cpg.Start_Position, asm_tumor, asm_tumor_mu, true); // Recall that for allele-specifc expresion, z is really the level of allele-specific expression, not the expression z score.
																																	   // TODO possibly add in normal expression
				}


			} // foreach CpG

		}

		static ASETools.GeneLocationsByNameAndChromosome geneLocationInformation;

		// One expression state for whole autosome
		static ASETools.RegionalExpressionState wholeAutosomeRegionalExpression = new ASETools.RegionalExpressionState();
		// Expression state for each chromosome, which will exclude the chromosome that the gene resides on
		static new Dictionary<string, ASETools.RegionalExpressionState> allButThisChromosomeAutosomalRegionalExpressionState = new Dictionary<string, ASETools.RegionalExpressionState>();   // "This chromosome" is the dictionary key
																																 // Expression state for each chromosome, which will include the chromsome that the gene resides on
		static ASETools.RegionalExpressionState[] perChromosomeRegionalExpressionState = new ASETools.RegionalExpressionState[ASETools.nHumanNuclearChromosomes];
		// per gene ASM
		static Dictionary<string, ASETools.GeneExpression> geneExpressions = new Dictionary<string, ASETools.GeneExpression>();


		static void Main(string[] args)
		{
			// Add cases with bisulfate data. We will use the bam and bisulfate files
			var bisulfiteCases = BisulfateCase.loadCases(ASETools.ASEConfirguation.bisulfiteCasesFilePathname);

			var configuration = ASETools.ASEConfirguation.loadFromFile(args);

			genome.load(configuration.indexDirectoryHg19);

			// load in gene information
			geneLocationInformation = new ASETools.GeneLocationsByNameAndChromosome(ASETools.readKnownGeneFile(ASETools.ASEConfirguation.hg19GeneLocationInformation));

			var cases = ASETools.Case.LoadCases(configuration.casesFilePathname);
			// filter out only bisulfite cases
			cases = cases.Where(case_ => bisulfiteCases.ContainsKey(case_.Value.case_id)).ToDictionary(x => x.Key, x => x.Value);

			foreach (var bCase in bisulfiteCases)
			{
				Console.WriteLine("Processing case " + bCase.Key);

				var bisulfate = bisulfiteCases[bCase.Key];

				var case_ = cases[bCase.Key];

				// reset allButThisChromosomeAutosomalRegionalExpressionState
				foreach (var geneEntry in geneLocationInformation.genesByName)
				{
					var chromosome = geneEntry.Value.chromosome;
					if (!allButThisChromosomeAutosomalRegionalExpressionState.ContainsKey(chromosome))
					{
						allButThisChromosomeAutosomalRegionalExpressionState.Add(chromosome, new ASETools.RegionalExpressionState());
					}
				}
				// reset perChromosomeRegionalExpressionState
				for (int whichChromosome = 0; whichChromosome < ASETools.nHumanNuclearChromosomes; whichChromosome++)
				{
					perChromosomeRegionalExpressionState[whichChromosome] = new ASETools.RegionalExpressionState();
				}

				// read in annotated variants
				if (bisulfate.realigned_selected_variants_filename == "" || bisulfate.realigned_maf_filename == "")
				{
					continue;
				}

				// variants and maflines realigned to hg19 genome
				var selectedVariants = ASETools.SelectedVariant.LoadFromFile(bisulfate.realigned_selected_variants_filename);
				var mafLines = ASETools.MAFLine.ReadFile(bisulfate.realigned_maf_filename, case_.maf_file_id, false);

				// read in cpgs
				List<ASETools.BedLine> cpgs = ASETools.BedLine.ReadFile(bisulfate.bed_tumor_filename);

				var indexedReadsFilename = ASETools.ASEConfirguation.bisulfiteDirectory + bisulfate.case_id + @"\" + bisulfate.bam_tumor_file_id;

				var variantsCount = 0;
				foreach (var selectedVariant in selectedVariants)
				{
					if (variantsCount % 500 == 0)
					{
						Console.WriteLine("Completed " + variantsCount + " variants.");
					}
					variantsCount++;

					// make sure the variant is not C
					if (selectedVariant.altBase == 'C' || selectedVariant.referenceBase == 'C')
					{
						continue;
					}

					// take all cpgs from bed file within 50 bp of this variant, but we also don't want CpG islands overlapping this variant
					var overlappingCpgs = cpgs.Where(r => r.Chromosome == selectedVariant.contig && Math.Abs(selectedVariant.locus - r.Start_Position) <= 50 && Math.Abs(selectedVariant.locus - r.Start_Position) > 5);

					ProcessCpGs(indexedReadsFilename, bisulfate.bam_tumor_filename + selectedVariant.getExtractedReadsExtension(), selectedVariant.contig, selectedVariant.locus,
						selectedVariant.referenceBase.ToString(),
						selectedVariant.altBase.ToString(),
						"SNP", // all of them are SNPs 
						overlappingCpgs);

				} // foreach variant

				var mafLineCount = 0;
				foreach (var mafLine in mafLines)
				{
					// If we do not have gene location information, ignore it
					if (!geneLocationInformation.genesByName.ContainsKey(mafLine.Hugo_Symbol))
					{
						Console.WriteLine("Could not find information for " + mafLine.Hugo_Symbol + " in list of known genes. Skipping.");
						continue;
					}

					// For debugging
					if (mafLineCount % 100 == 0)
					{
						Console.WriteLine("Completed " + mafLineCount + " mafLines.");
					}
					mafLineCount++;

					// make sure the variant is not C
					if (mafLine.Match_Norm_Seq_Allele1.Contains('C') || mafLine.Match_Norm_Seq_Allele2.Contains('C'))
					{
						continue;
					}

					// if Gene Symbol not yet in dictionary, add it
					if (!geneExpressions.ContainsKey(mafLine.Hugo_Symbol))
					{
						geneExpressions.Add(mafLine.Hugo_Symbol, new ASETools.GeneExpression(geneLocationInformation.genesByName[mafLine.Hugo_Symbol]));
					}

					// Increment the muation count by 1
					geneExpressions[mafLine.Hugo_Symbol].mutationCount++;

					// Take all cpgs from bed file within 50 bp of this variant, but we also don't want CpG islands overlapping this variant
					var overlappingCpgs = cpgs.Where(r => r.Chromosome == mafLine.Chromosome && Math.Abs(mafLine.Start_Position - r.Start_Position) <= 50 && Math.Abs(mafLine.Start_Position - r.Start_Position) > 5);

					ProcessCpGs(indexedReadsFilename, bisulfate.bam_tumor_filename + mafLine.getExtractedReadsExtension(), mafLine.Chromosome, mafLine.Start_Position, mafLine.Match_Norm_Seq_Allele1,
						mafLine.Tumor_Seq_Allele2, mafLine.Variant_Type, overlappingCpgs);


				} // foreach MAFLine

				//
				// Write the output file for ASM.
				//
				string directory = ASETools.GetDirectoryPathFromFullyQualifiedFilename(bisulfate.realigned_selected_variants_filename);
				string analysisId = ASETools.GetAnalysisIdFromPathname(bisulfate.realigned_selected_variants_filename);
				if ("" == directory || "" == analysisId)
				{
					Console.WriteLine("Couldn't parse input pathname, which is supposed to be absolute and include an analysis ID: " + bisulfate.realigned_selected_variants_filename);
					continue;
				}

				var outputFilename = directory + analysisId + (ASETools.bisulfiteAlleleSpecificMethylationExtension);
				Console.WriteLine("Saving ASM to file " + outputFilename);

				var outputFile = ASETools.CreateStreamWriterWithRetry(outputFilename);



				var allExpressions = new List<ASETools.GeneExpression>();
				foreach (var expressionEntry in geneExpressions)
				{
					allExpressions.Add(expressionEntry.Value);
				}

				allExpressions.Sort(ASETools.GeneExpression.CompareByGeneName);

				// TODO copy approach in ENM
				var fileHeader = "BisulfiteMethylation " + case_.case_id;
				var output = new ASETools.RegionalSignalFile(fileHeader,
					allExpressions,
					wholeAutosomeRegionalExpression, allButThisChromosomeAutosomalRegionalExpressionState,
					perChromosomeRegionalExpressionState);


				outputFile.Close();

				//////////////////
				// clear expression for next case
				//////////////////

				// One expression state for whole autosome
				wholeAutosomeRegionalExpression = new ASETools.RegionalExpressionState();
				// Expression state for each chromosome, which will exclude the chromosome that the gene resides on
				allButThisChromosomeAutosomalRegionalExpressionState = new Dictionary<string, ASETools.RegionalExpressionState>();   // "This chromosome" is the dictionary key
																																																	 // Expression state for each chromosome, which will include the chromsome that the gene resides on
				perChromosomeRegionalExpressionState = new ASETools.RegionalExpressionState[ASETools.nHumanNuclearChromosomes];
				// per gene ASM
				geneExpressions = new Dictionary<string, ASETools.GeneExpression>();

			} // foreach case

		} 
	}
}
