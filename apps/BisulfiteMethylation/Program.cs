using System;
using System.Collections.Generic;
using System.Linq;
using ASELib;

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


	class Program
	{
		static ASETools.Genome genome = new ASETools.Genome();

		static int nTotal = 0;
		static int nProcessed = 0;

		static Dictionary<string, List<Tuple<string, double>>> rawValuesTumor = new Dictionary<string, List<Tuple<string, double>>>();
		static Dictionary<string, List<Tuple<string, double>>> rawValuesNormal = new Dictionary<string, List<Tuple<string, double>>>();

		static void ProcessCpGs(ASETools.Case case_, bool forTumor, 
			string filename, 
			string file_id, 
			string contig, 
			int locus, 
			string reference_allele, 
			string alt_allele, 
			string variantType,
			IEnumerable<ASETools.BedLine> cpgs,
			ExpressionNearMutations.Program.RegionalSignal regionalSignal)
		{
			// make sure the variant is not C. If it is modified by busulfite, we can't determine which strands are from which allele
			if (alt_allele.Contains('C') || reference_allele.Contains('C'))
			{
				return;
			}

			// calculate the ASM for each CPG (it will be the same for all CpGs overlapping this site)
			foreach (var cpg in cpgs)
			{
				nTotal += 1;
				
				// count reads at CpG
				var methylationCounts = ASETools.MethylationCounts.ComputeMethylationReadCounts(filename, cpg.Chromosome,
					locus,
					reference_allele,
					alt_allele,
					variantType,
					cpg.Start_Position,
					file_id,
					genome);

				int chromosomeId = ASETools.ChromosomeNameToIndex(cpg.Chromosome);

				// TODO: possibly normalize by skew of normal dna in the first place?
				var nMatchingReferenceTumorDNA = methylationCounts.Item1.nMethylated + methylationCounts.Item1.nUnmethylated;
				var nMatchingVariantTumorDNA = methylationCounts.Item2.nMethylated + methylationCounts.Item2.nUnmethylated;

				// Make sure there are enough total reads, and that there is a valid mix of variant and reference DNA
				if (methylationCounts.Item1.nMethylated + methylationCounts.Item1.nUnmethylated < 4 ||
					methylationCounts.Item2.nMethylated + methylationCounts.Item2.nUnmethylated < 4 ||
					nMatchingReferenceTumorDNA * 3 < nMatchingVariantTumorDNA * 2 ||                // It's not more than 2/3 variant DNA
					nMatchingVariantTumorDNA * 3 < nMatchingReferenceTumorDNA * 2 ||
					methylationCounts.Item1.nNeither > 0 ||
					methylationCounts.Item2.nNeither > 0 && chromosomeId == ASETools.ChromosomeNameToIndex("chrX"))
				{
					continue;
				}

				nProcessed += 1;

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

				var interval = new GenomeBuild.Interval(cpg.Chromosome, cpg.Start_Position, cpg.End_Position, cpg.strand).ToString();
				// add to tumor raw values
				if (forTumor)
				{
					List<Tuple<string, double>> values;
					if (!rawValuesTumor.TryGetValue(interval, out values))
					{
						rawValuesTumor.Add(interval, new List<Tuple<string, double>>());
					}
					rawValuesTumor[interval].Add(new Tuple<string, double>(case_.case_id, asm_tumor));
				} else
				{
					List<Tuple<string, double>> values;
					if (!rawValuesNormal.TryGetValue(interval, out values))
					{
						rawValuesNormal.Add(interval, new List<Tuple<string, double>>());
					}
					rawValuesNormal[interval].Add(new Tuple<string, double>(case_.case_id, asm_tumor));
				}

				regionalSignal.Add(cpg.Chromosome, cpg.Start_Position, asm_tumor, asm_tumor_mu, false);

			} // foreach CpG


		}

		static ASETools.GeneLocationsByNameAndChromosome geneLocationInformation;

		static void ProcessFile(ASETools.BisulfateCase bisulfite, ASETools.Case case_, bool forTumor)
		{

			// read in annotated variants
			if (bisulfite.realigned_selected_variants_filename == "" || bisulfite.realigned_maf_filename == "")
			{
				return;
			}

			// variants and maflines realigned to hg19 genome
			var dir = ASETools.GetDirectoryFromPathname(bisulfite.realigned_maf_filename);
			var selectedVariants = ASETools.SelectedVariant.LoadFromFile(bisulfite.realigned_selected_variants_filename);
			var mafLines = ASETools.MAFLine.ReadFile(bisulfite.realigned_maf_filename, case_.maf_file_id, false);

			//// save variants and mafs as bed to read in IGV
			//var bedVariantsFilename = dir + @"\" + bisulfite.case_id + "_bedVariants.bed";
			//var bedVariants = selectedVariants.Select(r => new ASETools.BedLine(r.contig, r.locus - 1, r.locus, r.referenceBase + "_" + r.altBase, 1000, '.', r.locus, r.locus + 1));
			//ASETools.BedLine.writeFile(bedVariantsFilename, bedVariants.ToList());

			//var bedMAFFilename = dir + @"\" + bisulfite.case_id + "_bedMAFs.bed";
			//var bedMafs = mafLines.Select(r => new ASETools.BedLine(r.Chromosome, r.Start_Position - 1, r.Start_Position, r.Hugo_Symbol + "_" + r.Reference_Allele + "_" + r.Tumor_Seq_Allele2, 1000, '.', r.Start_Position, r.End_Positon));
			//ASETools.BedLine.writeFile(bedMAFFilename, bedMafs.ToList());

			// read in cpgs
			var bedFilename = bisulfite.bed_tumor_filename;
			var bamFilename = forTumor ? bisulfite.bam_tumor_filename : bisulfite.bam_normal_filename;

			if (bedFilename == "" || bamFilename == "" || bedFilename == null || bamFilename == null)
			{
				Console.WriteLine("No bisulfite data for " + bisulfite.case_id + ". Skipping...");
			}

			var geneExpressions = new Dictionary<string, ASETools.GeneExpression>();

			foreach (var mafLine in mafLines)
			{
				if (mafLine.Variant_Classification == "Silent")
				{
					continue;
				}

				if (!geneLocationInformation.genesByName.ContainsKey(mafLine.Hugo_Symbol))
				{
					//
					// Probably an inconsistent gene.  Skip it.
					//
					continue;
				}

				// if Gene Symbol not yet in dictionary, add it
				if (!geneExpressions.ContainsKey(mafLine.Hugo_Symbol))
				{
					geneExpressions.Add(mafLine.Hugo_Symbol, new ASETools.GeneExpression(geneLocationInformation.genesByName[mafLine.Hugo_Symbol]));
				}

				// Increment the muation count by 1
				geneExpressions[mafLine.Hugo_Symbol].mutationCount++;
			}

			// Stores methylation numbers
			var regionalSignal = new ExpressionNearMutations.Program.RegionalSignal(geneExpressions, geneLocationInformation, false);

			List<ASETools.BedLine> cpgs = ASETools.BedLine.ReadFile(bedFilename);

			string indexedReadsFilename;
			if (forTumor)
			{
				indexedReadsFilename = ASETools.ASEConfirguation.bisulfiteDirectory + bisulfite.case_id + @"\" + bisulfite.bam_tumor_file_id;
			}
			else
			{
				indexedReadsFilename = ASETools.ASEConfirguation.bisulfiteDirectory + bisulfite.case_id + @"\" + bisulfite.bam_normal_file_id;
			}

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


					ProcessCpGs(case_, forTumor, indexedReadsFilename, bamFilename + selectedVariant.getExtractedReadsExtension(), selectedVariant.contig, selectedVariant.locus,
					selectedVariant.referenceBase.ToString(),
					selectedVariant.altBase.ToString(),
					"SNP", // all of them are SNPs 
					overlappingCpgs,
					regionalSignal);

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

				ProcessCpGs(case_, forTumor, indexedReadsFilename, bamFilename + mafLine.getExtractedReadsExtension(), mafLine.Chromosome, mafLine.Start_Position, mafLine.Match_Norm_Seq_Allele1,
					mafLine.Tumor_Seq_Allele2, mafLine.Variant_Type, overlappingCpgs, regionalSignal);

			} // foreach MAFLine

			Console.WriteLine("Processed " + nProcessed + " CpGs out of " + nTotal);
			nTotal = 0;
			nProcessed = 0;

			//
			// Write the output file for ASM.
			//
			string directory = ASETools.GetDirectoryPathFromFullyQualifiedFilename(bisulfite.realigned_selected_variants_filename);
			string analysisId = ASETools.GetAnalysisIdFromPathname(bisulfite.realigned_selected_variants_filename);
			if ("" == directory || "" == analysisId)
			{
				Console.WriteLine("Couldn't parse input pathname, which is supposed to be absolute and include an analysis ID: " + bisulfite.realigned_selected_variants_filename);
				return;
			}

			var extention = forTumor ? ASETools.tumorBisulfiteAlleleSpecificMethylationExtension : ASETools.normalBisulfiteAlleleSpecificMethylationExtension;
			var outputFilename = directory + analysisId + extention;
			Console.WriteLine("Saving ASM to file " + outputFilename);

			var allExpressions = new List<ASETools.GeneExpression>();
			foreach (var expressionEntry in geneExpressions)
			{
				allExpressions.Add(expressionEntry.Value);
			}

			allExpressions.Sort(ASETools.GeneExpression.CompareByGeneName);

			var fileHeader = "BisulfiteMethylation " + case_.case_id;
			var output = new ASETools.RegionalSignalFile(fileHeader,
				allExpressions,
				regionalSignal.wholeAutosomeRegionalExpression,
				regionalSignal.allButThisChromosomeAutosomalRegionalExpressionState,
				regionalSignal.perChromosomeRegionalExpressionState);

			var columnSuffix = forTumor ? "tumor ASM" : "normal ASM";
			var printMu = false;
			output.WriteFile(outputFilename, printMu, 1, columnSuffix, true);

		}

		static void writeRawValues(string filename, bool isTumor)
		{
			var writer = ASETools.CreateStreamWriterWithRetry(filename);

			var values = isTumor ? rawValuesTumor : rawValuesNormal;

			// write header case ids
			var caseIds = values.SelectMany(r => r.Value.Select(k => k.Item1)).Distinct();
			writer.WriteLine("\t" + caseIds);

			foreach (var item in values)
			{
				// write interval
				writer.Write(item.Key);

				var rowIds = item.Value.Select(r => r.Item1).ToList();

				foreach (var id in caseIds)
				{
					int idx;
					if ((idx = rowIds.IndexOf(id)) >= 0)
					{
						// write ASM values for this case
						writer.Write("\t" + item.Value[idx]);
					} else {
						writer.Write("\tNA");
					}
				}

				writer.WriteLine();


			}

			writer.Close();
		}

		static int Main(string[] args)
		{
			// Add cases with bisulfite data. We will use the bam and bisulfite files
			var configuration = ASETools.ASEConfirguation.loadFromFile(args);

			var bisulfiteCases = ASETools.BisulfateCase.loadCases(ASETools.ASEConfirguation.bisulfiteCasesFilePathname);
			var cases = ASETools.Case.LoadCases(configuration.casesFilePathname);

			if (configuration.commandLineArgs.Count()  != 0)
			{
				Console.WriteLine("usage: BisulfiteAnalysis");
				return 1;
			}

			// load hg19 genome
			genome.load(configuration.indexDirectoryHg19);

			// load in gene information
			geneLocationInformation = new ASETools.GeneLocationsByNameAndChromosome(ASETools.readKnownGeneFile(ASETools.ASEConfirguation.hg19GeneLocationInformation));

			// filter out only bisulfite cases

			cases = cases.Where(case_ => bisulfiteCases.ContainsKey(case_.Value.case_id)).ToDictionary(x => x.Key, x => x.Value);

			foreach (var bCase in bisulfiteCases)
			{

				Console.WriteLine("Processing case " + bCase.Key);
				var bisulfite = bisulfiteCases[bCase.Key];

				var case_ = cases[bCase.Key];

				// Process Tumor methylation
				ProcessFile(bisulfite, case_, true);

				// Process Normal methylation
				ProcessFile(bisulfite, case_, false);

			} // foreach case

			// write out raw values
			var tumorOutput = ASETools.ASEConfirguation.bisulfiteDirectory + "tumorRawOutput.txt";
			var normalOutput = ASETools.ASEConfirguation.bisulfiteDirectory + "normalRawOutput.txt";
			
			writeRawValues(tumorOutput, true);
			writeRawValues(normalOutput, false);

			return 0;
		} 
	}
}
