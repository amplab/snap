using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using System.IO;
using ASELib;
using System.Diagnostics;
using System.Threading;

namespace ExpressionNearMutations
{
    class Program
    {

		static ASETools.GeneLocationsByNameAndChromosome geneLocationInformation;

		static void writeColumnNames(StreamWriter outputFile, bool forAlleleSpecificExpression, ASETools.Case case_, bool isTumor)
		{
			string columnSuffix = forAlleleSpecificExpression ? "(Tumor ase)" : "(Tumor z)";

			if (!isTumor)
			{
				columnSuffix = forAlleleSpecificExpression ? "(Normal ase)" : "(Normal z)";

			}

			string columnSuffix_mu = isTumor ? "(Tumor mu)" : "(Normal mu)";

			for (int sizeIndex = 0; sizeIndex < ASETools.GeneExpression.nRegionSizes; sizeIndex++)
			{
				outputFile.Write("\t" + ASETools.GeneExpression.regionSizeByRegionSizeIndex[sizeIndex] + columnSuffix);

			}

			outputFile.Write("\tWhole Autosome " + columnSuffix);

			if (!forAlleleSpecificExpression)
			{
				for (int sizeIndex = 0; sizeIndex < ASETools.GeneExpression.nRegionSizes; sizeIndex++)
				{
					outputFile.Write("\t" + ASETools.GeneExpression.regionSizeByRegionSizeIndex[sizeIndex] + columnSuffix_mu);
				}
				outputFile.Write("\tWhole Autosome " + columnSuffix_mu);
			}

			for (int sizeIndex = 0; sizeIndex < ASETools.GeneExpression.nRegionSizes; sizeIndex++)
			{
				outputFile.Write("\t" + ASETools.GeneExpression.regionSizeByRegionSizeIndex[sizeIndex] + " exclusive " + columnSuffix);

			}

			outputFile.Write("\tWhole Autosome exclusive " + columnSuffix);

			if (!forAlleleSpecificExpression)
			{
				for (int sizeIndex = 0; sizeIndex < ASETools.GeneExpression.nRegionSizes; sizeIndex++)
				{
					outputFile.Write("\t" + ASETools.GeneExpression.regionSizeByRegionSizeIndex[sizeIndex] + " exclusive " + columnSuffix_mu);

				}
				outputFile.Write("\tWhole Autosome exclusive " + columnSuffix_mu);

			}

			for (int whichChromosome = 0; whichChromosome < ASETools.nHumanNuclearChromosomes; whichChromosome++)
			{
				outputFile.Write("\t" + ASETools.ChromosomeIndexToName(whichChromosome, true) + columnSuffix);
			}

			if (!forAlleleSpecificExpression)
			{
				for (int whichChromosome = 0; whichChromosome < ASETools.nHumanNuclearChromosomes; whichChromosome++)
				{
					outputFile.Write("\t" + ASETools.ChromosomeIndexToName(whichChromosome, true) + columnSuffix_mu);

				}
			}
		}

		static void writeRow(StreamWriter outputFile, ASETools.GeneExpression allExpression,
			ASETools.RegionalExpressionState wholeAutosomeRegionalExpression,
			Dictionary<string, ASETools.RegionalExpressionState> allButThisChromosomeAutosomalRegionalExpressionState,
			ASETools.RegionalExpressionState[] perChromosomeRegionalExpressionState,
			bool forAlleleSpecificExpression, int minExamplesPerRegion)
		{
			// write tumor values
			for (int sizeIndex = 0; sizeIndex < ASETools.GeneExpression.nRegionSizes; sizeIndex++)
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

			if (wholeAutosomeRegionalExpression.nRegionsIncludedTumor >= minExamplesPerRegion)
			{
				outputFile.Write("\t" + wholeAutosomeRegionalExpression.totalTumorExpression / wholeAutosomeRegionalExpression.nRegionsIncludedTumor);
			}
			else
			{
				outputFile.Write("\t*");
			}

			if (!forAlleleSpecificExpression)
			{
				for (int sizeIndex = 0; sizeIndex < ASETools.GeneExpression.nRegionSizes; sizeIndex++)
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

				if (wholeAutosomeRegionalExpression.nRegionsIncludedTumor >= minExamplesPerRegion)
				{
					outputFile.Write("\t" + wholeAutosomeRegionalExpression.totalMeanTumorExpression / wholeAutosomeRegionalExpression.nRegionsIncludedTumor);
				}
				else
				{
					outputFile.Write("\t*");
				}
			}

			for (int sizeIndex = 0; sizeIndex < ASETools.GeneExpression.nRegionSizes; sizeIndex++)
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

			if (allButThisChromosomeAutosomalRegionalExpressionState[allExpression.geneLocationInfo.chromosome].nRegionsIncludedTumor >= minExamplesPerRegion)
			{
				outputFile.Write("\t" + allButThisChromosomeAutosomalRegionalExpressionState[allExpression.geneLocationInfo.chromosome].totalTumorExpression / allButThisChromosomeAutosomalRegionalExpressionState[allExpression.geneLocationInfo.chromosome].nRegionsIncludedTumor);
			}
			else
			{
				outputFile.Write("\t*");
			}

			if (!forAlleleSpecificExpression)
			{
				for (int sizeIndex = 0; sizeIndex < ASETools.GeneExpression.nRegionSizes; sizeIndex++)
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

				if (allButThisChromosomeAutosomalRegionalExpressionState[allExpression.geneLocationInfo.chromosome].nRegionsIncludedTumor >= minExamplesPerRegion)
				{
					outputFile.Write("\t" + allButThisChromosomeAutosomalRegionalExpressionState[allExpression.geneLocationInfo.chromosome].totalMeanTumorExpression / allButThisChromosomeAutosomalRegionalExpressionState[allExpression.geneLocationInfo.chromosome].nRegionsIncludedTumor);
				}
				else
				{
					outputFile.Write("\t*");
				}
			}

			for (int whichChromosome = 0; whichChromosome < ASETools.nHumanNuclearChromosomes; whichChromosome++)
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

			if (!forAlleleSpecificExpression)
			{
				for (int whichChromosome = 0; whichChromosome < ASETools.nHumanNuclearChromosomes; whichChromosome++)
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
			}


			// write normal values
			for (int sizeIndex = 0; sizeIndex < ASETools.GeneExpression.nRegionSizes; sizeIndex++)
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

			if (wholeAutosomeRegionalExpression.nRegionsIncludedNormal >= minExamplesPerRegion)
			{
				outputFile.Write("\t" + wholeAutosomeRegionalExpression.totalNormalExpression / wholeAutosomeRegionalExpression.nRegionsIncludedNormal);
			}
			else
			{
				outputFile.Write("\t*");
			}

			if (!forAlleleSpecificExpression)
			{
				for (int sizeIndex = 0; sizeIndex < ASETools.GeneExpression.nRegionSizes; sizeIndex++)
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


				if (wholeAutosomeRegionalExpression.nRegionsIncludedNormal >= minExamplesPerRegion)
				{
					outputFile.Write("\t" + wholeAutosomeRegionalExpression.totalMeanNormalExpression / wholeAutosomeRegionalExpression.nRegionsIncludedNormal);
				}
				else
				{
					outputFile.Write("\t*");
				}
			}

			for (int sizeIndex = 0; sizeIndex < ASETools.GeneExpression.nRegionSizes; sizeIndex++)
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

			if (allButThisChromosomeAutosomalRegionalExpressionState[allExpression.geneLocationInfo.chromosome].nRegionsIncludedNormal >= minExamplesPerRegion)
			{
				outputFile.Write("\t" + allButThisChromosomeAutosomalRegionalExpressionState[allExpression.geneLocationInfo.chromosome].totalNormalExpression / allButThisChromosomeAutosomalRegionalExpressionState[allExpression.geneLocationInfo.chromosome].nRegionsIncludedNormal);
			}
			else
			{
				outputFile.Write("\t*");
			}


			if (!forAlleleSpecificExpression)
			{
				for (int sizeIndex = 0; sizeIndex < ASETools.GeneExpression.nRegionSizes; sizeIndex++)
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


				if (allButThisChromosomeAutosomalRegionalExpressionState[allExpression.geneLocationInfo.chromosome].nRegionsIncludedNormal >= minExamplesPerRegion)
				{
					outputFile.Write("\t" + allButThisChromosomeAutosomalRegionalExpressionState[allExpression.geneLocationInfo.chromosome].totalMeanNormalExpression / allButThisChromosomeAutosomalRegionalExpressionState[allExpression.geneLocationInfo.chromosome].nRegionsIncludedNormal);
				}
				else
				{
					outputFile.Write("\t*");
				}
			}

			for (int whichChromosome = 0; whichChromosome < ASETools.nHumanNuclearChromosomes; whichChromosome++)
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

			if (!forAlleleSpecificExpression)
			{
				for (int whichChromosome = 0; whichChromosome < ASETools.nHumanNuclearChromosomes; whichChromosome++)
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

		class Temp
		{ // TODO rename
			public string chromosome;
			public int offset;

			public double z_tumor;
			public double z_normal;

			public double mu_tumor;
			public double mu_normal;

			public bool hasNormal;

			public Temp(string chromosome_, int offset_, double z_tumor_, double mu_tumor_, bool hasNormal_, double z_normal_ = 0, double mu_normal_ = 0)
			{
				this.chromosome = chromosome_;
				this.offset = offset_;
				this.z_tumor = z_tumor_;
				this.mu_tumor = mu_tumor_;

				this.hasNormal = hasNormal_;

				this.z_normal = z_normal_;
				this.mu_normal = mu_normal_;
			}
		}

		static void ProcessCases(List<ASETools.Case> casesToProcess, bool forAlleleSpecificExpression, int minExamplesPerRegion)
		{
            var timer = new Stopwatch();

			while (true)
			{
				ASETools.Case case_ = null;

				lock (casesToProcess)
				{
					if (casesToProcess.Count() == 0)
					{
						return;
					}

					case_ = casesToProcess[0];
					casesToProcess.RemoveAt(0);
				}

				timer.Reset();
				timer.Start();

				var inputFilename = forAlleleSpecificExpression ? case_.annotated_selected_variants_filename : case_.regional_expression_filename;

				if (inputFilename == "")
				{
					Console.WriteLine("Case " + case_.case_id + " doesn't have an input file yet.");
					continue;
				}

				// Load MAF file for this case
				var mafLines = ASETools.MAFLine.ReadFile(case_.extracted_maf_lines_filename, case_.maf_file_id, false);

				if (null == mafLines)
				{
					Console.WriteLine("Case " + case_.case_id + " failed to load extracted MAF lines.  Ignoring.");
					continue;
				}

				// dictionary of gene symbols
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

				List<ASETools.AnnotatedVariant> annotatedVariants = null;
				List<int> regionalExpressionLines = null;

				if (forAlleleSpecificExpression)
				{
					annotatedVariants = ASETools.AnnotatedVariant.readFile(case_.annotated_selected_variants_filename);
					if (annotatedVariants == null)
					{
						Console.WriteLine("Failed to read " + case_.annotated_selected_variants_filename + ".  Giving up on case " + case_.case_id);
						continue;
					}
				}
				else
				{
					// AM TODO
					//
					// Fill in something here when we have code for overall expression.
					//
				}

				// Variables storing expression state other than same-chromosome regional.

				// One expression state for whole autosome
				var wholeAutosomeRegionalExpression = new ASETools.RegionalExpressionState();
				// Expression state for each chromosome, which will exclude the chromosome that the gene resides on
				var allButThisChromosomeAutosomalRegionalExpressionState = new Dictionary<string, ASETools.RegionalExpressionState>();   // "This chromosome" is the dictionary key
																																// Expression state for each chromosome, which will include the chromsome that the gene resides on
				var perChromosomeRegionalExpressionState = new ASETools.RegionalExpressionState[ASETools.nHumanNuclearChromosomes];

				for (int whichChromosome = 0; whichChromosome < ASETools.nHumanNuclearChromosomes; whichChromosome++)
				{
					perChromosomeRegionalExpressionState[whichChromosome] = new ASETools.RegionalExpressionState();
				}

				foreach (var geneEntry in geneLocationInformation.genesByName)
				{
					var chromosome = geneEntry.Value.chromosome;
					if (!allButThisChromosomeAutosomalRegionalExpressionState.ContainsKey(chromosome))
					{
						allButThisChromosomeAutosomalRegionalExpressionState.Add(chromosome, new ASETools.RegionalExpressionState());
					}
				}

				var values = new List<Temp>();

				if (forAlleleSpecificExpression)
				{
					foreach (var annotatedVariant in annotatedVariants)
					{
						double alleleSpecificExpression;
						double alleleSpecificExpressionNormal = 0;
						bool hasNormal = false;

						if (geneLocationInformation.genesByChromosome.ContainsKey(annotatedVariant.contig) &&
							annotatedVariant.tumorDNAReadCounts.nMatchingReference + annotatedVariant.tumorDNAReadCounts.nMatchingAlt >= 10 &&
							annotatedVariant.tumorRNAReadCounts.nMatchingReference + annotatedVariant.tumorRNAReadCounts.nMatchingAlt >= 10 &&
							annotatedVariant.tumorDNAReadCounts.nMatchingReference * 3 >= annotatedVariant.tumorDNAReadCounts.nMatchingAlt * 2 &&
							annotatedVariant.tumorDNAReadCounts.nMatchingAlt * 3 >= annotatedVariant.tumorDNAReadCounts.nMatchingReference * 2)
						{
							double rnaFractionTumor = (double)annotatedVariant.tumorRNAReadCounts.nMatchingAlt / (annotatedVariant.tumorRNAReadCounts.nMatchingReference + annotatedVariant.tumorRNAReadCounts.nMatchingAlt);

							//
							// Now convert to the amount of allele-specific expression.  50% is no ASE, while 0 or 100% is 100% ASE.
							//
							alleleSpecificExpression = Math.Abs(rnaFractionTumor * 2.0 - 1.0);

							// If we have the normal DNA and RNA for this sample, compute the normal ASE
							if (annotatedVariant.normalRNAReadCounts != null && annotatedVariant.normalDNAReadCounts.nMatchingReference + annotatedVariant.normalDNAReadCounts.nMatchingAlt >= 10 &&   // We have at least 10 DNA reads
								annotatedVariant.normalRNAReadCounts.nMatchingReference + annotatedVariant.normalRNAReadCounts.nMatchingAlt >= 10 &&            // We have at least 10 RNA reads
								annotatedVariant.normalDNAReadCounts.nMatchingReference * 3 >= annotatedVariant.normalDNAReadCounts.nMatchingAlt * 2 &&         // It's not more than 2/3 variant DNA
								annotatedVariant.normalDNAReadCounts.nMatchingAlt * 3 >= annotatedVariant.normalDNAReadCounts.nMatchingReference * 2)
							{
								hasNormal = true;
								alleleSpecificExpressionNormal = Math.Abs(((double)annotatedVariant.normalRNAReadCounts.nMatchingAlt / (annotatedVariant.normalRNAReadCounts.nMatchingReference + annotatedVariant.normalRNAReadCounts.nMatchingAlt)) * 2.0 - 1.0);
							}

							// add to data
							// AM TODO naming
							values.Add(new Temp(annotatedVariant.contig, annotatedVariant.locus, alleleSpecificExpression, 0, hasNormal, alleleSpecificExpressionNormal, 0));
						}

					}
				}
				else
				{
					// AM TODO regional expression
				}

				foreach (var value in values)
				{
					if (ASETools.isChromosomeAutosomal(value.chromosome))
					{
						wholeAutosomeRegionalExpression.AddTumorExpression(value.z_tumor, value.mu_tumor);
						if (value.hasNormal)
							wholeAutosomeRegionalExpression.AddNormalExpression(value.z_normal, value.mu_normal);
						foreach (var entry in allButThisChromosomeAutosomalRegionalExpressionState)
						{
							if (entry.Key != value.chromosome)
							{
								entry.Value.AddTumorExpression(value.z_tumor, value.mu_tumor);
								if (value.hasNormal)
									entry.Value.AddNormalExpression(value.z_normal, value.mu_normal);
							}
						}
					}

					int chromosomeId = ASETools.ChromosomeNameToIndex(value.chromosome);
					if (chromosomeId != -1)
					{
						perChromosomeRegionalExpressionState[chromosomeId].AddTumorExpression(value.z_tumor, value.mu_tumor);

						if (value.hasNormal)
							perChromosomeRegionalExpressionState[chromosomeId].AddNormalExpression(value.z_normal, value.mu_normal);
					}
					foreach (var geneLocation in geneLocationInformation.genesByChromosome[value.chromosome])
					{
						if (!geneExpressions.ContainsKey(geneLocation.hugoSymbol))
						{
							geneExpressions.Add(geneLocation.hugoSymbol, new ASETools.GeneExpression(geneLocation));
						}
						geneExpressions[geneLocation.hugoSymbol].AddRegionalExpression(value.offset, value.z_tumor, value.mu_tumor, true); // Recall that for allele-specifc expresion, z is really the level of allele-specific expression, not the expression z score.

						if (value.hasNormal)
							geneExpressions[geneLocation.hugoSymbol].AddRegionalExpression(value.offset, value.z_normal, value.mu_normal, false); // Recall that for allele-specifc expresion, z is really the level of allele-specific expression, not the expression z score.
					}
				}

				// AM TODO separate normal and tumor files

				//
				// Write the output file.
				//
				string directory = ASETools.GetDirectoryPathFromFullyQualifiedFilename(inputFilename);
				string analysisId = ASETools.GetAnalysisIdFromPathname(inputFilename);
				if ("" == directory || "" == analysisId)
				{
					Console.WriteLine("Couldn't parse input pathname, which is supposed to be absolute and include an analysis ID: " + inputFilename);
					continue;
				}

				var normalOutputFilename = directory + analysisId + (forAlleleSpecificExpression ? ASETools.normalAlleleSpecificGeneExpressionExtension : ASETools.geneExpressionExtension);
				var tumorOutputFilename = directory + analysisId + (forAlleleSpecificExpression ? ASETools.tumorAlleleSpecificGeneExpressionExtension : ASETools.geneExpressionExtension);

				var columnSuffix = forAlleleSpecificExpression ? "ase" : "z";

				var allExpressions = new List<ASETools.GeneExpression>();
				foreach (var expressionEntry in geneExpressions)
				{
					allExpressions.Add(expressionEntry.Value);
				}

				allExpressions.Sort(ASETools.GeneExpression.CompareByGeneName);

				// AM TODO: save "ExpressionNearMutations v3.1 " as first line
				ASETools.ASEExpressionFile.WriteFile(tumorOutputFilename, allExpressions, wholeAutosomeRegionalExpression, allButThisChromosomeAutosomalRegionalExpressionState,
perChromosomeRegionalExpressionState, false, 1, case_, columnSuffix, true);

				// If normal RNA exists, save file
				if (values.Where(r => r.hasNormal).Count() > 0)
				{
					ASETools.ASEExpressionFile.WriteFile(normalOutputFilename, allExpressions, wholeAutosomeRegionalExpression, allButThisChromosomeAutosomalRegionalExpressionState,
perChromosomeRegionalExpressionState, false, 1, case_, columnSuffix, false);
				}

				timer.Stop();
				lock (casesToProcess)
				{
					var nRemaining = casesToProcess.Count();
					Console.WriteLine("Processed case " + case_.case_id + " in " + (timer.ElapsedMilliseconds + 500) / 1000 + "s.  " + nRemaining + " remain" + ((1 == nRemaining) ? "s" : "") + " queued.");
				}
			} // while true
        }

        static void PrintUsageMessage()
        {
            Console.WriteLine("usage: ExpressionNearMutations -a casesToProcess");
            Console.WriteLine("-a means to use allele-specific expression, as opposed to total expression.");
        }

       // static ASETools.ASEConfirguation configuration;

        static void Main(string[] args)
        {
			var timer = new Stopwatch();
			timer.Start();

			var configuration = ASETools.ASEConfirguation.loadFromFile(args);

			if (null == configuration)
			{
				Console.WriteLine("Giving up because we were unable to load configuration.");
				return;
			}

			// only allow flag for allele-specific expression or case ids
			if (configuration.commandLineArgs.Count() == 0 || configuration.commandLineArgs.Count() == 1 && configuration.commandLineArgs[0] == "-a")
            {
                PrintUsageMessage();
                return;
            }

			var cases = ASETools.Case.LoadCases(configuration.casesFilePathname);

			if (null == cases)
			{
				Console.WriteLine("Unable to load cases file " + configuration.casesFilePathname + ".  You must generate cases before running ExpressionNearMutations.");
			}

			// set ASE flag, if specified
			bool forAlleleSpecificExpression = configuration.commandLineArgs[0] == "-a";

			int argsConsumed = 0;
			if (forAlleleSpecificExpression)
			{
				argsConsumed = 1;
			}


			int minExamplesPerRegion = 1;   // If there are fewer than this, then ignore the region.

			//
			// Now build the map of mutations by gene.
			//
			timer.Reset();
            timer.Start();

			// Get information for current genome build
			geneLocationInformation = new ASETools.GeneLocationsByNameAndChromosome(ASETools.readKnownGeneFile(ASETools.ASEConfirguation.defaultGeneLocationInformation));

            timer.Stop();
            Console.WriteLine("Loaded mutations in " + geneLocationInformation.genesByName.Count() + " genes in " + (timer.ElapsedMilliseconds + 500) / 1000 + "s.");

			List<ASETools.Case> casesToProcess = new List<ASETools.Case>();

			for (int i = argsConsumed; i < configuration.commandLineArgs.Count(); i++)
			{
				if (!cases.ContainsKey(configuration.commandLineArgs[i]))
				{
					Console.WriteLine(configuration.commandLineArgs[i] + " does not appear to be a case ID.  Ignoring.");
				}
				else
				{
					casesToProcess.Add(cases[configuration.commandLineArgs[i]]);
				}
				
			}

			Console.WriteLine("Processing " + casesToProcess.Count() + " cases");

			//
			// Process the runs in parallel
			//
			timer.Reset();
			timer.Start();

			var threads = new List<Thread>();
			for (int i = 0; i < Environment.ProcessorCount; i++)
			{
				threads.Add(new Thread(() => ProcessCases(casesToProcess, forAlleleSpecificExpression, minExamplesPerRegion)));
			}

			threads.ForEach(t => t.Start());
			threads.ForEach(t => t.Join());

			timer.Stop();
            Console.WriteLine("Processed " + (configuration.commandLineArgs.Count() - (forAlleleSpecificExpression ? 1 : 0)) + " experiments in " + (timer.ElapsedMilliseconds + 500) / 1000 + " seconds");
        }
    }
}
