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

        class RegionalExpressionState
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

        class GeneExpression 
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

		static ASETools.GeneLocationsByNameAndChromosome geneLocationInformation;

		static void writeColumnNames(StreamWriter outputFile, bool forAlleleSpecificExpression, ASETools.Case case_, bool isTumor)
		{
			string columnSuffix = forAlleleSpecificExpression ? "(Tumor ase)" : "(Tumor z)";

			if (!isTumor)
			{
				columnSuffix = forAlleleSpecificExpression ? "(Normal ase)" : "(Normal z)";

			}

			string columnSuffix_mu = isTumor ? "(Tumor mu)" : "(Normal mu)";

			for (int sizeIndex = 0; sizeIndex < GeneExpression.nRegionSizes; sizeIndex++)
			{
				outputFile.Write("\t" + GeneExpression.regionSizeByRegionSizeIndex[sizeIndex] + columnSuffix);

			}

			outputFile.Write("\tWhole Autosome " + columnSuffix);

			if (!forAlleleSpecificExpression)
			{
				for (int sizeIndex = 0; sizeIndex < GeneExpression.nRegionSizes; sizeIndex++)
				{
					outputFile.Write("\t" + GeneExpression.regionSizeByRegionSizeIndex[sizeIndex] + columnSuffix_mu);
				}
				outputFile.Write("\tWhole Autosome " + columnSuffix_mu);
			}

			for (int sizeIndex = 0; sizeIndex < GeneExpression.nRegionSizes; sizeIndex++)
			{
				outputFile.Write("\t" + GeneExpression.regionSizeByRegionSizeIndex[sizeIndex] + " exclusive " + columnSuffix);

			}

			outputFile.Write("\tWhole Autosome exclusive " + columnSuffix);

			if (!forAlleleSpecificExpression)
			{
				for (int sizeIndex = 0; sizeIndex < GeneExpression.nRegionSizes; sizeIndex++)
				{
					outputFile.Write("\t" + GeneExpression.regionSizeByRegionSizeIndex[sizeIndex] + " exclusive " + columnSuffix_mu);

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

		static void writeRow(StreamWriter outputFile, GeneExpression allExpression,
			RegionalExpressionState wholeAutosomeRegionalExpression,
			Dictionary<string, RegionalExpressionState> allButThisChromosomeAutosomalRegionalExpressionState,
			RegionalExpressionState[] perChromosomeRegionalExpressionState,
			bool forAlleleSpecificExpression, int minExamplesPerRegion)
		{
				// write tumor values
				for (int sizeIndex = 0; sizeIndex < GeneExpression.nRegionSizes; sizeIndex++)
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
					for (int sizeIndex = 0; sizeIndex < GeneExpression.nRegionSizes; sizeIndex++)
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

				for (int sizeIndex = 0; sizeIndex < GeneExpression.nRegionSizes; sizeIndex++)
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
					for (int sizeIndex = 0; sizeIndex < GeneExpression.nRegionSizes; sizeIndex++)
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
			for (int sizeIndex = 0; sizeIndex < GeneExpression.nRegionSizes; sizeIndex++)
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
				for (int sizeIndex = 0; sizeIndex < GeneExpression.nRegionSizes; sizeIndex++)
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

			for (int sizeIndex = 0; sizeIndex < GeneExpression.nRegionSizes; sizeIndex++)
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
				for (int sizeIndex = 0; sizeIndex < GeneExpression.nRegionSizes; sizeIndex++)
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
                var geneExpressions = new Dictionary<string, GeneExpression>();    
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
                        geneExpressions.Add(mafLine.Hugo_Symbol, new GeneExpression(geneLocationInformation.genesByName[mafLine.Hugo_Symbol]));
                    }
 
					// Increment the muation count by 1
                    geneExpressions[mafLine.Hugo_Symbol].mutationCount++;
                }


                var reader = ASETools.CreateStreamReaderWithRetry(inputFilename);

                var headerLine = reader.ReadLine();
                if (null == headerLine)
                {
                    Console.WriteLine("Empty input file " + inputFilename);
                    continue;
                }

                string line;
                int lineNumber = 1;
                if (!forAlleleSpecificExpression)
                {
                    if (headerLine.Substring(0, 20) != "RegionalExpression v")
                    {
                        Console.WriteLine("Corrupt header line in file '" + inputFilename + "', line: " + headerLine);
                        continue;
                    }

                    if (headerLine.Substring(20, 2) != "3.")
                    {
                        Console.WriteLine("Unsupported version in file '" + inputFilename + "', header line: " + headerLine);
                        continue;
                    }
                    line = reader.ReadLine();   // The NumContigs line, which we just ignore
                    line = reader.ReadLine();   // The column header line, which we just ignore

                    if (null == line)
                    {
                        Console.WriteLine("Truncated file '" + inputFilename + "' ends after header line.");
                        continue;
                    }

                    lineNumber = 3;
                }
				// Variables storing expression state other than same-chromosome regional.

				// One expression state for whole autosome
                var wholeAutosomeRegionalExpression = new RegionalExpressionState();
				// Expression state for each chromosome, which will exclude the chromosome that the gene resides on
                var allButThisChromosomeAutosomalRegionalExpressionState = new Dictionary<string, RegionalExpressionState>();   // "This chromosome" is the dictionary key
				// Expression state for each chromosome, which will include the chromsome that the gene resides on
				var perChromosomeRegionalExpressionState = new RegionalExpressionState[ASETools.nHumanNuclearChromosomes];

                for (int whichChromosome = 0; whichChromosome < ASETools.nHumanNuclearChromosomes; whichChromosome++)
                {
                    perChromosomeRegionalExpressionState[whichChromosome] = new RegionalExpressionState();
                }

				foreach (var geneEntry in geneLocationInformation.genesByName)
                {
                    var chromosome = geneEntry.Value.chromosome;
                    if (!allButThisChromosomeAutosomalRegionalExpressionState.ContainsKey(chromosome))
                    {
                        allButThisChromosomeAutosomalRegionalExpressionState.Add(chromosome, new RegionalExpressionState());
                    }
                }

                bool seenDone = false;
                while (null != (line = reader.ReadLine()))
                {
                    lineNumber++;

                    if (seenDone)
                    {
                        Console.WriteLine("Saw data after **done** in file " + inputFilename + "', line " + lineNumber + ": " + line);
                        break;
                    }

                    if (line == "**done**")
                    {
                        seenDone = true;
                        continue;
                    }

 
                    string chromosome;
                    int offset;

                    // For allele-specific expression
                    double nMatchingReferenceTumorDNA = 0;
                    double nMatchingVariantTumorDNA = 0;
                    double nMatchingReferenceTumorRNA = 0;
                    double nMatchingVariantTumorRNA = 0;

					// normal allele specific expression
					double nMatchingReferenceNormalDNA = 0;
					double nMatchingVariantNormalDNA = 0;
					double nMatchingReferenceNormalRNA = 0;
					double nMatchingVariantNormalRNA = 0;
					bool normalRNAExists = false;

					// for regional expression
					double z_tumor = 0;
                    double mu_tumor = 0;

					double z_normal = 0;
					double mu_normal = 0;

					try {
                        if (forAlleleSpecificExpression) {
                            var alleleData = ASETools.AnnotatedVariant.fromText(line);

                            if (null == alleleData)
                            {
                                Console.WriteLine("Error parsing input line " + lineNumber + " in file " + inputFilename);
                                break;
                            }

                            chromosome = alleleData.contig;
                            offset = alleleData.locus;
                            nMatchingReferenceTumorDNA = alleleData.tumorDNAReadCounts.nMatchingReference;
                            nMatchingVariantTumorDNA = alleleData.tumorDNAReadCounts.nMatchingAlt;
                            nMatchingReferenceTumorRNA = alleleData.tumorRNAReadCounts.nMatchingReference;
                            nMatchingVariantTumorRNA = alleleData.tumorRNAReadCounts.nMatchingAlt;

							// Variables for normal ASE
							try
							{
								nMatchingReferenceNormalRNA = alleleData.normalRNAReadCounts.nMatchingReference;
								nMatchingVariantNormalRNA = alleleData.normalRNAReadCounts.nMatchingAlt;
								normalRNAExists = true;
							}
							catch (Exception) {
								// no normal RNA for this case. Skipping
							}

						}
						else {

                            var fields = line.Split('\t');
                            if (fields.Count() != 13)
                            {
                                Console.WriteLine("Badly formatted data line in file '" + inputFilename + "', line " + lineNumber + ": " + line);
                                break;
                            }

                            chromosome = fields[0]; // in chr form
                            offset = Convert.ToInt32(fields[1]);
                            z_tumor = Convert.ToDouble(fields[11]);
                            mu_tumor = Convert.ToDouble(fields[12]);

                            int nBasesExpressedWithBaselineExpression = Convert.ToInt32(fields[3]);
                            int nBasesUnexpressedWithBaselineExpression = Convert.ToInt32(fields[7]);

                            if (0 == nBasesExpressedWithBaselineExpression && 0 == nBasesUnexpressedWithBaselineExpression)
                            {
                                //
                                // No baseline expression for this region, skip it.
                                //
                                continue;
                            }
                        }
                    }
                    catch (FormatException)
                    {
                        Console.WriteLine("Format exception parsing data line in file '" + inputFilename + "', line " + lineNumber + ": " + line);
                        break;
                    }

					// Remove chr prefix from chromosome, if present
					if (!geneLocationInformation.genesByChromosome.ContainsKey(chromosome))
					{
						//
						// Try reversing the "chr" state of the chromosome.
						//

						if (chromosome.Count() > 3 && chromosome.Substring(0, 3) == "chr")
						{
							chromosome = chromosome.Substring(3);
						}
						else
						{
							chromosome = "chr" + chromosome;
						}
					}

					if (geneLocationInformation.genesByChromosome.ContainsKey(chromosome) && 
                        (!forAlleleSpecificExpression ||                                    // We only keep samples for allele specific expression if they meet certain criteria, to wit:
                            nMatchingReferenceTumorDNA + nMatchingVariantTumorDNA >= 10 &&            // We have at least 10 DNA reads
							nMatchingReferenceTumorRNA + nMatchingVariantTumorRNA >= 10 &&            // We have at least 10 RNA reads
							nMatchingReferenceTumorDNA * 3 >= nMatchingVariantTumorDNA * 2 &&         // It's not more than 2/3 variant DNA
							nMatchingVariantTumorDNA * 3 >= nMatchingReferenceTumorDNA * 2))          // It's not more than 2/3 reference DNA
                    {
                        if (forAlleleSpecificExpression)
                        {
                            double rnaFractionTumor = nMatchingVariantTumorRNA / (nMatchingReferenceTumorRNA + nMatchingVariantTumorRNA);

							//
							// Now convert to the amount of allele-specific expression.  50% is no ASE, while 0 or 100% is 100% ASE.
							//
							z_tumor = Math.Abs(rnaFractionTumor * 2.0 - 1.0); // Not really z, really alleleSpecificExpression
                            mu_tumor = 0;

							// If we have the normal DNA and RNA for this sample, compute the normal ASE
							if (normalRNAExists && nMatchingReferenceNormalDNA + nMatchingVariantNormalDNA >= 10 &&   // We have at least 10 DNA reads
							nMatchingReferenceNormalRNA + nMatchingVariantNormalRNA >= 10 &&            // We have at least 10 RNA reads
							nMatchingReferenceNormalDNA * 3 >= nMatchingVariantNormalDNA * 2 &&         // It's not more than 2/3 variant DNA
							nMatchingVariantNormalDNA * 3 >= nMatchingReferenceNormalDNA * 2)
							{
								double rnaFractionNormal = nMatchingVariantNormalRNA / (nMatchingReferenceNormalRNA + nMatchingVariantNormalRNA);

								// Convert to ASE
								z_normal = Math.Abs(rnaFractionNormal * 2.0 - 1.0); 
								mu_normal = 0;
							}
						}

                        if (ASETools.isChromosomeAutosomal(chromosome))
                        {
                            wholeAutosomeRegionalExpression.AddTumorExpression(z_tumor, mu_tumor);

							if (normalRNAExists)
								wholeAutosomeRegionalExpression.AddNormalExpression(z_normal, mu_normal);

							foreach (var entry in allButThisChromosomeAutosomalRegionalExpressionState) 
                            {
                                if (entry.Key != chromosome) {
                                    entry.Value.AddTumorExpression(z_tumor, mu_tumor);
									if (normalRNAExists)
										entry.Value.AddNormalExpression(z_normal, mu_normal);
									
								}
                            }
                        }

                        int chromosomeId = ASETools.ChromosomeNameToIndex(chromosome);
                        if (chromosomeId != -1)
                        {
                            perChromosomeRegionalExpressionState[chromosomeId].AddTumorExpression(z_tumor, mu_tumor);
							if (normalRNAExists)
								perChromosomeRegionalExpressionState[chromosomeId].AddNormalExpression(z_normal, mu_normal);
						}

                        foreach (var geneLocation in geneLocationInformation.genesByChromosome[chromosome])
                        {
                            if (!geneExpressions.ContainsKey(geneLocation.hugoSymbol))
                            {
                                geneExpressions.Add(geneLocation.hugoSymbol, new GeneExpression(geneLocation));
                            }

                            geneExpressions[geneLocation.hugoSymbol].AddRegionalExpression(offset, z_tumor, mu_tumor, true); // Recall that for allele-specifc expresion, z is really the level of allele-specific expression, not the expression z score.
							if (normalRNAExists)
								geneExpressions[geneLocation.hugoSymbol].AddRegionalExpression(offset, z_normal, mu_normal, false); // Recall that for allele-specifc expresion, z is really the level of allele-specific expression, not the expression z score.
						}
					}
                } // for each line in the input file

                if (!seenDone)
                {
                    Console.WriteLine("Truncated input file " + inputFilename);
                    continue;
                }

                //
                // Write the output file.
                //
                string directory = ASETools.GetDirectoryPathFromFullyQualifiedFilename(inputFilename);
                string analysisId = ASETools.GetAnalysisIdFromPathname(inputFilename);
                if ("" == directory || "" == analysisId) {
                    Console.WriteLine("Couldn't parse input pathname, which is supposed to be absolute and include an analysis ID: " + inputFilename);
                    continue;
                }

                var outputFilename = directory + analysisId + (forAlleleSpecificExpression ? ASETools.alleleSpecificGeneExpressionExtension : ASETools.geneExpressionExtension);

                var outputFile = ASETools.CreateStreamWriterWithRetry(outputFilename);

				outputFile.WriteLine("ExpressionNearMutations v3.1 " + case_.case_id + (forAlleleSpecificExpression ? " -a" : "")); // v3.1 uses ucsc gene locations
				outputFile.Write("Gene name\tnon-silent mutation count");

				writeColumnNames(outputFile, forAlleleSpecificExpression, case_, true);
				writeColumnNames(outputFile, forAlleleSpecificExpression, case_, false);

                outputFile.WriteLine();

                var allExpressions = new List<GeneExpression>();
                foreach (var expressionEntry in geneExpressions)
                {
                    allExpressions.Add(expressionEntry.Value);
                }

                allExpressions.Sort(GeneExpression.CompareByGeneName);

				for (int i = 0; i < allExpressions.Count(); i++)
				{
					outputFile.Write(ASETools.ConvertToExcelString(allExpressions[i].geneLocationInfo.hugoSymbol) + "\t" + allExpressions[i].mutationCount);

					writeRow(outputFile, allExpressions[i],
						wholeAutosomeRegionalExpression,
						allButThisChromosomeAutosomalRegionalExpressionState,
						perChromosomeRegionalExpressionState,
						forAlleleSpecificExpression,
						minExamplesPerRegion);
					//for (int sizeIndex = 0; sizeIndex < GeneExpression.nRegionSizes; sizeIndex++)
					//{
					//	if (allExpressions[i].regionalExpressionState[sizeIndex].nRegionsIncludedTumor >= minExamplesPerRegion)
					//	{
					//		outputFile.Write("\t" + allExpressions[i].regionalExpressionState[sizeIndex].totalTumorExpression / allExpressions[i].regionalExpressionState[sizeIndex].nRegionsIncludedTumor);

					//		// Write out normal RNA results, if present
					//		if (allExpressions[i].regionalExpressionState[sizeIndex].nRegionsIncludedNormal >= minExamplesPerRegion)
					//		{
					//			outputFile.Write("\t" + allExpressions[i].regionalExpressionState[sizeIndex].totalNormalExpression / allExpressions[i].regionalExpressionState[sizeIndex].nRegionsIncludedNormal);
					//		}
					//		else
					//		{
					//			outputFile.Write("\t*");
					//		}
					//	}
					//	else
					//	{
					//		outputFile.Write("\t*\t*");
					//	}
					//}

					//if (wholeAutosomeRegionalExpression.nRegionsIncludedTumor >= minExamplesPerRegion)
					//{
					//	outputFile.Write("\t" + wholeAutosomeRegionalExpression.totalTumorExpression / wholeAutosomeRegionalExpression.nRegionsIncludedTumor);

					//	if (wholeAutosomeRegionalExpression.nRegionsIncludedNormal >= minExamplesPerRegion)
					//	{
					//		outputFile.Write("\t" + wholeAutosomeRegionalExpression.totalNormalExpression / wholeAutosomeRegionalExpression.nRegionsIncludedNormal);
					//	}
					//	else
					//	{
					//		outputFile.Write("\t*");
					//	}
					//}
					//else
					//{
					//	outputFile.Write("\t*\t*");
					//}

					//if (!forAlleleSpecificExpression)
					//{
					//	for (int sizeIndex = 0; sizeIndex < GeneExpression.nRegionSizes; sizeIndex++)
					//	{
					//		if (allExpressions[i].regionalExpressionState[sizeIndex].nRegionsIncludedTumor >= minExamplesPerRegion)
					//		{
					//			outputFile.Write("\t" + allExpressions[i].regionalExpressionState[sizeIndex].totalMeanTumorExpression / allExpressions[i].regionalExpressionState[sizeIndex].nRegionsIncludedTumor);

					//			if (allExpressions[i].regionalExpressionState[sizeIndex].nRegionsIncludedNormal >= minExamplesPerRegion)
					//			{
					//				outputFile.Write("\t" + allExpressions[i].regionalExpressionState[sizeIndex].totalMeanNormalExpression / allExpressions[i].regionalExpressionState[sizeIndex].nRegionsIncludedNormal);
					//			}
					//			else
					//			{
					//				outputFile.Write("\t*");
					//			}
					//		}
					//		else
					//		{
					//			outputFile.Write("\t*\t*");
					//		}
					//	}

					//	if (wholeAutosomeRegionalExpression.nRegionsIncludedTumor >= minExamplesPerRegion)
					//	{
					//		outputFile.Write("\t" + wholeAutosomeRegionalExpression.totalMeanTumorExpression / wholeAutosomeRegionalExpression.nRegionsIncludedTumor);

					//		if (wholeAutosomeRegionalExpression.nRegionsIncludedNormal >= minExamplesPerRegion)
					//		{
					//			outputFile.Write("\t" + wholeAutosomeRegionalExpression.totalMeanNormalExpression / wholeAutosomeRegionalExpression.nRegionsIncludedNormal);
					//		}
					//		else
					//		{
					//			outputFile.Write("\t*");
					//		}
					//	}
					//	else
					//	{
					//		outputFile.Write("\t*\t*");
					//	}
					//}

					//for (int sizeIndex = 0; sizeIndex < GeneExpression.nRegionSizes; sizeIndex++)
					//{
					//	if (allExpressions[i].exclusiveRegionalExpressionState[sizeIndex].nRegionsIncludedTumor >= minExamplesPerRegion)
					//	{
					//		outputFile.Write("\t" + allExpressions[i].exclusiveRegionalExpressionState[sizeIndex].totalTumorExpression / allExpressions[i].exclusiveRegionalExpressionState[sizeIndex].nRegionsIncludedTumor);

					//		if (allExpressions[i].exclusiveRegionalExpressionState[sizeIndex].nRegionsIncludedNormal >= minExamplesPerRegion)
					//		{
					//			outputFile.Write("\t" + allExpressions[i].exclusiveRegionalExpressionState[sizeIndex].totalNormalExpression / allExpressions[i].exclusiveRegionalExpressionState[sizeIndex].nRegionsIncludedNormal);
					//		}
					//		else
					//		{
					//			outputFile.Write("\t*");
					//		}
					//	}
					//	else
					//	{
					//		outputFile.Write("\t*\t*");
					//	}
					//}

					//if (allButThisChromosomeAutosomalRegionalExpressionState[allExpressions[i].geneLocationInfo.chromosome].nRegionsIncludedTumor >= minExamplesPerRegion)
					//{
					//	outputFile.Write("\t" + allButThisChromosomeAutosomalRegionalExpressionState[allExpressions[i].geneLocationInfo.chromosome].totalTumorExpression / allButThisChromosomeAutosomalRegionalExpressionState[allExpressions[i].geneLocationInfo.chromosome].nRegionsIncludedTumor);

					//	if (allButThisChromosomeAutosomalRegionalExpressionState[allExpressions[i].geneLocationInfo.chromosome].nRegionsIncludedNormal >= minExamplesPerRegion)
					//	{
					//		outputFile.Write("\t" + allButThisChromosomeAutosomalRegionalExpressionState[allExpressions[i].geneLocationInfo.chromosome].totalNormalExpression / allButThisChromosomeAutosomalRegionalExpressionState[allExpressions[i].geneLocationInfo.chromosome].nRegionsIncludedNormal);
					//	}
					//	else
					//	{
					//		outputFile.Write("\t*");
					//	}

					//}
					//else
					//{
					//	outputFile.Write("\t*\t*");
					//}

					//if (!forAlleleSpecificExpression)
					//{
					//	for (int sizeIndex = 0; sizeIndex < GeneExpression.nRegionSizes; sizeIndex++)
					//	{
					//		if (allExpressions[i].exclusiveRegionalExpressionState[sizeIndex].nRegionsIncludedTumor >= minExamplesPerRegion)
					//		{
					//			outputFile.Write("\t" + allExpressions[i].exclusiveRegionalExpressionState[sizeIndex].totalMeanTumorExpression / allExpressions[i].exclusiveRegionalExpressionState[sizeIndex].nRegionsIncludedTumor);

					//			if (allExpressions[i].exclusiveRegionalExpressionState[sizeIndex].nRegionsIncludedNormal >= minExamplesPerRegion)
					//			{
					//				outputFile.Write("\t" + allExpressions[i].exclusiveRegionalExpressionState[sizeIndex].totalMeanNormalExpression / allExpressions[i].exclusiveRegionalExpressionState[sizeIndex].nRegionsIncludedNormal);
					//			}
					//			else
					//			{
					//				outputFile.Write("\t*");
					//			}
					//		}
					//		else
					//		{
					//			outputFile.Write("\t*\t*");
					//		}
					//	}

					//	if (allButThisChromosomeAutosomalRegionalExpressionState[allExpressions[i].geneLocationInfo.chromosome].nRegionsIncludedTumor >= minExamplesPerRegion)
					//	{
					//		outputFile.Write("\t" + allButThisChromosomeAutosomalRegionalExpressionState[allExpressions[i].geneLocationInfo.chromosome].totalMeanTumorExpression / allButThisChromosomeAutosomalRegionalExpressionState[allExpressions[i].geneLocationInfo.chromosome].nRegionsIncludedTumor);

					//		if (allButThisChromosomeAutosomalRegionalExpressionState[allExpressions[i].geneLocationInfo.chromosome].nRegionsIncludedNormal >= minExamplesPerRegion)
					//		{
					//			outputFile.Write("\t" + allButThisChromosomeAutosomalRegionalExpressionState[allExpressions[i].geneLocationInfo.chromosome].totalMeanNormalExpression / allButThisChromosomeAutosomalRegionalExpressionState[allExpressions[i].geneLocationInfo.chromosome].nRegionsIncludedNormal);
					//		}
					//		else
					//		{
					//			outputFile.Write("\t*");
					//		}
					//	}
					//	else
					//	{
					//		outputFile.Write("\t*\t*");
					//	}
					//}

					//for (int whichChromosome = 0; whichChromosome < ASETools.nHumanNuclearChromosomes; whichChromosome++)
					//{
					//	if (perChromosomeRegionalExpressionState[whichChromosome].nRegionsIncludedTumor >= minExamplesPerRegion)
					//	{
					//		outputFile.Write("\t" + perChromosomeRegionalExpressionState[whichChromosome].totalTumorExpression / perChromosomeRegionalExpressionState[whichChromosome].nRegionsIncludedTumor);

					//		if (perChromosomeRegionalExpressionState[whichChromosome].nRegionsIncludedNormal >= minExamplesPerRegion)
					//		{
					//			outputFile.Write("\t" + perChromosomeRegionalExpressionState[whichChromosome].totalNormalExpression / perChromosomeRegionalExpressionState[whichChromosome].nRegionsIncludedNormal);
					//		}
					//		else
					//		{
					//			outputFile.Write("\t*");
					//		}
					//	}
					//	else
					//	{
					//		outputFile.Write("\t*\t*");
					//	}
					//}

					//if (!forAlleleSpecificExpression)
					//{
					//	for (int whichChromosome = 0; whichChromosome < ASETools.nHumanNuclearChromosomes; whichChromosome++)
					//	{
					//		if (perChromosomeRegionalExpressionState[whichChromosome].nRegionsIncludedTumor >= minExamplesPerRegion)
					//		{
					//			outputFile.Write("\t" + perChromosomeRegionalExpressionState[whichChromosome].totalMeanTumorExpression / perChromosomeRegionalExpressionState[whichChromosome].nRegionsIncludedTumor);

					//			if (perChromosomeRegionalExpressionState[whichChromosome].nRegionsIncludedNormal >= minExamplesPerRegion)
					//			{
					//				outputFile.Write("\t" + perChromosomeRegionalExpressionState[whichChromosome].totalMeanNormalExpression / perChromosomeRegionalExpressionState[whichChromosome].nRegionsIncludedNormal);
					//			}
					//			else
					//			{
					//				outputFile.Write("\t*");
					//			}

					//		}
					//		else
					//		{
					//			outputFile.Write("\t*\t*");
					//		}
					//	}
					//}

					outputFile.WriteLine();
				} // for each gene

				outputFile.WriteLine("**done**");
                outputFile.Close();

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
			//timer.Reset();
			//timer.Start();

			//var threads = new List<Thread>();
			//for (int i = 0; i < Environment.ProcessorCount; i++)
			//{
			//	threads.Add(new Thread(() => ProcessCases(casesToProcess, forAlleleSpecificExpression, minExamplesPerRegion)));
			//}
			ProcessCases(casesToProcess, forAlleleSpecificExpression, minExamplesPerRegion);

			//threads.ForEach(t => t.Start());
   //         threads.ForEach(t => t.Join());

            timer.Stop();
            Console.WriteLine("Processed " + (configuration.commandLineArgs.Count() - (forAlleleSpecificExpression ? 1 : 0)) + " experiments in " + (timer.ElapsedMilliseconds + 500) / 1000 + " seconds");
        }
    }
}
