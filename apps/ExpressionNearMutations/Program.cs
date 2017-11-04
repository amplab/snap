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

		// Generic class to hold tumor/normal signal at a locus. Use for ASE, expression and ASM values.
		class RegionalSignal
		{
			// One expression state for whole autosome
			public ASETools.RegionalExpressionState wholeAutosomeRegionalExpression = new ASETools.RegionalExpressionState();
			// Expression state for each chromosome, which will exclude the chromosome that the gene resides on
			public Dictionary<string, ASETools.RegionalExpressionState> allButThisChromosomeAutosomalRegionalExpressionState = new Dictionary<string, ASETools.RegionalExpressionState>();   // "This chromosome" is the dictionary key
																																															 // Expression state for each chromosome, which will include the chromsome that the gene resides on
			public ASETools.RegionalExpressionState[] perChromosomeRegionalExpressionState = new ASETools.RegionalExpressionState[ASETools.nHumanNuclearChromosomes];

			public Dictionary<string, ASETools.GeneExpression> geneExpressions;

			public bool hasNormal;


		    public RegionalSignal(Dictionary<string, ASETools.GeneExpression> geneExpressions_, ASETools.GeneLocationsByNameAndChromosome geneLocationInformation, bool hasNormal_)
		    {
			    geneExpressions = geneExpressions_;
			    hasNormal = hasNormal_;

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
            } // RegionalSignal

            public void Add(string chromosome, int locus, double z_tumor, double mu_tumor, bool hasNormalValues, double z_normal = 0, double mu_normal = 0)
		    {
			    if (ASETools.isChromosomeAutosomal(chromosome))
			    {
				    // Autosome
				    wholeAutosomeRegionalExpression.AddTumorExpression(z_tumor, mu_tumor);
				    if (hasNormalValues)
					    wholeAutosomeRegionalExpression.AddNormalExpression(z_normal, mu_normal);
				    // Exclusive whole chromosomes
				    foreach (var entry in allButThisChromosomeAutosomalRegionalExpressionState)
				    {
					    if (entry.Key != chromosome)
					    {
						    entry.Value.AddTumorExpression(z_tumor, mu_tumor);
						    if (hasNormalValues)
							    entry.Value.AddNormalExpression(z_normal, mu_normal);
					    }
				    }
			    }

			    int chromosomeId = ASETools.ChromosomeNameToIndex(chromosome);
			    // Whole chromosomes
			    if (chromosomeId != -1)
			    {
				    perChromosomeRegionalExpressionState[chromosomeId].AddTumorExpression(z_tumor, mu_tumor);

				    if (hasNormalValues)
					    perChromosomeRegionalExpressionState[chromosomeId].AddNormalExpression(z_normal, mu_normal);
			    }
			    // Regional, same chromosome
			    foreach (var geneLocation in geneLocationInformation.genesByChromosome[chromosome])
			    {
				    if (!geneExpressions.ContainsKey(geneLocation.hugoSymbol))
				    {
					    geneExpressions.Add(geneLocation.hugoSymbol, new ASETools.GeneExpression(geneLocation));
				    }
				    geneExpressions[geneLocation.hugoSymbol].AddRegionalExpression(locus, z_tumor, mu_tumor, true); // Recall that for allele-specifc expresion, z is really the level of allele-specific expression, not the expression z score.

				    if (hasNormalValues)
					    geneExpressions[geneLocation.hugoSymbol].AddRegionalExpression(locus, z_normal, mu_normal, false); // Recall that for allele-specifc expresion, z is really the level of allele-specific expression, not the expression z score.
			    }
            } // Add
        } // RegionalSignal

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

				if (inputFilename == "" || case_.tumor_copy_number_filename == "")
				{
					Console.WriteLine("Case " + case_.case_id + " doesn't have an input or copy number file yet.");
					continue;
				}

                // Used to filter out annotated variants in regions with copy number variation
                // The cutoff we use here for choosing copy number variation is very loose (1.0).
                // Setting this value to 1.0 loses 10-20% of annotated variants, depending on the sample.
                var copyNumberVariation = ASETools.CopyNumberVariation.ReadFile(case_.tumor_copy_number_filename, case_.tumor_copy_number_file_id).Where(r => Math.Abs(r.Segment_Mean) > 1.0).ToList();
                List<ASETools.CopyNumberVariation> normalCopyNumberVariation = null;
                if (case_.normal_copy_number_filename != "")
                {
                    normalCopyNumberVariation = ASETools.CopyNumberVariation.ReadFile(case_.normal_copy_number_filename, case_.normal_copy_number_file_id).Where(r => Math.Abs(r.Segment_Mean) > 1.0).ToList();
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
				} // foreach maf line

				List<ASETools.AnnotatedVariant> annotatedVariants = null;
				List<ASETools.RegionalExpressionLine> regionalExpressionLines = null;

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
					regionalExpressionLines = ASETools.Region.readFile(inputFilename);

					if (regionalExpressionLines == null)
					{
						Console.WriteLine("Failed to read " + case_.regional_expression_filename + ". Giving up on case " + case_.case_id);
						continue;
					}
				}

				// Stores all regional information for expression values
				var hasNormal = forAlleleSpecificExpression ? annotatedVariants.Where(r => r.normalRNAReadCounts != null).Count() > 0: false;
				var expressionValues = new RegionalSignal(geneExpressions, geneLocationInformation, hasNormal);

				if (forAlleleSpecificExpression)
				{
					foreach (var annotatedVariant in annotatedVariants)
					{
						if (geneLocationInformation.genesByChromosome.ContainsKey(annotatedVariant.contig) &&
                            annotatedVariant.IsASECandidate(true, copyNumberVariation, configuration, perGeneASEMap, geneMap) && !annotatedVariant.somaticMutation)
						{
							// If we have the normal DNA and RNA for this sample, compute the normal ASE
							if (annotatedVariant.normalRNAReadCounts != null && annotatedVariant.IsASECandidate(false, normalCopyNumberVariation, configuration, perGeneASEMap, geneMap))
							{
								// add to data with normal
								expressionValues.Add(annotatedVariant.contig, annotatedVariant.locus, annotatedVariant.GetTumorAlleleSpecificExpression(), 0, true, annotatedVariant.GetNormalAlleleSpecificExpression(), 0);
							}
							else
                            {
								// add to data without normal
								expressionValues.Add(annotatedVariant.contig, annotatedVariant.locus, annotatedVariant.GetTumorAlleleSpecificExpression(), 0, false);
							}
						}
					}
				}
				else
				{
					foreach (var region in regionalExpressionLines)
					{
						if (0 == region.nbases_expressed_baseline && 0 == region.nbases_baseline_other_samples)
						{
							//
							// No baseline expression for this region, skip it.
							//
							continue;
						}

						if (geneLocationInformation.genesByChromosome.ContainsKey(region.contig))
						{
							expressionValues.Add(region.contig,  region.contig_offset, region.mean_z, region.mean_mu, false);
						}
					}
				} // not for ASE (i.e., total expression)

				//
				// Write the output file for tumor and normal, if data exists.
				//
				string directory = ASETools.GetDirectoryPathFromFullyQualifiedFilename(inputFilename);
				string analysisId = ASETools.GetAnalysisIdFromPathname(inputFilename);
				if ("" == directory || "" == analysisId)
				{
					Console.WriteLine("Couldn't parse input pathname, which is supposed to be absolute and include an analysis ID: " + inputFilename);
					continue;
				}

				var tumorOutputFilename = directory + analysisId + (forAlleleSpecificExpression ? ASETools.tumorAlleleSpecificGeneExpressionExtension : ASETools.geneExpressionExtension);

				var columnSuffix = forAlleleSpecificExpression ? "ase" : "z";

				var allExpressions = new List<ASETools.GeneExpression>();
				foreach (var expressionEntry in expressionValues.geneExpressions)
				{
					allExpressions.Add(expressionEntry.Value);
				}

				allExpressions.Sort(ASETools.GeneExpression.CompareByGeneName);

				var fileHeader = "ExpressionNearMutations v3.1 " + case_.case_id;
				var output = new ASETools.RegionalSignalFile(fileHeader,
					allExpressions, 
					expressionValues.wholeAutosomeRegionalExpression, 
					expressionValues.allButThisChromosomeAutosomalRegionalExpressionState,
					expressionValues.perChromosomeRegionalExpressionState);

				// ASE values do not have mu 
				var printMu = !forAlleleSpecificExpression;
				output.WriteFile(tumorOutputFilename, printMu, 1, columnSuffix, true);

				// If normal RNA exists, save file
				if (expressionValues.hasNormal)
				{
					var normalOutputFilename = directory + analysisId + (forAlleleSpecificExpression ? ASETools.normalAlleleSpecificGeneExpressionExtension : ASETools.geneExpressionExtension);
					output.WriteFile(normalOutputFilename, printMu, 1, columnSuffix, false);
				}

				timer.Stop();
				lock (casesToProcess)
				{
					var nRemaining = casesToProcess.Count();
					Console.WriteLine("Processed case " + case_.case_id + " in " + (timer.ElapsedMilliseconds + 500) / 1000 + "s.  " + nRemaining + " remain" + ((1 == nRemaining) ? "s" : "") + " queued.");
				}
			} // while true (i.e., for each case)
        } // ProcessCases

        static void PrintUsageMessage()
        {
            Console.WriteLine("usage: ExpressionNearMutations -a casesToProcess");
            Console.WriteLine("-a means to use allele-specific expression, as opposed to total expression.");
        }

        static ASETools.Configuration configuration;
        static ASETools.GeneMap geneMap;
        static Dictionary<string, ASETools.ASEMapPerGeneLine> perGeneASEMap;

        static void Main(string[] args)
        {
			var timer = new Stopwatch();
			timer.Start();

			configuration = ASETools.Configuration.loadFromFile(args);

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

            perGeneASEMap = ASETools.ASEMapPerGeneLine.ReadFromFileToDictionary(configuration.finalResultsDirectory + ASETools.PerGeneASEMapFilename);

            if (null == perGeneASEMap)
            {
                Console.WriteLine("You must first create the per-gene ASE map in " + configuration.finalResultsDirectory + ASETools.PerGeneASEMapFilename);
                return;
            }

            geneLocationInformation = new ASETools.GeneLocationsByNameAndChromosome(ASETools.readKnownGeneFile(configuration.geneLocationInformationFilename));
            geneMap = new ASETools.GeneMap(geneLocationInformation.genesByName);

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
			geneLocationInformation = new ASETools.GeneLocationsByNameAndChromosome(ASETools.readKnownGeneFile(configuration.geneLocationInformationFilename));

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
            Console.WriteLine("Processed " + (configuration.commandLineArgs.Count() - (forAlleleSpecificExpression ? 1 : 0)) + " cases in " + (timer.ElapsedMilliseconds + 500) / 1000 + " seconds");
        }
    }
}
