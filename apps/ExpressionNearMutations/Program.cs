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
    public class Program
    {


		static ASETools.GeneLocationsByNameAndChromosome geneLocationInformation;

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

				var annotatedVariants = ASETools.AnnotatedVariant.readFile(inputFilename);
				bool normalRNAExists = false;

				foreach (var alleleData in annotatedVariants)
                {
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

					// for regional expression
					double z_tumor = 0;
                    double mu_tumor = 0;

					double z_normal = 0;
					double mu_normal = 0;

                    if (forAlleleSpecificExpression) {

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

						throw new Exception("not implemented");
                        //var fields = line.Split('\t');
                        //if (fields.Count() != 13)
                        //{
                        //    Console.WriteLine("Badly formatted data line in file '" + inputFilename + "', line " + lineNumber + ": " + line);
                        //    break;
                        //}

                        //chromosome = fields[0]; // in chr form
                        //offset = Convert.ToInt32(fields[1]);
                        //z_tumor = Convert.ToDouble(fields[11]);
                        //mu_tumor = Convert.ToDouble(fields[12]);

                        //int nBasesExpressedWithBaselineExpression = Convert.ToInt32(fields[3]);
                        //int nBasesUnexpressedWithBaselineExpression = Convert.ToInt32(fields[7]);

                        //if (0 == nBasesExpressedWithBaselineExpression && 0 == nBasesUnexpressedWithBaselineExpression)
                        //{
                        //    //
                        //    // No baseline expression for this region, skip it.
                        //    //
                        //    continue;
                        //}
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
                                geneExpressions.Add(geneLocation.hugoSymbol, new ASETools.GeneExpression(geneLocation));
                            }

                            geneExpressions[geneLocation.hugoSymbol].AddRegionalExpression(offset, z_tumor, mu_tumor, true); // Recall that for allele-specifc expresion, z is really the level of allele-specific expression, not the expression z score.
							if (normalRNAExists)
								geneExpressions[geneLocation.hugoSymbol].AddRegionalExpression(offset, z_normal, mu_normal, false); // Recall that for allele-specifc expresion, z is really the level of allele-specific expression, not the expression z score.
						}
					}
                } // for each line in the input file

                //
                // Write the output file.
                //
                string directory = ASETools.GetDirectoryPathFromFullyQualifiedFilename(inputFilename);
                string analysisId = ASETools.GetAnalysisIdFromPathname(inputFilename);
                if ("" == directory || "" == analysisId) {
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



				ASETools.ASEExpressionFile.WriteFile(tumorOutputFilename, allExpressions, wholeAutosomeRegionalExpression, allButThisChromosomeAutosomalRegionalExpressionState,
perChromosomeRegionalExpressionState, false, 1, case_, columnSuffix, true);

				if (normalRNAExists)
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

			//threads.ForEach(t => t.Start());
			//threads.ForEach(t => t.Join());

			ProcessCases(casesToProcess, forAlleleSpecificExpression, minExamplesPerRegion);

			timer.Stop();
            Console.WriteLine("Processed " + (configuration.commandLineArgs.Count() - (forAlleleSpecificExpression ? 1 : 0)) + " experiments in " + (timer.ElapsedMilliseconds + 500) / 1000 + " seconds");
        }
    }
}
