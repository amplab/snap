using System;
using System.Diagnostics;
using System.Collections.Generic;
using System.Linq;
using System.Threading.Tasks;
using System.Threading;
using System.IO;
using ASELib;

namespace MethylationAnalysis
{
	class MethylationPoint : IComparer<MethylationPoint>
	{
		public string compositeREF;
		public double mValue;
		public bool isSingle;

		public int Compare(MethylationPoint a, MethylationPoint b)
		{
			return xCompare(a, b);
		}

		static public int xCompare(MethylationPoint a, MethylationPoint b)
		{
			if (a.mValue > b.mValue) return 1;
			if (a.mValue < b.mValue) return -1;
			return 0;
		}
	}



	class Program
	{

		static ASETools.GeneLocationsByNameAndChromosome geneLocationInformation;

		static void ProcessFile(ASETools.Case case_)
		{
			// select metadata whether processing normal or tumor methylation data
			var methylation_filename = case_.tumor_methylation_filename;
			var methylation_file_id = case_.tumor_methylation_file_id;

			// Load MAF file for this case
			var mafLines = ASETools.MAFLine.ReadFile(case_.extracted_maf_lines_filename, case_.maf_file_id, false);

			// Initialize objects to hold methylation values for regions, chr and autosome
			var wholeAutosomeRegionalExpression = new ASETools.RegionalExpressionState();
			var perChromosomeRegionalExpressionState = new ASETools.RegionalExpressionState[ASETools.nHumanNuclearChromosomes];
			var allButThisChromosomeAutosomalRegionalExpressionState = new Dictionary<string, ASETools.RegionalExpressionState>();   
		   // per gene methylation values over incremental regions
		   Dictionary <string, ASETools.GeneExpression> geneExpressions = new Dictionary<string, ASETools.GeneExpression>();

			// Initialization of per chromosome methylation values
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

			// 2. for each gene in the genelocations, count the number of mutations and get the methylation values for this gene. These will be 
			// used later.
			// // for each gene, count the methylation averages as you get farther and farther away from the gene (inc, chr, autosome)
			foreach (var geneLocation in geneLocationInformation.genesByName.Values)
			{
				if (!geneExpressions.ContainsKey(geneLocation.hugoSymbol))
				{
					geneExpressions.Add(geneLocation.hugoSymbol, new ASETools.GeneExpression(geneLocation));
				}

				// get mutation counts from MAFLines
				var mutationCount = mafLines.Where(r => r.Hugo_Symbol == geneLocation.hugoSymbol).Count();
				geneExpressions[geneLocation.hugoSymbol].mutationCount += mutationCount;
				
			}

			// read in tumor methylation file 
			var annotations = ASELib.ASETools.AnnotationLine.ReadFile(methylation_filename, methylation_file_id, false);

			foreach (var annotation in annotations)
			{
				if (ASETools.isChromosomeAutosomal(annotation.compositeREF.Chromosome))
				{
					wholeAutosomeRegionalExpression.AddTumorExpression(annotation.M_Value, 0);

					foreach (var entry in allButThisChromosomeAutosomalRegionalExpressionState)
					{
						if (entry.Key != annotation.compositeREF.Chromosome)
						{
							entry.Value.AddTumorExpression(annotation.M_Value, 0);

						}
					}
				}

				int chromosomeId = ASETools.ChromosomeNameToIndex(annotation.compositeREF.Chromosome);
				if (chromosomeId != -1)
				{
					perChromosomeRegionalExpressionState[chromosomeId].AddTumorExpression(annotation.M_Value, 0);
				}
				
				foreach (var geneLocation in geneLocationInformation.genesByChromosome[ASETools.chromosomeNameToNonChrForm(annotation.compositeREF.Chromosome)])
				{
					geneExpressions[geneLocation.hugoSymbol].AddRegionalExpression(annotation.compositeREF.Start, annotation.M_Value, 0, true);
				}
				// TODO what about max, min, ave just for some idea of where we are?

			}

			// Save values

			var allExpressions = new List<ASETools.GeneExpression>();
			foreach (var expressionEntry in geneExpressions)
			{
				allExpressions.Add(expressionEntry.Value);
			}

			allExpressions.Sort(ASETools.GeneExpression.CompareByGeneName);

			string directory = ASETools.GetDirectoryPathFromFullyQualifiedFilename(case_.tumor_methylation_filename);
			string analysisId = ASETools.GetAnalysisIdFromPathname(case_.tumor_methylation_filename);

			var outputFilename = directory + analysisId + (ASETools.regionalMethylationExtension);
			var outputFile = ASETools.CreateStreamWriterWithRetry(outputFilename);

			Console.WriteLine("saving to file " + outputFilename);

			var columnSuffix = "tumor M-val";
			//ExpressionNearMutations.Program.writeColumnNames(outputFile, true, case_, columnSuffix);

			for (int i = 0; i < allExpressions.Count(); i++)
			{
				outputFile.Write(ASETools.ConvertToExcelString(allExpressions[i].geneLocationInfo.hugoSymbol) + "\t" + allExpressions[i].mutationCount);

				//ExpressionNearMutations.Program.writeRow(outputFile, allExpressions[i],
				//	wholeAutosomeRegionalExpression,
				//	allButThisChromosomeAutosomalRegionalExpressionState,
				//	perChromosomeRegionalExpressionState,
				//	false,
				//	1);

				outputFile.WriteLine();
			} // for each gene

			// Next analysis (loading in these numbers). Get the expression values in this part of the pipeline.
			// run mann whitney for methylation scores for some vs no mutation for genes that we know are important in that cancer cell type

		}


		static void ProcessCases(List<ASETools.Case> cases, string gene = null)
		{
			while (true)
			{
				ASETools.Case case_;
				lock (cases)
				{
					if (cases.Count() == 0)
					{
						//
						// No more work, we're done.
						//
						return;
					}

					Console.WriteLine(cases.Count() + " remaining...");

					case_ = cases[0];

					cases.RemoveAt(0);
				}


				// verify methylation file exists
				if (case_.tumor_methylation_filename == "" || case_.extracted_maf_lines_filename == "" || case_.tumor_rna_filename == "")
				{
					Console.WriteLine("No tumor methylation data for case " + case_.case_id + ". Skipping...");
				}
				else
				{
					// process only tumor files
					ProcessFile(case_);
				}
			}
		}


		//private static void MethylationMannWhitney(List<KeyValuePair<string, ASETools.GeneLocationInfo>> compositeREFs_)
		//{

		//	outputLines = new List<ASETools.OutputLine>();

		//	while (true)
		//	{
		//		KeyValuePair<string, ASETools.GeneLocationInfo> compositeREF;

		//		lock (compositeREFs_)
		//		{
		//			if (compositeREFs_.Count() == 0)
		//			{
		//				//
		//				// No more work, we're done.
		//				//
		//				return;
		//			}

		//			Console.WriteLine(compositeREFs_.Count() + " remaining...");
		//			compositeREF = compositeREFs_[0];

		//			compositeREFs_.RemoveAt(0);
		//		}
		//		Console.WriteLine(Thread.CurrentThread.ToString() + ": " + compositeREF.Key);
		//		Dictionary<string, int> geneCounts = new Dictionary<string, int>();

		//		try
		//		{
		//			// get per gene information (tuples of file id, mutation count)
		//			geneCounts = mutationCounts[compositeREF.Value.hugoSymbol].ToDictionary(x => x.Item1, x => x.Item2);
		//		}
		//		catch (Exception)
		//		{
		//			// Gene does not have any mutations, and should not be included in analysis
		//			continue;

		//		}

		//		var methylationPoints = casesMethylation[compositeREF.Key].Select(r =>
		//		{
		//			var mutations = 0;
		//			try
		//			{
		//				mutations = geneCounts[r.Item1];
		//			}
		//			catch (Exception)
		//			{
		//				// no op
		//			}
		//			// Tuple of M-value, mutations
		//			return new Tuple<double, int>(r.Item2, mutations);
		//		}).Where(r => r.Item2 < 2).Select(r => {    // filter out more than 1 mutation and format to comparative MethylationPoints
		//			var m = new MethylationPoint();
		//			m.mValue = ASETools.AnnotationLine.betaToM(r.Item1);
		//			m.isSingle = r.Item2 == 0;
		//			m.compositeREF = compositeREF.Key;
		//			return m;
		//		}).ToList();

		//		// run MannWhitney test for this methylation point
		//		ASETools.MannWhitney<MethylationPoint>.GetValue getValue = new ASETools.MannWhitney<MethylationPoint>.GetValue(m => m.mValue);
		//		ASETools.MannWhitney<MethylationPoint>.WhichGroup whichGroup = new ASETools.MannWhitney<MethylationPoint>.WhichGroup(m => m.isSingle);

		//		double nSingle = 0;
		//		double nMultiple = 0;
		//		double p;
		//		bool reversed;
		//		double U;
		//		double z;

		//		bool enoughData;

		//		p = ASETools.MannWhitney<MethylationPoint>.ComputeMannWhitney(methylationPoints, methylationPoints[0], whichGroup, getValue, out enoughData, out reversed, out nSingle, out nMultiple, out U, out z);
		//		if (!enoughData)
		//		{
		//			continue;
		//		}

		//		var outputLine = new ASETools.OutputLine();
		//		outputLine.line = compositeREF.Value.hugoSymbol + "\t" + compositeREF.Key + "\t" + nSingle + "\t" + nMultiple + "\t" + U + "\t" + z;
		//		outputLine.p = p;
		//		outputLines.Add(outputLine);

		//	}
		//}

		//static void MethylationMannWhitney()
		//{
		//	// total number of REFs for Bonferroni correction
		//	int compositeREFCount = compositeREFs.Keys.Count();

		//	MethylationMannWhitney(compositeREFs.ToList());


		//	// open file
		//	var output = new StreamWriter(@"\\msr-genomics-0\d$\gdc\methyl_temp\MannWhitney_with450.txt");

		//	// write header
		//	output.WriteLine("Gene\tCompositeREF\tnSingle\tnMultiple\tU\tz\tp\tp_Bonferroni\t");

		//	var compositeCount = outputLines.Count();

		//	foreach (var outputLine in outputLines)
		//	{
		//		output.WriteLine(outputLine.line + "\t" + outputLine.p + "\t" + outputLine.p * compositeCount);
		//	}

		//	output.Close();
		//}

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

			if (configuration.commandLineArgs.Count() == 0)
			{
				PrintUsageMessage();
				return;
			}

			// lines to print out
			var lines = new List<string>();

			// dictionary by case ids, then dictionary of composite ref and ratio score for that ref
			var ratios = new Dictionary<string, Dictionary<string, double>>();
			var composites = new Dictionary<string, ASETools.CompositeREF>();

			var cases = ASETools.Case.LoadCases(configuration.casesFilePathname);

			var s = cases.Where(r => r.Value.tumor_methylation_filename.Contains("450") && r.Value.normal_methylation_filename != "" && r.Value.disease() == "brca").Select(r => r.Value);
			var count = 0;
			Console.WriteLine("Processing " + s.Count() + " cases...");
			
			foreach (var case_ in s)
			{
				ratios.Add(case_.case_id, new Dictionary<string, double>());

				var tumorMethyls = ASETools.AnnotationLine.ReadFile(case_.tumor_methylation_filename, case_.tumor_methylation_file_id, false).Where(r => r.compositeREF.Chromosome == "chr13");
				var normalMethyls = ASETools.AnnotationLine.ReadFile(case_.normal_methylation_filename, case_.normal_methylation_file_id, false).Where(r => r.compositeREF.Chromosome == "chr13")
					.ToDictionary(x => x.compositeREF.Composite_Element_REF, x => x);

				foreach (var tumorM in tumorMethyls)
				{
					// If composite REF not yet in dictionary, add it
					if (!composites.ContainsKey(tumorM.compositeREF.Composite_Element_REF))
					{
						composites.Add(tumorM.compositeREF.Composite_Element_REF, tumorM.compositeREF);
					}

					if (!normalMethyls.ContainsKey(tumorM.compositeREF.Composite_Element_REF))
					{
						Console.WriteLine("normal composite REF for " + tumorM.compositeREF.Composite_Element_REF + " does not exist. Skipping...");
						continue;
					}

					var ratio = tumorM.Beta_Value / normalMethyls[tumorM.compositeREF.Composite_Element_REF].Beta_Value;
					ratios[case_.case_id].Add(tumorM.compositeREF.Composite_Element_REF, ratio);
				}
				count += 1;
				Console.WriteLine(count + " cases done...");
			} // foreach case

			// Print final output
			List<ASETools.CompositeREF> sortedREFs = new List<ASETools.CompositeREF>();

			var grouped = composites.GroupBy(r => r.Value.Chromosome);

			foreach (var group in grouped)
			{
				var ordered = group.Select(r => r.Value).ToList().OrderBy(r => r.Start);

				foreach (var item in ordered)
				{
					sortedREFs.Add(item);
				}
			}


			var outputFilename = @"C:\Users\t-almorr\temp\DNMT1brcatest.txt";
			var outputFile = ASETools.CreateStreamWriterWithRetry(outputFilename);

			// write header

			foreach (var composite in sortedREFs)
			{
				//outputFile.Write("\t" + String.Join(",",composite.Gene_Symbol.Distinct()) + "-" + composite.Chromosome + ":" + 
				//	composite.Start + "-" + composite.End);
				outputFile.Write("\t" + composite.Start);
			}
			outputFile.WriteLine();

			foreach (var ratio in ratios)
			{
				outputFile.Write(ratio.Key);
				foreach (var composite in sortedREFs)
				{
					if (!ratio.Value.ContainsKey(composite.Composite_Element_REF))
					{
						outputFile.Write("\tNaN");
					} else
					{
						outputFile.Write("\t" + ratio.Value[composite.Composite_Element_REF]);
					}
				}

				outputFile.WriteLine();
			}

			outputFile.Close();

			return;


			if (null == cases)
			{
				Console.WriteLine("You must have selected cases before you can analyze methylation data.");
				return;
			}

			var selectedCases = new List<ASETools.Case>();
			foreach (var arg in configuration.commandLineArgs)
			{
				if (!cases.ContainsKey(arg))
				{
					Console.WriteLine(arg + " is not a valid case ID. Skipping...");
					continue;
				} else if (!cases[arg].tumor_methylation_filename.Contains("450"))
				{
					Console.WriteLine(arg + " does not have 450k methylation data. Skipping...");
					continue;
				}
				selectedCases.Add(cases[arg]);
			}

			// Get information for current genome build
			geneLocationInformation = new ASETools.GeneLocationsByNameAndChromosome(ASETools.readKnownGeneFile(ASETools.ASEConfirguation.defaultGeneLocationInformation));

			var threads = new List<Thread>();
			for (int i = 0; i < Environment.ProcessorCount; i++)
			{
				threads.Add(new Thread(() => ProcessCases(selectedCases)));
			}

			threads.ForEach(th => th.Start());
			threads.ForEach(th => th.Join());

			//MethylationMannWhitney();

			return;

		}
	}
}
