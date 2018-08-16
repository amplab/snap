using System;
using System.Diagnostics;
using System.Collections.Generic;
using System.Linq;
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

		// Distance from gene to process methylation values
		static int maxDistance = 100000;
		
		// keeps track of which diseases to evaluate
		static List<string> diseases = new List<string>();

		// Dictionary storing methylation REF values
		// key: Composite REF Id
		// Value: List of case ids and beta values for this composite Ref
		static Dictionary<string, List<Tuple<string, Double>>> casesMethylation = new Dictionary<string, List<Tuple<string, Double>>>();

		// Dictionary storing Composite REF information. 
		// Key: Composite REF Id
		// Value: Composite REF info
		static Dictionary<string, ASETools.CompositeREF> compositeREFs 
			= new Dictionary<string, ASETools.CompositeREF>();

		// Dictionary storing case data.
		// Key: case Id
		// Value: Tuple of disease, cell type (tumor or normal)
		static Dictionary<string, Tuple<string, string>> fileInfo = new Dictionary<string, Tuple<string, string>>();

		// dictionary storing mutation counts for each gene
		// Key: Hugo Symbol
		// Value: List of (case id, mutation count) tuples for this hugo symbol
		static Dictionary<string, List<Tuple<string, int>>> mutationCounts = new Dictionary<string, List<Tuple<string, int>>>();

		static void ProcessFile(ASETools.Case case_, bool isTumor, List<ASETools.GeneLocationInfo> genes)
		{
			// select metadata whether processing normal or tumor methylation data
			var label = isTumor ? "Tumor" : "Normal";
			var methylation_filename = isTumor ? case_.tumor_methylation_filename : case_.normal_methylation_filename;

			// read in tumor methylation file 
			var annotations = ASELib.ASETools.MethylationAnnotationLine.ReadFile(methylation_filename);

			fileInfo.Add(case_.case_id, new Tuple<string, string>(case_.disease(), label));

			foreach (var annotation in annotations)
			{
				// check if composite REF is close to gene
				var genesInRange = genes.Where(gene => gene.chromosome == annotation.compositeREF.Chromosome && Math.Abs(annotation.compositeREF.Start - gene.minLocus) < maxDistance ||
				Math.Abs(annotation.compositeREF.Start - gene.maxLocus) < maxDistance);

				if (genesInRange.Count() == 0)
				{
					continue;
				}

				lock (compositeREFs)
				{
					if (!compositeREFs.ContainsKey(annotation.compositeREF.Composite_Element_REF))
					{
						compositeREFs.Add(annotation.compositeREF.Composite_Element_REF, annotation.compositeREF);
					}
				}

				Tuple<String, Double> betaValueForCase = new Tuple<String, Double>(case_.case_id, annotation.Beta_Value);

				// add beta value 
				lock (casesMethylation)
				{
					if (!casesMethylation.ContainsKey(annotation.compositeREF.Composite_Element_REF))
					{
						// initialize file id to methylation
						casesMethylation.Add(annotation.compositeREF.Composite_Element_REF, new List<Tuple<string, double>>());
					}

					// append annotation to existing REF
					casesMethylation[annotation.compositeREF.Composite_Element_REF].Add(betaValueForCase);
				}
			}

			// read in maf file for counting mutations
			var mafLines = ASETools.MAFLine.ReadFile(case_.extracted_maf_lines_filename, case_.maf_file_id, false);

			var mutationDictionary = ASETools.MAFLine.GetMutationCounts(mafLines);

			foreach (var gene in genes)
			{
				int mutationCount;
				if (!mutationDictionary.TryGetValue(gene.hugoSymbol, out mutationCount))
				{
					mutationCount = 0;
				}

				lock (mutationCounts)
				{
					if (!mutationCounts.ContainsKey(gene.hugoSymbol))
					{
						mutationCounts.Add(gene.hugoSymbol, new List<Tuple<string, int>>());
					}

					mutationCounts[gene.hugoSymbol].Add(new Tuple<string, int>(case_.case_id, mutationCount));
				}
			}

		}

		static void ProcessCases(List<ASETools.Case> cases, List<ASETools.GeneLocationInfo> genes = null)
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

					if (cases.Count() % 100 == 0)
					{
						Console.WriteLine(cases.Count() + " remaining...");
					}

					case_ = cases[0];

					cases.RemoveAt(0);
				}


				// verify methylation file exists
				if (!case_.tumor_methylation_filename.Contains("HumanMethylation450") || case_.extracted_maf_lines_filename == "")
				{
					Console.WriteLine("No tumor methylation data for case " + case_.case_id + ". Skipping...");
					continue;
				}
				// process only tumor files
				ProcessFile(case_, true, genes);

			}
		}

		private static void MethylationMannWhitney(Dictionary<string, List<Tuple<string, Double>>> refsToProcess,
			List<ASETools.GeneLocationInfo> genes, string filename)
		{
			// Initialize file and write header
			var output = ASETools.CreateStreamWriterWithRetry(filename);
			WriteMannWhitneyHeader(output);

			var bonferroniCorrection = diseases.Count() * compositeREFs.Count();

			// for each gene in analysis
			foreach (var gene in genes)
			{
				// Dictionary of case id, mutation count
				Dictionary<string, int> geneCounts = new Dictionary<string, int>();

				try
				{
					// get per gene information (tuples of file id, mutation count)
					geneCounts = mutationCounts[gene.hugoSymbol].ToDictionary(x => x.Item1, x => x.Item2);
				}
				catch (Exception)
				{
					// Gene does not have any mutations, and should not be included in analysis
					continue;
				}

				// get all methylation points close within the range of this gene
				var compositeREFsForGene = compositeREFs.Where(r => Math.Abs(r.Value.Start - gene.minLocus) < maxDistance ||
						Math.Abs(r.Value.Start - gene.maxLocus) < maxDistance);
				var compositeREFIds = compositeREFsForGene.Select(r => r.Key);

				// get methylation points for all points close to this gene
				var methylationPoints = casesMethylation.Where(r => compositeREFIds.Contains(r.Key));

				// for each composite REF for this gene, run Mann Whitney and save values
				foreach (var methylationPoint in methylationPoints)
				{
					// get distance from gene to this composite REF
					var distanceFromGene = Math.Min(Math.Abs(compositeREFs[methylationPoint.Key].Start - gene.minLocus),
					Math.Abs(compositeREFs[methylationPoint.Key].Start - gene.maxLocus));

					var outputLine = new ASETools.OutputLine();
					outputLine.line = gene.hugoSymbol + "\t" + compositeREFs[methylationPoint.Key].Start + "\t" + 
						distanceFromGene + "\t" + methylationPoint.Key;

					// keeps track of number of diseases that got signifiant values
					var significantDiseases = 0;

					foreach (var disease in diseases)
					{
						var casesForDisease = fileInfo.Where(r => r.Value.Item1 == disease).Select(r => r.Key);

						// get cases that have this disease
						var pointsForDisease = methylationPoint.Value.Where(r => casesForDisease.Contains(r.Item1));

						var points = pointsForDisease.Select(r => {
								var mutations = geneCounts[r.Item1];
								return new Tuple<double, int>(r.Item2, mutations);
							}).Select(r => {   
							var m = new MethylationPoint();
							m.mValue = ASETools.MethylationAnnotationLine.betaToM(r.Item1);
							m.isSingle = r.Item2 == 1;
							m.compositeREF = methylationPoint.Key;
							return m;
							}).ToList();

						// run MannWhitney test for this methylation point
						ASETools.MannWhitney<MethylationPoint>.GetValue getValue = new ASETools.MannWhitney<MethylationPoint>.GetValue(m => m.mValue);
						ASETools.MannWhitney<MethylationPoint>.WhichGroup whichGroup = new ASETools.MannWhitney<MethylationPoint>.WhichGroup(m => m.isSingle);

						double nSingle = 0;
						double nMultiple = 0;
						double p;
						bool reversed;
						double U;
						double z;

						bool enoughData;

						var singleValues = points.Where(r => r.isSingle).Select(r => r.mValue).ToList();
						var multipleValues = points.Where(r => !r.isSingle).Select(r => r.mValue).ToList();

						var singleMean = ASETools.MeanOfList(singleValues);
						var multipleMean = ASETools.MeanOfList(multipleValues);

						var singleStdDev = ASETools.StandardDeviationOfList(singleValues.ToList());
						var multipleStdDev = ASETools.StandardDeviationOfList(multipleValues.ToList());

						p = ASETools.MannWhitney<MethylationPoint>.ComputeMannWhitney(points, points[0], whichGroup, getValue, out enoughData, out reversed, out nSingle, out nMultiple, out U, out z);

						// if this is a significant point for this gene, save it to file
						if (enoughData)
						{
							outputLine.line += "\t" + nSingle + "\t" + nMultiple + "\t"
								+ singleMean + "\t" + multipleMean
								+ "\t" + singleStdDev + "\t" + multipleStdDev 
								+ U + "\t" + z + "\t" 
								+ "\t" + p + "\t" + p * bonferroniCorrection;

							// if this disease is significant, record it
							if (p * bonferroniCorrection < 0.01)
								significantDiseases += 1;
						} // if enough data
						else
						{
							outputLine.line += "\t" + nSingle + "\t" + nMultiple + String.Concat(Enumerable.Repeat<string>("\tNaN", 8));
						}
					} // foreach disease

					if (significantDiseases > 0)
					{
						outputLine.line += "\t" + significantDiseases;
						output.WriteLine(outputLine.line);
					}

				} // foreach methylation point

			} // foreach gene

			output.Close();

		}

		// Writes out Mann Whitney results for methylation
		static void WriteMannWhitneyHeader(StreamWriter output)
		{

			// write header
			output.Write("HugoSymbol\tStart\tDistanceToGene\tCompositeREF\t");

			foreach (var disease in diseases)
			{
				output.Write(disease + "_nSingle\tnMultiple\tsingleMean\tmultipleMean\tsingleStdDev\tmultipleStdDev\t" 
					+ disease + "_U\tz\tp\tp_Bonferroni\t");
			}

			output.Write("SignificantDiseases");

			output.WriteLine();
		}

		static void printUsageMessage()
		{
			Console.WriteLine("MethylationAnalysis distance");
			Console.WriteLine("-distance specifies base pair range from gene to process methylation values");
		}

		// Program that computes differences in distributions for a list of genes
		// Runs Mann Whitney per disease on methylation points
		static void Main(string[] args)
		{
			var timer = new Stopwatch();
			timer.Start();

			var configuration = ASETools.Configuration.loadFromFile(args);
			var cases = ASETools.Case.LoadCases(configuration.casesFilePathname);

			if (null == cases)
			{
				Console.WriteLine("You must have selected cases before you can analyze methylation data.");
				return;
			}

			var selectedCases = cases.Select(kv => kv.Value).ToList();

			Console.WriteLine("Processing " + selectedCases.Count() + " cases");

			// get diseases from selected cases to process
			diseases = selectedCases.Select(r => r.disease()).Distinct().ToList();

			if (configuration.commandLineArgs.Count() < 1)
			{
				printUsageMessage();
				return;
			}

			// Get genes to process
			var knownGenes = ASETools.readKnownGeneFile(configuration.geneLocationInformationFilename);
			var knownGeneNames = knownGenes.Select(r => r.Value.hugoSymbol);

			var selectedGenes = new List<ASETools.GeneLocationInfo>();

			foreach (var arg in configuration.commandLineArgs)
			{
				if (knownGeneNames.Contains(arg))
				{
					selectedGenes.Add(knownGenes[arg]);
				}
			}

			// Get base pair distance from genes to process methylation values
			maxDistance = Convert.ToInt32(configuration.commandLineArgs[0]);

			// If no genes were selected, run all genes that are significant
			if (selectedGenes.Count() == 0)
			{
				var significantResults = ASETools.CreateStreamReaderWithRetry(configuration.finalResultsDirectory + ASETools.AlleleSpecificExpressionDistributionByMutationCountFilenameBase + "_bonferroni.txt");
				significantResults.ReadLine();

				string line;
				while ((line = significantResults.ReadLine()) != null)
				{
					var split = line.Split('\t');

					// for now, filter out significant genes that have low p-values at gene
					var zeroDistance = split[2] == "0Kbp";
					var geneName = ASETools.ConvertToNonExcelString(split[0]);
					
					if (split[3].ToLower() != "true" && split[3].ToLower() != "false")
					{
						continue;
					}
					var isSignificant = split[3] == "true";

					ASETools.GeneLocationInfo geneInformation;
					if (knownGenes.TryGetValue(geneName, out geneInformation) && isSignificant && !zeroDistance)
					{
						selectedGenes.Add(geneInformation);
					}
				}
			}

			// Process cases 
			var threads = new List<Thread>();
			ASETools.Shuffle<ASETools.Case>(selectedCases);
			for (int i = 0; i < Environment.ProcessorCount; i++)
			{
				threads.Add(new Thread(() => ProcessCases(selectedCases, selectedGenes)));
			}

			threads.ForEach(th => th.Start());
			threads.ForEach(th => th.Join());

			// Run Mann Whitney
			var mannWhitneyOutputFilename = configuration.finalResultsDirectory + "MannWhitneyMethylation_brca_sm.txt";

			MethylationMannWhitney(casesMethylation, selectedGenes, mannWhitneyOutputFilename);
		}
	}
}
