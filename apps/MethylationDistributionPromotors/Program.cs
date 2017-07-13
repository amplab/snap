using System;
using System.Collections.Generic;
using System.Linq;
using System.Threading;
using ASELib;
using System.IO;

namespace MethylationDistributionPromotors
{
	class Program
	{

		public class MethylationPoint : IComparer<MethylationPoint>
		{
			public string identifier;
			public double adjustedBValue;
			public double bValue;
			public bool hasOneMutation;

			public MethylationPoint(string identifier_, double bValue_, bool hasOneMutation_)
			{
				bValue = bValue_;
				identifier = identifier_;
				hasOneMutation = hasOneMutation_;

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


		// Dictionary of hugoSymbol, and list of methylation points
		static Dictionary<string, List<MethylationPoint>> elements = new Dictionary<string, List<MethylationPoint>>();

		static void ProcessCasesOnlyPromotors(List<ASETools.Case> cases, List<string> selectedHugoSymbols)
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

					case_ = cases[0];

					cases.RemoveAt(0);
				}

				// Read in methylation values and ASE values
				if (case_.tumor_allele_specific_gene_expression_filename == "" || case_.tumor_methylation_filename == "")
				{
					// no data. continue
					continue;
				}


				var methylationData = ASETools.AnnotationLine.ReadFile(case_.tumor_methylation_filename, case_.tumor_methylation_file_id, false);

				foreach (var hugoSymbol in selectedHugoSymbols)
				{
					// Check if the mutation counts for this gene and case id exist
					Dictionary<string, int> mutationCountsForThisGene;
					int mutationCount;
					if (mutationCounts.TryGetValue(hugoSymbol, out mutationCountsForThisGene))
					{
						if (!mutationCountsForThisGene.TryGetValue(case_.case_id, out mutationCount))
						{
							continue;
						}

					}
					else
					{
						continue;
					}

					var methylationForGene = methylationData.Where(r => 
						r.compositeREF.Gene_Symbol.Contains(ASETools.ConvertToNonExcelString(hugoSymbol)));

					// filter by promotors for this gene
					if (!allSites)
					{
						methylationForGene = methylationForGene.Where(r =>
						{
							var merged = r.compositeREF.Gene_Symbol.Zip(r.compositeREF.Position_to_TSS, (first, second) =>
								new Tuple<string, int>(first, second))
									.Where(t => t.Item1 == ASETools.ConvertToNonExcelString(hugoSymbol)).Select(f => f.Item2);
							return Math.Abs(merged.Average()) <= 2000;
						});
					}

					if (methylationForGene.Count() == 0)
					{
						continue;
					}

					var combinedScore = methylationForGene.Select(r => r.Beta_Value).Average();

					// Add to elements
					MethylationPoint x = new MethylationPoint(hugoSymbol, combinedScore, mutationCount == 1);

					lock (elements)
					{
						// If gene does not exist 
						List<MethylationPoint> element;
						if (!elements.TryGetValue(hugoSymbol, out element))
						{
							elements.Add(hugoSymbol, new List<MethylationPoint>());
						}
						elements[hugoSymbol].Add(x);
					}


				}

			}
		}

		static StreamWriter panCancerOutputFile;

		// stores distance for lowest pvalue
		static Dictionary<string, Tuple<string, double>> bestPValues = new Dictionary<string, Tuple<string, double>>();

		// Dictionary of hugoSymbol, then dictionary of (case_id, mutationCount)
		static Dictionary<string, Dictionary<string, int>> mutationCounts = new Dictionary<string, Dictionary<string, int>>();

		// Dictionary of hugoSymbol, then dictionare of (case_id, mutationCount)
		static Dictionary<string, Dictionary<string, double[]>> methylationValues = new Dictionary<string, Dictionary<string, double[]>>();

		static void readFinalValues(string filename)
		{

			var modValue = 11;

			// TODO move to bonferroni corrected values once available
			List<string[]> lines = ASETools.ReadAllLinesWithRetry(filename).Select(r => r.Split('\t')).ToList();

			var headers = lines[0];

			for (var i = 1; i < lines.Count(); i++)
			{
				var hugoSymbol = lines[i][0];

				// get p-values for 1 vs not 1

				var pvalues = lines[i].Select((value, index) => new Tuple<string, int>(value, index))
					.Where(r => (r.Item2 - 2) % modValue == 0)
					.Where(r => r.Item1 != "*").Select(r => new Tuple<double, int>(Convert.ToDouble(r.Item1), r.Item2));

				if (pvalues.Count() == 0)
				{
					continue;
				}

				var min = pvalues.Min(obj => obj.Item1);

				var bestPValue = pvalues.Where(r => r.Item1 == min).First();

				var label = headers[bestPValue.Item2];
				var distance = label.Substring(0, label.IndexOf(" 1 vs. not 1"));
				bestPValues.Add(hugoSymbol, new Tuple<string, double>(distance, min));
			}
		}

		static void loadMutationCounts(List<ASETools.Case> cases)
		{
			foreach (var case_ in cases)
			{
				if (case_.tumor_allele_specific_gene_expression_filename == "")
				{
					continue;
				}

				var aseFile = ASETools.RegionalSignalFile.ReadFile(case_.tumor_allele_specific_gene_expression_filename);

				foreach (var hugoSymbol in bestPValues.Keys)
				{
					double[] hugoData;

					if (!aseFile.Item1.TryGetValue(ASETools.ConvertToExcelString(hugoSymbol), out hugoData))
					{
						continue;
					}

					// get number of non-silent mutations
					var mutationCount = aseFile.Item3[ASETools.ConvertToExcelString(hugoSymbol)];
					if (!mutationCounts.ContainsKey(hugoSymbol))
					{
						mutationCounts.Add(hugoSymbol, new Dictionary<string, int>());
					}

					mutationCounts[hugoSymbol].Add(case_.case_id, mutationCount);
				}
			}
		}

		static bool allSites = true;

		static void Main(string[] args)
		{
			var configuration = ASETools.ASEConfirguation.loadFromFile(args);
			var cases = ASETools.Case.LoadCases(configuration.casesFilePathname);

			if (configuration.commandLineArgs.Count() > 0 && configuration.commandLineArgs[0] == "-p")
			{
				allSites = false;
			}

			string baseFileName = configuration.finalResultsDirectory +
				(allSites ? "methylationResults.all.txt" : "methylationResults.promotors.txt");

			panCancerOutputFile = ASETools.CreateStreamWriterWithRetry(baseFileName);

			var pValueFile = configuration.finalResultsDirectory + "AlleleSpecificExpressionDistributionByMutationCount.txt";
			readFinalValues(pValueFile);

			// Select only hugo symbols with low pvalues
			List<string> selectedHugoSymbols = bestPValues.Where(r => r.Value.Item2 < 1.0E-6).Select(r => r.Key).ToList();

			// load the mutations for each file
			loadMutationCounts(cases.Select(r => r.Value).ToList());

			ASETools.MannWhitney<MethylationPoint>.WhichGroup whichGroup = new ASETools.MannWhitney<MethylationPoint>.WhichGroup(m => m.hasOneMutation);
			ASETools.MannWhitney<MethylationPoint>.GetValue getValue = new ASETools.MannWhitney<MethylationPoint>.GetValue(x => x.adjustedBValue);
			bool twoTailed = true;

			var threads = new List<Thread>();
			var selectedCases = cases.Select(r => r.Value).ToList();

			for (int i = 0; i < Environment.ProcessorCount; i++)
			{
				threads.Add(new Thread(() => ProcessCasesOnlyPromotors(selectedCases, selectedHugoSymbols)));
			}

			threads.ForEach(th => th.Start());
			threads.ForEach(th => th.Join());

			// write header
			panCancerOutputFile.WriteLine("Gene Symbol\tbest pval\tbest location\t0 kb");

			foreach (var hugoSymbol in selectedHugoSymbols)
			{
				// dictionary of (case, mutationcount)
				Dictionary<string, int> mutationCountsForThisGene;
				if (!mutationCounts.TryGetValue(ASETools.ConvertToExcelString(hugoSymbol), out mutationCountsForThisGene))
				{
					continue;
				}

				var bestPValue = bestPValues[hugoSymbol];


				var hugoLine = hugoSymbol + "\t" + bestPValue.Item2 + "\t" + bestPValue.Item1;

				// If elements has no items, skip this hugo symbol
				List<MethylationPoint> element;
				if (!elements.TryGetValue(hugoSymbol, out element))
				{
					hugoLine += "\t*";
					panCancerOutputFile.WriteLine(hugoLine);
					continue;
				}

				bool reversed;
				bool enoughData;
				double nFirstGroup;
				double nSecondGroup;
				double U;
				double z;

				if (element.Count() > 0)
				{
					var p = ASETools.MannWhitney<MethylationPoint>.ComputeMannWhitney(element,
						element[0], whichGroup, getValue, out enoughData, out reversed,
						out nFirstGroup, out nSecondGroup, out U, out z, twoTailed, 20);
					if (!enoughData)
					{
						hugoLine += "\t*";
					}
					else
					{
						hugoLine += "\t" + p;
					}
				}
				else
				{
					hugoLine += "\t*";
				}


				panCancerOutputFile.WriteLine(hugoLine);
				panCancerOutputFile.Flush();

			} // foreach hugo symbol 

			panCancerOutputFile.Close();

		}
	}
}
