using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading;
using ASELib;
using System.IO;

namespace MethylationDistribution
{
	class Program
	{

		// keeps track of the tumor and normal methylation distributions
		static Dictionary<double, int> tumorMethylationDistribution = new Dictionary<double, int>();
		static Dictionary<double, int> normalMethylationDistribution = new Dictionary<double, int>();

		static Dictionary<string, List<ASETools.MethylationPoint[]>> elements = new Dictionary<string, List<ASETools.MethylationPoint[]>>();

		static void ProcessCases(List<ASETools.Case> cases, List<string> selectedHugoSymbols)
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
				if (case_.tumor_allele_specific_gene_expression_filename == "" || case_.tumor_regional_methylation_filename == "")
				{
					// no data. continue
					continue;
				}

				var regionalMethylationData = ASETools.RegionalSignalFile.ReadFile(case_.tumor_regional_methylation_filename);


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

					// TODO check for value
					double[] methylationForGene;
					if (!regionalMethylationData.Item1.TryGetValue(hugoSymbol, out methylationForGene))
					{
						continue;
					}

					// Add to elements
					ASETools.MethylationPoint[] array = methylationForGene.Select(r => new ASETools.MethylationPoint(hugoSymbol, r, mutationCount == 1))
						.ToArray();

					lock (elements)
					{
						// If gene does not exist 
						List<ASETools.MethylationPoint[]> element;
						if (!elements.TryGetValue(hugoSymbol, out element))
						{
							elements.Add(hugoSymbol, new List<ASETools.MethylationPoint[]>());
						}
						elements[hugoSymbol].Add(array);
					}
				}

			}
		}


		static StreamWriter panCancerOutputFile;

		// stores distance for lowest pvalue and first available ASE Value
		static Dictionary<string, Tuple<string, double, double>> bestPValues = new Dictionary<string, Tuple<string, double, double>>();

		// Dictionary of hugoSymbol, then dictionary of (case_id, mutationCount)
		static Dictionary<string, Dictionary<string, int>> mutationCounts = new Dictionary<string, Dictionary<string, int>>();

		// Dictionary of hugoSymbol, then dictionare of (case_id, mutationCount)
		static Dictionary<string, Dictionary<string, double[]>> methylationValues = new Dictionary<string, Dictionary<string, double[]>>();

		static void readFinalValues(string filename)
		{

			var pModValue = 11;
			var aseModValue = 7; // TODO

			// TODO move to bonferroni corrected values once available
			List<string[]> lines = ASETools.ReadAllLinesWithRetry(filename).Select(r => r.Split('\t')).ToList();

			var headers = lines[0];

			// Iterate through each line for a gene
			for (var i = 1; i < lines.Count(); i++)
			{
				var hugoSymbol = lines[i][0];

				// get p-values for 1 vs not 1
				var pvalues = lines[i].Select((value, index) => new Tuple<string, int>(value, index))
					.Where(r => (r.Item2 - 2) % pModValue == 0)
					.Where(r => r.Item1 != "*").Select(r => new Tuple<double, int>(Convert.ToDouble(r.Item1), r.Item2));

				// get first ASE value
				var firstAseValue = lines[i].Select((value, index) => new Tuple<string, int>(value, index))
					.Where(r => (r.Item2) % aseModValue == 0)
					.Where(r => r.Item1 != "*" || r.Item1.Trim() != "")
					.Select(r =>Convert.ToDouble(r.Item1)).First();

				if (pvalues.Count() == 0)
				{
					continue;
				}

				var min = pvalues.Min(obj => obj.Item1);

				var bestPValue = pvalues.Where(r => r.Item1 == min).First();

				var label = headers[bestPValue.Item2];
				var distance = label.Substring(0, label.IndexOf(" 1 vs. not 1"));
				bestPValues.Add(hugoSymbol, new Tuple<string, double, double>(distance, min, firstAseValue));
			}
		}

		static string loadMutationCounts(List<ASETools.Case> cases)
		{
			var header = "";
			foreach (var case_ in cases)
			{
				if (case_.tumor_allele_specific_gene_expression_filename == "")
				{
					continue;
				}

				var aseFile = ASETools.RegionalSignalFile.ReadFile(case_.tumor_allele_specific_gene_expression_filename);

				if (header.Length == 0)
				{
					header = String.Join("\t", aseFile.Item2);
				}

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
			return header;
		}

		static void Main(string[] args)
		{
			var configuration = ASETools.ASEConfirguation.loadFromFile(args);
			var cases = ASETools.Case.LoadCases(configuration.casesFilePathname);

			// Case 1: for any group of cases, compute the pvalues for adjusted beta value.
			// Run Mann Whitney for test for ASM values from adjusted beta value where the groups are 
			// ase vs no ase. (per gene) The purpose of this is to see if ASM significantly explains ASE for any locations.

			var regionsToProcess = 66;


			string baseFileName = configuration.finalResultsDirectory + "methylationResults.allSites.txt";

			panCancerOutputFile = ASETools.CreateStreamWriterWithRetry(baseFileName);

			var pValueFile = configuration.finalResultsDirectory + "AlleleSpecificExpressionDistributionByMutationCount.txt";
			readFinalValues(pValueFile);

			// Select only hugo symbols with low pvalues
			var selectedHugoSymbols = bestPValues.Where(r => r.Value.Item2 < 1.0E-6).Select(r => r.Key);

			// load the mutations for each file
			var header = loadMutationCounts(cases.Select(r => r.Value).ToList());

			ASETools.MannWhitney<ASETools.MethylationPoint>.WhichGroup whichGroup = new ASETools.MannWhitney<ASETools.MethylationPoint>.WhichGroup(m => m.hasOneMutation);
			ASETools.MannWhitney<ASETools.MethylationPoint>.GetValue getValue = new ASETools.MannWhitney<ASETools.MethylationPoint>.GetValue(x => x.adjustedBValue);
			bool twoTailed = true;

			// Preprocess cases by putting all required values in methylationValues
			var threads = new List<Thread>();
			var selectedCases = cases.Select(r => r.Value).ToList();

			for (int i = 0; i < Environment.ProcessorCount; i++)
			{
				threads.Add(new Thread(() => ProcessCases(selectedCases, selectedHugoSymbols.ToList())));
			}

			threads.ForEach(th => th.Start());
			threads.ForEach(th => th.Join());

			// write file header
			panCancerOutputFile.WriteLine("Gene Symbol\tbest pval\tbest location\tFirst ASE\t" + header);

			foreach (var hugoSymbol in selectedHugoSymbols)
			{
				Dictionary<string, int> mutationCountsForThisGene;
				if (!mutationCounts.TryGetValue(hugoSymbol, out mutationCountsForThisGene))
				{
					continue;
				}

				var bestPValue = bestPValues[hugoSymbol];

				// Write: best pvalue, location of best pvalue, ASE location at first spot
				var hugoLine = hugoSymbol + "\t" + bestPValue.Item2 + "\t" + bestPValue.Item1 + "\t" + bestPValue.Item3;

				// If elements has no items, skip this hugo symbol
				List<ASETools.MethylationPoint[]> element;
				if (!elements.TryGetValue(hugoSymbol, out element))
				{
					hugoLine += string.Concat(Enumerable.Repeat("\t*", regionsToProcess));
					panCancerOutputFile.WriteLine(hugoLine);
					continue;
				}

				// Iterate through all regions
				for (var i = 0; i < regionsToProcess; i++)
				{
					bool reversed;
					bool enoughData;
					double nFirstGroup;
					double nSecondGroup;
					double U;
					double z;

					var forRegion = element.Select(r => r[i]).Where(r => r.mValue > Double.NegativeInfinity).ToList(); // Filter by neg infinity
					if (forRegion.Count() > 0)
					{
						var p = ASETools.MannWhitney<ASETools.MethylationPoint>.ComputeMannWhitney(forRegion, 
							forRegion[0], whichGroup, getValue, out enoughData, out reversed, 
							out nFirstGroup, out nSecondGroup, out U, out z, twoTailed, 1);
						if (!enoughData)
						{
							hugoLine += "\t*";
						}
						else
						{
							hugoLine += "\t" + p;
						}
					}
					else {
						hugoLine += "\t*";
					}

				} // foreach region

				panCancerOutputFile.WriteLine(hugoLine);
				panCancerOutputFile.Flush();

			} // foreach hugo symbol 

			panCancerOutputFile.Close();


			// Case 2: Run t test, trying to assign distributions for methylation points. This is to only
			// look at 1 group (ie. ase) Thie purpose of this is to see what percentage of ASE significant
			// results can be explained by ASE. This should be written generically enough so it can
			// also be used to test no ASE for 0 mutations, testing what percentage of these cases fall into the
			// 'full methylation' category

		}
	}
}
