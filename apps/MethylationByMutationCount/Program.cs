using System;
using System.Collections.Generic;
using System.Linq;
using System.Diagnostics;
using ASELib;
using System.Threading;

// Case 2: Run t test, trying to assign distributions for methylation points. This is to only
// look at 1 group (ie. ase) Thie purpose of this is to see what percentage of ASE significant
// results can be explained by ASE. This should be written generically enough so it can
// also be used to test no ASE for 0 mutations, testing what percentage of these cases fall into the
// 'full methylation' category

namespace MethylationByMutationCount
{
	class Program
	{

		static void PrintUsageMessage()
		{
			Console.WriteLine("usage: MethylationByMutationCount hugoSymbols");
			Console.WriteLine("hugo symbols are separated by space");
		}

		// for each case id, map of composite refs and values
		static Dictionary<string, Dictionary<string, Dictionary<string, double>>> methylationValues = 
			new Dictionary<string, Dictionary<string, Dictionary<string, double>>>();

		static Tuple<Dictionary<string, Dictionary<string, int>>, Dictionary<string, Dictionary<string, double>>> aseValues =
			new Tuple<Dictionary<string, Dictionary<string, int>>, Dictionary<string, Dictionary<string, double>>>(new Dictionary<string, Dictionary<string, int>>(), new Dictionary<string, Dictionary<string, double>>());


		//// first is mutations, second is values
		//static Tuple<Dictionary<string, Dictionary<string, int>>, Dictionary<string, Dictionary<string, double>>> expressionValues =
		//	new Tuple<Dictionary<string, Dictionary<string, int>>, Dictionary<string, Dictionary<string, double>>>(new Dictionary<string, Dictionary<string, int>>(), new Dictionary<string, Dictionary<string, double>>());
		// dictionary of hugosymbol, dictionary of case id, FPKM
		static Dictionary<string, Dictionary<string, double>> expressionValues =
			new Dictionary<string, Dictionary<string, double>>();

		static void loadMethylation(List<ASETools.Case> cases, List<string> hugoSymbols)
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

				if (!case_.tumor_methylation_filename.Contains("HumanMethylation450") || case_.annotated_selected_variants_filename == "")
				{
					// no data. continue. We use the annotated selected variants to evaluate loss of heterozygosity at each gene
					continue;
				}

				var methylationData = ASETools.AnnotationLine.ReadFile(case_.tumor_methylation_filename, case_.tumor_methylation_file_id, false);

				foreach (var hugoSymbol in hugoSymbols)
				{
					lock (methylationValues)
					{
						Dictionary<string, Dictionary<string, double>> values;
						if (!methylationValues.TryGetValue(hugoSymbol, out values))
						{
							methylationValues.Add(hugoSymbol, new Dictionary<string, Dictionary<string, double>>());
						}
						methylationValues[hugoSymbol].Add(case_.case_id, new Dictionary<string, double>());
					}

					var methylationGenePoints = methylationData.Where(r => r.compositeREF.Gene_Symbol.Contains(ASETools.ConvertToNonExcelString(hugoSymbol)))
					.Where(r =>
					{
						var merged = r.compositeREF.Gene_Symbol.Zip(r.compositeREF.Position_to_TSS, (first, second) =>
							new Tuple<string, int>(first, second))
								.Where(t => t.Item1 == ASETools.ConvertToNonExcelString(hugoSymbol)).Select(f => f.Item2);
						return Math.Abs(merged.Average()) <= 2000;
					});

					foreach (var m in methylationGenePoints)
					{
						methylationValues[hugoSymbol][case_.case_id].Add(m.compositeREF.Composite_Element_REF, m.M_Value);
					}
				}

			}
		}


		static void loadRegionalSignal(List<ASETools.Case> cases, 
			List<string> hugoSymbols,
			Tuple<Dictionary<string, Dictionary<string, int>>, 
				Dictionary<string, Dictionary<string, double>>> values,
			Boolean forAlleleSpecificExpression)
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

				var filename = forAlleleSpecificExpression ? case_.tumor_allele_specific_gene_expression_filename : case_.gene_expression_filename;

				if (filename == "")
				{
					// no data. continue
					continue;
				}

				var regionalSignals = ASETools.RegionalSignalFile.ReadFile(filename, true, true);

				foreach (var hugoSymbol in hugoSymbols)
				{
					lock (values)
					{
						// make a spot for this gene in the dictionary of mutation counts
						Dictionary<string, int> mutationValue;
						if (!values.Item1.TryGetValue(hugoSymbol, out mutationValue))
						{
							values.Item1.Add(hugoSymbol, new Dictionary<string, int>());
						}

						// make a spot for this gene in the signal dictionary
						Dictionary<string, double> signalValue;
						if (!values.Item2.TryGetValue(hugoSymbol, out signalValue))
						{
							values.Item2.Add(hugoSymbol, new Dictionary<string, double>());
						}
					}
					double[] hugoData;
					if (regionalSignals.Item1.TryGetValue(ASETools.ConvertToNonExcelString(hugoSymbol), out hugoData))
					{
						lock (values.Item2)
						{
							// add regional signal at gene
							values.Item2[hugoSymbol].Add(case_.case_id, hugoData[1]);
						}
						lock (values.Item1)
						{
							// add mutation count for gene
							values.Item1[hugoSymbol].Add(case_.case_id, Convert.ToInt32(hugoData[0]));
						}
					}
				}
			}
		}


		static void loadFPKM(List<ASETools.Case> cases,
				List<string> hugoSymbols, Dictionary<string, Dictionary<string, double>> values)
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

				var filename = case_.tumor_fpkm_filename;

				if (filename == "")
				{
					// no data. continue
					continue;
				}

				var regionalSignals = ASETools.FPKMFile.ReadFile(filename);

				foreach (var hugoSymbol in hugoSymbols)
				{
					lock (values)
					{

						Dictionary<string, double> fpkm;
						if (!values.TryGetValue(hugoSymbol, out fpkm))
						{
							values.Add(hugoSymbol, new Dictionary<string, double>());
						}
					}
					double hugoData;
					if (regionalSignals.TryGetValue(ASETools.ConvertToNonExcelString(hugoSymbol), out hugoData))
					{
						lock (values)
						{

							values[hugoSymbol].Add(case_.case_id, hugoData);
						}
					}
				}
			}
		}

		static void Shuffle<T>(IList<T> list)
		{
			Random rng = new Random();

			int n = list.Count;
			while (n > 1)
			{
				n--;
				int k = rng.Next(n + 1);
				T value = list[k];
				list[k] = list[n];
				list[n] = value;
			}
		}

		static void Main(string[] args)
		{
			var timer = new Stopwatch();
			timer.Start();

			var configuration = ASETools.Configuration.loadFromFile(args);

			if (null == configuration)
			{
				Console.WriteLine("Giving up because we were unable to load configuration.");
				return;
			}

			// only allow flag for allele-specific expression or case ids
			if (configuration.commandLineArgs.Count() < 1)
			{
				PrintUsageMessage();
				return;
			}

			var cases = ASETools.Case.LoadCases(configuration.casesFilePathname).Select(r => r.Value).ToList();
			
			// shuffle cases to avoid contention among other tasks
			Shuffle<ASETools.Case>(cases);

			if (null == cases)
			{
				Console.WriteLine("Unable to load cases file " + configuration.casesFilePathname + ".  You must generate cases before running ExpressionNearMutations.");
			}

			List<string> hugoSymbols = new List<string>();
			foreach (var arg in configuration.commandLineArgs)
			{
				hugoSymbols.Add(arg);
			}

			// load ase values for all
			var aseCases = cases.ToList();
			Shuffle<ASETools.Case>(aseCases);

			// load ase
			var threads = new List<Thread>();
			for (int i = 0; i < Environment.ProcessorCount / 3; i++)
			{
				threads.Add(new Thread(() => loadRegionalSignal(aseCases, hugoSymbols, aseValues, true)));
			}

			var mutationCounts = aseValues.Item1;

			// load methylation
			var methylationCases = cases.ToList();
			Shuffle<ASETools.Case>(methylationCases);

			for (int i = 0; i < Environment.ProcessorCount / 3; i++)
			{
				threads.Add(new Thread(() => loadMethylation(methylationCases, hugoSymbols)));
			}

			// load expression
			var expCases = cases.ToList();
			Shuffle<ASETools.Case>(expCases);

			for (int i = 0; i < Environment.ProcessorCount / 3; i++)
			{
				threads.Add(new Thread(() => loadFPKM(expCases, hugoSymbols, expressionValues)));
			}

			threads.ForEach(t => t.Start());
			threads.ForEach(t => t.Join());

			foreach (var hugoSymbol in hugoSymbols)
			{
				var outputFilename = configuration.finalResultsDirectory + hugoSymbol + ".txt";

				var outputFile = ASETools.CreateStreamWriterWithRetry(outputFilename);

				// Get all composite ref names for this gene
				var compositeREFs = methylationValues[hugoSymbol].Values.SelectMany(r => r.Keys).Distinct();

				// write header
				outputFile.WriteLine("Id\tDisease\tMutationCount\tASEValue\tExpressionValue\t" + String.Join("\t", compositeREFs));

				foreach (var case_ in cases)
				{
					int mutationCount;
					if (!mutationCounts[hugoSymbol].TryGetValue(case_.case_id, out mutationCount))
					{
						continue;
					}

					// write case id
					outputFile.Write(case_.case_id + "\t");

					// write disease
					outputFile.Write(case_.disease() + "\t");

					// write mutation count
					outputFile.Write(mutationCount + "\t");

					// write ASE
					double value;
					if (aseValues.Item2[hugoSymbol].TryGetValue(case_.case_id, out value) && value > Double.NegativeInfinity)
					{
						outputFile.Write(value + "\t");

					}
					else
					{
						outputFile.Write("*\t");
					}

					// write expression
					if (expressionValues[hugoSymbol].TryGetValue(case_.case_id, out value) && value > Double.NegativeInfinity)
					{
						outputFile.Write(value + "\t");
					}
					else
					{
						outputFile.Write("*\t");
					}


					//// write methylation
					Dictionary<string, double> methylationValue;
					if (methylationValues[hugoSymbol].TryGetValue(case_.case_id, out methylationValue))
					{

						foreach (var cREF in compositeREFs)
						{
							double M_Value;
							if (methylationValue.TryGetValue(cREF, out M_Value))
							{
								outputFile.Write(M_Value + "\t");
							}
							else
							{
								outputFile.Write("*\t");
							}
						}

					}
					else
					{
						// Case no methylation data. write *'s. This can occur if there was no info for loss of heterozygousity
						outputFile.Write(String.Join("\t", compositeREFs.Select(r => "*")) + "\t");
					}

					outputFile.WriteLine();
				}

				outputFile.Close();

			} // foreach hugo symbol

		}
	}
}
