using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using System.Diagnostics;
using ASELib;
using MethylationDistributionPromotors;

// Case 2: Run t test, trying to assign distributions for methylation points. This is to only
// look at 1 group (ie. ase) Thie purpose of this is to see what percentage of ASE significant
// results can be explained by ASE. This should be written generically enough so it can
// also be used to test no ASE for 0 mutations, testing what percentage of these cases fall into the
// 'full methylation' category

namespace MethylationByMutationCount
{
	class Program
	{

		class MutationGroup
		{
			public double averageMethylation;
			public double averageASE;
			public double averageExpression;
		}


		static void PrintUsageMessage()
		{
			Console.WriteLine("usage: MethylationByMutationCount hugoSymbols");
			Console.WriteLine("hugo symbols are separated by space");
		}

		static Tuple<Dictionary<string, int>, Dictionary<string, double>> loadRegionalSignal(List<Tuple<string, string>> idsAndFilenames, string hugoSymbol)
		{
			var mutations = new Dictionary<string, int>();
			var values = new Dictionary<string, double>();

			while (true)
			{
				string id;
				string filename;
				lock (idsAndFilenames)
				{
					if (idsAndFilenames.Count() == 0)
					{
						//
						// No more work, we're done.
						//
						return new Tuple<Dictionary<string, int>, Dictionary<string, double>>(mutations, values);
					}
					id = idsAndFilenames[0].Item1;
					filename = idsAndFilenames[0].Item2;

					idsAndFilenames.RemoveAt(0);
				}

				if (filename == "")
				{
					// no data. continue
					continue;
				}

				var regionalSignals = ASETools.RegionalSignalFile.ReadFile(filename);

				double[] hugoData;
				if (regionalSignals.Item1.TryGetValue(ASETools.ConvertToExcelString(hugoSymbol), out hugoData))
				{
					values.Add(id, hugoData[0]);
					mutations.Add(id, regionalSignals.Item3[ASETools.ConvertToExcelString(hugoSymbol)]);
				}
			}
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
			if (configuration.commandLineArgs.Count() < 1)
			{
				PrintUsageMessage();
				return;
			}

			var cases = ASETools.Case.LoadCases(configuration.casesFilePathname);

			if (null == cases)
			{
				Console.WriteLine("Unable to load cases file " + configuration.casesFilePathname + ".  You must generate cases before running ExpressionNearMutations.");
			}

			var hugoSymbol = configuration.commandLineArgs[0];

			// load ase values for all
			var idsAndFilenames = cases.Select(r => 
				new Tuple<string, string>(r.Value.case_id, r.Value.tumor_allele_specific_gene_expression_filename)).ToList();
			var mutationsAndASE = loadRegionalSignal(idsAndFilenames, hugoSymbol);

			var mutationCounts = mutationsAndASE.Item1;
			var aseValues = mutationsAndASE.Item2;

			// load methylation
			idsAndFilenames = cases.Select(r => 
				new Tuple<string, string>(r.Value.case_id, r.Value.tumor_regional_methylation_filename)).ToList();
			Dictionary<string, double> methylationValues = loadRegionalSignal(idsAndFilenames, hugoSymbol).Item2;

			// load expression
			idsAndFilenames = cases.Select(r => 
				new Tuple<string, string>(r.Value.case_id, r.Value.gene_expression_filename)).ToList();
			Dictionary<string, double> expressionValues = loadRegionalSignal(idsAndFilenames, hugoSymbol).Item2;

			var outputFilename = configuration.finalResultsDirectory + hugoSymbol + ".txt";

			var outputFile = ASETools.CreateStreamWriterWithRetry(outputFilename);

			// write header
			outputFile.WriteLine("Mutation Count\tASE Value\tMethylation Value\tExpression Value");

			foreach (var case_ in cases)
			{
				// write mutation count
				int mutationCount;
				if (!mutationCounts.TryGetValue(case_.Value.case_id, out mutationCount))
				{
					continue;
				}

				outputFile.Write(mutationCounts[case_.Value.case_id] + "\t");

				// write ASE
				double value;
				if (aseValues.TryGetValue(case_.Value.case_id, out value) && value > Double.NegativeInfinity)
				{
					outputFile.Write(value + "\t");

				}
				else
				{
					outputFile.Write("*\t");
				}

				// write methylation
				if (methylationValues.TryGetValue(case_.Value.case_id, out value) && value > Double.NegativeInfinity)
				{
					outputFile.Write(value + "\t");
				}
				else
				{
					outputFile.Write("*\t");
				}


				if (expressionValues.TryGetValue(case_.Value.case_id, out value) && value > Double.NegativeInfinity)
				{
					outputFile.Write(value + "\t");
				}
				else
				{
					outputFile.Write("*\t");
				}

				outputFile.WriteLine();
			}

			outputFile.Close();

		}
	}
}
