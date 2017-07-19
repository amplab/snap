using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using ASELib;
using System.IO;

namespace FindSignificantGenes
{
	class Program
	{

		static void PrintUsageMessage()
		{
			Console.WriteLine("usage: FindSignificantGenes {-i}");
			Console.WriteLine("-i is inclusive, blank is exclusive");
		}

		static void Main(string[] args)
		{

			var configuration = ASETools.Configuration.loadFromFile(args);

			if (null == configuration)
			{
				Console.WriteLine("Giving up because we were unable to load configuration.");
				return;
			}

			bool inclusive = false;
			if (configuration.commandLineArgs.Count() == 1)
			{
				inclusive = configuration.commandLineArgs[0] == "-i";
			}
			else if (configuration.commandLineArgs.Count() > 1)
			{
				// error
				return;
			}


			var fileExtention = "AlleleSpecificExpressionDistributionByMutationCount";
			var diseases = new List<string>();

			string[] filePaths = Directory.GetFiles(configuration.finalResultsDirectory, fileExtention + "*.txt");

			var threshold = 1e-6;

			// dictionary of gene
			Dictionary<string, Dictionary<string, string>> pvalChart = new Dictionary<string, Dictionary<string, string>>();

			foreach (var filePath in filePaths)
			{
				var disease = filePath.Split(new string[] { fileExtention }, StringSplitOptions.None).Last();
				disease = disease.Split(new string[] { ".txt" }, StringSplitOptions.None).First();

				if (disease.Length == 0)
				{
					disease = "pan";

				}
				else
				{
					disease = disease.Substring(1);
				}

				diseases.Add(disease);

				// open file
				var results = ASETools.ExpressionResultsLine.readFile(filePath);


				// find all genes that have a pvalue below threshold
				List<ASETools.ExpressionResultsLine> significantGenes;
				if (inclusive)
				{
					significantGenes =
						results.Where(r => r.nonExclusiveResultsByRange.Where(e => e.oneVsNotOne < threshold && e.oneVsNotOne >= 0).Count() > 0)
						.ToList();
				}
				else {
					significantGenes = results.Where(r => r.exclusiveResultsByRange.Where(e => e.oneVsNotOne < threshold && e.oneVsNotOne >= 0).Count() > 0)
						.ToList();

				}

				foreach (var significantGene in significantGenes)
				{
					if (!pvalChart.ContainsKey(significantGene.hugo_symbol))
					{
						pvalChart.Add(significantGene.hugo_symbol, new Dictionary<string, string>());
					}

					double lowestPval;
					ASETools.SingleExpressionResult location;

					if (inclusive)
					{
						lowestPval = significantGene.nonExclusiveResultsByRange.Where(r => r.oneVsNotOne >= 0).Min(r => r.oneVsNotOne);
						location = significantGene.nonExclusiveResultsByRange.Where(r => r.oneVsNotOne == lowestPval).First();
					}
					else {
						lowestPval = significantGene.exclusiveResultsByRange.Where(r => r.oneVsNotOne >= 0).Min(r => r.oneVsNotOne);
						location = significantGene.exclusiveResultsByRange.Where(r => r.oneVsNotOne == lowestPval).First();
					}

					pvalChart[significantGene.hugo_symbol].Add(disease, lowestPval.ToString());
				}

			}


			// write out results to file
			var outputFilename = configuration.finalResultsDirectory + "overlappingGenes";
			outputFilename += inclusive? "_inclusive.txt" : "_exclusive.txt";

			Console.WriteLine("saving to " + outputFilename);

			var writer = ASETools.CreateStreamWriterWithRetry(outputFilename);

			// write header for diseases
			writer.Write("Gene_Symbol\t");
			writer.Write(String.Join("\t", diseases));
			writer.WriteLine();

			// foreach gene
			foreach (var gene in pvalChart)
			{
				writer.Write(gene.Key);

				foreach (var disease in diseases)
				{
					string value;
					if (gene.Value.TryGetValue(disease, out value))
					{
						writer.Write("\t" + value);
					}
					else
					{
						writer.Write("\t*");
					}
				}

				writer.WriteLine();
			} // foreach gene

			writer.Close();
		}
	}
}
