using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using ASELib;

namespace MethylationHeatmaps
{
	class Program
	{

		// map of binned ASE values, ASM values across distances
		// List is of int total counts, double is running sum of ASM values
		static Dictionary<int, List<Tuple<int, double>>> heatmap = new Dictionary<int, List<Tuple<int, double>>>();
		// List of heatmap xaxis labels
		static List<string> labels = new List<string>();

		// bins a value between 0 and 1
		static int binValue(double value, int binCount)
		{

			var x = value / (1.0 / binCount);
			var bin =  Convert.ToInt32(Math.Floor(x));

			// corner case for last bin
			if (bin == binCount)
				bin -= 1;
			return bin;
		}

		static void PrintUsageMessage()
		{
			Console.WriteLine("usage: MethylationHeatmaps {-a|-b} {-o}");
			Console.WriteLine("-a means 450k array data, -b means bisulfite data.");
		}

		static void ProcessMethylation(List<Tuple<string, string>> regionalFiles, int binCount)
		{

			while (true)
			{
				string aseFilename = null;
				string comparisonFilename = null;

				lock (regionalFiles)
				{
					if (regionalFiles.Count() == 0)
					{
						return;
					}

					aseFilename = regionalFiles[0].Item1;
					comparisonFilename = regionalFiles[0].Item2;
					regionalFiles.RemoveAt(0);
				}
				// Read in expression files
				var aseFile = new ASETools.ASEExpressionFile();
				aseFile.ReadFile(aseFilename);

				var comparisonFile = new ASETools.ASEExpressionFile();
				comparisonFile.ReadFile(comparisonFilename);

				foreach (var gene in aseFile.expressionMap)
				{
					var geneSymbol = gene.Key;
					var values = gene.Value;
					var aseValue = values[0];

					// we are only looking at ASE values at 0 distance
					if (aseValue >= 0)
					{

						// Process ASM values for this gene
						if (!comparisonFile.expressionMap.ContainsKey(geneSymbol))
						{
							// No data for this gene exists in the methylation set.
							continue;
						}

						var bin = binValue(aseValue, binCount);

						// Get data for ASM for all distances
						var asmData = comparisonFile.expressionMap[geneSymbol];

						if (labels.Count() == 0)
						{
							labels = comparisonFile.index;
						}

						lock (heatmap)
						{
							if (!heatmap.ContainsKey(bin))
							{
								// add key for this bin size
								heatmap.Add(bin, new List<Tuple<int, double>>());

								// Initialize areas for each distance from gene
								for (var i = 0; i < asmData.Count(); i++)
								{
									heatmap[bin].Add(new Tuple<int, double>(0, 0.0));
								}

							}
						}

						for (var i = 0; i < asmData.Count(); i++)
						{
							var asm = asmData[i];
							if (asm >= 0)
							{
								lock (heatmap)
								{
									// Update heatmap values for this bin and ASM distance
									heatmap[bin][i] = new Tuple<int, double>(heatmap[bin][i].Item1 + 1, heatmap[bin][i].Item2 + asm);
								}
							}
						} // foreach methylation site for this gene


					}

				} // foreach gene from the ASE results
			} // while true
		}

		static void ProcessBisulfite(Dictionary<string, ASETools.BisulfateCase> bisulfiteCases, Dictionary<string, ASETools.Case> cases, int binCount = 10)
		{
			var regionalFiles = new List<Tuple<string, string>>();

			foreach (var bCase in bisulfiteCases)
			{
				Console.WriteLine("Processing case " + bCase.Key);
				var bisulfateCase_ = bisulfiteCases[bCase.Key];
				var case_ = cases[bCase.Key];

				// TODO: save this somewhere in the cases file 
				string directory = ASETools.GetDirectoryPathFromFullyQualifiedFilename(bisulfateCase_.realigned_selected_variants_filename);
				string analysisId = ASETools.GetAnalysisIdFromPathname(bisulfateCase_.realigned_selected_variants_filename);
				if ("" == directory || "" == analysisId)
				{
					Console.WriteLine("Couldn't parse input pathname, which is supposed to be absolute and include an analysis ID: " + bisulfateCase_.realigned_selected_variants_filename);
					continue;
				}

				var asmFilename = directory + analysisId + (ASETools.bisulfiteAlleleSpecificMethylationExtension);

				regionalFiles.Add(new Tuple<string, string>(case_.tumor_allele_specific_gene_expression_filename, asmFilename));

			} // foreach bisulfite case

			ProcessMethylation(regionalFiles, binCount);

		}

		static void ProcessIlluminaArray(Dictionary<string, ASETools.Case> cases, int binCount = 10)
		{
			var regionalFiles = new List<Tuple<string, string>>();

			foreach (var caseEntry in cases)
			{
				Console.WriteLine("Processing case " + caseEntry.Key);
				var case_ = caseEntry.Value;

				regionalFiles.Add(new Tuple<string, string>(case_.tumor_allele_specific_gene_expression_filename, case_.tumor_regional_methylation_filename));

			}
			ProcessMethylation(regionalFiles, binCount);
		}

		static void Main(string[] args)
		{
			var configuration = ASETools.ASEConfirguation.loadFromFile(args);
			var cases = ASETools.Case.LoadCases(configuration.casesFilePathname);

			// only allow flag for allele-specific expression or case ids
			if (configuration.commandLineArgs.Count() < 1 && (configuration.commandLineArgs[0] == "-b" || configuration.commandLineArgs[0] == "-a"))
			{
				PrintUsageMessage();
				return;
			}

			var bisulfiteCases = ASETools.BisulfateCase.loadCases(ASETools.ASEConfirguation.bisulfiteCasesFilePathname);

			var bisulfite = configuration.commandLineArgs[0] == "-b";
			if (bisulfite)
			{
				// Add cases with bisulfate data. We will use the bam and bisulfate files

				// filter out only bisulfite cases
				cases = cases.Where(case_ => bisulfiteCases.ContainsKey(case_.Value.case_id)).ToDictionary(x => x.Key, x => x.Value);

				ProcessBisulfite(bisulfiteCases, cases);
			}
			else {
				cases = cases.Where(case_ => bisulfiteCases.ContainsKey(case_.Value.case_id)).ToDictionary(x => x.Key, x => x.Value);
				ProcessIlluminaArray(cases);
			}



			// get filename for heatmap
			var ext = bisulfite ? "_bisulfite.txt" : "_array.txt";
			var heatmapFilename = ASETools.ASEConfirguation.bisulfiteDirectory + @"heatmap" + ext;
			Console.WriteLine("Saving to " + heatmapFilename);

			var writer = ASETools.CreateStreamWriterWithRetry(heatmapFilename);

			// write header
			foreach (string label in labels)
			{
				var finalLabel = label;
				int index = label.LastIndexOf("(");
				if (index > 0)
					finalLabel = label.Substring(0, index);

				writer.Write('\t' + finalLabel);
			}
			writer.WriteLine();

			// write heatmap
			foreach (var bin in heatmap)
			{
				// write bin count
				writer.Write(bin.Key);

				// foreach distance, write distance count
				foreach (var value in bin.Value)
				{
					writer.Write("\t" + value.Item2 / value.Item1);
				}
				writer.WriteLine();
			}
			writer.Close();
		} // Main

	} // Process
}
