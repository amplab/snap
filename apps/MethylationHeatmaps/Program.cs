using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.IO;
using System.Threading.Tasks;
using ASELib;
using System.Diagnostics;
using System.Threading;

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
			Console.WriteLine("usage: MethylationHeatmaps {-array|-bilsulfite} {-t|-n} binCount {-a}");
			Console.WriteLine("-a means 450k array data, -b means bisulfite data.");
			Console.WriteLine("-t means tumor, -n means normal");
			Console.WriteLine("-a means allele specific expression. Absent means raw expression");
		}

		static void ProcessMethylation(List<Tuple<string, string>> regionalFiles, int binCount)
		{

			while (true)
			{
				string expressionFilename = null;
				string comparisonFilename = null;

				lock (regionalFiles)
				{
					if (regionalFiles.Count() == 0)
					{
						return;
					}

					if (regionalFiles.Count() % 500 == 0)
					{
						Console.WriteLine(regionalFiles.Count() + " cases left");
					}

					expressionFilename = regionalFiles[0].Item1;
					comparisonFilename = regionalFiles[0].Item2;
					regionalFiles.RemoveAt(0);
				}

				// Read in expression files
				var aseSignal = ASETools.RegionalSignalFile.ReadFile(expressionFilename);

				var comparisonSignal = ASETools.RegionalSignalFile.ReadFile(comparisonFilename);

				foreach (var gene in aseSignal.Item1)
				{
					var geneSymbol = gene.Key;
					var values = gene.Value;
					var aseValue = values[0];

					// we are only looking at ASE values at 0 distance
					if (aseValue >= 0)
					{

						// Process ASM values for this gene
						if (!comparisonSignal.Item1.ContainsKey(geneSymbol))
						{
							// No data for this gene exists in the methylation set.
							continue;
						}

						var bin = binValue(aseValue, binCount);

						// Get data for ASM for all distances
						var asmData = comparisonSignal.Item1[geneSymbol];

						if (labels.Count() == 0)
						{
							labels = comparisonSignal.Item2;
						}

						for (var i = 0; i < asmData.Count(); i++)
						{
							var asm = asmData[i];
							if (asm >= 0)
							{
								lock (heatmap[bin])
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

		static void ProcessBisulfite(Dictionary<string, ASETools.BisulfateCase> bisulfiteCases, Dictionary<string, ASETools.Case> cases, bool forTumor, int binCount)
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

				string asmFilename;
				string aseFilename;
				if (forTumor)
				{
					asmFilename = directory + analysisId + (ASETools.tumorBisulfiteAlleleSpecificMethylationExtension);
					aseFilename = case_.tumor_allele_specific_gene_expression_filename;
				}
				else
				{
					asmFilename = directory + analysisId + (ASETools.normalBisulfiteAlleleSpecificMethylationExtension);
					aseFilename = case_.normal_allele_specific_gene_expression_filename;
				}
				if (!File.Exists(asmFilename) || !File.Exists(aseFilename))
				{
					Console.WriteLine("ASM file " + asmFilename + " does not exist. Skipping...");
					continue;
				}

				regionalFiles.Add(new Tuple<string, string>(aseFilename, asmFilename));

			} // foreach bisulfite case


			var threads = new List<Thread>();
			for (int i = 0; i < Environment.ProcessorCount; i++)
			{
				threads.Add(new Thread(() => ProcessMethylation(regionalFiles, binCount)));
			}

			threads.ForEach(t => t.Start());
			threads.ForEach(t => t.Join());


		}

		static void ProcessIlluminaArray(Dictionary<string, ASETools.Case> cases, bool forTumor, int binCount, bool forAlleleSpecificExpression)
		{
			var regionalFiles = new List<Tuple<string, string>>();

			foreach (var caseEntry in cases)
			{
				var case_ = caseEntry.Value;

				if (forTumor)
				{
					if (forAlleleSpecificExpression)
					{
						if (case_.tumor_allele_specific_gene_expression_filename != "" && case_.tumor_regional_methylation_filename != "")
						{
							regionalFiles.Add(new Tuple<string, string>(case_.tumor_allele_specific_gene_expression_filename, case_.tumor_regional_methylation_filename));
						}
						else
						{
							Console.WriteLine("Case " + caseEntry.Key + " does not have required files. Skipping...");
						}
					}
					else
					{
						if (case_.gene_expression_filename != "" && case_.tumor_regional_methylation_filename != "")
						{
							regionalFiles.Add(new Tuple<string, string>(case_.gene_expression_filename, case_.tumor_regional_methylation_filename));
						}
						else
						{
							Console.WriteLine("Case " + caseEntry.Key + " does not have required files. Skipping...");
						}
					}
			
				}
				else
				{
					if (forAlleleSpecificExpression)
					{
						if (case_.normal_allele_specific_gene_expression_filename != "" && case_.normal_regional_methylation_filename != "")
						{
							regionalFiles.Add(new Tuple<string, string>(case_.normal_allele_specific_gene_expression_filename, case_.normal_regional_methylation_filename));
						}
						else
						{
							Console.WriteLine("Case " + caseEntry.Key + " does not have required files. Skipping...");
						}
					}
					else
					{
						throw new Exception("normal expression not yet implemented");
					}
				}

			}

			//var threads = new List<Thread>();
			//for (int i = 0; i < Environment.ProcessorCount; i++)
			//{
			//	threads.Add(new Thread(() => ProcessMethylation(regionalFiles, binCount)));
			//}

			//threads.ForEach(t => t.Start());
			//threads.ForEach(t => t.Join());
			ProcessMethylation(regionalFiles, binCount);
		}

		static void Main(string[] args)
		{
			var configuration = ASETools.ASEConfirguation.loadFromFile(args);
			var cases = ASETools.Case.LoadCases(configuration.casesFilePathname);

			// only allow flag for allele-specific expression or case ids
			if (configuration.commandLineArgs.Count() < 3 || configuration.commandLineArgs.Count() > 4 || (configuration.commandLineArgs[0] != "-bisulfite" && configuration.commandLineArgs[0] != "-array")
				 || (configuration.commandLineArgs[1] != "-t" && configuration.commandLineArgs[1] != "-n"))
			{
				PrintUsageMessage();
				return;
			}

			var bisulfite = configuration.commandLineArgs[0] == "-bisulfite";
			var forTumor = configuration.commandLineArgs[1] == "-t";
			bool forAlleleSpecificExpression = false;

			if (configuration.commandLineArgs.Count() == 4 && configuration.commandLineArgs[3] == "-a")
			{
				forAlleleSpecificExpression = true;
			}

			int binCount;
			try {
				binCount = Convert.ToInt32(configuration.commandLineArgs[2]);
			} catch (Exception)
			{
				Console.WriteLine("MethylationHeatmaps: invalid bin count");
				return;
			}

			var regionCount = 66; // TODO don't hardcode this

			// initialize heatmap
			for (var bin = 0; bin < binCount; bin++)
			{
				// add key for this bin size
				heatmap.Add(bin, new List<Tuple<int, double>>());

				// Initialize areas for each distance from gene
				for (var i = 0; i < regionCount; i++)
				{
					heatmap[bin].Add(new Tuple<int, double>(0, 0.0));
				}

			}

			if (bisulfite)
			{
				var bisulfiteCases = ASETools.BisulfateCase.loadCases(ASETools.ASEConfirguation.bisulfiteCasesFilePathname);
				// Add cases with bisulfate data. We will use the bam and bisulfate files

				// filter out only bisulfite cases
				cases = cases.Where(case_ => bisulfiteCases.ContainsKey(case_.Value.case_id)).ToDictionary(x => x.Key, x => x.Value);

				ProcessBisulfite(bisulfiteCases, cases, forTumor, binCount);
			}
			else {
				if (forTumor)
					cases = cases.Where(r => r.Value.tumor_methylation_filename.Contains("HumanMethylation450")).ToDictionary(x => x.Key, x => x.Value);
				else
					cases = cases.Where(r => r.Value.normal_methylation_filename.Contains("HumanMethylation450")).ToDictionary(x => x.Key, x => x.Value);
				var diseaseCounts = cases.Select(r => r.Value.disease()).GroupBy(r => r).Select(r => new Tuple<string, int>(r.Key, r.Count()));
				foreach (var dc in diseaseCounts)
				{
					Console.WriteLine(dc.Item1 + ": " + dc.Item2);
				}

				ProcessIlluminaArray(cases, forTumor, binCount, forAlleleSpecificExpression);
			}


			// get filename for heatmap
			var ext = bisulfite ? "_bisulfite.txt" : "_array.txt";
			var heatmapName = "heatmap" + (forTumor ? "_tumor" : "_normal");

			var heatmapFilename = ASETools.ASEConfirguation.bisulfiteDirectory  + heatmapName + ext;
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
