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
		static int binValue(double value, double binCount)
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

		static Dictionary<string, ASETools.GeneLocationInfo> knownGenes = 
			ASETools.readKnownGeneFile(ASETools.ASEConfirguation.defaultGeneLocationInformation);



		static void ProcessRawMethylation(List<Tuple<string, string>> regionalFiles, double binCount)
		{

			while (true)
			{
				string expressionFilename;
				string methylationFilename;

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
					methylationFilename = regionalFiles[0].Item2;
					regionalFiles.RemoveAt(0);
				}

				if (!methylationFilename.Contains("HumanMethylation450"))
				{
					continue;
				}

				// Read in expression files
				var aseSignal = ASETools.RegionalSignalFile.ReadFile(expressionFilename);

				// get dictionary for methylation values by chromosome
				Dictionary<string, List<ASETools.AnnotationLine>> methylationValues = 
					ASETools.AnnotationLine.ReadFile(methylationFilename, methylationFilename, false)
						.GroupBy(r => r.compositeREF.Chromosome)
						.ToDictionary(r => r.Key, t => t.Select(k => k).ToList());

				// For all genes with ASE for this sample, get the methylation values for all regions
				foreach (var gene in aseSignal.Item1)
				{
					var geneSymbol = gene.Key;
					var ASEValues = gene.Value;

					// make sure gene is verified and get location
					ASETools.GeneLocationInfo geneLocation;
					if (!knownGenes.TryGetValue(ASETools.ConvertToNonExcelString(geneSymbol), out geneLocation))
					{
						continue;
					}

					var aseValue = ASEValues[0]; // get ASE Value at 0 kb

					if (aseValue >= 0)
					{

						// initialize bin in heatmap, if necessary
						var bin = binValue(aseValue, binCount);
						lock (heatmap)
						{
							if (!heatmap.ContainsKey(bin))
							{
								heatmap.Add(bin, new List<Tuple<int, double>>());
								// Initialize areas for each distance from gene
								for (var k = 0; k < 20; k++) // TODO rm hardcoded 20
								{
									heatmap[bin].Add(new Tuple<int, double>(0, 0.0));
								}
							}
						}

						var start = geneLocation.isoforms[0].txStart;
						var end = geneLocation.isoforms[0].txEnd;

						double minDistance = start;
						double maxDistance = end;
						double lastMinDistance = start;
						double lastMaxDistance = end;

						// Foreach range
						// Because we have 0 (in the gene), this range is 2^(20 - 2) * 1000
						// here exclude the previous region size
						for (int i = 0; i < ASETools.GeneExpression.nRegionSizes; i++)
						{


							if (i == 0)
							{
								// special case: only look at gene
								minDistance = start;
								maxDistance = end;
								lastMinDistance = start;
								lastMaxDistance = end;
							}
							else
							{
								lastMinDistance = minDistance;
								lastMaxDistance = maxDistance;

								minDistance = Math.Max(0, start - Math.Pow(2, i - 1) * 1000);
								maxDistance = end + Math.Pow(2, i - 1) * 1000;
							}

							// find all methylation values that are < distance and > lastDistance. Take the mean.
							List<ASETools.AnnotationLine> methylationValuesForChromsome;
							if (!methylationValues.TryGetValue(geneLocation.chromosome, out methylationValuesForChromsome))
							{
								// chromosome could not be found
								continue;
							}

							List<ASETools.AnnotationLine> filteredMethylationValues;
							if (i == 0)
							{
								filteredMethylationValues = methylationValuesForChromsome.Where(r =>
								{
								// case if methylation point is "higher" than current gene
								return (r.compositeREF.Start > minDistance && r.compositeREF.Start < maxDistance);
								}).ToList();
							}
							else
							{
								filteredMethylationValues = methylationValuesForChromsome.Where(r =>
								{
									// case if methylation point is "higher" than current gene
									return ((r.compositeREF.Start > lastMaxDistance && r.compositeREF.Start < maxDistance) ||
									// ase if methylation point is "lower" than current gene
									(r.compositeREF.Start > minDistance && r.compositeREF.Start < lastMinDistance));
								}).ToList();
							}
						
							var totalMValues = filteredMethylationValues.Select(r => r.M_Value).ToList().Sum();

							lock (heatmap[bin])
							{
								// Update heatmap values for this bin and ASM distance
								heatmap[bin][i] = new Tuple<int, double>(heatmap[bin][i].Item1 + filteredMethylationValues.Count(),
									heatmap[bin][i].Item2 + totalMValues);
							}

						} // foreach distance from gene

					} // if ASM Value is valid

				} // foreach gene from the ASE results

			} // while true
		}

		static void ProcessMethylation(List<Tuple<string, string>> regionalFiles, double binCount)
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

						lock (heatmap)
						{
							if (!heatmap.ContainsKey(bin))
							{
								heatmap.Add(bin, new List<Tuple<int, double>>());
								// Initialize areas for each distance from gene
								for (var k = 0; k < 20; k++) // TODO rm hardcoded 20
								{
									heatmap[bin].Add(new Tuple<int, double>(0, 0.0));
								}
							}
						}

						// Get data for ASM for all distances
						var asmData = comparisonSignal.Item1[geneSymbol];

						if (labels.Count() == 0)
						{
							labels = comparisonSignal.Item2;
						}

						for (var i = 0; i < asmData.Count(); i++)
						{

							double asm = asmData[i];

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

		static void ProcessBisulfite(Dictionary<string, ASETools.BisulfateCase> bisulfiteCases, Dictionary<string, ASETools.Case> cases, bool forTumor, double binCount)
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

		static void ProcessIlluminaArray(Dictionary<string, ASETools.Case> cases, bool forTumor, double binCount, bool forAlleleSpecificExpression)
		{
			var regionalFiles = new List<Tuple<string, string>>();

			foreach (var caseEntry in cases)
			{
				var case_ = caseEntry.Value;

				if (forTumor)
				{
					if (forAlleleSpecificExpression)
					{
						if (case_.tumor_allele_specific_gene_expression_filename != "" && case_.tumor_methylation_filename != "")
						{
							regionalFiles.Add(new Tuple<string, string>(case_.tumor_allele_specific_gene_expression_filename, 
								case_.tumor_methylation_filename));
						}
						else
						{
							Console.WriteLine("Case " + caseEntry.Key + " does not have required files. Skipping...");
						}
					}
					else
					{
						if (case_.gene_expression_filename != "" && case_.tumor_methylation_filename != "")
						{
							regionalFiles.Add(new Tuple<string, string>(case_.gene_expression_filename, 
								case_.tumor_methylation_filename));
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
						if (case_.normal_allele_specific_gene_expression_filename != "" && case_.normal_methylation_filename != "")
						{
							regionalFiles.Add(new Tuple<string, string>(case_.normal_allele_specific_gene_expression_filename, 
								case_.normal_methylation_filename));
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

			var threads = new List<Thread>();
			for (int i = 0; i < Environment.ProcessorCount; i++)
			{
				threads.Add(new Thread(() => ProcessRawMethylation(regionalFiles, binCount)));
			}

			threads.ForEach(t => t.Start());
			threads.ForEach(t => t.Join());
		}

		static void Main(string[] args)
		{
			// filter out chr x and y
			knownGenes = knownGenes.Where(r => r.Value.chromosome != "chrX" && r.Value.chromosome != "chrY")
				.ToDictionary(x => x.Key, x => x.Value);


			var configuration = ASETools.ASEConfirguation.loadFromFile(args);
			var cases = ASETools.Case.LoadCases(configuration.casesFilePathname).ToDictionary(x => x.Key, x => x.Value);

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

			double binCount;
			try {
				binCount = Convert.ToDouble(configuration.commandLineArgs[2]);
			} catch (Exception)
			{
				Console.WriteLine("MethylationHeatmaps: invalid bin count");
				return;
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

				ProcessIlluminaArray(cases, forTumor, binCount, forAlleleSpecificExpression);
			}


			// get filename for heatmap
			var expExt = forAlleleSpecificExpression ? "_ase" : "_expression";
			var ext = bisulfite ? "_bisulfite.txt" : "_array.txt";
			var heatmapName = "heatmap" + (forTumor ? "_tumor" : "_normal");

			var heatmapFilename = ASETools.ASEConfirguation.bisulfiteDirectory + heatmapName + expExt + ext + ".test";
			Console.WriteLine("Saving to " + heatmapFilename);

			var writer = ASETools.CreateStreamWriterWithRetry(heatmapFilename);

			// write header
			for (int i = 0; i < ASETools.GeneExpression.nRegionSizes; i++)
			{
				
				var label = (i==0) ? 0 : Math.Pow(2, i - 1) * 1000;

				writer.Write("\t" + label + " kb");
			}
			writer.WriteLine();

			// write heatmap
			var sortedKeys = heatmap.Keys.ToList();
			sortedKeys.Sort();

			foreach (var key in sortedKeys)
			{
				var bin = heatmap[key];

				// write bin count
				writer.Write(key);

				// foreach distance, write distance count
				foreach (var value in bin)
				{
					writer.Write("\t" + value.Item2 / value.Item1);
				}
				writer.WriteLine();
			}
			writer.Close();
		} // Main

	} // Process
}
