using System;
using System.Collections.Generic;
using System.Linq;
using System.Threading;
using ASELib;

namespace MutationDistributions
{
	class Program
	{
		// Dictionary of (disease, bin counts)
		static Dictionary<string, int[]> mutationByDisease = new Dictionary<string, int[]>();

		// Dictionary of (disease, bin counts)
		static GermlineDictionary germlineByDisease;

		// Dictionary of (hugoSymbol, bin counts)
		static Dictionary<string, int[]> mutationByGene = new Dictionary<string, int[]>();

		// Dictionary of (hugoSymbol, bin counts)
		static GermlineDictionary germlineByGene;

		static ASETools.GeneLocationsByNameAndChromosome knownGenes;

		static string outDirectory;

		public class GermlineDictionary
		{
			public GermlineDictionary(int binCount_)
			{
				binCount = binCount_;

				germlineDictionary = new List<Dictionary<string, int[]>>();

				for (int i = 0; i < 5; i++)
				{
					germlineDictionary.Add(new Dictionary<string, int[]>());
				}
			}

			public int length()
			{
				return germlineDictionary.Count();
			}

			public Dictionary<string, int[]> get(int i)
			{
				return germlineDictionary[i];
			}

			int[] getKeys(int mutationCount, double mutationASEValue)
			{
				// 0 -> 0 mutations
				if (mutationCount == 0)
				{
					return new int[] {0};
				}
				// 1 -> 1 mutation
				// 2 -> 1 mutation with ref ASE (0 - 0.3)
				// 3 -> 1 mutation with no ASE (0.3 - 0.7)
				// 4 -> 1 mutation with mutant ASE (0.7 - 1)
				if (mutationASEValue <= referenceThreshold)
				{
					return new int[] { 1, 2};
				}
				else if (mutationASEValue > referenceThreshold && mutationASEValue < altThreshold)
				{
					return new int[] { 1, 3};
				}
				else
				{
					return new int[] { 1, 4};
				}

			}

			public void Add(string disease, int mutationCount, double germlineASE, double mutationASE)
			{
				var keys = getKeys(mutationCount, mutationASE);
				foreach (var key in keys)
				{
					lock (germlineDictionary)
					{
						if (!germlineDictionary[key].ContainsKey(disease))
						{
							germlineDictionary[key].Add(disease, new int[binCount]);
						}
					}

					int bin = binValue(germlineASE, binCount);

					lock (germlineDictionary[key][disease])
					{
						germlineDictionary[key][disease][bin] += 1;
					}
				}


			}

			// Dictionary of (hugoSymbol, bin counts)
			// 0 -> 0 mutations
			// 1 -> 1 mutation
			// 3 -> 1 mutation with no ASE (0.3 - 0.7)
			// 4 -> 1 mutation with mutant ASE (0.7 - 1)
			double referenceThreshold = 0.3;
			double altThreshold = 0.7;
			List<Dictionary<string, int[]>> germlineDictionary;
			int binCount;

		}
		


		// bins a value between 0 and 1

		static int binValue(double value, double binCount)

		{
			var x = value / (1.0 / binCount);
			var bin = Convert.ToInt32(Math.Floor(x));

			// corner case for last bin

			if (bin == binCount)
				bin -= 1;

			return bin;
		}
		static int casesProcessed = 0;

		static void ProcessCases(List<ASETools.Case> cases, int binCount)
		{
			while (true)
			{
				ASETools.Case case_;
				lock (cases)
				{
					if (cases.Count() == 0)
					{
						return;
					}

					case_ = cases[0];

					cases.RemoveAt(0);
					casesProcessed += 1;
				}

				if (casesProcessed % 1000 == 0 && casesProcessed > 2)
				{
					Console.WriteLine(cases.Count() + " cases remaining. Saving temporarily...");
				}

				if (case_.annotated_selected_variants_filename == "")
				{
					continue;
				}

				var disease = case_.disease();
				lock (mutationByDisease)
				{
					if (!mutationByDisease.ContainsKey(disease))
					{
						mutationByDisease.Add(disease, new int[binCount]);
					}
				}

				var annotatedVariants = ASETools.AnnotatedVariant.readFile(case_.annotated_selected_variants_filename);
				var somaticVariants = annotatedVariants.Where(r => r.somaticMutation);

				foreach (var annotatedVariant in annotatedVariants)
				{
					if (annotatedVariant.hasSufficientReads())
					{
						double rnaFractionTumor = annotatedVariant.getTumorAlleleSpecificExpression(false);

						int bin = binValue(rnaFractionTumor, binCount);

						// bin value by disease and gene. 
						// keep track of whether it is germline or somatic
						if (annotatedVariant.somaticMutation)
						{
							var gene = annotatedVariant.Hugo_symbol;

							lock (mutationByGene)
							{
								if (!mutationByGene.ContainsKey(gene))
								{
									mutationByGene.Add(gene, new int[binCount]);
								}
								mutationByGene[gene][bin] += 1;
							}

							lock (mutationByDisease[disease])
							{
								mutationByDisease[disease][bin] += 1;
							}
						}
						else
						{

							// check if this germline overlaps a gene AND this sample has a mutation at this gene
							var genesOverlappingGermline =
								knownGenes.genesByChromosome[annotatedVariant.contig].Where(r =>
								r.containsLocus(annotatedVariant.contig, annotatedVariant.locus));

							foreach (var gene in genesOverlappingGermline)
							{
								// how many mutations does this gene have?
								var mutations = somaticVariants.Where(r => r.Hugo_symbol == gene.hugoSymbol);
								var mutationCount = mutations.Count();

								double mutationASE;
								if (mutationCount > 0)
								{
									if (!mutations.First().hasSufficientReads())
									{
										continue;
									}
									mutationASE = mutations.First().getTumorAlleleSpecificExpression(false);
								}
								else
								{
									mutationASE = 0.0; // Doesnt matter, we won't use it
								}

								germlineByDisease.Add(disease, mutationCount, rnaFractionTumor, mutationASE);
								germlineByGene.Add(gene.hugoSymbol, mutationCount, rnaFractionTumor, mutationASE);

							}
						}
					}

				}

			}

		}

		static void writeBins(string directory, Dictionary<string, int[]> somaticVariants,
			GermlineDictionary germlineVariants)
		{
			foreach (var category in somaticVariants)
			{
				// hugo symbol or disease
				var identifier = category.Key;

				int[] somatic = category.Value;
				var binCount = somatic.Count();

				var filename = directory + identifier + ".txt";
				var outFile = ASETools.CreateStreamWriterWithRetry(filename);

				// write header
				outFile.WriteLine("Type\tBin\tCount");

				for (var i = 0; i < binCount; i++)
				{
					var value = somatic[i];
					var bin = (double)i / binCount;

					outFile.WriteLine("mutation\t" + bin + "\t" + value);
				}

				int[] germline;
				for(var i = 0; i < germlineVariants.length(); i++)
				{
					var germlineVariantsForMutationCount = germlineVariants.get(i);

					if (germlineVariantsForMutationCount.TryGetValue(identifier, out germline))
					{
						for (var j = 0; j < binCount; j++)
						{
							var value = germline[j];
							var bin = (double) j / binCount;

							outFile.WriteLine("germline_" + i + "\t" + bin + "\t" + value);
						}
					}
				}
	

				outFile.Close();
			}


		}

		static void PrintUsageMessage()
		{
			Console.WriteLine("Usage: MutationDistributions binCount");

		}
		static void Main(string[] args)
		{
			var configuration = ASETools.Configuration.loadFromFile(args);
			outDirectory = configuration.finalResultsDirectory + @"mutationDistributions\";


			var genes = ASETools.readKnownGeneFile(configuration.geneLocationInformationFilename)
					.Where(r => !r.Value.inconsistent);

			knownGenes = new ASETools.GeneLocationsByNameAndChromosome(genes.ToDictionary(x => x.Key, x => x.Value));

			if (null == configuration)
			{
				Console.WriteLine("Giving up because we were unable to load configuration.");
				return;
			}

			if (configuration.commandLineArgs.Count() != 1)
			{
				PrintUsageMessage();
				return;
			}
			var binCount = Convert.ToInt32(configuration.commandLineArgs[0]);

			germlineByDisease = new GermlineDictionary(binCount);
			germlineByGene = new GermlineDictionary(binCount);

			var cases = ASETools.Case.LoadCases(configuration.casesFilePathname).Select(r => r.Value).ToList();
			Console.WriteLine(cases.Count() + " cases to process");

			var threads = new List<Thread>();
			for (int i = 0; i < Environment.ProcessorCount; i++)
			{
				threads.Add(new Thread(() => ProcessCases(cases, binCount)));
			}

			threads.ForEach(th => th.Start());
			threads.ForEach(th => th.Join());


			// write per disease
			writeBins(outDirectory, mutationByDisease, germlineByDisease);

			// write per gene
			writeBins(outDirectory, mutationByGene, germlineByGene);
		}
	}
}
