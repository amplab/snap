using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using System.Threading;
using ASELib;

namespace RegionalSignalByDistance
{

	class Program
	{

		// stores metadata for each sample such as mutation count and disease
		// Dictionary<case id, Tuple<mutation count, disease>>
		static Dictionary<string, Tuple<int, string>> metaData = new Dictionary<string, Tuple<int, string>>();

		// stores average ase values for each bin size
		// Dictionary<case id, Dictionary<distance, Tuple<sum, count>>
		static Dictionary<string, Dictionary<int, Tuple<double, int>>> aseValues 
			= new Dictionary<string, Dictionary<int, Tuple<double, int>>>();

		static int binSize;

		public static void ProcessCases(List<ASETools.Case> cases, ASETools.GeneLocationInfo gene)
		{
			ASETools.Case case_;
			while (true)
			{
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

				if (case_.annotated_selected_variants_filename == "")
				{
					continue;
				}

				// Process case
				// Calculate ASE for all Variants and store by distance
				var annotatedVariants = ASETools.AnnotatedVariant.readFile(case_.annotated_selected_variants_filename)
					.Where(r => r.contig == gene.chromosome);

				var mutationCount = ASETools.MAFLine.ReadFile(case_.extracted_maf_lines_filename, case_.maf_file_id, false)
					.Where(r => r.Hugo_Symbol == gene.hugoSymbol).Count(); 

				// add information of disease type and mutation count
				lock (metaData)
				{
					metaData.Add(case_.case_id, new Tuple<int, string>(mutationCount, case_.disease()));
				}

				lock (aseValues)
				{
					aseValues.Add(case_.case_id, new Dictionary<int, Tuple<double, int>>());
				}

				foreach (var annotatedVariant in annotatedVariants)
				{
					// calculate ASE for this point
					if (annotatedVariant.tumorDNAReadCounts.nMatchingReference + annotatedVariant.tumorDNAReadCounts.nMatchingAlt >= 10 &&
						annotatedVariant.tumorRNAReadCounts.nMatchingReference + annotatedVariant.tumorRNAReadCounts.nMatchingAlt >= 10 &&
						annotatedVariant.tumorDNAReadCounts.nMatchingReference * 3 >= annotatedVariant.tumorDNAReadCounts.nMatchingAlt * 2 &&
						annotatedVariant.tumorDNAReadCounts.nMatchingAlt * 3 >= annotatedVariant.tumorDNAReadCounts.nMatchingReference * 2)
					{
						double rnaFractionTumor = (double)annotatedVariant.tumorRNAReadCounts.nMatchingAlt / (annotatedVariant.tumorRNAReadCounts.nMatchingReference + annotatedVariant.tumorRNAReadCounts.nMatchingAlt);

						var alleleSpecificExpression = Math.Abs(rnaFractionTumor * 2.0 - 1.0);

						// get bin as a function of binSize 
						int bin = annotatedVariant.locus - (annotatedVariant.locus % binSize); // TODO Verify 

						// bin distance by binsize and add to dictionary
						lock (aseValues[case_.case_id])
						{
							Tuple<double, int> runningSum;
							if (aseValues[case_.case_id].TryGetValue(bin, out runningSum))
							{
								aseValues[case_.case_id][bin] = new Tuple<double, int>(alleleSpecificExpression + runningSum.Item1,
									runningSum.Item2 + 1);
							} else
							{
								aseValues[case_.case_id].Add(bin, new Tuple<double, int>(alleleSpecificExpression, 1));
							}
						}

					}

				}

			}
		}

		static void PrintMessage()
		{
			Console.WriteLine("Usage: RegionalSignalByDistance hugoSymbol binSize");
		}


		static void Main(string[] args)
		{
			var configuration = ASETools.Configuration.loadFromFile(args);

			if (configuration.commandLineArgs.Count() != 2)
			{
				PrintMessage();
				return;
			}

			var hugoSymbol = configuration.commandLineArgs[0];
			binSize = Convert.ToInt32(configuration.commandLineArgs[1]);

			var knownGenes = ASETools.readKnownGeneFile(configuration.geneLocationInformationFilename)
				.Where(r => !r.Value.inconsistent).ToList(); 

			var matchingGenes = knownGenes.Where(r => r.Value.hugoSymbol == hugoSymbol);

			if (matchingGenes.Count() == 0)
			{
				Console.WriteLine("Invalid Hugo Symbol " + hugoSymbol);
				return;
			}

			var thisGene = matchingGenes.First();
			// get other genes on this chromosome
			var genesOnSameChr = knownGenes.Where(r => r.Value.chromosome == thisGene.Value.chromosome);

			var cases = ASETools.Case.LoadCases(configuration.casesFilePathname);

			var selectedCases = cases.Select(r => r.Value).ToList();
			var threads = new List<Thread>();
			for (int i = 0; i < Environment.ProcessorCount; i++)
			{
				threads.Add(new Thread(() => ProcessCases(selectedCases, thisGene.Value)));
			}

			threads.ForEach(t => t.Start());
			threads.ForEach(t => t.Join());


			// open file 
			var outputFilename = configuration.finalResultsDirectory + hugoSymbol + @"_alleleSpecificExpressionByDistance.txt";
			var writer = ASETools.CreateStreamWriterWithRetry(outputFilename);

			// write header
			writer.Write("MutationCount\t");
			writer.Write("Disease");

			// get max bin for this chromosome
			int maxBin = aseValues.Select(r => r.Value).SelectMany(r => r.Keys).Max();


			// write out genes
			for (var i = 0; i <= maxBin; i += binSize)
			{
				// if there are genes that overlap this bin, then print them out
				var overlappingGenes = genesOnSameChr.Where(r => r.Value.overlapsRange(r.Value.chromosome, i, i + binSize))
					.Select(r => r.Value.hugoSymbol);

				if (matchingGenes.Count() > 0)
				{
					writer.Write("\t" + String.Join(",", overlappingGenes));
				} else {
					writer.Write("\t");
				}

			}
			writer.WriteLine();

			// write out location. Two empty spots for Mutation and Disease
			writer.Write("\t");

			for (var i = 0; i <= maxBin; i += binSize)
			{
				writer.Write("\t" + i);
			}
			writer.WriteLine();

			foreach (var case_ in cases)
			{
				Tuple<int, string> meta;
				Dictionary<int, Tuple<double, int>> values;

				if (!metaData.TryGetValue(case_.Key, out meta) || !aseValues.TryGetValue(case_.Key, out values))
				{
					// no data for this case
					continue;
				}

				// write out this case
				writer.Write(meta.Item1 + "\t" + meta.Item2);

				for (var i = 0; i <= maxBin; i += binSize)
				{
					Tuple<double, int> binnedValue;
					if (values.TryGetValue(i, out binnedValue))
					{
						writer.Write("\t" + binnedValue.Item1 / binnedValue.Item2);
					}
					else {
						// no values for this bin
						writer.Write("\tNaN");
					}
				} // foreach bin
				writer.WriteLine();
			} // foreach case

			writer.Close();

		}
	}
}
