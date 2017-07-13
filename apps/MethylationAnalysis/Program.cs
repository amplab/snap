using System;
using System.Diagnostics;
using System.Collections.Generic;
using System.Linq;
using System.Threading.Tasks;
using System.Threading;
using System.IO;
using ASELib;
using ExpressionNearMutations;
using GenomeBuild;

namespace MethylationAnalysis

{

	class Program
	{
		static ASETools.GeneLocationsByNameAndChromosome geneLocationInformation;

		static Dictionary<string, List<GenomeBuild.Interval>> enhancers = new Dictionary<string, List<GenomeBuild.Interval>>();

		static void ProcessFile(ASETools.Case case_, bool forTumor)
		{
			// select metadata whether processing normal or tumor methylation data

			var methylation_filename = forTumor ? case_.tumor_methylation_filename : case_.normal_methylation_filename;
			var methylation_file_id = forTumor ? case_.tumor_methylation_file_id : case_.normal_methylation_file_id;

			// Load MAF file for this case
			var mafLines = ASETools.MAFLine.ReadFile(case_.extracted_maf_lines_filename, case_.maf_file_id, false);

		    Dictionary <string, ASETools.GeneExpression> geneExpressions = new Dictionary<string, ASETools.GeneExpression>();

			// 2. for each gene in the genelocations, count the number of mutations and get the methylation values for this gene. These will be 
			// used later.
			// for each gene, count the methylation averages as you get farther and farther away from the gene (inc, chr, autosome)
			foreach (var geneLocation in geneLocationInformation.genesByName.Values)
			{
				if (!geneExpressions.ContainsKey(geneLocation.hugoSymbol))
				{
					geneExpressions.Add(geneLocation.hugoSymbol, new ASETools.GeneExpression(geneLocation));
				}

				// get mutation counts from MAFLines
				var mutationCount = mafLines.Where(r => r.Hugo_Symbol == geneLocation.hugoSymbol).Count();
				geneExpressions[geneLocation.hugoSymbol].mutationCount += mutationCount;
			}
			var expressionValues = new ExpressionNearMutations.Program.RegionalSignal(geneExpressions, geneLocationInformation, false);

			// read in methylation file 
			var annotations = ASETools.AnnotationLine.ReadFile(methylation_filename, methylation_file_id, false);
			var initialCount = annotations.Count();
			var count = 0;

			var timer = new Stopwatch();
			timer.Start();

			// filter out annotations that are not by an enhancer or promotor
			var filteredAnnotations = new List<ASETools.AnnotationLine>();
			var filteredCount = 0;
			foreach (var r in annotations)
			{
				filteredCount++;
				List<GenomeBuild.Interval> chrEnhancers;
				enhancers.TryGetValue(ASETools.chromosomeNameToNonChrForm(r.compositeREF.Chromosome), out chrEnhancers);
				if (chrEnhancers != null)
				{
					if (r.compositeREF.Position_to_TSS.Average() <= 2000 ||
					chrEnhancers.Where(x => x.overlaps(new GenomeBuild.Interval(ASETools.chromosomeNameToNonChrForm(r.compositeREF.Chromosome), r.compositeREF.Start, r.compositeREF.End))).Count() > 0
					)
					{
						filteredAnnotations.Add(r);
					}

				}
				else if (r.compositeREF.Position_to_TSS.Average() <= 2000)
				{
					filteredAnnotations.Add(r);
				}
			}


			timer.Stop();

			Console.WriteLine("Lost " + (initialCount - filteredAnnotations.Count()) + " annotations for case" + case_.case_id);

			foreach (var annotation in filteredAnnotations)
			{
				count++;

				if (count % 100000 == 0 && false)
				{
					Console.WriteLine("processed " + count + "annotations out of " + filteredAnnotations.Count());
				}

				// get non chr form
				var chromosome = ASETools.chromosomeNameToNonChrForm(annotation.compositeREF.Chromosome);

				// filter out regions in promotor and enhancer regions
				expressionValues.Add(annotation.compositeREF.Chromosome, annotation.compositeREF.Start, annotation.M_Value, 0, false);

			}

			// Save values

			var allExpressions = new List<ASETools.GeneExpression>();
			foreach (var expressionEntry in expressionValues.geneExpressions)
			{
				allExpressions.Add(expressionEntry.Value);
			}

			allExpressions.Sort(ASETools.GeneExpression.CompareByGeneName);

			string directory = ASETools.GetDirectoryPathFromFullyQualifiedFilename(case_.extracted_maf_lines_filename);

			var extension = forTumor ? ASETools.tumorRegionalMethylationExtension : ASETools.normalRegionalMethylationExtension;
			var outputFilename = directory + case_.case_id + extension;
			Console.WriteLine("Saving to file " + outputFilename);

			var fileHeader = "Methylation Analysis " + case_.case_id;
			var output = new ASETools.RegionalSignalFile(fileHeader,
					allExpressions,
					expressionValues.wholeAutosomeRegionalExpression,
					expressionValues.allButThisChromosomeAutosomalRegionalExpressionState,
					expressionValues.perChromosomeRegionalExpressionState);

			var columnSuffix = forTumor ? "tumor M-val" : "normal M-val";
			var printMu = false;
			var isTumor = true;
			output.WriteFile(outputFilename, printMu, 1, columnSuffix, isTumor);
		}


		static void ProcessCases(List<ASETools.Case> cases)
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

					Console.WriteLine(cases.Count() + " remaining...");

					case_ = cases[0];

					cases.RemoveAt(0);
				}


				// verify methylation file exists for tumor data
				if (!case_.tumor_methylation_filename.Contains("HumanMethylation450") || case_.extracted_maf_lines_filename == "" || case_.tumor_rna_filename == "")
				{
					Console.WriteLine("No tumor methylation data for case " + case_.case_id + ". Skipping...");
				}
				else
				{
					ProcessFile(case_, true);
				}

				// verify methylation file exists for normal data
				if (!case_.normal_methylation_filename.Contains("HumanMethylation450") || case_.extracted_maf_lines_filename == "" || case_.normal_rna_filename == "")
				{
					Console.WriteLine("No normal methylation data for case " + case_.case_id + ". Skipping...");
				}
				else
				{
					ProcessFile(case_, false);
				}

			}
		}

		static void saveEnhancers()
		{
			Console.WriteLine("saving enhancers");

			// hg19 to hg38 build
			var build = new GenomeBuild.LiftOver();
			build.readChainFile(ASETools.ASEConfirguation.hg19Tohg38ChainFile);

			// read in enhancers
			string[] fileEntries = Directory.GetFiles(ASETools.ASEConfirguation.defaultBaseDirectory + "enhancers");
			List<GenomeBuild.Interval> enhancers = new List<GenomeBuild.Interval>();

			foreach (string fileName in fileEntries)
			{
				var bedLines = ASETools.BedLine.ReadFile(fileName)
					.Select(r => new GenomeBuild.Interval(r.Chromosome, r.Start_Position, r.End_Position)).ToList();
				bedLines.Sort();
				enhancers.AddRange(GenomeBuild.Interval.collapse(bedLines));
				Console.WriteLine("Finished file " + fileName);
			}
			Console.WriteLine(enhancers.Count() + " enhancers");

			enhancers = enhancers.SelectMany(r => build.mapCoordinates(r))
				.Select(r => new GenomeBuild.Interval(ASETools.chromosomeNameToNonChrForm(r.name), r.start, r.end)).ToList();
			Console.WriteLine(enhancers.Count() + "enhancers after running genome build");

			var collapsed = GenomeBuild.Interval.collapse(enhancers);
			Console.WriteLine(enhancers.Count() + "enhancers after running final collapse");

			var enhancerLines = collapsed.Select(r => new ASETools.BedLine(r.name, Convert.ToInt32(r.start), Convert.ToInt32(r.end), "enhancer", 1000, r.strand, r.start, r.end)).ToList();
			var enhancerFile = ASETools.ASEConfirguation.defaultBaseDirectory + @"enhancers\allEnhancers.bed";
			ASETools.BedLine.writeFile(enhancerFile, enhancerLines);
		}

		static void PrintUsageMessage()
		{
			Console.WriteLine("usage: MethylationAnalysis -450 casesToProcess");
		}

		static void Main(string[] args)
		{

			var timer = new Stopwatch();
			timer.Start();

			var configuration = ASETools.ASEConfirguation.loadFromFile(args);
			enhancers = ASETools.BedLine.ReadFile(ASETools.ASEConfirguation.defaultBaseDirectory + @"enhancers\allEnhancers.bed")
					.Select(r => new GenomeBuild.Interval(r.Chromosome, r.Start_Position, r.End_Position, r.strand)).GroupBy(r => r.name)
					.ToDictionary(x => x.Key, x => x.Select(r => r).ToList());

			if (configuration.commandLineArgs.Count() < 1)
			{
				PrintUsageMessage();
				return;
			}

			int argsConsumed = 0;
			bool onlyProcess450 = false;

			if (configuration.commandLineArgs[0] == "-450") {
				onlyProcess450 = true;
				argsConsumed += 1;
			}

			// lines to print out
			var lines = new List<string>();

			var cases = ASETools.Case.LoadCases(configuration.casesFilePathname);
			
			if (null == cases)
			{
				Console.WriteLine("You must have selected cases before you can analyze methylation data.");
				return;
			}

			var selectedCases = new List<ASETools.Case>();
			for (int i = argsConsumed; i < configuration.commandLineArgs.Count(); i++)
			{
				if (!cases.ContainsKey(configuration.commandLineArgs[i]))
				{
					Console.WriteLine(configuration.commandLineArgs[i] + " is not a valid case ID. Skipping...");
					continue;
				}
				selectedCases.Add(cases[configuration.commandLineArgs[i]]);
			}

			// Get information for current genome build in chr form
			geneLocationInformation = new ASETools.GeneLocationsByNameAndChromosome(ASETools.readKnownGeneFile(ASETools.ASEConfirguation.defaultGeneLocationInformation));

			var threads = new List<Thread>();
			for (int i = 0; i < Environment.ProcessorCount; i++)
			{
				threads.Add(new Thread(() => ProcessCases(selectedCases)));
			}

			threads.ForEach(th => th.Start());
			threads.ForEach(th => th.Join());


			return;

		}
	}
}
