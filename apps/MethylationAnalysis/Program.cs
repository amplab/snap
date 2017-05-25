using System;
using System.Diagnostics;
using System.Collections.Generic;
using System.Linq;
using System.Threading.Tasks;
using System.Threading;
using ASELib;

namespace MethylationAnalysis
{

	class Program
	{

		// Dictionary storing methylation REF, where key is Composite Ref and Value is (methylation_file_id, M value) tuple
		static Dictionary<string, List<Tuple<string, Double>>> casesMethylation = new Dictionary<string, List<Tuple<string, Double>>>();

		// Dictionary storing Composite REF information. Key is Composite REF, for which we store GeneLocationInfo tuples
		static Dictionary<string, ASETools.GeneLocationInfo> compositeREFs = new Dictionary<string, ASETools.GeneLocationInfo>();

		// Dictionary storing case data. Keys are methylation file ids and values are (disease, label) tuples
		static Dictionary<string, Tuple<string, string>> fileInfo = new Dictionary<string, Tuple<string, string>>();

		// dictionary of hugo symbol, (methylation file id, mutation counts) tuples
		static Dictionary<string, List<Tuple<string, int>>> mutationCounts = new Dictionary<string, List<Tuple<string, int>>>();

		static void ProcessFile(ASETools.Case case_, bool isTumor)
		{
			// select metadata whether processing normal or tumor methylation data
			var label = isTumor?"Tumor":"Normal";
			var methylation_filename = isTumor ? case_.tumor_methylation_filename : case_.normal_methylation_filename;
			var methylation_file_id = isTumor ? case_.tumor_methylation_file_id : case_.normal_methylation_file_id;

			// read in tumor methylation file 
			var annotations = ASELib.ASETools.AnnotationLine.ReadFile(methylation_filename, methylation_file_id, false);

			// add case to caseInfo
			if (fileInfo.ContainsKey(methylation_file_id))
			{
				Console.WriteLine("Error: already processed methylation file " + methylation_file_id + " for case " + case_.case_id);
				return;
			}

			fileInfo.Add(methylation_file_id, new Tuple<string, string>(case_.disease(), label));

			foreach (var annotation in annotations)
			{
				Tuple<String, Double> d = new Tuple<String, Double>(methylation_file_id, annotation.Beta_Value);

				lock (compositeREFs)
				{
					if (!compositeREFs.ContainsKey(annotation.Composite_Element_REF))
					{
						// initialize file id to methylation
						casesMethylation.Add(annotation.Composite_Element_REF, new List<Tuple<string, double>>());

						// Add to seen Composite REFs
						compositeREFs.Add(annotation.Composite_Element_REF, new ASETools.GeneLocationInfo());
						compositeREFs[annotation.Composite_Element_REF].chromosome = annotation.Chromsome;
						compositeREFs[annotation.Composite_Element_REF].minLocus = annotation.Start;
						compositeREFs[annotation.Composite_Element_REF].hugoSymbol = annotation.Gene_Symbol[0];
					}

					// append annotation to existing REF
					casesMethylation[annotation.Composite_Element_REF].Add(d);
				}
			}

			// read in maf file for counting mutations
			var mafLines = ASETools.MAFLine.ReadFile(case_.extracted_maf_lines_filename, case_.maf_file_id, false);

			// group MAFLines by hugo symbol to count mutations over symbol
			// TODO: this can be a function in ASETools
			// TODO: is this counting right?
			var caseMutations = mafLines.GroupBy(r => r.Hugo_Symbol)
				.Select(group => new
				{
					Hugo_Symbol = group.Key,
					Count = group.Count()
				});

			// traverse through all mutations and save to corresponding hugo symbols
			foreach (var i in caseMutations)
			{
				lock (mutationCounts)
				{
					if (!mutationCounts.ContainsKey(i.Hugo_Symbol))
					{
						mutationCounts.Add(i.Hugo_Symbol, new List<Tuple<string, int>>());
					}

					mutationCounts[i.Hugo_Symbol].Add(new Tuple<string, int>(methylation_file_id, i.Count));
				}

			}
		}

		static void ProcessCases(List<ASETools.Case> cases, Boolean only27)
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

				// verify methylation file exists
				if (case_.tumor_methylation_filename == "" || case_.extracted_maf_lines_filename == "")
				{
					Console.WriteLine("No tumor methylation data for case " + case_.case_id + ". Skipping...");
				}
				else
				{
					if (only27 && case_.tumor_methylation_filename.Contains("HumanMethylation27") || !only27)
						ProcessFile(case_, true);
					else
						Console.WriteLine("450k methylation for " + case_.tumor_methylation_filename + " Skipping.");

					if (case_.normal_methylation_filename != "")
					{
						Console.WriteLine("Both tumor and normal methylation data exists for case " + case_.case_id);
						if (only27 && case_.tumor_methylation_filename.Contains("HumanMethylation27") || !only27) {
							// no op. Do not process normal
							//	ProcessFile(case_, false);
						}
						else
							Console.WriteLine("450k methylation for " + case_.tumor_methylation_filename + " Skipping.");
					}

				}
			}
		}

		static int Main(string[] args)
		{
			var timer = new Stopwatch();
			timer.Start();

			var configuration = ASETools.ASEConfirguation.loadFromFile(args);
			var cases = ASETools.Case.LoadCases(configuration.casesFilePathname);

			if (null == cases)
			{
				Console.WriteLine("You must have selected cases before you can analyze methylation data.");
				return 1;
			}

			var targetDisease = "all";

			var selectedCases = cases.Select(kv => kv.Value).ToList();

			Console.WriteLine("Processing " + selectedCases.Count() + " cases");


			var threads = new List<Thread>();
			for (int i = 0; i < Environment.ProcessorCount; i++)
			{
				threads.Add(new Thread(() => ProcessCases(selectedCases, false)));
			}

			threads.ForEach(th => th.Start());
			threads.ForEach(th => th.Join());


			// Write file for B values and methylation counts
			var bValueFilename = @"C:\Users\t-almorr\temp\BVALUEMATRIX_" + targetDisease + ".txt";
			var mValueFilename = @"C:\Users\t-almorr\temp\MVALUEMATRIX_" + targetDisease + ".txt";

			var countMutationFilename = @"C:\Users\t-almorr\temp\COUNTMUTATIONMATRIX_" + targetDisease + ".txt";

			var bOutputFile = ASETools.CreateStreamWriterWithRetry(bValueFilename);
			var mOutputFile = ASETools.CreateStreamWriterWithRetry(mValueFilename);

			// Write header: all file ids
			var header = "Composite_REF\t" + String.Join("\t", fileInfo.Keys);
			bOutputFile.WriteLine(header);
			mOutputFile.WriteLine(header);

			// Stores final Composite REFs that were present in all samples
			var finalCompositeREFS = new List<string>();

			// write out caseMethylation lines to file
			foreach (KeyValuePair<string, List<Tuple<string, double>>> entry in casesMethylation)
			{
				// Filter out Gene Symbols that do not have data for all selected cases
				if (entry.Value.Count != fileInfo.Keys.Count) {
					Console.WriteLine("Don't have consistent values for gene " + compositeREFs[entry.Key].hugoSymbol + " for all files. Skipping...");
					continue;
				}

				// If this REF has data for all samples, add it to list of methylation sites to be processed
				finalCompositeREFS.Add(entry.Key);

				// write all Beta-values for this composite REF
				bOutputFile.WriteLine(entry.Key + "\t" + String.Join("\t", entry.Value.Select(r => r.Item2)));

				// write all m values for this composite REF
				mOutputFile.WriteLine(entry.Key + "\t" + String.Join("\t", entry.Value.Select(r => ASETools.AnnotationLine.betaToM(r.Item2))));

			}

			bOutputFile.Close();
			mOutputFile.Close();

			// TODO this is gonna be soo SLOW
			// write out mutations
			var countOutputFile = ASETools.CreateStreamWriterWithRetry(countMutationFilename);
			var countHeader = "Hugo_Symbol\t" + String.Join("\t", fileInfo.Keys );
			countOutputFile.WriteLine(header);

			// for each hugo symbol, get mutation counts for all cases
			foreach (KeyValuePair<string, List<Tuple<string, int>>> entry in mutationCounts)
			{

				// get case ids and set initial mutation counts to 0
				var countByCase = fileInfo.Keys.ToDictionary(x => x, x => 0);

				foreach (var v in entry.Value)
				{
					if (countByCase.ContainsKey(v.Item1))
					{
						countByCase[v.Item1] = v.Item2;
					}
				}
				
				// write all mutation counts for this Hugo Symbol
				countOutputFile.WriteLine(entry.Key + "\t" + String.Join("\t", countByCase.Values.ToList()));

			}

			countOutputFile.Close();


			// Write metadata for composite REFS
			var refFilename = @"C:\Users\t-almorr\temp\REFMETADATA_" + targetDisease + ".txt";
			var outputFile = ASETools.CreateStreamWriterWithRetry(refFilename);

			// Write header
			outputFile.WriteLine("Composite_REF\tChromosome\tPosition\tGene_Symbol");

			foreach(KeyValuePair<string, ASETools.GeneLocationInfo> entry in compositeREFs)
			{
				if (finalCompositeREFS.Contains(entry.Key))
				{
					outputFile.WriteLine(entry.Key + "\t" + entry.Value.chromosome + "\t" + entry.Value.minLocus + "\t" + entry.Value.hugoSymbol);
				}
			}

			outputFile.Close();

			// Write metadata for cases
			var caseFilename = @"C:\Users\t-almorr\temp\CASEMETADATA_" + targetDisease + ".txt";
			outputFile = ASETools.CreateStreamWriterWithRetry(caseFilename);

			// Write header
			outputFile.WriteLine("File_Id\tdisease\tType");

			foreach (KeyValuePair<string, Tuple<string, string>> entry in fileInfo)
			{
				outputFile.WriteLine(entry.Key + "\t" + entry.Value.Item1 + "\t" + entry.Value.Item2);
			}

			outputFile.Close();

			return 0;
		}
	}
}
