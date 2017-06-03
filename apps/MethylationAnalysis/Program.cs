using System;
using System.Diagnostics;
using System.Collections.Generic;
using System.Linq;
using System.Threading.Tasks;
using System.Threading;
using System.IO;
using ASELib;

namespace MethylationAnalysis
{
	class MethylationPoint : IComparer<MethylationPoint>
	{
		public string compositeREF;
		public double mValue;
		public bool isSingle;

		public int Compare(MethylationPoint a, MethylationPoint b)
		{
			return xCompare(a, b);
		}

		static public int xCompare(MethylationPoint a, MethylationPoint b)
		{
			if (a.mValue > b.mValue) return 1;
			if (a.mValue < b.mValue) return -1;
			return 0;
		}
	}

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

		// lines to write to file
		static List<ASETools.OutputLine> outputLines = new List<ASETools.OutputLine>();

		// saves a list of Composite REFs from methylation files that overlap both 27k and 450k
		public static void processCompositeREFs(ASETools.Case case_27k, ASETools.Case case_450k)
		{
			ProcessFile(case_27k, true);
			ProcessFile(case_450k, true);

			// get REFs that are found in both the 27k anbd 450k files
			var overlap = casesMethylation.Where(r => r.Value.Count() > 1).Select(r => r.Key).ToList();
			Console.WriteLine(overlap.First());

			var intersectingREFs = new List<KeyValuePair<string, ASETools.GeneLocationInfo>>();

			// add remaining CompositeREFs to list of merged
			foreach (var compositeREF in overlap)
			{
				intersectingREFs.Add(new KeyValuePair<string, ASETools.GeneLocationInfo>(compositeREF, compositeREFs[compositeREF]));
			}
				
		
			Console.WriteLine(intersectingREFs.First());

			// write intersecting refs to file
			var outputFile = ASETools.CreateStreamWriterWithRetry(ASETools.ASEConfirguation.methylationREFsFilename);

			// Write header
			outputFile.WriteLine("CompositeREF\tChromosome\tPosition\tHugo Symbol");

			foreach (KeyValuePair<string, ASETools.GeneLocationInfo> entry in intersectingREFs)
			{
				outputFile.WriteLine(entry.Key + "\t" + entry.Value.chromosome + "\t" + entry.Value.minLocus + "\t" + entry.Value.hugoSymbol);
			}

			outputFile.Close();

		}


		static void ProcessFile(ASETools.Case case_, bool isTumor, string gene = null)
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
				if (gene != null && annotation.Gene_Symbol[0] != gene) {
					continue;
				}

				// filter out regions in promotor regions
				//Console.WriteLine(annotation.Position_to_TSS.Average());
				//if (annotation.Position_to_TSS.Average() < -2000 || annotation.Position_to_TSS.Average() > 2000)
				//{
				//	// methylation site is probably not in promotor regions, so we dont really know what effect it is going to have
				//	continue;
				//}

				if (!compositeREFs.ContainsKey(annotation.Composite_Element_REF))
				{
					// Composite REF does not overlap REF intersection. Throw it away
					continue;
				}

				Tuple<String, Double> d = new Tuple<String, Double>(methylation_file_id, annotation.Beta_Value);

				lock (casesMethylation)
				{
					if (!casesMethylation.ContainsKey(annotation.Composite_Element_REF))
					{
						// initialize file id to methylation
						casesMethylation.Add(annotation.Composite_Element_REF, new List<Tuple<string, double>>());
					}

					// append annotation to existing REF
					casesMethylation[annotation.Composite_Element_REF].Add(d);
				}

			}

			// read in maf file for counting mutations
			var mafLines = ASETools.MAFLine.ReadFile(case_.extracted_maf_lines_filename, case_.maf_file_id, false);

			// group MAFLines by hugo symbol to count mutations over symbol
			// TODO: this can be a function in ASETools
			var caseMutations = mafLines.GroupBy(r => r.Hugo_Symbol)
				.Select(group => new
				{
					Hugo_Symbol = group.Key,
					Count = group.Count()
				});

			// traverse through all mutations and save to corresponding hugo symbols
			foreach (var i in caseMutations)
			{
				if (gene != null && i.Hugo_Symbol != gene) {
					continue;
				}
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

		static Dictionary<string, int> arrayOverlap = new Dictionary<string, int>();



		static void ProcessCases(List<ASETools.Case> cases, string gene = null)
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
				if (case_.tumor_methylation_filename == "" || case_.extracted_maf_lines_filename == "" || case_.tumor_methylation_filename.Contains("HumanMethylation450"))
				{
					Console.WriteLine("No tumor methylation data for case " + case_.case_id + ". Skipping...");
				}
				else
				{
					// process only tumor files
					ProcessFile(case_, true);

					if (case_.normal_methylation_filename != "")
					{
						ProcessFile(case_, false);
					}
				}
			}
		}


		private static void MethylationMannWhitney(List<KeyValuePair<string, ASETools.GeneLocationInfo>> compositeREFs_) {

			outputLines = new List<ASETools.OutputLine>();

			while (true)
			{
				KeyValuePair<string, ASETools.GeneLocationInfo> compositeREF;

				lock (compositeREFs_)
				{
					if (compositeREFs_.Count() == 0)
					{
						//
						// No more work, we're done.
						//
						return;
					}

					Console.WriteLine(compositeREFs_.Count() + " remaining...");
					compositeREF = compositeREFs_[0];

					compositeREFs_.RemoveAt(0);
				}
				Console.WriteLine(Thread.CurrentThread.ToString() + ": " + compositeREF.Key);
				Dictionary<string, int> geneCounts = new Dictionary<string, int>();

				try
				{
					// get per gene information (tuples of file id, mutation count)
					geneCounts = mutationCounts[compositeREF.Value.hugoSymbol].ToDictionary(x => x.Item1, x => x.Item2);
				}
				catch (Exception)
				{
					// Gene does not have any mutations, and should not be included in analysis
					continue;

				}

				var methylationPoints = casesMethylation[compositeREF.Key].Select(r =>
				{
					var mutations = 0;
					try
					{
						mutations = geneCounts[r.Item1];
					}
					catch (Exception)
					{
						// no op
					}
					// Tuple of M-value, mutations
					return new Tuple<double, int>(r.Item2, mutations);
				}).Where(r => r.Item2 < 2).Select(r => {    // filter out more than 1 mutation and format to comparative MethylationPoints
					var m = new MethylationPoint();
					m.mValue = ASETools.AnnotationLine.betaToM(r.Item1);
					m.isSingle = r.Item2 == 0;
					m.compositeREF = compositeREF.Key;
					return m;
				}).ToList();

				// run MannWhitney test for this methylation point
				ASETools.MannWhitney<MethylationPoint>.GetValue getValue = new ASETools.MannWhitney<MethylationPoint>.GetValue(m => m.mValue);
				ASETools.MannWhitney<MethylationPoint>.WhichGroup whichGroup = new ASETools.MannWhitney<MethylationPoint>.WhichGroup(m => m.isSingle);

				double nSingle = 0;
				double nMultiple = 0;
				double p;
				bool reversed;
				double U;
				double z;

				bool enoughData;

				p = ASETools.MannWhitney<MethylationPoint>.ComputeMannWhitney(methylationPoints, methylationPoints[0], whichGroup, getValue, out enoughData, out reversed, out nSingle, out nMultiple, out U, out z);
				if (!enoughData)
				{
					continue;
				}
				
				var outputLine = new ASETools.OutputLine();
				outputLine.line = compositeREF.Value.hugoSymbol + "\t" + compositeREF.Key + "\t" + nSingle + "\t" + nMultiple + "\t" + U + "\t" + z;
				outputLine.p = p;
				outputLines.Add(outputLine);

			}
		}

		static void MethylationMannWhitney()
		{
			// total number of REFs for Bonferroni correction
			int compositeREFCount = compositeREFs.Keys.Count();

			// TODO: this doesnt work
			//var threads = new List<Thread>();
			//for (int i = 0; i < Environment.ProcessorCount; i++)
			//{
			//	threads.Add(new Thread(() => MethylationMannWhitney(compositeREFs.ToList())));
			//}

			//threads.ForEach(th => th.Start());
			//threads.ForEach(th => th.Join());

			MethylationMannWhitney(compositeREFs.ToList());


			// open file
			// TODO remove hard coded
			var output = new StreamWriter(@"C:\Users\t-almorr\temp\MannWhitney_with450.txt");

			// write header
			output.WriteLine("Gene\tCompositeREF\tnSingle\tnMultiple\tU\tz\tp\tp_Bonferroni\t");

			var compositeCount = outputLines.Count();

			foreach (var outputLine in outputLines)
			{
				output.WriteLine(outputLine.line + "\t" + outputLine.p + "\t" + outputLine.p * compositeCount);
			}

			output.Close();
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

			var selectedCases = cases.Select(kv => kv.Value).ToList();

			// Compute file for all overlapping REFs
			// TODO move this to process manager
			// var case_27 = cases.Where(r => r.Value.tumor_methylation_filename.Contains("HumanMethylation27")).First().Value;
			// var case_450 = cases.Where(r => r.Value.tumor_methylation_filename.Contains("HumanMethylation450")).First().Value;
			// processCompositeREFs(case_27, case_450);

			// read in file of valid Composite REFs for methylation data
			compositeREFs = ASETools.CompositeREFInfoLine.ReadFile(ASETools.ASEConfirguation.methylationREFsFilename).ToDictionary(x => x.Key, x => x.Value);

			var threads = new List<Thread>();
			for (int i = 0; i < Environment.ProcessorCount; i++)
			{
				threads.Add(new Thread(() => ProcessCases(selectedCases)));
			}

			threads.ForEach(th => th.Start());
			threads.ForEach(th => th.Join());


			//MethylationMannWhitney();

			var targetDisease = "all";

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

			// Write metadata for cases
			var caseFilename = @"C:\Users\t-almorr\temp\CASEMETADATA_" + targetDisease + ".txt";
			var outputFile = ASETools.CreateStreamWriterWithRetry(caseFilename);

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
