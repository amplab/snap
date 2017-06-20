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

		// Dictionary storing methylation REF, where key is Composite Ref and Value is (case_id, M value) tuple
		static Dictionary<string, List<Tuple<string, Double>>> casesMethylation = new Dictionary<string, List<Tuple<string, Double>>>();

		// Dictionary storing Composite REF information. Key is Composite REF, for which we store GeneLocationInfo tuples
		static Dictionary<string, ASETools.GeneLocationInfo> compositeREFs = new Dictionary<string, ASETools.GeneLocationInfo>();

		// Dictionary storing case data. Keys are case ids and values are (disease, label, arrayType) tuples
		static Dictionary<string, Tuple<string, string, int>> fileInfo = new Dictionary<string, Tuple<string, string, int>>();

		// dictionary of hugo symbol, (case id, mutation counts) tuples
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
			outputFile.WriteLine("CompositeREF\tChromosome\tPosition\tHugo Symbol\tPositionToTSS");

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

			var arrayType = methylation_filename.Contains("HumanMethylation450") ? 1 : 0;

			// read in tumor methylation file 
			var annotations = ASELib.ASETools.AnnotationLine.ReadFile(methylation_filename, methylation_file_id, false);

			// add case to caseInfo
			if (fileInfo.ContainsKey(case_.case_id))
			{
				Console.WriteLine("Error: already processed methylation file " + methylation_file_id + " for case " + case_.case_id);
				return;
			}

			fileInfo.Add(case_.case_id, new Tuple<string, string, int>(case_.disease(), label, arrayType));

			foreach (var annotation in annotations)
			{
				if (gene != null && annotation.Gene_Symbol[0] != gene) {
					continue;
				}

				if (!compositeREFs.ContainsKey(annotation.Composite_Element_REF))
				{
					// Composite REF does not overlap REF intersection. Throw it away
					continue;
				}

				Tuple<String, Double> d = new Tuple<String, Double>(case_.case_id, annotation.Beta_Value);

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

					mutationCounts[i.Hugo_Symbol].Add(new Tuple<string, int>(case_.case_id, i.Count));
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
				if (case_.tumor_methylation_filename == "" || case_.extracted_maf_lines_filename == "")
				{
					Console.WriteLine("No tumor methylation data for case " + case_.case_id + ". Skipping...");
				}
				else
				{
					// process only tumor files
					ProcessFile(case_, true, gene);

					if (case_.normal_methylation_filename != "")
					{
						//ProcessFile(case_, false, gene);
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

			MethylationMannWhitney(compositeREFs.ToList());


			// open file
			var output = new StreamWriter(@"\\msr-genomics-0\d$\gdc\methyl_temp\MannWhitney_with450.txt");

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
			var bisulfiteCases = new List<string>();

			bisulfiteCases.Add("6fd72426-f6c8-47ca-a500-d5d3600b9b15");
			bisulfiteCases.Add("9d95a65b-e41d-4f93-92d4-99dce29ff40d");
			bisulfiteCases.Add("fba80122-d8b2-4d8d-a032-9767e8160f9f");
			bisulfiteCases.Add("c5c4a0a5-900d-483d-9282-475654d63265");

			// read in maf file for counting mutations
			foreach (var case_ in cases) {

				if (bisulfiteCases.Contains(case_.Value.case_id)) {
					var mafLines = ASETools.MAFLine.ReadFile(case_.Value. extracted_maf_lines_filename, case_.Value.maf_file_id, false);

					var mut = mafLines.Where(r => r.Hugo_Symbol == "TP53").ToList();

					Console.WriteLine("case " + case_.Value.case_id + " found tp53" + mafLines.Where(r => r.Hugo_Symbol == "TP53").Count() + "times");
				}

			}



			var selectedCases = cases.Select(kv => kv.Value).ToList();

			var methyl_dir = @"\\msr-genomics-0\d$\gdc\methyl_temp\";

			// read in file of valid Composite REFs for methylation data
			// This was for when we were using 27k and 450k data. for now, we will just use 450k data
			//compositeREFs = ASETools.CompositeREFInfoLine.ReadFile(ASETools.ASEConfirguation.methylationREFsFilename).ToDictionary(x => x.Key, x => x.Value);

			var threads = new List<Thread>();
			for (int i = 0; i < Environment.ProcessorCount; i++)
			{
				threads.Add(new Thread(() => ProcessCases(selectedCases)));
			}

			threads.ForEach(th => th.Start());
			threads.ForEach(th => th.Join());

			//MethylationMannWhitney();

			var targetDisease = "all";

			// Write file for B and M values and methylation counts
			var bValueFilename =  methyl_dir + "BVALUEMATRIX_" + targetDisease + ".txt";
			var mValueFilename = methyl_dir + "MVALUEMATRIX_" + targetDisease + ".txt";

			var countMutationFilename = methyl_dir + "COUNTMUTATIONMATRIX_" + targetDisease + ".txt";

			var bOutputFile = ASETools.CreateStreamWriterWithRetry(bValueFilename);
			var mOutputFile = ASETools.CreateStreamWriterWithRetry(mValueFilename);

			// Write header: all file ids
			var caseIds = fileInfo.Keys.ToList();
			caseIds.Sort();

			var header = "Composite_REF\t" + String.Join("\t", caseIds);
			bOutputFile.WriteLine(header);
			mOutputFile.WriteLine(header);

			// write out caseMethylation lines to file
			foreach (KeyValuePair<string, List<Tuple<string, double>>> entry in casesMethylation)
			{

				// write all Beta-values for this composite REF
				bOutputFile.Write(entry.Key);
				mOutputFile.Write(entry.Key);

				// dictionary of (case id, beta values) for this REF
				var casesForREF = entry.Value.ToDictionary(x => x.Item1, x => x.Item2);

				foreach (var case_id in caseIds)
				{
					// if case was found for this ref, put in value. Otherwise, assign 'NA'.
					if (casesForREF.ContainsKey(case_id)) {
						bOutputFile.Write("\t" + casesForREF[case_id]);
						mOutputFile.Write("\t" + ASETools.AnnotationLine.betaToM(casesForREF[case_id]));
					} else {
						bOutputFile.Write("\tna");
						mOutputFile.Write("\tna");
					}
				}

				bOutputFile.WriteLine();
				mOutputFile.WriteLine();

			}

			bOutputFile.Close();
			mOutputFile.Close();

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
			var caseFilename = methyl_dir + "CASEMETADATA_" + targetDisease + ".txt";
			var outputFile = ASETools.CreateStreamWriterWithRetry(caseFilename);

			// Write header
			outputFile.WriteLine("File_Id\tdisease\tType\tis450");

			foreach (KeyValuePair<string, Tuple<string, string, int>> entry in fileInfo)
			{
				outputFile.WriteLine(entry.Key + "\t" + entry.Value.Item1 + "\t" + entry.Value.Item2 + "\t" + entry.Value.Item3);
			}

			outputFile.Close();

			return 0;
		}
	}
}
