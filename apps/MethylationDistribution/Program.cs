using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading;
using ASELib;

namespace MethylationDistribution
{
	class Program
	{

		// keeps track of the tumor and normal methylation distributions
		static Dictionary<double, int> tumorMethylationDistribution = new Dictionary<double, int>();
		static Dictionary<double, int> normalMethylationDistribution = new Dictionary<double, int>();

		static void ProcessFile(ASETools.Case case_)
		{
			// Process tumor methylation data
			var annotations = ASETools.AnnotationLine.ReadFile(case_.tumor_methylation_filename, case_.tumor_methylation_file_id, false);

			foreach (var annotation in annotations)
			{
				var index = Math.Round(annotation.M_Value, 1); // 0.1, 0.2, ...
				int bin;
				lock (tumorMethylationDistribution)
				{
					tumorMethylationDistribution.TryGetValue(index, out bin);
					if (bin == 0)
					{
						tumorMethylationDistribution.Add(index, 1);
					}
					else
					{
						tumorMethylationDistribution[index] += 1;
					}
				}
			}

			// Process normal methylation data, if it exists
			if (case_.normal_methylation_filename.Contains("HumanMethylation450"))
			{
				// Process tumor methylation data
				annotations = ASETools.AnnotationLine.ReadFile(case_.normal_methylation_filename, case_.normal_methylation_file_id, false);

				foreach (var annotation in annotations)
				{
					var index = Math.Round(annotation.M_Value, 1); // 0.1, 0.2, ...
					int bin;
					lock (normalMethylationDistribution)
					{
						normalMethylationDistribution.TryGetValue(index, out bin);
						if (bin == 0)
						{
							normalMethylationDistribution.Add(index, 1);
						}
						else
						{
							normalMethylationDistribution[index] += 1;
						}
					}
				}
			}
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

					case_ = cases[0];

					cases.RemoveAt(0);
				}


				// verify methylation file exists for tumor data
				if (case_.tumor_methylation_filename.Contains("HumanMethylation450"))
				{
					ProcessFile(case_);
				}
			}

		}

		static void saveFile(string filename, bool forTumor)
		{
			// structure and save tumor file
			var min = forTumor ? tumorMethylationDistribution.Keys.Min() : normalMethylationDistribution.Keys.Min();
			var max = forTumor ? tumorMethylationDistribution.Keys.Max() : normalMethylationDistribution.Keys.Max();

			var output = ASETools.CreateStreamWriterWithRetry(filename);

			for (double i = min; i <= max; i += 0.1)
			{
				try
				{
					var count = tumorMethylationDistribution[i];
					output.WriteLine(i + "\t" + count);
				}
				catch (Exception)
				{
					output.WriteLine(i + "\t" + 0.0);
				}
			}

			output.Close();
		}



		static void Main(string[] args)
		{
			var configuration = ASETools.ASEConfirguation.loadFromFile(args);
			var cases = ASETools.Case.LoadCases(configuration.casesFilePathname);

			// Case 1: for any group of cases, compute the pvalues for adjusted beta value.
			// Run Mann Whitney for test for ASM values from adjusted beta value where the groups are 
			// ase vs no ase. (per gene) The purpose of this is to see if ASM significantly explains ASE anywhere



			// Case 2: Run t test, trying to assign distributions for methylation points. This is to only
			// look at 1 group (ie. ase) Thie purpose of this is to see what percentage of ASE significant
			// results can be explained by ASE. This should be written generically enough so it can
			// also be used to test no ASE for 0 mutations, testing what percentage of these cases fall into the
			// 'full methylation' category





			var threads = new List<Thread>();
			for (int i = 0; i < Environment.ProcessorCount; i++)
			{
				threads.Add(new Thread(() => ProcessCases(cases.Select(r => r.Value).ToList())));
			}


			threads.ForEach(th => th.Start());
			threads.ForEach(th => th.Join());

			saveFile(ASETools.ASEConfirguation.defaultBaseDirectory + "tumorMethylationDistribution.txt", true);
			saveFile(ASETools.ASEConfirguation.defaultBaseDirectory + "normalMethylationDistribution.txt", false);

		}
	}
}
