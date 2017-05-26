using System;
using System.Diagnostics;
using ASELib;

namespace MethylationAnalysis
{
	class Program
	{
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

			foreach (var caseEntry in cases)
			{
				var case_ = caseEntry.Value;

				// verify methylation file exists
				if (case_.tumor_methylation_filename == "")
				{
					Console.WriteLine("No methylation data for case " + case_.case_id + ". Skipping...");
				} else
				{
					if (case_.normal_methylation_filename != "") {
						Console.WriteLine("Both tumor and normal methylation data for " + case_.case_id);
					}

					var annotations = ASELib.ASETools.AnnotationLine.ReadFile(case_.tumor_methylation_filename, case_.tumor_methylation_file_id, false);
					Console.WriteLine(annotations[0].Beta_Value);
				}
				
			}

			return 0;
		}
	}
}
