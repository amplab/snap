using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using ASELib;

namespace WriteIGVScript
{
	class Program
	{
		static void PrintUsageMessage()
		{
			Console.WriteLine("Usage: WriteIGVScript.exe <case_ids> -g <Gene_name>");
			Console.WriteLine("Case_names is a list of case ids");
			Console.WriteLine("-p is optional position (gene_name or chrX:position)");
		}

		static void Main(string[] args)
		{
			var configuration = ASETools.Configuration.loadFromFile(args);
			var cases = ASETools.Case.LoadCases(configuration.casesFilePathname);

			if (null == cases)
			{
				Console.WriteLine("You must have selected cases before you can analyze methylation data.");
				return;
			}

			var selectedCases = new List<ASETools.Case>();
			if (configuration.commandLineArgs.Count() == 0)
			{
				PrintUsageMessage();
			}

			string position = "";

			for (var i = 0; i < configuration.commandLineArgs.Length; i++)
			{
				var arg = configuration.commandLineArgs[i];

				// getting optional position parameter. Should be last argument.
				if (arg == "-p")
				{
					position = configuration.commandLineArgs[i + 1];
					break;
				}

				ASETools.Case case_;
				if (cases.TryGetValue(arg, out case_))
				{
					selectedCases.Add(case_);
				}
				else
				{
					Console.WriteLine(arg + " is not a valid case id. Skipping ");
				}

			}

			var directory = ASETools.Configuration.defaultBaseDirectory + @"igv\";

			// write batch file
			Console.WriteLine("Writing batch files for " + selectedCases.Count() + " cases to " + directory);

			var ext = ".bat";
			foreach (var case_ in cases.Select(r => r.Value))
			{
				var filename = directory + case_.case_id + ext;
				ASETools.IGV.writeBatchFile(case_, filename, position);
			}

		}
	}
}
