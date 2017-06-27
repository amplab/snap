using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using ASELib;
using System.IO;
using System.Diagnostics;
using System.Threading;
using GenomeBuild;

namespace GenerateBisulfiteReadExtractionScript
{
	class Program
	{
		class CaseOutputFilenamePair
		{
			public CaseOutputFilenamePair(ASETools.Case case__, string outputFilename_)
			{
				case_ = case__;
				outputFilename = outputFilename_;
			}

			public readonly ASETools.Case case_;
			public readonly string outputFilename;
		}

		static int Main(string[] args)
		{

			var configuration = ASETools.ASEConfirguation.loadFromFile(args);

			if (configuration == null)
			{
				Console.WriteLine("Giving up because we were unable to load configuration");
				return 1;
			}

			if (configuration.commandLineArgs.Count() < 4 || configuration.commandLineArgs.Count() % 1 != 0 || configuration.commandLineArgs[0] != "-n" && configuration.commandLineArgs[0] != "-t")
			{
				Console.WriteLine("usage: GenerateReadExtractionScript <-n|-t> GenerateConsoldatedExtractedReadsPathname samtoolsPathname {(caseid outputfilename)}");
				return 1;
			}

			bool forTumor = configuration.commandLineArgs[0] == "-t";
			string generateConsoldatedExtractedReadsPathname = configuration.commandLineArgs[1];
			string samtoolsPathname = configuration.commandLineArgs[2];

			var cases = ASETools.Case.LoadCases(configuration.casesFilePathname);
			var bisulfiteCases = ASETools.BisulfateCase.loadCases(ASETools.ASEConfirguation.bisulfiteCasesFilePathname);

			var casesToProcess = new List<CaseOutputFilenamePair>();

			for (int i = 3; i< configuration.commandLineArgs.Count(); i ++)
			{
				if (!cases.ContainsKey(configuration.commandLineArgs[i]))
				{
					Console.WriteLine(configuration.commandLineArgs[i] + " does not appear to be a case id.  Ignoring.");
				}
				else
				{
					// filename path is based on bam id
					var filename = ASETools.ASEConfirguation.bisulfiteDirectory + configuration.commandLineArgs[i] + @"\" 
						+ bisulfiteCases[configuration.commandLineArgs[i]].bam_tumor_file_id;
					casesToProcess.Add(new CaseOutputFilenamePair(cases[configuration.commandLineArgs[i]], filename));
				}
			}

			//
			// Copy the generate consolodated extracted reads and samtools binaries to be local, so that we don't have a problem with overused remote
			// shares.
			//
			const string localSamtoolsPathname = @".\samtools.cmd";
			const string localGCERPathname = @".\GenerateConsolodatedExtractedReads.exe";
			bool copyWorked = false;
			for (int retryCount = 0; retryCount < 10; retryCount++)
			{
				try
				{
					File.Copy(generateConsoldatedExtractedReadsPathname, localGCERPathname, true);
					File.Copy(samtoolsPathname, localSamtoolsPathname, true);
					copyWorked = true;
					break;
				}
				catch
				{
					Console.WriteLine("Failed to copy binaries.  Sleeping 10s and retrying.");
					Thread.Sleep(10000);
				}
			}

			if (!copyWorked)
			{
				Console.WriteLine("Too many retries copying binaries.  Giving up.");
				File.Delete(localSamtoolsPathname);
				File.Delete(localGCERPathname);
				return 1;
			}


			foreach (var caseToProcess in casesToProcess)
			{
				var generatedScript = new List<string>();

				var case_ = caseToProcess.case_;
				var bisulfite_ = bisulfiteCases[case_.case_id];

				string bamFileId;
				string inputBamFilename;

				if (forTumor)
				{
					bamFileId = bisulfite_.bam_tumor_file_id;
					inputBamFilename = bisulfite_.bam_tumor_filename;
				}
				else
				{
					bamFileId = bisulfite_.bam_normal_file_id;
					inputBamFilename = bisulfite_.bam_normal_filename;
				}


				var selectedVariantsFilename = case_.selected_variants_filename;
				var extractedMAFLinesFilename = case_.extracted_maf_lines_filename;

				var mafLines = ASETools.MAFLine.ReadFile(extractedMAFLinesFilename, case_.maf_file_id, false);

				var build = new GenomeBuild.LiftOver();
				build.readChainFile(ASETools.ASEConfirguation.hg38Tohg19ChainFile); 

				// convert MAF file to hg19 build
				var MAFsHg19 = mafLines.Select(r =>
				   new Tuple<ASETools.MAFLine, List<GenomeBuild.Interval>>(r, build.mapCoordinates(new GenomeBuild.Interval(r.Chromosome, r.Start_Position, r.End_Positon))))
					.Where(r => r.Item2.Count() > 0).ToList(); // filter out mafLines that are not in the converted build

				mafLines = MAFsHg19.Select(r =>
				{
					r.Item1.Start_Position = Convert.ToInt32(r.Item2[0].start);
					r.Item1.End_Positon = Convert.ToInt32(r.Item2[0].end);
					return r.Item1;
				}).ToList();

				// save converted MAFLines
				var hg19FilenameMAF = ASETools.ASEConfirguation.bisulfiteDirectory + case_.case_id + @"\hg19" + case_.maf_filename.Split('\\').Last();
				hg19FilenameMAF = hg19FilenameMAF.Remove(hg19FilenameMAF.Length - ".gz".Length); // remove .gz extension
				ASETools.MAFLine.WriteToFile(hg19FilenameMAF, mafLines);

				// convert selected variants to hg19 build
				var selectedVariants = ASETools.SelectedVariant.LoadFromFile(selectedVariantsFilename);

				var convertedSelectedVariantsHg19 = selectedVariants.Select(r =>
				   new Tuple<ASETools.SelectedVariant, List<GenomeBuild.Interval>>(r, build.mapCoordinates(new GenomeBuild.Interval(r.contig, r.locus, r.locus))))
					.Where(r => r.Item2.Count() > 0).ToList(); // filter out mafLines that are not in the converted build

				selectedVariants = convertedSelectedVariantsHg19.Select(r =>
				{
					r.Item1.locus = Convert.ToInt32(r.Item2[0].start);
					return r.Item1;
				}).ToList();

				//  save file for new variants converted to hg19
				var hg19VariantsFilename = ASETools.ASEConfirguation.bisulfiteDirectory + case_.case_id + @"\hg19" + case_.selected_variants_filename.Split('\\').Last();
				var hg19VariantsFile = new StreamWriter(hg19VariantsFilename);
				hg19VariantsFile.WriteLine("SelectGermlineVariants v1.1 for input file " + case_.vcf_filename);      

				foreach (var selectedVariant in selectedVariants)
				{
					hg19VariantsFile.WriteLine(selectedVariant.contig + "\t" + selectedVariant.locus + "\t" + selectedVariant.contig + "\t" + selectedVariant.locus + "\t.\t" + selectedVariant.referenceBase + "\t" + selectedVariant.altBase);
				}

				hg19VariantsFile.WriteLine("**done**");
				hg19VariantsFile.Close();



				if (inputBamFilename == "")
				{
					Console.WriteLine("GenerateScriptsFromVariants: at least one input is missing for case " + case_.case_id);
					File.Delete(localSamtoolsPathname);
					File.Delete(localGCERPathname);
					return 1;
				}

				int nVariants = 0;

				foreach (var selectedVariant in selectedVariants)
				{

					generatedScript.Add("samtools view " + inputBamFilename + " " + selectedVariant.contig + ":" + Math.Max(1, selectedVariant.locus - 200) + "-" + (selectedVariant.locus + 10) +
						@" > " + inputBamFilename + selectedVariant.getExtractedReadsExtension());

					nVariants++;
				}

				if (null == mafLines)
				{
					return 1;   // ReadFile already printed an error message
				}

				foreach (var mafLine in mafLines)
				{
					generatedScript.Add("samtools view " + inputBamFilename + " " + mafLine.Chromosome + ":" + Math.Max(1, mafLine.Start_Position - 200) + "-" + (mafLine.End_Positon + 10) +
						@" > " + inputBamFilename + mafLine.getExtractedReadsExtension());
				}

				//
				// Create the output directory in case it doesn't exist.
				//
				Directory.CreateDirectory(ASETools.GetDirectoryFromPathname(caseToProcess.outputFilename));

				if (generatedScript.Count() == 0)
				{
					generatedScript.Add("This script intentionally left blank.");
				}

				//
				// And run GeneateConsolodatedExtractedReads on the file we created.
				//

				var startInfo = new ProcessStartInfo(localGCERPathname, " - " + caseToProcess.outputFilename + " 0 " + localSamtoolsPathname);
				startInfo.RedirectStandardInput = true;
				startInfo.UseShellExecute = false;
				var process = Process.Start(startInfo);

				foreach (var scriptLine in generatedScript)
				{
					process.StandardInput.WriteLine(scriptLine);
				}

				process.StandardInput.Close();

				process.WaitForExit();

				Console.WriteLine("Generate consolodated extracted reads for case id " + case_.case_id + " exited with code " + process.ExitCode);


				if (process.ExitCode != 0)
				{
					File.Delete(localSamtoolsPathname);
					File.Delete(localGCERPathname);

					File.Delete(caseToProcess.outputFilename);
					File.Delete(caseToProcess.outputFilename + ".index");
					return process.ExitCode;
				}
			} // Foreach case we're processing

			File.Delete(localSamtoolsPathname);
			File.Delete(localGCERPathname);

			return 0;
		}

	}
}
