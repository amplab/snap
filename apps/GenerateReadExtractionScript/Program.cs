using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using ASELib;
using System.IO;
using System.Diagnostics;
using System.Threading;

namespace GenerateScriptFromVariants
{
    class Program
    {
        static int Main(string[] args)
        {
            var generatedScript = new List<string>();

            var configuration = ASETools.ASEConfirguation.loadFromFile(args);

            if (configuration == null)
            {
                Console.WriteLine("Giving up because we were unable to load configuration");
                return 1;
            }

            if (configuration.commandLineArgs.Count() != 6 || configuration.commandLineArgs[1] != "-d" && configuration.commandLineArgs[1] != "-r" || configuration.commandLineArgs[2] != "-n" && configuration.commandLineArgs[2] != "-t")
            {
                Console.WriteLine("usage: GenerateReadExtractionScript case_id <-d|-r> <-n|-t> GenerateConsoldatedExtractedReadsPathname outputFileName samtoolsPathname");
                return 1;
            }

            var case_id = configuration.commandLineArgs[0];
            bool forDNA = configuration.commandLineArgs[1] == "-d";
            bool forTumor = configuration.commandLineArgs[2] == "-t";
            string generateConsoldatedExtractedReadsPathname = configuration.commandLineArgs[2];
            string outputFilename = configuration.commandLineArgs[3];
            string samtoolsPathname = configuration.commandLineArgs[4];



            var cases = ASETools.Case.LoadCases(configuration.casesFilePathname);

            if (!cases.ContainsKey(case_id))
            {
                Console.WriteLine("Unable to find case id " + case_id);
                return 1;
            }

            var case_ = cases[case_id];

            string bamFileId;
            string inputBamFilename;

            if (forDNA)
            {
                if (forTumor)
                {
                    bamFileId = case_.tumor_dna_file_id;
                    inputBamFilename = case_.tumor_dna_filename;
                } else
                {
                    bamFileId = case_.normal_dna_file_id;
                    inputBamFilename = case_.normal_dna_filename;
                }
            } else
            {
                if (forTumor)
                {
                    bamFileId = case_.normal_dna_file_id;
                    inputBamFilename = case_.normal_dna_filename;
                } else
                {
                    bamFileId = case_.normal_rna_file_id;
                    inputBamFilename = case_.normal_rna_filename;
                }
            }

            var selectedVariantsFilename = case_.selected_variants_filename;
            var extractedMAFLinesFilename = case_.extracted_maf_lines_filename;

            if (selectedVariantsFilename == "" || extractedMAFLinesFilename == "" || inputBamFilename == "")
            {
                Console.WriteLine("GenerateScriptsFromVariants: at least one input is missing for case " + case_id);
                return 1;
            }

            int nVariants = 0;

            var selectedVariants = ASETools.SelectedVariant.LoadFromFile(selectedVariantsFilename);

            foreach (var selectedVariant in selectedVariants)
            { 

                generatedScript.Add("samtools view " + inputBamFilename + " " + chromosomeName + ":" + Math.Max(1, snvPosition - 200) + "-" + (snvPosition + 10) +
                    @" > " + bamFileId + "-" + chromosomeName + "-" + snvPosition);

                nVariants++;
            }

            if (!seenDone && !failed)
            {
                Console.WriteLine("Selected variants file " + selectedVariantsFilename + " is truncated.");
                return 1;
            }

            //
            // Now do reads around mutations in the MAF.
            //
            var mafLines = ASETools.MAFLine.ReadFile(case_.extracted_maf_lines_filename, case_.maf_file_id, false);

            if (null == mafLines)
            {
                return 1;   // ReadFile already printed an error message
            }

            foreach (var mafLine in mafLines)
            {
                generatedScript.Add("samtools view " + inputBamFilename + " " + mafLine.Chromosome + ":" + Math.Max(1, mafLine.Start_Position - 200) + "-" + (mafLine.End_Positon + 10) +
                    @" > " + bamFileId + "-" + mafLine.getExtractedReadsExtension());
            }

            //
            // Create the output directory in case it doesn't exist.
            //
            Directory.CreateDirectory(ASETools.GetDirectoryFromPathname(outputFilename));

            if (generatedScript.Count() == 0)
            {
                generatedScript.Add("This script intentionally left blank.");
            }

            //
            // Copy the generate consolodated extracted reads and samtools binaries to be local, so that we don't have a problem with overused remote
            // shares.
            //
            const string localSamtoolsPathname = @".\samtools.exe";
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
                } catch {
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

            //
            // And run GeneateConsolodatedExtractedReads on the file we created.
            //

            var startInfo = new ProcessStartInfo(localGCERPathname, " - " + outputFilename + " 0 " + localSamtoolsPathname);
            startInfo.RedirectStandardInput = true;
            startInfo.UseShellExecute = false;
            var process = Process.Start(startInfo);

            foreach (var scriptLine in generatedScript)
            {
                process.StandardInput.WriteLine(scriptLine);
            }

            process.StandardInput.Close();

            process.WaitForExit();

            Console.WriteLine("Generate consolodated extracted reads exited with code " + process.ExitCode);

            if (process.ExitCode != 0)
            {
                File.Delete(outputFilename);
                File.Delete(outputFilename + ".index");
            }

            File.Delete(localSamtoolsPathname);
            File.Delete(localGCERPathname);

            return process.ExitCode;
        }
    }
}
