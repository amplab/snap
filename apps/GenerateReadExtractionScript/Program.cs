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

            if (configuration.commandLineArgs.Count() != 5 || configuration.commandLineArgs[1] != "-d" && configuration.commandLineArgs[1] != "-r")
            {
                Console.WriteLine("usage: GenerateReadExtractionScript case_id <-d|-r> GenerateConsoldatedExtractedReadsPathname outputFileName samtoolsPathname");
                return 1;
            }

            var case_id = configuration.commandLineArgs[0];
            bool forDNA = configuration.commandLineArgs[1] == "-d";
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

            var selectedVariantsFilename = case_.selected_variants_filename;
            var extractedMAFLinesFilename = case_.extracted_maf_lines_filename;
            var inputBamFilename = forDNA ? case_.tumor_dna_filename : case_.tumor_rna_filename;

            if (selectedVariantsFilename == "" || extractedMAFLinesFilename == "" || inputBamFilename == "")
            {
                Console.WriteLine("GenerateScriptsFromVariants: at least one input is missing for case " + case_id);
                return 1;
            }

            int nVariants = 0;

            var selectedVariantsInputFile = ASETools.CreateStreamReaderWithRetry(selectedVariantsFilename);
 
            string line = selectedVariantsInputFile.ReadLine();

            if (null == line)
            {
                Console.WriteLine("Empty selected variants file " + selectedVariantsFilename);
                return 1;
            }

            string headerLinePrefix = "SelectGermlineVariants v";
            if (line.Count() < headerLinePrefix.Count() + 2 || line.Substring(0,headerLinePrefix.Count()) != headerLinePrefix)
            {
                Console.WriteLine("Selected variants file " + selectedVariantsFilename + " contains invalid header line");
                return 1;
            }

            if (line.Substring(headerLinePrefix.Count(), 2) != "1.")
            {
                Console.WriteLine("Selected variants file " + selectedVariantsFilename + " is of the wrong version.");
                return 1;
            }


            bool seenDone = false;
            bool failed = false;

            while (null != (line = selectedVariantsInputFile.ReadLine()))
            {
                if (seenDone)
                {
                    Console.WriteLine("Selected variants file " + selectedVariantsFilename + " extends beyond **done**");
                    failed = true;
                    break;
                }

                if ("**done**" == line)
                {
                    seenDone = true;
                    continue;
                }

                var fields = line.Split('\t');
                if (fields.Count() < 3) {
                    Console.WriteLine("Bad line in selected variants file " + selectedVariantsFilename + ": " + line);
                    failed = true;
                    break;

                }

                string chromosomeName = fields[0];

                int snvPosition = 0;
                try {
                    snvPosition = Convert.ToInt32(fields[1]);
                } catch (FormatException) {
                    Console.WriteLine("Unparsable pos in selected variants file " + selectedVariantsFilename + ": " + line);
                    failed = true;
                    break;
                }

                generatedScript.Add("samtools view " + inputBamFilename + " " + chromosomeName + ":" + Math.Max(1, snvPosition - 200) + "-" + (snvPosition + 10) +
                    @" > " + (forDNA ? case_.tumor_dna_file_id : case_.tumor_rna_file_id) + "-" + chromosomeName + "-" + snvPosition);

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
                    @" > " + (forDNA ? case_.tumor_dna_file_id : case_.tumor_rna_file_id) + "-" + mafLine.Chromosome + "-" + Math.Max(1, mafLine.Start_Position - 200) + "-" + (mafLine.End_Positon + 10));
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
