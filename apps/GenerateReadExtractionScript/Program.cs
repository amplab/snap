﻿using System;
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
            var configuration = ASETools.Configuration.loadFromFile(args);

            if (configuration == null)
            {
                Console.WriteLine("Giving up because we were unable to load configuration");
                return 1;
            }

            if (configuration.commandLineArgs.Count() < 6 || configuration.commandLineArgs.Count() % 1 != 0 || configuration.commandLineArgs[0] != "-d" && configuration.commandLineArgs[0] != "-r" || configuration.commandLineArgs[1] != "-n" && configuration.commandLineArgs[1] != "-t")
            {
                Console.WriteLine("usage: GenerateReadExtractionScript <-d|-r> <-n|-t> GenerateConsoldatedExtractedReadsPathname samtoolsPathname {(caseid outputfilename)}");
                return 1;
            }

            bool forDNA = configuration.commandLineArgs[0] == "-d";
            bool forTumor = configuration.commandLineArgs[1] == "-t";
            string generateConsoldatedExtractedReadsPathname = configuration.commandLineArgs[2];
            string samtoolsPathname = configuration.commandLineArgs[3];

            var cases = ASETools.Case.LoadCases(configuration.casesFilePathname);

            var casesToProcess = new List<CaseOutputFilenamePair>();

            for (int i = 4; i + 1< configuration.commandLineArgs.Count(); i += 2)
            {
                if (!cases.ContainsKey(args[i]))
                {
                    Console.WriteLine(configuration.commandLineArgs[i] + " does not appear to be a case id.  Ignoring.");
                } else
                {
                    casesToProcess.Add(new CaseOutputFilenamePair(cases[configuration.commandLineArgs[i]], configuration.commandLineArgs[i+1]));
                }
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

                string bamFileId;
                string inputBamFilename;

                if (forDNA)
                {
                    if (forTumor)
                    {
                        bamFileId = case_.tumor_dna_file_id;
                        inputBamFilename = case_.tumor_dna_filename;
                    }
                    else
                    {
                        bamFileId = case_.normal_dna_file_id;
                        inputBamFilename = case_.normal_dna_filename;
                    }
                }
                else
                {
                    if (forTumor)
                    {
                        bamFileId = case_.tumor_rna_file_id;
                        inputBamFilename = case_.tumor_rna_filename;
                    }
                    else
                    {
                        bamFileId = case_.normal_rna_file_id;
                        inputBamFilename = case_.normal_rna_filename;
                    }
                }

                var selectedVariantsFilename = case_.selected_variants_filename;
                var extractedMAFLinesFilename = case_.extracted_maf_lines_filename;

                if (selectedVariantsFilename == "" || extractedMAFLinesFilename == "" || inputBamFilename == "")
                {
                    Console.WriteLine("GenerateScriptsFromVariants: at least one input is missing for case " + case_.case_id);
                    File.Delete(localSamtoolsPathname);
                    File.Delete(localGCERPathname);
                    return 1;
                }

                int nVariants = 0;

                var selectedVariants = ASETools.SelectedVariant.LoadFromFile(selectedVariantsFilename);

                foreach (var selectedVariant in selectedVariants)
                {

                    generatedScript.Add("samtools view " + inputBamFilename + " " + selectedVariant.contig + ":" + Math.Max(1, selectedVariant.locus - 200) + "-" + (selectedVariant.locus + 10) +
                        @" > " + bamFileId + selectedVariant.getExtractedReadsExtension());

                    nVariants++;
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
                        @" > " + bamFileId + mafLine.getExtractedReadsExtension());
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
