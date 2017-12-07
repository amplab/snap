using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using ASELib;
using System.IO;
using System.Diagnostics;

namespace ASEProcessManager
{
    //
    // The ASEProcessManager is responsible for driving the process of going from data stored in TCGA to a completed analysis.  There are two basic abstrations: Cases (patients with tumors
    // and apropriate data) and ProcessingStages, which are procedures that take a Case or set of Cases and run some analysis on it to produce more data.  Thus, the ProcessingStages
    // connect together to form the data flow graph for the system.  Each run of ASEProcessManager looks in the file system to find the state of the world and then generates scripts to
    // run processes that move things along toward a completed state.  These scripts will either download data from the Genome Data Commons or will read in existing data and do some processing
    // on it in order to produce an output file.
    //
    // So the overall way of running the complete experiment is to run ASEProcessManager, run the script that it produces, and repeat until ASEProcessManager says that all of the work is done.
    //
    class Program
    {
        const string scriptFilename = "ASENextSteps.cmd";
        const string linuxScriptFilename = "ASENextStepsLinux";
        const string downloadScriptFilename = "ASEDownload.cmd";

        static string jobAddString = "";

        //
        // A ProcessingStage examines the state of the world and if the prerequisites are there and the step isn't completed adds commands to the script to do some
        // work and advance the state of the world.
        //
        interface ProcessingStage
        {
            string GetStageName();
            bool NeedsCases();
            void EvaluateStage(StateOfTheWorld stateOfTheWorld, StreamWriter script, ASETools.RandomizingStreamWriter hpcScript, StreamWriter linuxScript, StreamWriter azureScript, out List<string> filesToDownload, out int nDone, out int nAddedToScript, out int nWaitingForPrerequisites);
            bool EvaluateDependencies(StateOfTheWorld stateOfTheWorld);
        }

        class MAFConfigurationProcessingStage : ProcessingStage 
        {
            public MAFConfigurationProcessingStage() { }

            public string GetStageName()
            {
                return "Generate MAF Configuration";
            }

            public bool NeedsCases() { return false; }

            public void EvaluateStage(StateOfTheWorld stateOfTheWorld, StreamWriter script, ASETools.RandomizingStreamWriter hpcScript, StreamWriter linuxScript, StreamWriter azureScript, out List<string> filesToDownload, out int nDone, out int nAddedToScript, out int nWaitingForPrerequisites)
            {
                filesToDownload = null;

                nWaitingForPrerequisites = 0; // This is the very first thing we do, there are never any prerequisites

                if (stateOfTheWorld.mafInfo != null)
                {
                    nDone = 1;
                    nAddedToScript = 0;
                    return;
                }

                script.WriteLine(stateOfTheWorld.configuration.binariesDirectory + "GenerateMAFConfiguration -configuration " + stateOfTheWorld.configuration.configuationFilePathname);

                nDone = 0;
                nAddedToScript = 1;
            }

            public bool EvaluateDependencies(StateOfTheWorld stateOfTheWorld)
            {
                //
                // This is the first stage to run, so there are no dependencies upstream of it.
                //
                return true;
            }
        }

        class GenerateCasesProcessingStage : ProcessingStage
        {
            public GenerateCasesProcessingStage() { }

            public string GetStageName()
            {
                return "Generate Cases";
            }

            public bool NeedsCases() { return false; }

            public void EvaluateStage(StateOfTheWorld stateOfTheWorld, StreamWriter script, ASETools.RandomizingStreamWriter hpcScript, StreamWriter linuxScript, StreamWriter azureScript, out List<string> filesToDownload, out int nDone, out int nAddedToScript, out int nWaitingForPrerequisites)
            {
                filesToDownload = null;

                if (stateOfTheWorld.cases != null)
                {
                    nDone = 1;
                    nAddedToScript = 0;
                    nWaitingForPrerequisites = 0;

                    return;
                }

                nDone = 0;

                if (stateOfTheWorld.mafInfo == null)
                {
                    nWaitingForPrerequisites = 1;
                    nAddedToScript = 0;

                    return;
                }

                //
                // See if we've downloaded all of the MAFs.
                //

                foreach (var mafEntry in stateOfTheWorld.mafInfo)
                {
                    if (!stateOfTheWorld.downloadedFiles.ContainsKey(mafEntry.Value.file_id))
                    {
                        if (null == filesToDownload)
                        {
                            filesToDownload = new List<string>();
                        }

                        filesToDownload.Add(mafEntry.Value.file_id);
                    }
                }

                if (null == filesToDownload)
                {
                    script.WriteLine(stateOfTheWorld.configuration.binariesDirectory + "GenerateCases -configuration " + stateOfTheWorld.configuration.configuationFilePathname);
                    nAddedToScript = 1;
                    nWaitingForPrerequisites = 0;
                }
                else
                {
                    nWaitingForPrerequisites = 1;
                    nAddedToScript = 0;
                }
            } // EvaluateStage

            public bool EvaluateDependencies(StateOfTheWorld stateOfTheWorld)
            {
                //
                // The cases file is updated every time we run, so it's kind of pointless to check to see if the mafs are newer than it.
                //
                return true;
            }

        } // GenerateCasesProcessingStage

        class AllcountProcesingStage : ProcessingStage
        {
            public AllcountProcesingStage() { }

            public string GetStageName()
            {
                return "Generate Allcount files";
            }

            public bool NeedsCases() { return true; }

            public void EvaluateStage(StateOfTheWorld stateOfTheWorld, StreamWriter script, ASETools.RandomizingStreamWriter hpcScript, StreamWriter linuxScript, StreamWriter azureScript, out List<string> filesToDownload, out int nDone, out int nAddedToScript, out int nWaitingForPrerequisites)
            {
                filesToDownload = null;
                nDone = 0;
                nAddedToScript = 0;
                nWaitingForPrerequisites = 0;
                filesToDownload = new List<string>();

                foreach (var caseEntry in stateOfTheWorld.cases)
                {
                    var case_ = caseEntry.Value;

                    HandleFile(stateOfTheWorld, case_.tumor_rna_file_id, case_.tumor_rna_file_bam_md5, case_.case_id, ASETools.DerivedFile.Type.TumorRNAAllcount,
                        ASETools.tumorRNAAllcountExtension, script, hpcScript, ref filesToDownload, ref nDone, ref nAddedToScript, ref nWaitingForPrerequisites);

                    HandleFile(stateOfTheWorld, case_.normal_dna_file_id, case_.normal_dna_file_bam_md5, case_.case_id, ASETools.DerivedFile.Type.NormalDNAAllcount,
                        ASETools.normalDNAAllcountExtension, script, hpcScript, ref filesToDownload, ref nDone, ref nAddedToScript, ref nWaitingForPrerequisites);

                    HandleFile(stateOfTheWorld, case_.tumor_dna_file_id, case_.tumor_dna_file_bam_md5, case_.case_id, ASETools.DerivedFile.Type.TumorDNAAllcount,
                        ASETools.tumorDNAAllcountExtension, script, hpcScript, ref filesToDownload, ref nDone, ref nAddedToScript, ref nWaitingForPrerequisites);

                    if (case_.normal_rna_file_id != "")
                    {
                        HandleFile(stateOfTheWorld, case_.normal_rna_file_id, case_.normal_rna_file_bam_md5, case_.case_id, ASETools.DerivedFile.Type.NormalRNAAllcount,
                            ASETools.normalRNAAllcountExtension, script, hpcScript, ref filesToDownload, ref nDone, ref nAddedToScript, ref nWaitingForPrerequisites);
                    }

                } // Foreach case
            }// EvaluateStage

            void HandleFile(StateOfTheWorld stateOfTheWorld, string file_id, string expectedMD5, string case_id, ASETools.DerivedFile.Type type, string extension, StreamWriter script, ASETools.RandomizingStreamWriter hpcScript, ref List<string> filesToDownload,ref int nDone, ref int nAddedToScript, ref int nWaitingForPrerequisites)
            {

                if (!stateOfTheWorld.downloadedFiles.ContainsKey(file_id))
                {
                    filesToDownload.Add(file_id);
                }
                else
                {
                    var downloadedFile = stateOfTheWorld.downloadedFiles[file_id];

                    if (!stateOfTheWorld.fileDownloadedAndVerified(file_id, expectedMD5))
                    {
                        nWaitingForPrerequisites++;
                    }
                    else if (stateOfTheWorld.containsDerivedFile(case_id, file_id, type))
                    {
                        if (stateOfTheWorld.getDrivedFile(case_id, file_id, type).fileinfo.Length < 200 * 1024)
                        {
                            Console.WriteLine("Suspiciously small allcount file of size " + stateOfTheWorld.getDrivedFile(case_id, file_id, type).fileinfo.Length + ": " + stateOfTheWorld.getDrivedFile(case_id, file_id, type).fileinfo.FullName);
                        }
                        nDone++;
                    }
                    else
                    {
                        nAddedToScript++;
                        string caseDirectory = ASETools.GetDirectoryFromPathname(stateOfTheWorld.downloadedFiles[file_id].fileInfo.FullName) + @"\..\..\" + stateOfTheWorld.configuration.derivedFilesDirectory + @"\" + case_id + @"\";
                        script.WriteLine("md " + caseDirectory);
                        script.WriteLine(stateOfTheWorld.configuration.binariesDirectory + "CountReadsCovering " + stateOfTheWorld.configuration.indexDirectory + " -a " + stateOfTheWorld.downloadedFiles[file_id].fileInfo.FullName + " - | gzip -9 > " +
                            caseDirectory + file_id + extension);

                        hpcScript.WriteLine(jobAddString + 
                            stateOfTheWorld.configuration.hpcBinariesDirectory + "MakeDirectoryAndCountReadsCovering.cmd " + caseDirectory + " " + stateOfTheWorld.configuration.hpcBinariesDirectory + " " +
                            stateOfTheWorld.configuration.hpcIndexDirectory + " " + stateOfTheWorld.downloadedFiles[file_id].fileInfo.FullName + " " + caseDirectory + file_id + extension);
                    }
                }
            } // HandleFile

            public bool EvaluateDependencies(StateOfTheWorld stateOFTheWorld) 
            {
                if (stateOFTheWorld.cases == null)
                {
                    return true;
                }

                bool allOK = true;
                foreach (var caseEntry in stateOFTheWorld.cases)
                {
                    var case_ = caseEntry.Value;
                    if (!stateOFTheWorld.containsDerivedFile(case_.case_id, case_.tumor_rna_file_id, ASETools.DerivedFile.Type.TumorRNAAllcount))
                    {
                        continue;
                    }

                    if (stateOFTheWorld.derivedFiles[case_.case_id].Where(x => x.type == ASETools.DerivedFile.Type.TumorRNAAllcount).Count() > 1)
                    {
                        Console.Write("More than one tumor RNA allcount file for case " + case_.case_id + ":");
                        foreach (var allcountFile in stateOFTheWorld.derivedFiles[case_.case_id].Where(x => x.type == ASETools.DerivedFile.Type.TumorRNAAllcount))
                        {
                            Console.Write(" " + allcountFile.fileinfo.FullName);
                        }
                        Console.WriteLine();
                        allOK = false;
                    }

                    var singleAllcountFile = stateOFTheWorld.derivedFiles[case_.case_id].Where(x => x.type == ASETools.DerivedFile.Type.TumorRNAAllcount).ToList()[0];

                    if (!stateOFTheWorld.downloadedFiles.ContainsKey(case_.tumor_rna_file_id))
                    {
                        Console.WriteLine("Allcount file " + singleAllcountFile.fileinfo.FullName + " exists, but the BAM from which it was generated does not");
                        allOK = false;
                    }
                }

                return allOK;
            }

        } // AllcountProcessingStage

        class DownloadProcessingStage : ProcessingStage
        {
            public DownloadProcessingStage() { }

            public string GetStageName()
            {
                return "Download";
            }
            public bool NeedsCases() { return true; }

            public void EvaluateStage(StateOfTheWorld stateOfTheWorld, StreamWriter script, ASETools.RandomizingStreamWriter hpcScript, StreamWriter linuxScript, StreamWriter azureScript, out List<string> filesToDownload, out int nDone, out int nAddedToScript, out int nWaitingForPrerequisites)
            {
                filesToDownload = new List<string>();
                nDone = 0;
                nAddedToScript = 0;
                nWaitingForPrerequisites = 0;

                foreach (var caseEntry in stateOfTheWorld.cases)
                {
                    var case_ = caseEntry.Value;

                    string[] idsToDownload = {case_.normal_dna_file_id, case_.tumor_dna_file_id, case_.normal_rna_file_id, case_.tumor_rna_file_id, case_.normal_methylation_file_id,
						case_.tumor_methylation_file_id, case_.normal_copy_number_file_id, case_.tumor_copy_number_file_id};

                    foreach (var id in idsToDownload) {
                        if (id != null && id != "" && !stateOfTheWorld.downloadedFiles.ContainsKey(id)) {
                            filesToDownload.Add(id);
                        }
                    }

                } // foreach case
            } // EvaluateStage

            public bool EvaluateDependencies(StateOfTheWorld stateOfTheWorld)
            {
                //
                // What we download depends on the cases selected, but we wouldn't re-download no matter what, since the data we get from the sevrer doesn't change.
                //
                return true;
            }
        } // DownloadProcessingStage

        class MD5ComputationProcessingStage : ProcessingStage
        {
            public MD5ComputationProcessingStage() { }

            public string GetStageName()
            {
                return "MD5 Computation";
            }

            public bool NeedsCases() { return true; }

            public void EvaluateStage(StateOfTheWorld stateOfTheWorld, StreamWriter script, ASETools.RandomizingStreamWriter hpcScript, StreamWriter linuxScript, StreamWriter azureScript, out List<string> filesToDownload, out int nDone, out int nAddedToScript, out int nWaitingForPrerequisites)
            {
                filesToDownload = null; // This stage never generates downloads
                nAddedToScript = 0;
                nDone = 0;

                nWaitingForPrerequisites = 0;

                foreach (var caseEntry in stateOfTheWorld.cases)
                {
                    var case_ = caseEntry.Value;

                    HandleFile(stateOfTheWorld, script, hpcScript, case_.tumor_rna_file_id, case_.tumor_rna_file_bam_md5, ref nDone, ref nAddedToScript, ref nWaitingForPrerequisites);
                    HandleFile(stateOfTheWorld, script, hpcScript, case_.normal_dna_file_id, case_.normal_dna_file_bam_md5, ref nDone, ref nAddedToScript, ref nWaitingForPrerequisites);
                    HandleFile(stateOfTheWorld, script, hpcScript, case_.tumor_dna_file_id, case_.tumor_dna_file_bam_md5, ref nDone, ref nAddedToScript, ref nWaitingForPrerequisites);
                    HandleFile(stateOfTheWorld, script, hpcScript, case_.normal_rna_file_id, case_.normal_rna_file_bam_md5, ref nDone, ref nAddedToScript, ref nWaitingForPrerequisites);
                    HandleFile(stateOfTheWorld, script, hpcScript, case_.tumor_methylation_file_id, case_.tumor_methylation_file_md5, ref nDone, ref nAddedToScript, ref nWaitingForPrerequisites);
        			HandleFile(stateOfTheWorld, script, hpcScript, case_.normal_methylation_file_id, case_.normal_methylation_file_md5, ref nDone, ref nAddedToScript, ref nWaitingForPrerequisites);
                    HandleFile(stateOfTheWorld, script, hpcScript, case_.tumor_copy_number_file_id, case_.tumor_copy_number_file_md5, ref nDone, ref nAddedToScript, ref nWaitingForPrerequisites);
					HandleFile(stateOfTheWorld, script, hpcScript, case_.normal_copy_number_file_id, case_.normal_copy_number_file_md5, ref nDone, ref nAddedToScript, ref nWaitingForPrerequisites);
				}
            } // EvaluateStage

            void HandleFile(StateOfTheWorld stateOfTheWorld, StreamWriter script, ASETools.RandomizingStreamWriter hpcScript, string fileId, string expectedMD5, ref int nDone, ref int nAddedToScript, ref int nWaitingForPrerequisites)
            {
                if (fileId == null || fileId == "")
                {
                    //
                    // There is no file at all (not just not downloaded), so it doesn't count in any of the counts.
                    //
                    return;
                }

                if (!stateOfTheWorld.downloadedFiles.ContainsKey(fileId) || null == expectedMD5 || "" == expectedMD5)
                {
                    nWaitingForPrerequisites++;
                    return;
                }

                var downloadedFile = stateOfTheWorld.downloadedFiles[fileId];

                if (downloadedFile.fileInfo.FullName.ToLower().EndsWith(".partial")) {
                    if (downloadedFile.fileInfo.LastWriteTime < DateTime.Now.AddDays(-1)) {
                        Console.WriteLine("Found partial download file that's more than a day old, it's probably abandoned and should be deleted: " + downloadedFile.fileInfo.FullName);
                    }
                    nWaitingForPrerequisites++;
                    return;
                }

                if (downloadedFile.storedMD5 != null && downloadedFile.storedMD5 != "")
                {
                    nDone++;

                    if (downloadedFile.storedMD5 != expectedMD5)
                    {
                        Console.WriteLine("MD5 checksum mismatch on file " + downloadedFile.fileInfo.FullName + " " + downloadedFile.storedMD5 + " != " + expectedMD5);
                    }

                    return;
                }

                script.WriteLine(stateOfTheWorld.configuration.binariesDirectory + "ComputeMD5 " + downloadedFile.fileInfo.FullName + " > " + downloadedFile.fileInfo.FullName + ".md5");
                hpcScript.WriteLine(jobAddString + stateOfTheWorld.configuration.hpcBinariesDirectory + "ComputeMD5IntoFile.cmd " +
                    stateOfTheWorld.configuration.hpcBinariesDirectory + " " + downloadedFile.fileInfo.FullName + " " + downloadedFile.fileInfo.FullName + ".md5");
                nAddedToScript++;
            }   // HandleFile

            public bool EvaluateDependencies(StateOfTheWorld stateOfTheWorld)
            {
                bool allOk = true;

                foreach (var fileEntry in stateOfTheWorld.downloadedFiles)
                {
                    var downloadedFile = fileEntry.Value;

                    if (downloadedFile.storedMD5 != null && downloadedFile.storedMD5 != "" && downloadedFile.md5FileInfo.LastWriteTime < downloadedFile.fileInfo.LastWriteTime)
                    {
                        Console.WriteLine("Downloaded file " + downloadedFile.fileInfo.FullName + " is newer than its md5 hash.");
                        allOk = false;
                    }
                }

                return allOk;
            } // EvaluateDependencies
        }

        class GermlineVariantCallingProcessingStage : ProcessingStage
        {
            public GermlineVariantCallingProcessingStage() { }

            public string GetStageName()
            {
                return "Germline Variant Calling";
            }

            public bool NeedsCases() { return true; }

            public void EvaluateStage(StateOfTheWorld stateOfTheWorld, StreamWriter script, ASETools.RandomizingStreamWriter hpcScript, StreamWriter linuxScript, StreamWriter azureScript, out List<string> filesToDownload, out int nDone, out int nAddedToScript, out int nWaitingForPrerequisites)
            {
                filesToDownload = new List<string>();
                nDone = 0;
                nAddedToScript = 0;
                nWaitingForPrerequisites = 0;

                foreach (var caseEntry in stateOfTheWorld.cases)
                {
                    var case_ = caseEntry.Value;

                    if (stateOfTheWorld.containsDerivedFile(case_.case_id, case_.normal_dna_file_id, ASETools.DerivedFile.Type.VCF))
                    {
                        nDone++;
                        continue;
                    }

                    //
                    // The Azure script downloads on the fly, so add every one that isn't done.
                    //
                    azureScript.Write("date\n");  // NB: we use Write and \n rather than WriteLine to avoid generating crlf text that would confuse Linux
                    azureScript.Write("rm -rf /mnt/downloaded_files/*\n"); // /mnt is a big but temporary filesystem on these azure instances
                    azureScript.Write("cd /mnt/downloaded_files\n");
                    azureScript.Write(@"~/gdc-client download --token-file ~/" + ASETools.GetFileNameFromPathname(stateOfTheWorld.configuration.accessTokenPathname) + " " + case_.normal_dna_file_id + "\n");
                    azureScript.Write("cd ~\n");
                    azureScript.Write("rm ~/x\n");    // We use a link from x to /mnt/downloaded_files/<download_directory> to make the command line shorter
                    azureScript.Write("ln -s /mnt/downloaded_files/" + case_.normal_dna_file_id + " ~/x\n");
                    azureScript.Write("cat ~/genomes/hg38-100k-regions | parallel -k -j `cat ~/ncores` \" freebayes --region {} --fasta-reference ~/genomes/hg38.fa ~/x/*.bam" +
                        " \" | ~/freebayes/vcflib/bin/vcffirstheader | ~/freebayes/vcflib/bin/vcfstreamsort -w 1000 | ~/freebayes/vcflib/bin/vcfuniq > ~/" +
                        case_.normal_dna_file_id + ASETools.vcfExtension + "\n");
                    azureScript.Write("if [ $? = 0 ]; then\n");
                    azureScript.Write("    mv " + case_.normal_dna_file_id + ASETools.vcfExtension + " ~/completed_vcfs/\n");
                    azureScript.Write("else\n");
                    azureScript.Write("    echo " + case_.normal_dna_file_id + " >> variant_calling_errors\n");
                    azureScript.Write("fi\n");
                    azureScript.Write("rm ~/" + case_.normal_dna_file_id + ASETools.vcfExtension + "\n");
                    azureScript.Write("rm -rf ~/downloaded_files/" + case_.normal_dna_file_id + "\n");
                    azureScript.Write("rm ~/x\n");

                    if (!stateOfTheWorld.fileDownloadedAndVerified(case_.normal_dna_file_id, case_.normal_dna_file_bam_md5))
                    {
                        nWaitingForPrerequisites++;
                        continue;
                    }

                    linuxScript.Write("date\n");    // NB: we use Write and \n rather than WriteLine to avoid generating crlf text that would confuse Linux
                    linuxScript.Write("cat ~/genomes/hg38-100k-regions | parallel -k -j `cat ~/ncores` \" freebayes --region {} --fasta-reference ~/genomes/hg38.fa " + 
                        ASETools.WindowsToLinuxPathname(stateOfTheWorld.downloadedFiles[case_.normal_dna_file_id].fileInfo.FullName) + 
                        " \" | ~/freebayes/vcflib/bin/vcffirstheader | ~/freebayes/vcflib/bin/vcfstreamsort -w 1000 | ~/freebayes/vcflib/bin/vcfuniq > " +
                        case_.normal_dna_file_id + ASETools.vcfExtension + "\n");
                    linuxScript.Write("if [ $? = 0 ]; then\n");
                    var outputDirectory = ASETools.WindowsToLinuxPathname(
                        ASETools.GetDirectoryPathFromFullyQualifiedFilename(stateOfTheWorld.downloadedFiles[case_.normal_dna_file_id].fileInfo.FullName) + @"..\..\" + stateOfTheWorld.configuration.derivedFilesDirectory + @"\" + case_.case_id + @"\"
                        );
                    linuxScript.Write(@"    mkdir " + outputDirectory + "\n");
                    linuxScript.Write(@"    cp " + case_.normal_dna_file_id + ASETools.vcfExtension + " " + outputDirectory + "\n");
                    linuxScript.Write("else\n");
                    linuxScript.Write(@"    echo " + case_.normal_dna_file_id + " >> variant_calling_errors\n");
                    linuxScript.Write("fi\n");
                    linuxScript.Write("rm " + case_.normal_dna_file_id + ASETools.vcfExtension + "\n");

                    nAddedToScript++;
                } // foreach case
            } // EvaluateStage

            public bool EvaluateDependencies(StateOfTheWorld stateOfTheWorld)
            {
                //
                // Skip this becuase we ran the variant calling on Azure, which downloaded the normal DNA files and then deleted them, 
                // so the local download of normal DNA may well be after the VCF was created.
                //

                return true;
            } // EvaluateDependencies

        }  // GermlineVariantCallingProcessingStage


		class AnnotateVariantsProcessingStage : ProcessingStage
		{

			public AnnotateVariantsProcessingStage() {}

			public string GetStageName() { return "Annotate Variants"; }

			public bool NeedsCases() { return true; }

			public void EvaluateStage(StateOfTheWorld stateOfTheWorld, StreamWriter script, ASETools.RandomizingStreamWriter hpcScript, StreamWriter linuxScript, StreamWriter azureScript, out List<string> filesToDownload, out int nDone, out int nAddedToScript, out int nWaitingForPrerequisites)
			{
				nDone = 0;
				nAddedToScript = 0;
				nWaitingForPrerequisites = 0;
				filesToDownload = null;

                string casesToProcess = "";

				foreach (var caseEntry in stateOfTheWorld.cases)
				{
					var case_ = caseEntry.Value;
					if (case_.annotated_selected_variants_filename != "")
					{
						nDone++;
						continue;
					}
					else if (case_.extracted_maf_lines_filename == "" || case_.selected_variants_filename == "" || case_.tumor_dna_reads_at_selected_variants_filename == "" || case_.tumor_dna_reads_at_selected_variants_index_filename == "" ||
					case_.tumor_rna_reads_at_selected_variants_filename == "" || case_.tumor_rna_reads_at_selected_variants_index_filename == "" || case_.normal_dna_reads_at_selected_variants_filename == "" || case_.normal_dna_reads_at_selected_variants_index_filename == "" ||
                    (case_.normal_rna_file_id != "" && (case_.normal_rna_reads_at_selected_variants_filename == "" || case_.normal_rna_reads_at_selected_variants_index_filename == "")))
					{
						nWaitingForPrerequisites++;
						continue;
					}
					else
					{
						nAddedToScript++;
					}

                    casesToProcess += case_.case_id + " ";  // Batch them both to allow for thread paralleleism within the program and also to reduce the number of (slow) job add commands to set up the script.

                    if (casesToProcess.Count() > 1000)
                    {
                        AddCasesToScripts(stateOfTheWorld, casesToProcess, script, hpcScript);
                        casesToProcess = "";
                    }
				} // foreach case

                if (casesToProcess != "")
                {
                    AddCasesToScripts(stateOfTheWorld, casesToProcess, script, hpcScript);
                }
			} // EvaluateStage

            void AddCasesToScripts(StateOfTheWorld stateOfTheWorld, string casesToProcess, StreamWriter script, ASETools.RandomizingStreamWriter hpcScript)
            {
                script.WriteLine(stateOfTheWorld.configuration.binariesDirectory + "AnnotateVariants.exe " + casesToProcess);
                hpcScript.WriteLine(jobAddString + stateOfTheWorld.configuration.hpcBinariesDirectory + "AnnotateVariants.exe " + casesToProcess);

            }

			public bool EvaluateDependencies(StateOfTheWorld stateOfTheWorld)
			{
				bool allOK = true;

				foreach (var caseEntry in stateOfTheWorld.cases)
				{
					var case_ = caseEntry.Value;

					if (case_.extracted_maf_lines_filename == "" || case_.selected_variants_filename == "" || case_.tumor_dna_reads_at_selected_variants_filename == "" || case_.tumor_dna_reads_at_selected_variants_index_filename == "" ||
					case_.tumor_rna_reads_at_selected_variants_filename == "" || case_.tumor_rna_reads_at_selected_variants_index_filename == "" || case_.normal_dna_reads_at_selected_variants_filename == "" || case_.normal_dna_reads_at_selected_variants_index_filename == "")
					{
						Console.WriteLine("Annotated variants file " + case_.annotated_selected_variants_filename + " exists, but dependencies do not.");
						allOK = false;
						continue;
					}

					if (case_.annotated_selected_variants_filename == "")
					{
						continue;
					}

					var annotatedVariantsWriteTime = new FileInfo(case_.annotated_selected_variants_filename).LastWriteTime;
					if (case_.annotated_selected_variants_filename == "")
					{
						Console.WriteLine("Annotated variants file " + case_.annotated_selected_variants_filename + " exists, but the precursor annotated selected variants file does not.");
						allOK = false;
						continue;
					}
				}
				return allOK;
			} // EvaluateDependencies

		} // AnnotateVariantsProcessingStage


		class MethylationProcessingStage : ProcessingStage
		{

			public MethylationProcessingStage() { }

			public string GetStageName() { return "Methylation"; }

			public bool NeedsCases() { return true; }

			public void EvaluateStage(StateOfTheWorld stateOfTheWorld, StreamWriter script, ASETools.RandomizingStreamWriter hpcScript, StreamWriter linuxScript, StreamWriter azureScript, out List<string> filesToDownload, out int nDone, out int nAddedToScript, out int nWaitingForPrerequisites)
			{
				nDone = 0;
				nAddedToScript = 0;
				nWaitingForPrerequisites = 0;
				filesToDownload = null;

				foreach (var caseEntry in stateOfTheWorld.cases)
				{
					var case_ = caseEntry.Value;
					if (case_.extracted_maf_lines_filename == "" || case_.tumor_methylation_filename == "")
					{
						nWaitingForPrerequisites++;
						continue;
					}
					else
					{
						nAddedToScript++;
					}

					script.WriteLine(stateOfTheWorld.configuration.binariesDirectory + "MethylationAnalysis.exe 1000000");
					hpcScript.WriteLine(jobAddString + stateOfTheWorld.configuration.hpcBinariesDirectory + "MethylationAnalysis.exe 1000000");
				}


			} // EvaluateStage

			public bool EvaluateDependencies(StateOfTheWorld stateOfTheWorld)
			{
				bool allOK = true;

				foreach (var caseEntry in stateOfTheWorld.cases)
				{
					var case_ = caseEntry.Value;

					if (case_.extracted_maf_lines_filename == "" || case_.tumor_methylation_filename == "")
					{
						Console.WriteLine("Regional methylation file " + case_.annotated_selected_variants_filename + " exists, but dependencies do not.");
						allOK = false;
						continue;
					}
				}
				return allOK;
			} // EvaluateDependencies

		} // MethylationProcessingStage


		class SelectVariantsProcessingStage : ProcessingStage
        {
            public SelectVariantsProcessingStage() { }

            public string GetStageName()
            {
                return "Select Germline Variants";
            }

            public bool NeedsCases() { return true; }

            public void EvaluateStage(StateOfTheWorld stateOfTheWorld, StreamWriter script, ASETools.RandomizingStreamWriter hpcScript, StreamWriter linuxScript, StreamWriter azureScript, out List<string> filesToDownload, out int nDone, out int nAddedToScript, out int nWaitingForPrerequisites)
            {
                nDone = 0;
                nAddedToScript = 0;
                filesToDownload = null;
                nWaitingForPrerequisites = 0;

                int nOnCurrentLine = 0;

                if (!File.Exists(stateOfTheWorld.configuration.redundantChromosomeRegionFilename))
                {
                    nWaitingForPrerequisites = 1;
                    return;
                }

                foreach (var caseEntry in stateOfTheWorld.cases)
                {
                    var case_ = caseEntry.Value;

                    if (case_.selected_variants_filename != "")
                    {
                        nDone++;
                        continue;
                    }

                    if (case_.vcf_filename == "" || case_.tumor_rna_allcount_filename == "" || case_.tumor_dna_allcount_filename == "" ||
                        case_.all_maf_lines_filename == "")
                    {
                        nWaitingForPrerequisites++;
                        continue;
                    }

                    nAddedToScript++;

                    if (nOnCurrentLine >= 100)
                    {
                        script.WriteLine();
                        hpcScript.WriteLine();
                        nOnCurrentLine = 0;
                    }

                    if (nOnCurrentLine == 0) 
                    {
                        script.Write(stateOfTheWorld.configuration.binariesDirectory + "SelectGermlineVariants.exe");
                        hpcScript.Write(jobAddString + stateOfTheWorld.configuration.hpcBinariesDirectory + "SelectGermlineVariants.exe");
                    }

                    script.Write(" " + case_.case_id);
                    hpcScript.Write(" " + case_.case_id);

                    nOnCurrentLine++;

                } // foreach case

                if (nOnCurrentLine > 0)
                {
                    script.WriteLine();
                    hpcScript.WriteLine();
                }
            } // EvaluateStage

            public bool EvaluateDependencies(StateOfTheWorld stateOfTheWorld)
            {
                bool allOK = true;

                foreach (var caseEntry in stateOfTheWorld.cases)
                {
                    var case_ = caseEntry.Value;

                    if (!stateOfTheWorld.containsDerivedFile(case_.case_id, case_.normal_dna_file_id, ASETools.DerivedFile.Type.SelectedVariants))
                    {
                        continue;
                    }

                    if (case_.vcf_filename == "" || case_.tumor_rna_allcount_filename == "" || case_.tumor_dna_allcount_filename == "")
                    {
                        Console.WriteLine(case_.selected_variants_filename + " depends on a file that is missing.");
                        allOK = false;
                        continue;
                    }

                    var selectedVariantsWriteTime = new FileInfo(case_.selected_variants_filename).LastWriteTime;
                    allOK &= checkOneDependency(case_.selected_variants_filename, selectedVariantsWriteTime, case_.vcf_filename);
                    allOK &= checkOneDependency(case_.selected_variants_filename, selectedVariantsWriteTime, case_.tumor_dna_allcount_filename);
                    allOK &= checkOneDependency(case_.selected_variants_filename, selectedVariantsWriteTime, case_.tumor_rna_allcount_filename);
                }

                return allOK;
            } // EvaluateDependencies

            bool checkOneDependency(string selectedVariantsFilename, DateTime selectedVariantsLastWriteTime, string sourceFilename)
            {
                if (selectedVariantsLastWriteTime < new FileInfo(sourceFilename).LastWriteTime)
                {
                    Console.WriteLine(selectedVariantsFilename + " is older than " + sourceFilename + ", upon which it depends.");
                    return false;
                }

                return true;
            } // checkOneDependency
        } // SelectVariantsProcessingStage


		class ExpressionDistributionProcessingStage : ProcessingStage
        {
            public ExpressionDistributionProcessingStage() { }

            public string GetStageName()
            {
                return "Per-disease mRNA expression distribution";
            }

            public bool NeedsCases() { return true; }

            public void EvaluateStage(StateOfTheWorld stateOfTheWorld, StreamWriter script, ASETools.RandomizingStreamWriter hpcScript, StreamWriter linuxScript, StreamWriter azureScript, out List<string> filesToDownload, out int nDone, out int nAddedToScript, out int nWaitingForPrerequisites)
            {
                nDone = 0;
                nAddedToScript = 0;
                nWaitingForPrerequisites = 0;
                filesToDownload = null;

                foreach (var disease in stateOfTheWorld.diseases)
                {
                    if (stateOfTheWorld.expressionFiles.ContainsKey(disease)) {
                        nDone++;
                    } else {
                        bool missingAny = false;
                        foreach (var caseEntry in stateOfTheWorld.cases.Where(x => x.Value.disease() == disease))
                        {
                            var case_ = caseEntry.Value;

                            if (!stateOfTheWorld.containsDerivedFile(case_.case_id, case_.tumor_rna_file_id, ASETools.DerivedFile.Type.TumorRNAAllcount) || 
                                !stateOfTheWorld.containsDerivedFile(case_.case_id, case_.tumor_rna_file_id, ASETools.DerivedFile.Type.TumorRNAMappedBaseCount))
                            {
                                nWaitingForPrerequisites++;
                                missingAny = true;
                                break;
                            }
                        }

                        if (missingAny)
                        {
                            continue;
                        }

                        string command = "ExpressionDistribution.exe " + stateOfTheWorld.configuration.casesFilePathname + " " +
                            stateOfTheWorld.configuration.expressionFilesDirectory + " " + ASETools.Case.ProjectColumn + " " + ASETools.Case.TumorRNAAllcountFilenameColumn + " " + ASETools.Case.TumorRNAMappedBaseCountColumn + " " + disease;

                        script.WriteLine(stateOfTheWorld.configuration.binariesDirectory + command);

                        hpcScript.WriteLine(jobAddString + stateOfTheWorld.configuration.hpcBinariesDirectory + command);
                        nAddedToScript++;
                    }
                } // foreach disease
            } // EvaluateStage

            public bool EvaluateDependencies(StateOfTheWorld stateOfTheWorld)
            {
                bool worked = true;

                foreach (var disease in stateOfTheWorld.diseases)
                {
                    if (!stateOfTheWorld.expressionFiles.ContainsKey(disease))
                    {
                        continue;
                    }

                    var expressionFileLastWriteTime = stateOfTheWorld.expressionFiles[disease].LastWriteTime;

                    foreach (var caseEntry in stateOfTheWorld.cases.Where(x => x.Value.disease() == disease))
                    {
                        var case_ = caseEntry.Value;

                        if (!stateOfTheWorld.containsDerivedFile(case_.case_id, case_.tumor_rna_file_id, ASETools.DerivedFile.Type.TumorRNAAllcount) ||
                            !stateOfTheWorld.containsDerivedFile(case_.case_id, case_.tumor_rna_file_id, ASETools.DerivedFile.Type.TumorRNAMappedBaseCount))
                        {
                            Console.WriteLine("Expression file " + stateOfTheWorld.expressionFiles[disease].FullName + " exists, but a prerequisite from case " + case_.case_id + " does not.");
                            worked = false;
                            break;   // Out of cases for this disease; keep checking other diseases
                        }

                        if (stateOfTheWorld.getDrivedFile(case_.case_id, case_.tumor_rna_file_id, ASETools.DerivedFile.Type.TumorRNAMappedBaseCount).fileinfo.LastWriteTime > expressionFileLastWriteTime)
                        {
                            Console.WriteLine("Expression file " + stateOfTheWorld.expressionFiles[disease].FullName + " is older than " +
                                stateOfTheWorld.getDrivedFile(case_.case_id, case_.tumor_rna_file_id, ASETools.DerivedFile.Type.TumorRNAMappedBaseCount).fileinfo.FullName + ", upon which it depends.");
                            worked = false;
                            break;  // Out of cases for this disease; keep checking other diseases
                        }

                        if (stateOfTheWorld.getDrivedFile(case_.case_id, case_.tumor_rna_file_id, ASETools.DerivedFile.Type.TumorRNAAllcount).fileinfo.LastWriteTime > expressionFileLastWriteTime)
                        {
                            Console.WriteLine("Expression file " + stateOfTheWorld.expressionFiles[disease].FullName + " is older than " +
                                stateOfTheWorld.getDrivedFile(case_.case_id, case_.tumor_rna_file_id, ASETools.DerivedFile.Type.TumorRNAAllcount).fileinfo.FullName + ", upon which it depends.");
                            worked = false;
                            break;  // Out of cases for this disease; keep checking other diseases
                        }
                    } // Foreach case
                } // foreach disease

                return worked;
            } // EvaluateDependencies
        } // ExpressionDistributionProcessingStage

        class ExtractMAFLinesProcessingStage : ProcessingStage
        {
            public string GetStageName() { return "Extract MAF Lines"; }

            public bool NeedsCases() { return true; }

            public void EvaluateStage(StateOfTheWorld stateOfTheWorld, StreamWriter script, ASETools.RandomizingStreamWriter hpcScript, StreamWriter linuxScript, StreamWriter azureScript, out List<string> filesToDownload, out int nDone, out int nAddedToScript, out int nWaitingForPrerequisites)
            {
                filesToDownload = null;

                nWaitingForPrerequisites = 0;
                nDone = 0;
                nAddedToScript = 0;

                foreach (var caseEntry in stateOfTheWorld.cases)
                {
                    var case_ = caseEntry.Value;

                    if (case_.maf_filename == null || case_.maf_filename == "")
                    {
                        nWaitingForPrerequisites++;
                    }
                    else if (case_.extracted_maf_lines_filename != "" && case_.all_maf_lines_filename != "")
                    {
                        nDone++;
                    }
                    else
                    {
                        nAddedToScript++;
                    }
                }

                if (nAddedToScript > 0)
                {
                    script.WriteLine(stateOfTheWorld.configuration.binariesDirectory + "ExtractMAFLines");
                    hpcScript.WriteLine(jobAddString + stateOfTheWorld.configuration.hpcBinariesDirectory + "ExtractMAFLines");
                }
            }

            public bool EvaluateDependencies(StateOfTheWorld stateOfTheWorld)
            {
                if (stateOfTheWorld.cases == null)
                {
                    return true;
                }

                bool allOK = true;
                foreach (var caseEntry in stateOfTheWorld.cases)
                {
                    var case_ = caseEntry.Value;

                    if (stateOfTheWorld.containsDerivedFile(case_.case_id, case_.case_id, ASETools.DerivedFile.Type.ExtractedMAFLines))
                    {
                        var derivedFile = stateOfTheWorld.derivedFiles[case_.case_id].Where(x => x.type == ASETools.DerivedFile.Type.ExtractedMAFLines).ToList()[0];

                        if (!stateOfTheWorld.downloadedFiles.ContainsKey(case_.maf_file_id)) {
                            Console.WriteLine("Case " + case_.case_id + " contains an extracted MAF lines file (" + derivedFile.fileinfo.FullName + "), but the corresponding MAF doesn't exist.");
                            allOK = false;
                            continue;
                        }

                        if (derivedFile.fileinfo.LastWriteTime < stateOfTheWorld.downloadedFiles[case_.maf_file_id].fileInfo.LastWriteTime) 
                        {
                            Console.WriteLine("Extracted MAF Lines file " + derivedFile.fileinfo.FullName + " is older than the MAF from which it's derived (" + stateOfTheWorld.downloadedFiles[case_.maf_file_id].fileInfo.FullName + ")");
                            allOK = false;
                        }
                    } // if the case has an extracted MAF Lines file
                } // foreach case

                return allOK;
            } // EvaluateDependencies

        } // ExtractMAFLinesProcessingStage

        class RegionalExpressionProcessingStage : ProcessingStage
        {
            public RegionalExpressionProcessingStage() { }

            public string GetStageName() { return "Regional Expression"; }

            public bool NeedsCases() { return true; }

            public void EvaluateStage(StateOfTheWorld stateOfTheWorld, StreamWriter script, ASETools.RandomizingStreamWriter hpcScript, StreamWriter linuxScript, StreamWriter azureScript, out List<string> filesToDownload, out int nDone, out int nAddedToScript, out int nWaitingForPrerequisites)
            {
                nDone = 0;
                nAddedToScript = 0;
                nWaitingForPrerequisites = 0;
                filesToDownload = null;

                var casesReadyToGoByDisease = new Dictionary<string, List<ASETools.Case>>();
                const int maxCasesPerCommandLine = 200;

                foreach (var caseEntry in stateOfTheWorld.cases)
                {
                    var case_ = caseEntry.Value;
                    if (case_.regional_expression_filename != "")
                    {
                        nDone++;
                    }
                    else if (!stateOfTheWorld.expressionFiles.ContainsKey(case_.disease()) || case_.tumor_rna_allcount_filename == "")
                    {
                        nWaitingForPrerequisites++;
                    }
                    else
                    {
                        nAddedToScript++;

                        if (!casesReadyToGoByDisease.ContainsKey(case_.disease()))
                        {
                            casesReadyToGoByDisease.Add(case_.disease(), new List<ASETools.Case>());
                        }

                        casesReadyToGoByDisease[case_.disease()].Add(case_);

                        if (casesReadyToGoByDisease[case_.disease()].Count() >= maxCasesPerCommandLine)
                        {
                            WriteScripts(stateOfTheWorld, casesReadyToGoByDisease[case_.disease()], script, hpcScript);
                            casesReadyToGoByDisease[case_.disease()] = new List<ASETools.Case>();
                        }
                    }
                } // foreach case

                foreach (var diseaseEntry in casesReadyToGoByDisease)
                {
                    if (diseaseEntry.Value.Count() > 0)
                    {
                        WriteScripts(stateOfTheWorld, diseaseEntry.Value, script, hpcScript);
                    }
                }
            } // EvaluateStage

            void WriteScripts(StateOfTheWorld stateOfTheWorld, List<ASETools.Case> cases, StreamWriter script, ASETools.RandomizingStreamWriter hpcScript)
            {
                script.Write(stateOfTheWorld.configuration.binariesDirectory + "RegionalExpression " + stateOfTheWorld.expressionFiles[cases[0].disease()].FullName + " " + stateOfTheWorld.configuration.regionalExpressionRegionSize + " ");
                hpcScript.Write(jobAddString + stateOfTheWorld.configuration.hpcBinariesDirectory + "RegionalExpression " + stateOfTheWorld.expressionFiles[cases[0].disease()].FullName + " " + stateOfTheWorld.configuration.regionalExpressionRegionSize + " ");
                foreach (var case_ in cases) {
                    script.Write(" " + case_.case_id);
                    hpcScript.Write(" " + case_.case_id);
                }
                script.WriteLine();
                hpcScript.WriteLine();
            } // WriteScripts

            public bool EvaluateDependencies(StateOfTheWorld stateOfTheWorld)
            {
                bool allOK = true;

                foreach (var caseEntry in stateOfTheWorld.cases)
                {
                    var case_ = caseEntry.Value;

                    if (case_.regional_expression_filename == "")
                    {
                        continue;
                    }

                    if (!stateOfTheWorld.expressionFiles.ContainsKey(case_.disease()))
                    {
                        Console.WriteLine("Missing per-disease expression file for " + case_.disease() + " even though regional expression file " + case_.regional_expression_filename + " exists.");
                        allOK = false;
                        continue;
                    }

                    var regionalExpressionWriteTime = new FileInfo(case_.regional_expression_filename).LastWriteTime;
                    if (stateOfTheWorld.expressionFiles[case_.disease()].LastWriteTime > regionalExpressionWriteTime)
                    {
                        Console.WriteLine("Regional expression file " + case_.regional_expression_filename + " is newer than the expression_ file on which it depends.");
                        allOK = false;
                        continue;
                    }

                    if (case_.tumor_rna_allcount_filename == "")
                    {
                        Console.WriteLine("Regional expression file " + case_.regional_expression_filename + " exists, but the precursor tumor rna allcount file does not.");
                        allOK = false;
                        continue;
                    }

                    if (new FileInfo(case_.tumor_rna_allcount_filename).LastWriteTime > regionalExpressionWriteTime)
                    {
                        Console.WriteLine("Regional expression file " + case_.regional_expression_filename + " is older than its tumor rna allcount file " + case_.tumor_rna_allcount_filename);
                        allOK = false;
                        continue;
                    }
                }

                return allOK;
            } // EvaluateDependencies

        } // RegionalExpressionProcessingStage

        class ExpressionNearMutationsProcessingStage : ProcessingStage
        {
			bool forAlleleSpecificExpression;

            public ExpressionNearMutationsProcessingStage(bool forAlleleSpecificExpression_) {
				forAlleleSpecificExpression = forAlleleSpecificExpression_;
			}

            public string GetStageName() { return "Expresssion Near Mutations"; }

            public bool NeedsCases() { return true; }

            public void EvaluateStage(StateOfTheWorld stateOfTheWorld, StreamWriter script, ASETools.RandomizingStreamWriter hpcScript, StreamWriter linuxScript, StreamWriter azureScript, out List<string> filesToDownload, out int nDone, out int nAddedToScript, out int nWaitingForPrerequisites)
            {
                nDone = 0;
                nAddedToScript = 0;
                nWaitingForPrerequisites = 0;
                filesToDownload = null;

                var currentCommandLine = "";

                foreach (var caseEntry in stateOfTheWorld.cases)
                {
                    var case_ = caseEntry.Value;
					if (forAlleleSpecificExpression)
					{
						if (case_.tumor_allele_specific_gene_expression_filename != "")
						{
							if (case_.normal_rna_filename == "")
							{
								nDone++;
								continue;
							}
							else // if normal rna exists and normal ASE was calculated
							{
								if (case_.normal_allele_specific_gene_expression_filename != "")
								{
									nDone++;
									continue;
								}
							}
						}
						else if (case_.maf_filename == "" || case_.annotated_selected_variants_filename == "" || case_.tumor_copy_number_filename == "")
						{
							nWaitingForPrerequisites++;
							continue;
						}
						else
						{
							nAddedToScript++;
						}
					}
					else
					{
						if (case_.gene_expression_filename != "")
						{
							nDone++;
							continue;
						}
						else if (case_.maf_filename == "" || case_.regional_expression_filename == ""  /* unfiltered counts, which we don't have a place for yet */)
						{
							nWaitingForPrerequisites++;
							continue;
						}
						else
						{
							nAddedToScript++;
						}
					}

                    if (currentCommandLine == "")
                    {
                        currentCommandLine = "ExpressionNearMutations.exe" + (forAlleleSpecificExpression ? " -a" : "");
                    }

                    currentCommandLine += " " + case_.case_id;

                    if (currentCommandLine.Count() > 5000)
                    {
                        script.WriteLine(stateOfTheWorld.configuration.binariesDirectory + currentCommandLine);
                        hpcScript.WriteLine(jobAddString + stateOfTheWorld.configuration.hpcBinariesDirectory + currentCommandLine);

                        currentCommandLine = "";
                    }

				}

                if (currentCommandLine != "")
                {
                    script.WriteLine(stateOfTheWorld.configuration.binariesDirectory + currentCommandLine);
                    hpcScript.WriteLine(jobAddString + stateOfTheWorld.configuration.hpcBinariesDirectory + currentCommandLine);
                }


            } // EvaluateStage

            public bool EvaluateDependencies(StateOfTheWorld stateOfTheWorld)
            {
				bool allOK = true;

				foreach (var caseEntry in stateOfTheWorld.cases)
				{
					var case_ = caseEntry.Value;

					if (case_.maf_filename == "")
					{
						Console.WriteLine("Gene expression file " + case_.gene_expression_filename + " exists, but the MAF file does not.");
						allOK = false;
						continue;
					}

					if (forAlleleSpecificExpression)
					{
						if (case_.tumor_allele_specific_gene_expression_filename == "" && case_.normal_allele_specific_gene_expression_filename == "")
						{
							continue;
						}

						if (case_.annotated_selected_variants_filename == "")
						{
							Console.WriteLine("AS gene expression file " + case_.tumor_allele_specific_gene_expression_filename + " exists, but the precursor annotated selected variants file does not.");
							allOK = false;
							continue;
						}

						var asGeneExpressionWriteTime = new FileInfo(case_.tumor_allele_specific_gene_expression_filename).LastWriteTime;
						if (new FileInfo(case_.annotated_selected_variants_filename).LastWriteTime > asGeneExpressionWriteTime)
						{
							Console.WriteLine("AS gene expression tumor file " + case_.tumor_allele_specific_gene_expression_filename + " is older than its regional file " + case_.regional_expression_filename);
							allOK = false;
							continue;
						}
						if (case_.normal_allele_specific_gene_expression_filename != "")
						{
							asGeneExpressionWriteTime = new FileInfo(case_.normal_allele_specific_gene_expression_filename).LastWriteTime;
							if (new FileInfo(case_.annotated_selected_variants_filename).LastWriteTime > asGeneExpressionWriteTime)
							{
								Console.WriteLine("AS gene expression normal file " + case_.normal_allele_specific_gene_expression_filename + " is older than its regional file " + case_.regional_expression_filename);
								allOK = false;
								continue;
							}
						}
					}
					else
					{
						if (case_.gene_expression_filename == "")
						{
							continue;
						}

						var geneExpressionWriteTime = new FileInfo(case_.gene_expression_filename).LastWriteTime;

						if (case_.regional_expression_filename == "")
						{
							Console.WriteLine("Gene expression file " + case_.gene_expression_filename + " exists, but the precursor regional expression file does not.");
							allOK = false;
							continue;
						}

						if (new FileInfo(case_.regional_expression_filename).LastWriteTime > geneExpressionWriteTime)
						{
							Console.WriteLine("Gene expression file " + case_.gene_expression_filename + " is older than its regional file " + case_.regional_expression_filename);
							allOK = false;
							continue;
						}
					}
				}
				return allOK;
			} // EvaluateDependencies

        } // ExpressionNearMutationsProcessingStage

        class ExtractReadsProcessingStage : ProcessingStage
        {
            public ExtractReadsProcessingStage() { }

            public string GetStageName() { return "Extract Reads"; }

            public bool NeedsCases() { return true; }

            void HandleFileAndType(StateOfTheWorld stateOfTheWorld, ASETools.Case case_, string flagsString, string readsAtSelectedVariantsFilename, string readsAtSelectedVariantsIndexFilename, string fileId, string md5Checksum, string inputFilename, string outputExtension,
                ref string outputSoFar, ref string outputSoFarHpc, StreamWriter script, ASETools.RandomizingStreamWriter hpcScript, ref int nDone, ref int nAddedToScript, ref int nWaitingForPrerequisites)
            {
                if (fileId == "")
                {
                    //
                    // Sometimes there is no normal RNA, which shows up as an empty fileID.
                    //
                    return;
                }

                if ((readsAtSelectedVariantsFilename == "") != (readsAtSelectedVariantsIndexFilename == ""))
                {
                    Console.WriteLine("Exactly one of reads at selected variants and reads at selected variants index files exists: " + readsAtSelectedVariantsFilename + " " + readsAtSelectedVariantsIndexFilename);
                    return;
                }

                if (readsAtSelectedVariantsFilename != "")
                {
                    if (false)  // This is off because it's slow
                    {
                        var writeTime1 = new FileInfo(readsAtSelectedVariantsFilename).LastWriteTime;
                        var writeTime2 = new FileInfo(readsAtSelectedVariantsIndexFilename).LastWriteTime;

                        if (writeTime1 > writeTime2.AddDays(1) || writeTime2 > writeTime1.AddDays(1))
                        {
                            Console.WriteLine("Extracted reads and index files differ by more than one day: " + readsAtSelectedVariantsFilename + " " + readsAtSelectedVariantsIndexFilename);
                        }
                    }

                    nDone++;
                    return;
                }

                if (!stateOfTheWorld.fileDownloadedAndVerified(fileId, md5Checksum) || case_.selected_variants_filename == "" || case_.extracted_maf_lines_filename == "")
                {
                    nWaitingForPrerequisites++;
                    return;
                }

                nAddedToScript++;

                if (outputSoFar == "")
                {
                    outputSoFar = stateOfTheWorld.configuration.binariesDirectory + "GenerateReadExtractionScript " + flagsString + " " + stateOfTheWorld.configuration.binariesDirectory + "GenerateConsolodatedExtractedReads.exe " +
                                stateOfTheWorld.configuration.binariesDirectory + "samtools.exe ";

                    outputSoFarHpc = jobAddString + stateOfTheWorld.configuration.hpcBinariesDirectory + "GenerateReadExtractionScript " + flagsString + " " + stateOfTheWorld.configuration.hpcBinariesDirectory + "GenerateConsolodatedExtractedReads.exe " +
                                stateOfTheWorld.configuration.hpcBinariesDirectory + "samtools.exe ";
                }

                string outputFilename = ASETools.GoUpFilesystemLevels(ASETools.GetDirectoryFromPathname(inputFilename), 2) + stateOfTheWorld.configuration.derivedFilesDirectory + @"\" + case_.case_id + @"\" + fileId + outputExtension + " ";

                outputSoFar += case_.case_id + " " + outputFilename;
                outputSoFarHpc += case_.case_id + " " + outputFilename;

                if (Math.Max(outputSoFar.Count(), outputSoFarHpc.Count()) > 2000)
                {
                    script.WriteLine(outputSoFar);
                    hpcScript.WriteLine(outputSoFarHpc);

                    outputSoFar = "";
                    outputSoFarHpc = "";
                }
            }

            public void EvaluateStage(StateOfTheWorld stateOfTheWorld, StreamWriter script, ASETools.RandomizingStreamWriter hpcScript, StreamWriter linuxScript, StreamWriter azureScript, out List<string> filesToDownload, out int nDone, out int nAddedToScript, out int nWaitingForPrerequisites)
            {
                filesToDownload = null;
                nDone = 0;
                nAddedToScript = 0;
                nWaitingForPrerequisites = 0;

                string tumorDNAOutput = "";
                string tumorDNAHpcOutput = "";
                string normalDNAOutput = "";
                string normalDNAHpcOutput = "";
                string tumorRNAOutput = "";
                string tumorRNAHpcOutput = "";
                string normalRNAOutput = "";
                string normalRNAHpcOutput = "";

                foreach (var caseEntry in stateOfTheWorld.cases)
                {
                    var case_ = caseEntry.Value;

                    HandleFileAndType(stateOfTheWorld, case_, "-d -t", case_.tumor_dna_reads_at_selected_variants_filename, case_.tumor_dna_reads_at_selected_variants_index_filename, case_.tumor_dna_file_id,  case_.tumor_dna_file_bam_md5,  case_.tumor_dna_filename,  ASETools.tumorDNAReadsAtSelectedVariantsExtension,  ref tumorDNAOutput,  ref tumorDNAHpcOutput,  script, hpcScript, ref nDone, ref nAddedToScript, ref nWaitingForPrerequisites);
                    HandleFileAndType(stateOfTheWorld, case_, "-d -n", case_.normal_dna_reads_at_selected_variants_filename, case_.normal_dna_reads_at_selected_variants_index_filename, case_.normal_dna_file_id, case_.normal_dna_file_bam_md5, case_.normal_dna_filename, ASETools.normalDNAReadsAtSelectedVariantsExtension, ref normalDNAOutput, ref normalDNAHpcOutput, script, hpcScript, ref nDone, ref nAddedToScript, ref nWaitingForPrerequisites);
                    HandleFileAndType(stateOfTheWorld, case_, "-r -t", case_.tumor_rna_reads_at_selected_variants_filename, case_.tumor_rna_reads_at_selected_variants_index_filename, case_.tumor_rna_file_id,  case_.tumor_rna_file_bam_md5,  case_.tumor_rna_filename,  ASETools.tumorRNAReadsAtSelectedVariantsExtension,  ref tumorRNAOutput,  ref tumorRNAHpcOutput,  script, hpcScript, ref nDone, ref nAddedToScript, ref nWaitingForPrerequisites);
                    HandleFileAndType(stateOfTheWorld, case_, "-r -n", case_.normal_rna_reads_at_selected_variants_filename, case_.normal_rna_reads_at_selected_variants_index_filename, case_.normal_rna_file_id, case_.normal_rna_file_bam_md5, case_.normal_rna_filename, ASETools.normalRNAReadsAtSelectedVariantsExtension, ref normalRNAOutput, ref normalRNAHpcOutput, script, hpcScript, ref nDone, ref nAddedToScript, ref nWaitingForPrerequisites);
                }

                string[] outputs = { tumorDNAOutput, normalDNAOutput, tumorRNAOutput, normalRNAOutput };
                foreach (var output in outputs)
                {
                    if (output != "")
                    {
                        script.WriteLine(output);
                    }
                }

                string[] hpcOutputs = { tumorDNAHpcOutput, normalDNAHpcOutput, tumorRNAHpcOutput, normalRNAHpcOutput };
                foreach (var hpcOutput in hpcOutputs)
                {
                    if (hpcOutput != "")
                    {
                        hpcScript.WriteLine(hpcOutput);
                    }
                }
            } // EvaluateStage


            public bool EvaluateDependencies(StateOfTheWorld stateOfTheWorld)
            {
                // Fill this in later
                return true;
            } // EvaluateDependencies
        } // ExtractReadsProcessingStage

        class SelectGenesProcessingStage : ProcessingStage
        {
            public SelectGenesProcessingStage() { }

            public string GetStageName()
            {
                return "Select Genes";
            }

            public bool NeedsCases() { return true; }

            public void EvaluateStage(StateOfTheWorld stateOfTheWorld, StreamWriter script, ASETools.RandomizingStreamWriter hpcScript, StreamWriter linuxScript, StreamWriter azureScript, out List<string> filesToDownload, out int nDone, out int nAddedToScript, out int nWaitingForPrerequisites)
            {
                nDone = 0;
                nAddedToScript = 0;
                nWaitingForPrerequisites = 0;
                filesToDownload = new List<string>();

                if (stateOfTheWorld.selectedGenes != null)
                {
                    nDone++;
                    return;
                }

                foreach (var caseEntry in stateOfTheWorld.cases)
                {
                    if (caseEntry.Value.extracted_maf_lines_filename == "")
                    {
                        nWaitingForPrerequisites++;
                        return;
                    }
                }

                nAddedToScript++;
                script.WriteLine(stateOfTheWorld.configuration.binariesDirectory + "SelectGenes.exe");
                hpcScript.WriteLine(jobAddString + stateOfTheWorld.configuration.hpcBinariesDirectory + "SelectGenes.exe");
            }

            public bool EvaluateDependencies(StateOfTheWorld stateOfTheWorld)
            {
                if (stateOfTheWorld.selectedGenes == null)
                {
                    return true;    // No output means no violated dependencies
                }

                var selectedGenesWriteTime = new FileInfo(stateOfTheWorld.configuration.selectedGenesFilename).LastWriteTime;

                foreach (var caseEntry in stateOfTheWorld.cases)
                {
                    var case_ = caseEntry.Value;

                    if (case_.extracted_maf_lines_filename == "" )
                    {
                        Console.WriteLine("Dependency violation: case " + case_.case_id + " does not have extracted MAF lines, but we have selected genes.");
                        return false;   // Don't worry about the others, we need to regenerate the selected genes file regardless.
                    }

                    if (new FileInfo (case_.extracted_maf_lines_filename).LastWriteTime > selectedGenesWriteTime)
                    {
                        Console.WriteLine("Dependency violation: the selected genes file is older than an extracted MAF lines file " + case_.extracted_maf_lines_filename);
                        return false;
                    }
                }

                return true;
            } // EvaluateDependencies
        } // SelectGenesProcessingStage

        class CountMappedBasesProcessingStage : ProcessingStage
        {
            public CountMappedBasesProcessingStage() { }

            public string GetStageName() { return "Count Mapped Bases"; }

            public bool NeedsCases() { return true; }

            public void EvaluateStage(StateOfTheWorld stateOfTheWorld, StreamWriter script, ASETools.RandomizingStreamWriter hpcScript, StreamWriter linuxScript, StreamWriter azureScript, out List<string> filesToDownload, out int nDone, out int nAddedToScript, out int nWaitingForPrerequisites)
            {
                nDone = 0;
                nAddedToScript = 0;
                nWaitingForPrerequisites = 0;
                filesToDownload = null;

                foreach (var caseEntry in stateOfTheWorld.cases)
                {
                    var case_ = caseEntry.Value;
                    HandleFile(stateOfTheWorld, script, hpcScript, case_.normal_dna_allcount_filename, case_.normal_dna_mapped_base_count_filename, ASETools.normalDNAMappedBaseCountExtension, ref nDone, ref nAddedToScript, ref nWaitingForPrerequisites);
                    HandleFile(stateOfTheWorld, script, hpcScript, case_.tumor_dna_allcount_filename, case_.tumor_dna_mapped_base_count_filename, ASETools.tumorDNAMappedBaseCountExtension, ref nDone, ref nAddedToScript, ref nWaitingForPrerequisites);
                    if (case_.normal_rna_file_id != "" && case_.normal_rna_file_id != null)
                    {
                        HandleFile(stateOfTheWorld, script, hpcScript, case_.normal_rna_allcount_filename, case_.normal_rna_mapped_base_count_filename, ASETools.normalRNAMappedBaseCountExtension, ref nDone, ref nAddedToScript, ref nWaitingForPrerequisites);
                    }
                    HandleFile(stateOfTheWorld, script, hpcScript, case_.tumor_rna_allcount_filename, case_.tumor_rna_mapped_base_count_filename, ASETools.tumorRNAMappedBaseCountExtension, ref nDone, ref nAddedToScript, ref nWaitingForPrerequisites);
                } // foreach case
            } // EvaluateStage

            void HandleFile(StateOfTheWorld stateOfTheWorld, StreamWriter script, ASETools.RandomizingStreamWriter hpcScript, string inputFilename, string existingOutputFilename, string outputExtension, ref int nDone, ref int nAddedToScript, ref int nWaitingForPrerequisites)
            {
                if (existingOutputFilename != "")
                {
                    nDone++;
                    return;
                }

                if (inputFilename == "")
                {
                    nWaitingForPrerequisites++;
                    return;
                }

                script.WriteLine(stateOfTheWorld.configuration.binariesDirectory + "CountMappedBases.exe " + inputFilename + " " + ASETools.GetDirectoryFromPathname(inputFilename) + @"\" + ASETools.GetFileIdFromPathname(inputFilename) + outputExtension);
                hpcScript.WriteLine(jobAddString + stateOfTheWorld.configuration.hpcBinariesDirectory + "CountMappedBases.exe " + inputFilename + " " + ASETools.GetDirectoryFromPathname(inputFilename) + @"\" + ASETools.GetFileIdFromPathname(inputFilename) + outputExtension);

                nAddedToScript++;
            } // HandleFile

            public bool EvaluateDependencies(StateOfTheWorld stateOfTheWorld)
            {
                bool worked = true;
                foreach (var caseEntry in stateOfTheWorld.cases)
                {
                    var case_ = caseEntry.Value;

                    worked &= HandleDependency(case_.normal_dna_mapped_base_count_filename, case_.normal_dna_allcount_filename);
                    worked &= HandleDependency(case_.tumor_dna_mapped_base_count_filename, case_.tumor_dna_allcount_filename);
                    worked &= HandleDependency(case_.normal_rna_mapped_base_count_filename, case_.normal_rna_allcount_filename);
                    worked &= HandleDependency(case_.tumor_rna_mapped_base_count_filename, case_.tumor_rna_allcount_filename);
                }

                return worked;
            } // EvaluateDependencies

            bool HandleDependency(string baseCountFilename, string allcountFilename)
            {
                if (baseCountFilename == "")
                {
                    return true;
                }

                if (allcountFilename == "")
                {
                    Console.WriteLine("Base count file " + baseCountFilename + " exists, but its prerequisite allcount file does not.");
                    return false;
                }

                if (new FileInfo(baseCountFilename).LastWriteTime < new FileInfo(allcountFilename).LastWriteTime)
                {
                    Console.WriteLine("Base count file " + baseCountFilename + " is older than its predecessor " + allcountFilename);
                    return false;
                }

                return true;
            } // HandleDependency

        } // CountMappedBasesProcessingStage

        class GenerateScatterGraphsProcessingStage : ProcessingStage
        {
            public GenerateScatterGraphsProcessingStage() { }

            public string GetStageName() { return "Generate Scatter Graphs"; }

            public bool NeedsCases() { return true; }

            public void EvaluateStage(StateOfTheWorld stateOfTheWorld, StreamWriter script, ASETools.RandomizingStreamWriter hpcScript, StreamWriter linuxScript, StreamWriter azureScript, out List<string> filesToDownload, out int nDone, out int nAddedToScript, out int nWaitingForPrerequisites)
            {
                filesToDownload = null;
                nDone = 0;
                nAddedToScript = 0;
                nWaitingForPrerequisites = 0;

                bool missingAnything = false;

                if (stateOfTheWorld.scatterGraphsSummaryFile == "")
                {
                    missingAnything = true;
                } 
                //
                // Don't check that there are files for the selected genes, since they might not be there due to too small n.
                //


                if (!missingAnything)
                {
                    nDone = 1;
                    return;
                }

                foreach (var caseEntry in stateOfTheWorld.cases)
                {
                    var case_ = caseEntry.Value;

                    if (case_.annotated_selected_variants_filename == "" || case_.tumor_dna_mapped_base_count_filename == "" || case_.tumor_rna_mapped_base_count_filename == "")
                    {
                        nWaitingForPrerequisites = 1;
                        return;
                    }
                }

                nAddedToScript = 1;

                script.WriteLine(stateOfTheWorld.configuration.binariesDirectory + "GenerateScatterGraphs.exe");
                hpcScript.WriteLine(jobAddString + stateOfTheWorld.configuration.hpcBinariesDirectory + "GenerateScatterGraphs.exe");
            } // EvaluateStage

            public bool EvaluateDependencies(StateOfTheWorld stateOfTheWorld)
            {
                DateTime oldestFile;
                if (stateOfTheWorld.scatterGraphsSummaryFile == "")
                {
                    return true;
                }
 
                oldestFile = new FileInfo(stateOfTheWorld.scatterGraphsSummaryFile).LastWriteTime;
                foreach (var selectedGene in stateOfTheWorld.selectedGenes)
                {
                    if (!stateOfTheWorld.scatterGraphsByHugoSymbol.ContainsKey(selectedGene.Hugo_Symbol))
                    {
                        return true;
                    }

                    var date = new FileInfo(stateOfTheWorld.scatterGraphsByHugoSymbol[selectedGene.Hugo_Symbol]).LastWriteTime;
                    if (date < oldestFile) {
                        oldestFile = date;
                    }
                }

                //
                // Now look at the dependencies.
                //
                foreach (var caseEntry in stateOfTheWorld.cases)
                {
                    var case_ = caseEntry.Value;

                    if (case_.annotated_selected_variants_filename == "" || new FileInfo(case_.annotated_selected_variants_filename).LastWriteTime > oldestFile)
                    {
                        Console.WriteLine("Case " + case_.case_id + " either doesn't have an annotated selected variants file, or it is newer than a gene scatter graph file that depends on it.");
                        return false;
                    }

                    if (case_.tumor_dna_mapped_base_count_filename == "" || new FileInfo(case_.tumor_dna_mapped_base_count_filename).LastWriteTime > oldestFile ||
                        case_.tumor_rna_mapped_base_count_filename == "" || new FileInfo(case_.tumor_rna_mapped_base_count_filename).LastWriteTime > oldestFile)
                    {
                        Console.WriteLine("Case " + case_.case_id + " either doesn't have a tumor mapped base count file, or it is newer than a gene scatter graph file that depends on it.");
                        return false;
                    }
                }

                return true;
            } // EvaluateDependencies

        } // GenerateScatterGraphsProcessingStage

        class GenerateAllLociProcessingStage : ProcessingStage
        {
            public GenerateAllLociProcessingStage() { }

            public string GetStageName() { return "Generate All Loci Files"; }

            public bool NeedsCases() { return false; }

            public void EvaluateStage(StateOfTheWorld stateOfTheWorld, StreamWriter script, ASETools.RandomizingStreamWriter hpcScript, StreamWriter linuxScript, StreamWriter azureScript, out List<string> filesToDownload, out int nDone, out int nAddedToScript, out int nWaitingForPrerequisites)
            {
                nDone = 0;
                nAddedToScript = 0;
                nWaitingForPrerequisites = 0;   // There are no prerequesites for this, so it will always stay zero
                filesToDownload = null;

                for (int chromosomeNumber = 1; chromosomeNumber <= ASETools.nHumanAutosomes; chromosomeNumber++)
                {
                    if (File.Exists(stateOfTheWorld.configuration.chromosomeMapsDirectory + "chr" + chromosomeNumber + ASETools.allLociExtension))
                    {
                        nDone++;
                    } else
                    {
                        script.WriteLine(stateOfTheWorld.configuration.binariesDirectory + "GenerateReadsForRepetitiveRegionDetection.exe " + stateOfTheWorld.configuration.indexDirectory + " " + chromosomeNumber + " " +
                            stateOfTheWorld.configuration.chromosomeMapsDirectory +"chr" + chromosomeNumber +
                            ASETools.allLociExtension + " 48");
                        nAddedToScript++;
                    }
                }
            }

            public bool EvaluateDependencies(StateOfTheWorld stateOfTheWorld) { return true; }   // Really no dependencies for this unless the reference genome changes
        } // GenerateAllLociProcessingStage

        class AlignAllLociProcessingStage : ProcessingStage
        {
            public AlignAllLociProcessingStage() { }

            public string GetStageName() { return "Align all loci"; }

            public bool NeedsCases() { return false; }

            public void EvaluateStage(StateOfTheWorld stateOfTheWorld, StreamWriter script, ASETools.RandomizingStreamWriter hpcScript, StreamWriter linuxScript, StreamWriter azureScript, out List<string> filesToDownload, out int nDone, out int nAddedToScript, out int nWaitingForPrerequisites)
            {
                nDone = 0;
                nAddedToScript = 0;
                nWaitingForPrerequisites = 0;
                filesToDownload = null;

                for (int chromosomeNumber = 1; chromosomeNumber <= ASETools.nHumanAutosomes; chromosomeNumber++)
                {
                    if (File.Exists(stateOfTheWorld.configuration.chromosomeMapsDirectory + "chr" + chromosomeNumber + ASETools.allLociAlignedExtension))
                    {
                        nDone++;
                    }
                    else if (!File.Exists(stateOfTheWorld.configuration.chromosomeMapsDirectory + "chr" + chromosomeNumber + ASETools.allLociExtension))
                    {
                        nWaitingForPrerequisites++;
                    }
                    else
                    {
                        script.WriteLine(stateOfTheWorld.configuration.binariesDirectory + "snap.exe single " + stateOfTheWorld.configuration.localIndexDirectory + " " + stateOfTheWorld.configuration.chromosomeMapsDirectory + "chr" + chromosomeNumber +
                            ASETools.allLociExtension + " -o " + stateOfTheWorld.configuration.chromosomeMapsDirectory +"chr" + chromosomeNumber + ASETools.allLociAlignedExtension + " -om 3 -omax 1 -x -D 3 -map -mrl 48");
                        nAddedToScript++;
                    }
                }
            }

            public bool EvaluateDependencies(StateOfTheWorld stateOfTheWorld) { return true; }   // Really no dependencies for this unless the reference genome changes, since the input files are effectively constants

        } // AlignAllLociProcessingStage

        class RepetitveRegionMapProcessingStage : ProcessingStage
        {
            public RepetitveRegionMapProcessingStage() { }

            public string GetStageName() { return "Make Repetitive Region Map"; }
            public bool NeedsCases() { return false; }
            public void EvaluateStage(StateOfTheWorld stateOfTheWorld, StreamWriter script, ASETools.RandomizingStreamWriter hpcScript, StreamWriter linuxScript, StreamWriter azureScript, out List<string> filesToDownload, out int nDone, out int nAddedToScript, out int nWaitingForPrerequisites)
            {
                nDone = 0;
                nAddedToScript = 0;
                nWaitingForPrerequisites = 0;
                filesToDownload = null;
                
                if (File.Exists(stateOfTheWorld.configuration.redundantChromosomeRegionFilename))
                {
                    nDone++;
                    return;
                }

                for (int chromosomeNumber = 1; chromosomeNumber <= ASETools.nHumanAutosomes; chromosomeNumber++)
                {
                    if (!File.Exists(stateOfTheWorld.configuration.chromosomeMapsDirectory + "chr" + chromosomeNumber + ASETools.allLociAlignedExtension))
                    {
                        nWaitingForPrerequisites++;
                        return;
                    }
                }

                script.WriteLine(stateOfTheWorld.configuration.binariesDirectory + "MakeRepetitiveRegionMap.exe");
                nAddedToScript++;
            }

            public bool EvaluateDependencies(StateOfTheWorld stateOfTheWorld) { return true; }   // Really no dependencies for this unless the reference genome changes, since the input files are effectively constants
        } // RepetitveRegionMapProcessingStage


        class FPKMProcessingStage : ProcessingStage
		{
			public FPKMProcessingStage() { }

			public string GetStageName() { return "Process FPKM Data"; }

			public bool NeedsCases() { return true; }

			public void EvaluateStage(StateOfTheWorld stateOfTheWorld, StreamWriter script, ASETools.RandomizingStreamWriter hpcScript, StreamWriter linuxScript, StreamWriter azureScript, out List<string> filesToDownload, out int nDone, out int nAddedToScript, out int nWaitingForPrerequisites)
			{
				filesToDownload = new List<string>();
				nDone = 0;
				nAddedToScript = 0;
				nWaitingForPrerequisites = 0;

				foreach (var caseEntry in stateOfTheWorld.cases)
				{
					var case_ = caseEntry.Value;

					if (case_.normal_fpkm_filename == "" && case_.normal_fpkm_file_id != "")
					{
						nWaitingForPrerequisites = 1;
						filesToDownload.Add(case_.normal_fpkm_file_id);
					}
					
					if (case_.tumor_fpkm_filename == "" && case_.tumor_fpkm_file_id != "")
					{
						nWaitingForPrerequisites = 1;
						filesToDownload.Add(case_.tumor_fpkm_file_id);
					}

					// nAddedToScript = 1; // Commented out until there's something to add to the script.
				}

			} // EvaluateStage

			public bool EvaluateDependencies(StateOfTheWorld stateOfTheWorld)
			{
				DateTime oldestFile;
				if (stateOfTheWorld.scatterGraphsSummaryFile == "")
				{
					return true;
				}

				oldestFile = new FileInfo(stateOfTheWorld.scatterGraphsSummaryFile).LastWriteTime;
				foreach (var selectedGene in stateOfTheWorld.selectedGenes)
				{
					if (!stateOfTheWorld.scatterGraphsByHugoSymbol.ContainsKey(selectedGene.Hugo_Symbol))
					{
						return true;
					}

					var date = new FileInfo(stateOfTheWorld.scatterGraphsByHugoSymbol[selectedGene.Hugo_Symbol]).LastWriteTime;
					if (date < oldestFile)
					{
						oldestFile = date;
					}
				}

				//
				// Now look at the dependencies.
				//
				foreach (var caseEntry in stateOfTheWorld.cases)
				{
					var case_ = caseEntry.Value;

					if (case_.annotated_selected_variants_filename == "" || new FileInfo(case_.annotated_selected_variants_filename).LastWriteTime > oldestFile)
					{
						Console.WriteLine("Case " + case_.case_id + " either doesn't have an annotated selected variants file, or it is newer than a gene scatter graph file that depends on it.");
						return false;
					}

					if (case_.tumor_dna_mapped_base_count_filename == "" || new FileInfo(case_.tumor_dna_mapped_base_count_filename).LastWriteTime > oldestFile ||
						case_.tumor_rna_mapped_base_count_filename == "" || new FileInfo(case_.tumor_rna_mapped_base_count_filename).LastWriteTime > oldestFile)
					{
						Console.WriteLine("Case " + case_.case_id + " either doesn't have a tumor mapped base count file, or it is newer than a gene scatter graph file that depends on it.");
						return false;
					}
				}

				return true;
			} // EvaluateDependencies

		} // FPKMProcessingStage

		//
		// This represents the state of the world.  Processing stages look at this state and generate actions to move us along.
		//
		class StateOfTheWorld
        {
            public StateOfTheWorld(ASETools.Configuration configuration_) 
            {
                configuration = configuration_;
            }

            public ASETools.Configuration configuration;
            public Dictionary<string, ASETools.DownloadedFile> downloadedFiles = null;
            public Dictionary<string, List<ASETools.DerivedFile>> derivedFiles = null;
            public Dictionary<string, FileInfo> expressionFiles = null;            
            public Dictionary<string, ASETools.MAFInfo> mafInfo = null;
            public Dictionary<string, ASETools.Case> cases = null;
            public List<string> diseases = null;
            public Dictionary<string, string> fileIdToCaseId = null;
            public Dictionary<string, long> fileSizesFromGDC = null;
            public List<ASETools.SelectedGene> selectedGenes = null;
            public string scatterGraphsSummaryFile = "";
            public Dictionary<string, string> scatterGraphsByHugoSymbol = new Dictionary<string, string>();

            public void DetermineTheStateOfTheWorld()
            {
                ASETools.ScanFilesystems(configuration, out downloadedFiles, out derivedFiles);

                if (File.Exists(configuration.selectedGenesFilename))
                {
                    selectedGenes = ASETools.SelectedGene.LoadFromFile(configuration.selectedGenesFilename);

                    foreach (var selectedGene in selectedGenes)
                    {
                        if (File.Exists(configuration.geneScatterGraphsDirectory + selectedGene.Hugo_Symbol + ".txt"))
                        {
                            scatterGraphsByHugoSymbol.Add(selectedGene.Hugo_Symbol, configuration.geneScatterGraphsDirectory + selectedGene.Hugo_Symbol + ".txt");
                        }
                    }
                }

                if (File.Exists(configuration.geneScatterGraphsDirectory + ASETools.scatterGraphsSummaryFilename))
                {
                    scatterGraphsSummaryFile = configuration.geneScatterGraphsDirectory + ASETools.scatterGraphsSummaryFilename;
                }

                mafInfo = ASETools.MAFInfo.LoadMAFManifest(configuration.mafManifestPathname);
                cases = ASETools.Case.LoadCases(configuration.casesFilePathname);

                if (null != cases)
                {
                    diseases = new List<string>();

                    foreach (var caseEntry in cases)
                    {
                        var case_ = caseEntry.Value;

                        if (!diseases.Contains(case_.disease()))
                        {
                            diseases.Add(case_.disease());
                        }
                    }

                    fileIdToCaseId = new Dictionary<string, string>();

                    foreach (var caseEntry in cases)
                    {
                        var case_ = caseEntry.Value;

                        fileIdToCaseId.Add(case_.tumor_dna_file_id, case_.case_id);
                        fileIdToCaseId.Add(case_.tumor_rna_file_id, case_.case_id);
                        fileIdToCaseId.Add(case_.normal_dna_file_id, case_.case_id);
                        if (null != case_.normal_rna_file_id && "" != case_.normal_rna_file_id)
                        {
                            fileIdToCaseId.Add(case_.normal_rna_file_id, case_.case_id);
                        }
                        if (null != case_.tumor_methylation_file_id && "" != case_.tumor_methylation_file_id)
                        {
                            fileIdToCaseId.Add(case_.tumor_methylation_file_id, case_.case_id);
                        }
						if (null != case_.normal_methylation_file_id && "" != case_.normal_methylation_file_id)
						{
							fileIdToCaseId.Add(case_.normal_methylation_file_id, case_.case_id);
						}
					}

                    //
                    // Check that the derived file cases are real cases.
                    //

                    foreach (var derivedFileCaseEntry in derivedFiles)
                    {
                        var caseId = derivedFileCaseEntry.Key;
                        var derivedFilesForThisCase = derivedFileCaseEntry.Value;

                        if (cases.ContainsKey(caseId))
                        {
                            continue;
                        }

                        Console.Write("There's a derived files directory for case id " + caseId + ", which isn't a known case.  It contains:");
                        foreach (var badDrivedFile in derivedFilesForThisCase)
                        {
                            Console.Write(" " + badDrivedFile.fileinfo.FullName);
                            if (fileIdToCaseId.ContainsKey(badDrivedFile.derived_from_file_id))
                            {
                                Console.WriteLine(" (derived from a fileID associated with case " + fileIdToCaseId[badDrivedFile.derived_from_file_id] + ")");
                            }
                        }
                        Console.WriteLine();
                    }


                    ASETools.Case.loadAllFileLocations(cases, downloadedFiles, derivedFiles);

                    int nNormalDNA = 0, nTumorDNA = 0, nNormalRNA = 0, nTumorRNA = 0, nMethylation = 0, nCopyNumber = 0;
                    ulong bytesNormalDNA = 0, bytesTumorDNA = 0, bytesNormalRNA = 0, bytesTumorRNA = 0, bytesMethylation = 0, bytesCopyNumber = 0;

                    foreach (var caseEntry in cases)
                    {
                        var case_ = caseEntry.Value;

                        if (downloadedFiles.ContainsKey(case_.normal_dna_file_id))
                        {
                            nNormalDNA++;
                            bytesNormalDNA += (ulong)downloadedFiles[case_.normal_dna_file_id].fileInfo.Length;
                        }

                        if (downloadedFiles.ContainsKey(case_.tumor_dna_file_id))
                        {
                            nTumorDNA++;
                            bytesTumorDNA += (ulong)downloadedFiles[case_.tumor_dna_file_id].fileInfo.Length;
                        }

                        if (downloadedFiles.ContainsKey(case_.normal_rna_file_id))
                        {
                            nNormalRNA++;
                            bytesNormalRNA += (ulong)downloadedFiles[case_.normal_rna_file_id].fileInfo.Length;
                        }

                        if (downloadedFiles.ContainsKey(case_.tumor_rna_file_id))
                        {
                            nTumorRNA++;
                            bytesTumorRNA += (ulong)downloadedFiles[case_.tumor_rna_file_id].fileInfo.Length;
                        }

                        if (downloadedFiles.ContainsKey(case_.normal_methylation_file_id))
                        {
                            nMethylation++;
                            bytesMethylation += (ulong)downloadedFiles[case_.normal_methylation_file_id].fileInfo.Length;
                        }

						if (downloadedFiles.ContainsKey(case_.tumor_methylation_file_id))
						{
							nMethylation++;
							bytesMethylation += (ulong)downloadedFiles[case_.tumor_methylation_file_id].fileInfo.Length;
						}

						if (downloadedFiles.ContainsKey(case_.normal_copy_number_file_id))
                        {
                            nCopyNumber++;
                            bytesCopyNumber += (ulong)downloadedFiles[case_.normal_copy_number_file_id].fileInfo.Length;
                        }

						if (downloadedFiles.ContainsKey(case_.tumor_copy_number_file_id))
						{
							nCopyNumber++;
							bytesCopyNumber += (ulong)downloadedFiles[case_.tumor_copy_number_file_id].fileInfo.Length;
						}
					} // foreach case


                    Console.WriteLine(nNormalDNA + "(" + ASETools.SizeToUnits(bytesNormalDNA) + "B) normal DNA, " + nTumorDNA + "(" + ASETools.SizeToUnits(bytesTumorDNA) + "B) tumor DNA, " +
                                      nNormalRNA + "(" + ASETools.SizeToUnits(bytesNormalRNA) + "B) normal RNA, " + nTumorRNA + "(" + ASETools.SizeToUnits(bytesTumorRNA) + "B) tumor RNA, " +
                                      nMethylation + "(" + ASETools.SizeToUnits(bytesMethylation) + "B) methylation, " + nCopyNumber + "(" + ASETools.SizeToUnits(bytesCopyNumber) + "B) copy number");

                } // If we loaded cases

                expressionFiles = new Dictionary<string, FileInfo>();

                if (Directory.Exists(configuration.expressionFilesDirectory))
                {
                    foreach (var filename in Directory.EnumerateFiles(configuration.expressionFilesDirectory, "expression_*")) {
                        var disease = filename.Substring(filename.LastIndexOf('_') + 1).ToLower();
                        if (!diseases.Contains(disease))
                        {
                            Console.WriteLine("Found expression file that doesn't seem to correspond to a disease: " + filename);
                        }
                        else
                        {
                            expressionFiles.Add(disease, new FileInfo(filename));
                        }
                    }
                }

                fileSizesFromGDC = new Dictionary<string, long>();

                foreach (var caseEntry in cases)
                {
                    var case_ = caseEntry.Value;

                    fileSizesFromGDC.Add(case_.normal_dna_file_id, case_.normal_dna_size);
                    fileSizesFromGDC.Add(case_.tumor_dna_file_id, case_.tumor_dna_size);
                    fileSizesFromGDC.Add(case_.tumor_rna_file_id, case_.tumor_rna_size);

                    if (case_.normal_rna_file_id != "")
                    {
                        fileSizesFromGDC.Add(case_.normal_rna_file_id, case_.normal_rna_size);
                    }

                    if (case_.tumor_methylation_file_id != "")
                    {
                        fileSizesFromGDC.Add(case_.tumor_methylation_file_id, case_.tumor_methylation_size);
                    }

					if (case_.normal_methylation_file_id != "")
					{
						fileSizesFromGDC.Add(case_.normal_methylation_file_id, case_.normal_methylation_size);
					}

					if (case_.tumor_copy_number_file_id != "")
                    {
                        fileSizesFromGDC.Add(case_.tumor_copy_number_file_id, case_.tumor_copy_number_size);
                    }

					if (case_.normal_copy_number_file_id != "")
					{
						fileSizesFromGDC.Add(case_.normal_copy_number_file_id, case_.normal_copy_number_size);
					}

					if (case_.tumor_fpkm_file_id != "")
					{
						fileSizesFromGDC.Add(case_.tumor_fpkm_file_id, case_.tumor_fpkm_size);
					}

					if (case_.normal_fpkm_file_id != "")
					{
						fileSizesFromGDC.Add(case_.normal_fpkm_file_id, case_.normal_fpkm_size);
					}
				}

            }

            public bool fileDownloadedAndVerified(string file_id, string expectedMD5)
            {
                return downloadedFiles.ContainsKey(file_id) && (null == expectedMD5 || "" == expectedMD5 || downloadedFiles[file_id].storedMD5 == expectedMD5);
            }

            public bool containsDerivedFile(string case_id, string derived_from_file_id, ASETools.DerivedFile.Type type)
            {
                return derivedFiles.ContainsKey(case_id) && derivedFiles[case_id].Where(x => x.derived_from_file_id == derived_from_file_id && x.type == type).Count() != 0;
            }

            public ASETools.DerivedFile getDrivedFile(string case_id, string derived_from_file_id, ASETools.DerivedFile.Type type)
            {
                if (!derivedFiles.ContainsKey(case_id))
                {
                    return null;
                }

                var set = derivedFiles[case_id].Where(x => x.derived_from_file_id == derived_from_file_id && x.type == type);

                if (set.Count() == 0)
                {
                    return null;
                }

                return set.ToList()[0];
            }
        } // StateOfTheWorld


        static void Main(string[] args)
        {
            var stopwatch = new Stopwatch();
            stopwatch.Start();

            var configuration = ASETools.Configuration.loadFromFile(args);

            if (null == configuration)
            {
                Console.WriteLine("Giving up because we were unable to load configuration.");
                return;
            }

            if (configuration.commandLineArgs.Count() > 1 || configuration.commandLineArgs.Count() == 1 && configuration.commandLineArgs[0] != "-d")
            {
                Console.WriteLine("usage: ASEProcessManager {-configuration configurationFilename} {-d}");
                Console.WriteLine("-d means to check dependencies.");
            }
            
            //
            // Delete any existing scripts.
            //
            File.Delete(configuration.scriptOutputDirectory + scriptFilename);
            File.Delete(configuration.scriptOutputDirectory + linuxScriptFilename);
            File.Delete(configuration.scriptOutputDirectory + downloadScriptFilename);
            if (configuration.hpcScriptFilename != "")
            {
                File.Delete(configuration.scriptOutputDirectory + configuration.hpcScriptFilename);
            }
            if (configuration.azureScriptFilename != "")
            {
                File.Delete(configuration.scriptOutputDirectory + configuration.azureScriptFilename);
            }

            bool checkDependencies = configuration.commandLineArgs.Count() >= 1 && configuration.commandLineArgs.Contains("-d");

            var stateOfTheWorld = new StateOfTheWorld(configuration);
            stateOfTheWorld.DetermineTheStateOfTheWorld();

            jobAddString = @"job add %1 /exclusive /numnodes:1-1 /scheduler:" + stateOfTheWorld.configuration.hpcScheduler + " ";
            
            Console.WriteLine();


            if (null != stateOfTheWorld.cases)
            {
                //
                // Rewrite the cases file, since we have just updated all of the file locations for downloaded and derived files.
                //
                ASETools.Case.SaveToFile(stateOfTheWorld.cases, configuration.casesFilePathname);
            }

            var script = ASETools.CreateStreamWriterWithRetry(configuration.scriptOutputDirectory + scriptFilename);

            if (configuration.completedVCFsDirectory != "")
            {
                //
                // Check to see if there are any completed VCFs (downloaded from Auzure) that need to be moved.
                //

                var casesByNormalDNAId = new Dictionary<string, ASETools.Case>();
                foreach (var caseEntry in stateOfTheWorld.cases) {
                    var case_ = caseEntry.Value;
                    casesByNormalDNAId.Add(case_.normal_dna_file_id, case_);
                }


                var vcfsToBeMoved = new List<string>();
                foreach (var completedVCF in Directory.EnumerateFiles(configuration.completedVCFsDirectory)) {
                    if (!completedVCF.EndsWith(ASETools.vcfExtension))
                    {
                        Console.WriteLine("Found non-VCF file in completed VCFs directory: " + completedVCF + ".  Ignoring.");
                        continue;
                    }

                    string fileId = ASETools.GetFileIdFromPathname(completedVCF);
                    if (!casesByNormalDNAId.ContainsKey(fileId)) {
                        Console.WriteLine("completed VCFs directory contains a file that doesn't seem to correspond to a normal DNA file id: " + completedVCF + ".  Ignoring.");
                        continue;
                    }

                    vcfsToBeMoved.Add(completedVCF);
                }

                if (vcfsToBeMoved.Count() > 0) {
                    var completedVCFsPathComponents = configuration.completedVCFsDirectory.Split('\\');

                    if (completedVCFsPathComponents.Count() < 2)
                    {
                        Console.WriteLine("The completed VCF directory in the configuration should be a pathname: " + configuration.completedVCFsDirectory);
                        return;
                    }


                    bool failed = false;
                    string [] dataPathComponents = null;
                    int completedComponentsToSkip = configuration.completedVCFsDirectory.EndsWith(@"\") ? 2 : 1;

                    foreach (var dataDirectory in configuration.dataDirectories)
                    {
                        dataPathComponents = dataDirectory.Split('\\');

                        failed = false;
                        for (int i = 0; i < completedVCFsPathComponents.Count() - completedComponentsToSkip; i++)
                        {
                            if (dataPathComponents[i] != completedVCFsPathComponents[i])
                            {
                                failed = true;
                                break;
                            }
                        }

                        if (!failed) {
                            break;
                        }
                    } // foreach data directory

                    if (failed) {
                        Console.WriteLine("Unable to find destination for completed VCFs (the completed VCFs directory doesn't share a parent with any data directory, and it must.)  Look at your configuration file.");
                        return;
                    }


                    string destinationDirectory = dataPathComponents[0];

                    for (int i = 1; i < completedVCFsPathComponents.Count() - completedComponentsToSkip; i++)
                    {
                        destinationDirectory += '\\' + dataPathComponents[i];
                    }

                    destinationDirectory += '\\' + configuration.derivedFilesDirectory + '\\';

                    foreach (var completedVCF in vcfsToBeMoved)
                    {
                        string normalDNAFileId = ASETools.GetFileIdFromPathname(completedVCF);

                        var case_ = casesByNormalDNAId[normalDNAFileId];
                        script.WriteLine("md " + destinationDirectory + case_.case_id);
                        script.WriteLine("mv " + completedVCF + " " + destinationDirectory + case_.case_id + @"\" + ASETools.GetFileNameFromPathname(completedVCF));
                    }

                    Console.WriteLine("Added " + vcfsToBeMoved.Count() + " vcfs to be moved from the completed_vcfs directory to their final locations.");
                }// If we had any completed VCFs to be moved.
            } // if we have a completed VCFs directory


            List<ProcessingStage> processingStages = new List<ProcessingStage>();

			var forAlleleSpecificExpression = true;

			processingStages.Add(new MAFConfigurationProcessingStage());
			processingStages.Add(new GenerateCasesProcessingStage());
			processingStages.Add(new AllcountProcesingStage());
			processingStages.Add(new DownloadProcessingStage());
			processingStages.Add(new MD5ComputationProcessingStage());
			processingStages.Add(new GermlineVariantCallingProcessingStage());
			processingStages.Add(new SelectVariantsProcessingStage());
			processingStages.Add(new AnnotateVariantsProcessingStage());
			processingStages.Add(new ExpressionDistributionProcessingStage());
			processingStages.Add(new ExtractMAFLinesProcessingStage());
			processingStages.Add(new RegionalExpressionProcessingStage());
			processingStages.Add(new ExpressionNearMutationsProcessingStage(forAlleleSpecificExpression));
			processingStages.Add(new ExtractReadsProcessingStage());
			processingStages.Add(new SelectGenesProcessingStage());
			processingStages.Add(new CountMappedBasesProcessingStage());
			processingStages.Add(new GenerateScatterGraphsProcessingStage());
			//processingStages.Add(new MethylationProcessingStage());
			processingStages.Add(new FPKMProcessingStage());
            processingStages.Add(new GenerateAllLociProcessingStage());
            processingStages.Add(new AlignAllLociProcessingStage());
            processingStages.Add(new RepetitveRegionMapProcessingStage());

            if (checkDependencies)
            {
                bool allDependenciesOK = true;
                foreach (var processingStage in processingStages)
                {
                    if (stateOfTheWorld.cases != null || !processingStage.NeedsCases())
                    {
                        allDependenciesOK &= processingStage.EvaluateDependencies(stateOfTheWorld);
                    }
                }

                if (!allDependenciesOK)
                {
                    Console.WriteLine("Not generating scripts because some dependencies have been violated.  Delete the stale generated files and rerun.");
                    return;
                }
            }

            ASETools.RandomizingStreamWriter hpcScript;
            StreamWriter azureScript;

            if (configuration.hpcScriptFilename == "")  // The empty string means not to generate an output.  We do this by making a Null stream.
            {
                hpcScript = new ASETools.RandomizingStreamWriter(new StreamWriter(Stream.Null));
            }
            else
            {
                hpcScript = new ASETools.RandomizingStreamWriter(ASETools.CreateStreamWriterWithRetry(configuration.scriptOutputDirectory + configuration.hpcScriptFilename));
            }


            if (configuration.azureScriptFilename == "")
            {
                azureScript = new StreamWriter(Stream.Null);
            }
            else
            {
                azureScript = ASETools.CreateStreamWriterWithRetry(configuration.scriptOutputDirectory + configuration.azureScriptFilename);
            }

            var linuxScript = ASETools.CreateStreamWriterWithRetry(configuration.scriptOutputDirectory + linuxScriptFilename);

            var allFilesToDownload = new List<string>();

            int longestStageName = 0;

            foreach (var processingStage in processingStages)
            {
                longestStageName = Math.Max(longestStageName, processingStage.GetStageName().Count());
            }

            const string stageNameHeader = "Stage Name";

            Console.Write(stageNameHeader);
            int paddingSize = Math.Max(0, longestStageName - stageNameHeader.Count());
            for (int i = 0; i < paddingSize; i++)
            {
                Console.Write(" ");
            }

            Console.WriteLine(" # Done  # Added  # Waiting  # Downloads");

            for (int i = 0; i < stageNameHeader.Count() + paddingSize; i++) {
                Console.Write("-");
            }
            Console.WriteLine(" ------  -------  ---------  -----------");



            foreach (var processingStage in processingStages)
            {
                int nDone;
                int nAddedToScript;
                int nWaitingForPrerequisites;
                List<string> stageFilesToDownload;

                if (stateOfTheWorld.cases != null || !processingStage.NeedsCases())
                {
                    processingStage.EvaluateStage(stateOfTheWorld, script, hpcScript, linuxScript, azureScript, out stageFilesToDownload, out nDone, out nAddedToScript, out nWaitingForPrerequisites);
                }
                else
                {
                    nDone = 0;
                    nAddedToScript = 0;
                    nWaitingForPrerequisites = 1;
                    stageFilesToDownload = null;
                }

                int nDownloadsRequested = 0;
                if (null != stageFilesToDownload)
                {
                    foreach (var file in stageFilesToDownload)
                    {
                        if (!allFilesToDownload.Contains(file))
                        {
                            nDownloadsRequested++;
                            allFilesToDownload.Add(file);
                        }
                    }
                }

                Console.WriteLine(String.Format("{0," + (stageNameHeader.Count() + paddingSize) + "}", processingStage.GetStageName()) + " " + String.Format("{0,6}", nDone) + " " + String.Format("{0,8}", nAddedToScript) + " " +
                    String.Format("{0,10}", nWaitingForPrerequisites) + " " + String.Format("{0,11}", nDownloadsRequested));
            } // foreach stage

            //
            // Now put downloads in their own script. They're separated out because they need to be run in one of the download
            // directories, and the user may want to split them across machines.
            //

            long bytesToDownload = 0;
            if (allFilesToDownload.Count() == 0)
            {
                File.Delete(configuration.scriptOutputDirectory + downloadScriptFilename);
            } 
            else
            {
                var downloadScript = ASETools.CreateStreamWriterWithRetry(configuration.scriptOutputDirectory + downloadScriptFilename);

                foreach (var file in allFilesToDownload)
                {
                    downloadScript.WriteLine(configuration.binariesDirectory + "gdc-client download --no-file-md5sum --token-file " + configuration.accessTokenPathname + " " + file);    // use the no MD5 sum option because we compute it ourselves later (and all a failed check does is print a message that we'll never see anyway)
                    bytesToDownload += stateOfTheWorld.fileSizesFromGDC[file];
                }

                downloadScript.Close();
            }

            script.Close();
            hpcScript.Close();
            linuxScript.Close();
            azureScript.Close();

            Console.WriteLine();
            Console.WriteLine("Downloading " + ASETools.SizeToUnits((ulong)bytesToDownload) + "B in " + allFilesToDownload.Count() + " files.");
            Console.WriteLine("ASEProcessManager took " + ASETools.ElapsedTimeInSeconds(stopwatch) + " and finished at " + DateTime.Now.ToLocalTime().ToString());
        }
    }
}
