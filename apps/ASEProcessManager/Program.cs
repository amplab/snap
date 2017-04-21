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

        //
        // A ProcessingStage examines the state of the world and if the prerequisites are there and the step isn't completed adds commands to the script to do some
        // work and advance the state of the world.
        //
        interface ProcessingStage
        {
            string GetStageName();
            void EvaluateStage(StateOfTheWorld stateOfTheWorld, StreamWriter script, StreamWriter hpcScript, StreamWriter linuxScript, out List<string> filesToDownload, out int nDone, out int nAddedToScript, out int nWaitingForPrerequisites);
            bool EvaluateDependencies(StateOfTheWorld stateOfTheWorld);
        }

        class MAFConfigurationProcessingStage : ProcessingStage 
        {
            public MAFConfigurationProcessingStage() { }

            public string GetStageName()
            {
                return "Generate MAF Configuration";
            }

            public void EvaluateStage(StateOfTheWorld stateOfTheWorld, StreamWriter script, StreamWriter hpcScript, StreamWriter linuxScript, out List<string> filesToDownload, out int nDone, out int nAddedToScript, out int nWaitingForPrerequisites)
            {
                filesToDownload = null;

                nWaitingForPrerequisites = 0; // This is the very first thing we do, there are never any prerequisites

                if (stateOfTheWorld.mafInfo != null)
                {
                    nDone = 1;
                    nAddedToScript = 0;
                    return;
                }

                script.WriteLine(stateOfTheWorld.configuration.binaryDirectory + "GenerateMAFConfiguration -configuration " + stateOfTheWorld.configuration.configuationFilePathname);

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

            public void EvaluateStage(StateOfTheWorld stateOfTheWorld, StreamWriter script, StreamWriter hpcScript, StreamWriter linuxScript, out List<string> filesToDownload, out int nDone, out int nAddedToScript, out int nWaitingForPrerequisites)
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
                    script.WriteLine(stateOfTheWorld.configuration.binaryDirectory + "GenerateCases -configuration " + stateOfTheWorld.configuration.configuationFilePathname);
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

            public void EvaluateStage(StateOfTheWorld stateOfTheWorld, StreamWriter script, StreamWriter hpcScript, StreamWriter linuxScript, out List<string> filesToDownload, out int nDone, out int nAddedToScript, out int nWaitingForPrerequisites)
            {
                filesToDownload = null;
                nDone = 0;
                nAddedToScript = 0;
                nWaitingForPrerequisites = 0;

                if (stateOfTheWorld.cases == null)
                {
                    nWaitingForPrerequisites = 1;
                    return;
                }

                foreach (var caseEntry in stateOfTheWorld.cases)
                {
                    var case_ = caseEntry.Value;

                    HandleFile(stateOfTheWorld, case_.tumor_rna_file_id, case_.tumor_rna_file_bam_md5, case_.case_id, ASETools.DerivedFile.Type.TumorRNAAllcount,
                        ASETools.allcountExtension, script, hpcScript, ref filesToDownload, ref nDone, ref nAddedToScript, ref nWaitingForPrerequisites);

                    HandleFile(stateOfTheWorld, case_.normal_dna_file_id, case_.normal_dna_file_bam_md5, case_.case_id, ASETools.DerivedFile.Type.NormalDNAAllcount,
                        ASETools.normalDNAAllcountExtension, script, hpcScript, ref filesToDownload, ref nDone, ref nAddedToScript, ref nWaitingForPrerequisites);

                    HandleFile(stateOfTheWorld, case_.tumor_dna_file_id, case_.tumor_dna_file_bam_md5, case_.case_id, ASETools.DerivedFile.Type.TumorDNAAllcount,
                        ASETools.tumorDNAAllcountExtension, script, hpcScript, ref filesToDownload, ref nDone, ref nAddedToScript, ref nWaitingForPrerequisites);

                } // Foreach case
            }// EvaluateStage

            void HandleFile(StateOfTheWorld stateOfTheWorld, string file_id, string expectedMD5, string case_id, ASETools.DerivedFile.Type type, string extension, StreamWriter script, StreamWriter hpcScript, ref List<string> filesToDownload,ref int nDone, ref int nAddedToScript, ref int nWaitingForPrerequisites)
            {

                if (!stateOfTheWorld.downloadedFiles.ContainsKey(file_id))
                {
                    if (null == filesToDownload)
                    {
                        filesToDownload = new List<string>();
                    }
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
                        nDone++;
                    }
                    else
                    {
                        nAddedToScript++;
                        string caseDirectory = ASETools.GetDirectoryFromPathname(stateOfTheWorld.downloadedFiles[file_id].fileInfo.FullName) + @"\..\..\" + stateOfTheWorld.configuration.derivedFilesDirectory + @"\" + case_id + @"\";
                        script.WriteLine("md " + caseDirectory);
                        script.WriteLine(stateOfTheWorld.configuration.binaryDirectory + "CountReadsCovering " + stateOfTheWorld.configuration.indexDirectory + " -a " + stateOfTheWorld.downloadedFiles[file_id].fileInfo.FullName + " - | gzip -9 > " +
                            caseDirectory + file_id + extension);

                        hpcScript.WriteLine(@"job add %1 /exclusive /numnodes:1-1 /scheduler:" + stateOfTheWorld.configuration.hpcScheduler + " " +
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

            public void EvaluateStage(StateOfTheWorld stateOfTheWorld, StreamWriter script, StreamWriter hpcScript, StreamWriter linuxScript, out List<string> filesToDownload, out int nDone, out int nAddedToScript, out int nWaitingForPrerequisites)
            {
                filesToDownload = new List<string>();
                nDone = 0;
                nAddedToScript = 0;
                nWaitingForPrerequisites = 0;

                if (stateOfTheWorld.cases == null)
                {
                    nWaitingForPrerequisites = 1;
                    return;
                }

                foreach (var caseEntry in stateOfTheWorld.cases)
                {
                    var case_ = caseEntry.Value;

                    string[] idsToDownload = {case_.normal_dna_file_id, case_.tumor_dna_file_id, case_.normal_rna_file_id, case_.tumor_rna_file_id, case_.methylation_file_id, case_.copy_number_file_id};

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

            public void EvaluateStage(StateOfTheWorld stateOfTheWorld, StreamWriter script, StreamWriter hpcScript, StreamWriter linuxScript, out List<string> filesToDownload, out int nDone, out int nAddedToScript, out int nWaitingForPrerequisites)
            {
                filesToDownload = null; // This stage never generates downloads
                nAddedToScript = 0;
                nDone = 0;

                if (stateOfTheWorld.cases == null)
                {
                    nWaitingForPrerequisites = 1;
                    return;
                }

                nWaitingForPrerequisites = 0;

                foreach (var caseEntry in stateOfTheWorld.cases)
                {
                    var case_ = caseEntry.Value;

                    HandleFile(stateOfTheWorld, script, hpcScript, case_.tumor_rna_file_id, case_.tumor_rna_file_bam_md5, ref nDone, ref nAddedToScript, ref nWaitingForPrerequisites);
                    HandleFile(stateOfTheWorld, script, hpcScript, case_.normal_dna_file_id, case_.normal_dna_file_bam_md5, ref nDone, ref nAddedToScript, ref nWaitingForPrerequisites);
                    HandleFile(stateOfTheWorld, script, hpcScript, case_.tumor_dna_file_id, case_.tumor_dna_file_bam_md5, ref nDone, ref nAddedToScript, ref nWaitingForPrerequisites);
                    HandleFile(stateOfTheWorld, script, hpcScript, case_.normal_rna_file_id, case_.normal_rna_file_bam_md5, ref nDone, ref nAddedToScript, ref nWaitingForPrerequisites);
                    if (case_.methylation_file_id != null && case_.methylation_file_id != "")
                    {
                        HandleFile(stateOfTheWorld, script, hpcScript, case_.methylation_file_id, case_.methylation_file_md5, ref nDone, ref nAddedToScript, ref nWaitingForPrerequisites);
                    }
                    if (case_.copy_number_file_id != null && case_.copy_number_file_id != "")
                    {
                        HandleFile(stateOfTheWorld, script, hpcScript, case_.copy_number_file_id, case_.copy_number_file_md5, ref nDone, ref nAddedToScript, ref nWaitingForPrerequisites);
                    }
                }
            } // EvaluateStage

            void HandleFile(StateOfTheWorld stateOfTheWorld, StreamWriter script, StreamWriter hpcScript, string fileId, string expectedMD5, ref int nDone, ref int nAddedToScript, ref int nWaitingForPrerequisites)
            {
                if (!stateOfTheWorld.downloadedFiles.ContainsKey(fileId) || null == expectedMD5 || "" == expectedMD5 || stateOfTheWorld.downloadedFiles[fileId].fileInfo.FullName.ToLower().Contains("partial"))
                {
                    nWaitingForPrerequisites++;
                    return;
                }

                var downloadedFile = stateOfTheWorld.downloadedFiles[fileId];

                if (downloadedFile.storedMD5 != null && downloadedFile.storedMD5 != "")
                {
                    nDone++;

                    if (downloadedFile.storedMD5 != expectedMD5)
                    {
                        Console.WriteLine("MD5 checksum mismatch on file " + downloadedFile.fileInfo.FullName + " " + downloadedFile.storedMD5 + " != " + expectedMD5);
                    }

                    return;
                }

                script.WriteLine(stateOfTheWorld.configuration.binaryDirectory + "ComputeMD5 " + downloadedFile.fileInfo.FullName + " > " + downloadedFile.fileInfo.FullName + ".md5");
                hpcScript.WriteLine(@"job add %1 /numnodes:1-1 /scheduler:" + stateOfTheWorld.configuration.hpcScheduler + " " + stateOfTheWorld.configuration.hpcBinariesDirectory + "ComputeMD5IntoFile.cmd " +
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

            public void EvaluateStage(StateOfTheWorld stateOfTheWorld, StreamWriter script, StreamWriter hpcScript, StreamWriter linuxScript, out List<string> filesToDownload, out int nDone, out int nAddedToScript, out int nWaitingForPrerequisites)
            {
                filesToDownload = new List<string>();
                nDone = 0;
                nAddedToScript = 0;
 
                if (stateOfTheWorld.cases == null)
                {
                    nWaitingForPrerequisites = 1;
                    return;
                }
                nWaitingForPrerequisites = 0;

                foreach (var caseEntry in stateOfTheWorld.cases)
                {
                    var case_ = caseEntry.Value;

                    if (!stateOfTheWorld.fileDownloadedAndVerified(case_.normal_dna_file_id, case_.normal_dna_file_bam_md5))
                    {
                        nWaitingForPrerequisites++;
                        continue;
                    }

                    if (stateOfTheWorld.containsDerivedFile(case_.case_id, case_.normal_dna_file_id, ASETools.DerivedFile.Type.VCF))
                    {
                        nDone++;
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
                // Placeholder.
                //
                return true;
            } // EvaluateDependencies

        }  // GermlineVariantCallingProcessingStage



        class SelectVariantsProcessingStage : ProcessingStage
        {
            public SelectVariantsProcessingStage() { }

            public string GetStageName()
            {
                return "Select Germline Variants";
            }

            public void EvaluateStage(StateOfTheWorld stateOfTheWorld, StreamWriter script, StreamWriter hpcScript, StreamWriter linuxScript, out List<string> filesToDownload, out int nDone, out int nAddedToScript, out int nWaitingForPrerequisites)
            {
                nDone = 0;
                nAddedToScript = 0;
                filesToDownload = null;

                if (stateOfTheWorld.cases == null)
                {
                    nWaitingForPrerequisites = 1;
                    return;
                }

                nWaitingForPrerequisites = 0;

                int nOnCurrentLine = 0;

                foreach (var caseEntry in stateOfTheWorld.cases)
                {
                    var case_ = caseEntry.Value;


                    if (case_.selected_variants_filename != null && case_.selected_variants_filename != "")
                    {
                        nDone++;
                    }

                    if (case_.vcf_filename == "" || case_.tumor_rna_allcount_filename == "" || case_.tumor_dna_allcount_filename == "")
                    {
                        nWaitingForPrerequisites++;
                        continue;
                    }

                    nAddedToScript++;

                    if (nOnCurrentLine >= 800)
                    {
                        script.WriteLine();
                        hpcScript.WriteLine();
                        nOnCurrentLine = 0;
                    }

                    if (nOnCurrentLine == 0) 
                    {
                        script.Write(stateOfTheWorld.configuration.binaryDirectory + "SelectGermlineVariants.exe");
                        hpcScript.Write(stateOfTheWorld.configuration.hpcBinariesDirectory + "SelectGermlineVariants.exe");
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
                //
                // Placeholder.
                //
                return true;
            } // EvaluateDependencies

        } // SelectVariantsProcessingStage

        class ExpressionDistributionProcessingStage : ProcessingStage
        {
            public ExpressionDistributionProcessingStage() { }

            public string GetStageName()
            {
                return "Per-disease mRNA expression distribution";
            }

            public void EvaluateStage(StateOfTheWorld stateOfTheWorld, StreamWriter script, StreamWriter hpcScript, StreamWriter linuxScript, out List<string> filesToDownload, out int nDone, out int nAddedToScript, out int nWaitingForPrerequisites)
            {
                nDone = 0;
                nAddedToScript = 0;
                nWaitingForPrerequisites = 0;
                filesToDownload = null;

                if (stateOfTheWorld.cases == null || stateOfTheWorld.diseases == null)
                {
                    nWaitingForPrerequisites = 1;
                    return;
                }

                bool missingAny = false;

                foreach (var disease in stateOfTheWorld.diseases)
                {
                    if (!stateOfTheWorld.expressionFiles.ContainsKey(disease))
                    {
                        missingAny = true;
                        break;
                    }
                }


                if (!missingAny)
                {
                    nDone = 1;
                    return;
                }

                foreach (var caseEntry in stateOfTheWorld.cases)
                {
                    var case_ = caseEntry.Value;

                    if (!stateOfTheWorld.containsDerivedFile(case_.case_id, case_.tumor_rna_file_id, ASETools.DerivedFile.Type.TumorRNAAllcount))
                    {
                        nWaitingForPrerequisites = 1;
                        return;
                    }
                }

                script.WriteLine(stateOfTheWorld.configuration.binaryDirectory + "ExpressionDistibution.exe " + stateOfTheWorld.configuration.casesFilePathname + " " +
                    stateOfTheWorld.configuration.expressionFilesDirectory);

                hpcScript.WriteLine(stateOfTheWorld.configuration.hpcBinariesDirectory + "ExpressionDistibution.exe " + stateOfTheWorld.configuration.casesFilePathname + " " +
                    stateOfTheWorld.configuration.expressionFilesDirectory);

                nAddedToScript++;
            }

            public bool EvaluateDependencies(StateOfTheWorld stateOfTheWorld)
            {
                //
                // Placeholder.
                //
                return true;
            } // EvaluateDependencies

        } // ExpressionDistributionProcessingStage

        class ExtractMAFLinesProcessingStage : ProcessingStage
        {
            public string GetStageName() { return "Extract MAF Lines"; }

            public void EvaluateStage(StateOfTheWorld stateOfTheWorld, StreamWriter script, StreamWriter hpcScript, StreamWriter linuxScript, out List<string> filesToDownload, out int nDone, out int nAddedToScript, out int nWaitingForPrerequisites)
            {
                filesToDownload = null;

                nWaitingForPrerequisites = 0;
                nDone = 0;
                nAddedToScript = 0;

                if (stateOfTheWorld.cases == null)
                {
                    nWaitingForPrerequisites = 1;
                    return;
                }


                foreach (var caseEntry in stateOfTheWorld.cases)
                {
                    var case_ = caseEntry.Value;

                    if (case_.maf_filename == null || case_.maf_filename == "")
                    {
                        nWaitingForPrerequisites++;
                    }
                    else if (case_.extracted_maf_lines != null && case_.extracted_maf_lines != "")
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
                    script.WriteLine(stateOfTheWorld.configuration.binaryDirectory + "ExtractMAFLines");
                    hpcScript.WriteLine(@"job add %1 /exclusive /numnodes:1-1 /scheduler:" + stateOfTheWorld.configuration.hpcScheduler + " " + stateOfTheWorld.configuration.hpcBinariesDirectory + "ExtractMAFLines");
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

            public string getStageName() { return "Regional Expression"; }

            public void EvaluateStage(StateOfTheWorld stateOfTheWorld, StreamWriter script, StreamWriter hpcScript, StreamWriter linuxScript, out List<string> filesToDownload, out int nDone, out int nAddedToScript, out int nWaitingForPrerequisites)
            {
                nDone = 0;
                nAddedToScript = 0;
                nWaitingForPrerequisites = 0;
                filesToDownload = null;

                if (stateOfTheWorld.cases == null)
                {
                    nWaitingForPrerequisites = 1;
                    return;
                }

                var casesReadyToGoByDisease = new Dictionary<string, List<ASETools.Case>>();
                const int maxCasesPerCommandLine = 800;

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

            void WriteScripts(StateOfTheWorld stateOfTheWorld, List<ASETools.Case> cases, StreamWriter script, StreamWriter hpcScript)
            {
                script.Write(stateOfTheWorld.configuration.binaryDirectory + "RegionalExpression " + stateOfTheWorld.expressionFiles[cases[0].disease()].FullName );
                hpcScript.Write(@"job add %1 /exclusive /numnodes:1-1 /scheduler:gcr " + stateOfTheWorld.configuration.hpcBinariesDirectory + "RegionalExpression " + stateOfTheWorld.expressionFiles[cases[0].disease()].FullName);
                foreach (var case_ in cases) {
                    script.Write(" " + case_.case_id);
                    hpcScript.Write(" " + case_.case_id);
                }
                script.WriteLine();
                hpcScript.WriteLine();
            }

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

        //
        // This represents the state of the world.  Processing stages look at this state and generate actions to move us along.
        //
        class StateOfTheWorld
        {
            public StateOfTheWorld(ASETools.ASEConfirguation configuration_) 
            {
                configuration = configuration_;
            }

            public ASETools.ASEConfirguation configuration;
            public Dictionary<string, ASETools.DownloadedFile> downloadedFiles = null;
            public Dictionary<string, List<ASETools.DerivedFile>> derivedFiles = null;
            public Dictionary<string, FileInfo> expressionFiles = null;            
            public Dictionary<string, ASETools.MAFInfo> mafInfo = null;
            public Dictionary<string, ASETools.Case> cases = null;
            public List<string> diseases = null;
            public Dictionary<string, string> fileIdToCaseId = null;
            public Dictionary<string, long> fileSizesFromGDC = null;

            public void DetermineTheStateOfTheWorld()
            {
                ASETools.ScanFilesystems(configuration, out downloadedFiles, out derivedFiles);

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
                        if (null != case_.methylation_file_id && "" != case_.methylation_file_id)
                        {
                            fileIdToCaseId.Add(case_.methylation_file_id, case_.case_id);
                        }
                    }


                    if (true)
                    {
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

                        if (downloadedFiles.ContainsKey(case_.methylation_file_id))
                        {
                            nMethylation++;
                            bytesMethylation += (ulong)downloadedFiles[case_.methylation_file_id].fileInfo.Length;
                        }

                        if (downloadedFiles.ContainsKey(case_.copy_number_file_id))
                        {
                            nCopyNumber++;
                            bytesCopyNumber += (ulong)downloadedFiles[case_.copy_number_file_id].fileInfo.Length;
                        }
                    } // foreach case

                    Console.WriteLine(nNormalDNA + "(" + ASETools.SizeToUnits(bytesNormalDNA) + "B) normal DNA, " + nTumorDNA + "(" + ASETools.SizeToUnits(bytesTumorDNA) + "B) tumor DNA, " +
                                      nNormalRNA + "(" + ASETools.SizeToUnits(bytesNormalRNA) + "B) normal RNA, " + nTumorRNA + "(" + ASETools.SizeToUnits(bytesTumorRNA) + "B) tumor RNA, " +
                                      nMethylation + "(" + ASETools.SizeToUnits(bytesMethylation) + "B) methylation, " + nCopyNumber + "(" + ASETools.SizeToUnits(bytesCopyNumber) + "B) copy number");

                } // If we loaded cases

                expressionFiles = new Dictionary<string, FileInfo>();

                if (Directory.Exists(configuration.expressionFilesDirectory))
                {
                    foreach (var filename in Directory.EnumerateFiles(configuration.expressionFilesDirectory + "expression_*"))
                    {
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

                    if (case_.methylation_file_id != "")
                    {
                        fileSizesFromGDC.Add(case_.methylation_file_id, case_.methylation_size);
                    }

                    if (case_.copy_number_file_id != "")
                    {
                        fileSizesFromGDC.Add(case_.copy_number_file_id, case_.copy_number_size);
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
        } // StateOfTheWorld


        static void Main(string[] args)
        {
            var stopwatch = new Stopwatch();
            stopwatch.Start();

            var configuration = ASETools.ASEConfirguation.loadFromFile(args);

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

            bool checkDependencies = configuration.commandLineArgs.Count() >= 1 && configuration.commandLineArgs.Contains("-d");

            var stateOfTheWorld = new StateOfTheWorld(configuration);
            stateOfTheWorld.DetermineTheStateOfTheWorld();
            
            Console.WriteLine();

            if (null != stateOfTheWorld.cases)
            {
                //
                // Rewrite the cases file, since we have just updated all of the file locations for downloaded and derived files.
                //
                ASETools.Case.SaveToFile(stateOfTheWorld.cases, configuration.casesFilePathname);
            }

            List<ProcessingStage> processingStages = new List<ProcessingStage>();

            processingStages.Add(new MAFConfigurationProcessingStage());
            processingStages.Add(new GenerateCasesProcessingStage());
            processingStages.Add(new AllcountProcesingStage());
            processingStages.Add(new DownloadProcessingStage());
            processingStages.Add(new MD5ComputationProcessingStage());
            processingStages.Add(new GermlineVariantCallingProcessingStage());
            processingStages.Add(new SelectVariantsProcessingStage());
            processingStages.Add(new ExpressionDistributionProcessingStage());
            processingStages.Add(new ExtractMAFLinesProcessingStage());

            if (checkDependencies)
            {
                bool allDependenciesOK = true;
                foreach (var processingStage in processingStages)
                {
                    allDependenciesOK &= processingStage.EvaluateDependencies(stateOfTheWorld);
                }

                if (!allDependenciesOK)
                {
                    Console.WriteLine("Not generating scripts because some dependencies have been violated.  Delete the stale generated files and rerun.");
                    return;
                }
            }

            var script = ASETools.CreateStreamWriterWithRetry("ASENextSteps.cmd");
            StreamWriter hpcScript;

            if (configuration.hpcScriptFilename == "")  // The empty string means not to generate an output.  We do this by making a Null stream.
            {
                hpcScript = new StreamWriter(Stream.Null);
            }
            else
            {
                hpcScript = ASETools.CreateStreamWriterWithRetry(configuration.hpcScriptFilename);
            }

            var linuxScript = ASETools.CreateStreamWriterWithRetry("ASENextStepsLinux");

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

                processingStage.EvaluateStage(stateOfTheWorld, script, hpcScript, linuxScript, out stageFilesToDownload, out nDone, out nAddedToScript, out nWaitingForPrerequisites);

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

                Console.WriteLine(String.Format("{0," + (stageNameHeader.Count() + paddingSize) + "}", processingStage.GetStageName()) + " " + String.Format("{0,6}", nDone) + " " + String.Format("{0,6}", nAddedToScript) + " " +
                    String.Format("{0,7}", nWaitingForPrerequisites) + " " + String.Format("{0,11}", nDownloadsRequested));
            } // foreach stage

            //
            // Now put downloads in their own script. They're separated out because they need to be run in one of the download
            // directories, and the user may want to split them across machines.
            //
            const string downloadScriptFilename = "ASEDownload.cmd";
            long bytesToDownload = 0;
            if (allFilesToDownload.Count() == 0)
            {
                File.Delete(downloadScriptFilename);
            } 
            else
            {
                var downloadScript = ASETools.CreateStreamWriterWithRetry(downloadScriptFilename);

                foreach (var file in allFilesToDownload)
                {
                    downloadScript.WriteLine(configuration.binaryDirectory + "gdc-client download --token-file " + configuration.accessTokenPathname + " " + file);
                    bytesToDownload += stateOfTheWorld.fileSizesFromGDC[file];
                }

                downloadScript.Close();
            }

            script.Close();
            hpcScript.Close();
            linuxScript.Close();

            Console.WriteLine();
            Console.WriteLine("Downloading " + ASETools.SizeToUnits((ulong)bytesToDownload) + "B in " + allFilesToDownload.Count() + " files.");
            Console.WriteLine("ASEProcessManager took " + ASETools.ElapsedTimeInSeconds(stopwatch) + " and finished at " + DateTime.Now.ToLocalTime().ToString());
        }
    }
}
