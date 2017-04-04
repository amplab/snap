using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using ASELib;
using System.IO;

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
        // A ProcessingStage takes some kind of 
        interface ProcessingStage
        {
            string GetStageName();
            void EvaluateStage(StateOfTheWorld stateOfTheWorld, StreamWriter script, StreamWriter hpcScript, StreamWriter linuxScript, out List<string> filesToDownload, out int nDone, out int nAddedToScript, out int nWaitingForPrerequisites);
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
            }
        }

        class RNAAllcountProcesingStage : ProcessingStage
        {
            public RNAAllcountProcesingStage() { }

            public string GetStageName()
            {
                return "Generate RNA Allcount files";
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
                    if (!stateOfTheWorld.downloadedFiles.ContainsKey(case_.tumor_rna_file_id))
                    {
                        if (null == filesToDownload)
                        {
                            filesToDownload = new List<string>();
                        }
                        filesToDownload.Add(case_.tumor_rna_file_id);
                    }
                    else
                    {
                        var downloadedFile = stateOfTheWorld.downloadedFiles[case_.tumor_rna_file_id];

                        if (stateOfTheWorld.fileDownloadedAndVerified(case_.tumor_rna_file_id, case_.tumor_rna_file_bam_md5))
                        {
                            nWaitingForPrerequisites++;
                        }
                        else if (stateOfTheWorld.containsDerivedFile(case_.case_id, case_.tumor_rna_file_id, ASETools.DerivedFile.Type.Allcount))
                        {
                            nDone++;
                        } 
                        else
                        {
                            nAddedToScript++;
                            string caseDirectory = ASETools.GetDirectoryFromPathname(stateOfTheWorld.downloadedFiles[case_.tumor_rna_file_id].fileInfo.FullName) + @"\..\..\" + stateOfTheWorld.configuration.derivedFilesDirectory + @"\" + case_.case_id + @"\";
                            script.WriteLine("md " + caseDirectory);
                            script.WriteLine(stateOfTheWorld.configuration.binaryDirectory + "CountReadsCovering " + stateOfTheWorld.configuration.indexDirectory + " -a " + stateOfTheWorld.downloadedFiles[case_.tumor_rna_file_id].fileInfo.FullName + " " +
                                caseDirectory + case_.tumor_rna_file_id + ".allcount.gz");

                            hpcScript.WriteLine(@"job add %1 /exclusive /numnodes:1-1 /scheduler:" + stateOfTheWorld.configuration.hpcScheduler + " " +
                                stateOfTheWorld.configuration.hpcBinariesDirectory + "MakeDirectoryAndCountReadsCovering.cmd " + caseDirectory + " " + stateOfTheWorld.configuration.hpcBinariesDirectory + " " +
                                stateOfTheWorld.configuration.hpcIndexDirectory + " " + stateOfTheWorld.downloadedFiles[case_.tumor_rna_file_id].fileInfo.FullName + " " + caseDirectory + case_.tumor_rna_file_id + ".allcount.gz");
                        }
                    }
                } // Foreach case
            }// EvaluateStage

        } // RNAAllcountProcessingStage

        class NormalDNAAndMethylationDownloadProcessingStage : ProcessingStage
        {
            public NormalDNAAndMethylationDownloadProcessingStage() { }

            public string GetStageName()
            {
                return "Normal DNA and and Methylation Download";
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
                    if (!stateOfTheWorld.downloadedFiles.ContainsKey(case_.normal_dna_file_id))
                    {
                        if (null == filesToDownload)
                        {
                            filesToDownload = new List<string>();
                        }
                        filesToDownload.Add(case_.normal_dna_file_id);
                    }


                    if (case_.methylation_file_id != null && case_.methylation_file_id != "" && !stateOfTheWorld.downloadedFiles.ContainsKey(case_.methylation_file_id))
                    {
                        if (null == filesToDownload)
                        {
                            filesToDownload = new List<string>();
                        }
                        filesToDownload.Add(case_.methylation_file_id);
                    }
                } // foreach case
            } // EvaluateStage
        } // NormalDNAAndMethylationDownloadProcessingStage

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
                    if (case_.methylation_file_id != null && case_.methylation_file_id != "")
                    {
                        // for some reason, these are all wrong, so ignore them for now.  HandleFile(stateOfTheWorld, script, case_.methylation_file_id, case_.tumor_rna_file_bam_md5, ref nDone, ref nAddedToScript, ref nWaitingForPrerequisites);
                    }
                }
            } // EvaluateStage

            void HandleFile(StateOfTheWorld stateOfTheWorld, StreamWriter script, StreamWriter hpcScript, string fileId, string expectedMD5, ref int nDone, ref int nAddedToScript, ref int nWaitingForPrerequisites)
            {
                if (!stateOfTheWorld.downloadedFiles.ContainsKey(fileId) || null == expectedMD5 || "" == expectedMD5)
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

                    if (!stateOfTheWorld.fileDownloadedAndVerified(case_.normal_dna_file_id, case_.normal_rna_file_bam_md5))
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
                    linuxScript.Write("cat ~/genomes/grch38-100k-regions | parallel -k -j `cat ~/ncores` \" freebayes --region {} --fasta-reference ~/genomes/grch38.fa " + 
                        ASETools.WindowsToLinuxPathname(stateOfTheWorld.downloadedFiles[case_.normal_dna_file_id].fileInfo.FullName) + 
                        " \" | ~/freebayes/vcflib/bin/vcffirstheader | ~/freebayes/vcflib/bin/vcfstreamsort -w 1000 | ~/freebayes/vcflib/bin/vcfuniq > " +
                        case_.normal_dna_file_id + ".vcf");
                    // and more stuff here.....


                }
            
            
            }
        }


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
            public Dictionary<string, ASETools.MAFInfo> mafInfo = null;
            public Dictionary<string, ASETools.Case> cases = null;

            public void DetermineTheStateOfTheWorld()
            {
                ASETools.ScanFilesystems(configuration, out downloadedFiles, out derivedFiles);

                mafInfo = ASETools.MAFInfo.LoadMAFManifest(configuration.mafManifestPathname);
                cases = ASETools.Case.LoadCases(configuration.casesFilePathname);

                ASETools.Case.loadAllFileLocations(cases, downloadedFiles, derivedFiles);
            }

            public bool fileDownloadedAndVerified(string file_id, string expectedMD5)
            {
                return downloadedFiles.ContainsKey(file_id) && (null == expectedMD5 || downloadedFiles[file_id].storedMD5 == expectedMD5);
            }

            public bool containsDerivedFile(string case_id, string derived_from_file_id, ASETools.DerivedFile.Type type)
            {
                return derivedFiles.ContainsKey(case_id) && derivedFiles[case_id].Where(x => x.derived_from_file_id == derived_from_file_id && x.type == type).Count() != 0;
            }
        }

 
        static void Main(string[] args)
        {
            var configuration = ASETools.ASEConfirguation.loadFromFile(args);

            if (null == configuration)
            {
                Console.WriteLine("Giving up because we were unable to load configuration.");
                return;
            }

            var stateOfTheWorld = new StateOfTheWorld(configuration);
            stateOfTheWorld.DetermineTheStateOfTheWorld();

            //
            // Rewrite the cases file, since we have just updated all of the file locations for downloaded and derived files.
            //
            ASETools.Case.SaveToFile(stateOfTheWorld.cases, configuration.casesFilePathname);

            List<ProcessingStage> processingStages = new List<ProcessingStage>();

            processingStages.Add(new MAFConfigurationProcessingStage());
            processingStages.Add(new GenerateCasesProcessingStage());
            processingStages.Add(new RNAAllcountProcesingStage());
            processingStages.Add(new NormalDNAAndMethylationDownloadProcessingStage());
            processingStages.Add(new MD5ComputationProcessingStage());
            processingStages.Add(new GermlineVariantCallingProcessingStage());

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

            var linuxScript = ASETools.CreateStreamWriterWithRetry("ASELinuxNextSteps");

            var allFilesToDownload = new List<string>();


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

                Console.WriteLine("Processing stage '" + processingStage.GetStageName() + "' had " + nDone + " done, added " + nAddedToScript + " to the script and had " + nWaitingForPrerequisites +
                    " blocked waiting for prerequisites and requested " + nDownloadsRequested + " files be downloaded.");
            }


            //
            // Now put downloads in their own script. They're separated out because they need to be run in one of the download
            // directories, and the user may want to split them across machines.
            //
            const string downloadScriptFilename = "ASEDownload.cmd";
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
                }

                downloadScript.Close();
            }

            script.Close();
            hpcScript.Close();
            linuxScript.Close();
 
        }
    }
}
