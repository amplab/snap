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
            void EvaluateStage(StateOfTheWorld stateOfTheWorld, StreamWriter script, out List<string> filesToDownload, out int nDone, out int nAddedToScript, out int nWaitingForPrerequisites);
        }

        class MAFConfigurationProcessingStage : ProcessingStage 
        {
            public MAFConfigurationProcessingStage() { }

            public string GetStageName()
            {
                return "Generate MAF Configuration";
            }

            public void EvaluateStage(StateOfTheWorld stateOfTheWorld, StreamWriter script, out List<string> filesToDownload, out int nDone, out int nAddedToScript, out int nWaitingForPrerequisites)
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

            public void EvaluateStage(StateOfTheWorld stateOfTheWorld, StreamWriter script, out List<string> filesToDownload, out int nDone, out int nAddedToScript, out int nWaitingForPrerequisites)
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
            public Dictionary<string, ASETools.MAFInfo> mafInfo = null;
            public Dictionary<string, ASETools.Case> cases = null;

            public void DetermineTheStateOfTheWorld()
            {
                downloadedFiles = ASETools.DownloadedFile.ScanFilesystems(configuration);

                mafInfo = ASETools.MAFInfo.LoadMAFManifest(configuration.mafManifestPathname);
                cases = ASETools.Case.LoadCases(configuration.casesFilePathname);
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

            List<ProcessingStage> processingStages = new List<ProcessingStage>();

            processingStages.Add(new MAFConfigurationProcessingStage());
            processingStages.Add(new GenerateCasesProcessingStage());

            var script = ASETools.CreateStreamWriterWithRetry("ASENextSteps.cmd");

            var allFilesToDownload = new List<string>();


            foreach (var processingStage in processingStages)
            {
                int nDone;
                int nAddedToScript;
                int nWaitingForPrerequisites;
                List<string> stageFilesToDownload;

                processingStage.EvaluateStage(stateOfTheWorld, script, out stageFilesToDownload, out nDone, out nAddedToScript, out nWaitingForPrerequisites);

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
 
        }
    }
}
