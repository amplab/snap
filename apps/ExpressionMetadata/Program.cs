using System;
using System.Collections;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using System.IO;
using System.Xml.Linq;
using System.Diagnostics;
using ExpressionLib;
using System.Threading;

//
// A bunch of aliases for string in order to make the code more readable.
//
using AnalysisID = System.String;
using ParticipantID = System.String;
using SampleID = System.String;
using Pathname = System.String;
using LibraryStrategy = System.String;
using Center = System.String;
using SampleType = System.String;
using TumorType = System.String;

namespace ExpressionMetadata
{
    class Program
    {
        static string baseDirectory = @"f:\temp\expression\";
        static string isoformDirectory = baseDirectory + @"isoform-files\";
        static List<string> hg18_likeReferences = new List<string>();

        public static void AddEntireTumorToMachine(Dictionary<ExpressionTools.AnalysisType, Pathname> tumorToMachineMapping, TumorType tumorType, Pathname machine) 
        {
            foreach (string libraryStrategy in ExpressionTools.LibraryStrategies)
            {
                var analysisType = new ExpressionTools.AnalysisType(tumorType, libraryStrategy, true);
                tumorToMachineMapping.Add(analysisType, machine);

                analysisType = new ExpressionTools.AnalysisType(tumorType, libraryStrategy, false);
                tumorToMachineMapping.Add(analysisType, machine);
            }
        }

        public class Machine
        {
            public Machine(string name_, int memoryInGB_, int diskInTB_, bool secondaryMachine_)
            {
                name = name_;
                memoryInGB = memoryInGB_;
                diskInTB = diskInTB_;
                secondaryMachine = secondaryMachine_;
            }

            public static void AddMachine(string name_, int memoryInGB_, int diskInTB_, bool secondaryMachine_ = false)
            {
                Machines.Add(name_, new Machine(name_, memoryInGB_, diskInTB_, secondaryMachine_));
            }

            public string name;
            public int memoryInGB;
            public int diskInTB;
            public bool secondaryMachine;
        }

        public static void InitializeMachines()
        {
            Machine.AddMachine("msr-srs-0", 8, 36, true);
            Machine.AddMachine("msr-srs-1", 8, 36);
            Machine.AddMachine("msr-srs-2", 8, 36);
            Machine.AddMachine("msr-srs-3", 8, 34, true);
            Machine.AddMachine("msr-srs-4", 8, 36);
            Machine.AddMachine("msr-srs-5", 16, 36);
            Machine.AddMachine("msr-srs-6", 16, 36);
            Machine.AddMachine("msr-srs-7", 48, 48);
            Machine.AddMachine("msr-srs-8", 48, 48);
            Machine.AddMachine("fds-k25-1", 48, 6);
            Machine.AddMachine("fds-k25-2", 48, 6);
            Machine.AddMachine("fds-k25-3", 48, 6);
            Machine.AddMachine("fds-k25-5", 48, 6);
            Machine.AddMachine("fds-k25-7", 48, 6);
            Machine.AddMachine("fds-k25-8", 48, 6);
            Machine.AddMachine("fds-k25-9", 48, 6);
            Machine.AddMachine("fds-k25-11", 48, 6, true);
            Machine.AddMachine("fds-k25-12", 48, 6);
            Machine.AddMachine("fds-k25-14", 48, 6, true);
            Machine.AddMachine("fds-k25-15", 48, 6);
            Machine.AddMachine("fds-k25-17", 48, 6);
            Machine.AddMachine("fds-k25-18", 48, 6);
            Machine.AddMachine("fds-k25-19", 48, 6);
            Machine.AddMachine("fds-k25-20", 48, 6);
            Machine.AddMachine("fds-k24-1", 12, 40, true);
            Machine.AddMachine("fds-k24-2", 12, 40, true);
            Machine.AddMachine("fds-k24-3", 12, 40, true);
            Machine.AddMachine("fds-k24-4", 12, 40, true);
            Machine.AddMachine("fds-k24-5", 12, 40, true);
            Machine.AddMachine("fds-k24-6", 12, 40, true);
            Machine.AddMachine("fds-k24-7", 12, 40, true);
            // There is no fds-k24-8.
            Machine.AddMachine("fds-k24-9", 12, 40, true);
            Machine.AddMachine("fds-k24-10", 12, 40, true);
            Machine.AddMachine("fds-k24-11", 12, 40, true);
            Machine.AddMachine("fds-k24-12", 12, 40, true);
            Machine.AddMachine("fds-k24-13", 12, 40, true);
            Machine.AddMachine("fds-k24-14", 12, 60, true);
            Machine.AddMachine("msr-genomics-0", 256, 240, true);
            Machine.AddMachine("msr-genomics-1", 256, 240, true);
        }

        public static Dictionary<string, Machine> Machines = new Dictionary<string, Machine>();

        public static Dictionary<ExpressionTools.AnalysisType, Pathname> GenerateTumorToMachineMapping()
        {
            var tumorToMachineMapping = new Dictionary<ExpressionTools.AnalysisType, Pathname>();

            AddEntireTumorToMachine(tumorToMachineMapping, "acc",  "fds-k25-2");
            AddEntireTumorToMachine(tumorToMachineMapping, "blca", "msr-srs-7");
            AddEntireTumorToMachine(tumorToMachineMapping, "chol", "msr-srs-6");
            AddEntireTumorToMachine(tumorToMachineMapping, "coad", "msr-srs-6");
            AddEntireTumorToMachine(tumorToMachineMapping, "dlbc", "msr-srs-7");
            AddEntireTumorToMachine(tumorToMachineMapping, "gbm",  "msr-srs-4");
            AddEntireTumorToMachine(tumorToMachineMapping, "hnsc", "msr-srs-0");
            AddEntireTumorToMachine(tumorToMachineMapping, "kich", "fds-k25-9");
            AddEntireTumorToMachine(tumorToMachineMapping, "kirc", "fds-k25-17");
            AddEntireTumorToMachine(tumorToMachineMapping, "kirp", "msr-srs-6");
            AddEntireTumorToMachine(tumorToMachineMapping, "laml", "msr-srs-8");
            AddEntireTumorToMachine(tumorToMachineMapping, "luad", "msr-srs-6");
            AddEntireTumorToMachine(tumorToMachineMapping, "meso", "fds-k25-15");
            AddEntireTumorToMachine(tumorToMachineMapping, "ov", "msr-srs-7");
            AddEntireTumorToMachine(tumorToMachineMapping, "paad", "fds-k25-5");
            AddEntireTumorToMachine(tumorToMachineMapping, "pcpg", "fds-k25-1");
            AddEntireTumorToMachine(tumorToMachineMapping, "prad", "msr-srs-2");
            AddEntireTumorToMachine(tumorToMachineMapping, "read", "fds-k25-7");
            AddEntireTumorToMachine(tumorToMachineMapping, "sarc", "msr-srs-4");
            AddEntireTumorToMachine(tumorToMachineMapping, "stad", "msr-srs-5");
            AddEntireTumorToMachine(tumorToMachineMapping, "tgct", "msr-srs-0");
            AddEntireTumorToMachine(tumorToMachineMapping, "thca", "msr-srs-4");
            AddEntireTumorToMachine(tumorToMachineMapping, "thym", "msr-srs-1");
            AddEntireTumorToMachine(tumorToMachineMapping, "ucec", "msr-srs-2");
            AddEntireTumorToMachine(tumorToMachineMapping, "ucs",  "msr-srs-4");
            AddEntireTumorToMachine(tumorToMachineMapping, "uvm",  "msr-srs-4");

            tumorToMachineMapping.Add(new ExpressionTools.AnalysisType("brca", "rna", true), "msr-srs-7");
            tumorToMachineMapping.Add(new ExpressionTools.AnalysisType("brca", "wgs", true), "msr-srs-7");
            tumorToMachineMapping.Add(new ExpressionTools.AnalysisType("brca", "wxs", true), "msr-srs-1");
            tumorToMachineMapping.Add(new ExpressionTools.AnalysisType("brca", "rna", false), "msr-srs-1");
            tumorToMachineMapping.Add(new ExpressionTools.AnalysisType("brca", "wgs", false), "msr-srs-1");
            tumorToMachineMapping.Add(new ExpressionTools.AnalysisType("brca", "wxs", false), "msr-srs-1");

            tumorToMachineMapping.Add(new ExpressionTools.AnalysisType("cesc", "rna", true), "msr-srs-7");
            tumorToMachineMapping.Add(new ExpressionTools.AnalysisType("cesc", "wgs", true), "msr-srs-7");
            tumorToMachineMapping.Add(new ExpressionTools.AnalysisType("cesc", "wxs", true), "msr-srs-7");
            tumorToMachineMapping.Add(new ExpressionTools.AnalysisType("cesc", "rna", false), "fds-k25-15");
            tumorToMachineMapping.Add(new ExpressionTools.AnalysisType("cesc", "wgs", false), "fds-k25-15");
            tumorToMachineMapping.Add(new ExpressionTools.AnalysisType("cesc", "wxs", false), "fds-k25-15");

            tumorToMachineMapping.Add(new ExpressionTools.AnalysisType("esca", "rna", true), "fds-k25-1");
            tumorToMachineMapping.Add(new ExpressionTools.AnalysisType("esca", "wgs", true), "fds-k25-1");
            tumorToMachineMapping.Add(new ExpressionTools.AnalysisType("esca", "wxs", true), "fds-k25-1");
            tumorToMachineMapping.Add(new ExpressionTools.AnalysisType("esca", "rna", false), "fds-k25-1");
            tumorToMachineMapping.Add(new ExpressionTools.AnalysisType("esca", "wgs", false), "fds-k25-1");
            tumorToMachineMapping.Add(new ExpressionTools.AnalysisType("esca", "wxs", false), "msr-srs-7");

            tumorToMachineMapping.Add(new ExpressionTools.AnalysisType("lgg", "rna", true), "fds-k25-14");
            tumorToMachineMapping.Add(new ExpressionTools.AnalysisType("lgg", "wgs", true), "fds-k25-14");
            tumorToMachineMapping.Add(new ExpressionTools.AnalysisType("lgg", "wxs", true), "fds-k25-14");
            tumorToMachineMapping.Add(new ExpressionTools.AnalysisType("lgg", "rna", false), "fds-k25-14");
            tumorToMachineMapping.Add(new ExpressionTools.AnalysisType("lgg", "wgs", false), "fds-k25-14");
            tumorToMachineMapping.Add(new ExpressionTools.AnalysisType("lgg", "wxs", false), "fds-k25-3");

            tumorToMachineMapping.Add(new ExpressionTools.AnalysisType("lihc", "rna", true), "msr-srs-2");
            tumorToMachineMapping.Add(new ExpressionTools.AnalysisType("lihc", "wgs", true), "msr-srs-2");
            tumorToMachineMapping.Add(new ExpressionTools.AnalysisType("lihc", "wxs", true), "msr-srs-2");
            tumorToMachineMapping.Add(new ExpressionTools.AnalysisType("lihc", "rna", false), "msr-srs-2");
            tumorToMachineMapping.Add(new ExpressionTools.AnalysisType("lihc", "wgs", false), "fds-k25-20");
            tumorToMachineMapping.Add(new ExpressionTools.AnalysisType("lihc", "wxs", false), "msr-srs-2");

            tumorToMachineMapping.Add(new ExpressionTools.AnalysisType("skcm", "rna", true), "fds-k25-8");
            tumorToMachineMapping.Add(new ExpressionTools.AnalysisType("skcm", "wgs", true), "fds-k25-8");
            tumorToMachineMapping.Add(new ExpressionTools.AnalysisType("skcm", "wxs", true), "fds-k25-8");
            tumorToMachineMapping.Add(new ExpressionTools.AnalysisType("skcm", "rna", false), "msr-srs-5");
            tumorToMachineMapping.Add(new ExpressionTools.AnalysisType("skcm", "wgs", false), "msr-srs-5");
            tumorToMachineMapping.Add(new ExpressionTools.AnalysisType("skcm", "wxs", false), "msr-srs-5");

            tumorToMachineMapping.Add(new ExpressionTools.AnalysisType("lusc", "rna", true), "fds-k25-18");
            tumorToMachineMapping.Add(new ExpressionTools.AnalysisType("lusc", "wgs", true), "fds-k25-18");
            tumorToMachineMapping.Add(new ExpressionTools.AnalysisType("lusc", "wxs", true), "fds-k25-18");
            tumorToMachineMapping.Add(new ExpressionTools.AnalysisType("lusc", "rna", false), "msr-srs-7");
            tumorToMachineMapping.Add(new ExpressionTools.AnalysisType("lusc", "wgs", false), "msr-srs-7");
            tumorToMachineMapping.Add(new ExpressionTools.AnalysisType("lusc", "wxs", false), "msr-srs-7");

             return tumorToMachineMapping;
        }


        static void GenerateHashCommandForFile(Pathname file, Dictionary<Pathname, StreamWriter> hashScripts)
        {
            //
            // The file is in \\ form.  Switch it to d: form.
            //

            string[] pathnameComponents = file.Split('\\');
            if (pathnameComponents.Count() < 4 || pathnameComponents[0] != "" || pathnameComponents[1] != "" || pathnameComponents[3].Count() != 2 || pathnameComponents[3].Substring(1) != "$")
            {
                Console.WriteLine("GenerateHashCommandForFile: expected \\ pathname, got " + file);
                return;
            }

            string machine = pathnameComponents[2];
            string localDirname = pathnameComponents[3].Substring(0, 1) + @":";
            string unqualifiedFilename = pathnameComponents[pathnameComponents.Count() - 1];

            for (int i = 4; i < pathnameComponents.Count() - 1; i++)
            {
                localDirname += @"\" + pathnameComponents[i];
            }

            if (!hashScripts.ContainsKey(machine))
            {
                hashScripts.Add(machine, new StreamWriter(baseDirectory + @"hashData-" + machine + ".cmd"));
            }

            hashScripts[machine].WriteLine(@"cd /d " + localDirname);
            hashScripts[machine].WriteLine(@"ComputeMD5 " + unqualifiedFilename + @"> " + unqualifiedFilename + ".md5");
        }

        public static Dictionary<AnalysisID, ExpressionTools.StoredBAM> LoadStoredBAMs(Dictionary<ExpressionTools.AnalysisType, Pathname> tumorToMachineMapping)
        {
            var storedBAMs = new Dictionary<AnalysisID, ExpressionTools.StoredBAM>();

            //
            // Directory structure is \\msr-srs-%n\d$\tcga\{rna,wgs,wxs}\{tumor, normal}\disease_abbr\analysis_id\*.{bam,bai}.  We need to call LoadStoredBAMsForDirectory on each
            // of the disease_abbr directories.  
            //
            var threads = new List<Thread>();
            foreach (var machine in Machines)
            {
                threads.Add(new Thread(() => ExpressionTools.LoadStoredBamsForMachine(@"\\" + machine.Value.name + @"\d$\tcga", storedBAMs)));
            }

            threads.Add(new Thread(() => ExpressionTools.LoadStoredBamsForMachine(@"\\msr-genomics-0\e$\tcga", storedBAMs)));
            threads.Add(new Thread(() => ExpressionTools.LoadStoredBamsForMachine(@"\\msr-genomics-1\e$\tcga", storedBAMs)));
            threads.Add(new Thread(() => ExpressionTools.LoadStoredBamsForMachine(@"\\bolosky\f$\tcga", storedBAMs)));

            threads.ForEach(t => t.Start());
            threads.ForEach(t => t.Join());

            var hashScripts = new Dictionary<Pathname, StreamWriter>();

            foreach (var entry in storedBAMs)
            {
                var storedBAM = entry.Value;

                if (storedBAM.bamHash == null && !storedBAM.bamHashBeingComputed)
                {
                    GenerateHashCommandForFile(storedBAM.bamInfo.FullName, hashScripts);
                }

                if (storedBAM.baiHash == null && !storedBAM.baiHashBeingComputed)
                {
                    GenerateHashCommandForFile(storedBAM.baiInfo.FullName, hashScripts);
                }
            }

            foreach (var entry in hashScripts)
            {
                entry.Value.Close();
            }

            return storedBAMs;
        }



        static void VerifyStoredBAMPaths(Dictionary<AnalysisID, ExpressionTools.StoredBAM> storedBAMs, Dictionary<AnalysisID, ExpressionTools.TCGARecord> tcgaRecords, Dictionary<ExpressionTools.AnalysisType, Pathname> tumorTypeToMachineMapping)
        {
            //
            // Make sure that all of the analyses are in the right place.  They should be in a path that looks like \\machine\d$\tcga\{disease}\{tumor, normal}\{RNA, WGS, WXS}\<analysis-id>\X.bam
            //

            const int diseaseComponent = 5;
            const int tumorNormalComponent = 6;
            const int libraryComponent = 7;

            int nRealigned = 0;
            int nChecksumsPresent = 0;
            int nChecksumsWrong = 0;

            foreach (var bamEntry in storedBAMs)
            {
                AnalysisID analysisID = bamEntry.Key;
                ExpressionTools.StoredBAM storedBAM = bamEntry.Value;
                string fullPath = storedBAM.bamInfo.FullName.ToLower();
                string[] pathComponents = fullPath.Split('\\');
                int pathLength = pathComponents.Count();    // Remember, this is one more than you think, because the first one is the empty string due to \\
                for (int i = 0; i < pathLength; i++)
                {
                    pathComponents[i] = pathComponents[i].ToLower();
                }

                if (!tcgaRecords.ContainsKey(analysisID)) {
                    Console.WriteLine("Stored BAM doesn't have accompanying TCGA record, pathname " + storedBAM.bamInfo.FullName);
                    continue;
                }

                var tcgaRecord = tcgaRecords[analysisID];

                if (storedBAM.bamHash != null && tcgaRecord.bamHash != null && tcgaRecord.bamHash != "")
                {
                    nChecksumsPresent++;
                    if (storedBAM.bamHash != tcgaRecord.bamHash)
                    {
                        Console.WriteLine("Checksum mismatch on file " + storedBAM.bamInfo.FullName + " (size " + storedBAM.bamInfo.Length + "), expected " + tcgaRecord.bamHash + " but got " + storedBAM.bamHash);
                        nChecksumsWrong++;
                    }
                }

                if (storedBAM.baiHash != null && storedBAM.baiHash != tcgaRecord.baiHash && tcgaRecord.bamHash != null && tcgaRecord.baiHash != "")
                {
                    Console.WriteLine("Checksum mismatch on file " + storedBAM.baiInfo.FullName + " (size " + storedBAM.bamInfo.Length + "), expected " + tcgaRecord.baiHash + " but got " + storedBAM.baiHash);
                    nChecksumsWrong++;
                }

                if (tcgaRecord.localRealign)
                {
                    nRealigned++;
                }

                if (pathLength < 8 || pathComponents[0] != "" || pathComponents[1] != "")
                {
                    Console.WriteLine("Unparsable directory path for stored BAM: " + fullPath);
                    continue;
                }

                if (pathLength != 10 || (pathComponents[3] != "d$" && pathComponents[3] != "e$" && pathComponents[3] != "f$") || pathComponents[4] != "tcga" || pathComponents[8] != analysisID ||
                    (pathComponents[libraryComponent] != "rna" && pathComponents[libraryComponent] != "wgs" && pathComponents[libraryComponent] != "wxs") || (pathComponents[tumorNormalComponent] != "tumor" && pathComponents[tumorNormalComponent] != "normal"))
                {
                    Console.WriteLine("Invalid path for analysis ID " + analysisID + " path " + fullPath);
                    continue;
                }

                //
                // Its shape is right.  Check to see if it's actually the right place.
                //
                if (tcgaRecord.library_strategy.ToLower() != pathComponents[libraryComponent] || tcgaRecord.tumorSample != (pathComponents[tumorNormalComponent] == "tumor") || tcgaRecord.disease_abbr != pathComponents[diseaseComponent])
                {
                    Console.WriteLine("Analysis " + analysisID + " is in the wrong directory.  It's at " + fullPath + " but it should be at " + tcgaRecord.disease_abbr + @"\" + (tcgaRecord.tumorSample ? "tumor" : "normal") + @"\" + tcgaRecord.library_strategy.ToLower());
                    continue;
                }

                string correctMachine = tumorTypeToMachineMapping[new ExpressionTools.AnalysisType(tcgaRecord)].ToLower();
                if (correctMachine != pathComponents[2] && (!Machines.ContainsKey(pathComponents[2]) || !Machines[pathComponents[2]].secondaryMachine)) 
                {
                    Console.WriteLine("BAM stored on wrong machine " + fullPath + " expected to be on " + correctMachine);
                    continue;
                }

                if (tcgaRecord.bamFileName.ToLower() != pathComponents[pathLength - 1].ToLower())
                {
                    Console.WriteLine("Unexpected BAM file " + fullPath + " different from expected file name " + tcgaRecord.bamFileName);
                    continue;
                }
            }

            Console.WriteLine("Verified " + nChecksumsPresent + " checksums on BAM files, of which " + nChecksumsWrong + " were wrong.");
        }


        public static string GenerateRealignmentFileString(string analysisID, ExpressionTools.TCGARecord tcgaRecord, string newRefassem)
        {
            string tumorOrNormal = tcgaRecord.tumorSample ? "tumor" : "normal";
            return
                analysisID +
                "*" + tcgaRecord.participant_id +
                "*" + analysisID + "-SNAP-realigned-" + tcgaRecord.analysis_id + "-" + tcgaRecord.disease_abbr + "-" + tumorOrNormal + ".bam" +
                "*" + tcgaRecord.lastModified.ToString() +
                "*" + tcgaRecord.uploadDate.ToString() +
                "*" + tcgaRecord.aliquot_id +
                "*" + newRefassem +
                "*" + tcgaRecord.disease_abbr +
                "*" + tcgaRecord.sampleType +
                "*" + tcgaRecord.library_strategy +
                "*" + tcgaRecord.center_name +
                "*" + tumorOrNormal +
                "*" + tcgaRecord.analysis_id;
        }

        public static void RebuildRealignmentAnalysesFromExperimentsFile(Dictionary<AnalysisID, ExpressionTools.TCGARecord> tcgaRecords)
        {
            StreamWriter realignFile = new StreamWriter(realignPathname);   // This overwrites the existing file (which probably got deleted or else you wouldn't be doing this).
            string [] experiments = File.ReadAllLines(baseDirectory + @"experiments.txt");

            foreach (var experiment in experiments) {
                string [] fields = experiment.Split('\t');
                if (fields.Count() == 0 || fields[0] == "disease_abbr")
                {
                    continue;
                }

                for (int whichField = 8; whichField <= 9; whichField++)
                {
                    string pathname = fields[whichField];
                    if (!pathname.Contains("-SNAP-realigned-"))
                    {
                        continue;
                    }

                    string[] pathComponents = pathname.Split('\\');
                    string filename = pathComponents[pathComponents.Count() - 1];
                    if (filename.Count() < 87)
                    {
                        Console.WriteLine("Bizarre filename in experiments file: " + pathname);
                        continue;
                    }

                    AnalysisID realignedAnalysisID = filename.Substring(0, 36);
                    AnalysisID realignedFromAnalysisID = filename.Substring(52, 36);

                    if (!tcgaRecords.ContainsKey(realignedFromAnalysisID))
                    {
                        Console.WriteLine("Couldn't find tcga record for source of realigned file " + pathname + " (realigned analysis ID: " + realignedFromAnalysisID + ")");
                        continue;
                    }

                    realignFile.WriteLine(GenerateRealignmentFileString(realignedAnalysisID, tcgaRecords[realignedFromAnalysisID], fields[1]));
                }
            }
            realignFile.Close();

        }
        public static string GenerateRealignmentFileString(ExpressionTools.TCGARecord tcgaRecord, string newRefassem)
        {
            string analysisID = Guid.NewGuid().ToString().ToLower();
            return GenerateRealignmentFileString(analysisID, tcgaRecord, newRefassem);
        }
        public static void GenerateRealignmentAnalyses(List<ExpressionTools.Experiment> experiments) 
        {
            StreamWriter realignedData = new StreamWriter(realignPathname, true); // This opens for append
            foreach (var experiment in experiments)
            {
                if (experiment.normalNeedsRealignment && !experiment.NormalDNAAnalysis.localRealign)
                {
                    realignedData.WriteLine(GenerateRealignmentFileString(experiment.NormalDNAAnalysis, experiment.TumorRNAAnalysis.refassemShortName));
                }

                if (experiment.tumorNeedsRealignment && !experiment.TumorDNAAnalysis.localRealign)
                {
                    realignedData.WriteLine(GenerateRealignmentFileString(experiment.TumorDNAAnalysis, experiment.TumorRNAAnalysis.refassemShortName));
                }
            }
            realignedData.Close();
        }

        public static Dictionary<SampleID, List<ExpressionTools.TCGARecord>> BuildAnalysesBySample(Dictionary<AnalysisID, ExpressionTools.TCGARecord> tcgaRecords)
        {
            var analysesBySample = new Dictionary<SampleID, List<ExpressionTools.TCGARecord>>();

            foreach (var tcgaEntry in tcgaRecords)
            {
                var entry = tcgaEntry.Value;

                if (!analysesBySample.ContainsKey(entry.aliquot_id))
                {
                    analysesBySample.Add(entry.aliquot_id, new List<ExpressionTools.TCGARecord>());
                }
                analysesBySample[entry.aliquot_id].Add(entry);
            }

            return analysesBySample;
        }




        public static Dictionary<string, string> GenerateCenterNameToDomainName()
        {
            //
            // This comes from the table at https://tcga-data.nci.nih.gov/datareports/codeTablesReport.htm?codeTable=center, then
            // modified by hand to correspond to the actual strings in the TCGA metadata and MAF files.
            //

            var mapping = new Dictionary<string, string>();

            mapping.Add("BI", "broad.mit.edu");
            mapping.Add("MSKCC", "mskcc.org");  // Memorial Sloan-Kettering
            mapping.Add("BCM", "hgsc.bcm.edu"); // Baylor college of medicine
            mapping.Add("WUSM", "genome.wustl.edu");
            mapping.Add("BCCAGSC", "bcgsc.ca"); // British Columbia Cancer Agency
            mapping.Add("UCSC", "ucsc.edu");

            return mapping;
        }


 
        public static Dictionary<string, string> InvertStringDictionary(Dictionary<string, string> input)
        {
            var output = new Dictionary<string, string>();

            foreach (var entry in input)
            {
                output.Add(entry.Value, entry.Key);
            }

            return output;
        }






        static void GenereateAdditionalTCGAMetadataGeneratingScript(Dictionary<AnalysisID, ExpressionTools.TCGARecord> tcgaRecords, Dictionary<AnalysisID, ExpressionTools.StoredBAM> storedBAMs, Dictionary<ExpressionTools.AnalysisType, Pathname> tumorToMachineMapping)
        {
            var globalScript = new StreamWriter(baseDirectory + @"extractTCGAAdditionalMetadata-all.cmd");
            globalScript.WriteLine(@"del \temp\analysisSummaries.txt");

            var writers = new Dictionary<Pathname, StreamWriter>();
            int nNeedingExtraction = 0;
            foreach (var entry in tcgaRecords)
            {
                var record = entry.Value;

                if (record.hasAdditionalData)
                {
                    continue;
                }

                if (record.localRealign)
                {
                    //
                    // We just copy these from their source.
                    //
                    if (record.realignSource.hasAdditionalData)
                    {
                        record.readSampleSize = record.realignSource.readSampleSize;
                        record.minReadLength = record.realignSource.minReadLength;
                        record.maxReadLength = record.realignSource.maxReadLength;
                        record.meanReadLength = record.realignSource.meanReadLength;
                        record.minGoodBases = record.realignSource.minGoodBases;
                        record.maxGoodBases = record.realignSource.maxGoodBases;
                        record.meanGoodBases = record.realignSource.meanGoodBases;
                        record.anyPaired = record.realignSource.anyPaired;
                        record.allPaired = record.realignSource.allPaired;
                        record.hasAdditionalData = true;
                    }

                    //
                    // In any case, we never generate them.
                    //
                    continue;
                }

                if (!storedBAMs.ContainsKey(record.analysis_id)) {
                    //
                    // Can't read it until we load it.
                    //
                    continue;
                }

                ExpressionTools.AnalysisType analysisType = new ExpressionTools.AnalysisType(record);
                Pathname machine = tumorToMachineMapping[analysisType];

                if (!writers.ContainsKey(machine))
                {
                    writers.Add(machine, new StreamWriter(baseDirectory + @"extractTCGAAdditionalMetadata-" + machine + ".cmd"));
                    writers[machine].WriteLine(@"del \temp\analysisSummaries.txt");
                }

                string line = @"SummarizeReads \sequence\indices\hg19-24 " + record.analysis_id + " " + storedBAMs[record.analysis_id].bamInfo.FullName + @" \temp\analysisSummaries.txt";
                writers[machine].WriteLine(line);
                globalScript.WriteLine(line);
                nNeedingExtraction++;
            }

            foreach (var entry in writers)
            {
                entry.Value.Close();
            }
            globalScript.Close();

            Console.WriteLine("" + nNeedingExtraction + " analyses require extraction of additional metadata (read length).");
        }



        public static void RecordDownloadAndRealign(ExpressionTools.TCGARecord record, ref int nRealigns, ref long nRealignBytes, ref int nDownloads, ref long nDownloadBytes)
        {
            if (record.storedBAM != null)
            {
                //
                // Got it, nothing to download or realign.
                //
                record.storedBAM.needed = true;
                if (record.realignSource != null && record.realignSource.storedBAM != null) {
                    record.realignSource.storedBAM.needed = true;
                }

                return;
            }

            if (record.localRealign)
            {
                nRealigns++;
                nRealignBytes += record.realignSource.totalFileSize;

                if (record.realignSource.storedBAM != null)
                {
                    record.realignSource.storedBAM.needed = true;
                }
                else
                {
                    nDownloads++;
                    nDownloadBytes += record.realignSource.totalFileSize;
                }
            }
        }

        public static List<ExpressionTools.Experiment> BuildExperiments(Dictionary<string, ExpressionTools.Participant> participants)
        {
            var experiments = new List<ExpressionTools.Experiment>();

            int nNormalDownloads = 0;
            long normalBytesToDownload = 0;
            int nTumorDownloads = 0;
            long tumorBytesToDownload = 0;
            int nRNADownloads = 0;
            long rnaBytesToDownload = 0;

            int nNormalRealigns = 0;
            int nTumorRealigns = 0;
            long nNormalRealignBytes = 0;
            long nTumorRealignBytes = 0;

            int nMutations = 0;

            //
            // We try to do an experiment for each participant.
            //
            foreach (var entry in participants)
            {
                var participant = entry.Value;
                string participantID = entry.Key;

                if (participant.mafs.Count() == 0) {
                    //
                    // No mutation calls for this person, just give up.
                    //
                    continue;
                }

                ExpressionTools.Experiment bestExperiment = null;
                int bestExperimentScore = -1;

                bool hg18Like = participant.mafs[0][0].NcbiBuild == "36";
                //
                // Iterate over all the possible RNA samples.  To be a candidate, we must have the same hg18-like status as the maf.
                // Beyond that, select the best normal/tumor possible.  At the end of the loop, see if what we've got is better than the
                // best experiment so far.
                //
                int bestScore = -1;
                ExpressionTools.TCGARecord bestAnalysis = null;
                foreach (var RNASampleEntry in participant.tumorSamples)
                {
                    var RNASample = RNASampleEntry.Value;

                    foreach (var RNAtcgaRecord in RNASample.RNA)
                    {
                        if (RNAtcgaRecord.refassemShortName == "unaligned")
                        {
                            continue;
                        }
                        string rnaIndex = ExpressionTools.RefassemToIndex(RNAtcgaRecord.refassemShortName);

                        var experiment = new ExpressionTools.Experiment();
                        experiment.participant = participant;
                        experiment.maf = participant.mafs[0];
                        experiment.disease_abbr = RNAtcgaRecord.disease_abbr;

                        if (hg18_likeReferences.Contains(RNAtcgaRecord.refassemShortName) != hg18Like)
                        {
                            continue;   // Never use the wrong reference class
                        }

                        experiment.TumorRNAAnalysis = RNAtcgaRecord;
                        //
                        // Next, find a normal sample.  Any will do, but prefer in order
                        // 0) Ones we've realigned
                        // 1) Ones aligned with the same reference as the RNA.
                        // 2) Ones we have downloaded
                        // 3) WGS Samples
                        //
                        // We do this by assigning points for each category and sorting.  We use 1 point for WGS, 2 for downloaded and 4 for right reference,
                        // so that it's just a straight select the most important one, go down on ties.
                        //
                        bestAnalysis = null;
                        bestScore = -1;
                        foreach (var normalSampleEntry in participant.normalSamples)
                        {
                            var normalSample = normalSampleEntry.Value;
                            foreach (var normalTCGARecord in normalSample.DNA)
                            {
                                int score = 0;
                                if (normalTCGARecord.localRealign)
                                {
                                    score += 8;
                                }
                                if (ExpressionTools.RefassemToIndex(normalTCGARecord.refassemShortName) == rnaIndex)
                                {
                                    score += 4;
                                }
                                if (normalTCGARecord.storedBAM != null)
                                {
                                    score += 2;
                                }
                                if (normalTCGARecord.isWGS)
                                {
                                    score += 1;
                                }

                                if (score > bestScore)
                                {
                                    bestScore = score;
                                    bestAnalysis = normalTCGARecord;
                                }
                            }
                        }

                        if (-1 == bestScore)
                        {
                            continue;
                        }

                        experiment.NormalDNAAnalysis = bestAnalysis;

                        //
                        // Find normal RNA.  Unlike tumor RNA and normal & tumor DNA, if we don't have this we're still willing to use this participant
                        //
                        bestAnalysis = null;
                        bestScore = -1;
                        foreach (var normalSampleEntry in participant.normalSamples)
                        {
                            var normalSample = normalSampleEntry.Value;
                            foreach (var normalTCGARecord in normalSample.RNA)
                            {
                                int score = 0;
                                if (ExpressionTools.RefassemToIndex(normalTCGARecord.refassemShortName) != rnaIndex)
                                {
                                    continue;   // We don't realign RNA, so if it doens't match, forget it.
                                }
                                score = 4;
                                if (normalTCGARecord.storedBAM != null)
                                {
                                    score += 2;
                                }
                                if (normalTCGARecord.isWGS)
                                {
                                    score += 1;
                                }

                                if (score > bestScore)
                                {
                                    bestScore = score;
                                    bestAnalysis = normalTCGARecord;
                                }
                            }
                        }

                        experiment.NormalRNAAnalysis = bestAnalysis;

                        //
                        // Finally, find Tumor DNA.  It must have an upload date within a week of the RNA sample.
                        // 1) It is aligned against the same reference as the tumor RNA.
                        // 2) WGS Samples
                        // 3) We have it downloaded.
                        //
                        bestAnalysis = null;
                        bestScore = -1;
                        foreach (var tumorSampleEntry in participant.tumorSamples)
                        {
                            var tumorSample = tumorSampleEntry.Value;

                            foreach (var tumorTCGARecord in tumorSample.DNA)
                            {
                                int score = 0;

                                if (ExpressionTools.RefassemToIndex(tumorTCGARecord.refassemShortName) == rnaIndex)
                                {
                                    score += 4;
                                }

                                if (tumorTCGARecord.isWGS)
                                {
                                    score += 2;
                                }

                                if (tumorTCGARecord.storedBAM != null)
                                {
                                    score += 1;
                                }

                                if (score > bestScore)
                                {
                                    bestScore = score;
                                    bestAnalysis = tumorTCGARecord;
                                }
                            } // for each tumor tcga record

                        } // for each tumor sample

                        if (null == bestAnalysis)
                        {
                            continue;
                        }

                        experiment.TumorDNAAnalysis = bestAnalysis;

                        int experimentScore = 0;

                        //
                        // Score this experiment.
                        // 1) Not needing to realign the tumor or normal samples
                        // 2) Not needing to download the samples
                        // 3) WGS
                        //
                        // 1 & 3 and be 0, 1 or 2, but 2 can be 0, 1, 2 or 3, so we use powers of four.
                        //

                        //
                        // Always realign the normals for now, FreeBayes seems to have problems with some of them as they are in TCGA.
                        //
                        experiment.normalNeedsRealignment = /*ExpressionTools.RefassemToIndex(experiment.NormalDNAAnalysis.refassemShortName) != rnaIndex*/ !experiment.NormalDNAAnalysis.localRealign ||experiment.NormalDNAAnalysis.localRealign && experiment.NormalDNAAnalysis.storedBAM == null;
                        experiment.tumorNeedsRealignment = ExpressionTools.RefassemToIndex(experiment.TumorDNAAnalysis.refassemShortName) != rnaIndex || experiment.TumorDNAAnalysis.localRealign && experiment.TumorDNAAnalysis.storedBAM == null;

                        if (experiment.NormalDNAAnalysis.refassemShortName == experiment.TumorRNAAnalysis.refassemShortName)
                        {
                            experimentScore += 16;
                        }

                        if (experiment.TumorDNAAnalysis.refassemShortName == experiment.TumorRNAAnalysis.refassemShortName)
                        {
                            experimentScore += 16;
                        }

                        if (experiment.NormalDNAAnalysis.storedBAM != null)
                        {
                            experimentScore += 4;
                        }

                        if (experiment.TumorRNAAnalysis.storedBAM != null)
                        {
                            experimentScore += 4;
                        }

                        if (experiment.TumorDNAAnalysis.storedBAM != null)
                        {
                            experimentScore += 4;
                        }

                        if (experiment.TumorDNAAnalysis.isWGS)
                        {
                            experimentScore += 1;
                        }

                        if (experiment.NormalDNAAnalysis.isWGS)
                        {
                            experimentScore += 1;
                        }

                        if (experimentScore > bestExperimentScore)
                        {
                            bestExperimentScore = experimentScore;
                            bestExperiment = experiment;
                        }
                    } // for each RNA tcga entry
                } // for each RNA tumor sample


                if (null == bestExperiment)
                {
                    continue;
                }

                experiments.Add(bestExperiment);

                RecordDownloadAndRealign(bestExperiment.NormalDNAAnalysis, ref nNormalRealigns, ref nNormalRealignBytes, ref nNormalDownloads, ref normalBytesToDownload);
                RecordDownloadAndRealign(bestExperiment.TumorDNAAnalysis, ref nTumorRealigns, ref nTumorRealignBytes, ref nTumorDownloads, ref tumorBytesToDownload);

                if (bestExperiment.TumorRNAAnalysis.storedBAM == null)
                {
                    nRNADownloads++;
                    rnaBytesToDownload += bestExperiment.TumorRNAAnalysis.totalFileSize;
                }
                else
                {
                    bestExperiment.TumorRNAAnalysis.storedBAM.needed = true;
                }

                if (bestExperiment.NormalRNAAnalysis != null && bestExperiment.NormalRNAAnalysis.storedBAM == null)
                {
                    nRNADownloads++;
                    rnaBytesToDownload += bestExperiment.NormalRNAAnalysis.totalFileSize;
                }

                nMutations += bestExperiment.maf.Count();

            }

            Console.WriteLine("Of " + participants.Count() + " participants we generated " + experiments.Count() + " experiments containing " + String.Format("{0:n0}",nMutations) + " mutatuions, for which we need to download " + (nNormalDownloads + nTumorDownloads + nRNADownloads) + 
                " analyses consiting of " + String.Format("{0:n0}", (normalBytesToDownload + tumorBytesToDownload + rnaBytesToDownload)) + " bytes and realign " + (nNormalRealigns + nTumorRealigns) + " analyses.");
            Console.WriteLine("Normal: " + nNormalDownloads + " downloads, " + String.Format("{0:n0}", normalBytesToDownload) + " bytes and " + nNormalRealigns + " realigns of " + String.Format("{0:n0}",nNormalRealignBytes) + " source bytes");
            Console.WriteLine("Tumor: " + nTumorDownloads + " downloads, " + String.Format("{0:n0}", tumorBytesToDownload) + " bytes and " + nTumorRealigns + " realigns of " + String.Format("{0:n0}",nTumorRealignBytes) + " source bytes");
            Console.WriteLine("RNA: " + nRNADownloads + " downloads, " + String.Format("{0:n0}", rnaBytesToDownload) + " bytes.");

            var countByDisease = new Dictionary<string, int>();
            var readyToGoByDisease = new Dictionary<string, int>();
            var readyToGoExceptForNormalByDisease = new Dictionary<string, int>();
            foreach (var experiment in experiments)
            {
                string disease = experiment.NormalDNAAnalysis.disease_abbr;
                if (!countByDisease.ContainsKey(disease)) {
                    countByDisease.Add(disease, 0);
                    readyToGoByDisease.Add(disease, 0);
                    readyToGoExceptForNormalByDisease.Add(disease, 0);
                }
                countByDisease[disease]++;
                if (experiment.NormalDNAAnalysis.storedBAM != null && experiment.TumorRNAAnalysis.storedBAM != null && experiment.TumorDNAAnalysis.storedBAM != null && !experiment.normalNeedsRealignment && !experiment.tumorNeedsRealignment) {
                    readyToGoByDisease[disease]++;
                }
                if (experiment.TumorRNAAnalysis.storedBAM != null && experiment.TumorDNAAnalysis.storedBAM != null && !experiment.tumorNeedsRealignment)
                {
                    readyToGoExceptForNormalByDisease[disease]++;
                }
            }

            foreach (var entry in countByDisease) {
                if (entry.Value != readyToGoByDisease[entry.Key])
                {
                    Console.WriteLine(entry.Key + ":" + entry.Value + " (" + readyToGoByDisease[entry.Key] + " ready to go, " + readyToGoExceptForNormalByDisease[entry.Key] + " if you don't count normal)");
                }
            }

            return experiments;
        }

        public static void AddDownloadToScript(ExpressionTools.TCGARecord tcgaRecord, Dictionary<string, StreamWriter> downloadScripts, Dictionary<ExpressionTools.AnalysisType, Pathname> tumorToMachineMapping, Dictionary<string, long> downloadAmounts)
        {
            string machine = tumorToMachineMapping[new ExpressionTools.AnalysisType(tcgaRecord)];
            if (!downloadScripts.ContainsKey(machine))
            {
                downloadScripts.Add(machine, new StreamWriter(baseDirectory + @"loadFromTCGA-" + machine + ".cmd"));
                downloadAmounts.Add(machine, 0);
            }

            downloadScripts[machine].WriteLine(@"cd /d d:\tcga\" + tcgaRecord.disease_abbr + (tcgaRecord.tumorSample ? @"\tumor\" : @"\normal\") + tcgaRecord.library_strategy);
            downloadScripts[machine].WriteLine(@"gtdownload -c c:\bolosky\ravip-cghub.key -k 30 " + tcgaRecord.analysis_id);

            downloadAmounts[machine] += tcgaRecord.totalFileSize;

            if (-1 != tcgaRecord.bamFileName.IndexOf("HOLD_QC"))
            {
                Console.WriteLine("Warning: generated download command for analysis " + tcgaRecord.analysis_id + " that has bam file name " + tcgaRecord.bamFileName + ".  Perhaps you should exclude it.");
            }
        }

        public static void GenerateDownloadScripts(List<ExpressionTools.Experiment> experiments, Dictionary<ExpressionTools.AnalysisType, Pathname> tumorToMachineMapping)
        {
            Console.WriteLine();
            var downloadScripts = new Dictionary<string, StreamWriter>();
            var downloadAmounts = new Dictionary<string, long>();
            foreach (var experiment in experiments)
            {
                if (experiment.TumorRNAAnalysis.storedBAM == null && !experiment.TumorRNAAnalysis.localRealign) {
                    AddDownloadToScript(experiment.TumorRNAAnalysis, downloadScripts, tumorToMachineMapping, downloadAmounts);
                }

                if (experiment.NormalRNAAnalysis != null && experiment.NormalRNAAnalysis.storedBAM == null && !experiment.NormalRNAAnalysis.localRealign)
                {
                    AddDownloadToScript(experiment.NormalRNAAnalysis, downloadScripts, tumorToMachineMapping, downloadAmounts);
                }

                if (experiment.TumorDNAAnalysis.storedBAM == null)
                {
                    if (experiment.TumorDNAAnalysis.localRealign)
                    {
                        if (experiment.TumorDNAAnalysis.realignSource.storedBAM == null)
                        {
                            AddDownloadToScript(experiment.TumorDNAAnalysis.realignSource, downloadScripts, tumorToMachineMapping, downloadAmounts);
                        }
                    }
                    else
                    {
                        AddDownloadToScript(experiment.TumorDNAAnalysis, downloadScripts, tumorToMachineMapping, downloadAmounts);
                    }
                }

                if (experiment.NormalDNAAnalysis.storedBAM == null)
                {
                    if (experiment.NormalDNAAnalysis.localRealign)
                    {
                        if (experiment.NormalDNAAnalysis.realignSource.storedBAM == null)
                        {
                            AddDownloadToScript(experiment.NormalDNAAnalysis.realignSource, downloadScripts, tumorToMachineMapping, downloadAmounts);
                        }
                    }
                    else
                    {
                        AddDownloadToScript(experiment.NormalDNAAnalysis, downloadScripts, tumorToMachineMapping, downloadAmounts);
                    }
                }
            }

            foreach (var entry in downloadScripts)
            {
                Console.WriteLine(entry.Key + " needs to download " + String.Format("{0:n0}", downloadAmounts[entry.Key]) + " bytes.");
                entry.Value.Close();
            }
        }


        static public void GenerateUnneededLists(Dictionary<AnalysisID, ExpressionTools.StoredBAM> storedBAMs)
        {
            var perMachineLists = new Dictionary<string, List<ExpressionTools.StoredBAM>>();
            var perMachineDeleteBytes = new Dictionary<string, long>();

            long totalNeededBytes = 0;
            int nNeededBAMs = 0;
            int nUnneededBAMs = 0;
            long totalUnneededBytes = 0;

            foreach (var bamEntry in storedBAMs) {
                var storedBAM = bamEntry.Value;

                long size = storedBAM.bamInfo.Length + storedBAM.baiInfo.Length;

                if (storedBAM.needed)
                {
                    nNeededBAMs++;
                    totalNeededBytes += size;
                }
                else
                {
                    nUnneededBAMs++;
                    totalUnneededBytes += size;

                    string machine = storedBAM.bamInfo.FullName.Split('\\')[2];
                    if (!perMachineDeleteBytes.ContainsKey(machine))
                    {
                        perMachineDeleteBytes.Add(machine, 0);
                        perMachineLists.Add(machine, new List<ExpressionTools.StoredBAM>());
                    }

                    perMachineDeleteBytes[machine] += size;
                    perMachineLists[machine].Add(storedBAM);
                }

            }

            Console.WriteLine("Need " + nNeededBAMs + " bams consisting of " + String.Format("{0:n0}", totalNeededBytes) + " bytes; do not need " + nUnneededBAMs + " consisting of " + String.Format("{0:n0}", totalUnneededBytes));

            var deleteLog = new StreamWriter(baseDirectory + @"deletable-files.txt");

            foreach (var machineEntry in perMachineDeleteBytes)
            {
                string machine = machineEntry.Key;
                Console.WriteLine(machine + " doesn't need " + String.Format("{0:n0}", perMachineDeleteBytes[machine]) + " bytes that it's storing.");
                foreach (var storedBAM in perMachineLists[machine])
                {
                    deleteLog.WriteLine("" + (storedBAM.bamInfo.Length + storedBAM.baiInfo.Length) + "\t" + machine + "\t" + storedBAM.bamInfo.FullName);
                }
            }
            deleteLog.Close();
        }


        public static void AddSingleRealignmentToScript(ExpressionTools.TCGARecord record, StreamWriter script, Dictionary<AnalysisID, ExpressionTools.TCGARecord> tcgaRecords, Dictionary<AnalysisID, ExpressionTools.StoredBAM> storedBAMs,
            Dictionary<ExpressionTools.AnalysisType, Pathname> tumorToMachineMapping, bool local, bool bigMem, bool cluster)
        {
            long lastChunkOfGuid = long.Parse(record.analysis_id.Substring(24), System.Globalization.NumberStyles.HexNumber);
            var analysisType = new ExpressionTools.AnalysisType(record);
            var realignSource = tcgaRecords[record.realignedFrom];
            if (realignSource.storedBAM == null || !realignSource.hasAdditionalData)
            {
                //
                // First we have to download the source and get its additional data (read length, etc.), then we can realign it.
                //
                return;
            }

            int seedLength;
            if (realignSource.meanGoodBases / 3 < 13)
            {
                seedLength = 13;
            }
            else if (realignSource.meanGoodBases / 3 > 32)
            {
                seedLength = 32;
            }
            else
            {
                seedLength = (int)(realignSource.meanGoodBases / 3 + .99);
            }

            string tumorOrNormal = (record.tumorSample ? "tumor" : "normal");
            string outputFileName = record.analysis_id + "-SNAP-realigned-" + realignSource.analysis_id + "-" + record.disease_abbr + "-" + tumorOrNormal + ".bam";
            string destinationDirectory;

            if (!cluster && (local || bigMem))// bigMem machines are the msr-genomics ones, and they have tons of storage.  Just hold them here.
            {
                destinationDirectory = @"d:\tcga\" + record.disease_abbr + @"\" + tumorOrNormal + @"\" + record.library_strategy + @"\" + record.analysis_id + @"\";
            } else {
                //
                // Just write directly to msr-genomics-0 or 1, d or e.
                //
                bool useGenomics0 = lastChunkOfGuid % 2 == 0;
                bool useD = (lastChunkOfGuid / 2) % 3 > 0;  // Favor d over e 2 to 1.

                if (cluster && (lastChunkOfGuid / 1000 ) % 4 >= 1 && false)  // 0% of the time write to the non-genomics location
                {
                    destinationDirectory = @"\\" + tumorToMachineMapping[analysisType] + @"\d$\tcga\" + record.disease_abbr + @"\" + tumorOrNormal + @"\" + record.library_strategy + @"\" + record.analysis_id + @"\";
                }
                else
                {
                    destinationDirectory = @"\\" + (useGenomics0 ? "msr-genomics-0" : "msr-genomics-1") + @"\" + (useD ? "d" : "e") + @"$\tcga\" + record.disease_abbr + @"\" + tumorOrNormal + @"\" + record.library_strategy + @"\" + record.analysis_id + @"\";
                }
            }

            if (cluster || local || bigMem) 
            {
                outputFileName = destinationDirectory + outputFileName;
            }

            if (local || bigMem || cluster)
            {
                script.WriteLine(@"md " + destinationDirectory);
            }

            if (cluster)
            {
                script.Write(@"job add %1 /exclusive /numnodes:1-1 /scheduler:gcr ");
            }

            ExpressionTools.UseAnalysis(record.realignedFrom, null, storedBAMs);

            string indexDrive = "";

            if (bigMem)
            {
                indexDrive = @"d:";
            }
            else if (cluster)
            {
                bool useGenomics0 = lastChunkOfGuid % 2 == 0;
                indexDrive = @"\\msr-genomics-" + (useGenomics0 ? "0" : "1") + @"\d$";
                script.Write(@"\\gcr\scratch\b99\bolosky\");
            }

            script.Write("snap ");
            if (realignSource.allPaired)
            {
                script.Write(@"paired " + indexDrive + @"\sequence\indices\");
            }
            else
            {
                script.Write(@"single " + indexDrive + @"\sequence\indices\");
            }

            script.Write(record.refassemShortName + "-");

            script.Write("" + seedLength + " " + realignSource.storedBAM.bamInfo.FullName + " -o " + outputFileName);
            if (realignSource.allPaired)
            {
                script.Write(" -s 0 1000");
            }
            script.Write(" -di -so -sm " + ((bigMem || cluster) ? "60" : "5 -kts"));
            if (bigMem)
            {
                script.Write(" -sid .");
            }

            if (cluster)
            {
                script.Write(@" -sid d:\scratch\bolosky\ -pre");
            }

            script.WriteLine(" -lp -map -pc -pro -mrl " + Math.Min(50, seedLength * 2));

            if (bigMem)
            {
                script.WriteLine(@"del .\*.tmp");
            }
            else if (!cluster)
            {
                script.WriteLine("del " + outputFileName + ".tmp"); // In case the alignment failed for whatever reason, don't copy the intermediate file, just get rid of it
            }

            if (!local && !bigMem && !cluster)
            {
                script.WriteLine(@"md " + destinationDirectory);
                script.WriteLine("copy " + outputFileName + "* " + destinationDirectory);
                script.WriteLine(@"del " + outputFileName + "*");
            }
        }

        public static void GenerateRealignmentScript(
            List<ExpressionTools.Experiment>                    experiments, 
            Dictionary<AnalysisID, ExpressionTools.TCGARecord>  tcgaRecords, 
            Dictionary<AnalysisID, ExpressionTools.StoredBAM>   storedBAMs, 
            Dictionary<ExpressionTools.AnalysisType, Pathname>  tumorToMachineMapping)
        {
            //
            // Find any experiments that rely on a realignment, and generate the SNAP command to do the realignment.
            //
            StreamWriter realignNormalScript = null;
            StreamWriter realignTumorScript = null;
            StreamWriter realignBigMemScript = null;
            StreamWriter realignClusterScript = null;
            var perMachineScripts = new Dictionary<Pathname, StreamWriter>();

            var perMachineOutput = new Dictionary<Pathname, long>();

            const long limitForBigMem = (long)50 * 1024 * 1024 * 1024;
            const long limitFor48 = (long)25 * 1024 * 1024 * 1024;

            foreach (var experiment in experiments)
            {
                List<ExpressionTools.TCGARecord> analyses = new List<ExpressionTools.TCGARecord>();

                if (experiment.NormalDNAAnalysis.localRealign && experiment.NormalDNAAnalysis.storedBAM == null)
                {
                    analyses.Add(experiment.NormalDNAAnalysis);
                }
                if (experiment.TumorDNAAnalysis.localRealign && experiment.TumorDNAAnalysis.storedBAM == null)
                {
                    analyses.Add(experiment.TumorDNAAnalysis);
                }

                foreach (var analysis in analyses)
                {
                    if (analysis.realignSource.storedBAM == null)
                    {
                        //
                        // Waiting for the input data.
                        //
                        continue;
                    }

                    StreamWriter script;
                    bool local;
                    bool bigMem;

                    if (analysis.realignSource.storedBAM.totalSize >= limitForBigMem)
                    {
                        if (realignBigMemScript == null)
                        {
                            realignBigMemScript = new StreamWriter(baseDirectory + @"realignBigMem.cmd");
                        }
                        script = realignBigMemScript;
                        local = false;
                        bigMem = true;
                    }
                    else
                    {
                        var analysisType = new ExpressionTools.AnalysisType(analysis);
                        int localMem = Machines[tumorToMachineMapping[analysisType]].memoryInGB;

                        bigMem = false;

                        if (localMem < 48 || localMem == 48 && analysis.realignSource.storedBAM.totalSize >= limitFor48 || true /* don't generate local anymore*/)
                        {
                            local = false;
                            if (analysis.tumorSample)
                            {
                                if (realignTumorScript == null)
                                {
                                    realignTumorScript = new StreamWriter(baseDirectory + @"realignTumor.cmd");
                                }
                                script = realignTumorScript;
                            }
                            else
                            {
                                if (realignNormalScript == null)
                                {
                                    realignNormalScript = new StreamWriter(baseDirectory + @"realignNormal.cmd");
                                }
                                script = realignNormalScript;
                            }
                        }
                        else
                        {
                            if (!perMachineScripts.ContainsKey(tumorToMachineMapping[analysisType]))
                            {
                                perMachineScripts.Add(tumorToMachineMapping[analysisType], new StreamWriter(@"f:\temp\expression\realign-" + tumorToMachineMapping[analysisType] + ".cmd"));
                            }

                            script = perMachineScripts[tumorToMachineMapping[analysisType]];
                            local = true;
                        }
                    }

                    Pathname destMachine = tumorToMachineMapping[new ExpressionTools.AnalysisType(analysis)];
                    if (!perMachineOutput.ContainsKey(destMachine))
                    {
                        perMachineOutput.Add(destMachine, 0);
                    }
                    perMachineOutput[destMachine] += analysis.realignSource.totalFileSize;
                    AddSingleRealignmentToScript(analysis, script, tcgaRecords, storedBAMs, tumorToMachineMapping, local, bigMem, false);
                    if (realignClusterScript == null)
                    {
                        var createJobScript = new StreamWriter(baseDirectory + "createjob.cmd");
                        createJobScript.WriteLine(@"job new /emailaddress:bolosky@microsoft.com /nodegroup:B99,ExpressQ /exclusive:true /failontaskfailure:false /jobname:align12 /memorypernode:100000 /notifyoncompletion:true /numnodes:1-10 /runtime:12:00 /scheduler:gcr /jobtemplate:ExpressQ /estimatedprocessmemory:100000");
                        createJobScript.Close();
                        realignClusterScript = new StreamWriter(baseDirectory + "realignCluster.cmd");
                    }
                    AddSingleRealignmentToScript(analysis, realignClusterScript, tcgaRecords, storedBAMs, tumorToMachineMapping, false, false, true);
                }
             }

            if (realignNormalScript != null)
            {
                realignNormalScript.Close();
            }

            if (realignTumorScript != null)
            {
                realignTumorScript.Close();
            }

            if (realignBigMemScript != null)
            {
                realignBigMemScript.Close();
            }

            foreach (var entry in perMachineScripts)
            {
                entry.Value.Close();
            }

            if (realignClusterScript != null)
            {
                //realignClusterScript.WriteLine("job submit /id:%1 /scheduler:gcr");
                realignClusterScript.Close();
            }

            foreach (var entry in perMachineOutput)
            {
                Console.WriteLine("Realign for " + entry.Key + " will consume about " + String.Format("{0:n0}", entry.Value) + " bytes.");
            }
        }

        public class DiseaseReferencePair
        {
            public DiseaseReferencePair(string disease_, string reference_)
            {
                disease = disease_;
                reference = reference_;
            }
            public string disease;
            public string reference;
        }

        static List<string> refassemWithoutChr = new List<string>();
        static Dictionary<string, string> refassemChrStateSwitched = new Dictionary<string, string>();
        static Dictionary<string, string> refassemTranslation = new Dictionary<string, string>();        // refassems in the input that are equivalent to ones we actually have


        static public void GenerateExtractionScripts(List<ExpressionTools.Experiment> experiments)
        {
            Console.WriteLine("Generating extraction scripts for " + experiments.Count() + " experiments");
            var stopwatch = Stopwatch.StartNew();
            var overallStopwatch = Stopwatch.StartNew();
            var extractRNAInputs = new Dictionary<string, StreamWriter>();  // Maps reference->input file
            var extractDNAInputs = new Dictionary<string, StreamWriter>();  // Maps reference->input file
            StreamWriter extractRNAScript = new StreamWriter(baseDirectory + @"extractRNAScript.cmd");
            StreamWriter extractTumorDNAScript = new StreamWriter(baseDirectory + @"extractDNASamples.cmd");
            var perMachineExtractScripts = new Dictionary<Pathname, StreamWriter>();

            var experimentSets = new Dictionary<DiseaseReferencePair, List<ExpressionTools.MAFRecord>>();
            int nExperimentsProcessed = 0;
            int nExperimentsWithVCFProcessed = 0;
            foreach (var experiment in experiments)
            {
                if (stopwatch.ElapsedMilliseconds >= 60000)
                {
                    Console.WriteLine("" + overallStopwatch.ElapsedMilliseconds / 1000 + "s: Processed " + nExperimentsProcessed + "/" + experiments.Count() + " experiments, of which " + nExperimentsWithVCFProcessed + " had VCFs");
                    stopwatch.Restart();
                }
                nExperimentsProcessed++;

                var diseaseAndReference = new DiseaseReferencePair(experiment.NormalDNAAnalysis.disease_abbr, experiment.TumorRNAAnalysis.refassemShortName);
                if (!experimentSets.ContainsKey(diseaseAndReference))
                {
                    experimentSets.Add(diseaseAndReference, new List<ExpressionTools.MAFRecord>());
                }

                string refassem = experiment.TumorDNAAnalysis.refassemShortName.ToLower();

                if (!extractDNAInputs.ContainsKey(refassem))
                {
                    extractDNAInputs.Add(refassem, new StreamWriter(baseDirectory + @"tumorDNAInput-" + refassem + ".txt"));
                    extractRNAInputs.Add(refassem, new StreamWriter(baseDirectory + @"tumorRNAInput-" + refassem + ".txt"));
                }

                ExpressionTools.VCF vcf = null;
                if (experiment.NormalDNAAnalysis.storedBAM != null && experiment.NormalDNAAnalysis.storedBAM.vcfInfo != null && false)
                {
                    vcf = new ExpressionTools.VCF(experiment.NormalDNAAnalysis.storedBAM.vcfInfo.FullName);
                    nExperimentsWithVCFProcessed++;
                }

                int maxDistanceBack = 200;
                int maxDistanceForward = 10;


                foreach (var mafRecord in experiment.maf)
                {
                    if (experiment.TumorDNAAnalysis.storedBAM != null)
                    {
                        //
                        // Write out all the variant lines that are close to this.  Note that this works correctly if vcf is null.
                        //
                        var enumerator = new ExpressionTools.VCFEnumerator(vcf, mafRecord.Chrom, mafRecord.Start_position - maxDistanceBack, mafRecord.Start_position + maxDistanceForward);
                        foreach (var variant in enumerator)
                        {
                            extractDNAInputs[refassem].WriteLine("variant\t" + variant.chromosome + "\t" + variant.pos + "\t" + variant.reference + "\t" + variant.alt + "\t" + variant.qual);
                        }

                        string tumorDNAFileName = diseaseAndReference.disease + @"\" + experiment.TumorDNAAnalysis.analysis_id + "-" +
                            mafRecord.Hugo_symbol + "-chr" + mafRecord.Chrom + "-" + Math.Max(1, mafRecord.Start_position - 200) + "-" +
                           (mafRecord.Start_position + 10);

                        extractDNAInputs[refassem].WriteLine(tumorDNAFileName + "\t" + mafRecord.entire_maf_line);
                        string line = "samtools view " + experiment.TumorDNAAnalysis.storedBAM.bamInfo.FullName +
                            " " + ExpressionTools.ChromPrefixFromRefassemChromosomeAndBam(experiment.TumorDNAAnalysis.storedBAM, mafRecord.Chrom, refassem) + mafRecord.Chrom.ToUpper() + ":" + Math.Max(1, mafRecord.Start_position - 200) + "-" + (mafRecord.Start_position + 10) + @" > " +
                                                             tumorDNAFileName;
                        extractTumorDNAScript.WriteLine(line);
                    }

                    if (experiment.TumorRNAAnalysis.storedBAM != null)
                    {
                        //
                        // Write out all the variant lines that are close to this.  Note that this works correctly if vcf is null.
                        //
                        var enumerator = new ExpressionTools.VCFEnumerator(vcf, mafRecord.Chrom, mafRecord.Start_position - maxDistanceBack, mafRecord.Start_position + maxDistanceForward);
                        foreach (var variant in enumerator)
                        {
                            extractRNAInputs[refassem].WriteLine("variant\t" + variant.chromosome + "\t" + variant.pos + "\t" + variant.reference + "\t" + variant.alt + "\t" + variant.qual);
                        }

                        string tumorRNAFileName = diseaseAndReference.disease + @"\" + experiment.TumorRNAAnalysis.analysis_id + "-" +
                            mafRecord.Hugo_symbol + "-chr" + mafRecord.Chrom + "-" + Math.Max(1, mafRecord.Start_position - 200) + "-" +
                           (mafRecord.Start_position + 10) + "-RNA";

                        extractRNAInputs[refassem].WriteLine(tumorRNAFileName + "\t" + mafRecord.entire_maf_line);
                        string line = "samtools view " + experiment.TumorRNAAnalysis.storedBAM.bamInfo.FullName +
                            " " + ExpressionTools.ChromPrefixFromRefassemChromosomeAndBam(experiment.TumorRNAAnalysis.storedBAM, mafRecord.Chrom, refassem) + mafRecord.Chrom.ToUpper() + ":" + Math.Max(1, mafRecord.Start_position - 200) + "-" + (mafRecord.Start_position + 10) + @" > " +
                                                             tumorRNAFileName;
                        extractRNAScript.WriteLine(line);
                    }
                }
            }

            Console.WriteLine("Processed " + nExperimentsProcessed + " experiments with " + nExperimentsWithVCFProcessed + " VCFs");
            extractRNAScript.Close();
            extractTumorDNAScript.Close();

            foreach (var entry in extractDNAInputs)
            {
                extractDNAInputs[entry.Key].Close();
                extractRNAInputs[entry.Key].Close();
            }

            foreach (var entry in perMachineExtractScripts)
            {
                entry.Value.Close();
            }


        }

        public static void GenerateLAMLMoves(Dictionary<AnalysisID, ExpressionTools.TCGARecord> tcgaRecords)
        {
            StreamWriter script = new StreamWriter(baseDirectory + @"moveLAML.cmd");

            foreach (var entry in tcgaRecords)
            {
                var record = entry.Value;
                if (record.disease_abbr != "laml")
                {
                    continue;
                }

                string destination;
                if (record.tumorSample) {
                    destination = @"tumor\";
                } else {
                    destination = @"normal\";
                }
                destination += record.library_strategy + @"\";
                script.WriteLine(@"mv " + record.analysis_id + " " + destination + record.analysis_id);
            }

            script.Close();

        }




        public static void AddAnalysisToMakeIsoformScript(ExpressionTools.TCGARecord analysis, Dictionary<Pathname, StreamWriter> scripts, string disease_abbr)
        {
            if (analysis.storedBAM != null && analysis.storedBAM.isoformCount == null)
            {
                if (!scripts.ContainsKey(analysis.storedBAM.machineName))
                {
                    scripts.Add(analysis.storedBAM.machineName, new StreamWriter(baseDirectory + @"CountIsoforms-" + analysis.storedBAM.machineName + ".cmd"));
                }

                scripts[analysis.storedBAM.machineName].WriteLine(@"cd /d " + analysis.storedBAM.driveLetter + @":\tcga\" + disease_abbr + @"\tumor\" + analysis.library_strategy + @"\" + analysis.analysis_id);
                scripts[analysis.storedBAM.machineName].WriteLine(@"d:\tcga\CountReadsCovering d:\sequence\indices\" + analysis.refassemShortName + @"-24 d:\sequence\gene_data\knownGene-hg" +
                    (hg18_likeReferences.Contains(analysis.refassemShortName) ? "18" : "19") + ".txt " + analysis.storedBAM.bamInfo.FullName + @" " + analysis.analysis_id + "-isoforms.txt");
            }

        }

        public static void GenerateMakeIsoformsScripts(List<ExpressionTools.Experiment> experiments)
        {
            var scripts = new Dictionary<Pathname, StreamWriter>();

            foreach (var experiment in experiments) {
                AddAnalysisToMakeIsoformScript(experiment.TumorRNAAnalysis, scripts, experiment.disease_abbr);
                AddAnalysisToMakeIsoformScript(experiment.TumorDNAAnalysis, scripts, experiment.disease_abbr);
            }

            foreach (var entry in scripts) {
                entry.Value.Close();
            }
        }


        public static int MedianOfSortedList(List<int> list) {
            int n = list.Count();
            if (n == 0) return 0;

            if (n % 2 == 1)
            {
                //
                // For an odd length list, just take the middle element.
                //
                return list[n / 2]; // Imagine a + 1 because of the roundoff and a -1 because of being 0-based
            }

            //
            // For an even length list, take the mean of the two elements in the center (rounding up, because that's the convention for .5)
            //
            return (list[n / 2] + list[n / 2 - 1] + 1) / 2;
        }

        public static Dictionary<string, int> GenerateSizeOfProteinCodingRegionMap()
        {
            //
            // This maps hugo symbol->count of bases.  Of course, most genes have multiple isoforms.  This just uses the one that's
            // listed as the ucsc symbol in the genenames.org database.
            //

            var ucscToSize = new Dictionary<string, int>();

            string[] ucscLines = File.ReadAllLines(@"f:\sequence\gene_data\knownGene-hg19.txt");
            for (int i = 1; i < ucscLines.Count(); i++ ) // start at 1 to skip the header line
            {
                var line = ucscLines[i];
                string[] fields = line.Split('\t');
                string[] exonStarts = fields[8].Split(',');
                string[] exonEnds = fields[9].Split(',');

                if (exonStarts.Count() != exonEnds.Count())
                {
                    Console.WriteLine("Non-matching exon start and end count on line " + (i + 1));
                    continue;
                }

                int totalSize = 0;
                for (int exonNumber = 0; exonNumber < exonStarts.Count() - 1; exonNumber++) // stop at Count-1 because these have trailing commas, and so the last field will always be empty
                {
                    int exonStart = Convert.ToInt32(exonStarts[exonNumber]);
                    int exonEnd = Convert.ToInt32(exonEnds[exonNumber]);

                    if (exonStart >= exonEnd)
                    {
                        Console.WriteLine("Zero or negative sized exon, line " + (i + 1) + " exon number " + (exonNumber + 1));
                        continue;
                    }

                    totalSize += exonEnd - exonStart;
                }

                ucscToSize.Add(fields[0], totalSize);
            }

            var geneToSize = new Dictionary<string, int>();

            string[] hgnc_lines = File.ReadAllLines(@"f:\sequence\gene_data\hgnc_complete_set.txt");

            for (int i = 1; i < hgnc_lines.Count(); i++) // start at 1 to skip header line
            {
                string[] fields = hgnc_lines[i].Split('\t');

                if (ucscToSize.ContainsKey(fields[21]))
                {
                    geneToSize.Add(fields[1], ucscToSize[fields[21]]);
                }
            }


                return geneToSize;
        }

        public static void GenerateMutationDistributionByGene(List<ExpressionTools.Experiment> experiments)
        {
            //
            // This determines the distribution of mutation counts based on whether each gene has 0, 1 or more than 1 mutation in each experiment.  It also normalizes them
            // by the size of the protein coding region.
            //
            var singleMutations = new Dictionary<string, List<int>>();

            //
            // First, determine all of the genes by running through all of the mafs and adding them to the singleMutations list.
            //
            foreach (var experiment in experiments)
            {
                foreach (var mafRecord in experiment.maf)
                {
                    if (!singleMutations.ContainsKey(mafRecord.Hugo_symbol))
                    {
                        singleMutations.Add(mafRecord.Hugo_symbol, new List<int>());
                    }
                }
            }

            //
            // Now make empty lists for the no mutations and multipleMutations groups.
            //
            var multipleMutations = new Dictionary<string, List<int>>();
            var noMutations = new Dictionary<string, List<int>>();
            foreach (var entry in singleMutations) {
                multipleMutations.Add(entry.Key, new List<int>());
                noMutations.Add(entry.Key, new List<int>());
            }

            //
            // Now go through each experiment, and decide how many mutations it has for each gene, and add it to the appropriate list.
            //
            foreach (var experiment in experiments) {
                int nMutations = experiment.maf.Count();
                foreach (var entry in singleMutations) {
                    string gene = entry.Key;

                    int count = experiment.maf.Where(mafRecord => mafRecord.Hugo_symbol == gene).Count();
                    if (count == 0) {
                        noMutations[gene].Add(nMutations);
                    } else if (count == 1) {
                        singleMutations[gene].Add(nMutations);
                    } else {
                        multipleMutations[gene].Add(nMutations);
                    }
                }
            }

            //
            // Finally, generate the output file.
            //

            var outputFile = new StreamWriter(baseDirectory + @"MutationDistribution.txt");
            outputFile.WriteLine("Gene\tnNone\tnOne\tnMultiple\tMutationsNone\tMutationsOne\tMutationsMultiple\tMedianNone\tMedianOne\tMedianMultiple\tMeanNone\tMeanOne\tMeanMultiple");

            foreach (var entry in singleMutations) {
                string gene = entry.Key;
                int nNone = noMutations[gene].Count();
                int nSingle = singleMutations[gene].Count();
                int nMultiple = multipleMutations[gene].Count();

                int countNone = noMutations[gene].Sum();
                int countSingle = singleMutations[gene].Sum();
                int countMultiple = multipleMutations[gene].Sum();

                noMutations[gene].Sort();
                singleMutations[gene].Sort();
                multipleMutations[gene].Sort();

                double meanNone, meanSingle, meanMutiple;

                if (nNone == 0) {
                    meanNone = 0;
                } else {
                    meanNone = ((double)countNone)/nNone;
                }

                if (nSingle == 0) {
                    meanSingle = 0;
                } else{
                    meanSingle = ((double)countSingle)/nSingle;
                }

                if (nMultiple == 0) {
                    meanMutiple = 0;
                } else {
                    meanMutiple = ((double)countMultiple)/nMultiple;
                }

                outputFile.WriteLine(gene + "\t" + 
                    nNone + "\t" + nSingle + "\t" + nMultiple + "\t" + 
                    countNone + "\t" + countSingle + "\t" + countMultiple + "\t" + 
                    MedianOfSortedList(noMutations[gene]) + "\t" + MedianOfSortedList(singleMutations[gene]) + "\t" + MedianOfSortedList(multipleMutations[gene]) + "\t" +
                    meanNone + "\t" + meanSingle + "\t" + meanMutiple);
            }

            outputFile.Close();
        }

        class GeneAndIsoform
        {
            public GeneAndIsoform(string gene_, string isoform)
            {
                gene = gene_;
                isoforms.Add(isoform);
            }

            public string gene;
            public List<string> isoforms = new List<string>();
        }

        public static void GenerateIsoformFileListsByDiseaseAndMutation(List<ExpressionTools.Experiment> experiments)
        {
            var genesToConsider = new List<GeneAndIsoform>();
            genesToConsider.Add(new GeneAndIsoform("TP53", "uc010cni"));
            genesToConsider.Add(new GeneAndIsoform("CDKN2A", "uc010miu"));
            genesToConsider.Add(new GeneAndIsoform("HUWE1", "uc004dsp"));
            genesToConsider.Add(new GeneAndIsoform("KEAP1", "uc002mor"));
            genesToConsider.Add(new GeneAndIsoform("IDH1", "uc002vcs"));
            genesToConsider.Add(new GeneAndIsoform("KRAS", "uc001rgp"));
            genesToConsider.Add(new GeneAndIsoform("EGFR", "uc022adn"));
            genesToConsider.Add(new GeneAndIsoform("PTEN", "uc001kfb"));
            genesToConsider.Add(new GeneAndIsoform("FAT1", "uc003izf"));
            genesToConsider.Add(new GeneAndIsoform("BRAF", "uc003vwc"));


            var generatorScript = new StreamWriter(isoformDirectory + "generate_isoform_usage.cmd");
            generatorScript.WriteLine(@"del isoform_counts.txt");

            var mutantStreamsByGene = new Dictionary<string, Dictionary<TumorType, StreamWriter>>();   // maps gene -> mapping of disease->stream
            var nonmutantStreamsByGene = new Dictionary<string, Dictionary<TumorType, StreamWriter>>();   // maps gene -> mapping of disease->stream

            var genesAndIsoformsFile = new StreamWriter(baseDirectory + @"genes_and_isoforms.txt");
            foreach (var gene in genesToConsider)
            {
                foreach (var isoform in gene.isoforms)
                {
                    genesAndIsoformsFile.WriteLine(gene.gene + "\t" + isoform);
                }
                mutantStreamsByGene.Add(gene.gene, new Dictionary<TumorType, StreamWriter>());
                nonmutantStreamsByGene.Add(gene.gene, new Dictionary<TumorType, StreamWriter>());
            }

            genesAndIsoformsFile.Close();

            foreach (var experiment in experiments)
            {
                if (experiment.TumorRNAAnalysis.storedBAM == null || experiment.TumorRNAAnalysis.storedBAM.isoformCount == null) {
                    // This experiment doesn't yet have what we need.  Skip it.
                    continue;
                }

                foreach (var gene in genesToConsider) {
                    int countOfNonSilentMutations = experiment.maf.Where(maf => maf.Hugo_symbol == gene.gene && maf.Variant_classification.ToLower() != "silent").Count();
                    Dictionary<string, Dictionary<TumorType, StreamWriter>> mapping;
                    string suffix;
                    if (countOfNonSilentMutations > 0)
                    {
                        mapping = mutantStreamsByGene;
                        suffix = "-mutant.txt";
                    } else {
                        mapping = nonmutantStreamsByGene;
                        suffix = "-non-mutant.txt";
                    }

                    if (!mapping[gene.gene].ContainsKey(experiment.disease_abbr))
                    {
                        string isoformFilesFileName = isoformDirectory + @"isoform_files-" + experiment.disease_abbr + "-" + gene.gene + suffix;
                        mapping[gene.gene].Add(experiment.disease_abbr, new StreamWriter(isoformFilesFileName));
                        foreach (var isoform in gene.isoforms)
                        {
                            generatorScript.WriteLine(@"findstr /f:" + isoformFilesFileName + " " + isoform + @">> isoform_counts.txt");
                        }
                    }
                    mapping[gene.gene][experiment.disease_abbr].WriteLine(experiment.TumorRNAAnalysis.storedBAM.isoformCount);
                } // foreach gene
            } // foreach experiment

            //
            // Now close all of the streams.  Wouldn't it be better if this just happened automatically at program exit time?
            //
            var mappings = new List<Dictionary<string, Dictionary<TumorType, StreamWriter>>>();
            mappings.Add(mutantStreamsByGene);
            mappings.Add(nonmutantStreamsByGene);
            foreach (var mapping in mappings)
            {
                foreach (var gene in mapping)
                {
                    foreach (var stream in gene.Value)
                    {
                        stream.Value.Close();
                    }
                }
            }
            generatorScript.Close();

        }


        public static void GenerateVariantCallingScript(List<ExpressionTools.Experiment> experiments)
        {
            int nSamplesToVariantCall = 0;
            int nSamplesAlreadyCalled = 0;
            int nSamplesNotYetReadyToCall = 0;
            int nSkipped = 0;

            var referencesNeedingAlternateChrState = new List<string>();

            StreamWriter variantCallScript = null;

            foreach (var experiment in experiments)
            {
                string refassem = experiment.NormalDNAAnalysis.refassemShortName.ToLower();
                if (refassemTranslation.ContainsKey(refassem))
                {
                    refassem = refassemTranslation[refassem].ToLower();
                }

                if (experiment.NormalDNAAnalysis.storedBAM == null)
                {
                    nSamplesNotYetReadyToCall++;
                }
                else if (experiment.NormalDNAAnalysis.storedBAM.vcfInfo != null) 
                {
                    nSamplesAlreadyCalled++;
                }
                else if (!experiment.NormalDNAAnalysis.localRealign && (!experiment.NormalDNAAnalysis.storedBAM.chrStateKnown ||
                    refassemWithoutChr.Contains(refassem) == experiment.NormalDNAAnalysis.storedBAM.usesChr && !refassemChrStateSwitched.ContainsKey(refassem)))
                {
                    if (experiment.NormalDNAAnalysis.storedBAM.chrStateKnown && !referencesNeedingAlternateChrState.Contains(refassem))
                    {
                        referencesNeedingAlternateChrState.Add(refassem);
                    }
                    nSkipped++;
                }
                else
                {
                    string referenceToUse;

                    if (!experiment.NormalDNAAnalysis.localRealign)
                    {
                        //
                        // We must have chrStateKnown or an alternate, because we checked for it above
                        //
                        if (refassemWithoutChr.Contains(refassem) == experiment.NormalDNAAnalysis.storedBAM.usesChr) // Subtle: reverse the chr state if the reference not-chr state is the same as the bam chr state.
                        {
                            if (refassemChrStateSwitched.ContainsKey(refassem))
                            {
                                referenceToUse = refassemChrStateSwitched[refassem];
                            } else {
                                Console.WriteLine("Got to unexpected place with ref " + refassem);
                                nSkipped++; // For now, just skip it.
                                referenceToUse = "This shouldn't happen";   // Because we picked it off above in the nSkipped case.  Need this to convince the compiler, though
                            }
                        } else {
                            referenceToUse = refassem;
                        }
                    }
                    else
                    {
                        referenceToUse = refassem;
                    }
                    nSamplesToVariantCall++;

                    if (variantCallScript == null)
                    {
                        variantCallScript = new StreamWriter(baseDirectory + @"callVariants");
                    }

                    variantCallScript.Write("date\n");
                    variantCallScript.Write(@"cat ~/genomes/" + referenceToUse + "-100k-regions | parallel -k -j `cat ~/ncores` \" freebayes --region {} --fasta-reference ~/genomes/" + 
                        referenceToUse + ".fa " +
                        ExpressionTools.WindowsToLinuxPathname(experiment.NormalDNAAnalysis.storedBAM.bamInfo.FullName) + " \" | ~/freebayes/vcflib/bin/vcffirstheader | ~/freebayes/vcflib/bin/vcfstreamsort -w 1000 | ~/freebayes/vcflib/bin/vcfuniq > " +
                        experiment.NormalDNAAnalysis.analysis_id + ".vcf\n");
                    variantCallScript.Write("if [ $? = 0 ]; then\n");
                    variantCallScript.Write(@"    cp " + experiment.NormalDNAAnalysis.analysis_id + ".vcf " + 
                        ExpressionTools.WindowsToLinuxPathname(ExpressionTools.GetDirectoryPathFromFullyQualifiedFilename(experiment.NormalDNAAnalysis.storedBAM.bamInfo.FullName)) + "\n");
                    variantCallScript.Write("else\n");
                    variantCallScript.Write(@"    echo " + experiment.NormalDNAAnalysis.analysis_id + " >> variant_calling_errors\n");
                    variantCallScript.Write("fi\n");
                    variantCallScript.Write("rm " + experiment.NormalDNAAnalysis.analysis_id + ".vcf\n");
                }
            }

            if (variantCallScript != null)
            {
                variantCallScript.Close();
            }

            Console.WriteLine("" + nSamplesAlreadyCalled + " have been variant called, " + nSamplesNotYetReadyToCall + " are waiting for other things, " + nSkipped + " are skipped, and " + nSamplesToVariantCall + " are ready to call.");
            if (referencesNeedingAlternateChrState.Count() != 0)
            {
                Console.Write("References needing alternate chr state: ");
                foreach (var reference in referencesNeedingAlternateChrState)
                {
                    Console.Write(reference + " ");
                }
                Console.WriteLine();
            }
        }

        static public void DumpSampleToParticipantIDMap( Dictionary<SampleID, ParticipantID> sampleToParticipantIdMap)
        {
            StreamWriter outputFile = new StreamWriter(baseDirectory + "sampleToParticipantIdMap.txt");
            outputFile.WriteLine("SampleID\tParticipantID");
            foreach (var entry in sampleToParticipantIdMap)
            {
                outputFile.WriteLine(entry.Key + "\t" + entry.Value);
            }
            outputFile.Close();
        }

        static public void GenerateAllcountScript(Dictionary<AnalysisID, ExpressionTools.TCGARecord> tcgaRecords)
        {
            int nWithAllcount = 0;
            int nAllcountsToGenerate = 0;
            int nAwaitingDownload = 0;
            int indexMachine = 0;
            StreamWriter allcountScript = null;
            foreach (var entry in tcgaRecords)
            {
                var tcgaRecord = entry.Value;
                AnalysisID analysisID = entry.Key;

                if (tcgaRecord.library_strategy == "RNA")
                {
                    if (tcgaRecord.storedBAM == null)
                    {
                        nAwaitingDownload++;
                    }
                    else if (tcgaRecord.storedBAM.allCountInfo != null)
                    {
                        nWithAllcount++;
                    }
                    else
                    {
                        nAllcountsToGenerate++;
                        if (allcountScript == null)
                        {
                            allcountScript = new StreamWriter(baseDirectory + "allcountCluster.cmd");
                            StreamWriter jobScript = new StreamWriter(baseDirectory + "createAllcountJob.cmd");
                            jobScript.WriteLine("job new /emailaddress:bolosky@microsoft.com /nodegroup:B99,ExpressQ /exclusive:true /failontaskfailure:false /jobname:allcount /memorypernode:10000 /notifyoncompletion:true /numnodes:1-40 /runtime:12:00 /scheduler:gcr /jobtemplate:ExpressQ /estimatedprocessmemory:10000");
                            jobScript.Close();
                        }
                        allcountScript.WriteLine(@"job add %1 /exclusive /numnodes:1-1 /scheduler:gcr " + @"\\gcr\scratch\b99\bolosky\countAll.cmd \\msr-genomics-" + indexMachine + @"\d$\sequence\indices\" + tcgaRecord.refassemShortName + "-24 " + tcgaRecord.storedBAM.bamInfo.FullName + " " +
                            tcgaRecord.storedBAM.bamInfo.DirectoryName + @"\" + analysisID + ".allcount.gz");
                        indexMachine = 1 - indexMachine;
                    }
                }
            } // foreach TCGA record
            if (allcountScript != null)
            {
                allcountScript.Close();
            }
            Console.WriteLine("Found " + nWithAllcount + " RNA analyses with allcounts, generated a script to build " + nAllcountsToGenerate + " more, and " + nAwaitingDownload + " are awaiting downloads");
        }

        static void GenerateListOfAllcountFiles(Dictionary<AnalysisID, ExpressionTools.StoredBAM> storedBAMs)
        {
            StreamWriter output = new StreamWriter(baseDirectory + "allcount_files.txt");

            foreach (var entry in storedBAMs)
            {
                var storedBAM = entry.Value;

                if (storedBAM.allCountInfo != null)
                {
                    output.WriteLine(storedBAM.allCountInfo.FullName);
                }
            }

            output.Close();
        }

        delegate ExpressionTools.TCGARecord RecordExtractor(ExpressionTools.Experiment experiment);
        static void GenerateTP53CountScript(RecordExtractor recordExtractor, List<ExpressionTools.Experiment> experiments, StreamWriter script)
         {
             foreach (var experiment in experiments)
             {
                 if (recordExtractor(experiment) != null && recordExtractor(experiment).storedBAM != null)
                 {
                     script.Write("samtools view " + recordExtractor(experiment).storedBAM.bamInfo.FullName + " ");
                     if (recordExtractor(experiment).refassemShortName.ToLower() == "hg19")
                     {
                         script.Write("chr");
                     }

                     if (recordExtractor(experiment).refassemShortName.ToLower() == "ncbi36_bccagsc_variant")
                     {
                         script.Write("17:7519125-7519126");   // hg18 coordinates (this is the only refassem with Tumor RNA aligned to any hg18 reference)
                     }
                     else
                     {
                         script.Write("17:7578400-7578401");   // hg19 coordinates
                     }

                     script.WriteLine(" > " + recordExtractor(experiment).analysis_id + ".tp53locus-reads.txt");
                 }
             }
         }

        static void GenerateTP53CountScripts(List<ExpressionTools.Experiment> experiments)
        {
            var script = new StreamWriter(baseDirectory + "generateTP53Count.cmd");
            GenerateTP53CountScript(e => e.TumorRNAAnalysis, experiments, script);
            script.Close();

            script = new StreamWriter(baseDirectory + "generateTP53Count-normal.cmd");
            GenerateTP53CountScript(e => e.NormalRNAAnalysis, experiments, script);
            script.Close();
        }

        public const string realignPathname = @"f:\sequence\reads\tcga\realigns.txt";

        static void GenerateListOfBamsNeedingChrState(Dictionary<AnalysisID, ExpressionTools.StoredBAM> storedBAMs)
        {
            int n = 0;
            int nNeedingChrState = 0;
            var writer = new StreamWriter(baseDirectory + "BamsNeedingChrState");
            foreach (var entry in storedBAMs)
            {
                var storedBAM = entry.Value;

                n++;
                if (!storedBAM.chrStateKnown)
                {
                    writer.WriteLine(storedBAM.bamInfo.FullName);
                    nNeedingChrState++;
                }
            }

            writer.Close();
            Console.WriteLine("" + nNeedingChrState + " of " + n + " bams need to have their chr state measured.");
        }

        static void AddRegionalExpressionProgramRunToScript(StreamWriter script, List<ParticipantID> participants, string disease_abbr)
        {
            if (participants.Count() == 0) {
                return;
            }
            string line = /*@"job add %1 /exclusive /numnodes:1-1 /scheduler:gcr" +*/ @"RegionalExpression f:\sequence\reads\tcga\expression\expression_" + disease_abbr + " 1000 ";
            foreach (var participant in participants) {
                line = line + participant + " ";
            }
            script.WriteLine(line);
        }

        static void GenerateRegionalExpressionScripts(List<ExpressionTools.Experiment> experiments)
        {
            int nWithRegionalExpression = 0;
            int nInScript = 0;
            int nNeedingPrecursors = 0;
            var perDiseaseLines = new Dictionary<string, List<ParticipantID>>(); // Maps disease->list of analyses waiting to be run.

            const int nDiseasesPerProgramRun = 800; // The command line is limited to 32K and GUIDs are 36 characters plus one for the space.  800 * 37 < 30K, so there's some headroom for the program name & args

            string createJobScriptFilename = baseDirectory + "createRegionalExpressionjob.cmd";
            string runJobScriptFilename = baseDirectory + @"runRegionalExpression.cmd";

            StreamWriter script = null;

            foreach (var experiment in experiments)
            {
                if (experiment.TumorRNAAnalysis.storedBAM == null)
                {
                    nNeedingPrecursors++;
                }
                else if (experiment.TumorRNAAnalysis.storedBAM.regionalExpressionInfo != null)
                {
                    nWithRegionalExpression++;
                }
                else
                {
                    nInScript++;

                    if (script == null)
                    {
                        var createJobScript = new StreamWriter(createJobScriptFilename);
                        createJobScript.WriteLine(@"job new /emailaddress:bolosky@microsoft.com /nodegroup:B99,ExpressQ /exclusive:true /failontaskfailure:false /jobname:regionalExpression /memorypernode:100000 /notifyoncompletion:true /numnodes:1-20 /runtime:12:00 /scheduler:gcr /jobtemplate:ExpressQ /estimatedprocessmemory:100000");
                        createJobScript.Close();
                        script = new StreamWriter(runJobScriptFilename);
                    }

                    string disease_abbr;
                    if (experiment.disease_abbr == "ov" && experiment.TumorRNAAnalysis.refassemShortName == "ncbi36_bccagsc_variant")
                    {
                        disease_abbr = "ov_hg18";    // Special case, since this is aligned against two different references
                    } else {
                        disease_abbr = experiment.disease_abbr;
                    }

                    if (!perDiseaseLines.ContainsKey(disease_abbr)) {
                        perDiseaseLines.Add(disease_abbr, new List<ParticipantID>());
                    }

                    perDiseaseLines[disease_abbr].Add(experiment.participant.participantId);
                    if (perDiseaseLines[disease_abbr].Count() >= nDiseasesPerProgramRun) {
                        AddRegionalExpressionProgramRunToScript(script, perDiseaseLines[disease_abbr], disease_abbr);
                        perDiseaseLines[disease_abbr] = new List<ParticipantID>();
                    }
                }
            } // foreach experiment

            if (null != script) {
                foreach (var entry in perDiseaseLines) {
                    AddRegionalExpressionProgramRunToScript(script, entry.Value, entry.Key);
                }
                script.Close();
            }
            else
            {
                //
                // Get rid of any stale scripts from previous runs.  File.Delete does nothing if the file doesn't exist.
                //
                File.Delete(createJobScriptFilename);
                File.Delete(runJobScriptFilename);
            }

            Console.WriteLine(nWithRegionalExpression + " analyses have regional expression, generated a script to make " + nInScript + " more and " + nNeedingPrecursors + " need to have their RNA data downloaded first.");
        }

        static void GenerateSelectedVariantsScript(List<ExpressionTools.Experiment> experiments)
        {
            int nWithSelectedVariants = 0;
            int nNeedingPrecursors = 0;
            int nReadyToGo = 0;

            string scriptFilename = baseDirectory + "generateSelectedVariants.cmd";
            File.Delete(scriptFilename);

            StreamWriter scriptFile = null;
            int nOnCurrentLine = 0;

            foreach (var experiment in experiments)
            {
                if (experiment.NormalDNAAnalysis.storedBAM == null || experiment.NormalDNAAnalysis.storedBAM.vcfInfo == null) {
                    nNeedingPrecursors++;
                }
                else if (experiment.NormalDNAAnalysis.storedBAM.selectedVariantsInfo != null)
                {
                    nWithSelectedVariants++;
                }
                else
                {
                    if (null == scriptFile)
                    {
                        scriptFile = new StreamWriter(scriptFilename);
                    }

                    if (nOnCurrentLine >= 800)
                    {
                        scriptFile.WriteLine();
                        nOnCurrentLine = 0;
                    }

                    if (0 == nOnCurrentLine)
                    {
                        scriptFile.Write("SelectGermlineVariants");
                    }

                    scriptFile.Write(" " + experiment.participant.participantId);
                    nOnCurrentLine++;
                    nReadyToGo++;
                }
            }

            if (null != scriptFile)
            {
                scriptFile.WriteLine();
                scriptFile.Close();
            }

            Console.WriteLine("" + nWithSelectedVariants + " have selected variants, " + nNeedingPrecursors + " aren't ready to select and " + nReadyToGo + " added to script to select.");
        }

        static void GenerateGeneExpressionScripts(List<ExpressionTools.Experiment> experiments)
        {
            int nWithGeneExpression = 0;
            int nInScript = 0;
            int nNeedingPrecursors = 0;
            var readyToGo = new List<ParticipantID>();

            string createJobScriptFilename = baseDirectory + "createGeneExpressionJob.cmd";
            string scheduleJobScriptFilename = baseDirectory + "scheduleGeneExpressionJob.cmd";

            foreach (var experiment in experiments)
            {
                if (experiment.TumorRNAAnalysis.storedBAM == null || experiment.TumorRNAAnalysis.storedBAM.regionalExpressionInfo == null)
                {
                    nNeedingPrecursors++;
                }
                else if (experiment.TumorRNAAnalysis.storedBAM.geneExpressionInfo != null)
                {
                    nWithGeneExpression++;
                }
                else
                {
                    readyToGo.Add(experiment.participant.participantId);
                }
            }

            nInScript = readyToGo.Count();
            if (0 != nInScript) {
                var createJobScript = new StreamWriter(createJobScriptFilename);
                createJobScript.WriteLine(@"job new /emailaddress:bolosky@microsoft.com /nodegroup:B99,ExpressQ /exclusive:true /failontaskfailure:false /jobname:geneExpression /memorypernode:32000 /notifyoncompletion:true /numnodes:1-20 /runtime:12:00 /scheduler:gcr /jobtemplate:ExpressQ /estimatedprocessmemory:20000");
                createJobScript.Close();

                int nPerJob = 40;
                int nInThisJob = 0;
                string jobHeader = @"job add %1 /exclusive /numnodes:1-1 /scheduler:gcr \\gcr\scratch\b99\bolosky\ExpressionNearMutations ";
                string thisJob = "";

                StreamWriter scheduleJobScript = new StreamWriter(scheduleJobScriptFilename);

                foreach (var participant in readyToGo)
                {
                    if (nInThisJob >= nPerJob) {
                        scheduleJobScript.WriteLine(jobHeader + thisJob);
                        thisJob = "";
                        nInThisJob = 0;
                    }

                    thisJob += participant + " ";
                    nInThisJob++;
                }

                if (0 != nInThisJob) {
                    scheduleJobScript.WriteLine(jobHeader + thisJob);
                }

                scheduleJobScript.Close();

            }
            else
            {
                //
                // Get rid of any stale scripts left over from a previous run.  File.Delete does nothing if the file doesn't exist.
                //
                File.Delete(createJobScriptFilename);
                File.Delete(scheduleJobScriptFilename);
            }

            Console.WriteLine("" + nInScript + " are ready to have gene expression computed, " + nNeedingPrecursors + " still need preliminary work done and " + nWithGeneExpression + " are done.");
        }


        static void CountTP53Mutations(List<ExpressionTools.Experiment> experiments)
        {
            var counts = new Dictionary<int, int>();
            counts.Add(0, 0);
            counts.Add(1, 0);

            int nBiggerThanOne = 0;

            foreach (var experiment in experiments)
            {
                int count = 0;
                foreach (var mutation in experiment.participant.mafs[0])
                {
                    if (mutation.Hugo_symbol.ToLower() == "tp53")
                    {
                        count++;
                    }
                }

                if (!counts.ContainsKey(count))
                {
                    counts.Add(count, 0);
                }

                counts[count]++;
                if (count > 1)
                {
                    nBiggerThanOne++;
                }
            }

            Console.WriteLine("" + counts[0] + " have no TP53 mutations, " +counts[1] + " have exactly one and " + nBiggerThanOne + " have more than one.");
            foreach (var count in counts)
            {
                Console.WriteLine("" + count.Key + ": " + count.Value);
            }

        }
        static void Main(string[] args)
        {
            var timer = new Stopwatch();
            timer.Start();

            hg18_likeReferences.Add("NCBI36_BCCAGSC_variant".ToLower());
            hg18_likeReferences.Add("NCBI36_BCM_variant".ToLower());
            hg18_likeReferences.Add("NCBI36_WUGSC_variant".ToLower());
            hg18_likeReferences.Add("hg18");
            hg18_likeReferences.Add("HG18_Broad_variant".ToLower());

            refassemWithoutChr.Add("NCBI36_BCCAGSC_variant".ToLower());
            refassemWithoutChr.Add("GRCH37-lite".ToLower());
            refassemWithoutChr.Add("hg19_broad_variant");
            refassemWithoutChr.Add("HS37D5".ToLower());
            refassemWithoutChr.Add("ncbi36.54");
            refassemWithoutChr.Add("NCBI36_WUGSC_variant".ToLower());
            refassemWithoutChr.Add("grch37-lite-+-hpv_redux-build");

            refassemChrStateSwitched.Add("hg19", "hg19-no-chr");
            refassemChrStateSwitched.Add("hs37d5", "hs37d5-with-chr");
            refassemChrStateSwitched.Add("grch37", "grch37-no-chr");
            refassemChrStateSwitched.Add("grch37-lite", "grch37-lite-with-chr");
            refassemChrStateSwitched.Add("hg19_broad_variant", "hg19_broad_variant-with-chr");
            refassemChrStateSwitched.Add("NCBI36_BCM_variant".ToLower(), "NCBI36_BCM_variant-no-chr".ToLower());
            refassemChrStateSwitched.Add("NCBI36_WUGSC_variant".ToLower(), "NCBI36_WUGSC_variant-with-chr".ToLower());
            refassemTranslation.Add("grch37-lite_wugsc_variant_1", "grch37-lite");


            InitializeMachines();

            //var geneToSize = GenerateSizeOfProteinCodingRegionMap();

            List<AnalysisID> excludedAnalyses = ExpressionTools.LoadExcludedAnalyses();

            var tumorToMachineMapping = GenerateTumorToMachineMapping();
            var storedBAMs = LoadStoredBAMs(tumorToMachineMapping);
            ExpressionTools.LoadChrStateFile( @"f:\sequence\reads\tcga\chrState", storedBAMs);
            GenerateListOfBamsNeedingChrState(storedBAMs);
            GenerateListOfAllcountFiles(storedBAMs);
            var tcgaRecords = ExpressionTools.LoadTCGARecords(storedBAMs, excludedAnalyses, @"f:\sequence\Reads\tcga-all.xml");
            ExpressionTools.LoadTCGARecordsForLocalRealigns(tcgaRecords, storedBAMs, realignPathname);
            ExpressionTools.LoadTCGAAdditionalMetadata(tcgaRecords);
            GenereateAdditionalTCGAMetadataGeneratingScript(tcgaRecords, storedBAMs, tumorToMachineMapping);
            VerifyStoredBAMPaths(storedBAMs, tcgaRecords, tumorToMachineMapping);
            GenerateAllcountScript(tcgaRecords);

            var sampleToParticipantIDMap = ExpressionTools.CreateSampleToParticipantIDMap(tcgaRecords);
            DumpSampleToParticipantIDMap(sampleToParticipantIDMap);
            var centerNameToDomainNameMap = GenerateCenterNameToDomainName();
            var domainNameToCenterNameMap = InvertStringDictionary(centerNameToDomainNameMap);

            //
            // Now build up participant data.
            //
            Dictionary<string, ExpressionTools.Sample> allSamples;
            var participants = ExpressionTools.BuildParticipantData(tcgaRecords, out allSamples);

            //
            // Load all the MAF files.
            //
            //var mutations = new Dictionary<string, Dictionary<string, List<ExpressionTools.MAFRecord>>>();
            ExpressionTools.AddAllMAFFilesToParticipants(participants, sampleToParticipantIDMap);


            //BuildMAFsByGene(participants);
            //AddCountsToMAFs(participants, allSamples);

            var experiments = BuildExperiments(participants);
            //CountTP53Mutations(experiments);
            ExpressionTools.DumpExperimentsToFile(experiments, baseDirectory + "experiments.txt");
            GenerateSelectedVariantsScript(experiments);
            GenerateRegionalExpressionScripts(experiments);
            GenerateGeneExpressionScripts(experiments);
            GenerateTP53CountScripts(experiments);
            GenerateIsoformFileListsByDiseaseAndMutation(experiments);
            //GenerateMutationDistributionByGene(experiments);
            GenerateMakeIsoformsScripts(experiments);
            GenerateRealignmentScript(experiments, tcgaRecords, storedBAMs, tumorToMachineMapping);
            GenerateUnneededLists(storedBAMs);
            GenerateRealignmentAnalyses(experiments);
            GenerateDownloadScripts(experiments, tumorToMachineMapping);

            GenerateExtractionScripts(experiments);
            GenerateVariantCallingScript(experiments);
            //GenerateLAMLMoves(tcgaRecords);

            timer.Stop();
            Console.WriteLine("Expression metadata took " + (timer.ElapsedMilliseconds + 500) / 1000 + "s.");
        }
    }
}
