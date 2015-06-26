using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using System.IO;
using System.Xml.Linq;

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
        static List<string> hg18_likeReferences = new List<string>();

        public struct AnalysisType
        {
            public AnalysisType(TumorType tumorType_, LibraryStrategy libraryStrategy_, bool isNormal_)
            {
                tumorType = tumorType_.ToLower();
                libraryStrategy = libraryStrategy_.ToLower();
                isNormal = isNormal_;
            }

            public AnalysisType(TCGARecord tcgaRecord)
            {
                tumorType = tcgaRecord.disease_abbr.ToLower();
                libraryStrategy = tcgaRecord.library_strategy.ToLower();
                isNormal = !tcgaRecord.tumorSample;
            }

            public TumorType   tumorType;
            public LibraryStrategy libraryStrategy;
            public bool isNormal;
        };

        static List<LibraryStrategy> LibraryStrategies = new List<LibraryStrategy>();

        public static void AddEntireTumorToMachine(Dictionary<AnalysisType, Pathname> tumorToMachineMapping, TumorType tumorType, Pathname machine) 
        {
            foreach (string libraryStrategy in LibraryStrategies)
            {
                var analysisType = new AnalysisType(tumorType, libraryStrategy, true);
                tumorToMachineMapping.Add(analysisType, machine);

                analysisType = new AnalysisType(tumorType, libraryStrategy, false);
                tumorToMachineMapping.Add(analysisType, machine);
            }
        }

        public class Machine
        {
            public Machine(string name_, int memoryInGB_, int diskInTB_)
            {
                name = name_;
                memoryInGB = memoryInGB_;
                diskInTB = diskInTB_;
            }

            public static void AddMachine(string name_, int memoryInGB_, int diskInTB_)
            {
                Machines.Add(name_, new Machine(name_, memoryInGB_, diskInTB_));
            }

            public string name;
            public int memoryInGB;
            public int diskInTB;
        }

        public static void InitializeMachines()
        {
            Machine.AddMachine("msr-srs-0", 8, 36);
            Machine.AddMachine("msr-srs-1", 8, 36);
            Machine.AddMachine("msr-srs-2", 8, 36);
            Machine.AddMachine("msr-srs-3", 8, 36);
            Machine.AddMachine("msr-srs-4", 8, 36);
            Machine.AddMachine("msr-srs-5", 8, 36);
            Machine.AddMachine("msr-srs-6", 8, 36);
            Machine.AddMachine("msr-srs-7", 48, 48);
            Machine.AddMachine("msr-srs-8", 48, 36);
            Machine.AddMachine("fds-326-k25-1", 48, 6);
            Machine.AddMachine("fds-326-k25-2", 48, 6);
            Machine.AddMachine("fds-326-k25-3", 48, 6);
            Machine.AddMachine("fds-326-k25-4", 48, 6);
            Machine.AddMachine("fds-326-k25-5", 48, 6);
            Machine.AddMachine("fds-326-k25-6", 48, 6);
            Machine.AddMachine("fds-326-k25-7", 48, 6);
            Machine.AddMachine("fds-326-k25-8", 48, 6);
            Machine.AddMachine("fds-326-k25-9", 48, 6);
            Machine.AddMachine("fds-326-k25-11", 48, 6);
            Machine.AddMachine("fds-326-k25-12", 48, 6);
            Machine.AddMachine("fds-326-k25-14", 48, 6);
            Machine.AddMachine("fds-326-k25-15", 48, 6);
            Machine.AddMachine("fds-326-k25-16", 48, 6);
            Machine.AddMachine("fds-326-k25-17", 48, 6);
            Machine.AddMachine("fds-326-k25-18", 48, 6);
            Machine.AddMachine("fds-326-k25-20", 48, 6);
        }

        public static Dictionary<string, Machine> Machines = new Dictionary<string, Machine>();

        public static Dictionary<AnalysisType, Pathname> GenerateTumorToMachineMapping()
        {
            var tumorToMachineMapping = new Dictionary<AnalysisType, Pathname>();

            AddEntireTumorToMachine(tumorToMachineMapping, "acc", "fds-326-k25-2");
            AddEntireTumorToMachine(tumorToMachineMapping, "blca", "msr-srs-3");
            AddEntireTumorToMachine(tumorToMachineMapping, "chol", "msr-srs-6");
            AddEntireTumorToMachine(tumorToMachineMapping, "gbm", "msr-srs-4");
            AddEntireTumorToMachine(tumorToMachineMapping, "hnsc", "msr-srs-0");
            AddEntireTumorToMachine(tumorToMachineMapping, "kich", "fds-326-k25-9");
            AddEntireTumorToMachine(tumorToMachineMapping, "kirp", "msr-srs-6");
            AddEntireTumorToMachine(tumorToMachineMapping, "laml", "msr-srs-8");
            AddEntireTumorToMachine(tumorToMachineMapping, "luad", "msr-srs-6");
            AddEntireTumorToMachine(tumorToMachineMapping, "ov", "msr-srs-7");
            AddEntireTumorToMachine(tumorToMachineMapping, "paad", "fds-326-k25-5");
            AddEntireTumorToMachine(tumorToMachineMapping, "pcpg", "fds-326-k25-1");
            AddEntireTumorToMachine(tumorToMachineMapping, "prad", "msr-srs-2");
            AddEntireTumorToMachine(tumorToMachineMapping, "read", "fds-326-k25-7");
            AddEntireTumorToMachine(tumorToMachineMapping, "stad", "msr-srs-5");
            AddEntireTumorToMachine(tumorToMachineMapping, "thca", "msr-srs-4");
            AddEntireTumorToMachine(tumorToMachineMapping, "ucec", "msr-srs-2");
            AddEntireTumorToMachine(tumorToMachineMapping, "ucs", "msr-srs-4");
            AddEntireTumorToMachine(tumorToMachineMapping, "uvm", "msr-srs-4");

            tumorToMachineMapping.Add(new AnalysisType("brca", "rna", true), "fds-326-k25-4");
            tumorToMachineMapping.Add(new AnalysisType("brca", "wgs", true), "fds-326-k25-4");
            tumorToMachineMapping.Add(new AnalysisType("brca", "wxs", true), "msr-srs-1");
            tumorToMachineMapping.Add(new AnalysisType("brca", "rna", false), "msr-srs-1");
            tumorToMachineMapping.Add(new AnalysisType("brca", "wgs", false), "msr-srs-1");
            tumorToMachineMapping.Add(new AnalysisType("brca", "wxs", false), "msr-srs-1");

            tumorToMachineMapping.Add(new AnalysisType("cesc", "rna", true), "msr-srs-7");
            tumorToMachineMapping.Add(new AnalysisType("cesc", "wgs", true), "msr-srs-7");
            tumorToMachineMapping.Add(new AnalysisType("cesc", "wxs", true), "msr-srs-7");
            tumorToMachineMapping.Add(new AnalysisType("cesc", "rna", false), "fds-326-k25-15");
            tumorToMachineMapping.Add(new AnalysisType("cesc", "wgs", false), "fds-326-k25-15");
            tumorToMachineMapping.Add(new AnalysisType("cesc", "wxs", false), "fds-326-k25-15");

            tumorToMachineMapping.Add(new AnalysisType("esca", "rna", true), "fds-326-k25-11");
            tumorToMachineMapping.Add(new AnalysisType("esca", "wgs", true), "fds-326-k25-11");
            tumorToMachineMapping.Add(new AnalysisType("esca", "wxs", true), "fds-326-k25-11");
            tumorToMachineMapping.Add(new AnalysisType("esca", "rna", false), "fds-326-k25-11");
            tumorToMachineMapping.Add(new AnalysisType("esca", "wgs", false), "fds-326-k25-11");
            tumorToMachineMapping.Add(new AnalysisType("esca", "wxs", false), "msr-srs-7");

            tumorToMachineMapping.Add(new AnalysisType("lgg", "rna", true), "msr-srs-3");
            tumorToMachineMapping.Add(new AnalysisType("lgg", "wgs", true), "msr-srs-3");
            tumorToMachineMapping.Add(new AnalysisType("lgg", "wxs", true), "msr-srs-3");
            tumorToMachineMapping.Add(new AnalysisType("lgg", "rna", false), "fds-326-k25-14");
            tumorToMachineMapping.Add(new AnalysisType("lgg", "wgs", false), "fds-326-k25-14");
            tumorToMachineMapping.Add(new AnalysisType("lgg", "wxs", false), "fds-326-k25-3");

            tumorToMachineMapping.Add(new AnalysisType("lihc", "rna", true), "fds-326-k25-6");
            tumorToMachineMapping.Add(new AnalysisType("lihc", "wgs", true), "fds-326-k25-6");
            tumorToMachineMapping.Add(new AnalysisType("lihc", "wxs", true), "fds-326-k25-6");
            tumorToMachineMapping.Add(new AnalysisType("lihc", "rna", false), "msr-srs-2");
            tumorToMachineMapping.Add(new AnalysisType("lihc", "wgs", false), "fds-326-k25-20");
            tumorToMachineMapping.Add(new AnalysisType("lihc", "wxs", false), "msr-srs-2");

            tumorToMachineMapping.Add(new AnalysisType("skcm", "rna", true), "fds-326-k25-8");
            tumorToMachineMapping.Add(new AnalysisType("skcm", "wgs", true), "fds-326-k25-8");
            tumorToMachineMapping.Add(new AnalysisType("skcm", "wxs", true), "fds-326-k25-8");
            tumorToMachineMapping.Add(new AnalysisType("skcm", "rna", false), "msr-srs-5");
            tumorToMachineMapping.Add(new AnalysisType("skcm", "wgs", false), "msr-srs-5");
            tumorToMachineMapping.Add(new AnalysisType("skcm", "wxs", false), "msr-srs-5");

            tumorToMachineMapping.Add(new AnalysisType("lusc", "rna", true), "fds-326-k25-18");
            tumorToMachineMapping.Add(new AnalysisType("lusc", "wgs", true), "fds-326-k25-18");
            tumorToMachineMapping.Add(new AnalysisType("lusc", "wxs", true), "fds-326-k25-18");
            tumorToMachineMapping.Add(new AnalysisType("lusc", "rna", false), "fds-326-k25-16");
            tumorToMachineMapping.Add(new AnalysisType("lusc", "wgs", false), "fds-326-k25-19");
            tumorToMachineMapping.Add(new AnalysisType("lusc", "wxs", false), "fds-326-k25-16");

             return tumorToMachineMapping;
        }

        public class StoredBAM
        {
            public AnalysisID analysisID;
            public Pathname directoryName;
            public FileInfo bamInfo = null;
            public FileInfo baiInfo = null;
            public long totalSize;
            public bool needed = false; // Do we need this as input for some experiment?
            public bool isRealingned = false;   // Was this realined, or is it straight from TCGA
        }


        static void LoadStoredBAMsForDirectory(Pathname directory, Dictionary<AnalysisID, StoredBAM> storedBAMs) 
        {
            if (!Directory.Exists(directory))
            {
                Console.WriteLine("Unable to find expected tcga data directory " + directory);
                return;
            }

            foreach (var subdir in Directory.GetDirectories(directory))
            {
                var pathnameComponents = subdir.Split('\\');
                AnalysisID analysisID = pathnameComponents[pathnameComponents.Count() - 1].ToLower();
                if (analysisID.Count() == 36)  // The length of a GUID string
                {
                    var storedBAM = new StoredBAM();
                    storedBAM.analysisID = analysisID;

                    bool hasTarFile = false;    // Some RNA analyses include a tarred FASTQ.  If we see one, don't generate a warning for not having BAM files.

                    foreach (var file in Directory.GetFiles(subdir))
                    {
                        if (file.Count() > 4 && file.Substring(file.Count() - 4) == ".bam")
                        {
                            if (null != storedBAM.bamInfo)
                            {
                                Console.WriteLine("Saw multiple BAM files in the same analysis directory " + file + " and " + storedBAM.bamInfo.FullName);
                            }
                            storedBAM.bamInfo = new FileInfo(file);
                        }
                        else if (file.Count() > 4 && file.Substring(file.Count() - 4) == ".bai")
                        {
                            if (null != storedBAM.baiInfo)
                            {
                                Console.WriteLine("Saw multiple BAI files in the same analysis directory " + file + " and " + storedBAM.baiInfo.FullName);
                            }
                            storedBAM.baiInfo = new FileInfo(file);
                        } else if (file.Count() > 4 && file.Substring(file.Count() -4) == ".tar") {
                            hasTarFile = true;
                        }
                        else
                        {
                            Console.WriteLine("Saw unexpected file " + file);
                        }
                    }

                    if (storedBAM.baiInfo == null || storedBAM.bamInfo == null)
                    {
                        if (!hasTarFile) {
                            Console.WriteLine("Analysis directory " + subdir + " doesn't contain both .bam and .bai files");
                        }
                        continue;
                    }

                    storedBAM.totalSize = storedBAM.bamInfo.Length + storedBAM.baiInfo.Length;
                    storedBAM.directoryName = subdir;

                    if (storedBAMs.ContainsKey(analysisID))
                    {
                        Console.WriteLine("Duplicate stored BAM, " + subdir + " and " + storedBAMs[analysisID].directoryName);
                    }
                    else
                    {
                        storedBAMs.Add(analysisID, storedBAM);
                    }
                } else if (analysisID.Contains(".partial")) {
                    Console.WriteLine("Partial download at " + subdir);
                }
            }
        }

        public static Dictionary<AnalysisID, StoredBAM> LoadStoredBAMs(Dictionary<AnalysisType, Pathname> tumorToMachineMapping)
        {
            var storedBAMs = new Dictionary<AnalysisID, StoredBAM>();

            //
            // Directory structure is \\msr-srs-%n\d$\tcga\{rna,wgs,wxs}\{tumor, normal}\disease_abbr\analysis_id\*.{bam,bai}.  We need to call LoadStoredBAMsForDirectory on each
            // of the disease_abbr directories.  In addition, there are files on \\erg00\r$\tcga\laml.
            //


            foreach (var entry in tumorToMachineMapping)
            {
                LoadStoredBAMsForDirectory(@"\\" + entry.Value + @"\d$\tcga\" + entry.Key.tumorType + @"\" + (entry.Key.isNormal ? "normal" : "tumor") + @"\" + entry.Key.libraryStrategy, storedBAMs);
            }

            LoadStoredBAMsForDirectory(@"\\erg00\r$\tcga\laml", storedBAMs);

            return storedBAMs;
        }

        public class TCGARecord
        {
            public AnalysisID analysis_id;
            public ParticipantID participant_id;
            public Pathname bamFileName;
            public Pathname baiFileName;
            public DateTime lastModified;
            public DateTime uploadDate;
            public SampleID aliquot_id;
            public string refassemShortName;
            public TumorType disease_abbr;
            public long totalFileSize;
            public SampleType sampleType;
            public LibraryStrategy library_strategy;
            public string center_name;
            public bool tumorSample;
            public bool localRealign;
            public AnalysisID realignedFrom = null;
            public TCGARecord realignSource = null;
            public bool isWGS;

            public bool hasAdditionalData;  // Do we have the extra metadata that's shown below?
            public int readSampleSize;
            public int minReadLength;
            public int maxReadLength;
            public float meanReadLength;
            public int minGoodBases;
            public int maxGoodBases;
            public float meanGoodBases;
            public bool anyPaired;
            public bool allPaired;


            public string chromPrefix;
            public StoredBAM storedBAM;

            public string GetContainingDirectory()
            {
                return @"d:\tcga\" + library_strategy + (tumorSample ? @"\tumor\" : @"\normal\") + disease_abbr;
            }
        }

        public static Dictionary<AnalysisID, TCGARecord> LoadTCGARecords(Dictionary<AnalysisID, StoredBAM> storedBAMs)
        {
            var tcgaRecords = new Dictionary<AnalysisID, TCGARecord>();

            XDocument tcga_all = XDocument.Load(@"f:\sequence\Reads\tcga-all.xml");

            var samplesQuery = from Result in tcga_all.Element("ResultSet").Elements("Result")
                               select Result;

            foreach (var result in samplesQuery)
            {
                TCGARecord record = new TCGARecord();

                SampleType sample_type = result.Element("sample_type").Value;
                if (sample_type.Count() < 2) {
                    continue;
                }

                if (sample_type.Count() > 2) {
                    if (sample_type.Count() != 9 || sample_type.Substring(0, 7) != "TARGET_") {
                        Console.WriteLine("Unknown sample type " + sample_type);
                        continue;
                    }

                    sample_type = sample_type.Substring(7,2);
                }

                int sampleTypeInt = Convert.ToInt32(sample_type);
                if (sampleTypeInt <= 0 || sampleTypeInt > 19)
                {
                    continue;
                }

                LibraryStrategy library_strategy = result.Element("library_strategy").Value;

                if (library_strategy != "RNA-Seq" && library_strategy != "WXS" && library_strategy != "WGS")
                {
                    // miRNA or validation or something else we don't care about.
                    continue;
                }

                if (library_strategy == "RNA-Seq")
                {
                    library_strategy = "RNA";   // Just makes my life easier.
                }

                record.analysis_id = result.Element("analysis_id").Value.ToLower();
                record.aliquot_id = result.Element("aliquot_id").Value.ToLower();
                record.center_name = result.Element("center_name").Value.ToLower();
                record.disease_abbr = result.Element("disease_abbr").Value.ToLower();
                record.lastModified = Convert.ToDateTime(result.Element("last_modified").Value);
                record.library_strategy = library_strategy;
                record.participant_id = result.Element("participant_id").Value.ToLower();
                record.refassemShortName = result.Element("refassem_short_name").Value.ToLower();
                record.sampleType = sample_type;
                record.uploadDate = Convert.ToDateTime(result.Element("upload_date").Value);
                record.tumorSample = sampleTypeInt >= 1 && sampleTypeInt <= 9;
                record.isWGS = record.library_strategy.ToLower() == "wgs";

                record.bamFileName = "";
                record.baiFileName = "";
                record.totalFileSize = 0;
                foreach (var file in result.Element("files").Elements("file"))
                {
                    var fileNameElement = file.Element("filename");
                    Pathname fileName = fileNameElement.Value;
                    if (fileName.Length > 4 && ".bam" == fileName.Substring(fileName.Length - 4))
                    {
                        record.bamFileName = fileName;
                    }
                    else if (fileName.Length > 4 && ".bai" == fileName.Substring(fileName.Length - 4))
                    {
                        record.baiFileName = fileName;
                    }
                    record.totalFileSize += Convert.ToInt64(file.Element("filesize").Value);
                }
                record.localRealign = false;

                string refassem = record.refassemShortName;
                if (refassem == "NCBI36_BCCAGSC_variant" || refassem == "GRCh37-lite" || refassem == "GRCh37" || refassem == "HG19_Broad_variant" || refassem == "GRCh37-lite_WUGSC_variant_1" || refassem == "GRCh37-lite_WUGSC_variant_2"
                    || refassem == "GRCh37-lite-+-HPV_Redux-build" || refassem == "HS37D5" || refassem == "NCBI36_WUGSC_variant")
                {
                    record.chromPrefix = "";
                }
                else
                {
                    record.chromPrefix = " chr";
                }

                if (storedBAMs.ContainsKey(record.analysis_id)) {
                    record.storedBAM = storedBAMs[record.analysis_id];
                    if (record.storedBAM.bamInfo.Length + record.storedBAM.baiInfo.Length != record.totalFileSize)
                    {
                        Console.WriteLine("File size mismatch for " + record.storedBAM.bamInfo.FullName + " should have total size " + record.totalFileSize + " but has " + record.storedBAM.bamInfo.Length + record.storedBAM.baiInfo.Length);
                    }
                } else {
                    record.storedBAM = null;
                }

                tcgaRecords.Add(record.analysis_id, record);
            }

            return tcgaRecords;

        }

        public static void LoadTCGAAdditionalMetadata(Dictionary<AnalysisID, TCGARecord> tcgaRecords)
        {    
            string[] additionalMetadata = File.ReadAllLines(@"f:\sequence\reads\tcga\tcgaAdditionalMetadata.txt");
            foreach (var line in additionalMetadata)
            {
                string[] fields = line.Split('\t'); // Format is 0:analysisID 1:filename 2:readsSampled 3:minLength 4:maxLength 5:totalLength 6:minGood 7:maxGood 8:totalGoodBases 9:anyPaired 10:allPaired

                if (!tcgaRecords.ContainsKey(fields[0]))
                {
                    Console.WriteLine("Found additional metadata for unknown analysis " + line);
                }
                else
                {
                    var record = tcgaRecords[fields[0]];
                    record.hasAdditionalData = true;
                    record.readSampleSize = Convert.ToInt32(fields[2]);
                    record.minReadLength = Convert.ToInt32(fields[3]);
                    record.maxReadLength = Convert.ToInt32(fields[4]);
                    record.meanReadLength = Convert.ToSingle(fields[5]) / (float)record.readSampleSize;
                    record.minGoodBases = Convert.ToInt32(fields[6]);
                    record.maxGoodBases = Convert.ToInt32(fields[7]);
                    record.meanGoodBases = Convert.ToSingle(fields[8]) / (float)record.readSampleSize;
                    record.anyPaired = fields[9] == "1";
                    record.allPaired = fields[10] == "1";
                }
            }

        }

        static void VerifyStoredBAMPaths(Dictionary<AnalysisID, StoredBAM> storedBAMs, Dictionary<AnalysisID, TCGARecord> tcgaRecords, Dictionary<AnalysisType, Pathname> tumorTypeToMachineMapping)
        {
            //
            // Make sure that all of the analyses are in the right place.  They should be in a path that looks like \\machine\d$\tcga\{disease}\{tumor, normal}\{RNA, WGS, WXS}\<analysis-id>\X.bam
            // Or \\erg00\r$\laml\<analysis-id>\X.bam
            //

            const int diseaseComponent = 5;
            const int tumorNormalComponent = 6;
            const int libraryComponent = 7;

            int nRealigned = 0;

            foreach (var bamEntry in storedBAMs)
            {
                AnalysisID analysisID = bamEntry.Key;
                StoredBAM storedBAM = bamEntry.Value;
                string fullPath = storedBAM.bamInfo.FullName.ToLower();
                string[] pathComponents = fullPath.Split('\\');
                int pathLength = pathComponents.Count();    // Remember, this is one more than you think, because the first one is the empty string due to \\
                for (int i = 0; i < pathLength; i++)
                {
                    pathComponents[i] = pathComponents[i].ToLower();
                }

                if (!tcgaRecords.ContainsKey(analysisID)) {
                    //
                    // Special case it if it's stored on erg00, which has a bunch of uncategorized mirna and stuff.
                    //
                    if (!(pathLength == 8 && pathComponents[2] == "erg00"))
                    {
                        Console.WriteLine("Stored BAM doesn't have accompanying TCGA record, pathname " + storedBAM.bamInfo.FullName);
                    }
                    continue;
                }

                TCGARecord tcgaRecord = tcgaRecords[analysisID];

                if (tcgaRecord.localRealign)
                {
                    nRealigned++;
                }

                if (pathLength < 8 || pathComponents[0] != "" || pathComponents[1] != "")
                {
                    Console.WriteLine("Unparsable directory path for stored BAM: " + fullPath);
                    continue;
                }

                if (pathComponents[2] == "erg00")
                {
                    if (pathLength != 8 || pathComponents[3] != "r$" || pathComponents[4] != "tcga" || pathComponents[5] != "laml" || pathComponents[6] != analysisID)
                    {
                        Console.WriteLine("Invalid erg00 path for analysis ID " + analysisID + " path " + fullPath);
                        continue;
                    }

                    if (tcgaRecord.disease_abbr != "laml")
                    {
                        Console.WriteLine("Non-AML sample stored on erg00, analysisID " + analysisID + " tumor type " + tcgaRecord.disease_abbr);
                        continue;
                    }
                }
                else
                {
                    if (pathLength != 10 || pathComponents[3] != "d$" || pathComponents[4] != "tcga" || pathComponents[8] != analysisID ||
                        (pathComponents[libraryComponent] != "rna" && pathComponents[libraryComponent] != "wgs" && pathComponents[libraryComponent] != "wxs") || (pathComponents[tumorNormalComponent] != "tumor" && pathComponents[tumorNormalComponent] != "normal"))
                    {
                        Console.WriteLine("Invalid non-erg00 path for analysis ID " + analysisID + " path " + fullPath);
                        continue;
                    }

                    //
                    // Its shape is right.  Check to see if it's actually the right place.
                    //
                    if (tcgaRecord.library_strategy.ToLower() != pathComponents[libraryComponent])
                    {
                        Console.WriteLine("Analysis " + analysisID + " has library strategy " + tcgaRecord.library_strategy + " but is stored at " + fullPath);
                        continue;
                    }

                    if (tcgaRecord.tumorSample != (pathComponents[tumorNormalComponent] == "tumor"))
                    {
                        Console.WriteLine("Analysis " + analysisID + " in wrong tumor/normal directory " + fullPath);
                        continue;
                    }

                    if (tcgaRecord.disease_abbr != pathComponents[diseaseComponent])
                    {
                        Console.WriteLine("Sample in wrong tumor type.  Should be " + tcgaRecord.disease_abbr + " but it's in " + fullPath);
                        continue;
                    }

                    string correctMachine = tumorTypeToMachineMapping[new AnalysisType(tcgaRecord)].ToLower();
                    if (correctMachine != pathComponents[2])
                    {
                        Console.WriteLine("BAM stored on wrong machine " + fullPath + " expected to be on " + correctMachine);
                        continue;
                    }
                }

                if (tcgaRecord.bamFileName.ToLower() != pathComponents[pathLength - 1].ToLower())
                {
                    Console.WriteLine("Unexpected BAM file " + fullPath + " different from expected file name " + tcgaRecord.bamFileName);
                    continue;
                }
            }

        }

        const string realignPathname = @"f:\sequence\reads\tcga\realigns.txt";
        public static void LoadTCGARecordsForLocalRealigns(Dictionary<AnalysisID, TCGARecord> tcgaRecords, Dictionary<AnalysisID, StoredBAM> storedBAMs) 
        {
            if (!File.Exists(realignPathname))
            {
                return;
            }
            var realigns = File.ReadAllLines(realignPathname);

            //
            // Format is the list of fields from a TCGARecord separated by *s.
            //
            foreach (var realign in realigns)
            {
                var fields = realign.Split('*');
                var record = new TCGARecord();
                record.analysis_id = fields[0].ToLower();
                record.participant_id = fields[1].ToLower();
                record.bamFileName = fields[2];
                record.baiFileName = fields[2] + ".bai";
                record.lastModified = Convert.ToDateTime(fields[3]);
                record.uploadDate = Convert.ToDateTime(fields[4]);
                record.aliquot_id = fields[5].ToLower();
                record.refassemShortName = fields[6].ToLower();
                record.disease_abbr = fields[7].ToLower();
                record.totalFileSize = 0;   // Don't know for the locally generated ones; have to get this from storedBAMs
                record.sampleType = fields[8];
                record.library_strategy = fields[9];
                record.center_name = fields[10].ToLower();
                record.tumorSample = fields[11] == "tumor";
                record.localRealign = true;
                record.realignedFrom = fields[12];
                record.realignSource = tcgaRecords[record.realignedFrom];

                UseAnalysis(record.realignedFrom, null, storedBAMs);

                if (storedBAMs.ContainsKey(record.analysis_id))
                {
                    record.storedBAM = storedBAMs[record.analysis_id];
                }

                tcgaRecords.Add(record.analysis_id, record);
            }

        }

        public static string GenerateRealignmentFileString(TCGARecord tcgaRecord, string newRefassem)
        {
            string analysisID = Guid.NewGuid().ToString().ToLower();
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
        public static void GenerateRealignmentAnalyses(List<Experiment> experiments) 
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

        public static Dictionary<SampleID, List<TCGARecord>> BuildAnalysesBySample(Dictionary<AnalysisID, TCGARecord> tcgaRecords)
        {
            var analysesBySample = new Dictionary<SampleID, List<TCGARecord>>();

            foreach (var tcgaEntry in tcgaRecords)
            {
                TCGARecord entry = tcgaEntry.Value;

                if (!analysesBySample.ContainsKey(entry.aliquot_id))
                {
                    analysesBySample.Add(entry.aliquot_id, new List<TCGARecord>());
                }
                analysesBySample[entry.aliquot_id].Add(entry);
            }

            return analysesBySample;
        }

        public class MAFRecord
        {
            public string entire_maf_line;  // The raw line.

            public string Hugo_symbol;
            public string center;
            public string NcbiBuild;
            public string Chrom;
            public int Start_position;
            public int End_position;
            public string Variant_classification;
            public string Variant_type;
            public string Reference_allele;
            public string Tumor_seq_allele_1;
            public string Tumor_seq_allele_2;
            public string Tumor_sample_barcode;
            public string Matched_norm_sample_barcode;
            public SampleID Tumor_sample_uuid;
            public SampleID Matched_normal_sample_uuid;
            public Pathname MAFFileName;
        }

        static void AddMAFFileToParticipants(Pathname inputFilename, Dictionary<ParticipantID, Participant> participants, Dictionary<SampleID, ParticipantID> sampleToParticipantIDMap)
        {
            string[] inputLines = File.ReadAllLines(inputFilename);

            int tooShortLines = 0;

            // Do this in two phases.  First, parse all the MAF entries in the file into lists based on participant ID.  Then, add the lists to the participants.
            // Since none of the expriments we do care about mitochondrial mutations, we just filter them out here.
            //
            var mafs = new Dictionary<ParticipantID, List<MAFRecord>>();
 
            for (int i = 1; i < inputLines.Count(); i++)    // Line 0 is a header, skip it.
            {
                string[] fields = inputLines[i].Split('\t');
                if (fields.Count() < 37)
                {
                    tooShortLines++;
                    continue;
                }

                var center = fields[2];

                SampleID tumorSampleId = fields[32].ToLower();
                SampleID normalSampleID = fields[33].ToLower();

                if (fields[4].ToLower() == "mt" || fields[4].ToLower() == "chrm"  || fields[4].ToLower() == "m")
                {
                    continue;   // Just ignore mitochondrial mutations
                }

                if (!sampleToParticipantIDMap.ContainsKey(tumorSampleId) || !sampleToParticipantIDMap.ContainsKey(normalSampleID))
                {
                    //Console.WriteLine("Found MAF entry with no sample record in TCGA, tumor ID " + tumorSampleId + " normal ID " + normalSampleID);
                    // This seems to happen pretty often.  Just ignore the MAF line.
                    continue;
                }

                if (sampleToParticipantIDMap[tumorSampleId] != sampleToParticipantIDMap[normalSampleID])
                {
                    Console.WriteLine("MAF entry contains samples from different participants.");
                    return;
                }
                ParticipantID participantID = sampleToParticipantIDMap[tumorSampleId];

                MAFRecord mafRecord = new MAFRecord();
                mafRecord.entire_maf_line = inputLines[i];

                mafRecord.Hugo_symbol = fields[0];
                mafRecord.center = fields[2];
                mafRecord.NcbiBuild = fields[3];
                mafRecord.Chrom = fields[4];
                mafRecord.Start_position = Convert.ToInt32(fields[5]);
                mafRecord.End_position = Convert.ToInt32(fields[6]);
                mafRecord.Variant_classification = fields[8];
                mafRecord.Variant_type = fields[9];
                mafRecord.Reference_allele = fields[10];
                mafRecord.Tumor_seq_allele_1 = fields[11];
                mafRecord.Tumor_seq_allele_2 = fields[12];
                mafRecord.Tumor_sample_barcode = fields[15];
                mafRecord.Matched_norm_sample_barcode = fields[16];
                mafRecord.Tumor_sample_uuid = tumorSampleId;
                mafRecord.Matched_normal_sample_uuid = normalSampleID;
                mafRecord.MAFFileName = inputFilename;

                if (!mafs.ContainsKey(participantID))
                {
                    mafs.Add(participantID, new List<MAFRecord>());
                }

                mafs[participantID].Add(mafRecord);
            }

            //
            // Now go through what we've built up and add them to the participants.
            //
            foreach (var mafEntry in mafs)
            {
                ParticipantID participantID = mafEntry.Key;
                Participant participant = participants[participantID];
                participant.mafs.Add(mafEntry.Value);
            }

        }

        static Dictionary<SampleID, ParticipantID> CreateSampleToParticipantIDMap(Dictionary<AnalysisID, TCGARecord> tcgaEntries)
        {
            var sampleToParticipantIDMap = new Dictionary<SampleID, ParticipantID>();

            foreach (var entry in tcgaEntries)
            {
                string sampleID = entry.Value.aliquot_id;
                string participantID = entry.Value.participant_id;
                if (sampleToParticipantIDMap.ContainsKey(sampleID))
                {
                    if (sampleToParticipantIDMap[sampleID] != participantID)
                    {
                        Console.WriteLine("Multiple mapping for sampleID -> participant ID");
                    }
                }
                else 
                {
                    sampleToParticipantIDMap.Add(entry.Value.aliquot_id, entry.Value.participant_id);
                }
            }

            return sampleToParticipantIDMap;
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

        public class Sample
        {
            public SampleID sampleId;
            public bool isTumor;
            public ParticipantID participantID;
            public List<TCGARecord> RNA = new List<TCGARecord>();
            public List<TCGARecord> DNA = new List<TCGARecord>();
        }

        public class Participant
        {
            public string participantId;
            public List<List<MAFRecord>> mafs = new List<List<MAFRecord>>();
            public Dictionary<SampleID, Sample> tumorSamples = new Dictionary<SampleID, Sample>();   // Maps sampleID->sample for this participant
            public Dictionary<SampleID, Sample> normalSamples = new Dictionary<SampleID, Sample>();  // Maps sampleID->sample for this participant
            public Dictionary<string, List<TCGARecord>> tumorRNAByReference = new Dictionary<string, List<TCGARecord>>();
            public Dictionary<string, List<TCGARecord>> normalRNAByReference = new Dictionary<string, List<TCGARecord>>();
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

        public static Dictionary<string, Participant> BuildParticipantData(Dictionary<string, TCGARecord> tcgaEntries, out Dictionary<string, Sample> allSamples)
        {
            var participants = new Dictionary<string, Participant>();

            allSamples = new Dictionary<string, Sample>();

            foreach (var tcgaEntry in tcgaEntries)
            {
                string library_strategy = tcgaEntry.Value.library_strategy;
                string participant_id = tcgaEntry.Value.participant_id;

                Participant participant;
                if (!participants.ContainsKey(participant_id))
                {
                    participant = new Participant();
                    participant.participantId = participant_id;
                    participants.Add(participant_id, participant);
                }
                else
                {
                    participant = participants[participant_id];
                }

                string sample_id = tcgaEntry.Value.aliquot_id;
                Sample sample;
                Dictionary<string, Sample> samples;
                if (tcgaEntry.Value.tumorSample)
                {
                    samples = participant.tumorSamples;
                }
                else
                {
                    samples = participant.normalSamples;
                }

                if (!samples.ContainsKey(sample_id))
                {
                    sample = new Sample();
                    sample.sampleId = sample_id;
                    sample.isTumor = tcgaEntry.Value.tumorSample;
                    sample.participantID = participant_id;
                    samples.Add(sample_id, sample);
                    allSamples.Add(sample_id, sample);
                }
                else
                {
                    sample = samples[sample_id];
                    if (sample.isTumor != tcgaEntry.Value.tumorSample) {
                        Console.WriteLine("Sample " + sample_id + " appears to be both tumor and normal");
                        //
                        // Just skip it.
                        //
                        continue;
                    }
                }


                if (library_strategy == "WXS" || library_strategy == "WGS")
                {
                    sample.DNA.Add(tcgaEntry.Value);
                }
                else if (library_strategy == "RNA")
                {
                    sample.RNA.Add(tcgaEntry.Value);
                    Dictionary<string, List<TCGARecord>> rnaRecords;    // Either tumor or normal
                    if (tcgaEntry.Value.tumorSample)
                    {
                        rnaRecords = participant.tumorRNAByReference;
                    }
                    else
                    {
                        rnaRecords = participant.normalRNAByReference;
                    }

                    if (!rnaRecords.ContainsKey(tcgaEntry.Value.refassemShortName))
                    {
                        rnaRecords.Add(tcgaEntry.Value.refassemShortName, new List<TCGARecord>());
                    }

                    rnaRecords[tcgaEntry.Value.refassemShortName].Add(tcgaEntry.Value);
                } 
                else
                {
                    Console.WriteLine("Unknown library strategy (this should have been caught on import) " + library_strategy + " for sample " + sample_id);
                    continue;
                }
            }


            return participants;
        }

        public static void UseAnalysis(AnalysisID analysisID, List<AnalysisID> neededFiles, Dictionary<AnalysisID, StoredBAM> storedBAMs)
        {
            if (storedBAMs.ContainsKey(analysisID))
            {
                storedBAMs[analysisID].needed = true;
            }
            else if (neededFiles != null && !neededFiles.Contains(analysisID))
            {
                neededFiles.Add(analysisID);
            }
        }

        public class Experiment1Run
        {
            public MAFRecord mafRecord;
            public StoredBAM rnaBAM;
            public StoredBAM normalBAM;
            public TCGARecord rnaRecord;
            public TCGARecord normalRecord;
        }

        public class Experiment1Participant // The inputs for all the experiment 1 runs for a particular participant
        {
            public Experiment1Participant ()
            {
                for (int silent = 0; silent < 2; silent++)
                    for (int sex = 0; sex < 2; sex++)
                        for (int multipleMutations = 0; multipleMutations < 2; multipleMutations++)
                            runs[silent, sex, multipleMutations] = new List<Experiment1Run>();
            }
            public Participant participant;
            public List<Experiment1Run>[, ,] runs = new List<Experiment1Run>[2, 2, 2];    // Dimensions are slient/non-silent, autosome/sex (ignore MT), single/multiple mutations in this gene
        }


        public static string RefassemToIndex(string refassem)
        {
            string refassem_l = refassem.ToLower();
            if (refassem_l == "hg19")
            {
                return @"hg19-20";
            }

            if (refassem_l == "grch37-lite" || refassem_l == "GRCh37-lite_WUGSC_variant_1".ToLower() || refassem_l == "GRCh37-lite_WUGSC_variant_2".ToLower()
            ||refassem_l == "GRCh37-lite-+-HPV_Redux-build".ToLower())
            {
                return @"GRCh37-lite-20";
            }

            if (refassem_l == "NCBI36_BCCAGSC_variant".ToLower() || refassem_l == "NCBI36_BCM_variant".ToLower() || refassem_l == "NCBI36_WUGSC_variant".ToLower())
            {
                return @"NCBI36.54-20";
            }

            if (refassem_l == "grch37" || refassem_l == "GRCh37_BI_Variant".ToLower() || refassem_l == "hs37d5" /* this one's a guess*/)
            {
                return @"GRCh37-20";
            }

            if (refassem_l == "hg18" || refassem_l == "hg18_broad_variant" /*dubious*/)
            {
                return @"hg18-20";
            }

            if (refassem_l == "hg19_broad_variant")
            {
                return @"hg19_broad_variant-20";
            }

            if (refassem_l == "unaligned")
            {
                return @"hg19-20";  // Doesn't really matter, so just use the most generic one.
            }

            Console.WriteLine("RefassemToIndex: unknown refassem " + refassem);

            return "UnknownRefassem";
        }
        public static string RefassemToIndexPath(string refassem)
        {
            return @"\sequence\indices\" + RefassemToIndex(refassem); 

        }

        static void GenereateAdditionalTCGAMetadataGeneratingScript(Dictionary<AnalysisID, TCGARecord> tcgaRecords, Dictionary<AnalysisID, StoredBAM> storedBAMs, Dictionary<AnalysisType, Pathname> tumorToMachineMapping)
        {
            var writers = new Dictionary<Pathname, StreamWriter>();
            int nNeedingExtraction = 0;
            foreach (var entry in tcgaRecords)
            {
                TCGARecord record = entry.Value;

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

                AnalysisType analysisType = new AnalysisType(record);
                Pathname machine = tumorToMachineMapping[analysisType];

                if (!writers.ContainsKey(machine))
                {
                    writers.Add(machine, new StreamWriter(@"f:\temp\expression\extractTCGAAdditionalMetadata-" + machine + ".cmd"));
                    writers[machine].WriteLine(@"del \temp\analysisSummaries.txt");
                }


                writers[machine].WriteLine(@"SummarizeReads \sequence\indices\hg19-24 " + record.analysis_id + " " + storedBAMs[record.analysis_id].bamInfo.FullName + @" \temp\analysisSummaries.txt");
                nNeedingExtraction++;
            }

            foreach (var entry in writers)
            {
                entry.Value.Close();
            }

            Console.WriteLine("" + nNeedingExtraction + " analyses require extraction of additional metadata (read length).");
        }


        public class Experiment
        {
            public Participant participant;
            public TCGARecord TumorRNAAnalysis;
            public TCGARecord TumorDNAAnalysis;
            public TCGARecord NormalDNAAnalysis;
            public List<MAFRecord> maf;
            public bool normalNeedsRealignment = false;
            public bool tumorNeedsRealignment = false;
        }

        public static void RecordDownloadAndRealign(TCGARecord record, ref int nRealigns, ref long nRealignBytes, ref int nDownloads, ref long nDownloadBytes)
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

        public static List<Experiment> BuildExperiments(Dictionary<string, Participant> participants)
        {
            var experiments = new List<Experiment>();

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

            int nAMLbjb = 0;    // mostly just to set a breakpoint

            //
            // We try to do an experiment for each participant.
            //
            foreach (var entry in participants)
            {
                Participant participant = entry.Value;
                string participantID = entry.Key;

                if (participant.mafs.Count() == 0) {
                    //
                    // No mutation calls for this person, just give up.
                    //
                    continue;
                }

                Experiment bestExperiment = null;
                int bestExperimentScore = -1;

                bool hg18Like = participant.mafs[0][0].NcbiBuild == "36";
                //
                // Iterate over all the possible RNA samples.  To be a candidate, we must have the same hg18-like status as the maf.
                // Beyond that, select the best normal/tumor possible.  At the end of the loop, see if what we've got is better than the
                // best experiment so far.
                //
                int bestScore = -1;
                TCGARecord bestAnalysis = null;
                foreach (var RNASampleEntry in participant.tumorSamples)
                {
                    Sample RNASample = RNASampleEntry.Value;

                    foreach (var RNAtcgaRecord in RNASample.RNA)
                    {
                        if (RNAtcgaRecord.refassemShortName == "unaligned")
                        {
                            continue;
                        }
                        string rnaIndex = RefassemToIndex(RNAtcgaRecord.refassemShortName);

                        Experiment experiment = new Experiment();
                        experiment.participant = participant;
                        experiment.maf = participant.mafs[0];

                        if (RNAtcgaRecord.disease_abbr == "laml")
                        {
                            nAMLbjb++;
                        }

                        if (hg18_likeReferences.Contains(RNAtcgaRecord.refassemShortName) != hg18Like)
                        {
                            continue;   // Never use the wrong reference class
                        }

                        experiment.TumorRNAAnalysis = RNAtcgaRecord;
                        //
                        // Next, find a normal sample.  Any will do, but prefer in order
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
                            Sample normalSample = normalSampleEntry.Value;
                            foreach (var normalTCGARecord in normalSample.DNA)
                            {
                                int score = 0;
                                if (RefassemToIndex(normalTCGARecord.refassemShortName) == rnaIndex)
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
                        // Finally, find Tumor DNA.  It must have an upload date within a week of the RNA sample.
                        // 1) It is aligned against the same reference as the tumor RNA.
                        // 2) WGS Samples
                        // 3) We have it downloaded.
                        //
                        bestAnalysis = null;
                        bestScore = -1;
                        foreach (var tumorSampleEntry in participant.tumorSamples)
                        {
                            Sample tumorSample = tumorSampleEntry.Value;

                            foreach (var tumorTCGARecord in tumorSample.DNA)
                            {
                                int score = 0;

                                if (RefassemToIndex(tumorTCGARecord.refassemShortName) == rnaIndex)
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
                        experiment.normalNeedsRealignment = RefassemToIndex(experiment.NormalDNAAnalysis.refassemShortName) != rnaIndex || experiment.NormalDNAAnalysis.localRealign && experiment.NormalDNAAnalysis.storedBAM == null;
                        experiment.tumorNeedsRealignment = RefassemToIndex(experiment.TumorDNAAnalysis.refassemShortName) != rnaIndex || experiment.TumorDNAAnalysis.localRealign && experiment.TumorDNAAnalysis.storedBAM == null;

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
                Console.WriteLine(entry.Key + ":" + entry.Value + " (" + readyToGoByDisease[entry.Key] + " ready to go, " + readyToGoExceptForNormalByDisease[entry.Key] + " if you don't count normal)");
            }

            return experiments;
        }

        public static void AddDownloadToScript(TCGARecord tcgaRecord, Dictionary<string, StreamWriter> downloadScripts,  Dictionary<AnalysisType, Pathname> tumorToMachineMapping)
        {
            string machine = tumorToMachineMapping[new AnalysisType(tcgaRecord)];
            if (!downloadScripts.ContainsKey(machine))
            {
                downloadScripts.Add(machine, new StreamWriter(@"f:\temp\expression\loadFromTCGA-" + machine + ".cmd"));
            }

            downloadScripts[machine].WriteLine(@"cd /d d:\tcga\" + tcgaRecord.disease_abbr + (tcgaRecord.tumorSample ? @"\tumor\" : @"\normal\") + tcgaRecord.library_strategy);
            downloadScripts[machine].WriteLine(@"gtdownload -c c:\bolosky\ravip-cghub.key " + tcgaRecord.analysis_id);
        }

        public static void GenerateDownloadScripts(List<Experiment> experiments,  Dictionary<AnalysisType, Pathname> tumorToMachineMapping)
        {
            var downloadScripts = new Dictionary<string, StreamWriter>();
            foreach (var experiment in experiments)
            {
                if (experiment.TumorRNAAnalysis.storedBAM == null && !experiment.TumorRNAAnalysis.localRealign) {
                    AddDownloadToScript(experiment.TumorRNAAnalysis, downloadScripts, tumorToMachineMapping);
                }

                if (experiment.TumorDNAAnalysis.storedBAM == null && !experiment.TumorDNAAnalysis.localRealign)
                {
                    AddDownloadToScript(experiment.TumorDNAAnalysis, downloadScripts, tumorToMachineMapping);
                }

                if (experiment.NormalDNAAnalysis.storedBAM == null && !experiment.NormalDNAAnalysis.localRealign)
                {
                    AddDownloadToScript(experiment.NormalDNAAnalysis, downloadScripts, tumorToMachineMapping);
                }
            }

            foreach (var entry in downloadScripts)
            {
                entry.Value.Close();
            }
        }


        static public void GenerateUnneededLists(Dictionary<AnalysisID, StoredBAM> storedBAMs)
        {
            var perMachineLists = new Dictionary<string, List<StoredBAM>>();
            var perMachineDeleteBytes = new Dictionary<string, long>();

            long totalNeededBytes = 0;
            int nNeededBAMs = 0;
            int nUnneededBAMs = 0;
            long totalUnneededBytes = 0;

            foreach (var bamEntry in storedBAMs) {
                StoredBAM storedBAM = bamEntry.Value;

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
                        perMachineLists.Add(machine, new List<StoredBAM>());
                    }

                    perMachineDeleteBytes[machine] += size;
                    perMachineLists[machine].Add(storedBAM);
                }

            }

            Console.WriteLine("Need " + nNeededBAMs + " bams consisting of " + String.Format("{0:n0}", totalNeededBytes) + " bytes; do not need " + nUnneededBAMs + " consisting of " + String.Format("{0:n0}", totalUnneededBytes));

            var deleteLog = new StreamWriter(@"f:\temp\expression\deletable-files.txt");

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

        public static void AddSingleRealignmentToScript(TCGARecord record, StreamWriter script, Dictionary<AnalysisID, TCGARecord> tcgaRecords, Dictionary<AnalysisID, StoredBAM> storedBAMs,
            Dictionary<AnalysisType, Pathname> tumorToMachineMapping, bool local, bool bigMem)
        {
            AnalysisType analysisType = new AnalysisType(record);
            TCGARecord realignSource = tcgaRecords[record.realignedFrom];
            if (realignSource.storedBAM == null || !realignSource.hasAdditionalData)
            {
                //
                // First we have to download the source and get its additional data (read length, etc.), then we can realign it.
                //
                return;
            }

            string tumorOrNormal = (record.tumorSample ? "tumor" : "normal");
            string outputFileName = record.analysis_id + "-SNAP-realigned-" + realignSource.analysis_id + "-" + record.disease_abbr + "-" + tumorOrNormal + ".bam";
            string destinationDirectory;

            if (local)
            {
                destinationDirectory = @"d:\tcga\" + record.disease_abbr + @"\" + tumorOrNormal + @"\" + record.library_strategy + @"\" + record.analysis_id + @"\";
                outputFileName = destinationDirectory + outputFileName;
                script.WriteLine(@"md " + destinationDirectory);
            } else {
                destinationDirectory = @"\\" + tumorToMachineMapping[analysisType] + @"\d$\tcga\" + record.disease_abbr + @"\" + tumorOrNormal + @"\" + record.library_strategy + @"\" + record.analysis_id + @"\";
            }

            UseAnalysis(record.realignedFrom, null, storedBAMs);

            script.Write("snap ");
            if (realignSource.anyPaired)
            {
                script.Write(@"paired \sequence\indices\");
            }
            else
            {
                script.Write(@"single \sequence\indices\");
            }

            script.Write(record.refassemShortName + "-");
            int seedLength;
            if (realignSource.meanGoodBases * 3 < 16)
            {
                seedLength = 16;
            }
            else if (realignSource.meanGoodBases * 3 > 24)
            {
                seedLength = 24;
            }
            else
            {
                seedLength = (int)(realignSource.meanGoodBases * 3 + .99);
            }

            script.Write("" + seedLength + " " + realignSource.storedBAM.bamInfo.FullName + " -o " + outputFileName);
            if (realignSource.anyPaired)
            {
                script.Write(" -s 0 1000");
            }
            if (bigMem)
            {
                script.Write(" -so -sm 100 -xf 5");
            }
            else
            {
                script.Write(" -so -sm 5 -di -xf 1.5");
            }
            script.WriteLine(" -lp -map -pc -mrl " + Math.Min(50, seedLength * 2));

            if (!local)
            {
                script.WriteLine(@"md " + destinationDirectory);
                script.WriteLine("copy " + outputFileName + "* " + destinationDirectory);
                script.WriteLine(@"del " + outputFileName + "*");
            }
        }

        public static void GenerateRealignmentScript(List<Experiment> experiments, Dictionary<AnalysisID, TCGARecord> tcgaRecords, Dictionary<AnalysisID, StoredBAM> storedBAMs, Dictionary<AnalysisType, Pathname> tumorToMachineMapping)
        {
            //
            // Find any experiments that rely on a realignment, and generate the SNAP command to do the realignment.
            //
            var realignNormalScript = new StreamWriter(@"f:\temp\expression\realignNormal.cmd");
            var realignTumorScript = new StreamWriter(@"f:\temp\expression\realignTumor.cmd");
            var realignBigMemScript = new StreamWriter(@"f:\temp\expression\realignBigMem.cmd");
            var perMachineScripts = new Dictionary<Pathname, StreamWriter>();

            const long limitForBigMem = (long)50 * 1024 * 1024 * 1024;
            const long limitFor48 = (long)25 * 1024 * 1024 * 1024;

            foreach (var experiment in experiments)
            {
                List<TCGARecord> analyses = new List<TCGARecord>();

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
                        script = realignBigMemScript;
                        local = false;
                        bigMem = true;
                    }
                    else
                    {
                        var analysisType = new AnalysisType(analysis);
                        int localMem = Machines[tumorToMachineMapping[analysisType]].memoryInGB;

                        bigMem = false;

                        if (localMem < 48 || localMem == 48 && analysis.realignSource.storedBAM.totalSize >= limitFor48)
                        {
                            local = false;
                            if (analysis.tumorSample)
                            {
                                script = realignTumorScript;
                            }
                            else
                            {
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

                    AddSingleRealignmentToScript(analysis, script, tcgaRecords, storedBAMs, tumorToMachineMapping, local, bigMem);
                }
             }

            realignNormalScript.Close();
            realignTumorScript.Close();
            realignBigMemScript.Close();

            foreach (var entry in perMachineScripts)
            {
                entry.Value.Close();
            }
        }

        static public void GenerateExtractionScripts(List<Experiment> experiments)
        {
            StreamWriter extractScript = new StreamWriter(@"f:\temp\expression\extract_genetic_data.cmd");
            StreamWriter extractRNAInput = new StreamWriter(@"f:\temp\expression\extractRNAInput.txt");
            StreamWriter extractTumorDNAInput = new StreamWriter(@"f:\temp\expression\extractTumorDNAInput.txt");
            foreach (var experiment in experiments)
            {
                foreach (MAFRecord mafRecord in experiment.maf)
                {
                    if (mafRecord.Chrom == "MT" || mafRecord.Chrom == "M")
                    {
                        //
                        // Just totally ignore mitochondrial mutations.
                        //
                        continue;
                    }
                }
            }
        }
  
        static void Main(string[] args)
        {
            hg18_likeReferences.Add("NCBI36_BCCAGSC_variant".ToLower());
            hg18_likeReferences.Add("NCBI36_BCM_variant".ToLower());
            hg18_likeReferences.Add("NCBI36_WUGSC_variant".ToLower());
            hg18_likeReferences.Add("hg18");
            hg18_likeReferences.Add("HG18_Broad_variant".ToLower());

            LibraryStrategies.Add("wgs");
            LibraryStrategies.Add("wxs");
            LibraryStrategies.Add("rna");

            InitializeMachines();

            var tumorToMachineMapping = GenerateTumorToMachineMapping();
            Dictionary<AnalysisID, StoredBAM> storedBAMs = LoadStoredBAMs(tumorToMachineMapping);
            var tcgaRecords = LoadTCGARecords(storedBAMs);
            LoadTCGARecordsForLocalRealigns(tcgaRecords, storedBAMs);
            LoadTCGAAdditionalMetadata(tcgaRecords);
            GenereateAdditionalTCGAMetadataGeneratingScript(tcgaRecords, storedBAMs, tumorToMachineMapping);
            VerifyStoredBAMPaths(storedBAMs, tcgaRecords, tumorToMachineMapping);

            // off while it's running. GenereateAdditionalTCGAMetadataGeneratingScript(tcgaRecords, storedBAMs);
            
            var sampleToParticipantIDMap = CreateSampleToParticipantIDMap(tcgaRecords);
            var centerNameToDomainNameMap = GenerateCenterNameToDomainName();
            var domainNameToCenterNameMap = InvertStringDictionary(centerNameToDomainNameMap);

            //
            // Now build up participant data.
            //
            Dictionary<string, Sample> allSamples;
            var participants = BuildParticipantData(tcgaRecords, out allSamples);

            //
            // Load all the MAF files.
            //
            var mutations = new Dictionary<string, Dictionary<string, List<MAFRecord>>>();
            foreach (var file in Directory.GetFiles(@"f:\sequence\Reads\tcga\mafs\"))
            {
                AddMAFFileToParticipants(file, participants, sampleToParticipantIDMap);
            }

            var experiments = BuildExperiments(participants);
            GenerateRealignmentScript(experiments, tcgaRecords, storedBAMs, tumorToMachineMapping);
            GenerateUnneededLists(storedBAMs);
            GenerateRealignmentAnalyses(experiments);
            GenerateDownloadScripts(experiments, tumorToMachineMapping);

            GenerateExtractionScripts(experiments);

        }
    }
}
