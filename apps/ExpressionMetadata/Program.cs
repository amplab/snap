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
        static string baseDirectory = @"f:\temp\expression\";
        static string isoformDirectory = baseDirectory + @"isoform-files\";
        static List<string> hg18_likeReferences = new List<string>();

        static string ShareFromPathname(string pathname)
        {
            //
            // A pathname is of the form \\computer\share\...
            // Get up to share, discard the rest.
            //
            string[] fields = pathname.Split('\\');
            if (fields.Count() < 4 || fields[0] != "" || fields[1] != "")
            {
                Console.WriteLine("ShareFromPathname: can't parse pathname " + pathname);
                return "";
            }

            string share = "";
            for (int i = 1; i < 4; i++) // Skip 0 because we know it's the empty string, and it makes the \ part easier.
            {
                share += @"\" + fields[i];
            }

            return share;
        }

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
        static List<string> TumorAndNormal = new List<string>();

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
            Machine.AddMachine("msr-srs-0", 8, 36);
            Machine.AddMachine("msr-srs-1", 8, 36);
            Machine.AddMachine("msr-srs-2", 8, 36);
            Machine.AddMachine("msr-srs-3", 8, 36);
            Machine.AddMachine("msr-srs-4", 8, 36);
            Machine.AddMachine("msr-srs-5", 16, 36);
            Machine.AddMachine("msr-srs-6", 16, 36);
            Machine.AddMachine("msr-srs-7", 48, 48);
            Machine.AddMachine("msr-srs-8", 48, 36);
            Machine.AddMachine("fds-k25-1", 48, 6);
            Machine.AddMachine("fds-k25-2", 48, 6);
            Machine.AddMachine("fds-k25-3", 48, 6);
            Machine.AddMachine("fds-k25-5", 48, 6);
            Machine.AddMachine("fds-k25-7", 48, 6);
            Machine.AddMachine("fds-k25-8", 48, 6);
            Machine.AddMachine("fds-k25-9", 48, 6);
            Machine.AddMachine("fds-k25-14", 48, 6);
            Machine.AddMachine("fds-k25-15", 48, 6);
            Machine.AddMachine("fds-k25-18", 48, 6);
            Machine.AddMachine("fds-k25-20", 48, 6);
            Machine.AddMachine("msr-genomics-0", 256, 240, true);
            Machine.AddMachine("msr-genomics-1", 256, 240, true);
        }

        public static Dictionary<string, Machine> Machines = new Dictionary<string, Machine>();

        public static Dictionary<AnalysisType, Pathname> GenerateTumorToMachineMapping()
        {
            var tumorToMachineMapping = new Dictionary<AnalysisType, Pathname>();

            AddEntireTumorToMachine(tumorToMachineMapping, "acc",  "fds-k25-2");
            AddEntireTumorToMachine(tumorToMachineMapping, "blca", "msr-srs-7");
            AddEntireTumorToMachine(tumorToMachineMapping, "chol", "msr-srs-6");
            AddEntireTumorToMachine(tumorToMachineMapping, "coad", "msr-srs-6");
            AddEntireTumorToMachine(tumorToMachineMapping, "dlbc", "msr-srs-7");
            AddEntireTumorToMachine(tumorToMachineMapping, "gbm",  "msr-srs-4");
            AddEntireTumorToMachine(tumorToMachineMapping, "hnsc", "msr-srs-0");
            AddEntireTumorToMachine(tumorToMachineMapping, "kich", "fds-k25-9");
            AddEntireTumorToMachine(tumorToMachineMapping, "kirc", "msr-srs-3");
            AddEntireTumorToMachine(tumorToMachineMapping, "kirp", "msr-srs-6");
            AddEntireTumorToMachine(tumorToMachineMapping, "laml", "msr-srs-8");
            AddEntireTumorToMachine(tumorToMachineMapping, "luad", "msr-srs-6");
            AddEntireTumorToMachine(tumorToMachineMapping, "meso", "msr-srs-3");
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

            tumorToMachineMapping.Add(new AnalysisType("brca", "rna", true), "msr-srs-7");
            tumorToMachineMapping.Add(new AnalysisType("brca", "wgs", true), "msr-srs-7");
            tumorToMachineMapping.Add(new AnalysisType("brca", "wxs", true), "msr-srs-1");
            tumorToMachineMapping.Add(new AnalysisType("brca", "rna", false), "msr-srs-1");
            tumorToMachineMapping.Add(new AnalysisType("brca", "wgs", false), "msr-srs-1");
            tumorToMachineMapping.Add(new AnalysisType("brca", "wxs", false), "msr-srs-1");

            tumorToMachineMapping.Add(new AnalysisType("cesc", "rna", true), "msr-srs-7");
            tumorToMachineMapping.Add(new AnalysisType("cesc", "wgs", true), "msr-srs-7");
            tumorToMachineMapping.Add(new AnalysisType("cesc", "wxs", true), "msr-srs-7");
            tumorToMachineMapping.Add(new AnalysisType("cesc", "rna", false), "fds-k25-15");
            tumorToMachineMapping.Add(new AnalysisType("cesc", "wgs", false), "fds-k25-15");
            tumorToMachineMapping.Add(new AnalysisType("cesc", "wxs", false), "fds-k25-15");

            tumorToMachineMapping.Add(new AnalysisType("esca", "rna", true), "fds-k25-1");
            tumorToMachineMapping.Add(new AnalysisType("esca", "wgs", true), "fds-k25-1");
            tumorToMachineMapping.Add(new AnalysisType("esca", "wxs", true), "fds-k25-1");
            tumorToMachineMapping.Add(new AnalysisType("esca", "rna", false), "fds-k25-1");
            tumorToMachineMapping.Add(new AnalysisType("esca", "wgs", false), "fds-k25-1");
            tumorToMachineMapping.Add(new AnalysisType("esca", "wxs", false), "msr-srs-7");

            tumorToMachineMapping.Add(new AnalysisType("lgg", "rna", true), "msr-srs-3");
            tumorToMachineMapping.Add(new AnalysisType("lgg", "wgs", true), "msr-srs-3");
            tumorToMachineMapping.Add(new AnalysisType("lgg", "wxs", true), "msr-srs-3");
            tumorToMachineMapping.Add(new AnalysisType("lgg", "rna", false), "fds-k25-14");
            tumorToMachineMapping.Add(new AnalysisType("lgg", "wgs", false), "fds-k25-14");
            tumorToMachineMapping.Add(new AnalysisType("lgg", "wxs", false), "fds-k25-3");

            tumorToMachineMapping.Add(new AnalysisType("lihc", "rna", true), "msr-srs-2");
            tumorToMachineMapping.Add(new AnalysisType("lihc", "wgs", true), "msr-srs-2");
            tumorToMachineMapping.Add(new AnalysisType("lihc", "wxs", true), "msr-srs-2");
            tumorToMachineMapping.Add(new AnalysisType("lihc", "rna", false), "msr-srs-2");
            tumorToMachineMapping.Add(new AnalysisType("lihc", "wgs", false), "fds-k25-20");
            tumorToMachineMapping.Add(new AnalysisType("lihc", "wxs", false), "msr-srs-2");

            tumorToMachineMapping.Add(new AnalysisType("skcm", "rna", true), "fds-k25-8");
            tumorToMachineMapping.Add(new AnalysisType("skcm", "wgs", true), "fds-k25-8");
            tumorToMachineMapping.Add(new AnalysisType("skcm", "wxs", true), "fds-k25-8");
            tumorToMachineMapping.Add(new AnalysisType("skcm", "rna", false), "msr-srs-5");
            tumorToMachineMapping.Add(new AnalysisType("skcm", "wgs", false), "msr-srs-5");
            tumorToMachineMapping.Add(new AnalysisType("skcm", "wxs", false), "msr-srs-5");

            tumorToMachineMapping.Add(new AnalysisType("lusc", "rna", true), "fds-k25-18");
            tumorToMachineMapping.Add(new AnalysisType("lusc", "wgs", true), "fds-k25-18");
            tumorToMachineMapping.Add(new AnalysisType("lusc", "wxs", true), "fds-k25-18");
            tumorToMachineMapping.Add(new AnalysisType("lusc", "rna", false), "msr-srs-7");
            tumorToMachineMapping.Add(new AnalysisType("lusc", "wgs", false), "msr-srs-7");
            tumorToMachineMapping.Add(new AnalysisType("lusc", "wxs", false), "msr-srs-7");

             return tumorToMachineMapping;
        }

        public class StoredBAM
        {
            public AnalysisID analysisID;
            public Pathname directoryName;
            public string machineName;
            public char driveLetter;
            public FileInfo bamInfo = null;
            public FileInfo baiInfo = null;
            public FileInfo vcfInfo = null;
            public FileInfo allCountInfo = null;
            public long totalSize;
            public bool needed = false; // Do we need this as input for some experiment?
            public bool isRealingned = false;   // Was this realined, or is it straight from TCGA
            public string bamHash = null;  // MD5 of the BAM file
            public string baiHash = null;  // MD5 of the BAI file
            public string isoformCount = null;  // File with the count of isoform references
        }

        static public string GetDirectoryPathFromFullyQualifiedFilename(string filename) 
        {
            if (filename.Count() < 1 || filename[0] != '\\') {
                Console.WriteLine("GetDirectoryPathFromFullyQualifiedFilename -- invalid input " + filename);
                return "";
            }

            string[] pathnameComponents = filename.Split('\\');
            if (pathnameComponents.Count() < 4 || pathnameComponents[0] != "" || pathnameComponents[1] != "" || pathnameComponents[3].Count() != 2 || pathnameComponents[3].Substring(1) != "$")
            {
                Console.WriteLine("GetDirectoryPathFromFullyQualifiedFilename: expected \\ pathname, got " + filename);
                return "";
            }

            string dirName = @"";

            for (int i = 0; i < pathnameComponents.Count() - 1; i++)
            {
                dirName += pathnameComponents[i] + @"\";
            }


            return dirName;
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
        static StreamWriter AllIsoformsFile;
        static void LoadStoredBAMsForDirectory(Pathname directory, Dictionary<AnalysisID, StoredBAM> storedBAMs, Dictionary<Pathname, StreamWriter> hashScripts) 
        {
            if (!Directory.Exists(directory))
            {
                Console.WriteLine("Unable to find expected tcga data directory " + directory);
                return;
            }

            foreach (var subdir in Directory.GetDirectories(directory))
            {
                var pathnameComponents = subdir.Split('\\');
                if (pathnameComponents.Count() < 4)
                {
                    Console.WriteLine("LoadStoredBAMsForDirectory: bizarre pathname " + subdir);
                    continue;
                }

                AnalysisID analysisID = pathnameComponents[pathnameComponents.Count() - 1].ToLower();
                if (analysisID.Count() == 36)  // The length of a GUID string
                {
                    var storedBAM = new StoredBAM();
                    storedBAM.analysisID = analysisID;
                    storedBAM.machineName = pathnameComponents[2].ToLower();
                    storedBAM.driveLetter = char.ToLower(pathnameComponents[3][0]);

                    bool hasTarFile = false;    // Some RNA analyses include a tarred FASTQ.  If we see one, don't generate a warning for not having BAM files.

                    foreach (var file in Directory.GetFiles(subdir))
                    {
                        if (file.Count() > 4 && file.Substring(file.Count() - 4) == ".bam")
                        {
                            if (-1 != file.IndexOf("HOLD_QC"))
                            {
                                Console.WriteLine("File " + file + " appears to be QC pending; you may want to exclude it.");
                            }
                            if (null != storedBAM.bamInfo)
                            {
                                Console.WriteLine("Saw multiple BAM files in the same analysis directory " + file + " and " + storedBAM.bamInfo.FullName);
                            }
                            storedBAM.bamInfo = new FileInfo(file);
                            string md5Pathname = file + ".md5";
                            if (File.Exists(md5Pathname))
                            {
                                try
                                {
                                    storedBAM.bamHash = File.ReadAllText(md5Pathname);
                                    if (storedBAM.bamHash == "")
                                    {
                                        storedBAM.bamHash = null;   // Probably it got interrupted while being computed, just ignore it and it'll get run again.
                                    }
                                }
                                catch (IOException e)
                                {
                                    //
                                    // It's probably a sharing violation because we're computing hashes.  Just warn and go on.
                                    //
                                    Console.WriteLine("IO Exception reading m5 file " + md5Pathname + ".  If you're not computing that hash now, maybe something's wrong; error: " + e.Message);
                                }
                            }
                            else
                            {
                                GenerateHashCommandForFile(file, hashScripts);
                            }
                        }
                        else if (file.Count() > 4 && file.Substring(file.Count() - 4) == ".bai")
                        {
                            if (null != storedBAM.baiInfo)
                            {
                                Console.WriteLine("Saw multiple BAI files in the same analysis directory " + file + " and " + storedBAM.baiInfo.FullName);
                            }
                            storedBAM.baiInfo = new FileInfo(file);
                            string md5Pathname = file + ".md5";

                            if (File.Exists(md5Pathname))
                            {
                                try {
                                    storedBAM.baiHash = File.ReadAllText(md5Pathname);
                                    if (storedBAM.baiHash == "")
                                    {
                                        storedBAM.baiHash = null;   // Probably it got interrupted while being computed, just ignore it and it'll get run again.
                                    }
                                }
                                catch (IOException e)
                                {
                                    //
                                    // It's probably a sharing violation because we're computing hashes.  Just warn and go on.
                                    //
                                    Console.WriteLine("IO Exception reading m5 file " + md5Pathname + ".  If you're not computing that hash now, maybe something's wrong; error: " + e.Message);
                                }
                            }
                            else
                            {
                                GenerateHashCommandForFile(file, hashScripts);
                            }
                        }
                        else if (file.Count() > 4 && file.Substring(file.Count() - 4) == ".tar")
                        {
                            hasTarFile = true;
                        } else if (file.Count() > 3 && file.Substring(file.Count() - 4) == ".md5") {
                            // Do nothing; checksum files are expected.
                        } else if (file == subdir + @"\" + analysisID + @".vcf") {
                            storedBAM.vcfInfo = new FileInfo(file);
                        }
                        else if (file.Count() > 12 && file.Substring(file.Count() - 13).ToLower() == "-isoforms.txt")
                        {
                            string expectedFilename = analysisID + "-isoforms.txt";
                            if (file.Count() < expectedFilename.Count() || file.Substring(file.Count() - expectedFilename.Count()) != expectedFilename)
                            {
                                Console.WriteLine("Incorrect isoform file " + file);
                            } else {
                                storedBAM.isoformCount = file;
                                AllIsoformsFile.WriteLine(file);
                            }
                        }
                        else if (file.Count() > 12 && file.Substring(file.Count() - 12).ToLower() == ".allcount.gz")
                        {
                            string expectedFilename = analysisID + ".allcount.gz";
                            if (file.Count() < expectedFilename.Count() || file.Substring(file.Count() - expectedFilename.Count()) != expectedFilename)
                            {
                                Console.WriteLine("Incorrect allcount file " + file);
                            }
                            else
                            {
                                storedBAM.allCountInfo = new FileInfo(file);
                                if (storedBAM.allCountInfo.Length < 1 * 1024 * 1024)
                                {
                                    Console.WriteLine("Unusually small allCount file " + file + ", " + storedBAM.allCountInfo.Length);
                                }
                            }
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

        public static void LoadStoredBamsForMachine(Pathname directory, Dictionary<AnalysisID, StoredBAM> storedBAMs, Dictionary<Pathname, StreamWriter> hashScripts)
        {
            //
            // Load BAMs from a machine that offloads storage from the main machines.  These machines are allowed to store data
            // in places that's other than its primary storage location.  We just work though the directory hierarchy and
            // load each lowest level directory individually.  "directory" should point at the directory containing the
            // disease subdirectories, i.e., \\msr-genomics-0\d$\tcga.
            //

            if (!Directory.Exists(directory))
            {
                Console.WriteLine("Unable to find expected secondary machine directory " + directory);
                return;
            }

            Console.WriteLine("Loading bams for " + directory);

            foreach (var diseaseDir in Directory.GetDirectories(directory)) 
            {
                foreach (var tumorOrNormal in TumorAndNormal)
                {
                    if (!Directory.Exists(diseaseDir + @"\" + tumorOrNormal))
                    {
                        continue;
                    }

                    foreach (var libraryStrategy in LibraryStrategies)
                    {
                        var sampleDir = diseaseDir + @"\" + tumorOrNormal + @"\" + libraryStrategy;
                        if (Directory.Exists(sampleDir))
                        {
                            LoadStoredBAMsForDirectory(sampleDir, storedBAMs, hashScripts);
                        }
                    }
                }
            }
        }

        public static Dictionary<AnalysisID, StoredBAM> LoadStoredBAMs(Dictionary<AnalysisType, Pathname> tumorToMachineMapping)
        {
            var storedBAMs = new Dictionary<AnalysisID, StoredBAM>();
            var hashScripts = new Dictionary<Pathname, StreamWriter>();

            //
            // Directory structure is \\msr-srs-%n\d$\tcga\{rna,wgs,wxs}\{tumor, normal}\disease_abbr\analysis_id\*.{bam,bai}.  We need to call LoadStoredBAMsForDirectory on each
            // of the disease_abbr directories.  
            //
            foreach (var machine in Machines)
            {
                LoadStoredBamsForMachine(@"\\" + machine.Value.name + @"\d$\tcga", storedBAMs, hashScripts);
            }

            LoadStoredBamsForMachine(@"\\msr-genomics-0\e$\tcga", storedBAMs, hashScripts);
            LoadStoredBamsForMachine(@"\\msr-genomics-1\e$\tcga", storedBAMs, hashScripts);
            LoadStoredBamsForMachine(@"\\bolosky\f$\tcga", storedBAMs, hashScripts);

            foreach (var entry in hashScripts)
            {
                entry.Value.Close();
            }

            return storedBAMs;
        }

        public class TCGARecord
        {
            public AnalysisID analysis_id;
            public ParticipantID participant_id;
            public Pathname bamFileName;
            public string bamHash = null;
            public Pathname baiFileName;
            public string baiHash = null;
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

        public static Dictionary<AnalysisID, TCGARecord> LoadTCGARecords(Dictionary<AnalysisID, StoredBAM> storedBAMs, List<AnalysisID> excludedAnalyses)
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

                if (excludedAnalyses.Contains(record.analysis_id)) {
                    continue;   // Just pretend we didn't see it.
                }

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
                        record.bamHash = file.Element("checksum").Value;
                    }
                    else if (fileName.Length > 4 && ".bai" == fileName.Substring(fileName.Length - 4))
                    {
                        record.baiFileName = fileName;
                        record.baiHash = file.Element("checksum").Value;
                    }
                    record.totalFileSize += Convert.ToInt64(file.Element("filesize").Value);
                }

                if (record.bamFileName == "" || record.baiFileName == "")
                {
                    //
                    // This might be a broken analysis with a name like foo.bam_HOLD_QC_PENDING or something.  Ignore it.
                    //
                    continue;
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
            if (!File.Exists(@"f:\sequence\reads\tcga\tcgaAdditionalMetadata.txt"))
            {
                return;
            }
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
                StoredBAM storedBAM = bamEntry.Value;
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

                TCGARecord tcgaRecord = tcgaRecords[analysisID];

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

                string correctMachine = tumorTypeToMachineMapping[new AnalysisType(tcgaRecord)].ToLower();
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
                if (!tcgaRecords.ContainsKey(record.realignedFrom))
                {
                    Console.WriteLine("Realignment analysis " + record.analysis_id + " is realigned from an unknown source; most likely something we excluded after we generated this realignment.  Maybe you want to delete it from realigns.txt");
                    continue;
                }
                record.realignSource = tcgaRecords[record.realignedFrom];

                UseAnalysis(record.realignedFrom, null, storedBAMs);

                if (storedBAMs.ContainsKey(record.analysis_id))
                {
                    record.storedBAM = storedBAMs[record.analysis_id];
                    if (record.realignSource.storedBAM != null && record.realignSource.storedBAM.bamInfo.Length  > record.storedBAM.bamInfo.Length * 2)
                    {
                        Console.WriteLine("Suspiciously small realigned file " + record.storedBAM.bamInfo.FullName + " size " + record.storedBAM.bamInfo.Length + 
                            " while source " + record.realignSource.storedBAM.bamInfo.FullName + " is " + record.realignSource.storedBAM.bamInfo.Length);
                    }
                }

                tcgaRecords.Add(record.analysis_id, record);
            }

        }

        public static string GenerateRealignmentFileString(string analysisID, TCGARecord tcgaRecord, string newRefassem)
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

        public static void RebuildRealignmentAnalysesFromExperimentsFile(Dictionary<AnalysisID, TCGARecord> tcgaRecords)
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
        public static string GenerateRealignmentFileString(TCGARecord tcgaRecord, string newRefassem) {
            string analysisID = Guid.NewGuid().ToString().ToLower();
            return GenerateRealignmentFileString(analysisID, tcgaRecord, newRefassem);
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

            public bool hasDNACounts = false;
            public int DNALineNumber;
            public int nDNAMatchingNormal;
            public int nDNAMatchingTumor;
            public int nDNAMatchingNeither;
            public int nDNAMatchingBoth;

            public bool hasRNACounts = false;
            public int RNALineNumber;
            public int nRNAMatchingNormal;
            public int nRNAMatchingTumor;
            public int nRNAMatchingNeither;
            public int nRNAMatchingBoth;
        }

        static void AddMAFFileToParticipants(Pathname inputFilename, Dictionary<ParticipantID, Participant> participants, Dictionary<SampleID, ParticipantID> sampleToParticipantIDMap)
        {
            string[] inputLines = File.ReadAllLines(inputFilename);

            int tooShortLines = 0;

            // Do this in two phases.  First, parse all the MAF entries in the file into lists based on participant ID.  Then, add the lists to the participants.
            // Since none of the expriments we do care about mitochondrial mutations, we just filter them out here.
            //
            var mafs = new Dictionary<ParticipantID, Dictionary<string, List<MAFRecord>>>();
 
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
                mafRecord.NcbiBuild = fields[3].ToLower();
                mafRecord.Chrom = fields[4].ToLower();
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
                    mafs.Add(participantID, new Dictionary<string,List<MAFRecord>>());
                }

                string key = mafRecord.Tumor_sample_uuid + mafRecord.Matched_normal_sample_uuid;
                if (!mafs[participantID].ContainsKey(key))
                {
                    mafs[participantID].Add(key, new List<MAFRecord>());
                }

                mafs[participantID][key].Add(mafRecord);
            }

            //
            // Now go through what we've built up and add them to the participants.
            //
            foreach (var participantEntry in mafs)
            {
                ParticipantID participantID = participantEntry.Key;
                Participant participant = participants[participantID];
                foreach (var mafEntry in participantEntry.Value)
                {
                    participant.mafs.Add(mafEntry.Value);
                }
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
            public Dictionary<string, List<MAFRecord>> mafSets = new Dictionary<string, List<MAFRecord>>();
            public Dictionary<string, List<MAFRecord>> mafsByGene = new Dictionary<TumorType, List<MAFRecord>>();
            public Dictionary<SampleID, Sample> tumorSamples = new Dictionary<SampleID, Sample>();   // Maps sampleID->sample for this participant
            public Dictionary<SampleID, Sample> normalSamples = new Dictionary<SampleID, Sample>();  // Maps sampleID->sample for this participant
            public Dictionary<string, List<TCGARecord>> tumorRNAByReference = new Dictionary<string, List<TCGARecord>>();
            public Dictionary<string, List<TCGARecord>> normalRNAByReference = new Dictionary<string, List<TCGARecord>>();
            public bool metadataPresent = false;
            public string gender = "";
            public string race = "";
            public string ethnicity = "";
            public string historyOfOtherMalignancy = "";
            public string tumorStatus = "";
            public string vitalStatus = "";
            public string occupationPrimary = "";
            public string treatmentOutcomeFirstCourse = "";
            public int daysToBirth = 0; // This is negative
            public int daysToDeath = -1;    // -1 means not applicable
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

            //
            // Now add in the per-participant clinical metadata.
            //
            foreach (var clinical_file in Directory.EnumerateFiles(@"f:\sequence\reads\tcga\clinical"))
            {
                ProcessClinicalFile(participants, clinical_file);
            }

            //
            // Print out some stats.
            //
            int nMale = 0; 
            int nFemale = 0;
            int nAlive = 0;
            int nDead = 0;
            var raceCounts = new Dictionary<string, int>();

            foreach (var entry in participants)
            {
                Participant participant = entry.Value;
                if (participant.gender == "male")
                {
                    nMale++;
                }
                else if (participant.gender == "female")
                {
                    nFemale++;
                }

                if (participant.vitalStatus == "alive")
                {
                    nAlive++;
                }
                else if (participant.vitalStatus == "dead") 
                {
                    nDead++;
                }

                if (!raceCounts.ContainsKey(participant.race))
                {
                    raceCounts.Add(participant.race, 0);
                }
                raceCounts[participant.race]++;
            }

            Console.WriteLine("Found " + participants.Count() + " participants, " + participants.Where(x => x.Value.metadataPresent).Count() + " with metadata, of whom " + nMale + " are male, " + nFemale + " are female, " + nAlive + " are alive and " + nDead + " are dead.");

            Console.WriteLine("Races:");
            foreach (var entry in raceCounts)
            {
                Console.WriteLine(entry.Key + "\t" + entry.Value);
            }

            return participants;
        }

        public static void ProcessClinicalFile(Dictionary<ParticipantID, Participant> participants, string filename) 
        {

            //
            // These files don't all have the same headers.  Scan through the first line to find the column numbers we care about.
            //
            var fieldExtractorGenerators = new Dictionary<string, Func<int, Action<Participant, string []>>> ();
            fieldExtractorGenerators.Add("gender", i => ((p, fields) => p.gender = fields[i].ToLower()));
            fieldExtractorGenerators.Add("race", i => ((p, fields) => p.race = fields[i].ToLower()));
            fieldExtractorGenerators.Add("ethnicity", i => ((p, fields) => p.ethnicity = fields[i].ToLower()));
            fieldExtractorGenerators.Add("history_other_malignancy", i=> ((p, fields) => p.historyOfOtherMalignancy = fields[i].ToLower()));
            fieldExtractorGenerators.Add("tumor_status", i => ((p, fields) => p.tumorStatus = fields[i].ToLower()));
            fieldExtractorGenerators.Add("vital_status", i=> ((p, fields) => p.vitalStatus = fields[i].ToLower()));
            fieldExtractorGenerators.Add("occupation_primary", i=> ((p, fields) => p.occupationPrimary = fields[i].ToLower()));
            fieldExtractorGenerators.Add("treatment_outcome_first_course", i => ((p, fields) => p.treatmentOutcomeFirstCourse = fields[i].ToLower()));
            fieldExtractorGenerators.Add("days_to_birth", i => ((p, fields) => 
                {
                    try {
                        p.daysToBirth = Convert.ToInt32(fields[i]);
                    } catch (FormatException) {
                        p.daysToBirth = 0;
                    }
                }));

            fieldExtractorGenerators.Add("days_to_death", i => ((p, fields) => 
                {
                    try {
                        p.daysToDeath = Convert.ToInt32(fields[i]);
                    } catch (FormatException) {
                        p.daysToDeath = -1;
                    }
                }));

            string [] lines = File.ReadAllLines(filename);
            var fieldExtractors = new List<Action<Participant, string []>>();
            string [] headers = lines[0].Split('\t');
            int participantIdField = -1;
            int maxFieldNeeded = 0;
            for (int i =0; i < headers.Count(); i++) {
                if (headers[i] == "bcr_patient_uuid") {
                    participantIdField = i;
                    maxFieldNeeded = i;
                } 
                else if (fieldExtractorGenerators.ContainsKey(headers[i])) 
                {
                    fieldExtractors.Add(fieldExtractorGenerators[headers[i]](i));
                    maxFieldNeeded = i;
                }
            }

            if (participantIdField == -1) {
                Console.WriteLine("Couldn't find participant ID field in clinical data file " + filename);
                return;
            }

            //
            // Now process each of the actual lines.
            //
            for(int i = 2; i < lines.Count(); i++) {    // Skip the first three, they're headers
                string [] fields = lines[i].Split('\t');
                if (fields.Count() <= maxFieldNeeded || !participants.ContainsKey(fields[participantIdField].ToLower())) {
                    continue;
                }

                Participant participant = participants[fields[participantIdField].ToLower()];
                participant.metadataPresent = true;
                foreach (var extractor in fieldExtractors)
                {
                    extractor(participant, fields);
                }
            }
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
            var globalScript = new StreamWriter(baseDirectory + @"extractTCGAAdditionalMetadata-all.cmd");
            globalScript.WriteLine(@"del \temp\analysisSummaries.txt");

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


        public class Experiment
        {
            public Participant participant;
            public string disease_abbr;
            public TCGARecord TumorRNAAnalysis;
            public TCGARecord TumorDNAAnalysis;
            public TCGARecord NormalDNAAnalysis;
            public TCGARecord NormalRNAAnalysis;    // This is optional, and filled in if it exists, but we'll generate the experiment even if it doesn't.
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
                        experiment.disease_abbr = RNAtcgaRecord.disease_abbr;

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
                        // Find normal RNA.  Unlike tumor RNA and normal & tumor DNA, if we don't have this we're still willing to use this participant
                        //
                        bestAnalysis = null;
                        bestScore = -1;
                        foreach (var normalSampleEntry in participant.normalSamples)
                        {
                            Sample normalSample = normalSampleEntry.Value;
                            foreach (var normalTCGARecord in normalSample.RNA)
                            {
                                int score = 0;
                                if (RefassemToIndex(normalTCGARecord.refassemShortName) != rnaIndex)
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
                Console.WriteLine(entry.Key + ":" + entry.Value + " (" + readyToGoByDisease[entry.Key] + " ready to go, " + readyToGoExceptForNormalByDisease[entry.Key] + " if you don't count normal)");
            }

            return experiments;
        }

        public static void AddDownloadToScript(TCGARecord tcgaRecord, Dictionary<string, StreamWriter> downloadScripts,  Dictionary<AnalysisType, Pathname> tumorToMachineMapping, Dictionary<string, long> downloadAmounts)
        {
            string machine = tumorToMachineMapping[new AnalysisType(tcgaRecord)];
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

        public static void GenerateDownloadScripts(List<Experiment> experiments,  Dictionary<AnalysisType, Pathname> tumorToMachineMapping)
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


        public static void AddSingleRealignmentToScript(TCGARecord record, StreamWriter script, Dictionary<AnalysisID, TCGARecord> tcgaRecords, Dictionary<AnalysisID, StoredBAM> storedBAMs,
            Dictionary<AnalysisType, Pathname> tumorToMachineMapping, bool local, bool bigMem, bool cluster)
        {
            long lastChunkOfGuid = long.Parse(record.analysis_id.Substring(24), System.Globalization.NumberStyles.HexNumber);
            AnalysisType analysisType = new AnalysisType(record);
            TCGARecord realignSource = tcgaRecords[record.realignedFrom];
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

                if (cluster && (lastChunkOfGuid / 1000 ) % 4 >= 1)  // 50% of the time write to the non-genomics location
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

            UseAnalysis(record.realignedFrom, null, storedBAMs);

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

        public static void GenerateRealignmentScript(List<Experiment> experiments, Dictionary<AnalysisID, TCGARecord> tcgaRecords, Dictionary<AnalysisID, StoredBAM> storedBAMs, Dictionary<AnalysisType, Pathname> tumorToMachineMapping)
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
                        var analysisType = new AnalysisType(analysis);
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

                    Pathname destMachine = tumorToMachineMapping[new AnalysisType(analysis)];
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
                realignClusterScript.WriteLine("job submit /id:%1 /scheduler:gcr");
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

        static public void GenerateExtractionScripts(List<Experiment> experiments)
        {
            var extractRNAInputs = new Dictionary<string, StreamWriter>();  // Maps reference->input file
            var extractDNAInputs = new Dictionary<string, StreamWriter>();  // Maps reference->input file
            StreamWriter extractRNAScript = new StreamWriter(baseDirectory + @"extractRNAScript.cmd");
            StreamWriter extractTumorDNAScript = new StreamWriter(baseDirectory + @"extractDNASamples.cmd");
            var perMachineExtractScripts = new Dictionary<Pathname, StreamWriter>();

            const int nLinesPerPartialExtractionScript = 250000;
            int whichRNAPartialExtractionScript = -1;
            int whichDNAPartialExtractionScript = -1;
            int nLinesInRNAPartialScript = nLinesPerPartialExtractionScript;    // So we'll start a new one right away
            int nLinesInDNAPartialScript = nLinesPerPartialExtractionScript;    // So we'll start a new one right away
            StreamWriter extractRNAPartialScript = null;
            StreamWriter extractDNAPartialScript = null;

            var experimentSets = new Dictionary<DiseaseReferencePair, List<MAFRecord>>();
            foreach (var experiment in experiments)
            {
                var diseaseAndReference = new DiseaseReferencePair(experiment.NormalDNAAnalysis.disease_abbr, experiment.TumorRNAAnalysis.refassemShortName);
                if (!experimentSets.ContainsKey(diseaseAndReference))
                {
                    experimentSets.Add(diseaseAndReference, new List<MAFRecord>());
                }

                string chromPrefix;
                string refassem = experiment.TumorDNAAnalysis.refassemShortName.ToLower();
                if (refassem == "NCBI36_BCCAGSC_variant".ToLower() || refassem == "grch37-lite" || refassem == "grch37" || refassem == "hg19_broad_variant" || refassem == "grch37-lite_wugsc_variant_1" || refassem == "grch37-lite_wugsc_variant_2"
                    || refassem == "grch37-lite-+-hpv_redux-build" || refassem == "HS37D5".ToLower() || refassem == "NCBI36_WUGSC_variant".ToLower())
                {
                    chromPrefix = " ";
                }
                else
                {
                    chromPrefix = " chr";
                }

                if (!extractDNAInputs.ContainsKey(refassem))
                {
                    extractDNAInputs.Add(refassem, new StreamWriter(baseDirectory + @"tumorDNAInput-" + refassem + ".txt"));
                    extractRNAInputs.Add(refassem, new StreamWriter(baseDirectory + @"tumorRNAInput-" + refassem + ".txt"));
                }

                 foreach (MAFRecord mafRecord in experiment.maf)
                {
                    if (experiment.TumorDNAAnalysis.storedBAM != null)
                    {
                        string tumorDNAFileName = diseaseAndReference.disease + @"\" + experiment.TumorDNAAnalysis.analysis_id + "-" +
                            mafRecord.Hugo_symbol + "-chr" + mafRecord.Chrom + "-" + Math.Max(1, mafRecord.Start_position - 200) + "-" +
                           (mafRecord.Start_position + 10);

                        extractDNAInputs[refassem].WriteLine(tumorDNAFileName + "\t" + mafRecord.entire_maf_line);
                        string line = "samtools view " + experiment.TumorDNAAnalysis.storedBAM.bamInfo.FullName + chromPrefix + mafRecord.Chrom.ToUpper() + ":" + Math.Max(1, mafRecord.Start_position - 200) + "-" + (mafRecord.Start_position + 10) + @" > " +
                                                             tumorDNAFileName;
                        extractTumorDNAScript.WriteLine(line);
                        if (nLinesInDNAPartialScript >= nLinesPerPartialExtractionScript)
                        {
                            if (null != extractDNAPartialScript)
                            {
                                extractDNAPartialScript.Close();
                            }
                            whichDNAPartialExtractionScript++;
                            extractDNAPartialScript = new StreamWriter(baseDirectory + @"extractDNAScript-" + whichDNAPartialExtractionScript + ".cmd");
                            nLinesInDNAPartialScript = 0;
                        }
                        extractDNAPartialScript.WriteLine(line);
                        nLinesInDNAPartialScript++;
                    }

                    if (experiment.TumorRNAAnalysis.storedBAM != null)
                    {
                        string tumorRNAFileName = diseaseAndReference.disease + @"\" + experiment.TumorRNAAnalysis.analysis_id + "-" +
                            mafRecord.Hugo_symbol + "-chr" + mafRecord.Chrom + "-" + Math.Max(1, mafRecord.Start_position - 200) + "-" +
                           (mafRecord.Start_position + 10) + "-RNA";

                        extractRNAInputs[refassem].WriteLine(tumorRNAFileName + "\t" + mafRecord.entire_maf_line);
                        string line = "samtools view " + experiment.TumorRNAAnalysis.storedBAM.bamInfo.FullName + chromPrefix + mafRecord.Chrom.ToUpper() + ":" + Math.Max(1, mafRecord.Start_position - 200) + "-" + (mafRecord.Start_position + 10) + @" > " +
                                                             tumorRNAFileName;
                        extractRNAScript.WriteLine(line);
                        if (nLinesInRNAPartialScript >= nLinesPerPartialExtractionScript)
                        {
                            if (null != extractRNAPartialScript)
                            {
                                extractRNAPartialScript.Close();
                            }
                            whichRNAPartialExtractionScript++;
                            extractRNAPartialScript = new StreamWriter(baseDirectory + @"extractRNAScript-" + whichRNAPartialExtractionScript + ".cmd");
                            nLinesInRNAPartialScript = 0;
                        }
                        extractRNAPartialScript.WriteLine(line);
                        nLinesInRNAPartialScript++;
                    }
                }
            }

            extractRNAScript.Close();
            extractTumorDNAScript.Close();
            if (null != extractRNAPartialScript)
            {
                extractRNAPartialScript.Close();
            }

            if (null != extractDNAPartialScript)
            {
                extractDNAPartialScript.Close();
            }

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

        public static void GenerateLAMLMoves(Dictionary<AnalysisID, TCGARecord> tcgaRecords)
        {
            StreamWriter script = new StreamWriter(baseDirectory + @"moveLAML.cmd");

            foreach (var entry in tcgaRecords)
            {
                TCGARecord record = entry.Value;
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

        public static List<AnalysisID> LoadExcludedAnalyses() 
        {
            var excludedAnalyses = new List<AnalysisID>();

            string [] analyses = File.ReadAllLines(@"f:\sequence\reads\tcga\excluded_analyses.txt");

            foreach (var analysis in analyses) {
                excludedAnalyses.Add(analysis.ToLower());
            }

            return excludedAnalyses;
        }

        public static void AddOneTypeOfCountsToMAFs(string[] counts, Dictionary<string, Participant> participants, Dictionary<string, Sample> allSamples, bool isDNA)
        {
            int nDuplicateLines = 0;
            int nNoSample = 0;
            for (int lineNumber = 1; lineNumber <= counts.Count(); lineNumber++ )   // Use 1-based to correspond to how people usually use file lines
            {
                string count = counts[lineNumber-1]; // -1 is because we're 1 based
                string[] fields = count.Split('\t');
                if (fields.Count() != 43)
                {
                    Console.WriteLine("Badly formatted dna count line " + count);
                    continue;
                }

                //
                // Get the analysis ID, which is contained in the filename in field 1.
                //
                string[] filenameFields = fields[1].Split('\\');
                if (filenameFields.Count() != 2 || filenameFields[1].Count() < 36)
                {
                    Console.WriteLine("Can't parse filename field in dna record " + count);
                    continue;
                }
                AnalysisID analysisID = filenameFields[1].Substring(0, 36);
                string gene = fields[2];

                if (!allSamples.ContainsKey(fields[34]) || !allSamples.ContainsKey(fields[35]))
                {
                    //Console.WriteLine("Unable to find tumor and/or normal sample for DNA line " + count);
                    nNoSample++;
                    continue;
                }

                Sample tumorSample = allSamples[fields[34]];
                Sample normalSample = allSamples[fields[35]];

                if (tumorSample.participantID != normalSample.participantID)
                {
                    Console.WriteLine("Samples have different participant IDs " + fields[35] + "(" + tumorSample.participantID + ") and " + fields[36] + "(" + normalSample.participantID);
                    continue;
                }

                Participant participant = participants[tumorSample.participantID];

                //
                // Now find the actual MAFRecord.
                //
                if (!participant.mafsByGene.ContainsKey(gene))
                {
                    Console.WriteLine("Participant " + participant.participantId + " doesn't seem to have any mutations for gene " + gene + " though it has a MAF line");
                    continue;
                }

                MAFRecord mafRecord = null;
                foreach (var mafCandidate in participant.mafsByGene[gene])
                {
                    if (mafCandidate.Start_position == Convert.ToInt32(fields[7]) && mafCandidate.End_position == Convert.ToInt32(fields[8]))
                    {
                        mafRecord = mafCandidate;
                        break;
                    }
                }
                if (null == mafRecord)
                {
                    Console.WriteLine("Couldn't find matching MAF record for count line " + count);
                    continue;
                }

                if (isDNA)
                {
                    if (mafRecord.hasDNACounts)
                    {
                        //
                        // For reasons that I don't entirely understand, some mutations are listed more than once in the MAF files, and so wind up more than once here.
                        //
                        if (mafRecord.nDNAMatchingNormal != Convert.ToInt32(fields[39]) || mafRecord.nDNAMatchingTumor != Convert.ToInt32(fields[40]) ||
                            mafRecord.nDNAMatchingNeither != Convert.ToInt32(fields[41]) || mafRecord.nDNAMatchingBoth != Convert.ToInt32(fields[42]))
                        {
                            Console.WriteLine("Found DNA count record for already filled in MAF with different values (line numbers " + mafRecord.DNALineNumber + " and " + lineNumber + " line: " + count);
                        }
                        nDuplicateLines++;
                        continue;
                    }
                    mafRecord.hasDNACounts = true;
                    mafRecord.DNALineNumber = lineNumber;
                    mafRecord.nDNAMatchingNormal = Convert.ToInt32(fields[39]);
                    mafRecord.nDNAMatchingTumor = Convert.ToInt32(fields[40]);
                    mafRecord.nDNAMatchingNeither = Convert.ToInt32(fields[41]);
                    mafRecord.nDNAMatchingBoth = Convert.ToInt32(fields[42]);
                }
                else
                {
                    if (mafRecord.hasRNACounts)
                    {
                        if (mafRecord.nRNAMatchingNormal != Convert.ToInt32(fields[39]) || mafRecord.nRNAMatchingTumor != Convert.ToInt32(fields[40]) ||
                            mafRecord.nRNAMatchingNeither != Convert.ToInt32(fields[41]) || mafRecord.nRNAMatchingBoth != Convert.ToInt32(fields[42]))
                        {
                            Console.WriteLine("Found RNA count record for already filled in MAF with different values (line numbers " + mafRecord.RNALineNumber + " and " + lineNumber + " line: " + count);
                        }
                        nDuplicateLines++;
                        continue;
                    }
                    mafRecord.hasRNACounts = true;
                    mafRecord.RNALineNumber = lineNumber;
                    mafRecord.nRNAMatchingNormal = Convert.ToInt32(fields[39]);
                    mafRecord.nRNAMatchingTumor = Convert.ToInt32(fields[40]);
                    mafRecord.nRNAMatchingNeither = Convert.ToInt32(fields[41]);
                    mafRecord.nRNAMatchingBoth = Convert.ToInt32(fields[42]);
                }
            }

            if (0 != nDuplicateLines || 0 != nNoSample)
            {
                Console.WriteLine("Found " + nDuplicateLines + " duplicate and  " + nNoSample + " no sample " + (isDNA ? "DNA" : "RNA") + " counts, presumably because of duplicate lines in the MAF files from TCGA.");
            }

        }


        public static void AddCountsToMAFs(Dictionary<string, Participant> participants, Dictionary<string, Sample> allSamples)
        {
            string[] dnaCounts = File.ReadAllLines(@"f:\sequence\Reads\tcga\tumor_dna.txt");    // These are the counts of tumor/normal/both/neither for each MAF line for which we could get data
            AddOneTypeOfCountsToMAFs(dnaCounts, participants, allSamples, true);

            dnaCounts = null;
            string[] rnaCounts = File.ReadAllLines(@"f:\sequence\Reads\tcga\tumor_rna.txt");
            AddOneTypeOfCountsToMAFs(rnaCounts, participants, allSamples, false);
        }

        public static void BuildMAFsByGene(Dictionary<string, Participant> participants)
        {
            foreach (var entry in participants)
            {
                var participant = entry.Value;
                if (participant.mafs.Count() == 0)
                {
                    continue;
                }

                foreach (MAFRecord mafRecord in participant.mafs[0])   // We only ever use mafs[0]; other ones are just different mutation calls for the same participant, which we ignore
                {
                    if (!participant.mafsByGene.ContainsKey(mafRecord.Hugo_symbol))
                    {
                        participant.mafsByGene.Add(mafRecord.Hugo_symbol, new List<MAFRecord>());
                    }

                    participant.mafsByGene[mafRecord.Hugo_symbol].Add(mafRecord);
                }

            }
        }

        public static void AddAnalysisToMakeIsoformScript(TCGARecord analysis, Dictionary<Pathname, StreamWriter> scripts, string disease_abbr)
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

        public static void GenerateMakeIsoformsScripts(List<Experiment> experiments)
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

        static void DumpExperiments(List<Experiment> experiments) {
            StreamWriter outputFile = new StreamWriter(baseDirectory + @"experiments.txt");

            outputFile.WriteLine("disease_abbr\treference\tparticipantID\tTumorRNAAnalysis\tTumorDNAAnalysis\tNormalDNAAnalysis\tNormalRNAAnalysis\ttumorRNAPathname\ttumorDNAPathname\t" +
                "normalDNAPathname\tNormalRNAPathname\tVCFPathname\tgender\tdaysToBirth\tdaysToDeath\tOrigTumorDNAAliquotID\tTumorAllcountFile\tNormalAllcountFile");

            foreach (Experiment experiment in experiments) {
                outputFile.Write(experiment.disease_abbr + "\t" + experiment.TumorRNAAnalysis.refassemShortName + "\t" + experiment.participant.participantId + "\t" + experiment.TumorRNAAnalysis.analysis_id + "\t" + experiment.TumorDNAAnalysis.analysis_id + "\t" + experiment.NormalDNAAnalysis.analysis_id + "\t");

                if (null != experiment.NormalRNAAnalysis)
                {
                    outputFile.Write(experiment.NormalRNAAnalysis.analysis_id + "\t");
                }
                else
                {
                    outputFile.Write("\t");
                }


                if (experiment.TumorRNAAnalysis.storedBAM != null)
                {
                    outputFile.Write(experiment.TumorRNAAnalysis.storedBAM.bamInfo.FullName + "\t");
                }
                else
                {
                    outputFile.Write("\t");
                }

                if (experiment.TumorDNAAnalysis.storedBAM != null)
                {
                    outputFile.Write(experiment.TumorDNAAnalysis.storedBAM.bamInfo.FullName + "\t");
                }
                else
                {
                    outputFile.Write("\t");
                }

                if (experiment.NormalDNAAnalysis.storedBAM != null)
                {
                    outputFile.Write(experiment.NormalDNAAnalysis.storedBAM.bamInfo.FullName + "\t");
                }
                else
                {
                    outputFile.Write("\t");
                }

                if (experiment.NormalRNAAnalysis != null && experiment.NormalRNAAnalysis.storedBAM != null)
                {
                    outputFile.Write(experiment.NormalRNAAnalysis.storedBAM.bamInfo.FullName + "\t");
                }
                else
                {
                    outputFile.Write("\t");
                }

                if (experiment.NormalDNAAnalysis.storedBAM != null && experiment.NormalDNAAnalysis.storedBAM.vcfInfo != null)
                {
                    outputFile.Write(experiment.NormalDNAAnalysis.storedBAM.vcfInfo.FullName + "\t");
                }
                else
                {
                    outputFile.Write("\t");
                }

                if (experiment.participant == null || !experiment.participant.metadataPresent)
                {
                    outputFile.Write("\t\t\t");
                } 
                else
                {
                    outputFile.Write(experiment.participant.gender + "\t");
                    if (experiment.participant.daysToBirth != 0)
                    {
                        outputFile.Write("" + experiment.participant.daysToBirth + "\t");
                    }
                    else
                    {
                        outputFile.Write("\t");
                    }

                    if (experiment.participant.daysToDeath >= 0)
                    {
                        outputFile.Write("" + experiment.participant.daysToDeath + "\t");
                    }
                    else
                    {
                        outputFile.Write("\t");
                    }
                }

                if (experiment.TumorDNAAnalysis.realignSource != null)
                {
                    outputFile.Write(experiment.TumorDNAAnalysis.realignSource.aliquot_id + "\t");
                }
                else
                {
                    outputFile.Write(experiment.TumorDNAAnalysis.aliquot_id + "\t");
                }

                if (experiment.TumorRNAAnalysis.storedBAM != null && experiment.TumorRNAAnalysis.storedBAM.allCountInfo != null)
                {
                    outputFile.Write(experiment.TumorRNAAnalysis.storedBAM.allCountInfo.FullName + "\t");
                }
                else
                {
                    outputFile.Write("\t");
                }

                if (experiment.NormalRNAAnalysis != null && experiment.NormalRNAAnalysis.storedBAM != null && experiment.NormalRNAAnalysis.storedBAM.allCountInfo != null)
                {
                    outputFile.WriteLine(experiment.NormalRNAAnalysis.storedBAM.allCountInfo.FullName + "\t");
                }
                else
                {
                    outputFile.WriteLine("\t");
                }

            }

            outputFile.Close();
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

        public static void GenerateMutationDistributionByGene(List<Experiment> experiments)
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
                foreach (MAFRecord mafRecord in experiment.maf)
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

        public static void GenerateIsoformFileListsByDiseaseAndMutation(List<Experiment> experiments)
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

        public static string WindowsToLinuxPathname(string windowsPathname)
        {
            string[] components = windowsPathname.Split('\\');


            if (components.Count() < 5)
            {
                Console.WriteLine("WindowsToLinuxPathname: unable to parse " + windowsPathname);
                return @"/dev/null";
            }

            string linuxPathname = @"/mnt/" + components[2] + @"/" + components[3][0];  // components[3][0] is the drive letter, and strips off the $

            for (int i = 4; i < components.Count(); i++)
            {
                linuxPathname += @"/" + components[i];
            }


            return linuxPathname;
        }

        public static void GenerateVariantCallingScript(List<Experiment> experiments)
        {
            int nSamplesToVariantCall = 0;
            int nSamplesAlreadyCalled = 0;
            int nSamplesNotYetReadyToCall = 0;
            int nSkipped = 0;

            StreamWriter variantCallScript = null;

            foreach (var experiment in experiments)
            {
                if (experiment.NormalDNAAnalysis.storedBAM == null)
                {
                    nSamplesNotYetReadyToCall++;
                }
                else if (experiment.NormalDNAAnalysis.storedBAM.vcfInfo != null) 
                {
                    nSamplesAlreadyCalled++;
                } 
                else if (!experiment.NormalDNAAnalysis.localRealign)
                {
                    nSkipped++; // First pass, just do the ones we locally realigned, because there's some confusion matching the chromosome names with the references for the ones in TCGA
                }
                else
                {
                    nSamplesToVariantCall++;

                    if (variantCallScript == null)
                    {
                        variantCallScript = new StreamWriter(baseDirectory + @"callVariants");
                    }

                    variantCallScript.Write("date\n");
                    variantCallScript.Write(@"cat ~/genomes/" + experiment.NormalDNAAnalysis.refassemShortName + "-100k-regions | parallel -k -j `cat ~/ncores` \" freebayes --region {} --fasta-reference ~/genomes/" + experiment.NormalDNAAnalysis.refassemShortName + ".fa " +
                        WindowsToLinuxPathname(experiment.NormalDNAAnalysis.storedBAM.bamInfo.FullName) + " \" | ~/freebayes/vcflib/bin/vcffirstheader | ~/freebayes/vcflib/bin/vcfstreamsort -w 1000 | ~/freebayes/vcflib/bin/vcfuniq > " +
                        experiment.NormalDNAAnalysis.analysis_id + ".vcf\n");
                    variantCallScript.Write("if [ $? = 0 ]; then\n");
                    variantCallScript.Write(@"    cp " + experiment.NormalDNAAnalysis.analysis_id + ".vcf " + WindowsToLinuxPathname(GetDirectoryPathFromFullyQualifiedFilename(experiment.NormalDNAAnalysis.storedBAM.bamInfo.FullName)) + "\n");
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

            Console.WriteLine("" + nSamplesAlreadyCalled + " have been variant called, " + nSamplesNotYetReadyToCall + " are waiting for other things, " + nSkipped + " were skipped, and " + nSamplesToVariantCall + " are ready to call.");
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

        static public void GenerateAllcountScript( Dictionary<AnalysisID, TCGARecord> tcgaRecords)
        {
            int nWithAllcount = 0;
            int nAllcountsToGenerate = 0;
            int nAwaitingDownload = 0;
            int indexMachine = 0;
            StreamWriter allcountScript = null;
            foreach (var entry in tcgaRecords)
            {
                TCGARecord tcgaRecord = entry.Value;
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

        static void GenerateListOfAllcountFiles(Dictionary<AnalysisID, StoredBAM> storedBAMs)
        {
            StreamWriter output = new StreamWriter(baseDirectory + "allcount_files.txt");

            foreach (var entry in storedBAMs)
            {
                StoredBAM storedBAM = entry.Value;

                if (storedBAM.allCountInfo != null)
                {
                    output.WriteLine(storedBAM.allCountInfo.FullName);
                }
            }

            output.Close();
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

            TumorAndNormal.Add("tumor");
            TumorAndNormal.Add("normal");

            InitializeMachines();

            //var geneToSize = GenerateSizeOfProteinCodingRegionMap();

            AllIsoformsFile = new StreamWriter(baseDirectory + @"isoform-files\all_isoform_files.txt");

            List<AnalysisID> excludedAnalyses = LoadExcludedAnalyses();

            var tumorToMachineMapping = GenerateTumorToMachineMapping();
            Dictionary<AnalysisID, StoredBAM> storedBAMs = LoadStoredBAMs(tumorToMachineMapping);
            AllIsoformsFile.Close();
            GenerateListOfAllcountFiles(storedBAMs);
            var tcgaRecords = LoadTCGARecords(storedBAMs, excludedAnalyses);
            LoadTCGARecordsForLocalRealigns(tcgaRecords, storedBAMs);
//RebuildRealignmentAnalysesFromExperimentsFile(tcgaRecords); return;
            LoadTCGAAdditionalMetadata(tcgaRecords);
            GenereateAdditionalTCGAMetadataGeneratingScript(tcgaRecords, storedBAMs, tumorToMachineMapping);
            VerifyStoredBAMPaths(storedBAMs, tcgaRecords, tumorToMachineMapping);
            GenerateAllcountScript(tcgaRecords);

            var sampleToParticipantIDMap = CreateSampleToParticipantIDMap(tcgaRecords);
            DumpSampleToParticipantIDMap(sampleToParticipantIDMap);
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

            //BuildMAFsByGene(participants);
            //AddCountsToMAFs(participants, allSamples);

            var experiments = BuildExperiments(participants);
            DumpExperiments(experiments);
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
        }
    }
}
