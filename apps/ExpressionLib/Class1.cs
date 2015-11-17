using System;
using System.Collections;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using System.IO;
using System.Xml.Linq;
using System.IO.Compression;


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

namespace ExpressionLib
{
    public class ExpressionTools
    {
        public static string ShareFromPathname(string pathname)
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

        public class Variant : IComparable<Variant>
        {
            public static Variant GenerateFromVCF(string line)
            {
                var variant = new Variant();

                var fields = line.Split('\t');

                try
                {
                    variant.chromosome = fields[0].ToLower();
                    variant.pos = Convert.ToInt32(fields[1]);
                    variant.reference = fields[3];
                    variant.alt = fields[4];
                    variant.qual = Convert.ToDouble(fields[5]);
                    variant.line = line;
                }
                catch (FormatException)
                {
                    return null;
                }

                return variant;
            }

            public int CompareTo(Variant other)
            {
                if (other == null) return 1;
                return pos.CompareTo(other.pos);
            }

            public string chromosome;
            public int pos;
            public string reference;    // it's "ref" in the format spec, but that's a C# keyword
            public string alt;
            public double qual;
            public string line;
        }


        public class VCF
        {
            public Dictionary<string, List<Variant>> variants = new Dictionary<string, List<Variant>>();

            public VCF(string filename)
            {
                var reader = new StreamReader(filename);

                string line;

                while (null != (line = reader.ReadLine()))
                {
                    if (line.Count() == 0 || line[0] == '#')
                    {
                        continue;   // Not a variant line
                    }

                    Variant variant = Variant.GenerateFromVCF(line);
                    if (null == variant)
                    {
                        Console.WriteLine("VCF file " + filename + " contains an unparsable line " + line);
                        break;
                    }
                    if (!variants.ContainsKey(variant.chromosome))
                    {
                        variants.Add(variant.chromosome, new List<Variant>());
                    }

                    List<Variant> lov = variants[variant.chromosome];
                    int lovSize = lov.Count();

                    if (lovSize != 0 && lov[lovSize - 1].pos > variant.pos)
                    {
                        Console.WriteLine("vcf " + filename + " has out-of-order line " + line);
                        return;
                    }

                    lov.Add(variant);
                }

                reader.Close();
            }
        }

        public class VCFEnumerator : IEnumerator<Variant>, IEnumerable<Variant>
        {
            public VCFEnumerator(VCF vcf, string chromosome, int startPos, int _endPos)
            {
                endPos = _endPos;
                current = null;

                chromosome = chromosome.ToLower();

                if (vcf == null)
                {
                    variants = null;
                    return;
                }

                if (!vcf.variants.ContainsKey(chromosome))
                {
                    chromosome = "chr" + chromosome;  // Trty prepending "chr". The lack of standards is so annoying!

                    if (!vcf.variants.ContainsKey(chromosome))
                    {
                        variants = null;
                        return;
                    }

                }

                variants = vcf.variants[chromosome];
                numVariants = variants.Count();

                var dummyVariant = new Variant();
                dummyVariant.pos = startPos;
                currentVapIndex = variants.BinarySearch(dummyVariant);

                if (currentVapIndex < 0)
                {
                    currentVapIndex = ~currentVapIndex; // Negative means it wasn't found exactly, and returns the index of the next greater element, or of the size of the list if it's bigger than everything
                }
                else
                {
                    //
                    // Since we may have equal valued elements in variants, and BinarySearch just returns any one of them, walk backward until we hit something different or run
                    // off the beginning of the list.
                    //
                    while (currentVapIndex > 0 && variants[currentVapIndex - 1].pos == startPos)
                    {
                        currentVapIndex--;
                    }
                }

                firstVapIndex = currentVapIndex;
            }

            public bool MoveNext()
            {
                if (null == variants)
                {
                    return false;
                }

                if (null != current)
                {
                    currentVapIndex++;
                }

                if (currentVapIndex >= numVariants || variants[currentVapIndex].pos > endPos)
                {
                    return false;
                }

                current = variants[currentVapIndex];

                return true;
            }

            public void Reset()
            {
                current = null;
                currentVapIndex = firstVapIndex;
            }

            public IEnumerator<Variant> GetEnumerator()
            {
                return this;    // Yick.
            }

            IEnumerator IEnumerable.GetEnumerator()
            {
                return this;
            }

            public Variant current;
            public Variant Current { get { return current; } }
            object System.Collections.IEnumerator.Current { get { return Current; } }


            public void Dispose() { }

            List<Variant> variants;
            int numVariants;
            int firstVapIndex;
            int currentVapIndex;
            int endPos;

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

            public TumorType tumorType;
            public LibraryStrategy libraryStrategy;
            public bool isNormal;
        }
        public class StoredBAM
        {
            public AnalysisID analysisID = null;
            public Pathname directoryName = null;
            public string machineName = null;
            public char driveLetter = ' ';
            public FileInfo bamInfo = null;
            public FileInfo baiInfo = null;
            public FileInfo vcfInfo = null;
            public FileInfo allCountInfo = null;
            public long totalSize = 0;
            public bool needed = false; // Do we need this as input for some experiment?
            public bool isRealingned = false;   // Was this realined, or is it straight from TCGA
            public bool bamHashBeingComputed = false;   // Did we get an IO error reading the bam hash, meaning it's being computed now?
            public string bamHash = null;  // MD5 of the BAM file
            public bool baiHashBeingComputed = false;   // Did we get an IO error reading the bai hash, meaning it's being computed now?
            public string baiHash = null;  // MD5 of the BAI file
            public string isoformCount = null;  // File with the count of isoform references
            public bool chrStateKnown = false;  // Do we know whether this analysis uses "chr" before chromosome names (i.e, is usesChr set)
            public bool usesChr; // Does this analysis use "chr" before chromosome names (you can only check this is chrStateKnown is set)
        }

        static public string GetDirectoryPathFromFullyQualifiedFilename(string filename)
        {
            if (filename.Count() < 1 || filename[0] != '\\')
            {
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
                                    storedBAM.bamHashBeingComputed = true;
                                    Console.WriteLine("IO Exception reading m5 file " + md5Pathname + ".  If you're not computing that hash now, maybe something's wrong; error: " + e.Message);
                                }
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
                                try
                                {
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
                                    storedBAM.baiHashBeingComputed = true;
                                    Console.WriteLine("IO Exception reading m5 file " + md5Pathname + ".  If you're not computing that hash now, maybe something's wrong; error: " + e.Message);
                                }
                            }
                        }
                        else if (file.Count() > 4 && file.Substring(file.Count() - 4) == ".tar")
                        {
                            hasTarFile = true;
                        }
                        else if (file.Count() > 3 && file.Substring(file.Count() - 4) == ".md5")
                        {
                            // Do nothing; checksum files are expected.
                        }
                        else if (file == subdir + @"\" + analysisID + @".vcf")
                        {
                            storedBAM.vcfInfo = new FileInfo(file);
                        }
                        else if (file.Count() > 12 && file.Substring(file.Count() - 13).ToLower() == "-isoforms.txt")
                        {
                            string expectedFilename = analysisID + "-isoforms.txt";
                            if (file.Count() < expectedFilename.Count() || file.Substring(file.Count() - expectedFilename.Count()) != expectedFilename)
                            {
                                Console.WriteLine("Incorrect isoform file " + file);
                            }
                            else
                            {
                                storedBAM.isoformCount = file;
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
                        if (!hasTarFile)
                        {
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
                }
                else if (analysisID.Contains(".partial"))
                {
                    Console.WriteLine("Partial download at " + subdir);
                }
            }
        }

        public static string [] TumorAndNormal = {"tumor", "normal"};
        public static string[] LibraryStrategies = { "wgs", "wxs", "rna" };


        public static void LoadStoredBamsForMachine(Pathname directory, Dictionary<AnalysisID, StoredBAM> storedBAMs)
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
                            LoadStoredBAMsForDirectory(sampleDir, storedBAMs);
                        }
                    }
                }
            }
        }

        public static void LoadChrStateFile(Pathname filename, Dictionary<AnalysisID, StoredBAM> storedBAMs)
        {
            if (!File.Exists(filename))
            {
                Console.WriteLine("Couldn't find chr state file " + filename + ", assuming you haven't computed any chr states yet");
                return;
            }
            var reader = new StreamReader(filename);

            string line;
            while (null != (line = reader.ReadLine()))
            {
                var fields = line.Split('\t');
                if (fields.Count() != 2)
                {
                    Console.WriteLine("Chr state file " + filename + " has bad line " + line);
                    return;
                }

                string analysisID = fields[0];
                bool usesChr;
                if (fields[1] == "chr")
                {
                    usesChr = true;
                }
                else if (fields[1] == "noChr")
                {
                    usesChr = false;
                }
                else
                {
                    Console.WriteLine("Couldn't parse chr state for line " + line + " from file " + filename);
                    return;
                }

                if (storedBAMs.ContainsKey(analysisID))
                {
                    storedBAMs[analysisID].chrStateKnown = true;
                    storedBAMs[analysisID].usesChr = usesChr;
                }
            }

            reader.Close();
        }


        public class TCGARecord
        {
            public AnalysisID analysis_id = null;
            public ParticipantID participant_id = null;
            public Pathname bamFileName = null;
            public string bamHash = null;
            public Pathname baiFileName = null;
            public string baiHash = null;
            public DateTime lastModified;
            public DateTime uploadDate;
            public SampleID aliquot_id = null;
            public string refassemShortName = null;
            public TumorType disease_abbr = null;
            public long totalFileSize = 0;
            public SampleType sampleType;
            public LibraryStrategy library_strategy = null;
            public string center_name = null;
            public bool tumorSample = false;
            public bool localRealign = false;
            public AnalysisID realignedFrom = null;
            public TCGARecord realignSource = null;
            public bool isWGS = false;

            public bool hasAdditionalData = false;  // Do we have the extra metadata that's shown below?
            public int readSampleSize = 0;
            public int minReadLength = 0;
            public int maxReadLength = 0;
            public float meanReadLength = 0;
            public int minGoodBases = 0;
            public int maxGoodBases = 0;
            public float meanGoodBases = 0;
            public bool anyPaired = false;
            public bool allPaired = false;


            public string chromPrefix = null;
            public StoredBAM storedBAM = null;

            public string GetContainingDirectory()
            {
                return @"d:\tcga\" + library_strategy + (tumorSample ? @"\tumor\" : @"\normal\") + disease_abbr;
            }
        } // TCGARecord

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

        public static void DumpExperimentsToFile(List<Experiment> experiments, string filename)
        {
            StreamWriter outputFile = new StreamWriter(filename);

            outputFile.WriteLine("disease_abbr\treference\tparticipantID\tTumorRNAAnalysis\tTumorDNAAnalysis\tNormalDNAAnalysis\tNormalRNAAnalysis\ttumorRNAPathname\ttumorDNAPathname\t" +
                "normalDNAPathname\tNormalRNAPathname\tVCFPathname\tgender\tdaysToBirth\tdaysToDeath\tOrigTumorDNAAliquotID\tTumorAllcountFile\tNormalAllcountFile\tmafFile");

            foreach (Experiment experiment in experiments)
            {
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
                    outputFile.Write(experiment.NormalRNAAnalysis.storedBAM.allCountInfo.FullName + "\t");
                }
                else
                {
                    outputFile.Write("\t");
                }

                outputFile.Write(experiment.participant.mafFile + "\t");

                outputFile.WriteLine();
            }

            outputFile.Close();
        } // DumpExperimentsToFile

        public static List<Experiment> LoadExperimentsFromFile(string filename, Dictionary<ParticipantID, Participant> participants, Dictionary<AnalysisID, TCGARecord> tcgaRecords)
        {
            var experiments = new List<Experiment>();

            var reader = new StreamReader(filename);

            string line;
            reader.ReadLine();  // Skip the header line
            while (null != (line = reader.ReadLine())) {
                var fields = line.Split('\t');

                var experiment = new Experiment();
                experiment.disease_abbr = fields[0];
                experiment.participant = participants[fields[2]];
                experiment.TumorRNAAnalysis = tcgaRecords[fields[3]];
                experiment.TumorDNAAnalysis = tcgaRecords[fields[4]];
                experiment.NormalDNAAnalysis = tcgaRecords[fields[5]];
                if (fields[6] == "") {
                    experiment.NormalRNAAnalysis = null;
                } else {
                    experiment.NormalRNAAnalysis = tcgaRecords[fields[6]];
                }
                experiment.maf = experiment.participant.mafs[0];

                experiments.Add(experiment);
            }

            return experiments;
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
        } // MAFRecord

        public class Sample
        {
            public SampleID sampleId;
            public bool isTumor;
            public ParticipantID participantID;
            public List<TCGARecord> RNA = new List<TCGARecord>();
            public List<TCGARecord> DNA = new List<TCGARecord>();
        } // Sample


        public class Participant
        {
            public string participantId;
            public List<List<MAFRecord>> mafs = new List<List<MAFRecord>>();
            public string mafFile = "";
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
        } // Participant

        public static void AddMAFFileToParticipants(Pathname inputFilename, Dictionary<ParticipantID, Participant> participants, Dictionary<SampleID, ParticipantID> sampleToParticipantIDMap)
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

                if (fields[4].ToLower() == "mt" || fields[4].ToLower() == "chrm" || fields[4].ToLower() == "m")
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
                    mafs.Add(participantID, new Dictionary<string, List<MAFRecord>>());
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
                    if (participant.mafs.Count() == 0)
                    {
                        participant.mafFile = inputFilename;    // Just keep the first one, it's all we use anyway.
                    }
                    participant.mafs.Add(mafEntry.Value);
                }
            }

        } // AddMAFFileToParticipants

        public static void AddAllMAFFilesToParticipants(Dictionary<ParticipantID, Participant> participants, Dictionary<SampleID, ParticipantID> sampleToParticipantIDMap)
        {
            foreach (var file in Directory.GetFiles(@"f:\sequence\Reads\tcga\mafs\"))
            {
                AddMAFFileToParticipants(file, participants, sampleToParticipantIDMap);
            }
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
                    if (sample.isTumor != tcgaEntry.Value.tumorSample)
                    {
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

            return participants;
        } // BuildParticipantData

        public static void ProcessClinicalFile(Dictionary<ParticipantID, Participant> participants, string filename)
        {

            //
            // These files don't all have the same headers.  Scan through the first line to find the column numbers we care about.
            //
            var fieldExtractorGenerators = new Dictionary<string, Func<int, Action<Participant, string[]>>>();
            fieldExtractorGenerators.Add("gender", i => ((p, fields) => p.gender = fields[i].ToLower()));
            fieldExtractorGenerators.Add("race", i => ((p, fields) => p.race = fields[i].ToLower()));
            fieldExtractorGenerators.Add("ethnicity", i => ((p, fields) => p.ethnicity = fields[i].ToLower()));
            fieldExtractorGenerators.Add("history_other_malignancy", i => ((p, fields) => p.historyOfOtherMalignancy = fields[i].ToLower()));
            fieldExtractorGenerators.Add("tumor_status", i => ((p, fields) => p.tumorStatus = fields[i].ToLower()));
            fieldExtractorGenerators.Add("vital_status", i => ((p, fields) => p.vitalStatus = fields[i].ToLower()));
            fieldExtractorGenerators.Add("occupation_primary", i => ((p, fields) => p.occupationPrimary = fields[i].ToLower()));
            fieldExtractorGenerators.Add("treatment_outcome_first_course", i => ((p, fields) => p.treatmentOutcomeFirstCourse = fields[i].ToLower()));
            fieldExtractorGenerators.Add("days_to_birth", i => ((p, fields) =>
            {
                try
                {
                    p.daysToBirth = Convert.ToInt32(fields[i]);
                }
                catch (FormatException)
                {
                    p.daysToBirth = 0;
                }
            }));

            fieldExtractorGenerators.Add("days_to_death", i => ((p, fields) =>
            {
                try
                {
                    p.daysToDeath = Convert.ToInt32(fields[i]);
                }
                catch (FormatException)
                {
                    p.daysToDeath = -1;
                }
            }));

            string[] lines = File.ReadAllLines(filename);
            var fieldExtractors = new List<Action<Participant, string[]>>();
            string[] headers = lines[0].Split('\t');
            int participantIdField = -1;
            int maxFieldNeeded = 0;
            for (int i = 0; i < headers.Count(); i++)
            {
                if (headers[i] == "bcr_patient_uuid")
                {
                    participantIdField = i;
                    maxFieldNeeded = i;
                }
                else if (fieldExtractorGenerators.ContainsKey(headers[i]))
                {
                    fieldExtractors.Add(fieldExtractorGenerators[headers[i]](i));
                    maxFieldNeeded = i;
                }
            }

            if (participantIdField == -1)
            {
                Console.WriteLine("Couldn't find participant ID field in clinical data file " + filename);
                return;
            }

            //
            // Now process each of the actual lines.
            //
            for (int i = 2; i < lines.Count(); i++)
            {    // Skip the first three, they're headers
                string[] fields = lines[i].Split('\t');
                if (fields.Count() <= maxFieldNeeded || !participants.ContainsKey(fields[participantIdField].ToLower()))
                {
                    continue;
                }

                Participant participant = participants[fields[participantIdField].ToLower()];
                participant.metadataPresent = true;
                foreach (var extractor in fieldExtractors)
                {
                    extractor(participant, fields);
                }
            }
        } // ProcessClinicalFile

        public static string RefassemToIndex(string refassem)
        {
            string refassem_l = refassem.ToLower();
            if (refassem_l == "hg19")
            {
                return @"hg19-20";
            }

            if (refassem_l == "grch37-lite" || refassem_l == "GRCh37-lite_WUGSC_variant_1".ToLower() || refassem_l == "GRCh37-lite_WUGSC_variant_2".ToLower()
            || refassem_l == "GRCh37-lite-+-HPV_Redux-build".ToLower())
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
        static void AddOneTypeOfCountsToMAFs(string[] counts, Dictionary<string, Participant> participants, Dictionary<string, Sample> allSamples, bool isDNA)
        {
            int nDuplicateLines = 0;
            int nNoSample = 0;
            for (int lineNumber = 1; lineNumber <= counts.Count(); lineNumber++)   // Use 1-based to correspond to how people usually use file lines
            {
                string count = counts[lineNumber - 1]; // -1 is because we're 1 based
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

        } // AddOneTypeOfCountsToMAFs


        public static void AddCountsToMAFs(Dictionary<string, Participant> participants, Dictionary<string, Sample> allSamples)
        {
            string[] dnaCounts = File.ReadAllLines(@"f:\sequence\Reads\tcga\tumor_dna.txt");    // These are the counts of tumor/normal/both/neither for each MAF line for which we could get data
            AddOneTypeOfCountsToMAFs(dnaCounts, participants, allSamples, true);

            dnaCounts = null;
            string[] rnaCounts = File.ReadAllLines(@"f:\sequence\Reads\tcga\tumor_rna.txt");
            AddOneTypeOfCountsToMAFs(rnaCounts, participants, allSamples, false);
        } // AddCountsToMAFs

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
        } // BuildMAFsByGene
        public static Dictionary<AnalysisID, ExpressionTools.TCGARecord> LoadTCGARecords(Dictionary<AnalysisID, ExpressionTools.StoredBAM> storedBAMs, List<AnalysisID> excludedAnalyses, string filename)
        {
            var tcgaRecords = new Dictionary<AnalysisID, ExpressionTools.TCGARecord>();

            XDocument tcga_all = XDocument.Load(filename);

            var samplesQuery = from Result in tcga_all.Element("ResultSet").Elements("Result")
                               select Result;

            foreach (var result in samplesQuery)
            {
                var record = new ExpressionTools.TCGARecord();

                SampleType sample_type = result.Element("sample_type").Value;
                if (sample_type.Count() < 2)
                {
                    continue;
                }

                if (sample_type.Count() > 2)
                {
                    if (sample_type.Count() != 9 || sample_type.Substring(0, 7) != "TARGET_")
                    {
                        Console.WriteLine("Unknown sample type " + sample_type);
                        continue;
                    }

                    sample_type = sample_type.Substring(7, 2);
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

                if (excludedAnalyses != null && excludedAnalyses.Contains(record.analysis_id))
                {
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

                if (storedBAMs != null && storedBAMs.ContainsKey(record.analysis_id))
                {
                    record.storedBAM = storedBAMs[record.analysis_id];
                    if (record.storedBAM.bamInfo.Length + record.storedBAM.baiInfo.Length != record.totalFileSize)
                    {
                        Console.WriteLine("File size mismatch for " + record.storedBAM.bamInfo.FullName + " should have total size " + record.totalFileSize + " but has " + record.storedBAM.bamInfo.Length + record.storedBAM.baiInfo.Length);
                    }
                }
                else
                {
                    record.storedBAM = null;
                }

                tcgaRecords.Add(record.analysis_id, record);
            }

            return tcgaRecords;

        }

        public static void LoadTCGAAdditionalMetadata(Dictionary<AnalysisID, ExpressionTools.TCGARecord> tcgaRecords)
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

        } // LoadTCGAAdditionalMetadata


        public static void LoadTCGARecordsForLocalRealigns(Dictionary<AnalysisID, ExpressionTools.TCGARecord> tcgaRecords, Dictionary<AnalysisID, ExpressionTools.StoredBAM> storedBAMs, string filename)
        {
            if (!File.Exists(filename))
            {
                return;
            }
            var realigns = File.ReadAllLines(filename);

            //
            // Format is the list of fields from a TCGARecord separated by *s.
            //
            foreach (var realign in realigns)
            {
                var fields = realign.Split('*');
                var record = new ExpressionTools.TCGARecord();
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

                if (storedBAMs != null && storedBAMs.ContainsKey(record.analysis_id))
                {
                    record.storedBAM = storedBAMs[record.analysis_id];
                    if (record.realignSource.storedBAM != null && record.realignSource.storedBAM.bamInfo.Length > record.storedBAM.bamInfo.Length * 2)
                    {
                        Console.WriteLine("Suspiciously small realigned file " + record.storedBAM.bamInfo.FullName + " size " + record.storedBAM.bamInfo.Length +
                            " while source " + record.realignSource.storedBAM.bamInfo.FullName + " is " + record.realignSource.storedBAM.bamInfo.Length);
                    }
                }

                tcgaRecords.Add(record.analysis_id, record);
            }
        } // LoadTCGARecordsForLocalRealigns

        public static void UseAnalysis(AnalysisID analysisID, List<AnalysisID> neededFiles, Dictionary<AnalysisID, ExpressionTools.StoredBAM> storedBAMs)
        {
            if (storedBAMs != null && storedBAMs.ContainsKey(analysisID))
            {
                storedBAMs[analysisID].needed = true;
            }
            else if (neededFiles != null && !neededFiles.Contains(analysisID))
            {
                neededFiles.Add(analysisID);
            }
        } // UseAnalysis

        public static Dictionary<SampleID, ParticipantID> CreateSampleToParticipantIDMap(Dictionary<AnalysisID, ExpressionTools.TCGARecord> tcgaEntries)
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
        } // CreateSampleToParticipantIDMap

        public class GeneMap
        {
            public class Region : IComparable<Region>
            {
                public Region(int start_, int end_, string gene)
                {
                    start = start_;
                    end = end_;
                    genes = new List<string>();
                    genes.Add(gene);
                }

                public Region(int start_)
                {
                    start = start_; // Version for use in comparison
                }

                public Region(Region oldRegion)
                {
                    start = oldRegion.start;
                    end = oldRegion.end;
                    if (oldRegion.genes != null)
                    {
                        genes = new List<string>();

                        foreach (var gene in oldRegion.genes)
                        {
                            genes.Add(gene);
                        }
                    }
                }

                public int CompareTo(Region other)
                {
                    return start.CompareTo(other.start);
                }

                public int start;
                public int end;
                public List<string> genes = null;
            }

            public class Chromosome
            {
                public Chromosome(string name_)
                {
                    name = name_;
                }

                public void AddMapping(string gene, int start, int end) // NB: This code is pretty much the same as the C++ code in CountReadsCovering.cpp.
                {
                    var comparisonRegion = new Region(start);
                    int coveredUpTo = start;
                    while (coveredUpTo < end)
                    {
                        comparisonRegion.start = coveredUpTo;

                        int index = regions.BinarySearch(comparisonRegion);
                        if (index < 0)
                        {
                            //
                            // The List<T>.BinarySearch routine returns the bitwise complement of the index of the next greater element
                            // if there is no exact match.  Fix that.
                            //
                            index = ~index;

                            //
                            // Furthermore, we want the first that's less than or equal to our index, so back up one.
                            //
                            index--;     // NB: This may now be -1 if we were at the beginning of the list
                        }

                        Region region;
                        if (index < 0)
                        {
                            region = null;
                        }
                        else
                        {
                            region = regions[index];
                        }

                        if (region == null || region.end <= coveredUpTo)
                        {
                            //
                            // There's no region covering the first part of what we need.  Create a new one.
                            //
                            Region nextRegion;
                            if (index < regionsSize - 1)
                            {
                                nextRegion = regions[index + 1];
                            }
                            else
                            {
                                nextRegion = null;
                            }

                            int thisEnd;
                            if (null != nextRegion && nextRegion.start < end)
                            {
                                thisEnd = nextRegion.start;
                            }
                            else
                            {
                                thisEnd = end;
                            }

                            regions.Insert(index+1, new Region(coveredUpTo, thisEnd, gene));
                            regionsSize++;
                            coveredUpTo = thisEnd;
                        }
                        else
                        {
                            //
                            // The region before us covers up our start.
                            //
                            if (region.start < coveredUpTo)
                            {
                                //
                                // It has part hanging off of the beginning.  This should only happen if we're the first part of an exon.  Chop it in two, but don't
                                // add us to the second half; instead that'll happen the next time around the loop.
                                //
                                var newRegion = new Region(region);
                                region.end = coveredUpTo;
                                newRegion.start = coveredUpTo;
                                regions.Insert(index+1, newRegion);
                                regionsSize++;
                            }
                            else
                            {
                                //
                                // It starts exactly where we do.
                                //
                                if (region.end > end)
                                {
                                    //
                                    // But ends beyond where we do.  Split it and add ourself to the first half.
                                    //
                                    Region newRegion = new Region(region);
                                    newRegion.start = end;
                                    region.end = end;
                                    if (!region.genes.Contains(gene))
                                    {
                                        region.genes.Add(gene);
                                    }
                                    regions.Insert(index+1, newRegion);
                                    regionsSize++;
                                }
                                else
                                {
                                    //
                                    // Just add ourself to it.  It ends before us or exactly where we end.
                                    //
                                    if (!region.genes.Contains(gene))
                                    {
                                        region.genes.Add(gene);
                                    }
                                }
                            }
                            coveredUpTo = region.end;
                        } // else we weren't filling in the previously unoccupied first part


                    } // While we have more of this new exon to add

                } // AddMapping

                public readonly string name;
                public List<Region> regions = new List<Region>();
                public int regionsSize = 0; // So we don't have to use Count() all the time.
            }


            public GeneMap(string gencodeFilename)
            {
                var reader = new StreamReader(new GZipStream(new StreamReader(gencodeFilename).BaseStream, CompressionMode.Decompress));

                string line;
                while (null != (line = reader.ReadLine()))
                {
                    if (line.Count() == 0 || line[0] == '#') {
                        //
                        // Blank or comment/header line.  Skip.
                        //
                        continue;
                    }

                    var fields = line.Split('\t');
                    if (fields.Count() != 9)
                    {
                        Console.WriteLine("Malformed gencode line " + line);
                        continue;
                    }

                    string chromosome = fields[0];

                    if (!chromosomes.ContainsKey(chromosome))
                    {
                        chromosomes.Add(chromosome, new Chromosome(chromosome));
                    }

                    int start = Convert.ToInt32(fields[3]);
                    int end = Convert.ToInt32(fields[4]);

                    string mess = fields[8];
                    var messFields = mess.Split(';');
                    if (messFields.Count() < 5)
                    {
                        Console.WriteLine("too few fields in giant blob at end of gencode line " + line);
                        continue;
                    }

                    var geneNameFields = messFields[4].Split('\"'); // That's splitting by double quote.  Field should have the form gene_name "<gene name>"
                    if (geneNameFields.Count() != 3 || geneNameFields[0] != " gene_name " || geneNameFields[2] != "")
                    {
                        Console.WriteLine("Unparsable gene_name field in genecode line " + line);
                        continue;
                    }

                    string geneName = geneNameFields[1];

                    chromosomes[chromosome].AddMapping(geneName, start, end);
                }

                reader.Close();
            }

            public List<string> GetGenesCovering(string chromosomeName, int location)
            {
                chromosomeName = chromosomeName.ToLower();
                if (chromosomeName == "m" || chromosomeName == "mt" || chromosomeName == "chrm" || chromosomeName == "chrmt") {
                    //
                    // Just say no to mitochondrial genes.
                    //
                    return emptyList;
                }


                if (!chromosomes.ContainsKey(chromosomeName))
                {
                    //
                    // Try adding a leading "chr" if it doesn't have one
                    //
                    if (chromosomeName.Count() < 4 || chromosomeName.Substring(0, 3) != "chr")
                    {
                        if (chromosomes.ContainsKey("chr" + chromosomeName))
                        {
                            chromosomeName = "chr" + chromosomeName;
                        }
                    }

                    return emptyList;
                }

                Chromosome chromosome = chromosomes[chromosomeName];

                matchingRegion.start = location;
                int index = chromosome.regions.BinarySearch(matchingRegion);
                if (index < 0)
                {
                    //
                    // This means it wasn't found exactly.  index is the bitwise complement of the index of the next larger item.
                    // Since we want the next lower one, complement it and then subtract one (which might leave it negative).
                    //
                    index = ~index;
                    index -= 1;

                    if (index < 0)
                    {
                        return emptyList;
                    }
                }

                if (chromosome.regions[index].end <= location)
                {
                    return emptyList;
                }

                return chromosome.regions[index].genes;
            }

            Dictionary<string, Chromosome> chromosomes = new Dictionary<string, Chromosome>();
            readonly List<string> emptyList = new List<string>(); // We use this to return an answer when there's no matching chromosome or no gene mapped to that offset.
            Region matchingRegion = new Region(0);
        }




    } // ExpressionTools

}
