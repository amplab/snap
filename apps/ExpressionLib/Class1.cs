using System;
using System.Collections;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using System.IO;
using System.Xml.Linq;
using System.IO.Compression;
using System.Diagnostics;
using MathNet.Numerics;
using System.Threading;
using System.Runtime.InteropServices;

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

        public const int nHumanNuclearChromosomes = 24;   // 1-22, X and Y.

        public static int ChromosomeNameToIndex(string chromosomeName)
        {
            var name = chromosomeName.ToLower();

            for (int i = 1; i <= 22; i++)
            {
                var comp = Convert.ToString(i);
                if (name == comp || name == "chr" + comp) return i;
            }

            //
            // Use 0 and 23 for the sex chromosomes, so the array indices correspond to the chromosome names.
            //
            if (name == "chrx" || name == "x") return 0;
            if (name == "chry" || name == "y") return 23;

            return -1;
        }

        public static string ChromosomeIndexToName(int index, bool useChr = false)
        {
            if (index >= 1 && index <= 22)
            {
                return (useChr ? "chr" : "") + index;
            }

            if (index == 0) return (useChr ? "chr" : "") + "X";

            if (index == 23) return (useChr ? "chr" : "") + "Y";

            return "Invalid chromosome index: " + index;
        }
        public static string ChromPrefixFromRefassem(string refassem)
        {
            refassem = refassem.ToLower();
            if (refassem == "NCBI36_BCCAGSC_variant".ToLower() || refassem == "grch37-lite" || refassem == "grch37" || refassem == "hg19_broad_variant" || refassem == "grch37-lite_wugsc_variant_1" || refassem == "grch37-lite_wugsc_variant_2"
                || refassem == "grch37-lite-+-hpv_redux-build" || refassem == "HS37D5".ToLower() || refassem == "NCBI36_WUGSC_variant".ToLower())
            {
               return "";
            }
            else
            {
                return "chr";
            }

        }

        public static string ChromPrefixFromRefassemChromosomeAndBam(StoredBAM bam, string chromosome, string refassem)
        {
            if (chromosome.Count() > 2 && chromosome.Substring(0, 2).ToLower() == "gl")
            {
                return "";  // So we don't generate chrGL000...
            }
            if (null == bam || !bam.chrStateKnown) {
                return ChromPrefixFromRefassem(refassem);
            }

            if (bam.usesChr)
            {
                return "chr";
            }

            return "";
        }

        public static string chromosomeNameToNonChrForm(string rawName)
        {
            if (rawName.Substring(0, 3).ToLower() != "chr")
            {
                return rawName.ToLower();
            }

            return rawName.Substring(3).ToLower();
        }

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

        public static string GetFileNameFromPathname(Pathname pathname, bool excludeExtension = false)
        {
            string lastComponent;
            if (pathname.LastIndexOf('\\') == -1)
            {
                lastComponent = pathname;
            }
            else
            {
                lastComponent = pathname.Substring(pathname.LastIndexOf('\\') + 1);
            }

            if (!excludeExtension || lastComponent.LastIndexOf('.') == -1)
            {
                return lastComponent;
            }

            return lastComponent.Substring(0, lastComponent.LastIndexOf('.'));
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
            public FileInfo selectedVariantsInfo = null;
            public FileInfo allCountInfo = null;
            public FileInfo regionalExpressionInfo = null;
            public FileInfo geneExpressionInfo = null;
            public FileInfo readsAtSelectedVariantsInfo = null;
            public FileInfo readsAtSelectedVariantsIndexInfo = null;
            public FileInfo annotatedSelectedVariantsInfo = null;
            public FileInfo alleleSpecificGeneExpressionInfo = null;
            public FileInfo geneCoverageInfo = null;
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

            public ulong DiskSpaceUsed()
            {
                ulong totalSpace = 0;

                if (bamInfo != null) totalSpace += (ulong)bamInfo.Length;
                if (baiInfo != null) totalSpace += (ulong)baiInfo.Length;
                if (vcfInfo != null) totalSpace += (ulong)vcfInfo.Length;
                if (selectedVariantsInfo != null) totalSpace += (ulong)selectedVariantsInfo.Length;
                if (allCountInfo != null) totalSpace += (ulong)allCountInfo.Length;
                if (regionalExpressionInfo != null) totalSpace += (ulong)regionalExpressionInfo.Length;
                if (geneExpressionInfo != null) totalSpace += (ulong)geneExpressionInfo.Length;
                if (readsAtSelectedVariantsInfo != null) totalSpace += (ulong)readsAtSelectedVariantsInfo.Length;
                if (readsAtSelectedVariantsIndexInfo != null) totalSpace += (ulong)readsAtSelectedVariantsIndexInfo.Length;
                if (annotatedSelectedVariantsInfo != null) totalSpace += (ulong)annotatedSelectedVariantsInfo.Length;
                if (alleleSpecificGeneExpressionInfo != null) totalSpace += (ulong)alleleSpecificGeneExpressionInfo.Length;

                return totalSpace;
            }
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

        //
        // For filenames of the from ...\analysis-id\filename
        //
        static public string GetAnalysisIdFromPathname(string pathname)
        {
            if (!pathname.Contains('\\'))
            {
                Console.WriteLine("ExpressionTools.GetAnalysisIdFromPathname: invalid pathname input " + pathname);
                return "";
            }

            string directoryName;
            if (pathname.LastIndexOf('\\') != pathname.Count())
            {
                //
                // Does not end in a trailing \.
                //
                directoryName = pathname.Substring(0, pathname.LastIndexOf('\\'));
            }
            else
            {
                //
                // Ends in a trailing backslash, just trim it off.
                //
                directoryName = pathname.Substring(0, pathname.Count() - 1);
            }

            if (directoryName.Contains('\\'))
            {
                directoryName = directoryName.Substring(directoryName.LastIndexOf('\\') + 1);
            }

            if (directoryName.Count() != 36)
            {
                Console.WriteLine("ExpressionTools.GetAnalysisIdFromPathname: invalid pathname input does not include analysis id as last directory " + pathname);
                return "";
            }


            return directoryName;
        }

        public const string bamExtension = ".bam";
        public const string baiExtension = ".bai";
        public const string tarExtension = ".tar";
        public const string md5Extension = ".md5";
        public const string vcfExtension = ".vcf";
        public const string isoformsExtension = "-isoforms.txt";
        public const string allcountsExtension = ".allcount.gz";
        public const string regionalExpressionExtension = ".regional_expression.txt";
        public const string geneExpressionExtension = ".gene_expression.txt";
        public const string selectedVariantsExtension = ".selectedVariants";
        public const string annotatedSelectedVariantsExtension = ".annotatedSelectedVariants";
        public const string readsAtSelectedVariantsExtension = ".reads-at-selected-variants";
        public const string readsAtSelectedVariantsIndexExtension = readsAtSelectedVariantsExtension + ".index";
        public const string alleleSpecificGeneExpressionExtension = ".allele_specific_gene_expression";
        public const string geneCoverageExtension = ".gene_coverage.txt";

        static readonly int bamExtensionLength = bamExtension.Count();
        static readonly int baiExtensionLength = baiExtension.Count();
        static readonly int tarExtensionLength = tarExtension.Count();
        static readonly int md5ExtensionLength = md5Extension.Count();
        static readonly int vcfExtensionLength = vcfExtension.Count();
        static readonly int isoformsExtensionLength = isoformsExtension.Count();
        static readonly int allcountsExtensionLength = allcountsExtension.Count();
        static readonly int regionalExpressionExtensionLength = regionalExpressionExtension.Count();
        static readonly int geneExpressionExtensionLength = geneExpressionExtension.Count();
        static readonly int selectedVariantsExtensionLength = selectedVariantsExtension.Count();
        static readonly int annotatedSelectedVariantsExtensionLength = annotatedSelectedVariantsExtension.Count();
        static readonly int readsAtSelectedVariantsExtensionLength = readsAtSelectedVariantsExtension.Count();
        static readonly int readsAtSelectedVariantsIndexExtensionLength = readsAtSelectedVariantsIndexExtension.Count();
        static readonly int alleleSpecificGeneExpressionExtensionLength = alleleSpecificGeneExpressionExtension.Count();
        static readonly int geneCoverageExtensionLength = geneCoverageExtension.Count();

        static public bool CheckFileForDone(string filename, bool compressedFile) // Check that the given file is a text file who's last line is **done**
        {

            StreamReader reader = null;

            try
            {
                if (!compressedFile)
                {
                    var fileStream = new FileStream(filename, FileMode.Open);
                    if (fileStream.Length > 10 * 1024) {
                        fileStream.Seek(-10 * 1024, SeekOrigin.End);
                    }
                    reader = new StreamReader(fileStream);
                }
                else
                {
                    reader = new StreamReader(new GZipStream(new StreamReader(filename).BaseStream, CompressionMode.Decompress));
                    
                }

                bool seenDone = false;

                string line;
                while (null != (line = reader.ReadLine())) {
                    if (seenDone)
                    {
                        Console.WriteLine("File " + filename + " continues beyond **done**");
                        return false;
                    }

                    if ("**done**" == line || "**done**\t\t" == line) {
                        seenDone = true;
                    }
                }

                if (!seenDone)
                {
                    Console.WriteLine("File " + filename + " is truncated.");
                    return false;
                }
            }
            catch (IOException)
            {
                Console.WriteLine("Unable to to open file " + filename + " to check to see if it's truncated.  It's probably being written now, so assuming it's OK");
                return true;
            }

            return true;
        }


        static int LoadStoredBAMsForDirectory(Pathname directory, Dictionary<AnalysisID, StoredBAM> storedBAMs)
        {
            int nBamsLoaded = 0;

            if (!Directory.Exists(directory))
            {
                Console.WriteLine("Unable to find expected tcga data directory " + directory);
                return nBamsLoaded;
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
                        if (file.Count() > bamExtensionLength && file.Substring(file.Count() - bamExtensionLength) == bamExtension)
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
                            string md5Pathname = file + md5Extension;
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
                        else if (file.Count() > baiExtensionLength && file.Substring(file.Count() - baiExtensionLength) == baiExtension)
                        {
                            if (null != storedBAM.baiInfo)
                            {
                                Console.WriteLine("Saw multiple BAI files in the same analysis directory " + file + " and " + storedBAM.baiInfo.FullName);
                            }
                            storedBAM.baiInfo = new FileInfo(file);
                            string md5Pathname = file + md5Extension;

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
                        else if (file.Count() > tarExtensionLength && file.Substring(file.Count() - tarExtensionLength) == tarExtension)
                        {
                            hasTarFile = true;
                        }
                        else if (file.Count() > md5ExtensionLength && file.Substring(file.Count() - md5ExtensionLength) == md5Extension)
                        {
                            // Do nothing; checksum files are expected.
                        }
                        else if (file == subdir + @"\" + analysisID + vcfExtension)
                        {
                            storedBAM.vcfInfo = new FileInfo(file);
                            if (storedBAM.vcfInfo.Length < 1024 * 1024)
                            {
                                // This isn't my format and doesn't have a **done** terminator, so I can't check it, just report.

                                Console.WriteLine("Unusually small vcf file: " + file + ", " + storedBAM.vcfInfo.Length);
                            }
                        }
                        else if (file == subdir + @"\" + analysisID + selectedVariantsExtension)
                        {
                            storedBAM.selectedVariantsInfo = new FileInfo(file);

                            if (storedBAM.selectedVariantsInfo.Length < 1024 * 1024)
                            {
                                if (!CheckFileForDone(file, false))
                                {
                                    Console.WriteLine("Truncated (or empty) selected variants file " + file);
                                    storedBAM.selectedVariantsInfo = null;
                                }
                            }
                        }
                        else if (file == subdir + @"\" + analysisID + annotatedSelectedVariantsExtension)
                        {
                            storedBAM.annotatedSelectedVariantsInfo = new FileInfo(file);
                            if (storedBAM.annotatedSelectedVariantsInfo.Length < 1024 * 1024 && !CheckFileForDone(file, false))
                            {
                                Console.WriteLine("Truncated annotated selected variants file " + file);
                                storedBAM.annotatedSelectedVariantsInfo = null;
                            }
                        }
                        else if (file.Count() > isoformsExtensionLength && file.Substring(file.Count() - isoformsExtensionLength).ToLower() == isoformsExtension)
                        {
                            string expectedFilename = analysisID + isoformsExtension;
                            if (file.Count() < expectedFilename.Count() || file.Substring(file.Count() - expectedFilename.Count()) != expectedFilename)
                            {
                                Console.WriteLine("Incorrect isoform file " + file);
                            }
                            else
                            {
                                storedBAM.isoformCount = file;
                            }
                        }
                        else if (file.Count() > allcountsExtensionLength && file.Substring(file.Count() - allcountsExtensionLength).ToLower() == allcountsExtension)
                        {
                            string expectedFilename = analysisID + allcountsExtension;
                            if (file.Count() < expectedFilename.Count() || file.Substring(file.Count() - expectedFilename.Count()) != expectedFilename)
                            {
                                Console.WriteLine("Incorrect allcount file " + file);
                            }
                            else
                            {
                                storedBAM.allCountInfo = new FileInfo(file);
                                if (storedBAM.allCountInfo.Length < 1 * 1024 * 1024 && !CheckFileForDone(file, true))
                                {
                                    storedBAM.allCountInfo = null;
                                }
                            }
                        }
                        else if (file.Count() > regionalExpressionExtensionLength && file.Substring(file.Count() - regionalExpressionExtensionLength).ToLower() == regionalExpressionExtension)
                        {
                            string expectedFileName = analysisID + regionalExpressionExtension;
                            if (file.Count() < expectedFileName.Count() || file.Substring(file.Count() - expectedFileName.Count()) != expectedFileName)
                            {
                                Console.WriteLine("Incorrect regional expression file " + file);
                            }
                            else
                            {
                                storedBAM.regionalExpressionInfo = new FileInfo(file);
                                if (storedBAM.regionalExpressionInfo.Length < 1 * 1024 * 1024 && !CheckFileForDone(file, false))
                                {
                                    Console.WriteLine("Truncated small regional expression file " + file + ", " + storedBAM.regionalExpressionInfo.Length);
                                    storedBAM.regionalExpressionInfo = null;
                                }
                            }
                        }
                        else if (file.Count() > geneExpressionExtensionLength && file.Substring(file.Count() - geneExpressionExtensionLength).ToLower() == geneExpressionExtension)
                        {
                            string expectedFileName = analysisID + geneExpressionExtension;
                            if (file.Count() < expectedFileName.Count() || file.Substring(file.Count() - expectedFileName.Count()).ToLower() != expectedFileName.ToLower())
                            {
                                Console.WriteLine("Incorrect gene expression file " + file + ", expected filename " + expectedFileName);
                            }
                            else
                            {
                                storedBAM.geneExpressionInfo = new FileInfo(file);
                                if (/*storedBAM.geneExpressionInfo.Length < 1 * 1024 * 1024 &&*/ !CheckFileForDone(file, false))
                                {
                                    Console.WriteLine("Truncated gene expression file " + file + ", " + storedBAM.geneExpressionInfo.Length);
                                    storedBAM.geneExpressionInfo = null;
                                }
                            }
                        }
                        else if (file.Count() > alleleSpecificGeneExpressionExtensionLength && file.Substring(file.Count() - alleleSpecificGeneExpressionExtensionLength).ToLower() == alleleSpecificGeneExpressionExtension)
                        {
                            string expectedFileName = analysisID + alleleSpecificGeneExpressionExtension;
                            if (file.Count() < expectedFileName.Count() || file.Substring(file.Count() - expectedFileName.Count()) != expectedFileName)
                            {
                                Console.WriteLine("Incorrect allele-specific gene expression file " + file);
                            }
                            else
                            {
                                storedBAM.alleleSpecificGeneExpressionInfo = new FileInfo(file);
                                if (storedBAM.alleleSpecificGeneExpressionInfo.Length < 1 * 1024 * 1024)
                                {
                                    //
                                    // Check it to see if it's truncated.
                                    //
                                    if (!CheckFileForDone(file, false)) 
                                    {
                                        Console.WriteLine("Truncated allele-specific gene expression file " + file + ", size " + storedBAM.alleleSpecificGeneExpressionInfo.Length);
                                        storedBAM.alleleSpecificGeneExpressionInfo = null;
                                    }
                                }
                            }
                        }
                        else if (file.Count() > readsAtSelectedVariantsExtensionLength && file.Substring(file.Count() - readsAtSelectedVariantsExtensionLength).ToLower() == readsAtSelectedVariantsExtension)
                        {
                            string expectedFileName = analysisID + readsAtSelectedVariantsExtension;
                            if (file.Count() < expectedFileName.Count() || file.Substring(file.Count() - expectedFileName.Count()) != expectedFileName)
                            {
                                Console.WriteLine("Incorrect reads at expected variants file " + file);
                            }
                            else
                            {
                                storedBAM.readsAtSelectedVariantsInfo = new FileInfo(file);
                            }
                        }
                        else if (file.Count() > readsAtSelectedVariantsIndexExtensionLength && file.Substring(file.Count() - readsAtSelectedVariantsIndexExtensionLength).ToLower() == readsAtSelectedVariantsIndexExtension)
                        {
                            string expectedFileName = analysisID + readsAtSelectedVariantsIndexExtension;
                            if (file.Count() < expectedFileName.Count() || file.Substring(file.Count() - expectedFileName.Count()) != expectedFileName)
                            {
                                Console.WriteLine("Incorrect reads at expected variants index file " + file);
                            }
                            else
                            {
                                storedBAM.readsAtSelectedVariantsIndexInfo = new FileInfo(file);
                                if (storedBAM.readsAtSelectedVariantsIndexInfo.Length < 1024 * 1024 && !CheckFileForDone(file, false))
                                {
                                    Console.WriteLine("File " + file + " is truncated.");
                                    storedBAM.readsAtSelectedVariantsIndexInfo = null;
                                }
                            }
                        }
                        else if (file.Count() > geneCoverageExtensionLength && file.Substring(file.Count() - geneCoverageExtensionLength).ToLower() == geneCoverageExtension)
                        {
                            string expectedFileName = analysisID + geneCoverageExtension;
                            if (file.Count() < expectedFileName.Count() || file.Substring(file.Count() - expectedFileName.Count()) != expectedFileName)
                            {
                                Console.WriteLine("Incorrect gene coverage file " + file);
                            }
                            else
                            {
                                storedBAM.geneCoverageInfo = new FileInfo(file);
                                if (storedBAM.geneCoverageInfo.Length < 1024 * 1024 && !CheckFileForDone(file, false))
                                {
                                    Console.WriteLine("File " + file + " is truncated.");
                                    storedBAM.geneCoverageInfo = null;
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

                    if ((storedBAM.readsAtSelectedVariantsInfo == null) != (storedBAM.readsAtSelectedVariantsIndexInfo == null))
                    {
                        Console.WriteLine("Analysis directory " + subdir + " has exactly one of the .reads-at-selected-variants and reads-at-selected-variants.index files.");
                        storedBAM.readsAtSelectedVariantsInfo = null;
                        storedBAM.readsAtSelectedVariantsIndexInfo = null;
                    }

                    storedBAM.totalSize = storedBAM.bamInfo.Length + storedBAM.baiInfo.Length;
                    storedBAM.directoryName = subdir;

                    lock (storedBAMs) {
                        if (storedBAMs.ContainsKey(analysisID))
                        {
                            Console.WriteLine("Duplicate stored BAM, " + subdir + " and " + storedBAMs[analysisID].directoryName);
                        }
                        else
                        {
                            storedBAMs.Add(analysisID, storedBAM);
                            nBamsLoaded++;
                        }
                    }
                }
                else if (analysisID.Contains(".partial"))
                {
                    Console.WriteLine("Partial download at " + subdir);
                }
            }

            return nBamsLoaded;
        }

        public static string [] TumorAndNormal = {"tumor", "normal"};
        public static string[] LibraryStrategies = { "wgs", "wxs", "rna" };


        public static void LoadStoredBamsForMachine(Pathname directory, Dictionary<AnalysisID, StoredBAM> storedBAMs, ref ulong totalFreeDiskSpace)
        {
            //
            // Load BAMs from a machine that offloads storage from the main machines.  These machines are allowed to store data
            // in places that's other than its primary storage location.  We just work though the directory hierarchy and
            // load each lowest level directory individually.  "directory" should point at the directory containing the
            // disease subdirectories, i.e., \\msr-genomics-0\d$\tcga.
            //

            int nBamsLoaded = 0;
            var stopwatch = new Stopwatch();
            stopwatch.Start();

            if (!Directory.Exists(directory))
            {
                Console.WriteLine("Unable to find expected secondary machine directory " + directory);
                return;
            }

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
                            nBamsLoaded += LoadStoredBAMsForDirectory(sampleDir, storedBAMs);
                        }
                    }
                }
            }

            stopwatch.Stop();

            lock (storedBAMs)
            {
                var freeDiskSpace = GetVolumeFreeSpace(directory);
                totalFreeDiskSpace += freeDiskSpace;
                Console.WriteLine("Loaded " + nBamsLoaded + " bams for " + directory + " in " + (stopwatch.ElapsedMilliseconds + 500) / 1000 + "s, " + (nBamsLoaded * 1000 / (stopwatch.ElapsedMilliseconds + 1)) + " bams/s, "
                    + SizeToUnits(freeDiskSpace) + "B free");
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
                line = line.ToLower();
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
                else if (fields[1] == "nochr")
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
            public Pathname allcountFileName = null;
            public Pathname regionalExpressionFileName = null;
            public Pathname geneExpressionFileName = null;
            public Pathname alleleSpecificGeneExpressionFileName = null;
            public Pathname selectedVariantsFileName = null;
            public Pathname readsAtSelectedVariantsFileName = null;
            public Pathname annotatedSelectedVariantsFileName = null;
            public Pathname geneCoverageFileName = null;
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
            public Participant participant = null;
            public string disease_abbr;
            public TCGARecord TumorRNAAnalysis;
            public TCGARecord TumorDNAAnalysis;
            public TCGARecord NormalDNAAnalysis;
            public TCGARecord NormalRNAAnalysis = null;    // This is optional, and filled in if it exists, but we'll generate the experiment even if it doesn't.
            public List<MAFRecord> maf;
            public bool normalNeedsRealignment = false;
            public bool tumorNeedsRealignment = false;
        }

        public const int DiseaseAbbrFieldNumber = 0;
        public const int ReferenceFieldNumber = 1;
        public const int ParticipantIDFieldNumber = 2;
        public const int TumorRNAAnalysisFieldNumber = 3;
        public const int TumorDNAAnalysisFieldNumber = 4;
        public const int NormalDNAAnalysisFieldNumber = 5;
        public const int NormalRNAAnalysisFieldNumber = 6;
        public const int TumorRNAPathnameFieldNumber = 7;
        public const int TumorDNAPathnameFieldNumber = 8;
        public const int NormalDNAPathnameFieldNumber = 9;
        public const int NormalRNAPathnameFieldNumber = 10;
        public const int VCFPathnameFieldNumber = 11;
        public const int GenderFieldNumber = 12;
        public const int DaysToBirthFieldNumber = 13;
        public const int DaysToDeathFieldNumber = 14;
        public const int OrigTumorDNAAliquotIDFieldNumber = 15;
        public const int TumorRNAAllcountFileFieldNumber = 16;
        public const int NormalRNAAllcountFileFieldNumber = 17;
        public const int MAFFileFieldNumber = 18;
        public const int RegionallExpressionFileFieldNumber = 19;
        public const int GeneExpressionFileFieldNumber = 20;
        public const int SelectedVariantsFileFieldNumber = 21;
        public const int ReadsAtSelectedVariantsDNAFileFieldNumber = 22;
        public const int ReadsAtSelectedVariantsRNAFileFieldNumber = 23;
        public const int AnnotatedSelectedVariantsFieldNumber = 24;
        public const int SourceNormalDNAFileFieldNumber = 25;
        public const int SourceTumorDNAFileFieldNumber = 26;
        public const int TumorDNAAllcountFileFieldNumber = 27;
        public const int AlleleSpecificGeneExpressionFileFieldNumber = 28;
        public const int ReadsAtSelectedVariantsNormalRNAFileFieldNumber = 29;
        public const int TumorDNAGeneCoverageFileFieldNumber = 30;

        public static void DumpExperimentsToFile(List<Experiment> experiments, string filename)
        {
            StreamWriter outputFile = null;

            bool threw;
            do
            {
                threw = false;
                try
                {
                    outputFile = new StreamWriter(filename);
                }
                catch (IOException)
                {
                    Console.WriteLine("IOException opening file " + filename + " for write.  Do you perhaps have it open in Excel?  Sleeping and retrying.");
                    Thread.Sleep(10000);
                    threw = true;
                }
            } while (threw);

            outputFile.WriteLine("disease_abbr\treference\tparticipantID\tTumorRNAAnalysis\tTumorDNAAnalysis\tNormalDNAAnalysis\tNormalRNAAnalysis\ttumorRNAPathname\ttumorDNAPathname\t" +
                "normalDNAPathname\tNormalRNAPathname\tVCFPathname\tgender\tdaysToBirth\tdaysToDeath\tOrigTumorDNAAliquotID\tTumorRNAAllcountFile\tNormalRNAAllcountFile\tmafFile\t" +
                "RegionalExpressionFilename\tGeneExpressionFilename\tSelectedVariantsFilename\tReadsAtSelectedVariantsDNAFilename\tReadsAtSelectedVariantsRNAFilename\tAnotatedSelectedVariantsFile\t" +
                "SourceNormalDNA\tSourceTumorDNA\tTumorDNAAllcountFilename\tAlleleSpecifcGeneExpressionFile\tNormal RNA ReadsAtSelectedVariants\tGene Coverage");

            foreach (Experiment experiment in experiments)
            {
                // 0-5
                outputFile.Write(experiment.disease_abbr + "\t" + experiment.TumorRNAAnalysis.refassemShortName + "\t" + experiment.participant.participantId + "\t" + experiment.TumorRNAAnalysis.analysis_id + "\t" + experiment.TumorDNAAnalysis.analysis_id + "\t" + experiment.NormalDNAAnalysis.analysis_id + "\t");

                // 6
                if (null != experiment.NormalRNAAnalysis)
                {
                    outputFile.Write(experiment.NormalRNAAnalysis.analysis_id + "\t");
                }
                else
                {
                    outputFile.Write("\t");
                }


                // 7
                if (experiment.TumorRNAAnalysis.storedBAM != null)
                {
                    outputFile.Write(experiment.TumorRNAAnalysis.storedBAM.bamInfo.FullName + "\t");
                }
                else
                {
                    outputFile.Write("\t");
                }

                // 8
                if (experiment.TumorDNAAnalysis.storedBAM != null)
                {
                    outputFile.Write(experiment.TumorDNAAnalysis.storedBAM.bamInfo.FullName + "\t");
                }
                else
                {
                    outputFile.Write("\t");
                }

                // 9
                if (experiment.NormalDNAAnalysis.storedBAM != null)
                {
                    outputFile.Write(experiment.NormalDNAAnalysis.storedBAM.bamInfo.FullName + "\t");
                }
                else
                {
                    outputFile.Write("\t");
                }

                // 10
                if (experiment.NormalRNAAnalysis != null && experiment.NormalRNAAnalysis.storedBAM != null)
                {
                    outputFile.Write(experiment.NormalRNAAnalysis.storedBAM.bamInfo.FullName + "\t");
                }
                else
                {
                    outputFile.Write("\t");
                }

                // 11
                if (experiment.NormalDNAAnalysis.storedBAM != null && experiment.NormalDNAAnalysis.storedBAM.vcfInfo != null)
                {
                    outputFile.Write(experiment.NormalDNAAnalysis.storedBAM.vcfInfo.FullName + "\t");
                }
                else
                {
                    outputFile.Write("\t");
                }

                // 12 - 14
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

                // 15
                if (experiment.TumorDNAAnalysis.realignSource != null)
                {
                    outputFile.Write(experiment.TumorDNAAnalysis.realignSource.aliquot_id + "\t");
                }
                else
                {
                    outputFile.Write(experiment.TumorDNAAnalysis.aliquot_id + "\t");
                }

                // 16
                if (experiment.TumorRNAAnalysis.storedBAM != null && experiment.TumorRNAAnalysis.storedBAM.allCountInfo != null)
                {
                    outputFile.Write(experiment.TumorRNAAnalysis.storedBAM.allCountInfo.FullName + "\t");
                }
                else
                {
                    outputFile.Write("\t");
                }

                // 17
                if (experiment.NormalRNAAnalysis != null && experiment.NormalRNAAnalysis.storedBAM != null && experiment.NormalRNAAnalysis.storedBAM.allCountInfo != null)
                {
                    outputFile.Write(experiment.NormalRNAAnalysis.storedBAM.allCountInfo.FullName + "\t");
                }
                else
                {
                    outputFile.Write("\t");
                }

                // 18
                outputFile.Write(experiment.participant.mafFile + "\t");

                // 19
                if (experiment.TumorRNAAnalysis.regionalExpressionFileName == null)
                {
                    outputFile.Write("\t");
                }
                else
                {
                    outputFile.Write(experiment.TumorRNAAnalysis.regionalExpressionFileName + "\t");
                }

                // 20
                if (experiment.TumorRNAAnalysis.geneExpressionFileName == null)
                {
                    outputFile.Write("\t");
                }
                else
                {
                    outputFile.Write(experiment.TumorRNAAnalysis.geneExpressionFileName + "\t");
                }

                // 21
                if (experiment.NormalDNAAnalysis.storedBAM != null && experiment.NormalDNAAnalysis.storedBAM.selectedVariantsInfo != null)
                {
                    outputFile.Write(experiment.NormalDNAAnalysis.storedBAM.selectedVariantsInfo.FullName + "\t");
                }
                else
                {
                    outputFile.Write("\t");
                }

                // 22
                if (experiment.TumorDNAAnalysis.storedBAM != null && experiment.TumorDNAAnalysis.storedBAM.readsAtSelectedVariantsInfo != null)
                {
                    outputFile.Write(experiment.TumorDNAAnalysis.storedBAM.readsAtSelectedVariantsInfo.FullName + "\t");
                }
                else
                {
                    outputFile.Write("\t");
                }

                // 23
                if (experiment.TumorRNAAnalysis.storedBAM != null && experiment.TumorRNAAnalysis.storedBAM.readsAtSelectedVariantsInfo != null)
                {
                    outputFile.Write(experiment.TumorRNAAnalysis.storedBAM.readsAtSelectedVariantsInfo.FullName + "\t");
                }
                else
                {
                    outputFile.Write("\t");
                }

                // 24
                if (experiment.NormalDNAAnalysis.storedBAM != null && experiment.NormalDNAAnalysis.storedBAM.annotatedSelectedVariantsInfo != null)
                {
                    outputFile.Write(experiment.NormalDNAAnalysis.storedBAM.annotatedSelectedVariantsInfo.FullName + "\t");
                }
                else
                {
                    outputFile.Write("\t");
                }

                // 25
                if (experiment.NormalDNAAnalysis.realignSource != null && experiment.NormalDNAAnalysis.realignSource.storedBAM != null && experiment.NormalDNAAnalysis.realignSource.storedBAM.bamInfo != null)
                {
                    outputFile.Write(experiment.NormalDNAAnalysis.realignSource.storedBAM.bamInfo.FullName + "\t");
                } else {
                    outputFile.Write("\t");
                }

                // 26
                if (experiment.TumorDNAAnalysis.realignSource != null && experiment.TumorDNAAnalysis.realignSource.storedBAM != null && experiment.TumorDNAAnalysis.realignSource.storedBAM.bamInfo != null)
                {
                    outputFile.Write(experiment.TumorDNAAnalysis.realignSource.storedBAM.bamInfo.FullName + "\t");
                } else {
                    outputFile.Write("\t");
                }

                // 27
                if (experiment.TumorDNAAnalysis.storedBAM != null && experiment.TumorDNAAnalysis.storedBAM.allCountInfo != null)
                {
                    outputFile.Write(experiment.TumorDNAAnalysis.storedBAM.allCountInfo.FullName + "\t");
                }
                else
                {
                    outputFile.Write("\t");
                }

                // 28
                if (experiment.NormalDNAAnalysis.storedBAM == null || experiment.NormalDNAAnalysis.storedBAM.alleleSpecificGeneExpressionInfo == null)
                {
                    outputFile.Write("\t");
                }
                else
                {
                    outputFile.Write(experiment.NormalDNAAnalysis.storedBAM.alleleSpecificGeneExpressionInfo.FullName + "\t");
                }

                // 28
                if (experiment.NormalRNAAnalysis == null || experiment.NormalRNAAnalysis.storedBAM == null || experiment.NormalDNAAnalysis.storedBAM.readsAtSelectedVariantsInfo == null)
                {
                    outputFile.Write("\t");
                }
                else
                {
                    outputFile.Write(experiment.NormalDNAAnalysis.storedBAM.readsAtSelectedVariantsInfo.FullName + "\t");
                }

                // 29
                if (experiment.TumorDNAAnalysis == null || experiment.TumorDNAAnalysis.storedBAM == null || experiment.TumorDNAAnalysis.geneCoverageFileName == null)
                {
                    outputFile.Write("\t");
                }
                else
                {
                    outputFile.Write(experiment.TumorDNAAnalysis.geneCoverageFileName + "\t");
                }

                outputFile.WriteLine();
            }

            outputFile.Close();
        } // DumpExperimentsToFile

        public static List<Experiment> LoadExperimentsFromFile(string filename, Dictionary<ParticipantID, Participant> participants, Dictionary<AnalysisID, TCGARecord> tcgaRecords)
        {
            var experiments = new List<Experiment>();

            StreamReader reader = null;
            bool threw = false;
            do
            {
                threw = false;
                try
                {
                    reader = new StreamReader(filename);
                }
                catch (IOException)
                {
                    Console.WriteLine("IOException opening " + filename + ".  Most likely it's open in another program, probably Excel.  Sleeping 10s and retrying.");
                    Thread.Sleep(10000);
                    threw = true;
                }
            } while (threw);

            string line;
            reader.ReadLine();  // Skip the header line
            while (null != (line = reader.ReadLine())) {
                var fields = line.Split('\t');

                if (fields.Count() < 28)
                {
                    Console.WriteLine("LoadExperimentsFromFile::line doesn't have enough fields.  Perhaps it's from an old format or damaged experiments.txt file, line " + line + ", filename " + filename);
                    return null;
                }

                var experiment = new Experiment();
                experiment.disease_abbr = fields[DiseaseAbbrFieldNumber];
                if (null != participants)
                {
                    experiment.participant = participants[fields[ParticipantIDFieldNumber]];
                }
                experiment.TumorRNAAnalysis = tcgaRecords[fields[TumorRNAAnalysisFieldNumber]];
                experiment.TumorRNAAnalysis.bamFileName = fields[TumorRNAPathnameFieldNumber];
                experiment.TumorDNAAnalysis = tcgaRecords[fields[TumorDNAAnalysisFieldNumber]];
                experiment.TumorDNAAnalysis.bamFileName = fields[TumorDNAPathnameFieldNumber];
                experiment.NormalDNAAnalysis = tcgaRecords[fields[NormalDNAAnalysisFieldNumber]];
                experiment.NormalDNAAnalysis.bamFileName = fields[NormalDNAPathnameFieldNumber];
                if (fields[NormalRNAAnalysisFieldNumber] == "")
                {
                    experiment.NormalRNAAnalysis = null;
                } else {
                    experiment.NormalRNAAnalysis = tcgaRecords[fields[NormalRNAAnalysisFieldNumber]];
                    experiment.NormalRNAAnalysis.bamFileName = fields[NormalRNAPathnameFieldNumber];
                }
                if (null != experiment.participant && experiment.participant.mafs.Count() > 0)
                {
                    experiment.maf = experiment.participant.mafs[0];
                }
                else
                {
                    experiment.maf = null;
                }

                if (null != experiment.NormalRNAAnalysis)
                {
                    experiment.NormalRNAAnalysis.allcountFileName = fields[NormalRNAAllcountFileFieldNumber];
                }
                experiment.TumorRNAAnalysis.allcountFileName = fields[TumorRNAAllcountFileFieldNumber];
                experiment.TumorRNAAnalysis.regionalExpressionFileName = fields[RegionallExpressionFileFieldNumber];
                experiment.TumorRNAAnalysis.geneExpressionFileName = fields[GeneExpressionFileFieldNumber];
                experiment.NormalDNAAnalysis.selectedVariantsFileName = fields[SelectedVariantsFileFieldNumber];
                experiment.TumorDNAAnalysis.readsAtSelectedVariantsFileName = fields[ReadsAtSelectedVariantsDNAFileFieldNumber];
                experiment.TumorRNAAnalysis.readsAtSelectedVariantsFileName = fields[ReadsAtSelectedVariantsRNAFileFieldNumber];
                experiment.NormalDNAAnalysis.annotatedSelectedVariantsFileName = fields[AnnotatedSelectedVariantsFieldNumber];
                experiment.TumorDNAAnalysis.allcountFileName = fields[TumorDNAAllcountFileFieldNumber];
                experiment.NormalDNAAnalysis.alleleSpecificGeneExpressionFileName = fields[AlleleSpecificGeneExpressionFileFieldNumber];
                if (experiment.NormalRNAAnalysis != null)
                {
                    experiment.NormalRNAAnalysis.readsAtSelectedVariantsFileName = fields[ReadsAtSelectedVariantsNormalRNAFileFieldNumber];
                }
                experiment.TumorDNAAnalysis.geneCoverageFileName = fields[TumorDNAGeneCoverageFileFieldNumber];
                experiments.Add(experiment);
            }

            reader.Close();
            return experiments;
        }

        public class MAFRecord : IComparable<MAFRecord>
        {
            public string entire_maf_line;  // The raw line.

            public string ReferenceClass()
            {
                if (NcbiBuild == "36")
                {
                    return "hg18";
                }
                else
                {
                    return "hg19";
                }
            }

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

            public int CompareTo(MAFRecord peer) {
                if (Chrom != peer.Chrom)
                {
                    return Chrom.CompareTo(peer.Chrom);
                }

                return Start_position.CompareTo(peer.Start_position);
            }
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
            string[] inputLines = ReadAllLinesWithRetry(inputFilename);

            int tooShortLines = 0;

            //
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

                mafRecord.Hugo_symbol = ConvertToNonExcelString(fields[0]);
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

        public static void AddAllMAFFilesToParticipants(Dictionary<ParticipantID, Participant> participants, Dictionary<SampleID, ParticipantID> sampleToParticipantIDMap, string path = @"f:\sequence\Reads\tcga\mafs\")
        {
            foreach (var file in Directory.GetFiles(path))
            {
                AddMAFFileToParticipants(file, participants, sampleToParticipantIDMap);
            }
        }


        public static Dictionary<ParticipantID, Participant> BuildParticipantData(Dictionary<string, TCGARecord> tcgaEntries, out Dictionary<string, Sample> allSamples, string clinicalDirectory = @"f:\sequence\reads\tcga\clinical")
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
            foreach (var clinical_file in Directory.EnumerateFiles(clinicalDirectory))
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
        public static Dictionary<AnalysisID, ExpressionTools.TCGARecord> LoadTCGARecords(Dictionary<AnalysisID, ExpressionTools.StoredBAM> storedBAMs, List<AnalysisID> excludedAnalyses, string filename = @"f:\sequence\reads\tcga-all.xml")
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

                    if (record.storedBAM.allCountInfo != null)
                    {
                        record.allcountFileName = record.storedBAM.allCountInfo.FullName;
                    }

                    if (record.storedBAM.regionalExpressionInfo != null)
                    {
                        record.regionalExpressionFileName = record.storedBAM.regionalExpressionInfo.FullName;
                    }

                    if (record.storedBAM.geneExpressionInfo != null)
                    {
                        record.geneExpressionFileName = record.storedBAM.geneExpressionInfo.FullName;
                    }

                    if (record.storedBAM.selectedVariantsInfo != null) {
                        record.selectedVariantsFileName = record.storedBAM.selectedVariantsInfo.FullName;
                    }

                    if (record.storedBAM.readsAtSelectedVariantsInfo != null) {
                        record.readsAtSelectedVariantsFileName = record.storedBAM.readsAtSelectedVariantsInfo.FullName;
                    }

                    if (record.storedBAM.annotatedSelectedVariantsInfo != null) {
                        record.annotatedSelectedVariantsFileName = record.storedBAM.annotatedSelectedVariantsInfo.FullName;
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

        public static void LoadTCGAAdditionalMetadata(Dictionary<AnalysisID, ExpressionTools.TCGARecord> tcgaRecords, string additonalMetadataFilename = @"f:\sequence\reads\tcga\tcgaAdditionalMetadata.txt")
        {
            if (!File.Exists(additonalMetadataFilename))
            {
                return;
            }
            string[] additionalMetadata = File.ReadAllLines(additonalMetadataFilename);
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

        public static List<AnalysisID> LoadExcludedAnalyses(string filename = @"f:\sequence\reads\tcga\excluded_analyses.txt")
        {
            var excludedAnalyses = new List<AnalysisID>();

            string[] analyses = File.ReadAllLines(filename);

            foreach (var analysis in analyses)
            {
                excludedAnalyses.Add(analysis.ToLower());
            }

            return excludedAnalyses;
        }

        public static void LoadTCGARecordsForLocalRealigns(Dictionary<AnalysisID, ExpressionTools.TCGARecord> tcgaRecords, Dictionary<AnalysisID, ExpressionTools.StoredBAM> storedBAMs, string filename = @"f:\sequence\reads\tcga\realigns.txt")
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
                        else
                        {
                            return emptyList;
                        }
                    }
                }
                else
                {
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


        public struct MeanAndStdDev
        {
            public double mean;
            public double stddev;
        }

        //
        // Maps chromosomeName -> (offset -> MeanAndStdDev)
        public static Dictionary<string, Dictionary<int, MeanAndStdDev>> LoadExpressionFile(string filename)
        {
            var expression = new Dictionary<string, Dictionary<int, ExpressionTools.MeanAndStdDev>>();

            StreamReader reader = new StreamReader(filename);

            string currentChromosome = null;
            Dictionary<int, ExpressionTools.MeanAndStdDev> currentMapping = null;

            string statLine;
            while (null != (statLine = reader.ReadLine()))
            {
                string[] statLineFields = statLine.Split('\t');
                if (statLineFields.Count() == 1)
                {
                    //
                    // It's a chromosome header.
                    //
                    currentChromosome = statLine.ToLower();
                    if (currentChromosome.Substring(0, 3) == "chr") // Clip off the leading "chr," if it exists
                    {
                        currentChromosome = currentChromosome.Substring(3);
                    }
                    currentMapping = new Dictionary<int, MeanAndStdDev>();
                    expression.Add(currentChromosome, currentMapping);
                    continue;
                }

                MeanAndStdDev stats = new MeanAndStdDev();
                int chromosomeOffset = Convert.ToInt32(statLineFields[0]);
                stats.mean = Convert.ToDouble(statLineFields[2]);
                stats.stddev = Convert.ToDouble(statLineFields[3]);

                if (stats.mean > 0 && stats.stddev > 0)
                {
                    currentMapping.Add(chromosomeOffset, stats);
                }
            }

            return expression;
        }

        public static Dictionary<string, Experiment> BuildParticipantToExperimentMapping(List<Experiment> experiments)
        {
            var mapping = new Dictionary<string, Experiment>();

            foreach (var experiment in experiments)
            {
                mapping.Add(experiment.participant.participantId, experiment);
            }

            return mapping;
        }

        public static double MeanOfList(List<double> values)
        {
            return values.Sum() / values.Count();
        }

        public static double StandardDeviationOfList(List<double> values)
        {
            var mean = MeanOfList(values);

            double variance = 0;
            foreach (var value in values)
            {
                var difference = value - mean;
                variance = difference * difference;
            }

            return Math.Sqrt(variance);
        }

        public class MannWhitney<T>
        {
            public delegate bool WhichGroup(T element);
            public delegate double GetValue(T element);
            public static double ComputeMannWhitney(List<T> elements, IComparer<T> comparer, WhichGroup whichGroup, GetValue getValue, out bool enoughData, out bool reversed, 
                out double nFirstGroup, out double nSecondGroup, out double U, out double z, bool twoTailed = true, int minGroupSize = 1)
            {
                elements.Sort(comparer);

                reversed = false;

                double RfirstGroup = 0; // The rank sum for the first group
                nFirstGroup = 0;
                nSecondGroup = 0;
                U = 0;
                z = 0;
                int n = elements.Count();

                for (int i = 0; i < n; i++)
                {
                    if (whichGroup(elements[i]))
                    {
                        int cumulativeR = n - i;
                        int nTied = 1;
                        //
                        // Now add in adjascent indices if there's a tie.  For ties, we use the mean of all the indices in the tied region for each of them (regardless of whether they're single or multiple).
                        //
                        for (int j = i - 1; j >= 0 && getValue(elements[j]) == getValue(elements[i]); j--)
                        {
                            cumulativeR += n - j;
                            nTied++;
                        }

                        for (int j = i + 1; j < elements.Count() && getValue(elements[j]) == getValue(elements[i]); j++)
                        {
                            cumulativeR += n - j;
                            nTied++;
                        }


                        RfirstGroup += cumulativeR / nTied;
                        nFirstGroup++;
                    }
                    else
                    {
                        nSecondGroup++;
                    }
                }

                if (nFirstGroup < minGroupSize || nSecondGroup < minGroupSize)
                {
                    //
                    // Not enough data, reject this gene
                    //
                    enoughData = false;
                    return -1;
                }

                U = (double)RfirstGroup - (double)nFirstGroup * (nFirstGroup + 1.0) / 2.0;

                z = (U - nFirstGroup * nSecondGroup / 2) / Math.Sqrt(nSecondGroup * nFirstGroup * (nSecondGroup + nFirstGroup + 1) / 12);

                double p = MathNet.Numerics.Distributions.Normal.CDF(0, 1.0, z);

                if (twoTailed)
                {
                    //
                    // We're doing two-tailed, so we need to see if this is on the other end
                    if (p > 0.5)
                    {
                        p = 1.0 - p;
                        reversed = true;
                    }

                    //
                    // And then multiply by two (because this is a two-tailed test).
                    //
                    p *= 2.0;
                }

                enoughData = true;
                return p;
            }
        }

        static public List<string> GetListOfDiseases(List<Experiment> experiments)
        {
            var listOfDiseases = new List<string>();

            foreach (var experiment in experiments)
            {
                if (!listOfDiseases.Contains(experiment.disease_abbr))
                {
                    listOfDiseases.Add(experiment.disease_abbr);
                }
            }

            return listOfDiseases;
        }

        //
        // Code for GetVolumeFreeSpace adapted from http://stackoverflow.com/questions/14465187/get-available-disk-free-space-for-a-given-path-on-windows
        //
        [DllImport("kernel32.dll", SetLastError = true, CharSet = CharSet.Auto)]
        [return: MarshalAs(UnmanagedType.Bool)]
        static extern bool GetDiskFreeSpaceEx(string lpDirectoryName,
           out ulong lpFreeBytesAvailable,
           out ulong lpTotalNumberOfBytes,
           out ulong lpTotalNumberOfFreeBytes);

        public static ulong GetVolumeFreeSpace(string pathname)
        {
            ulong FreeBytesAvailable;
            ulong TotalNumberOfBytes;
            ulong TotalNumberOfFreeBytes;

            if (!GetDiskFreeSpaceEx(pathname,
                                        out FreeBytesAvailable,
                                        out TotalNumberOfBytes,
                                        out TotalNumberOfFreeBytes))
            {
                throw new System.ComponentModel.Win32Exception();
            }

            return FreeBytesAvailable;
        }

        public static string SizeToUnits(ulong size)
        {
            if (size < 1024) return "" + size;

            if (size < 1024 * 1024) return "" + ((size + 512) / 1024) + "K";

            if (size < 1024 * 1024 * 1024) return "" + ((size + 512 * 1024) / (1024 * 1024)) + "M";

            if (size < (ulong)1024 * 1024 * 1024 * 1024) return "" + ((size + 512 * 1024 * 1024) / (1024 * 1024 * 1024)) + "G";

            if (size < (ulong)1024 * 1024 * 1024 * 1024 * 1024) return "" + ((size + (ulong)512 * 1024 * 1204 * 1024) / ((ulong)1024 * 1024 * 1024 * 1024)) + "T";

            return "" + ((size + (ulong)512 * 1024 * 1204 * 1024 * 1024) / ((ulong)1024 * 1024 * 1024 * 1024 * 1024)) + "P";
        }

        //
        // The BeatAML dataset doesn't have TCGA records or GUIDs for its participants and analyses.  Generate them together with filename mappings.
        //
        const string BeatAMLBaseDirectory = @"\\msr-genomics-1\e$\BeatAML\";
        const string BeatAMLGuidsFile = @"f:\sequence\reads\tcga\BeatAMLGuids";

        public static string BeatAMLFilenameToSampleName(string filename) 
        {
            //
            // First strip off the directory, if it's there.
            //
            string workingString;

            if (filename.Contains('\\'))
            {
                workingString = filename.Substring(filename.LastIndexOf('\\') + 1).ToLower();
            }
            else
            {
                workingString = filename.ToLower();
            }

            //
            // Now strip off anything after the second _.
            //

            if (!workingString.Contains('_'))
            {
                throw new FormatException();
            }

            var indexOfFirstUnderscore = workingString.IndexOf('_');
            if (workingString.Substring(indexOfFirstUnderscore + 1).Contains('_'))
            {
                workingString = workingString.Substring(0, workingString.Substring(indexOfFirstUnderscore + 1).IndexOf('_') + indexOfFirstUnderscore + 1);
            }

            if (workingString.Count() != 15 || workingString.Substring(0, 7) != "sample_" || workingString[9] != '-')
            {
                throw new FormatException();
            }

            return workingString;
        }

        public class BeatAMLSample
        {
            public BeatAMLSample(string sampleName_) {
                sampleName = sampleName_;
            }

            public string sampleName;
            public string normalDNAGuid = null;
            public string tumorDNAGuid = null;
            public string tumorRNAGuid = null;
        }

        static public void GenerateBeatAMLGuids(Dictionary<string, BeatAMLSample> extantGuids)
        {
            bool anyChanged = false;
            foreach (string filename in Directory.GetFiles(BeatAMLBaseDirectory + @"seqcap\bam\" , "*_Skin.realign.bam"))
            {
                var sampleName = BeatAMLFilenameToSampleName(filename);

                if (!extantGuids.ContainsKey(sampleName))
                {
                    extantGuids.Add(sampleName, new BeatAMLSample(sampleName));
                    anyChanged = true;
                }

                if (extantGuids[sampleName].normalDNAGuid == null)
                {
                    extantGuids[sampleName].normalDNAGuid = Guid.NewGuid().ToString();
                    anyChanged = true;
                }
            }

            foreach (string filename in Directory.GetFiles(BeatAMLBaseDirectory + @"seqcap\bam\" , "*_AML.realign.bam"))
            {
                var sampleName = BeatAMLFilenameToSampleName(filename);

                if (!extantGuids.ContainsKey(sampleName))
                {
                    extantGuids.Add(sampleName, new BeatAMLSample(sampleName));
                    anyChanged = true;
                }

                if (extantGuids[sampleName].tumorDNAGuid == null)
                {
                    extantGuids[sampleName].tumorDNAGuid = Guid.NewGuid().ToString();
                    anyChanged = true;
                }
            }

            foreach (string filename in Directory.GetFiles(BeatAMLBaseDirectory + @"rnaseq\bam-sorted\", "*_sorted.bam"))
            {
                var sampleName = BeatAMLFilenameToSampleName(filename);

                if (!extantGuids.ContainsKey(sampleName))
                {
                    extantGuids.Add(sampleName, new BeatAMLSample(sampleName));
                    anyChanged = true;
                }

                if (extantGuids[sampleName].tumorRNAGuid == null)
                {
                    extantGuids[sampleName].tumorRNAGuid = Guid.NewGuid().ToString();
                    anyChanged = true;
                }
            }

            if (anyChanged)
            {
                var outputFile = new StreamWriter(BeatAMLGuidsFile);

                foreach (var guidSet in extantGuids)
                {
                    outputFile.Write(guidSet.Key + "\t");

                    if (guidSet.Value.normalDNAGuid != null)
                    {
                        outputFile.Write(guidSet.Value.normalDNAGuid);
                    }
                    outputFile.Write("\t");

                    if (guidSet.Value.tumorDNAGuid != null)
                    {
                        outputFile.Write(guidSet.Value.tumorDNAGuid);
                    }
                    outputFile.Write("\t");

                    if (guidSet.Value.tumorRNAGuid != null)
                    {
                        outputFile.Write(guidSet.Value.tumorRNAGuid);
                    }
                    outputFile.WriteLine();
                }

                outputFile.Close();
            }
        }


        static public void LoadStateFromExperimentsFile(out List<AnalysisID> excludedAnalyses, out Dictionary<ParticipantID, TCGARecord> tcgaRecords, out Dictionary<SampleID, ParticipantID> sampleToParticipantIDMap,
            out Dictionary<ParticipantID, Participant> participants, out List<Experiment> experiments, out Dictionary<string, ExpressionTools.Sample> allSamples)
        {
            excludedAnalyses = ExpressionTools.LoadExcludedAnalyses(@"\\gcr\scratch\b99\bolosky\excluded_analyses.txt");

            tcgaRecords = ExpressionTools.LoadTCGARecords(null /* stored BAMs*/, excludedAnalyses, @"\\gcr\scratch\b99\bolosky\tcga-all.xml");
            ExpressionTools.LoadTCGARecordsForLocalRealigns(tcgaRecords, null, @"\\gcr\scratch\b99\bolosky\realigns.txt");
            ExpressionTools.LoadTCGAAdditionalMetadata(tcgaRecords, @"\\gcr\scratch\b99\bolosky\tcgaAdditionalMetadata.txt");

            sampleToParticipantIDMap = ExpressionTools.CreateSampleToParticipantIDMap(tcgaRecords);

            participants = ExpressionTools.BuildParticipantData(tcgaRecords, out allSamples, @"\\gcr\scratch\b99\bolosky\clinical");

            ExpressionTools.AddAllMAFFilesToParticipants(participants, sampleToParticipantIDMap, @"\\gcr\scratch\b99\bolosky\mafs");

            experiments = ExpressionTools.LoadExperimentsFromFile(@"\\gcr\scratch\b99\bolosky\experiments.txt", participants, tcgaRecords);
        }

        static public bool isChromosomeSex(string chromosome)
        {
            string lowerChromosome = chromosome.ToLower();
            return lowerChromosome == "x" || lowerChromosome == "chrx" || lowerChromosome == "y" || lowerChromosome == "chry";
        }

        static public bool isChromosomeMitochondrial(string chromosome)
        {
            string lowerChromosome = chromosome.ToLower();
            return lowerChromosome == "m" || lowerChromosome == "mt" || lowerChromosome == "chrm" || lowerChromosome == "chrmt";
        }

        static public bool isChromosomeAutosomal(string chromosome)
        {
            return !isChromosomeSex(chromosome) && !isChromosomeMitochondrial(chromosome);
        }

        public class Gene
        {

            public Gene(string hugo_symbol_, string chromosome_, int offset)
            {
                hugo_symbol = hugo_symbol_;
                chromosome = chromosome_;
                minOffset = offset;
                maxOffset = offset;

                sex = isChromosomeSex(chromosome);
                mitochondrial = isChromosomeMitochondrial(chromosome);
                autosomal = isChromosomeAutosomal(chromosome);
            }

            public string hugo_symbol;  // The gene name
            public string chromosome;
            public int minOffset;
            public int maxOffset;
            public readonly bool autosomal;
            public readonly bool sex;
            public readonly bool mitochondrial;

            public bool inconsistent = false;
        }

        public class MutationMap
        {
            public static int largestAllowedGene = 2100000; // At little bigger than DMD, the largest human gene
            public MutationMap() { }

            public void AddMutation(string hugo_symbol, string chromosome, int offset)
            {
                if (!genes.ContainsKey(hugo_symbol))
                {
                    genes.Add(hugo_symbol, new Gene(hugo_symbol, chromosome, offset));
                }

                var gene = genes[hugo_symbol];

                if (!genes[hugo_symbol].inconsistent)
                {
                    if (gene.chromosome != chromosome)
                    {
                        Console.WriteLine("Gene " + hugo_symbol + " occurs on (at least) two different chromosomes: " + chromosome + " and " + gene.chromosome + ".  Ignoring gene.");
                        gene.inconsistent = true;
                    }
                    else
                    {
                        gene.minOffset = Math.Min(gene.minOffset, offset);
                        gene.maxOffset = Math.Max(gene.maxOffset, offset);
                        if (gene.maxOffset - gene.minOffset > largestAllowedGene)
                        {
                            Console.WriteLine("Gene " + hugo_symbol + " has too big a range of mutation offsets: " + chromosome + ":" + gene.minOffset + "-" + gene.maxOffset);
                            gene.inconsistent = true;
                        }
                    }
                }
            }

            public void DoneAddingMutations()
            {
                foreach (var geneEntry in genes)
                {
                    var gene = geneEntry.Value;

                    if (!gene.inconsistent)
                    {
                        if (!genesByChromosome.ContainsKey(gene.chromosome))
                        {
                            genesByChromosome.Add(gene.chromosome, new List<Gene>());
                        }
                        genesByChromosome[gene.chromosome].Add(gene);
                        genesByName.Add(gene.hugo_symbol, gene);
                    }
                }
            }

            public int Count()
            {
                return genes.Count();
            }

            Dictionary<string, Gene> genes = new Dictionary<string, Gene>();
            public Dictionary<string, List<Gene>> genesByChromosome = new Dictionary<string, List<Gene>>();
            public Dictionary<string, Gene> genesByName = new Dictionary<string, Gene>();
        }

        static public Dictionary<string, MutationMap> GenerateMutationMapFromExperiments(List<Experiment> experiments, Dictionary<string, ExpressionTools.Experiment> experimentsByParticipant)
        {

            var mutations = new Dictionary<string, ExpressionTools.MutationMap>();  // Maps reference type (hg18 or hg19) to MutationMap
            mutations.Add("hg18", new ExpressionTools.MutationMap());
            mutations.Add("hg19", new ExpressionTools.MutationMap());


                
            var badHugoSymbols = new List<string>();    // These are corruped by Excel.  They're in the MAFs that I downloaded, and I'm removing them by hand.
            badHugoSymbols.Add("1-Mar");
            badHugoSymbols.Add("1-Dec");
            badHugoSymbols.Add("1-Sep");
            badHugoSymbols.Add("10-Mar");
            badHugoSymbols.Add("11-Mar");
            badHugoSymbols.Add("12-Sep");
            badHugoSymbols.Add("14-Sep");
            badHugoSymbols.Add("1SEPT4");
            badHugoSymbols.Add("2-Mar");
            badHugoSymbols.Add("2-Sep");
            badHugoSymbols.Add("3-Mar");
            badHugoSymbols.Add("3-Sep");
            badHugoSymbols.Add("4-Mar");
            badHugoSymbols.Add("4-Sep");
            badHugoSymbols.Add("5-Sep");
            badHugoSymbols.Add("6-Mar");
            badHugoSymbols.Add("6-Sep");
            badHugoSymbols.Add("7-Mar");
            badHugoSymbols.Add("7-Sep");
            badHugoSymbols.Add("8-Mar");
            badHugoSymbols.Add("9-Sep");

            int nMutations = 0;
            foreach (var experiment in experiments)
            {
                experimentsByParticipant.Add(experiment.participant.participantId, experiment);

                foreach (ExpressionTools.MAFRecord maf in experiment.maf)
                {
                    nMutations++;

                    string chromosome = chromosomeNameToNonChrForm(maf.Chrom);

                    if (badHugoSymbols.Contains(maf.Hugo_symbol))
                    {
                        Console.WriteLine("Bad hugo symbol " + maf.Hugo_symbol);
                    }
                    mutations[maf.ReferenceClass()].AddMutation(maf.Hugo_symbol, chromosome, maf.Start_position);
                }
            }

            foreach (var referenceClass in mutations)
            {
                referenceClass.Value.DoneAddingMutations();
            }

            return mutations;
        }

        static public string switchChrState(string inputContigName)
        {
            if (inputContigName.Count() < 3 || inputContigName.Substring(0, 3).ToLower() != "chr")
            {
                return "chr" + inputContigName;
            }

            return inputContigName.Substring(3);
        }


        public class AllcountReader
        {
            public AllcountReader(string filename_) 
            {
                filename = filename_;
            }


            public bool openFile(out long mappedHQNuclearReads, out int out_numContigs)
            {
                if (null != allcountReader)
                {
                    Console.WriteLine("ExpressionTools.AllcountReader.openFile(): file is already open.  You can only do this once per allcount file.");
                }

                allcountReader = null;
                out_numContigs = numContigs = 0;
                contigs = null;
                mappedHQNuclearReads = 0;

                try
                {
                    compressedStreamReader = CreateCompressedStreamReaderWithRetry(filename);
                    allcountReader = new StreamReader(new GZipStream(compressedStreamReader.BaseStream, CompressionMode.Decompress));
                }
                catch (IOException)
                {
                    Console.WriteLine("IOException opening allcount file " + filename);
                    return false;
                }

                //
                // The format is two header lines like:
                //      CountReadsCovering v1.1 \\msr-genomics-0\d$\sequence\indices\grch37-24 -a \\msr-genomics-0\e$\tcga\hnsc\tumor\rna\0415c9cc-48e6-46e7-9076-06b6b56bb4be\PCAWG.e7697e88-96df-4c96-ad8d-6c54df95d29b.STAR.v1.bam -
                //      83303682 mapped high quality nuclear reads, 9891899 low quality reads, 0 mitochondrial reads, 100740272 total reads
                // followed by a blank line followed by
                //      NumContigs: 93
                //      ContigName\tLength
                //      <ContigName\tLength> x numContigs
                // Then for each contig:
                //      >ContigName
                // And within each contig one of three line types:
                //      offset\tcount
                // Which is an offset in the contig (in HEX) followed by the count of high quality reads that map there (also in hex) or:
                //      xnumber
                // which is the number of consecutive loci with the same number of HQ reads mapped as the previous locus (x is literal, count is in hex) or:
                //      count
                // which is the new count of reads mapped for the next locus (in hex)
                //
                // And the final line in the file is
                //      **done**
                // 
                //

                var line = allcountReader.ReadLine();

                const string headerBeginning = "CountReadsCovering v";

                if (null == line || line.Count() < headerBeginning.Count() + 1 || line.Substring(0, headerBeginning.Count()) != headerBeginning)
                {
                    Console.WriteLine("Empty or corrupt allcount file " + filename);
                    return false;
                }

                if (line[headerBeginning.Count()] != '1')
                {
                    Console.WriteLine("Unsupported major version of allcount file " + filename + ".  Header line: " + line);
                    return false;
                }

                line = allcountReader.ReadLine();
                if (null == line)
                {
                    Console.WriteLine("Corrupt or tuncated allcount file " + filename);
                    return false;
                }

                var fields = line.Split(' ');
                if (fields.Count() != 16)
                {
                    Console.WriteLine("Corrupt or tuncated allcount file " + filename + ".  Second line has " + fields.Count() + " fields: " + line);
                    return false;
                }

                try
                {
                    mappedHQNuclearReads = Convert.ToInt64(fields[0]);
                }
                catch (FormatException)
                {
                    Console.WriteLine("Format exception parsing mapped HQ read count for file " + filename + " from line: " + line);
                    return false;
                }

                line = allcountReader.ReadLine();   // The blank line

                line = allcountReader.ReadLine();
                if (null == line)
                {
                    Console.WriteLine("Allcount file truncated before contig count: " + filename);
                    return false;
                }

                const string numContigsLineBeginning = "NumContigs: ";
                if (line.Count() < numContigsLineBeginning.Count() + 1 || line.Substring(0, numContigsLineBeginning.Count()) != numContigsLineBeginning)
                {
                    Console.WriteLine("Malformed NumContigs line in " + filename + ": " + line);
                    return false;
                }

                try
                {
                    out_numContigs = numContigs = Convert.ToInt32(line.Substring(numContigsLineBeginning.Count()));
                }
                catch (FormatException)
                {
                    Console.WriteLine("Couldn't parse NumContigs line in file " + filename + ": " + line);
                    return false;
                }

                if (numContigs < 1)
                {
                    Console.WriteLine("Invalid numContigs in " + filename + ": " + line);
                    return false;
                }

                line = allcountReader.ReadLine();   // The header line for the contigs.

                contigs = new Contig[numContigs];

                int whichContig;

                for (whichContig = 0; whichContig < numContigs; whichContig++)
                {
                    contigs[whichContig] = new Contig();

                    line = allcountReader.ReadLine();
                    if (null == line)
                    {
                        Console.WriteLine("File truncated in contig list " + filename);
                        return false;
                    }

                    fields = line.Split('\t');
                    if (fields.Count() != 2)
                    {
                        Console.WriteLine("Incorrect contig line format in file " + filename + ": " + line);
                        return false;
                    }

                    contigs[whichContig].name = fields[0].ToLower();
                    try
                    {
                        contigs[whichContig].length = Convert.ToInt64(fields[1]);
                    }
                    catch (FormatException)
                    {
                        Console.WriteLine("Incorrect contig line format in file " + filename + ": " + line);
                        return false;
                    }
                } // for all expected contigs

                return true;
            }

            public delegate void ProcessBase(string strippedContigName, int location, int currentMappedReadCount);

            public bool ReadAllcountFile(ProcessBase processBase)
            {
                int currentOffset = -1;
                int currentMappedReadCount = -1;
                int whichContig = -1;

                bool sawDone = false;
                string strippedContigName = "";

                string line;

                while (null != (line = allcountReader.ReadLine()))
                {
                    if (sawDone)
                    {
                        Console.WriteLine("File " + filename + " continues after **done** line: " + line);
                        return false;
                    }

                    if ("**done**" == line)
                    {
                        sawDone = true;
                        continue;
                    }

                    if (line.Count() == 0)
                    {
                        Console.WriteLine("Unexpected blank line in " + filename);
                        return false;
                    }

                    if (line[0] == '>')
                    {
                        whichContig++;
                        if (whichContig >= numContigs)
                        {
                            Console.WriteLine("Saw too many contigs in " + filename + ": " + line);
                            return false;
                        }

                        if (line.Substring(1).ToLower() != contigs[whichContig].name)
                        {
                            Console.WriteLine("Unexpected contig in " + filename + ".  Expected " + contigs[whichContig].name + ", got ", line.Substring(1));
                            return false;
                        }

                        if (line.Count() > 4 && line.Substring(1, 3).ToLower() == "chr")
                        {
                            strippedContigName = line.Substring(4).ToLower();
                        }
                        else
                        {
                            strippedContigName = line.Substring(1).ToLower();
                        }


                        currentOffset = -1;
                        currentMappedReadCount = -1;
                        continue;
                    }

                    if (-1 == whichContig)
                    {
                        Console.WriteLine("Expected contig line after list of contigs, got " + line);
                        return false;
                    }

                    var fields = line.Split('\t');
                    if (fields.Count() == 1)
                    {
                        //
                        // Either xRepeatCount or newCount
                        //
                        if (line[0] == 'x')
                        {
                            if (currentMappedReadCount < 0 || currentOffset < 0)
                            {
                                Console.WriteLine("Got unexpected x line " + line);
                                return false;
                            }

                            int repeatCount;
                            try
                            {
                                repeatCount = Convert.ToInt32(line.Substring(1), 16); // 16 means the string is in hex
                                if (repeatCount <= 1)
                                {
                                    Console.WriteLine("Bogus repeat count " + line);
                                    return false;
                                }
                            }
                            catch (SystemException ex)
                            {
                                if (ex is FormatException || ex is ArgumentOutOfRangeException) {
                                    Console.WriteLine("Format exception processing x line " + line);
                                    return false;
                                } else {
                                    throw ex;
                                }
                            }
                            for (; repeatCount > 1; repeatCount--)  // > 1 because this count includes the locus that specified the mapped read count (the previous non-x line) and we already emitted that one.
                            {
                                currentOffset++;
                                processBase(strippedContigName, currentOffset, currentMappedReadCount);
                            }
                        }
                        else
                        {
                            //
                            // A new mapped read count line
                            //
                            if (currentOffset <= 0)
                            {
                                Console.WriteLine("Got unexpected mapped read count line " + line);
                                return false;
                            }
                            try
                            {
                                currentMappedReadCount = Convert.ToInt32(line, 16); // 16 means the string is in hex
                                if (currentMappedReadCount <= 0)
                                {
                                    Console.WriteLine("Bogus current count " + line);
                                    return false;
                                }
                            }
                            catch (SystemException ex)
                            {
                                if (ex is FormatException || ex is ArgumentOutOfRangeException)
                                {
                                    Console.WriteLine("Format exception processing count line " + line);
                                    return false;
                                }
                                else
                                {
                                    throw ex;
                                }
                            }

                            currentOffset++;
                            processBase(strippedContigName, currentOffset, currentMappedReadCount);
                        }

                        continue;
                    }

                    if (fields.Count() != 2)
                    {
                        Console.WriteLine("Saw too many fields in line " + line);
                        return false;
                    }

                    //
                    // An offset + count line
                    //
                    try
                    {  
                        currentOffset = Convert.ToInt32(fields[0], 16); // 16 means the string is in hex
                        currentMappedReadCount = Convert.ToInt32(fields[1], 16); // 16 means the string is in hex

                        if (currentOffset <= 0 || currentMappedReadCount <= 0)
                        {
                            Console.WriteLine("Bogus offset + count line " + line);
                            return false;
                        }
                    }
                    catch (SystemException ex)
                    {
                        if (ex is FormatException || ex is ArgumentOutOfRangeException)
                        {
                            Console.WriteLine("Unable to parse offset + count line " + line);
                            return false;
                        }
                        else
                        {
                            throw ex;
                        }
                    }

                    processBase(strippedContigName, currentOffset, currentMappedReadCount);
                }

                if (!sawDone)
                {
                    Console.WriteLine("Truncated allcount file " + filename);
                    return false;
                }

                return true;
            }

            public void Close()
            {
                if (null != allcountReader)
                {
                    allcountReader.Close();
                    compressedStreamReader.Close();
                }
            }

            StreamReader compressedStreamReader = null;
            StreamReader allcountReader = null;
            int numContigs = 0;
            Contig[] contigs = null;

            public readonly string filename;


            class Contig
            {
                public string name = "";
                public long length = -1;
            }

        }

        //
        // Excel has the bad habit of taking strings that looks like dates and converting them into actual dates.  This is a problem for some gene names ("MARCH1"), which is
        // bad enough that it's actually made into the literature (not to mention the mafs that I imported).  This takes a string and converts it into a format that Excel
        // will not munge.
        //
        public static string ConvertToExcelString(string input)
        {
            if (input.Count() < 3 || input.Substring(0, 2) != "=\"" || input[input.Count() - 1] != '\"')
            {
                return "=\"" + input + "\"";
            }

            //
            // It's already in excel format, just keep it.
            //
            return input;
        }

        public static string ConvertToNonExcelString(string input)
        {
            if (input.Count() < 3 || input.Substring(0, 2) != "=\"" || input[input.Count() - 1] != '\"')
            {
                //
                // It's not an excel string, just return it.
                //
                return input;
            }

            return input.Substring(2, input.Count() - 3);
        }

        public static StreamWriter CreateStreamWriterWithRetry(string filename)
        {
            while (true)
            {
                try
                {
                    var writer = new StreamWriter(filename);
                    return writer;
                }
                catch (IOException)
                {
                    Console.WriteLine("IOException opening " + filename + " for write.  Sleeping and retrying.");
                    Thread.Sleep(10 * 1000);
                }
            }
        }

        public static StreamReader CreateStreamReaderWithRetry(string filename)
        {
            while (true)
            {
                try
                {
                    var reader = new StreamReader(filename);
                    return reader;
                }
                catch (IOException)
                {
                    Console.WriteLine("IOException opening " + filename + " for read.  Sleeping and retrying.");
                    Thread.Sleep(10 * 1000);
                }
            }
        }

        public static StreamReader CreateCompressedStreamReaderWithRetry(string filename)
        {
            return new StreamReader(new GZipStream(CreateStreamReaderWithRetry(filename).BaseStream, CompressionMode.Decompress));
        }

        public static string[] ReadAllLinesWithRetry(string filename)
        {
            while (true)
            {
                try
                {
                    var lines = File.ReadAllLines(filename);
                    return lines;
                }
                catch (IOException)
                {
                    Console.WriteLine("IOException reading " + filename + ".  Sleeping and retrying.");
                    Thread.Sleep(10 * 1000);
                }
            }
        }


        public class HistogramResultLine
        {
            public string minValue;
            public int count = 0;
            public double total = 0;    // The sum of all of the values in this line
            public double pdfValue = 0;
            public double cdfValue = 0;
        }

        public class Histogram
        {
            public Histogram() { }

            public void addValue(double value)
            {
                values.Add(value);
            }

            public double min() // Throws ArgumentNullExcpetion if no data jas been added
            {
                return values.Min();
            }

            public double max()
            {
                return values.Max();
            }


            public int count()
            {
                return values.Count();
            }
            public HistogramResultLine [] ComputeHistogram(double min, double max, double increment, string format = "G")
            {
                int nBuckets = (int)((max - min) / increment) + 1;  // +1 is for "more"

                var result = new HistogramResultLine[nBuckets];

                double x = min;
                for (int i = 0; i < nBuckets - 1; i++)
                {
                    result[i] = new HistogramResultLine();
                    result[i].minValue = x.ToString(format);
                    x += increment;
                }
                result[nBuckets - 1] = new HistogramResultLine();
                result[nBuckets - 1].minValue = "More";

                foreach (var value in values)
                {
                    int whichBucket;
                    if (value > max)
                    {
                        whichBucket = nBuckets - 1;
                    } else {
                        whichBucket = (int)((value - min) / increment);
                    }

                    result[whichBucket].count++;
                    result[whichBucket].total += value;
                }

                int overallCount = values.Count();
                int runningCount = 0;
                for (int whichBucket = 0; whichBucket < nBuckets; whichBucket++)
                {
                    runningCount += result[whichBucket].count;

                    result[whichBucket].pdfValue = ((double)result[whichBucket].count) / overallCount;
                    result[whichBucket].cdfValue = ((double)runningCount) / overallCount;
                }

                return result;  
            }

            List<double> values = new List<double>();
        } // Histogram

        public class GeneScatterGraphLine
        {
            public string RNAFile;
            public string DNAFile;
            public string Hugo_Symbol;
            public string Entrez_Gene_Id;
            public string Center;
            public string NCBI_Build;
            public string Chromosome;
            public string Start_Position;
            public string End_Position;
            public string Strand;
            public string Variant_Classification;
            public string Variant_Type;
            public string Reference_Allele;
            public string Tumor_Seq_Allele_1;
            public string Tumor_Seq_Allele_2;
            public string dbSNP_RS;
            public string dbSNP_Val_Status;
            public string Tumor_Sample_Barcode;
            public string Matched_Norm_Sample_Barcode;
            public string Match_Norm_Seq_Allele1;
            public string Match_Norm_Seq_Allele2;
            public string Tumor_Validation_Allele1;
            public string Tumor_Validation_Allele2;
            public string Match_Norm_Validation_Allele1;
            public string Match_Norm_Validation_Allele2;
            public string Verification_Status;
            public string Validation_Status;
            public string Mutation_Status;
            public string Sequencing_Phase;
            public string Sequence_Source;
            public string Validation_Method;
            public string Score;
            public string BAM_File;
            public string Sequencer;
            public string Tumor_Sample_UUID;
            public string Matched_Norm_Sample_UUID;
            public string File_Name;
            public string Archive_Name;
            public int Line_Number;
            public int n_DNA_Matching_Reference;
            public int n_DNA_Matching_Tumor;
            public int n_DNA_Matching_Neither;
            public int n_DNA_Matching_Both;
            public int n_RNA_Matching_Reference;
            public int n_RNA_Matching_Tumor;
            public int n_RNA_Matching_Neither;
            public int n_RNA_Matching_Both;
            public double tumorDNAFraction;
            public double tumorRNAFraction;
            public double tumorDNAMultiple;
            public double tumorRNAMultiple;
            public double tumorDNARatio;
            public double tumorRNARatio;
            public string FractionOfMeanExpression; // This doesn't seem to be filled in, hence string
            public string zOfmeanExpression;
            public double ratioOfRatios;
            public bool IsSingle;
            public string CancerType;
            public string gender;
            public double zTumor;
            public double zNormal;
            public double z2Tumor;
            public double z2Normal;
            public double percentMeanTumor; // Expressed as a fraction (i.e., 100% -> 1.0)
            public double percentMeanNormal;// Expressed as a fraction (i.e., 100% -> 1.0)

            static GeneScatterGraphLine fromLine(string inputLine)
            {
                var newLine = new GeneScatterGraphLine();

                var fields = inputLine.Split('\t');

                if (fields.Count() < 65) {
                    Console.WriteLine("GeneScatterGraphEntry.fromLine: too few fields " + fields.Count() + " in input line " + inputLine);
                    return null;
                }

                try
                {
                    newLine.RNAFile = fields[0];
                    newLine.DNAFile = fields[1];
                    newLine.Hugo_Symbol = ConvertToNonExcelString(fields[2]);
                    newLine.Entrez_Gene_Id = fields[3];
                    newLine.Center = fields[4];
                    newLine.NCBI_Build = fields[5];
                    newLine.Chromosome = fields[6];
                    newLine.Start_Position = fields[7];
                    newLine.End_Position = fields[8];
                    newLine.Strand = fields[9];
                    newLine.Variant_Classification = fields[10];
                    newLine.Variant_Type = fields[11];
                    newLine.Reference_Allele = fields[12];
                    newLine.Tumor_Seq_Allele_1 = fields[13];
                    newLine.Tumor_Seq_Allele_2 = fields[14];
                    newLine.dbSNP_RS = fields[15];
                    newLine.dbSNP_Val_Status = fields[16];
                    newLine.Tumor_Sample_Barcode = fields[17];
                    newLine.Matched_Norm_Sample_Barcode = fields[18];
                    newLine.Match_Norm_Seq_Allele1 = fields[19];
                    newLine.Match_Norm_Seq_Allele2 = fields[20];
                    newLine.Tumor_Validation_Allele1 = fields[21];
                    newLine.Tumor_Validation_Allele2 = fields[22];
                    newLine.Match_Norm_Validation_Allele1 = fields[23];
                    newLine.Match_Norm_Validation_Allele2 = fields[24];
                    newLine.Verification_Status = fields[25];
                    newLine.Validation_Status = fields[26];
                    newLine.Mutation_Status = fields[27];
                    newLine.Sequencing_Phase = fields[28];
                    newLine.Sequence_Source = fields[29];
                    newLine.Validation_Method = fields[30];
                    newLine.Score = fields[31];
                    newLine.BAM_File = fields[32];
                    newLine.Sequencer = fields[33];
                    newLine.Tumor_Sample_UUID = fields[34];
                    newLine.Matched_Norm_Sample_UUID = fields[35];
                    newLine.File_Name = fields[36];
                    newLine.Archive_Name = fields[37];
                    newLine.Line_Number = Convert.ToInt32(fields[38]);
                    newLine.n_DNA_Matching_Reference = Convert.ToInt32(fields[39]);
                    newLine.n_DNA_Matching_Tumor = Convert.ToInt32(fields[40]);
                    newLine.n_DNA_Matching_Neither = Convert.ToInt32(fields[41]);
                    newLine.n_DNA_Matching_Both = Convert.ToInt32(fields[42]);
                    newLine.n_RNA_Matching_Reference = Convert.ToInt32(fields[43]);
                    newLine.n_RNA_Matching_Tumor = Convert.ToInt32(fields[44]);
                    newLine.n_RNA_Matching_Neither = Convert.ToInt32(fields[45]);
                    newLine.n_RNA_Matching_Both = Convert.ToInt32(fields[46]);
                    newLine.tumorDNAFraction = Convert.ToDouble(fields[47]);
                    newLine.tumorRNAFraction = Convert.ToDouble(fields[48]);
                    newLine.tumorDNAMultiple = Convert.ToDouble(fields[49]);
                    newLine.tumorRNAMultiple = Convert.ToDouble(fields[50]);
                    newLine.tumorDNARatio = Convert.ToDouble(fields[51]);
                    newLine.tumorRNARatio = Convert.ToDouble(fields[52]);
                    newLine.FractionOfMeanExpression = fields[53];
                    newLine.zOfmeanExpression = fields[54];
                    newLine.ratioOfRatios = Convert.ToDouble(fields[55]);
                    newLine.IsSingle = Convert.ToBoolean(fields[56]);
                    newLine.CancerType = fields[57];
                    newLine.gender = fields[58];
                    newLine.zTumor = Convert.ToDouble(fields[59]);
                    newLine.zNormal = Convert.ToDouble(fields[60]);
                    newLine.z2Tumor = Convert.ToDouble(fields[61]);
                    newLine.z2Normal = Convert.ToDouble(fields[62]);
                    newLine.percentMeanTumor = Convert.ToDouble(fields[63]) / 100.0; // Expressed as a fraction (i.e., 100% -> 1.0 = fields[];
                    newLine.percentMeanNormal = Convert.ToDouble(fields[64]) / 100.0;// Expressed as a fraction (i.e., 100% -> 1.0 = fields[];
                }
                catch (FormatException)
                {
                    Console.WriteLine("GeneScatterGraphEntry.fromLine: error parsing line " + inputLine);
                    return null;
                }

                return newLine;
            } // fromLine
        }

        public static void LoadAllGeneScatterGraphEntries(out List<GeneScatterGraphLine> geneScatterGraphEntries, string directoryName = @"f:\temp\gene_scatter_graphs")
        {
            geneScatterGraphEntries = new List<GeneScatterGraphLine>();

            foreach (var filename in Directory.EnumerateFiles(directoryName, "*.txt")) {
                if (filename.Count() == 0 || filename[0] == '_')
                {
                    continue;   // Summary file like _MannWhitney rather than a gene file
                }

                var lines = ReadAllLinesWithRetry(filename);
                for (int i = 1 ; i < lines.Count(); i++) {

                }
            }

        }

        public class Exon
        {
            public Exon(string startString, string endString)
            {
                start = Convert.ToInt32(startString);
                end = Convert.ToInt32(endString);
            }

            public int start;
            public int end;
        }

        public class Isoform    // These come from the knownGene file
        {
            public string ucscId;
            public string chromosome;   // The non-chr version
            public string strand;
            public int txStart;
            public int txEnd;
            public int cdsStart;
            public int cdsEnd;
            public string proteinId;
            public string alignId;
            public Exon [] exons;

            static Isoform fromFileLine(string fileLine)
            {
                var fields = fileLine.Split('\t');
                if (fields.Count() != 12)
                {
                    Console.WriteLine("Isoform.fromFileLine: line had wrong number of fields, " + fields.Count() + " != 12");
                    return null;
                }

                var isoform = new Isoform();

                isoform.ucscId = ConvertToNonExcelString(fields[0]);
                isoform.chromosome = chromosomeNameToNonChrForm(ConvertToNonExcelString(fields[1]));
                isoform.strand = ConvertToNonExcelString(fields[2]);

                int exonCount;
                try
                {
                    isoform.txStart = Convert.ToInt32(ConvertToNonExcelString(fields[3]));
                    isoform.txEnd = Convert.ToInt32(ConvertToNonExcelString(fields[4]));

                    if (isoform.txEnd <= isoform.txStart)
                    {
                        Console.WriteLine("Isoform.fromFileLine: warning: isoform " + isoform.ucscId + " has empty or negative transcription region " + fileLine);
                    }

                    isoform.cdsStart = Convert.ToInt32(ConvertToNonExcelString(fields[5]));
                    isoform.cdsEnd = Convert.ToInt32(ConvertToNonExcelString(fields[6]));
                    exonCount = Convert.ToInt32(ConvertToNonExcelString(fields[7]));


                    var exonStartStrings = ConvertToNonExcelString(fields[8]).Split(',');
                    var exonEndStrings = ConvertToNonExcelString(fields[9]).Split(',');

                    //
                    // They have trailing commas, so they should have one more field than there are exons, and also should have an empty string
                    // as their last element.
                    //
                    if (exonStartStrings.Count() != exonCount + 1 || exonEndStrings.Count() != exonCount + 1 || exonStartStrings[exonCount] != "" || exonEndStrings[exonCount] != "")
                    {
                        Console.WriteLine("Isoform.fromFileLine: Bad exon start/end: " + fileLine);
                    }

                    isoform.exons = new Exon[exonCount];
                    for (int i = 0; i < exonCount; i++)
                    {
                        isoform.exons[i] = new Exon(exonStartStrings[i], exonEndStrings[i]);
                    }


                }
                catch (FormatException)
                {
                    Console.WriteLine("Isoform.fromFileLine: Format exception parsing a numeric field in line: " + fileLine);
                    return null;
                }

                isoform.proteinId = fields[10];
                isoform.alignId = fields[11];

                return isoform;
            }

            public static Dictionary<string, Isoform> readKnownGeneFile(string fileName)
            {
                var file = CreateStreamReaderWithRetry(fileName);

                var result = new Dictionary<string, Isoform>();

                file.ReadLine();    // Skip the header
                string line;

                while (null != (line = file.ReadLine()))
                {
                    var isoform = Isoform.fromFileLine(line);
                    result.Add(isoform.ucscId, isoform);
                }

                file.Close();

                return result;
            }
        }

            

        public static Dictionary<string, List<Isoform>> loadIsoformGroupMapFromFile(string filename, Dictionary<string, Isoform> isoforms)    // This loads the known isoforms file.
        {
            var file = CreateStreamReaderWithRetry(filename);

            var result = new Dictionary<string, List<Isoform>>();

            file.ReadLine();    // Header

            string line;
            int currentClusterId = 0;
            List<Isoform> currentCluster = new List<Isoform>();

            while (null != (line = file.ReadLine()))
            {
                var fields = line.Split('\t');
                if (fields.Count() != 2)
                {
                    Console.WriteLine("loadIsoformGroupMapFromFile: wrong field count in input file line " + line);
                    return null;
                }

                var clusterId = Convert.ToInt32(ConvertToNonExcelString(fields[0]));

                if (clusterId != currentClusterId)
                {
                    foreach (var isoform in currentCluster)
                    {
                        result.Add(isoform.ucscId, currentCluster); // Maps one isoform to the whole cluster (gene)
                    }

                    currentCluster = new List<Isoform>();
                }

                currentClusterId = clusterId;

                var isoformId = ConvertToNonExcelString(fields[1]);

                if (!isoforms.ContainsKey(isoformId))
                {
                    Console.WriteLine("Isoform group file contains isoform ID " + isoformId + " that's not in the isoform dictionary.  Ignoring.");
                    continue;
                }

                currentCluster.Add(isoforms[isoformId]);
            }

            file.Close();

            return result;
        }

        static string UcscNameToStrippedUcscName(string ucscIdWithVersion)
        {
            if (-1 == ucscIdWithVersion.LastIndexOf('.'))
            {
                Console.WriteLine("UcscNameToStrippedUcscName: ucsc ID doesn't contain a dot: " + ucscIdWithVersion);
                return ucscIdWithVersion;
            }

            return ucscIdWithVersion.Substring(0, ucscIdWithVersion.LastIndexOf('.'));
        }

        public static Dictionary<string, List<Isoform>> LoadUcscNameToIsoformsMap(string knownGenesFilename, string knownIsoformsFilename, bool stripVersionFromUcscName)  // Produces a map from non-version ucsc id to isoform groups
        {
            var isoforms = Isoform.readKnownGeneFile(knownGenesFilename);
            var groupMap = loadIsoformGroupMapFromFile(knownIsoformsFilename, isoforms);

            var result = new Dictionary<string, List<Isoform>>();

            //
            // Chop off the version numbers in the keys, and build a map from the chopped versions to all of the isoforms.
            //
            foreach (var entry in groupMap)
            {
                var ucscIdWithVersion = entry.Key;
                var ucscIdToUse = stripVersionFromUcscName ? UcscNameToStrippedUcscName(ucscIdWithVersion) : ucscIdWithVersion;

                if (!result.ContainsKey(ucscIdToUse))
                {
                    result.Add(ucscIdToUse, new List<Isoform>());
                }

                result[ucscIdToUse].AddRange(entry.Value);
            }


            return result;
        }

        public static Dictionary<string, GeneInformation> LoadGeneInformationFromGeneWithProtein(string geneWithProteinProductFilename, string knownGenesFilename, string knownIsoformsFilename)
        {
            var isoformsByNonVersionName = LoadUcscNameToIsoformsMap(knownGenesFilename, knownIsoformsFilename, true);

            var genesByHugoId = GeneInformation.LoadFromFile(geneWithProteinProductFilename);

            foreach (var geneEntry in genesByHugoId)
            {
                var gene = geneEntry.Value;

                if (gene.ucsc_id == "")
                {
                    continue;   // No ucsc_id means no idea about the gene.
                }

                var strippedUcscId = UcscNameToStrippedUcscName(gene.ucsc_id);

                if (!isoformsByNonVersionName.ContainsKey(strippedUcscId) || isoformsByNonVersionName[strippedUcscId].Count() == 0)
                {
                    Console.WriteLine("Hugo symbol " + gene.symbol + " with ucsc id " + gene.ucsc_id + " doens't appear in known genes file or has no isoforms.");
                    continue;
                }

                var isoformList = isoformsByNonVersionName[strippedUcscId];
                var chromosome = isoformList[0].chromosome;
                int minLocus = isoformList[0].txStart;
                int maxLocus = isoformList[0].txEnd;

                for (int i = 1; i < isoformList.Count(); i++)
                {
                    if (chromosome != isoformList[i].chromosome)
                    {
                        Console.WriteLine("LoadGeneInformation: Warning: different isoforms of " + gene.symbol + " are on different chromosomes: " + chromosome + " and " + isoformList[i].chromosome);
                    }
                    else
                    {
                        minLocus = Math.Min(minLocus, isoformList[i].txStart);
                        maxLocus = Math.Max(maxLocus, isoformList[i].txEnd);
                    }
                }

                gene.chromosome = chromosome;
                gene.minLocus = minLocus;
                gene.maxLocus = maxLocus;
                gene.isoforms = isoformList;
            }

            return genesByHugoId;
        }

        public class GeneLocationInfo {
            public string hugoSymbol;
            public string chromosome;   // in non-chr form
            public int minLocus;
            public int maxLocus;

            public bool inconsistent = false;

            public List<Isoform> isoforms = new List<Isoform>();
        }

        public static Dictionary<string, GeneLocationInfo> LoadGeneLocationInfo(string knownGenesFilename, string kgXrefFilename)
        {
            var isoforms = Isoform.readKnownGeneFile(knownGenesFilename);

            var result = new Dictionary<string, GeneLocationInfo>();

            var kgXrefFile = CreateStreamReaderWithRetry(kgXrefFilename);

            string line;
            while (null != (line = kgXrefFile.ReadLine()))
            {
                var fields = line.Split('\t');
                if (fields.Count() != 10 && fields.Count() != 8)
                {
                    Console.WriteLine("LoadGeneLocationInfo: wrong number of fields in kgXref file " + fields.Count() + " != 10 or 8 " + ": " + line);
                    continue;
                }

                var ucscId = ConvertToNonExcelString(fields[0]);

                if (!isoforms.ContainsKey(ucscId)) {
                    Console.WriteLine("LoadGeneLocationInfo: found unknown isoform: " + line);
                    continue;
                }

                var isoform = isoforms[ucscId];

                var hugoSymbol = ConvertToNonExcelString(fields[4]);
                if (!result.ContainsKey(hugoSymbol))
                {
                    result.Add(hugoSymbol, new GeneLocationInfo());
                    result[hugoSymbol].hugoSymbol = hugoSymbol;
                    result[hugoSymbol].chromosome = isoform.chromosome;
                    result[hugoSymbol].minLocus = isoform.txStart;
                    result[hugoSymbol].maxLocus = isoform.txEnd;
                }
                else
                {
                    if (result[hugoSymbol].chromosome != isoform.chromosome && !result[hugoSymbol].inconsistent)
                    {
                        Console.WriteLine("LoadGeneLocationInfo: mismatched chromsome for " + hugoSymbol + ": " + result[hugoSymbol].chromosome + " != " + isoform.chromosome);
                        result[hugoSymbol].inconsistent = true;
                        continue;
                    }

                    result[hugoSymbol].minLocus = Math.Min(result[hugoSymbol].minLocus, isoform.txStart);
                    result[hugoSymbol].maxLocus = Math.Max(result[hugoSymbol].maxLocus, isoform.txEnd);
                }

                result[hugoSymbol].isoforms.Add(isoform);
            }

            kgXrefFile.Close();
            return result;
        }

        public class GeneLocationsByNameAndChromosome
        {
            public GeneLocationsByNameAndChromosome(Dictionary<string, GeneLocationInfo> genesByName_)
            {
                genesByName = genesByName_;

                foreach (var geneEntry in genesByName)
                {
                    var gene = geneEntry.Value;
                    if (!genesByChromosome.ContainsKey(gene.chromosome))
                    {
                        genesByChromosome.Add(gene.chromosome, new List<GeneLocationInfo>());
                    }

                    genesByChromosome[gene.chromosome].Add(gene);
                }
            }

            public Dictionary<string, GeneLocationInfo> genesByName;
            public Dictionary<string, List<GeneLocationInfo>> genesByChromosome = new Dictionary<string ,List<GeneLocationInfo>>();    // chromosome is in non-chr form.
        }

        public class ExonicMap
        {
            public ExonicMap(GeneLocationsByNameAndChromosome geneMap)
            {
                foreach (var entry in geneMap.genesByChromosome)
                {
                    chromsomeMapSizes.Add(entry.Key, 0);
                }

                foreach (var geneEntry in geneMap.genesByName)
                {
                    var gene = geneEntry.Value;

                    foreach (var isoform in gene.isoforms)
                    {
                        foreach (var exon in isoform.exons)
                        {
                            chromsomeMapSizes[gene.chromosome] = Math.Max(chromsomeMapSizes[gene.chromosome], exon.end);
                        }
                    }
                }

                //
                // We now know the largest base in any exon in each chromosome.  Allocate the maps.
                //

                foreach (var chromsomeEntry in chromsomeMapSizes)
                {
                    var chromosome = chromsomeEntry.Key;
                    var highestBaseInAnExon = chromsomeEntry.Value;
                    map.Add(chromosome, new bool[highestBaseInAnExon + 1]);  // +1 becuase the genome isn't 0 based.  Biologists, what can you do?
                    for (int i = 0; i <= highestBaseInAnExon; i++)
                    {
                        map[chromosome][i] = false;
                    }
                }

                //
                // Now run through all the exons again and populate the map.
                //
                foreach (var geneEntry in geneMap.genesByName)
                {
                    var gene = geneEntry.Value;

                    foreach (var isoform in gene.isoforms)
                    {
                        foreach (var exon in isoform.exons)
                        {
                            for (int i = exon.start; i <= exon.end; i++)
                            {
                                if (!map[gene.chromosome][i])
                                {
                                    nExonicLoci++;
                                    map[gene.chromosome][i] = true;
                                }
                            }
                        }
                    }
                }

            }

            public bool isLocationInAnExon(string chromosome, int location)
            {
                return map.ContainsKey(chromosome) && chromsomeMapSizes[chromosome] <= location && map[chromosome][location];
            }

            Dictionary<string, int> chromsomeMapSizes = new Dictionary<string,int>();
            Dictionary<string, bool[]> map = new Dictionary<string,bool[]>(); // Maps non-chr chromosome name->array of bools.  And yes, this could be done more space efficiently, but it's peanuts compared to a lot of stuff we do.

            int nExonicLoci = 0;

            public int getNExonicLoci() { return nExonicLoci; }
        }

        public class GeneInformation
        {
            public GeneInformation(string gene_with_protein_product_line)
            {
                var fields = gene_with_protein_product_line.Split('\t');

                if (fields.Count() != 48)
                {
                    Console.WriteLine("GeneInformation.GeneInfomation: wrong field count in input line, " + fields.Count() + " != 48");
                    return;
                }

                HgncId = ConvertToNonExcelString(fields[0]);
                symbol = ConvertToNonExcelString(fields[1]);
                name = ConvertToNonExcelString(fields[2]);
                locus_group = ConvertToNonExcelString(fields[3]);
                locus_type = ConvertToNonExcelString(fields[4]);
                status = ConvertToNonExcelString(fields[5]);
                location = ConvertToNonExcelString(fields[6]);
                location_sortable = ConvertToNonExcelString(fields[7]);
                alias_symbol = ConvertToNonExcelString(fields[8]);
                alias_name = ConvertToNonExcelString(fields[9]);
                prev_symbol = ConvertToNonExcelString(fields[10]);
                prev_name = ConvertToNonExcelString(fields[11]);
                gene_family = ConvertToNonExcelString(fields[12]);
                gene_family_id = ConvertToNonExcelString(fields[13]);
                date_approved_reserved = ConvertToNonExcelString(fields[14]);
                date_symbol_changed = ConvertToNonExcelString(fields[15]);
                date_name_changed = ConvertToNonExcelString(fields[16]);
                date_modified = ConvertToNonExcelString(fields[17]);
                entrez_id = ConvertToNonExcelString(fields[18]);
                ensembl_gene_id = ConvertToNonExcelString(fields[19]);
                vega_id = ConvertToNonExcelString(fields[20]);
                ucsc_id = ConvertToNonExcelString(fields[21]);
                ena = ConvertToNonExcelString(fields[22]);
                refseq_accession = ConvertToNonExcelString(fields[23]);
                ccds_id = ConvertToNonExcelString(fields[24]);
                uniprot_ids = ConvertToNonExcelString(fields[25]);
                pubmed_id = ConvertToNonExcelString(fields[26]);
                mgd_id = ConvertToNonExcelString(fields[27]);
                rgd_id = ConvertToNonExcelString(fields[28]);
                lsdb = ConvertToNonExcelString(fields[29]);
                cosmic = ConvertToNonExcelString(fields[30]);
                omim_id = ConvertToNonExcelString(fields[31]);
                mirbase = ConvertToNonExcelString(fields[32]);
                homeodb = ConvertToNonExcelString(fields[33]);
                snornabase = ConvertToNonExcelString(fields[34]);
                bioparadigms_slc = ConvertToNonExcelString(fields[35]);
                orphanet = ConvertToNonExcelString(fields[36]);
                pseudogene_dot_org = ConvertToNonExcelString(fields[37]);
                horde_id = ConvertToNonExcelString(fields[38]);
                merops = ConvertToNonExcelString(fields[39]);
                imgt = ConvertToNonExcelString(fields[40]);
                iuphar = ConvertToNonExcelString(fields[41]);
                kznf_gene_catalog = ConvertToNonExcelString(fields[42]);
                mamit_trnadb = ConvertToNonExcelString(fields[43]);
                cd = ConvertToNonExcelString(fields[44]);
                lncrnadb = ConvertToNonExcelString(fields[45]);
                enzyme_id = ConvertToNonExcelString(fields[46]);
                intermediate_filament_db = ConvertToNonExcelString(fields[47]);
            }
            public static Dictionary<string, GeneInformation> LoadFromFile(string filename)
            {
                var genes = new Dictionary<string, GeneInformation>();

                var reader = CreateStreamReaderWithRetry(filename);

                reader.ReadLine();  // Skip the header

                string line;
                while (null != (line = reader.ReadLine()))
                {
                    var gene = new GeneInformation(line);

                    genes.Add(gene.symbol, gene);
                }

                reader.Close();

                return genes;
            }

            public string HgncId;
            public string symbol;
            public string name;
            public string locus_group;
            public string locus_type;
            public string status;
            public string location;
            public string location_sortable;
            public string alias_symbol;
            public string alias_name;
            public string prev_symbol;
            public string prev_name;
            public string gene_family;
            public string gene_family_id;
            public string date_approved_reserved;
            public string date_symbol_changed;
            public string date_name_changed;
            public string date_modified;
            public string entrez_id;
            public string ensembl_gene_id;
            public string vega_id;
            public string ucsc_id;
            public string ena;
            public string refseq_accession;
            public string ccds_id;
            public string uniprot_ids;
            public string pubmed_id;
            public string mgd_id;
            public string rgd_id;
            public string lsdb;
            public string cosmic;
            public string omim_id;
            public string mirbase;
            public string homeodb;
            public string snornabase;
            public string bioparadigms_slc;
            public string orphanet;
            public string pseudogene_dot_org;
            public string horde_id;
            public string merops;
            public string imgt;
            public string iuphar;
            public string kznf_gene_catalog;
            public string mamit_trnadb;
            public string cd;
            public string lncrnadb;
            public string enzyme_id;
            public string intermediate_filament_db;

            public List<Isoform> isoforms = new List<Isoform>();
            public string chromosome = "Unknown";               // In the non-chr form
            public int minLocus = 0;                            // Min locus of any isoform
            public int maxLocus = 0;                            // Max locus of any isoform

        } // GeneInformation

    } // ExpressionTools

}
