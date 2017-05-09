using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using System.IO;
using System.Threading;
using System.IO.Compression;
using System.Net;
using System.Web;
using System.Runtime.Serialization;
using System.Runtime.Serialization.Json;
using System.Diagnostics;
using System.Runtime.InteropServices;

namespace ASELib
{
    public class ASETools
    {
        public const string urlPrefix = @"https://gdc-api.nci.nih.gov/";

        public const int GuidStringLength = 36;

        //
        // A Case is a person with TCGA data.
        //
        public class Case
        {
            //
            // Mandatory metadata that's created by GenerateCases.
            //
            public string case_id;
            public string normal_dna_file_id;
            public string tumor_dna_file_id; 
            public string normal_rna_file_id = "";    // This is optional
            public string tumor_rna_file_id;
            public string maf_file_id;
            public string methylation_file_id;
            public string copy_number_file_id = "";
            public string project_id;   // This is TCGA_<DiseaseType> for TCGA.
            public List<string> sample_ids = new List<string>();

            //
            // Pathnames for downloaded files.
            //
            public string normal_dna_filename = "";
            public string tumor_dna_filename = "";
            public string normal_rna_filename = "";
            public string tumor_rna_filename = "";
            public string methylation_filename = "";
            public string copy_number_filename = "";

            //
            // Sizes for downloaded files.
            //
            public long normal_dna_size = 0;
            public long tumor_dna_size = 0;
            public long normal_rna_size = 0;
            public long tumor_rna_size = 0;
            public long methylation_size = 0;
            public long copy_number_size = 0;

            //
            // Pathnames for derived files.
            //
            public string normal_dna_allcount_filename = "";
            public string tumor_dna_allcount_filename = "";
            public string normal_rna_allcount_filename = "";
            public string tumor_rna_allcount_filename = "";
            public string maf_filename = "";
            public string regional_expression_filename = "";
            public string gene_expression_filename = "";
            public string selected_variants_filename = "";
            public string dna_reads_at_selected_variants_filename = "";
            public string dna_reads_at_selected_variants_index_filename = "";
            public string rna_reads_at_selected_variants_filename = "";
            public string rna_reads_at_selected_variants_index_filename = "";
            public string annotated_selected_variants_filename = "";
            public string allele_specific_gene_expression_filename = "";
            public string tumor_dna_gene_coverage_filname = "";
            public string vcf_filename = "";
            public string extracted_maf_lines_filename = "";

            //
            // Checksums for downloaded files. The tumor DNA BAMs aren't included here
            // because they're BAM sliced and so don't have server-side md5s.
            //
            public string normal_rna_file_bam_md5 = "";
            public string normal_rna_file_bai_md5 = "";
            public string tumor_rna_file_bam_md5 = "";
            public string tumor_rna_file_bai_md5 = "";
            public string normal_dna_file_bam_md5 = "";
            public string normal_dna_file_bai_md5 = "";
            public string tumor_dna_file_bam_md5 = "";
            public string tumor_dna_file_bai_md5 = "";
            public string methylation_file_md5 = "";
            public string copy_number_file_md5 = "";

            //
            // The column numbers from the cases file for these fields.  They're used by C++ programs, which don't have access to the HeaderizedFile class,
            // so so just take the column numbers on the command line.  (That seemed better than compiling file format knowledge into them.)
            //
            static public int ProjectColumn = -1;
            static public int TumorRNAAllcountFilenameColumn = -1;

            public string disease()
            {
                if (project_id.Contains('-'))
                {
                    return project_id.Substring(project_id.LastIndexOf('-') + 1).ToLower();
                }

                if (project_id.Contains('_'))
                {
                    return project_id.Substring(project_id.LastIndexOf('_') + 1).ToLower();
                }

                return project_id.ToLower();
            }

            public static Dictionary<string, Case> LoadCases(string inputFilename)
            {
                if (!File.Exists(inputFilename))
                {
                    return null;   // Nothing to load because we haven't generated a cases file yet.
                }

                var wantedFields = new List<string>();

                foreach (var info in AllFields)
                {
                    wantedFields.Add(info.columnName);
                }

                var inputFile = CreateStreamReaderWithRetry(inputFilename);
                if (null == inputFile)
                {
                    return null;
                }

                var headerizedFile = new HeaderizedFile<Case>(inputFile, false, true, "", wantedFields);

                List<Case> listOfCases;
                Dictionary<string, int> fieldMappings;
                if (!headerizedFile.ParseFile(fromSaveFileLine, out listOfCases, out fieldMappings)) {
                    inputFile.Close();
                    return null;
                }

                ProjectColumn = fieldMappings["Project ID"];
                TumorRNAAllcountFilenameColumn = fieldMappings["Tumor RNA Allcount Filename"];

                inputFile.Close();

                var cases = new Dictionary<string, Case>();
                foreach (var case_ in listOfCases) {
                    cases.Add(case_.case_id, case_);
                }
 
                return cases;
            } // LoadCases

            public delegate string ColumnGetter(Case case_);
            public delegate void ColumnSetter(Case case_, string value);
            public delegate string ExpectedIdGetter(Case case_);
            public class FieldInformation
            {
                public FieldInformation(string columnName_, ColumnGetter getter_, ColumnSetter setter_, DerivedFile.Type type_, string extension_, ExpectedIdGetter idGetter_)
                {
                    columnName = columnName_;
                    getter = getter_;
                    setter = setter_;

                    //
                    // type, extension and id getter apply only to derived files.  For other fields, use the other constructor.
                    //
                    type = type_;
                    extension = extension_;
                    idGetter = idGetter_;
                }

                public FieldInformation(string columnName_, ColumnGetter getter_, ColumnSetter setter_)
                {
                    columnName = columnName_;
                    getter = getter_;
                    setter = setter_;
                }

                public string getValue(Case case_)
                {
                    return getter(case_);
                }

                public void setValue(Case case_, string value)
                {
                    setter(case_, value);
                }

                public string getExpectedId(Case case_)
                {
                    return idGetter(case_);
                }
 
                public readonly string columnName;
                ColumnGetter getter;
                ColumnSetter setter;
                ExpectedIdGetter idGetter = c => "";
                public readonly DerivedFile.Type type = DerivedFile.Type.Unknown;
                public readonly string extension = "";
            } // FieldInformation 

            public static FieldInformation[] AllFields = 
            {
                new FieldInformation("Case ID",                                         c => c.case_id, (c, v) => c.case_id = v),
                new FieldInformation("Normal DNA File ID",                              c => c.normal_dna_file_id, (c, v) => c.normal_dna_file_id = v),
                new FieldInformation("Tumor DNA File ID",                               c => c.tumor_dna_file_id, (c,v) => c.tumor_dna_file_id = v),
                new FieldInformation("Normal RNA File ID",                              c => c.normal_rna_file_id, (c,v) => c.normal_rna_file_id = v),
                new FieldInformation("Tumor RNA File ID",                               c => c.tumor_rna_file_id, (c,v) => c.tumor_rna_file_id = v),
                new FieldInformation("MAF File ID",                                     c => c.maf_file_id, (c,v) => c.maf_file_id = v),
                new FieldInformation("Methylation File ID",                             c => c.methylation_file_id, (c,v) => c.methylation_file_id = v),
                new FieldInformation("Copy Number File ID",                             c => c.copy_number_file_id, (c,v) => c.copy_number_file_id = v),
                new FieldInformation("Project ID",                                      c => c.project_id, (c,v) => c.project_id = v),
                new FieldInformation("Sample IDs",                                      c => c.sampleIdsInCommaSeparatedList(), (c,v) => c.sample_ids = v.Split(',').ToList()),

                new FieldInformation("Normal DNA Filename",                             c => c.normal_dna_filename, (c,v) => c.normal_dna_filename = v),
                new FieldInformation("Tumor DNA Filename",                              c => c.tumor_dna_filename, (c,v) => c.tumor_dna_filename = v),
                new FieldInformation("Normal RNA Filename",                             c => c.normal_rna_filename, (c,v) => c.normal_rna_filename = v),
                new FieldInformation("Tumor RNA Filename",                              c => c.tumor_rna_filename, (c,v) => c.tumor_rna_filename = v),
                new FieldInformation("Methylation Filename",                            c => c.methylation_filename, (c,v) => c.methylation_filename = v),
                new FieldInformation("Copy Number Filename",                            c => c.copy_number_filename, (c,v) => c.copy_number_filename = v),
                new FieldInformation("MAF Filename",                                    c => c.maf_filename, (c,v) => c.maf_filename = v),

                new FieldInformation("Normal DNA Size",                                 c => Convert.ToString(c.normal_dna_size), (c,v) => c.normal_dna_size = LongFromString(v)),
                new FieldInformation("Tumor DNA Size",                                  c => Convert.ToString(c.tumor_dna_size), (c,v) => c.tumor_dna_size = LongFromString(v)),
                new FieldInformation("Normal RNA Size",                                 c => Convert.ToString(c.normal_rna_size), (c,v) => c.normal_rna_size = LongFromString(v)),
                new FieldInformation("Tumor RNA Size",                                  c => Convert.ToString(c.tumor_rna_size), (c,v) => c.tumor_rna_size = LongFromString(v)),
                new FieldInformation("Methylation Size",                                c => Convert.ToString(c.methylation_size), (c,v) => c.methylation_size = LongFromString(v)),
                new FieldInformation("Copy Number Size",                                c => Convert.ToString(c.copy_number_size), (c,v) => c.copy_number_size = LongFromString(v)),

                new FieldInformation("Normal DNA Allcount Filename",                    c => c.normal_dna_allcount_filename, (c,v) => c.normal_dna_allcount_filename = v, DerivedFile.Type.NormalDNAAllcount, normalDNAAllcountExtension, c => c.normal_dna_file_id),
                new FieldInformation("Tumor DNA Allcount Filename",                     c => c.tumor_dna_allcount_filename, (c,v) => c.tumor_dna_allcount_filename = v, DerivedFile.Type.TumorDNAAllcount, tumorDNAAllcountExtension, c => c.tumor_dna_file_id),
                new FieldInformation("Normal RNA Allcount Filename",                    c => c.normal_rna_allcount_filename, (c,v) => c.normal_rna_allcount_filename = v, DerivedFile.Type.NormalRNAAllcount, normalRNAAllcountExtension, c => c.normal_rna_file_id),
                new FieldInformation("Tumor RNA Allcount Filename",                     c => c.tumor_rna_allcount_filename, (c,v) => c.tumor_rna_allcount_filename = v, DerivedFile.Type.TumorRNAAllcount, tumorRNAAllcountExtension, c => c.tumor_rna_file_id),
                new FieldInformation("Regional Expression Filename",                    c => c.regional_expression_filename, (c,v) => c.regional_expression_filename = v, DerivedFile.Type.RegionalExpression, regionalExpressionExtension, c => c.tumor_rna_file_id),
                new FieldInformation("Gene Expression Filename",                        c => c.gene_expression_filename, (c,v) => c.gene_expression_filename = v, DerivedFile.Type.GeneExpression, geneExpressionExtension, c => c.tumor_rna_file_id),
                new FieldInformation("Selected Variants Filename",                      c => c.selected_variants_filename, (c,v) => c.selected_variants_filename = v, DerivedFile.Type.SelectedVariants, selectedVariantsExtension, c => c.normal_dna_file_id),
                new FieldInformation("DNA Reads At Selected Variants Filename",         c => c.dna_reads_at_selected_variants_filename, (c,v) => c.dna_reads_at_selected_variants_filename = v, DerivedFile.Type.DNAReadsAtSelectedVariants, dnaReadsAtSelectedVariantsExtension, c => c.tumor_dna_file_id),
                new FieldInformation("DNA Reads At Selected Variants Index Filename",   c => c.dna_reads_at_selected_variants_index_filename, (c,v) => c.dna_reads_at_selected_variants_index_filename = v, DerivedFile.Type.DNAReadsAtSelectedVariantsIndex, dnaReadsAtSelectedVariantsIndexExtension, c => c.tumor_dna_file_id),
                new FieldInformation("RNA Reads At Selected Variants Filename",         c => c.rna_reads_at_selected_variants_filename, (c,v) => c.rna_reads_at_selected_variants_filename = v, DerivedFile.Type.RNAReadsAtSelectedVariants, rnaReadsAtSelectedVariantsExtension, c => c.tumor_rna_file_id),
                new FieldInformation("RNA Reads At Selected Variants Index Filename",   c => c.rna_reads_at_selected_variants_index_filename, (c,v) => c.rna_reads_at_selected_variants_index_filename = v, DerivedFile.Type.RNAReadsAtSelectedVariantsIndex, rnaReadsAtSelectedVariantsIndexExtension, c => c.tumor_rna_file_id),
                new FieldInformation("Annotated Selected Variants Filename",            c => c.annotated_selected_variants_filename, (c,v) => c.annotated_selected_variants_filename = v, DerivedFile.Type.AnnotatedSelectedVariants, annotatedSelectedVariantsExtension, c => c.case_id),
                new FieldInformation("Allele Specific Gene Expression Filename",        c => c.allele_specific_gene_expression_filename, (c,v) => c.allele_specific_gene_expression_filename = v, DerivedFile.Type.AlleleSpecificGeneExpression, alleleSpecificGeneExpressionExtension, c => c.case_id),
                new FieldInformation("Tumor DNA Gene Coverage Filename",                c => c.tumor_dna_gene_coverage_filname, (c,v) => c.tumor_dna_gene_coverage_filname = v, DerivedFile.Type.TumorDNAGeneCoverage, tumorDNAGeneCoverageExtension, c => c.tumor_dna_file_id),
                new FieldInformation("VCF Filename",                                    c => c.vcf_filename, (c,v) => c.vcf_filename = v, DerivedFile.Type.VCF, vcfExtension, c => c.normal_dna_file_id),
                new FieldInformation("Extracted MAF Lines Filename",                    c => c.extracted_maf_lines_filename, (c,v) => c.extracted_maf_lines_filename = v, DerivedFile.Type.ExtractedMAFLines, extractedMAFLinesExtension, c => c.case_id),

                new FieldInformation("Normal RNA BAM MD5",                              c => c.normal_rna_file_bam_md5, (c,v) => c.normal_rna_file_bam_md5 = v),
                new FieldInformation("Normal RNA BAI MD5",                              c => c.normal_rna_file_bai_md5, (c,v) => c.normal_rna_file_bai_md5 = v),
                new FieldInformation("Tumor RNA BAM MD5",                               c => c.tumor_rna_file_bam_md5, (c,v) => c.tumor_rna_file_bam_md5 = v),
                new FieldInformation("Tumor RNA BAI MD5",                               c => c.tumor_rna_file_bai_md5, (c,v) => c.tumor_rna_file_bai_md5 = v),
                new FieldInformation("Normal DNA BAM MD5",                              c => c.normal_dna_file_bam_md5, (c,v) => c.normal_dna_file_bam_md5 = v),
                new FieldInformation("Normal DNA BAI MD5",                              c => c.normal_dna_file_bai_md5, (c,v) => c.normal_dna_file_bai_md5 = v),
                new FieldInformation("Tumor DNA BAM MD5",                               c => c.tumor_dna_file_bam_md5, (c,v) => c.tumor_dna_file_bam_md5 = v),
                new FieldInformation("Tumor DNA BAI MD5",                               c => c.tumor_dna_file_bai_md5, (c,v) => c.tumor_dna_file_bai_md5 = v),
                new FieldInformation("Methylation MD5",                                 c => c.methylation_file_md5, (c,v) => c.methylation_file_md5 = v),
                new FieldInformation("Copy Number MD5",                                 c => c.copy_number_file_md5, (c,v) => c.copy_number_file_md5 = v),
            }; // fieldInformation

  
            static public Case fromSaveFileLine(Dictionary<string, int> fieldMappings, string[] fields) 
            {
                var case_ = new Case();

                foreach (var info in AllFields)
                {
                    info.setValue(case_, fields[fieldMappings[info.columnName]]);
                }
 
                return case_;
            } // fromSaveFile

            string sampleIdsInCommaSeparatedList()
            {
                if (sample_ids.Count() == 0) {
                    return "";
                }

                 var retVal = sample_ids[0];

                for (int i = 1; i < sample_ids.Count(); i++) {
                    retVal += "," + sample_ids[i];
                }

                return retVal;

            }   // sampleIdsInCommaSeparatedList    

            public string GenerateLine()
            {

                string value = AllFields[0].getValue(this);

                for (int i = 1; i < AllFields.Count(); i++)
                {
                    value += "\t" + AllFields[i].getValue(this);
                }

                return value;

            } // GenerateLine

            public static void SaveToFile(Dictionary<string, Case> cases, string filename)
            {
                var file = CreateStreamWriterWithRetry(filename);

                //
                // Header
                //
                for (int i = 0; i < AllFields.Count() - 1; i++)
                {
                    file.Write(AllFields[i].columnName + "\t");
                }

                file.WriteLine(AllFields[AllFields.Count() - 1].columnName);

                //
                // Cases
                //
                foreach (var case_ in cases)
                {
                    file.WriteLine(case_.Value.GenerateLine());
                }
                file.WriteLine("**done**");
                file.Close();
            } // SaveToFile

            public static void loadAllFileLocations(Dictionary<string, Case> cases, Dictionary<string, DownloadedFile> downloadedFiles, Dictionary<string, List<DerivedFile>> derivedFiles)
            {
                foreach (var caseEntry in cases)
                {
                    caseEntry.Value.loadFileLocations(downloadedFiles, derivedFiles);
                }
            }

            //
            // Fill in all the derived file locations if they exist.
            //
            public void loadFileLocations(Dictionary<string, DownloadedFile> downloadedFiles, Dictionary<string, List<DerivedFile>> derivedFiles)
            {

                if (downloadedFiles.ContainsKey(normal_dna_file_id)) {
                    normal_dna_filename = downloadedFiles[normal_dna_file_id].fileInfo.FullName;
                }
                else
                {
                    normal_dna_filename = "";
                }

                if (downloadedFiles.ContainsKey(tumor_dna_file_id))
                {
                  tumor_dna_filename = downloadedFiles[tumor_dna_file_id].fileInfo.FullName;
                }
                else
                {
                    tumor_dna_filename = "";
                }

                if (normal_rna_file_id != "" && downloadedFiles.ContainsKey(normal_rna_file_id))
                {
                  normal_rna_filename = downloadedFiles[normal_rna_file_id].fileInfo.FullName;
                }
                else
                {
                    normal_rna_filename = "";
                }

                if (downloadedFiles.ContainsKey(tumor_rna_file_id))
                {
                  tumor_rna_filename = downloadedFiles[tumor_rna_file_id].fileInfo.FullName;
                }
                else
                {
                    tumor_rna_filename = "";
                }

                if (downloadedFiles.ContainsKey(maf_file_id))
                {
                  maf_filename  = downloadedFiles[maf_file_id].fileInfo.FullName;
                }
                else
                {
                    maf_filename = "";
                }

                if (methylation_filename != "" && downloadedFiles.ContainsKey(methylation_file_id))
                {
                  methylation_filename  = downloadedFiles[methylation_file_id].fileInfo.FullName;
                }
                else
                {
                    methylation_filename = "";
                }

                if (copy_number_filename != "" && downloadedFiles.ContainsKey(copy_number_file_id))
                {
                    copy_number_filename = downloadedFiles[copy_number_file_id].fileInfo.FullName;
                }
                else
                {
                    copy_number_filename = "";
                }

                if (!derivedFiles.ContainsKey(case_id))
                {
                    tumor_rna_allcount_filename = "";
                    allele_specific_gene_expression_filename = "";
                    annotated_selected_variants_filename = "";
                    dna_reads_at_selected_variants_filename = "";
                    dna_reads_at_selected_variants_index_filename = "";
                    gene_expression_filename = "";
                    dna_reads_at_selected_variants_filename = "";
                    regional_expression_filename = "";
                    rna_reads_at_selected_variants_filename = "";
                    rna_reads_at_selected_variants_index_filename = "";
                    selected_variants_filename = "";
                    vcf_filename = "";
                    return;
                }

                Dictionary<DerivedFile.Type, List<DerivedFile>> derivedFilesByType = new Dictionary<DerivedFile.Type,List<DerivedFile>>();
                var derivedFilesForThisCase = derivedFiles[case_id];
                foreach (var type in (DerivedFile.Type[])Enum.GetValues(typeof(DerivedFile.Type))) {
                    derivedFilesByType.Add(type, derivedFilesForThisCase.Where(x => x.type == type).ToList());
                    if (type != DerivedFile.Type.Unknown && derivedFilesByType[type].Count() > 1) {
                        Console.Write("Found more than one file of type " + type + " for case " + case_id + ":");
                        foreach (var file in derivedFilesByType[type]) {
                            Console.Write(" " + file.fileinfo.FullName);
                        }
                        Console.WriteLine();
                    }
                }

                foreach (var field in AllFields)
                {
                    if (field.type == DerivedFile.Type.Unknown)
                    {
                        continue;
                    }

                    if (derivedFilesByType[field.type].Count() > 0)
                    {
                        var filename = derivedFilesByType[field.type][0].fileinfo.FullName;
                        field.setValue(this, filename);
                        if (GetFileIdFromPathname(filename) != field.getExpectedId(this))
                        {
                            Console.WriteLine("Found derived file with unexpected file ID: " + filename + ", expected file id: " + field.getExpectedId(this));
                            derivedFilesByType[field.type] = new List<DerivedFile>();   // Just blow it away and don't use it.
                        }
                    }
                    else
                    {
                        field.setValue(this, "");
                    }
                }
            } // loadFileLocations
        } // Case

        public static long LongFromString(string String)
        {
            if (String == null || String == "") return 0;

            try
            {
                return Convert.ToInt64(String);
            }
            catch (FormatException)
            {
                Console.WriteLine("LongFromString: string '" + String + "' didn't parse.");
                return 0;
            }
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
                catch (IOException e)
                {
                    if (e is FileNotFoundException)
                    {
                        return null;
                    }
                    Console.WriteLine("IOException opening " + filename + " for read.  Sleeping and retrying.");
                    Thread.Sleep(10 * 1000);
                }
            }
        }

        public static StreamReader CreateCompressedStreamReaderWithRetry(string filename)
        {
            var innerReader = CreateStreamReaderWithRetry(filename);
            if (innerReader == null) {
                return null;
            }
            return new StreamReader(new GZipStream(innerReader.BaseStream, CompressionMode.Decompress));
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

        public class ASEConfirguation
        {
            public const string defaultBaseDirectory = @"\\msr-genomics-0\d$\gdc\";
            public const string defaultConfigurationFilePathame = defaultBaseDirectory + "configuration.txt";

            public string accessTokenPathname = defaultBaseDirectory +  @"access_token.txt";
            public List<string> dataDirectories = new List<string>();
            public string mafManifestPathname = defaultBaseDirectory + "mafManifest.txt";
            public string mutationCaller = "mutect";
            public List<string> programNames = new List<string>();
            public string binariesDirectory = defaultBaseDirectory + @"bin\";
            public string configuationFilePathname = defaultConfigurationFilePathame;
            public string casesFilePathname = defaultBaseDirectory + "cases.txt";
            public string indexDirectory = @"d:\gdc\indices\hg38-20";
            public string derivedFilesDirectory = "derived_files";    // This is relative to each download directory
            public string hpcScriptFilename = "";    // The empty string says to black hole this script
            public string hpcScheduler = "gcr";
            public string hpcBinariesDirectory = @"\\gcr\scratch\b99\bolosky\";
            public string hpcIndexDirectory = @"\\msr-genomics-0\d$\gdc\indices\hg38-20";
            public string azureScriptFilename = ""; // The empty string says to black hole this script
            public string expressionFilesDirectory = @"\\msr-genomics-0\d$\gdc\expression\";
            public string completedVCFsDirectory = "";  // Where completed VCFs from Azure are dropped
            public int regionalExpressionRegionSize = 1000;

            public string[] commandLineArgs = null;    // The args excluding -configuration <filename>

            ASEConfirguation()
            {
                programNames.Add("TCGA");   // The default value
                dataDirectories.Add(defaultBaseDirectory + @"downloaded_files\");
            }

            //
            // Parse the args to find the configuration file pathname, and then load from that path (or the default if it's not present).
            //
            public static ASEConfirguation loadFromFile(string [] args) 
            {
                string pathname = @"\\msr-genomics-0\d$\gdc\configuration.txt";

                var nonConsumedArgs = new List<string>();
                for (int i = 0; i < args.Count(); i++) {
                    if (args[i] == "-configuration") {
                        if (i >= args.Count() - 1) {
                            Console.WriteLine("-configuation can't be the last parameter, it needs to be followed by a file path.  Ignoring.");
                        } else {
                            pathname = args[i+1];
                            i++;
                        }
                    }
                    else
                    {
                        nonConsumedArgs.Add(args[i]);
                    }
                }

                var retVal = new ASEConfirguation();
                retVal.commandLineArgs = nonConsumedArgs.ToArray();

                if (!File.Exists(pathname))
                {
                    //
                    // No config file means to use the default.
                    //
                    return retVal;
                }

                retVal.configuationFilePathname = pathname; // Don't load this from the file, just keep track of where we got it from.

                var lines = ReadAllLinesWithRetry(pathname);
                bool seenProgramNames = false;
                bool seenDataDirectories = false;

                foreach (var line in lines) 
                {
                    var fields = line.Split('\t');
                    if (fields.Count() != 2) {
                        Console.WriteLine("ASEConfiguration.loadFromFile: configuration file " + pathname + " contains a line that doesn't have exactly two tab separated fields: " + line + ".  Ignoring.");
                        continue;
                    }

                    var type = fields[0].ToLower();
                    if (type == "access token") {
                        retVal.accessTokenPathname = fields[1];
                    } else if (type == "data directory") {
                        if (!seenDataDirectories)   // Don't keep the default
                        {
                            retVal.dataDirectories = new List<string>();
                            seenDataDirectories = true;
                        }
                        retVal.dataDirectories.Add(fields[1]);
                    } else if (type == "maf manifest") {
                        retVal.mafManifestPathname = fields[1];
                    } else if (type == "mutation caller") {
                        retVal.mutationCaller = fields[1];
                    } else if (type == "program name") {
                        if (!seenProgramNames)  // If we've got the default value, override it.
                        {
                            retVal.programNames = new List<string>();
                            seenProgramNames = true;
                        }
                        retVal.programNames.Add(fields[1]);
                    } else if (type == "binary directory") {
                        retVal.binariesDirectory = fields[1];
                    } else if (type == "cases") {
                        retVal.casesFilePathname = fields[1];
                    } else if (type == "index directory") {
                        retVal.indexDirectory = fields[1];
                    } else if (type == "derived files") {
                        retVal.derivedFilesDirectory = fields[1];
                    } else if (type == "hpc script name") {
                        retVal.hpcScriptFilename = fields[1];
                    } else if (type == "hpc scheduler") {
                        retVal.hpcScheduler = fields[1];
                    } else if (type == "hpc binaries directory") {
                        retVal.hpcBinariesDirectory = fields[1];
                    } else if (type == "hpc index directory") {
                        retVal.hpcIndexDirectory = fields[1];
                    } else if (type == "expression files directory") {
                        retVal.expressionFilesDirectory = fields[1];
                    } else if (type == "azure script name") {
                        retVal.azureScriptFilename = fields[1];
                    } else if (type == "completed vcfs directory") {
                        retVal.completedVCFsDirectory = fields[1];
                    } else if (type == "regional expression region size") {
                        retVal.regionalExpressionRegionSize = Convert.ToInt32(fields[1]);   // should probably wrap this in a try
                    } else {
                        Console.WriteLine("ASEConfiguration.loadFromFile: configuration file " + pathname + " contains a line with an unknown configuration parameter type: " + line + ".  Ignoring.");
                        continue;
                    }
                }

                return retVal;
            }

        }

        public static System.Net.WebClient getWebClient()
        {
            ServicePointManager.SecurityProtocol = SecurityProtocolType.Tls12;  // NIH uses TLS v 1.2, and for some reason the auto negotiate doesn't pick that up, so we'll just set it explicitly

            var webClient = new WebClient();

            //
            // Check the status of gdc to make sure that we're in the right version.
            //

            var statusSerializer = new DataContractJsonSerializer(typeof(ASETools.Status));
            ASETools.Status status = (ASETools.Status)statusSerializer.ReadObject(new MemoryStream(webClient.DownloadData(ASETools.urlPrefix + "status")));

            if (status.status != "OK")
            {
                Console.WriteLine("GDC returned not OK status of " + status.status + ", aborting");
                return null;
            }

            if (status.version != "1")
            {
                Console.WriteLine("GDC returned a version number that we don't understand (" + status.version + ", you probably need to update this program.");
                return null;
            }

            return webClient;
        }

        [DataContract]
        public class Status
        {
            [DataMember]
            public string commit = "";  // Default values are to void having Visual Studio generate a warning, since it can't see that they're set by the JSON serializer

            [DataMember]
            public string status = "";

            [DataMember]
            public string tag = "";

            [DataMember]
            public string version = "";
        }

        [DataContract]
        public class GDCPagination
        {
            [DataMember]
            public int count = 0;

            [DataMember]
            public string sort = "";

            [DataMember]
            public int from = 0;

            [DataMember]
            public int page = 0;

            [DataMember]
            public int total = 0;

            [DataMember]
            public int pages = 0;

            [DataMember]
            public string size = "";

            static public GDCPagination extractFromString(string inputString)
            {
                int paginationIndex = inputString.IndexOf("\"pagination\":");

                if (-1 == paginationIndex)
                {
                    return null;
                }

                int openCurlyBraceIndex = paginationIndex;

                int size = inputString.Count();

                while (openCurlyBraceIndex < size && inputString[openCurlyBraceIndex] != '{')
                {
                    openCurlyBraceIndex++;
                }

                if (openCurlyBraceIndex >= size)
                {
                    return null;
                }

                int currentIndex = openCurlyBraceIndex + 1;
                int openCurlyBraceCount = 1;

                while (currentIndex < size && openCurlyBraceCount > 0)
                {
                    if (inputString[currentIndex] == '}')
                    {
                        openCurlyBraceCount--;
                    }
                    else if (inputString[currentIndex] == '{')
                    {
                        openCurlyBraceCount++;
                    }
                    currentIndex++;
                }

                if (0 != openCurlyBraceCount)
                {
                    throw new FormatException();
                }


                var paginationString = inputString.Substring(openCurlyBraceIndex, currentIndex - openCurlyBraceIndex);
                var serializer = new DataContractJsonSerializer(typeof(GDCPagination));
                return (GDCPagination)serializer.ReadObject(new MemoryStream(Encoding.ASCII.GetBytes(paginationString)));
            }
        }

        [DataContract]
        public class GDCCaseProject
        {
            [DataMember]
            public string project_id;
        }

        [DataContract]
        public class GDCCase
        {
            [DataMember]
            public string[] sample_ids = { };

            [DataMember]
            public string[] portion_ids = { };

            [DataMember]
            public string updated_datetime = "";

            [DataMember]
            public string created_datetime = "";

//            [DataMember]
//            public string[] submitter_aliquot_ids = { };

//            [DataMember]
//            public string[] submitter_portion_ids = { };

//            [DataMember]
//            public string[] submitter_analyte_ids = { };

//            [DataMember]
//            public string[] analyte_ids = { };

            [DataMember]
            public string submitter_id = "";

            [DataMember]
            public string case_id = "";

            [DataMember]
            public string state = "";

//            [DataMember]
//            public string[] aliquot_ids = { };

//            [DataMember]
//            public string[] slide_ids = { };

//            [DataMember]
//            public string[] submitter_sample_ids = { };

//            [DataMember]
//            public string[] submitter_slide_ids = { };

            [DataMember]
            public GDCCaseProject project = null;

            public const string fields = //"sample_ids,portion_ids,updated_datetime,created_datetime,submitter_aliquot_ids,submitter_portion_ids,submitter_analyte_ids,analyte_ids,submitter_id,case_id,state,aliquot_ids,slide_ids,submitter_sample_ids,submitter_slide_ids,project.project_id";
                "sample_ids,portion_ids,updated_datetime,created_datetime,submitter_id,case_id,state,project.project_id";

            public string debugString()
            {
                var retVal = "sample IDs (" + sample_ids.Count() + "): {";
                foreach (var sampleId in sample_ids)
                {
                    retVal += sampleId + " ";
                }
                retVal += "} updated datetime: " + updated_datetime + ", created datetime:  " + created_datetime + ", submitter_id: " + submitter_id + ", case_id: " + case_id + ", state: " + state + ", project: " + project.project_id;

                return retVal;
            }
        }

        [DataContract]
        public class GDCHits<containedClass>
        {
            [DataMember]
            public containedClass[] hits = { };
        }

        [DataContract]
        public class GDCData<containedClass>
        {
            [DataMember]
            public GDCHits<containedClass> data = null;
        }

         [DataContract]
        public class GDCSamples
        {
            [DataMember]
            public string sample_type_id = "";

            [DataMember]
            public string sample_type = "";

            [DataMember]
            public string sample_id = "";
        }

        [DataContract]
        public class GDCFileCases
        {
            [DataMember]
            public GDCSamples[] samples = { };
        }

        [DataContract]
        public class GDCFileAnalysis
        {
            [DataMember]
            public string workflow_link = "";
        }

        [DataContract]
        public class GDCFile
        {
            [DataMember]
            public string data_type = "";

            [DataMember]
            public string updated_datetime = "";

            [DataMember]
            public string created_datetime = "";

            [DataMember]
            public string file_name = "";

            [DataMember]
            public string md5sum = "";

            [DataMember]
            public string data_format = "";

            [DataMember]
            public string[] acl = { };

            [DataMember]
            public string access = "";

            [DataMember]
            public string platform = "";

            [DataMember]
            public string state = "";

            [DataMember]
            public string file_id = "";

            [DataMember]
            public string data_category = "";

            [DataMember]
            public long file_size = 0;

            [DataMember]
            public string submitter_id = "";

            [DataMember]
            public string type = "";

            [DataMember]
            public string file_state = "";

            [DataMember]
            public string experimental_strategy = "";

            [DataMember]
            public GDCFileCases[] cases = { };

            [DataMember]
            public GDCFileAnalysis analysis = null;

            public static GDCFile selectNewestUpdated(GDCFile one, GDCFile two)
            {
                if (null == one) return two;
                if (null == two) return one;

                if (Convert.ToDateTime(one.updated_datetime) < Convert.ToDateTime(two.updated_datetime))
                {
                    return two;
                }
                else
                {
                    return one;
                }
            }
        } // GDCFile


        public static string generateFilterList(List<string> itemsToAdd)
        {
            var outputString = "[";

            bool addedOneAlready = false;
            foreach (var item in itemsToAdd)
            {
                if (addedOneAlready)
                {
                    outputString += ",";
                }
                else
                {
                    addedOneAlready = true;
                }
                outputString += "\"" + item + "\"";
            }

            outputString += "]";

            return outputString;
        }

        public static string ElapsedTimeInSeconds(Stopwatch stopwatch)
        {
            return "" + (stopwatch.ElapsedMilliseconds + 500) / 1000 + "s";
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

        public static string GetFileNameFromPathname(string pathname, bool excludeExtension = false)
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

        public static string GetDirectoryFromPathname(string pathname)
        {
            if (!pathname.Contains('\\'))
            {
                throw new FormatException();
            }

            return pathname.Substring(0, pathname.LastIndexOf('\\'));
        }

        public static string GetFileIdFromPathname(string pathname) 
        {
            string filename = GetFileNameFromPathname(pathname);

            if (filename.Count() < FileIdLength)
            {
                throw new FormatException();
            }

            return filename.Substring(0, FileIdLength);
        }

        static readonly int FileIdLength = "01c5f902-ec3c-4eb7-9e38-1d29ae6ab959".Count();

        public class DownloadedFile
        {
            public readonly string file_id;
            public readonly FileInfo fileInfo;
            public readonly string storedMD5;    // This is the null string for files for which we don't know the sum (i.e., they haven't yet been computed or we don't know the right answer).
            public readonly FileInfo md5FileInfo;

            public DownloadedFile(string file_id_, string pathname, string storedMd5Sum_, string md5FilePathname)
            {
                file_id = file_id_;
                storedMD5 = storedMd5Sum_;

                fileInfo = new FileInfo(pathname);

                if (md5FilePathname == null || md5FilePathname == "")
                {
                    md5FileInfo = null;
                }
                else
                {
                    md5FileInfo = new FileInfo(md5FilePathname);
                }
            }
        } // DownloadedFile


        public const string tumorRNAAllcountExtension = ".allcount.gz";
        public const string normalRNAAllcountExtension = ".normal_rna_allcount.gz";
        public const string normalDNAAllcountExtension = ".normal_dna_allcount.gz";
        public const string tumorDNAAllcountExtension = ".tumor_dna_allcount.gz";
        public const string selectedVariantsExtension = ".selectedVariants";
        public const string annotatedSelectedVariantsExtension = ".annotatedSeletedVariants";
        public const string regionalExpressionExtension = ".regional_expression.txt";
        public const string geneExpressionExtension = ".gene_expression.txt";
        public const string alleleSpecificGeneExpressionExtension = ".allele-specific_gene_expression.txt";
        public const string rnaReadsAtSelectedVariantsExtension = ".rna-reads-at-selected-variants.txt";
        public const string rnaReadsAtSelectedVariantsIndexExtension = ".rna-reads-at-selected-variants.txt.index";
        public const string dnaReadsAtSelectedVariantsExtension = ".dna-reads-at-selected-variants.txt";
        public const string dnaReadsAtSelectedVariantsIndexExtension = ".dna-reads-at-selected-variants.txt.index";
        public const string vcfExtension = ".vcf";
        public const string extractedMAFLinesExtension = ".extracted_maf_lines.txt";
        public const string tumorDNAGeneCoverageExtension = ".tumor_dna_gene_coverage.txt";

        public class DerivedFile
        {
            public readonly string derived_from_file_id;
            public readonly Type type;
            public readonly FileInfo fileinfo;
            public readonly string case_id;

            public DerivedFile(string pathname, string case_id_)
            {
                fileinfo = new FileInfo(pathname);
                type = Type.Unknown;

                var filename = GetFileNameFromPathname(pathname).ToLower();
                case_id = case_id_;

                if (filename.Count() <= FileIdLength + 1) {
                    Console.WriteLine("Derived file has too short of a name: " + pathname);
                    derived_from_file_id = "";
                    return;
                }

                derived_from_file_id = filename.Substring(0,FileIdLength);

                var extension = filename.Substring(FileIdLength).ToLower();

                foreach (var field in Case.AllFields)
                {
                    if (field.type != Type.Unknown && extension == field.extension.ToLower())
                    {
                        type = field.type;
                    }
                }

                if (Type.Unknown == type)
                {
                    Console.WriteLine("Derived file with unknown extension: " + pathname);
                }
            } // DerviedFile.DerivedFile


            public enum Type { Unknown, NormalRNAAllcount, TumorRNAAllcount, NormalDNAAllcount, TumorDNAAllcount, RegionalExpression, GeneExpression, TumorDNAGeneCoverage,
            SelectedVariants, DNAReadsAtSelectedVariants, DNAReadsAtSelectedVariantsIndex, RNAReadsAtSelectedVariants, RNAReadsAtSelectedVariantsIndex, AnnotatedSelectedVariants, AlleleSpecificGeneExpression, VCF, ExtractedMAFLines
            };
        } // DerivedFile

        class ScanFilesystemState
        {
            public ulong totalFreeBytes = 0;
            public ulong totalBytesInDownloadedFiles = 0;
            public ulong totalBytesInDerivedFiles = 0;
            public int nDownloadedFiles = 0;
            public int nDerivedFiles = 0;

            public List<DownloadedFile> downloadedFiles = new List<DownloadedFile>();
            public List<DerivedFile> derivedFiles = new List<DerivedFile>();
        }

        static void ScanOneFilesystem(ASEConfirguation configuration, string downloadedFilesDirectory, ScanFilesystemState state, Stopwatch stopwatch, int directoryFieldLength) 
        {
            if (!Directory.Exists(downloadedFilesDirectory))
            {
                Console.WriteLine("Warning: can't find download directory " + downloadedFilesDirectory + ".  Ignoring.");
                return;
            }

            int nDownloadedFiles = 0;
            int nDerivedFiles = 0;
            ulong totalBytesInDownloadedFiles = 0;
            ulong totalBytesInDerivedFiles = 0;

            var downloadedFiles = new List<DownloadedFile>();
            var derivedFiles = new List<DerivedFile>();

            foreach (var subdir in Directory.EnumerateDirectories(downloadedFilesDirectory))
            {
                var file_id = GetFileNameFromPathname(subdir).ToLower();    // The directory is the same as the file id.

                if (file_id.Count() != GuidStringLength)
                {
                    lock (state)
                    {
                        Console.WriteLine("Found subdirectory of a downloaded files dierectory whose name doesn't appear to be a guid, ignoring: " + subdir);
                    }
                    continue;
                }

                //
                // Look through the subdirectory to find the downloaded file and also any .md5 versions of it.
                //
                string candidatePathname = null;
                string md5Pathname = null;
                foreach (var pathname in Directory.EnumerateFiles(subdir))
                {
                    var filename = GetFileNameFromPathname(pathname).ToLower();
                    if (filename == "annotations.txt")
                    {
                        continue;
                    }

                    if (filename.Count() > 4 && filename.Substring(filename.Count() - 4) == ".bai")
                    {
                        continue;
                    }


                    if (filename.Count() > 4 && filename.Substring(filename.Count() - 4) == ".md5")
                    {
                        if (null != md5Pathname)
                        {
                            Console.WriteLine("Directory " + subdir + " contains more than one file with extension .md5");
                            continue;
                        }
                        md5Pathname = pathname;
                        continue;
                    }

                    if (null != candidatePathname)
                    {
                        Console.WriteLine("Found more than one file candidate in " + subdir);
                        continue;
                    }

                    candidatePathname = pathname;
                }

                if (candidatePathname == null && md5Pathname != null)
                {
                    Console.WriteLine("Found md5 file in directory without accompanying downloaded file (?): " + md5Pathname);
                    continue;
                }

                string md5Value = "";
                if (md5Pathname != null)
                {
                    if (md5Pathname.ToLower() != candidatePathname.ToLower() + ".md5")
                    {
                        Console.WriteLine("md5 file has wrong basename: " + md5Pathname);
                    }
                    else
                    {
                        try
                        {
                            var lines = File.ReadAllLines(md5Pathname); // Don't use the with retry version, because we don't want to get stuck behind a worker hashing a file.
                            if (lines.Count() != 1 || lines[0].Count() != 32)
                            {
                                Console.WriteLine("md5 file " + md5Pathname + " has a non-md5 content.  Ignoring.");
                            }
                            else
                            {
                                md5Value = lines[0];
                            }
                        }
                        catch (IOException)
                        {
                            Console.WriteLine("Skipping md5 file " + md5Pathname + " because of an IO error.  It's most likely being created now.");
                        }
                    }
                }


                if (null == candidatePathname)
                {
                    Console.WriteLine("Couldn't find downloaded file in directory " + subdir);
                }
                else
                {
                    nDownloadedFiles++;
                    var downloadedFile = new DownloadedFile(file_id, candidatePathname, md5Value, md5Pathname);
                    downloadedFiles.Add(downloadedFile);

                    totalBytesInDownloadedFiles += (ulong)downloadedFile.fileInfo.Length;
                }
            } // foreach subdir

            string derivedFilesDirectory = downloadedFilesDirectory + @"..\" + configuration.derivedFilesDirectory;
            foreach (var derivedCase in Directory.EnumerateDirectories(derivedFilesDirectory))
            {
                var caseId = ASETools.GetFileNameFromPathname(derivedCase).ToLower();

                if (caseId.Count() != GuidStringLength)
                {
                    lock (state)
                    {
                        Console.WriteLine("Found subdirectory of a derived files dierectory whose name doesn't appear to be a guid, ignoring: " + derivedCase);
                    }
                    continue;
                }
 
                foreach (var derivedFilePathname in Directory.EnumerateFiles(derivedCase))
                {
                    var derivedFile = new DerivedFile(derivedFilePathname, caseId);
                    nDerivedFiles++;
                    totalBytesInDerivedFiles += (ulong)derivedFile.fileinfo.Length;
                    derivedFiles.Add(derivedFile);
                }
            }

            ulong freeBytesAvailable, totalBytes, totalNumberOfFreeBytes;
            GetDiskFreeSpaceEx(downloadedFilesDirectory, out freeBytesAvailable, out totalBytes, out totalNumberOfFreeBytes);

            lock (state)
            {
                Console.WriteLine(String.Format("{0," + directoryFieldLength + "}", downloadedFilesDirectory) + " " + String.Format("{0,16}", "" + nDownloadedFiles + " (" + SizeToUnits(totalBytesInDownloadedFiles) + "B)") + " " +
                    String.Format("{0,13}", "" + nDerivedFiles + " (" + SizeToUnits(totalBytesInDerivedFiles) + "B)") + " " + String.Format("{0,10}", SizeToUnits(freeBytesAvailable) +"B") + " " +
                    String.Format("{0,9}", ElapsedTimeInSeconds(stopwatch)));

                state.nDownloadedFiles += nDownloadedFiles;
                state.nDerivedFiles += nDerivedFiles;
                state.totalBytesInDownloadedFiles += totalBytesInDownloadedFiles;
                state.totalBytesInDerivedFiles += totalBytesInDerivedFiles;
                state.totalFreeBytes += freeBytesAvailable;
                state.downloadedFiles.AddRange(downloadedFiles);
                state.derivedFiles.AddRange(derivedFiles);
            }
        }

        //
        // Grovel the file system(s) to look for downloaded and derived files.  Returns a dictionary that maps from file_id -> DownloadedFile object.
        //
        public static void ScanFilesystems(ASEConfirguation configuration, out Dictionary<string, DownloadedFile> downloadedFiles, out Dictionary<string, List<DerivedFile>> derivedFiles)
        {
            var stopwatch = new Stopwatch();
            stopwatch.Start();

            int longestDirectoryNameLength = 0;

            foreach (var directory in configuration.dataDirectories)
            {
                longestDirectoryNameLength = Math.Max(longestDirectoryNameLength, directory.Count());
            }

            const string directoryHeader = "Directory";
            int paddingLength = Math.Max(0, longestDirectoryNameLength - directoryHeader.Count());

            Console.Write("Directory");
            for (int i = 0; i < paddingLength; i++)
            {
                Console.Write(" ");
            }

            Console.WriteLine(" Downloaded Files Derived Files Free Size Scan Time");
            for (int i = 0; i < directoryHeader.Count() + paddingLength; i++ )
            {
                Console.Write("-");
            }
            Console.WriteLine(" ---------------- ------------- ---------- ---------");

            var state = new ScanFilesystemState();
            var threads = new List<Thread>();

            foreach (var directory in configuration.dataDirectories)
            {
                threads.Add(new Thread(() => ScanOneFilesystem(configuration, directory, state, stopwatch, paddingLength + directoryHeader.Count())));
            } // foreach data directory

            threads.ForEach(t => t.Start());
            threads.ForEach(t => t.Join());

            downloadedFiles = new Dictionary<string, DownloadedFile>();
            derivedFiles = new Dictionary<string, List<DerivedFile>>();

            foreach (var downloadedFile in state.downloadedFiles)
            {
                if (downloadedFiles.ContainsKey(downloadedFile.file_id))
                {
                    Console.WriteLine("Found multiple directories with the downloaded file id " + downloadedFile.file_id + ": " + downloadedFile.fileInfo.FullName + " and " + downloadedFiles[downloadedFile.file_id].fileInfo.FullName);
                    continue;
                }
                downloadedFiles.Add(downloadedFile.file_id, downloadedFile);
            }

            foreach (var derivedFile in state.derivedFiles)
            {
                if (!derivedFiles.ContainsKey(derivedFile.case_id))
                {
                    derivedFiles.Add(derivedFile.case_id, new List<DerivedFile>());
                }

                derivedFiles[derivedFile.case_id].Add(derivedFile);
            }

            Console.WriteLine("Scanned " + configuration.dataDirectories.Count() + " data directories in " + ElapsedTimeInSeconds(stopwatch) + ", containing " + state.nDownloadedFiles + " (" + SizeToUnits(state.totalBytesInDownloadedFiles) +
                "B) downloaded and " + state.nDerivedFiles + " (" + SizeToUnits(state.totalBytesInDerivedFiles) + "B) derived files.  " + SizeToUnits(state.totalFreeBytes) + "B remain free.");


        } // ScanFilesystems

        public class MAFInfo
        {
            public MAFInfo(string file_id_, string md5Sum_, string filename_)
            {
                file_id = file_id_;
                md5Sum = md5Sum_;
                filename = filename_;
            }
            public static MAFInfo fromSaveFileLine(string saveFileLine)
            {
                var fields = saveFileLine.Split('\t');

                if (fields.Count() != 3)
                {
                    Console.WriteLine("Incorrect number of fields in maf configuration file line: " + saveFileLine);
                    return null;
                }

                return new MAFInfo(fields[0], fields[1], fields[2]);
            }

            public static Dictionary<string, MAFInfo> LoadMAFManifest(string filename)
            {
                if (!File.Exists(filename))
                {
                    return null;
                }

                var lines = ReadAllLinesWithRetry(filename);

                if (lines.Count() < 2)
                {
                    Console.WriteLine("MAF manifest file " + filename + " has too few lines.  Ignoring.");
                    return null;
                }

                var headerPrefix = "MAF Manifest v1.0 generated at ";
                if (lines[0].Count() < headerPrefix.Count() || lines[0].Substring(0,headerPrefix.Count()) != headerPrefix) {
                    Console.WriteLine("Corrupt or unrecognized version in maf manifest header, ignoring: " + lines[0]);
                    return null;
                }

                var retVal = new Dictionary<string, MAFInfo>();

                if (lines[lines.Count() -1] != "**done**") {
                    Console.WriteLine("maf manifest file " + filename + " does not end with **done**, and so is probably truncated.  Ignoring.");
                    return null;
                }

                for (int i = 1; i < lines.Count() - 1; i++) {
                    var mafInfo = MAFInfo.fromSaveFileLine(lines[i]);

                    if (null == mafInfo) {
                        Console.WriteLine("failed to parse maf manifest line " + lines[i] + ".  Ignoring configuration file.");
                        return null;
                    }

                    if (retVal.ContainsKey(mafInfo.file_id)) {
                        Console.WriteLine("Duplicate file id " + mafInfo.file_id + " in maf manifest file " + filename);
                        return null;
                    }

                    retVal.Add(mafInfo.file_id, mafInfo);
                }

                return retVal;
            }

            public string file_id;
            public string md5Sum;
            public string filename;
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

        public class HeaderizedFile<outputType>
        {
            public delegate outputType Parse(Dictionary<string, int> fieldMappings, string[] fields);
            public HeaderizedFile(StreamReader inputFile_, bool hasVersion_, bool hasDone_, string expectedVersion_, List<string> wantedFields_)
            {
                inputFile = inputFile_;
                hasVersion = hasVersion_;
                hasDone = hasDone_;
                expectedVersion = expectedVersion_;
                wantedFields = wantedFields_;
            }

            //
            // Just eat the fieldMappings output the clunky way, since out parameters can't have default values.
            //
            public bool ParseFile(Parse parser, out List<outputType> result)
            {
                Dictionary<string, int> fieldMappings;

                return ParseFile(parser, out result, out fieldMappings);
            }

            public bool ParseFile(Parse parser, out List<outputType> result, out Dictionary<string, int> fieldMappings_out)
            {
                fieldMappings_out = null;

                if (hasVersion)
                {
                    var versionString = inputFile.ReadLine();

                    if (versionString == null || expectedVersion != null && versionString != expectedVersion)
                    {
                        result = null;
                        return false;
                    }
                }

                var header = inputFile.ReadLine();
                if (null == header)
                {
                    result = null;
                    return false;
                }

                var columns = header.Split('\t');
                var fieldMappings = new Dictionary<string, int>();
                int maxNeededField = -1;

                for (int i = 0; i < columns.Count(); i++)
                {
                    if (wantedFields.Contains(columns[i]))
                    {
                        if (fieldMappings.ContainsKey(columns[i]))
                        {
                            Console.WriteLine("Duplicate needed column in headerized file (or code bug or something): " + columns[i]);
                            result = null;
                            return false;
                        }
                        
                        fieldMappings.Add(columns[i], i);
                        maxNeededField = i;
                    }
                }

                if (fieldMappings.Count() != wantedFields.Count())
                {
                    var missingColumns = new List<string>();
                    foreach (var wantedField in wantedFields)
                    {
                        if (!fieldMappings.ContainsKey(wantedField))
                        {
                            missingColumns.Add(wantedField);
                        }
                    }

                    if (fieldMappings.Count() + missingColumns.Count() != wantedFields.Count())
                    {
                        Console.WriteLine("Got the wrong number of missing fields.  Code bug.");
                        result = null;
                        return false;

                    }
                    Console.Write("Headerized file: missing columns:");
                    foreach (var missingColumn in missingColumns)
                    {
                        Console.Write(" " + missingColumn);
                        fieldMappings.Add(missingColumn, maxNeededField + 1);
                    }
                    Console.WriteLine(".  Filling in with the empty string.");
                    hasMissingFields = true;

                }

                string inputLine;
                bool sawDone = false;
                result = new List<outputType>();
                while (null != (inputLine = inputFile.ReadLine())) {
                    if (sawDone) {
                        Console.WriteLine("HeaderizedFile: Saw data after **done**");
                        result = null;
                        return false;
                    }

                    if ("**done**" == inputLine) {
                        sawDone = true;
                        continue;
                    }

                    var fields = inputLine.Split('\t');
                    if (fields.Count() <= maxNeededField)
                    {
                        Console.WriteLine("HeaderizedFile.Parse: input line didn't include a needed field " + inputLine);
                        result = null;
                        return false;
                    }
                    else if (hasMissingFields && fields.Count() <= maxNeededField + 1)
                    {
                        var extendedFields = new string[maxNeededField + 2];
                        for (int i = 0; i <= maxNeededField; i++)
                        {
                            extendedFields[i] = fields[i];
                        }
                        fields = extendedFields;
                    }

                    if (hasMissingFields)
                    {
                        fields[maxNeededField + 1] = "";    // This is for all missing fields
                    }

                    result.Add(parser(fieldMappings, fields));
                }

                if (hasDone && !sawDone)
                {
                    Console.WriteLine("HeaderizedFile.Parse: missing **done**");
                    result = null;
                    return false;
                }
                else if (!hasDone && sawDone)
                {
                    Console.WriteLine("Saw unepected **done**.  Ignoring.");
                    result = null;
                    return false;
                }


                fieldMappings_out = fieldMappings;
                return true;
            } // ParseFile

            StreamReader inputFile;
            bool hasVersion;
            bool hasDone;
            string expectedVersion;
            List<string> wantedFields;
            bool hasMissingFields = false;
        } // HeaderizedFile

        public class MAFLine
        {
            public readonly string Hugo_Symbol;
            public readonly string NCBI_Build;
            public readonly string Chromosome;
            public readonly int Start_Position;
            public readonly int End_Positon;
            public readonly string Variant_Classification;
            public readonly string Variant_Type;
            public readonly string Reference_Allele;
            public readonly string Tumor_Seq_Allele1;
            public readonly string Tumor_Seq_Allele2;
            public readonly string Match_Norm_Seq_Allele1;
            public readonly string Match_Norm_Seq_Allele2;
            public readonly string Tumor_Sample_UUID;
            public readonly string Matched_Norm_Sample_UUID;
            public readonly string tumor_bam_uuid;
            public readonly string normal_bam_uuid;

            public readonly string maf_file_id; // Of the MAF file

            MAFLine(string Hugo_Symbol_, 
             string NCBI_Build_,
             string Chromosome_,
             int Start_Position_,
             int End_Positon_,
             string Variant_Classification_,
             string Variant_Type_,
             string Reference_Allele_,
             string Tumor_Seq_Allele1_,
             string Tumor_Seq_Allele2_,
             string Match_Norm_Seq_Allele1_,
             string Match_Norm_Seq_Allele2_,
             string Tumor_Sample_UUID_,
             string Matched_Norm_Sample_UUID_,
             string tumor_bam_uuid_,
             string normal_bam_uuid_,
             string maf_file_id_)
            {
                Hugo_Symbol = Hugo_Symbol_;
                NCBI_Build = NCBI_Build_;
                Chromosome = Chromosome_;
                Start_Position = Start_Position_;
                End_Positon = End_Positon_;
                Variant_Classification = Variant_Classification_;
                Variant_Type = Variant_Type_;
                Reference_Allele = Reference_Allele_;
                Tumor_Seq_Allele1 = Tumor_Seq_Allele1_;
                Tumor_Seq_Allele2_ = Tumor_Seq_Allele2;
                Match_Norm_Seq_Allele1 = Match_Norm_Seq_Allele1_;
                Match_Norm_Seq_Allele2 = Match_Norm_Seq_Allele2_;
                Tumor_Sample_UUID = Tumor_Sample_UUID_;
                Matched_Norm_Sample_UUID = Matched_Norm_Sample_UUID_;
                tumor_bam_uuid = tumor_bam_uuid_;
                normal_bam_uuid = normal_bam_uuid_;
                maf_file_id = maf_file_id_;
            }

            static MAFLine ParseLine(Dictionary<string, int> fieldMappings, string[] fields, string maf_file_id)
            {
                return new MAFLine(
                    fields[fieldMappings["Hugo_Symbol"]],
                    fields[fieldMappings["NCBI_Build"]],
                    fields[fieldMappings["Chromosome"]],
                    Convert.ToInt32(fields[fieldMappings["Start_Position"]]),
                    Convert.ToInt32(fields[fieldMappings["End_Position"]]),
                    fields[fieldMappings["Variant_Classification"]],
                    fields[fieldMappings["Variant_Type"]],
                    fields[fieldMappings["Reference_Allele"]],
                    fields[fieldMappings["Tumor_Seq_Allele1"]],
                    fields[fieldMappings["Tumor_Seq_Allele2"]],
                    fields[fieldMappings["Match_Norm_Seq_Allele1"]],
                    fields[fieldMappings["Match_Norm_Seq_Allele2"]],
                    fields[fieldMappings["Tumor_Sample_UUID"]],
                    fields[fieldMappings["Matched_Norm_Sample_UUID"]],
                    fields[fieldMappings["tumor_bam_uuid"]],
                    fields[fieldMappings["normal_bam_uuid"]],
                    maf_file_id
                    );
            }

            static public List<MAFLine> ReadFile(string filename, string file_id, bool fileHasVersion)
            {
                StreamReader inputFile;

                if (filename.Count() > 2 && filename.Substring(filename.Count() - 3, 3) == ".gz")
                {
                    inputFile = CreateCompressedStreamReaderWithRetry(filename);
                }
                else
                {
                    inputFile = CreateStreamReaderWithRetry(filename);
                }

                var neededFields = new List<string>();
                neededFields.Add("Hugo_Symbol");
                neededFields.Add("NCBI_Build");
                neededFields.Add("Chromosome");
                neededFields.Add("Start_Position");
                neededFields.Add("End_Position");
                neededFields.Add("Variant_Classification");
                neededFields.Add("Variant_Type");
                neededFields.Add("Reference_Allele");
                neededFields.Add("Tumor_Seq_Allele1");
                neededFields.Add("Tumor_Seq_Allele2");
                neededFields.Add("Match_Norm_Seq_Allele1");
                neededFields.Add("Match_Norm_Seq_Allele2");
                neededFields.Add("Tumor_Sample_UUID");
                neededFields.Add("Matched_Norm_Sample_UUID");
                neededFields.Add("tumor_bam_uuid");
                neededFields.Add("normal_bam_uuid");

                var headerizedFile = new HeaderizedFile<MAFLine>(inputFile, fileHasVersion, !fileHasVersion, "#version 2.4", neededFields);

                List<MAFLine> result;

                if (!headerizedFile.ParseFile((a, b) => ParseLine(a, b, file_id), out result))
                {
                    Console.WriteLine("Error reading MAF File " + filename);
                    return null;
                }

                inputFile.Close();

                return result;
            } // ReadFile

            public static void WriteHeaderLine(StreamWriter outputStream)
            {
                outputStream.WriteLine("Hugo_Symbol\tNCBI_Build\tChromosome\tStart_Position\tEnd_Position\tVariant_Classification\tVariant_Type\tReference_Allele\tTumor_Seq_Allele1\tTumor_Seq_Allele2\tMatch_Norm_Seq_Allele1\tMatch_Norm_Seq_Allele2\tTumor_Sample_UUID\tMatched_Norm_Sample_UUID\ttumor_bam_uuid\tnormal_bam_uuid");
            }

            public void WriteToStream(StreamWriter output)
            {
                output.WriteLine(
                    Hugo_Symbol                 + "\t" +
                    NCBI_Build                  + "\t" +
                    Chromosome                  + "\t" +
                    Start_Position              + "\t" +
                    End_Positon                 + "\t" +
                    Variant_Classification      + "\t" +
                    Reference_Allele            + "\t" +
                    Reference_Allele            + "\t" +
                    Tumor_Seq_Allele1           + "\t" +
                    Tumor_Seq_Allele2           + "\t" +
                    Match_Norm_Seq_Allele1      + "\t" +
                    Match_Norm_Seq_Allele2      + "\t" +
                    Tumor_Sample_UUID           + "\t" +
                    Matched_Norm_Sample_UUID    + "\t" +
                    tumor_bam_uuid              + "\t" +
                    normal_bam_uuid
                    );
            }

            public static void WriteToFile(string filename, List<MAFLine> linesToWrite)
            {
                var output = CreateStreamWriterWithRetry(filename);

                WriteHeaderLine(output);

                foreach (var mafLine in linesToWrite)
                {
                    mafLine.WriteToStream(output);                    
                }

                output.WriteLine("**done**");
                output.Close();
            }

        } // MAFLine

        public static string SizeToUnits(ulong size)
        {
            if (size < 1024) return "" + size;

            if (size < 1024 * 1024) return "" + ((size + 512) / 1024) + "K";

            if (size < 1024 * 1024 * 1024) return "" + ((size + 512 * 1024) / (1024 * 1024)) + "M";

            if (size < (ulong)1024 * 1024 * 1024 * 1024) return "" + ((size + 512 * 1024 * 1024) / (1024 * 1024 * 1024)) + "G";

            if (size < (ulong)1024 * 1024 * 1024 * 1024 * 1024) return "" + ((size + (ulong)512 * 1024 * 1204 * 1024) / ((ulong)1024 * 1024 * 1024 * 1024)) + "T";

            return "" + ((size + (ulong)512 * 1024 * 1204 * 1024 * 1024) / ((ulong)1024 * 1024 * 1024 * 1024 * 1024)) + "P";
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
                    compressedStreamReader = CreateStreamReaderWithRetry(filename);
                    allcountReader = new StreamReader(new GZipStream(compressedStreamReader.BaseStream, CompressionMode.Decompress));


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
                } // try
                catch (Exception e)
                {
                    if (e is IOException || e is InvalidDataException)
                    {
                        Console.WriteLine("IOException or InvalidDataException opening allcount file " + filename);
                        return false;
                    }
                    else
                    {
                        throw e;
                    }
                }

                return true;
            } // openFile

            public delegate void ProcessBase(string strippedContigName, int location, int currentMappedReadCount);

            public bool ReadAllcountFile(ProcessBase processBase)
            {
                int currentOffset = -1;
                int currentMappedReadCount = -1;
                int whichContig = -1;

                bool sawDone = false;
                string contigName = "";

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


                        contigName = line.Substring(1).ToLower();

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
                                if (ex is FormatException || ex is ArgumentOutOfRangeException)
                                {
                                    Console.WriteLine("Format exception processing x line " + line);
                                    return false;
                                }
                                else
                                {
                                    throw ex;
                                }
                            }
                            for (; repeatCount > 1; repeatCount--)  // > 1 because this count includes the locus that specified the mapped read count (the previous non-x line) and we already emitted that one.
                            {
                                currentOffset++;
                                processBase(contigName, currentOffset, currentMappedReadCount);
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
                            processBase(contigName, currentOffset, currentMappedReadCount);
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

                    processBase(contigName, currentOffset, currentMappedReadCount);
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

        } // AllcountReader

        public static string DownloadStringWithRetry(WebClient webclient, string address)
        {
            while (true)
            {
                try
                {
                    return webclient.DownloadString(address);
                }
                catch (WebException e)
                {
                    Console.WriteLine("DownloadStringWithRetry: caught WebException: " + e.Message + ".  Pausing 10 seconds and retrying.");
                    Thread.Sleep(10000);
                }
            }
        } // DownloadStringWithRetry

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
        } // WindowsToLinuxPathname

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
        } // GetDirectoryPathFromFullyQualifiedFilename

        class MAFLoadStatus
        {
            public int nToLoad;
            public int nLoaded = 0;
            public int nFailed = 0;
        } // MAFLoadStatus

        static void ReadMafFileAndAppendToList(string filename, string file_id, List<ASETools.MAFLine> allMAFLines, MAFLoadStatus loadStatus)
        {
            var mafLinesForThisFile = ASETools.MAFLine.ReadFile(filename, file_id, true);

            lock (loadStatus)
            {
                loadStatus.nLoaded++;

                Console.Write(loadStatus.nLoaded + "/" + loadStatus.nToLoad + " ");

                if (null == mafLinesForThisFile)
                {
                    loadStatus.nFailed++;
                    Console.WriteLine("MAF load of " + filename + " failed.");
                }

                allMAFLines.AddRange(mafLinesForThisFile);
            }
        } // ReadMafFileAndAppendToList

        public static List<MAFLine> LoadMAFs(ASEConfirguation configuration, Dictionary<string, DownloadedFile> downloadedFiles, out Dictionary<string, List<MAFLine>> byTumorSampleId)
        {
            var stopwatch = new Stopwatch();
            stopwatch.Start();

            byTumorSampleId = null;

            if (configuration.mafManifestPathname == null || configuration.mafManifestPathname == "")
            {
                Console.WriteLine("The MAF manifest must exist in order to run GenerateCases.  Go back and run the process manager to generate the script to make the maf manifest.");
                return null;
            }

            var mafManifest = ASETools.MAFInfo.LoadMAFManifest(configuration.mafManifestPathname);
            if (null == mafManifest)
            {
                Console.WriteLine("Unable to load maf manifest " + configuration.mafManifestPathname);
                return null;
            }

 
            foreach (var mafEntry in mafManifest)
            {
                var mafInfo = mafEntry.Value;

                if (!downloadedFiles.ContainsKey(mafInfo.file_id))
                {
                    Console.WriteLine("Missing MAF file " + mafInfo.file_id + ".  Download it (which should be in the script generated by ASEProcessManager) and rerun.");
                    return null;
                }
            }

            Console.Write("Loading MAF entries..");
            var allMAFLines = new List<MAFLine>();

            var threads = new List<Thread>();

            var loadStatus = new MAFLoadStatus();
            loadStatus.nToLoad = mafManifest.Count();


            foreach (var mafEntry in mafManifest)
            {
                var mafInfo = mafEntry.Value;

                threads.Add(new Thread(() => ReadMafFileAndAppendToList(downloadedFiles[mafInfo.file_id].fileInfo.FullName, mafInfo.file_id, allMAFLines, loadStatus)));
            }
 

            threads.ForEach(t => t.Start());
            threads.ForEach(t => t.Join());

            if (loadStatus.nFailed > 0)
            {
                Console.WriteLine("Giving up due to failed MAF load(s).");
                return null;
            }

            byTumorSampleId = new Dictionary<string,List<MAFLine>>();

            foreach (MAFLine mafLine in allMAFLines) {
                if (!byTumorSampleId.ContainsKey(mafLine.Tumor_Sample_UUID))
                {
                    byTumorSampleId.Add(mafLine.Tumor_Sample_UUID, new List<MAFLine>());
                }

                byTumorSampleId[mafLine.Tumor_Sample_UUID].Add(mafLine);
            }

            Console.WriteLine(ElapsedTimeInSeconds(stopwatch) + " to load " + allMAFLines.Count() + " mutations.");

            return allMAFLines;
        } // LoadMAFs

        public struct MeanAndStdDev
        {
            public double mean;
            public double stddev;
        }

        //
        // Maps chromosomeName -> (offset -> MeanAndStdDev)
        //
        public static Dictionary<string, Dictionary<int, MeanAndStdDev>> LoadExpressionFile(string filename)
        {
            var expression = new Dictionary<string, Dictionary<int, MeanAndStdDev>>();

            StreamReader reader = new StreamReader(filename);

            string currentChromosome = null;
            Dictionary<int, MeanAndStdDev> currentMapping = null;

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
        } // LoadExpressionFile

    } // ASETools
}
