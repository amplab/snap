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
            public string project_id;   // This is TCGA_<DiseaseType> for TCGA.

            //
            // Pathnames for downloaded files.
            //
            public string normal_dna_filename = "";
            public string tumor_dna_filename = "";
            public string normal_rna_filename = "";
            public string tumor_rna_filename = "";
            public string methylation_filename = "";

            //
            // Pathnames for derived files.
            //
            //public string tumor_dna_allcount_filename = "";
            public string tumor_rna_allcount_filename = "";
            public string maf_filename = "";
            public string regional_expression_filename = "";
            public string gene_expression_filename = "";
            public string selected_variants_filename = "";
            public string dna_reads_at_selected_variants_filename = "";
            public string rna_reads_at_selected_variants_filename = "";
            public string annotated_selected_variants_filename = "";
            public string allele_specific_gene_expression_filename = "";
            public string tumor_dna_gene_coverage_filname = "";
            public string vcf_filename = "";

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
            public string methylation_file_md5 = "";

            public static Dictionary<string, Case> LoadCases(string inputFilename)
            {
                if (!File.Exists(inputFilename))
                {
                    return null;   // Nothing to load because we haven't generated a cases file yet.
                }

                var wantedFields = new List<string>();
                wantedFields.Add("Case ID");
                wantedFields.Add("Normal DNA File ID");
                wantedFields.Add("Tumor DNA File ID");
                wantedFields.Add("Normal RNA File ID");
                wantedFields.Add("Tumor RNA File ID");
                wantedFields.Add("MAF File ID");
                wantedFields.Add("Methylation File ID");
                wantedFields.Add("Project ID");
                wantedFields.Add("Normal DNA Filename");
                wantedFields.Add("Tumor DNA Filename");
                wantedFields.Add("Normal RNA Filename");
                wantedFields.Add("Tumor RNA Filename");
                wantedFields.Add("VCF Filename");
                wantedFields.Add("Methylation Filename");
                //wantedFields.Add("Tumor DNA Allcount Filename"); // We're not doing this
                wantedFields.Add("Tumor RNA Allcount Filename");
                wantedFields.Add("MAF Filename");
                wantedFields.Add("Regional Expression Filename");
                wantedFields.Add("Gene Expression Filename");
                wantedFields.Add("Selected Variants Filename");
                wantedFields.Add("DNA Reads At Selected Variants Filename");
                wantedFields.Add("RNA Reads At Selected Variants Filename");
                wantedFields.Add("Annotated Selected Variants Filename");
                wantedFields.Add("Allele Specific Gene Expression Filename");
                wantedFields.Add("Tumor DNA Gene Coverage Filename");
                wantedFields.Add("Normal RNA BAM MD5");
                wantedFields.Add("Normal RNA BAI MD5");
                wantedFields.Add("Tumor RNA BAM MD5");
                wantedFields.Add("Tumor RNA BAI MD5");
                wantedFields.Add("Normal DNA BAM MD5");
                wantedFields.Add("Normal DNA BAI MD5");
                wantedFields.Add("Methylation MD5");

                var inputFile = CreateStreamReaderWithRetry(inputFilename);
                var headerizedFile = new HeaderizedFile<Case>(inputFile, false, true, "", wantedFields);

                List<Case> listOfCases;
                if (!headerizedFile.ParseFile(fromSaveFileLine, out listOfCases)) {
                    inputFile.Close();
                    return null;
                }

                inputFile.Close();

                var cases = new Dictionary<string, Case>();
                foreach (var case_ in listOfCases) {
                    cases.Add(case_.case_id, case_);
                }
 
                return cases;
            } // LoadCases

            static public Case fromSaveFileLine(Dictionary<string, int> fieldMappings, string[] fields)
            {
                var case_ = new Case();

                case_.case_id = fields[fieldMappings["Case ID"]];
                case_.normal_dna_file_id = fields[fieldMappings["Normal DNA File ID"]];
                case_.tumor_dna_file_id = fields[fieldMappings["Tumor DNA File ID"]];
                case_.normal_rna_file_id = fields[fieldMappings["Normal RNA File ID"]];
                case_.tumor_rna_file_id = fields[fieldMappings["Tumor RNA File ID"]];
                case_.maf_file_id = fields[fieldMappings["MAF File ID"]];
                case_.methylation_file_id = fields[fieldMappings["Methylation File ID"]];
                case_.project_id = fields[fieldMappings["Project ID"]];
                case_.normal_dna_filename = fields[fieldMappings["Normal DNA Filename"]];
                case_.tumor_dna_filename = fields[fieldMappings["Tumor DNA Filename"]];
                case_.normal_rna_filename = fields[fieldMappings["Normal RNA Filename"]];
                case_.tumor_rna_filename = fields[fieldMappings["Tumor RNA Filename"]];
                case_.vcf_filename = fields[fieldMappings["VCF Filename"]];
                case_.methylation_filename = fields[fieldMappings["Methylation Filename"]];
                //case_.tumor_dna_allcount_filename = fields[fieldMappings["Tumor DNA Allcount Filename"]];
                case_.tumor_rna_allcount_filename = fields[fieldMappings["Tumor RNA Allcount Filename"]];
                case_.maf_filename = fields[fieldMappings["MAF Filename"]];
                case_.regional_expression_filename = fields[fieldMappings["Regional Expression Filename"]];
                case_.gene_expression_filename = fields[fieldMappings["Gene Expression Filename"]];
                case_.selected_variants_filename = fields[fieldMappings["Selected Variants Filename"]];
                case_.dna_reads_at_selected_variants_filename = fields[fieldMappings["DNA Reads At Selected Variants Filename"]];
                case_.rna_reads_at_selected_variants_filename = fields[fieldMappings["RNA Reads At Selected Variants Filename"]];
                case_.annotated_selected_variants_filename = fields[fieldMappings["Annotated Selected Variants Filename"]];
                case_.allele_specific_gene_expression_filename = fields[fieldMappings["Allele Specific Gene Expression Filename"]];
                case_.tumor_dna_gene_coverage_filname = fields[fieldMappings["Tumor DNA Gene Coverage Filename"]];
                case_.normal_rna_file_bam_md5 = fields[fieldMappings["Normal RNA BAM MD5"]];
                case_.normal_rna_file_bai_md5 = fields[fieldMappings["Normal RNA BAI MD5"]];
                case_.tumor_rna_file_bam_md5 = fields[fieldMappings["Tumor RNA BAM MD5"]];
                case_.tumor_rna_file_bai_md5 = fields[fieldMappings["Tumor RNA BAI MD5"]];
                case_.normal_dna_file_bam_md5 = fields[fieldMappings["Normal DNA BAM MD5"]];
                case_.normal_dna_file_bai_md5 = fields[fieldMappings["Normal DNA BAI MD5"]];
                case_.methylation_file_md5 = fields[fieldMappings["Methylation MD5"]];
 
                return case_;
            } // fromSaveFile

            public string GenerateLine()
            {
                return case_id + "\t" + normal_dna_file_id + "\t" + tumor_dna_file_id + "\t" + normal_rna_file_id + "\t" + tumor_rna_file_id + "\t" + maf_file_id + "\t" + methylation_file_id + "\t" + project_id + "\t" +
                    normal_dna_filename + "\t" + tumor_dna_filename + "\t" + normal_rna_filename + "\t" + tumor_rna_filename  + "\t" + methylation_filename + /*"\t" + tumor_dna_allcount_filename +*/ "\t" + tumor_rna_allcount_filename + "\t" +
                    maf_filename + "\t" + regional_expression_filename + "\t" + gene_expression_filename + "\t" +
                    selected_variants_filename + "\t" + dna_reads_at_selected_variants_filename + "\t" + rna_reads_at_selected_variants_filename + "\t" + annotated_selected_variants_filename + "\t" +
                    allele_specific_gene_expression_filename + "\t" + tumor_dna_gene_coverage_filname + "\t" + "\t" + vcf_filename +
                    normal_rna_file_bam_md5 + "\t" + normal_rna_file_bai_md5 + "\t" + tumor_rna_file_bam_md5 + "\t" + tumor_rna_file_bai_md5 + "\t" + normal_dna_file_bam_md5 + "\t" + normal_dna_file_bai_md5 + "\t" +  methylation_file_md5;
            } // GenerateLine

            public static void SaveToFile(Dictionary<string, Case> cases, string filename)
            {
                var file = CreateStreamWriterWithRetry(filename);
                file.WriteLine(
                    "Case ID" + "\t" +
                    "Normal DNA File ID" + "\t" +
                    "Tumor DNA File ID" + "\t" +
                    "Normal RNA File ID" + "\t" +
                    "Tumor RNA File ID" + "\t" +
                    "MAF File ID" + "\t" +
                    "Methylation File ID" + "\t" +
                    "Project ID" + "\t" +
                    "Normal DNA Filename" + "\t" +
                    "Tumor DNA Filename" + "\t" +
                    "Normal RNA Filename" + "\t" +
                    "Tumor RNA Filename" + "\t" +
                    "Methylation Filename" + "\t" +
                    //"Tumor DNA Allcount Filename" + "\t" +
                    "Tumor RNA Allcount Filename" + "\t" +
                    "MAF Filename" + "\t" +
                    "Regional Expression Filename" + "\t" +
                    "Gene Expression Filename" + "\t" +
                    "Selected Variants Filename" + "\t" +
                    "DNA Reads At Selected Variants Filename" + "\t" +
                    "RNA Reads At Selected Variants Filename" + "\t" +
                    "Annotated Selected Variants Filename" + "\t" +
                    "Allele Specific Gene Expression Filename" + "\t" +
                    "Tumor DNA Gene Coverage Filename" + "\t" +
                    "VCF Filename" + "\t" +
                    "Normal RNA BAM MD5" + "\t" +
                    "Normal RNA BAI MD5" + "\t" +
                    "Tumor RNA BAM MD5" + "\t" +
                    "Tumor RNA BAI MD5" + "\t" +
                    "Normal DNA BAM MD5" + "\t" +
                    "Normal DNA BAI MD5" + "\t" +
                    "Methylation MD5"
                    );

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

                if (downloadedFiles.ContainsKey(tumor_dna_file_id))
                {
                  tumor_dna_filename = downloadedFiles[tumor_dna_file_id].fileInfo.FullName;
                }

                if (normal_rna_file_id != "" && downloadedFiles.ContainsKey(normal_rna_file_id))
                {
                  normal_rna_filename = downloadedFiles[normal_rna_file_id].fileInfo.FullName;
                }

                if (downloadedFiles.ContainsKey(tumor_rna_file_id))
                {
                  tumor_rna_filename = downloadedFiles[tumor_rna_file_id].fileInfo.FullName;
                }

                if (downloadedFiles.ContainsKey(maf_file_id))
                {
                  maf_filename  = downloadedFiles[maf_file_id].fileInfo.FullName;
                }

                if (methylation_filename != "" && downloadedFiles.ContainsKey(methylation_file_id))
                {
                  methylation_filename  = downloadedFiles[methylation_file_id].fileInfo.FullName;
                }

                if (!derivedFiles.ContainsKey(case_id))
                {
                    return;
                }

                Dictionary<DerivedFile.Type, List<DerivedFile>> derivedFilesByType = new Dictionary<DerivedFile.Type,List<DerivedFile>>();
                foreach (var type in (DerivedFile.Type[])Enum.GetValues(typeof(DerivedFile.Type))) {
                    derivedFilesByType.Add(type, derivedFiles[case_id].Where(x => x.type == type).ToList());
                    if (derivedFilesByType[type].Count() > 1) {
                        Console.Write("Found more than one file of type " + type + " for case " + case_id + ":");
                        foreach (var file in derivedFilesByType[type]) {
                            Console.Write(" " + file.fileinfo.FullName);
                        }
                        Console.WriteLine();
                    }
                }

                if (derivedFilesByType[DerivedFile.Type.Allcount].Count() > 0)
                {
                    tumor_rna_allcount_filename = derivedFilesByType[DerivedFile.Type.Allcount][0].fileinfo.FullName;
                }

                if (derivedFilesByType[DerivedFile.Type.AlleleSpecificGeneExpression].Count() > 0)
                {
                    allele_specific_gene_expression_filename = derivedFilesByType[DerivedFile.Type.AlleleSpecificGeneExpression][0].fileinfo.FullName;
                }

                if (derivedFilesByType[DerivedFile.Type.AnnotatedSelectedVariants].Count() > 0)
                {
                    annotated_selected_variants_filename = derivedFilesByType[DerivedFile.Type.AnnotatedSelectedVariants][0].fileinfo.FullName;
                }

                if (derivedFilesByType[DerivedFile.Type.DNAReadsAtSelectedVariants].Count() > 0)
                {
                    dna_reads_at_selected_variants_filename = derivedFilesByType[DerivedFile.Type.DNAReadsAtSelectedVariants][0].fileinfo.FullName;
                }

                if (derivedFilesByType[DerivedFile.Type.GeneExpression].Count() > 0)
                {
                    gene_expression_filename = derivedFilesByType[DerivedFile.Type.GeneExpression][0].fileinfo.FullName;
                }

                if (derivedFilesByType[DerivedFile.Type.DNAReadsAtSelectedVariants].Count() > 0)
                {
                    dna_reads_at_selected_variants_filename = derivedFilesByType[DerivedFile.Type.DNAReadsAtSelectedVariants][0].fileinfo.FullName;
                }

                if (derivedFilesByType[DerivedFile.Type.RegionalExpression].Count() > 0)
                {
                    regional_expression_filename = derivedFilesByType[DerivedFile.Type.RegionalExpression][0].fileinfo.FullName;
                }

                if (derivedFilesByType[DerivedFile.Type.RNAReadsAtSelectedVariants].Count() > 0)
                {
                    rna_reads_at_selected_variants_filename = derivedFilesByType[DerivedFile.Type.RNAReadsAtSelectedVariants][0].fileinfo.FullName;
                }

                if (derivedFilesByType[DerivedFile.Type.SelectedVariants].Count() > 0)
                {
                    selected_variants_filename = derivedFilesByType[DerivedFile.Type.SelectedVariants][0].fileinfo.FullName;
                }

            } // loadFileLocations
        } // Case

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

        public class ASEConfirguation
        {
            public const string defaultBaseDirectory = @"\\msr-genomics-0\d$\gdc\";
            public const string defaultConfigurationFilePathame = defaultBaseDirectory + "configuration.txt";

            public string accessTokenPathname = defaultBaseDirectory +  @"access_token.txt";
            public List<string> dataDirectories = new List<string>();
            public string mafManifestPathname = defaultBaseDirectory + "mafManifest.txt";
            public string mutationCaller = "mutect";
            public List<string> programNames = new List<string>();
            public string binaryDirectory = defaultBaseDirectory + @"bin\";
            public string configuationFilePathname = defaultConfigurationFilePathame;
            public string casesFilePathname = defaultBaseDirectory + "cases.txt";
            public string indexDirectory = @"d:\gdc\indices\hg38-20";
            public string derivedFilesDirectory = "derived_files";    // This is relative to each download directory
            public string hpcScriptFilename = "";    // The empty string says to black hole this script
            public string hpcScheduler = "gcr";
            public string hpcBinariesDirectory = @"\\gcr\scratch\b99\bolosky\";
            public string hpcIndexDirectory = @"\\msr-genomics-0\d$\gdc\indices\hg38-20";

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

                for (int i = 0; i < args.Count(); i++) {
                    if (args[i] == "-configuration") {
                        if (i >= args.Count() - 1) {
                            Console.WriteLine("-configuation can't be the last parameter, it needs to be followed by a file path.  Ignoring.");
                        } else {
                            pathname = args[i+1];
                        }
                    }
                }

                var retVal = new ASEConfirguation();

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
                        retVal.binaryDirectory = fields[1];
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

        static readonly int FileIdLength = "01c5f902-ec3c-4eb7-9e38-1d29ae6ab959".Count();

        public class DownloadedFile
        {
            public readonly string file_id;
            public readonly FileInfo fileInfo;
            public readonly string storedMD5;    // This is the null string for files for which we don't know the sum (i.e., they haven't yet been computed or we don't know the right answer).

            public DownloadedFile(string file_id_, string pathname, string storedMd5Sum_)
            {
                file_id = file_id_;
                storedMD5 = storedMd5Sum_;

                fileInfo = new FileInfo(pathname);
            }
        } // DownloadedFile


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

                var extension = filename.Substring(FileIdLength);
 
                if (extension == ".allcount.gz") {
                    type = Type.Allcount;
                }
                else if (extension == ".selectedVariants")
                {
                    type = Type.SelectedVariants;
                } else if (extension == ".regional_expression.txt") {
                    type = Type.RegionalExpression;
                } else if (extension == ".gene_expression.txt") {
                    type = Type.GeneExpression;
                } else if (extension == ".rna-reads-at-selected-variants.txt") {
                    type = Type.RNAReadsAtSelectedVariants;
                }
                else if (extension == ".dna-reads-at-selected-variants.txt")
                {
                    type = Type.DNAReadsAtSelectedVariants;
                }
                else if (extension == ".vcf")
                {
                    type = Type.VCF;
                }
                else
                {
                    Console.WriteLine("Derived file with unknown extension: " + pathname);
                }
            } // DerviedFile.DerivedFile


            public enum Type { Unknown, Allcount, RegionalExpression, GeneExpression, SelectedVariants, DNAReadsAtSelectedVariants, RNAReadsAtSelectedVariants, AnnotatedSelectedVariants, AlleleSpecificGeneExpression, VCF };
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

        static void ScanOneFilesystem(ASEConfirguation configuration, string downloadedFilesDirectory, ScanFilesystemState state) 
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
                        Console.WriteLine("Found subdirectory a downloaded files dierectory whose name doesn't appear to be a guid, ignoring: " + subdir);
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
                    var downloadedFile = new DownloadedFile(file_id, candidatePathname, md5Value);
                    downloadedFiles.Add(downloadedFile);

                    totalBytesInDownloadedFiles += (ulong)downloadedFile.fileInfo.Length;
                }
            } // foreach subdir

            string derivedFilesDirectory = downloadedFilesDirectory + @"..\" + configuration.derivedFilesDirectory;
            foreach (var derivedCase in Directory.EnumerateDirectories(derivedFilesDirectory))
            {
                var caseId = ASETools.GetFileNameFromPathname(derivedCase).ToLower();
 
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
                Console.WriteLine("Data directory " + downloadedFilesDirectory + " has " + SizeToUnits(freeBytesAvailable) + "B free, with " + nDownloadedFiles + " (" + SizeToUnits(totalBytesInDownloadedFiles) + "B) downloaded files and " + nDerivedFiles + " (" +
                    SizeToUnits(totalBytesInDerivedFiles) + "B) derived files.");

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

            var state = new ScanFilesystemState();
            var threads = new List<Thread>();

            foreach (var directory in configuration.dataDirectories)
            {
                threads.Add(new Thread(() => ScanOneFilesystem(configuration, directory, state)));
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

            public bool ParseFile(Parse parser, out List<outputType> result)
            {
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
                    Console.WriteLine("Headerized file: missing some columns.");
                    result = null;
                    return false;
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

                return true;
            } // ParseFile

            StreamReader inputFile;
            bool hasVersion;
            bool hasDone;
            string expectedVersion;
            List<string> wantedFields;
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
            public readonly string file_id; // Of the MAF file

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
             string file_id_)
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
                file_id = file_id_;
            }

            static MAFLine ParseLine(Dictionary<string, int> fieldMappings, string[] fields, string file_id)
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
                    file_id
                    );
            }

            static public List<MAFLine> ReadFile(string filename, string file_id)
            {
                var compressedStreamReader = CreateStreamReaderWithRetry(filename);

                var inputFile = new StreamReader(new GZipStream(compressedStreamReader.BaseStream, CompressionMode.Decompress));

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

                var headerizedFile = new HeaderizedFile<MAFLine>(inputFile, true, false, "#version 2.4", neededFields);

                List<MAFLine> result;

                if (!headerizedFile.ParseFile((a, b) => ParseLine(a, b, file_id), out result))
                {
                    Console.WriteLine("Error reading MAF File " + filename);
                    return null;
                }

                inputFile.Close();
                compressedStreamReader.Close();

                return result;
            } // ReadFile

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


    } // ASETools
}
