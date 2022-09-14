using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using ASELib;
using System.Threading;
using System.Diagnostics;
using System.IO;

namespace CheckDone
{
    class Program
    {

        static Dictionary<string, List<string>> FilenamesByDataDirectory = new Dictionary<string, List<string>>();
        static ASETools.Configuration configuration;
        static int totalFiles = 0;
        static int nFilesProcessed = 0;

        static void ProcessOneDataDirectory(List<string> filenames)
        {
            const int nBytesToCheck = 10; // **done**\n\r
            const int nBytesToCheckIndex = 12; // **done**\t\t\n\r
            var lastBitOfFile = new byte[nBytesToCheck];
            var lastBitOfFileIndex = new byte[nBytesToCheckIndex];

            foreach (var filename in filenames)
            {
                lock (FilenamesByDataDirectory)
                {
                    nFilesProcessed++;
                    if (nFilesProcessed % nPerDot == 0) Console.Write(".");
                }

                if (filename.ToLower().EndsWith(".gz"))
                {
                    //
                    // We can't seek in gzipped files, so just read the whole thing.
                    //
                    var reader = ASETools.CreateCompressedStreamReaderWithRetry(filename);

                    string line;
                    bool seenDone = false;
                    while (null != (line = reader.ReadLine()))
                    {
                        if (seenDone)
                        {
                            Console.WriteLine(filename + " continues after **done**");
                            break;
                        }

                        if (line == "**done**")
                        {
                            seenDone = true;
                        }
                    } // while we have an input line

                    if (!seenDone)
                    {
                        Console.WriteLine(filename + " is truncated.");
                    }
                } else
                {
                    //
                    // Not compressed, so just seek to the end of the file.  This skips checking for multiple **done** lines, but
                    // that's not really a problem I can forsee happening.
                    //
                    FileStream filestream = null;
                    int retryCount = 0;
                    while (true)
                    {
                        try
                        {
                            filestream = new FileStream(filename, FileMode.Open);
                            break;
                        }
                        catch (Exception e)
                        {
                            if (e is IOException && (e.HResult & 0xffff) == 32 /* sharing violation */ && retryCount < 3)
                            {
                                retryCount++;
                                Thread.Sleep(1000);
                            }
                            else
                            {
                                Console.WriteLine("Exception opening " + filename + ", ignoring.  " + e.GetType() + ": " + e.Message + " retry count " + retryCount + ", HRESULT " + (e.HResult & 0xffff));
                                break;
                            }
                        }
                    }

                    if (filestream == null)
                    {
                        continue;
                    }

                    var index = filename.EndsWith(".index");

                    int bytesToCheckThisFile = index ? nBytesToCheckIndex : nBytesToCheck;
                    byte[] buffer = index ? lastBitOfFileIndex : lastBitOfFile;

                    if (filestream.Length < bytesToCheckThisFile)
                    {
                        Console.WriteLine(filename + " is truncated.");
                        continue;
                    } else if (filestream.Length == bytesToCheckThisFile)
                    {
                        Console.WriteLine(filename + " appears to be nothing but the **done** line.");
                    }

                    filestream.Position = filestream.Length - bytesToCheckThisFile;
                    if (bytesToCheckThisFile != filestream.Read(buffer, 0, bytesToCheckThisFile))
                    {
                        Console.WriteLine("Error reading " + filename);
                        continue;
                    }

                    if ((!index && (
                            buffer[0] != '*' ||
                            buffer[1] != '*' ||
                            buffer[2] != 'd' ||
                            buffer[3] != 'o' ||
                            buffer[4] != 'n' ||
                            buffer[5] != 'e' ||
                            buffer[6] != '*' ||
                            buffer[7] != '*' ||
                            buffer[8] != '\r' ||
                            buffer[9] != '\n'))     ||
                        (index && 
                            (buffer[0] != '*' ||    // **done**\t\t\r\n
                            buffer[1] != '*' ||
                            buffer[2] != 'd' ||
                            buffer[3] != 'o' ||
                            buffer[4] != 'n' ||
                            buffer[5] != 'e' ||
                            buffer[6] != '*' ||
                            buffer[7] != '*' ||
                            buffer[8] != '\t' ||
                            buffer[9] != '\t' ||
                            buffer[10] != '\r' ||
                            buffer[11] != '\n') &&

                            (buffer[2] != '*' || // **done**\r\n
                            buffer[3] != '*' ||
                            buffer[4] != 'd' ||
                            buffer[5] != 'o' ||
                            buffer[6] != 'n' ||
                            buffer[7] != 'e' ||
                            buffer[8] != '*' ||
                            buffer[9] != '*' ||
                            buffer[10] != '\r' ||
                            buffer[11] != '\n')
                        )
                    )
                    {
                        Console.WriteLine(filename + " is truncated.");
                    }
                } // If it's not gzipped
            } // foreach filename
        } // ProcessOneDataDirectory


        static int nPerDot;
        static void Main(string[] args)
        {
            var timer = new Stopwatch();
            timer.Start();

            configuration = ASETools.Configuration.loadFromFile(args);

            bool skipAllcount = configuration.commandLineArgs.Contains("-a");

            string[] fileTypesToUse = null;

            if (configuration.commandLineArgs.Any(x => x != "-a"))
            {
                fileTypesToUse = configuration.commandLineArgs.Where(x => x != "-a").Select(x => x.ToLower()).ToArray();
            }

            var cases = ASETools.Case.LoadCases(configuration.casesFilePathname);

            if (null == cases)
            {
                Console.WriteLine("Unable to load cases.  You must generate it before running this tool.");
                return;
            }

            foreach (var caseEntry in cases)
            {
                var case_ = caseEntry.Value;

                if (!skipAllcount)
                {
                    HandleFilename(case_.normal_dna_allcount_filename, "NormalDNAAllcount", fileTypesToUse);
                    HandleFilename(case_.normal_rna_allcount_filename, "NormalRNAAllcount", fileTypesToUse);
                    HandleFilename(case_.tumor_dna_allcount_filename, "TumorDNAAllcount", fileTypesToUse);
                    HandleFilename(case_.tumor_rna_allcount_filename, "TumorRNAAllcount", fileTypesToUse);
                }
                HandleFilename(case_.regional_expression_filename, "RegionalExpression", fileTypesToUse);
                HandleFilename(case_.gene_expression_filename, "GeneExpression", fileTypesToUse);
                HandleFilename(case_.tentative_selected_variants_filename, "TentativeSelectedVariants", fileTypesToUse);
                HandleFilename(case_.normal_dna_reads_at_tentative_selected_variants_index_filename, "NormalDNAReadsAtSelectedVariants", fileTypesToUse);
                HandleFilename(case_.normal_rna_reads_at_tentative_selected_variants_index_filename, "NormalRNAReadsAtSelectedVariants", fileTypesToUse);
                HandleFilename(case_.tumor_dna_reads_at_tentative_selected_variants_index_filename, "TumorDNAReadsAtSelectedVariants", fileTypesToUse);
                HandleFilename(case_.tumor_rna_reads_at_tentative_selected_variants_index_filename, "TumorRNAReadsAtSelectedVariants", fileTypesToUse);
                HandleFilename(case_.tentative_annotated_selected_variants_filename, "TentativeAnnotatedSelectedVariants", fileTypesToUse);
                HandleFilename(case_.annotated_selected_variants_filename, "AnnotatedSelectedVariants", fileTypesToUse);
                HandleFilename(case_.normal_allele_specific_gene_expression_filename, "NormalAlleleSpecificExpression", fileTypesToUse);
                HandleFilename(case_.tumor_allele_specific_gene_expression_filename, "TumorAlleleSpecificExpression", fileTypesToUse);
                HandleFilename(case_.tumor_dna_gene_coverage_filname, "TumorDNAGeneCoverage", fileTypesToUse);
                HandleFilename(case_.extracted_maf_lines_filename, "ExtractedMAFLines", fileTypesToUse);
                HandleFilename(case_.all_maf_lines_filename, "AllMAFLines", fileTypesToUse);
                HandleFilename(case_.normal_dna_mapped_base_count_filename, "NormalDNAMappedBaseCount", fileTypesToUse);
                HandleFilename(case_.tumor_dna_mapped_base_count_filename, "TumorDNAMappedBaseCount", fileTypesToUse);
                HandleFilename(case_.normal_rna_mapped_base_count_filename, "NormalRNAMappedBaseCount", fileTypesToUse);
                HandleFilename(case_.tumor_rna_mapped_base_count_filename, "TumorRNAMappedBaseCount", fileTypesToUse);
                HandleFilename(case_.selected_variant_counts_by_gene_filename, "SelectedVariantCountsByGene", fileTypesToUse);
                HandleFilename(case_.selected_regulatory_maf_filename, "SelectedRegulatoryMAF", fileTypesToUse);
                HandleFilename(case_.annotated_regulatory_regions_filename, "AnnotatedRegulatoryRegions", fileTypesToUse);
                HandleFilename(case_.regulatory_mutations_near_mutations_filename, "RegulatoryMutationsNearMutations", fileTypesToUse);
                HandleFilename(case_.annotated_geneHancer_lines_filename, "AnnotatedGeneHancerLines", fileTypesToUse);
                HandleFilename(case_.expression_by_gene_filename, "ExpressionByGene", fileTypesToUse);
                HandleFilename(case_.isoform_read_counts_filename, "IsoformReadCounts", fileTypesToUse);
                HandleFilename(case_.case_metadata_filename, "CaseMetadata", fileTypesToUse);
                HandleFilename(case_.tentative_asv_without_cnvs_filename, "TentativeASVsWithoutCNVs", fileTypesToUse);
                HandleFilename(case_.variant_phasing_filename, "VariantPhasing", fileTypesToUse);
                HandleFilename(case_.vcf_statistics_filename, "VCFStatistics", fileTypesToUse);
                HandleFilename(case_.read_statictics_filename, "ReadStatictics", fileTypesToUse);
                HandleFilename(case_.gene_expression_fraction_filename, "GeneExpressionFraction", fileTypesToUse);

                foreach (var tumor in ASETools.BothBools)
                {
                    foreach (var variantCaller in ASETools.EnumUtil.GetValues<ASETools.VariantCaller>())
                    {
                        HandleFilename(case_.perVariantCaller[variantCaller][tumor].venn_filename, "Venn", fileTypesToUse);
                    }
                }
            }

            HandleFilenameIfExists(configuration.finalResultsDirectory + ASETools.SingleReadPhasingFilename, "SingleReadPhasing", fileTypesToUse);

            var expressionByChromosomeMap = ASETools.ExpressionDistributionByChromosomeMap.LoadFromFile(configuration.expression_distribution_by_chromosome_map_filename);
            foreach (var file in expressionByChromosomeMap.allFiles())
            {
                HandleFilename(file, "ExpressionDistribution", fileTypesToUse);
            }

            foreach (var disease in ASETools.GetListOfDiseases(cases))
            {
                // these don't have **done** for some reason.  They should! HandleFilenameIfExists(configuration.expressionFilesDirectory + ASETools.Expression_filename_base + disease, "Expression", fileTypesToUse);

                foreach (var chromosome in ASETools.chromosomes)
                {
                    HandleFilenameIfExists(configuration.geneScatterGraphsLinesWithPercentilesDirectory + ASETools.GeneScatterGraphLinesWithPercentilesPrefix + disease + "_" + ASETools.chromosomeNameToNonChrForm(chromosome), "ScatterGraphLinesWithPercentiles", fileTypesToUse);
                }
            }

            ASETools.PrintMessageAndNumberBar("Processing", "files in " + FilenamesByDataDirectory.Count() + " data directories", totalFiles, out nPerDot);

            //
            // Run one thread per data directory, since this is likely IO bound by the server disks rather than
            // by the local processors.
            //
            var threads = new List<Thread>();
            foreach (var dataDirectoryEntry in FilenamesByDataDirectory)
            {
                threads.Add(new Thread(() => ProcessOneDataDirectory(dataDirectoryEntry.Value)));
            }

            threads.ForEach(t => t.Start());
            threads.ForEach(t => t.Join());

            Console.WriteLine();
            Console.WriteLine("Processed " + totalFiles + " files in " + ASETools.ElapsedTimeInSeconds(timer));
        } // Main

        static void HandleFilenameIfExists(string filename, string fileType, string [] fileTypesToUse)
        {
            if (File.Exists(filename))
            {
                HandleFilename(filename, fileType, fileTypesToUse);
            }
        }

        static void HandleFilename(string filename, string fileType, string [] fileTypesToUse)
        {
            if (filename == "" || fileTypesToUse != null && !fileTypesToUse.Contains(fileType.ToLower()))
            {
                return;
            }

            var dataDirectory = ASETools.GetDataDirectoryFromFilename(filename, configuration);

            if (!FilenamesByDataDirectory.ContainsKey(dataDirectory))
            {
                FilenamesByDataDirectory.Add(dataDirectory, new List<string>());
            }

            FilenamesByDataDirectory[dataDirectory].Add(filename);
            totalFiles++;
        }
    }
}
