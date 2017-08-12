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
                    if (nFilesProcessed % 100 == 0) Console.Write(".");
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
                    FileStream filestream;
                    try
                    {
                        filestream = new FileStream(filename, FileMode.Open);
                    } catch
                    {
                        Console.WriteLine("Exception opening " + filename + ", ignoring.");
                        continue;
                    }

                    int bytesToCheckThisFile = nBytesToCheck;
                    byte[] buffer = lastBitOfFile;

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

                    if (buffer[0] != '*' ||
                        buffer[1] != '*' ||
                        buffer[2] != 'd' ||
                        buffer[3] != 'o' ||
                        buffer[4] != 'n' ||
                        buffer[5] != 'e' ||
                        buffer[6] != '*' ||
                        buffer[7] != '*' ||
                        buffer[8] != '\r' ||
                        buffer[9] != '\n')
                    {
                        Console.WriteLine(filename + " is truncated.");
                    }
                } // If it's not gzipped
            } // foreach filename
        } // ProcessOneDataDirectory

        static void Main(string[] args)
        {
            var timer = new Stopwatch();
            timer.Start();

            configuration = ASETools.Configuration.loadFromFile(args);

            if (configuration.commandLineArgs.Count() != 0 && (configuration.commandLineArgs.Count() != 1 || configuration.commandLineArgs[0] != "-a"))
            {
                Console.WriteLine("usage: CheckDone {-configuration configurationFileName} {-a}");
                Console.WriteLine("-a means to skip allcount files (which is where most of the time goes).");
                return;
            }

            bool skipAllcount = configuration.commandLineArgs.Contains("-a");

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
                    HandleFilename(case_.normal_dna_allcount_filename);
                    HandleFilename(case_.normal_rna_allcount_filename);
                    HandleFilename(case_.tumor_dna_allcount_filename);
                    HandleFilename(case_.tumor_rna_allcount_filename);
                }
                HandleFilename(case_.regional_expression_filename);
                HandleFilename(case_.gene_expression_filename);
                HandleFilename(case_.selected_variants_filename);
                HandleFilename(case_.normal_dna_reads_at_selected_variants_index_filename);
                HandleFilename(case_.normal_rna_reads_at_selected_variants_index_filename);
                HandleFilename(case_.tumor_dna_reads_at_selected_variants_index_filename);
                HandleFilename(case_.tumor_rna_reads_at_selected_variants_index_filename);
                HandleFilename(case_.annotated_selected_variants_filename);
                HandleFilename(case_.tumor_dna_gene_coverage_filname);
                HandleFilename(case_.extracted_maf_lines_filename);
                HandleFilename(case_.normal_dna_mapped_base_count_filename);
                HandleFilename(case_.tumor_dna_mapped_base_count_filename);
                HandleFilename(case_.normal_rna_mapped_base_count_filename);
                HandleFilename(case_.tumor_rna_mapped_base_count_filename);
                HandleFilename(case_.selected_variant_counts_by_gene_filename);
            }

            Console.Write("Processing " + totalFiles + " files in " + FilenamesByDataDirectory.Count() + " data directories (one dot/hundred): ");

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

            Console.WriteLine("Processed " + totalFiles + " files in " + ASETools.ElapsedTimeInSeconds(timer));
        } // Main

        static void HandleFilename(string filename)
        {
            if (filename == "")
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
