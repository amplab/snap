using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using ASELib;
using System.IO;
using System.Diagnostics;

namespace MakeCaseMetadata
{
    class Program
    {
        static ASETools.CommonData commonData;

        static void Main(string[] args)
        {
            commonData = ASETools.CommonData.LoadCommonData(args);
            if (null == commonData)
            {
                return;
            }

            var caseIdsToRun = commonData.configuration.commandLineArgs.ToList();
            var casesToRun = commonData.listOfCases.Where(_ => caseIdsToRun.Contains(_.case_id)).ToList();

            var casesWithoutEnoughData = casesToRun.Where(_ => _.normal_dna_filename == "" || _.tumor_dna_filename == "" || _.tumor_rna_filename == "" || _.normal_dna_reads_at_tentative_selected_variants_filename == "" ||
                                                _.tumor_dna_reads_at_tentative_selected_variants_filename == "" || _.tumor_rna_reads_at_tentative_selected_variants_filename == "" ||
                                                (_.normal_rna_file_id != "" && (_.normal_rna_filename == "" || _.normal_rna_reads_at_tentative_selected_variants_filename == ""))).ToList();

            if (casesWithoutEnoughData.Count() != 0)
            {
                Console.Write("Missing some input data for cases:");
                casesWithoutEnoughData.ForEach(_ => Console.Write(" " + _.case_id));
                Console.WriteLine();
                Console.WriteLine("Skipping them.");

                casesToRun = casesToRun.Where(_ => !casesWithoutEnoughData.Contains(_)).ToList();
            }
            
            int nPerDot;
            ASETools.PrintMessageAndNumberBar("Processing", "cases", casesToRun.Count(), out nPerDot);
            var threading = new ASETools.WorkerThreadHelper<ASETools.Case, int>(casesToRun, HandleOneCase, null, null, nPerDot);
            threading.run();


            Console.WriteLine(ASETools.ElapsedTimeInSeconds(commonData.timer));
        } // Main

        static void HandleOneCase(ASETools.Case case_, int state)
        {
            var bamMetadata = new Dictionary<bool, Dictionary<bool, ASETools.BAMMetadata>>();
            var sample = new Dictionary<bool, Dictionary<bool, string>>();
            foreach (var tumor in ASETools.BothBools)
            {
                bamMetadata.Add(tumor, new Dictionary<bool, ASETools.BAMMetadata>());
                sample.Add(tumor, new Dictionary<bool, string>());

                foreach (var dna in ASETools.BothBools)
                {
                    bamMetadata[tumor].Add(dna, GetBamMetadataFromFile(case_.getReadsAtTentativeSelectedVariantsFilename(tumor, dna)));

                    sample[tumor].Add(dna, GetSampleFromFile(case_.getDownloadedReadsFilename(tumor, dna)));
                } // dna
            } // tumor

            var outputFilename = ASETools.GetDirectoryFromPathname(case_.normal_dna_reads_at_tentative_selected_variants_filename) + @"\" + case_.case_id + ASETools.caseMetadataExtension;
            var outputFile = ASETools.CreateStreamWriterWithRetry(outputFilename);
            if (null == outputFile)
            {
                Console.WriteLine("Unable to open output file " + outputFilename);
                return;
            }

            outputFile.Write("Case ID");
            foreach (var tumor in ASETools.BothBools)
            {
                foreach (var dna in ASETools.BothBools)
                {
                    var specifier = ASETools.CaseMetadata.getSpecifier(tumor, dna);
                    outputFile.Write("\t" + specifier + "is paired");
                    outputFile.Write("\t" + specifier + "min read length");
                    outputFile.Write("\t" + specifier + "max read length");
                    outputFile.Write("\t" + specifier + "mean read length");
                    outputFile.Write("\t" + specifier + "median read length");
                    outputFile.Write("\t" + specifier + "min insert");
                    outputFile.Write("\t" + specifier + "max insert");
                    outputFile.Write("\t" + specifier + "mean insert");
                    outputFile.Write("\t" + specifier + "median insert");
                    outputFile.Write("\t" + specifier + "Sample");
                }
            }
            outputFile.WriteLine();

            outputFile.Write(case_.case_id);
            foreach (var tumor in ASETools.BothBools)
            {
                foreach (var dna in ASETools.BothBools)
                {
                    outputFile.Write("\t" + bamMetadata[tumor][dna].isPaired);
                    outputFile.Write("\t" + bamMetadata[tumor][dna].minReadLength);
                    outputFile.Write("\t" + bamMetadata[tumor][dna].maxReadLength);
                    outputFile.Write("\t" + bamMetadata[tumor][dna].meanReadLength);
                    outputFile.Write("\t" + bamMetadata[tumor][dna].medianReadLength);
                    outputFile.Write("\t" + bamMetadata[tumor][dna].minInsert);
                    outputFile.Write("\t" + bamMetadata[tumor][dna].maxInsert);
                    outputFile.Write("\t" + bamMetadata[tumor][dna].meanInsert);
                    outputFile.Write("\t" + bamMetadata[tumor][dna].medianInsert);
                    outputFile.Write("\t" + sample[tumor][dna]);
                } // dna
            } // tumor

            outputFile.WriteLine();
            outputFile.WriteLine("**done**");
            outputFile.Close();
        } // HandleOneCase

        static ASETools.BAMMetadata GetBamMetadataFromFile(string filename)
        {
            var retVal = new ASETools.BAMMetadata(false, int.MaxValue, -1, double.NegativeInfinity, double.NegativeInfinity, int.MaxValue, -1, double.NegativeInfinity, double.NegativeInfinity);

            if (filename == "")
            {
                return retVal;
            }

            var reader = new ASETools.ConsolodatedFileReader();

            if (!reader.open(filename))
            {
                throw new Exception("Unable to open reads file (or its index): " + filename);
            }

            var readLengths = new List<int>(); // we keep them all, so that we can find the median to return.  
            var insertLengths = new List<int>();    

            var allSubfiles = reader.allSubfiles();

            foreach (var chromosome in ASETools.chromosomes)
            {
                var subfilesToRead = allSubfiles.Where(_ => _.Contains(chromosome)).Take(10).ToList();

                foreach (var subfileName in subfilesToRead)
                {
                    var subfile = reader.getSubfile(subfileName);

                    var samLines = ASETools.SAMLine.ReadFromFile(subfile);
                    subfile.Close();

                    foreach (var samLine in samLines)
                    {
                        retVal.isPaired |= samLine.isPaired();
                        int readLen = samLine.seq.Length;
                        readLengths.Add(readLen);
                        retVal.minReadLength = Math.Min(retVal.minReadLength, readLen);
                        retVal.maxReadLength = Math.Max(retVal.maxReadLength, readLen);
                        if (samLine.isPaired() && !samLine.isUnmapped() && !samLine.isNextUnmapped() && samLine.pos != 0 && samLine.pnext != 0 && (samLine.rnext == "=" || samLine.rnext == samLine.rname))
                        {
                            insertLengths.Add(Math.Abs(samLine.pos - samLine.pnext));
                        }
                    } // sam lines
                } // subfiles
            } // chromosome

            var nInsertLengths = insertLengths.Count();

            if (nInsertLengths >= 10)
            {
                insertLengths.Sort();
                retVal.minInsert = insertLengths[0];
                retVal.maxInsert = insertLengths[nInsertLengths - 1];

                retVal.medianInsert = insertLengths.Median();
                retVal.meanInsert = (double)insertLengths.Select(_ => (long)_).Sum() / insertLengths.Count();
            }
            else
            {
                retVal.isPaired = false;    // Not paired enough to count
                retVal.medianInsert = double.NegativeInfinity;
                retVal.meanInsert = double.NegativeInfinity;
            }

            var nReadLengths = readLengths.Count();
            if (nReadLengths > 10)
            {
                retVal.medianReadLength = readLengths.Median();
                retVal.meanReadLength = readLengths.Select(_ => (double)_).Sum() / nReadLengths;
            }

            return retVal;
        } // GetBamMetadataFromFile

        static string GetSampleFromFile(string filename)
        {
            if (filename == "")
            {
                return "";
            }

            var rgLines = ASETools.RunProcessAndGetOutput(commonData.configuration.binariesDirectory + "samtools.exe", "view -H " + filename).Where(_ => _.StartsWith("@RG"));

            string sample = "";

            foreach (var rgLine in rgLines)
            {
                var fields = rgLine.Split('\t');

                foreach (var smField in fields.Where(_ => _.StartsWith("SM:")))
                {
                    if (sample != "" && sample != smField.Substring(3))
                    {
                        Console.WriteLine(filename + " has more than one sample: " + sample + " and " + smField.Substring(3));
                    }
                    sample = smField.Substring(3);
                } // smField
            } // rgLine

            if (sample == "")
            {
                Console.WriteLine(filename + " has no sample");
            }

            return sample;
        } // GetSampleFromFile
    }
}
