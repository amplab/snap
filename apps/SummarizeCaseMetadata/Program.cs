using System;
using System.CodeDom;
using System.Collections.Generic;
using System.ComponentModel;
using System.ComponentModel.Design.Serialization;
using System.Linq;
using System.Runtime.Versioning;
using System.Text;
using System.Threading.Tasks;
using ASELib;

namespace SummarizeCaseMetadata
{
    class Program
    {
        static void Main(string[] args)
        {
            var commonData = ASETools.CommonData.LoadCommonData(args);
            if (null == commonData)
            {
                return;
            }

            if (commonData.listOfCases.Any(_ => _.case_metadata_filename == ""))
            {
                Console.WriteLine("Not all cases have metadata.");
                return;
            }

            //
            // All of these map tumor->dna->value
            //
            var nPaired = new Dictionary<bool, Dictionary<bool, int>>(); // fraction is determined by dividing by cases, except for normal RNA, where you need to count
            var minReadLengthHistogram = new Dictionary<bool, Dictionary<bool, ASETools.PreBucketedHistogram>>();
            var maxReadLengthHistogram = new Dictionary<bool, Dictionary<bool, ASETools.PreBucketedHistogram>>();
            var meanReadLengthHistogram = new Dictionary<bool, Dictionary<bool, ASETools.PreBucketedHistogram>>();
            var medianReadLengthHistogram = new Dictionary<bool, Dictionary<bool, ASETools.PreBucketedHistogram>>();
            var minInsert = new Dictionary<bool, Dictionary<bool, ASETools.PreBucketedHistogram>>();
            var maxInsert = new Dictionary<bool, Dictionary<bool, ASETools.PreBucketedHistogram>>();
            var meanInsert = new Dictionary<bool, Dictionary<bool, ASETools.PreBucketedHistogram>>();
            var medianInsert = new Dictionary<bool, Dictionary<bool, ASETools.PreBucketedHistogram>>();

            foreach (var tumor in ASETools.BothBools)
            {
                nPaired.Add(tumor, new Dictionary<bool, int>());
                minReadLengthHistogram.Add(tumor, new Dictionary<bool, ASETools.PreBucketedHistogram>());
                maxReadLengthHistogram.Add(tumor, new Dictionary<bool, ASETools.PreBucketedHistogram>());
                meanReadLengthHistogram.Add(tumor, new Dictionary<bool, ASETools.PreBucketedHistogram>());
                medianReadLengthHistogram.Add(tumor, new Dictionary<bool, ASETools.PreBucketedHistogram>());
                minInsert.Add(tumor, new Dictionary<bool, ASETools.PreBucketedHistogram>());
                maxInsert.Add(tumor, new Dictionary<bool, ASETools.PreBucketedHistogram>());
                meanInsert.Add(tumor, new Dictionary<bool, ASETools.PreBucketedHistogram>());
                medianInsert.Add(tumor, new Dictionary<bool, ASETools.PreBucketedHistogram>());

                foreach (var dna in ASETools.BothBools)
                {
                    nPaired[tumor].Add(dna, 0);
                    minReadLengthHistogram[tumor].Add(dna, new ASETools.PreBucketedHistogram(0, 200, 1));
                    maxReadLengthHistogram[tumor].Add(dna, new ASETools.PreBucketedHistogram(0, 200, 1));
                    meanReadLengthHistogram[tumor].Add(dna, new ASETools.PreBucketedHistogram(0, 200, 1));
                    medianReadLengthHistogram[tumor].Add(dna, new ASETools.PreBucketedHistogram(0, 200, 1));
                    minInsert[tumor].Add(dna, new ASETools.PreBucketedHistogram(0, 2000, 10));
                    maxInsert[tumor].Add(dna, new ASETools.PreBucketedHistogram(0, 2000, 10));
                    meanInsert[tumor].Add(dna, new ASETools.PreBucketedHistogram(0, 2000, 10));
                    medianInsert[tumor].Add(dna, new ASETools.PreBucketedHistogram(0, 2000, 10));
                } // dna
            } // tumor

            int nPerDot;
            ASETools.PrintMessageAndNumberBar("Processing", "cases", commonData.listOfCases.Count(), out nPerDot);
            int nProcessed = 0;

            foreach (var case_ in commonData.listOfCases)
            {
                var metadata = ASETools.CaseMetadata.ReadFromFile(case_.case_metadata_filename);
                if (null == metadata)
                {
                    Console.WriteLine("Unable to load case metadata for " + case_.case_id);
                    return;
                }

                foreach (var tumor in ASETools.BothBools)
                {
                    foreach (var dna in ASETools.BothBools)
                    {
                        if (!tumor && !dna && case_.normal_rna_file_id == "")
                        {
                            continue;   // This one doesn't exist.
                        }

                        var bamData = metadata.getBAMMetadata(tumor, dna);

                        if (bamData.isPaired)
                        {
                            nPaired[tumor][dna]++;
                        }

                        minReadLengthHistogram[tumor][dna].addValue(bamData.minReadLength);
                        maxReadLengthHistogram[tumor] [dna].addValue(bamData.maxReadLength);
                        meanReadLengthHistogram[tumor] [dna].addValue(bamData.meanReadLength);
                        medianReadLengthHistogram[tumor] [dna].addValue(bamData.medianReadLength);

                        if (bamData.isPaired)
                        {
                            minInsert[tumor][dna].addValue(bamData.minInsert);
                            maxInsert[tumor][dna].addValue(bamData.maxInsert);
                            meanInsert[tumor][dna].addValue(bamData.meanInsert);
                            medianInsert[tumor][dna].addValue(bamData.medianInsert);
                        }
                    } // dna
                } // tumor

                nProcessed++;
                if (nProcessed % nPerDot == 0)
                {
                    Console.Write(".");
                }
            } // foreach case

            Console.WriteLine();

            var ouputFilename = commonData.configuration.finalResultsDirectory + ASETools.CaseMetadataSummaryFilename;
            var outputFile = ASETools.CreateStreamWriterWithRetry(ouputFilename);
            if (null == outputFile)
            {
                Console.WriteLine("Unable to open output file " + ouputFilename);
                return;
            }


            foreach (var tumor in ASETools.BothBools)
            {
                foreach (var dna in ASETools.BothBools)
                {
                    int nCases = (!tumor && !dna) ? commonData.listOfCases.Where(_ => _.normal_rna_file_id != "").Count() : commonData.listOfCases.Count();
                    outputFile.WriteLine(ASETools.tumorToString[tumor] + " " + ASETools.dnaToString[dna] + " nPaired " + nPaired[tumor][dna] + " of " + nCases + ", fraction " + (double)nPaired[tumor][dna] / nCases);
                    outputFile.WriteLine(ASETools.tumorToString[tumor] + " " + ASETools.dnaToString[dna] + " distribution of min read length");
                    minReadLengthHistogram[tumor][dna].WriteHistogram(outputFile);
                    outputFile.WriteLine(ASETools.tumorToString[tumor] + " " + ASETools.dnaToString[dna] + " distribution of max read length");
                    maxReadLengthHistogram[tumor] [dna].WriteHistogram(outputFile);
                    outputFile.WriteLine(ASETools.tumorToString[tumor] + " " + ASETools.dnaToString[dna] + " distribution of mean read length");
                    meanReadLengthHistogram[tumor] [dna].WriteHistogram(outputFile);
                    outputFile.WriteLine(ASETools.tumorToString[tumor] + " " + ASETools.dnaToString[dna] + " distribution of median read length");
                    medianReadLengthHistogram[tumor] [dna].WriteHistogram(outputFile);
                    outputFile.WriteLine(ASETools.tumorToString[tumor] + " " + ASETools.dnaToString[dna] + " distribution of min insert length");
                    minInsert[tumor] [dna].WriteHistogram(outputFile);
                    outputFile.WriteLine(ASETools.tumorToString[tumor] + " " + ASETools.dnaToString[dna] + " distribution of max insert length");
                    maxInsert[tumor] [dna].WriteHistogram(outputFile);
                    outputFile.WriteLine(ASETools.tumorToString[tumor] + " " + ASETools.dnaToString[dna] + " distribution of mean insert length");
                    meanInsert[tumor] [dna].WriteHistogram(outputFile);
                    outputFile.WriteLine(ASETools.tumorToString[tumor] + " " + ASETools.dnaToString[dna] + " distribution of median insert length");
                    medianInsert[tumor] [dna].WriteHistogram(outputFile);
                } // dna
            } // tumor

            outputFile.WriteLine("**done**");
            outputFile.Close();

            Console.WriteLine("Processed " + commonData.listOfCases.Count() + " cases in " + ASETools.ElapsedTimeInSeconds(commonData.timer));
        } // Main
    }
}
