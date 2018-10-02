using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using ASELib;
using System.IO;

namespace ReadLengthDistributions
{
    class Program
    {
        class PerThreadState
        {
            public Dictionary<int,Dictionary<bool, Dictionary<bool, ASETools.PreBucketedHistogram>>> histograms = new Dictionary<int,Dictionary<bool, Dictionary<bool, ASETools.PreBucketedHistogram>>>();  // {0 min, 1 average, 2 max} ->tumor->dna->histogram

            public PerThreadState()
            {
                for (int mam = 0; mam <= 2; mam++)
                {
                    histograms.Add(mam, new Dictionary<bool, Dictionary<bool, ASETools.PreBucketedHistogram>>());
                    foreach (var tumor in ASETools.BothBools)
                    {
                        histograms[mam].Add(tumor, new Dictionary<bool, ASETools.PreBucketedHistogram>());
                        foreach (var dna in ASETools.BothBools)
                        {
                            histograms[mam][tumor].Add(dna, new ASETools.PreBucketedHistogram(0, 200, 1));
                        } // dna
                    } // mutant
                } // mam
            } // ctor

            public void merge(PerThreadState peer)
            {
                for (int mam = 0; mam <= 2; mam++)
                {
                    foreach (var tumor in ASETools.BothBools)
                    {
                        foreach (var dna in ASETools.BothBools)
                        {
                            histograms[mam][tumor][dna].merge(peer.histograms[mam][tumor][dna]);
                        } // dna
                    } // mutant
                } // mam
            } // merge
        } // PerThreadState

        static PerThreadState globalState = new PerThreadState();
        static ASETools.CommonData commonData;

        static void Main(string[] args)
        {
            commonData = ASETools.CommonData.LoadCommonData(args);

            if (null == commonData)
            {
                return;
            }

            if (commonData.cases.Select(_ => _.Value).Any(_ => _.normal_dna_reads_at_selected_variants_filename == "" || _.tumor_dna_reads_at_selected_variants_filename == "" || _.tumor_rna_reads_at_selected_variants_filename == "" ||
                                                          _.selected_variants_filename == ""))
            {
                Console.WriteLine("Some cases are missing reads at selected variants or selected variants files.  Make them and try again.");
                return;
            }

            var threading = new ASETools.WorkerThreadHelper<ASETools.Case, PerThreadState>(commonData.cases.Select(_ => _.Value).ToList(), HandleOneCase, FinishUp, null, 100);

            Console.WriteLine("Processing " + commonData.cases.Count() + " cases, one dot/100 cases.");
            ASETools.PrintNumberBar(commonData.cases.Count() / 100);

            threading.run();

            Console.WriteLine();

            var outputFilename = commonData.configuration.finalResultsDirectory + ASETools.ReadLengthHistogramFilename;
            var outputFile = ASETools.CreateStreamWriterWithRetry(outputFilename);

            if (null == outputFile)
            {
                Console.WriteLine("Unable to open " + outputFilename);
                return;
            }

            for (int mam = 0; mam <= 2; mam++)
            {
                foreach (var tumor in ASETools.BothBools)
                {
                    foreach (var dna in ASETools.BothBools)
                    {
                        if (0 == mam)
                        {
                            outputFile.Write("Min");
                        } else if (1 == mam)
                        {
                            outputFile.Write("Mean");
                        } else
                        {
                            outputFile.Write("Max");
                        }

                        outputFile.WriteLine(" read depth distribution for tumor: " + tumor + " dna: " + dna);
                        outputFile.WriteLine(ASETools.HistogramResultLine.Header());
                        globalState.histograms[mam][tumor][dna].ComputeHistogram().ToList().ForEach(_ => outputFile.WriteLine(_.ToString()));
                        outputFile.WriteLine();
                    } // dna
                } // tumor
            } // mam

            outputFile.WriteLine("**done**");
            outputFile.Close();

            Console.WriteLine("Elapsed time " + ASETools.ElapsedTimeInSeconds(commonData.timer));
        } // Main

        static void HandleOneCase(ASETools.Case case_, PerThreadState perThreadState)
        {
            var selectedVariants = ASETools.SelectedVariant.LoadFromFile(case_.selected_variants_filename);

            foreach (var tumor in ASETools.BothBools)
            {
                foreach (var dna in ASETools.BothBools)
                {
                    if (case_.getReadsAtSelectedVariantsFilename(tumor, dna) == "")
                    {
                        //
                        // This is normal DNA for a case that doesn't have it.
                        //
                        continue;
                    }

                    int nVariantsChecked = 0;
                    int nReadsChecked = 0;
                    int minLength = 200;
                    int maxLength = 0;
                    int totalLength = 0;

                    var consolodatedFile = new ASETools.ConsolodatedFileReader();

                    if (!consolodatedFile.open(case_.getReadsAtSelectedVariantsFilename(tumor, dna)))
                    {
                        Console.WriteLine("Unable to open reads at selected variants file " + case_.getReadsAtSelectedVariantsFilename(tumor, dna));
                        continue;
                    }

                    foreach (var variant in selectedVariants)
                    {
                        var subfileName = case_.getDownloadedReadsFileId(tumor, dna) + variant.getExtractedReadsExtension();
                        var subfileReader = consolodatedFile.getSubfile(subfileName);
                        if (null == subfileReader)
                        {
                            Console.WriteLine("ComputeReadCounts: no subfile named " + subfileName);
                            continue;
                        }

                        string line;
                        while (null != (line = subfileReader.ReadLine()))
                        {
                            ASETools.SAMLine samLine;

                            try
                            {
                                samLine = new ASETools.SAMLine(line);
                            }
                            catch (Exception e)
                            {
                                Console.WriteLine("Unable to parse sam line in extracted reads main file + " + case_.getReadsAtSelectedVariantsFilename(tumor, dna) + " subfile " + subfileName + ": " + line);
                                throw e;
                            }

                            var length = samLine.nonclippedBases;
                            minLength = Math.Min(minLength, length);
                            maxLength = Math.Max(maxLength, length);
                            totalLength += length;
                            nReadsChecked++;
                        } // for each read
                        nVariantsChecked++;

                        subfileReader.Close();

                        if (nReadsChecked > 200 && nVariantsChecked > 10)
                        {
                            perThreadState.histograms[0][tumor][dna].addValue(minLength);
                            perThreadState.histograms[1][tumor][dna].addValue(totalLength / nReadsChecked); // Rounds down, but whatever
                            perThreadState.histograms[2][tumor][dna].addValue(maxLength);
                            break;
                        }
                    } // foreach variant
                } // dna
            } // mutant
        } // HandleOneCase

        static void FinishUp(PerThreadState perThreadState)
        {
            lock (globalState)
            {
                globalState.merge(perThreadState);
            }
        } // FinishUp
    }
}
