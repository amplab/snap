using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using ASELib;
using System.IO;
using System.Diagnostics;

namespace MakeSignificant012Graphs
{
    class Program
    {

        static List<ResultToReport> resultsToReport = new List<ResultToReport>();
        static int nSignificantResults;

        const int nHistogramsPerResult = 4; // 0, 1, > 1, 1 silent.  
        const double rangeStep = 0.01;

        class ResultToReport : IComparable<ResultToReport>
        {
            public ResultToReport(ASETools.BonferroniCorrectedASEDistributionLine line_, ASETools.GeneLocationsByNameAndChromosome geneLocationInformation, int whichResult, bool significant_, string overrideRange = "")
            {
                line = line_;
                significant = significant_;

                if (overrideRange != "")
                {
                    range = overrideRange;  // This is a hack to allow looking at genes/ranges that aren't significant
                }
                else
                {
                    range = line.significantAtArray[whichResult];
                }


                var histogramNameTrailer = " mutations in " + line.hugo_symbol + " at range index " + range;
                histograms[0] = new ASETools.PreBucketedHistogram(0, 1, rangeStep, "0 " + histogramNameTrailer);

                for (int i = 0; i < nHistogramsPerResult; i++)
                {
                    histograms[i] = new ASETools.PreBucketedHistogram(0, 1, rangeStep, "" + ((i < 2) ? i.ToString() : ((i == 2) ? ">1" : "1 silent")) + histogramNameTrailer);
                }

                var gene = geneLocationInformation.genesByName[line.hugo_symbol];
                chromosome = gene.chromosome;

                //
                // The range is a number from 0 to 20 possibly followed by an E.  0 means in gene, 20 means whole autosome,
                // x = 1-19 means a range of 1000 * 2^(x-1) on either side of the gene.  E means exclusive, which for 1 - 19 means
                // also not in the next smaller range, and for 20 means chromosomes other than the one with the gene.
                //

                int rangeIndex;
                bool exclusive;
                if (exclusive = range.EndsWith("E"))
                {
                    rangeIndex = Convert.ToInt32(range.Substring(0, range.Count() - 1));
                } else
                {
                    rangeIndex = Convert.ToInt32(range);
                }


                if (rangeIndex == 0)
                {
                    minLocus[0] = gene.minLocus;
                    maxLocus[0] = gene.maxLocus;
                    minLocus[1] = maxLocus[1] = -1;
                } else if (rangeIndex == 20)
                {
                    countsIfInRange = false;
                    if (exclusive)
                    {
                        minLocus[0] = 0;
                        maxLocus[0] = int.MaxValue;
                    } else
                    {
                        minLocus[0] = maxLocus[0] = -1;
                    }
                    minLocus[1] = maxLocus[1] = -1;
                } else
                {
                    int rangeSize = 1000 * (1 << (rangeIndex - 1));

                    minLocus[0] = Math.Max(0, gene.minLocus - rangeSize);

                    if (exclusive)
                    {
                        maxLocus[1] = gene.maxLocus + rangeSize;
                        if (rangeSize == 1)
                        {
                            maxLocus[0] = gene.minLocus;
                            minLocus[1] = gene.maxLocus;
                        }
                        else
                        {
                            maxLocus[0] = gene.minLocus - rangeSize / 2;
                            minLocus[1] = gene.maxLocus + rangeSize / 2;
                        }
                    }
                    else
                    {
                        maxLocus[0] = gene.maxLocus + rangeSize;

                        minLocus[1] = maxLocus[1] = -1;
                    }
                }
            } // SignificantResult ctor

            public int CompareTo(ResultToReport peer)
            {
                if (peer.line.hugo_symbol != line.hugo_symbol)
                {
                    return line.hugo_symbol.CompareTo(peer.line.hugo_symbol);
                }

                bool peerExclusive = peer.range.EndsWith("E");
                bool exclusive = range.EndsWith("E");

                if (exclusive != peerExclusive)
                {
                    if (exclusive)
                    {
                        return 1;
                    }

                    return -1;
                }

                return ASETools.RangeValueFromRange(range).CompareTo(ASETools.RangeValueFromRange(peer.range));
            } // CompareTo

            public readonly ASETools.BonferroniCorrectedASEDistributionLine line;
            public readonly string chromosome;
            public readonly int[] minLocus = new int[2];
            public readonly int[] maxLocus = new int [2];
            public readonly bool countsIfInRange = true;
            public readonly string range;
            public readonly bool significant;

            public ASETools.PreBucketedHistogram[] histograms = new ASETools.PreBucketedHistogram[nHistogramsPerResult];
        }

        class PerThreadState {
            public PerThreadState()
            {
                histograms = new ASETools.PreBucketedHistogram[nSignificantResults, nHistogramsPerResult];
                for (int i = 0; i < nSignificantResults; i++)
                {
                    for (int j = 0; j < nHistogramsPerResult; j++)
                    {
                        histograms[i, j] = new ASETools.PreBucketedHistogram(0, 1, 0.01);
                    }
                }
            }
            public ASETools.PreBucketedHistogram [,]histograms;
        }

        static void FinishUp(PerThreadState state)
        {
            lock (resultsToReport)
            {
                for (int i = 0; i < nSignificantResults; i++)
                {
                    for (int j = 0; j < nHistogramsPerResult; j++)
                    {
                        resultsToReport[i].histograms[j].merge(state.histograms[i, j]);
                    }
                }
            }
        } // FinishUp

        static void HandleOneCase(ASETools.Case case_, PerThreadState state, List<ASETools.AnnotatedVariant> annotatedSelectedVariants, ASETools.ASECorrection ASECorrection, Dictionary<bool, List<ASETools.CopyNumberVariation>> copyNumber)
        {

            var somaticVariants = annotatedSelectedVariants.Where(x => x.somaticMutation && !x.isSilent()).ToList();
            var germlineVariants = annotatedSelectedVariants.Where(x => x.somaticMutation == false && x.IsASECandidate(true, copyNumber, configuration, perGeneASEMap, geneMap)).ToList();

            for (int i = 0; i < nSignificantResults; i++)
            {
                var result = resultsToReport[i];
                int nValid = 0;
                double totalASE = 0;

                foreach (var variant in germlineVariants)
                {
                    var inRange = variant.contig == result.chromosome && (variant.locus >= result.minLocus[0] && variant.locus <= result.maxLocus[0] || variant.locus >= result.minLocus[1] && variant.locus <= result.maxLocus[1]);

                    if (inRange == result.countsIfInRange)
                    {
                        nValid++;
                        totalASE += variant.GetAlleleSpecificExpression(true, aseCorrection);
                    }
                }

                if (nValid > 0)
                {
                    int mutationCountIndex = ASETools.ZeroOneMany(somaticVariants.Where(x => x.Hugo_symbol == result.line.hugo_symbol).Count());
                    result.histograms[mutationCountIndex].addValue(totalASE / nValid);

                    // Handle the one silent mutation count case.
                    if (annotatedSelectedVariants.Where(x => x.somaticMutation && x.isSilent() && x.Hugo_symbol == result.line.hugo_symbol).Count() == 1)
                    {
                        result.histograms[3].addValue(totalASE / nValid);
                    }
                }
            } // for each significant result
        } // HandleOneCase

        static ASETools.Configuration configuration;
        static Dictionary<string, ASETools.ASEMapPerGeneLine> perGeneASEMap;
        static ASETools.ASECorrection aseCorrection;
        static ASETools.GeneLocationsByNameAndChromosome geneLocationInformation;
        static ASETools.GeneMap geneMap;

        static void Main(string[] args)
        {
            var timer = new Stopwatch();
            timer.Start();

            configuration = ASETools.Configuration.loadFromFile(args);
            if (null == configuration)
            {
                Console.WriteLine("Unable to load configuration.");
                return;
            }

            perGeneASEMap = ASETools.ASEMapPerGeneLine.ReadFromFileToDictionary(configuration.finalResultsDirectory + ASETools.PerGeneASEMapFilename);

            if (null == perGeneASEMap)
            {
                Console.WriteLine("You must first create the per-gene ASE map in " + configuration.finalResultsDirectory + ASETools.PerGeneASEMapFilename);
                return;
            }

            aseCorrection = ASETools.ASECorrection.LoadFromFile(configuration.finalResultsDirectory + ASETools.ASECorrectionFilename);
            if (null == aseCorrection)
            {
                Console.WriteLine("Unable to load ASE correction");
                return;
            }

            geneLocationInformation = new ASETools.GeneLocationsByNameAndChromosome(ASETools.readKnownGeneFile(configuration.geneLocationInformationFilename));
            geneMap = new ASETools.GeneMap(geneLocationInformation.genesByName);


            var bonferroniLines = ASETools.BonferroniCorrectedASEDistributionLine.readFromFile(configuration.finalResultsDirectory + "AlleleSpecificExpressionDistributionByMutationCount_bonferroni.txt");


            if (null == bonferroniLines)
            {
                Console.WriteLine("Unable to load Bonferroni lines.");
                return;
            }

            string[] genesAlwaysToReport = { "APC", "BRCA1", "BRCA2", "CD95", "ST5", "YPEL3", "ST7", "ST14",    // Tumor suppressors
                "EGFR", "VEGFR", "MYC", "HER2", "PDGF", "IDH1", "IDH2", // Oncogenes 
                "ATCB", "HK1"                                  // controls
            };

            foreach (var bonferroniLine in bonferroniLines)
            {
                if (bonferroniLine.significantAtArray.Count() > 0 && bonferroniLine.significant01 || genesAlwaysToReport.Contains(bonferroniLine.hugo_symbol))
                {
                    bool saw5 = false;
                    bool saw11 = false;
                    bool saw0 = false;

                    for (int i = 0; i < bonferroniLine.significantAtArray.Count(); i++)
                    {
                        if (bonferroniLine.significantAtArray[i] != "" && (i == 0 || bonferroniLine.significantAtArray[i - 1] != bonferroniLine.significantAtArray[i])) // The equality check is because there can be duplicates if 1 vs. not 1 and 1 vs. many are both significant
                        {
                            saw5 |= bonferroniLine.significantAtArray[i] == "5";
                            saw11 |= bonferroniLine.significantAtArray[i] == "11";
                            saw0 |= bonferroniLine.significantAtArray[i] == "0";
                            resultsToReport.Add(new ResultToReport(bonferroniLine, geneLocationInformation, i, true));
                        }
                    }

                    if (!saw5)
                    {
                        resultsToReport.Add(new ResultToReport(bonferroniLine, geneLocationInformation, 0, false, "5"));
                    }

                    if (!saw11)
                    {
                        resultsToReport.Add(new ResultToReport(bonferroniLine, geneLocationInformation, 0, false, "11"));
                    }

                    if (!saw0)
                    {
                        resultsToReport.Add(new ResultToReport(bonferroniLine, geneLocationInformation, 0, false, "0"));
                    }
                } // If there are any significant results for this gene.

            }

            resultsToReport.Sort();

            nSignificantResults = resultsToReport.Count();

            var cases = ASETools.Case.LoadCases(configuration.casesFilePathname);

            if (null == cases)
            {
                Console.WriteLine("Unable to load cases.");
                return;
            }

            var casesToProcess = cases.Select(x => x.Value).Where(x => x.annotated_selected_variants_filename != "" && x.tumor_copy_number_filename != "").ToList();

            Console.Write("Processing " + casesToProcess.Count() + " of " + cases.Count() + " cases and " + resultsToReport.Count() + " significant results, 1 dot/100 cases: ");

            var threading = new ASETools.ASVThreadingHelper<PerThreadState>(casesToProcess, aseCorrection, (x, y) => true, HandleOneCase, FinishUp, null, 100);
            threading.run();

            Console.WriteLine("Took " + ASETools.ElapsedTimeInSeconds(timer));

            var nLinesPerHistogram = resultsToReport[0].histograms[0].ComputeHistogram().Count();

            foreach (var resultToReport in resultsToReport)
            {
                var outputFile = ASETools.CreateStreamWriterWithRetry(ASETools.Configuration.zero_one_two_directory + resultToReport.line.hugo_symbol + "-" + resultToReport.range + (resultToReport.significant ? "" : "I") + ".txt");

                //
                // Pull the three CDF columns at the top of the file next to one another, since that's what we make into graphs.
                //
                outputFile.WriteLine("CDFs for 0, 1 and >1 mutations and 1 silent mutation.");
                ASETools.HistogramResultLine[][] histogramLines = new ASETools.HistogramResultLine[nHistogramsPerResult][];
                for (int i = 0; i < nHistogramsPerResult; i++)
                {
                    histogramLines[i] = resultToReport.histograms[i].ComputeHistogram();
                }

                outputFile.WriteLine("minValue\t0 mutations (n = " + histogramLines[0].Select(x => x.count).Sum() + ")\t1 mutation (n = " + histogramLines[1].Select(x => x.count).Sum() + ")\t> 1 mutation (n = " + histogramLines[2].Select(x => x.count).Sum() + ")" +
                    "\t1 silent mutation (n = " + histogramLines[3].Select(x => x.count).Sum() + ")");
                for (int whichLine = 0; whichLine < nLinesPerHistogram; whichLine++)
                {
                    if (whichLine == nLinesPerHistogram - 1)
                    {
                        outputFile.Write("1");  // This replaces "More," which in turn eliminates one more step in making the graphs
                    }
                    else
                    {
                        outputFile.Write(histogramLines[0][whichLine].minValue);
                    }
                    for (int i = 0; i < nHistogramsPerResult; i++)
                    {
                        outputFile.Write("\t" + histogramLines[i][whichLine].cdfValue);
                    }
                    outputFile.WriteLine();
                }
                outputFile.WriteLine();


                for (int i = 0; i < nHistogramsPerResult; i++)
                {
                    outputFile.WriteLine(resultToReport.histograms[i].name);
                    outputFile.WriteLine(ASETools.HistogramResultLine.Header());
                    resultToReport.histograms[i].ComputeHistogram().ToList().ForEach(x => outputFile.WriteLine(x));
                    outputFile.WriteLine();
                }
                outputFile.WriteLine("**done**");
                outputFile.Close();
            }


            //
            // Now do per-gene files for genes with more than one significant result other than whole autosome.  CountsIfInRange is false iff whole autosome.
            //
            var multipleResults = resultsToReport.Where(x => x.countsIfInRange).GroupByToDict(x => x.line.hugo_symbol).Where(x => x.Value.Count() > 1);

            foreach (var gene in multipleResults)
            {
                var lines = gene.Value;
                var hugoSymbol = gene.Key;
                var outputFile = ASETools.CreateStreamWriterWithRetry(ASETools.Configuration.zero_one_two_directory + hugoSymbol + "-all.txt");

                for (int i = 0; i < nHistogramsPerResult; i++)
                {
                    if (i == nHistogramsPerResult -1)
                    {
                        outputFile.WriteLine(">" + (i - 1) + " mutations");
                    } else
                    {
                        outputFile.WriteLine("" + i + " mutation" + ((i == 0) ? "s" : ""));
                    }

                    var histograms = lines.Select(x => x.histograms[i].ComputeHistogram()).ToList();

                    outputFile.Write("minValue");
                    for (int whichLine = 0; whichLine < lines.Count(); whichLine++)
                    {
                        var line = lines[whichLine];
                        outputFile.Write("\t" + ASETools.RangeToDescriptiveString(line.range) + " (n= " + histograms[whichLine].Select(x => x.count).Sum() + ")");
                    }
                    outputFile.WriteLine();

                    for (int whichBucket = 0; whichBucket < nLinesPerHistogram; whichBucket++) {
                        outputFile.Write((histograms[0][whichBucket].minValue == "More" ? "1" : histograms[0][whichBucket].minValue));

                        foreach (var histogram in histograms)
                        {
                            outputFile.Write("\t" + histogram[whichBucket].cdfValue);
                        }
                        outputFile.WriteLine();
                    }

                    outputFile.WriteLine();
                }

                outputFile.WriteLine("**done**");
                outputFile.Close();
            }
        }
    }
}
