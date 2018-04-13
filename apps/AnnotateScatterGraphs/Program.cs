using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using System.Diagnostics;
using System.IO;
using ASELib;

namespace AnnotateScatterGraphs
{
    class Program
    {
        static ASETools.Configuration configuration;

        static void Main(string[] args)
        {
            var timer = new Stopwatch();
            timer.Start();

            configuration = ASETools.Configuration.loadFromFile(args);

            if (configuration == null)
            {
                Console.WriteLine("Unable to load configuration.");
                return;
            }

            var subTimer = new Stopwatch();
            subTimer.Start();
            Console.Write("Loading cases...");
            var cases = ASETools.Case.LoadCases(configuration.casesFilePathname).Select(x => x.Value).ToList();
            if (cases == null)
            {
                Console.WriteLine("Unable to load cases.");
                return;
            }
            Console.WriteLine(ASETools.ElapsedTimeInSeconds(subTimer));

            subTimer.Stop();
            subTimer.Reset();
            subTimer.Start();
            Console.Write("Loading RNA mapped base counts...");
            var tumorRNAMappedBaseCounts = new Dictionary<string, ASETools.MappedBaseCount>();
            foreach (var case_ in cases)
            {
                if (case_.tumor_rna_mapped_base_count_filename != "")
                {
                    tumorRNAMappedBaseCounts.Add(case_.case_id, ASETools.MappedBaseCount.readFromFile(case_.tumor_rna_mapped_base_count_filename));
                }
            }
            Console.WriteLine(ASETools.ElapsedTimeInSeconds(subTimer));



            subTimer.Stop();
            subTimer.Reset();
            subTimer.Start();
            Console.Write("Loading scatter graph lines...");
            var scatterGraphLines = ASETools.GeneScatterGraphLine.LoadAllGeneScatterGraphLines(configuration.geneScatterGraphsDirectory, false, "*");
            Console.WriteLine(ASETools.ElapsedTimeInSeconds(subTimer));

            if (scatterGraphLines == null)
            {
                Console.WriteLine("Unable to load scatter graph lines");
                return;
            }

            scatterGraphLines.Sort((x, y) => x.disease.CompareTo(y.disease));   // Sorted by disease so we don't have to load the expression distributions more than once.
            var nScatterGraphLines = scatterGraphLines.Count();

            ASETools.ExpressionDistributionMap expressionDistributionMap = null;
            string loadedDisease = "";

            for (int i = 0; i < nScatterGraphLines; i++)
            {
                var scatterGraphLine = scatterGraphLines[i];

                if (scatterGraphLine.disease != loadedDisease)
                {
                    var filename = ASETools.Configuration.expression_distribution_directory + ASETools.Expression_distribution_filename_base + scatterGraphLine.disease;
                    Console.Write("Loading expression map for " + scatterGraphLine.disease + " from " + filename + "...");
                    var loadTimer = new Stopwatch();
                    loadTimer.Start();

                    expressionDistributionMap = null;   // So GC can get it while we're loading the next one
                    expressionDistributionMap = new ASETools.ExpressionDistributionMap(filename);
                    Console.WriteLine(ASETools.ElapsedTimeInSeconds(loadTimer));
                    loadedDisease = scatterGraphLine.disease;
                }

                if (expressionDistributionMap.map.ContainsKey(scatterGraphLine.Chromosome) && expressionDistributionMap.map[scatterGraphLine.Chromosome].ContainsKey(scatterGraphLine.Start_Position)) {
                    var percentiles = expressionDistributionMap.map[scatterGraphLine.Chromosome][scatterGraphLine.Start_Position];
                    if (scatterGraphLine.tumorRNAReadCounts != null && tumorRNAMappedBaseCounts.ContainsKey(scatterGraphLine.case_id))
                    {
                        scatterGraphLine.tumorRNAFracAltPercentile = new double[11];
                        scatterGraphLine.tumorRNAFracRefPercentile = new double[11];
                        scatterGraphLine.tumorRNAFracAllPercentile = new double[11];

                        for (int j = 0; j < 11; j++)
                        {
                            if (percentiles.getPercentile(j * 10) == 0)
                            {
                                scatterGraphLine.tumorRNAFracAltPercentile[j] = double.NegativeInfinity;
                                scatterGraphLine.tumorRNAFracRefPercentile[j] = double.NegativeInfinity;
                                scatterGraphLine.tumorRNAFracAllPercentile[j] = double.NegativeInfinity;
                            } else
                            {
                                scatterGraphLine.tumorRNAFracAltPercentile[j] = ((double)scatterGraphLine.tumorRNAReadCounts.nMatchingAlt / tumorRNAMappedBaseCounts[scatterGraphLine.case_id].mappedBaseCount) / percentiles.getPercentile(j * 10);
                                scatterGraphLine.tumorRNAFracRefPercentile[j] = ((double)scatterGraphLine.tumorRNAReadCounts.nMatchingReference / tumorRNAMappedBaseCounts[scatterGraphLine.case_id].mappedBaseCount) / percentiles.getPercentile(j * 10);
                                scatterGraphLine.tumorRNAFracAllPercentile[j] = (((double)scatterGraphLine.tumorRNAReadCounts.totalReads()) / tumorRNAMappedBaseCounts[scatterGraphLine.case_id].mappedBaseCount) / percentiles.getPercentile(j * 10);
                            }
                        }
                    }
                }
            } // for each scatter graph line.


            var scatterGraphLinesByGene = scatterGraphLines.GroupByToDict(x => x.Hugo_Symbol);

            Console.Write("Generating output for " + scatterGraphLinesByGene.Count() + " genes...");
            subTimer.Stop();
            subTimer.Reset();
            subTimer.Start();

            var combinedHistogramsFile = ASETools.CreateStreamWriterWithRetry(configuration.geneScatterGraphsDirectory + ASETools.annotated_scatter_graphs_histogram_filename);

            var genesWithData = new List<string>();

            foreach (var geneEntry in scatterGraphLinesByGene)
            {
                genesWithData.Add(geneEntry.Key);
            }

            genesWithData.Sort();   // So they're in alphabetical order in the histograms file.

            foreach (var hugo_symbol in genesWithData) {
                var lines = scatterGraphLinesByGene[hugo_symbol];
                StreamWriter outputFile = null;

                lines.Sort((x, y) => (x.tumorRNAFracRefPercentile == null) || y.tumorRNAFracRefPercentile == null ? 0 : ((x.MultipleMutationsInThisGene == y.MultipleMutationsInThisGene) ? x.tumorRNAFracRefPercentile[5].CompareTo(y.tumorRNAFracRefPercentile[5]) : x.MultipleMutationsInThisGene.CompareTo(y.MultipleMutationsInThisGene)));

                var refHistograms = new Dictionary<bool, ASETools.PreBucketedHistogram>();
                var altHistograms = new Dictionary<bool, ASETools.PreBucketedHistogram>();
                var allHistograms = new Dictionary<bool, ASETools.PreBucketedHistogram>();

                foreach (var multiple in ASETools.BothBools)
                {
                    refHistograms.Add(multiple, new ASETools.PreBucketedHistogram(0, 2.02, 0.02));
                    altHistograms.Add(multiple, new ASETools.PreBucketedHistogram(0, 2.02, 0.02));
                    allHistograms.Add(multiple, new ASETools.PreBucketedHistogram(0, 2.02, 0.02));
                }

                foreach (var line in lines)
                {

                    if (line.tumorRNAFracAltPercentile == null)
                    {
                        continue;
                    }

                    if (outputFile == null)
                    {
                        outputFile = ASETools.CreateStreamWriterWithRetry(configuration.geneScatterGraphsDirectory + hugo_symbol + ASETools.annotated_scatter_graph_filename_extension);
                        if (outputFile == null)
                        {
                            Console.WriteLine("Unable to open output file for " + hugo_symbol);
                            continue;
                        }
                        outputFile.WriteLine(ASETools.GeneScatterGraphLine.annotatedHeaderLine);
                    }

                    outputFile.Write(line.rawInputLine);
                    for (int i = 0; i < 11; i++)
                    {
                        if (line.tumorRNAFracAltPercentile[i] == double.NegativeInfinity)
                        {
                            outputFile.Write("\t*");
                        }
                        else
                        {
                            outputFile.Write("\t" + line.tumorRNAFracRefPercentile[i]);
                            if (5 == i)
                            {
                                refHistograms[line.MultipleMutationsInThisGene].addValue(line.tumorRNAFracRefPercentile[5] * 2);
                            }
                        }
                    }
                    for (int i = 0; i < 11; i++)
                    {
                        if (line.tumorRNAFracAltPercentile[i] == double.NegativeInfinity)
                        {
                            outputFile.Write("\t*");
                        }
                        else
                        {
                            outputFile.Write("\t" + line.tumorRNAFracAltPercentile[i]);
                            if (5 == i)
                            {
                                altHistograms[line.MultipleMutationsInThisGene].addValue(line.tumorRNAFracAltPercentile[5] * 2);
                            }
                        }
                    }
                    for (int i = 0; i < 11; i++)
                    {
                        if (line.tumorRNAFracAltPercentile[i] == double.NegativeInfinity)
                        {
                            outputFile.Write("\t*");
                        }
                        else
                        {
                            outputFile.Write("\t" + line.tumorRNAFracAllPercentile[i]);
                            if (5 == i)
                            {
                                allHistograms[line.MultipleMutationsInThisGene].addValue(line.tumorRNAFracAllPercentile[5]);
                            }
                        }
                    }
                    outputFile.WriteLine();
                } // foreach line

                if (refHistograms[true].count() != 0 || refHistograms[false].count() != 0) 
                {
                    combinedHistogramsFile.WriteLine(hugo_symbol + "\twild-type one mutation (n = " + refHistograms[false].count() + ")\tmutant one mutation (n = " + altHistograms[false].count() + ")\tall one mutation (n = " + allHistograms[false].count() + ")\twild-type > 1 mutation (n = " +
                        refHistograms[true].count() + ")\tmutant > 1 mutation (n = " + altHistograms[true].count() + ")\tall > 1 mutation (n = " + allHistograms[true].count() + ")");

                    var cdfs = new List<double>[6];

                    cdfs[0] = refHistograms[false].ComputeHistogram().Select(x => x.cdfValue).ToList();
                    cdfs[1] = altHistograms[false].ComputeHistogram().Select(x => x.cdfValue).ToList();
                    cdfs[2] = allHistograms[false].ComputeHistogram().Select(x => x.cdfValue).ToList();
                    cdfs[3] = refHistograms[true].ComputeHistogram().Select(x => x.cdfValue).ToList();
                    cdfs[4] = altHistograms[true].ComputeHistogram().Select(x => x.cdfValue).ToList();
                    cdfs[5] = allHistograms[true].ComputeHistogram().Select(x => x.cdfValue).ToList();

                    for (int i =0; i < cdfs[0].Count(); i++)
                    {
                        if (i == cdfs[0].Count() - 1)
                        {
                            combinedHistogramsFile.Write("More");
                        } else
                        {
                            combinedHistogramsFile.Write(0.02 * i);
                        }

                        for (int j = 0; j < 6; j++)
                        {
                            combinedHistogramsFile.Write("\t" + cdfs[j][i]);
                        }
                        combinedHistogramsFile.WriteLine();
                    }

                    combinedHistogramsFile.WriteLine();
                } // If we had anything for histograms

                if (outputFile != null)
                {
                    outputFile.WriteLine("**done**");
                    outputFile.Close();
                }
            } // foreach gene

            combinedHistogramsFile.WriteLine("**done**");
            combinedHistogramsFile.Close();

            Console.WriteLine(ASETools.ElapsedTimeInSeconds(subTimer));

            Console.WriteLine("Overall run time " + ASETools.ElapsedTimeInSeconds(timer));
        }
    }
}
