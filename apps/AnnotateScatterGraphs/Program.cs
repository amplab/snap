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
#if false
    public class DiseaseAndChromosome : IComparable<DiseaseAndChromosome>
    {
        public readonly string disease;
        public readonly string chromosome;   // in short form

        public DiseaseAndChromosome(string disease_, string chromosome_)
        {
            disease = disease_;
            chromosome = chromosome_;
        }

        public int CompareTo(DiseaseAndChromosome peer)
        {
            if (disease != peer.disease)
            {
                return disease.CompareTo(peer.disease);
            }

            return chromosome.CompareTo(peer.chromosome);
        } // CompareTo


    } // DiseaseAndChromosome


    static class Extensions
    {
        public static DiseaseAndChromosome getDiseaseAndChromosome(this ASETools.GeneScatterGraphLine geneScatterGraphLine)
        {
            return new DiseaseAndChromosome(geneScatterGraphLine.disease, geneScatterGraphLine.Chromosome);
        }
    }
#endif // false

    class Program
    {
        static ASETools.Configuration configuration;

        class DataPoint : IComparer<DataPoint>
        {
            public readonly double value;
            public readonly bool alt;   // Or else ref
            public DataPoint(double value_, bool alt_)
            {
                value = value_;
                alt = alt_;
            }

            public int Compare(DataPoint x, DataPoint y)
            {
                return x.value.CompareTo(y.value);
            }
        }

        class HugoPAndMedians
        {
            public readonly string hugo_symbol;
            public readonly double p;
            public readonly double refMedian;
            public readonly double altMedian;

            public HugoPAndMedians(string hugo, double p_, double refMedian_, double altMedian_)
            {
                hugo_symbol = hugo;
                p = p_;
                refMedian = refMedian_;
                altMedian = altMedian_;
            }
        }

        static ASETools.CommonData commonData;
        static List<ASETools.GeneScatterGraphLine> scatterGraphLines;
        static Dictionary<string, ASETools.MappedBaseCount> tumorRNAMappedBaseCounts;

        static void Main(string[] args)
        {
            commonData = ASETools.CommonData.LoadCommonData(args);
            if (null == commonData)
            {
                return;
            }

            configuration = commonData.configuration;

            var subTimer = new Stopwatch();
            subTimer.Start();
            var cases = commonData.listOfCases;

            if (commonData.diseases.Any(disease => ASETools.chromosomes.Any(chromosome => !File.Exists(commonData.configuration.geneScatterGraphsLinesWithPercentilesDirectory + ASETools.GeneScatterGraphLinesWithPercentilesPrefix + disease + "_" + ASETools.chromosomeNameToNonChrForm(chromosome)))))
            {
                Console.WriteLine("Not all of the gene scatter lines by disease and chromsome (the percentile ones) exist.  Create them and try again.");
                return;
            }

            Console.Write("Loading RNA mapped base counts...");
            tumorRNAMappedBaseCounts = new Dictionary<string, ASETools.MappedBaseCount>();
            foreach (var case_ in cases)
            {
                if (case_.tumor_rna_mapped_base_count_filename != "")
                {
                    tumorRNAMappedBaseCounts.Add(case_.case_id, ASETools.MappedBaseCount.readFromFile(case_.tumor_rna_mapped_base_count_filename));
                }
            }
            Console.WriteLine(ASETools.ElapsedTimeInSeconds(subTimer));

            var perGeneASEMap = commonData.perGeneASEMap;

            if (null == perGeneASEMap)
            {
                Console.WriteLine("You must first create the per-gene ASE map in " + configuration.finalResultsDirectory + ASETools.PerGeneASEMapFilename);
                return;
            }

            var geneLocationInformation = commonData.geneLocationInformation;
            var geneMap = commonData.geneMap;

            subTimer.Stop();
            subTimer.Reset();
            subTimer.Start();
            Console.Write("Loading copy number files...");
            var copyNumberByCase = new Dictionary<string, List<ASETools.CopyNumberVariation>>();
            foreach (var case_ in cases)
            {
                copyNumberByCase.Add(case_.case_id, ASETools.CopyNumberVariation.ReadFile(case_.tumor_copy_number_filename));
            }
            Console.WriteLine(ASETools.ElapsedTimeInSeconds(subTimer));


            ASETools.ASERepetitiveRegionMap repetitiveRegionMap = null; // Don't use this exclusion criterion here.

            subTimer.Stop();
            subTimer.Reset();
            subTimer.Start();
            Console.Write("Loading scatter graph lines...");
            scatterGraphLines = ASETools.GeneScatterGraphLine.LoadAllLinesFromPercentilesDirectory(configuration.geneScatterGraphsLinesWithPercentilesDirectory);
            Console.WriteLine(ASETools.ElapsedTimeInSeconds(subTimer));

            if (scatterGraphLines == null)
            {
                Console.WriteLine("Unable to load scatter graph lines");
                return;
            }

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

            int nCandidates = 0;

            List<HugoPAndMedians> results = new List<HugoPAndMedians>();

            foreach (var hugo_symbol in genesWithData) {
                var lines = scatterGraphLinesByGene[hugo_symbol];
                StreamWriter outputFile = null;

                lines.Sort((x, y) => CompareScatterGraphLines(x, y, perGeneASEMap, geneMap, copyNumberByCase, repetitiveRegionMap));                

                var refHistograms = new Dictionary<bool, ASETools.Histogram>();
                var altHistograms = new Dictionary<bool, ASETools.Histogram>();
                var allHistograms = new Dictionary<bool, ASETools.Histogram>();

                var singleMutationDataPoints = new List<DataPoint>();

                foreach (var multiple in ASETools.BothBools)
                {
                    refHistograms.Add(multiple, new ASETools.Histogram());
                    altHistograms.Add(multiple, new ASETools.Histogram());
                    allHistograms.Add(multiple, new ASETools.Histogram());

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

                    string whyNot;
                    var aseCandidate = line.isASECandidate(out whyNot, copyNumberByCase[line.case_id], configuration, perGeneASEMap, geneMap, repetitiveRegionMap);

                    string whyNotIgnoringNMD = whyNot;
                    var aseCandidateIgnoringNMD = aseCandidate || line.isASECandidate(out whyNotIgnoringNMD, copyNumberByCase[line.case_id], configuration, perGeneASEMap, geneMap, repetitiveRegionMap, true);

                    outputFile.WriteLine("\t" + aseCandidate + "\t" + whyNot + "\t" + whyNotIgnoringNMD);

                    refHistograms[line.MultipleMutationsInThisGene].addValue(line.tumorRNAFracRefPercentile[5]);
                    if (!line.MultipleMutationsInThisGene)
                    {
                        singleMutationDataPoints.Add(new DataPoint(line.tumorRNAFracRefPercentile[5], false));
                    }

                    altHistograms[line.MultipleMutationsInThisGene].addValue(line.tumorRNAFracAltPercentile[5]);
                    if (!line.MultipleMutationsInThisGene)
                    {
                        singleMutationDataPoints.Add(new DataPoint(line.tumorRNAFracAltPercentile[5], true));
                    }

                    allHistograms[line.MultipleMutationsInThisGene].addValue(line.tumorRNAFracAllPercentile[5]);
                } // foreach line


                if ((refHistograms[true].count() >= 10 || refHistograms[false].count() >= 10) && singleMutationDataPoints.Count() > 0)
                {
                    bool enoughData, reversed;
                    double nFirstGroup, nSecondGroup, U, z;
                    var mw = ASETools.MannWhitney<DataPoint>.ComputeMannWhitney(singleMutationDataPoints, singleMutationDataPoints[0], x => x.alt, x => x.value, out enoughData, out reversed, out nFirstGroup, out nSecondGroup, out U, out z);

                    combinedHistogramsFile.Write(hugo_symbol);
                    if (enoughData)
                    {
                        results.Add(new HugoPAndMedians(hugo_symbol, mw, refHistograms[false].median(), altHistograms[false].median()));
                        combinedHistogramsFile.Write(" uncorrected p = " + mw);
                        nCandidates++;
                    }

                    combinedHistogramsFile.WriteLine("\twild-type one mutation (n = " + refHistograms[false].count() + ")\tmutant one mutation (n = " + altHistograms[false].count() + ")\tall one mutation (n = " + allHistograms[false].count() + ")\twild-type > 1 mutation (n = " +
                        refHistograms[true].count() + ")\tmutant > 1 mutation (n = " + altHistograms[true].count() + ")\tall > 1 mutation (n = " + allHistograms[true].count() + ")");

                    const int nHistograms = 6;
                    var histograms = new ASETools.Histogram[nHistograms];
                    histograms[0] = refHistograms[false];
                    histograms[1] = altHistograms[false];
                    histograms[2] = allHistograms[false];
                    histograms[3] = refHistograms[true];
                    histograms[4] = altHistograms[true];
                    histograms[5] = allHistograms[true];


                    var cdfs = new List<double>[nHistograms];
                    for (int i = 0; i < nHistograms; i++)
                    {
                        cdfs[i] = histograms[i].ComputeHistogram(0, 3.02, .02).Select(x => x.cdfValue).ToList();
                    }

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

                    //
                    // Now write out the R format lines to allow easily making R boxplots.  We pad the shorter ones out with "NA" so they're all equal length
                    // so that we can put them in an R data.frame.
                    //
                    int maxLength = refHistograms.Select(_ => _.Value.getValues().Count()).Max();   // The ref is always longest, since it includes some Nonsense Mediated Decay that's otherwise excluded.

                    var perGeneRawDataFilename = configuration.geneScatterGraphsDirectory + hugo_symbol + ASETools.raw_median_data_extension;

                    double biggestValue = 0;
                    for (int i = 0; i < nHistograms; i++)
                    {
                        if (histograms[i].count() > 0)
                        {
                            biggestValue = Math.Max(biggestValue, histograms[i].max());
                        }
                    }

                    combinedHistogramsFile.WriteLine(hugo_symbol + " <- read.table(\"" + perGeneRawDataFilename.Replace('\\', '/') + "\", sep=\"\\t\", header=TRUE)");    // switch from backslashes to forward, since that's what R wants
                    combinedHistogramsFile.WriteLine("boxplot(" + hugo_symbol + ", names=c(\"Ref 1\", \"Mut 1\", \"All 1\", \"Ref >1\", \"Mut >1\", \"All >1\"), ylim=c(0," + Math.Min(3, biggestValue) + 
                        "), ylab=\"Expression (multiple of median)\")");

                    combinedHistogramsFile.WriteLine(); // Space between genes

                    var perGeneRawDataFile = ASETools.CreateStreamWriterWithRetry(perGeneRawDataFilename);
                    perGeneRawDataFile.WriteLine("Ref 1 mutation\tAlt 1 mutation\tAll 1 mutation\tRef > 1 mutation\tAlt > 1 mutation\tAll > 1 mutation");
                    var ref_ = refHistograms[false].getValues().EnumerateWithCommas(maxLength, "").Split(',');
                    var alt = altHistograms[false].getValues().EnumerateWithCommas(maxLength, "").Split(',');
                    var all = allHistograms[false].getValues().EnumerateWithCommas(maxLength, "").Split(',');
                    var ref_many = refHistograms[true].getValues().EnumerateWithCommas(maxLength, "").Split(',');
                    var alt_many = altHistograms[true].getValues().EnumerateWithCommas(maxLength, "").Split(',');
                    var all_many = allHistograms[true].getValues().EnumerateWithCommas(maxLength, "").Split(',');

                    for (int i = 0; i < ref_.Count(); i++)
                    {
                        perGeneRawDataFile.WriteLine(ref_[i] + "\t" + alt[i] + "\t" + all[i] + "\t" + ref_many[i] + "\t" + alt_many[i] + "\t" + all_many[i]);
                    }
                    perGeneRawDataFile.Close();
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

            results.Sort((x, y) => x.p.CompareTo(y.p));
            for (int i = 0; i < results.Count(); i++)
            {
                if (results[i].p * nCandidates <= 0.01)
                {
                    Console.WriteLine(results[i].hugo_symbol + " has ref median of " + results[i].refMedian + " and alt median of " + results[i].altMedian + " corrected p value of " + results[i].p * nCandidates);
                }
            }

            Console.WriteLine("nNonsenseMediatedDecay = " + ASETools.nNonsenseMediatedDecay + ", nNoReadCounts = " + ASETools.nNoReadCounts + ", nBadReadCounts = " + ASETools.nBadReadCounts + ", nBadGene = " + ASETools.nBadGene + ", nNoCopyNumber = " + ASETools.nNoCopyNumber + 
                ", nBadCopyNumber = " + ASETools.nBadCopyNumber + ", nRepetitive = " + ASETools.nRepetitive + ", nCandidate = " + ASETools.nCandidate);
            Console.WriteLine("Total of " + nCandidates + " Mann-Whitney computations by this program");

            Console.WriteLine("Overall run time " + ASETools.ElapsedTimeInSeconds(commonData.timer));
        } // Main

        static int CompareScatterGraphLines(ASETools.GeneScatterGraphLine a, ASETools.GeneScatterGraphLine b, Dictionary<string, ASETools.ASEMapPerGeneLine> perGeneASEMap, ASETools.GeneMap geneMap, 
            Dictionary<string, List<ASETools.CopyNumberVariation>> copyNumberByCase, ASETools.ASERepetitiveRegionMap repetitiveRegionMap)
        {
            var aIsASECanddiate = a.isASECandidate(copyNumberByCase[a.case_id], configuration, perGeneASEMap, geneMap, repetitiveRegionMap);
            var bIsASECanddiate = b.isASECandidate(copyNumberByCase[b.case_id], configuration, perGeneASEMap, geneMap, repetitiveRegionMap);

            if (aIsASECanddiate != bIsASECanddiate)
            {
                return aIsASECanddiate ? -1 : 1;
            }

            if (a.tumorRNAFracAllPercentile == null || b.tumorRNAFracAllPercentile == null)
            {
                return 0;
            }

            if (a.MultipleMutationsInThisGene != b.MultipleMutationsInThisGene)
            {
                return a.MultipleMutationsInThisGene.CompareTo(b.MultipleMutationsInThisGene);
            }

            return a.tumorRNAFracAllPercentile[5].CompareTo(b.tumorRNAFracAllPercentile[5]);
        } // CompareScatterGraphLines
    }
}
