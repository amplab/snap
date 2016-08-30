using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using System.IO;
using ExpressionLib;
using System.Diagnostics;

namespace ExpressionNearMutations
{
    class Program
    {

        class RegionalExpressionState
        {
            public int nTumors = 0;
            public double minExpression = 100000;
            public double maxExpression = -100000;
            public double totalExpression = 0;
        }
        class Gene
        {
            static Gene()
            {
                regionSizeByRegionSizeIndex[0] = 0;
                regionSizeByRegionSizeIndex[1] = 1000;
                for (int i = 2; i < nRegionSizes; i++)
                {
                    regionSizeByRegionSizeIndex[i] = regionSizeByRegionSizeIndex[i - 1] * 2;
                }
            }
            public Gene(string hugo_symbol_, string chromosome_, int offset)
            {
                hugo_symbol = hugo_symbol_;
                chromosome = chromosome_;
                minOffset = offset;
                maxOffset = offset;

                for (int nMutations = 0; nMutations < 3; nMutations++)
                {
                    for (int size = 0; size < nRegionSizes; size++)
                    {
                        regionalExpressionState[nMutations, size] = new RegionalExpressionState();
                    }
                }
            }

            public void AddRegionalExpression(int nMutations, int chromosomeOffset, double z)
            {
                int distance;
                if (chromosomeOffset >= minOffset && chromosomeOffset <= maxOffset)
                {
                    distance = 0;
                }
                else if (chromosomeOffset < minOffset)
                {
                    distance = minOffset - chromosomeOffset;
                }
                else
                {
                    distance = chromosomeOffset - maxOffset;
                }

                for (int i = nRegionSizes - 1; i >= 0; i--)
                {
                    RegionalExpressionState state = regionalExpressionState[nMutations, i];

                    if (regionSizeByRegionSizeIndex[i] < distance)
                    {
                        break;
                    }

                    state.nTumors++;
                    state.totalExpression += z;
                    state.minExpression = Math.Min(state.minExpression, z);
                    state.maxExpression = Math.Max(state.maxExpression, z);
                }
            }

            public string hugo_symbol;  // The gene name
            public string chromosome;
            public int minOffset;
            public int maxOffset;

            public bool inconsistent = false;

            public const int nRegionSizes = 20;    // Because we have 0 (in the gene), this range is 2^(20 - 2) * 1000 = 262 Mbases on either side, i.e., the entire chromosome
            public static readonly int[] regionSizeByRegionSizeIndex = new int[nRegionSizes];

            public RegionalExpressionState[,] regionalExpressionState = new RegionalExpressionState[3, nRegionSizes]; // First dimension is nMutations (0, 1, more than 1) and second is log2(regionSize) - 1


        }

        class MutationMap
        {
            public static int largestAllowedGene = 2100000; // At little bigger than DMD, the largest human gene
            public MutationMap() {}

            public void AddMutation(string hugo_symbol, string chromosome, int offset)
            {
                if (!genes.ContainsKey(hugo_symbol))
                {
                    genes.Add(hugo_symbol, new Gene(hugo_symbol, chromosome, offset));
                }

                var gene = genes[hugo_symbol];

                if (!genes[hugo_symbol].inconsistent) {
                    if (gene.chromosome != chromosome) {
                        Console.WriteLine("Gene " + hugo_symbol + " occurs on (at least) two different chromosomes: " + chromosome + " and " + gene.chromosome + ".  Ignoring gene.");
                        gene.inconsistent = true;
                    }  else {
                        gene.minOffset = Math.Min(gene.minOffset, offset);
                        gene.maxOffset = Math.Max(gene.maxOffset, offset);
                        if (gene.maxOffset - gene.minOffset > largestAllowedGene) {
                            Console.WriteLine("Gene " + hugo_symbol + " has too big a range of mutation offsets: " + chromosome + ":" + gene.minOffset + "-" + gene.maxOffset);
                            gene.inconsistent = true;
                        }
                    }
                }
            }

            public void DoneAddingMutations()
            {
                foreach (var geneEntry in genes) {
                    var gene = geneEntry.Value;

                    if (!gene.inconsistent)
                    {
                        if (!genesByChromosome.ContainsKey(gene.chromosome))
                        {
                            genesByChromosome.Add(gene.chromosome, new List<Gene>());
                        }
                        genesByChromosome[gene.chromosome].Add(gene);
                    }
                }
            }

            public int Count()
            {
                return genes.Count();
            }

            Dictionary<string, Gene> genes = new Dictionary<string, Gene>();
            public Dictionary<string, List<Gene>> genesByChromosome = new Dictionary<string, List<Gene>>();
        }

        static void Main(string[] args)
        {

            Stopwatch timer = new Stopwatch();
            timer.Start();

            var excludedAnalyses = ExpressionTools.LoadExcludedAnalyses();

            var tcgaRecords = ExpressionTools.LoadTCGARecords(null /* stored BAMs*/, excludedAnalyses, @"f:\sequence\Reads\tcga-all.xml");
            ExpressionTools.LoadTCGARecordsForLocalRealigns(tcgaRecords, null, @"f:\sequence\reads\tcga\realigns.txt");
            ExpressionTools.LoadTCGAAdditionalMetadata(tcgaRecords);

            var sampleToParticipantIDMap = ExpressionTools.CreateSampleToParticipantIDMap(tcgaRecords);

            Dictionary<string, ExpressionTools.Sample> allSamples;
            var participants = ExpressionTools.BuildParticipantData(tcgaRecords, out allSamples);

            ExpressionTools.AddAllMAFFilesToParticipants(participants, sampleToParticipantIDMap);

            var experiments = ExpressionTools.LoadExperimentsFromFile(@"f:\temp\expression\experiments.txt", participants, tcgaRecords);

            timer.Stop();
            Console.WriteLine("Loaded " + experiments.Count() + " experiments with MAFs in " + (timer.ElapsedMilliseconds + 500) / 1000 + "s");

            //
            // Now build the map of mutations by gene.
            //

            timer.Reset();
            timer.Start();
            var mutations = new Dictionary<string, MutationMap>();  // Maps reference type (hg18 or hg19) to MutationMap
            mutations.Add("hg18", new MutationMap());
            mutations.Add("hg19", new MutationMap());

            int nMutations = 0;
            foreach (var experiment in experiments)
            {
                foreach (ExpressionTools.MAFRecord maf in experiment.maf)
                {
                    nMutations++;

                    string chromosome;
                    if (maf.Chrom.Count() > 3 && maf.Chrom.Substring(0, 3) == "chr")
                    {
                        chromosome = maf.Chrom.Substring(3);
                    }
                    else
                    {
                        chromosome = maf.Chrom;
                    }

                    mutations[maf.ReferenceClass()].AddMutation(maf.Hugo_symbol, chromosome, maf.Start_position);
                }
            }

            foreach (var referenceClass in mutations)
            {
                referenceClass.Value.DoneAddingMutations();
            }

            timer.Stop();
            Console.WriteLine("Loaded " + nMutations + " mutations in " + mutations["hg19"].Count() + " genes in " + (timer.ElapsedMilliseconds + 500) / 1000 + "s.");

            //
            // Run through each regional expression file and look at the regions around each mutation, categorizing it by whether the 
            // tumor has zero, one or more than one mutation in that gene.
            //
            int nWithoutRegionalExpression = 0;

            timer.Start();
            int nExperimentsProcessed = 0;
            long elapsedMillisecondsAtLastPrint = 0;
            foreach (var experiment in experiments)
            {
                if (experiment.TumorRNAAnalysis.regionalExpressionFileName == "")
                {
                    nWithoutRegionalExpression++;
                    continue;
                }

                var mutationsForThisReference = mutations[experiment.maf[0].ReferenceClass()];

                var genesWithMutations = new Dictionary<string, int>();    // hugo_symbol->{1,2} for exactly one and more than one mutation.
                foreach (var maf in experiment.maf)
                {
                    if (maf.Variant_classification == "Silent") {
                        continue;
                    }

                    if (genesWithMutations.ContainsKey(maf.Hugo_symbol))
                    {
                        genesWithMutations[maf.Hugo_symbol] = 2;
                    }
                    else
                    {
                        genesWithMutations.Add(maf.Hugo_symbol, 1);
                    }
                }

                var reader = new StreamReader(experiment.TumorRNAAnalysis.regionalExpressionFileName);

                var headerLine = reader.ReadLine();
                if (null == headerLine)
                {
                    Console.WriteLine("Empty regional expression file " + experiment.TumorRNAAnalysis.regionalExpressionFileName);
                    nWithoutRegionalExpression++;
                    continue;
                }

                if (headerLine.Count() < 23)
                {
                    Console.WriteLine("Truncated regional expression header line in file '" + experiment.TumorRNAAnalysis + "', line: ", headerLine);
                    nWithoutRegionalExpression++;
                    continue;
                }

                if (headerLine.Substring(0, 20) != "RegionalExpression v")
                {
                    Console.WriteLine("Corrupt header line in file '" + experiment.TumorRNAAnalysis.regionalExpressionFileName + "', line: " + headerLine);
                    nWithoutRegionalExpression++;
                    continue;
                }

                if (headerLine.Substring(20, 1) != "2")
                {
                    Console.WriteLine("Unsupported version in file '" + experiment.TumorRNAAnalysis.regionalExpressionFileName + "', header line: " + headerLine);
                    nWithoutRegionalExpression++;
                    continue;
                }

                var line = reader.ReadLine();   // The NumContigs line, which we just ignore
                line = reader.ReadLine();   // The column header line, which we just ignore
                if (null == line)
                {
                    Console.WriteLine("Truncated file '" + experiment.TumorRNAAnalysis.regionalExpressionFileName + "' ends after header line.");
                    nWithoutRegionalExpression++;
                    continue;
                }

                bool seenDone = false;
                int lineNumber = 3;
                while (null != (line = reader.ReadLine()))
                {
                    lineNumber++;

                    if (seenDone)
                    {
                        Console.WriteLine("Saw data after **done** in file " + experiment.TumorRNAAnalysis.regionalExpressionFileName + "', line " + lineNumber + ": " + line);
                        break;
                    }

                    if (line == "**done**")
                    {
                        seenDone = true;
                        continue;
                    }

                    var fields = line.Split('\t');
                    if (fields.Count() != 12)
                    {
                        Console.WriteLine("Badly formatted data line in file '" + experiment.TumorRNAAnalysis.regionalExpressionFileName + "', line " + lineNumber + ": " + line);
                        break;
                    }

                    string chromosome;
                    int offset;
                    double z;

                    try
                    {
                        chromosome = fields[0];
                        offset = Convert.ToInt32(fields[1]);
                        z = Convert.ToDouble(fields[11]);
                    }
                    catch (FormatException)
                    {
                        Console.WriteLine("Format exception parsing data line in file '" + experiment.TumorRNAAnalysis.regionalExpressionFileName + "', line " + lineNumber + ": " + line);
                        nWithoutRegionalExpression++;
                        break;
                    }

                    if (mutationsForThisReference.genesByChromosome.ContainsKey(chromosome))
                    {
                        foreach (var gene in mutationsForThisReference.genesByChromosome[chromosome])
                        {
                            if (!genesWithMutations.ContainsKey(gene.hugo_symbol))
                            {
                                genesWithMutations.Add(gene.hugo_symbol, 0);
                            }

                            gene.AddRegionalExpression(genesWithMutations[gene.hugo_symbol], offset, z);
                        }
                    }
                }

                if (!seenDone)
                {
                    Console.WriteLine("Truncated regional expression file " + experiment.TumorRNAAnalysis.regionalExpressionFileName);
                }

                nExperimentsProcessed++;
                if (timer.ElapsedMilliseconds - elapsedMillisecondsAtLastPrint > 60000)
                {
                    elapsedMillisecondsAtLastPrint = timer.ElapsedMilliseconds;
                    Console.WriteLine("Processed " + nExperimentsProcessed + " experiments in " + (elapsedMillisecondsAtLastPrint + 30000) / 60000 + "m.  " + elapsedMillisecondsAtLastPrint / 1000 / nExperimentsProcessed + " s/experiment");
                }
             } // foreach experiment

            //
            // Now dump out the results.
            //
            var outputFile = new StreamWriter(@"f:\temp\expression\regional_expression_by_gene.txt");
            outputFile.WriteLine("ExpressionNearMutations v1.0 ");
            outputFile.Write("Reference class\tgene\tchromosome\tminOffset\tmaxOffset\tnMutations");

            var fullDumpFile = new StreamWriter(@"f:\temp\expression\reginal_expression_all.txt");
            fullDumpFile.WriteLine("ExpressionNearMutations v1.0");
            fullDumpFile.Write("Reference class\tparticipantId\tgene\tchromosome\tminOffset\tmaxOffset\tnMutations");


            for (int regionSizeIndex = 0; regionSizeIndex < Gene.nRegionSizes; regionSizeIndex++)
            {
                outputFile.Write("\t" + Gene.regionSizeByRegionSizeIndex[regionSizeIndex]);
                fullDumpFile.Write("\t" + Gene.regionSizeByRegionSizeIndex[regionSizeIndex]);
            }
            outputFile.WriteLine("\tnTumors");
            fullDumpFile.WriteLine();

            foreach (var referenceClass in mutations)
            {
                foreach (var chromosomeGenes in referenceClass.Value.genesByChromosome)
                {
                    foreach (var gene in chromosomeGenes.Value)
                    {
                        for (nMutations = 0; nMutations < 3; nMutations++)
                        {
                            outputFile.Write(referenceClass.Key + "\t\"" + gene.hugo_symbol + "\"\t" + gene.chromosome + "\t" + gene.minOffset + "\t" + gene.maxOffset + "\t");
                            if (2 == nMutations)
                            {
                                outputFile.Write(">1");
                            }
                            else
                            {
                                outputFile.Write(nMutations);
                            }

                            for (int regionSizeIndex = 0; regionSizeIndex < Gene.nRegionSizes; regionSizeIndex++)
                            {
                                if (gene.regionalExpressionState[nMutations, regionSizeIndex].nTumors == 0)
                                {
                                    outputFile.Write("\t*");
                                }
                                else
                                {
                                    outputFile.Write("\t" + gene.regionalExpressionState[nMutations, regionSizeIndex].totalExpression / gene.regionalExpressionState[nMutations, regionSizeIndex].nTumors);
                                }
                            } // for each region size
                            outputFile.WriteLine();
                        } // for each nMutations
                    } // foreach gene
                } // foreach chromosome
            } // foreach reference class

            outputFile.WriteLine("**done**");
            outputFile.Close();

            fullDumpFile.WriteLine("**done**");
            fullDumpFile.Close();
        }
    }
}
