using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using System.IO;
using ExpressionLib;
using System.Diagnostics;
using System.Threading;

namespace ExpressionNearMutations
{
    class Program
    {

        class RegionalExpressionState
        {
            public int nRegionsIncluded = 0;
            public double minExpression = 100000;
            public double maxExpression = -100000;
            public double totalExpression = 0;

            public double minMeanExpression = 100000;
            public double maxMeanExpression = -1;
            public double totalMeanExpression = 0;

            public void AddExpression(double z, double mu)
            {
                nRegionsIncluded++;
                totalExpression += z;
                minExpression = Math.Min(minExpression, z);
                maxExpression = Math.Max(maxExpression, z);

                totalMeanExpression += mu;
                minMeanExpression = Math.Min(minMeanExpression, mu);
                maxMeanExpression = Math.Max(maxMeanExpression, mu);
            }
        }
        class Gene
        {
 
            public Gene(string hugo_symbol_, string chromosome_, int offset)
            {
                hugo_symbol = hugo_symbol_;
                chromosome = chromosome_;
                minOffset = offset;
                maxOffset = offset;
            }

            public string hugo_symbol;  // The gene name
            public string chromosome;
            public int minOffset;
            public int maxOffset;

            public bool inconsistent = false;
        }

        class GeneExpression 
        {
            static GeneExpression()
            {
                regionSizeByRegionSizeIndex[0] = 0;
                regionSizeByRegionSizeIndex[1] = 1000;
                for (int i = 2; i < nRegionSizes; i++)
                {
                    regionSizeByRegionSizeIndex[i] = regionSizeByRegionSizeIndex[i - 1] * 2;
                }

                comparer = StringComparer.OrdinalIgnoreCase;
            }

            public GeneExpression(Gene gene_) 
            {
                gene = gene_;

                for (int sizeIndex = 0; sizeIndex < nRegionSizes; sizeIndex++)
                {
                    regionalExpressionState[sizeIndex] = new RegionalExpressionState();
                }
            }

            public void AddRegionalExpression(int chromosomeOffset, double z, double mu)
            {
                int distance;
                if (chromosomeOffset >= gene.minOffset && chromosomeOffset <= gene.maxOffset)
                {
                    distance = 0;
                }
                else if (chromosomeOffset < gene.minOffset)
                {
                    distance = gene.minOffset - chromosomeOffset;
                }
                else
                {
                    distance = chromosomeOffset - gene.maxOffset;
                }


                if (0 == distance)
                {
                    regionalExpressionState[0].AddExpression(z, mu);
                }
                else
                {
                    for (int sizeIndex = nRegionSizes - 1; sizeIndex > 0; sizeIndex--)  // Don't do 0, so we exclude the gene from the surronding region
                    {
                        if (regionSizeByRegionSizeIndex[sizeIndex] < distance)
                        {
                            break;
                        }

                        regionalExpressionState[sizeIndex].AddExpression(z, mu);
                    }
                }
            }

            public static int CompareByGeneName(GeneExpression a, GeneExpression b)
            {
                return comparer.Compare(a.gene.hugo_symbol, b.gene.hugo_symbol);
            }

            public const int nRegionSizes = 20;    // Because we have 0 (in the gene), this range is 2^(20 - 2) * 1000 = 262 Mbases on either side, i.e., the entire chromosome
            public static readonly int[] regionSizeByRegionSizeIndex = new int[nRegionSizes];

            public RegionalExpressionState[] regionalExpressionState = new RegionalExpressionState[nRegionSizes]; // Dimension is log2(regionSize) - 1

            public Gene gene;
            public int mutationCount = 0;
            public static StringComparer comparer;
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
                        genesByName.Add(gene.hugo_symbol, gene);
                    }
                }
            }

            public int Count()
            {
                return genes.Count();
            }

            Dictionary<string, Gene> genes = new Dictionary<string, Gene>();
            public Dictionary<string, List<Gene>> genesByChromosome = new Dictionary<string, List<Gene>>();
            public Dictionary<string, Gene> genesByName = new Dictionary<string, Gene>();
        }

        static Dictionary<string, MutationMap> mutations;

        static void ProcessParticipants(List<string> participantsToProcess)
        {
            var timer = new Stopwatch();

            string participantId;

            while (true)
            {
                lock (participantsToProcess)
                {
                    if (participantsToProcess.Count() == 0)
                    {
                        return;
                    }

                    participantId = participantsToProcess[0];
                    participantsToProcess.RemoveAt(0);
                }

                timer.Reset();
                timer.Start();

                if (!experimentsByParticipant.ContainsKey(participantId))
                {
                    Console.WriteLine("Couldn't find experiment for participant ID " + participantId);
                    continue;
                }

                var experiment = experimentsByParticipant[participantId];

                if (experiment.TumorRNAAnalysis.regionalExpressionFileName == "")
                {
                    Console.WriteLine("Participant " + participantId + " doesn't have a regional expression file yet.");
                    continue;
                }

                var mutationsForThisReference = mutations[experiment.maf[0].ReferenceClass()];

                var geneExpressions = new Dictionary<string, GeneExpression>();    
                foreach (var maf in experiment.maf)
                {
                    if (maf.Variant_classification == "Silent")
                    {
                        continue;
                    }

                    if (!mutationsForThisReference.genesByName.ContainsKey(maf.Hugo_symbol)) {
                        //
                        // Probably an inconsistent gene.  Skip it.
                        //
                        continue;
                    }

                    if (!geneExpressions.ContainsKey(maf.Hugo_symbol))
                    {
                        geneExpressions.Add(maf.Hugo_symbol, new GeneExpression(mutationsForThisReference.genesByName[maf.Hugo_symbol]));
 
                    }
 
                    geneExpressions[maf.Hugo_symbol].mutationCount++;
                }

                var reader = new StreamReader(experiment.TumorRNAAnalysis.regionalExpressionFileName);

                var headerLine = reader.ReadLine();
                if (null == headerLine)
                {
                    Console.WriteLine("Empty regional expression file " + experiment.TumorRNAAnalysis.regionalExpressionFileName);
                    continue;
                }

                if (headerLine.Count() < 23)
                {
                    Console.WriteLine("Truncated regional expression header line in file '" + experiment.TumorRNAAnalysis + "', line: ", headerLine);
                    continue;
                }

                if (headerLine.Substring(0, 20) != "RegionalExpression v")
                {
                    Console.WriteLine("Corrupt header line in file '" + experiment.TumorRNAAnalysis.regionalExpressionFileName + "', line: " + headerLine);
                    continue;
                }

                if (headerLine.Substring(20, 1) != "3")
                {
                    Console.WriteLine("Unsupported version in file '" + experiment.TumorRNAAnalysis.regionalExpressionFileName + "', header line: " + headerLine);
                    continue;
                }

                var line = reader.ReadLine();   // The NumContigs line, which we just ignore
                line = reader.ReadLine();   // The column header line, which we just ignore
                if (null == line)
                {
                    Console.WriteLine("Truncated file '" + experiment.TumorRNAAnalysis.regionalExpressionFileName + "' ends after header line.");
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
                    if (fields.Count() != 13)
                    {
                        Console.WriteLine("Badly formatted data line in file '" + experiment.TumorRNAAnalysis.regionalExpressionFileName + "', line " + lineNumber + ": " + line);
                        break;
                    }

                    string chromosome;
                    int offset;
                    double z;
                    double mu;

                    try
                    {
                        chromosome = fields[0];
                        offset = Convert.ToInt32(fields[1]);
                        z = Convert.ToDouble(fields[11]);
                        mu = Convert.ToDouble(fields[12]);

                        int nBasesExpressedWithBaselineExpression = Convert.ToInt32(fields[3]);
                        int nBasesUnexpressedWithBaselineExpression = Convert.ToInt32(fields[7]);

                        if (0 == nBasesExpressedWithBaselineExpression && 0 == nBasesUnexpressedWithBaselineExpression)
                        {
                            //
                            // No baseline expression for this region, skip it.
                            //
                            continue;
                        }
                    }
                    catch (FormatException)
                    {
                        Console.WriteLine("Format exception parsing data line in file '" + experiment.TumorRNAAnalysis.regionalExpressionFileName + "', line " + lineNumber + ": " + line);
                        break;
                    }

                    if (mutationsForThisReference.genesByChromosome.ContainsKey(chromosome))
                    {
                        foreach (var gene in mutationsForThisReference.genesByChromosome[chromosome])
                        {
                            if (!geneExpressions.ContainsKey(gene.hugo_symbol))
                            {
                                geneExpressions.Add(gene.hugo_symbol, new GeneExpression(gene));
                            }

                            geneExpressions[gene.hugo_symbol].AddRegionalExpression(offset, z, mu);
                        }
                    }
                }

                if (!seenDone)
                {
                    Console.WriteLine("Truncated regional expression file " + experiment.TumorRNAAnalysis.regionalExpressionFileName);
                    continue;
                }

                //
                // Write the output file.
                //
                int indexOfLastSlash = experiment.TumorRNAAnalysis.regionalExpressionFileName.LastIndexOf('\\');
                if (-1 == indexOfLastSlash)
                {
                    Console.WriteLine("Couldn't find a backslash in regional expression pathname, which is supposed to be absolute: " + experiment.TumorRNAAnalysis.regionalExpressionFileName);
                    continue;
                }

                string directory = experiment.TumorRNAAnalysis.regionalExpressionFileName.Substring(0, indexOfLastSlash + 1);  // Includes trailing backslash
                var outputFilename = directory + experiment.TumorRNAAnalysis.analysis_id + ".gene_expression.txt";

                var outputFile = new StreamWriter(outputFilename);

                outputFile.WriteLine("ExpressionNearMutations v2.0 " + participantId);
                outputFile.Write("Gene name\tnon-silent mutation count");
                for (int sizeIndex = 0; sizeIndex < GeneExpression.nRegionSizes; sizeIndex++)
                {
                    outputFile.Write("\t" + GeneExpression.regionSizeByRegionSizeIndex[sizeIndex] + "(z)");
                }

                for (int sizeIndex = 0; sizeIndex < GeneExpression.nRegionSizes; sizeIndex++)
                {
                    outputFile.Write("\t" + GeneExpression.regionSizeByRegionSizeIndex[sizeIndex] + "(mu)");
                }

                outputFile.WriteLine();

                var allExpressions = new List<GeneExpression>();
                foreach (var expressionEntry in geneExpressions)
                {
                    allExpressions.Add(expressionEntry.Value);
                }

                allExpressions.Sort(GeneExpression.CompareByGeneName);

                for (int i = 0; i < allExpressions.Count(); i++)
                {
                    outputFile.Write(allExpressions[i].gene.hugo_symbol + "\t" + allExpressions[i].mutationCount);

                    for (int sizeIndex = 0; sizeIndex < GeneExpression.nRegionSizes; sizeIndex++)
                    {
                        if (allExpressions[i].regionalExpressionState[sizeIndex].nRegionsIncluded != 0)
                        {
                            outputFile.Write("\t" + allExpressions[i].regionalExpressionState[sizeIndex].totalExpression / allExpressions[i].regionalExpressionState[sizeIndex].nRegionsIncluded);
                        }
                        else
                        {
                            outputFile.Write("\t*");
                        }
                    }

                    for (int sizeIndex = 0; sizeIndex < GeneExpression.nRegionSizes; sizeIndex++)
                    {
                        if (allExpressions[i].regionalExpressionState[sizeIndex].nRegionsIncluded != 0)
                        {
                            outputFile.Write("\t" + allExpressions[i].regionalExpressionState[sizeIndex].totalMeanExpression / allExpressions[i].regionalExpressionState[sizeIndex].nRegionsIncluded);
                        }
                        else
                        {
                            outputFile.Write("\t*");
                        }
                    }
                    outputFile.WriteLine();
                }

                outputFile.WriteLine("**done**");
                outputFile.Close();

                timer.Stop();
                lock (participantsToProcess)
                {
                    var nRemaining = participantsToProcess.Count();
                    Console.WriteLine("Processed participant " + participantId + " in " + (timer.ElapsedMilliseconds + 500) / 1000 + "s.  " + nRemaining + " remain" + ((1 == nRemaining) ? "s" : "") + " queued.");
                }
            }
        }

        static List<ExpressionTools.Experiment> experiments;
        static Dictionary<string, ExpressionTools.Participant> participants;
        static Dictionary<string, ExpressionTools.Experiment> experimentsByParticipant = new Dictionary<string, ExpressionTools.Experiment>();
 
        static void Main(string[] args)
        {
            if (args.Count() == 0)
            {
                Console.WriteLine("usage: ExpressionNearMutations <participantIdsToProcess>");
                return;
            }

            Stopwatch timer = new Stopwatch();
            timer.Start();

            var excludedAnalyses = ExpressionTools.LoadExcludedAnalyses(@"\\gcr\scratch\b99\bolosky\excluded_analyses.txt");

            var tcgaRecords = ExpressionTools.LoadTCGARecords(null /* stored BAMs*/, excludedAnalyses, @"\\gcr\scratch\b99\bolosky\tcga-all.xml");
            ExpressionTools.LoadTCGARecordsForLocalRealigns(tcgaRecords, null, @"\\gcr\scratch\b99\bolosky\realigns.txt");
            ExpressionTools.LoadTCGAAdditionalMetadata(tcgaRecords, @"\\gcr\scratch\b99\bolosky\tcgaAdditionalMetadata.txt");

            var sampleToParticipantIDMap = ExpressionTools.CreateSampleToParticipantIDMap(tcgaRecords);

            Dictionary<string, ExpressionTools.Sample> allSamples;
            participants = ExpressionTools.BuildParticipantData(tcgaRecords, out allSamples, @"\\gcr\scratch\b99\bolosky\clinical");

            ExpressionTools.AddAllMAFFilesToParticipants(participants, sampleToParticipantIDMap, @"\\gcr\scratch\b99\bolosky\mafs");

            experiments = ExpressionTools.LoadExperimentsFromFile(@"\\gcr\scratch\b99\bolosky\experiments.txt", participants, tcgaRecords);

            timer.Stop();
            Console.WriteLine("Loaded " + experiments.Count() + " experiments with MAFs in " + (timer.ElapsedMilliseconds + 500) / 1000 + "s");

            //
            // Now build the map of mutations by gene.
            //

            timer.Reset();
            timer.Start();

            mutations = new Dictionary<string, MutationMap>();  // Maps reference type (hg18 or hg19) to MutationMap
            mutations.Add("hg18", new MutationMap());
            mutations.Add("hg19", new MutationMap());

            var badHugoSymbols = new List<string>();    // These are corruped by Excel.  They're in the MAFs that I downloaded, and I'm removing them by hand.
            badHugoSymbols.Add("1-Mar");
            badHugoSymbols.Add("1-Dec");
            badHugoSymbols.Add("1-Sep");
            badHugoSymbols.Add("10-Mar");
            badHugoSymbols.Add("11-Mar");
            badHugoSymbols.Add("12-Sep");
            badHugoSymbols.Add("14-Sep");
            badHugoSymbols.Add("1SEPT4");
            badHugoSymbols.Add("2-Mar");
            badHugoSymbols.Add("2-Sep");
            badHugoSymbols.Add("3-Mar");
            badHugoSymbols.Add("3-Sep");
            badHugoSymbols.Add("4-Mar");
            badHugoSymbols.Add("4-Sep");
            badHugoSymbols.Add("5-Sep");
            badHugoSymbols.Add("6-Mar");
            badHugoSymbols.Add("6-Sep");
            badHugoSymbols.Add("7-Mar");
            badHugoSymbols.Add("7-Sep");
            badHugoSymbols.Add("8-Mar");
            badHugoSymbols.Add("9-Sep");

            int nMutations = 0;
            foreach (var experiment in experiments)
            {
                experimentsByParticipant.Add(experiment.participant.participantId, experiment);

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

                    if (badHugoSymbols.Contains(maf.Hugo_symbol))
                    {
                        Console.WriteLine("Bad hugo symbol " + maf.Hugo_symbol);
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

            var participantsToProcess = new List<string>();
            foreach (var arg in args)
            {
                participantsToProcess.Add(arg);
            }

            //
            // Process the runs in parallel
            //
            int totalNumberOfExperiments = experiments.Count();
            timer.Reset();
            timer.Start();

            var threads = new List<Thread>();
            for (int i = 0; i < Environment.ProcessorCount; i++)
            {
                threads.Add(new Thread(() => ProcessParticipants(participantsToProcess)));
            }

            threads.ForEach(t => t.Start());
            threads.ForEach(t => t.Join());

            timer.Stop();
            Console.WriteLine("Processed " + args.Count() + " experiments in " + (timer.ElapsedMilliseconds + 500) / 1000 + " seconds");
        }
    }
}
