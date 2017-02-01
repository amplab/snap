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


        class GeneExpression 
        {
            static GeneExpression() // This is a static initializer that runs once at program start time, it's not a constructor.
            {
                regionSizeByRegionSizeIndex[0] = 0;
                regionSizeByRegionSizeIndex[1] = 1000;
                for (int i = 2; i < nRegionSizes; i++)
                {
                    regionSizeByRegionSizeIndex[i] = regionSizeByRegionSizeIndex[i - 1] * 2;
                }

                comparer = StringComparer.OrdinalIgnoreCase;
            }

            public GeneExpression(ExpressionTools.GeneLocationInfo gene_) 
            {
                geneLocationInfo = gene_;

                for (int sizeIndex = 0; sizeIndex < nRegionSizes; sizeIndex++)
                {
                    regionalExpressionState[sizeIndex] = new RegionalExpressionState();
                    exclusiveRegionalExpressionState[sizeIndex] = new RegionalExpressionState();
                }
            }

            public void AddRegionalExpression(int chromosomeOffset, double z, double mu)
            {
                int distance;
                if (chromosomeOffset >= geneLocationInfo.minLocus && chromosomeOffset <= geneLocationInfo.maxLocus)
                {
                    distance = 0;
                }
                else if (chromosomeOffset < geneLocationInfo.minLocus)
                {
                    distance = geneLocationInfo.minLocus - chromosomeOffset;
                }
                else
                {
                    distance = chromosomeOffset - geneLocationInfo.maxLocus;
                }

                for (int sizeIndex = nRegionSizes - 1; sizeIndex >= 0; sizeIndex--)
                {
                    if (regionSizeByRegionSizeIndex[sizeIndex] < distance)
                    {
                        if (sizeIndex != nRegionSizes - 1)
                        {
                            exclusiveRegionalExpressionState[sizeIndex + 1].AddExpression(z, mu);
                            break;
                        }
                    }
                    regionalExpressionState[sizeIndex].AddExpression(z, mu);
                }

                if (0 == distance)  // Have to special case this, since exclusive gets added when we're one smaller, and there is nothing smaller than sizeIndex 0.
                {
                    exclusiveRegionalExpressionState[0].AddExpression(z, mu);
                }
            }

 
            public static int CompareByGeneName(GeneExpression a, GeneExpression b)
            {
                return comparer.Compare(a.geneLocationInfo.hugoSymbol, b.geneLocationInfo.hugoSymbol);
            }

            public const int nRegionSizes = 20;    // Because we have 0 (in the gene), this range is 2^(20 - 2) * 1000 = 262 Mbases on either side, i.e., the entire chromosome
            public static readonly int[] regionSizeByRegionSizeIndex = new int[nRegionSizes];

            public RegionalExpressionState[] regionalExpressionState = new RegionalExpressionState[nRegionSizes]; // Dimension is log2(regionSize) - 1
            public RegionalExpressionState[] exclusiveRegionalExpressionState = new RegionalExpressionState[nRegionSizes];  // Expression in this region but not closer, so from log2(regionSize - 1) to log2(regionSize) - 1.  The zero element is the same as regionalExpressionState

            public ExpressionTools.GeneLocationInfo geneLocationInfo;
            public int mutationCount = 0;
            public static StringComparer comparer;
        }



        //static Dictionary<string, ExpressionTools.MutationMap> mutations;

        static Dictionary<string, ExpressionTools.GeneLocationsByNameAndChromosome> geneLocationInformation;

        static void ProcessParticipants(List<string> participantsToProcess, bool forAlleleSpecificExpression, int minExamplesPerRegion)
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

                var inputFilename = forAlleleSpecificExpression ? experiment.NormalDNAAnalysis.annotatedSelectedVariantsFileName : experiment.TumorRNAAnalysis.regionalExpressionFileName;

                if (inputFilename == "")
                {
                    Console.WriteLine("Participant " + participantId + " doesn't have an input file yet.");
                    continue;
                }

                //var mutationsForThisReference = mutations[experiment.maf[0].ReferenceClass()];
                var geneLocationsForThisReference = geneLocationInformation[experiment.maf[0].ReferenceClass()];

                var geneExpressions = new Dictionary<string, GeneExpression>();    
                foreach (var maf in experiment.maf)
                {
                    if (maf.Variant_classification == "Silent")
                    {
                        continue;
                    }

                    if (!geneLocationsForThisReference.genesByName.ContainsKey(maf.Hugo_symbol))
                    {
                        //
                        // Probably an inconsistent gene.  Skip it.
                        //
                        continue;
                    }

                    if (!geneExpressions.ContainsKey(maf.Hugo_symbol))
                    {
                        geneExpressions.Add(maf.Hugo_symbol, new GeneExpression(geneLocationsForThisReference.genesByName[maf.Hugo_symbol]));
                    }
 
                    geneExpressions[maf.Hugo_symbol].mutationCount++;
                }

                var reader = ExpressionTools.CreateStreamReaderWithRetry(inputFilename);

                var headerLine = reader.ReadLine();
                if (null == headerLine)
                {
                    Console.WriteLine("Empty input file " + inputFilename);
                    continue;
                }

                string line;
                int lineNumber = 1;
                if (!forAlleleSpecificExpression)
                {
                    if (headerLine.Substring(0, 20) != "RegionalExpression v")
                    {
                        Console.WriteLine("Corrupt header line in file '" + inputFilename + "', line: " + headerLine);
                        continue;
                    }

                    if (headerLine.Substring(20, 2) != "3.")
                    {
                        Console.WriteLine("Unsupported version in file '" + inputFilename + "', header line: " + headerLine);
                        continue;
                    }
                    line = reader.ReadLine();   // The NumContigs line, which we just ignore
                    line = reader.ReadLine();   // The column header line, which we just ignore

                    if (null == line)
                    {
                        Console.WriteLine("Truncated file '" + inputFilename + "' ends after header line.");
                        continue;
                    }

                    lineNumber = 3;
                }

                var wholeAutosomeRegionalExpression = new RegionalExpressionState();
                var allButThisChromosomeAutosomalRegionalExpressionState = new Dictionary<string, RegionalExpressionState>();   // "This chromosome" is the dictionary key
                var perChromosomeRegionalExpressionState = new RegionalExpressionState[ExpressionTools.nHumanNuclearChromosomes];

                for (int whichChromosome = 0; whichChromosome < ExpressionTools.nHumanNuclearChromosomes; whichChromosome++)
                {
                    perChromosomeRegionalExpressionState[whichChromosome] = new RegionalExpressionState();
                }

                foreach (var geneEntry in geneLocationsForThisReference.genesByName)
                {
                    var chromosome = geneEntry.Value.chromosome;
                    if (!allButThisChromosomeAutosomalRegionalExpressionState.ContainsKey(chromosome))
                    {
                        allButThisChromosomeAutosomalRegionalExpressionState.Add(chromosome, new RegionalExpressionState());
                    }
                }
               
                bool seenDone = false;
                while (null != (line = reader.ReadLine()))
                {
                    lineNumber++;

                    if (seenDone)
                    {
                        Console.WriteLine("Saw data after **done** in file " + inputFilename + "', line " + lineNumber + ": " + line);
                        break;
                    }

                    if (line == "**done**")
                    {
                        seenDone = true;
                        continue;
                    }

 
                    string chromosome;
                    int offset;

                    // For allele-specific expression
                    double nMatchingReferenceDNA = 0;
                    double nMatchingVariantDNA = 0;
                    double nMatchingReferenceRNA = 0;
                    double nMatchingVariantRNA = 0;

                    // for regional expression
                    double z = 0;
                    double mu = 0;

                    try {
                        if (forAlleleSpecificExpression) {
                            var alleleData = ExpressionTools.AnnotatedSelectedVariantLine.fromText(line);

                            if (null == alleleData)
                            {
                                Console.WriteLine("Error parsing input line " + lineNumber + " in file " + inputFilename);
                                break;
                            }

                            chromosome = alleleData.contig;
                            offset = alleleData.loc;
                            nMatchingReferenceDNA = alleleData.nMatchingReferenceDNA;
                            nMatchingVariantRNA = alleleData.nMatchingVariantDNA;
                            nMatchingReferenceRNA = alleleData.nMatchingReferenceRNA;
                            nMatchingVariantRNA = alleleData.nMatchingVariantRNA;

                            if (!geneLocationsForThisReference.genesByChromosome.ContainsKey(chromosome))
                            {
                                //
                                // Try reversing the "chr" state of the chromosome.
                                //

                                if (chromosome.Count() > 3 && chromosome.Substring(0, 3) == "chr")
                                {
                                    chromosome = chromosome.Substring(3);
                                }
                                else
                                {
                                    chromosome = "chr" + chromosome;
                                }
                            }
                        } else {

                            var fields = line.Split('\t');
                            if (fields.Count() != 13)
                            {
                                Console.WriteLine("Badly formatted data line in file '" + inputFilename + "', line " + lineNumber + ": " + line);
                                break;
                            }

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
                    }
                    catch (FormatException)
                    {
                        Console.WriteLine("Format exception parsing data line in file '" + inputFilename + "', line " + lineNumber + ": " + line);
                        break;
                    }

                    if (geneLocationsForThisReference.genesByChromosome.ContainsKey(chromosome) && 
                        (!forAlleleSpecificExpression ||                                    // We only keep samples for allele specific expression if they meet certain criteria, to wit:
                            nMatchingReferenceDNA + nMatchingVariantDNA >= 10 &&            // We have at least 10 DNA reads
                            nMatchingReferenceRNA + nMatchingVariantRNA >= 10 &&            // We have at least 10 RNA reads
                            nMatchingReferenceDNA * 3 >= nMatchingVariantDNA * 2 &&         // It's not more than 2/3 variant DNA
                            nMatchingVariantDNA * 3 >= nMatchingReferenceDNA * 2))          // It's not more than 2/3 reference DNA
                    {
                        if (forAlleleSpecificExpression)
                        {
                            double rnaFraction = nMatchingVariantRNA / (nMatchingReferenceRNA + nMatchingVariantRNA);

                            //
                            // Now convert to the amount of allele-specific expression.  50% is no ASE, while 0 or 100% is 100% ASE.
                            //
                            z = Math.Abs(rnaFraction * 2.0 - 1.0); // Not really z, really alleleSpecificExpression
                            mu = 0;
                        }

                        if (ExpressionTools.isChromosomeAutosomal(chromosome))
                        {
                            wholeAutosomeRegionalExpression.AddExpression(z, mu);

                            foreach (var entry in allButThisChromosomeAutosomalRegionalExpressionState) 
                            {
                                if (entry.Key != chromosome) {
                                    entry.Value.AddExpression(z, mu);
                                }
                            }
                        }

                        int chromosomeId = ExpressionTools.ChromosomeNameToIndex(chromosome);
                        if (chromosomeId != -1)
                        {
                            perChromosomeRegionalExpressionState[chromosomeId].AddExpression(z, mu);
                        }

                        foreach (var geneLocation in geneLocationsForThisReference.genesByChromosome[chromosome])
                        {
                            if (!geneExpressions.ContainsKey(geneLocation.hugoSymbol))
                            {
                                geneExpressions.Add(geneLocation.hugoSymbol, new GeneExpression(geneLocation));
                            }

                            geneExpressions[geneLocation.hugoSymbol].AddRegionalExpression(offset, z, mu); // Recall that for allele-specifc expresion, z is really the level of allele-specific expression, not the expression z score.
                        }
                    }
                } // for each line in the input file

                if (!seenDone)
                {
                    Console.WriteLine("Truncated input file " + inputFilename);
                    continue;
                }

                //
                // Write the output file.
                //
                string directory = ExpressionTools.GetDirectoryPathFromFullyQualifiedFilename(inputFilename);
                string analysisId = ExpressionTools.GetAnalysisIdFromPathname(inputFilename);
                if ("" == directory || "" == analysisId) {
                    Console.WriteLine("Couldn't parse input pathname, which is supposed to be absolute and include an analysis ID: " + inputFilename);
                    continue;
                }

                var outputFilename = directory + analysisId + (forAlleleSpecificExpression ? ExpressionTools.alleleSpecificGeneExpressionExtension : ExpressionTools.geneExpressionExtension);

                var outputFile = ExpressionTools.CreateStreamWriterWithRetry(outputFilename);

                outputFile.WriteLine("ExpressionNearMutations v3.1 " + participantId + (forAlleleSpecificExpression ? " -a" : "")); // v3.1 uses ucsc gene locations
                outputFile.Write("Gene name\tnon-silent mutation count");
                for (int sizeIndex = 0; sizeIndex < GeneExpression.nRegionSizes; sizeIndex++)
                {
                    outputFile.Write("\t" + GeneExpression.regionSizeByRegionSizeIndex[sizeIndex] + (forAlleleSpecificExpression ? "(ase)" : "(z)"));
                }

                outputFile.Write("\tWhole Autosome " + (forAlleleSpecificExpression ? "(ase)" : "(z)"));

                if (!forAlleleSpecificExpression) {
                    for (int sizeIndex = 0; sizeIndex < GeneExpression.nRegionSizes; sizeIndex++)
                    {
                        outputFile.Write("\t" + GeneExpression.regionSizeByRegionSizeIndex[sizeIndex] + "(mu)");
                    }
                    outputFile.Write("\tWhole Autosome (mu)");
                }

                for (int sizeIndex = 0; sizeIndex < GeneExpression.nRegionSizes; sizeIndex++)
                {
                    outputFile.Write("\t" + GeneExpression.regionSizeByRegionSizeIndex[sizeIndex] + " exclusive " + (forAlleleSpecificExpression ? "(ase)" : "(z)"));
                }

                outputFile.Write("\tWhole Autosome exclusive " + (forAlleleSpecificExpression ? "(ase)" : "(z)"));

                if (!forAlleleSpecificExpression)
                {
                    for (int sizeIndex = 0; sizeIndex < GeneExpression.nRegionSizes; sizeIndex++)
                    {
                        outputFile.Write("\t" + GeneExpression.regionSizeByRegionSizeIndex[sizeIndex] + " exclusive (mu)");
                    }
                    outputFile.Write("\tWhole Autosome exclusive (mu)");
                }

                for (int whichChromosome = 0; whichChromosome < ExpressionTools.nHumanNuclearChromosomes; whichChromosome++)
                {
                    outputFile.Write("\t" + ExpressionTools.ChromosomeIndexToName(whichChromosome, true));
                }

                if (!forAlleleSpecificExpression)
                {
                    for (int whichChromosome = 0; whichChromosome < ExpressionTools.nHumanNuclearChromosomes; whichChromosome++)
                    {
                        outputFile.Write("\t" + ExpressionTools.ChromosomeIndexToName(whichChromosome, true) + " mu");
                    }
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
                    outputFile.Write(ExpressionTools.ConvertToExcelString(allExpressions[i].geneLocationInfo.hugoSymbol) + "\t" + allExpressions[i].mutationCount);

                    for (int sizeIndex = 0; sizeIndex < GeneExpression.nRegionSizes; sizeIndex++)
                    {
                        if (allExpressions[i].regionalExpressionState[sizeIndex].nRegionsIncluded >= minExamplesPerRegion)
                        {
                            outputFile.Write("\t" + allExpressions[i].regionalExpressionState[sizeIndex].totalExpression / allExpressions[i].regionalExpressionState[sizeIndex].nRegionsIncluded);
                        }
                        else
                        {
                            outputFile.Write("\t*");
                        }
                    }

                    if (wholeAutosomeRegionalExpression.nRegionsIncluded >= minExamplesPerRegion)
                    {
                        outputFile.Write("\t" + wholeAutosomeRegionalExpression.totalExpression / wholeAutosomeRegionalExpression.nRegionsIncluded);
                    }
                    else
                    {
                        outputFile.Write("\t*");
                    }

                    if (!forAlleleSpecificExpression) {
                        for (int sizeIndex = 0; sizeIndex < GeneExpression.nRegionSizes; sizeIndex++)
                        {
                            if (allExpressions[i].regionalExpressionState[sizeIndex].nRegionsIncluded >= minExamplesPerRegion)
                            {
                                outputFile.Write("\t" + allExpressions[i].regionalExpressionState[sizeIndex].totalMeanExpression / allExpressions[i].regionalExpressionState[sizeIndex].nRegionsIncluded);
                            }
                            else
                            {
                                outputFile.Write("\t*");
                            }
                        }

                        if (wholeAutosomeRegionalExpression.nRegionsIncluded >= minExamplesPerRegion)
                        {
                            outputFile.Write("\t" + wholeAutosomeRegionalExpression.totalMeanExpression / wholeAutosomeRegionalExpression.nRegionsIncluded);
                        }
                        else
                        {
                            outputFile.Write("\t*");
                        }
                    }

                    for (int sizeIndex = 0; sizeIndex < GeneExpression.nRegionSizes; sizeIndex++)
                    {
                        if (allExpressions[i].exclusiveRegionalExpressionState[sizeIndex].nRegionsIncluded >= minExamplesPerRegion)
                        {
                            outputFile.Write("\t" + allExpressions[i].exclusiveRegionalExpressionState[sizeIndex].totalExpression / allExpressions[i].exclusiveRegionalExpressionState[sizeIndex].nRegionsIncluded);
                        }
                        else
                        {
                            outputFile.Write("\t*");
                        }
                    }

                    if (allButThisChromosomeAutosomalRegionalExpressionState[allExpressions[i].geneLocationInfo.chromosome].nRegionsIncluded >= minExamplesPerRegion)
                    {
                        outputFile.Write("\t" + allButThisChromosomeAutosomalRegionalExpressionState[allExpressions[i].geneLocationInfo.chromosome].totalExpression / allButThisChromosomeAutosomalRegionalExpressionState[allExpressions[i].geneLocationInfo.chromosome].nRegionsIncluded);
                    }
                    else
                    {
                        outputFile.Write("\t*");
                    }

                    if (!forAlleleSpecificExpression)
                    {
                        for (int sizeIndex = 0; sizeIndex < GeneExpression.nRegionSizes; sizeIndex++)
                        {
                            if (allExpressions[i].exclusiveRegionalExpressionState[sizeIndex].nRegionsIncluded >= minExamplesPerRegion)
                            {
                                outputFile.Write("\t" + allExpressions[i].exclusiveRegionalExpressionState[sizeIndex].totalMeanExpression / allExpressions[i].exclusiveRegionalExpressionState[sizeIndex].nRegionsIncluded);
                            }
                            else
                            {
                                outputFile.Write("\t*");
                            }
                        }

                        if (allButThisChromosomeAutosomalRegionalExpressionState[allExpressions[i].geneLocationInfo.chromosome].nRegionsIncluded >= minExamplesPerRegion)
                        {
                            outputFile.Write("\t" + allButThisChromosomeAutosomalRegionalExpressionState[allExpressions[i].geneLocationInfo.chromosome].totalMeanExpression / allButThisChromosomeAutosomalRegionalExpressionState[allExpressions[i].geneLocationInfo.chromosome].nRegionsIncluded);
                        }
                        else
                        {
                            outputFile.Write("\t*");
                        }
                    }

                    for (int whichChromosome = 0; whichChromosome < ExpressionTools.nHumanNuclearChromosomes; whichChromosome++)
                    {
                        if (perChromosomeRegionalExpressionState[whichChromosome].nRegionsIncluded >= minExamplesPerRegion)
                        {
                            outputFile.Write("\t" + perChromosomeRegionalExpressionState[whichChromosome].totalExpression / perChromosomeRegionalExpressionState[whichChromosome].nRegionsIncluded);
                        }
                        else
                        {
                            outputFile.Write("\t*");
                        }
                    }

                    if (!forAlleleSpecificExpression)
                    {
                        for (int whichChromosome = 0; whichChromosome < ExpressionTools.nHumanNuclearChromosomes; whichChromosome++)
                        {
                            if (perChromosomeRegionalExpressionState[whichChromosome].nRegionsIncluded >= minExamplesPerRegion)
                            {
                                outputFile.Write("\t" + perChromosomeRegionalExpressionState[whichChromosome].totalMeanExpression / perChromosomeRegionalExpressionState[whichChromosome].nRegionsIncluded);
                            }
                            else
                            {
                                outputFile.Write("\t*");
                            }
                        }                        
                    }

                    outputFile.WriteLine();
                } // for each gene

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

        static void PrintUsageMessage()
        {
            Console.WriteLine("usage: ExpressionNearMutations {-a} {-f fileContainingParticipantIds | <participantIdsToProcess>}");
            Console.WriteLine("-a means to use allele-specific expression rather than total expression.");
        }

        static void Main(string[] args)
        {
            if (args.Count() == 0 || args.Count() == 1 && args[0] == "-a")
            {
                PrintUsageMessage();
                return;
            }

            int minExamplesPerRegion;   // If there are fewer than this, then ignore the region.

            bool forAlleleSpecificExpression = args[0] == "-a";

            int argsConsumed = 0;
            if (forAlleleSpecificExpression) {
                argsConsumed = 1;
            }

            if (args[argsConsumed] == "-f" && args.Count() != argsConsumed + 2)
            {
                PrintUsageMessage();
                return;
            }

            if (forAlleleSpecificExpression)
            {
                minExamplesPerRegion = 1;
            }
            else
            {
                minExamplesPerRegion = 1;
            }

            Stopwatch timer = new Stopwatch();
            timer.Start();

            List<string> excludedAnalyses;
            Dictionary<string, ExpressionTools.TCGARecord> tcgaRecords;
            Dictionary<string, string> sampleToParticipantIDMap;
            Dictionary<string, ExpressionTools.Sample> allSamples;

            ExpressionTools.LoadStateFromExperimentsFile(out excludedAnalyses, out tcgaRecords, out sampleToParticipantIDMap, out participants, out experiments, out allSamples);

            timer.Stop();
            Console.WriteLine("Loaded " + experiments.Count() + " experiments with MAFs in " + (timer.ElapsedMilliseconds + 500) / 1000 + "s");

            //
            // Now build the map of mutations by gene.
            //

            timer.Reset();
            timer.Start();

            //mutations = ExpressionTools.GenerateMutationMapFromExperiments(experiments, experimentsByParticipant);

            foreach (var experiment in experiments) // This got done in GenerateMutationMapFromExperiments, but now we need to do it here.
            {
                experimentsByParticipant.Add(experiment.participant.participantId, experiment);
            }

            geneLocationInformation = new Dictionary<string, ExpressionTools.GeneLocationsByNameAndChromosome>(); // // Maps reference type (hg18 or hg19) to map from hugo symbol to location info
            geneLocationInformation.Add("hg18", new ExpressionTools.GeneLocationsByNameAndChromosome(ExpressionTools.LoadGeneLocationInfo(@"\\gcr\scratch\b99\bolosky\knownGene-hg18.txt", @"\\gcr\scratch\b99\bolosky\kgXref-hg18.txt")));
            geneLocationInformation.Add("hg19", new ExpressionTools.GeneLocationsByNameAndChromosome(ExpressionTools.LoadGeneLocationInfo(@"\\gcr\scratch\b99\bolosky\knownGene-hg19.txt", @"\\gcr\scratch\b99\bolosky\kgXref-hg19.txt")));

            timer.Stop();
            Console.WriteLine("Loaded mutations in " + geneLocationInformation["hg19"].genesByName.Count() + " genes in " + (timer.ElapsedMilliseconds + 500) / 1000 + "s.");

            List<string> participantsToProcess;

            if (args[argsConsumed] == "-f")
            {
                participantsToProcess = ExpressionTools.ReadAllLinesWithRetry(args[argsConsumed + 1]).ToList();
            }
            else
            {
                participantsToProcess = new List<string>();
                for (int i = argsConsumed; i < args.Count(); i++)
                {
                    participantsToProcess.Add(args[i]);
                }
            }

            //
            // Process the runs in parallel
            //
            int totalNumberOfExperiments = experiments.Count();
            timer.Reset();
            timer.Start();

            var threads = new List<Thread>();
            for (int i = 0; i < Environment.ProcessorCount/* 1*/; i++)
            {
                threads.Add(new Thread(() => ProcessParticipants(participantsToProcess, forAlleleSpecificExpression, minExamplesPerRegion)));
            }

            threads.ForEach(t => t.Start());
            threads.ForEach(t => t.Join());

            timer.Stop();
            Console.WriteLine("Processed " + (args.Count() - (forAlleleSpecificExpression ? 1 : 0)) + " experiments in " + (timer.ElapsedMilliseconds + 500) / 1000 + " seconds");
        }
    }
}
