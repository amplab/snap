using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using System.IO;
using ASELib;
using System.Diagnostics;
using System.Threading;

namespace ExpressionByMutationCount
{
    class Program
    {
        class ExpressionInstance : IComparer<ExpressionInstance>
        {
            public ExpressionInstance(string tumorType_, int nMutations_, double z_, string case_id_)
            {
                tumorType = tumorType_;
                nMutations = nMutations_;
                z = z_;
                case_id = case_id_;
            }

            public readonly int nMutations;
            public readonly double z;       // Or mu, as appropriate
            public readonly string tumorType;
            public readonly string case_id;

            public int Compare(ExpressionInstance a, ExpressionInstance b)
            {
                if (a.z > b.z) return 1;
                if (a.z < b.z) return -1;
                return 0;
            }
        }

        class Mutation
        {
            public Mutation(bool lowDNA_, bool highDNA_, bool truncatingMutation_, bool insufficientDNA_)
            {
                lowDNA = lowDNA_;
                highDNA = highDNA_;
                truncatingMutation = truncatingMutation_;
                insufficientDNA = insufficientDNA_;
            }

            public void Update(bool lowDNA_, bool highDNA_, bool truncatingMutation_, bool insufficientDNA_)
            {
                lowDNA |= lowDNA_;
                highDNA |= highDNA_;
                truncatingMutation |= truncatingMutation_;
                insufficientDNA |= insufficientDNA_;
                multipleMutations = true;
            }

            public bool lowDNA ;                 // DNA ratio < 1/3; a minor subclone
            public bool highDNA;                // DNA ratio > 2/3; a possible loss of heterozygosity
            public bool truncatingMutation;     // A nonsense or splice site mutation that would reduce the size of the mRNA
            public bool insufficientDNA;

            public bool multipleMutations = false;
        }

        class GeneState
        {
            public GeneState(string hugo_symbol_, List<ASETools.GeneScatterGraphLine> scatterGraphLines_)
            {
                hugo_symbol = hugo_symbol_;
                scatterGraphLines = scatterGraphLines_;

                foreach (var wholeAutosome in ASETools.BothBools)
                {
                    expressionByRange.Add(wholeAutosome, new Dictionary<bool, Dictionary<bool, Dictionary<int, List<ExpressionInstance>>>>());

                    foreach (var usingMu in ASETools.BothBools)
                    {
                        expressionByRange[wholeAutosome].Add(usingMu, new Dictionary<bool, Dictionary<int, List<ExpressionInstance>>>());

                        foreach (var exclusive in ASETools.BothBools)
                        {
                            expressionByRange[wholeAutosome][usingMu].Add(exclusive, new Dictionary<int, List<ExpressionInstance>>());
                        }
                    }
                }

                foreach (var usingMu in ASETools.BothBools)
                {
                    perChromosomeExpression.Add(usingMu, new List<ExpressionInstance>[ASETools.nHumanNuclearChromosomes]);
                }
            }

            public void AddExpression(string tumorType, int range, int nMutations, double z /* or mu, as approproate*/, string case_id, bool usingMu, bool exclusive)
            {
                bool wholeAutosome = range == ASETools.nRegions - 1;

                if (!expressionByRange[wholeAutosome][usingMu][exclusive].ContainsKey(range))
                {
                    expressionByRange[wholeAutosome][usingMu][exclusive].Add(range, new List<ExpressionInstance>());
                }

                expressionByRange[wholeAutosome][usingMu][exclusive][range].Add(new ExpressionInstance(tumorType, nMutations, z, case_id));
            }

            public void AddWholeAutosomeExpression(string tumorType, int nMutations, double z /* or mu, as approproate*/, string case_id, bool usingMu, bool exclusive)
            {
                AddExpression(tumorType, ASETools.nRegions - 1, nMutations, z, case_id, usingMu, exclusive);
            }

            public void AddPerChromosomeExpression(string tumorType, int nMutations, double z /*or mu as appropriate*/, string case_id, bool usingMu, int whichChromosome)
            {
                perChromosomeExpression[usingMu][whichChromosome].Add(new ExpressionInstance(tumorType, nMutations, z, case_id));
            }

            public bool loadPerCaseState(Dictionary<string, string> sampleToParticipantIDMap, Dictionary<string, ASETools.Case> experimentsByRNAAnalysisID, ASETools.Configuration configuration)
            {

                var geneScatterPlotLines = ASETools.GeneScatterGraphLine.LoadAllGeneScatterGraphEntries(configuration.geneScatterGraphsDirectory, true, hugo_symbol);
                if (null == geneScatterPlotLines || geneScatterPlotLines.Count() == 0)
                {
                    return false;
                }

                foreach (var geneScatterPlotLine in geneScatterPlotLines)
                {

                    if (perCaseState.ContainsKey(geneScatterPlotLine.case_id))
                    {
                        perCaseState[geneScatterPlotLine.case_id].Update(
                            geneScatterPlotLine.tumorDNAFraction < 1.0 / 3.0, geneScatterPlotLine.tumorDNAFraction > 2.0 / 3.0, geneScatterPlotLine.Variant_Classification == "Nonsense_Mutation" || geneScatterPlotLine.Variant_Classification == "Splice_Site", 
                            geneScatterPlotLine.tumorDNAReadCounts.nMatchingAlt + geneScatterPlotLine.tumorDNAReadCounts.nMatchingReference < 10);
                    }
                    else
                    {
                        perCaseState.Add(geneScatterPlotLine.case_id,
                            new Mutation(geneScatterPlotLine.tumorDNAFraction < 1.0 / 3.0, geneScatterPlotLine.tumorDNAFraction > 2.0 / 3.0,
                                geneScatterPlotLine.Variant_Classification == "Nonsense_Mutation" || geneScatterPlotLine.Variant_Classification == "Splice_Site", geneScatterPlotLine.tumorDNAReadCounts.nMatchingAlt + geneScatterPlotLine.tumorDNAReadCounts.nMatchingReference < 10));
                    }
                }

                return true;
            }

            public string hugo_symbol;
            public Dictionary<bool, List<ExpressionInstance>[]> perChromosomeExpression = new Dictionary<bool, List<ExpressionInstance>[]>();   // First index is usingMu, then chromsome index
            public Dictionary<string, Mutation> perCaseState = new Dictionary<string, Mutation>();

            //
            // A single data structure to look up all of the expression by range dictionaries.  Order of bools is wholeAutosome, usingMu, exclusive.
            //
            public Dictionary<bool, Dictionary<bool, Dictionary<bool, Dictionary<int, List<ExpressionInstance>>>>> expressionByRange = new Dictionary<bool,Dictionary<bool, Dictionary<bool, Dictionary<int, List<ExpressionInstance>>>>>();

            public readonly List<ASETools.GeneScatterGraphLine> scatterGraphLines;
            public List<string> perGeneLines = new List<string>();

            public int nInputFilesPastThisGene = 0;
        }

        static List<ExpressionInstance> FilterUnusualMutantsIfRequested(List<ExpressionInstance> instances, GeneState geneState, bool excludeUnusualMutants)
        {
            if (!excludeUnusualMutants)
            {
                return instances;
            }

            var retVal = new List<ExpressionInstance>();

            foreach (var instance in instances)
            {
                if (geneState.perCaseState.ContainsKey(instance.case_id))
                {
                    var perExperimentState = geneState.perCaseState[instance.case_id];
                    if (perExperimentState.truncatingMutation || perExperimentState.highDNA || perExperimentState.lowDNA)
                    {
                        //
                        // This is to be excluded, skip it.
                        //
                        continue;
                    }
                }
                retVal.Add(instance);
            }

            return retVal;
            
        }

        static void ComputeMannWhitneyAndPrint(StreamWriter outputFile, List<ExpressionInstance> instances, ASETools.MannWhitney<ExpressionInstance>.WhichGroup selector, bool twoTailed, GeneState geneState, bool excludeUnusualMutants)
        {
            bool enoughData;
            double p = -1;

            if (instances.Count() != 0)
            {
                bool reversed;
                double nFirstGroup;
                double nSecondGroup;
                double U;
                double z;

                var instancesToUse = FilterUnusualMutantsIfRequested(instances, geneState, excludeUnusualMutants);

                p = ASETools.MannWhitney<ExpressionInstance>.ComputeMannWhitney(
                        instancesToUse,
                        instances[0], selector, getValue, out enoughData, out reversed, out nFirstGroup, out nSecondGroup, out U, out z, twoTailed, 10);
            }
            else
            {
                enoughData = false;
            }

            if (enoughData)
            {
                outputFile.Write("\t" + p);
            }
            else
            {
                outputFile.Write("\t*");
            }
        }

        static void ComputeNMeanAndStandardDeviationAndPrint(StreamWriter outputFile, List<ExpressionInstance> instances, GeneState geneState, bool excludeUnusualMutants)
        {
            var values = FilterUnusualMutantsIfRequested(instances, geneState, excludeUnusualMutants).ConvertAll(ExtractZFromExpressionInstance);

            outputFile.Write("\t" + values.Count());
            if (values.Count() == 0)
            {
                outputFile.Write("\t*\t*");
                return;
            }

            outputFile.Write("\t" + ASETools.MeanOfList(values) + "\t" + ASETools.StandardDeviationOfList(values));
        }

        static ASETools.MannWhitney<ExpressionInstance>.GetValue getValue = new ASETools.MannWhitney<ExpressionInstance>.GetValue(x => x.z);

        static void WriteHeaderGroup(string regionName, StreamWriter outputFile, string muString) // Writes out the header for the complete set of measurements for one region
        {
            outputFile.Write(/*"\t" + groupId + " 0 vs. 1" + muString + "\t" + groupId + " 0 vs. not zero" + muString +*/ "\t" + regionName + " 1 vs. many" + muString + "\t" + regionName + " 1 vs. not 1" + muString);
            string[] mutationSets = { "0", "1", ">1" };

            for (int i = 0; i < mutationSets.Count(); i++)
            {
                outputFile.Write("\t" + regionName + " " + mutationSets[i] + " mutation" + muString + " N\t" + regionName + " " + mutationSets[i] + " mutation " + muString + " mean\t" + regionName + " " + mutationSets[i] + " mutation " + muString + " stdDev");
            }
        }

        static void WriteFileHeader(StreamWriter outputFile, bool forAlleleSpecificExpression, bool perChromosome)
        {
            outputFile.Write("Hugo Symbol");
            for (int exclusive = 0; exclusive < 2; exclusive++)
            {
                for (int mu = 0; mu < (forAlleleSpecificExpression ? 1 : 2); mu++)
                {
                    int width = 0;
                    string muString = ((mu == 0) ? "" : " mu") + (exclusive == 0 ? "" : " exclusive");
                    for (int i = 0; i < ASETools.nRegions - 1 /* -1 is because we don't do whole autosome here*/; i++)
                    {
                        WriteHeaderGroup(width + "Kbp", outputFile, muString);

                        if (0 == width)
                        {
                            width = 1;
                        }
                        else
                        {
                            width *= 2;
                        }
                    }

                    WriteHeaderGroup("whole autosome", outputFile, muString);
                }
            }

            if (perChromosome)
            {
                for (int mu = 0; mu < (forAlleleSpecificExpression ? 1 : 2); mu++)
                {
                    string muString = (mu == 0) ? "" : " mu";

                    for (int whichChromosome = 0; whichChromosome < ASETools.nHumanNuclearChromosomes; whichChromosome++)
                    {
                        var name = ASETools.ChromosomeIndexToName(whichChromosome, true);

                        WriteHeaderGroup(ASETools.ChromosomeIndexToName(whichChromosome, true), outputFile, muString);
                    }
                }
            }

            outputFile.WriteLine("\tnTumorsExcluded\tnZero\tnOne\tnMore");
        }

        static void WriteCounts(StreamWriter outputFile, GeneState geneToProcess, string disease, Dictionary<string, ASETools.Case> cases)
        {
            int nExcluded = 0;
            int nZero = 0;
            int nOne = 0;
            int nMore = 0;

            if (!geneToProcess.expressionByRange[true][false][false].ContainsKey(0))
            {
                outputFile.WriteLine("\t*\t*\t*\t*");
                return;
            }

            foreach (var expression in geneToProcess.expressionByRange[true][false][false][0])
            {
                if (disease != "" && cases[expression.case_id].disease() != disease)
                {
                    continue;
                }

                if (geneToProcess.perCaseState.ContainsKey(expression.case_id))
                {
                    var state = geneToProcess.perCaseState[expression.case_id];
                    if (state.highDNA || state.lowDNA || state.truncatingMutation)
                    {
                        //
                        // Ignore excluded ones.
                        //
                        nExcluded++;
                        continue;
                    }
                }

                if (expression.nMutations == 0)
                {
                    nZero++;
                }
                else if (expression.nMutations == 1)
                {
                    nOne++;
                } else {
                    nMore++;
                }
            }

            outputFile.WriteLine("\t" + nExcluded + "\t" + nZero + "\t" + nOne + "\t" + nMore);
        }

        static ASETools.MannWhitney<ExpressionInstance>.WhichGroup isZero = new ASETools.MannWhitney<ExpressionInstance>.WhichGroup(x => x.nMutations == 0);
        static ASETools.MannWhitney<ExpressionInstance>.WhichGroup isOne = new ASETools.MannWhitney<ExpressionInstance>.WhichGroup(x => x.nMutations == 1);
        static ASETools.MannWhitney<ExpressionInstance>.WhichGroup isNotOne = new ASETools.MannWhitney<ExpressionInstance>.WhichGroup(x => x.nMutations != 1);
        static ASETools.MannWhitney<ExpressionInstance>.WhichGroup isMoreThanOne = new ASETools.MannWhitney<ExpressionInstance>.WhichGroup(x => x.nMutations > 1);

        static double ExtractZFromExpressionInstance(ExpressionInstance expressionInstance)
        {
            return expressionInstance.z;
        }

        static void WriteMannWhitneyToFiles(StreamWriter panCancerOutputFile, Dictionary<string, StreamWriter> outputFilesByDisease, List<ExpressionInstance> expressionInstances, GeneState geneToProcess, bool hasEnough)
        {
            if (hasEnough)
            {

                foreach (var perDiseaseOutputFileEntry in outputFilesByDisease)
                {
                    var perDiseaseOutputFile = perDiseaseOutputFileEntry.Value;
                    var disease = perDiseaseOutputFileEntry.Key;

                    //ComputeMannWhitneyAndPrint(perDiseaseOutputFile, expressionInstances.Where(x => x.nMutations < 2 && x.tumorType == disease).ToList(), isZero, true, geneToProcess, true);
                    //ComputeMannWhitneyAndPrint(perDiseaseOutputFile, expressionInstances.Where(x => x.tumorType == disease).ToList(), isZero, true, geneToProcess, true);
                    ComputeMannWhitneyAndPrint(perDiseaseOutputFile, expressionInstances.Where(x => x.nMutations != 0 && x.tumorType == disease).ToList(), isOne, true, geneToProcess, true);
                    ComputeMannWhitneyAndPrint(perDiseaseOutputFile, expressionInstances.Where(x => x.tumorType == disease).ToList(), isOne, false, geneToProcess, true);

                    ComputeNMeanAndStandardDeviationAndPrint(perDiseaseOutputFile, expressionInstances.Where(x => x.nMutations == 0 && x.tumorType == disease).ToList(), geneToProcess, true);
                    ComputeNMeanAndStandardDeviationAndPrint(perDiseaseOutputFile, expressionInstances.Where(x => x.nMutations == 1 && x.tumorType == disease).ToList(), geneToProcess, true);
                    ComputeNMeanAndStandardDeviationAndPrint(perDiseaseOutputFile, expressionInstances.Where(x => x.nMutations >1 && x.tumorType == disease).ToList(), geneToProcess, true);
                }

                //ComputeMannWhitneyAndPrint(panCancerOutputFile, expressionInstances.Where(x => x.nMutations < 2).ToList(), isZero, true, geneToProcess, true);
                //ComputeMannWhitneyAndPrint(panCancerOutputFile, expressionInstances, isZero, true, geneToProcess, true);
                ComputeMannWhitneyAndPrint(panCancerOutputFile, expressionInstances.Where(x => x.nMutations != 0).ToList(), isOne, false, geneToProcess, true);
                ComputeMannWhitneyAndPrint(panCancerOutputFile, expressionInstances, isOne, false, geneToProcess, true);

                ComputeNMeanAndStandardDeviationAndPrint(panCancerOutputFile, expressionInstances.Where(x => x.nMutations == 0).ToList(), geneToProcess, true);
                ComputeNMeanAndStandardDeviationAndPrint(panCancerOutputFile, expressionInstances.Where(x => x.nMutations == 1).ToList(), geneToProcess, true);
                ComputeNMeanAndStandardDeviationAndPrint(panCancerOutputFile, expressionInstances.Where(x => x.nMutations > 1).ToList(), geneToProcess, true);
            }
            else
            {
                foreach (var perDiseaseOutputFileEntry in outputFilesByDisease)
                {
                    perDiseaseOutputFileEntry.Value.Write("\t*\t*\t*\t*\t*\t*\t*\t*\t*\t*\t*");
                }
                panCancerOutputFile.Write("\t*\t*\t*\t*\t*\t*\t*\t*\t*\t*\t*");
            }
        }

        static void PrintUsage()
        {
            Console.WriteLine("usage: ExpressionByMutationCount {-a} {-c}");
            Console.WriteLine("-a means allele-specific");
            Console.WriteLine("-c means do per-chromosome");
        }


        static void Main(string[] args)
        {
            var timer = new Stopwatch();
            timer.Start();

            var configuration = ASETools.Configuration.loadFromFile(args);
            if (null == configuration)
            {
                Console.WriteLine("Unable to load configuration.");
            }

            bool forAlleleSpecificExpression = false;
            bool perChromosome = false;

            foreach (var arg in configuration.commandLineArgs)
            {
                if (arg == "-a")
                {
                    forAlleleSpecificExpression = true;
                }
                else if (arg == "-c")
                {
                    perChromosome = true;
                }
                else
                {
                    PrintUsage();
                    return;
                }
            }

            var comparer = StringComparer.OrdinalIgnoreCase;

            var cases = ASETools.Case.LoadCases(configuration.casesFilePathname);
            if (null == cases)
            {
                Console.WriteLine("You must generate cases first.");
                return;
            }

            var selectedGenes = ASETools.SelectedGene.LoadFromFile(configuration.selectedGenesFilename);
            if (null == selectedGenes)
            {
                Console.WriteLine("Must first select genes.");
                return;
            }

            int missingCount = cases.Where(caseEntry => (forAlleleSpecificExpression ? caseEntry.Value.tumor_allele_specific_gene_expression_filename : caseEntry.Value.gene_expression_filename) == "").Count();



            var casesToProcess = new List<ASETools.Case>();
            foreach (var caseEntry in cases)
            {
                var case_ = caseEntry.Value;

                if (forAlleleSpecificExpression && case_.tumor_allele_specific_gene_expression_filename != "" ||
                    !forAlleleSpecificExpression && case_.gene_expression_filename != "")
                {
                    casesToProcess.Add(case_);
                }
            }

            var genesToProcess = new Dictionary<string, GeneState>();

            int nGenesSkipped = 0;
            foreach (var selectedGene in selectedGenes)
            {
                //
                // Load the unfiltered scatter graph for this gene into memory.
                //

                var geneScatterPlotLines = ASETools.GeneScatterGraphLine.LoadAllGeneScatterGraphEntries(configuration.geneScatterGraphsDirectory, true, selectedGene.Hugo_Symbol);
                if (geneScatterPlotLines.Count() == 0)
                {
                    //
                    // Probably not enough tumors to make the cut, skip this gene.
                    //
                    nGenesSkipped++;
                    continue;
                }

               genesToProcess.Add(selectedGene.Hugo_Symbol, new GeneState(selectedGene.Hugo_Symbol, geneScatterPlotLines));
            }

            Console.WriteLine("Loaded " + cases.Count() + " cases, of which " + missingCount + " are missing gene expression, and " + genesToProcess.Count() + " genes in " + ASETools.ElapsedTimeInSeconds(timer));

            timer.Reset();
            timer.Start();

            //
            // Run through all of the {AlleleSpecific} expression files and add them to the per-gene state.  Do this in parallel, since there are a lot of them.
            //
            var threads = new List<Thread>();
            for (int i = 0; i < Environment.ProcessorCount * 2 /* *2 because this is often IO bound*/; i++)
            {
                threads.Add(new Thread(() => ProcessRegionalExpressionFile(forAlleleSpecificExpression, casesToProcess, genesToProcess)));
            }

            Console.Write("Loading expression files (1 dot/100): ");
            threads.ForEach(t => t.Start());
            threads.ForEach(t => t.Join());

            Console.WriteLine();
            Console.WriteLine("Loaded and processed expression files in " + ASETools.ElapsedTimeInSeconds(timer));

Console.WriteLine("lock: " + lockMilliseconds + "ms; read: " + readMilliseconds + "ms, process: " + processMilliseconds + "ms, 1: " + process1Milliseconds + "ms, 2: " + process2Milliseconds + "ms, 3: " + process3Milliseconds + "ms.");

            string baseFileName = configuration.finalResultsDirectory + (forAlleleSpecificExpression ? "AlleleSpecific" : "") + "ExpressionDistributionByMutationCount";

            var panCancerOutputFile = ASETools.CreateStreamWriterWithRetry(baseFileName + ".txt");
            WriteFileHeader(panCancerOutputFile, forAlleleSpecificExpression, perChromosome);
 
            var outputFilesByDisease = new Dictionary<string, StreamWriter>();

            var listOfDiseases = ASETools.GetListOfDiseases(cases);

            foreach (var disease in listOfDiseases)
            {
                outputFilesByDisease.Add(disease, ASETools.CreateStreamWriterWithRetry(baseFileName + "_" + disease + ".txt"));
                WriteFileHeader(outputFilesByDisease[disease], forAlleleSpecificExpression, perChromosome);
            }
            
            timer.Reset();
            timer.Start();

            //
            // Now generate output for all of the genes.
            //
            foreach (var geneToProcessEntry in genesToProcess)
            {
                var geneToProcess = geneToProcessEntry.Value;

                var perGeneLinesFile = StreamWriter.Null; /* ASETools.CreateStreamWriterWithRetry(configuration.regionalExpressionDirectory + geneToProcess.hugo_symbol.ToLower() + (forAlleleSpecificExpression ? "_allele_specific" : "") + "_lines.txt");*/
                perGeneLinesFile.Write("Case ID\tdisease abbr.\tHugo Symbol\tMutation Count");

                foreach (var exclusive in ASETools.BothBools)
                {
                    ulong width = 1000;
                    perGeneLinesFile.Write("\t0" + (exclusive ? " exclusive" : ""));

                    for (int i = 1; i < ASETools.nRegions - 1 /* -1 because we don't do whole autosome here*/; i++)
                    {
                        perGeneLinesFile.Write("\t" + ASETools.SizeToUnits(width) + "b" + (exclusive ? " exclusive" : ""));
                        width *= 2;
                    }

                    perGeneLinesFile.Write("\tWhole Autosome" + (exclusive ? " exclusive" : ""));

                    if (!forAlleleSpecificExpression)
                    {

                        perGeneLinesFile.Write("\t0" + (exclusive ? " exclusive" : "") + " mean");

                        width = 1000;
                        for (int i = 1; i < ASETools.nRegions - 1 /* -1 because we don't do whole autosome here*/; i++)
                        {
                            perGeneLinesFile.Write("\t" + ASETools.SizeToUnits(width) + "b mean" + (exclusive ? " exclusive" : ""));
                            width *= 2;
                        }

                        perGeneLinesFile.Write("\tWhole Autosome mean " + (exclusive ? " exclusive" : ""));
                    }
                } // exclusive

                //
                // Even if we're not doing per-chromosome, write the headers, since they're in the input file that's being copied here.
                //
                for (int whichChromosome = 0; whichChromosome < ASETools.nHumanNuclearChromosomes; whichChromosome++)
                {
                    perGeneLinesFile.Write("\t" + ASETools.ChromosomeIndexToName(whichChromosome, true));
                }

                if (!forAlleleSpecificExpression)
                {
                    for (int whichChromosome = 0; whichChromosome < ASETools.nHumanNuclearChromosomes; whichChromosome++)
                    {
                        perGeneLinesFile.Write("\t" + ASETools.ChromosomeIndexToName(whichChromosome, true) + " mean");
                    }
                }

                perGeneLinesFile.WriteLine();
                foreach (var line in geneToProcess.perGeneLines)
                {
                    perGeneLinesFile.WriteLine(line);
                }
                perGeneLinesFile.Close();

                //
                // Compute and write out the results.
                //

                panCancerOutputFile.Write(ASETools.ConvertToExcelString(geneToProcess.hugo_symbol));
                foreach (var perDiseaseOutputFileEntry in outputFilesByDisease)
                {
                    perDiseaseOutputFileEntry.Value.Write(ASETools.ConvertToExcelString(geneToProcess.hugo_symbol));
                }

                foreach (var exclusive in ASETools.BothBools)
                {
                    bool[] musToUse = { false };
                    if (!forAlleleSpecificExpression)
                    {
                        musToUse = ASETools.BothBools;
                    }
                    foreach (var mu in musToUse)  
                    {
                        Dictionary<int, List<ExpressionInstance>> byRange;
                        List<ExpressionInstance> wholeAutosome;

                        byRange = geneToProcess.expressionByRange[false][mu][exclusive];

                        if (byRange.Count() == 0)
                        {
                            continue;
                        }
                            
                        wholeAutosome = geneToProcess.expressionByRange[true][mu][exclusive][ASETools.nRegions - 1];

                        for (int i = 0; i < ASETools.nRegions - 1 /* -1 because we don't do whole autosome here*/; i++)
                        {
                            if (byRange.ContainsKey(i))
                            {
                                WriteMannWhitneyToFiles(panCancerOutputFile, outputFilesByDisease, byRange[i], geneToProcess, true);
                            }
                            else
                            {
                                WriteMannWhitneyToFiles(panCancerOutputFile, outputFilesByDisease, null, geneToProcess, false);
                            }
                        } // For all widths

                        WriteMannWhitneyToFiles(panCancerOutputFile, outputFilesByDisease, wholeAutosome, geneToProcess, wholeAutosome.Count() > 0);
                    } // mu
                } // exclusive

                if (perChromosome)
                {
                    for (int whichChromosome = 0; whichChromosome < ASETools.nHumanNuclearChromosomes; whichChromosome++)
                    {
                        WriteMannWhitneyToFiles(panCancerOutputFile, outputFilesByDisease, geneToProcess.perChromosomeExpression[false][whichChromosome],
                            geneToProcess, geneToProcess.perChromosomeExpression[false][whichChromosome].Count() > 0);
                    }

                    if (!forAlleleSpecificExpression)
                    {
                        for (int whichChromosome = 0; whichChromosome < ASETools.nHumanNuclearChromosomes; whichChromosome++)
                        {
                            WriteMannWhitneyToFiles(panCancerOutputFile, outputFilesByDisease, geneToProcess.perChromosomeExpression[true][whichChromosome],
                                geneToProcess, geneToProcess.perChromosomeExpression[true][whichChromosome].Count() > 0);
                        }
                    }
                }

                WriteCounts(panCancerOutputFile, geneToProcess, "", cases);
                foreach (var perDiseaseOutputFileEntry in outputFilesByDisease)
                {
                    WriteCounts(perDiseaseOutputFileEntry.Value, geneToProcess, perDiseaseOutputFileEntry.Key, cases);
                }
            }  // foreach gene



            panCancerOutputFile.Close();

            foreach (var perDiseaseOutputFileEntry in outputFilesByDisease)
            {
                perDiseaseOutputFileEntry.Value.Close();
            }
        } // Main

        static int nLoaded = 0;

        static long lockMilliseconds = 0;
        static long readMilliseconds = 0;
        static long processMilliseconds = 0;
        static long process1Milliseconds = 0;
        static long process2Milliseconds = 0;
        static long process3Milliseconds = 0;

        static void ProcessRegionalExpressionFile(bool forAlleleSpecificExpression, List<ASETools.Case>cases, Dictionary<string, GeneState> genesToProcessInput)
        {
            //
            // Randomize our list of genes, so that we don't wind up with a lock convoy with all of the threads processing them in the same order.
            //

            Random rand = new Random();
            var genesToProcess = genesToProcessInput.OrderBy(c => rand.Next()).Select(c => c.Value).ToList();

            var lockTimer = new Stopwatch();
            var readTimer = new Stopwatch();
            var processTimer = new Stopwatch();
            var process1Timer = new Stopwatch();
            var process2Timer = new Stopwatch();
            var process3Timer = new Stopwatch();

            while (true)
            {
                ASETools.Case case_;
                lock(cases)
                {
                    if (cases.Count() == 0) {
                        lockMilliseconds += lockTimer.ElapsedMilliseconds;
                        readMilliseconds += readTimer.ElapsedMilliseconds;
                        processMilliseconds += processTimer.ElapsedMilliseconds;
                        process1Milliseconds += process1Timer.ElapsedMilliseconds;
                        process2Milliseconds += process2Timer.ElapsedMilliseconds;
                        process3Milliseconds += process3Timer.ElapsedMilliseconds;
                        return;
                    }

                    case_ = cases[0];
                    cases.RemoveAt(0);
                }

                readTimer.Start();
                var regionalSignalFile = ASETools.RegionalSignalFile.ReadFile(forAlleleSpecificExpression ? case_.tumor_allele_specific_gene_expression_filename : case_.gene_expression_filename, true, true);
                readTimer.Stop();

                var expressionForThisCase = new Dictionary<string, ASETools.AlleleSpecificSignal>();
                foreach (var regionalSignalEntry in regionalSignalFile.Item1)
                {
                    expressionForThisCase.Add(regionalSignalEntry.Key, new ASETools.AlleleSpecificSignal(regionalSignalEntry.Key, regionalSignalEntry.Value));
                }

                var disease = case_.disease();

                foreach (var geneToProcess in genesToProcess)
                {

                    if (!expressionForThisCase.ContainsKey(geneToProcess.hugo_symbol))
                    {
                        continue;
                    }

                    var thisExpression = expressionForThisCase[geneToProcess.hugo_symbol];

                    lockTimer.Start();
                    lock (geneToProcess)
                    {
                        lockTimer.Stop();
                        processTimer.Start();

                        process1Timer.Start();
                        //geneToProcess.perGeneLines.Add(case_.case_id + "\t" + disease + "\t" + thisExpression.OutputString());
                        process1Timer.Stop();

                        for (int i = 0; i < ASETools.nRegions; i++)
                        {
                            process2Timer.Start();
                            if (thisExpression.nonExclusive[i] != double.NegativeInfinity)
                            {
                                geneToProcess.AddExpression(disease, i, thisExpression.mutationCount, thisExpression.nonExclusive[i], case_.case_id, false, false);
                            }
                            process2Timer.Stop();

                            process3Timer.Start();
                            if (thisExpression.exclusive[i] != double.NegativeInfinity)
                            {
                                geneToProcess.AddExpression(disease, i, thisExpression.mutationCount, thisExpression.exclusive[i], case_.case_id, false, true);
                            }
                            process3Timer.Stop();
                        }
                        processTimer.Stop();
                    } // lock
                }

                lock (cases) {
                    nLoaded++;
                    if (nLoaded % 100 == 0) Console.Write(".");
                }
            } // while true
        } // LoadRegionalSignalWorker
    }
}
