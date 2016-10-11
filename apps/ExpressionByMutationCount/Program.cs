using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using System.IO;
using ExpressionLib;
using System.Diagnostics;
using System.Threading;

namespace ExpressionByMutationCount
{
    class Program
    {
        class GeneExpressionFile {
            public GeneExpressionFile(string filename_, string tumorType_, string participantID_, bool forAlleleSpecificExpression_) {
                filename = filename_;
                tumorType = tumorType_;
                participantID = participantID_;
                forAlleleSpecificExpression = forAlleleSpecificExpression_;

                int headerPrefixLength = headerPrefix.Count();

                reader = new StreamReader(filename);

                var line = reader.ReadLine();

                if (line == null || line.Count() < headerPrefixLength || line.Substring(0, headerPrefixLength) != headerPrefix || line.Substring(headerPrefixLength, 2) != "2.")
                {
                    Console.WriteLine("Truncated corrupt or wrong version gene expression file " + filename);
                    seenDone = true;
                    return;
                }

                line = reader.ReadLine();
                if (null == line)
                {
                    Console.WriteLine("Truncated gene expresion file at header line " + filename);
                    seenDone = true;
                    return;
                }
            }

            public readonly string filename;
            public readonly string tumorType;
            public readonly string participantID;
            public StreamReader reader = null;
            public bool seenDone = false;
            public string currentHugoSymbol = null;
            public int currentNMutations = 0;
            public double[] currentZ = new double[nWidths];
            public bool[] currentZValid = new bool[nWidths];
            public double[] currentMu = new double[nWidths];
            public bool[] currentMuValid = new bool[nWidths];
            public const int nWidths = 20;
            const string headerPrefix = "ExpressionNearMutations v";
            public string currentLine = null;
            bool forAlleleSpecificExpression;

            System.Random random = new System.Random();

            public bool GetNextLine()
            {
                //
                // We queue up a bunch of lines so that we're not having the system cache manager read all 8K files at once,
                // resulting in an insane incast that can take down the net of the machine running this program.  We randomize the
                // amount we read to avoid having all the readers read at the same time, resulting in the same incast we were 
                // trying to avoid.
                //
                if (queuedLines.Count() == 0)
                {
                    int nExtraLinesToQueue = random.Next(minLinesToQueue);
                    for (int i = 0; i < minLinesToQueue + nExtraLinesToQueue; i++)
                    {
                        var lineToQueue = reader.ReadLine();

                        if (lineToQueue == null)
                        {
                            break;
                        }

                        queuedLines.Add(lineToQueue);
                    }
                }

                if (queuedLines.Count() == 0)
                {
                    if (!seenDone)
                    {
                        Console.WriteLine("File " + filename + " is truncated.");
                    }

                    return false;
                }

                string line = queuedLines[0];
                queuedLines.RemoveAt(0);

                currentLine = line;

                if (seenDone)
                {
                    Console.WriteLine("File " + filename + " extends after **done**");
                    return false;
                }

                if ("**done**" == line)
                {
                    seenDone = true;
                    return false;
                }

                var fields = line.Split('\t');
                if (fields.Count() != (forAlleleSpecificExpression ? 1 : 2) * nWidths + 2)
                {
                    Console.WriteLine("File " + filename + " has unparsable line " + line);
                    return false;
                }

                currentHugoSymbol = fields[0];

                try
                {
                    currentNMutations = Convert.ToInt32(fields[1]);
                    for (int i = 0; i < nWidths; i++)
                    {
                        if (fields[i + 2] == "*")
                        {
                            currentZValid[i] = false;
                        }
                        else
                        {
                            currentZ[i] = Convert.ToDouble(fields[i + 2]);
                            currentZValid[i] = true;
                        }

                        if (forAlleleSpecificExpression || fields[i + nWidths + 2] == "*")
                        {
                            currentMuValid[i] = false;
                        }
                        else
                        {
                            currentMu[i] = Convert.ToDouble(fields[i + nWidths + 2]);
                            currentMuValid[i] = true;
                        }
                    }
                }
                catch (FormatException)
                {
                    Console.WriteLine("Corrupt data line in " + filename + ": " + line);
                    return false;
                }

                return true;
            }

            List<string> queuedLines = new List<string>();
            const int minLinesToQueue = 1000;
        }

        class ExpressionInstance : IComparer<ExpressionInstance>
        {
            public ExpressionInstance(string tumorType_, int nMutations_, double z_, string participantID_)
            {
                tumorType = tumorType_;
                nMutations = nMutations_;
                z = z_;
                participantID = participantID_;
            }

            public readonly int nMutations;
            public readonly double z;       // Or mu, as appropriate
            public readonly string tumorType;
            public readonly string participantID;

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
            public GeneState(string hugo_symbol_)
            {
                hugo_symbol = hugo_symbol_;
            }

            public void AddExpression(string tumorType, int range, int nMutations, double z /* or mu, as approproate*/, string participantID, bool usingMu)
            {
                Dictionary<int, List<ExpressionInstance>> byRange;
                if (usingMu) {
                    byRange = meanExpressionByRange;
                } else {
                    byRange = expressionByRange;
                }

                if (!byRange.ContainsKey(range))
                {
                    byRange.Add(range, new List<ExpressionInstance>());
                }
                byRange[range].Add(new ExpressionInstance(tumorType, nMutations, z, participantID));
            }

            public bool loadPerExperimentState(Dictionary<string, string> sampleToParticipantIDMap, Dictionary<string, ExpressionTools.Experiment> experimentsByRNAAnalysisID)
            {
                string scatterPlotFilename = @"f:\temp\gene_mutations_with_counts\" + hugo_symbol + "_unfiltered_counts.txt";
                if (!File.Exists(scatterPlotFilename))
                {
                    return false;
                }

                StreamReader reader = null;
                bool threw = false;

                do
                {
                    threw = false;
                    try
                    {
                        reader = new StreamReader(scatterPlotFilename);
                    }
                    catch (IOException)
                    {
                        Console.WriteLine("IOException opening file " + scatterPlotFilename + " perhaps it's open in another app (like Excel).  Sleeping 10s and retrying.");
                        threw = true;
                        Thread.Sleep(10000);
                    }
                } while (threw);

                reader.ReadLine(); // Skip the header

                string line = null;

                while (null != (line = reader.ReadLine())) {
                    var fields = line.Split('\t');

                    int nDNAMatchingReference = Convert.ToInt32(fields[39]);
                    int nDNAMatchingTumor = Convert.ToInt32(fields[40]);
                    int nDNAMatchingNeither = Convert.ToInt32(fields[41]);
                    bool insufficientDNA = nDNAMatchingReference + nDNAMatchingTumor + nDNAMatchingNeither < 30;

                    double tumorDNAFraction;
                    if (insufficientDNA) {
                        tumorDNAFraction = 0.5;
                    } else {
                        tumorDNAFraction = (double)nDNAMatchingTumor / ((double)(nDNAMatchingReference + nDNAMatchingTumor + nDNAMatchingNeither));
                    }
                    
                    string variantClassification = fields[10];

                    //
                    // Because I cleverly didn't put the participant ID in scatter graph file, we just parse the tumor RNA analysis ID out of the
                    // filename, which is of format <disease_abbr>\<analysis_id>-<gene>-<chromosome>-<begin-offset>-<end-offset>-RNA
                    //
                    string tumorRNAAnalysisID = fields[0].Substring(fields[0].IndexOf('\\') + 1, 36);
                    if (!experimentsByRNAAnalysisID.ContainsKey(tumorRNAAnalysisID))
                    {
                        Console.WriteLine("Can't find experiment for line " + line);
                        continue;
                    }

                    string participantId = experimentsByRNAAnalysisID[tumorRNAAnalysisID].participant.participantId;
                    if (perExperimentState.ContainsKey(participantId)) {
                        perExperimentState[participantId].Update(tumorDNAFraction < 1.0 / 3.0, tumorDNAFraction > 2.0 / 3.0, variantClassification == "Nonsense_Mutation" || variantClassification == "Splice_Site", insufficientDNA);
                    } else {
                        perExperimentState.Add(participantId, new Mutation(tumorDNAFraction < 1.0 / 3.0, tumorDNAFraction > 2.0 / 3.0, variantClassification == "Nonsense_Mutation" || variantClassification == "Splice_Site", insufficientDNA));
                    }
                }

                return true;
            }

            public string hugo_symbol;
            public Dictionary<int, List<ExpressionInstance>> expressionByRange = new Dictionary<int, List<ExpressionInstance>>();
            public Dictionary<int, List<ExpressionInstance>> meanExpressionByRange = new Dictionary<int, List<ExpressionInstance>>();   // i.e., using mu not z
            public Dictionary<string, Mutation> perExperimentState = new Dictionary<string, Mutation>();

            public int nInputFilesPastThisGene = 0;
        }

        static void ComputeMannWhitneyAndPrint(StreamWriter outputFile, List<ExpressionInstance> instances, ExpressionTools.MannWhitney<ExpressionInstance>.WhichGroup selector, bool twoTailed, GeneState geneState, bool excludeUnusualMutants)
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

                List<ExpressionInstance> instancesToUse;

                if (excludeUnusualMutants)
                {
                    instancesToUse = new List<ExpressionInstance>();

                    foreach (var instance in instances)
                    {
                        if (geneState.perExperimentState.ContainsKey(instance.participantID))
                        {
                            var perExperimentState = geneState.perExperimentState[instance.participantID];
                            if (!(perExperimentState.truncatingMutation || perExperimentState.highDNA || perExperimentState.lowDNA))
                            {
                                instancesToUse.Add(instance);
                            }
                        }
                        else
                        {
                            instancesToUse.Add(instance);
                        }
                    }
                }
                else
                {
                    instancesToUse = instances;
                }

                p = ExpressionTools.MannWhitney<ExpressionInstance>.ComputeMannWhitney(
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

        static ExpressionTools.MannWhitney<ExpressionInstance>.GetValue getValue = new ExpressionTools.MannWhitney<ExpressionInstance>.GetValue(x => x.z);

        static void WriteFileHeader(StreamWriter outputFile, bool forAlleleSpecificExpression)
        {
            outputFile.Write("Hugo Symbol");
            int width = 0;
            for (int mu = 0; mu < (forAlleleSpecificExpression ? 1 : 2); mu++)
            {
                string muString = (mu == 0) ? "" : " mu";
                for (int i = 0; i < GeneExpressionFile.nWidths; i++)
                {
                    outputFile.Write("\t" + width + "Kbp 0 vs. 1" + muString + "\t" + width + "Kbp 0 vs. not zero" + muString + "\t" + width + "Kbp 1.vs many" + muString + "\t" + width + "Kbp 1 vs. not 1" + muString);

                    if (0 == width)
                    {
                        width = 1;
                    }
                    else
                    {
                        width *= 2;
                    }
                }
            }

            outputFile.WriteLine("\tnTumorsExcluded\tnZero\tnOne\tnMore");
        }

        static void WriteCounts(StreamWriter outputFile, GeneState geneToProcess, string disease, Dictionary<string, ExpressionTools.Experiment> experimentsByParticipantID)
        {
            int nExcluded = 0;
            int nZero = 0;
            int nOne = 0;
            int nMore = 0;

            if (!geneToProcess.expressionByRange.ContainsKey(0))
            {
                outputFile.WriteLine("\t*\t*\t*\t*");
                return;
            }

            foreach (var expression in geneToProcess.expressionByRange[0])
            {
                if (disease != "" && experimentsByParticipantID[expression.participantID].disease_abbr != disease)
                {
                    continue;
                }

                if (geneToProcess.perExperimentState.ContainsKey(expression.participantID))
                {
                    var state = geneToProcess.perExperimentState[expression.participantID];
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


        static void Main(string[] args)
        {
            if (args.Count() > 1 || args.Count() == 1 && args[0] != "-a")
            {
                Console.WriteLine("usage: ExpressionByMutationCount {-a}");
                return;
            }

            bool forAlleleSpecificExpression = args.Count() == 1;

            Stopwatch timer = new Stopwatch();
            timer.Start();

            var comparer = StringComparer.OrdinalIgnoreCase;

            var excludedAnalyses = ExpressionTools.LoadExcludedAnalyses();

            var tcgaRecords = ExpressionTools.LoadTCGARecords(null /* stored BAMs*/, excludedAnalyses);
            ExpressionTools.LoadTCGARecordsForLocalRealigns(tcgaRecords, null);
            ExpressionTools.LoadTCGAAdditionalMetadata(tcgaRecords);

            var sampleToParticipantIDMap = ExpressionTools.CreateSampleToParticipantIDMap(tcgaRecords);

            Dictionary<string, ExpressionTools.Sample> allSamples;
            var participants = ExpressionTools.BuildParticipantData(tcgaRecords, out allSamples);

            var experiments = ExpressionTools.LoadExperimentsFromFile(@"f:\temp\expression\experiments.txt", participants, tcgaRecords);

            var experimentsByRNAAnalysisID = new Dictionary<string, ExpressionTools.Experiment>();
            var experimentsByParticipantId = new Dictionary<string, ExpressionTools.Experiment>();
            foreach (var experiment in experiments)
            {
                experimentsByRNAAnalysisID.Add(experiment.TumorRNAAnalysis.analysis_id, experiment);
                experimentsByParticipantId.Add(experiment.participant.participantId, experiment);
            }

            timer.Stop();
            Console.WriteLine("Loaded " + experiments.Count() + " experiments in " + (timer.ElapsedMilliseconds + 500) / 1000 + "s, of which " + 
                experiments.Where(e => e.TumorRNAAnalysis == null || e.TumorRNAAnalysis.geneExpressionFileName == null || e.TumorRNAAnalysis.geneExpressionFileName == "").Count() + " are missing gene expression.");

            var inputFiles = new List<GeneExpressionFile>();

            timer.Reset();
            timer.Start();
            foreach (var experiment in experiments)
            {
                string filename = forAlleleSpecificExpression ? experiment.NormalDNAAnalysis.alleleSpecificGeneExpressionFileName : experiment.TumorRNAAnalysis.geneExpressionFileName;
                if (filename != null && filename != "")
                {
                    var file = new GeneExpressionFile(filename, experiment.disease_abbr, experiment.participant.participantId, forAlleleSpecificExpression);
  
                    inputFiles.Add(file);
                }
            }

            //
            // Run through all of the input files and get the first line.
            //
            foreach (var inputFile in inputFiles)
            {
                if (!inputFile.GetNextLine())
                {
                    Console.WriteLine("Input file " + inputFile.filename + " appears to be empty (?)");
                }
            }

            timer.Stop();
            Console.WriteLine("Opened " + inputFiles.Count() + " files and read their first lines in " + (timer.ElapsedMilliseconds + 500) / 1000 + "s");

            string baseFileName = @"f:\temp\expression\" + (forAlleleSpecificExpression ? "AlleleSpecific" : "") + "ExpressionDistributionByMutationCount";

            var panCancerOutputFile = new StreamWriter(baseFileName + ".txt");
            WriteFileHeader(panCancerOutputFile, forAlleleSpecificExpression);
 
            var outputFilesByDisease = new Dictionary<string, StreamWriter>();

            var listOfDiseases = ExpressionTools.GetListOfDiseases(experiments);

///*BJB - off for now*/ listOfDiseases = new List<string>();

            foreach (var disease in listOfDiseases)
            {
                outputFilesByDisease.Add(disease, new StreamWriter(baseFileName + "_" + disease + ".txt"));
                WriteFileHeader(outputFilesByDisease[disease], forAlleleSpecificExpression);
            }

            ExpressionTools.MannWhitney<ExpressionInstance>.WhichGroup isZero = new ExpressionTools.MannWhitney<ExpressionInstance>.WhichGroup(x => x.nMutations == 0);
            ExpressionTools.MannWhitney<ExpressionInstance>.WhichGroup isOne = new ExpressionTools.MannWhitney<ExpressionInstance>.WhichGroup(x => x.nMutations == 1);
            ExpressionTools.MannWhitney<ExpressionInstance>.WhichGroup isNotOne = new ExpressionTools.MannWhitney<ExpressionInstance>.WhichGroup(x => x.nMutations != 1);
            ExpressionTools.MannWhitney<ExpressionInstance>.WhichGroup isMoreThanOne = new ExpressionTools.MannWhitney<ExpressionInstance>.WhichGroup(x => x.nMutations > 1);

            char currentFirstLetter = 'A';
            int nGenesWithCurrentFirstLetter = 0;
            timer.Reset();
            timer.Start();

            int nGenesProcessed = 0;
            int nGenesSkipped = 0;
            var fiveMinuteTimer = new Stopwatch();
            fiveMinuteTimer.Start();
            long lastFiveMinutePrint = 0;

            while (inputFiles.Count() > 0)
            {

                string tooBigHugoSymbol = "ZZZZZZZZZZ";   // Larger than the larges gene name, which is ZZZ3
                string nextHugoSymbol = tooBigHugoSymbol;

                foreach (var inputFile in inputFiles)
                {
                    if (comparer.Compare(inputFile.currentHugoSymbol, nextHugoSymbol) < 0)
                    {
                        nextHugoSymbol = inputFile.currentHugoSymbol;
                    }
                }

                if (nextHugoSymbol == tooBigHugoSymbol)
                {
                    Console.WriteLine("Some input files seem to have a hugo symbol bigger than " + tooBigHugoSymbol);
                    break;
                }

                if (Char.ToUpper(nextHugoSymbol[0]) != currentFirstLetter)
                {
                    timer.Stop();
                    Console.WriteLine("Processed " + nGenesWithCurrentFirstLetter + " genes starting with " + currentFirstLetter + " in " + (timer.ElapsedMilliseconds + 500) / 1000 + "s");
                    timer.Reset();
                    timer.Start();
                    nGenesWithCurrentFirstLetter = 0;
                    currentFirstLetter = nextHugoSymbol[0];
                }
                nGenesWithCurrentFirstLetter++;

                var completedInputFiles = new List<GeneExpressionFile>();
                var geneToProcess = new GeneState(nextHugoSymbol);

                var perGeneLines = new List<string>();

                foreach (var inputFile in inputFiles.Where(f => f.currentHugoSymbol == nextHugoSymbol))
                {
                    perGeneLines.Add(inputFile.participantID + "\t" + inputFile.currentLine);

                    for (int i = 0 ; i < GeneExpressionFile.nWidths; i++) {
                        if (inputFile.currentZValid[i])
                        {
                            geneToProcess.AddExpression(inputFile.tumorType, i, inputFile.currentNMutations, inputFile.currentZ[i], inputFile.participantID, false);
                        }

                        if (!forAlleleSpecificExpression && inputFile.currentMuValid[i]) {
                            geneToProcess.AddExpression(inputFile.tumorType, i, inputFile.currentNMutations, inputFile.currentMu[i], inputFile.participantID, true);
                        }
                    }

                    if (!inputFile.GetNextLine())
                    {
                        completedInputFiles.Add(inputFile);
                    }
                }


                if (!geneToProcess.loadPerExperimentState(sampleToParticipantIDMap, experimentsByRNAAnalysisID))
                {
                    //
                    // There wasn't a scatter graphs file for it, probably because it had too few mutations.  Skip it.
                    //
                    nGenesSkipped++;
                }
                else
                {

                    var perGeneLinesFile = new StreamWriter(@"f:\temp\expression\RegionalExpressionByGene\" + nextHugoSymbol.ToLower() + (forAlleleSpecificExpression ? "_allele_specific" : "") + "_lines.txt");
                    perGeneLinesFile.Write("ParticipantID\tHugo Symbol\tMutation Count\t0");
                    ulong width = 1000;
                    for (int i = 1; i < GeneExpressionFile.nWidths; i++)
                    {
                        perGeneLinesFile.Write("\t" + ExpressionTools.SizeToUnits(width) + "b");
                        width *= 2;
                    }

                    foreach (var line in perGeneLines)
                    {
                        perGeneLinesFile.WriteLine(line);
                    }
                    perGeneLinesFile.Close();


                    //
                    // Compute and write out the results.
                    //

                    panCancerOutputFile.Write(nextHugoSymbol);
                    foreach (var perDiseaseOutputFileEntry in outputFilesByDisease)
                    {
                        perDiseaseOutputFileEntry.Value.Write(nextHugoSymbol);
                    }

                    for (int mu = 0; mu < (forAlleleSpecificExpression ? 1 : 2); mu++)  // shouldn't we be able to do this for a bool somehow?
                    {
                        Dictionary<int, List<ExpressionInstance>> byRange;
                        if (mu == 0)
                        {
                            byRange = geneToProcess.expressionByRange;
                        }
                        else
                        {
                            byRange = geneToProcess.meanExpressionByRange;
                        }

                        for (int i = 0; i < GeneExpressionFile.nWidths; i++)
                        {
                            if (byRange.ContainsKey(i))
                            {
                                foreach (var perDiseaseOutputFileEntry in outputFilesByDisease)
                                {
                                    var perDiseaseOutputFile = perDiseaseOutputFileEntry.Value;
                                    var disease = perDiseaseOutputFileEntry.Key;

                                    ComputeMannWhitneyAndPrint(perDiseaseOutputFile, byRange[i].Where(x => x.nMutations < 2 && x.tumorType == disease).ToList(), isZero, true, geneToProcess, true);
                                    ComputeMannWhitneyAndPrint(perDiseaseOutputFile, byRange[i].Where(x => x.tumorType == disease).ToList(), isZero, true, geneToProcess, true);
                                    ComputeMannWhitneyAndPrint(perDiseaseOutputFile, byRange[i].Where(x => x.nMutations != 0 && x.tumorType == disease).ToList(), isOne, false, geneToProcess, true);
                                    ComputeMannWhitneyAndPrint(perDiseaseOutputFile, byRange[i].Where(x => x.tumorType == disease).ToList(), isOne, false, geneToProcess, true);
                                }

                                ComputeMannWhitneyAndPrint(panCancerOutputFile, byRange[i].Where(x => x.nMutations < 2).ToList(), isZero, true, geneToProcess, true);
                                ComputeMannWhitneyAndPrint(panCancerOutputFile, byRange[i], isZero, true, geneToProcess, true);
                                ComputeMannWhitneyAndPrint(panCancerOutputFile, byRange[i].Where(x => x.nMutations != 0).ToList(), isOne, false, geneToProcess, true);
                                ComputeMannWhitneyAndPrint(panCancerOutputFile, byRange[i], isOne, false, geneToProcess, true);
                            }
                            else
                            {
                                foreach (var perDiseaseOutputFileEntry in outputFilesByDisease)
                                {
                                    perDiseaseOutputFileEntry.Value.WriteLine("\t*\t*\t*\t*");
                                }
                                panCancerOutputFile.Write("\t*\t*\t*\t*");
                            }
                        }
                    }

                    WriteCounts(panCancerOutputFile, geneToProcess, "", experimentsByParticipantId);
                    foreach (var perDiseaseOutputFileEntry in outputFilesByDisease)
                    {
                        WriteCounts(perDiseaseOutputFileEntry.Value, geneToProcess, perDiseaseOutputFileEntry.Key, experimentsByParticipantId);
                    }
                }

                foreach (var inputFile in completedInputFiles)
                {
                    inputFiles.Remove(inputFile);
                }

                nGenesProcessed++;
                var now = fiveMinuteTimer.ElapsedMilliseconds;
                if (now - lastFiveMinutePrint > 5 * 60 * 1000)
                {
                    Console.WriteLine("Processed " + nGenesProcessed + " (of which " + nGenesSkipped + " were skipped) in " + (now + 30000) / 60000 + "m");
                    lastFiveMinutePrint = now;
                }
            }
            panCancerOutputFile.Close();

            foreach (var perDiseaseOutputFileEntry in outputFilesByDisease)
            {
                perDiseaseOutputFileEntry.Value.Close();
            }
        }
    }
}
