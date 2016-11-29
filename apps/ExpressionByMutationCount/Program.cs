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

                reader = ExpressionTools.CreateStreamReaderWithRetry(filename);

                var line = reader.ReadLine();

                if (line == null || line.Count() < headerPrefixLength || line.Substring(0, headerPrefixLength) != headerPrefix || line.Substring(headerPrefixLength, 4) != "3.0 ") // v 3.0 adds whole autosome
                {
                    Console.WriteLine("Truncated corrupt or wrong version gene expression file " + filename);
                    seenDone = true;
                    return;
                }

                //
                // Now read the second header line (the one with the column headers).
                //
                line = reader.ReadLine();
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
            public double[] currentZExclusive = new double[nWidths];
            public bool[] currentZExclusiveValid = new bool[nWidths];
            public double[] currentMuExclusive = new double[nWidths];
            public bool[] currentMuExclusiveValid = new bool[nWidths];
            public double[] perChromosomeZ = new double[ExpressionTools.nHumanNuclearChromosomes];
            public bool[] perChromosomeZValid = new bool[ExpressionTools.nHumanNuclearChromosomes];
            public double[] perChromosomeMu = new double[ExpressionTools.nHumanNuclearChromosomes];
            public bool[] perChromosomeMuValid = new bool[ExpressionTools.nHumanNuclearChromosomes];
            public double wholeAutosomeZ = 0;
            public bool wholeAutosomeZValid = false;
            public double wholeAutosomeZExclusive = 0;
            public bool wholeAutosomeZExclusiveValid = false; 
            public double wholeAutosomeMu = 0;
            public bool wholeAutosomeMuValid = false;
            public double wholeAutosomeMuExclusive = 0;
            public bool wholeAutosomeMuExclusiveValid = false;
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
                if (fields.Count() != 2 * (forAlleleSpecificExpression ? 1 : 2) * (nWidths + 1) + 2 + 
                    (forAlleleSpecificExpression ? 1 : 2) * ExpressionTools.nHumanNuclearChromosomes)
                {
                    Console.WriteLine("File " + filename + " has unparsable line " + line);
                    return false;
                }

                currentHugoSymbol = ExpressionTools.ConvertToNonExcelString(fields[0]);

                try
                {
                    currentNMutations = Convert.ToInt32(fields[1]);
                    int exclusiveOffset = forAlleleSpecificExpression ? (nWidths + 1) : 2 * (nWidths + 1);  // +1 on nWidths is because of the whole autosome column
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

                        if (fields[i + exclusiveOffset + 2] == "*")
                        {
                            currentZExclusiveValid[i] = false;
                        }
                        else
                        {
                            currentZExclusive[i] = Convert.ToDouble(fields[i + exclusiveOffset + 2]);
                            currentZExclusiveValid[i] = true;
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

                        if (forAlleleSpecificExpression || fields[i + nWidths + 1 + exclusiveOffset + 2] == "*")
                        {
                            currentMuExclusiveValid[i] = false;
                        }
                        else
                        {
                            currentMuExclusive[i] = Convert.ToDouble(fields[i + nWidths + 1 + exclusiveOffset + 2]);
                            currentMuExclusiveValid[i] = true;
                        }
                    }

                    if (fields[nWidths + 2] == "*")
                    {
                        wholeAutosomeZValid = false;
                    }
                    else
                    {
                        wholeAutosomeZ = Convert.ToDouble(fields[nWidths + 2]);
                        wholeAutosomeZValid = true;
                    }

                    if (forAlleleSpecificExpression || fields[2 * nWidths + 4] == "*")
                    {
                        wholeAutosomeMuValid = false;
                    }
                    else
                    {
                        wholeAutosomeMu = Convert.ToDouble(fields[2 * nWidths + 4]);
                        wholeAutosomeMuValid = true;
                    }

                    if (fields[nWidths + 2 + exclusiveOffset] == "*")
                    {
                        wholeAutosomeZExclusiveValid = false;
                    }
                    else
                    {
                        wholeAutosomeZExclusiveValid = true;
                        wholeAutosomeZExclusive = Convert.ToDouble(fields[nWidths + 2 + exclusiveOffset]);
                    }

                    if (forAlleleSpecificExpression || fields[2 * nWidths + exclusiveOffset + 3] == "*")
                    {
                        wholeAutosomeMuExclusiveValid = false;
                    }
                    else
                    {
                        wholeAutosomeMuExclusive = Convert.ToDouble(fields[2 * nWidths + exclusiveOffset + 3]);
                        wholeAutosomeMuExclusiveValid = true;
                    }

                    //
                    // Now do the per-chromosome ones.  There is no exclusive, but there is a mu iff !forAlleleSpecificExpression.
                    //
                    int perChromosomeOffset = forAlleleSpecificExpression ? nWidths + 2 + exclusiveOffset + 1 : 2 * nWidths + exclusiveOffset + 4;
                    int perChromosomeMuOffset = perChromosomeOffset + ExpressionTools.nHumanNuclearChromosomes;
                    for (int whichChromosome = 0; whichChromosome < ExpressionTools.nHumanNuclearChromosomes; whichChromosome++)
                    {
                        if (fields[whichChromosome + perChromosomeOffset] == "*")
                        {
                            perChromosomeZValid[whichChromosome] = false;
                        }
                        else
                        {
                            perChromosomeZValid[whichChromosome] = true;
                            perChromosomeZ[whichChromosome] = Convert.ToDouble(fields[whichChromosome + perChromosomeOffset]);
                        }

                        if (forAlleleSpecificExpression || fields[whichChromosome + perChromosomeMuOffset] == "*")
                        {
                            perChromosomeMuValid[whichChromosome] = false;
                        }
                        else
                        {
                            perChromosomeMuValid[whichChromosome] = true;
                            perChromosomeMu[whichChromosome] = Convert.ToDouble(fields[whichChromosome + perChromosomeMuOffset]);
                        }
                    } // For each chromosome
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

                for (int i = 0; i < ExpressionTools.nHumanNuclearChromosomes; i++)
                {
                    perChromosomeExpression[i] = new List<ExpressionInstance>();
                    perChromosomeMeanExpression[i] = new List<ExpressionInstance>();
                }
            }

            public void AddExpression(string tumorType, int range, int nMutations, double z /* or mu, as approproate*/, string participantID, bool usingMu, bool exclusive)
            {
                Dictionary<int, List<ExpressionInstance>> byRange;
                if (usingMu) {
                    if (exclusive)
                    {
                        byRange = exclusiveMeanExpressionByRange;
                    }
                    else
                    {
                        byRange = meanExpressionByRange;
                    }
                } else {
                    if (exclusive) {
                        byRange = exclusiveExpressionByRange;
                    }
                    else
                    {
                        byRange = expressionByRange;
                    }
                }

                if (!byRange.ContainsKey(range))
                {
                    byRange.Add(range, new List<ExpressionInstance>());
                }
                byRange[range].Add(new ExpressionInstance(tumorType, nMutations, z, participantID));
            }

            public void AddWholeAutosomeExpression(string tumorType, int nMutations, double z /* or mu, as approproate*/, string participantID, bool usingMu, bool exclusive)
            {
                var newEi = new ExpressionInstance(tumorType, nMutations, z, participantID);
                if (usingMu)
                {
                    if (exclusive)
                    {
                        wholeAutosomeMeanExpressionExclusive.Add(newEi);
                    }
                    else
                    {
                        wholeAutosomeMeanExpression.Add(newEi);
                    }
                }
                else
                {
                    if (exclusive)
                    {
                        wholeAutosomeExpressionExclusive.Add(newEi);
                    }
                    else
                    {
                        wholeAutosomeExpression.Add(newEi);
                    }
                }
            }

            public void AddPerChromosomeExpression(string tumorType, int nMutations, double z /*or mu as appropriate*/, string participantID, bool usingMu, int whichChromosome)
            {
                var newEi = new ExpressionInstance(tumorType, nMutations, z, participantID);

                if (usingMu)
                {
                    perChromosomeMeanExpression[whichChromosome].Add(newEi);
                }
                else
                {
                    perChromosomeExpression[whichChromosome].Add(newEi);
                }
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
                        reader = ExpressionTools.CreateStreamReaderWithRetry(scatterPlotFilename);
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
            public Dictionary<int, List<ExpressionInstance>> exclusiveExpressionByRange = new Dictionary<int, List<ExpressionInstance>>();
            public Dictionary<int, List<ExpressionInstance>> exclusiveMeanExpressionByRange = new Dictionary<int, List<ExpressionInstance>>();   // i.e., using mu not z
            public List<ExpressionInstance> wholeAutosomeExpression = new List<ExpressionInstance>();
            public List<ExpressionInstance> wholeAutosomeExpressionExclusive = new List<ExpressionInstance>();
            public List<ExpressionInstance> wholeAutosomeMeanExpression = new List<ExpressionInstance>();
            public List<ExpressionInstance> wholeAutosomeMeanExpressionExclusive = new List<ExpressionInstance>();
            public List<ExpressionInstance>[] perChromosomeExpression = new List<ExpressionInstance>[ExpressionTools.nHumanNuclearChromosomes];
            public List<ExpressionInstance>[] perChromosomeMeanExpression = new List<ExpressionInstance>[ExpressionTools.nHumanNuclearChromosomes];
            public Dictionary<string, Mutation> perExperimentState = new Dictionary<string, Mutation>();

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
                if (geneState.perExperimentState.ContainsKey(instance.participantID))
                {
                    var perExperimentState = geneState.perExperimentState[instance.participantID];
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

                var instancesToUse = FilterUnusualMutantsIfRequested(instances, geneState, excludeUnusualMutants);

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

        static void ComputeNMeanAndStandardDeviationAndPrint(StreamWriter outputFile, List<ExpressionInstance> instances, GeneState geneState, bool excludeUnusualMutants)
        {
            var values = FilterUnusualMutantsIfRequested(instances, geneState, excludeUnusualMutants).ConvertAll(ExtractZFromExpressionInstance);

            outputFile.Write("\t" + values.Count());
            if (values.Count() == 0)
            {
                outputFile.Write("\t*\t*");
                return;
            }

            outputFile.Write("\t" + ExpressionTools.MeanOfList(values) + "\t" + ExpressionTools.StandardDeviationOfList(values));
        }

        static ExpressionTools.MannWhitney<ExpressionInstance>.GetValue getValue = new ExpressionTools.MannWhitney<ExpressionInstance>.GetValue(x => x.z);

        static void WriteHeaderGroup(string regionName, StreamWriter outputFile, string muString) // Writes out the header for the complete set of measurements for one region
        {
            outputFile.Write(/*"\t" + groupId + " 0 vs. 1" + muString + "\t" + groupId + " 0 vs. not zero" + muString +*/ "\t" + regionName + "Kbp 1 vs. many" + muString + "\t" + regionName + " 1 vs. not 1" + muString);
            string[] mutationSets = { "0", "1", ">1" };

            for (int i = 0; i < mutationSets.Count(); i++)
            {
                outputFile.Write("\t" + regionName + " " + mutationSets[i] + " mutation" + muString + " N\t" + regionName + " " + mutationSets[i] + " mutation " + muString + " mean\t" + regionName + " " + mutationSets[i] + " mutation " + muString + " stdDev");
            }
        }

        static void WriteFileHeader(StreamWriter outputFile, bool forAlleleSpecificExpression)
        {
            outputFile.Write("Hugo Symbol");
            for (int exclusive = 0; exclusive < 2; exclusive++)
            {
                for (int mu = 0; mu < (forAlleleSpecificExpression ? 1 : 2); mu++)
                {
                    int width = 0;
                    string muString = ((mu == 0) ? "" : " mu") + (exclusive == 0 ? "" : " exclusive");
                    for (int i = 0; i < GeneExpressionFile.nWidths; i++)
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

            for (int mu = 0; mu < (forAlleleSpecificExpression ? 1 : 2); mu++)
            {
                string muString = (mu == 0) ? "" : " mu";

                for (int whichChromosome = 0; whichChromosome < ExpressionTools.nHumanNuclearChromosomes; whichChromosome++)
                {
                    var name = ExpressionTools.ChromosomeIndexToName(whichChromosome, true);

                    WriteHeaderGroup(ExpressionTools.ChromosomeIndexToName(whichChromosome, true), outputFile, muString);
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

            foreach (var expression in geneToProcess.wholeAutosomeExpression)
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

        static ExpressionTools.MannWhitney<ExpressionInstance>.WhichGroup isZero = new ExpressionTools.MannWhitney<ExpressionInstance>.WhichGroup(x => x.nMutations == 0);
        static ExpressionTools.MannWhitney<ExpressionInstance>.WhichGroup isOne = new ExpressionTools.MannWhitney<ExpressionInstance>.WhichGroup(x => x.nMutations == 1);
        static ExpressionTools.MannWhitney<ExpressionInstance>.WhichGroup isNotOne = new ExpressionTools.MannWhitney<ExpressionInstance>.WhichGroup(x => x.nMutations != 1);
        static ExpressionTools.MannWhitney<ExpressionInstance>.WhichGroup isMoreThanOne = new ExpressionTools.MannWhitney<ExpressionInstance>.WhichGroup(x => x.nMutations > 1);

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
            int missingCount;
            if (!forAlleleSpecificExpression) {
                missingCount = experiments.Where(e => e.TumorRNAAnalysis == null || e.TumorRNAAnalysis.geneExpressionFileName == null || e.TumorRNAAnalysis.geneExpressionFileName == "").Count();
            } else {
                missingCount = experiments.Where(e => e.NormalDNAAnalysis == null || e.NormalDNAAnalysis.alleleSpecificGeneExpressionFileName == null || e.NormalDNAAnalysis.alleleSpecificGeneExpressionFileName == "").Count();
            }

            Console.WriteLine("Loaded " + experiments.Count() + " experiments in " + (timer.ElapsedMilliseconds + 500) / 1000 + "s, of which " + missingCount + " are missing gene expression.");
            
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

            var panCancerOutputFile = ExpressionTools.CreateStreamWriterWithRetry(baseFileName + ".txt");
            WriteFileHeader(panCancerOutputFile, forAlleleSpecificExpression);
 
            var outputFilesByDisease = new Dictionary<string, StreamWriter>();

            var listOfDiseases = ExpressionTools.GetListOfDiseases(experiments);

            foreach (var disease in listOfDiseases)
            {
                outputFilesByDisease.Add(disease, ExpressionTools.CreateStreamWriterWithRetry(baseFileName + "_" + disease + ".txt"));
                WriteFileHeader(outputFilesByDisease[disease], forAlleleSpecificExpression);
            }
            
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

                if (!geneToProcess.loadPerExperimentState(sampleToParticipantIDMap, experimentsByRNAAnalysisID))
                {
                    //
                    // There wasn't a scatter graphs file for it, probably because it had too few mutations.  Skip it.
                    //
                    nGenesSkipped++;

                    //
                    // Eat up the lines for the skipped genes, and remember any files that have hit EOF.
                    //
                    foreach (var inputFile in inputFiles.Where(f => f.currentHugoSymbol == nextHugoSymbol))
                    {
                        if (!inputFile.GetNextLine())
                        {
                            completedInputFiles.Add(inputFile);
                        }
                    }

                }
                else
                {
                    foreach (var inputFile in inputFiles.Where(f => f.currentHugoSymbol == nextHugoSymbol))
                    {
                        perGeneLines.Add(inputFile.participantID + "\t" + experimentsByParticipantId[inputFile.participantID].disease_abbr + "\t" + inputFile.currentLine);

                        for (int i = 0 ; i < GeneExpressionFile.nWidths; i++) {
                            if (inputFile.currentZValid[i])
                            {
                                geneToProcess.AddExpression(inputFile.tumorType, i, inputFile.currentNMutations, inputFile.currentZ[i], inputFile.participantID, false, false);
                            }

                            if (inputFile.currentZExclusiveValid[i])
                            {
                                geneToProcess.AddExpression(inputFile.tumorType, i, inputFile.currentNMutations, inputFile.currentZExclusive[i], inputFile.participantID, false, true);
                            }

                            if (!forAlleleSpecificExpression && inputFile.currentMuValid[i]) {
                                geneToProcess.AddExpression(inputFile.tumorType, i, inputFile.currentNMutations, inputFile.currentMu[i], inputFile.participantID, true, false);
                            }

                            if (!forAlleleSpecificExpression && inputFile.currentMuExclusiveValid[i])
                            {
                                geneToProcess.AddExpression(inputFile.tumorType, i, inputFile.currentNMutations, inputFile.currentMuExclusive[i], inputFile.participantID, true, true);
                            }
                        }

                        if (inputFile.wholeAutosomeZValid) {
                            geneToProcess.AddWholeAutosomeExpression(inputFile.tumorType, inputFile.currentNMutations, inputFile.wholeAutosomeZ, inputFile.participantID, false, false);
                        }

                        if (inputFile.wholeAutosomeZExclusiveValid) {
                            geneToProcess.AddWholeAutosomeExpression(inputFile.tumorType, inputFile.currentNMutations, inputFile.wholeAutosomeZExclusive, inputFile.participantID, false, true);
                        }
                        
                        if (inputFile.wholeAutosomeMuValid) {
                            geneToProcess.AddWholeAutosomeExpression(inputFile.tumorType, inputFile.currentNMutations, inputFile.wholeAutosomeMu, inputFile.participantID, true, false);
                        }

                        if (inputFile.wholeAutosomeMuExclusiveValid) {
                            geneToProcess.AddWholeAutosomeExpression(inputFile.tumorType, inputFile.currentNMutations, inputFile.wholeAutosomeMuExclusive, inputFile.participantID, true, true);
                        }

                        for (int whichChromosome = 0; whichChromosome < ExpressionTools.nHumanNuclearChromosomes; whichChromosome++)
                        {
                            if (inputFile.perChromosomeZValid[whichChromosome]) 
                            {
                                geneToProcess.AddPerChromosomeExpression(inputFile.tumorType, inputFile.currentNMutations, inputFile.perChromosomeZ[whichChromosome], inputFile.participantID, false, whichChromosome);
                            }                            
                            
                            if (inputFile.perChromosomeMuValid[whichChromosome]) 
                            {
                                geneToProcess.AddPerChromosomeExpression(inputFile.tumorType, inputFile.currentNMutations, inputFile.perChromosomeMu[whichChromosome], inputFile.participantID, true, whichChromosome);
                            }
                        }

                        if (!inputFile.GetNextLine())
                        {
                            completedInputFiles.Add(inputFile);
                        }
                    }

                    var perGeneLinesFile = ExpressionTools.CreateStreamWriterWithRetry(@"f:\temp\expression\RegionalExpressionByGene\" + nextHugoSymbol.ToLower() + (forAlleleSpecificExpression ? "_allele_specific" : "") + "_lines.txt");
                    perGeneLinesFile.Write("ParticipantID\tdisease abbr.\tHugo Symbol\tMutation Count");
                    for (int exclusive = 0; exclusive < 2; exclusive++)
                    {
                        ulong width = 1000;
                        perGeneLinesFile.Write("\t0" + (exclusive == 1 ? " exclusive" : ""));

                        for (int i = 1; i < GeneExpressionFile.nWidths; i++)
                        {
                            perGeneLinesFile.Write("\t" + ExpressionTools.SizeToUnits(width) + "b" + ((exclusive == 1) ? " exclusive" : ""));
                            width *= 2;
                        }

                        perGeneLinesFile.Write("\tWhole Autosome" + ((exclusive == 1) ? " exclusive" : ""));

                        if (!forAlleleSpecificExpression)
                        {
                            for (int i = 1; i < GeneExpressionFile.nWidths; i++)
                            {
                                perGeneLinesFile.Write("\t" + ExpressionTools.SizeToUnits(width) + "b mean" + ((exclusive == 1) ? " exclusive" : ""));
                                width *= 2;
                            }

                            perGeneLinesFile.Write("\tWhole Autosome mean " + ((exclusive == 1) ? " exclusive" : ""));
                        }
                    }

                    for (int whichChromosome = 0; whichChromosome < ExpressionTools.nHumanNuclearChromosomes; whichChromosome++)
                    {
                        perGeneLinesFile.Write("\t" + ExpressionTools.ChromosomeIndexToName(whichChromosome, true));
                    }

                    if (!forAlleleSpecificExpression) 
                    {
                        for (int whichChromosome = 0; whichChromosome < ExpressionTools.nHumanNuclearChromosomes; whichChromosome++)
                        {
                            perGeneLinesFile.Write("\t" + ExpressionTools.ChromosomeIndexToName(whichChromosome, true) + " mean");
                        }
                    }

                    perGeneLinesFile.WriteLine();

                    foreach (var line in perGeneLines)
                    {
                        perGeneLinesFile.WriteLine(line);
                    }
                    perGeneLinesFile.Close();


                    //
                    // Compute and write out the results.
                    //

                    panCancerOutputFile.Write(ExpressionTools.ConvertToExcelString(nextHugoSymbol));
                    foreach (var perDiseaseOutputFileEntry in outputFilesByDisease)
                    {
                        perDiseaseOutputFileEntry.Value.Write(ExpressionTools.ConvertToExcelString(nextHugoSymbol));
                    }

                    for (var exclusive = 0; exclusive < 2; exclusive++)
                    {
                        for (int mu = 0; mu < (forAlleleSpecificExpression ? 1 : 2); mu++)  // shouldn't we be able to do this for a bool somehow?
                        {
                            Dictionary<int, List<ExpressionInstance>> byRange;
                            List<ExpressionInstance> wholeAutosome;
                            if (mu == 0)
                            {
                                if (exclusive == 0)
                                {
                                    byRange = geneToProcess.expressionByRange;
                                    wholeAutosome = geneToProcess.wholeAutosomeExpression;
                                }
                                else
                                {
                                    byRange = geneToProcess.exclusiveExpressionByRange;
                                    wholeAutosome = geneToProcess.wholeAutosomeExpressionExclusive;
                                }
                            }
                            else
                            {
                                if (exclusive == 0)
                                {
                                    byRange = geneToProcess.meanExpressionByRange;
                                    wholeAutosome = geneToProcess.wholeAutosomeMeanExpression;
                                }
                                else
                                {
                                    byRange = geneToProcess.exclusiveMeanExpressionByRange;
                                    wholeAutosome = geneToProcess.wholeAutosomeMeanExpressionExclusive;
                                }
                            }

                            for (int i = 0; i < GeneExpressionFile.nWidths; i++)
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

                    for (int whichChromosome = 0; whichChromosome < ExpressionTools.nHumanNuclearChromosomes; whichChromosome++)
                    {
                        WriteMannWhitneyToFiles(panCancerOutputFile, outputFilesByDisease, geneToProcess.perChromosomeExpression[whichChromosome],
                            geneToProcess, geneToProcess.perChromosomeExpression[whichChromosome].Count() > 0);
                    }

                    if (!forAlleleSpecificExpression)
                    {
                        for (int whichChromosome = 0; whichChromosome < ExpressionTools.nHumanNuclearChromosomes; whichChromosome++)
                        {
                            WriteMannWhitneyToFiles(panCancerOutputFile, outputFilesByDisease, geneToProcess.perChromosomeMeanExpression[whichChromosome],
                                geneToProcess, geneToProcess.perChromosomeExpression[whichChromosome].Count() > 0);
                        }
                    }

                    WriteCounts(panCancerOutputFile, geneToProcess, "", experimentsByParticipantId);
                    foreach (var perDiseaseOutputFileEntry in outputFilesByDisease)
                    {
                        WriteCounts(perDiseaseOutputFileEntry.Value, geneToProcess, perDiseaseOutputFileEntry.Key, experimentsByParticipantId);
                    }
                } // If we're not skipping this gene

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
