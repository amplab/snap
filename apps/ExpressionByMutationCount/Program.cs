using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using System.IO;
using ExpressionLib;
using System.Diagnostics;

namespace ExpressionByMutationCount
{
    class Program
    {
        class GeneExpressionFile {
            public GeneExpressionFile(string filename_, string tumorType_) {
                filename = filename_;
                tumorType = tumorType_;

                int headerPrefixLength = headerPrefix.Count();

                reader = new StreamReader(filename);

                var line = reader.ReadLine();

                if (line == null || line.Count() < headerPrefixLength || line.Substring(0, headerPrefixLength) != headerPrefix || line.Substring(headerPrefixLength, 2) != "1.")
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
            public StreamReader reader = null;
            public bool seenDone = false;
            public string currentHugoSymbol = null;
            public int currentNMutations = 0;
            public double[] currentZ = new double[nWidths];
            public bool[] currentZValid = new bool[nWidths];
            public const int nWidths = 20;
            const string headerPrefix = "ExpressionNearMutations v";

            public bool GetNextLine()
            {
                string line = reader.ReadLine();

                if (null == line)
                {
                    if (!seenDone)
                    {
                        Console.WriteLine("File " + filename + " is truncated.");
                    }

                    return false;
                }

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
                if (fields.Count() != nWidths + 2)
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
                    }
                }
                catch (FormatException)
                {
                    Console.WriteLine("Corrupt data line in " + filename + ": " + line);
                    return false;
                }

                return true;
            }
        }

        class ExpressionInstance : IComparer<ExpressionInstance>
        {
            public ExpressionInstance(string tumorType_, int nMutations_, double z_)
            {
                tumorType = tumorType_;
                nMutations = nMutations_;
                z = z_;
            }

            public readonly int nMutations;
            public readonly double z;
            public readonly string tumorType;

            public int Compare(ExpressionInstance a, ExpressionInstance b)
            {
                if (a.z > b.z) return 1;
                if (a.z < b.z) return -1;
                return 0;
            }
        }

        class GeneState
        {
            public GeneState(string hugo_symbol_)
            {
                hugo_symbol = hugo_symbol_;
            }

            public void AddExpression(string tumorType, int range, int nMutations, double z)
            {
                if (!expressionByRange.ContainsKey(range))
                {
                    expressionByRange.Add(range, new List<ExpressionInstance>());
                }
                expressionByRange[range].Add(new ExpressionInstance(tumorType, nMutations, z));
            }

            public string hugo_symbol;
            public Dictionary<int, List<ExpressionInstance>> expressionByRange = new Dictionary<int, List<ExpressionInstance>>();

            public int nInputFilesPastThisGene = 0;
        }

        static void ComputeMannWhitneyAndPrint(StreamWriter outputFile, List<ExpressionInstance> instances, ExpressionTools.MannWhitney<ExpressionInstance>.WhichGroup selector, bool twoTailed)
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
                p = ExpressionTools.MannWhitney<ExpressionInstance>.ComputeMannWhitney(
                        instances,
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

        static void WriteFileHeader(StreamWriter outputFile)
        {
            outputFile.Write("Hugo Symbol");
            int width = 0;
            for (int i = 0; i < GeneExpressionFile.nWidths; i++)
            {
                outputFile.Write("\t" + width + "Kbp 0 vs. 1\t" + width + "Kbp 0 vs. not zero\t" + width + "Kbp 1.vs many\t" + width + "Kbp 1 vs. not 1");

                if (0 == width)
                {
                    width = 1;
                }
                else
                {
                    width *= 2;
                }
            }
            outputFile.WriteLine();
        }


        static void Main(string[] args)
        {
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

            timer.Stop();
            Console.WriteLine("Loaded " + experiments.Count() + " experiments in " + (timer.ElapsedMilliseconds + 500) / 1000 + "s, of which " + 
                experiments.Where(e => e.TumorRNAAnalysis == null || e.TumorRNAAnalysis.geneExpressionFileName == null || e.TumorRNAAnalysis.geneExpressionFileName == "").Count() + " are missing gene expression.");

            var inputFiles = new List<GeneExpressionFile>();

            timer.Reset();
            timer.Start();
            foreach (var experiment in experiments)
            {
                if (experiment.TumorRNAAnalysis != null && experiment.TumorRNAAnalysis.geneExpressionFileName != null && experiment.TumorRNAAnalysis.geneExpressionFileName != "")
                {
                    var file = new GeneExpressionFile(experiment.TumorRNAAnalysis.geneExpressionFileName, experiment.disease_abbr);
  
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

            string baseFileName = @"f:\temp\expression\ExpressionNearMutations";

            var panCancerOutputFile = new StreamWriter(baseFileName + ".txt");
            WriteFileHeader(panCancerOutputFile);
 
            var outputFilesByDisease = new Dictionary<string, StreamWriter>();

            var listOfDiseases = ExpressionTools.GetListOfDiseases(experiments);

            foreach (var disease in listOfDiseases)
            {
                outputFilesByDisease.Add(disease, new StreamWriter(baseFileName + "_" + disease + ".txt"));
                WriteFileHeader(outputFilesByDisease[disease]);
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
                foreach (var inputFile in inputFiles.Where(f => f.currentHugoSymbol == nextHugoSymbol))
                {
                    for (int i = 0 ; i < GeneExpressionFile.nWidths; i++) {
                        if (inputFile.currentZValid[i])
                        {
                            geneToProcess.AddExpression(inputFile.tumorType, i, inputFile.currentNMutations, inputFile.currentZ[i]);
                        }
                    }

                    if (!inputFile.GetNextLine())
                    {
                        completedInputFiles.Add(inputFile);
                    }
                 }

                //
                // Compute and write out the results.
                //

                panCancerOutputFile.Write(nextHugoSymbol);
                for (int i = 0; i < GeneExpressionFile.nWidths; i++)
                {
                    if (geneToProcess.expressionByRange.ContainsKey(i))
                    {
                        foreach (var perDiseaseOutputFileEntry in outputFilesByDisease) {
                            var perDiseaseOutputFile = perDiseaseOutputFileEntry.Value;
                            var disease = perDiseaseOutputFileEntry.Key;

                            ComputeMannWhitneyAndPrint(perDiseaseOutputFile, geneToProcess.expressionByRange[i].Where(x => x.nMutations < 2 && x.tumorType == disease).ToList(), isZero, true);
                            ComputeMannWhitneyAndPrint(perDiseaseOutputFile, geneToProcess.expressionByRange[i].Where(x => x.tumorType == disease).ToList(), isZero, true);
                            ComputeMannWhitneyAndPrint(perDiseaseOutputFile, geneToProcess.expressionByRange[i].Where(x => x.nMutations != 0 && x.tumorType == disease).ToList(), isOne, false);
                            ComputeMannWhitneyAndPrint(perDiseaseOutputFile, geneToProcess.expressionByRange[i].Where(x => x.tumorType == disease).ToList(), isNotOne, false);
                        }

                        ComputeMannWhitneyAndPrint(panCancerOutputFile, geneToProcess.expressionByRange[i].Where(x => x.nMutations < 2).ToList(), isZero, true);
                        ComputeMannWhitneyAndPrint(panCancerOutputFile, geneToProcess.expressionByRange[i], isZero, true);
                        ComputeMannWhitneyAndPrint(panCancerOutputFile, geneToProcess.expressionByRange[i].Where(x => x.nMutations != 0).ToList(), isOne, false);
                        ComputeMannWhitneyAndPrint(panCancerOutputFile, geneToProcess.expressionByRange[i], isNotOne, false);
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

                panCancerOutputFile.WriteLine();

                foreach (var inputFile in completedInputFiles)
                {
                    inputFiles.Remove(inputFile);
                }

                nGenesProcessed++;
                var now = fiveMinuteTimer.ElapsedMilliseconds;
                if (now - lastFiveMinutePrint > 5 * 60 * 1000)
                {
                    Console.WriteLine("Processed " + nGenesProcessed + " in " + (now + 30000) / 60000 + "m");
                    lastFiveMinutePrint = now;
                }
            }
            panCancerOutputFile.Close();
        }
    }
}
