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

        class ExpressionInstance
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

            public static int Compare(ExpressionInstance a, ExpressionInstance b)
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

            var outputFile = new StreamWriter(@"f:\temp\expression\ExpressionNearMutations.txt");
            outputFile.Write("Hugo Symbol");
            int width = 0;
            for (int i = 0; i < GeneExpressionFile.nWidths; i++)
            {
                outputFile.Write("\t" + width + "Kbp 0 vs. 1\t" + width + "Kbp 0 vs. many\t" + width + "Kbp 1.vs many\t");

                if (0 == width) {
                    width = 1;
                } else {
                    width *= 2;
                }
            }


            ExpressionTools.MannWhitney<ExpressionInstance>.GetValue getValue = new ExpressionTools.MannWhitney<ExpressionInstance>.GetValue(x => x.z);
            ExpressionTools.MannWhitney<ExpressionInstance>.WhichGroup isZero = new ExpressionTools.MannWhitney<ExpressionInstance>.WhichGroup(x => x.nMutations == 0);
            ExpressionTools.MannWhitney<ExpressionInstance>.WhichGroup isOne = new ExpressionTools.MannWhitney<ExpressionInstance>.WhichGroup(x => x.nMutations == 1);
            ExpressionTools.MannWhitney<ExpressionInstance>.WhichGroup isMoreThanOne = new ExpressionTools.MannWhitney<ExpressionInstance>.WhichGroup(x => x.nMutations > 1);

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

                var completedInputFiles = new List<GeneExpressionFile>();
                var geneToProcess = new GeneState(nextHugoSymbol);
                foreach (var inputFile in inputFiles.Where(f => f.currentHugoSymbol == nextHugoSymbol))
                {
                    for (int i = 0 ; i < GeneExpressionFile.nWidths; i++) {
                        if (inputFile.currentZValid[i])
                        {
                            geneToProcess.AddExpression(inputFile.tumorType, i - 2, inputFile.currentNMutations, inputFile.currentZ[i]);
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

                for (int i = 0; i < GeneExpressionFile.nWidths; i++)
                {
                    bool enoughData;
                    bool reversed;
                    int nFirstGroup;
                    int nSecondGroup;
                    double U;
                    double z;
                    ExpressionTools.MannWhitney<ExpressionInstance>.ComputeMannWhitney(geneToProcess.expressionByRange[i].Where(x => x.nMutations < 2), ExpressionInstance.Compare, isOne, getValue, out enoughData, out reversed, out nFirstGroup, out nSecondGroup, out U, out z);
                }

                foreach (var inputFile in completedInputFiles)
                {
                    inputFiles.Remove(inputFile);
                }

                outputFile.Close();
            }
        }
    }
}
