using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using ExpressionLib;
using System.IO;

namespace NormalExpression
{
    class Program
    {
        class Result
        {
            public string Hugo_symbol;
            public int nAboveOrEqual = 0;
            public int nBelow = 0;
            public double totalMu = 0;
            public double totalZ = 0;
            public double uncorrectedP;
        }

        static void Main(string[] args)
        {
            var excludedAnalyses = ExpressionTools.LoadExcludedAnalyses();

            var tcgaRecords = ExpressionTools.LoadTCGARecords(null /* stored BAMs*/, excludedAnalyses);
            ExpressionTools.LoadTCGARecordsForLocalRealigns(tcgaRecords, null);

            Dictionary<string, ExpressionTools.Sample> allSamples;
            var participants = ExpressionTools.BuildParticipantData(tcgaRecords, out allSamples);

            var experiments = ExpressionTools.LoadExperimentsFromFile(@"f:\temp\expression\experiments.txt", participants, tcgaRecords);

            var experimentsByRNAAnalysisID = new Dictionary<string, ExpressionTools.Experiment>();
            foreach (var experiment in experiments)
            {
                experimentsByRNAAnalysisID.Add(experiment.TumorRNAAnalysis.analysis_id, experiment);
            }

            var geneScatterGraphLines = ExpressionTools.GeneScatterGraphLine.LoadAllGeneScatterGraphEntries(ExpressionTools.geneScatterGraphsDirectory, false, experimentsByRNAAnalysisID, "*");

            var scatterGraphLinesByGene = new Dictionary<string, List<ExpressionTools.GeneScatterGraphLine>>();

            foreach (var line in geneScatterGraphLines)
            {
                if (!scatterGraphLinesByGene.ContainsKey(line.Hugo_Symbol))
                {
                    scatterGraphLinesByGene.Add(line.Hugo_Symbol, new List<ExpressionTools.GeneScatterGraphLine>());
                }

                scatterGraphLinesByGene[line.Hugo_Symbol].Add(line);
            }

            var results = new List<Result>();

            foreach (var geneEntry in scatterGraphLinesByGene) 
            {
                var result = new Result();
                result.Hugo_symbol = geneEntry.Key;

                foreach (var line in geneEntry.Value) 
                {
                    var totalDNA = line.n_DNA_Matching_Tumor + line.n_DNA_Matching_Reference;
                    if (totalDNA < 10 || line.tumorDNARatio < .4 || line.tumorDNARatio > .6 || !line.zKnown) {
                        //
                        // DNA level too far off or no data.
                        //
                        continue;
                    }

                    double normalizedZNormal = line.
                    if (normalizedZNormal)
                }
            }

            var outputFile = ExpressionTools.CreateStreamWriterWithRetry(@"f:\temp\expression\normal_expression.txt");

            outputFile.WriteLine("Hugo_symbol")
        }
    }
}
