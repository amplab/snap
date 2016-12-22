using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using System.IO;
using ExpressionLib;

namespace ClassifyGenes
{
    class Program
    {
        static void Main(string[] args)
        {
            List<string> excludedAnalyses;
            Dictionary<string, ExpressionTools.TCGARecord> tcgaRecords;
            Dictionary<string, string> sampleToParticipantIDMap;
            Dictionary<string, ExpressionTools.Participant> participants;
            List<ExpressionTools.Experiment> experiments;
            Dictionary<string, ExpressionTools.Sample> allSamples;

            ExpressionTools.LoadStateFromExperimentsFile(out excludedAnalyses, out tcgaRecords, out sampleToParticipantIDMap, out participants, out experiments, out allSamples);

            var geneMap = new ExpressionTools.GeneLocationsByNameAndChromosome(ExpressionTools.LoadGeneLocationInfo(@"\\gcr\scratch\b99\bolosky\knownGene-hg19.txt", @"\\gcr\scratch\b99\bolosky\kgXref-hg19.txt"));


        }
    }
}
