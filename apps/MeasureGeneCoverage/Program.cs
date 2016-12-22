using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using System.Threading;
using System.IO;
using ExpressionLib;

namespace MeasureGeneCoverage
{
    class Program
    {
        static Dictionary<string, ExpressionTools.GeneLocationsByNameAndChromosome> geneLocationInformation;    // Map of the genes for each reference class (hg18 or hg19)
        static List<string> queuedParticipantIDs;
        static Dictionary<string, ExpressionTools.Experiment> experimentsByParticipantId = new Dictionary<string, ExpressionTools.Experiment>();
        static Dictionary<string, ExpressionTools.ExonicMap> exonicMaps = new Dictionary<string,ExpressionTools.ExonicMap>();

        static void ProcessParticipants()
        {
            while (true)
            {
                string participantId;
                lock (queuedParticipantIDs)
                {
                    if (queuedParticipantIDs.Count() == 0)
                    {
                        return;
                    }

                    participantId = queuedParticipantIDs[0];
                    queuedParticipantIDs.RemoveAt(0);
                }

                if (!experimentsByParticipantId.ContainsKey(participantId))
                {
                    Console.WriteLine("Can't find experiment for participant ID " + participantId);
                    continue;
                }

                var experiment = experimentsByParticipantId[participantId];

                if (experiment.TumorDNAAnalysis.allcountFileName == null || experiment.TumorDNAAnalysis.allcountFileName == "")
                {
                    Console.WriteLine("participant " + participantId + " doesn't have a tumor DNA allcount file");
                    continue;
                }

                var reader = ExpressionTools.CreateCompressedStreamReaderWithRetry(experiment.TumorDNAAnalysis.allcountFileName);

                var referenceClass = (experiment.TumorDNAAnalysis.refassemShortName == "hg18" || experiment.TumorDNAAnalysis.refassemShortName == "hg18_broad_variant") ? "hg18" : "hg19";
                var exonMap = exonicMaps[referenceClass];
                var geneMap = geneLocationInformation[referenceClass];





                reader.Close();

            }
        }

        static void Main(string[] args)
        {
            if (args.Count() == 0)
            {
                Console.WriteLine("usage: MeasureGeneCoverage <participant id>");
                return;
            }

            List<string> excludedAnalyses;
            Dictionary<string, ExpressionTools.Participant> participants;
            Dictionary<string, ExpressionTools.TCGARecord> tcgaRecords;
            Dictionary<string, string> sampleToParticipantIDMap;

            List<ExpressionTools.Experiment> experiments;
            Dictionary<string, ExpressionTools.Sample> allSamples;

            ExpressionTools.LoadStateFromExperimentsFile(out excludedAnalyses, out tcgaRecords, out sampleToParticipantIDMap, out participants, out experiments, out allSamples);

            foreach (var experiment in experiments)
            {
                experimentsByParticipantId.Add(experiment.participant.participantId, experiment);
            }

            geneLocationInformation = new Dictionary<string, ExpressionTools.GeneLocationsByNameAndChromosome>(); // // Maps reference class (hg18 or hg19) to map from hugo symbol to location info
            geneLocationInformation.Add("hg18", new ExpressionTools.GeneLocationsByNameAndChromosome(ExpressionTools.LoadGeneLocationInfo(@"\\gcr\scratch\b99\bolosky\knownGene-hg18.txt", @"\\gcr\scratch\b99\bolosky\kgXref-hg18.txt")));
            geneLocationInformation.Add("hg19", new ExpressionTools.GeneLocationsByNameAndChromosome(ExpressionTools.LoadGeneLocationInfo(@"\\gcr\scratch\b99\bolosky\knownGene-hg19.txt", @"\\gcr\scratch\b99\bolosky\kgXref-hg19.txt")));

            foreach (var entry in geneLocationInformation)
            {
                exonicMaps.Add(entry.Key, new ExpressionTools.ExonicMap(entry.Value));  // Make an exoic map for each reference class
            }

            queuedParticipantIDs = args.ToList();

            var threads = new List<Thread>();
            for (int i = 0; i < Environment.ProcessorCount/* 1*/; i++)
            {
                threads.Add(new Thread(() => ProcessParticipants()));
            }

            threads.ForEach(t => t.Start());
            threads.ForEach(t => t.Join());
        }
    }
}
