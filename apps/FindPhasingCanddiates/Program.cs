using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using ExpressionLib;
using System.IO;
using System.Diagnostics;
using System.Threading;

namespace FindPhasingCanddiates
{
    class Program
    {
        static StreamWriter outputFile;
        static List<ExpressionTools.Experiment> experiments;
        static List<ExpressionTools.AnnotatedSelectedVariantLine> annotatedSelectedVariants;
        static Dictionary<string, List<ExpressionTools.GeneScatterGraphLine>> scatterGraphEntriesByParticipant;
        static int nExperimentsProcessed = 0;

        static int nMissingAnnotatedSelectedVariants = 0;

        static void WorkerThread()
        {
            int nUnreportedExperimentsProcessed = 0;
            var timer = new Stopwatch();
            timer.Start();

            while (true)
            {
                ExpressionTools.Experiment experiment;

                lock (experiments)
                {
                    if (experiments.Count() == 0)
                    {
                        break;
                    }

                    experiment = experiments[0];
                    experiments.RemoveAt(0);
                }

                nUnreportedExperimentsProcessed++;

                if (timer.ElapsedMilliseconds >= 10000)
                {
                    lock (experiments)
                    {
                        nExperimentsProcessed += nUnreportedExperimentsProcessed;
                    }

                    nUnreportedExperimentsProcessed = 0;
                    timer.Restart();
                }

                if (experiment.NormalDNAAnalysis.annotatedSelectedVariantsFileName == null || experiment.NormalDNAAnalysis.annotatedSelectedVariantsFileName == "")
                {
                    lock (experiments)
                    {
                        nMissingAnnotatedSelectedVariants++;
                    }
                    continue;
                }

                annotatedSelectedVariants = ExpressionTools.AnnotatedSelectedVariantLine.readFile(experiment.NormalDNAAnalysis.annotatedSelectedVariantsFileName);
                if (null == annotatedSelectedVariants)
                {
                    Console.WriteLine("Error reading annotated selected variants file " + experiment.NormalDNAAnalysis.annotatedSelectedVariantsFileName);
                    nMissingAnnotatedSelectedVariants++;
                    continue;
                }

                foreach (var variant in annotatedSelectedVariants)
                {
                    if (variant.nMatchingReferenceRNA + variant.nMatchingVariantRNA < 10 || variant.nMatchingReferenceDNA < 10 || variant.nMatchingVariantDNA < 10)
                    {
                        //
                        // Not enough reads to consider this variant.
                        //
                        continue;
                    }
                    //
                    // See if this variant is close to any mutations.
                    //
                    var chromosome = ExpressionTools.chromosomeNameToNonChrForm(variant.contig);
                    int rangeToConsider;

                    if (experiment.TumorDNAAnalysis.anyPaired && experiment.TumorDNAAnalysis.medianPairedReadGap != -1)
                    {
                        rangeToConsider = (int)(experiment.TumorDNAAnalysis.meanGoodBases + experiment.TumorDNAAnalysis.medianPairedReadGap) * 3;
                    }
                    else
                    {
                        rangeToConsider = (int)experiment.TumorDNAAnalysis.meanGoodBases * 3;
                    }

                    if (scatterGraphEntriesByParticipant.ContainsKey(experiment.participant.participantId))
                    {
                        foreach (var mutation in scatterGraphEntriesByParticipant[experiment.participant.participantId])
                        {
                            if (mutation.Chromosome == chromosome && Math.Abs(mutation.Start_Position - variant.loc) < rangeToConsider)
                            {
                                lock (outputFile)
                                {
                                    outputFile.WriteLine(experiment.participant.participantId + "\t" + ExpressionTools.ConvertToExcelString(mutation.Hugo_Symbol) + "\t" + (mutation.IsSingle ? "Single" : "Multiple") + "\t" +
                                        chromosome + "\t" + mutation.Start_Position + "\t" + mutation.Tumor_Seq_Allele_1 + "\t" + mutation.Match_Norm_Seq_Allele1 + "\t" +
                                        mutation.n_DNA_Matching_Tumor + "\t" + mutation.n_DNA_Matching_Reference + "\t" + mutation.n_RNA_Matching_Tumor + "\t" + mutation.n_RNA_Matching_Reference + "\t" +
                                        variant.loc + "\t" + variant.alt + "\t" + variant.Ref + "\t" + variant.nMatchingVariantDNA + "\t" + variant.nMatchingReferenceDNA + "\t" + variant.nMatchingVariantRNA + "\t" + variant.nMatchingReferenceRNA + "\t" +
                                        experiment.TumorDNAAnalysis.bamFileName + "\t" + experiment.TumorRNAAnalysis.bamFileName + "\t" + experiment.NormalDNAAnalysis.bamFileName);
                                }
                            }
                        }
                    }
                }
            }

            lock (experiments)
            {
                nExperimentsProcessed += nUnreportedExperimentsProcessed;
            }
        }

        static void Main(string[] args)
        {
            List<string> excludedAnalyses;
            Dictionary<string, ExpressionTools.TCGARecord> tcgaRecords;
            Dictionary<string, string> sampleToParticipantIDMap;
            Dictionary<string, ExpressionTools.Sample> allSamples;
            Dictionary<string, ExpressionTools.Participant> participants;

            var timer = new Stopwatch();
            timer.Start();

            ExpressionTools.LoadStateFromExperimentsFile(out excludedAnalyses, out tcgaRecords, out sampleToParticipantIDMap, out participants, out experiments, out allSamples);

            var experimentsByRNAAnalysisId = new Dictionary<string, ExpressionTools.Experiment>();
            foreach (var experiment in experiments) {
                if (experiment.TumorRNAAnalysis != null) {
                    experimentsByRNAAnalysisId.Add(experiment.TumorRNAAnalysis.analysis_id, experiment);
                }
            }

            var allScatterGraphEntries = ExpressionTools.GeneScatterGraphLine.LoadAllGeneScatterGraphEntries(ExpressionTools.geneScatterGraphsDirectory, false, experimentsByRNAAnalysisId, "*");

            scatterGraphEntriesByParticipant = new Dictionary<string,List<ExpressionTools.GeneScatterGraphLine>>(); 
            foreach (var scatterGraphEntry in allScatterGraphEntries)
            {
                scatterGraphEntry.Chromosome = ExpressionTools.chromosomeNameToNonChrForm(scatterGraphEntry.Chromosome);

                if (!scatterGraphEntriesByParticipant.ContainsKey(scatterGraphEntry.participantId)) {
                    scatterGraphEntriesByParticipant.Add(scatterGraphEntry.participantId, new List<ExpressionTools.GeneScatterGraphLine>());
                }

                scatterGraphEntriesByParticipant[scatterGraphEntry.participantId].Add(scatterGraphEntry);
            }

            Console.WriteLine("Loaded state in " + (timer.ElapsedMilliseconds + 500) / 1000 + "s");

            timer.Restart();

            //
            // Now run through the experiments and find candidates.
            //

            outputFile = ExpressionTools.CreateStreamWriterWithRetry(@"f:\temp\expression\phasingCandidates.txt");

            outputFile.WriteLine("ParticipantID\tGene\tIsSingle\tChromosome\tMutation Start Location\tMutatant Allele\tMutant Wild Type\tMutant DNA Count\tWild Type DNA Count\tMutant RNA Count\tWild Type RNA Count" + 
                "\tVariant Location\tVariant Allele 1\tVariant Allele 2\tVariant Allele 1 DNA count\tVariant Allele 2 DNA count\t" +
                "Variant Allele 1 RNA Count\tVariant Allele 2 RNA Count\tTumor DNA BAM\tTumor RNA Bam\tNormal DNA Bam");


            var threads = new List<Thread>();
            for (int i = 0; i < Environment.ProcessorCount; i++)
            {
                threads.Add(new Thread(() => WorkerThread()));
            }

            threads.ForEach(t => t.Start());
            threads.ForEach(t => t.Join());

            outputFile.WriteLine("**done**");
            outputFile.Close();

            Console.WriteLine("Processed " + experiments.Count() + " experiments in " + (timer.ElapsedMilliseconds + 500) / 1000 + "s, of which " + nMissingAnnotatedSelectedVariants + " had problems with their annotated selected variants files.");
        }
    }
}
