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

        class ReadCoverageState
        {
            public ReadCoverageState(ExpressionTools.ExonicMap exonicMap_)
            {
                exonicMap = exonicMap_;
            }

            public Dictionary<string, Dictionary<int, int>> coverageCounts = new Dictionary<string, Dictionary<int, int>>();   // chromosome->locus->count.  Dictionary<> since exons are only 1%
            ExpressionTools.ExonicMap exonicMap;
            public long totalCoverageOnAllExonicLoci = 0;

            public void processBase(string strippedContigName, int location, int currentMappedReadCount)
            {
                if (!exonicMap.isLocationInAnExon(strippedContigName, location))
                {
                    //
                    // Not in an exon, so not interesting.
                    //
                    return;
                }

                if (!coverageCounts.ContainsKey(strippedContigName))
                {
                    coverageCounts.Add(strippedContigName, new Dictionary<int, int>());
                }

                coverageCounts[strippedContigName].Add(location, currentMappedReadCount);
                totalCoverageOnAllExonicLoci += currentMappedReadCount;
            }
        }

        class IsoformWithCounts
        {
            public IsoformWithCounts(ExpressionTools.Isoform isoform_)
            {
                isoform = isoform_;
            }

            public ExpressionTools.Isoform isoform;
            public int nBases = 0;
            public long totalBaseDepth = 0;
        }

        class GeneWithIsoformsAndCounts
        {
            public GeneWithIsoformsAndCounts(ExpressionTools.GeneLocationInfo geneLocationInfo_)
            {
                geneLocationInfo = geneLocationInfo_;
            }

            public ExpressionTools.GeneLocationInfo geneLocationInfo;
            public List<IsoformWithCounts> isoforms = new List<IsoformWithCounts>();

            public int totalBasesCovered = 0;
            public int totalBaseCoverage = 0;
        }

        static void ProcessParticipants()
        {
            ExpressionTools.AllcountReader reader = null;
            while (true)
            {
                if (null != reader)
                {
                    reader.Close();
                }

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

                string referenceClass = (experiment.TumorDNAAnalysis.refassemShortName == "hg18" || experiment.TumorDNAAnalysis.refassemShortName == "hg18_broad_variant") ?
                    "hg18" : "hg19";

                reader = new ExpressionTools.AllcountReader(experiment.TumorDNAAnalysis.allcountFileName);

                long mappedHQNuclearReads;
                int numContigs;
                if (!reader.openFile(out mappedHQNuclearReads, out numContigs))
                {
                    Console.WriteLine("Error in allcount file for participant " + participantId + ", it's probably corrupt.");
                    continue;
                }

                var readCoverageState = new ReadCoverageState(exonicMaps[referenceClass]);

                ExpressionTools.AllcountReader.ProcessBase processBase = (x, y, z) => readCoverageState.processBase(x, y, z);

                if (!reader.ReadAllcountFile(processBase))
                {
                    Console.WriteLine("Error reading allcount file for participant id " + participantId + ", it's probably corrupt.");
                    continue;
                }

                //
                // We now know the total coverage depth of all exonic loci.  For each gene go through its exons and report its average depth.
                //
                var genes = new List<GeneWithIsoformsAndCounts>();
                var wholeExomeSeenBases = new Dictionary<string, Dictionary<int, bool>>();
                int wholeExomeNBases = 0;
                long wholeExomeTotalCoverage = 0;

                foreach (var geneEntry in geneLocationInformation[referenceClass].genesByName)
                {
                    var gene = geneEntry.Value;
                    var seenBases = new Dictionary<int, bool>();    // Keeps track of bases we've already seen so that we don't multiply count regions in more than one exon

                    var geneWithIsoformsAndCounts = new GeneWithIsoformsAndCounts(gene);
                    genes.Add(geneWithIsoformsAndCounts);
                    foreach (var isoform in gene.isoforms)
                    {
                        var isoformWithCounts = new IsoformWithCounts(isoform);
                        geneWithIsoformsAndCounts.isoforms.Add(isoformWithCounts);

                        foreach (var exon in isoform.exons)
                        {
                            isoformWithCounts.nBases += exon.end - exon.start + 1;  // +1 because the range is inclusive

                            for (int locus = exon.start; locus <= exon.end; locus++)
                            {

                                int coverageAtThisLocus;
                                if (readCoverageState.coverageCounts.ContainsKey(isoform.chromosome) &&
                                        readCoverageState.coverageCounts[isoform.chromosome].ContainsKey(locus))
                                {
                                    coverageAtThisLocus = readCoverageState.coverageCounts[isoform.chromosome][locus];
                                } else {
                                    coverageAtThisLocus = 0;
                                }

                                isoformWithCounts.totalBaseDepth += coverageAtThisLocus;

                                if (!seenBases.ContainsKey(locus))  // If we haven't seen this base in some other isoform of this gene
                                {
                                    seenBases.Add(locus, true);

                                    geneWithIsoformsAndCounts.totalBasesCovered++;
                                    geneWithIsoformsAndCounts.totalBaseCoverage += coverageAtThisLocus;
 
                                    if (!wholeExomeSeenBases.ContainsKey(gene.chromosome)) {
                                        wholeExomeSeenBases.Add(gene.chromosome, new Dictionary<int,bool>());
                                    }

                                    if (!wholeExomeSeenBases[gene.chromosome].ContainsKey(locus)) { // If we haven't seen this locus in some other gene or isoform
                                        wholeExomeSeenBases[gene.chromosome].Add(locus, true);

                                        wholeExomeNBases++;
                                        wholeExomeTotalCoverage += coverageAtThisLocus;
                                    }
                                } // if we've seen this base in this gene
                            }
                        }
                    } // foreach isoform
                } // foreach gene

                genes.Sort((x, y) => x.geneLocationInfo.hugoSymbol.CompareTo(y.geneLocationInfo.hugoSymbol));   // Sort them alphabetically

                var outputFile = ExpressionTools.CreateStreamWriterWithRetry(
                    ExpressionTools.GetDirectoryPathFromFullyQualifiedFilename(experiment.TumorDNAAnalysis.allcountFileName) + '\\' +
                    ExpressionTools.GetAnalysisIdFromPathname(experiment.TumorDNAAnalysis.allcountFileName) + ".gene_coverage.txt");
                outputFile.WriteLine("MeasureGeneCoverage v1.0\tnBasesInExome: " + wholeExomeNBases + "\tTotalCoverageInExome: " + wholeExomeTotalCoverage);
                outputFile.WriteLine("Hugo Symbol\tTotal Bases Covered\tTotal Base Coverage\tnIsoforms\tucscId\tnBasesInIsoform\ttotalCoverageInIsoform");
                     
                for (int i = 0; i < genes.Count(); i++)
                {
                    outputFile.Write(genes[i].geneLocationInfo.hugoSymbol + "\t" + genes[i].totalBasesCovered + "\t" + genes[i].totalBaseCoverage + "\t" +
                        genes[i].isoforms.Count());

                    foreach (var isoform in genes[i].isoforms)
                    {
                        outputFile.Write("\t" + isoform.isoform.ucscId + "\t" + isoform.nBases + "\t" + isoform.totalBaseDepth);
                    } // foreach isoform
                    outputFile.WriteLine();
                }
                outputFile.WriteLine("**done**");
                outputFile.Close();
            } // while (true) (the outer loop over pariticpants)
        } // ProcessParticipants

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
