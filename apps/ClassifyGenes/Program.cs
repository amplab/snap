using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using System.IO;
using ExpressionLib;
using System.Diagnostics;

namespace ClassifyGenes
{
    class Program
    {

        class PerGeneState
        {
            public PerGeneState(string hugoSymbol_)
            {
                hugoSymbol = hugoSymbol_;
            }


            public string hugoSymbol;
            //
            // The classification counts.  These are mutatally exclusive catagories, with a tumor being slotted into the first one that fits.
            //
            public int nSex = 0;                                // Sex chromosome genes
            public int nDNACoverageOff = 0;                     // Tumors with DNA coverage not close enough to overall average
            public int nZUnknown = 0;                           // z score unknown for this mutation
            public int nLowTumorDNA = 0;                        // Tumors with too little tumor DNA (minor subclones)
            public int nHighTumorDNA = 0;                       // Tumors with too little normal DNA (loss of heterozygozity or copy number variations)
            public int nDoubleMutant = 0;                       // Tumors with two (or more) mutations in this gene.
            public int nOneMutationPlusLackOfExpression = 0;    // Tumors that lost one through mutation and one through epigenetic loss (as determined at the mutation site).
            public int nOneMutationWithExpression = 0;          // Tumors with exactly one muation, but with reasonable expression of the wild type allele
            public int nZeroMutationsWithLowExpression = 0;     // Tumors with no mutations, but no more than 10% of mean overall expression across all isoforms
            public int nZeroMutationsWithSomeExpression = 0;    // Everything else

            public int overallZeroMutations = 0;                // No mutations regardless of coverage
            public int overallOneMutation = 0;                  // One mutation, regardless of coverage
            // No special field for multiple mutations, since nDoubleMutant covers them
        }

        static string ComputeFraction(int numerator, int denominator)
        {
            if (0 == denominator) return "-";

            return Convert.ToString((double)numerator / (double)denominator);
        }

        static void Main(string[] args)
        {
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

            var allScatterGraphLines = ExpressionTools.GeneScatterGraphLine.LoadAllGeneScatterGraphEntries(ExpressionTools.geneScatterGraphsDirectory, false, experimentsByRNAAnalysisID, "*");

            var genesToConsider = new Dictionary<string, PerGeneState>();

            foreach (var scatterGraphLine in allScatterGraphLines)
            {
                if (!genesToConsider.ContainsKey(scatterGraphLine.Hugo_Symbol)) // We only have scatter graph lines for the genes that interest us (the ones with enough mutations in the data set)
                {
                    genesToConsider.Add(scatterGraphLine.Hugo_Symbol, new PerGeneState(scatterGraphLine.Hugo_Symbol));
                }
            }

            var scatterGraphLinesByParticipant = new Dictionary<string, List<ExpressionTools.GeneScatterGraphLine>>();

            foreach (var experiment in experiments)
            {
                scatterGraphLinesByParticipant.Add(experiment.participant.participantId, new List<ExpressionTools.GeneScatterGraphLine>());
            }

            foreach (var line in allScatterGraphLines) 
            {
                scatterGraphLinesByParticipant[line.participantId].Add(line);
            }

            Console.WriteLine("Total of " + allScatterGraphLines.Count() + " gene scatter graph lines loaded.");

            //
            // Now run through all of the experiments and classify all of their genes that are in genesToConsider.
            //
            int nExperimentsProcessed = 0;
            int nExperimentsWithMissingPrecursors = 0;

            var timer = new Stopwatch();
            timer.Start();

            var lastPrint = timer.ElapsedMilliseconds;
            var nExperimentsSkippedForLowDNACoverage = 0;


            foreach (var experiment in experiments)
            {
                nExperimentsProcessed++;

                if (timer.ElapsedMilliseconds - lastPrint > 60000)
                {
                    Console.WriteLine("Processed " + nExperimentsProcessed + " experiments in " + (timer.ElapsedMilliseconds + 500) / 1000 + "s");
                    lastPrint = timer.ElapsedMilliseconds;
                }

                if (experiment.TumorDNAAnalysis == null || experiment.TumorDNAAnalysis.geneCoverageFileName == null || experiment.TumorDNAAnalysis.geneCoverageFileName == "")
                {
                    nExperimentsWithMissingPrecursors++;
                    continue;
                }

                var geneCoverageFile = ExpressionTools.CreateStreamReaderWithRetry(experiment.TumorDNAAnalysis.geneCoverageFileName);
                var headerLine = geneCoverageFile.ReadLine();
                var fields = headerLine.Split('\t');

                const string exomeSizeString = "nBasesInExome: ";
                const string exomeCoverageString = "TotalCoverageInExome: ";

                if (fields.Count() != 3 || fields[0] != "MeasureGeneCoverage v1.0" || 
                    fields[1].Count() <= exomeSizeString.Count() || fields[1].Substring(0, exomeSizeString.Count()) != exomeSizeString ||
                    fields[2].Count() <= exomeCoverageString.Count() || fields[2].Substring(0, exomeCoverageString.Count()) != exomeCoverageString)
                {
                    Console.WriteLine("Corrupt or wrong version gene coverage file " + experiment.TumorDNAAnalysis.geneCoverageFileName + ", header line: " + headerLine + ", skipping");
                    continue;
                }

                try
                {
                    double nBasesInExome = 0;
                    double totalExomeCoverage = 0;

                    nBasesInExome = Convert.ToInt64(fields[1].Substring(exomeSizeString.Count()));
                    totalExomeCoverage = Convert.ToInt64(fields[2].Substring(exomeCoverageString.Count()));

                    double overallCoverage = totalExomeCoverage / nBasesInExome;

                    if (overallCoverage < 5)
                    {
                        //Console.WriteLine("Participant " + experiment.participant.participantId + " doesn't have enough Tumor DNA coverage to proceed. " + overallCoverage + " < 5.  Skipping.");
                        nExperimentsSkippedForLowDNACoverage++;
                        continue;
                    }

                    geneCoverageFile.ReadLine();    // Skip the column header
                    bool seenDone = false;

                    string line;
                    while (null != (line = geneCoverageFile.ReadLine()))
                    {
                        if (seenDone)
                        {
                            Console.WriteLine("File " + experiment.TumorDNAAnalysis.geneCoverageFileName + " continues after **done**.");
                            break;
                        }

                        if ("**done**" == line)
                        {
                            seenDone = true;
                            continue;
                        }

                        fields = line.Split('\t');

                        if (fields.Count() < 3)
                        {
                            Console.WriteLine("Not enough fields in input line " + line + " in file " + experiment.TumorDNAAnalysis.geneCoverageFileName);
                            break;
                        }

                        string hugo_symbol = ExpressionTools.ConvertToNonExcelString(fields[0]).ToUpper();
                        if (!genesToConsider.ContainsKey(hugo_symbol))
                        {
                            //
                            // This isn't a gene we care about, skip it.
                            //
                            continue;
                        }

                        double nBasesInGeneExome = Convert.ToInt64(fields[1]);
                        double coverageInGeneExome = Convert.ToInt64(fields[2]);

                        double coverage = coverageInGeneExome / nBasesInGeneExome;

                        var mutations =  scatterGraphLinesByParticipant[experiment.participant.participantId].Where(x => x.Hugo_Symbol.ToUpper() == hugo_symbol && x.participantId == experiment.participant.participantId && x.Variant_Classification.ToLower() != "silent").ToList();

                        if (mutations.Count() == 0)
                        {
                            genesToConsider[hugo_symbol].overallZeroMutations++;
                        }
                        else if (mutations.Count() == 1)
                        {
                            genesToConsider[hugo_symbol].overallOneMutation++;

                        }
                        else
                        {
                            genesToConsider[hugo_symbol].nDoubleMutant++;
                            continue;
                        }

                        if (mutations.Count() == 1)
                        {
                            var mutation = mutations[0];

                            int nDNATotal = mutation.n_DNA_Matching_Neither + mutation.n_DNA_Matching_Reference + mutation.n_DNA_Matching_Tumor;
                            int nRNATotal = mutation.n_RNA_Matching_Neither + mutation.n_RNA_Matching_Reference + mutation.n_RNA_Matching_Tumor;

                            if (nDNATotal < 10)
                            {
                                genesToConsider[hugo_symbol].nDNACoverageOff++;
                            } 
                            else if ((double)mutation.n_DNA_Matching_Tumor / nDNATotal >= .6) {
                                genesToConsider[hugo_symbol].nHighTumorDNA++;
                            } 
                            else if ((double)mutation.n_DNA_Matching_Tumor / nDNATotal < .1) {
                                genesToConsider[hugo_symbol].nLowTumorDNA++;
                            } 
                            else if (!mutation.zKnown) {
                                genesToConsider[hugo_symbol].nZUnknown++;
                            } 
                            else
                            {
                                if (mutation.percentMeanNormal <= 0.25)
                                {
                                    genesToConsider[hugo_symbol].nOneMutationPlusLackOfExpression++;
                                }
                                else
                                {
                                    genesToConsider[hugo_symbol].nOneMutationWithExpression++;
                                }
                            }
                        }
                        else
                        {
                            if (coverage * 1.5 < overallCoverage || overallCoverage * 1.5 < coverage)
                            {
                                genesToConsider[hugo_symbol].nDNACoverageOff++;
                            }
                            else
                            {
                                // Figure out expression over whole gene (get from "lines" files).  For now, just dump them in one bucket.
                                genesToConsider[hugo_symbol].nZeroMutationsWithSomeExpression++;
                            }

                        }
                    }
                }
                catch (FormatException)
                {
                    Console.WriteLine("Format exception processing input file " + experiment.TumorDNAAnalysis.geneCoverageFileName + ", skipping.");
                    continue;
                }

                geneCoverageFile.Close();

            }

            //
            // And write the output.
            //
            var outputFile = ExpressionTools.CreateStreamWriterWithRetry(@"f:\temp\expression\gene_classification.txt");
            outputFile.WriteLine("Hugo Symbol\tDNA Coverage low\tz Unknown\tlow tumor DNA\thigh tumor DNA\t% high tumor DNA\t>1 mutation\t% >1 mutation\t1 mutation, one loss\t% 1 mutation, one loss\tsingle mutation\t% single mutation\tlow overall expression\t% low overall expression\tnone of the above\t% none of the above\tZero mutations\tOne mutation\tMultiple mutations");

            foreach (var geneEntry in genesToConsider)
            {
                var hugo_symbol = geneEntry.Key;
                var gene = geneEntry.Value;

                int totalNonAbnormal = gene.nHighTumorDNA + gene.nDoubleMutant + gene.nOneMutationPlusLackOfExpression + gene.nOneMutationWithExpression + gene.nZeroMutationsWithLowExpression + gene.nZeroMutationsWithSomeExpression;

                outputFile.WriteLine(
                    ExpressionTools.ConvertToExcelString(hugo_symbol) +
                    "\t" + (gene.nDNACoverageOff + nExperimentsSkippedForLowDNACoverage) +
                    "\t" + gene.nZUnknown +
                    "\t" + gene.nLowTumorDNA +
                    "\t" + gene.nHighTumorDNA + "\t" + ComputeFraction(gene.nHighTumorDNA, totalNonAbnormal) +
                    "\t" + gene.nDoubleMutant + "\t" + ComputeFraction(gene.nDoubleMutant, totalNonAbnormal) +
                    "\t" + gene.nOneMutationPlusLackOfExpression + "\t" + ComputeFraction(gene.nOneMutationPlusLackOfExpression, totalNonAbnormal) +
                    "\t" + gene.nOneMutationWithExpression + "\t" + ComputeFraction(gene.nOneMutationWithExpression, totalNonAbnormal) +
                    "\t" + gene.nZeroMutationsWithLowExpression + "\t" + ComputeFraction(gene.nZeroMutationsWithLowExpression, totalNonAbnormal) +
                    "\t" + gene.nZeroMutationsWithSomeExpression + "\t" + ComputeFraction(gene.nZeroMutationsWithSomeExpression, totalNonAbnormal) +
                    "\t" + gene.overallZeroMutations +
                    "\t" + gene.overallOneMutation +
                    "\t" + gene.nDoubleMutant +
                    "");
            }

            outputFile.Close();

            Console.WriteLine("Processed " + experiments.Count() + " experiments in " + (timer.ElapsedMilliseconds + 500) / 1000 + "s, of which " + nExperimentsSkippedForLowDNACoverage + " were skipped due to low DNA coverage.");
        }
    }
}
