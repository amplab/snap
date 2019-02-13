using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using ASELib;
using System.IO;
using System.Diagnostics;

namespace SingleReadPairPhasing
{
    class Program
    {

        static ASETools.CommonData commonData;

        class Phase
        {
            public enum WhichPhase { inPhase, outOfPhase, other};
            public Dictionary<WhichPhase, int> counts = new Dictionary<WhichPhase, int>();
            public Dictionary<WhichPhase, List<string>> casesByPhase = new Dictionary<WhichPhase, List<string>>();

            static List<WhichPhase> allPhases = new List<WhichPhase>();
            static Phase()
            {
                allPhases.Add(WhichPhase.inPhase);
                allPhases.Add(WhichPhase.outOfPhase);
                allPhases.Add(WhichPhase.other);
            }

            public Phase()
            {
                foreach (var phase in allPhases)
                {
                    counts.Add(phase, 0);
                    casesByPhase.Add(phase, new List<string>());
                }
            }

            public void mergeWith(Phase peer)
            {
                foreach (var phase in allPhases)
                {
                    counts[phase] += peer.counts[phase];
                    casesByPhase[phase].AddRange(peer.casesByPhase[phase]);
                }
            } // mergeWith

            public int totalCases()
            {
                return allPhases.Sum(_ => counts[_]);
            }

            public void addCase(string caseId, WhichPhase phase)
            {
                counts[phase]++;
                casesByPhase[phase].Add(caseId);
            }
        } // Phase

        class SingleGenePhasing
        {
            public Dictionary<int, Phase> phasesByRange = new Dictionary<int, Phase>();

            public void mergeWith(SingleGenePhasing peer)
            {
                foreach (var entry in peer.phasesByRange)
                {
                    var range = entry.Key;
                    var phase = entry.Value;

                    if (!phasesByRange.ContainsKey(range))
                    {
                        phasesByRange.Add(range, phase);
                    } else
                    {
                        phasesByRange[range].mergeWith(phase);
                    }
                }
            } // mergeWith

            public int count() // total cases phased 
            {
                return phasesByRange.Select(_ => _.Value.totalCases()).Sum();
            }

            public void counts(out int totalInPhase, out int totalOutOfPhase, out int totalOther)
            {
                totalInPhase = phasesByRange.Select(_ => _.Value.counts[Phase.WhichPhase.inPhase]).Sum();
                totalOutOfPhase = phasesByRange.Select(_ => _.Value.counts[Phase.WhichPhase.outOfPhase]).Sum();
                totalOther = phasesByRange.Select(_ => _.Value.counts[Phase.WhichPhase.other]).Sum();
            }
        } // SingleGenePhasing

        class AllGenesPhasing
        {
            public Dictionary<string, SingleGenePhasing> allGenes = new Dictionary<string, SingleGenePhasing>();

            public void mergeWith(AllGenesPhasing peer)
            {
                foreach (var entry in peer.allGenes)
                {
                    var hugo_symbol = entry.Key;
                    var singleGenePhasing = entry.Value;

                    if (!allGenes.ContainsKey(hugo_symbol))
                    {
                        allGenes.Add(hugo_symbol, singleGenePhasing);
                    } else
                    {
                        allGenes[hugo_symbol].mergeWith(singleGenePhasing);
                    }
                }
            } // mergeWith
        } // AllGenesPhasing

        static AllGenesPhasing globalResults = new AllGenesPhasing();

        static void Main(string[] args)
        {
            commonData = ASETools.CommonData.LoadCommonData(args);
            if (null == commonData)
            {
                return;
            }

            if (commonData.listOfCases.Any(_ => _.tumor_dna_reads_at_tentative_selected_variants_filename == "" || _.tumor_rna_reads_at_tentative_selected_variants_filename == "" || _.extracted_maf_lines_filename == ""))
            {
                Console.WriteLine("Not all cases have tumor DNA and RNA reads at selected variants files.");
                return;
            }

            var outputFile = ASETools.CreateStreamWriterWithRetry(commonData.configuration.finalResultsDirectory + ASETools.SingleReadPhasingFilename);
            if (null == outputFile)
            {
                Console.WriteLine("Unable to open output file.");
                return;
            }

            int nPerDot;

            ASETools.PrintMessageAndNumberBar("Processing", "cases", commonData.listOfCases.Count(), out nPerDot);

            var threading = new ASETools.WorkerThreadHelper<ASETools.Case, AllGenesPhasing>(commonData.listOfCases, HandleOneCase, FinishUp, null, nPerDot);
            threading.run();

            Console.WriteLine();

            var phaseHistogram = new ASETools.PreBucketedHistogram(0, 1.01, 0.01);

            foreach (var entry in globalResults.allGenes)
            {
                var hugo_symbol = entry.Key;
                var singleGenePhasing = entry.Value;
                int totalCases = singleGenePhasing.count();



                int inPhase, outOfPhase, other;
                singleGenePhasing.counts(out inPhase, out outOfPhase, out other);

                if (inPhase + outOfPhase < /*BJB 10*/2)
                {
                    //
                    // Not enough cases to phase.
                    //
                    continue;
                }

                phaseHistogram.addValue((double)inPhase / (inPhase + outOfPhase));

                outputFile.WriteLine(hugo_symbol + "\t" + inPhase + "(" + ASETools.ratioToPercentage(inPhase, totalCases) + "%) in phase, " + outOfPhase + "(" + ASETools.ratioToPercentage(outOfPhase, totalCases) + "%) out of phase, " +
                    other + "(" + ASETools.ratioToPercentage(other, totalCases) + "%) other.");
                int maxDistance = singleGenePhasing.phasesByRange.Select(_ => _.Key).Max();

                outputFile.WriteLine("Top Of Range\tin phase\tout of phase\tother\t% in phase\t% out of phase\t%other");
                const int rangeGroupSize = 50;
                for (int rangeGroup = rangeGroupSize; rangeGroup < ((maxDistance + rangeGroupSize - 1) / rangeGroupSize) * rangeGroupSize; rangeGroup += rangeGroupSize)
                {
                    var inRangeGroup = singleGenePhasing.phasesByRange.Where(_ => _.Key > rangeGroup - rangeGroupSize && _.Key <= rangeGroup).ToList();
                    if (inRangeGroup.Count() == 0)
                    {
                        continue;
                    }

                    var totalCasesInRangeGroup = inRangeGroup.Select(_ => _.Value.totalCases()).Sum();
                    outputFile.WriteLine(rangeGroup + "\t" + inRangeGroup.Select(_ => _.Value.counts[Phase.WhichPhase.inPhase]).Sum() + "\t" + inRangeGroup.Select(_ => _.Value.counts[Phase.WhichPhase.outOfPhase]).Sum() + "\t" + inRangeGroup.Select(_ => _.Value.counts[Phase.WhichPhase.other]).Sum() + "\t" +
                        ASETools.ratioToPercentage(inRangeGroup.Select(_ => _.Value.counts[Phase.WhichPhase.inPhase]).Sum(), totalCasesInRangeGroup) + "%\t" + ASETools.ratioToPercentage(inRangeGroup.Select(_ => _.Value.counts[Phase.WhichPhase.outOfPhase]).Sum(), totalCasesInRangeGroup) + "%\t" +
                        ASETools.ratioToPercentage(inRangeGroup.Select(_ => _.Value.counts[Phase.WhichPhase.other]).Sum(), totalCasesInRangeGroup) + "%");
                } // for each range group

                writeCases(outputFile, singleGenePhasing, Phase.WhichPhase.inPhase);
                writeCases(outputFile, singleGenePhasing, Phase.WhichPhase.outOfPhase);
                outputFile.WriteLine();
            } // for each gene

            outputFile.WriteLine("Histogram of genes by fraction in phase");
            outputFile.WriteLine(ASETools.HistogramResultLine.Header());
            phaseHistogram.ComputeHistogram().ToList().ForEach(_ => outputFile.WriteLine(_));

            outputFile.WriteLine("**done**");
            outputFile.Close();

            Console.WriteLine("Run time " + ASETools.ElapsedTimeInSeconds(commonData.timer));
        } // Main

        static void writeCases(StreamWriter outputFile, SingleGenePhasing singleGenePhasing, Phase.WhichPhase whichPhase)
        {
            var caseIds = new List<string>();

            singleGenePhasing.phasesByRange.ToList().ForEach(_ => caseIds.AddRange(_.Value.casesByPhase[whichPhase]));

            if (caseIds.Count() == 0)
            {
                return;
            }

            outputFile.Write(whichPhase + ":");
            caseIds.ForEach(_ => outputFile.Write("\t" + _));
            outputFile.WriteLine();
        }

        static void HandleOneCase(ASETools.Case case_, AllGenesPhasing state)
        {
            var tumorDNAExtractedReadsFile = new ASETools.ConsolodatedFileReader();
            if (!tumorDNAExtractedReadsFile.open(case_.tumor_dna_reads_at_tentative_selected_variants_filename))
            {
                throw new Exception("Unable to open tumor DNA extracted reads file " + case_.tumor_dna_reads_at_tentative_selected_variants_filename);
            }

            var tumorRNAExtractedReadsFile = new ASETools.ConsolodatedFileReader();
            if (!tumorRNAExtractedReadsFile.open(case_.tumor_rna_reads_at_tentative_selected_variants_filename))
            {
                throw new Exception("Unable to open tumor RNA extracted reads file " + case_.tumor_rna_reads_at_tentative_selected_variants_filename);
            }


            var readers = new Dictionary<string, ASETools.ConsolodatedFileReader>();
            readers.Add(case_.tumor_dna_file_id, tumorDNAExtractedReadsFile);
            readers.Add(case_.tumor_rna_file_id, tumorRNAExtractedReadsFile);

            var mafLines = ASETools.MAFLine.ReadFileOnlyOnePerLocus(case_.extracted_maf_lines_filename, case_.case_id, false);

            var mutationsByGene = mafLines.Where(_ => !_.IsSilent() && _.t_ref_count + _.t_alt_count >= 10 && (double)_.t_ref_count / (_.t_ref_count + _.t_alt_count) >= 0.4 && (double)_.t_ref_count / (_.t_ref_count + _.t_alt_count) <= 0.6).GroupByToDict(_ => _.Hugo_Symbol);
            foreach (var mutatedGeneEntry in mutationsByGene)
            {
                var hugo_symbol = mutatedGeneEntry.Key;
                var mutations = mutatedGeneEntry.Value;
                int nMutations = mutations.Count();

                if (nMutations != 2 || mutations.Any(_ => _.Variant_Type != "SNP")) // Could theoretically handle > 2 and indels, but this is just a cheezy first cut to see if we need to run a real phasing tool.
                {
                    continue;
                }

                if (Math.Abs(mutations[0].Start_Position - mutations[1].Start_Position) < 5) // Ignore really close ones because they're probably the same mutation called as two in the MAF.  Assumes nMutations == 2
                {
                    //
                    // Two mutations at the same locus are one mutation.
                    //
                    continue;
                }

                List<ASETools.SAMLine> allReads = new List<ASETools.SAMLine>();

                foreach (var mutation in mutations)
                {
                    foreach (var readsFileEntry in readers) {
                        var fileId = readsFileEntry.Key;
                        var readsFile = readsFileEntry.Value;

                        var subfile = readsFile.getSubfile(fileId  + mutation.getExtractedReadsExtension());
                        if (subfile == null)
                        {
                            throw new Exception("Unable to open subfile for case " + case_.case_id + ", subfile name " + fileId + mutation.getExtractedReadsExtension());
                        }

                        allReads.AddRange(ASETools.SAMLine.ReadFromFile(subfile).Where(_ => !_.isUnmapped() && _.mapq >= 10 && !_.isSecondaryAlignment()));

                        subfile.Close();
                    } // reader
                } // mutation

                var phaseCounts = new Phase(); // With this one, we're counting read sets, not mutations.  It will be used to generate a single mutation entry for this gene eventually.
                int nInPhaseRef = 0;
                int nInPhaseAlt = 0;

                var readSets = allReads.GroupByToDict(_ => _.qname);    // This will group together both paired-end reads (which have the same qnames) and duplicates, but the duplicates aren't an issue.

                foreach (var readSet in readSets.Select(_ => _.Value))
                {
                    //
                    //
                    int[] refCount = new int[mutations.Count()];
                    int[] altCount = new int[mutations.Count()];
                    int[] otherCount = new int[mutations.Count()];

                    for (int whichMutation = 0; whichMutation < nMutations; whichMutation++)
                    {
                        var mutation = mutations[whichMutation];
                        foreach (var read in readSet)
                        {
                            if (ASETools.chromosomeNameToNonChrForm(read.rname) != ASETools.chromosomeNameToNonChrForm(mutation.Chromosome))
                            {
                                continue;
                            }

                            if (read.mappedBases.ContainsKey(mutation.Start_Position))
                            {
                                var readBase = read.mappedBases[mutation.Start_Position];
                                if (readBase == mutation.Reference_Allele[0])
                                {
                                    refCount[whichMutation]++;
                                } else if (readBase == mutation.Tumor_Seq_Allele2[0]) // For whatever reason, Tumor_Seq_Allele2 is the alt for heterozygous loci, not Tumor_Seq_Allele1
                                {
                                    altCount[whichMutation]++;
                                } else
                                {
                                    otherCount[whichMutation]++;
                                }
                            } // if this read maps this base
                        } // for each read in the read set (paired ends)
                    } // for each mutation of this gene in this patient

                    int minReadCount = Enumerable.Range(0, nMutations).Select(_ => refCount[_] + altCount[_] + otherCount[_]).Min();

                    if (minReadCount == 0)
                    {
                        // Not all mutations covered by this read, skip it.
                        continue;
                    }

                    if (Enumerable.Range(0, nMutations).Select(_ => (double)Math.Max(refCount[_], altCount[_]) / (refCount[_] + altCount[_] + otherCount[_])).Min() < 0.9)
                    {
                        //
                        // At least one of the mutations is floating (not at least 90% one allele), so this is an other.  Since this is essentially a pair of reads (with maybe some duplicates), this really means that either a read
                        // shows something other than ref or alt, or both reads of a paired end cover a mutation and have opposite values (which indicates a sequencing error, since it should be the same R/DNA).
                        //
                        phaseCounts.addCase(case_.case_id, Phase.WhichPhase.other);
                    }
                    else
                    {
                        bool[] mutationIsRef = new bool[nMutations];
                        for (int i = 0; i < nMutations; i++)
                        {
                            mutationIsRef[i] = refCount[i] > altCount[i];
                        }

                        //
                        // Here we assume nMutations == 2 (which we forced way above, but didn't afterward assume anywhere else)
                        //
                        if (mutationIsRef[0] == mutationIsRef[1])
                        {
                            phaseCounts.addCase(case_.case_id, Phase.WhichPhase.inPhase);
                            if (mutationIsRef[0])
                            {
                                nInPhaseRef++;
                            } else
                            {
                                nInPhaseAlt++;
                            }
                        } else
                        {
                            phaseCounts.addCase(case_.case_id, Phase.WhichPhase.outOfPhase);
                        }
                    }
                } // for each set of reads (essentially, two paired-end reads)

                if (phaseCounts.totalCases() < 5)
                {
                    //
                    // Too little read depth.
                    //
                    continue;
                }

                if (!state.allGenes.ContainsKey(hugo_symbol))
                {
                    state.allGenes.Add(hugo_symbol, new SingleGenePhasing());
                }

                // Again assuming nMutations == 2
                int mutationsRange = Math.Abs(mutations[0].Start_Position - mutations[1].Start_Position); // Distance between the two mutations
                if (!state.allGenes[hugo_symbol].phasesByRange.ContainsKey(mutationsRange))
                {
                    state.allGenes[hugo_symbol].phasesByRange.Add(mutationsRange, new Phase());
                }

                if ((double)phaseCounts.counts[Phase.WhichPhase.inPhase] >= phaseCounts.totalCases() * 0.9)
                {
                    double altFraction = (double)nInPhaseAlt / (nInPhaseAlt + nInPhaseRef);
                    if (altFraction < 0.4 || altFraction > 0.6)
                    {
                        //
                        // Too many one or the other is suspicious.  Bounce it.
                        //
                        state.allGenes[hugo_symbol].phasesByRange[mutationsRange].addCase(case_.case_id, Phase.WhichPhase.other);
                    }
                    else
                    {
                        state.allGenes[hugo_symbol].phasesByRange[mutationsRange].addCase(case_.case_id, Phase.WhichPhase.inPhase);
                    }
                } else if ((double)phaseCounts.counts[Phase.WhichPhase.outOfPhase] >= phaseCounts.totalCases() * 0.9)
                {
                    state.allGenes[hugo_symbol].phasesByRange[mutationsRange].addCase(case_.case_id, Phase.WhichPhase.outOfPhase);
                } else
                {
                    state.allGenes[hugo_symbol].phasesByRange[mutationsRange].addCase(case_.case_id, Phase.WhichPhase.other);
                }

            } // gene
        } // HandleOneCase

        static void FinishUp(AllGenesPhasing state)
        {
            lock (globalResults)
            {
                globalResults.mergeWith(state);
            }
        } // FinishUp
    }
}
