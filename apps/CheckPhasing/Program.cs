using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using System.IO;
using System.Diagnostics;
using ASELib;

namespace CheckPhasing
{
    class Program
    {
        static ASETools.CommonData commonData;
        static ASETools.ASERepetitiveRegionMap repetitiveRegionMap;

        static int minReadsToDecideConsistency = 3;  // This maybe should be a configuration parameter

        static void Main(string[] args)
        {
            commonData = ASETools.CommonData.LoadCommonData(args);

            if (null == commonData)
            {
                return; // LoadCommonData printed an error message
            }

            if (commonData.configuration.commandLineArgs.Any(caseId => !commonData.listOfCases.Select(_ => _.case_id).Contains(caseId)))
            {
                Console.WriteLine("At least one of the command line args wasn't a case ID.");
                return;
            }

            repetitiveRegionMap = ASETools.ASERepetitiveRegionMap.loadFromFile(commonData.configuration.redundantChromosomeRegionFilename);
            if (null == repetitiveRegionMap)
            {
                Console.WriteLine("Unable to load repetitive region map from " + commonData.configuration.redundantChromosomeRegionFilename);
                return;
            }

            var casesToProcess = commonData.listOfCases.Where(_ => commonData.configuration.commandLineArgs.Contains(_.case_id)).ToList();

            int nPerDot;
            ASETools.PrintMessageAndNumberBar("Processing", "cases", casesToProcess.Count(), out nPerDot);

            var threading = new ASETools.WorkerThreadHelper<ASETools.Case, WorkerThreadState>(casesToProcess, HandleOneCase, FinishUp, null, nPerDot);

            threading.run();

            Console.WriteLine(ASETools.ElapsedTimeInSeconds(commonData.timer));

        } // Main

        class AlleleBalance
        {
            public int[,] counts = new int[3, 3];   // [ref = 0 alt = 1, other = 2 for first variant, ref = 0, alt = 1, other = 2 for second variant] value is count for this
            public const int ref_allele = 0;
            public const int alt_allele = 1;
            public const int other_allele = 2;
        } // AlleleBalance

        class DnaRnaAlleleBalance
        {
            public Dictionary<bool, AlleleBalance>  balanceByNucleicAcidType = new Dictionary<bool, AlleleBalance>();
            public const bool DNA = true;
            public const bool RNA = true;
        } // DnaRnaAlleleBalance

        class WorkerThreadState
        {
            public Dictionary<string, List<DnaRnaAlleleBalance>> alleleBalanceByGene = new Dictionary<string, List<DnaRnaAlleleBalance>>();
        }

        static Dictionary<string, List<DnaRnaAlleleBalance>> alleleBalanceByGene = new Dictionary<string, List<DnaRnaAlleleBalance>>();

        static void FinishUp(WorkerThreadState state)
        {
            lock (alleleBalanceByGene)
            {
                foreach (var hugo_symbol in state.alleleBalanceByGene.Select(_ => _.Key))
                {
                    if (!alleleBalanceByGene.ContainsKey(hugo_symbol))
                    {
                        alleleBalanceByGene.Add(hugo_symbol, state.alleleBalanceByGene[hugo_symbol]);
                    } else
                    {
                        alleleBalanceByGene[hugo_symbol].AddRange(state.alleleBalanceByGene[hugo_symbol]);
                    }
                } // hugo symbol
            } // lock
        } // FinishUp


        static void HandleOneCase(ASETools.Case case_, WorkerThreadState state)
        {
            var tentativeAnnotatedVariants = ASETools.AnnotatedVariant.readFile(case_.tentative_annotated_selected_variants_filename).Where(_ => _.variantType == "SNP").ToList();  // We use the tentative ones because the selected ones discard germline variants that are too close to one another
            if (null == tentativeAnnotatedVariants)
            {
                Console.WriteLine("Unable to load tentative annotated variants for case " + case_.case_id + ".  Skipping.");
                return;
            }

            var tumorDNAExtractedReadsFile = new ASETools.ConsolodatedFileReader();
            if (!tumorDNAExtractedReadsFile.open(case_.tumor_dna_reads_at_tentative_selected_variants_filename))
            {
                Console.WriteLine("Unable to open tumor DNA extracted reads file " + case_.tumor_dna_reads_at_tentative_selected_variants_filename + " for case " + case_.case_id + ".  Skipping.");
                return;
            }

            var tumorRNAExtractedReadsFile = new ASETools.ConsolodatedFileReader();
            if (!tumorRNAExtractedReadsFile.open(case_.tumor_rna_reads_at_tentative_selected_variants_filename))
            {
                Console.WriteLine("Unable to open tumor RNA extracted reads file " + case_.tumor_rna_reads_at_tentative_selected_variants_filename + " for case " + case_.case_id + ".  Skipping.");
                return;
            }

            var readers = new Dictionary<bool, ASETools.ConsolodatedFileReader>();
            readers.Add(true, tumorDNAExtractedReadsFile);
            readers.Add(false, tumorRNAExtractedReadsFile);

            var copyNumber = ASETools.CopyNumberVariation.ReadBothFiles(case_);
            if (null == copyNumber)
            {
                Console.WriteLine("Unable to read copy number for case " + case_.case_id + ".  Skipping.");
                return;
            }

            var readsForGivenVariant = new Dictionary<string, Dictionary<bool, List<ASETools.SAMLine>>>();    // Save reads we've already parsed so we don't have to do it more than once.



            var outputLines = new List<string>();


            foreach (var variant in tentativeAnnotatedVariants)
            {
                if (!variant.IsASECandidate(true, copyNumber, commonData) || !ASETools.isChromosomeAutosomal(variant.contig) || 
                    repetitiveRegionMap.isCloseToRepetitiveRegion(ASETools.ChromosomeNameToIndex(variant.contig), variant.locus, commonData.configuration.minDistanceFromRepetitiveRegion))
                {
                    //
                    // This isn't exactly the same test we use to check for selecting variants. There, we first do a test of the confidence interval of the VAFs for variants here or close to here, and only
                    // use this map if that is inconclusive.  But this is close enough for us.
                    //
                    continue;
                }

                //
                // See if there are other variants that would pair with this one, meaning that they're within distanceFromVariantToKeepReads.  We only select ones that have a higher locus than this one
                // to avoid selecting the same pairs twice, once for each end.
                //
                var matchingVariants = tentativeAnnotatedVariants.Where(_ => _.contig == variant.contig && _.locus > variant.locus && _.locus <= variant.locus + 1500).ToList();

                foreach (var matchingVariant in matchingVariants)
                {
                    var variants = new ASETools.AnnotatedVariant[2];
                    variants[0] = variant;
                    variants[1] = matchingVariant;

                    var reads = new Dictionary<bool, List<ASETools.SAMLine>>();    // dna -> reads
                    ASETools.BothBools.ToList().ForEach(_ => reads.Add(_, new List<ASETools.SAMLine>()));
                    bool loadFailed = false;

                    for (int i = 0; i < 2; i++)
                    {
                        var variantToLoad = variants[i];
                        var key = variantToLoad.contig + ":" + variantToLoad.locus;
                        if (!readsForGivenVariant.ContainsKey(variantToLoad.contig + ":" + variantToLoad.locus))
                        {
                            readsForGivenVariant.Add(key, new Dictionary<bool, List<ASETools.SAMLine>>());
                            foreach (var dna in ASETools.BothBools)
                            {
                                var subfileName = (dna ? case_.tumor_dna_file_id : case_.tumor_rna_file_id) + variantToLoad.getExtractedReadsExtension();
                                if (readers[dna].isSubfileTooBigToRead(subfileName))
                                {
                                    readsForGivenVariant.Remove(key);
                                    loadFailed = true;
                                    break;
                                }

                                var subfileReader = readers[dna].getSubfile(subfileName);

                                if (subfileReader == null)
                                {
                                    Console.WriteLine("Unable to load subfile reader for " + (dna ? case_.tumor_dna_file_id : case_.tumor_rna_file_id) + variantToLoad.getExtractedReadsExtension());
                                    return;
                                }
                                readsForGivenVariant[key].Add(dna, ASETools.SAMLine.ReadFromFile(subfileReader).Where(_ => !_.isUnmapped() && _.mapq >= 10 && !_.isSecondaryAlignment() && _.rname == variant.contig).ToList()); // Skip unmapped, secondary and low-confidence reads

                                subfileReader.Close();
                            } // dna/rna
                        } // if we don't have it cached

                        if (loadFailed)
                        {
                            break;
                        }

                        foreach (var dna in ASETools.BothBools)
                        {
                            reads[dna].AddRange(readsForGivenVariant[key][dna]);
                        }
                    } // for each variant

                    if (loadFailed)
                    {
                        continue;
                    }

                    var inPhaseRef = new Dictionary<bool, int>();
                    ASETools.BothBools.ToList().ForEach(dna => inPhaseRef.Add(dna, 0));

                    var inPhaseAlt = new Dictionary<bool, int>();
                    ASETools.BothBools.ToList().ForEach(dna => inPhaseAlt.Add(dna, 0));

                    var outOfPhaseRA = new Dictionary<bool, int>();
                    ASETools.BothBools.ToList().ForEach(dna => outOfPhaseRA.Add(dna, 0));

                    var outOfPhaseAR = new Dictionary<bool, int>();
                    ASETools.BothBools.ToList().ForEach(dna => outOfPhaseAR.Add(dna, 0));

                    var weird = new Dictionary<bool, int>();
                    ASETools.BothBools.ToList().ForEach(dna => weird.Add(dna, 0));

                    foreach (var dna in ASETools.BothBools)
                    {
                        var readsById = reads[dna].GroupByToDict(_ => _.qname);

                        foreach (var id in readsById.Select(_ => _.Key))
                        {
                            var mapsRefThisPair = new int[2];
                            var mapsAltThisPair = new int[2];
                            var mapsOtherThisPair = new int[2];

                            for (var whichVariant = 0; whichVariant < 2; whichVariant++)
                            {
                                var locus = variants[whichVariant].locus;
                                foreach (var read in readsById[id])
                                {
                                    if (read.mappedBases.ContainsKey(locus)) 
                                    {
                                        if (read.mappedBases[locus].ToString() == variants[whichVariant].reference_allele)
                                        {
                                            mapsRefThisPair[whichVariant]++;
                                        } else if (read.mappedBases[locus].ToString() == variants[whichVariant].alt_allele)
                                        {
                                            mapsAltThisPair[whichVariant]++;
                                        } else
                                        {
                                            mapsOtherThisPair[whichVariant]++;
                                        }
                                    }
                                } // foreach read in the set
                            } // foreach variant

                            // Recall that Enumerable.Range has parameters (start, count), so Enumerable.Range(0,2) is {0, 1}.

                            if (Enumerable.Range(0,2).Any(whichVariant => mapsRefThisPair[whichVariant] + mapsAltThisPair[whichVariant] + mapsOtherThisPair[whichVariant] == 0))
                            {
                                //
                                // No coverage on at least one variant, nothing to see here.
                                //
                                continue;
                            }

                            if (Enumerable.Range(0, 2).Any(whichVariant => mapsOtherThisPair[whichVariant] > 0 || mapsRefThisPair[whichVariant] > 0 && mapsAltThisPair[whichVariant] > 0))    // We have an other or a both on a single variant.  Weird.
                            {
                                weird[dna]++;
                            } else if (mapsRefThisPair[0] > 0 && mapsRefThisPair[1] > 0)
                            {
                                inPhaseRef[dna]++;
                            } else if (mapsAltThisPair[0] > 0 && mapsAltThisPair[1] > 0)
                            {
                                inPhaseAlt[dna]++;
                            } else if (mapsRefThisPair[0] > 0 && mapsAltThisPair[1] > 0)
                            {
                                outOfPhaseRA[dna]++;
                            } else if (mapsAltThisPair[0] > 0 && mapsRefThisPair[1] > 0)
                            {
                                outOfPhaseAR[dna]++;
                            } else
                            {
                                throw new Exception("Logic error -- this shouldn't be possible.");
                            }
                        } // for each read ID
                    } // dna/rna

                    //
                    // Now decide whether this variant pair is in phase, out of phase, or both, and whether that is supported by evidence from both phases and from both DNA & RNA.
                    //
                    if (ASETools.BothBools.Select(dna => inPhaseRef[dna] + inPhaseAlt[dna] + outOfPhaseAR[dna] + outOfPhaseRA[dna]).Sum() >= minReadsToDecideConsistency)
                    {
                        // Contig\tlocus 0\tlocus 1\tSomatic 0\tSomatic 1\tRef 0\tAlt 0\tRef 1\tAlt 1\tDNA ref ref\tDNA alt alt\tDNA ref alt\tDNA alt ref\tDNA weird\tRNA ref ref\tRNA alt alt\tRNA ref alt\tRNA alt ref\tRNA weird\tDNA consistent\tRNA consistent\tDNA consistent with RNA
                        var outputLine = variant.contig + "\t" + variant.locus + "\t" + matchingVariant.locus + "\t" + variant.somaticMutation + "\t" + matchingVariant.somaticMutation + "\t";
                        outputLine += variant.reference_allele + "\t" + variant.alt_allele + "\t" + matchingVariant.reference_allele + "\t" + matchingVariant.alt_allele;

                        foreach (var dna in ASETools.BothBools)
                        {
                            outputLine += "\t" + inPhaseRef[dna] + "\t" + inPhaseAlt[dna] + "\t" + outOfPhaseRA[dna] + "\t" + outOfPhaseAR[dna] + "\t" + weird[dna];
                        }

                        outputLine += "\t" + isConsistent(true, inPhaseRef, inPhaseAlt, outOfPhaseRA, outOfPhaseAR, weird) + "\t" + isConsistent(false, inPhaseRef, inPhaseAlt, outOfPhaseRA, outOfPhaseAR, weird);
                        outputLine += "\t" + DNAandRNAConsistent(inPhaseRef, inPhaseAlt, outOfPhaseRA, outOfPhaseAR, weird);

                        outputLines.Add(outputLine);
                    }
                } // foreach matching variant for the first chosen variant

                if (readsForGivenVariant.ContainsKey(variant.contig + ":" + variant.locus))
                {
                    //
                    // Since we go sequentially, we'll never need this again.  And they can get very big.  Dump it from the cache
                    // so that gc can get it.
                    //
                    readsForGivenVariant.Remove(variant.contig + ":" + variant.locus);
                }
            } // foreach variant for this case


            var outputFilename = ASETools.GetDirectoryFromPathname(case_.tentative_annotated_selected_variants_filename) + @"\" + case_.case_id + ASETools.variantPhasingExtension;
            var outputFile = ASETools.CreateStreamWriterWithRetry(outputFilename);
            if (null == outputFile)
            {
                Console.WriteLine("Unable to open " + outputFilename);
                return;
            }

            outputFile.WriteLine("Contig\tlocus 0\tlocus 1\tSomatic 0\tSomatic 1\tRef 0\tAlt 0\tRef 1\tAlt 1\tDNA ref ref\tDNA alt alt\tDNA ref alt\tDNA alt ref\tDNA weird\tRNA ref ref\tRNA alt alt\tRNA ref alt\tRNA alt ref\tRNA weird" +
                "\tDNA consistent\tRNA consistent\tDNA consistent with RNA");
            outputLines.ForEach(_ => outputFile.WriteLine(_));

            outputFile.WriteLine("**done**");
            outputFile.Close();
        } // HandleOneCase

        static bool isConsistent(bool dna, Dictionary<bool, int> inPhaseRef, Dictionary<bool, int> inPhaseAlt, Dictionary<bool, int> outOfPhaseRA, Dictionary<bool, int> outOfPhaseAR, Dictionary<bool, int> weird)
        {
            int totalMeasurements = inPhaseRef[dna] + inPhaseAlt[dna] + outOfPhaseRA[dna] + outOfPhaseRA[dna] + weird[dna];
            if (inPhaseRef[dna] + outOfPhaseRA[dna] < minReadsToDecideConsistency || inPhaseAlt[dna] + outOfPhaseAR[dna] < minReadsToDecideConsistency || 
                inPhaseRef[dna] + outOfPhaseAR[dna] < minReadsToDecideConsistency || inPhaseAlt[dna] + outOfPhaseRA[dna] < minReadsToDecideConsistency) // We have at least minReadsToDecideConsistency cases in the ref and alt of at least one variant
            {
                return true;    // Nothing to see here
            }

            double inPhaseFraction = ((double)inPhaseRef[dna] + inPhaseAlt[dna]) / totalMeasurements;
            return (inPhaseFraction <= 0.1 || inPhaseFraction >= 0.9) && (double)weird[dna] / totalMeasurements < 0.1;
        } // isConsistent

        static bool DNAandRNAConsistent(Dictionary<bool, int> inPhaseRef, Dictionary<bool, int> inPhaseAlt, Dictionary<bool, int> outOfPhaseRA, Dictionary<bool, int> outOfPhaseAR, Dictionary<bool, int> weird)
        {
            if (ASETools.BothBools.Any(dna => !isConsistent(dna, inPhaseRef, inPhaseAlt, outOfPhaseRA, outOfPhaseAR, weird)))
            {
                return false;
            }

            if (ASETools.BothBools.Any(dna => inPhaseRef[dna] + outOfPhaseRA[dna] < minReadsToDecideConsistency || inPhaseAlt[dna] + outOfPhaseAR[dna] < minReadsToDecideConsistency ||
                inPhaseRef[dna] + outOfPhaseAR[dna] < minReadsToDecideConsistency || inPhaseAlt[dna] + outOfPhaseRA[dna] < minReadsToDecideConsistency))
            {
                return true;    // Anything is consistent with too little measurement.
            }

            var dnaIsInPhase = inPhaseRef[true] + inPhaseAlt[true] > outOfPhaseAR[true] + outOfPhaseRA[true];
            var rnaIsInPhase = inPhaseRef[false] + inPhaseAlt[false] > outOfPhaseAR[false] + outOfPhaseRA[false];

            return dnaIsInPhase == rnaIsInPhase;
        } // DNAandRNAConsistent

    } // Program
} // Namespace
