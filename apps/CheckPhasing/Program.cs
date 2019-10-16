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

            var outputFilename = commonData.configuration.finalResultsDirectory + ASETools.PhasingForNearbyVariantsFilename;
            var outputFile = ASETools.CreateStreamWriterWithRetry(outputFilename);
            if (null == outputFile)
            {
                Console.WriteLine("Unable to open outputfile " + outputFilename);
                return;
            }

            outputFile.WriteLine("**done**");
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
                var matchingVariants = tentativeAnnotatedVariants.Where(_ => _.contig == variant.contig && _.locus > variant.locus && _.locus <= variant.locus + commonData.configuration.distanceFromVariantToKeepReads);

                foreach (var matchingVariant in matchingVariants)
                {
                    var variants = new ASETools.AnnotatedVariant[2];
                    variants[0] = variant;
                    variants[1] = matchingVariant;

                    var reads = new Dictionary<bool, List<ASETools.SAMLine>>();    // dna -> reads

                    for (int i = 0; i < 2; i++)
                    {
                        var variantToLoad = variants[i];
                        var key = variantToLoad.contig + ":" + variantToLoad.locus;
                        if (!readsForGivenVariant.ContainsKey(variantToLoad.contig + ":" + variantToLoad.locus))
                        {
                            readsForGivenVariant.Add(key, new Dictionary<bool, List<ASETools.SAMLine>>());
                            foreach (var dna in ASETools.BothBools)
                            {
                                var subfileReader = readers[dna].getSubfile((dna ? case_.tumor_dna_file_id : case_.tumor_rna_file_id) + variantToLoad.getExtractedReadsExtension());
                                if (subfileReader == null)
                                {
                                    Console.WriteLine("Unable to load subfile reader for " + (dna ? case_.tumor_dna_file_id : case_.tumor_rna_file_id) + variantToLoad.getExtractedReadsExtension());
                                    return;
                                }
                                readsForGivenVariant[key].Add(dna, ASETools.SAMLine.ReadFromFile(subfileReader).Where(_ => !_.isUnmapped() && _.mapq >= 10 && !_.isSecondaryAlignment() && _.rname == variant.contig).ToList()); // Skip unmapped, secondary and low-confidence reads

                                subfileReader.Close();
                            } // dna/rna
                        } // if we don't have it cached

                        foreach (var dna in ASETools.BothBools)
                        {
                            reads[dna].AddRange(readsForGivenVariant[key][dna]);
                        }
                    } // for each variant


                    foreach (var dna in ASETools.BothBools)
                    {
                        var readsById = reads[dna].GroupByToDict(_ => _.qname);

                        var mapsRef = new int[2];
                        var mapsAlt = new int[2];
                        var mapsOther = new int[2];

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
                        }
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


        } // HandleOneCase
    } // Program
} // Namespace
