using ASELib;
using System;
using System.Collections.Generic;
using System.Diagnostics;
using System.Linq;

namespace ChooseAnnotatedVariants
{
    class Program
    {
        static ASETools.Configuration configuration;
        static ASETools.ASERepetitiveRegionMap repetitiveRegionMap;
        static PerThreadState globalState = new PerThreadState();
        static ASETools.GeneLocationsByNameAndChromosome geneLocationInformation;
        static ASETools.GeneMap geneMap;

        static void Main(string[] args)
        {
            Stopwatch timer = new Stopwatch();
            timer.Start();

            configuration = ASETools.Configuration.loadFromFile(args);
            if (null == configuration)
            {
                Console.WriteLine("Unable to load configuration.");
                return;
            }

            geneLocationInformation = new ASETools.GeneLocationsByNameAndChromosome(ASETools.readKnownGeneFile(configuration.geneLocationInformationFilename));
            geneMap = new ASETools.GeneMap(geneLocationInformation.genesByName);

            var cases = ASETools.Case.LoadCases(configuration.casesFilePathname);
            if (null == cases)
            {
                Console.WriteLine("Unable to load cases");
                return;
            }

            repetitiveRegionMap = ASETools.ASERepetitiveRegionMap.loadFromFile(configuration.redundantChromosomeRegionFilename);
            if (null == repetitiveRegionMap)
            {
                Console.WriteLine("Unable to load repetitive region map.");
                return;
            }

            var listOfCases = cases.Select(_ => _.Value).Where(_ => _.tentative_annotated_selected_variants_filename != "").ToList();
            if (listOfCases.Count() != cases.Count())
            {
                Console.WriteLine("At least one case is missing tentative annotated selected variants.");
                //off for debugging. return;
            }

            int nPerDot;
            ASETools.PrintMessageAndNumberBar("Loading tentative variants for", "cases", listOfCases.Count(), out nPerDot);
            var threading = new ASETools.WorkerThreadHelper<ASETools.Case, PerThreadState>(listOfCases, HandleOneCase, FinishUp, null, nPerDot);
            threading.run();
            Console.WriteLine();

            var summaryFilename = configuration.finalResultsDirectory + ASETools.variantSelectionSummaryFilename;
            var summaryFile = ASETools.CreateStreamWriterWithRetry(summaryFilename);
            if (null == summaryFile)
            {
                Console.WriteLine("Unable to open " + summaryFilename + ".  Skipping.");
            } else
            {
                var depthHistogram = new ASETools.PreBucketedHistogram(0, 1000, 1);
                var nearbyHistogram = new ASETools.PreBucketedHistogram(0, 1000, 1);

                foreach (var chr in globalState.variantMap.Select(_ => _.Key).ToList())
                {
                    foreach (var roundedLocus in globalState.variantMap[chr].Select(_ => _.Key).ToList())
                    {
                        foreach (var locus in globalState.variantMap[chr][roundedLocus].Select(_ => _.Key))
                        {
                            depthHistogram.addValue(globalState.variantMap[chr][roundedLocus][locus].Count());

                            int actualDistance;
                            int nNearby = globalState.getVariantsInProximityToALocus(chr, locus, configuration.maxProximityForReferenceBiasCheck, int.MaxValue, out actualDistance).Count();
                            nearbyHistogram.addValue(nNearby);
                        } // locus
                    } // rounded locus
                } // chr

                summaryFile.WriteLine("Distribution of tentative annotated variants by count at the same locus.");
                summaryFile.WriteLine(ASETools.HistogramResultLine.Header());
                depthHistogram.ComputeHistogram().ToList().ForEach(_ => summaryFile.WriteLine(_.ToString()));

                summaryFile.WriteLine("Distribution of tentative annotated variants by count within " + configuration.maxProximityForReferenceBiasCheck + " of variant.");
                summaryFile.WriteLine(ASETools.HistogramResultLine.Header());
                nearbyHistogram.ComputeHistogram().ToList().ForEach(_ => summaryFile.WriteLine(_.ToString()));

                summaryFile.WriteLine("**done**");
            }

            Console.WriteLine(ASETools.ElapsedTimeInSeconds(timer));

            ASETools.PrintMessageAndNumberBar("Selecting and writing variants for", "cases", listOfCases.Count(), out nPerDot);
            var threading2 = new ASETools.WorkerThreadHelper<ASETools.Case, SelectVariantsPerThreadState>(listOfCases, SelectVariantsForOneCase, FinishUpSelectVariants, null, nPerDot);
            threading2.run();

            Console.WriteLine();
            Console.WriteLine("Processed " + listOfCases.Count() + " cases with " + globalCounts.nTentativeVariants + " tentative variants in " + ASETools.ElapsedTimeInSeconds(timer));

            Console.WriteLine("Accepted " + globalCounts.nIncluded + " germline and " + globalCounts.nSomatic + " somatic variants (" + globalCounts.nIncludedForNoObservedBias + " for lack of observed bias in either direction), rejected " + globalCounts.nExcludedByIsASE +
                " because of IsASE, " + globalCounts.nExcludedForObservedBias + " because of observed bias, " + globalCounts.nExcludedByMappedReads +
                " because of repetitive mapped reads, and " + globalCounts.nExcludedByProximityToOtherVariants + " because of proximity to other variants.");

            Console.WriteLine("Acceptance rate " + (double)globalCounts.nIncluded / globalCounts.nTentativeVariants);
        }

        const int treeGranularity = 100;


        class PerThreadState
        {
            public Dictionary<string, Dictionary<int, Dictionary<int, List<double>>>> variantMap = new Dictionary<string, Dictionary<int, Dictionary<int, List<double>>>>(); // chr->locus rounded down to tree granularity->locus->VAFs

            public PerThreadState()
            {
                for (int chr = 1; chr <= ASETools.nHumanAutosomes; chr++)
                {
                    variantMap.Add("chr" + chr, new Dictionary<int, Dictionary<int, List<double>>>());
                }

                variantMap.Add("chrX", new Dictionary<int, Dictionary<int, List<double>>>());
                variantMap.Add("chrY", new Dictionary<int, Dictionary<int, List<double>>>());
            }

            public List<double> getVariantsInProximityToALocus(string chr, int locus, int maxDistance, int minNeeded, out int actualDistance)
            {
                if (maxDistance > treeGranularity)
                {
                    throw new Exception("getVariantsInProximityToALocus: looking too far afield: " + maxDistance + " > " + treeGranularity);
                }


                var retVal = new List<double>();

                int roundedLocus = (locus / treeGranularity) * treeGranularity;

                // Special case distance 0, since there's only one of them.
                if (variantMap[chr].ContainsKey(roundedLocus) && variantMap[chr][roundedLocus].ContainsKey(locus))
                {
                    retVal.AddRange(variantMap[chr][roundedLocus][locus]);
                    if (retVal.Count() >= minNeeded)
                    {
                        actualDistance = 0;
                        return retVal;
                    }
                }

                for (int distance = 1; distance <= maxDistance; distance++)
                {
                    for (int multiplier = -1; multiplier <= 1; multiplier += 2) // fancy way of doing -1, 1
                    {
                        int candidateLocus = locus + distance * multiplier;
                        int roundedCandidateLocus = (candidateLocus / treeGranularity) * treeGranularity;

                        if (variantMap[chr].ContainsKey(roundedCandidateLocus) && variantMap[chr][roundedCandidateLocus].ContainsKey(candidateLocus))
                        {
                            retVal.AddRange(variantMap[chr][roundedCandidateLocus][candidateLocus]);
                        }
                    } // multiplier

                    if (retVal.Count() >= minNeeded)
                    {
                        actualDistance = distance;
                        return retVal;
                    }
                } // distance

                actualDistance = maxDistance;
                return retVal;
            }
        }

        static void FinishUp(PerThreadState perThreadState)
        {
            lock (globalState)
            {
                foreach (var chr in perThreadState.variantMap.Select(_ => _.Key).ToList())
                {
                    foreach (var roundedLocus in perThreadState.variantMap[chr].Select(_ => _.Key).ToList())
                    {
                        if (!globalState.variantMap[chr].ContainsKey(roundedLocus))
                        {
                            globalState.variantMap[chr].Add(roundedLocus, new Dictionary<int, List<double>>());
                        }

                        foreach (var locus in perThreadState.variantMap[chr][roundedLocus].Select(_ => _.Key).ToList())
                        {
                            if (!globalState.variantMap[chr][roundedLocus].ContainsKey(locus))
                            {
                                globalState.variantMap[chr][roundedLocus].Add(locus, new List<double>());
                            }

                            globalState.variantMap[chr][roundedLocus][locus].AddRange(perThreadState.variantMap[chr][roundedLocus][locus]);
                        } // locus
                    } // rounded locus
                } // chr
            } // lock

            perThreadState.variantMap = null;
        } // FinishUp

        static void HandleOneCase(ASETools.Case case_, PerThreadState perThreadState)
        {
            var annotatedVariants = ASETools.AnnotatedVariant.readFile(case_.tentative_annotated_selected_variants_filename);

            if (annotatedVariants == null)
            {
                Console.WriteLine("Unable to load tentative annotated variants from " + case_.tentative_annotated_selected_variants_filename);
                return;
            }

            var copyNumber = ASETools.CopyNumberVariation.ReadBothFiles(case_);
            if (null == copyNumber)
            {
                Console.WriteLine("Unable to load copy number for case " + case_.case_id);
                return;
            }

            foreach (var annotatedVariant in annotatedVariants)
            {
                string whyNot;
                if (!annotatedVariant.IsASECandidate(out whyNot, true, copyNumber, configuration, null, geneMap))
                {
                    continue;
                }

                int roundedLocus = (annotatedVariant.locus / treeGranularity) * treeGranularity;

                if (!perThreadState.variantMap[annotatedVariant.contig].ContainsKey(roundedLocus))
                {
                    perThreadState.variantMap[annotatedVariant.contig].Add(roundedLocus, new Dictionary<int, List<double>>());
                }

                if (!perThreadState.variantMap[annotatedVariant.contig][roundedLocus].ContainsKey(annotatedVariant.locus))
                {
                    perThreadState.variantMap[annotatedVariant.contig][roundedLocus].Add(annotatedVariant.locus, new List<double>());
                }
                perThreadState.variantMap[annotatedVariant.contig][roundedLocus][annotatedVariant.locus].Add(annotatedVariant.GetTumorAltAlleleFraction());
            } // variant
        } // HandleOneCase

        class SelectVariantsPerThreadState
        {
            public int nTentativeVariants = 0;
            public int nExcludedByIsASE = 0;
            public int nIncludedForNoObservedBias = 0;
            public int nExcludedForObservedBias = 0;
            public int nExcludedByMappedReads = 0;
            public int nExcludedByProximityToOtherVariants = 0;
            public int nIncluded = 0;
            public int nSomatic = 0;

            public void merge(SelectVariantsPerThreadState peer)
            {
                nTentativeVariants += peer.nTentativeVariants;
                nExcludedByIsASE += peer.nExcludedByIsASE;
                nIncludedForNoObservedBias += peer.nIncludedForNoObservedBias;
                nExcludedForObservedBias += peer.nExcludedForObservedBias;
                nExcludedByMappedReads += peer.nExcludedByMappedReads;
                nExcludedByProximityToOtherVariants += peer.nExcludedByProximityToOtherVariants;
                nIncluded += peer.nIncluded;
                nSomatic += peer.nSomatic;
            }
        }

        static SelectVariantsPerThreadState globalCounts = new SelectVariantsPerThreadState();

        static void SelectVariantsForOneCase(ASETools.Case case_, SelectVariantsPerThreadState state)
        {
            //
            // Now go back through the tentative variants and select the ones that we want to actually keep.  We throw them
            // out on three bases: being in reference-biased regions; being otherwise not ASE-quality (for instance for copy number
            // variations or DNA read counts being off) and being too close to other selected variants.  The other conditions were already
            // checked when the variants were initially selected.
            //

            var tentativeAnnotatedSelectedVariants = ASETools.AnnotatedVariant.readFile(case_.tentative_annotated_selected_variants_filename);

            if (tentativeAnnotatedSelectedVariants == null)
            {
                Console.WriteLine("Unable to read tentative annotated selected variants from file " + case_.tentative_annotated_selected_variants_filename);
                return;
            }

            var copyNumber = ASETools.CopyNumberVariation.ReadBothFiles(case_);
            if (null == copyNumber)
            {
                Console.WriteLine("Unable to read copy number for case " + case_.case_id);
                return;
            }

            state.nTentativeVariants += tentativeAnnotatedSelectedVariants.Count();

            var annotatedSelectedVariants = tentativeAnnotatedSelectedVariants.Where(_ => _.somaticMutation).ToList();  // The final result.  We take all somatic mutations, so we start with them.

            //
            // Filter the ones with bad read counts, DNA imbalance or copy number variants.  We alrady put the somatic mutations in the final list, so cull them
            // here for the rest of the processing.
            //
            var aseCandidates = tentativeAnnotatedSelectedVariants.Where(_ => _.IsASECandidate(true, copyNumber, configuration, null, geneMap) && !_.somaticMutation).ToList();
            state.nExcludedByIsASE += tentativeAnnotatedSelectedVariants.Where(_ => !_.somaticMutation).Count() - aseCandidates.Count();

            var afterRepetitiveFilter = new List<ASETools.AnnotatedVariant>();
            foreach (var tentativeVariant in aseCandidates)
            {
                //
                // Figure out if it's in a repetitive region.  We do this in two ways.  First, we compute the 95% confidence interval of the mean of the VAF for all samples at this location (or all samples within 25 bases
                // of this location) and see if that fits entirely within the range of 0.4-0.6.  If not, then we look at the repetitive regions map that was generated by mapping the reference back against itself.
                //
                bool passedIntervalTest = false;
                bool failedIntervalTest = false;

                for (int genomeWidth = 0; genomeWidth <= configuration.maxProximityForReferenceBiasCheck; genomeWidth += configuration.maxProximityForReferenceBiasCheck) // fancy way of having genomeWidth be 0 then max
                {
                    int actualDistance;
                    double mean, confidenceIntervalWidth;
                    var samples = globalState.getVariantsInProximityToALocus(tentativeVariant.contig, tentativeVariant.locus, genomeWidth, configuration.minDepthForReferenceBiasCheck, out actualDistance);
                    if (samples.Count() < configuration.minDepthForReferenceBiasCheck)
                    {
                        continue;
                    }

                    ASETools.ComputeConfidenceInterval(samples, configuration.repetitiveRegionConfidence, out mean, out confidenceIntervalWidth);
                    if (mean - confidenceIntervalWidth >= 0.40 && mean + confidenceIntervalWidth <= 0.60)
                    {
                        passedIntervalTest = true;
                        break;
                    } else if (mean + confidenceIntervalWidth < 0.40 || mean - confidenceIntervalWidth > 0.60)
                    {
                        //
                        // The 95% confidence interval is entirely outside the 0.4-0.6 range.  Exclude this variant.
                        //
                        failedIntervalTest = true;
                        break;
                    }
                } // genomeWidth

                if (passedIntervalTest)
                {
                    state.nIncludedForNoObservedBias++;
                    afterRepetitiveFilter.Add(tentativeVariant);
                }
                else
                {
                    if (failedIntervalTest)
                    {
                        state.nExcludedForObservedBias++;
                    } 
                    else if (repetitiveRegionMap.isCloseToRepetitiveRegion(ASETools.ChromosomeNameToIndex(tentativeVariant.contig), tentativeVariant.locus,
                        configuration.minDistanceFromRepetitiveRegion))
                    {
                        state.nExcludedByMappedReads++;
                    }
                    else
                    {
                        afterRepetitiveFilter.Add(tentativeVariant);
                    }
                } // !passedIntervalTest
            } // tentative variant

            var acceptedTree = new ASETools.AVLTree<ASETools.AnnotatedVariant>(); 

            //
            // Now weed out the ones that are too close to others that are still candidates, and write the remainder to the ASV file.  Just be greedy; we could
            // optimize a little and get a handful of extra ASV's, but in real life the mean distance between them is more than 1000x the min, so mostly it
            // won't matter much.
            //

            foreach (var candidate in afterRepetitiveFilter)
            {
                ASETools.AnnotatedVariant nearby;
                if (acceptedTree.FindFirstGreaterThanOrEqualTo(candidate, out nearby) && nearby.contig == candidate.contig && Math.Abs(nearby.locus - candidate.locus) < configuration.minDistanceBetweenGermlineVariants ||
                    acceptedTree.FindFirstLessThanOrEqualTo   (candidate, out nearby) && nearby.contig == candidate.contig && Math.Abs(nearby.locus - candidate.locus) < configuration.minDistanceBetweenGermlineVariants)
                {
                    state.nExcludedByProximityToOtherVariants++;
                }
                else
                {
                    acceptedTree.Insert(candidate);
                    annotatedSelectedVariants.Add(candidate);
                    state.nIncluded++;
                }
            }

            ASETools.AnnotatedVariant.writeFile(ASETools.GetDirectoryFromPathname(case_.tentative_annotated_selected_variants_filename) + @"\" +
                case_.case_id + ASETools.annotatedSelectedVariantsExtension, annotatedSelectedVariants);

            
        } // SelectVariantsForOneCase

        static void FinishUpSelectVariants(SelectVariantsPerThreadState state)
        {
            lock (globalCounts)
            {
                globalCounts.merge(state);
            }
        } // FinishUpSelectVariants
    }
}
