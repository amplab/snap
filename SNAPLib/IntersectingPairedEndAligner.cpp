/*++

Module Name:

    IntersectingPairedEndAligner.cpp

Abstract:

    A paired-end aligner based on set intersections to narrow down possible candidate locations.

Authors:

    Bill Bolosky, February, 2013

Environment:

    User mode service.

Revision History:

--*/

#include "stdafx.h"
#include "IntersectingPairedEndAligner.h"
#include "SeedSequencer.h"
#include "mapq.h"
#include "exit.h"
#include "Error.h"
#include "BigAlloc.h"
#include "AlignerOptions.h"

#ifdef  _DEBUG
extern bool _DumpAlignments;    // From BaseAligner.cpp
#endif  // _DEBUG

IntersectingPairedEndAligner::IntersectingPairedEndAligner(
        GenomeIndex  *index_,
        unsigned      maxReadSize_,
        unsigned      maxHits_,
        unsigned      maxK_,
        unsigned      numSeedsFromCommandLine_,
        double        seedCoverage_,
        unsigned      minSpacing_,                 // Minimum distance to allow between the two ends.
        unsigned      maxSpacing_,                 // Maximum distance to allow between the two ends.
        unsigned      maxBigHits_,
        unsigned      extraSearchDepth_,
        unsigned      maxCandidatePoolSize,
        int           maxSecondaryAlignmentsPerContig_,
        BigAllocator  *allocator,
        bool          noUkkonen_,
        bool          noOrderedEvaluation_,
		bool          noTruncation_) :
    index(index_), maxReadSize(maxReadSize_), maxHits(maxHits_), maxK(maxK_), numSeedsFromCommandLine(__min(MAX_MAX_SEEDS,numSeedsFromCommandLine_)), minSpacing(minSpacing_), maxSpacing(maxSpacing_),
	landauVishkin(NULL), reverseLandauVishkin(NULL), maxBigHits(maxBigHits_), seedCoverage(seedCoverage_),
    extraSearchDepth(extraSearchDepth_), nLocationsScored(0), noUkkonen(noUkkonen_), noOrderedEvaluation(noOrderedEvaluation_), noTruncation(noTruncation_), 
    maxSecondaryAlignmentsPerContig(maxSecondaryAlignmentsPerContig_)
{
    doesGenomeIndexHave64BitLocations = index->doesGenomeIndexHave64BitLocations();

    unsigned maxSeedsToUse;
    if (0 != numSeedsFromCommandLine) {
        maxSeedsToUse = numSeedsFromCommandLine;
    } else {
        maxSeedsToUse = (unsigned)(maxReadSize * seedCoverage / index->getSeedLength());
    }
    allocateDynamicMemory(allocator, maxReadSize, maxBigHits, maxSeedsToUse, maxK, extraSearchDepth, maxCandidatePoolSize, maxSecondaryAlignmentsPerContig);

    rcTranslationTable['A'] = 'T';
    rcTranslationTable['G'] = 'C';
    rcTranslationTable['C'] = 'G';
    rcTranslationTable['T'] = 'A';
    rcTranslationTable['N'] = 'N';

    for (unsigned i = 0; i < 256; i++) {
        nTable[i] = 0;
    }

    nTable['N'] = 1;

    seedLen = index->getSeedLength();

    genome = index->getGenome();
    genomeSize = genome->getCountOfBases();
}

IntersectingPairedEndAligner::~IntersectingPairedEndAligner()
{
}

    size_t
IntersectingPairedEndAligner::getBigAllocatorReservation(GenomeIndex * index, unsigned maxBigHitsToConsider, unsigned maxReadSize, unsigned seedLen, unsigned numSeedsFromCommandLine,
                                                         double seedCoverage, unsigned maxEditDistanceToConsider, unsigned maxExtraSearchDepth, unsigned maxCandidatePoolSize,
                                                         int maxSecondaryAlignmentsPerContig)
{
    unsigned maxSeedsToUse;
    if (0 != numSeedsFromCommandLine) {
        maxSeedsToUse = numSeedsFromCommandLine;
    } else {
        maxSeedsToUse = (unsigned)(maxReadSize * seedCoverage / index->getSeedLength());
    }
    CountingBigAllocator countingAllocator;
    {
        IntersectingPairedEndAligner aligner; // This has to be in a nested scope so its destructor is called before that of the countingAllocator
        aligner.index = index;
        aligner.doesGenomeIndexHave64BitLocations = index->doesGenomeIndexHave64BitLocations();

        aligner.allocateDynamicMemory(&countingAllocator, maxReadSize, maxBigHitsToConsider, maxSeedsToUse, maxEditDistanceToConsider, maxExtraSearchDepth, maxCandidatePoolSize,
            maxSecondaryAlignmentsPerContig);
        return sizeof(aligner) + countingAllocator.getMemoryUsed();
    }
}

    void
IntersectingPairedEndAligner::allocateDynamicMemory(BigAllocator *allocator, unsigned maxReadSize, unsigned maxBigHitsToConsider, unsigned maxSeedsToUse,
                                                    unsigned maxEditDistanceToConsider, unsigned maxExtraSearchDepth, unsigned maxCandidatePoolSize,
                                                    int maxSecondaryAlignmentsPerContig)
{
    seedUsed = (BYTE *) allocator->allocate(100 + (maxReadSize + 7) / 8);

    for (unsigned whichRead = 0; whichRead < NUM_READS_PER_PAIR; whichRead++) {
        rcReadData[whichRead] = (char *)allocator->allocate(maxReadSize);
        rcReadQuality[whichRead] = (char *)allocator->allocate(maxReadSize);

        for (Direction dir = 0; dir < NUM_DIRECTIONS; dir++) {
            reversedRead[whichRead][dir] = (char *)allocator->allocate(maxReadSize);
            hashTableHitSets[whichRead][dir] =(HashTableHitSet *)allocator->allocate(sizeof(HashTableHitSet)); /*new HashTableHitSet();*/
            hashTableHitSets[whichRead][dir]->firstInit(maxSeedsToUse, maxMergeDistance, allocator, doesGenomeIndexHave64BitLocations);
        }
    }

    scoringCandidatePoolSize = min(maxCandidatePoolSize, maxBigHitsToConsider * maxSeedsToUse * NUM_READS_PER_PAIR);

    scoringCandidates = (ScoringCandidate **) allocator->allocate(sizeof(ScoringCandidate *) * (maxEditDistanceToConsider + maxExtraSearchDepth + 1));  //+1 is for 0.
    scoringCandidatePool = (ScoringCandidate *)allocator->allocate(sizeof(ScoringCandidate) * scoringCandidatePoolSize);

    for (unsigned i = 0; i < NUM_READS_PER_PAIR; i++) {
        scoringMateCandidates[i] = (ScoringMateCandidate *) allocator->allocate(sizeof(ScoringMateCandidate) * scoringCandidatePoolSize / NUM_READS_PER_PAIR);
    }

    mergeAnchorPoolSize = scoringCandidatePoolSize;
    mergeAnchorPool = (MergeAnchor *)allocator->allocate(sizeof(MergeAnchor) * mergeAnchorPoolSize);

    if (maxSecondaryAlignmentsPerContig > 0) {
        size_t size = sizeof(*hitsPerContigCounts) * index->getGenome()->getNumContigs();
        hitsPerContigCounts = (HitsPerContigCounts *)allocator->allocate(size);
        memset(hitsPerContigCounts, 0, size);
        contigCountEpoch = 0;
    } else {
        hitsPerContigCounts = NULL;
    }
}

    void
IntersectingPairedEndAligner::align(
        Read                  *read0,
        Read                  *read1,
        PairedAlignmentResult *result,
        int                    maxEditDistanceForSecondaryResults,
        int                    secondaryResultBufferSize,
        int                   *nSecondaryResults,
        PairedAlignmentResult *secondaryResults,             // The caller passes in a buffer of secondaryResultBufferSize and it's filled in by align()
        int                    singleSecondaryBufferSize,
        int                    maxSecondaryResultsToReturn,
        int                   *nSingleEndSecondaryResultsForFirstRead,
        int                   *nSingleEndSecondaryResultsForSecondRead,
        SingleAlignmentResult *singleEndSecondaryResults     // Single-end secondary alignments for when the paired-end alignment didn't work properly
        )
{
    result->nLVCalls = 0;
    result->nSmallHits = 0;

    *nSecondaryResults = 0;
    *nSingleEndSecondaryResultsForFirstRead = 0;
    *nSingleEndSecondaryResultsForSecondRead = 0;

    int maxSeeds;
    if (numSeedsFromCommandLine != 0) {
        maxSeeds = (int)numSeedsFromCommandLine;
    } else {
        maxSeeds = (int)(max(read0->getDataLength(), read1->getDataLength()) * seedCoverage / index->getSeedLength());
    }

#ifdef  _DEBUG
    if (_DumpAlignments) {
        printf("\nIntersectingAligner aligning reads '%*.s' and '%.*s' with data '%.*s' and '%.*s'\n", read0->getIdLength(), read0->getId(), read1->getIdLength(), read1->getId(), read0->getDataLength(), read0->getData(), read1->getDataLength(), read1->getData());
    }
#endif  // _DEBUG

    lowestFreeScoringCandidatePoolEntry = 0;
    for (unsigned k = 0; k <= maxK + extraSearchDepth; k++) {
        scoringCandidates[k] = NULL;
    }

    for (unsigned i = 0; i < NUM_SET_PAIRS; i++) {
        lowestFreeScoringMateCandidate[i] = 0;
    }
    firstFreeMergeAnchor = 0;

    Read rcReads[NUM_READS_PER_PAIR];

    GenomeLocation bestResultGenomeLocation[NUM_READS_PER_PAIR];
    Direction bestResultDirection[NUM_READS_PER_PAIR];
    unsigned bestResultScore[NUM_READS_PER_PAIR];
    unsigned popularSeedsSkipped[NUM_READS_PER_PAIR];

    reads[0][FORWARD] = read0;
    reads[1][FORWARD] = read1;

    //
    // Don't bother if one or both reads are too short.  The minimum read length here is the seed length, but usually there's a longer
	// minimum enforced by our called
    //
    if (read0->getDataLength() < seedLen || read1->getDataLength() < seedLen) {
         return;
    }

    //
    // Build the RC reads.
    //
    unsigned countOfNs = 0;

    for (unsigned whichRead = 0; whichRead < NUM_READS_PER_PAIR; whichRead++) {
        Read *read = reads[whichRead][FORWARD];
        readLen[whichRead] = read->getDataLength();
        popularSeedsSkipped[whichRead] = 0;
        countOfHashTableLookups[whichRead] = 0;
#if 0
        hitLocations[whichRead]->clear();
        mateHitLocations[whichRead]->clear();
#endif // 0

        for (Direction dir = FORWARD; dir < NUM_DIRECTIONS; dir++) {
            totalHashTableHits[whichRead][dir] = 0;
            largestHashTableHit[whichRead][dir] = 0;
            hashTableHitSets[whichRead][dir]->init();
        }

        if (readLen[whichRead] > maxReadSize) {
            WriteErrorMessage("IntersectingPairedEndAligner:: got too big read (%d > %d)\n"
                              "Change MAX_READ_LENTH at the beginning of Read.h and recompile.\n", readLen[whichRead], maxReadSize);
            soft_exit(1);
        }

        for (unsigned i = 0; i < reads[whichRead][FORWARD]->getDataLength(); i++) {
            rcReadData[whichRead][i] = rcTranslationTable[read->getData()[readLen[whichRead] - i - 1]];
            rcReadQuality[whichRead][i] = read->getQuality()[readLen[whichRead] - i - 1];
            countOfNs += nTable[read->getData()[i]];
        }
        reads[whichRead][RC] = &rcReads[whichRead];
        reads[whichRead][RC]->init(read->getId(), read->getIdLength(), rcReadData[whichRead], rcReadQuality[whichRead], read->getDataLength());
    }

    if (countOfNs > maxK) {
        return;
    }

    //
    // Build the reverse data for both reads in both directions for the backwards LV to use.
    //
    for (unsigned whichRead = 0; whichRead < NUM_READS_PER_PAIR; whichRead++) {
        for (Direction dir = 0; dir < NUM_DIRECTIONS; dir++) {
            Read *read = reads[whichRead][dir];

            for (unsigned i = 0; i < read->getDataLength(); i++) {
                reversedRead[whichRead][dir][i] = read->getData()[read->getDataLength() - i - 1];
            }
        }
    }

    unsigned thisPassSeedsNotSkipped[NUM_READS_PER_PAIR][NUM_DIRECTIONS] = {{0,0}, {0,0}};

    //
    // Initialize the member variables that are effectively stack locals, but are in the object
    // to avoid having to pass them to score.
    //
    double probabilityOfBestPair = 0;
    localBestPairProbability[0] = 0;
    localBestPairProbability[1] = 0;
    double probabilityOfAllPairs = 0;
    unsigned bestPairScore = 65536;
    unsigned scoreLimit = maxK + extraSearchDepth;

    //
    // Phase 1: do the hash table lookups for each of the seeds for each of the reads and add them to the hit sets.
    //
    for (unsigned whichRead = 0; whichRead < NUM_READS_PER_PAIR; whichRead++) {
        int nextSeedToTest = 0;
        unsigned wrapCount = 0;
        int nPossibleSeeds = (int)readLen[whichRead] - seedLen + 1;
        memset(seedUsed, 0, (__max(readLen[0], readLen[1]) + 7) / 8);
        bool beginsDisjointHitSet[NUM_DIRECTIONS] = {true, true};

        while (countOfHashTableLookups[whichRead] < nPossibleSeeds && countOfHashTableLookups[whichRead] < maxSeeds) {
            if (nextSeedToTest >= nPossibleSeeds) {
                wrapCount++;
				beginsDisjointHitSet[FORWARD] = beginsDisjointHitSet[RC] = true;
                if (wrapCount >= seedLen) {
                    //
                    // There aren't enough valid seeds in this read to reach our target.
                    //
                    break;
                }
                nextSeedToTest = GetWrappedNextSeedToTest(seedLen, wrapCount);
            }


            while (nextSeedToTest < nPossibleSeeds && IsSeedUsed(nextSeedToTest)) {
                //
                // This seed is already used.  Try the next one.
                //
                nextSeedToTest++;
            }

            if (nextSeedToTest >= nPossibleSeeds) {
                //
                // Unusable seeds have pushed us past the end of the read.  Go back around the outer loop so we wrap properly.
                //
                continue;
            }

            SetSeedUsed(nextSeedToTest);

            if (!Seed::DoesTextRepresentASeed(reads[whichRead][FORWARD]->getData() + nextSeedToTest, seedLen)) {
                //
                // It's got Ns in it, so just skip it.
                //
                nextSeedToTest++;
                continue;
            }

            Seed seed(reads[whichRead][FORWARD]->getData() + nextSeedToTest, seedLen);
            //
            // Find all instances of this seed in the genome.
            //
            _int64 nHits[NUM_DIRECTIONS];
            const GenomeLocation *hits[NUM_DIRECTIONS];
            const unsigned *hits32[NUM_DIRECTIONS];

            if (doesGenomeIndexHave64BitLocations) {
                index->lookupSeed(seed, &nHits[FORWARD], &hits[FORWARD], &nHits[RC], &hits[RC], 
                            hashTableHitSets[whichRead][FORWARD]->getNextSingletonLocation(), hashTableHitSets[whichRead][RC]->getNextSingletonLocation());
            } else {
                index->lookupSeed32(seed, &nHits[FORWARD], &hits32[FORWARD], &nHits[RC], &hits32[RC]);
            }

            countOfHashTableLookups[whichRead]++;
            for (Direction dir = FORWARD; dir < NUM_DIRECTIONS; dir++) {
                int offset;
                if (dir == FORWARD) {
                    offset = nextSeedToTest;
                } else {
                    offset = readLen[whichRead] - seedLen - nextSeedToTest;
                }
                if (nHits[dir] < maxBigHits) {
                    totalHashTableHits[whichRead][dir] += nHits[dir];
                    if (doesGenomeIndexHave64BitLocations) {
                        hashTableHitSets[whichRead][dir]->recordLookup(offset, nHits[dir], hits[dir], beginsDisjointHitSet[dir]);
                    } else {
                        hashTableHitSets[whichRead][dir]->recordLookup(offset, nHits[dir], hits32[dir], beginsDisjointHitSet[dir]);
                    }
                    beginsDisjointHitSet[dir]= false;
                } else {
                    popularSeedsSkipped[whichRead]++;
                }
            }

            //
            // If we don't have enough seeds left to reach the end of the read, space out the seeds more-or-less evenly.
            //
            if ((maxSeeds - countOfHashTableLookups[whichRead] + 1) * (int)seedLen + nextSeedToTest < nPossibleSeeds) {
                _ASSERT((nPossibleSeeds - nextSeedToTest - 1) / (maxSeeds - countOfHashTableLookups[whichRead] + 1) >= (int)seedLen);
                nextSeedToTest += (nPossibleSeeds - nextSeedToTest - 1) / (maxSeeds - countOfHashTableLookups[whichRead] + 1);
                _ASSERT(nextSeedToTest < nPossibleSeeds);   // We haven't run off the end of the read.
            } else {
                nextSeedToTest += seedLen;
            }
        } // while we need to lookup seeds for this read
    } // for each read

    readWithMoreHits = totalHashTableHits[0][FORWARD] + totalHashTableHits[0][RC] > totalHashTableHits[1][FORWARD] + totalHashTableHits[1][RC] ? 0 : 1;
    readWithFewerHits = 1 - readWithMoreHits;

#ifdef  _DEBUG
    if (_DumpAlignments) {
        printf("Read 0 has %d hits, read 1 has %d hits\n", totalHashTableHits[0][FORWARD] + totalHashTableHits[0][RC], totalHashTableHits[1][FORWARD] + totalHashTableHits[1][RC]);
    }
#endif  // _DEBUG

    Direction setPairDirection[NUM_SET_PAIRS][NUM_READS_PER_PAIR] = {{FORWARD, RC}, {RC, FORWARD}};


    //
    // Phase 2: find all possible candidates and add them to candidate lists (for the reads with fewer and more hits).
    //
    unsigned maxUsedBestPossibleScoreList = 0;

    for (unsigned whichSetPair = 0; whichSetPair < NUM_SET_PAIRS; whichSetPair++) {
        HashTableHitSet *setPair[NUM_READS_PER_PAIR];

        if (whichSetPair == 0) {
            setPair[0] = hashTableHitSets[0][FORWARD];
            setPair[1] = hashTableHitSets[1][RC];
        } else {
            setPair[0] = hashTableHitSets[0][RC];
            setPair[1] = hashTableHitSets[1][FORWARD];
        }


        unsigned            lastSeedOffsetForReadWithFewerHits;
        GenomeLocation      lastGenomeLocationForReadWithFewerHits;
        GenomeLocation      lastGenomeLocationForReadWithMoreHits;
        unsigned            lastSeedOffsetForReadWithMoreHits;

        bool                outOfMoreHitsLocations = false;

        //
        // Seed the intersection state by doing a first lookup.
        //
        if (setPair[readWithFewerHits]->getFirstHit(&lastGenomeLocationForReadWithFewerHits, &lastSeedOffsetForReadWithFewerHits)) {
            //
            // No hits in this direction.
            //
            continue;   // The outer loop over set pairs.
        }

        lastGenomeLocationForReadWithMoreHits = InvalidGenomeLocation;

        //
        // Loop over the candidates in for the read with more hits.  At the top of the loop, we have a candidate but don't know if it has
        // a mate.  Each pass through the loop considers a single hit on the read with fewer hits.
        //
        for (;;) {

            //
            // Loop invariant: lastGenomeLocationForReadWithFewerHits is the highest genome offset that has not been considered.
            // lastGenomeLocationForReadWithMoreHits is also the highest genome offset on that side that has not been
            // considered (or is InvalidGenomeLocation), but higher ones within the appropriate range might already be in scoringMateCandidates.
            // We go once through this loop for each
            //

            if (lastGenomeLocationForReadWithMoreHits > lastGenomeLocationForReadWithFewerHits + maxSpacing) {
                //
                // The more hits side is too high to be a mate candidate for the fewer hits side.  Move it down to the largest
                // location that's not too high.
                //
                if (!setPair[readWithMoreHits]->getNextHitLessThanOrEqualTo(lastGenomeLocationForReadWithFewerHits + maxSpacing,
                                                                             &lastGenomeLocationForReadWithMoreHits, &lastSeedOffsetForReadWithMoreHits)) {
                    break;  // End of all of the mates.  We're done with this set pair.
                }
            }

            if ((lastGenomeLocationForReadWithMoreHits + maxSpacing < lastGenomeLocationForReadWithFewerHits || outOfMoreHitsLocations) &&
                (0 == lowestFreeScoringMateCandidate[whichSetPair] ||
                !genomeLocationIsWithin(scoringMateCandidates[whichSetPair][lowestFreeScoringMateCandidate[whichSetPair]-1].readWithMoreHitsGenomeLocation, lastGenomeLocationForReadWithFewerHits, maxSpacing))) {
                //
                // No mates for the hit on the read with fewer hits.  Skip to the next candidate.
                //
                if (outOfMoreHitsLocations) {
                    //
                    // Nothing left on the more hits side, we're done with this set pair.
                    //
                    break;
                }

                if (!setPair[readWithFewerHits]->getNextHitLessThanOrEqualTo(lastGenomeLocationForReadWithMoreHits + maxSpacing, &lastGenomeLocationForReadWithFewerHits,
                                                        &lastSeedOffsetForReadWithFewerHits)) {
                    //
                    // No more candidates on the read with fewer hits side.  We're done with this set pair.
                    //
                    break;
                }
                continue;
            }

            //
            // Add all of the mate candidates for this fewer side hit.
            //

            GenomeLocation previousMoreHitsLocation = lastGenomeLocationForReadWithMoreHits;
            while (lastGenomeLocationForReadWithMoreHits + maxSpacing >= lastGenomeLocationForReadWithFewerHits && !outOfMoreHitsLocations) {
				unsigned bestPossibleScoreForReadWithMoreHits;
				if (noTruncation) {
					bestPossibleScoreForReadWithMoreHits = 0;
				} else {
					bestPossibleScoreForReadWithMoreHits = setPair[readWithMoreHits]->computeBestPossibleScoreForCurrentHit();
				} 

                if (lowestFreeScoringMateCandidate[whichSetPair] >= scoringCandidatePoolSize / NUM_READS_PER_PAIR) {
                    WriteErrorMessage("Ran out of scoring candidate pool entries.  Perhaps trying with a larger value of -mcp will help.\n");
                    soft_exit(1);
                }
                scoringMateCandidates[whichSetPair][lowestFreeScoringMateCandidate[whichSetPair]].init(
                                lastGenomeLocationForReadWithMoreHits, bestPossibleScoreForReadWithMoreHits, lastSeedOffsetForReadWithMoreHits);

#ifdef _DEBUG
                if (_DumpAlignments) {
                    printf("SetPair %d, added more hits candidate %d at genome location %u, bestPossibleScore %d, seedOffset %d\n",
                            whichSetPair, lowestFreeScoringMateCandidate[whichSetPair], lastGenomeLocationForReadWithMoreHits,
                            bestPossibleScoreForReadWithMoreHits,
                            lastSeedOffsetForReadWithMoreHits);
                }
#endif // _DEBUG

                lowestFreeScoringMateCandidate[whichSetPair]++;

                previousMoreHitsLocation = lastGenomeLocationForReadWithMoreHits;

                if (!setPair[readWithMoreHits]->getNextLowerHit(&lastGenomeLocationForReadWithMoreHits, &lastSeedOffsetForReadWithMoreHits)) {
                    lastGenomeLocationForReadWithMoreHits = 0;
                    outOfMoreHitsLocations = true;
                    break; // out of the loop looking for candidates on the more hits side.
                }
            }

            //
            // And finally add the hit from the fewer hit side.  To compute its best possible score, we need to look at all of the mates; we couldn't do it in the
            // loop immediately above because some of them might have already been in the mate list from a different, nearby fewer hit location.
            //
			unsigned bestPossibleScoreForReadWithFewerHits;
			
			if (noTruncation) {
				bestPossibleScoreForReadWithFewerHits = 0;
			} else {
				bestPossibleScoreForReadWithFewerHits = setPair[readWithFewerHits]->computeBestPossibleScoreForCurrentHit();
			}

            unsigned lowestBestPossibleScoreOfAnyPossibleMate = maxK + extraSearchDepth;
            for (int i = lowestFreeScoringMateCandidate[whichSetPair] - 1; i >= 0; i--) {
                if (scoringMateCandidates[whichSetPair][i].readWithMoreHitsGenomeLocation > lastGenomeLocationForReadWithFewerHits + maxSpacing) {
                    break;
                }
                lowestBestPossibleScoreOfAnyPossibleMate = __min(lowestBestPossibleScoreOfAnyPossibleMate, scoringMateCandidates[whichSetPair][i].bestPossibleScore);
            }

            if (lowestBestPossibleScoreOfAnyPossibleMate + bestPossibleScoreForReadWithFewerHits <= maxK + extraSearchDepth) {
                //
                // There's a set of ends that we can't prove doesn't have too large of a score.  Allocate a fewer hit candidate and stick it in the
                // correct weight list.
                //
                if (lowestFreeScoringCandidatePoolEntry >= scoringCandidatePoolSize) {
                    WriteErrorMessage("Ran out of scoring candidate pool entries.  Perhaps rerunning with a larger value of -mcp will help.\n");
                    soft_exit(1);
                }

                //
                // If we have noOrderedEvaluation set, just stick everything on list 0, regardless of what it really is.  This will cause us to
                // evaluate the candidates in more-or-less inverse genome order.
                //
                unsigned bestPossibleScore = noOrderedEvaluation ? 0 : lowestBestPossibleScoreOfAnyPossibleMate + bestPossibleScoreForReadWithFewerHits;

                scoringCandidatePool[lowestFreeScoringCandidatePoolEntry].init(lastGenomeLocationForReadWithFewerHits, whichSetPair, lowestFreeScoringMateCandidate[whichSetPair] - 1,
                                                                                lastSeedOffsetForReadWithFewerHits, bestPossibleScoreForReadWithFewerHits,
                                                                                scoringCandidates[bestPossibleScore]);


                scoringCandidates[bestPossibleScore] = &scoringCandidatePool[lowestFreeScoringCandidatePoolEntry];

#ifdef _DEBUG
                if (_DumpAlignments) {
                    printf("SetPair %d, added fewer hits candidate %d at genome location %u, bestPossibleScore %d, seedOffset %d\n",
                            whichSetPair, lowestFreeScoringCandidatePoolEntry, lastGenomeLocationForReadWithFewerHits,
                            lowestBestPossibleScoreOfAnyPossibleMate + bestPossibleScoreForReadWithFewerHits,
                            lastSeedOffsetForReadWithFewerHits);
                }
#endif // _DEBUG

                lowestFreeScoringCandidatePoolEntry++;
                maxUsedBestPossibleScoreList = max(maxUsedBestPossibleScoreList, bestPossibleScore);
            }

            if (!setPair[readWithFewerHits]->getNextLowerHit(&lastGenomeLocationForReadWithFewerHits, &lastSeedOffsetForReadWithFewerHits)) {
                break;
            }
        }
    } // For each set pair


    //
    // Phase 3: score and merge the candidates we've found.
    //
    unsigned currentBestPossibleScoreList = 0;
    scoreLimit = maxK + extraSearchDepth;
    //
    // Loop until we've scored all of the candidates, or proven that what's left must have too high of a score to be interesting.
    //
    while (currentBestPossibleScoreList <= maxUsedBestPossibleScoreList && currentBestPossibleScoreList <= scoreLimit) {
        if (scoringCandidates[currentBestPossibleScoreList] == NULL) {
            //
            // No more candidates on this list.  Skip to the next one.
            //
            currentBestPossibleScoreList++;
            continue;
        }


        //
        // Grab the first candidate on the highest list and score it.
        //
        ScoringCandidate *candidate = scoringCandidates[currentBestPossibleScoreList];

        unsigned fewerEndScore;
        double fewerEndMatchProbability;
        int fewerEndGenomeLocationOffset;

        scoreLocation(readWithFewerHits, setPairDirection[candidate->whichSetPair][readWithFewerHits], candidate->readWithFewerHitsGenomeLocation,
            candidate->seedOffset, scoreLimit, &fewerEndScore, &fewerEndMatchProbability, &fewerEndGenomeLocationOffset);

        _ASSERT(-1 == fewerEndScore || fewerEndScore >= candidate->bestPossibleScore);

#ifdef _DEBUG
        if (_DumpAlignments) {
            printf("Scored fewer end candidate %d, set pair %d, read %d, location %u, seed offset %d, score limit %d, score %d, offset %d\n", (int)(candidate - scoringCandidatePool),
                candidate->whichSetPair, readWithFewerHits, candidate->readWithFewerHitsGenomeLocation, candidate->seedOffset,
                scoreLimit, fewerEndScore, fewerEndGenomeLocationOffset);
        }
#endif // DEBUG

        if (fewerEndScore != -1) {
            //
            // Find and score mates.  The index in scoringMateCandidateIndex is the lowest mate (i.e., the highest index number).
            //
            unsigned mateIndex = candidate->scoringMateCandidateIndex;

            for (;;) {

                ScoringMateCandidate *mate = &scoringMateCandidates[candidate->whichSetPair][mateIndex];
                _ASSERT(genomeLocationIsWithin(mate->readWithMoreHitsGenomeLocation, candidate->readWithFewerHitsGenomeLocation, maxSpacing));
                if (!genomeLocationIsWithin(mate->readWithMoreHitsGenomeLocation, candidate->readWithFewerHitsGenomeLocation, minSpacing) && mate->bestPossibleScore <= scoreLimit - fewerEndScore) {
                    //
                    // It's within the range and not necessarily too poor of a match.  Consider it.
                    //

                    //
                    // If we haven't yet scored this mate, or we've scored it and not gotten an answer, but had a higher score limit than we'd
                    // use now, score it.
                    //
                    if (mate->score == -2 || mate->score == -1 && mate->scoreLimit < scoreLimit - fewerEndScore) {
                        scoreLocation(readWithMoreHits, setPairDirection[candidate->whichSetPair][readWithMoreHits], mate->readWithMoreHitsGenomeLocation,
                            mate->seedOffset, scoreLimit - fewerEndScore, &mate->score, &mate->matchProbability,
                            &mate->genomeOffset);
#ifdef _DEBUG
                        if (_DumpAlignments) {
                            printf("Scored mate candidate %d, set pair %d, read %d, location %u, seed offset %d, score limit %d, score %d, offset %d\n",
                                (int)(mate - scoringMateCandidates[candidate->whichSetPair]), candidate->whichSetPair, readWithMoreHits, mate->readWithMoreHitsGenomeLocation,
                                mate->seedOffset, scoreLimit - fewerEndScore, mate->score, mate->genomeOffset);
                        }
#endif // _DEBUG

                        _ASSERT(-1 == mate->score || mate->score >= mate->bestPossibleScore);

                        mate->scoreLimit = scoreLimit - fewerEndScore;
                    }

                    if (mate->score != -1) {
                        double pairProbability = mate->matchProbability * fewerEndMatchProbability;
                        unsigned pairScore = mate->score + fewerEndScore;
                        //
                        // See if this should be ignored as a merge, or if we need to back out a previously scored location
                        // because it's a worse version of this location.
                        //
                        MergeAnchor *mergeAnchor = candidate->mergeAnchor;

                        if (NULL == mergeAnchor) {
                            //
                            // Look up and down the array of candidates to see if we have possible merge candidates.
                            //
                            for (ScoringCandidate *mergeCandidate = candidate - 1;
                                        mergeCandidate >= scoringCandidatePool &&
                                        genomeLocationIsWithin(mergeCandidate->readWithFewerHitsGenomeLocation, candidate->readWithFewerHitsGenomeLocation + fewerEndGenomeLocationOffset, 50) &&
                                        mergeCandidate->whichSetPair == candidate->whichSetPair;
                                        mergeCandidate--) {

                                if (mergeCandidate->mergeAnchor != NULL) {
                                    candidate->mergeAnchor = mergeAnchor = mergeCandidate->mergeAnchor;
                                    break;
                                }
                            }

                            if (NULL == mergeAnchor) {
                                for (ScoringCandidate *mergeCandidate = candidate + 1;
                                            mergeCandidate < scoringCandidatePool + lowestFreeScoringCandidatePoolEntry &&
                                            genomeLocationIsWithin(mergeCandidate->readWithFewerHitsGenomeLocation, candidate->readWithFewerHitsGenomeLocation + fewerEndGenomeLocationOffset, 50) &&
                                            mergeCandidate->whichSetPair == candidate->whichSetPair;
                                            mergeCandidate--) {

                                    if (mergeCandidate->mergeAnchor != NULL) {
                                        candidate->mergeAnchor = mergeAnchor = mergeCandidate->mergeAnchor;
                                        break;
                                    }
                                }
                            }
                        }

                        bool merged;

                        double oldPairProbability;

                        if (NULL == mergeAnchor) {
                            if (firstFreeMergeAnchor >= mergeAnchorPoolSize) {
                                WriteErrorMessage("Ran out of merge anchor pool entries.  Perhaps rerunning with a larger value of -mcp will help\n");
                                soft_exit(1);
                            }

                            mergeAnchor = &mergeAnchorPool[firstFreeMergeAnchor];

                            firstFreeMergeAnchor++;

                            mergeAnchor->init(mate->readWithMoreHitsGenomeLocation + mate->genomeOffset, candidate->readWithFewerHitsGenomeLocation + fewerEndGenomeLocationOffset,
                                pairProbability, pairScore);

                            merged = false;
                            oldPairProbability = 0;
                            candidate->mergeAnchor = mergeAnchor;
                        } else {
                            merged = mergeAnchor->checkMerge(mate->readWithMoreHitsGenomeLocation + mate->genomeOffset, candidate->readWithFewerHitsGenomeLocation + fewerEndGenomeLocationOffset,
                                pairProbability, pairScore, &oldPairProbability);
                        }

                        if (!merged) {
                            //
                            // Back out the probability of the old match that we're merged with, if any.  The max
                            // is necessary because a + b - b is not necessarily a in floating point.  If there
                            // was no merge, the oldPairProbability is 0.
                            //
                            probabilityOfAllPairs = __max(0, probabilityOfAllPairs - oldPairProbability);

                            bool isBestHit = false;

                            if (pairScore <= maxK && (pairScore < bestPairScore ||
                                (pairScore == bestPairScore && pairProbability > probabilityOfBestPair))) {
                                //
                                // A new best hit.
                                //
                                if (maxEditDistanceForSecondaryResults != -1 && (unsigned)maxEditDistanceForSecondaryResults >= pairScore - bestPairScore) {
                                    //
                                    // Move the old best to be a secondary alignment.  This won't happen on the first time we get a valid alignment,
                                    // because bestPairScore is initialized to be very large.
                                    //
                                    //
                                    if (*nSecondaryResults >= secondaryResultBufferSize) {
                                        WriteErrorMessage("IntersectingPairedEndAligner::align(): out of secondary result buffer\n");
                                        soft_exit(1);
                                    }

                                    PairedAlignmentResult *result = &secondaryResults[*nSecondaryResults];
                                    result->alignedAsPair = true;
                                    result->fromAlignTogether = true;

                                    for (int r = 0; r < NUM_READS_PER_PAIR; r++) {
                                        result->direction[r] = bestResultDirection[r];
                                        result->location[r] = bestResultGenomeLocation[r];
                                        result->mapq[r] = 0;
                                        result->score[r] = bestResultScore[r];
                                        result->status[r] = MultipleHits;
                                    }
 
                                    (*nSecondaryResults)++;

                                }
                                bestPairScore = pairScore;
                                probabilityOfBestPair = pairProbability;
                                bestResultGenomeLocation[readWithFewerHits] = candidate->readWithFewerHitsGenomeLocation + fewerEndGenomeLocationOffset;
                                bestResultGenomeLocation[readWithMoreHits] = mate->readWithMoreHitsGenomeLocation + mate->genomeOffset;
                                bestResultScore[readWithFewerHits] = fewerEndScore;
                                bestResultScore[readWithMoreHits] = mate->score;
                                bestResultDirection[readWithFewerHits] = setPairDirection[candidate->whichSetPair][readWithFewerHits];
                                bestResultDirection[readWithMoreHits] = setPairDirection[candidate->whichSetPair][readWithMoreHits];

                                if (!noUkkonen) {
                                    scoreLimit = bestPairScore + extraSearchDepth;
                                }

                                isBestHit = true;
                            } else {
                                if (maxEditDistanceForSecondaryResults != -1 && (unsigned)maxEditDistanceForSecondaryResults >= pairScore - bestPairScore) {
                                    //
                                    // A secondary result to save.
                                    //
                                    if (*nSecondaryResults >= secondaryResultBufferSize) {
                                        WriteErrorMessage("IntersectingPairedEndAligner::align(): out of secondary result buffer.  Read ID %.*s\n", read0->getIdLength(), read0->getId());
                                        soft_exit(1);
                                    }

                                    PairedAlignmentResult *result = &secondaryResults[*nSecondaryResults];
                                    result->alignedAsPair = true;
                                    result->direction[readWithMoreHits] = setPairDirection[candidate->whichSetPair][readWithMoreHits];
                                    result->direction[readWithFewerHits] = setPairDirection[candidate->whichSetPair][readWithFewerHits];
                                    result->fromAlignTogether = true;
                                    result->location[readWithMoreHits] = mate->readWithMoreHitsGenomeLocation + mate->genomeOffset;
                                    result->location[readWithFewerHits] = candidate->readWithFewerHitsGenomeLocation + fewerEndGenomeLocationOffset;
                                    result->mapq[0] = result->mapq[1] = 0;
                                    result->score[readWithMoreHits] = mate->score;
                                    result->score[readWithFewerHits] = fewerEndScore;
                                    result->status[readWithFewerHits] = result->status[readWithMoreHits] = MultipleHits;

                                    (*nSecondaryResults)++;
                                }
                            }

                            probabilityOfAllPairs += pairProbability;
    #ifdef  _DEBUG
                            if (_DumpAlignments) {
                                printf("Added %e (= %e * %e) @ (%u, %u), giving new probability of all pairs %e, score %d = %d + %d%s\n",
                                    pairProbability, mate->matchProbability , fewerEndMatchProbability,
                                    candidate->readWithFewerHitsGenomeLocation + fewerEndGenomeLocationOffset, mate->readWithMoreHitsGenomeLocation + mate->genomeOffset,
                                    probabilityOfAllPairs,
                                    pairScore, fewerEndScore, mate->score, isBestHit ? " New best hit" : "");
                            }
    #endif  // _DEBUG

                            if (probabilityOfAllPairs >= 4.9 && -1 == maxEditDistanceForSecondaryResults) {
                                //
                                // Nothing will rescue us from a 0 MAPQ, so just stop looking.
                                //
                                goto doneScoring;
                            }
                        }
                    }// if the mate has a non -1 score
                }

                if (mateIndex == 0 || !genomeLocationIsWithin(scoringMateCandidates[candidate->whichSetPair][mateIndex-1].readWithMoreHitsGenomeLocation, candidate->readWithFewerHitsGenomeLocation, maxSpacing)) {
                    //
                    // Out of mate candidates.
                    //
                    break;
                }

                mateIndex--;
            }
        }

        //
        // Remove us from the head of the list and proceed to the next candidate to score.
        //
        scoringCandidates[currentBestPossibleScoreList] = candidate->scoreListNext;
     }

doneScoring:

    if (bestPairScore == 65536) {
        //
        // Found nothing.
        //
        for (unsigned whichRead = 0; whichRead < NUM_READS_PER_PAIR; whichRead++) {
            result->location[whichRead] = InvalidGenomeLocation;
            result->mapq[whichRead] = 0;
            result->score[whichRead] = -1;
            result->status[whichRead] = NotFound;
#ifdef  _DEBUG
            if (_DumpAlignments) {
                printf("No sufficiently good pairs found.\n");
            }
#endif  // DEBUG
        }
    } else {
        for (unsigned whichRead = 0; whichRead < NUM_READS_PER_PAIR; whichRead++) {
            result->location[whichRead] = bestResultGenomeLocation[whichRead];
            result->direction[whichRead] = bestResultDirection[whichRead];
            result->mapq[whichRead] = computeMAPQ(probabilityOfAllPairs, probabilityOfBestPair, bestResultScore[whichRead], popularSeedsSkipped[0] + popularSeedsSkipped[1]);
            result->status[whichRead] = result->mapq[whichRead] > MAPQ_LIMIT_FOR_SINGLE_HIT ? SingleHit : MultipleHits;
            result->score[whichRead] = bestResultScore[whichRead];
        }
#ifdef  _DEBUG
            if (_DumpAlignments) {
                printf("Returned %u %s %u %s with MAPQ %d and %d, probability of all pairs %e, probability of best pair %e\n",
                    result->location[0], result->direction[0] == RC ? "RC" : "", result->location[1], result->direction[1] == RC ? "RC" : "", result->mapq[0], result->mapq[1],
                    probabilityOfAllPairs, probabilityOfBestPair);
            }
#endif  // DEBUG
    }

    //
    // Get rid of any secondary results that are too far away from the best score.  (NB: the rest of the code in align() is very similar to BaseAligner::finalizeSecondaryResults.  Sorry)
    //
    int i = 0;
    while (i < *nSecondaryResults) {
        if ((int)(secondaryResults[i].score[0] + secondaryResults[i].score[1]) > (int)bestPairScore + maxEditDistanceForSecondaryResults) {
            secondaryResults[i] = secondaryResults[(*nSecondaryResults) - 1];
            (*nSecondaryResults)--;
        } else {
            i++;
        }
    }

    //
    // Now check to see if there are too many for any particular contig.
    //
    if (maxSecondaryAlignmentsPerContig > 0 && result->status[0] != NotFound) {
        //
        // Run through the results and count the number of results per contig, to see if any of them are too big.
        // First, record the primary result.
        //

        bool anyContigHasTooManyResults = false;
        contigCountEpoch++;

        int primaryContigNum = genome->getContigNumAtLocation(result->location[0]);
        hitsPerContigCounts[primaryContigNum].hits = 1;
        hitsPerContigCounts[primaryContigNum].epoch = contigCountEpoch;
        

        for (i = 0; i < *nSecondaryResults; i++) {
            int contigNum = genome->getContigNumAtLocation(secondaryResults[i].location[0]);    // We know they're on the same contig, so either will do
            if (hitsPerContigCounts[contigNum].epoch != contigCountEpoch) {
                hitsPerContigCounts[contigNum].epoch = contigCountEpoch;
                hitsPerContigCounts[contigNum].hits = 0;
            }

            hitsPerContigCounts[contigNum].hits++;
            if (hitsPerContigCounts[contigNum].hits > maxSecondaryAlignmentsPerContig) {
                anyContigHasTooManyResults = true;
                break;
            }
        }

        if (anyContigHasTooManyResults) {
            //
            // Just sort them all, in order of contig then hit depth.
            //
            qsort(secondaryResults, *nSecondaryResults, sizeof(*secondaryResults), PairedAlignmentResult::compareByContigAndScore);

            //
            // Now run through and eliminate any contigs with too many hits.  We can't use the same trick at the first loop above, because the
            // counting here relies on the results being sorted.  So, instead, we just copy them as we go.
            //
            int currentContigNum = -1;
            int currentContigCount = 0;
            int destResult = 0;

            for (int sourceResult = 0; sourceResult < *nSecondaryResults; sourceResult++) {
                int contigNum = genome->getContigNumAtLocation(secondaryResults[sourceResult].location[0]);
                if (contigNum != currentContigNum) {
                    currentContigNum = contigNum;
                    currentContigCount = (contigNum == primaryContigNum) ? 1 : 0;
                }

                currentContigCount++;

                if (currentContigCount <= maxSecondaryAlignmentsPerContig) {
                    //
                    // Keep it.  If we don't get here, then we don't copy the result and
                    // don't increment destResult.  And yes, this will sometimes copy a
                    // result over itself.  That's harmless.
                    //
                    secondaryResults[destResult] = secondaryResults[sourceResult];
                    destResult++;
                }
            } // for each source result
            *nSecondaryResults = destResult;
        }
    } // if we're limiting by contig


    if (*nSecondaryResults > maxSecondaryResultsToReturn) {
        qsort(secondaryResults, *nSecondaryResults, sizeof(*secondaryResults), PairedAlignmentResult::compareByScore);
        *nSecondaryResults = maxSecondaryResultsToReturn;   // Just truncate it
    }
}

    void
IntersectingPairedEndAligner::scoreLocation(
    unsigned             whichRead,
    Direction            direction,
    GenomeLocation       genomeLocation,
    unsigned             seedOffset,
    unsigned             scoreLimit,
    unsigned            *score,
    double              *matchProbability,
    int                 *genomeLocationOffset)
{
    nLocationsScored++;

    Read *readToScore = reads[whichRead][direction];
    unsigned readDataLength = readToScore->getDataLength();
    GenomeDistance genomeDataLength = readDataLength + MAX_K; // Leave extra space in case the read has deletions
    const char *data = genome->getSubstring(genomeLocation, genomeDataLength);

#if		0 // This only happens when genomeLocation is in the padding, which can lead to no good.  Just say no.
    if (NULL == data) {
        //
        // We're up against the end of a contig.  Reduce the extra space enough that it isn't too
        // long.  We're willing to reduce it to less than the length of a read, because the read could
        // butt up against the end of the contig and have insertions in it.
        //
        const Genome::Contig *contig = genome->getContigAtLocation(genomeLocation);

        GenomeLocation endLocation;
        if (genomeLocation + readDataLength + MAX_K >= genome->getCountOfBases()) {
            endLocation = genome->getCountOfBases();
        } else {
            const Genome::Contig *nextContig = genome->getContigAtLocation(genomeLocation + readDataLength + MAX_K);
            _ASSERT(NULL != contig && contig->beginningLocation <= genomeLocation && contig != nextContig);

            endLocation = nextContig->beginningLocation;
        }
        genomeDataLength = endLocation - genomeLocation - 1;
        if (genomeDataLength >= readDataLength - MAX_K) {
            data = genome->getSubstring(genomeLocation, genomeDataLength);
            _ASSERT(NULL != data);
        }
    }

#endif // 0 This only happens when genomeLocation is in the padding, which can lead to no good.  Just say no.

    if (NULL == data) {
        *score = -1;
        *matchProbability = 0;
        return;
    }


    // Compute the distance separately in the forward and backward directions from the seed, to allow
    // arbitrary offsets at both the start and end but not have to pay the cost of exploring all start
    // shifts in BoundedStringDistance
    double matchProb1, matchProb2;
    int score1, score2;
    // First, do the forward direction from where the seed aligns to past of it
    int readLen = readToScore->getDataLength();
    int seedLen = index->getSeedLength();
    int tailStart = seedOffset + seedLen;

    _ASSERT(!memcmp(data+seedOffset, readToScore->getData() + seedOffset, seedLen));    // that the seed actually matches

    int textLen;
    if (genomeDataLength - tailStart > INT32_MAX) {
        textLen = INT32_MAX;
    } else {
        textLen = (int)(genomeDataLength - tailStart);
    }
    score1 = landauVishkin->computeEditDistance(data + tailStart, textLen, readToScore->getData() + tailStart, readToScore->getQuality() + tailStart, readLen - tailStart,
        scoreLimit, &matchProb1);
    if (score1 == -1) {
        *score = -1;
    } else {
        // The tail of the read matched; now let's reverse the reference genome data and match the head
        int limitLeft = scoreLimit - score1;
        score2 = reverseLandauVishkin->computeEditDistance(data + seedOffset, seedOffset + MAX_K, reversedRead[whichRead][direction] + readLen - seedOffset,
                                                                    reads[whichRead][OppositeDirection(direction)]->getQuality() + readLen - seedOffset, seedOffset, limitLeft, &matchProb2, genomeLocationOffset);

        if (score2 == -1) {
            *score = -1;
        } else {
            *score = score1 + score2;
            _ASSERT(*score <= scoreLimit);
            // Map probabilities for substrings can be multiplied, but make sure to count seed too
            *matchProbability = matchProb1 * matchProb2 * pow(1 - SNP_PROB, seedLen);
        }
    }

    if (*score == -1) {
        *matchProbability = 0;
    }
}

    void
 IntersectingPairedEndAligner::HashTableHitSet::firstInit(unsigned maxSeeds_, unsigned maxMergeDistance_, BigAllocator *allocator, bool doesGenomeIndexHave64BitLocations_)
 {
    maxSeeds = maxSeeds_;
    maxMergeDistance = maxMergeDistance_;
    doesGenomeIndexHave64BitLocations = doesGenomeIndexHave64BitLocations_;
    nLookupsUsed = 0;
    if (doesGenomeIndexHave64BitLocations) {
        lookups64 = (HashTableLookup<GenomeLocation> *)allocator->allocate(sizeof(HashTableLookup<GenomeLocation>) * maxSeeds);
        lookups32 = NULL;
    } else {
        lookups32 = (HashTableLookup<unsigned> *)allocator->allocate(sizeof(HashTableLookup<unsigned>) * maxSeeds);
        lookups64 = NULL;
    }
    disjointHitSets = (DisjointHitSet *)allocator->allocate(sizeof(DisjointHitSet) * maxSeeds);
 }
    void
IntersectingPairedEndAligner::HashTableHitSet::init()
{
    nLookupsUsed = 0;
    currentDisjointHitSet = -1;
    if (doesGenomeIndexHave64BitLocations) {
        lookupListHead64->nextLookupWithRemainingMembers = lookupListHead64->prevLookupWithRemainingMembers = lookupListHead64;
        lookupListHead32->nextLookupWithRemainingMembers = lookupListHead32->prevLookupWithRemainingMembers = NULL;
    } else {
        lookupListHead32->nextLookupWithRemainingMembers = lookupListHead32->prevLookupWithRemainingMembers = lookupListHead32;
        lookupListHead64->nextLookupWithRemainingMembers = lookupListHead64->prevLookupWithRemainingMembers = NULL;
    }
}


//
// I apologize for this, but I had to do two versions of recordLookup, one for the 32 bit and one for the 64 bit version.  The options were
// copying the code or doing a macro with the types as parameters.  I chose macro, so you get ugly but unlikely to accidentally diverge.
// At least it's just isolated to the HashTableHitSet class.
//

#define RL(lookups, glType, lookupListHead)                                                                                                                 \
    void                                                                                                                                                    \
IntersectingPairedEndAligner::HashTableHitSet::recordLookup(unsigned seedOffset, _int64 nHits, const glType *hits, bool beginsDisjointHitSet)               \
{                                                                                                                                                           \
    _ASSERT(nLookupsUsed < maxSeeds);                                                                                                                       \
    if (beginsDisjointHitSet) {                                                                                                                             \
        currentDisjointHitSet++;                                                                                                                            \
        _ASSERT(currentDisjointHitSet < (int)maxSeeds);                                                                                                     \
        disjointHitSets[currentDisjointHitSet].countOfExhaustedHits = 0;                                                                                    \
    }                                                                                                                                                       \
                                                                                                                                                            \
    if (0 == nHits) {                                                                                                                                       \
        disjointHitSets[currentDisjointHitSet].countOfExhaustedHits++;                                                                                      \
    } else {                                                                                                                                                \
        _ASSERT(currentDisjointHitSet != -1);    /* Essentially that beginsDisjointHitSet is set for the first recordLookup call */                         \
        lookups[nLookupsUsed].currentHitForIntersection = 0;                                                                                                \
        lookups[nLookupsUsed].hits = hits;                                                                                                                  \
        lookups[nLookupsUsed].nHits = nHits;                                                                                                                \
        lookups[nLookupsUsed].seedOffset = seedOffset;                                                                                                      \
        lookups[nLookupsUsed].whichDisjointHitSet = currentDisjointHitSet;                                                                                  \
                                                                                                                                                            \
        /* Trim off any hits that are smaller than seedOffset, since they are clearly meaningless. */                                                       \
                                                                                                                                                            \
        while (lookups[nLookupsUsed].nHits > 0 && lookups[nLookupsUsed].hits[lookups[nLookupsUsed].nHits - 1] < lookups[nLookupsUsed].seedOffset) {         \
            lookups[nLookupsUsed].nHits--;                                                                                                                  \
        }                                                                                                                                                   \
                                                                                                                                                            \
        /* Add this lookup into the non-empty lookup list. */                                                                                               \
                                                                                                                                                            \
        lookups[nLookupsUsed].prevLookupWithRemainingMembers = lookupListHead;                                                                              \
        lookups[nLookupsUsed].nextLookupWithRemainingMembers = lookupListHead->nextLookupWithRemainingMembers;                                              \
        lookups[nLookupsUsed].prevLookupWithRemainingMembers->nextLookupWithRemainingMembers =                                                              \
            lookups[nLookupsUsed].nextLookupWithRemainingMembers->prevLookupWithRemainingMembers = &lookups[nLookupsUsed];                                  \
                                                                                                                                                            \
        if (doAlignerPrefetch) {                                                                                                                            \
            _mm_prefetch((const char *)&lookups[nLookupsUsed].hits[lookups[nLookupsUsed].nHits / 2], _MM_HINT_T2);                                          \
        }                                                                                                                                                   \
                                                                                                                                                            \
        nLookupsUsed++;                                                                                                                                     \
    }                                                                                                                                                       \
}

RL(lookups32, unsigned, lookupListHead32)
RL(lookups64, GenomeLocation, lookupListHead64)

#undef RL


	unsigned
IntersectingPairedEndAligner::HashTableHitSet::computeBestPossibleScoreForCurrentHit()
{
 	//
	// Now compute the best possible score for the hit.  This is the largest number of misses in any disjoint hit set.
	//
    for (int i = 0; i <= currentDisjointHitSet; i++) {
        disjointHitSets[i].missCount = disjointHitSets[i].countOfExhaustedHits;
    }

    //
    // Another macro.  Sorry again.
    //
#define loop(glType, lookupListHead)                                                                                                                                \
	for (HashTableLookup<glType> *lookup = lookupListHead->nextLookupWithRemainingMembers; lookup != lookupListHead;                                                \
         lookup = lookup->nextLookupWithRemainingMembers) {                                                                                                         \
                                                                                                                                                                    \
		if (!(lookup->currentHitForIntersection != lookup->nHits &&                                                                                                 \
				genomeLocationIsWithin(lookup->hits[lookup->currentHitForIntersection], mostRecentLocationReturned + lookup->seedOffset,  maxMergeDistance) ||      \
			lookup->currentHitForIntersection != 0 &&                                                                                                               \
				genomeLocationIsWithin(lookup->hits[lookup->currentHitForIntersection-1], mostRecentLocationReturned + lookup->seedOffset,  maxMergeDistance))) {   \
                                                                                                                                                                    \
			/* This one was not close enough. */                                                                                                                    \
                                                                                                                                                                    \
			disjointHitSets[lookup->whichDisjointHitSet].missCount++;                                                                                               \
		}                                                                                                                                                           \
	}

    if (doesGenomeIndexHave64BitLocations) {
        loop(GenomeLocation, lookupListHead64);
    } else {
        loop(unsigned, lookupListHead32);
    }
#undef loop

    unsigned bestPossibleScoreSoFar = 0;
    for (int i = 0; i <= currentDisjointHitSet; i++) {
        bestPossibleScoreSoFar = max(bestPossibleScoreSoFar, disjointHitSets[i].missCount);
    }

	return bestPossibleScoreSoFar;
}

	bool
IntersectingPairedEndAligner::HashTableHitSet::getNextHitLessThanOrEqualTo(GenomeLocation maxGenomeLocationToFind, GenomeLocation *actualGenomeLocationFound, unsigned *seedOffsetFound)
{

    bool anyFound = false;
    GenomeLocation bestLocationFound = 0;
    for (unsigned i = 0; i < nLookupsUsed; i++) {
        //
        // Binary search from the current starting offset to either the right place or the end.
        //
        _int64 limit[2];
        GenomeLocation maxGenomeLocationToFindThisSeed;

        if (doesGenomeIndexHave64BitLocations) {
            limit[0] = (_int64)lookups64[i].currentHitForIntersection;
            limit[1] = (_int64)lookups64[i].nHits - 1;
            maxGenomeLocationToFindThisSeed = maxGenomeLocationToFind + lookups64[i].seedOffset;
        } else {
            limit[0] = (_int64)lookups32[i].currentHitForIntersection;
            limit[1] = (_int64)lookups32[i].nHits - 1;
            maxGenomeLocationToFindThisSeed = maxGenomeLocationToFind + lookups32[i].seedOffset;
        }
 
        while (limit[0] <= limit[1]) {
            _int64 probe = (limit[0] + limit[1]) / 2;
            if (doAlignerPrefetch) { // not clear this helps.  We're probably not far enough ahead.
                if (doesGenomeIndexHave64BitLocations) {
                    _mm_prefetch((const char *)&lookups64[i].hits[(limit[0] + probe) / 2 - 1], _MM_HINT_T2);
                    _mm_prefetch((const char *)&lookups64[i].hits[(limit[1] + probe) / 2 + 1], _MM_HINT_T2);
                } else {
                    _mm_prefetch((const char *)&lookups32[i].hits[(limit[0] + probe) / 2 - 1], _MM_HINT_T2);
                    _mm_prefetch((const char *)&lookups32[i].hits[(limit[1] + probe) / 2 + 1], _MM_HINT_T2);
                }
            }
            //
            // Recall that the hit sets are sorted from largest to smallest, so the strange looking logic is actually right.
            // We're evaluating the expression "lookups[i].hits[probe] <= maxGenomeOffsetToFindThisSeed && (probe == 0 || lookups[i].hits[probe-1] > maxGenomeOffsetToFindThisSeed)"
            // It's written in this strange way just so the profile tool will show us where the time's going.
            //
            GenomeLocation probeHit;
            GenomeLocation probeMinusOneHit;
            unsigned seedOffset;
            if (doesGenomeIndexHave64BitLocations) {
                probeHit = lookups64[i].hits[probe];
                probeMinusOneHit = lookups64[i].hits[probe-1];
                seedOffset = lookups64[i].seedOffset;
            } else {
                probeHit = lookups32[i].hits[probe];
                probeMinusOneHit = lookups32[i].hits[probe-1];
                seedOffset = lookups32[i].seedOffset;
            }
            unsigned clause1 =  probeHit <= maxGenomeLocationToFindThisSeed;
            unsigned clause2 = probe == 0;

            if (clause1 && (clause2 || probeMinusOneHit > maxGenomeLocationToFindThisSeed)) {
                if (probeHit - seedOffset > bestLocationFound) {
					anyFound = true;
                    mostRecentLocationReturned = *actualGenomeLocationFound = bestLocationFound = probeHit - seedOffset;
                    *seedOffsetFound = seedOffset;
                }

                if (doesGenomeIndexHave64BitLocations) {
                    lookups64[i].currentHitForIntersection = probe;
                } else {
                    lookups32[i].currentHitForIntersection = probe;
                }
                break;
            }

            if (probeHit > maxGenomeLocationToFindThisSeed) {   // Recode this without the if to avoid the hard-to-predict branch.
                limit[0] = probe + 1;
            } else {
                limit[1] = probe - 1;
            }
        } // While we're looking

        if (limit[0] > limit[1]) {
            // We're done with this lookup.
            if (doesGenomeIndexHave64BitLocations) {
                lookups64[i].currentHitForIntersection = lookups64[i].nHits;    
            } else {
                lookups32[i].currentHitForIntersection = lookups32[i].nHits;
            }
        }
    } // For each lookup

    _ASSERT(!anyFound || *actualGenomeLocationFound <= maxGenomeLocationToFind);

    return anyFound;
}


    bool
IntersectingPairedEndAligner::HashTableHitSet::getFirstHit(GenomeLocation *genomeLocation, unsigned *seedOffsetFound)
{
    bool anyFound = false;
    *genomeLocation = 0;

    //
    // Yet another macro.  This makes me want to write in a better language sometimes.  But then it would be too slow.  :-(
    //

#define LOOP(lookups)                                                                                                                       \
    for (unsigned i = 0; i < nLookupsUsed; i++) {                                                                                           \
        if (lookups[i].nHits > 0 && lookups[i].hits[0] - lookups[i].seedOffset > GenomeLocationAsInt64(*genomeLocation)) {                  \
            mostRecentLocationReturned = *genomeLocation = lookups[i].hits[0] - lookups[i].seedOffset;                                      \
            *seedOffsetFound = lookups[i].seedOffset;                                                                                       \
            anyFound = true;                                                                                                                \
        }                                                                                                                                   \
    }

    if (doesGenomeIndexHave64BitLocations) {
        LOOP(lookups64);
    } else {
        LOOP(lookups32);
    }

#undef LOOP

	return !anyFound;
}

    bool
IntersectingPairedEndAligner::HashTableHitSet::getNextLowerHit(GenomeLocation *genomeLocation, unsigned *seedOffsetFound)
{
    //
    // Look through all of the lookups and find the one with the highest location smaller than the current one.
    //
    GenomeLocation foundLocation = 0;
    bool anyFound = false;

    //
    // Run through the lookups pushing up any that are at the most recently returned
    //

    for (unsigned i = 0; i < nLookupsUsed; i++) {
        _int64 *currentHitForIntersection;
        _int64 nHits;
        GenomeLocation hitLocation;
        unsigned seedOffset;

        //
        // A macro to initialize stuff that we need to avoid a bigger macro later.
        //
#define initVars(lookups)                                                                                               \
        currentHitForIntersection = &lookups[i].currentHitForIntersection;                                              \
        nHits = lookups[i].nHits;                                                                                       \
        seedOffset = lookups[i].seedOffset;                                                                             \
        if (nHits != *currentHitForIntersection) {                                                                      \
            hitLocation = lookups[i].hits[*currentHitForIntersection];                                                  \
        }


        if (doesGenomeIndexHave64BitLocations) {
            initVars(lookups64);
        } else {
            initVars(lookups32);
        }
#undef  initVars

        _ASSERT(*currentHitForIntersection == nHits || hitLocation - seedOffset <= mostRecentLocationReturned || hitLocation < seedOffset);

        if (*currentHitForIntersection != nHits && hitLocation - seedOffset == mostRecentLocationReturned) {
            (*currentHitForIntersection)++;
            if (*currentHitForIntersection == nHits) {
                continue;
            }
            if (doesGenomeIndexHave64BitLocations) {
                hitLocation = lookups64[i].hits[*currentHitForIntersection];
            } else {
                hitLocation = lookups32[i].hits[*currentHitForIntersection];
            }
        }

        if (*currentHitForIntersection != nHits) {
            if (foundLocation < hitLocation - seedOffset && // found location is OK
                hitLocation >= seedOffset) // found location isn't too small to push us before the beginning of the genome
            {
                *genomeLocation = foundLocation = hitLocation - seedOffset;
                *seedOffsetFound = seedOffset;
                anyFound = true;
            }
        }
    }

    if (anyFound) {
        mostRecentLocationReturned = foundLocation;
    }

    return anyFound;
}

            bool
IntersectingPairedEndAligner::MergeAnchor::checkMerge(GenomeLocation newMoreHitLocation, GenomeLocation newFewerHitLocation, double newMatchProbability, int newPairScore,
                        double *oldMatchProbability)
{
    if (locationForReadWithMoreHits == InvalidGenomeLocation || !doesRangeMatch(newMoreHitLocation, newFewerHitLocation)) {
        //
        // No merge.  Remember the new one.
        //
        locationForReadWithMoreHits = newMoreHitLocation;
        locationForReadWithFewerHits = newFewerHitLocation;
        matchProbability = newMatchProbability;
        pairScore = newPairScore;
        *oldMatchProbability = 0.0;
        return false;
    }  else {
        //
        // Within merge distance.  Keep the better score (or if they're tied the better match probability).
        //
        if (newPairScore < pairScore || newPairScore == pairScore && newMatchProbability > matchProbability) {
#ifdef _DEBUG
            if (_DumpAlignments) {
                printf("Merge replacement at anchor (%u, %u), loc (%u, %u), old match prob %e, new match prob %e, old pair score %d, new pair score %d\n",
                    locationForReadWithMoreHits, locationForReadWithFewerHits, newMoreHitLocation, newFewerHitLocation,
                    matchProbability, newMatchProbability, pairScore, newPairScore);
            }
#endif // DEBUG

            *oldMatchProbability = matchProbability;
            matchProbability = newMatchProbability;
            pairScore = newPairScore;
            return false;
        } else {
            //
            // The new one should just be ignored.
            //
#ifdef _DEBUG
            if (_DumpAlignments) {
                printf("Merged at anchor (%u, %u), loc (%u, %u), old match prob %e, new match prob %e, old pair score %d, new pair score %d\n",
                    locationForReadWithMoreHits, locationForReadWithFewerHits, newMoreHitLocation, newFewerHitLocation,
                    matchProbability, newMatchProbability, pairScore, newPairScore);
            }
#endif // DEBUG
            return true;
        }
    }

    _ASSERT(!"NOTREACHED");
}

const unsigned IntersectingPairedEndAligner::maxMergeDistance = 31;
