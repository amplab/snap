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
        unsigned      maxKForIndels_,
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
		bool          noTruncation_,
        bool          useAffineGap_,    
        bool          ignoreAlignmentAdjustmentsForOm_,
		bool		  altAwareness_,
        unsigned      maxScoreGapToPreferNonAltAlignment_,
        unsigned      matchReward_,
        unsigned      subPenalty_,
        unsigned      gapOpenPenalty_,
        unsigned      gapExtendPenalty_) :
    index(index_), maxReadSize(maxReadSize_), maxHits(maxHits_), maxK(maxK_), maxKForIndels(maxKForIndels_), numSeedsFromCommandLine(__min(MAX_MAX_SEEDS,numSeedsFromCommandLine_)), minSpacing(minSpacing_), maxSpacing(maxSpacing_),
	landauVishkin(NULL), reverseLandauVishkin(NULL), maxBigHits(maxBigHits_), seedCoverage(seedCoverage_),
    extraSearchDepth(extraSearchDepth_), nLocationsScored(0), noUkkonen(noUkkonen_), noOrderedEvaluation(noOrderedEvaluation_), noTruncation(noTruncation_), useAffineGap(useAffineGap_),
    maxSecondaryAlignmentsPerContig(maxSecondaryAlignmentsPerContig_), alignmentAdjuster(index->getGenome()), ignoreAlignmentAdjustmentsForOm(ignoreAlignmentAdjustmentsForOm_), altAwareness(altAwareness_),
    maxScoreGapToPreferNonAltAlignment(maxScoreGapToPreferNonAltAlignment_), matchReward(matchReward_), subPenalty(subPenalty_), gapOpenPenalty(gapOpenPenalty_), gapExtendPenalty(gapExtendPenalty_)
{
    doesGenomeIndexHave64BitLocations = index->doesGenomeIndexHave64BitLocations();

    unsigned maxSeedsToUse;
    if (0 != numSeedsFromCommandLine) {
        maxSeedsToUse = numSeedsFromCommandLine;
    } else {
        maxSeedsToUse = (unsigned)(maxReadSize * seedCoverage / index->getSeedLength());
    }
    allocateDynamicMemory(allocator, maxReadSize, maxBigHits, maxSeedsToUse, MAX_K, extraSearchDepth, maxCandidatePoolSize, maxSecondaryAlignmentsPerContig);

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
    seedUsed = (BYTE *) allocator->allocate(100 + ((size_t)maxReadSize + 7) / 8);

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

    scoringCandidates = (ScoringCandidate **) allocator->allocate(sizeof(ScoringCandidate *) * ((size_t)maxEditDistanceToConsider + maxExtraSearchDepth + 1));  //+1 is for 0.
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

    bool
IntersectingPairedEndAligner::align(
        Read                  *read0,
        Read                  *read1,
        PairedAlignmentResult* result,
        PairedAlignmentResult* firstALTResult,
        int                    maxEditDistanceForSecondaryResults,
        _int64                 secondaryResultBufferSize,
        _int64                *nSecondaryResults,
        PairedAlignmentResult *secondaryResults,             // The caller passes in a buffer of secondaryResultBufferSize and it's filled in by align()
        _int64                 singleSecondaryBufferSize,
        _int64                 maxSecondaryResultsToReturn,
        _int64                *nSingleEndSecondaryResultsForFirstRead,
        _int64                *nSingleEndSecondaryResultsForSecondRead,
        SingleAlignmentResult *singleEndSecondaryResults,     // Single-end secondary alignments for when the paired-end alignment didn't work properly
        _int64				  maxLVCandidatesForAffineGapBufferSize,
        _int64				  *nLVCandidatesForAffineGap,
        PairedAlignmentResult *lvCandidatesForAffineGap, // Landau-Vishkin candidates that need to be rescored using affine gap
        int                   maxK_
	)
{
    maxK = maxK_;
    if (!useAffineGap) {
        //
        // Version with no affine gap scoring
        //
        return alignLandauVishkin(read0, read1, result, firstALTResult, maxEditDistanceForSecondaryResults, secondaryResultBufferSize,
	        nSecondaryResults, secondaryResults, singleSecondaryBufferSize, maxSecondaryResultsToReturn, nSingleEndSecondaryResultsForFirstRead, nSingleEndSecondaryResultsForSecondRead,
	        singleEndSecondaryResults, maxLVCandidatesForAffineGapBufferSize, nLVCandidatesForAffineGap, lvCandidatesForAffineGap);
    }
    else {
        //
        // Perform seeding, set intersection, LV alignment and identify promising candidates for affine gap scoring
        //
        bool fitInSecondaryBuffer = alignLandauVishkin(read0, read1, result, firstALTResult, maxEditDistanceForSecondaryResults, secondaryResultBufferSize,
	        nSecondaryResults, secondaryResults, singleSecondaryBufferSize, maxSecondaryResultsToReturn, nSingleEndSecondaryResultsForFirstRead, nSingleEndSecondaryResultsForSecondRead,
	        singleEndSecondaryResults, maxLVCandidatesForAffineGapBufferSize, nLVCandidatesForAffineGap, lvCandidatesForAffineGap);

        if (*nLVCandidatesForAffineGap > maxLVCandidatesForAffineGapBufferSize) {
	        *nLVCandidatesForAffineGap = maxLVCandidatesForAffineGapBufferSize + 1;
	        return false;
        }

        if (!fitInSecondaryBuffer) {
            return false;
        }

        //
        // Perform affine gap scoring for promising candidates to get best scoring hit
        //
        fitInSecondaryBuffer = alignAffineGap(read0, read1, result, firstALTResult, maxEditDistanceForSecondaryResults, secondaryResultBufferSize,
	        nSecondaryResults, secondaryResults, singleSecondaryBufferSize, maxSecondaryResultsToReturn, nSingleEndSecondaryResultsForFirstRead, nSingleEndSecondaryResultsForSecondRead,
	        singleEndSecondaryResults, maxLVCandidatesForAffineGapBufferSize, nLVCandidatesForAffineGap, lvCandidatesForAffineGap);

        if (!fitInSecondaryBuffer) {
            return false;
        }

        return true;
    }
}

bool
    IntersectingPairedEndAligner::alignLandauVishkin(
        Read* read0,
        Read* read1,
        PairedAlignmentResult* result,
        PairedAlignmentResult* firstALTResult,
        int                    maxEditDistanceForSecondaryResults,
        _int64                 secondaryResultBufferSize,
        _int64* nSecondaryResults,
        PairedAlignmentResult* secondaryResults,             // The caller passes in a buffer of secondaryResultBufferSize and it's filled in by align()
        _int64                 singleSecondaryBufferSize,
        _int64                 maxSecondaryResultsToReturn,
        _int64* nSingleEndSecondaryResultsForFirstRead,
        _int64* nSingleEndSecondaryResultsForSecondRead,
        SingleAlignmentResult* singleEndSecondaryResults,     // Single-end secondary alignments for when the paired-end alignment didn't work properly
        _int64                 maxLVCandidatesForAffineGapBufferSize,
        _int64* nLVCandidatesForAffineGap,
        PairedAlignmentResult* lvCandidatesForAffineGap
)
{
#if INSTRUMENTATION_FOR_PAPER
    _int64 startTime = timeInNanos();
    _int64 nScored = 0;
    _int64 setIntersectionSize = 0;
    bool bestCandidateScoredFirst = false;
#endif // INSTRUMENTATION_FOR_PAPER

    firstALTResult->status[0] = firstALTResult->status[1] = NotFound;

    result->nLVCalls = 0;
    result->nSmallHits = 0;
    result->clippingForReadAdjustment[0] = result->clippingForReadAdjustment[1] = 0;
    result->usedAffineGapScoring[0] = result->usedAffineGapScoring[1] = false;
    result->basesClippedBefore[0] = result->basesClippedBefore[1] = 0;
    result->basesClippedAfter[0] = result->basesClippedAfter[1] = 0;
    result->agScore[0] = result->agScore[1] = 0;

    *nSecondaryResults = 0;
    *nSingleEndSecondaryResultsForFirstRead = 0;
    *nSingleEndSecondaryResultsForSecondRead = 0;
    *nLVCandidatesForAffineGap = 0;

    int maxSeeds;
    if (numSeedsFromCommandLine != 0) {
        maxSeeds = (int)numSeedsFromCommandLine;
    }
    else {
        maxSeeds = (int)(max(read0->getDataLength(), read1->getDataLength()) * seedCoverage / index->getSeedLength());
    }

#ifdef  _DEBUG
    if (_DumpAlignments) {
        printf("\nIntersectingAligner aligning reads '%*.s' and '%.*s' with data '%.*s' and '%.*s'\n", read0->getIdLength(), read0->getId(), read1->getIdLength(), read1->getId(), read0->getDataLength(), read0->getData(), read1->getDataLength(), read1->getData());
    }
#endif  // _DEBUG

    lowestFreeScoringCandidatePoolEntry = 0;
    for (int k = 0; k <= maxK + extraSearchDepth; k++) {
        scoringCandidates[k] = NULL;
    }

    for (unsigned i = 0; i < NUM_SET_PAIRS; i++) {
        lowestFreeScoringMateCandidate[i] = 0;
    }
    firstFreeMergeAnchor = 0;

    Read rcReads[NUM_READS_PER_PAIR];

    ScoreSet scoresForAllAlignments;
    ScoreSet scoresForNonAltAlignments;

    unsigned popularSeedsSkipped[NUM_READS_PER_PAIR];

    reads[0][FORWARD] = read0;
    reads[1][FORWARD] = read1;

    //
    // Don't bother if one or both reads are too short.  The minimum read length here is the seed length, but usually there's a longer
    // minimum enforced by our caller
    //
    if (read0->getDataLength() < seedLen || read1->getDataLength() < seedLen) {
        return true;
    }

    //
    // Build the RC reads.
    //
    unsigned countOfNs = 0;

    for (unsigned whichRead = 0; whichRead < NUM_READS_PER_PAIR; whichRead++) {
        Read* read = reads[whichRead][FORWARD];
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

    if ((int)countOfNs > maxK) {
        return true;
    }

    //
    // Build the reverse data for both reads in both directions for the backwards LV to use.
    //
    for (unsigned whichRead = 0; whichRead < NUM_READS_PER_PAIR; whichRead++) {
        for (Direction dir = 0; dir < NUM_DIRECTIONS; dir++) {
            Read* read = reads[whichRead][dir];

            for (unsigned i = 0; i < read->getDataLength(); i++) {
                reversedRead[whichRead][dir][i] = read->getData()[read->getDataLength() - i - 1];
            }
        }
    }

    unsigned thisPassSeedsNotSkipped[NUM_READS_PER_PAIR][NUM_DIRECTIONS] = { {0,0}, {0,0} };

    //
    // Initialize the member variables that are effectively stack locals, but are in the object
    // to avoid having to pass them to score.
    //

    localBestPairProbability[0] = 0;
    localBestPairProbability[1] = 0;

    //
    // Phase 1: do the hash table lookups for each of the seeds for each of the reads and add them to the hit sets.
    //

    for (unsigned whichRead = 0; whichRead < NUM_READS_PER_PAIR; whichRead++) {
        int nextSeedToTest = 0;
        unsigned wrapCount = 0;
        int nPossibleSeeds = (int)readLen[whichRead] - seedLen + 1;
        memset(seedUsed, 0, (__max(readLen[0], readLen[1]) + 7) / 8);
        bool beginsDisjointHitSet[NUM_DIRECTIONS] = { true, true };

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
            const GenomeLocation* hits[NUM_DIRECTIONS];
            const unsigned* hits32[NUM_DIRECTIONS];

            if (doesGenomeIndexHave64BitLocations) {
                index->lookupSeed(seed, &nHits[FORWARD], &hits[FORWARD], &nHits[RC], &hits[RC],
                    hashTableHitSets[whichRead][FORWARD]->getNextSingletonLocation(), hashTableHitSets[whichRead][RC]->getNextSingletonLocation());
            }
            else {
                index->lookupSeed32(seed, &nHits[FORWARD], &hits32[FORWARD], &nHits[RC], &hits32[RC]);
            }

            countOfHashTableLookups[whichRead]++;
            for (Direction dir = FORWARD; dir < NUM_DIRECTIONS; dir++) {
                int offset;
                if (dir == FORWARD) {
                    offset = nextSeedToTest;
                }
                else {
                    offset = readLen[whichRead] - seedLen - nextSeedToTest;
                }

                if (nHits[dir] < maxBigHits) {
                    totalHashTableHits[whichRead][dir] += nHits[dir];
                    if (doesGenomeIndexHave64BitLocations) {
                        hashTableHitSets[whichRead][dir]->recordLookup(offset, nHits[dir], hits[dir], beginsDisjointHitSet[dir]);
                    }
                    else {
                        hashTableHitSets[whichRead][dir]->recordLookup(offset, nHits[dir], hits32[dir], beginsDisjointHitSet[dir]);
                    }
                    beginsDisjointHitSet[dir] = false;
                }
                else {
                    popularSeedsSkipped[whichRead]++;
                }
            } // for each direction

            //
            // If we don't have enough seeds left to reach the end of the read, space out the seeds more-or-less evenly.
            //
            if ((maxSeeds - countOfHashTableLookups[whichRead] + 1) * (int)seedLen + nextSeedToTest < nPossibleSeeds) {
                _ASSERT((nPossibleSeeds - nextSeedToTest - 1) / (maxSeeds - countOfHashTableLookups[whichRead] + 1) >= (int)seedLen);
                nextSeedToTest += (nPossibleSeeds - nextSeedToTest - 1) / (maxSeeds - countOfHashTableLookups[whichRead] + 1);
                _ASSERT(nextSeedToTest < nPossibleSeeds);   // We haven't run off the end of the read.
            }
            else {
                nextSeedToTest += seedLen;
            }
        } // while we need to lookup seeds for this read
    } // for each read

#if INSTRUMENTATION_FOR_PAPER
    int hashTableHits[NUM_READS_PER_PAIR] = { hashTableHitSets[0][FORWARD]->getNumDistinctHitLocations(0) + hashTableHitSets[0][RC]->getNumDistinctHitLocations(0),
                                                hashTableHitSets[1][FORWARD]->getNumDistinctHitLocations(0) + hashTableHitSets[1][RC]->getNumDistinctHitLocations(0) };

    int log2HashTableHits[NUM_READS_PER_PAIR] = { __min(cheezyLogBase2(hashTableHits[0]), MAX_HIT_SIZE_LOG_2),__min(cheezyLogBase2(hashTableHits[1]), MAX_HIT_SIZE_LOG_2) };

    // = { __min(cheezyLogBase2(totalHashTableHits[0][FORWARD] + totalHashTableHits[0][RC]), MAX_HIT_SIZE_LOG_2), __min(cheezyLogBase2(totalHashTableHits[1][FORWARD] + totalHashTableHits[1][RC]), MAX_HIT_SIZE_LOG_2) };

#endif // INSTRUMENTATION_FOR_PAPER

    readWithMoreHits = totalHashTableHits[0][FORWARD] + totalHashTableHits[0][RC] > totalHashTableHits[1][FORWARD] + totalHashTableHits[1][RC] ? 0 : 1;
    readWithFewerHits = 1 - readWithMoreHits;

#ifdef  _DEBUG
    if (_DumpAlignments) {
        printf("Read 0 has %lld hits, read 1 has %lld hits\n", totalHashTableHits[0][FORWARD] + totalHashTableHits[0][RC], totalHashTableHits[1][FORWARD] + totalHashTableHits[1][RC]);
    }
#endif  // _DEBUG

    Direction setPairDirection[NUM_SET_PAIRS][NUM_READS_PER_PAIR] = {{FORWARD, RC}, {RC, FORWARD}};

    //
    // Phase 2: find all possible candidates and add them to candidate lists (for the reads with fewer and more hits).
    //
    int maxUsedBestPossibleScoreList = 0;

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
            // We go once through this loop for each.  Because the index is ordered high to low, we also go in that direction.
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
            // Add all of the mate candidates for the fewer side hit.
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
                    printf("SetPair %d, added more hits candidate %d at genome location %s:%llu, bestPossibleScore %d, seedOffset %d\n",
                            whichSetPair, lowestFreeScoringMateCandidate[whichSetPair], 
                            genome->getContigAtLocation(lastGenomeLocationForReadWithMoreHits)->name,
                            lastGenomeLocationForReadWithMoreHits - genome->getContigAtLocation(lastGenomeLocationForReadWithMoreHits)->beginningLocation,
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
			int bestPossibleScoreForReadWithFewerHits;
			
			if (noTruncation) {
				bestPossibleScoreForReadWithFewerHits = 0;
			} else {
				bestPossibleScoreForReadWithFewerHits = setPair[readWithFewerHits]->computeBestPossibleScoreForCurrentHit();
			}

            int lowestBestPossibleScoreOfAnyPossibleMate = maxK + extraSearchDepth;
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
                int bestPossibleScore = noOrderedEvaluation ? 0 : lowestBestPossibleScoreOfAnyPossibleMate + bestPossibleScoreForReadWithFewerHits;

                scoringCandidatePool[lowestFreeScoringCandidatePoolEntry].init(lastGenomeLocationForReadWithFewerHits, whichSetPair, lowestFreeScoringMateCandidate[whichSetPair] - 1,
                                                                                lastSeedOffsetForReadWithFewerHits, bestPossibleScoreForReadWithFewerHits,
                                                                                scoringCandidates[bestPossibleScore]);


                scoringCandidates[bestPossibleScore] = &scoringCandidatePool[lowestFreeScoringCandidatePoolEntry];

#if     INSTRUMENTATION_FOR_PAPER
                setIntersectionSize++;
#endif // INSTRUMENTATION_FOR_PAPER

#ifdef _DEBUG
                if (_DumpAlignments) {
                    printf("SetPair %d, added fewer hits candidate %d at genome location %s:%llu, bestPossibleScore %d, seedOffset %d\n",
                            whichSetPair, lowestFreeScoringCandidatePoolEntry, 
                            genome->getContigAtLocation(lastGenomeLocationForReadWithFewerHits)->name, lastGenomeLocationForReadWithFewerHits - genome->getContigAtLocation(lastGenomeLocationForReadWithFewerHits)->beginningLocation,
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
        } // forever (the loop that does the intersection walk)
    } // For each set pair

    //
    // Phase 2a: mark any scoring candidates that are near enough to one another to be indels so that we can increase the score limit
    // when considering them.
    //
    for (int whichSetPair = 0; whichSetPair < NUM_SET_PAIRS; whichSetPair++) {
        //
        // For every candidate, we're trying to find the candidate that's the farthest away from it in genome
        // space that's still within maxDistForIndels and then mark that distance in the candidate.  We do that
        // by keeping two indices into the array, a bottom and top (relative to the array; the actual genome
        // locations go the other way).  Each pass through the loop we move one or the other.  We move up top if
        // bottom is one less or they're within maxDistForIndels of one another.  Otherwise we move up bottom.
        //
        int bottom = 0;
        int top = 1;
        while (top < lowestFreeScoringMateCandidate[whichSetPair]) {
            _ASSERT(bottom < top);
            _ASSERT(scoringMateCandidates[whichSetPair][bottom].readWithMoreHitsGenomeLocation > scoringMateCandidates[whichSetPair][top].readWithMoreHitsGenomeLocation);  // Remember, they're in backward order.
            GenomeDistance spread = DistanceBetweenGenomeLocations(scoringMateCandidates[whichSetPair][bottom].readWithMoreHitsGenomeLocation, scoringMateCandidates[whichSetPair][top].readWithMoreHitsGenomeLocation);
            if (spread < maxKForIndels) {
                scoringMateCandidates[whichSetPair][bottom].largestBigIndelDetected = __max(spread, scoringMateCandidates[whichSetPair][bottom].largestBigIndelDetected);
                scoringMateCandidates[whichSetPair][top].largestBigIndelDetected = __max(spread, scoringMateCandidates[whichSetPair][top].largestBigIndelDetected);
#if _DEBUG
                if (_DumpAlignments) {
                    fprintf(stderr,"Set largest big indel detected to %d, %d for set pair %d mate candidates %d and %d\n", 
                            scoringMateCandidates[whichSetPair][bottom].largestBigIndelDetected, scoringMateCandidates[whichSetPair][top].largestBigIndelDetected,
                            whichSetPair, bottom, top);
                }
#endif // DEBUG
                top++;
            } else if (bottom < top - 1) {
                //
                // Move up bottom since it will still be less than top.
                //
                bottom++;
            } else {
                //
                // They're next to each other and still too far apart.  Move up both.
                //
                bottom++;
                top++;
            }
        } // While we're still considering candidates in one set pair.
    } // for each set pair

    //
    // Now do the fewer end candidates.  This works the same way as the loop above (except they're not segregated by set pair).
    //
    { // a scope for us to declare locals in.
        int bottom = 0;
        int top = 1;

        while (top < lowestFreeScoringCandidatePoolEntry) {
            _ASSERT(bottom < top);

            if (scoringCandidatePool[bottom].whichSetPair != scoringCandidatePool[top].whichSetPair) {
                bottom = top;
                top = top + 1;
                continue;
            }

            _ASSERT(scoringCandidatePool[bottom].readWithFewerHitsGenomeLocation > scoringCandidatePool[top].readWithFewerHitsGenomeLocation); // They're in backward order

            GenomeDistance spread = DistanceBetweenGenomeLocations(scoringCandidatePool[bottom].readWithFewerHitsGenomeLocation, scoringCandidatePool[top].readWithFewerHitsGenomeLocation);
            if (spread < maxKForIndels) {
                scoringCandidatePool[bottom].largestBigIndelDetected = __max(spread, scoringCandidatePool[bottom].largestBigIndelDetected);
                scoringCandidatePool[top].largestBigIndelDetected = __max(spread, scoringCandidatePool[top].largestBigIndelDetected);
#if _DEBUG
                if (_DumpAlignments) {
                    fprintf(stderr, "Set largest big indel to %d, %d for fewer end candidates %d and %d\n",
                            scoringCandidatePool[bottom].largestBigIndelDetected, scoringCandidatePool[top].largestBigIndelDetected, bottom, top);
                }
#endif // DEBUG
                top++;
            } else if (bottom < top - 1) {
                bottom++;
            } else {
                bottom++;
                top++;
            }
        } // While we're still looking.
    } // Scope for fewer candidate possible indel detection.  

    //
    // Phase 3: score and merge the candidates we've found using Laundau-Vishkin (edit distance, not affine gap).
    //
    int currentBestPossibleScoreList = 0;

    //
    // Loop until we've scored all of the candidates, or proven that what's left must have too high of a score to be interesting.
    // 
    //
    while (currentBestPossibleScoreList <= maxUsedBestPossibleScoreList && 
            currentBestPossibleScoreList <= extraSearchDepth + min(maxK, max(   // Never look for worse than our worst interesting score
                                                                           min(scoresForAllAlignments.bestPairScore, scoresForNonAltAlignments.bestPairScore - maxScoreGapToPreferNonAltAlignment),   // Worst we care about for ALT
                                                                           min(scoresForAllAlignments.bestPairScore + maxScoreGapToPreferNonAltAlignment, scoresForNonAltAlignments.bestPairScore)))) // And for non-ALT
    {
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

        int fewerEndScore;
        double fewerEndMatchProbability;
        int fewerEndGenomeLocationOffset;

        bool nonALTAlignment = (!altAwareness) || !genome->isGenomeLocationALT(candidate->readWithFewerHitsGenomeLocation);

        int scoreLimit = computeScoreLimit(nonALTAlignment, &scoresForAllAlignments, &scoresForNonAltAlignments, candidate->largestBigIndelDetected);

        if (currentBestPossibleScoreList > scoreLimit) {
            //
            // Remove us from the head of the list and proceed to the next candidate to score.  We can get here because now we know ALT/non-ALT, which have different limits.
            //
            scoringCandidates[currentBestPossibleScoreList] = candidate->scoreListNext;
            continue;
        }

        scoreLocation(readWithFewerHits, setPairDirection[candidate->whichSetPair][readWithFewerHits], candidate->readWithFewerHitsGenomeLocation,
            candidate->seedOffset, scoreLimit, &fewerEndScore, &fewerEndMatchProbability, &fewerEndGenomeLocationOffset, &candidate->usedAffineGapScoring,
            &candidate->basesClippedBefore, &candidate->basesClippedAfter, &candidate->agScore, &candidate->lvIndels);

#if INSTRUMENTATION_FOR_PAPER
        nScored++;
#endif // INSTRUMENTATION_FOR_PAPER

        candidate->matchProbability = fewerEndMatchProbability;

        _ASSERT(ScoreAboveLimit == fewerEndScore || fewerEndScore >= candidate->bestPossibleScore);

#ifdef _DEBUG
        if (_DumpAlignments) {
            printf("Scored fewer end candidate %d, set pair %d, read %d, location %s:%llu, seed offset %d, score limit %d, score %d, offset %d, agScore %d, matchProb %e\n", 
                (int)(candidate - scoringCandidatePool),
                candidate->whichSetPair, readWithFewerHits, 
                genome->getContigAtLocation(candidate->readWithFewerHitsGenomeLocation)->name, 
                candidate->readWithFewerHitsGenomeLocation - genome->getContigAtLocation(candidate->readWithFewerHitsGenomeLocation)->beginningLocation,
                candidate->seedOffset,
                scoreLimit, fewerEndScore, fewerEndGenomeLocationOffset, candidate->agScore, fewerEndMatchProbability);
        }
#endif // DEBUG

        if (fewerEndScore != ScoreAboveLimit) {
            //
            // Find and score mates.  The index in scoringMateCandidateIndex is the lowest mate (i.e., the highest index number).
            //
            unsigned mateIndex = candidate->scoringMateCandidateIndex;

            for (;;) {

                ScoringMateCandidate *mate = &scoringMateCandidates[candidate->whichSetPair][mateIndex];
                scoreLimit = computeScoreLimit(nonALTAlignment, &scoresForAllAlignments, &scoresForNonAltAlignments, mate->largestBigIndelDetected);

                _ASSERT(genomeLocationIsWithin(mate->readWithMoreHitsGenomeLocation, candidate->readWithFewerHitsGenomeLocation, maxSpacing));
                if (!genomeLocationIsWithin(mate->readWithMoreHitsGenomeLocation, candidate->readWithFewerHitsGenomeLocation, minSpacing) && mate->bestPossibleScore <= scoreLimit - fewerEndScore) {
                    //
                    // It's within the range and not necessarily too poor of a match.  Consider it.
                    //

                    //
                    // If we haven't yet scored this mate, or we've scored it and not gotten an answer, but had a higher score limit than we'd
                    // use now, score it.
                    //
                    if (mate->score == ScoringMateCandidate::LocationNotYetScored || (mate->score == ScoreAboveLimit && mate->scoreLimit < scoreLimit - fewerEndScore)) {
	                    scoreLocation(readWithMoreHits, setPairDirection[candidate->whichSetPair][readWithMoreHits], GenomeLocationAsInt64(mate->readWithMoreHitsGenomeLocation),
                            mate->seedOffset, scoreLimit - fewerEndScore, &mate->score, &mate->matchProbability,
                            &mate->genomeOffset, &mate->usedAffineGapScoring, &mate->basesClippedBefore, &mate->basesClippedAfter, &mate->agScore, &mate->lvIndels);
#ifdef _DEBUG
                        if (_DumpAlignments) {
                            printf("Scored mate candidate %d, set pair %d, read %d, location %s:%llu, seed offset %d, score limit %d, score %d, offset %d, agScore %d, matchProb %e\n",
                                (int)(mate - scoringMateCandidates[candidate->whichSetPair]), candidate->whichSetPair, readWithMoreHits, 
                                genome->getContigAtLocation(mate->readWithMoreHitsGenomeLocation)->name,
                                mate->readWithMoreHitsGenomeLocation - genome->getContigAtLocation(mate->readWithMoreHitsGenomeLocation)->beginningLocation,
                                mate->seedOffset, scoreLimit - fewerEndScore, mate->score, mate->genomeOffset, mate->agScore, mate->matchProbability);
                        }
#endif // _DEBUG

#if INSTRUMENTATION_FOR_PAPER
                        nScored++;
#endif // INSTRUMENTATION_FOR_PAPER

                        _ASSERT(ScoreAboveLimit == mate->score || mate->score >= mate->bestPossibleScore);

                        mate->scoreLimit = scoreLimit - fewerEndScore;
                    }

                    if (mate->score != ScoreAboveLimit && fewerEndScore + mate->score <= scoreLimit) { // We need to check to see that we're below scoreLimit because we may have scored this earlier when scoreLimit was higher.
                        double pairProbability = mate->matchProbability * fewerEndMatchProbability;
                        int pairScore = mate->score + fewerEndScore;
                        int pairAGScore = mate->agScore + candidate->agScore;
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
                                            mergeCandidate++) {

                                    if (mergeCandidate->mergeAnchor != NULL) {
                                        candidate->mergeAnchor = mergeAnchor = mergeCandidate->mergeAnchor;
                                        break;
                                    }
                                }
                            }
                        }

                        bool eliminatedByMerge; // Did we merge away this result.  If this is false, we may still have merged away a previous result.

                        double oldPairProbability;

                        bool mergeReplacement = false; // Did we replace the anchor with the new candidate ?

                        if (NULL == mergeAnchor) {
                            if (firstFreeMergeAnchor >= mergeAnchorPoolSize) {
                                WriteErrorMessage("Ran out of merge anchor pool entries.  Perhaps rerunning with a larger value of -mcp will help\n");
                                soft_exit(1);
                            }

                            mergeAnchor = &mergeAnchorPool[firstFreeMergeAnchor];

                            firstFreeMergeAnchor++;

                            mergeAnchor->init(mate->readWithMoreHitsGenomeLocation + mate->genomeOffset, candidate->readWithFewerHitsGenomeLocation + fewerEndGenomeLocationOffset,
                                pairProbability, pairScore, pairAGScore);

                            eliminatedByMerge = false;
                            oldPairProbability = 0;
                            candidate->mergeAnchor = mergeAnchor;
                        } else {
                            eliminatedByMerge = mergeAnchor->checkMerge(mate->readWithMoreHitsGenomeLocation + mate->genomeOffset, candidate->readWithFewerHitsGenomeLocation + fewerEndGenomeLocationOffset,
                                pairProbability, pairScore, pairAGScore, &oldPairProbability, &mergeReplacement);
                        }

                        if (!eliminatedByMerge) {
                            //
                            // Back out the probability of the old match that we're merged with, if any.  The max
                            // is necessary because a + b - b is not necessarily a in floating point.  If there
                            // was no merge, the oldPairProbability is 0. 
                            //

                            scoresForAllAlignments.updateProbabilityOfAllPairs(oldPairProbability);
                            if (nonALTAlignment) {
                                scoresForNonAltAlignments.updateProbabilityOfAllPairs(oldPairProbability);
                            }

                            if (pairProbability > scoresForAllAlignments.probabilityOfBestPair && maxEditDistanceForSecondaryResults != -1 && maxEditDistanceForSecondaryResults >= scoresForAllAlignments.bestPairScore - pairScore) {
                                //
                                // Move the old best to be a secondary alignment.  This won't happen on the first time we get a valid alignment,
                                // because bestPairScore is initialized to be very large.
                                //
                                //
                                if (*nSecondaryResults >= secondaryResultBufferSize) {
                                    *nSecondaryResults = secondaryResultBufferSize + 1;
                                    return false;
                                }

                                PairedAlignmentResult *secondaryResult = &secondaryResults[*nSecondaryResults];
                                secondaryResult->alignedAsPair = true;

                                for (int r = 0; r < NUM_READS_PER_PAIR; r++) {
                                    secondaryResult->direction[r] = scoresForAllAlignments.bestResultDirection[r];
                                    secondaryResult->location[r] = scoresForAllAlignments.bestResultGenomeLocation[r];
                                    secondaryResult->origLocation[r] = scoresForAllAlignments.bestResultOrigGenomeLocation[r];
                                    secondaryResult->mapq[r] = 0;
                                    secondaryResult->score[r] = scoresForAllAlignments.bestResultScore[r];
                                    secondaryResult->status[r] = MultipleHits;
                                    secondaryResult->usedAffineGapScoring[r] = scoresForAllAlignments.bestResultUsedAffineGapScoring[r];
                                    secondaryResult->basesClippedBefore[r] = scoresForAllAlignments.bestResultBasesClippedBefore[r];
                                    secondaryResult->basesClippedAfter[r] = scoresForAllAlignments.bestResultBasesClippedAfter[r];
                                    secondaryResult->agScore[r] = scoresForAllAlignments.bestResultAGScore[r];
                                    secondaryResult->seedOffset[r] = scoresForAllAlignments.bestResultSeedOffset[r];
                                    secondaryResult->popularSeedsSkipped[r] = popularSeedsSkipped[r];
                                    secondaryResult->lvIndels[r] = scoresForAllAlignments.bestResultLVIndels[r];
                                    secondaryResult->matchProbability[r] = scoresForAllAlignments.bestResultMatchProbability[r];
                                }
 
                                (*nSecondaryResults)++;

                            } // If we're saving the old best score as a secondary result

                            if (!mergeReplacement && (pairProbability > scoresForAllAlignments.probabilityOfBestPair) && (maxLVCandidatesForAffineGapBufferSize > 0) && (extraSearchDepth >= scoresForAllAlignments.bestPairScore - pairScore)) {
                                //
                                // This is close enough that scoring it with affine gap scoring might make it be the best result.  Save it for possible consideration in pase 4.
                                //
                                if (*nLVCandidatesForAffineGap >= maxLVCandidatesForAffineGapBufferSize) {
                                    *nLVCandidatesForAffineGap = maxLVCandidatesForAffineGapBufferSize + 1;
                                    return false;
                                }
                                PairedAlignmentResult *agResult = &lvCandidatesForAffineGap[*nLVCandidatesForAffineGap];
                                agResult->alignedAsPair = true;

                                for (int r = 0; r < NUM_READS_PER_PAIR; r++) {
                                    agResult->direction[r] = scoresForAllAlignments.bestResultDirection[r];
                                    agResult->location[r] = scoresForAllAlignments.bestResultGenomeLocation[r];
                                    agResult->origLocation[r] = scoresForAllAlignments.bestResultOrigGenomeLocation[r];
                                    agResult->mapq[r] = 0;
                                    agResult->score[r] = scoresForAllAlignments.bestResultScore[r];
                                    agResult->status[r] = MultipleHits;
                                    agResult->usedAffineGapScoring[r] = scoresForAllAlignments.bestResultUsedAffineGapScoring[r];
                                    agResult->basesClippedBefore[r] = scoresForAllAlignments.bestResultBasesClippedBefore[r];
                                    agResult->basesClippedAfter[r] = scoresForAllAlignments.bestResultBasesClippedAfter[r];
                                    agResult->agScore[r] = scoresForAllAlignments.bestResultAGScore[r];
                                    agResult->seedOffset[r] = scoresForAllAlignments.bestResultSeedOffset[r];
                                    agResult->popularSeedsSkipped[r] = popularSeedsSkipped[r];
                                    agResult->lvIndels[r] = scoresForAllAlignments.bestResultLVIndels[r];
                                    agResult->matchProbability[r] = scoresForAllAlignments.bestResultMatchProbability[r];
                                }

                                (*nLVCandidatesForAffineGap)++;
                            }

                            if (nonALTAlignment) {
                                scoresForNonAltAlignments.updateBestHitIfNeeded(pairScore, pairAGScore, pairProbability, fewerEndScore, readWithMoreHits, fewerEndGenomeLocationOffset, candidate, mate);
                            }

                            bool updatedBestScore = scoresForAllAlignments.updateBestHitIfNeeded(pairScore, pairAGScore, pairProbability, fewerEndScore, readWithMoreHits, fewerEndGenomeLocationOffset, candidate, mate);
#if INSTRUMENTATION_FOR_PAPER
                            if (updatedBestScore) {
                                bestCandidateScoredFirst = nScored == 2;
                            }
#endif // INSTRUMENTATION_FOR_PAPER


                            // scoreLimit = computeScoreLimit(nonALTAlignment, &scoresForAllAlignments, &scoresForNonAltAlignments);
                            
                            if ((!updatedBestScore) && maxEditDistanceForSecondaryResults != -1 && pairScore <= maxK && maxEditDistanceForSecondaryResults >= pairScore - scoresForAllAlignments.bestPairScore) {
                                
                                //
                                // A secondary result to save.
                                //
                                if (*nSecondaryResults >= secondaryResultBufferSize) {
                                    *nSecondaryResults = secondaryResultBufferSize + 1;
                                    return false;
                                }

                                PairedAlignmentResult *result = &secondaryResults[*nSecondaryResults];
                                result->alignedAsPair = true;
                                result->direction[readWithMoreHits] = setPairDirection[candidate->whichSetPair][readWithMoreHits];
                                result->direction[readWithFewerHits] = setPairDirection[candidate->whichSetPair][readWithFewerHits];
                                result->location[readWithMoreHits] = mate->readWithMoreHitsGenomeLocation + mate->genomeOffset;
                                result->location[readWithFewerHits] = candidate->readWithFewerHitsGenomeLocation + fewerEndGenomeLocationOffset;
                                result->origLocation[readWithMoreHits] = mate->readWithMoreHitsGenomeLocation;
                                result->origLocation[readWithFewerHits] = candidate->readWithFewerHitsGenomeLocation;
                                result->mapq[0] = result->mapq[1] = 0;
                                result->score[readWithMoreHits] = mate->score;
                                result->score[readWithFewerHits] = fewerEndScore;
                                result->status[readWithFewerHits] = result->status[readWithMoreHits] = MultipleHits;
                                result->usedAffineGapScoring[readWithMoreHits] = mate->usedAffineGapScoring;
                                result->usedAffineGapScoring[readWithFewerHits] = candidate->usedAffineGapScoring;
                                result->basesClippedBefore[readWithFewerHits] = candidate->basesClippedBefore;
                                result->basesClippedAfter[readWithFewerHits] = candidate->basesClippedAfter;
                                result->basesClippedBefore[readWithMoreHits] = mate->basesClippedBefore;
                                result->basesClippedAfter[readWithMoreHits] = mate->basesClippedAfter;
                                result->agScore[readWithMoreHits] = mate->agScore;
                                result->agScore[readWithFewerHits] = candidate->agScore;
                                result->seedOffset[readWithMoreHits] = mate->seedOffset;
                                result->seedOffset[readWithFewerHits] = candidate->seedOffset;
                                result->lvIndels[readWithMoreHits] = mate->lvIndels;
                                result->lvIndels[readWithFewerHits] = candidate->lvIndels;
                                result->matchProbability[readWithMoreHits] = mate->matchProbability;
                                result->matchProbability[readWithFewerHits] = candidate->matchProbability;
                                result->popularSeedsSkipped[readWithMoreHits] = popularSeedsSkipped[readWithMoreHits];
                                result->popularSeedsSkipped[readWithFewerHits] = popularSeedsSkipped[readWithFewerHits];

                                (*nSecondaryResults)++;
                            }


                            if ((!updatedBestScore) && maxLVCandidatesForAffineGapBufferSize > 0 && pairScore <= maxK && (extraSearchDepth >= pairScore - scoresForAllAlignments.bestPairScore)) {

                                if (*nLVCandidatesForAffineGap >= maxLVCandidatesForAffineGapBufferSize) {
                                    *nLVCandidatesForAffineGap = maxLVCandidatesForAffineGapBufferSize + 1;
                                    return false;
                                }

                                PairedAlignmentResult *result = &lvCandidatesForAffineGap[*nLVCandidatesForAffineGap];
                                result->alignedAsPair = true;
                                result->direction[readWithMoreHits] = setPairDirection[candidate->whichSetPair][readWithMoreHits];
                                result->direction[readWithFewerHits] = setPairDirection[candidate->whichSetPair][readWithFewerHits];
                                result->location[readWithMoreHits] = mate->readWithMoreHitsGenomeLocation + mate->genomeOffset;
                                result->location[readWithFewerHits] = candidate->readWithFewerHitsGenomeLocation + fewerEndGenomeLocationOffset;
                                result->origLocation[readWithMoreHits] = mate->readWithMoreHitsGenomeLocation;
                                result->origLocation[readWithFewerHits] = candidate->readWithFewerHitsGenomeLocation;
                                result->mapq[0] = result->mapq[1] = 0;
                                result->score[readWithMoreHits] = mate->score;
                                result->score[readWithFewerHits] = fewerEndScore;
                                result->status[readWithFewerHits] = result->status[readWithMoreHits] = MultipleHits;
                                result->usedAffineGapScoring[readWithMoreHits] = mate->usedAffineGapScoring;
                                result->usedAffineGapScoring[readWithFewerHits] = candidate->usedAffineGapScoring;
                                result->basesClippedBefore[readWithFewerHits] = candidate->basesClippedBefore;
                                result->basesClippedAfter[readWithFewerHits] = candidate->basesClippedAfter;
                                result->basesClippedBefore[readWithMoreHits] = mate->basesClippedBefore;
                                result->basesClippedAfter[readWithMoreHits] = mate->basesClippedAfter;
                                result->agScore[readWithMoreHits] = mate->agScore;
                                result->agScore[readWithFewerHits] = candidate->agScore;
                                result->seedOffset[readWithMoreHits] = mate->seedOffset;
                                result->seedOffset[readWithFewerHits] = candidate->seedOffset;
                                result->lvIndels[readWithMoreHits] = mate->lvIndels;
                                result->lvIndels[readWithFewerHits] = candidate->lvIndels;
                                result->matchProbability[readWithMoreHits] = mate->matchProbability;
                                result->matchProbability[readWithFewerHits] = candidate->matchProbability;
                                result->popularSeedsSkipped[readWithMoreHits] = popularSeedsSkipped[readWithMoreHits];
                                result->popularSeedsSkipped[readWithFewerHits] = popularSeedsSkipped[readWithFewerHits];

                                (*nLVCandidatesForAffineGap)++;
                            }
                            
    #ifdef  _DEBUG
                            if (_DumpAlignments) {
                                printf("Added %e (= %e * %e) @ (%s:%llu, %s:%llu), giving new probability of all pairs %e, score %d = %d + %d, agScore %d = %d + %d%s\n",
                                    pairProbability, mate->matchProbability , fewerEndMatchProbability,
                                    genome->getContigAtLocation(candidate->readWithFewerHitsGenomeLocation.location + fewerEndGenomeLocationOffset)->name,
                                    (candidate->readWithFewerHitsGenomeLocation + fewerEndGenomeLocationOffset) - genome->getContigAtLocation(candidate->readWithFewerHitsGenomeLocation.location + fewerEndGenomeLocationOffset)->beginningLocation,
                                    genome->getContigAtLocation(mate->readWithMoreHitsGenomeLocation + mate->genomeOffset)->name,
                                    (mate->readWithMoreHitsGenomeLocation + mate->genomeOffset) - genome->getContigAtLocation(mate->readWithMoreHitsGenomeLocation.location + mate->genomeOffset)->beginningLocation,
                                    scoresForNonAltAlignments.probabilityOfAllPairs,
                                    pairScore, fewerEndScore, mate->score, candidate->agScore + mate->agScore, candidate->agScore, mate->agScore, updatedBestScore ? " New best hit" : "");
                            }
    #endif  // _DEBUG

                            if ((altAwareness ? scoresForNonAltAlignments.probabilityOfAllPairs : scoresForAllAlignments.probabilityOfAllPairs) >= 4.9 && -1 == maxEditDistanceForSecondaryResults) {
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

     ScoreSet* scoreSetToEmit;
     if ((!altAwareness) || scoresForNonAltAlignments.bestPairScore > scoresForAllAlignments.bestPairScore + maxScoreGapToPreferNonAltAlignment) {
         scoreSetToEmit = &scoresForAllAlignments;
     } else {
         scoreSetToEmit = &scoresForNonAltAlignments;
     }

    if (scoreSetToEmit->bestPairScore == TooBigScoreValue) {
        //
        // Found nothing.
        //
        for (unsigned whichRead = 0; whichRead < NUM_READS_PER_PAIR; whichRead++) {
            result->location[whichRead] = InvalidGenomeLocation;
			result->origLocation[whichRead] = InvalidGenomeLocation;
            result->mapq[whichRead] = 0;
            result->score[whichRead] = ScoreAboveLimit;
            result->status[whichRead] = NotFound;
            result->clippingForReadAdjustment[whichRead] = 0;
            result->usedAffineGapScoring[whichRead] = false;
            result->basesClippedBefore[whichRead] = 0;
            result->basesClippedAfter[whichRead] = 0;
            result->agScore[whichRead] = ScoreAboveLimit;
            result->seedOffset[whichRead] = 0;
            result->lvIndels[whichRead] = 0;
            result->popularSeedsSkipped[whichRead] = popularSeedsSkipped[whichRead];
            result->matchProbability[whichRead] = 0.0;

            firstALTResult->status[whichRead] = NotFound;
#ifdef  _DEBUG
            if (_DumpAlignments) {
                printf("No sufficiently good pairs found.\n");
            }
#endif  // DEBUG
        }
        result->probabilityAllPairs = 0.0;

        
    } else {
        scoreSetToEmit->fillInResult(result, popularSeedsSkipped);
        if (altAwareness && scoreSetToEmit == &scoresForNonAltAlignments &&
            (scoresForAllAlignments.bestResultGenomeLocation[0] != scoresForNonAltAlignments.bestResultGenomeLocation[0] || 
             scoresForAllAlignments.bestResultGenomeLocation[1] != scoresForNonAltAlignments.bestResultGenomeLocation[1]))

        {
            _ASSERT(genome->isGenomeLocationALT(scoresForAllAlignments.bestResultGenomeLocation[0]));
            scoresForAllAlignments.fillInResult(firstALTResult, popularSeedsSkipped);
            for (int whichRead = 0; whichRead < NUM_READS_PER_PAIR; whichRead++)
            {
                firstALTResult->supplementary[whichRead] = true;
            }
        } else {
            for (int whichRead = 0; whichRead < NUM_READS_PER_PAIR; whichRead++)
            {
                firstALTResult->status[whichRead] = NotFound;
            }
        }
#ifdef  _DEBUG
            if (_DumpAlignments) {
                printf("Returned %s:%llu %s %s:%llu %s with MAPQ %d and %d, probability of all pairs %e, probability of best pair %e, pair score %d\n",
                    genome->getContigAtLocation(result->location[0])->name, result->location[0] - genome->getContigAtLocation(result->location[0])->beginningLocation,
                    result->direction[0] == RC ? "RC" : "", 
                    genome->getContigAtLocation(result->location[1])->name, result->location[1] - genome->getContigAtLocation(result->location[1])->beginningLocation, 
                    result->direction[1] == RC ? "RC" : "", result->mapq[0], result->mapq[1], scoreSetToEmit->probabilityOfAllPairs, scoreSetToEmit->probabilityOfBestPair,
                    scoreSetToEmit->bestPairScore);

                if (firstALTResult->status[0] != NotFound) {
                    printf("Returned first ALT Result %s:%llu %s %s:%llu %s with MAPQ %d and %d, probability of all pairs %e, probability of best pair %e, pair score %d\n",
                        genome->getContigAtLocation(firstALTResult->location[0])->name, firstALTResult->location[0] - genome->getContigAtLocation(firstALTResult->location[0])->beginningLocation,
                        firstALTResult->direction[0] == RC ? "RC" : "",
                        genome->getContigAtLocation(firstALTResult->location[1])->name, firstALTResult->location[1] - genome->getContigAtLocation(firstALTResult->location[1])->beginningLocation,
                        firstALTResult->direction[1] == RC ? "RC" : "", firstALTResult->mapq[0], firstALTResult->mapq[1], scoresForAllAlignments.probabilityOfAllPairs, scoresForAllAlignments.probabilityOfBestPair,
                        scoresForAllAlignments.bestPairScore);
                } // If we're also returning an ALT result
            }
#endif  // DEBUG
    }

    //
    // Get rid of any secondary results that are too far away from the best score.  (NB: the rest of the code in align() is very similar to BaseAligner::finalizeSecondaryResults.  Sorry)
    //


    Read *inputReads[2] = { read0, read1 };
    for (int whichRead = 0; whichRead < NUM_READS_PER_PAIR; whichRead++) {
        result->scorePriorToClipping[whichRead] = result->score[whichRead];
    }

    if (!ignoreAlignmentAdjustmentsForOm) {
        //
        // Start by adjusting the alignments.
        //
        alignmentAdjuster.AdjustAlignments(inputReads, result);
        if (result->status[0] != NotFound && result->status[1] != NotFound && !ignoreAlignmentAdjustmentsForOm) {
            scoreSetToEmit->bestPairScore = result->score[0] + result->score[1];
        }

        for (int i = 0; i < *nSecondaryResults; i++) {
            for (int whichRead = 0; whichRead < NUM_READS_PER_PAIR; whichRead++) {
                secondaryResults[i].scorePriorToClipping[whichRead] = secondaryResults[i].score[whichRead];
            }
            alignmentAdjuster.AdjustAlignments(inputReads, &secondaryResults[i]);
            if (secondaryResults[i].status[0] != NotFound && secondaryResults[i].status[1] != NotFound && !ignoreAlignmentAdjustmentsForOm) {
                scoreSetToEmit->bestPairScore = __min(scoreSetToEmit->bestPairScore, secondaryResults[i].score[0] + secondaryResults[i].score[1]);
            }
        }
    } else {
        for (int i = 0; i < *nSecondaryResults; i++) {
            for (int whichRead = 0; whichRead < NUM_READS_PER_PAIR; whichRead++) {
                secondaryResults[i].scorePriorToClipping[whichRead] = secondaryResults[i].score[whichRead];
            }
        }
    }

    int i = 0;
    while (i < *nSecondaryResults) {
        if ((int)(secondaryResults[i].score[0] + secondaryResults[i].score[1]) >(int)scoreSetToEmit->bestPairScore + maxEditDistanceForSecondaryResults ||
            secondaryResults[i].status[0] == NotFound || secondaryResults[i].status[1] == NotFound) {

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

#if INSTRUMENTATION_FOR_PAPER
    _int64 runTime = timeInNanos() - startTime;
    if (runTime >= 0) { // Really don't understand why timeInNanos() sometimes produces garbage, but it does.
        InterlockedAdd64AndReturnNewValue(&g_alignmentCountByHitCountsOfEachSeed[log2HashTableHits[0]][log2HashTableHits[1]], 1);
        InterlockedAdd64AndReturnNewValue(&g_alignmentTimeByHitCountsOfEachSeed[log2HashTableHits[0]][log2HashTableHits[1]], runTime);
        InterlockedAdd64AndReturnNewValue(&g_scoreCountByHitCountsOfEachSeed[log2HashTableHits[0]][log2HashTableHits[1]], nScored);
        InterlockedAdd64AndReturnNewValue(&g_setIntersectionSizeByHitCountsOfEachSeed[log2HashTableHits[0]][log2HashTableHits[1]], setIntersectionSize);
        if (__min(hashTableHits[0], hashTableHits[1]) > 0) {
            InterlockedAdd64AndReturnNewValue(&g_100xtotalRatioOfSetIntersectionSizeToSmallerSeedHitCountByCountsOfEachSeed[log2HashTableHits[0]][log2HashTableHits[1]], 
                (100 * setIntersectionSize) / __min(hashTableHits[0], hashTableHits[1]));
        }
        InterlockedAdd64AndReturnNewValue(&g_totalSizeOfSmallerHitSet, __min(hashTableHits[0], hashTableHits[1]));
        InterlockedAdd64AndReturnNewValue(&g_totalSizeOfSetIntersection, setIntersectionSize);
    }

    if (nScored > 2) {
        InterlockedAdd64AndReturnNewValue(&g_alignmentsWithMoreThanOneCandidate, 1);
        if (bestCandidateScoredFirst) {
            InterlockedAdd64AndReturnNewValue(&g_alignmentsWithMoreThanOneCandidateWhereTheBestCandidateIsScoredFirst, 1);
        }
    }
#endif // INSTRUMENTATION_FOR_PAPER


    return true;
}

    bool
IntersectingPairedEndAligner::alignAffineGap(
        Read                  *read0,
        Read                  *read1,
        PairedAlignmentResult* result,
        PairedAlignmentResult* firstALTResult,
        int                    maxEditDistanceForSecondaryResults,
        _int64                 secondaryResultBufferSize,
        _int64                *nSecondaryResults,
        PairedAlignmentResult *secondaryResults,             // The caller passes in a buffer of secondaryResultBufferSize and it's filled in by align()
        _int64                 singleSecondaryBufferSize,
        _int64                 maxSecondaryResultsToReturn,
        _int64                *nSingleEndSecondaryResultsForFirstRead,
        _int64                *nSingleEndSecondaryResultsForSecondRead,
        SingleAlignmentResult *singleEndSecondaryResults,     // Single-end secondary alignments for when the paired-end alignment didn't work properly
        _int64                 maxLVCandidatesForAffineGapBufferSize,
        _int64                *nLVCandidatesForAffineGap,
        PairedAlignmentResult *lvCandidatesForAffineGap
    )
{
    //
    // Phase 4: Re-score candidates that need to be using affine gap scoring, and change the result if necessary.
    //

    if (result->status[0] == NotFound || result->status[1] == NotFound) {
        return true;
    }

    //
    // We rebuild the RC reads and extract the reference again before scoring candidates using affine gap
    //
    Read rcReads[NUM_READS_PER_PAIR];

    reads[0][FORWARD] = read0;
    reads[1][FORWARD] = read1;

    //
    // Don't bother if one or both reads are too short.  The minimum read length here is the seed length, but usually there's a longer
    // minimum enforced by our caller
    //
    if (read0->getDataLength() < seedLen || read1->getDataLength() < seedLen) {
        return true;
    }

    //
    // Build the RC reads.
    //
    unsigned countOfNs = 0;

    for (unsigned whichRead = 0; whichRead < NUM_READS_PER_PAIR; whichRead++) {
        Read *read = reads[whichRead][FORWARD];
        readLen[whichRead] = read->getDataLength();

        for (unsigned i = 0; i < reads[whichRead][FORWARD]->getDataLength(); i++) {
            rcReadData[whichRead][i] = rcTranslationTable[read->getData()[readLen[whichRead] - i - 1]];
            rcReadQuality[whichRead][i] = read->getQuality()[readLen[whichRead] - i - 1];
            countOfNs += nTable[read->getData()[i]];
        }

        reads[whichRead][RC] = &rcReads[whichRead];
        reads[whichRead][RC]->init(read->getId(), read->getIdLength(), rcReadData[whichRead], rcReadQuality[whichRead], read->getDataLength());
    }

    if ((int)countOfNs > maxK) {
        return true;
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

    _ASSERT(maxLVCandidatesForAffineGapBufferSize > 0);

    int maxKForSameAlignment = gapOpenPenalty / (subPenalty - gapExtendPenalty);
    int bestPairScore = result->score[0] + result->score[1];
    int scoreLimit, scoreLimitALT;
    scoreLimit = scoreLimitALT = maxKForIndels + extraSearchDepth;  // XXX maxKForIndels is too much!
    int genomeOffset[NUM_READS_PER_PAIR] = { 0, 0 };
    bool skipAffineGap[NUM_READS_PER_PAIR] = { false, false };

    //
    // Keep track of old bestPairProbability as this is used in updating the new match probability after affine gap scoring
    //
    double oldPairProbabilityBestResult = result->matchProbability[0] * result->matchProbability[1];
    double oldPairProbabilityBestResultALT = (firstALTResult->status[0] != NotFound) ? firstALTResult->matchProbability[0] * firstALTResult->matchProbability[1] : 0.0;

    for (int r = 0; r < NUM_READS_PER_PAIR; r++) {
        if (result->score[r] > maxKForSameAlignment) {
            //
            // Use affine gap scoring to determine if bases need to be clipped
            //
            result->usedAffineGapScoring[r] = true;
            scoreLocationWithAffineGap(r, result->direction[r], result->origLocation[r],
	            result->seedOffset[r], scoreLimit, &result->score[r], &result->matchProbability[r],
	            &genomeOffset[r], &result->basesClippedBefore[r], &result->basesClippedAfter[r], &result->agScore[r]);
#if _DEBUG
            if (_DumpAlignments) {
                fprintf(stderr, "Affine gap scored read %d at %s:%llu result %d\n", 
                        r, genome->getContigAtLocation(result->origLocation[r])->name, result->origLocation[r] - genome->getContigAtLocation(result->origLocation[r])->beginningLocation,
                        result->agScore[r]
                        );
            }
#endif  // _DEBUG

            if (result->score[r] != ScoreAboveLimit) {
                result->location[r] = result->origLocation[r] + genomeOffset[r];
                scoreLimit -= result->score[r];
            } else {
                result->status[r] = NotFound;
            }

            //
            // Use affine gap scoring for ALT result if it was computed in Phase 3
            //
            if (firstALTResult->status[r] != NotFound) {
                if (firstALTResult->score[r] > maxKForSameAlignment) { // affine gap may produce a better alignment
                    firstALTResult->usedAffineGapScoring[r] = true;
                    scoreLocationWithAffineGap(r, firstALTResult->direction[r], firstALTResult->origLocation[r],
                        firstALTResult->seedOffset[r], scoreLimitALT, &firstALTResult->score[r], &firstALTResult->matchProbability[r],
                        &genomeOffset[r], &firstALTResult->basesClippedBefore[r], &firstALTResult->basesClippedAfter[r], &firstALTResult->agScore[r]);

                    if (firstALTResult->score[r] != ScoreAboveLimit) {
                        firstALTResult->location[r] = firstALTResult->origLocation[r] + genomeOffset[r];
                        scoreLimitALT -= firstALTResult->score[r];
                    }
                    else {
                        firstALTResult->status[r] = NotFound;
                    }
                }
            }
        } else {
            //
            // Skip affine gap scoring for reads we know that LV and affine gap will agree on the alignment
            //
            result->usedAffineGapScoring[r] = false;
            skipAffineGap[r] = true;
        }
    }

    if (result->status[0] == NotFound || result->status[1] == NotFound) {
        //
        // Found nothing from the paired-end aligner if one of the reads in the pair is unmapped.
        //
        for (unsigned whichRead = 0; whichRead < NUM_READS_PER_PAIR; whichRead++) {
            result->location[whichRead] = InvalidGenomeLocation;
            result->origLocation[whichRead] = InvalidGenomeLocation;
            result->mapq[whichRead] = 0;
            result->score[whichRead] = ScoreAboveLimit;
            result->status[whichRead] = NotFound;
            result->clippingForReadAdjustment[whichRead] = 0;
            result->usedAffineGapScoring[whichRead] = false;
            result->basesClippedBefore[whichRead] = 0;
            result->basesClippedAfter[whichRead] = 0;
            result->agScore[whichRead] = ScoreAboveLimit;
            result->seedOffset[whichRead] = 0;
            result->lvIndels[whichRead] = 0;
            result->matchProbability[whichRead] = 0.0;

            firstALTResult->status[whichRead] = NotFound;
        #ifdef  _DEBUG
            if (_DumpAlignments) {
                printf("No sufficiently good pairs found.\n");
            }
        #endif  // DEBUG
        }
        result->probabilityAllPairs = 0.0;
        return true;
    }

    ScoreSet scoresForAllAlignments;
    ScoreSet scoresForNonAltAlignments;

    //
    // In the beginning we only have the best alignment result in the score set.
    // It is important to initialize the score set here and not before affine gap scoring, since only affine gap does clipping of alignments
    //
    bool nonALTBestAlignment = (!altAwareness) || !genome->isGenomeLocationALT(result->location[0]);
    if (firstALTResult->status[0] != NotFound) { // best result is an ALT result
        scoresForAllAlignments.init(firstALTResult);
    }
    else {
        scoresForAllAlignments.init(result);
    }
    if (nonALTBestAlignment) {
        scoresForNonAltAlignments.init(result);
    }

    //
    // Update match probability for reads rescored with affine gap
    //
    if (!skipAffineGap[0] || !skipAffineGap[1]) {
        double newPairProbability = result->matchProbability[0] * result->matchProbability[1];
        if (firstALTResult->status[0] != NotFound) { // best result is an ALT result
            double newPairProbabilityALT = firstALTResult->matchProbability[0] * firstALTResult->matchProbability[1];
            scoresForAllAlignments.updateProbabilityOfAllPairs(oldPairProbabilityBestResultALT);
            scoresForAllAlignments.updateProbabilityOfBestPair(newPairProbabilityALT);
        }
        else {
            scoresForAllAlignments.updateProbabilityOfAllPairs(oldPairProbabilityBestResult);
            scoresForAllAlignments.updateProbabilityOfBestPair(newPairProbability);
        }
        if (nonALTBestAlignment) {
            scoresForNonAltAlignments.updateProbabilityOfAllPairs(oldPairProbabilityBestResult);
            scoresForNonAltAlignments.updateProbabilityOfBestPair(newPairProbability);
        }
    }

    //
    // Evaluate LV candidates with affine gap scoring
    //
    if ((*nLVCandidatesForAffineGap > 0) && (!skipAffineGap[0] || !skipAffineGap[1])) {

        //
        // Reset score limit
        //
        scoreLimit = bestPairScore + extraSearchDepth;

        //
        // We sort all all promising LV candidates and score them with affine gap starting with the best one
        //
        qsort(lvCandidatesForAffineGap, *nLVCandidatesForAffineGap, sizeof(*lvCandidatesForAffineGap), PairedAlignmentResult::compareByScore);

        for (int i = 0; i < *nLVCandidatesForAffineGap; i++) {

            PairedAlignmentResult* lvResult = &lvCandidatesForAffineGap[i];
            int lvPairScore = lvResult->score[0] + lvResult->score[1];
            int lvPairIndels = lvResult->lvIndels[0] + lvResult->lvIndels[1];

            //
            // Use the maximum scoreLimit if we expect the read could have alignments with large indels
            //
            if ((lvPairScore > bestPairScore + extraSearchDepth) && (lvPairIndels > 1)) {
                scoreLimit = maxK + extraSearchDepth;
            }

            if ((lvPairScore <= bestPairScore + extraSearchDepth) || (lvPairIndels > 1)) {
                _ASSERT(lvResult->status[0] != NotFound && lvResult->status[1] != NotFound);
                bool nonALTAlignment = (!altAwareness) || !genome->isGenomeLocationALT(lvResult->location[0]);
                double oldPairProbability = lvResult->matchProbability[0] * lvResult->matchProbability[1];

                if (!skipAffineGap[0]) {
                    //
                    // Score first read with affine gap
                    //
                    lvResult->usedAffineGapScoring[0] = true;
                    scoreLocationWithAffineGap(0, lvResult->direction[0], lvResult->origLocation[0],
	                    lvResult->seedOffset[0], scoreLimit, &lvResult->score[0], &lvResult->matchProbability[0],
	                    &genomeOffset[0], &lvResult->basesClippedBefore[0], &lvResult->basesClippedAfter[0], &lvResult->agScore[0]);
                }

                if (lvResult->score[0] != ScoreAboveLimit) {
                    lvResult->location[0] = lvResult->origLocation[0] + genomeOffset[0];

                    if (!skipAffineGap[1]) {
                        //
                        // Score mate with affine gap
                        //
                        lvResult->usedAffineGapScoring[1] = true;
                        scoreLocationWithAffineGap(1, lvResult->direction[1], lvResult->origLocation[1],
	                        lvResult->seedOffset[1], scoreLimit - lvResult->score[0], &lvResult->score[1], &lvResult->matchProbability[1],
	                        &genomeOffset[1], &lvResult->basesClippedBefore[1], &lvResult->basesClippedAfter[1], &lvResult->agScore[1]);
                    } // if !skipAffineGap[1]

                    if (lvResult->score[1] != ScoreAboveLimit) {
                        lvResult->location[1] = lvResult->origLocation[1] + genomeOffset[1];
                        double pairProbability = lvResult->matchProbability[0] * lvResult->matchProbability[1];
                        int pairScore = lvResult->score[0] + lvResult->score[1];
                        int pairAGScore = lvResult->agScore[0] + lvResult->agScore[1];

                        //
                        // Update match probabilities for read pair and best hit if better
                        //
                        scoresForAllAlignments.updateProbabilityOfAllPairs(oldPairProbability);
                        bool updatedBestScore = scoresForAllAlignments.updateBestHitIfNeeded(pairScore, pairAGScore, pairProbability, lvResult);
                        if (nonALTAlignment) {
	                        scoresForNonAltAlignments.updateProbabilityOfAllPairs(oldPairProbability);
	                        scoresForNonAltAlignments.updateBestHitIfNeeded(pairScore, pairAGScore, pairProbability, lvResult);
                        } // nonALTAlignment

                        //
                        // Update scoreLimit so that we only look for alignments extraSearchDepth worse than the best
                        //
                        scoreLimit = computeScoreLimit(nonALTAlignment, &scoresForAllAlignments, &scoresForNonAltAlignments, 0);
                    } // lvResult->score[1] != ScoreAboveLimit
                } // lvResult->score[0] != ScoreAboveLimit
            } // If we want to score this candidate with affine gap
        } // for each candidate from LV
    }

    //
    // Emit the final result (i.e., ALT/non-ALT best result and first ALT result, if any)
    //
    ScoreSet* scoreSetToEmit;
    if ((!altAwareness) || scoresForNonAltAlignments.bestPairScore > scoresForAllAlignments.bestPairScore + maxScoreGapToPreferNonAltAlignment) {
        scoreSetToEmit = &scoresForAllAlignments;
    }
    else {
        scoreSetToEmit = &scoresForNonAltAlignments;
    }

    scoreSetToEmit->fillInResult(result, result->popularSeedsSkipped);
    if (altAwareness && scoreSetToEmit == &scoresForNonAltAlignments &&
	    (scoresForAllAlignments.bestResultGenomeLocation[0] != scoresForNonAltAlignments.bestResultGenomeLocation[0] ||
		    scoresForAllAlignments.bestResultGenomeLocation[1] != scoresForNonAltAlignments.bestResultGenomeLocation[1]))

    {
        _ASSERT(genome->isGenomeLocationALT(scoresForAllAlignments.bestResultGenomeLocation[0]));
        scoresForAllAlignments.fillInResult(firstALTResult, firstALTResult->popularSeedsSkipped);
        for (int whichRead = 0; whichRead < NUM_READS_PER_PAIR; whichRead++)
        {
            firstALTResult->supplementary[whichRead] = true;
        }
    }
    else {
        for (int whichRead = 0; whichRead < NUM_READS_PER_PAIR; whichRead++)
        {
            firstALTResult->status[whichRead] = NotFound;
        }
    }

    return true;
}

	void
IntersectingPairedEndAligner::scoreLocationWithAffineGap(
	unsigned             whichRead,
	Direction            direction,
	GenomeLocation       genomeLocation,
	unsigned             seedOffset,
	int                  scoreLimit,
	int                 *score,
	double              *matchProbability,
	int                 *genomeLocationOffset,
	int                 *basesClippedBefore,
	int                 *basesClippedAfter,
	int                 *agScore
	)
{
    Read *readToScore = reads[whichRead][direction];
    unsigned readDataLength = readToScore->getDataLength();
    GenomeDistance genomeDataLength = readDataLength + MAX_K; // Leave extra space in case the read has deletions
    const char *data = genome->getSubstring(genomeLocation, genomeDataLength);

    *genomeLocationOffset = 0;

    if (NULL == data) {
        *score = ScoreAboveLimit;
        *matchProbability = 0;
        *genomeLocationOffset = 0;
        *agScore = ScoreAboveLimit;
        return;
    }

    *basesClippedBefore = 0;
    *basesClippedAfter = 0;

    double matchProb1 = 1.0, matchProb2 = 1.0;
    int score1 = 0, score2 = 0; // edit distance
    // First, do the forward direction from where the seed aligns to past of it
    int readLen = readToScore->getDataLength();
    int seedLen = index->getSeedLength();
    int tailStart = seedOffset + seedLen;
    int agScore1 = seedLen, agScore2 = 0; // affine gap scores

    _ASSERT(!memcmp(data+seedOffset, readToScore->getData() + seedOffset, seedLen));    // that the seed actually matches

    int textLen;
    if (genomeDataLength - tailStart > INT32_MAX) {
        textLen = INT32_MAX;
    } else {
        textLen = (int)(genomeDataLength - tailStart);
    }

    if (tailStart != readLen) {
        int patternLen = readLen - tailStart;
        //
        // Try banded affine-gap when pattern is long and band needed is small
        //
        if (patternLen >= (3 * (2 * (int)scoreLimit + 1))) {
            agScore1 = affineGap->computeScoreBanded(data + tailStart,
                textLen,
                readToScore->getData() + tailStart,
                readToScore->getQuality() + tailStart,
                readLen - tailStart,
                scoreLimit,
                readLen,
                NULL,
                basesClippedAfter,
                &score1,
                &matchProb1);
        }
        else {
            agScore1 = affineGap->computeScore(data + tailStart,
                textLen,
                readToScore->getData() + tailStart,
                readToScore->getQuality() + tailStart,
                readLen - tailStart,
                scoreLimit,
                readLen,
                NULL,
                basesClippedAfter,
                &score1,
                &matchProb1);
        }

        agScore1 += (seedLen - readLen);
    }
    if (score1 != ScoreAboveLimit) {
        if (seedOffset != 0) {
            int limitLeft = scoreLimit - score1;
            int patternLen = seedOffset;
            //
            // Try banded affine-gap when pattern is long and band needed is small
            //
            if (patternLen >= (3 * (2 * limitLeft + 1))) {
                agScore2 = reverseAffineGap->computeScoreBanded(data + seedOffset,
	                seedOffset + limitLeft,
	                reversedRead[whichRead][direction] + readLen - seedOffset,
	                reads[whichRead][OppositeDirection(direction)]->getQuality() + readLen - seedOffset,
	                seedOffset,
	                limitLeft,
	                readLen, // FIXME: Assumes the rest of the read matches perfectly
	                genomeLocationOffset,
	                basesClippedBefore,
	                &score2,
	                &matchProb2);
            }
            else {
                agScore2 = reverseAffineGap->computeScore(data + seedOffset,
	                seedOffset + limitLeft,
	                reversedRead[whichRead][direction] + readLen - seedOffset,
	                reads[whichRead][OppositeDirection(direction)]->getQuality() + readLen - seedOffset,
	                seedOffset,
	                limitLeft,
	                readLen,
	                genomeLocationOffset,
	                basesClippedBefore,
	                &score2,
	                &matchProb2);
            }

            agScore2 -= (readLen);

            if (score2 == ScoreAboveLimit) {
	            *score = ScoreAboveLimit;
	            *genomeLocationOffset = 0;
	            *agScore = -1;
            }
        }
    }
    else {
        *score = ScoreAboveLimit;
        *genomeLocationOffset = 0;
        *agScore = -1;
    }

    if (score1 != ScoreAboveLimit && score2 != ScoreAboveLimit) {
        *score = score1 + score2;
        // _ASSERT(*score <= scoreLimit);
        // Map probabilities for substrings can be multiplied, but make sure to count seed too
        *matchProbability = matchProb1 * matchProb2 * pow(1 - SNP_PROB, seedLen);

        *agScore = agScore1 + agScore2;
    }
    else {
        *score = ScoreAboveLimit;
        *agScore = -1;
        *matchProbability = 0.0;
    }
}

    void
IntersectingPairedEndAligner::scoreLocation(
    unsigned             whichRead,
    Direction            direction,
    GenomeLocation       genomeLocation,
    unsigned             seedOffset,
    int                  scoreLimit,
    int                 *score,
    double              *matchProbability,
    int                 *genomeLocationOffset,
    bool                *usedAffineGapScoring,
    int                 *basesClippedBefore,
    int                 *basesClippedAfter,
    int                 *agScore,
    int                 *totalIndelsLV)
{
    nLocationsScored++;

    if (noUkkonen) {
        scoreLimit = maxK + extraSearchDepth;
    }

    Read *readToScore = reads[whichRead][direction];
    unsigned readDataLength = readToScore->getDataLength();
    GenomeDistance genomeDataLength = readDataLength + MAX_K; // Leave extra space in case the read has deletions
    const char *data = genome->getSubstring(genomeLocation, genomeDataLength);

    *genomeLocationOffset = 0;

    if (NULL == data) {
        *score = ScoreAboveLimit;
        *matchProbability = 0;
        *genomeLocationOffset = 0;
        *agScore = ScoreAboveLimit;
        return;
    }

    *basesClippedBefore = 0;
    *basesClippedAfter = 0;

    // Compute the distance separately in the forward and backward directions from the seed, to allow
    // arbitrary offsets at both the start and end but not have to pay the cost of exploring all start
    // shifts in BoundedStringDistance
    double matchProb1 = 1.0, matchProb2 = 1.0;
    int score1 = 0, score2 = 0; // edit distance
    // First, do the forward direction from where the seed aligns to past of it
    int readLen = readToScore->getDataLength();
    int seedLen = index->getSeedLength();
    int tailStart = seedOffset + seedLen;
    int agScore1 = seedLen, agScore2 = 0; // affine gap scores

    _ASSERT(!memcmp(data+seedOffset, readToScore->getData() + seedOffset, seedLen));    // that the seed actually matches

    int textLen;
    if (genomeDataLength - tailStart > INT32_MAX) {
        textLen = INT32_MAX;
    } else {
        textLen = (int)(genomeDataLength - tailStart);
    }

    int totalIndels1 = 0, totalIndels2 = 0;

    score1 = landauVishkin->computeEditDistance(data + tailStart, textLen, readToScore->getData() + tailStart, readToScore->getQuality() + tailStart, readLen - tailStart,
        scoreLimit, &matchProb1, NULL, &totalIndels1);

    agScore1 = (seedLen + readLen - tailStart - score1) * matchReward - score1 * subPenalty;

    if (score1 != ScoreAboveLimit) {
        int limitLeft = scoreLimit - score1;
        score2 = reverseLandauVishkin->computeEditDistance(data + seedOffset, seedOffset + MAX_K, reversedRead[whichRead][direction] + readLen - seedOffset,
            reads[whichRead][OppositeDirection(direction)]->getQuality() + readLen - seedOffset, seedOffset, limitLeft, &matchProb2, genomeLocationOffset, &totalIndels2);

        agScore2 = (seedOffset - score2) * matchReward - score2 * subPenalty;
    }

    if (score1 != ScoreAboveLimit && score2 != ScoreAboveLimit) {
        *score = score1 + score2;
        _ASSERT(*score <= scoreLimit);
        // Map probabilities for substrings can be multiplied, but make sure to count seed too
        *matchProbability = matchProb1 * matchProb2 * pow(1 - SNP_PROB, seedLen);

        *agScore = agScore1 + agScore2;

        *totalIndelsLV = totalIndels1 + totalIndels2;
    }
    else {
        *score = ScoreAboveLimit;
        *agScore = -1;
        *matchProbability = 0.0;
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
                        int newPairAGScore, double *oldMatchProbability, bool *mergeReplacement)
{
    if (locationForReadWithMoreHits == InvalidGenomeLocation || !doesRangeMatch(newMoreHitLocation, newFewerHitLocation)) {
        //
        // No merge.  Remember the new one.
        //
        locationForReadWithMoreHits = newMoreHitLocation;
        locationForReadWithFewerHits = newFewerHitLocation;
        matchProbability = newMatchProbability;
        pairScore = newPairScore;
        pairAGScore = newPairAGScore;
        *oldMatchProbability = 0.0;
        *mergeReplacement = false;
        return false;
    }  else {
        //
        // Within merge distance.  Keep the better score (or if they're tied the better match probability).
        //
        if (newPairAGScore > pairAGScore || (newPairAGScore == pairAGScore && newMatchProbability > matchProbability)) {
#ifdef _DEBUG
            if (_DumpAlignments) {
                printf("Merge replacement at anchor (%llu, %llu), loc (%llu, %llu), old match prob %e, new match prob %e, old pair score %d, new pair score %d\n",
                    locationForReadWithMoreHits.location, locationForReadWithFewerHits.location, newMoreHitLocation.location, newFewerHitLocation.location,
                    matchProbability, newMatchProbability, pairScore, newPairScore);
            }
#endif // DEBUG

            * oldMatchProbability = matchProbability;
            matchProbability = newMatchProbability;
            pairScore = newPairScore;
            pairAGScore = newPairAGScore;
            *mergeReplacement = true;
            return false;
        } else {
            //
            // The new one should just be ignored.
            //
#ifdef _DEBUG
            if (_DumpAlignments) {
                printf("Merged at anchor (%llu, %llu), loc (%llu, %llu), old match prob %e, new match prob %e, old pair score %d, new pair score %d\n",
                    locationForReadWithMoreHits.location, locationForReadWithFewerHits.location, newMoreHitLocation.location, newFewerHitLocation.location,
                    matchProbability, newMatchProbability, pairScore, newPairScore);
            }
#endif // DEBUG
            *mergeReplacement = false;
            return true;
        }
    }

    _ASSERT(!"NOTREACHED");
}

void IntersectingPairedEndAligner::ScoreSet::updateProbabilityOfAllPairs(double oldPairProbability) 
{
    probabilityOfAllPairs = __max(0, probabilityOfAllPairs - oldPairProbability);
}

bool IntersectingPairedEndAligner::ScoreSet::updateBestHitIfNeeded(int pairScore, int pairAGScore, double pairProbability, int fewerEndScore, int readWithMoreHits, GenomeDistance fewerEndGenomeLocationOffset, ScoringCandidate* candidate, ScoringMateCandidate* mate)
{
    probabilityOfAllPairs += pairProbability;
    int readWithFewerHits = 1 - readWithMoreHits;

    if (pairAGScore > bestPairAGScore || (pairAGScore == bestPairAGScore && pairProbability > probabilityOfBestPair)) {
        bestPairScore = pairScore;
        bestPairAGScore = pairAGScore;
        probabilityOfBestPair = pairProbability;
        bestResultGenomeLocation[readWithFewerHits] = candidate->readWithFewerHitsGenomeLocation + fewerEndGenomeLocationOffset;
        bestResultGenomeLocation[readWithMoreHits] = mate->readWithMoreHitsGenomeLocation + mate->genomeOffset;
        bestResultOrigGenomeLocation[readWithFewerHits] = candidate->readWithFewerHitsGenomeLocation;
        bestResultOrigGenomeLocation[readWithMoreHits] = mate->readWithMoreHitsGenomeLocation;
        bestResultScore[readWithFewerHits] = fewerEndScore;
        bestResultScore[readWithMoreHits] = mate->score;
        bestResultDirection[readWithFewerHits] = setPairDirection[candidate->whichSetPair][readWithFewerHits];
        bestResultDirection[readWithMoreHits] = setPairDirection[candidate->whichSetPair][readWithMoreHits];
        bestResultUsedAffineGapScoring[readWithFewerHits] = candidate->usedAffineGapScoring;
        bestResultUsedAffineGapScoring[readWithMoreHits] = mate->usedAffineGapScoring;
        bestResultBasesClippedBefore[readWithFewerHits] = candidate->basesClippedBefore;
        bestResultBasesClippedAfter[readWithFewerHits] = candidate->basesClippedAfter;
        bestResultBasesClippedBefore[readWithMoreHits] = mate->basesClippedBefore;
        bestResultBasesClippedAfter[readWithMoreHits] = mate->basesClippedAfter;
        bestResultAGScore[readWithFewerHits] = candidate->agScore;
        bestResultAGScore[readWithMoreHits] = mate->agScore;
        bestResultSeedOffset[readWithFewerHits] = candidate->seedOffset;
        bestResultSeedOffset[readWithMoreHits] = mate->seedOffset;
        bestResultMatchProbability[readWithFewerHits] = candidate->matchProbability;
        bestResultMatchProbability[readWithMoreHits] = mate->matchProbability;
        bestResultLVIndels[readWithFewerHits] = candidate->lvIndels;
        bestResultLVIndels[readWithMoreHits] = mate->lvIndels;

        return true;
    } else {
        return false;
    }
} // updateBestHitIfNeeded

bool IntersectingPairedEndAligner::ScoreSet::updateBestHitIfNeeded(int pairScore, int pairAGScore, double pairProbability, PairedAlignmentResult* newResult)
{
    probabilityOfAllPairs += pairProbability;

    if (pairAGScore > bestPairAGScore || (pairAGScore == bestPairAGScore && pairProbability > probabilityOfBestPair)) {
        bestPairScore = pairScore;
        bestPairAGScore = pairAGScore;
        probabilityOfBestPair = pairProbability;
        for (int r = 0; r < NUM_READS_PER_PAIR; r++) {
            bestResultGenomeLocation[r] = newResult->location[r];
            bestResultOrigGenomeLocation[r] = newResult->origLocation[r];
            bestResultScore[r] = newResult->score[r];
            bestResultDirection[r] = newResult->direction[r];
            bestResultUsedAffineGapScoring[r] = newResult->usedAffineGapScoring[r];
            bestResultBasesClippedBefore[r] = newResult->basesClippedBefore[r];
            bestResultBasesClippedAfter[r] = newResult->basesClippedAfter[r];
            bestResultAGScore[r] = newResult->agScore[r];
            bestResultSeedOffset[r] = newResult->seedOffset[r];
            bestResultMatchProbability[r] = newResult->matchProbability[r];
            bestResultLVIndels[r] = newResult->lvIndels[r];
        }
        return true;
    }
    else {
        return false;
    }
} // updateBestHitIfNeeded


void IntersectingPairedEndAligner::ScoreSet::fillInResult(PairedAlignmentResult* result, unsigned *popularSeedsSkipped) 
{
    for (unsigned whichRead = 0; whichRead < NUM_READS_PER_PAIR; whichRead++) {
        result->location[whichRead] = bestResultGenomeLocation[whichRead];
        result->origLocation[whichRead] = bestResultOrigGenomeLocation[whichRead];
        result->direction[whichRead] = bestResultDirection[whichRead];
        result->mapq[whichRead] = computeMAPQ(probabilityOfAllPairs, probabilityOfBestPair, bestResultScore[whichRead], popularSeedsSkipped[0] + popularSeedsSkipped[1]);
        result->status[whichRead] = result->mapq[whichRead] > MAPQ_LIMIT_FOR_SINGLE_HIT ? SingleHit : MultipleHits;
        result->score[whichRead] = bestResultScore[whichRead];
        result->clippingForReadAdjustment[whichRead] = 0;
        result->usedAffineGapScoring[whichRead] = bestResultUsedAffineGapScoring[whichRead];
        result->basesClippedBefore[whichRead] = bestResultBasesClippedBefore[whichRead];
        result->basesClippedAfter[whichRead] = bestResultBasesClippedAfter[whichRead];
        result->agScore[whichRead] = bestResultAGScore[whichRead];
        result->seedOffset[whichRead] = bestResultSeedOffset[whichRead];
        result->lvIndels[whichRead] = bestResultLVIndels[whichRead];
        result->matchProbability[whichRead] = bestResultMatchProbability[whichRead];
    } // for each read in the pair
    result->probabilityAllPairs = probabilityOfAllPairs;
} // fillInResult

const unsigned IntersectingPairedEndAligner::maxMergeDistance = 31;

int IntersectingPairedEndAligner::computeScoreLimit(bool nonALTAlignment, const ScoreSet * scoresForAllAlignments, const ScoreSet * scoresForNonAltAlignments, GenomeDistance maxBigIndelSeen)
{
    if (nonALTAlignment) {
        //
        // For a non-ALT alignment to matter, it must be no worse than maxScoreGapToPreferNonAltAlignment of the best ALT alignment and at least as good as the best non-ALT alignment.
        //
        return extraSearchDepth + min(maxK + maxBigIndelSeen, min(scoresForAllAlignments->bestPairScore + maxScoreGapToPreferNonAltAlignment, scoresForNonAltAlignments->bestPairScore));
    } else {
        //
        // For an ALT alignment to matter, it has to be at least maxScoreGapToPreferNonAltAlignment better than the best non-ALT alignment, and better than the best ALT alignment.
        //
        return extraSearchDepth + min(maxK + maxBigIndelSeen, min(scoresForAllAlignments->bestPairScore, scoresForNonAltAlignments->bestPairScore - maxScoreGapToPreferNonAltAlignment));
    }
}

#if INSTRUMENTATION_FOR_PAPER
    int
IntersectingPairedEndAligner::HashTableHitSet::getNumDistinctHitLocations(unsigned maxK)
{
    int nextSeed[MAX_MAX_SEEDS];
    int nDistinctHits = 0;

    GenomeLocation currentAnchor = InvalidGenomeLocation;
    int nSetsRemaining = nLookupsUsed;

    for (int i = 0; i < nLookupsUsed; i++) {
        nextSeed[i] = 0;
    }

#define BODY(lookups)                                                                                       \
    while (true) {                                                                                          \
        int setForThisPass = -1;                                                                            \
        GenomeLocation locationForThisPass = 0;                                                             \
        for (int whichSet = 0; whichSet < nLookupsUsed; whichSet++) {                                       \
            if (lookups[whichSet].nHits == nextSeed[whichSet]) {                                            \
                nextSeed[whichSet]++;                                                                       \
                nSetsRemaining--;                                                                           \
            } else if (lookups[whichSet].nHits > nextSeed[whichSet] &&                                      \
                ((GenomeLocation)lookups[whichSet].hits[nextSeed[whichSet]] -                               \
                    lookups[whichSet].seedOffset) > locationForThisPass) {                                  \
                setForThisPass = whichSet;                                                                  \
                locationForThisPass = lookups[whichSet].hits[nextSeed[whichSet]] -                          \
                    lookups[whichSet].seedOffset;                                                           \
            }                                                                                               \
        }                                                                                                   \
                                                                                                            \
        _ASSERT(nSetsRemaining == 0 || setForThisPass != -1);                                               \
        if (nSetsRemaining == 0) {                                                                          \
            return nDistinctHits;                                                                           \
        }                                                                                                   \
                                                                                                            \
        if (currentAnchor - locationForThisPass > maxK) {                                                   \
            nDistinctHits++;                                                                                \
            currentAnchor = locationForThisPass;                                                            \
        }                                                                                                   \
                                                                                                            \
        nextSeed[setForThisPass]++;                                                                         \
    }

    if (doesGenomeIndexHave64BitLocations) {
        BODY(lookups64);
    } else {
        while (true) {            
            int setForThisPass = -1;                                                                            
            GenomeLocation locationForThisPass = 0;                                         
            for (int whichSet = 0; whichSet < nLookupsUsed; whichSet++) {
                     if (lookups32[whichSet].nHits == nextSeed[whichSet]) {                            
                            nSetsRemaining--;         
                            nextSeed[whichSet]++;
                    } else if (lookups32[whichSet].nHits > nextSeed[whichSet] &&
                        ((GenomeLocation)lookups32[whichSet].hits[nextSeed[whichSet]] - lookups32[whichSet].seedOffset) > locationForThisPass) {
                            
                            setForThisPass = whichSet;                                                                  
                            locationForThisPass = lookups32[whichSet].hits[nextSeed[whichSet]] - lookups32[whichSet].seedOffset;
                    }                                                                                               
            }                                                                                                   
                    
            _ASSERT(nSetsRemaining == 0 || setForThisPass != -1);                                               
            if (nSetsRemaining == 0) {                                    
                return nDistinctHits;                                                                           
            }                                                                                                   
                                   
            if (currentAnchor - locationForThisPass > maxK) {
                nDistinctHits++;                                                                                
                currentAnchor = locationForThisPass;
            }                                                                                                   

            nextSeed[setForThisPass]++;                                                                         
        }
    }
#undef BODY

    _ASSERT(!"NOTREACHED");
    return -1;
}
#endif // INSTRUMENTATION_FOR_PAPER
