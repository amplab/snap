/*++

Module Name:

    BaseAligner.cpp

Abstract:

    SNAP genome aligner

Authors:

    Bill Bolosky, August, 2011

Environment:

    User mode service.

    This class is NOT thread safe.  It's the caller's responsibility to ensure that
    at most one thread uses an instance at any time.

Revision History:

    Adapted from Matei Zaharia's Scala implementation.

--*/

#include "stdafx.h"
#include "BaseAligner.h"
#include "Compat.h"
#include "LandauVishkin.h"
#include "BigAlloc.h"
#include "mapq.h"
#include "SeedSequencer.h"
#include "exit.h"
#include "AlignerOptions.h"
#include "Error.h"

using std::min;

#ifdef TRACE_ALIGNER    // If you turn this on, then stdout writing won't work.
#define TRACE printf
#else
#define TRACE(...) {}
#endif

BaseAligner::BaseAligner(
    GenomeIndex    *i_genomeIndex,
    unsigned        i_maxHitsToConsider,
    unsigned        i_maxK,
    unsigned        i_maxReadSize,
    unsigned        i_maxSeedsToUseFromCommandLine,
    double          i_maxSeedCoverage,
    unsigned        i_minWeightToCheck,
    unsigned        i_extraSearchDepth,
    bool            i_noUkkonen,
    bool            i_noOrderedEvaluation,
	bool			i_noTruncation,
    int             i_maxSecondaryAlignmentsPerContig,
    LandauVishkin<1>*i_landauVishkin,
    LandauVishkin<-1>*i_reverseLandauVishkin,
    AlignerStats   *i_stats,
    BigAllocator   *allocator) :
        genomeIndex(i_genomeIndex), maxHitsToConsider(i_maxHitsToConsider), maxK(i_maxK),
        maxReadSize(i_maxReadSize), maxSeedsToUseFromCommandLine(i_maxSeedsToUseFromCommandLine),
        maxSeedCoverage(i_maxSeedCoverage), readId(-1), extraSearchDepth(i_extraSearchDepth),
        explorePopularSeeds(false), stopOnFirstHit(false), stats(i_stats), 
        noUkkonen(i_noUkkonen), noOrderedEvaluation(i_noOrderedEvaluation), noTruncation(i_noTruncation),
		minWeightToCheck(max(1u, i_minWeightToCheck)), maxSecondaryAlignmentsPerContig(i_maxSecondaryAlignmentsPerContig)
/*++

Routine Description:

    Constructor for the BaseAligner class.  Aligners align reads against an indexed genome.

Arguments:

    i_genomeIndex       - The index against which to do the alignments
    i_maxHitsToConsider - The maximum number of hits to use from a seed lookup.  Any lookups that return more
                          than this are ignored.
    i_maxK              - The largest string difference to consider for any comparison.
    i_maxReadSize       - Bound on the number of bases in any read.  There's no reason to make it tight, it just affects a little memory allocation.
    i_maxSeedsToUse     - The maximum number of seeds to use when aligning any read (not counting ones ignored because they resulted in too many
                          hits).  Once we've looked up this many seeds, we just score what we've got.
    i_maxSeedCoverage   - The maximum number of seeds to use expressed as readSize/seedSize
    i_extraSearchDepth  - How deeply beyond bestScore do we search?
    i_noUkkonen         - Don't use Ukkonen's algorithm (i.e., don't reduce the max edit distance depth as we score candidates)
    i_noOrderedEvaluation-Don't order evaluating the reads by the hit count in order to drive down the max edit distance more quickly
	i_noTruncation       - Don't truncate searches based on count of disjoint seed misses
    i_maxSecondaryAlignmentsPerContig - Maximum secondary alignments per contig; -1 means don't limit this
    i_landauVishkin     - an externally supplied LandauVishkin string edit distance object.  This is useful if we're expecting repeated computations and use the LV cache.
    i_reverseLandauVishkin - the same for the reverse direction.
    i_stats             - an object into which we report out statistics
    allocator           - an allocator that's used to allocate our local memory.  This is useful for TLB optimization.  If this is supplied, the caller
                          is responsible for deallocation, we'll not deallocate any dynamic memory in our destructor.

 --*/
{
    hadBigAllocator = allocator != NULL;

    nHashTableLookups = 0;
    nLocationsScored = 0;
    nHitsIgnoredBecauseOfTooHighPopularity = 0;
    nReadsIgnoredBecauseOfTooManyNs = 0;
    nIndelsMerged = 0;

    genome = genomeIndex->getGenome();
    seedLen = genomeIndex->getSeedLength();
    doesGenomeIndexHave64BitLocations = genomeIndex->doesGenomeIndexHave64BitLocations();

    probDistance = new ProbabilityDistance(SNP_PROB, GAP_OPEN_PROB, GAP_EXTEND_PROB);  // Match Mason

    if ((i_landauVishkin == NULL) != (i_reverseLandauVishkin == NULL)) {
        WriteErrorMessage("Must supply both or neither of forward & reverse Landau-Vishkin objects.  You tried exactly one.\n");
        soft_exit(1);
    }

    if (i_landauVishkin == NULL) {
        if (allocator) {
            landauVishkin = new (allocator) LandauVishkin<>;
            reverseLandauVishkin = new (allocator) LandauVishkin<-1>;
        } else {
            landauVishkin = new LandauVishkin<>;
            reverseLandauVishkin = new LandauVishkin<-1>;
        }
        ownLandauVishkin = true;
    } else {
        landauVishkin = i_landauVishkin;
        reverseLandauVishkin = i_reverseLandauVishkin;
        ownLandauVishkin = false;
    }

    unsigned maxSeedsToUse;
    if (0 != maxSeedsToUseFromCommandLine) {
        maxSeedsToUse = maxSeedsToUseFromCommandLine;
    } else {
        maxSeedsToUse = (int)(maxSeedCoverage * maxReadSize / genomeIndex->getSeedLength());
    }

    numWeightLists = maxSeedsToUse + 1;

    candidateHashTablesSize = (maxHitsToConsider * maxSeedsToUse * 3)/2;    // *1.5 for hash table slack
    hashTableElementPoolSize = maxHitsToConsider * maxSeedsToUse * 2 ;   // *2 for RC

    if (allocator) {
        rcReadData = (char *)allocator->allocate(sizeof(char) * maxReadSize * 2); // The *2 is to allocte space for the quality string
    } else {
        rcReadData = (char *)BigAlloc(sizeof(char) * maxReadSize * 2); // The *2 is to allocte space for the quality string
    }

    rcReadQuality = rcReadData + maxReadSize;

    if (allocator) {
        reversedRead[FORWARD] = (char *)allocator->allocate(sizeof(char) * maxReadSize * 4 + 2 * MAX_K); // Times 4 to also hold RC version and genome data (+2MAX_K is for genome data)
    } else {
        reversedRead[FORWARD] = (char *)BigAlloc(sizeof(char) * maxReadSize * 4 + 2 * MAX_K); // Times 4 to also hold RC version and genome data (+2MAX_K is for genome data)
    }

    rcReadData = (char *)BigAlloc(sizeof(char) * maxReadSize);

    // treat everything but ACTG like N
    for (unsigned i = 0; i < 256; i++) {
        nTable[i] = 1;
        rcTranslationTable[i] = 'N';
    }
    reversedRead[RC] = reversedRead[FORWARD] + maxReadSize;

    rcTranslationTable['A'] = 'T';
    rcTranslationTable['G'] = 'C';
    rcTranslationTable['C'] = 'G';
    rcTranslationTable['T'] = 'A';
    rcTranslationTable['N'] = 'N';

    memset(nTable, 0, sizeof(nTable));

    nTable['N'] = 1;

    if (allocator) {
        seedUsed = (BYTE *)allocator->allocate((sizeof(BYTE) * (maxReadSize + 7 + 128) / 8));    // +128 to make sure it extends at both
    } else {
        seedUsed = (BYTE *)BigAlloc((sizeof(BYTE) * (maxReadSize + 7 + 128) / 8));    // +128 to make sure it extends at both
    }

    seedUsedAsAllocated = seedUsed; // Save the pointer for the delete.
    seedUsed += 8;  // This moves the pointer up an _int64, so we now have the appropriate before buffer.

    nUsedHashTableElements = 0;

    if (allocator) {
        candidateHashTable[FORWARD] = (HashTableAnchor *)allocator->allocate(sizeof(HashTableAnchor) * candidateHashTablesSize);
        candidateHashTable[RC] = (HashTableAnchor *)allocator->allocate(sizeof(HashTableAnchor) * candidateHashTablesSize);
        weightLists = (HashTableElement *)allocator->allocate(sizeof(HashTableElement) * numWeightLists);
        hashTableElementPool = (HashTableElement *)allocator->allocate(sizeof(HashTableElement) * hashTableElementPoolSize); // Allocte last, because it's biggest and usually unused.  This puts all of the commonly used stuff into one large page.
        hitCountByExtraSearchDepth = (unsigned *)allocator->allocate(sizeof(*hitCountByExtraSearchDepth) * extraSearchDepth);
        if (maxSecondaryAlignmentsPerContig > 0) {
            hitsPerContigCounts = (HitsPerContigCounts *)allocator->allocate(sizeof(*hitsPerContigCounts) * genome->getNumContigs());
            memset(hitsPerContigCounts, 0, sizeof(*hitsPerContigCounts) * genome->getNumContigs());
        } else {
            hitsPerContigCounts = NULL;
        }
    } else {
        candidateHashTable[FORWARD] = (HashTableAnchor *)BigAlloc(sizeof(HashTableAnchor) * candidateHashTablesSize);
        candidateHashTable[RC] = (HashTableAnchor *)BigAlloc(sizeof(HashTableAnchor) * candidateHashTablesSize);
        weightLists = (HashTableElement *)BigAlloc(sizeof(HashTableElement) * numWeightLists);
        hashTableElementPool = (HashTableElement *)BigAlloc(sizeof(HashTableElement) * hashTableElementPoolSize);
        hitCountByExtraSearchDepth = (unsigned *)BigAlloc(sizeof(*hitCountByExtraSearchDepth) * extraSearchDepth);
        if (maxSecondaryAlignmentsPerContig > 0) {
            hitsPerContigCounts = (HitsPerContigCounts *)BigAlloc(sizeof(*hitsPerContigCounts) * genome->getNumContigs());
            memset(hitsPerContigCounts, 0, sizeof(*hitsPerContigCounts) * genome->getNumContigs());
        }
        else {
            hitsPerContigCounts = NULL;
        }
    }

    for (unsigned i = 0; i < hashTableElementPoolSize; i++) {
        hashTableElementPool[i].init();
    }

    for (unsigned i = 0; i < maxSeedsToUse + 1; i++) {
        weightLists[i].init();
    }

    for (Direction rc = 0; rc < NUM_DIRECTIONS; rc++) {
        memset(candidateHashTable[rc],0,sizeof(HashTableAnchor) * candidateHashTablesSize);
    }
    hashTableEpoch = 0;

 
}


#ifdef  _DEBUG
bool _DumpAlignments = false;
#endif  // _DEBUG

    void
BaseAligner::AlignRead(
        Read                    *inputRead,
        SingleAlignmentResult   *primaryResult,
        int                      maxEditDistanceForSecondaryResults,
        int                      secondaryResultBufferSize,
        int                     *nSecondaryResults,
        int                      maxSecondaryResults,
        SingleAlignmentResult   *secondaryResults             // The caller passes in a buffer of secondaryResultBufferSize and it's filled in by AlignRead()
    )
/*++

Routine Description:

    Align a particular read, possibly constraining the search around a given location.

Arguments:

    read                                - the read to align
    primaryResult                       - the best alignment result found
    maxEditDistanceForSecondaryResults  - How much worse than the primary result should we look?
    secondaryResultBufferSize           - the size of the secondaryResults buffer.  If provided, it must be at least maxK * maxSeeds * 2.
    nRescondaryResults                  - returns the number of secondary results found
    maxSecondaryResults                 - limit the number of secondary results to this
    secondaryResults                    - returns the secondary results


Return Value:

    true if there was enough space in secondaryResults, false otherwise

--*/
{   
    memset(hitCountByExtraSearchDepth, 0, sizeof(*hitCountByExtraSearchDepth) * extraSearchDepth);

    if (NULL != nSecondaryResults) {
        *nSecondaryResults = 0;
    }

    firstPassSeedsNotSkipped[FORWARD] = firstPassSeedsNotSkipped[RC] = 0;
    smallestSkippedSeed[FORWARD] = smallestSkippedSeed[RC] = 0x8fffffffffffffff;
    highestWeightListChecked = 0;

    unsigned maxSeedsToUse;
    if (0 != maxSeedsToUseFromCommandLine) {
        maxSeedsToUse = maxSeedsToUseFromCommandLine;
    } else {
        maxSeedsToUse = (int)(2 * maxSeedCoverage * inputRead->getDataLength() / genomeIndex->getSeedLength()); // 2x is for FORWARD/RC
    }

    primaryResult->location = InvalidGenomeLocation; // Value to return if we don't find a location.
    primaryResult->direction = FORWARD;              // So we deterministically print the read forward in this case.
    primaryResult->score = UnusedScoreValue;
    primaryResult->status = NotFound;

    unsigned lookupsThisRun = 0;

    popularSeedsSkipped = 0;

    //
    // A bitvector for used seeds, indexed on the starting location of the seed within the read.
    //
    if (inputRead->getDataLength() > maxReadSize) {
        WriteErrorMessage("BaseAligner:: got too big read (%d > %d)\n" 
                          "Increase MAX_READ_LENGTH at the beginning of Read.h and recompile\n", inputRead->getDataLength(), maxReadSize);
        soft_exit(1);
    }

    if ((int)inputRead->getDataLength() < seedLen) {
        //
        // Too short to have any seeds, it's hopeless.
        // No need to finalize secondary results, since we don't have any.
        //
        return;
    }

#ifdef TRACE_ALIGNER
    printf("Aligning read '%.*s':\n%.*s\n%.*s\n", inputRead->getIdLength(), inputRead->getId(), inputRead->getDataLength(), inputRead->getData(),
            inputRead->getDataLength(), inputRead->getQuality());
#endif

#ifdef  _DEBUG
    if (_DumpAlignments) {
        printf("BaseAligner: aligning read ID '%.*s', data '%.*s'\n", inputRead->getIdLength(), inputRead->getId(), inputRead->getDataLength(), inputRead->getData());
    }
#endif  // _DEBUG

    //
    // Clear out the seed used array.
    //
    memset(seedUsed, 0, (inputRead->getDataLength() + 7) / 8);

    unsigned readLen = inputRead->getDataLength();
    const char *readData = inputRead->getData();
    const char *readQuality = inputRead->getQuality();
    unsigned countOfNs = 0;
    for (unsigned i = 0; i < readLen; i++) {
        char baseByte = readData[i];
        char complement = rcTranslationTable[baseByte];
        rcReadData[readLen - i - 1] = complement;
        rcReadQuality[readLen - i - 1] = readQuality[i];
        reversedRead[FORWARD][readLen - i - 1] = baseByte;
        reversedRead[RC][i] = complement;
        countOfNs += nTable[baseByte];
    }

    if (countOfNs > maxK) {
        nReadsIgnoredBecauseOfTooManyNs++;
        // No need to finalize secondary results, since we don't have any.
        return;
    }

    //
    // Block off any seeds that would contain an N.
    //
    if (countOfNs > 0) {
        int minSeedToConsiderNing = 0; // In English, any word can be verbed. Including, apparently, "N."
        for (int i = 0; i < (int) readLen; i++) {
            if (BASE_VALUE[readData[i]] > 3) {
                int limit = __min(i + seedLen - 1, readLen-1);
                for (int j = __max(minSeedToConsiderNing, i - (int) seedLen + 1); j <= limit; j++) {
                    SetSeedUsed(j);
                }
                minSeedToConsiderNing = limit+1;
                if (minSeedToConsiderNing >= (int) readLen)
                    break;
            }
        }
    }

    Read reverseComplimentRead;
    Read *read[NUM_DIRECTIONS];
    read[FORWARD] = inputRead;
    read[RC] = &reverseComplimentRead;
    read[RC]->init(NULL, 0, rcReadData, rcReadQuality, readLen);

    clearCandidates();

    //
    // Initialize the bases table, which represents which bases we've checked.
    // We have readSize - seeds size + 1 possible seeds.
    //
    unsigned nPossibleSeeds = readLen - seedLen + 1;
    TRACE("nPossibleSeeds: %d\n", nPossibleSeeds);

    unsigned nextSeedToTest = 0;
    unsigned wrapCount = 0;
    lowestPossibleScoreOfAnyUnseenLocation[FORWARD] = lowestPossibleScoreOfAnyUnseenLocation[RC] = 0;
    mostSeedsContainingAnyParticularBase[FORWARD] = mostSeedsContainingAnyParticularBase[RC] = 1;  // Instead of tracking this for real, we're just conservative and use wrapCount+1.  It's faster.
    bestScore = UnusedScoreValue;
    secondBestScore = UnusedScoreValue;
    nSeedsApplied[FORWARD] = nSeedsApplied[RC] = 0;
    lvScores = 0;
    lvScoresAfterBestFound = 0;
    probabilityOfAllCandidates = 0.0;
    probabilityOfBestCandidate = 0.0;

    scoreLimit = maxK + extraSearchDepth; // For MAPQ computation

    while (nSeedsApplied[FORWARD] + nSeedsApplied[RC] < maxSeedsToUse) {
        //
        // Choose the next seed to use.  Choose the first one that isn't used
        //
        if (nextSeedToTest >= nPossibleSeeds) {
            //
            // We're wrapping.  We want to space the seeds out as much as possible, so if we had
            // a seed length of 20 we'd want to take 0, 10, 5, 15, 2, 7, 12, 17.  To make the computation
            // fast, we use use a table lookup.
            //
            wrapCount++;
            if (wrapCount >= seedLen) {
                //
                // We tried all possible seeds without matching or even getting enough seeds to
                // exceed our seed count.  Do the best we can with what we have.
                //
#ifdef TRACE_ALIGNER
                printf(stderr, "Calling score with force=true because we wrapped around enough\n");
#endif
                score(
                    true,
                    read,
                    primaryResult,
                    maxEditDistanceForSecondaryResults,
                    secondaryResultBufferSize,
                    nSecondaryResults,
                    secondaryResults);

#ifdef  _DEBUG
                if (_DumpAlignments) printf("\tFinal result score %d MAPQ %d (%e probability of best candidate, %e probability of all candidates)  at %u\n", 
                                            primaryResult->score, primaryResult->mapq, probabilityOfBestCandidate, probabilityOfAllCandidates, primaryResult->location);
#endif  // _DEBUG
                finalizeSecondaryResults(*primaryResult, nSecondaryResults, secondaryResults, maxSecondaryResults, maxEditDistanceForSecondaryResults, bestScore);
                return;
            }
            nextSeedToTest = GetWrappedNextSeedToTest(seedLen, wrapCount);

            mostSeedsContainingAnyParticularBase[FORWARD] = mostSeedsContainingAnyParticularBase[RC] = wrapCount + 1;
        }

        while (nextSeedToTest < nPossibleSeeds && IsSeedUsed(nextSeedToTest)) {
            //
            // This seed is already used.  Try the next one.
            //
            TRACE("Skipping due to IsSeedUsed\n");
            nextSeedToTest++;
        }

        if (nextSeedToTest >= nPossibleSeeds) {
            //
            // Unusable seeds have pushed us past the end of the read.  Go back around the outer loop so we wrap properly.
            //
            TRACE("Eek, we're past the end of the read\n");
            continue;
        }

        SetSeedUsed(nextSeedToTest);

        if (!Seed::DoesTextRepresentASeed(read[FORWARD]->getData() + nextSeedToTest, seedLen)) {
            continue;
        }

        Seed seed(read[FORWARD]->getData() + nextSeedToTest, seedLen);

        _int64        nHits[NUM_DIRECTIONS];                // Number of times this seed hits in the genome
        const GenomeLocation  *hits[NUM_DIRECTIONS];        // The actual hits (of size nHits)
        GenomeLocation singletonHits[NUM_DIRECTIONS];       // Storage for single hits (this is required for 64 bit genome indices, since they might use fewer than 8 bytes internally)

        const unsigned *hits32[NUM_DIRECTIONS];

        if (doesGenomeIndexHave64BitLocations) {
            genomeIndex->lookupSeed(seed, &nHits[FORWARD], &hits[FORWARD], &nHits[RC], &hits[RC], &singletonHits[FORWARD], &singletonHits[RC]);
        } else {
            genomeIndex->lookupSeed32(seed, &nHits[FORWARD], &hits32[FORWARD], &nHits[RC], &hits32[RC]);
        }

        nHashTableLookups++;
        lookupsThisRun++;


#ifdef  _DEBUG
        if (_DumpAlignments) {
            printf("\tSeed offset %2d, %4d hits, %4d rcHits.", nextSeedToTest, nHits[0], nHits[1]);
            for (int rc = 0; rc < 2; rc++) {
                for (unsigned i = 0; i < __min(nHits[rc], 5); i++) {
                    printf(" %sHit at %9llu.", rc == 1 ? "RC " : "", doesGenomeIndexHave64BitLocations ? hits[rc][i] : (_int64)hits32[rc][i]);
                }
            }
            printf("\n");
        }
#endif  // _DEUBG

#ifdef TRACE_ALIGNER
        printf("Looked up seed %.*s (offset %d): hits=%u, rchits=%u\n",
                seedLen, inputRead->getData() + nextSeedToTest, nextSeedToTest, nHits[0], nHits[1]);
        for (int rc = 0; rc < 2; rc++) {
            if (nHits[rc] <= maxHitsToConsider) {
                printf("%sHits:", rc == 1 ? "RC " : "");
                for (unsigned i = 0; i < nHits[rc]; i++)
                    printf(" %u", hits[rc][i]);
                printf("\n");
            }
        }
#endif

        bool appliedEitherSeed = false;

        for (Direction direction = 0; direction < NUM_DIRECTIONS; direction++) {
            if (nHits[direction] > maxHitsToConsider && !explorePopularSeeds) {
                //
                // This seed is matching too many places.  Just pretend we never looked and keep going.
                //
                nHitsIgnoredBecauseOfTooHighPopularity++;
                popularSeedsSkipped++;
                smallestSkippedSeed[direction] = __min(nHits[direction], smallestSkippedSeed[direction]);
            } else {
                if (0 == wrapCount) {
                    firstPassSeedsNotSkipped[direction]++;
                }

                //
                // Update the candidates list with any hits from this seed.  If lowest possible score of any unseen location is
                // more than best_score + confDiff then we know that if this location is newly seen then its location won't ever be a
                // winner, and we can ignore it.
                //

                unsigned offset;
                if (direction == FORWARD) {
                    offset = nextSeedToTest;
                } else {
                    //
                    // The RC seed is at offset ReadSize - SeedSize - seed offset in the RC seed.
                    //
                    // To see why, imagine that you had a read that looked like 0123456 (where the digits
                    // represented some particular bases, and digit' is the base's complement). Then the
                    // RC of that read is 6'5'4'3'2'1'.  So, when we look up the hits for the seed at
                    // offset 0 in the forward read (i.e. 012 assuming a seed size of 3) then the index
                    // will also return the results for the seed's reverse complement, i.e., 3'2'1'.
                    // This happens as the last seed in the RC read.
                    //
                    offset = readLen - seedLen - nextSeedToTest;
                }

                const unsigned prefetchDepth = 30;
                _int64 limit = min(nHits[direction], (_int64)maxHitsToConsider) + prefetchDepth;
                for (unsigned iBase = 0 ; iBase < limit; iBase += prefetchDepth) {
                    //
                    // This works in two phases: we launch prefetches for a group of hash table lines,
                    // then we do all of the inserts, and then repeat.
                    //

		            _int64 innerLimit = min((_int64)iBase + prefetchDepth, min(nHits[direction], (_int64)maxHitsToConsider));
                    if (doAlignerPrefetch) {
                        for (unsigned i = iBase; i < innerLimit; i++) {
                            if (doesGenomeIndexHave64BitLocations) {
                                prefetchHashTableBucket(GenomeLocationAsInt64(hits[direction][i]) - offset, direction);
                            } else {
                                prefetchHashTableBucket(hits32[direction][i] - offset, direction);
                            }
                        }
                    }

                    for (unsigned i = iBase; i < innerLimit; i++) {
                        //
                        // Find the genome location where the beginning of the read would hit, given a match on this seed.
                        //
                        GenomeLocation genomeLocationOfThisHit;
                        if (doesGenomeIndexHave64BitLocations) {
                            genomeLocationOfThisHit = hits[direction][i] - offset;
                        } else {
                            genomeLocationOfThisHit = hits32[direction][i] - offset;
                        }

                        Candidate *candidate = NULL;
                        HashTableElement *hashTableElement;

                        findCandidate(genomeLocationOfThisHit, direction, &candidate, &hashTableElement);

                        if (NULL != hashTableElement) {
                            if (!noOrderedEvaluation) {     // If noOrderedEvaluation, just leave them all on the one-hit weight list so they get evaluated in whatever order
                                incrementWeight(hashTableElement);
                            }
                            candidate->seedOffset = offset;
                            _ASSERT((unsigned)candidate->seedOffset <= readLen - seedLen);
                        } else if (lowestPossibleScoreOfAnyUnseenLocation[direction] <= scoreLimit || noTruncation) {
                            _ASSERT(offset <= readLen - seedLen);
                            allocateNewCandidate(genomeLocationOfThisHit, direction, lowestPossibleScoreOfAnyUnseenLocation[direction],
                                    offset, &candidate, &hashTableElement);
                        }
                    }
                }
                nSeedsApplied[direction]++;
                appliedEitherSeed = true;
            } // not too popular
        }   // directions

#if 1
        nextSeedToTest += seedLen;
#else   // 0

        //
        // If we don't have enough seeds left to reach the end of the read, space out the seeds more-or-less evenly.
        //
        if ((maxSeedsToUse - (nSeedsApplied[FORWARD] + nSeedsApplied[RC]) + 1) * seedLen + nextSeedToTest < nPossibleSeeds) {
            _ASSERT((nPossibleSeeds + nextSeedToTest) / (maxSeedsToUse - (nSeedsApplied[FORWARD] + nSeedsApplied[RC]) + 1) > seedLen);
            nextSeedToTest += (nPossibleSeeds + nextSeedToTest) / (maxSeedsToUse - (nSeedsApplied[FORWARD] + nSeedsApplied[RC]) + 1);
        } else {
            nextSeedToTest += seedLen;
        }

#endif // 0


        if (appliedEitherSeed) {
            //
            // And finally, try scoring.
            //
            if (score(
                    false,
                    read,
                    primaryResult,
                    maxEditDistanceForSecondaryResults,
                    secondaryResultBufferSize,
                    nSecondaryResults,
                    secondaryResults)) {

#ifdef  _DEBUG
                if (_DumpAlignments) printf("\tFinal result score %d MAPQ %d at %u\n", primaryResult->score, primaryResult->mapq, primaryResult->location);
#endif  // _DEBUG

                finalizeSecondaryResults(*primaryResult, nSecondaryResults, secondaryResults, maxSecondaryResults, maxEditDistanceForSecondaryResults, bestScore);
                return;
            }
        }
    }

    //
    // Do the best with what we've got.
    //
#ifdef TRACE_ALIGNER
    printf("Calling score with force=true because we ran out of seeds\n");
#endif
    score(
        true,
        read,
        primaryResult,
        maxEditDistanceForSecondaryResults,
        secondaryResultBufferSize,
        nSecondaryResults,
        secondaryResults);

#ifdef  _DEBUG
    if (_DumpAlignments) printf("\tFinal result score %d MAPQ %d (%e probability of best candidate, %e probability of all candidates) at %u\n", primaryResult->score, primaryResult->mapq, probabilityOfBestCandidate, probabilityOfAllCandidates, primaryResult->location);
#endif  // _DEBUG

    finalizeSecondaryResults(*primaryResult, nSecondaryResults, secondaryResults, maxSecondaryResults, maxEditDistanceForSecondaryResults, bestScore);
    return;
}

    bool
BaseAligner::score(
        bool                     forceResult,
        Read                    *read[NUM_DIRECTIONS],
        SingleAlignmentResult   *primaryResult,
        int                      maxEditDistanceForSecondaryResults,
        int                      secondaryResultBufferSize,
        int                     *nSecondaryResults,
        SingleAlignmentResult   *secondaryResults)
/*++

Routine Description:

    Make progress in scoring a possibly partial alignment.  This is a private method of the BaseAligner class that's used
    only by AlignRead.

    It does a number of things.  First, it computes the lowest possible score of any unseen location.  This is useful
    because once we have a scored hit that's more than confDiff better than all unseen locations, there's no need to
    lookup more of them, we can just score what we've got and be sure that the answer is right (unless errors have
    pushed the read to be closer to a wrong location than to the correct one, in which case it's hopeless).

    It then decides whether it should score a location, and if so what one to score.  It chooses the unscored
    location that's got the highest weight (i.e., appeared in the most hash table lookups), since that's most
    likely to be the winner.  If there are multiple candidates with the same best weight, it breaks the tie using
    the best possible score for the candidates (which is determined when they're first hit).  Remaining ties are
    broken arbitrarily.

    It merges indels with scored candidates.  If there's an insertion or deletion in the read, then we'll get very
    close but unequal results out of the hash table lookup for parts of the read on opposite sides of the
    insertion or deletion.  This throws out the one with the worse score.

    It then figures out if we have a definitive answer, and says what that is.

Arguments:

    forceResult                             - should we generate an answer even if it's not definitive?
    read                                    - the read we're aligning in both directions
    result                                  - returns the result if we reach one
    singleHitGenomeLocation                 - returns the location in the genome if we return a single hit
    hitDirection                            - if we return a single hit, indicates its direction
    candidates                              - in/out the array of candidates that have hit and possibly been scored
    mapq                                    - returns the map quality if we've reached a final result
    secondary                               - returns secondary alignment locations & directions (optional)

Return Value:

    true iff we've reached a result.  When called with forceResult, we'll always return true.

--*/
{
#ifdef TRACE_ALIGNER
    printf("score() called with force=%d nsa=%d nrcsa=%d best=%u bestloc=%u 2nd=%u\n",
        forceResult, nSeedsApplied[FORWARD], nSeedsApplied[RC], bestScore, bestScoreGenomeLocation, secondBestScore);
    //printf("Candidates:\n");
    //for (int i = 0; i < nCandidates; i++) {
    //    Candidate* c = candidates + i;
    //    printf("  loc=%u rc=%d weight=%u minps=%u scored=%d score=%u r=%u-%u\n",
    //        c->genomeLocation, c->isRC, c->weight, c->minPossibleScore, c->scored,
    //        c->score, c->minRange, c->maxRange);
    //}
    //printf("\n\n");
#endif

    if (0 == mostSeedsContainingAnyParticularBase[FORWARD] && 0 == mostSeedsContainingAnyParticularBase[RC]) {
        //
        // The only way we can get here is if we've tried all of the seeds that we're willing
        // to try and every one of them generated too many hits to process.  Give up.
        //
        _ASSERT(forceResult);
        primaryResult->status = NotFound;
        primaryResult->mapq = 0;
        return true;
    }

    //
    // Recompute lowestPossibleScore.
    //
    for (Direction direction = 0; direction < 2; direction++) {
        if (0 != mostSeedsContainingAnyParticularBase[direction]) {
            lowestPossibleScoreOfAnyUnseenLocation[direction] =
                __max(lowestPossibleScoreOfAnyUnseenLocation[direction],
                      nSeedsApplied[direction] / mostSeedsContainingAnyParticularBase[direction]);
        }
    }

#ifdef TRACE_ALIGNER
    printf("Lowest possible scores for unseen locations: %d (fwd), %d (RC)\n",
        lowestPossibleScoreOfAnyUnseenLocation[FORWARD],
        lowestPossibleScoreOfAnyUnseenLocation[RC]);
#endif

    unsigned weightListToCheck = highestUsedWeightList;

    do {
        //
        // Grab the next element to score, and score it.
        //

        while (weightListToCheck > 0 && weightLists[weightListToCheck].weightNext == &weightLists[weightListToCheck]) {
            weightListToCheck--;
            highestUsedWeightList = weightListToCheck;
        }

        if ((__min(lowestPossibleScoreOfAnyUnseenLocation[FORWARD],lowestPossibleScoreOfAnyUnseenLocation[RC]) > scoreLimit && !noTruncation) || forceResult) {
            if (weightListToCheck< minWeightToCheck) {
                //
                // We've scored all live candidates and excluded all non-candidates, or we've checked enough that we've hit the cutoff.  We have our
                // answer.
                //
                primaryResult->score = bestScore;
                if (bestScore <= maxK) {
                    primaryResult->location = bestScoreGenomeLocation;
                    primaryResult->mapq = computeMAPQ(probabilityOfAllCandidates, probabilityOfBestCandidate, bestScore, popularSeedsSkipped);
                    if (primaryResult->mapq >= MAPQ_LIMIT_FOR_SINGLE_HIT) {
                        primaryResult->status = SingleHit;
                    } else {
                        primaryResult->status = MultipleHits;
                    }
                    return true;
                } else {
                    primaryResult->status = NotFound;
                    primaryResult->mapq = 0;
                    return true;
                }
            }
            //
            // Nothing that we haven't already looked up can possibly be the answer.  Score what we've got and exit.
            //
            forceResult = true;
        } else if (weightListToCheck == 0) {
            //
            // No candidates, look for more.
            //
            return false;
        }

        HashTableElement *elementToScore = weightLists[weightListToCheck].weightNext;
        _ASSERT(!elementToScore->allExtantCandidatesScored);
        _ASSERT(elementToScore->candidatesUsed != 0);
        _ASSERT(elementToScore != &weightLists[weightListToCheck]);

        if (doAlignerPrefetch) {
            //
            // Our prefetch pipeline is one loop out we get the genome data for the next loop, and two loops out we get the element to score.
            //
            _mm_prefetch((const char *)(elementToScore->weightNext->weightNext), _MM_HINT_T2);   // prefetch the next element, it's likely to be the next thing we score.
            genome->prefetchData(elementToScore->weightNext->baseGenomeLocation);
        }

        if (elementToScore->lowestPossibleScore <= scoreLimit) {

            unsigned long candidateIndexToScore;
            _uint64 candidatesMask = elementToScore->candidatesUsed;
            while (_BitScanForward64(&candidateIndexToScore,candidatesMask)) {
                _uint64 candidateBit = ((_uint64)1 << candidateIndexToScore);
                candidatesMask &= ~candidateBit;
                if ((elementToScore->candidatesScored & candidateBit) != 0) {
                    // Already scored it, or marked it as scored due to using ProbabilityDistance
                    continue;
                }

                bool anyNearbyCandidatesAlreadyScored = elementToScore->candidatesScored != 0;

                elementToScore->candidatesScored |= candidateBit;
                _ASSERT(candidateIndexToScore < hashTableElementSize);
                Candidate *candidateToScore = &elementToScore->candidates[candidateIndexToScore];

                GenomeLocation genomeLocation = elementToScore->baseGenomeLocation + candidateIndexToScore;
                GenomeLocation elementGenomeLocation = genomeLocation;    // This is the genome location prior to any adjustments for indels

                //
                // We're about to run edit distance computation on the genome.  Launch a prefetch for it
                // so that it's in cache when we do (or at least on the way).
                //
                if (doAlignerPrefetch) {
                    genomeIndex->prefetchGenomeData(genomeLocation);
                }

                unsigned score = -1;
                double matchProbability = 0;
                unsigned readDataLength = read[elementToScore->direction]->getDataLength();
                GenomeDistance genomeDataLength = readDataLength + MAX_K; // Leave extra space in case the read has deletions
                const char *data = genome->getSubstring(genomeLocation, genomeDataLength);

#if 0 // This only happens when we're in the padding region, and genomeLocations there just lead to problems.  Just say no.
                if (NULL == data) {
                    //
                    // We're up against the end of a chromosome.  Reduce the extra space enough that it isn't too
                    // long.  We're willing to reduce it to less than the length of a read, because the read could
                    // butt up against the end of the chromosome and have insertions in it.
                    //
                    const Genome::Contig *contig = genome->getContigAtLocation(genomeLocation);

                    if (contig != NULL) {
                        GenomeLocation endLocation;
                        if (genomeLocation + readDataLength + MAX_K >= GenomeLocation(0) + genome->getCountOfBases()) {
                            endLocation = GenomeLocation(0) + genome->getCountOfBases();
                        } else {
                            const Genome::Contig *nextContig = genome->getNextContigAfterLocation(genomeLocation);
                            _ASSERT(contig->beginningLocation <= genomeLocation && contig != nextContig);

                            endLocation = nextContig->beginningLocation;
                        }
                        genomeDataLength = endLocation - genomeLocation - 1;
                        if (genomeDataLength >= readDataLength - MAX_K) {
                            data = genome->getSubstring(genomeLocation, genomeDataLength);
                            _ASSERT(NULL != data);
                        }
                    }
                }

#endif // 0
                if (data != NULL) {
                    Read *readToScore = read[elementToScore->direction];

                    _ASSERT(candidateToScore->seedOffset + seedLen <= readToScore->getDataLength());

                    //
                    // Compute the distance separately in the forward and backward directions from the seed, to allow
                    // arbitrary offsets at both the start and end.
                    //
                    double matchProb1, matchProb2;
                    int score1, score2;
                    // First, do the forward direction from where the seed aligns to past of it
                    int readLen = readToScore->getDataLength();
                    int seedLen = genomeIndex->getSeedLength();
                    int seedOffset = candidateToScore->seedOffset; // Since the data is reversed
                    int tailStart = seedOffset + seedLen;

                    _ASSERT(!memcmp(data+seedOffset, readToScore->getData() + seedOffset, seedLen));

                    int textLen = (int)__min(genomeDataLength - tailStart, 0x7ffffff0);
                    score1 = landauVishkin->computeEditDistance(data + tailStart, textLen, readToScore->getData() + tailStart, readToScore->getQuality() + tailStart, readLen - tailStart,
                        scoreLimit, &matchProb1);

                    if (score1 == -1) {
                        score = -1;
                    } else {
                        // The tail of the read matched; now let's reverse match the reference genome and the head
                        int limitLeft = scoreLimit - score1;
                        int genomeLocationOffset;
                        score2 = reverseLandauVishkin->computeEditDistance(data + seedOffset, seedOffset + MAX_K, reversedRead[elementToScore->direction] + readLen - seedOffset,
                                                                                    read[OppositeDirection(elementToScore->direction)]->getQuality() + readLen - seedOffset, seedOffset, limitLeft, &matchProb2,
                                                                                    &genomeLocationOffset);

                        if (score2 == -1) {
                            score = -1;
                        } else {
                            score = score1 + score2;
                            // Map probabilities for substrings can be multiplied, but make sure to count seed too
                            matchProbability = matchProb1 * matchProb2 * pow(1 - SNP_PROB, seedLen);

                            //
                            // Adjust the genome location based on any indels that we found.
                            //
                            genomeLocation += genomeLocationOffset;

                            //
                            // We could mark as scored anything in between the old and new genome offsets, but it's probably not worth the effort since this is
                            // so rare and all it would do is same time.
                            //
                        }
                    }
                } else { // if we had genome data to compare against
                    matchProbability = 0;
                }
#ifdef TRACE_ALIGNER
                printf("Computing distance at %u (RC) with limit %d: %d (prob %g)\n",
                        genomeLocation, scoreLimit, score, matchProbability);
#endif


#ifdef  _DEBUG
                if (_DumpAlignments) printf("Scored %9u weight %2d limit %d, result %2d %s\n", genomeLocation, elementToScore->weight, scoreLimit, score, elementToScore->direction ? "RC" : "");
#endif  // _DEBUG

                candidateToScore->score = score;

                nLocationsScored++;
                lvScores++;
                lvScoresAfterBestFound++;

                //
                // Handle the special case where we just scored a different offset for a region that's already been scored.  This can happen when
                // there are indels in the read.  In this case, we want to treat them as a single aignment, not two different ones (which would
                // cause us to lose confidence in the alignment, since they're probably both pretty good).
                //
                if (anyNearbyCandidatesAlreadyScored) {
                    if (elementToScore->bestScore < score || elementToScore->bestScore == score && matchProbability <= elementToScore->matchProbabilityForBestScore) {
                        //
                        // This is a no better mapping than something nearby that we already tried.  Just ignore it.
                        //
                        continue;
                    }
                } else {
                    _ASSERT(elementToScore->matchProbabilityForBestScore == 0.0);
                }

                elementToScore->bestScoreGenomeLocation = genomeLocation;

                //
                // Look up the hash table element that's closest to the genomeLocation but that doesn't
                // contain it, to check if this location is already scored.
                //
                // We do this computation in a strange way in order to avoid generating a branch instruction that
                // the processor's branch predictor will get wrong half of the time.  Think about it like this:
                // The genome location lies in a bucket of size hashTableElementSize.  Its offset in the bucket
                // is genomeLocation % hashTableElementSize.  If we take that quantity and integer divide it by
                // hashTableElementSize / 2, we get 0 if it's in the first half and 1 if it's in the second.  Double that and subtract
                // one, and you're at the right place with no branches.
                //
                HashTableElement *nearbyElement;
                GenomeLocation nearbyGenomeLocation;
                if (-1 != score) {
                    nearbyGenomeLocation = elementGenomeLocation + (2*(GenomeLocationAsInt64(elementGenomeLocation) % hashTableElementSize / (hashTableElementSize/2)) - 1) * (hashTableElementSize/2);
                    _ASSERT((GenomeLocationAsInt64(elementGenomeLocation) % hashTableElementSize >= (hashTableElementSize/2) ? elementGenomeLocation + (hashTableElementSize/2) : elementGenomeLocation - (hashTableElementSize/2)) == nearbyGenomeLocation);   // Assert that the logic in the above comment is right.

                    findElement(nearbyGenomeLocation, elementToScore->direction, &nearbyElement);
                } else {
                    nearbyElement = NULL;
                }

                if (NULL != nearbyElement && nearbyElement->candidatesScored != 0) {
                    //
                    // Just because there's a "nearby" element doesn't mean it's really within the maxMergeDist.  Check that now.
                    //
                    if (!genomeLocationIsWithin(genomeLocation, nearbyElement->bestScoreGenomeLocation, maxMergeDist)) {

                        //
                        // There's a nearby element, but its best score is too far away to merge.  Forget it.
                        //
                        nearbyElement = NULL;
                    }

                    if (NULL != nearbyElement) {
                        if (nearbyElement->bestScore < score || nearbyElement->bestScore == score && nearbyElement->matchProbabilityForBestScore >= matchProbability) {
                            //
                            // Again, this no better than something nearby we already tried.  Give up.
                            //
                            continue;
                        }
                        anyNearbyCandidatesAlreadyScored = true;
                        probabilityOfAllCandidates = __max(0.0, probabilityOfAllCandidates - nearbyElement->matchProbabilityForBestScore);
                        nearbyElement->matchProbabilityForBestScore = 0;    // keeps us from backing it out twice
                    }
                }

                probabilityOfAllCandidates = __max(0.0, probabilityOfAllCandidates - elementToScore->matchProbabilityForBestScore); // need the max due to floating point lossage.
                probabilityOfAllCandidates += matchProbability; // Don't combine this with the previous line, it introduces floating point unhappiness.
                elementToScore->matchProbabilityForBestScore = matchProbability;
                elementToScore->bestScore = score;

                if (bestScore > score ||
                    (bestScore == score && matchProbability > probabilityOfBestCandidate)) {

                    //
                    // We have a new best score.  The old best score becomes the second best score, unless this is the same as the best or second best score
                    //

                    if ((secondBestScore == UnusedScoreValue || !(secondBestScoreGenomeLocation + maxMergeDist > genomeLocation && secondBestScoreGenomeLocation < genomeLocation + maxMergeDist)) &&
                        (bestScore == UnusedScoreValue || !(bestScoreGenomeLocation + maxMergeDist > genomeLocation && bestScoreGenomeLocation < genomeLocation + maxMergeDist)) &&
                        (!anyNearbyCandidatesAlreadyScored || (GenomeLocationAsInt64(bestScoreGenomeLocation) / maxMergeDist != GenomeLocationAsInt64(genomeLocation) / maxMergeDist &&
                                                               GenomeLocationAsInt64(secondBestScoreGenomeLocation) / maxMergeDist != GenomeLocationAsInt64(genomeLocation) / maxMergeDist))) {
                            secondBestScore = bestScore;
                            secondBestScoreGenomeLocation = bestScoreGenomeLocation;
                            secondBestScoreDirection = primaryResult->direction;
                    }

                    //
                    // If we're tracking secondary alignments, put the old best score in as a new secondary alignment
                    //
                    if (NULL != secondaryResults && (int)(bestScore - score) <= maxEditDistanceForSecondaryResults) { // bestScore is initialized to UnusedScoreValue, which is large, so this won't fire if this is the first candidate
                        if (secondaryResultBufferSize <= *nSecondaryResults) {
                            WriteErrorMessage("Out of secondary result buffer in BaseAliner::score(), which shouldn't be possible");
                            soft_exit(1);
                        }

                        SingleAlignmentResult *result = &secondaryResults[*nSecondaryResults];
                        result->direction = primaryResult->direction;
                        result->location = bestScoreGenomeLocation;
                        result->mapq = 0;
                        result->score = bestScore;
                        result->status = MultipleHits;

                        _ASSERT(result->score != -1);

                        (*nSecondaryResults)++;
                    }

                    bestScore = score;
                    probabilityOfBestCandidate = matchProbability;
                    _ASSERT(probabilityOfBestCandidate <= probabilityOfAllCandidates);
                    bestScoreGenomeLocation = genomeLocation;
                    primaryResult->location = bestScoreGenomeLocation;
                    primaryResult->score = bestScore;
                    primaryResult->direction = elementToScore->direction;

                    lvScoresAfterBestFound = 0;
                } else {
                    if (secondBestScore > score) {
                        //
                        // A new second best.
                        //
                        secondBestScore = score;
                        secondBestScoreGenomeLocation = genomeLocation;
                        secondBestScoreDirection = elementToScore->direction;
                    }

                    //
                    // If this is close enough, record it as a secondary alignment.
                    //
                    if (-1 != maxEditDistanceForSecondaryResults && NULL != secondaryResults && (int)(bestScore - score) <= maxEditDistanceForSecondaryResults && score != -1) {
                         if (secondaryResultBufferSize <= *nSecondaryResults) {
                            WriteErrorMessage("Out of secondary result buffer in BaseAliner::score(), which shouldn't be possible");
                            soft_exit(1);
                        }

                        SingleAlignmentResult *result = &secondaryResults[*nSecondaryResults];
                        result->direction = elementToScore->direction;
                        result->location = genomeLocation;
                        result->mapq = 0;
                        result->score = score;
                        result->status = MultipleHits;

                        _ASSERT(result->score != -1);

                        (*nSecondaryResults)++;
                    }
                }

                if (stopOnFirstHit && bestScore <= maxK) {
                    // The user just wanted to find reads that match the database within some distance, but doesn't
                    // care about the best alignment. Stop now but mark the result as MultipleHits because we're not
                    // confident that it's the best one.  We don't support mapq in this secnario, because we haven't
                    // explored enough to compute it.
                    primaryResult->status = MultipleHits;
                    primaryResult->mapq = 0;
                    return true;
                }

                // Update scoreLimit since we may have improved bestScore or secondBestScore
                if (!noUkkonen) {   // If we've turned off Ukkonen, then don't drop the score limit, just leave it at maxK + extraSearchDepth always
                    scoreLimit = min(bestScore, maxK) + extraSearchDepth;
                } else {
                    _ASSERT(scoreLimit == maxK + extraSearchDepth);
                }
            }   // While candidates exist in the element
        }   // If the element could possibly affect the result

        //
        // Remove the element from the weight list.
        //
        elementToScore->allExtantCandidatesScored = true;
        elementToScore->weightNext->weightPrev = elementToScore->weightPrev;
        elementToScore->weightPrev->weightNext = elementToScore->weightNext;
        elementToScore->weightNext = elementToScore->weightPrev = elementToScore;

    } while (forceResult);

    return false;
}


    void
BaseAligner::prefetchHashTableBucket(GenomeLocation genomeLocation, Direction direction)
{
    HashTableAnchor *hashTable = candidateHashTable[direction];

    _uint64 lowOrderGenomeLocation;
    _uint64 highOrderGenomeLocation;

    decomposeGenomeLocation(genomeLocation, &highOrderGenomeLocation, &lowOrderGenomeLocation);

    _uint64 hashTableIndex = hash(highOrderGenomeLocation) % candidateHashTablesSize;

    _mm_prefetch((const char *)&hashTable[hashTableIndex], _MM_HINT_T2);
}

    bool
BaseAligner::findElement(
    GenomeLocation   genomeLocation,
    Direction        direction,
    HashTableElement **hashTableElement)
{
    HashTableAnchor *hashTable = candidateHashTable[direction];

    _uint64 lowOrderGenomeLocation;
    _uint64 highOrderGenomeLocation;

    decomposeGenomeLocation(genomeLocation, &highOrderGenomeLocation, &lowOrderGenomeLocation);

    _uint64 hashTableIndex = hash(highOrderGenomeLocation) % candidateHashTablesSize;
    HashTableAnchor *anchor = &hashTable[hashTableIndex];
    if (anchor->epoch != hashTableEpoch) {
        //
        // It's empty.
        //
        *hashTableElement = NULL;
        return false;
    }

    HashTableElement *lookedUpElement = anchor->element;
    while (NULL != lookedUpElement && lookedUpElement->baseGenomeLocation != highOrderGenomeLocation) {
        lookedUpElement = lookedUpElement->next;
    }
    *hashTableElement = lookedUpElement;
    return lookedUpElement != NULL;
}


    void
BaseAligner::findCandidate(
    GenomeLocation   genomeLocation,
    Direction        direction,
    Candidate        **candidate,
    HashTableElement **hashTableElement)
/*++

Routine Description:

    Find a candidate in the hash table, optionally allocating it if it doesn't exist (but the element does).

Arguments:

    genomeLocation - the location of the candidate we'd like to look up
    candidate - The candidate that was found or created
    hashTableElement - the hashTableElement for the candidate that was found.
    allocateNew - if this doesn't already exist, should we allocate it?

--*/
{
    _uint64 lowOrderGenomeLocation;
 
    decomposeGenomeLocation(genomeLocation, NULL, &lowOrderGenomeLocation);
    if (!findElement(genomeLocation, direction, hashTableElement)) {
        *hashTableElement = NULL;
        *candidate = NULL;
        return;
    }

    _uint64 bitForThisCandidate = (_uint64)1 << lowOrderGenomeLocation;

    *candidate = &(*hashTableElement)->candidates[lowOrderGenomeLocation];

    (*hashTableElement)->allExtantCandidatesScored = (*hashTableElement)->allExtantCandidatesScored && ((*hashTableElement)->candidatesUsed & bitForThisCandidate);
    (*hashTableElement)->candidatesUsed |= bitForThisCandidate;

}

bool doAlignerPrefetch = true;

    void
BaseAligner::allocateNewCandidate(
    GenomeLocation      genomeLocation,
    Direction           direction,
    unsigned            lowestPossibleScore,
    int                 seedOffset,
    Candidate **        candidate,
    HashTableElement ** hashTableElement)
/*++

Routine Description:

Arguments:

Return Value:

--*/
{
    HashTableAnchor *hashTable = candidateHashTable[direction];

    _uint64 lowOrderGenomeLocation;
    _uint64 highOrderGenomeLocation;

    decomposeGenomeLocation(genomeLocation, &highOrderGenomeLocation, &lowOrderGenomeLocation);

    unsigned hashTableIndex = hash(highOrderGenomeLocation) % candidateHashTablesSize;

    HashTableAnchor *anchor = &hashTable[hashTableIndex];
    if (doAlignerPrefetch) {
        _mm_prefetch((const char *)anchor, _MM_HINT_T2);    // Prefetch our anchor.  We don't have enough computation to completely hide the prefetch, but at least we get some for free here.
    }
    HashTableElement *element;

#if     DBG
    element = hashTable[hashTableIndex].element;

    while (anchor->epoch == hashTableEpoch && NULL != element && element->genomeLocation != highOrderGenomeLocation) {
        element = element->next;
    }
    _ASSERT(NULL == element || anchor->epoch != hashTableEpoch);
#endif  // DBG

    _ASSERT(nUsedHashTableElements < hashTableElementPoolSize);
    element = &hashTableElementPool[nUsedHashTableElements];
    nUsedHashTableElements++;

    if (doAlignerPrefetch) {
        //
        // Fetch the next candidate so we don't cache miss next time around.
        //
        _mm_prefetch((const char *)&hashTableElementPool[nUsedHashTableElements], _MM_HINT_T2);
    }

    element->candidatesUsed = (_uint64)1 << lowOrderGenomeLocation;
    element->candidatesScored = 0;
    element->lowestPossibleScore = lowestPossibleScore;
    element->direction = direction;
    element->weight = 1;
    element->baseGenomeLocation = highOrderGenomeLocation;
    element->bestScore = UnusedScoreValue;
    element->allExtantCandidatesScored = false;
    element->matchProbabilityForBestScore = 0;

    //
    // And insert it at the end of weight list 1.
    //
    element->weightNext = &weightLists[1];
    element->weightPrev = weightLists[1].weightPrev;
    element->weightNext->weightPrev = element;
    element->weightPrev->weightNext = element;

    *candidate = &element->candidates[lowOrderGenomeLocation];
    (*candidate)->seedOffset = seedOffset;
    *hashTableElement = element;

    highestUsedWeightList = __max(highestUsedWeightList,(unsigned)1);

    if (anchor->epoch == hashTableEpoch) {
        element->next = anchor->element;
    } else {
        anchor->epoch = hashTableEpoch;
        element->next = NULL;
    }
    anchor->element = element;

}

BaseAligner::~BaseAligner()
/*++

Routine Description:

Arguments:

Return Value:

--*/
{
    delete probDistance;

    if (hadBigAllocator) {
        //
        // Since these got allocated with the alloator rather than new, we want to call
        // their destructors without freeing their memory (which is the responsibility of
        // the owner of the allocator).
        //
        if (ownLandauVishkin) {
            if (NULL != landauVishkin) {
                landauVishkin->~LandauVishkin();
            }
            if (NULL != reverseLandauVishkin) {
                reverseLandauVishkin->~LandauVishkin();
            }
        }
    } else {

        if (ownLandauVishkin) {
            if (NULL != landauVishkin) {
                delete landauVishkin;
            }

            if (NULL != reverseLandauVishkin) {
                delete reverseLandauVishkin;
            }
        }

        BigDealloc(rcReadData);
        rcReadData = NULL;

        BigDealloc(reversedRead[FORWARD]);
        reversedRead[FORWARD] = NULL;
        reversedRead[RC] = NULL;

        BigDealloc(seedUsedAsAllocated);
        seedUsed = NULL;

        BigDealloc(candidateHashTable[FORWARD]);
        candidateHashTable[FORWARD] = NULL;

        BigDealloc(candidateHashTable[RC]);
        candidateHashTable[RC] = NULL;

        BigDealloc(weightLists);
        weightLists = NULL;

        BigDealloc(hashTableElementPool);
        hashTableElementPool = NULL;

        if (NULL != hitsPerContigCounts) {
            BigDealloc(hitsPerContigCounts);
            hitsPerContigCounts = NULL;
        }
    }
}

BaseAligner::HashTableElement::HashTableElement()
{
    init();
}

    void
BaseAligner::HashTableElement::init()
{
    weightNext = NULL;
    weightPrev = NULL;
    next = NULL;
    candidatesUsed = 0;
    baseGenomeLocation = 0;
    weight = 0;
    lowestPossibleScore = UnusedScoreValue;
    bestScore = UnusedScoreValue;
    direction = FORWARD;
    allExtantCandidatesScored = false;
    matchProbabilityForBestScore = 0;
}

    void
BaseAligner::Candidate::init()
{
    score = UnusedScoreValue;
}

    void
BaseAligner::clearCandidates() {
    hashTableEpoch++;
    nUsedHashTableElements = 0;
    highestUsedWeightList = 0;
    for (unsigned i = 1; i < numWeightLists; i++) {
        weightLists[i].weightNext = weightLists[i].weightPrev = &weightLists[i];
    }
}

    void
BaseAligner::incrementWeight(HashTableElement *element)
{
    if (element->allExtantCandidatesScored) {
        //
        // It's already scored, so it shouldn't be on a weight list.
        //
        _ASSERT(element->weightNext == element);
        _ASSERT(element->weightPrev == element);
        return;
    }
    //
    // It's possible to have elements with weight > maxSeedsToUse.  This
    // happens when a single seed occurs more than once within a particular
    // element (imagine an element with bases ATATATATATATATAT..., it will
    // match the appropriate seed at offset 0, 2, 4, etc.)  If that happens,
    // just don't let the weight get too big.
    //
    if (element->weight >= numWeightLists - 1) {
        return;
    }

    //
    // Remove it from its existing list.
    //
    element->weightNext->weightPrev = element->weightPrev;
    element->weightPrev->weightNext = element->weightNext;

    element->weight++;
    highestUsedWeightList = __max(highestUsedWeightList,element->weight);

    //
    // And insert it at the tail of the new list.
    //
    element->weightNext = &weightLists[element->weight];
    element->weightPrev = weightLists[element->weight].weightPrev;
    element->weightNext->weightPrev = element;
    element->weightPrev->weightNext = element;
}

    size_t
BaseAligner::getBigAllocatorReservation(GenomeIndex *index, bool ownLandauVishkin, unsigned maxHitsToConsider, unsigned maxReadSize,
                unsigned seedLen, unsigned numSeedsFromCommandLine, double seedCoverage, int maxSecondaryAlignmentsPerContig)
{
    unsigned maxSeedsToUse;
    if (0 != numSeedsFromCommandLine) {
        maxSeedsToUse = numSeedsFromCommandLine;
    } else {
        maxSeedsToUse = (unsigned)(maxReadSize * seedCoverage / seedLen);
    }
    size_t candidateHashTablesSize = (maxHitsToConsider * maxSeedsToUse * 3)/2;    // *1.5 for hash table slack
    size_t hashTableElementPoolSize = maxHitsToConsider * maxSeedsToUse * 2 ;   // *2 for RC
    size_t contigCounters;
    if (maxSecondaryAlignmentsPerContig > 0) {
        contigCounters = sizeof(HitsPerContigCounts)* index->getGenome()->getNumContigs();
    } else {
        contigCounters = 0;
    }

    return
        contigCounters                                                  +
        sizeof(_uint64) * 14                                            + // allow for alignment
        sizeof(BaseAligner)                                             + // our own member variables
        (ownLandauVishkin ?
            LandauVishkin<>::getBigAllocatorReservation() +
            LandauVishkin<-1>::getBigAllocatorReservation() : 0)        + // our LandauVishkin objects
        sizeof(char) * maxReadSize * 2                                  + // rcReadData
        sizeof(char) * maxReadSize * 4 + 2 * MAX_K                      + // reversed read (both)
        sizeof(BYTE) * (maxReadSize + 7 + 128) / 8                      + // seed used
        sizeof(HashTableElement) * hashTableElementPoolSize             + // hash table element pool
        sizeof(HashTableAnchor) * candidateHashTablesSize * 2           + // candidate hash table (both)
        sizeof(HashTableElement) * (maxSeedsToUse + 1);                   // weight lists
}

    void 
BaseAligner::finalizeSecondaryResults(
    SingleAlignmentResult    primaryResult,
    int                     *nSecondaryResults,                     // in/out
    SingleAlignmentResult   *secondaryResults,
    int                      maxSecondaryResults,
    int                      maxEditDistanceForSecondaryResults,
    int                      bestScore)
{
    //
    // There's no guarantee that the results are actually within the bound; the aligner records anything that's
    // within the bound when it's scored, but if we subsequently found a better fit, then it may no longer be
    // close enough.  Get rid of those now.
    //
    // NB: This code is very similar to code at the end of IntersectingPairedEndAligner::align().  Sorry.
    //

    int worstScoreToKeep = min((int)maxK, bestScore + maxEditDistanceForSecondaryResults);

    int i = 0;
    while (i < *nSecondaryResults) {
        if (secondaryResults[i].score > worstScoreToKeep) {
            //
            // This one is too bad to keep.  Move the last one from the array here and decrement the
            // count.  Don't move up i, because the one we just moved in may also be too
            // bad.
            //
            secondaryResults[i] = secondaryResults[(*nSecondaryResults)-1];
            (*nSecondaryResults)--;
        } else {
            i++;
        }
    }

    if (maxSecondaryAlignmentsPerContig > 0 && primaryResult.status != NotFound) {
        //
        // Run through the results and count the number of results per contig, to see if any of them are too big.
        //

        bool anyContigHasTooManyResults = false;

        int primaryResultContigNum = genome->getContigNumAtLocation(primaryResult.location);
        hitsPerContigCounts[primaryResultContigNum].hits = 1;
        hitsPerContigCounts[primaryResultContigNum].epoch = hashTableEpoch;

        for (i = 0; i < *nSecondaryResults; i++) {
            int contigNum = genome->getContigNumAtLocation(secondaryResults[i].location);
            if (hitsPerContigCounts[contigNum].epoch != hashTableEpoch) {
                hitsPerContigCounts[contigNum].epoch = hashTableEpoch;
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
            qsort(secondaryResults, *nSecondaryResults, sizeof(*secondaryResults), SingleAlignmentResult::compareByContigAndScore);

            //
            // Now run through and eliminate any contigs with too many hits.  We can't use the same trick at the first loop above, because the
            // counting here relies on the results being sorted.  So, instead, we just copy them as we go.
            //
            int currentContigNum = -1;
            int currentContigCount = 0;
            int destResult = 0;

            for (int sourceResult = 0; sourceResult < *nSecondaryResults; sourceResult++) {
                int contigNum = genome->getContigNumAtLocation(secondaryResults[sourceResult].location);
                if (contigNum != currentContigNum) {
                    currentContigNum = contigNum;
                    currentContigCount = (contigNum == primaryResultContigNum) ? 1 : 0;
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
    } // if maxSecondaryAlignmentsPerContig > 0

    if (*nSecondaryResults > maxSecondaryResults) {
        qsort(secondaryResults, *nSecondaryResults, sizeof(*secondaryResults), SingleAlignmentResult::compareByScore);
        *nSecondaryResults = maxSecondaryResults;   // Just truncate it
    }
}

    unsigned 
BaseAligner::getMaxSecondaryResults(unsigned maxSeedsToUse, double maxSeedCoverage, unsigned maxReadSize, unsigned maxHits, unsigned seedLength)
{

    if (0 != maxSeedsToUse) {
        return maxHits * maxSeedsToUse * NUM_DIRECTIONS; // Can't have more alignments than total possible hits 
    } else {
        return (unsigned)((maxSeedCoverage * maxReadSize + seedLength) / seedLength) * maxHits * NUM_DIRECTIONS;
    }
}
