/*++

Module Name:

    BaseAligner.cpp

Abstract:

   Single-end aligner

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

// #define TRACE_ALIGNER 1
#define EXACT_DISJOINT_MISS_COUNT 1

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
    bool            i_useAffineGap,
    bool            i_ignoreAlignmentAdjustmentsForOm,
	bool            i_altAwareness,
    bool            i_emitALTAlignments,
    int             i_maxScoreGapToPreferNonAltAlignment,
	int             i_maxSecondaryAlignmentsPerContig,
    LandauVishkin<1>*i_landauVishkin,
    LandauVishkin<-1>*i_reverseLandauVishkin,
    unsigned             i_matchReward,
    unsigned             i_subPenalty,
    unsigned             i_gapOpenPenalty,
    unsigned             i_gapExtendPenalty,
    AlignerStats   *i_stats,
    BigAllocator   *allocator) :
        genomeIndex(i_genomeIndex), maxHitsToConsider(i_maxHitsToConsider), maxK(i_maxK),
        maxReadSize(i_maxReadSize), maxSeedsToUseFromCommandLine(i_maxSeedsToUseFromCommandLine),
        maxSeedCoverage(i_maxSeedCoverage), readId(-1), extraSearchDepth(i_extraSearchDepth),
        explorePopularSeeds(false), stopOnFirstHit(false), stats(i_stats), 
        noUkkonen(i_noUkkonen), noOrderedEvaluation(i_noOrderedEvaluation), noTruncation(i_noTruncation),
		useAffineGap(i_useAffineGap), matchReward(i_matchReward), subPenalty(i_subPenalty), 
        gapOpenPenalty(i_gapOpenPenalty), gapExtendPenalty(i_gapExtendPenalty),
        minWeightToCheck(max(1u, i_minWeightToCheck)), maxSecondaryAlignmentsPerContig(i_maxSecondaryAlignmentsPerContig),
        alignmentAdjuster(i_genomeIndex->getGenome()), ignoreAlignmentAdjustmentsForOm(i_ignoreAlignmentAdjustmentsForOm),
		altAwareness(i_altAwareness), emitALTAlignments(i_emitALTAlignments),
        maxScoreGapToPreferNonAltAlignment(i_maxScoreGapToPreferNonAltAlignment)
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
    i_useAffineGap      - Use affine gap scoring for seed extension
    i_ignoreAlignmentAdjustmentsForOm - When a read score is adjusted because of soft clipping for being near the end of a contig, don't use the adjusted score when computing what to keep for -om
    i_maxSecondaryAlignmentsPerContig - Maximum secondary alignments per contig; -1 means don't limit this
    i_landauVishkin     - an externally supplied LandauVishkin string edit distance object.  This is useful if we're expecting repeated computations and use the LV cache.
    i_reverseLandauVishkin - the same for the reverse direction.
    i_matchReward       - affine gap score for a match
    i_subPenalty        - affine gap score for a substitution
    i_gapOpenPenalty    - affine gap cost for opening a gap (indel)
    i_gapExtendPenalty  - affine gap cost for extending a gap (indel)
	i_altAwareness      - treat reads mapped to ALT contigs differently than normal ones
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

    if (i_subPenalty > (i_gapOpenPenalty + i_gapExtendPenalty)) {
        WriteErrorMessage("Substitutions may be penalized too high to be seen in alignments. Make sure subPenalty < gapOpen + gapExtend\n");
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

    if (allocator) {
        // affineGap = new (allocator) AffineGap<>(i_matchReward, i_subPenalty, i_gapOpenPenalty, i_gapExtendPenalty);
        // reverseAffineGap = new (allocator) AffineGap<-1>(i_matchReward, i_subPenalty, i_gapOpenPenalty, i_gapExtendPenalty);
        affineGap = new (allocator) AffineGapVectorized<>(i_matchReward, i_subPenalty, i_gapOpenPenalty, i_gapExtendPenalty);
        reverseAffineGap = new (allocator) AffineGapVectorized<-1>(i_matchReward, i_subPenalty, i_gapOpenPenalty, i_gapExtendPenalty);
    } else {
        // affineGap = new AffineGap<>(i_matchReward, i_subPenalty, i_gapOpenPenalty, i_gapExtendPenalty);
        // reverseAffineGap = new AffineGap<-1>(i_matchReward, i_subPenalty, i_gapOpenPenalty, i_gapExtendPenalty);
        affineGap = new AffineGapVectorized<>(i_matchReward, i_subPenalty, i_gapOpenPenalty, i_gapExtendPenalty); // This is a bad idea, it'll result in false sharing in the single-end aligner.  Use BigAlloc().
        reverseAffineGap = new AffineGapVectorized<-1>(i_matchReward, i_subPenalty, i_gapOpenPenalty, i_gapExtendPenalty);
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
        seedUsed = (BYTE *)allocator->allocate((sizeof(BYTE) * ((_int64)maxReadSize + 7 + 128) / 8));    // +128 to make sure it extends at both
    } else {
        seedUsed = (BYTE *)BigAlloc((sizeof(BYTE) * ((_int64)maxReadSize + 7 + 128) / 8));    // +128 to make sure it extends at both
    }

    seedUsedAsAllocated = seedUsed; // Save the pointer for the delete.
    seedUsed += 8;  // This moves the pointer up an _int64, so we now have the appropriate before buffer.

    nUsedHashTableElements = 0;

    if (allocator) {
        candidateHashTable[FORWARD] = (HashTableAnchor *)allocator->allocate(sizeof(HashTableAnchor) * candidateHashTablesSize);
        candidateHashTable[RC] = (HashTableAnchor *)allocator->allocate(sizeof(HashTableAnchor) * candidateHashTablesSize);
        weightLists = (HashTableElement *)allocator->allocate(sizeof(HashTableElement) * numWeightLists);
        hashTableElementPool = (HashTableElement *)allocator->allocate(sizeof(HashTableElement) * hashTableElementPoolSize); // Allocate last, because it's biggest and usually unused.  This puts all of the commonly used stuff into one large page.
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
        } else {
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
} // BaseAligner::BaseAligner


#ifdef  _DEBUG
bool _DumpAlignments = false;
#endif  // _DEBUG

    bool
BaseAligner::AlignRead(
        Read                    *inputRead,
        SingleAlignmentResult   *primaryResult,
        SingleAlignmentResult   *firstALTResult,
        int                      maxEditDistanceForSecondaryResults,
        _int64                   secondaryResultBufferSize,
        _int64                  *nSecondaryResults,
        _int64                   maxSecondaryResults,
        SingleAlignmentResult   *secondaryResults,             // The caller passes in a buffer of secondaryResultBufferSize and it's filled in by AlignRead()
        _int64                   maxCandidatesForAffineGapBufferSize,
        _int64                  *nCandidatesForAffineGap,
        SingleAlignmentResult   *candidatesForAffineGap, // Alignment candidates that need to be rescored using affine gap
        bool                     useHamming
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
#if _DEBUG
    const size_t genomeLocationBufferSize = 200;
    char genomeLocationBuffer[genomeLocationBufferSize];
#endif // _DEBUG

    bool overflowedSecondaryResultsBuffer = false;
    memset(hitCountByExtraSearchDepth, 0, sizeof(*hitCountByExtraSearchDepth) * extraSearchDepth);

    if (NULL != nSecondaryResults) {
        *nSecondaryResults = 0;
    }

    firstPassSeedsNotSkipped[FORWARD] = firstPassSeedsNotSkipped[RC] = 0;
    highestWeightListChecked = 0;

    scoresForAllAlignments.bestScore = scoresForNonAltAlignments.bestScore = TooBigScoreValue;

    unsigned maxSeedsToUse;
    if (0 != maxSeedsToUseFromCommandLine) {
        maxSeedsToUse = maxSeedsToUseFromCommandLine;
    } else {
        maxSeedsToUse = (int)(NUM_DIRECTIONS * maxSeedCoverage * inputRead->getDataLength() / genomeIndex->getSeedLength()); 
    }

    primaryResult->location = InvalidGenomeLocation; // Value to return if we don't find a location.
    primaryResult->direction = FORWARD;              // So we deterministically print the read forward in this case.
    primaryResult->score = UnusedScoreValue;
    primaryResult->status = NotFound;
    primaryResult->clippingForReadAdjustment = 0;
    primaryResult->usedAffineGapScoring = false;
    primaryResult->basesClippedBefore = 0;
    primaryResult->basesClippedAfter = 0;
    primaryResult->agScore = 0;
    primaryResult->seedOffset = 0;
    primaryResult->supplementary = false;

    unsigned lookupsThisRun = 0;

    popularSeedsSkipped = 0;
    nAddedToHashTable = 0;

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
        return true;
    }

#ifdef TRACE_ALIGNER
    printf("Aligning read '%.*s':\n%.*s\n%.*s\n", inputRead->getIdLength(), inputRead->getId(), inputRead->getDataLength(), inputRead->getData(),
            inputRead->getDataLength(), inputRead->getQuality());
#endif

#ifdef  _DEBUG
    if (_DumpAlignments) {
        printf("BaseAligner: aligning read ID '%.*s', data '%.*s' %s\n", inputRead->getIdLength(), inputRead->getId(), inputRead->getDataLength(), inputRead->getData(), useHamming ? "Hamming" : "");
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
        return true;
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
    wrapCount = 0;
    lowestPossibleScoreOfAnyUnseenLocation[FORWARD] = lowestPossibleScoreOfAnyUnseenLocation[RC] = 0;
    currRoundLowestPossibleScoreOfAnyUnseenLocation[FORWARD] = currRoundLowestPossibleScoreOfAnyUnseenLocation[RC] = 0;
    mostSeedsContainingAnyParticularBase[FORWARD] = mostSeedsContainingAnyParticularBase[RC] = 1;  // Instead of tracking this for real, we're just conservative and use wrapCount+1.  It's faster.

    scoresForAllAlignments.init();
    if (altAwareness) {
        scoresForNonAltAlignments.init();
    }

    nSeedsApplied[FORWARD] = nSeedsApplied[RC] = 0;
    lvScores = 0;
    lvScoresAfterBestFound = 0;

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
                printf("Calling score with force=true because we wrapped around enough\n");
#endif

                score(
                    true,
                    read,
                    primaryResult,
                    firstALTResult,
                    maxEditDistanceForSecondaryResults,
                    secondaryResultBufferSize,
                    nSecondaryResults,
                    secondaryResults,
                    &overflowedSecondaryResultsBuffer,
                    maxCandidatesForAffineGapBufferSize,
                    nCandidatesForAffineGap,
                    candidatesForAffineGap,
                    useHamming);

#ifdef  _DEBUG
                if (_DumpAlignments) printf("Final result score %d MAPQ %d (%e probability of best candidate, %e probability of all candidates, non ALT-aware)  at %s:%llu\n\n", 
                                            primaryResult->score, primaryResult->mapq, scoresForAllAlignments.probabilityOfBestCandidate, scoresForAllAlignments.probabilityOfAllCandidates, 
                                            genome->getContigAtLocation(primaryResult->location)->name, primaryResult->location - genome->getContigAtLocation(primaryResult->location)->beginningLocation);
#endif  // _DEBUG
                if (overflowedSecondaryResultsBuffer) {
                    return false;
                }

                finalizeSecondaryResults(read[FORWARD], primaryResult, nSecondaryResults, secondaryResults, maxSecondaryResults, maxEditDistanceForSecondaryResults, primaryResult->score);
                return true;
            }
            nextSeedToTest = GetWrappedNextSeedToTest(seedLen, wrapCount);

            mostSeedsContainingAnyParticularBase[FORWARD] = mostSeedsContainingAnyParticularBase[RC] = wrapCount + 1;

            currRoundLowestPossibleScoreOfAnyUnseenLocation[FORWARD] = currRoundLowestPossibleScoreOfAnyUnseenLocation[RC] = 0;
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
            printf("\tSeed offset %2d, %4lld hits, %4lld rcHits.", nextSeedToTest, nHits[0], nHits[1]);
            for (int rc = 0; rc < 2; rc++) {
                for (unsigned i = 0; i < __min(nHits[rc], 2); i++) {
                    printf(" %sHit at %s.", rc == 1 ? "RC " : "", genome->genomeLocationInStringForm(doesGenomeIndexHave64BitLocations ? hits[rc][i].location : (_int64)hits32[rc][i], genomeLocationBuffer, genomeLocationBufferSize));
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
                    printf(" %9llu", doesGenomeIndexHave64BitLocations ? hits[rc][i].location : (_int64)hits32[rc][i]);
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
                    // will also return the results for the seed's reverse complement, i.e., 2'1'0'.
                    // This happens as the last seed in the RC read.
                    //
                    offset = readLen - seedLen - nextSeedToTest;
                }

                const unsigned prefetchDepth = 6;

                //
                // We keep prefetches outstanding prefetchDepth ahead of where we are.  Start by launching the first
                // prefetchDepth of them, and then each time we process a hit launch the prefetch for prefetchDepth farther
                // along, assuming it exists.
                //

                if (doAlignerPrefetch) {
                    _int64 prefetchLimit = __min(prefetchDepth, __min(nHits[direction], (_int64)maxHitsToConsider));
                    for (int prefetchIndex = 0; prefetchIndex < prefetchLimit; prefetchIndex++) {
                        if (doesGenomeIndexHave64BitLocations) {
                            prefetchHashTableBucket(GenomeLocationAsInt64(hits[direction][prefetchIndex]) - offset, direction);
                        }
                        else {
                            prefetchHashTableBucket(hits32[direction][prefetchIndex] - offset, direction);
                        }
                    }
                }

                _int64 limit = min(nHits[direction], (_int64)maxHitsToConsider);
 

                for (unsigned i = 0; i < limit; i++) {
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

                    bool candidateIsALT = altAwareness && genome->isGenomeLocationALT(genomeLocationOfThisHit);

                    if (NULL != hashTableElement) {
                        if (!noOrderedEvaluation) {     // If noOrderedEvaluation, just leave them all on the one-hit weight list so they get evaluated in whatever order
                            incrementWeight(hashTableElement);
                        }
                        candidate->seedOffset = offset;
                        _ASSERT((unsigned)candidate->seedOffset <= readLen - seedLen);
                    } else if (lowestPossibleScoreOfAnyUnseenLocation[direction] <= scoreLimit(candidateIsALT) || noTruncation) {
                            _ASSERT(offset <= readLen - seedLen);
                            allocateNewCandidate(genomeLocationOfThisHit, direction, lowestPossibleScoreOfAnyUnseenLocation[direction],
                                offset, &candidate, &hashTableElement);
                            nAddedToHashTable++;
                    }

                    if (doAlignerPrefetch && (_int64)i + prefetchDepth < limit) {
                        if (doesGenomeIndexHave64BitLocations) {
                            prefetchHashTableBucket(GenomeLocationAsInt64(hits[direction][i + prefetchDepth]) - offset, direction);
                        }
                        else {
                            prefetchHashTableBucket(hits32[direction][i + prefetchDepth] - offset, direction);
                        }
                    }
                }

                nSeedsApplied[direction]++;
                currRoundLowestPossibleScoreOfAnyUnseenLocation[direction]++;
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
                firstALTResult,
                maxEditDistanceForSecondaryResults,
                secondaryResultBufferSize,
                nSecondaryResults,
                secondaryResults,
                &overflowedSecondaryResultsBuffer,
                maxCandidatesForAffineGapBufferSize,
                nCandidatesForAffineGap,
                candidatesForAffineGap,
                useHamming)) {

#ifdef  _DEBUG
                if (_DumpAlignments) printf("Final result score %d MAPQ %d at %s:%llu\n", primaryResult->score, primaryResult->mapq,
                    genome->getContigAtLocation(primaryResult->location)->name, primaryResult->location - genome->getContigAtLocation(primaryResult->location)->beginningLocation);
#endif  // _DEBUG
                if (overflowedSecondaryResultsBuffer) {
                    return false;
                }
                finalizeSecondaryResults(read[FORWARD], primaryResult, nSecondaryResults, secondaryResults, maxSecondaryResults, maxEditDistanceForSecondaryResults, primaryResult->score);
                return true;
            } // If score says we have a difinitive answer
        } // If we applied a seed, and so something's changed.
    } // While we're still applying seeds

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
        firstALTResult,
        maxEditDistanceForSecondaryResults,
        secondaryResultBufferSize,
        nSecondaryResults,
        secondaryResults,
        &overflowedSecondaryResultsBuffer,
        maxCandidatesForAffineGapBufferSize,
        nCandidatesForAffineGap,
        candidatesForAffineGap,
        useHamming);

#ifdef  _DEBUG
    if (_DumpAlignments) printf("Final result score %d MAPQ %d (%e probability of best candidate, %e probability of all candidates non ALT-aware) at %s:%llu\n", 
        primaryResult->score, primaryResult->mapq, scoresForAllAlignments.probabilityOfBestCandidate, scoresForAllAlignments.probabilityOfAllCandidates, 
        genome->getContigAtLocation(primaryResult->location)->name, primaryResult->location - genome->getContigAtLocation(primaryResult->location)->beginningLocation);
    if (_DumpAlignments && firstALTResult->status != NotFound) printf("Emitting ALT result score %d MAPQ %d at %llu\n",
        firstALTResult->score, firstALTResult->mapq, firstALTResult->location.location);
#endif  // _DEBUG

    if (overflowedSecondaryResultsBuffer) {
        return false;
    }

    finalizeSecondaryResults(read[FORWARD], primaryResult, nSecondaryResults, secondaryResults, maxSecondaryResults, maxEditDistanceForSecondaryResults, primaryResult->score);
    return true;
}

    void
BaseAligner::scoreLocationWithAffineGap(
    Read* reads[NUM_DIRECTIONS],
    Direction            direction,
    GenomeLocation       genomeLocation,
    unsigned             seedOffset,
    int                  scoreLimit,
    int* score,
    double* matchProbability,
    int* genomeLocationOffset,
    int* basesClippedBefore,
    int* basesClippedAfter,
    int* agScore
)
{
    Read* readToScore = reads[direction];
    unsigned readDataLength = readToScore->getDataLength();
    GenomeDistance genomeDataLength = (GenomeDistance)readDataLength + MAX_K; // Leave extra space in case the read has deletions
    const char* data = genome->getSubstring(genomeLocation, genomeDataLength);

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
    int tailStart = seedOffset + seedLen;
    int agScore1 = seedLen, agScore2 = 0; // affine gap scores

    _ASSERT(!memcmp(data + seedOffset, readToScore->getData() + seedOffset, seedLen));    // that the seed actually matches


    int textLen;
    if (genomeDataLength - tailStart > INT32_MAX) {
        textLen = INT32_MAX;
    }
    else {
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
                    reversedRead[direction] + readLen - seedOffset,
                    reads[OppositeDirection(direction)]->getQuality() + readLen - seedOffset,
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
                    reversedRead[direction] + readLen - seedOffset,
                    reads[OppositeDirection(direction)]->getQuality() + readLen - seedOffset,
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

    bool
BaseAligner::score(
        bool                     forceResult,
        Read                    *read[NUM_DIRECTIONS],
        SingleAlignmentResult   *primaryResult,
        SingleAlignmentResult   *firstALTResult,                        // This only gets filled in if there's a good ALT result that's not the primary result && emitALTAlignments
        int                      maxEditDistanceForSecondaryResults,
        _int64                   secondaryResultBufferSize,
        _int64                  *nSecondaryResults,
        SingleAlignmentResult   *secondaryResults,
        bool                    *overflowedSecondaryBuffer,
        _int64                   maxCandidatesForAffineGapBufferSize,
        _int64                  *nCandidatesForAffineGap,
        SingleAlignmentResult   *candidatesForAffineGap, // Alignment candidates that need to be rescored using affine gap
        bool                     useHamming)
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
#if _DEBUG
    const size_t genomeLocationBufferSize = 200;
    char genomeLocationBuffer[genomeLocationBufferSize];
#endif // _DEBUG

    *overflowedSecondaryBuffer = false;
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
    for (Direction direction = 0; direction < NUM_DIRECTIONS; direction++) {
        if (0 != mostSeedsContainingAnyParticularBase[direction]) {
#ifdef EXACT_DISJOINT_MISS_COUNT
            lowestPossibleScoreOfAnyUnseenLocation[direction] =
                __max(lowestPossibleScoreOfAnyUnseenLocation[direction],
                    currRoundLowestPossibleScoreOfAnyUnseenLocation[direction]);
#else
            lowestPossibleScoreOfAnyUnseenLocation[direction] =
                __max(lowestPossibleScoreOfAnyUnseenLocation[direction],
                      nSeedsApplied[direction] / mostSeedsContainingAnyParticularBase[direction]);
#endif
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


        if ((__min(lowestPossibleScoreOfAnyUnseenLocation[FORWARD],lowestPossibleScoreOfAnyUnseenLocation[RC]) > max(scoreLimit(true), scoreLimit(false)) && !noTruncation) || forceResult) {
            if (weightListToCheck < minWeightToCheck) {
                //
                // We've scored all live candidates and excluded all non-candidates, or we've checked enough that we've hit the cutoff.  We have our
                // answer.
                //
                ScoreSet *scoreSetOfFinalResult;
                if (!altAwareness || scoresForNonAltAlignments.bestScore > scoresForAllAlignments.bestScore + maxScoreGapToPreferNonAltAlignment) {
                    scoreSetOfFinalResult = &scoresForAllAlignments;
                    firstALTResult->status = NotFound;
                } else {
                    scoreSetOfFinalResult = &scoresForNonAltAlignments;
                    if (emitALTAlignments && scoresForAllAlignments.bestScore <= scoresForNonAltAlignments.bestScore && scoresForAllAlignments.bestScoreGenomeLocation != scoresForNonAltAlignments.bestScoreGenomeLocation) {
                        scoresForAllAlignments.fillInSingleAlignmentResult(firstALTResult, popularSeedsSkipped);
                    } else {
                        firstALTResult->status = NotFound;
                    }
                }

                primaryResult->score = scoreSetOfFinalResult->bestScore;
                if (scoreSetOfFinalResult->bestScore <= maxK || (useHamming && scoreSetOfFinalResult->bestScore != UnusedScoreValue)) {
                    scoreSetOfFinalResult->fillInSingleAlignmentResult(primaryResult, popularSeedsSkipped);
                    primaryResult->supplementary = false;
                    return true;
                } else {
                    primaryResult->status = NotFound;
                    primaryResult->mapq = 0;
                    return true;
                }
            } // If we don't still have weight lists to check

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

        int scoreLimitForThisElement = scoreLimit(altAwareness && genome->isGenomeLocationALT(elementToScore->baseGenomeLocation)); // All nearby genome locations are either ALT or non-ALT, so it's OK to be close here
        if (elementToScore->lowestPossibleScore <= scoreLimitForThisElement) {

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
                GenomeLocation origGenomeLocation = genomeLocation;
                GenomeLocation elementGenomeLocation = genomeLocation;    // This is the genome location prior to any adjustments for indels

                bool genomeLocationIsNonALT = (!altAwareness) || !genome->isGenomeLocationALT(genomeLocation);

                //
                // We're about to run edit distance computation on the genome.  Launch a prefetch for it
                // so that it's in cache when we do (or at least on the way).
                //
                if (doAlignerPrefetch) {
                    genomeIndex->prefetchGenomeData(genomeLocation);
                }

                unsigned score = ScoreAboveLimit;
                double matchProbability = 0.0;
                unsigned readDataLength = read[elementToScore->direction]->getDataLength();
                GenomeDistance genomeDataLength = (GenomeDistance)readDataLength + MAX_K; // Leave extra space in case the read has deletions
                const char *data = genome->getSubstring(genomeLocation, genomeDataLength);
                bool usedAffineGapScoring = false;
                int basesClippedBefore = 0;
                int basesClippedAfter = 0;
                int agScore = -1;
                int scoreGapless = -1;

                if (data != NULL) {
                    Read *readToScore = read[elementToScore->direction];

                    _ASSERT(candidateToScore->seedOffset + seedLen <= readToScore->getDataLength());

                    //
                    // Compute the distance separately in the forward and backward directions from the seed, to allow
                    // arbitrary offsets at both the start and end.
                    //
                    double matchProb1 = 1.0, matchProb2 = 1.0;
                    int score1 = 0, score2 = 0;
                    // First, do the forward direction from where the seed aligns to past of it
                    int readLen = readToScore->getDataLength();
                    int seedLen = genomeIndex->getSeedLength();
                    int seedOffset = candidateToScore->seedOffset; // Since the data is reversed
                    int tailStart = seedOffset + seedLen;
                    int agScore1 = seedLen, agScore2 = 0;

                    // Compute maxK for which edit distance and affine gap scoring report the same alignment given a set of scoring parameters
                    // GapOpenPenalty + k.GapExtendPenalty >= k * SubPenalty
                    int maxKForSameAlignment = gapOpenPenalty / (subPenalty - gapExtendPenalty);

                    int totalIndels = 0;
                    int genomeLocationOffset = 0;

                    _ASSERT(!memcmp(data+seedOffset, readToScore->getData() + seedOffset, seedLen));
                    int textLen = (int)__min(genomeDataLength - tailStart, 0x7ffffff0);

                    // Try gapless scoring to see if we can align the read after clipping
                    int score1Gapless = 0, score2Gapless = 0; // gapless scores are only for the unclipped portions

                    if (!useHamming) {
                        score1 = landauVishkin->computeEditDistance(data + tailStart, textLen, readToScore->getData() + tailStart, readToScore->getQuality() + tailStart, readLen - tailStart,
                            scoreLimitForThisElement, &matchProb1, NULL, &totalIndels);

                        agScore1 = (seedLen + readLen - tailStart - score1) * matchReward - score1 * subPenalty;

                        if (score1 != ScoreAboveLimit) {
                            // The tail of the read matched; now let's reverse match the reference genome and the head
                            int limitLeft = scoreLimitForThisElement - score1;
                            totalIndels = 0;
                            score2 = reverseLandauVishkin->computeEditDistance(data + seedOffset, seedOffset + MAX_K, reversedRead[elementToScore->direction] + readLen - seedOffset,
                                read[OppositeDirection(elementToScore->direction)]->getQuality() + readLen - seedOffset, seedOffset, limitLeft, &matchProb2,
                                &genomeLocationOffset, &totalIndels);

                            agScore2 = (seedOffset - score2) * matchReward - score2 * subPenalty;
                        }
                    } else {
                        if (tailStart != readLen) {
                            agScore1 = affineGap->computeGaplessScore(data + tailStart, textLen, readToScore->getData() + tailStart, readToScore->getQuality() + tailStart, readLen - tailStart,
                                readLen, scoreLimitForThisElement, &score1, NULL, NULL, &matchProb1, &score1Gapless);
                            agScore1 += (seedLen - readLen);
                        }
                        if (score1Gapless != ScoreAboveLimit) {
                            int limitLeft = scoreLimitForThisElement - score1Gapless;
                            if (seedOffset != 0) {
                                agScore2 = reverseAffineGap->computeGaplessScore(data + seedOffset, seedOffset + MAX_K, reversedRead[elementToScore->direction] + readLen - seedOffset,
                                    read[OppositeDirection(elementToScore->direction)]->getQuality() + readLen - seedOffset, seedOffset, readLen, limitLeft, &score2, &genomeLocationOffset, NULL, &matchProb2, &score2Gapless);

                                agScore2 -= (readLen);

                                if (score2Gapless == ScoreAboveLimit) {
                                    score = ScoreAboveLimit;
                                    genomeLocationOffset = 0;
                                    agScore = ScoreAboveLimit;
                                }
                            }
                        }
                    }

                    if (!useHamming && (score1 != ScoreAboveLimit && score2 != ScoreAboveLimit)) {
                        // Check if affine gap must be called
                        if (useAffineGap && ((score1 + score2) > maxKForSameAlignment) && elementToScore->lowestPossibleScore <= scoresForAllAlignments.bestScore) {
                            score1 = 0;  score2 = 0;  agScore1 = seedLen; agScore2 = 0;
                            usedAffineGapScoring = true;
                            if (tailStart != readLen) {
                                int patternLen = readLen - tailStart;
                                //
                                // Try banded affine-gap when pattern is long and band needed is small
                                //
                                if (patternLen >= (3 * (2 * scoreLimitForThisElement + 1))) {
                                    agScore1 = affineGap->computeScoreBanded(data + tailStart,
                                        textLen,
                                        readToScore->getData() + tailStart,
                                        readToScore->getQuality() + tailStart,
                                        readLen - tailStart,
                                        scoreLimitForThisElement,
                                        readLen,
                                        NULL,
                                        &basesClippedAfter,
                                        &score1,
                                        &matchProb1);
                                }
                                else {
                                    agScore1 = affineGap->computeScore(data + tailStart,
                                        textLen,
                                        readToScore->getData() + tailStart,
                                        readToScore->getQuality() + tailStart,
                                        readLen - tailStart,
                                        scoreLimitForThisElement,
                                        readLen,
                                        NULL,
                                        &basesClippedAfter,
                                        &score1,
                                        &matchProb1);
                                }

                                agScore1 += (seedLen - readLen);
                            }

                            if (score1 != ScoreAboveLimit) {
                                if (seedOffset != 0) {
                                    int limitLeft = scoreLimitForThisElement - score1;
                                    int patternLen = seedOffset;
                                    //
                                    // Try banded affine-gap when pattern is long and band needed is small
                                    //
                                    if (patternLen >= (3 * (2 * limitLeft + 1))) {
                                        agScore2 = reverseAffineGap->computeScoreBanded(data + seedOffset,
                                            seedOffset + limitLeft,
                                            reversedRead[elementToScore->direction] + readLen - seedOffset,
                                            read[OppositeDirection(elementToScore->direction)]->getQuality() + readLen - seedOffset,
                                            seedOffset,
                                            limitLeft,
                                            readLen,
                                            &genomeLocationOffset,
                                            &basesClippedBefore,
                                            &score2,
                                            &matchProb2);
                                    }
                                    else {
                                        agScore2 = reverseAffineGap->computeScore(data + seedOffset,
                                            seedOffset + limitLeft,
                                            reversedRead[elementToScore->direction] + readLen - seedOffset,
                                            read[OppositeDirection(elementToScore->direction)]->getQuality() + readLen - seedOffset,
                                            seedOffset,
                                            limitLeft,
                                            readLen,
                                            &genomeLocationOffset,
                                            &basesClippedBefore,
                                            &score2,
                                            &matchProb2);
                                    }

                                    agScore2 -= (readLen);

                                    if (score2 == ScoreAboveLimit) {
                                        score = ScoreAboveLimit;
                                        agScore = -1;
                                    }
                                }
                            } else {
                                score = ScoreAboveLimit;
                                agScore = -1;
                            }
                        }
                    }

                    if ((score1 != ScoreAboveLimit && score2 != ScoreAboveLimit) || (useHamming && score1Gapless != ScoreAboveLimit && score2Gapless != ScoreAboveLimit)) {
                        score = score1 + score2;
                        scoreGapless = score1Gapless + score2Gapless;
                        // Map probabilities for substrings can be multiplied, but make sure to count seed too
                        matchProbability = matchProb1 * matchProb2 * pow(1 - SNP_PROB, seedLen);

                        //
                        // Adjust the genome location based on any indels that we found.
                        //
                        genomeLocation += genomeLocationOffset;

                        agScore = agScore1 + agScore2;

                        //
                        // We could mark as scored anything in between the old and new genome offsets, but it's probably not worth the effort since this is
                        // so rare and all it would do is same time.
                        //
                    } else {
                        score = ScoreAboveLimit;
                        agScore = ScoreAboveLimit;
                        matchProbability = 0.0;
                    }
                } else { // if we had genome data to compare against
                    matchProbability = 0.0;
                }
#ifdef TRACE_ALIGNER
                printf("Computing distance at %u (RC) with limit %d: %d (prob %g)\n",
                        genomeLocation, scoreLimit, score, matchProbability);
#endif


#ifdef  _DEBUG
                if (_DumpAlignments) printf("\t\t%cScored %s weight %2d limit %d, result %2d %s, agScore %d, usedAffine %d, gaplessScore %d, usedHamming %d, matchProb %g. %d added to hash table\n",
                    score != ScoreAboveLimit ? '*' : ' ',
                    genome->genomeLocationInStringForm(genomeLocation.location, genomeLocationBuffer, genomeLocationBufferSize), elementToScore->weight, scoreLimitForThisElement, score, 
                    (elementToScore->direction ? "RC" : ""), agScore, usedAffineGapScoring, scoreGapless, useHamming, matchProbability,
                    nAddedToHashTable);
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
                    //
                    // Match probability is a better indicator of "goodness" of alignment for Hamming distance scoring
                    //
                    if (useHamming && matchProbability <= elementToScore->matchProbabilityForBestScore) {
                        continue;
                    }
                    if (elementToScore->bestScore < score || (elementToScore->bestScore == score && matchProbability <= elementToScore->matchProbabilityForBestScore)) {
                    // if (matchProbability <= elementToScore->matchProbabilityForBestScore) {
						//
                        // This is a no better mapping than something nearby that we already tried.  Just ignore it.
                        //
                        continue;
                    }
                } else {
                    _ASSERT(elementToScore->matchProbabilityForBestScore == 0.0);
                }

                elementToScore->bestScoreGenomeLocation = genomeLocation;
                elementToScore->bestScoreOrigGenomeLocation = origGenomeLocation;
                elementToScore->usedAffineGapScoring = usedAffineGapScoring;
                elementToScore->basesClippedBefore = basesClippedBefore;
                elementToScore->basesClippedAfter = basesClippedAfter;
                elementToScore->agScore = agScore;
                elementToScore->seedOffset = candidateToScore->seedOffset;

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
                if (ScoreAboveLimit != score && score < 2) {
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
                    } else {
                        if (useHamming && nearbyElement->matchProbabilityForBestScore >= matchProbability) {
                            continue;
                        }
                        if (nearbyElement->bestScore < score || (nearbyElement->bestScore == score && nearbyElement->matchProbabilityForBestScore >= matchProbability)) {
                           //
                            // Again, this is no better than something nearby we already tried.  Give up.
                            //
                            continue;
                        }
                        
                        scoresForAllAlignments.updateProbabilitiesForNearbyMatch(nearbyElement->matchProbabilityForBestScore);
                        if (genomeLocationIsNonALT) {
                            scoresForNonAltAlignments.updateProbabilitiesForNearbyMatch(nearbyElement->matchProbabilityForBestScore);
                        }
                        anyNearbyCandidatesAlreadyScored = true;
                        nearbyElement->matchProbabilityForBestScore = 0;    // keeps us from backing it out twice
                    }
                }

                scoresForAllAlignments.updateProbabilitiesForNewMatch(matchProbability, elementToScore->matchProbabilityForBestScore);
                if (genomeLocationIsNonALT) {
                    scoresForNonAltAlignments.updateProbabilitiesForNewMatch(matchProbability, elementToScore->matchProbabilityForBestScore);
                }

                elementToScore->matchProbabilityForBestScore = matchProbability;
                elementToScore->bestScore = score;

                if (useHamming) {
                    scoresForAllAlignments.updateBestScore(genomeLocation, origGenomeLocation, score, usedAffineGapScoring, agScore, matchProbability, lvScoresAfterBestFound, elementToScore,
                                                secondaryResults, nSecondaryResults, secondaryResultBufferSize,
                                                anyNearbyCandidatesAlreadyScored, maxEditDistanceForSecondaryResults, overflowedSecondaryBuffer,
                                                maxCandidatesForAffineGapBufferSize, nCandidatesForAffineGap, candidatesForAffineGap, extraSearchDepth);
                }
                else {
                    scoresForAllAlignments.updateBestScore(genomeLocation, origGenomeLocation, score, usedAffineGapScoring, agScore, matchProbability, lvScoresAfterBestFound, elementToScore,
                        secondaryResults, nSecondaryResults, secondaryResultBufferSize,
                        anyNearbyCandidatesAlreadyScored, maxEditDistanceForSecondaryResults, overflowedSecondaryBuffer,
                        0, NULL, NULL, extraSearchDepth);
                }

                if (*overflowedSecondaryBuffer) {
                    return false;
                }

                if (NULL != candidatesForAffineGap && *nCandidatesForAffineGap >= maxCandidatesForAffineGapBufferSize) {
                    *nCandidatesForAffineGap = maxCandidatesForAffineGapBufferSize + 1;
                    return false;
                }

                if (genomeLocationIsNonALT) {
                    //
                    // Don't update secondary results here; we don't exclude ALT alignments from them, only from the primary result.
                    //
                    scoresForNonAltAlignments.updateBestScore(genomeLocation, origGenomeLocation, score, usedAffineGapScoring, agScore, matchProbability, lvScoresAfterBestFound, elementToScore,
                                                    NULL, 0, 0, anyNearbyCandidatesAlreadyScored, -1, NULL, 0, NULL, NULL, extraSearchDepth);

                }
                                        

                if (stopOnFirstHit && ((scoresForAllAlignments.bestScore <= maxK) || (useHamming && scoresForAllAlignments.bestScore != UnusedScoreValue))) {
                    // The user just wanted to find reads that match the database within some distance, but doesn't
                    // care about the best alignment. Stop now but mark the result as MultipleHits because we're not
                    // confident that it's the best one.  We don't support mapq in this secnario, because we haven't
                    // explored enough to compute it.
                    primaryResult->status = MultipleHits;
                    primaryResult->mapq = 0;
                    return true;
                }

                // Taken from intersecting paired-end aligner.
                // Assuming maximum probability among unseen candidates is 1 and MAPQ < 1, find probability of
                // all candidates for which we can terminate early without exploring any more MAPQ < 1 alignments
                // i.e., -10 log10(1 - 1/x) < 1 
                // i.e.,  x > 4.89 ~ 4.9
                if ((altAwareness ? scoresForNonAltAlignments.probabilityOfAllCandidates : scoresForAllAlignments.probabilityOfAllCandidates) >= 4.9 && -1 == maxEditDistanceForSecondaryResults) {
                    //
                    // nothing will rescue us from a 0 mapq, so just stop looking.
                    //
                    (altAwareness ? scoresForNonAltAlignments : scoresForAllAlignments).fillInSingleAlignmentResult(primaryResult, popularSeedsSkipped);
                    firstALTResult->status = NotFound;
                    return true;
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

    bool
BaseAligner::alignAffineGap(
        Read* inputRead,
        SingleAlignmentResult* result,
        SingleAlignmentResult* firstALTResult,
        _int64                 nCandidatesForAffineGap,
        SingleAlignmentResult* candidatesForAffineGap // Alignment candidates that need to be rescored using affine gap
    )
{
    if (result->status == NotFound) {
        return true;
    }

    unsigned readLen = inputRead->getDataLength();
    const char* readData = inputRead->getData();
    const char* readQuality = inputRead->getQuality();
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
        return true;
    }

    Read reverseComplementRead;
    Read* read[NUM_DIRECTIONS];
    read[FORWARD] = inputRead;
    read[RC] = &reverseComplementRead;
    read[RC]->init(inputRead->getId(), inputRead->getIdLength(), rcReadData, rcReadQuality, readLen);

    int bestScore = result->score;
    int scoreLimitForCandidate, scoreLimitForCandidateALT;

    scoreLimitForCandidate = scoreLimitForCandidateALT = MAX_K - 1;

    int genomeOffset = 0;
    bool skipAffineGap = false;

    //
    // Keep track of old bestPairProbability as this is used in updating the new match probability after affine gap scoring
    //
    double oldProbabilityBestResult = result->matchProbability;
    double oldProbabilityBestResultALT = (firstALTResult->status != NotFound) ? firstALTResult->matchProbability : 0.0;

    int maxKForSameAlignment = gapOpenPenalty / (subPenalty - gapExtendPenalty);

    result->usedAffineGapScoring = false;
    if (result->score > maxKForSameAlignment) {
        result->usedAffineGapScoring = true;
        scoreLocationWithAffineGap(read, result->direction, result->origLocation,
            result->seedOffset, scoreLimitForCandidate, &result->score, &result->matchProbability,
            &genomeOffset, &result->basesClippedBefore, &result->basesClippedAfter, &result->agScore);

        if (result->score != ScoreAboveLimit) {
            result->location = result->origLocation + genomeOffset;
        } else {
            result->status = NotFound;
        }

#if _DEBUG
        if (_DumpAlignments) {
            fprintf(stderr, "Affine gap scored read at %s:%llu score %d, agScore %d\n",
                genome->getContigAtLocation(result->location)->name, result->location - genome->getContigAtLocation(result->location)->beginningLocation,
                result->score, result->agScore
            );
        }
#endif  // _DEBUG

    } else {
        skipAffineGap = true;
    }

    //
    // Use affine gap scoring for ALT result if it was computed in Phase 3
    //
    if (firstALTResult->status != NotFound) {
        if (firstALTResult->score > maxKForSameAlignment) { // affine gap may produce a better alignment
            firstALTResult->usedAffineGapScoring = true;
            scoreLocationWithAffineGap(read, firstALTResult->direction, firstALTResult->origLocation,
                firstALTResult->seedOffset, scoreLimitForCandidateALT, &firstALTResult->score, &firstALTResult->matchProbability,
                &genomeOffset, &firstALTResult->basesClippedBefore, &firstALTResult->basesClippedAfter, &firstALTResult->agScore);

            if (firstALTResult->score != ScoreAboveLimit) {
                firstALTResult->location = firstALTResult->origLocation + genomeOffset;
            }
            else {
                firstALTResult->status = NotFound;
            }

#if _DEBUG
            if (_DumpAlignments) {
                fprintf(stderr, "Affine gap scored read ALT at %s:%llu score %d, agScore %d\n",
                    genome->getContigAtLocation(firstALTResult->location)->name, firstALTResult->location - genome->getContigAtLocation(firstALTResult->location)->beginningLocation,
                    firstALTResult->score, firstALTResult->agScore
                );
            }
#endif  // _DEBUG

        }
    }

    if (result->status == NotFound || result->score > MAX_K - 1) {
        //
        // Found nothing from the aligner if read is unmapped.
        //
        result->location = InvalidGenomeLocation;
        result->mapq = 0;
        result->score = ScoreAboveLimit;
        result->status = NotFound;
        result->clippingForReadAdjustment = 0;
        result->usedAffineGapScoring = false;
        result->basesClippedBefore = 0;
        result->basesClippedAfter = 0;
        result->agScore = ScoreAboveLimit;
        result->seedOffset = 0;
        result->matchProbability = 0.0;

        firstALTResult->status = NotFound;
#ifdef  _DEBUG
            if (_DumpAlignments) {
                printf("Affine: No sufficiently good alignment found.\n");
            }
#endif  // DEBUG
        return true;
    }

    ScoreSet scoresForAllAlignments;
    ScoreSet scoresForNonAltAlignments;

    //
    // In the beginning we only have the best alignment result in the score set.
    // It is important to initialize the score set here and not before affine gap scoring, since only affine gap does clipping of alignments
    //
    bool nonALTAlignment = (!altAwareness) || !genome->isGenomeLocationALT(result->location);
    scoresForAllAlignments.init(result);
    bool altBestAlignment = false;
    if (firstALTResult->status != NotFound) {
        altBestAlignment = scoresForAllAlignments.updateBestScore(firstALTResult);
    }
    if (nonALTAlignment) {
        scoresForNonAltAlignments.init(result);
    }

    //
    // Update match probability for reads rescored with affine gap
    //
    if (!skipAffineGap) {
        double newProbability = result->matchProbability;
        if (altBestAlignment) { // best result is an ALT result
            double newProbabilityALT = firstALTResult->matchProbability;
            scoresForAllAlignments.updateProbabilityOfAllMatches(oldProbabilityBestResultALT);
            scoresForAllAlignments.updateProbabilityOfBestMatch(newProbabilityALT);
        }
        else {
            scoresForAllAlignments.updateProbabilityOfAllMatches(oldProbabilityBestResult);
            scoresForAllAlignments.updateProbabilityOfBestMatch(newProbability);
        }
        if (nonALTAlignment) {
            scoresForNonAltAlignments.updateProbabilityOfAllMatches(oldProbabilityBestResult);
            scoresForNonAltAlignments.updateProbabilityOfBestMatch(newProbability);
        }
    }

    //
    // Evaluate other alignment candidates with affine gap scoring
    //
    if (nCandidatesForAffineGap > 0 && !skipAffineGap) {

        //
        // Reset score limit
        //
        scoreLimitForCandidate = min(maxK, bestScore) + extraSearchDepth;

        //
        // We sort all promising alignment candidates and score them with affine gap starting with the best one
        //
        qsort(candidatesForAffineGap, nCandidatesForAffineGap, sizeof(*candidatesForAffineGap), SingleAlignmentResult::compareByScore);

        for (int i = 0; i < nCandidatesForAffineGap; i++) {

            SingleAlignmentResult* candidate = &candidatesForAffineGap[i];
            bool nonALTAlignment = (!altAwareness) || !genome->isGenomeLocationALT(candidate->location);
            double oldProbability = candidate->matchProbability;

            if (!skipAffineGap) {
                candidate->usedAffineGapScoring = true;
                scoreLocationWithAffineGap(read, candidate->direction, candidate->origLocation,
                    candidate->seedOffset, scoreLimitForCandidate, &candidate->score, &candidate->matchProbability,
                    &genomeOffset, &candidate->basesClippedBefore, &candidate->basesClippedAfter, &candidate->agScore);

                if (candidate->score != ScoreAboveLimit && (candidate->score <= MAX_K - 1)) {
                    candidate->location = candidate->origLocation + genomeOffset;

                    //
                    // Do not lower MAPQ if we get the same alignment again
                    //
                    if (result->location == candidate->location) {
                        continue;
                    }

                    //
                    // Update match probabilities for read if better
                    //
                    scoresForAllAlignments.updateProbabilityOfAllMatches(oldProbability);
                    scoresForAllAlignments.updateBestScore(candidate);
                    if (nonALTAlignment) {
                        scoresForNonAltAlignments.updateProbabilityOfAllMatches(oldProbability);
                        scoresForNonAltAlignments.updateBestScore(candidate);
                    } // nonALTAlignment

                    //
                    // Update scoreLimit so that we only look for alignments extraSearchDepth worse than the best
                    //
                    scoreLimitForCandidate = scoreLimit(altAwareness && !nonALTAlignment);
                } // candidate->score != ScoreAboveLimit
            } // If we want to score this candidate with affine gap
        } // for each candidate
    }

    //
    // Emit the final result (i.e., ALT/non-ALT best result and first ALT result, if any)
    //
    ScoreSet* scoreSetToEmit;
    if ((!altAwareness) || scoresForNonAltAlignments.bestScore > scoresForAllAlignments.bestScore + maxScoreGapToPreferNonAltAlignment) {
        scoreSetToEmit = &scoresForAllAlignments;
    }
    else {
        scoreSetToEmit = &scoresForNonAltAlignments;
    }

    scoreSetToEmit->fillInSingleAlignmentResult(result, result->popularSeedsSkipped);
    if (altAwareness && scoreSetToEmit == &scoresForNonAltAlignments &&
        (scoresForAllAlignments.bestScoreGenomeLocation != scoresForNonAltAlignments.bestScoreGenomeLocation))

    {
        _ASSERT(genome->isGenomeLocationALT(scoresForAllAlignments.bestScoreGenomeLocation));
        scoresForAllAlignments.fillInSingleAlignmentResult(firstALTResult, firstALTResult->popularSeedsSkipped);

        firstALTResult->supplementary = true;
    }
    else {
        firstALTResult->status = NotFound;
    }

#ifdef  _DEBUG
    if (_DumpAlignments) printf("Final affine gap result score %d agScore %d MAPQ %d (%e probability of best candidate, %e probability of all candidates, non ALT-aware)  at %s:%llu\n\n",
        result->score, result->agScore, result->mapq, scoresForAllAlignments.probabilityOfBestCandidate, scoresForAllAlignments.probabilityOfAllCandidates,
        genome->getContigAtLocation(result->location)->name, result->location - genome->getContigAtLocation(result->location)->beginningLocation);
#endif  // _DEBUG

    return true;
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
    element->usedAffineGapScoring = false;
    element->basesClippedBefore = 0;
    element->basesClippedAfter = 0;
    element->agScore = 0;

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

        if (NULL != affineGap) {
            // affineGap->~AffineGap();
            affineGap->~AffineGapVectorized();
        }
        if (NULL != reverseAffineGap) {
            // reverseAffineGap->~AffineGap();
            reverseAffineGap->~AffineGapVectorized();
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


        if (NULL != affineGap) {
            delete affineGap;
        }
        if (NULL != reverseAffineGap) {
            delete reverseAffineGap;
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
    bestAGScore = -1;
    direction = FORWARD;
    allExtantCandidatesScored = false;
    matchProbabilityForBestScore = 0;
    usedAffineGapScoring = false;
    basesClippedBefore = 0;
    basesClippedAfter = 0;
    agScore = 0;
    seedOffset = 0;
}

BaseAligner::ScoreSet::ScoreSet() 
{
    init();
}

void BaseAligner::ScoreSet::init()
{
    bestScore = UnusedScoreValue;
    bestScoreGenomeLocation = InvalidGenomeLocation;
    bestScoreOrigGenomeLocation = InvalidGenomeLocation;
    bestScoreDirection = FORWARD;
    bestScoreUsedAffineGapScoring = false;
    bestScoreBasesClippedBefore = 0;
    bestScoreBasesClippedAfter = 0;
    bestScoreAGScore = -1;
    bestScoreSeedOffset = 0;
    bestScoreMatchProbability = 0.0;

    probabilityOfAllCandidates = 0;
    probabilityOfBestCandidate = 0;
}

void BaseAligner::ScoreSet::init(SingleAlignmentResult* result) {
    bestScore = result->score;
    bestScoreGenomeLocation = result->location;
    bestScoreOrigGenomeLocation = result->origLocation;
    bestScoreDirection = result->direction;
    bestScoreUsedAffineGapScoring = result->usedAffineGapScoring;
    bestScoreBasesClippedBefore = result->basesClippedBefore;
    bestScoreBasesClippedAfter = result->basesClippedAfter;
    bestScoreAGScore = result->agScore;
    bestScoreSeedOffset = result->seedOffset;
    bestScoreMatchProbability = result->matchProbability;

    probabilityOfAllCandidates = result->probabilityAllCandidates;
    probabilityOfBestCandidate = result->matchProbability;
}

void BaseAligner::ScoreSet::updateProbabilitiesForNearbyMatch(double probabilityOfMatchBeingReplaced)
{
    probabilityOfAllCandidates = __max(0.0, probabilityOfAllCandidates - probabilityOfMatchBeingReplaced);
}

void BaseAligner::ScoreSet::updateProbabilitiesForNewMatch(double newProbability, double matchProbabilityOfNearbyMatch)
{
    probabilityOfAllCandidates = __max(0.0, probabilityOfAllCandidates - matchProbabilityOfNearbyMatch); // need the max due to floating point lossage.
    probabilityOfAllCandidates += newProbability; // Don't combine this with the previous line, it introduces floating point unhappiness.
}

void BaseAligner::ScoreSet::updateBestScore(
                                GenomeLocation genomeLocation,
                                GenomeLocation origGenomeLocation,
                                unsigned score, 
                                bool useAffineGap, 
                                int agScore, 
                                double matchProbability, 
                                unsigned int &lvScoresAfterBestFound,
                                BaseAligner::HashTableElement* elementToScore,
                                SingleAlignmentResult* secondaryResults, 
                                _int64* nSecondaryResults, 
                                _int64 secondaryResultBufferSize,
                                bool anyNearbyCandidatesAlreadyScored,
                                int maxEditDistanceForSecondaryResults, 
                                bool* overflowedSecondaryBuffer,
                                _int64 maxCandidatesForAffineGapBufferSize,
                                _int64* nCandidatesForAffineGap,
                                SingleAlignmentResult* candidatesForAffineGap,
                                unsigned extraSearchDepth)
{
    bool seenNewBestScore = (agScore > bestScoreAGScore) || (bestScoreAGScore == agScore && matchProbability > probabilityOfBestCandidate);

    if (seenNewBestScore) {
        if (bestScore >= score) {
            //
            // If we're tracking secondary alignments, put the old best score in as a new secondary alignment
            //
            if (NULL != secondaryResults && (int)(bestScore - score) <= maxEditDistanceForSecondaryResults) { // bestScore is initialized to UnusedScoreValue, which is large, so this won't fire if this is the first candidate
                if (secondaryResultBufferSize <= *nSecondaryResults) {
                    *overflowedSecondaryBuffer = true;
                    return;
                }

                SingleAlignmentResult* result = &secondaryResults[*nSecondaryResults];
                result->direction = bestScoreDirection;
                result->location = bestScoreGenomeLocation;
                result->origLocation = bestScoreOrigGenomeLocation;
                result->mapq = 0;
                result->score = bestScore;
                result->status = MultipleHits;
                result->clippingForReadAdjustment = 0;
                result->usedAffineGapScoring = bestScoreUsedAffineGapScoring;
                result->basesClippedBefore = bestScoreBasesClippedBefore;
                result->basesClippedAfter = bestScoreBasesClippedAfter;
                result->agScore = bestScoreAGScore;
                result->matchProbability = bestScoreMatchProbability;
                result->seedOffset = bestScoreSeedOffset;

                _ASSERT(result->score != ScoreAboveLimit);

                (*nSecondaryResults)++;
            }

            if (NULL != candidatesForAffineGap && bestScore >= score && (int)(bestScore - score) <= extraSearchDepth) { // bestScore is initialized to UnusedScoreValue, which is large, so this won't fire if this is the first candidate
                if (*nCandidatesForAffineGap >= maxCandidatesForAffineGapBufferSize) {
                    *nCandidatesForAffineGap = maxCandidatesForAffineGapBufferSize + 1;
                    return;
                }

                SingleAlignmentResult* result = &candidatesForAffineGap[*nCandidatesForAffineGap];
                result->direction = bestScoreDirection;
                result->location = bestScoreGenomeLocation;
                result->origLocation = bestScoreOrigGenomeLocation;
                result->mapq = 0;
                result->score = bestScore;
                result->status = MultipleHits;
                result->clippingForReadAdjustment = 0;
                result->usedAffineGapScoring = bestScoreUsedAffineGapScoring;
                result->basesClippedBefore = bestScoreBasesClippedBefore;
                result->basesClippedAfter = bestScoreBasesClippedAfter;
                result->agScore = bestScoreAGScore;
                result->matchProbability = bestScoreMatchProbability;
                result->seedOffset = bestScoreSeedOffset;

                _ASSERT(result->score != ScoreAboveLimit);

                (*nCandidatesForAffineGap)++;
            }
        }

        bestScore = score;
        bestScoreAGScore = agScore;
        probabilityOfBestCandidate = matchProbability;
        _ASSERT(probabilityOfBestCandidate <= probabilityOfAllCandidates);
        bestScoreGenomeLocation = genomeLocation;
        bestScoreOrigGenomeLocation = origGenomeLocation;
        bestScoreDirection = elementToScore->direction;
        bestScoreUsedAffineGapScoring = elementToScore->usedAffineGapScoring;
        bestScoreBasesClippedBefore = elementToScore->basesClippedBefore;
        bestScoreBasesClippedAfter = elementToScore->basesClippedAfter;
        bestScoreSeedOffset = elementToScore->seedOffset;
        bestScoreMatchProbability = elementToScore->matchProbabilityForBestScore;

        lvScoresAfterBestFound = 0;
    } else {
        //
        // If this is close enough, record it as a secondary alignment.
        //
        if (-1 != maxEditDistanceForSecondaryResults && NULL != secondaryResults && (int)(bestScore - score) <= maxEditDistanceForSecondaryResults && score != ScoreAboveLimit && bestScore >= score) {
            if (secondaryResultBufferSize <= *nSecondaryResults) {
                *overflowedSecondaryBuffer = true;
                return;
            }

            SingleAlignmentResult* result = &secondaryResults[*nSecondaryResults];
            result->direction = elementToScore->direction;
            result->location = genomeLocation;
            result->origLocation = origGenomeLocation;
            result->mapq = 0;
            result->score = score;
            result->status = MultipleHits;
            result->clippingForReadAdjustment = 0;
            result->usedAffineGapScoring = elementToScore->usedAffineGapScoring;
            result->basesClippedBefore = elementToScore->basesClippedBefore;
            result->basesClippedAfter = elementToScore->basesClippedAfter;
            result->agScore = elementToScore->agScore;
            result->seedOffset = elementToScore->seedOffset;
            result->matchProbability = elementToScore->matchProbabilityForBestScore;

            _ASSERT(result->score != ScoreAboveLimit);

            (*nSecondaryResults)++;
        }

        if (NULL != candidatesForAffineGap && (int)(bestScore - score) <= extraSearchDepth && score != ScoreAboveLimit && bestScore >= score) {
            if (*nCandidatesForAffineGap >= maxCandidatesForAffineGapBufferSize) {
                *nCandidatesForAffineGap = maxCandidatesForAffineGapBufferSize + 1;
                return;
            }

            SingleAlignmentResult* result = &candidatesForAffineGap[*nCandidatesForAffineGap];
            result->direction = elementToScore->direction;
            result->location = genomeLocation;
            result->origLocation = origGenomeLocation;
            result->mapq = 0;
            result->score = score;
            result->status = MultipleHits;
            result->clippingForReadAdjustment = 0;
            result->usedAffineGapScoring = elementToScore->usedAffineGapScoring;
            result->basesClippedBefore = elementToScore->basesClippedBefore;
            result->basesClippedAfter = elementToScore->basesClippedAfter;
            result->agScore = elementToScore->agScore;
            result->seedOffset = elementToScore->seedOffset;
            result->matchProbability = elementToScore->matchProbabilityForBestScore;

            _ASSERT(result->score != ScoreAboveLimit);

            (*nCandidatesForAffineGap)++;
        }
    }
} // updateBestAndSecondBestScores

void BaseAligner::ScoreSet::fillInSingleAlignmentResult(SingleAlignmentResult* result, int popularSeedsSkipped) {
    result->agScore = bestScoreAGScore;
    result->basesClippedAfter = bestScoreBasesClippedAfter;
    result->basesClippedBefore = bestScoreBasesClippedBefore;
    result->clippingForReadAdjustment = 0;  // This isn't filled in by align()
    result->direction = bestScoreDirection;
    result->location = bestScoreGenomeLocation;
    result->origLocation = bestScoreOrigGenomeLocation;
    result->mapq = computeMAPQ(probabilityOfAllCandidates, probabilityOfBestCandidate, bestScore, popularSeedsSkipped);
    result->score = bestScore;
    result->usedAffineGapScoring = bestScoreUsedAffineGapScoring;
    result->seedOffset = bestScoreSeedOffset;
    result->matchProbability = bestScoreMatchProbability;
    result->popularSeedsSkipped = popularSeedsSkipped;

    if (result->mapq >= MAPQ_LIMIT_FOR_SINGLE_HIT) {
        result->status = SingleHit;
    }
    else {
        result->status = MultipleHits;
    }
    result->probabilityAllCandidates = probabilityOfAllCandidates;
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
                unsigned seedLen, unsigned numSeedsFromCommandLine, double seedCoverage, int maxSecondaryAlignmentsPerContig, unsigned extraSearchDepth)
{
    unsigned maxSeedsToUse;
    if (0 != numSeedsFromCommandLine) {
        maxSeedsToUse = numSeedsFromCommandLine;
    } else {
        maxSeedsToUse = (unsigned)(maxReadSize * seedCoverage / seedLen);
    }
    size_t candidateHashTablesSize = (maxHitsToConsider * maxSeedsToUse * 3)/2;    // *1.5 for hash table slack
    size_t hashTableElementPoolSize = (_int64)maxHitsToConsider * maxSeedsToUse * 2 ;   // *2 for RC
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
        // AffineGap<>::getBigAllocatorReservation()                       + 
        // AffineGap<-1>::getBigAllocatorReservation()                     + // our AffineGap objects
        AffineGapVectorized<>::getBigAllocatorReservation()             + 
        AffineGapVectorized<-1>::getBigAllocatorReservation()           + // our AffineGap objects
        sizeof(char) * maxReadSize * 2                                  + // rcReadData
        sizeof(char) * maxReadSize * 4 + 2 * MAX_K                      + // reversed read (both)
        sizeof(BYTE) * ((_int64)maxReadSize + 7 + 128) / 8              + // seed used
        sizeof(HashTableElement) * hashTableElementPoolSize             + // hash table element pool
        sizeof(HashTableAnchor) * candidateHashTablesSize * 2           + // candidate hash table (both)
        sizeof(HashTableElement) * ((_int64)maxSeedsToUse + 1)          + // weight lists
        sizeof(unsigned) * extraSearchDepth;                              // hitCountByExtraSearchDepth
}

    void 
BaseAligner::finalizeSecondaryResults(
    Read                    *read,
    SingleAlignmentResult   *primaryResult,
    _int64                  *nSecondaryResults,                     // in/out
    SingleAlignmentResult   *secondaryResults,
    _int64                   maxSecondaryResults,
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

    _ASSERT(bestScore == primaryResult->score);

    primaryResult->scorePriorToClipping = primaryResult->score;

    if (!ignoreAlignmentAdjustmentsForOm) {
        //
        // Start by adjusting the alignments for the primary and secondary reads, since that can affect their score
        // and hence whether they should be kept.
        //
        alignmentAdjuster.AdjustAlignment(read, primaryResult);
        if (primaryResult->status != NotFound) {
            bestScore = primaryResult->score;
        } else {
            bestScore = TooBigScoreValue;
        }

        for (int i = 0; i < *nSecondaryResults; i++) {
            secondaryResults[i].scorePriorToClipping = secondaryResults[i].score;
            alignmentAdjuster.AdjustAlignment(read, &secondaryResults[i]);
            if (secondaryResults[i].status != NotFound) {
                bestScore = __min(bestScore, secondaryResults[i].score);
            }
        }
    }

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
            if (ignoreAlignmentAdjustmentsForOm) {
                secondaryResults[i].scorePriorToClipping = secondaryResults[i].score;
            }

            secondaryResults[i].supplementary = altAwareness && genome->isGenomeLocationALT(secondaryResults[i].location);
            i++;
        }
    }

    if (maxSecondaryAlignmentsPerContig > 0 && primaryResult->status != NotFound) {
        //
        // Run through the results and count the number of results per contig, to see if any of them are too big.
        //

        bool anyContigHasTooManyResults = false;

        int primaryResultContigNum = genome->getContigNumAtLocation(primaryResult->location);
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
} // BaseAligner::finalizeSecondaryResults

    int
BaseAligner::scoreLimit(bool forALT) 
{
    if (noUkkonen) {
        return maxK + extraSearchDepth; // We're testing the value of truncating our searches by not doing so.
    }

    if (forALT) {
        return __min(maxK + extraSearchDepth, __min(scoresForAllAlignments.bestScore + extraSearchDepth, scoresForNonAltAlignments.bestScore - maxScoreGapToPreferNonAltAlignment));
    }

    return __min(maxK + extraSearchDepth, __min(scoresForAllAlignments.bestScore + maxScoreGapToPreferNonAltAlignment, scoresForNonAltAlignments.bestScore + extraSearchDepth));
} // BaseAligner::scoreLimit
