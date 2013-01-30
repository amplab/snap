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
#include "Hamming.h"

using std::min;

#ifdef TRACE_ALIGNER
#define TRACE printf
#else
#define TRACE(...) {}
#endif


BaseAligner::BaseAligner(
    GenomeIndex    *i_genomeIndex, 
    unsigned        i_confDiff, 
    unsigned        i_maxHitsToConsider, 
    unsigned        i_maxK,
    unsigned        i_maxReadSize,
    unsigned        i_maxSeedsToUse,
    unsigned        i_adaptiveConfDiffThreshold,
    LandauVishkin<>*i_landauVishkin,
    SimilarityMap  *i_similarityMap,
    AlignerStats   *i_stats,
    BigAllocator   *allocator) : 
        genomeIndex(i_genomeIndex), confDiff(i_confDiff), maxHitsToConsider(i_maxHitsToConsider), maxK(i_maxK), 
        maxReadSize(i_maxReadSize), maxSeedsToUse(i_maxSeedsToUse), readId(-1),
        adaptiveConfDiffThreshold(i_adaptiveConfDiffThreshold),
        similarityMap(i_similarityMap), explorePopularSeeds(false), stopOnFirstHit(false), stats(i_stats)
/*++

Routine Description:

    Constructor for the BaseAligner class.  Aligners align reads against an indexed genome.

Arguments:

    i_genomeIndex       - The index against which to do the alignments
    i_confDiff          - The string difference between two matches necessary to believe they're really different
    i_maxHitsToConsider - The maximum number of hits to use from a seed lookup.  Any lookups that return more
                          than this are ignored.
    i_maxK              - The largest string difference to consider for any comparison.
    i_maxReadSize       - Bound on the number of bases in any read.  There's no reason to make it tight, it just affects a little memory allocation.
    i_maxSeedsToUse     - The maximum number of seeds to use when aligning any read (not counting ones ignored because they resulted in too many
                          hits).  Once we've looked up this many seeds, we just score what we've got.
    i_adaptiveConfDiffThreshold - the number of hash table hits larger than maxHitsToConsider beyond which we effectively increase confDiff
    i_landauVishkin     - an externally supplied LandauVishkin string edit distance object.  This is useful if we're expecting repeated computations and use the LV cache.
    i_similarityMap     - a similarity map used to handle repetitive regions.  Optional.
    i_stats             - an object into which we report out statistics
    allocator           - an allocator that's used to allocate our local memory.  This is useful for TLB optimization.  If this is supplied, the caller
                          is responsible for deallocation, we'll not deallocate any dynamic memory in our destructor.
 
 --*/
{
#if     MAINTAIN_HISTOGRAMS
    lvHistogram = new Histogram(20, true);
    lookupHistogram = new Histogram(20, true);
    lvHistogramForMulti = new Histogram(20, true);
    lvCountWhenBestFound = new Histogram(maxHitsToConsider*4,false);
#endif  // MAINTAIN_HISTOGRAMS

    hadBigAllocator = allocator != NULL;

    nHashTableLookups = 0;
    nLocationsScored = 0;
    nHitsIgnoredBecauseOfTooHighPopularity = 0;
    nReadsIgnoredBecauseOfTooManyNs = 0;
    nIndelsMerged = 0;

    genome = genomeIndex->getGenome();
    seedLen = genomeIndex->getSeedLength();

    probDistance = new ProbabilityDistance(SNP_PROB, GAP_OPEN_PROB, GAP_EXTEND_PROB);  // Match Mason

#if defined(USE_BOUNDED_STRING_DISTANCE)
    boundedStringDist = new BoundedStringDistance<>(2, 2, SNP_PROB, GAP_OPEN_PROB, GAP_EXTEND_PROB);
    landauVishkin = NULL;
    reverseLandauVishkin = NULL;
#else // Landau-Vishkin

    if (i_landauVishkin == NULL) {
        if (allocator) {
            landauVishkin = new (allocator) LandauVishkin<>;
        } else {
            landauVishkin = new LandauVishkin<>;
        }
        ownLandauVishkin = true;
    } else {
        landauVishkin = i_landauVishkin;
        ownLandauVishkin = false;
    }

    if (allocator) {
        reverseLandauVishkin = new (allocator) LandauVishkin<-1>;
    } else {
        reverseLandauVishkin = new LandauVishkin<-1>;
    }
#endif  // BSD or LV

    unsigned nCandidates = __min(maxHitsToConsider * (maxReadSize - seedLen + 1) * 2, 200000);  // *2 is for reverse complement
    candidateHashTablesSize = (maxHitsToConsider * maxSeedsToUse * 3)/2;    // *1.5 for hash table slack
    hashTableElementPoolSize = maxHitsToConsider * maxSeedsToUse * 2 ;   // *2 for RC
   
    if (allocator) {
        candidates = (Candidate *)allocator->allocate(sizeof(Candidate) * nCandidates);
    } else {
        candidates = (Candidate *)BigAlloc(sizeof(Candidate) * nCandidates);
    }
    for (unsigned i = 0 ; i < nCandidates; i++) {
        candidates[i].init();
    }

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

    reversedRead[RC] = reversedRead[FORWARD] + maxReadSize;
    reversedGenomeData = reversedRead[FORWARD] + 3 * maxReadSize + MAX_K;   // WARNING: This buffer extends in both directions from the pointer.
    

    rcTranslationTable['A'] = 'T';
    rcTranslationTable['G'] = 'C';
    rcTranslationTable['C'] = 'G';
    rcTranslationTable['T'] = 'A';
    rcTranslationTable['N'] = 'N';

    for (unsigned i = 0; i < 256; i++) {
        nTable[i] = 0;
    }

    nTable['N'] = 1;

    if (allocator) {
        seedUsed = (BYTE *)allocator->allocate((sizeof(BYTE) * (maxReadSize + 7 + 128) / 8));    // +128 to make sure it extends at both 
    } else {
        seedUsed = (BYTE *)BigAlloc((sizeof(BYTE) * (maxReadSize + 7 + 128) / 8));    // +128 to make sure it extends at both  
    }

    seedUsedAsAllocated = seedUsed; // Save the pointer for the delete.
    seedUsed += 8;  // This moves the pointer up an _int64, so we now have the appropriate before buffer.

    nUsedHashTableElements = 0;

    /*
    if (maxK*2 > 64) {
        fprintf(stderr,"The HashTableElement struct uses a 64 bit bitvector for used candidates.  As a result, you can't 
        exit(1);
    }
    */

    if (allocator) {
        candidateHashTable[FORWARD] = (HashTableAnchor *)allocator->allocate(sizeof(HashTableAnchor) * candidateHashTablesSize);
        candidateHashTable[RC] = (HashTableAnchor *)allocator->allocate(sizeof(HashTableAnchor) * candidateHashTablesSize);
        weightLists = (HashTableElement *)allocator->allocate(sizeof(HashTableElement) * (maxSeedsToUse + 1));
        hashTableElementPool = (HashTableElement *)allocator->allocate(sizeof(HashTableElement) * hashTableElementPoolSize); // Allocte last, because it's biggest and usually unused.  This puts all of the commonly used stuff into one large page.
    } else {
        candidateHashTable[FORWARD] = (HashTableAnchor *)BigAlloc(sizeof(HashTableAnchor) * candidateHashTablesSize);
        candidateHashTable[RC] = (HashTableAnchor *)BigAlloc(sizeof(HashTableAnchor) * candidateHashTablesSize);
        weightLists = (HashTableElement *)BigAlloc(sizeof(HashTableElement) * (maxSeedsToUse + 1));
        hashTableElementPool = (HashTableElement *)BigAlloc(sizeof(HashTableElement) * hashTableElementPoolSize);
    }

      if (NULL != stats) {
        int onTheStack = 42;
        stats->threadEntry->candidateEntries = hashTableElementPool;
        stats->threadEntry->alignerObject = this;
        stats->threadEntry->hashAnchor[0] = candidateHashTable[0];
        stats->threadEntry->hashAnchor[1] = candidateHashTable[1];
        stats->threadEntry->stackPointer = &onTheStack;
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

    AlignmentResult
BaseAligner::AlignRead(Read *inputRead, unsigned *genomeLocation, Direction *hitDirection, int *finalScore, int *mapq)
{
    return AlignRead(inputRead, genomeLocation, hitDirection, finalScore, mapq, 0, 0, false, 0, NULL, NULL, NULL, NULL);
}

    AlignmentResult
BaseAligner::AlignRead(
    Read      *inputRead,
    unsigned  *genomeLocation,
    Direction *hitDirection,
    int       *finalScore,
    int       *mapq,
    unsigned   searchRadius,
    unsigned   searchLocation,
    Direction  searchDirection)
{
    return AlignRead(inputRead, genomeLocation, hitDirection, finalScore, mapq,
        searchRadius, searchLocation, searchDirection, 0, NULL, NULL, NULL, NULL);
}

#ifdef  _DEBUG
bool _DumpAlignments = false;
#endif  // _DEBUG

    AlignmentResult
BaseAligner::AlignRead(
    Read      *inputRead,
    unsigned  *genomeLocation,
    Direction *hitDirection,
    int       *finalScore,
    int       *mapq,
    unsigned   searchRadius,
    unsigned   searchLocation,
    Direction  searchDirection,
    int        maxHitsToGet,
    int       *multiHitsFound,
    unsigned  *multiHitLocations,
    Direction *multiHitDirections,
    int       *multiHitScores)  
/*++

Routine Description:

    Align a particular read, possibly constraining the search around a given location.

Arguments:

    inputRead           - the read to align
    genomeLocation      - if this aligned to a SingleHit, the 0-based offset in the genome that this hit.  The aligner treats the entire
                          genome as a single string, even though it's really a set of chrosomes.  It just makes the code simpler.
                          The caller can convert to chromosome+offset by looking up the piece boudaries using the Genome object.
    hitDirection        - the aligner tries to align both the given read and its reverse complement.  If we found a SingleHit this
                          is set to indicate whether that hit was on the forward or reverse complement.
    finalScore          - if a single or confident hit, this is the score of the hit (i.e., the LV distance from the read)
    mapq                - return the mapping quality
    searchRadius        - if non-zero, constrain the search to this distance around searchLocation, in orientation searchRC
    searchLocation      - location to search around if searchRadius is given
    searchDirection     - whether to search in forward or reverse complement orientation if searchRadius is given
    maxHitsToGet        - if greater than 0, output up to this many hits within confDiff of the best (if any) in multiHitLocation
                          writing their count in multiHitsFound, instead of returning MultipleHits immediately
    multiHitsFound      - output parameter for number of alternative hits found if maxHitsToGet is true
    multiHitLocations   - output parameter for locations of alternative hits found if maxHitsToGet is true

Return Value:

    ConfidentHit, SingleHit, MultiHit or NotFound depending on how the alignment went.

--*/
{
    //
    // mapq and finalScore are optional parameters.  Rather than checking all over the code to see if
    // they're null, we just point them at locals if they're not passed in.
    //
    int unusedMapq;
    int unusedFinalScore;
    firstPassSeedsNotSkipped[FORWARD] = firstPassSeedsNotSkipped[RC] = 0;
    smallestSkippedSeed[FORWARD] = smallestSkippedSeed[RC] = 0xffffffff;
    biggestClusterScored = 1;
    highestWeightListChecked = 0;
    usedHammingThisAlignment = false;

    if (NULL == mapq) {
        mapq = &unusedMapq;
    }

    if (NULL == finalScore) {
        finalScore = &unusedFinalScore;
    }

    *genomeLocation = 0xFFFFFFFF; // Value to return if we don't find a location.
    *hitDirection = FORWARD;             // So we deterministically print the read forward in this case.
    *finalScore = UnusedScoreValue;

    unsigned lookupsThisRun = 0;

    popularSeedsSkipped = 0;

    // Range of genome locations to search in.
    unsigned minLocation = 0;
    unsigned maxLocation = 0xFFFFFFFF;
    if (searchRadius != 0) {
        minLocation = (searchLocation > searchRadius) ? searchLocation - searchRadius : 0;
        maxLocation = (searchLocation < 0xFFFFFFFF - searchRadius) ? searchLocation + searchRadius : 0xFFFFFFFF;
    }

    // If asked to return the locations of multiple hits, make sure we get sensible counts and results
    if (maxHitsToGet > 0) {
        memset(hitCount, 0, MAX_K * sizeof(unsigned));
        *multiHitsFound = 0;
    }

    AlignmentResult finalResult;

    //
    // A bitvector for used seeds, indexed on the starting location of the seed within the read.
    //
    if (inputRead->getDataLength() > maxReadSize) {
        fprintf(stderr,"BaseAligner:: got too big read (%d > %d)", inputRead->getDataLength(),maxReadSize);
        exit(1);
    }

    if ((int)inputRead->getDataLength() < seedLen) {
        //
        // Too short to have any seeds, it's hopeless.
        //
        return NotFound;
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
        return NotFound;
    }

    //
    // Block off any seeds that would contain an N.
    //
    //if (countOfNs > 0) {
    //    int minSeedToConsiderNing = 0; // In English, any word can be verbed. Including, apparently, "N."
    //    for (int i = 0; i < (int) readLen; i++) {
    //        if (readData[i] == 'N') {
    //            int limit = __min(i + seedLen - 1, readLen-1);
    //            for (int j = __max(minSeedToConsiderNing, i - (int) seedLen + 1); j <= limit; j++) {
    //                SetSeedUsed(j);
    //            }
    //            minSeedToConsiderNing = limit+1;
    //            if (minSeedToConsiderNing >= (int) readLen)
    //                break;
    //        }
    //    }
    //}

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

    if (maxHitsToGet > 0) {
        scoreLimit = maxK + 3;
    } else {
        scoreLimit = maxK + 2; // For MAPQ computation
        //scoreLimit = maxK + confDiff - 1;
    }

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
                    &finalResult,
                    finalScore,
                    genomeLocation,
                    hitDirection,
                    candidates,
                    maxHitsToGet,
                    mapq);
#if     MAINTAIN_HISTOGRAMS
                lvHistogram->addToCount(lvScores,lvScores);
                lookupHistogram->addToCount(lookupsThisRun,lookupsThisRun);
                if (MultipleHits == finalResult) {
                    lvHistogramForMulti->addToCount(lvScores,lvScores);
                } else if (SingleHit == finalResult) {
                    lvCountWhenBestFound->addToCount(lvScores-lvScoresAfterBestFound);
                }
#endif  // MAINTAIN_HISTOGRAMS

                fillHitsFound(maxHitsToGet, multiHitsFound,
                              multiHitLocations, multiHitDirections, multiHitScores);
#ifdef  _DEBUG
                if (_DumpAlignments) printf("\tFinal result score %d MAPQ %d at %u\n", *finalScore, *mapq, *genomeLocation);
#endif  // _DEBUG
                return finalResult;
            }
            nextSeedToTest = getWrappedNextSeedToTest(wrapCount);

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

        unsigned        nHits[NUM_DIRECTIONS];      // Number of times this seed hits in the genome
        const unsigned  *hits[NUM_DIRECTIONS];      // The actual hits (of size nHits)
    
        unsigned minSeedLoc = (minLocation < readLen ? 0 : minLocation - readLen);
        unsigned maxSeedLoc = (maxLocation > 0xFFFFFFFF - readLen ? 0xFFFFFFFF : maxLocation + readLen);
        genomeIndex->lookupSeed(seed, minSeedLoc, maxSeedLoc, &nHits[0], &hits[0], &nHits[1], &hits[1]);
        nHashTableLookups++;
        lookupsThisRun++;


#ifdef  _DEBUG
        if (_DumpAlignments) {
            printf("\tSeed offset %2d, %4d hits, %4d rcHits.", nextSeedToTest, nHits[0], nHits[1]);
            for (int rc = 0; rc < 2; rc++) {
                for (unsigned i = 0; i < __min(nHits[rc], 5); i++) {
                    printf(" %sHit at %9u.", rc == 1 ? "RC " : "", hits[rc][i]);
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
            if (searchRadius != 0 && searchDirection != direction) {
                //
                // We're looking only for hits in the other direction.
                //
                continue;
            }

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

                for (unsigned i = 0 ; i < min(nHits[direction], maxHitsToConsider); i++) {
                    //
                    // Find the genome location where the beginning of the read would hit, given a match on this seed.
                    //

                    unsigned genomeLocationOfThisHit = hits[direction][i] - offset;
                    if (genomeLocationOfThisHit < minLocation || genomeLocationOfThisHit > maxLocation) {
                        continue;
                    }
    
                    Candidate *candidate = NULL;
                    HashTableElement *hashTableElement;

                    findCandidate(genomeLocationOfThisHit, direction, &candidate, &hashTableElement);
                    if (NULL != hashTableElement) {
                        incrementWeight(hashTableElement);
                        candidate->seedOffset = offset;
                        _ASSERT((unsigned)candidate->seedOffset <= readLen - seedLen);
                    } else if (lowestPossibleScoreOfAnyUnseenLocation[direction] <= scoreLimit) {
                        _ASSERT(offset <= readLen - seedLen);
                        allocateNewCandidate(genomeLocationOfThisHit, direction, lowestPossibleScoreOfAnyUnseenLocation[direction],
                                offset, &candidate, &hashTableElement);
                    }
                }
                nSeedsApplied[direction]++;
                appliedEitherSeed = true;
            } // not too popular
        }   // directions 

        //
        // Move us along.
        //
        nextSeedToTest += seedLen;

        if (appliedEitherSeed) {
            //
            // And finally, try scoring.
            //
            if (score(  false,
                        read,
                        &finalResult,
                        finalScore,
                        genomeLocation,
                        hitDirection,
                        candidates,
                        maxHitsToGet,
                        mapq)) {
#if     MAINTAIN_HISTOGRAMS
                lvHistogram->addToCount(lvScores,lvScores);
                lookupHistogram->addToCount(lookupsThisRun,lookupsThisRun);
                if (MultipleHits == finalResult) {
                    lvHistogramForMulti->addToCount(lvScores,lvScores);
                } else if (SingleHit == finalResult) {
                    lvCountWhenBestFound->addToCount(lvScores-lvScoresAfterBestFound);
                }
#endif  // MAINTAIN_HISTOGRAMS

                fillHitsFound(maxHitsToGet, multiHitsFound,
                              multiHitLocations, multiHitDirections, multiHitScores);
#ifdef  _DEBUG
                if (_DumpAlignments) printf("\tFinal result score %d MAPQ %d at %u\n", *finalScore, *mapq, *genomeLocation);
#endif  // _DEBUG
                return finalResult;
            }
        }
    }

    //
    // Do the best with what we've got.
    //
#ifdef TRACE_ALIGNER
    printf("Calling score with force=true because we ran out of seeds\n");
#endif
    score(  true,
            read,
            &finalResult,
            finalScore,
            genomeLocation,
            hitDirection,
            candidates,
            maxHitsToGet,
            mapq);
#if     MAINTAIN_HISTOGRAMS
        lvHistogram->addToCount(lvScores,lvScores);
        lookupHistogram->addToCount(lookupsThisRun,lookupsThisRun);
        if (MultipleHits == finalResult) {
            lvHistogramForMulti->addToCount(lvScores,lvScores);
        } else if (SingleHit == finalResult) {
            lvCountWhenBestFound->addToCount(lvScores-lvScoresAfterBestFound);
        }
#endif  // MAINTAIN_HISTOGRAMS
    
    fillHitsFound(maxHitsToGet, multiHitsFound,
                  multiHitLocations, multiHitDirections, multiHitScores);
#ifdef  _DEBUG
    if (_DumpAlignments) printf("\tFinal result score %d MAPQ %d at %u\n", *finalScore, *mapq, *genomeLocation);
#endif  // _DEBUG

    return finalResult;
}

    bool
BaseAligner::score(
    bool             forceResult,
    Read            *read[NUM_DIRECTIONS],
    AlignmentResult *result,
    int             *finalScore,
    unsigned        *singleHitGenomeLocation,
    Direction       *hitDirection,
    Candidate       *candidates,
    unsigned         maxHitsToGet,
    int             *mapq)
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
    maxHitsToGet                            - whether the user wants to get multiple non-confident hits in the MultiHit case
    mapq                                    - returns the map quality if we've reached a final result

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
        // to try and every one of them generated too many hits to process.  Declare
        // a multi hit and give up.
        //
        _ASSERT(forceResult);
        *result = MultipleHits;
        *mapq = 0;
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

    // Return early if we haven't tried enough seeds
    //if (__max(nSeedsApplied, nRCSeedsApplied) < 4 && !forceResult) {
    //    return false;
    //}

    unsigned weightListToCheck = highestUsedWeightList;

    do {
        HashTableElement *elementToScore;

        _ASSERT(weightListToCheck <= maxSeedsToUse);
        //
        // Grab the next element to score, and score it.
        //

        while (weightListToCheck > 0 && weightLists[weightListToCheck].weightNext == &weightLists[weightListToCheck]) {
            weightListToCheck--;
            highestUsedWeightList = weightListToCheck;
        }

        _ASSERT(weightListToCheck <= maxSeedsToUse);

        if (__min(lowestPossibleScoreOfAnyUnseenLocation[FORWARD],lowestPossibleScoreOfAnyUnseenLocation[RC]) > scoreLimit || forceResult) {
            if (weightListToCheck == 0 /*|| weightListToCheck == 1 && forceResult && lvScores > 20*/) {
                //
                // We've scored all live candidates and excluded all non-candidates, or we've checked enough that we've hit the cutoff.  We have our
                // answer.
                //
                int realConfDiff = confDiff + (popularSeedsSkipped >= adaptiveConfDiffThreshold ? 1 : 0);
                *finalScore = bestScore;
                if (bestScore + realConfDiff <= secondBestScore && bestScore <= maxK) {
                    if (popularSeedsSkipped == 0 && !forceResult) {
                        *result = CertainHit;
                    } else {
                        *result = SingleHit;
                    }
                    *singleHitGenomeLocation = bestScoreGenomeLocation;
                    *mapq = computeMAPQ(probabilityOfAllCandidates, probabilityOfBestCandidate, bestScore, firstPassSeedsNotSkipped,  
                                            smallestSkippedSeed, bestScoreGenomeLocation, popularSeedsSkipped, similarityMap, biggestClusterScored, usedHammingThisAlignment);
                    return true;
                } else if (bestScore > maxK) {
                    // If none of our seeds was below the popularity threshold, report this as MultipleHits; otherwise,
                    // report it as NotFound
                    *result = (nSeedsApplied[FORWARD] == 0 && nSeedsApplied[RC] == 0) ? MultipleHits : NotFound;
                    *mapq = 0;
                    return true;
                } else {
                    *result = MultipleHits;
                    *mapq = computeMAPQ(probabilityOfAllCandidates, probabilityOfBestCandidate, bestScore, firstPassSeedsNotSkipped,  
                                        smallestSkippedSeed, bestScoreGenomeLocation, popularSeedsSkipped, similarityMap, biggestClusterScored, usedHammingThisAlignment);
                    return true;
                }
            }
            //
            // Nothing that we haven't already looked up can possibly be the answer.  Score what we've got and exit.
            //
            forceResult = true;
        } else if (weightListToCheck == 0) {
            _ASSERT(!forceResult);  // This assert is here beacause I'm deleting some code that I believe is dead (and that otherwise looked wrong).  --bb
            //
            // No candidates, look for more.
            //
            return false;
        }

        elementToScore = weightLists[weightListToCheck].weightNext;
        _ASSERT(!elementToScore->allExtantCandidatesScored);
        _ASSERT(elementToScore->candidatesUsed != 0);
        _ASSERT(elementToScore != &weightLists[weightListToCheck]);

        if (elementToScore->lowestPossibleScore <= scoreLimit) {

            unsigned long candidateIndexToScore;
            _uint64 candidatesMask = elementToScore->candidatesUsed;
            while (_BitScanForward64(&candidateIndexToScore,candidatesMask)) {
#ifdef  TIME_STRING_DISTANCE
                _int64 timeInBSD = 0;
#endif  // TIME_STRING_DISTANCE
                bool usedHamming;
                _uint64 candidateBit = ((_uint64)1 << candidateIndexToScore);
                candidatesMask &= ~candidateBit;
                if ((elementToScore->candidatesScored & candidateBit) != 0) {
                    // Already scored it, or marked it as scored due to using ProbabilityDistance
                    continue;
                }

                bool anyNearbyCandidatesAlreadyScored = elementToScore->candidatesScored != 0;

                elementToScore->candidatesScored |= candidateBit;
                Candidate *candidateToScore = &elementToScore->candidates[candidateIndexToScore];
 
                unsigned genomeLocation = elementToScore->baseGenomeLocation + candidateIndexToScore;

                //
                // We're about to run edit distance computations on the genome.  Launch a prefetch for it
                // so that it's in cache when we do (or at least on the way).
                //
                if (doAlignerPrefetch) {
                    genomeIndex->prefetchGenomeData(genomeLocation);
                }

    
                unsigned score = -1;
                double matchProbability;
                unsigned readDataLength = read[elementToScore->direction]->getDataLength();
                unsigned genomeDataLength = readDataLength + MAX_K; // Leave extra space in case the read has deletions
                const char *data = genome->getSubstring(genomeLocation, genomeDataLength);
                if (NULL == data) {
                    //
                    // We're up against the end of a chromosome.  Reduce the extra space enough that it isn't too
                    // long.  We're willing to reduce it to less than the length of a read, because the read could
                    // but up against the end of the chromosome and have insertions in it.
                    //
                    const Genome::Piece *piece = genome->getPieceAtLocation(genomeLocation);
                    
                    unsigned endOffset;
                    if (genomeLocation + readDataLength + MAX_K >= genome->getCountOfBases()) {
                        endOffset = genome->getCountOfBases();
                    } else {
                        const Genome::Piece *nextPiece = genome->getPieceAtLocation(genomeLocation + readDataLength + MAX_K);
                        _ASSERT(NULL != piece && piece->beginningOffset <= genomeLocation && piece != nextPiece);

                        endOffset = nextPiece->beginningOffset;
                    }
                    genomeDataLength = endOffset - genomeLocation - 1;
                    if (genomeDataLength >= readDataLength - MAX_K) {
                        data = genome->getSubstring(genomeLocation, genomeDataLength);
                        _ASSERT(NULL != data);
                    }
                }

                if (data != NULL) {
                    Read *readToScore = read[elementToScore->direction];

                    if (false && 0 == candidatesMask /*&& NULL == nearbyElement*/ && forceResult && elementToScore->weight > 1) {
                        //
                        // No evidence of indels near here and we're done with hash table lookups.  Use Hamming distance.
                        //
#ifdef  TIME_STRING_DISTANCE
                        _int64 startTime = timeInNanos();
#endif  // TIME_STRING_DISTANCE
                        score = ComputeHammingDistance(readToScore->getData(), data, readToScore->getDataLength(), readToScore->getQuality(), scoreLimit, &matchProbability);
#ifdef  TIME_STRING_DISTANCE
                        if (NULL != stats) {
                            stats->hammingNanos += timeInNanos() - startTime;
                            stats->hammingCount++;
                        }
#endif  // TIME_STRING_DISTANCE
                        usedHamming = true;
                        usedHammingThisAlignment = true;
                    } else {
                        usedHamming = false;
                        _ASSERT(candidateToScore->seedOffset + seedLen <= readToScore->getDataLength());
#if defined(USE_PROBABILITY_DISTANCE)
                        score = probDistance->compute(data, readToBeScored->getData(), readToBeScored->getQuality(),
                                readToBeScored->getDataLength(), 6, 12, &matchProbability);
 
                        // Since we're allowing a startShift in ProbabilityDistance, mark nearby locations as scored too
                        for (int shift = 1; shift <= 6; shift++) {
                            elementToScore->candidatesScored |= (candidateBit << shift);
                            elementToScore->candidatesScored |= (candidateBit >> shift);
                        }
                        // Update the biggest cluster scored if this is a sufficiently likely cluster
                        //if (similarityMap != NULL && matchProbability >= __max(probabilityOfBestCandidate * 0.01, 1e-30)) {
                        //    biggestClusterScored = __max(biggestClusterScored, similarityMap->getNumClusterMembers(genomeLocation));
                        //}
#elif defined(USE_BOUNDED_STRING_DISTANCE)
                        int maxStartShift = __min(scoreLimit, 7);
                        for (int shift = 1; shift <= maxStartShift; shift++) {
                            elementToScore->candidatesScored |= (candidateBit << shift);
                            elementToScore->candidatesScored |= (candidateBit >> shift);
                        }
                        //score = boundedStringDist->compute(data, rcRead->getData(), rcRead->getDataLength(),
                        //        maxStartShift, scoreLimit, &matchProbability, rcRead->getQuality());
                        
                        
                        // Compute the distance separately in the forward and backward directions from the seed, to allow
                        // arbitrary offsets at both the start and end but not have to pay the cost of exploring all start
                        // shifts in BoundedStringDistance
                        double matchProb1, matchProb2;
                        int score1, score2;
                        // First, do the forward direction from where the seed aligns to past of it
                        int readLen = readToScore->getDataLength();
                        int seedLen = genomeIndex->getSeedLength();
                        int seedOffset = candidateToScore->seedOffset; // Since the data is reversed
                        int tailStart = seedOffset + seedLen;
#ifdef  TIME_STRING_DISTANCE
                        _int64 bsdStartTime = timeInNanos();
#endif  // TIME_STRING_DISTANCE
                        score1 = boundedStringDist->compute(data + tailStart, readToScore->getData() + tailStart,
                                readLen - tailStart, 0, scoreLimit, &matchProb1, readToScore->getQuality() + tailStart);

                        if (score1 == -1) {
                            score = -1;
                        } else {
                            // The tail of the read matched; now let's reverse the reference genome data and match the head
                            int limitLeft = scoreLimit - score1;

                            for (int i = -MAX_K; i <= seedOffset + MAX_K; i++) {
                                reversedGenomeData[i] = data[seedOffset - 1 - i];
                            }
 
                            // Note that we use the opposite direction read for the quality, since it's the reverse of our direction's
                            score2 = boundedStringDist->compute(reversedGenomeData, reversedRead[elementToScore->direction] + readLen - seedOffset,
                                   seedOffset, 0, limitLeft, &matchProb2, read[OppositeDirection(elementToScore->direction)]->getQuality() + readLen - seedOffset);
                            // TODO: Use endShift to mark scored candidates more correctly and to report real location in SAM

                            if (score2 == -1) {
                                score = -1;
                            } else {
                                score = score1 + score2;
                                // Map probabilities for substrings can be multiplied, but make sure to count seed too
                                matchProbability = matchProb1 * matchProb2 * pow(1 - SNP_PROB, seedLen);
                            }
                        }
#ifdef  TIME_STRING_DISTANCE
                        timeInBSD = timeInNanos() - bsdStartTime;
#endif  // TIME_STRING_DISTANCE
                        
                        if (score != -1) {
                            if (similarityMap != NULL) {
                                biggestClusterScored = __max(biggestClusterScored,
                                        similarityMap->getNumClusterMembers(genomeLocation));
                            }
                        } else {
                            matchProbability = 0;
                        }
#else   // Landau-Vishkin
                        int maxStartShift = __min(scoreLimit, 7);
                        for (int shift = 1; shift <= maxStartShift; shift++) {
                            elementToScore->candidatesScored |= (candidateBit << shift);
                            elementToScore->candidatesScored |= (candidateBit >> shift);
                        }
                        
                        // Compute the distance separately in the forward and backward directions from the seed, to allow
                        // arbitrary offsets at both the start and end but not have to pay the cost of exploring all start
                        // shifts in BoundedStringDistance
                        double matchProb1, matchProb2;
                        int score1, score2;
                        // First, do the forward direction from where the seed aligns to past of it
                        int readLen = readToScore->getDataLength();
                        int seedLen = genomeIndex->getSeedLength();
                        int seedOffset = candidateToScore->seedOffset; // Since the data is reversed
                        int tailStart = seedOffset + seedLen;

                        score1 = landauVishkin->computeEditDistance(data + tailStart, genomeDataLength - tailStart, readToScore->getData() + tailStart, readToScore->getQuality() + tailStart, readLen - tailStart,
                            scoreLimit, &matchProb1);
                        if (score1 == -1) {
                            score = -1;
                        } else {
                            // The tail of the read matched; now let's reverse the reference genome data and match the head
                            int limitLeft = scoreLimit - score1;
                            score2 = reverseLandauVishkin->computeEditDistance(data + seedOffset, seedOffset + MAX_K, reversedRead[elementToScore->direction] + readLen - seedOffset,
                                                                                        read[OppositeDirection(elementToScore->direction)]->getQuality() + readLen - seedOffset, seedOffset, limitLeft, &matchProb2);

                            if (score2 == -1) {
                                score = -1;
                            } else {
                                score = score1 + score2;
                                // Map probabilities for substrings can be multiplied, but make sure to count seed too
                                matchProbability = matchProb1 * matchProb2 * pow(1 - SNP_PROB, seedLen);
                            }
                        }

                        if (score != -1) {
                            if (similarityMap != NULL) {
                                biggestClusterScored = __max(biggestClusterScored,
                                        similarityMap->getNumClusterMembers(genomeLocation));
                            }
                        } else {
                            matchProbability = 0;
                        }

#endif
                    } // If we used Hamming
                } // if we had genome data to compare against
#ifdef TRACE_ALIGNER
                printf("Computing distance at %u (RC) with limit %d: %d (prob %g)\n",
                        genomeLocation, scoreLimit, score, matchProbability);
#endif


#ifdef  _DEBUG
                if (_DumpAlignments) printf("Scored %9u weight %2d limit %d, result %2d %s\n", genomeLocation, elementToScore->weight, scoreLimit, score, elementToScore->direction ? "RC" : "");
#endif  // _DEBUG
                
                if (maxHitsToGet > 0 && score != -1 && hitCount[score] < maxHitsToGet) {
                    // Remember the location of this hit because we don't have enough at this distance
                    hitLocations[score][hitCount[score]] = genomeLocation;
                    hitDirections[score][hitCount[score]] = elementToScore->direction;
                    hitCount[score]++;
                }
                
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

                //
                // Look up the hash table element that's closest to the genomeLocation but that doesn't
                // contain it, to check if this location is already scored.
                //
                // We do this computation in a strange way in order to avoid generating a branch instruction that
                // the processor's branch predictor will get wrong half of the time.  Think about it like this:
                // The genome location lies in a bucket of size 2 * maxMergeDist.  Its offset in the bucket
                // is genomeLocation % (2 * maxMergeDist).  If we take that quantity and integer divide it by
                // maxMergeDist, we get 0 if it's in the first half and 1 if it's in the second.  Double that and subtract
                // one, and you're at the right place with no branches.
                //
                HashTableElement *nearbyElement;
                if (-1 != score) {
                    unsigned nearbyGenomeLocation = genomeLocation + (2*(genomeLocation % (2 * maxMergeDist) / maxMergeDist) - 1) * maxMergeDist;
                    _ASSERT((genomeLocation % (2 * maxMergeDist) >= maxMergeDist ? genomeLocation + maxMergeDist : genomeLocation - maxMergeDist) == nearbyGenomeLocation);   // Assert that the logic in the above comment is right.

                    findElement(nearbyGenomeLocation, elementToScore->direction, &nearbyElement);
                } else {
                    nearbyElement = NULL;
                }

                if (NULL != nearbyElement && nearbyElement->candidatesScored != 0) {
                    if (nearbyElement->bestScore < score || nearbyElement->bestScore == score && nearbyElement->matchProbabilityForBestScore >= matchProbability) {
                        //
                        // Again, this no better than something nearby we already tried.  Give up.
                        //
#ifdef  TIME_STRING_DISTANCE
                        if (NULL != stats && !usedHamming) {
                            stats->nanosTimeInBSD[1][score==-1 ? 0:1] += timeInBSD;
                            stats->BSDCounts[1][score==-1 ? 0:1]++;
                        }
#endif  // TIME_STRING_DISTANCE
                        continue;
                    }
                    anyNearbyCandidatesAlreadyScored = true;
                    probabilityOfAllCandidates = __max(0.0, probabilityOfAllCandidates - nearbyElement->matchProbabilityForBestScore);
                    nearbyElement->matchProbabilityForBestScore = 0;    // keeps us from backing it out twice
                }

                probabilityOfAllCandidates = __max(0.0, probabilityOfAllCandidates - elementToScore->matchProbabilityForBestScore); // need the max due to floating point lossage.
                probabilityOfAllCandidates += matchProbability; // Don't combine this with the previous line, it introduces floating point unhappiness.
                elementToScore->matchProbabilityForBestScore = matchProbability;
                elementToScore->bestScore = score;
#ifdef  TIME_STRING_DISTANCE
                if (NULL != stats && !usedHamming) {
                    stats->nanosTimeInBSD[anyNearbyCandidatesAlreadyScored ? 1 : 0][score==-1 ? 0:1] += timeInBSD;
                    stats->BSDCounts[anyNearbyCandidatesAlreadyScored ? 1 : 0][score==-1 ? 0:1]++;
                }
#endif  // TIME_STRING_DISTANCE
    
                if (bestScore > score ||
                    bestScore == score && matchProbability > probabilityOfBestCandidate) {

                    //
                    // We have a new best score.  The old best score becomes the second best score, unless this is the same as the best or second best score
                    //                    

                    if ((secondBestScore == UnusedScoreValue || !(secondBestScoreGenomeLocation + maxK > genomeLocation && secondBestScoreGenomeLocation < genomeLocation + maxK)) &&
                        (bestScore == UnusedScoreValue || !(bestScoreGenomeLocation + maxK > genomeLocation && bestScoreGenomeLocation < genomeLocation + maxK)) &&
                        (!anyNearbyCandidatesAlreadyScored || (bestScoreGenomeLocation / (2 * maxMergeDist) != genomeLocation / (2 * maxMergeDist) &&
                                                               secondBestScoreGenomeLocation / (2 * maxMergeDist) != genomeLocation / (2 * maxMergeDist)))) {
                            secondBestScore = bestScore;
                            secondBestScoreGenomeLocation = bestScoreGenomeLocation;
                            secondBestScoreDirection = *hitDirection;
                    }
  
                    bestScore = score;
                    probabilityOfBestCandidate = matchProbability;
                    _ASSERT(probabilityOfBestCandidate <= probabilityOfAllCandidates);
                    bestScoreGenomeLocation = genomeLocation;
                    *singleHitGenomeLocation = bestScoreGenomeLocation;
                    *finalScore = bestScore;
                    *hitDirection = elementToScore->direction;
    
                    lvScoresAfterBestFound = 0;
    
                } else if (secondBestScore > score) {
                    //
                    // A new second best.
                    //
                    secondBestScore = score;
                    secondBestScoreGenomeLocation = genomeLocation;
                    secondBestScoreDirection = elementToScore->direction;
                }
 
                if (stopOnFirstHit && bestScore <= maxK) {
                    // The user just wanted to find reads that match the database within some distance, but doesn't
                    // care about the best alignment. Stop now but mark the result as MultipleHits because we're not
                    // confident that it's the best one.  We don't support mapq in this secnario, because we haven't
                    // explored enough to compute it.
                    *result = MultipleHits;
                    *mapq = 0;
                    return true;
                }

                // Update scoreLimit since we may have improved bestScore or secondBestScore
                if (maxHitsToGet > 0) {
                    scoreLimit = min(bestScore, maxK) + 3;
                } else {
                    scoreLimit = min(bestScore, maxK) + 2;
                    //scoreLimit = min(bestScore, maxK) + confDiff - 1 + (popularSeedsSkipped >= adaptiveConfDiffThreshold ? 1 : 0);
                    //if (bestScore != UnusedScoreValue && secondBestScore < bestScore + confDiff) {
                    //    // Since secondBestScore already means that our best location so far won't be a SingleHit,
                    //    // we really just care about finding better locations than that. However, still check for
                    //    // scores up to bestScore - 1 since those might become a new secondBestScore.
                    //    if (bestScore == 0) {
                    //        scoreLimit = 0;
                    //    } else {
                    //        scoreLimit = bestScore - 1;
                    //    }
                    //}
                    //// Make sure that the score limit is 
                    //// always search for something at least one worse than the best we found to drive the denominator of the MAPQ computation
                    //scoreLimit = __min(__max(scoreLimit,bestScore+2),maxK); 

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
        
        // Check whether we know we'll return MAPQ=0 for this read (NOTE: doesn't work right and doesn't give much speed)
        //if (bestScore != UnusedScoreValue) {
        //   double bestPossibleUnseenScore = __min(firstPassSeedsNotSkipped, firstPassRCSeedsNotSkipped);
        //   double bestPossibleUnseenMatchProb = 0.001 * pow(0.5, __max(0, bestPossibleUnseenScore-1));
        //   double bestPossibleProb = __max(bestPossibleUnseenMatchProb, probabilityOfBestCandidate);
        //   if (probabilityOfAllCandidates - probabilityOfBestCandidate > 0.9 * bestPossibleProb) {
        //       *result = MultipleHits;
        //       *singleHitGenomeLocation = bestScoreGenomeLocation;
        //       *finalScore = bestScore;
        //       *mapq = 0;
        //       return true;
        //   }
        //}
    } while (forceResult);

    return false;
}

    void 
BaseAligner::fillHitsFound(
    unsigned    maxHitsToGet, 
    int        *multiHitsFound, 
    unsigned   *multiHitLocations,
    Direction  *multiHitDirections,
    int        *multiHitScores)
/*++

Routine Description:

    Return up to maxHitsToGet best hits found for the user in the given array and count parameters, as long
    as maxHitsToGet >= 1. We also need the popularSeedsSkipped parameter to set confDiff adaptively.

--*/
{
    if (maxHitsToGet > 0) {
        unsigned realConfDiff = confDiff + (popularSeedsSkipped > adaptiveConfDiffThreshold ? 1 : 0);
        *multiHitsFound = 0;
        int firstDist = 0; // Distance of the best hit
        while (firstDist < MAX_K && hitCount[firstDist] == 0) {
            firstDist++;
        }
        for (int dist = firstDist; dist < min(firstDist + 4, MAX_K); dist++) {
            for (unsigned i = 0; i < hitCount[dist]; i++) {
                multiHitLocations[*multiHitsFound] = hitLocations[dist][i];
                multiHitDirections[*multiHitsFound] = hitDirections[dist][i];
                multiHitScores[*multiHitsFound] = dist;
                *multiHitsFound += 1;
                if (*multiHitsFound == maxHitsToGet) {
                    return;
                }
            }
        }
    }
}

    bool
BaseAligner::findElement(
    unsigned         genomeLocation, 
    Direction        direction,
    HashTableElement **hashTableElement)
{
    HashTableAnchor *hashTable = candidateHashTable[direction];

    unsigned lowOrderGenomeLocation = genomeLocation % (2 * maxMergeDist);
    unsigned highOrderGenomeLocation = genomeLocation - lowOrderGenomeLocation;

    unsigned hashTableIndex = hash(highOrderGenomeLocation) % candidateHashTablesSize;
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
    unsigned         genomeLocation, 
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
    unsigned lowOrderGenomeLocation = genomeLocation % (2 * maxMergeDist);

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
    unsigned    genomeLocation,
    Direction   direction,
    unsigned    lowestPossibleScore,
    int         seedOffset,
    Candidate **candidate,
    HashTableElement **hashTableElement)
/*++

Routine Description:

Arguments:

Return Value:

--*/
{
    HashTableAnchor *hashTable = candidateHashTable[direction];

    unsigned lowOrderGenomeLocation = genomeLocation % (2 * maxMergeDist);
    unsigned highOrderGenomeLocation = genomeLocation - lowOrderGenomeLocation;

    unsigned hashTableIndex = hash(highOrderGenomeLocation) % candidateHashTablesSize;
    HashTableAnchor *anchor = &hashTable[hashTableIndex];
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

    if (anchor->epoch == hashTableEpoch) {
        element->next = anchor->element;
    } else {
        anchor->epoch = hashTableEpoch;
        element->next = NULL;
    }
    anchor->element = element;

    element->candidatesUsed = 0;
    element->candidatesScored = 0;
    element->lowestPossibleScore = lowestPossibleScore;
    element->direction = direction;
    element->candidatesUsed = (_uint64)1 << lowOrderGenomeLocation;
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
}

#if     DBG
    bool
BaseAligner::AssertCandidatesAreSorted()
/*++

Routine Description:

Arguments:

Return Value:

--*/
{
    for (int i = 1; i < nCandidates; i++) {
        _ASSERT(candidates[i].genomeLocation <= candidates[i-1].genomeLocation);
    }
}
#endif  // DBG

BaseAligner::~BaseAligner()
/*++

Routine Description:

Arguments:

Return Value:

--*/
{
#if     MAINTAIN_HISTOGRAMS
    delete lvHistogram;
    delete lookupHistogram;
    delete lvHistogramForMulti;
    delete lvCountWhenBestFound;
#endif  // MAINTAIN_HISTOGRAMS

    delete probDistance;

#if     defined(USE_BOUNDED_STRING_DISTANCE)
    delete boundedStringDist;
#endif  // !bsd   

    if (hadBigAllocator) {
        //
        // Since these got allocated with the alloator rather than new, we want to call
        // their destructors without freeing their memory (which is the responsibility of
        // the owner of the alllocator).
        //
        if (ownLandauVishkin && NULL != landauVishkin) {
            landauVishkin->~LandauVishkin();
        }
        if (NULL != reverseLandauVishkin) {
            reverseLandauVishkin->~LandauVishkin();
        }
    } else {

        if (ownLandauVishkin && NULL != landauVishkin) {
            delete landauVishkin;
        }

        if (NULL != reverseLandauVishkin) {
            delete reverseLandauVishkin;
        }

        BigDealloc(candidates);
        candidates = NULL;

        BigDealloc(rcReadData);
        rcReadData = NULL;

        BigDealloc(reversedRead[FORWARD]);
        reversedRead[FORWARD] = NULL;
        reversedRead[RC] = NULL;

        BigDealloc(seedUsed);
        seedUsed = NULL;

        BigDealloc(candidateHashTable[FORWARD]);
        candidateHashTable[FORWARD] = NULL;
                
        BigDealloc(candidateHashTable[RC]);
        candidateHashTable[RC] = NULL;

        BigDealloc(weightLists);
        weightLists = NULL;

        BigDealloc(hashTableElementPool);
        hashTableElementPool = NULL;
    }
}

    void
BaseAligner::ComputeHitDistribution(
        Read        *read,
        unsigned     correctGenomeLocation,
        Direction    correctHitDirection,
        unsigned    *hitCountBySeed[NUM_DIRECTIONS],
        unsigned    *nSeedsApplied[NUM_DIRECTIONS],
        unsigned    *hitsContainingCorrectLocation)
{
    nSeedsApplied[FORWARD] = nSeedsApplied[RC] = 0;
 
    for (unsigned i = 0; i < maxSeedsToUse; i++) {
        hitCountBySeed[FORWARD][i] = hitCountBySeed[RC][i] = 0;
        hitsContainingCorrectLocation[i] = 0;
    }

    //
    // A bitvector for used seeds, indexed on the starting location of the seed within the read.
    //
    if (read->getDataLength() > maxReadSize) {
        fprintf(stderr,"BaseAligner:: got too big read (%d > %d)", read->getDataLength(),maxReadSize);
        exit(1);
    }

    if (read->getDataLength() < seedLen) {
        //
        // Too short to have any seeds, it's hopeless.
        //
        return;
    }

    //
    // Clear out the seed used array.
    //
    memset(seedUsed,0,(read->getDataLength() + 7) / 8);

    // Compute the reversed complement and count of Ns
    int readLen = (int)read->getDataLength();
    const char *readData = read->getData();
    unsigned countOfNs = 0;
    for (int i = 0; i < readLen; i++) {
        char baseByte = readData[i];
        char complement = rcTranslationTable[baseByte];
        rcReadData[readLen - i - 1] = complement;
        countOfNs += nTable[baseByte];
    }

    if (countOfNs > maxK) {
        return;
    }

    //
    // Block off any seeds that would contain an N.
    //
    //if (countOfNs > 0) {
    //    int minSeedToConsiderNing = 0; // In English, any word can be verbed. Including, apparently, "N."
    //    for (int i = 0; i < readLen; i++) {
    //        if (readData[i] == 'N') {
    //            int limit = __max(i + (int)seedLen - 1, readLen-1);
    //            for (int j = __max(minSeedToConsiderNing, i - (int)seedLen + 1); j <= limit; j++) {
    //                SetSeedUsed(j);
    //            }
    //            minSeedToConsiderNing = limit+1;
    //        }
    //    }
    //}

    Read rcRead[1];
    rcRead->init(NULL,0,rcReadData,NULL,read->getDataLength());

    unsigned nPossibleSeeds = read->getDataLength() - seedLen + 1;

    unsigned nextSeedToTest = 0;
    unsigned wrapCount = 0;

    while (*nSeedsApplied[FORWARD] + *nSeedsApplied[RC] < maxSeedsToUse) {
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
                return;
            }
            nextSeedToTest = getWrappedNextSeedToTest(wrapCount);
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

        Seed seed(read->getData() + nextSeedToTest, seedLen);

        const unsigned  *hits[NUM_DIRECTIONS];         // The actual hits (of size nHits)
    
        genomeIndex->lookupSeed(seed, &hitCountBySeed[FORWARD][*nSeedsApplied[FORWARD]], &hits[FORWARD], &hitCountBySeed[RC][*nSeedsApplied[RC]], &hits[RC]);

        SetSeedUsed(nextSeedToTest);

        for (Direction direction = 0; direction < NUM_DIRECTIONS; direction++) {
            if (hitCountBySeed[direction][*nSeedsApplied[direction]] <= maxHitsToConsider) {
                if (correctHitDirection == direction) {
                    //
                    // See if the correct answer is in this hit.
                    //
                    _int64 min = 0;
                    _int64 max = ((_int64)hitCountBySeed[direction][*nSeedsApplied[direction]]) - 1;

                    while (min <= max) {
                        _int64 probe = (min + max) / 2;
                        _ASSERT(probe >= min && probe <= max);

                        if (hits[direction][probe] == correctGenomeLocation + nextSeedToTest) {
                            hitsContainingCorrectLocation[*nSeedsApplied[direction]] = hitCountBySeed[direction][*nSeedsApplied[direction]];
                            break;
                        }

                        //
                        // Slice off half of the region.  Recall that the hits in the table are sorted backwards, with
                        // hits[0] > hits[1], etc.
                        //
                        if (hits[direction][probe] > correctGenomeLocation) {
                            min = probe+1;
                        } else {
                            max = probe-1;
                        }
                    }
                }
                nSeedsApplied++;
            }
        } // for each direction

        nextSeedToTest += seedLen;
    }

    return;
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
    for (unsigned i = 1; i <= maxSeedsToUse; i++) {
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
    if (element->weight >= maxSeedsToUse) {
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
BaseAligner::getBigAllocatorReservation(bool ownLandauVishkin, unsigned maxHitsToConsider, unsigned maxReadSize, unsigned seedLen, unsigned maxSeedsToUse)
{
    unsigned nCandidates = __min(maxHitsToConsider * (maxReadSize - seedLen + 1) * 2, 200000);  // *2 is for reverse complement
    size_t candidateHashTablesSize = (maxHitsToConsider * maxSeedsToUse * 3)/2;    // *1.5 for hash table slack
    size_t hashTableElementPoolSize = maxHitsToConsider * maxSeedsToUse * 2 ;   // *2 for RC

    return
        sizeof(BaseAligner)                                         + // our own member variables
        (ownLandauVishkin ? 
            LandauVishkin<>::getBigAllocatorReservation() +
            LandauVishkin<-1>::getBigAllocatorReservation() : 0)    + // our LandauVishkin objects
        sizeof(Candidate) * nCandidates                             + // candidates
        sizeof(char) * maxReadSize * 2                              + // rcReadData
        sizeof(char) * maxReadSize * 4 + 2 * MAX_K                  + // reversed read (both)
        sizeof(BYTE) * (maxReadSize + 7 + 128) / 8                  + // seed used
        sizeof(HashTableElement) * hashTableElementPoolSize         + // hash table element pool
        sizeof(HashTableAnchor) * candidateHashTablesSize * 2       + // candidate hash table (both)
        sizeof(HashTableElement) * (maxSeedsToUse + 1);               // weight lists
}