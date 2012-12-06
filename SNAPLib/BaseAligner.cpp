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

using std::min;


BaseAligner::BaseAligner(
    GenomeIndex    *i_genomeIndex, 
    unsigned        i_confDiff, 
    unsigned        i_maxHitsToConsider, 
    unsigned        i_maxK,
    unsigned        i_maxReadSize,
    unsigned        i_maxSeedsToUse,
    unsigned        i_lvCutoff,
    unsigned        i_adaptiveConfDiffThreshold,
    LandauVishkin  *i_landauVishkin) : 
        genomeIndex(i_genomeIndex), confDiff(i_confDiff), maxHitsToConsider(i_maxHitsToConsider), maxK(i_maxK), 
        maxReadSize(i_maxReadSize), maxSeedsToUse(i_maxSeedsToUse), lvCutoff(i_lvCutoff), readId(-1),
        adaptiveConfDiffThreshold(i_adaptiveConfDiffThreshold), landauVishkin(i_landauVishkin),
	explorePopularSeeds(false), stopOnFirstHit(false)
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

--*/
{
#if     MAINTAIN_HISTOGRAMS
    lvHistogram = new Histogram(20, true);
    lookupHistogram = new Histogram(20, true);
    lvHistogramForMulti = new Histogram(20, true);
    lvCountWhenBestFound = new Histogram(maxHitsToConsider*4,false);
#endif  // MAINTAIN_HISTOGRAMS

    nHashTableLookups = 0;
    nLocationsScored = 0;
    nHitsIgnoredBecauseOfTooHighPopularity = 0;
    nReadsIgnoredBecauseOfTooManyNs = 0;
    nIndelsMerged = 0;

    genome = genomeIndex->getGenome();
    seedLen = genomeIndex->getSeedLength();

    bsd = new BoundedStringDistance<>(3);

    if (i_landauVishkin == NULL) {
        landauVishkin = new LandauVishkin;
        ownLandauVishkin = true;
    } else {
        ownLandauVishkin = false;
    }

    unsigned nCandidates = maxHitsToConsider * (maxReadSize - seedLen + 1) * 2;  // *2 is for reverse complement
    candidates = (Candidate *)BigAlloc(sizeof(Candidate) * nCandidates);
    for (unsigned i = 0 ; i < nCandidates; i++) {
        candidates[i].init();
    }

    rcReadData = (char *)BigAlloc(sizeof(char) * maxReadSize);

    rcTranslationTable['A'] = 'T';
    rcTranslationTable['G'] = 'C';
    rcTranslationTable['C'] = 'G';
    rcTranslationTable['T'] = 'A';
    rcTranslationTable['N'] = 'N';

    for (unsigned i = 0; i < 256; i++) {
        nTable[i] = 0;
    }

    nTable['N'] = 1;

    seedUsed = (BYTE *)BigAlloc(sizeof(BYTE) * (maxReadSize + 7 + 128) / 8);    // +128 to make sure it extends at both before and after at least an _int64

    seedUsedAsAllocated = seedUsed; // Save the pointer for the delete.
    seedUsed += 8;  // This moves the pointer up an _int64, so we now have the appropriate before buffer.

    nUsedHashTableElements = 0;
    hashTableElementPoolSize = maxHitsToConsider * maxSeedsToUse * 2 ;   // *2 for RC

    /*
    if (maxK*2 > 64) {
        fprintf(stderr,"The HashTableElement struct uses a 64 bit bitvector for used candidates.  As a result, you can't have maxK > 32.  Not that you'd want to anyway.\n");
        exit(1);
    }
    */

    candidateHashTablesSize = (maxHitsToConsider * maxSeedsToUse * 3)/2;    // *1.5 for hash table slack

    hashTableElementPool = (HashTableElement *)BigAlloc(sizeof(HashTableElement) * hashTableElementPoolSize);
    candidateHashTable[0] = (HashTableAnchor *)BigAlloc(sizeof(HashTableAnchor) * candidateHashTablesSize);
    candidateHashTable[1] = (HashTableAnchor *)BigAlloc(sizeof(HashTableAnchor) * candidateHashTablesSize);
    weightLists = (HashTableElement *)BigAlloc(sizeof(HashTableElement) * (maxSeedsToUse + 1));

    for (unsigned i = 0; i < hashTableElementPoolSize; i++) {
        hashTableElementPool[i].init();
    }

    for (unsigned i = 0; i < maxSeedsToUse + 1; i++) {
        weightLists[i].init();
    }

    memset(candidateHashTable[0],0,sizeof(HashTableAnchor) * candidateHashTablesSize);
    memset(candidateHashTable[1],0,sizeof(HashTableAnchor) * candidateHashTablesSize);
    hashTableEpoch = 0;
}

    AlignmentResult
BaseAligner::AlignRead(Read *read, unsigned *genomeLocation, bool *hitIsRC, int *finalScore)
{
    return AlignRead(read, genomeLocation, hitIsRC, finalScore, 0, 0, false, 0, NULL, NULL, NULL, NULL);
}

    AlignmentResult
BaseAligner::AlignRead(
    Read      *read,
    unsigned  *genomeLocation,
    bool      *hitIsRC,
    int       *finalScore,
    unsigned   searchRadius,
    unsigned   searchLocation,
    bool       searchRC)
{
    return AlignRead(read, genomeLocation, hitIsRC, finalScore,
        searchRadius, searchLocation, searchRC, 0, NULL, NULL, NULL, NULL);
}

    AlignmentResult
BaseAligner::AlignRead(
    Read      *read,
    unsigned  *genomeLocation,
    bool      *hitIsRC,
    int       *finalScore,
    unsigned   searchRadius,
    unsigned   searchLocation,
    bool       searchRC,
    int        maxHitsToGet,
    int       *multiHitsFound,
    unsigned  *multiHitLocations,
    bool      *multiHitRCs,
    int       *multiHitScores)  
/*++

Routine Description:

    Align a particular read, possibly constraining the search around a given location.

Arguments:

    read                - the read to align
    genomeLocation      - if this aligned to a SingleHit, the 0-based offset in the genome that this hit.  The aligner treats the entire
                          genome as a single string, even though it's really a set of chrosomes.  It just makes the code simpler.
                          The caller can convert to chromosome+offset by looking up the piece boudaries using the Genome object.
    hitIsRC             - the aligner tries to align both the given read and its reverse complement.  If we found a SingleHit this
                          bool is set to indicate whether that hit was on the reverse complement.
    finalScore          - if a single or confident hit, this is the score of the hit (i.e., the LV distance from the read)
    searchRadius        - if non-zero, constrain the search to this distance around searchLocation, in orientation searchRC
    searchLocation      - location to search around if searchRadius is given
    searchRC            - whether to search in reverse complement orientation if searchRadius is given
    maxHitsToGet        - if greater than 0, output up to this many hits within confDiff of the best (if any) in multiHitLocations,
                          writing their count in multiHitsFound, instead of returning MultipleHits immediately
    multiHitsFound      - output parameter for number of alternative hits found if maxHitsToGet is true
    multiHitLocations   - output parameter for locations of alternative hits found if maxHitsToGet is true

Return Value:

    ConfidentHit, SingleHit, MultiHit or NotFound depending on how the alignment went.

--*/
{
    *genomeLocation = 0xFFFFFFFF; // Value to return if we don't find a location.
    *hitIsRC = false;             // So we deterministically print the read forward in this case.
    if (finalScore != NULL) {
        *finalScore = UnusedScoreValue;
    }

    unsigned lookupsThisRun = 0;
    unsigned lvScores = 0;
    unsigned lvScoresAfterBestFound = 0;
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
    if (read->getDataLength() > maxReadSize) {
        fprintf(stderr,"BaseAligner:: got too big read (%d > %d)", read->getDataLength(),maxReadSize);
        exit(1);
    }

    if ((int)read->getDataLength() < seedLen) {
        //
        // Too short to have any seeds, it's hopeless.
        //
        return NotFound;
    }

    //
    // Clear out the seed used array.
    //
    memset(seedUsed, 0, (read->getDataLength() + 7) / 8);

    unsigned readLen = read->getDataLength();
    const char *readData = read->getData();
    unsigned countOfNs = 0;
    for (unsigned i = 0; i < readLen; i++) {
        char baseByte = readData[i];
        rcReadData[readLen - i - 1] = rcTranslationTable[baseByte];
        countOfNs += nTable[baseByte];
    }

    if (countOfNs > maxK) {
        nReadsIgnoredBecauseOfTooManyNs++;
        return NotFound;
    }

    //
    // Block off any seeds that would contain an N.
    //
    if (countOfNs > 0) {
        int minSeedToConsiderNing = 0; // In English, any word can be verbed. Including, apparently, "N."
        for (int i = 0; i < (int) readLen; i++) {
            if (readData[i] == 'N') {
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

    Read rcRead[1];
    rcRead->init(NULL, 0, rcReadData, NULL, readLen);

    clearCandidates();

    //
    // Initialize the bases table, which represents which bases we've checked.
    // We have readSize - seeds size + 1 possible seeds.
    //
    unsigned nPossibleSeeds = readLen - seedLen + 1;

    unsigned nextSeedToTest = 0;
    unsigned wrapCount = 0;
    unsigned lowestPossibleScoreOfAnyUnseenLocation = 0;
    unsigned lowestPossibleRCScoreOfAnyUnseenLocation = 0;
    unsigned mostSeedsContainingAnyParticularBase = 1;  // Instead of tracking this for real, we're just conservative and use wrapCount+1.  It's faster.
    unsigned mostRCSeedsContainingAnyParticularBase = 1;// ditto
    unsigned bestScore = UnusedScoreValue;
    unsigned bestScoreGenomeLocation;
    unsigned secondBestScore = UnusedScoreValue;
    unsigned secondBestScoreGenomeOffset;
    bool     secondBestScoreIsRC;
    unsigned nSeedsApplied = 0;
    unsigned nRCSeedsApplied = 0;
    unsigned scoreLimit;
    if (maxHitsToGet > 0) {
        scoreLimit = maxK + 3;
    } else {
       scoreLimit = maxK + confDiff - 1;
    }

    while (nSeedsApplied + nRCSeedsApplied < maxSeedsToUse) {
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
                    rcRead,
                    &finalResult,
                    finalScore,
                    genomeLocation,
                    hitIsRC,
                    nSeedsApplied,
                    nRCSeedsApplied,
                    mostSeedsContainingAnyParticularBase,
                    mostRCSeedsContainingAnyParticularBase,
                    lowestPossibleScoreOfAnyUnseenLocation,
                    lowestPossibleRCScoreOfAnyUnseenLocation,
                    candidates,
                    bestScore,
                    bestScoreGenomeLocation,
                    secondBestScore,
                    secondBestScoreGenomeOffset,
                    secondBestScoreIsRC,
                    scoreLimit,
                    lvScores,
                    lvScoresAfterBestFound,
                    maxHitsToGet);
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
                              multiHitLocations, multiHitRCs, multiHitScores);
                return finalResult;
            }
            nextSeedToTest = getWrappedNextSeedToTest(wrapCount);

            mostSeedsContainingAnyParticularBase = wrapCount + 1;
            mostRCSeedsContainingAnyParticularBase = wrapCount + 1;
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

        unsigned        nHits;         // Number of times this seed hits in the genome
        const unsigned  *hits;         // The actual hits (of size nHits)
        unsigned        nRCHits;       // Number of hits for the seed's reverse complement
        const unsigned  *rcHits;       // The actual reverse complement hits
    
        unsigned minSeedLoc = (minLocation < readLen ? 0 : minLocation - readLen);
        unsigned maxSeedLoc = (maxLocation > 0xFFFFFFFF - readLen ? 0xFFFFFFFF : maxLocation + readLen);
        genomeIndex->lookupSeed(seed, minSeedLoc, maxSeedLoc, &nHits, &hits, &nRCHits, &rcHits);
        nHashTableLookups++;
        lookupsThisRun++;

#ifdef TRACE_ALIGNER
        printf("Looked up seed %llx: hits=%u, rchits=%u\n", seed.getBases(), nHits, nRCHits);
        printf("Hits:");
        for (int i = 0; i < nHits; i++)
            printf(" %u", hits[i]);
        printf("\n");
        printf("RC hits:");
        for (int i = 0; i < nRCHits; i++)
            printf(" %u", rcHits[i]);
        printf("\n");
#endif

        SetSeedUsed(nextSeedToTest);

        bool appliedEitherSeed = false;

        //
        // Look for the seed in the forward direction, unless we were asked to search only the RC.
        //
        if (searchRadius == 0 || searchRC == false) {
            if (nHits > maxHitsToConsider && !explorePopularSeeds) {
                //
                // This seed is matching too many places.  Just pretend we never looked and keep going.
                //
                nHitsIgnoredBecauseOfTooHighPopularity++;
                popularSeedsSkipped++;
            } else {
                //
                // Update the candidates list with any hits from this seed.  If lowest possible score of any unseen location is
                // more than best_score + confDiff then we know that if this location is newly seen then its location won't ever be a
                // winner, and we can ignore it.
                //
                for (unsigned i = 0 ; i < min(nHits, maxHitsToConsider); i++) {
                    //
                    // Find the genome location where the beginning of the read would hit, given a match on this seed.
                    //
                    unsigned genomeLocationOfThisHit = hits[i] - nextSeedToTest;
                    if (genomeLocationOfThisHit < minLocation || genomeLocationOfThisHit > maxLocation)
                        continue;
    
                    Candidate *candidate;
                    HashTableElement *hashTableElement;

                    findCandidate(genomeLocationOfThisHit,false,&candidate,&hashTableElement);
                    if (NULL != hashTableElement) {
                        incrementWeight(hashTableElement);
                    } else if (lowestPossibleScoreOfAnyUnseenLocation <= scoreLimit) {
                        allocateNewCandidate(genomeLocationOfThisHit,false, lowestPossibleScoreOfAnyUnseenLocation, &candidate, &hashTableElement);
                    }
                }
                nSeedsApplied++;
                appliedEitherSeed = true;
            }
        }

        //
        // Look for the seed in the RC direction, unless we were asked to search only forwards.
        //
        if (searchRadius == 0 || searchRC == true) {
            if (nRCHits > maxHitsToConsider && !explorePopularSeeds) {
                //
                // This seed is matching too many places.  Just pretend we never looked and keep going.
                //
                nHitsIgnoredBecauseOfTooHighPopularity++;
                popularSeedsSkipped++;
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
                unsigned rcOffset = readLen - seedLen - nextSeedToTest;
    
                //
                // And now update the candidates list with any hits from this seed.  If lowest possible score of any unseen RC location is
                // more than best_score + confDiff then we know that if this location is newly seen then its location won't ever be a
                // winner, and we can ignore it.
                //
                for (unsigned i = 0 ; i < min(nRCHits, maxHitsToConsider); i++) {
                    //
                    // Find the genome location where the beginning of the read would hit, given a match on this seed.
                    //
                    unsigned genomeLocationOfThisHit = rcHits[i] - rcOffset;
                    if (genomeLocationOfThisHit < minLocation || genomeLocationOfThisHit > maxLocation)
                        continue;
    
                    Candidate *candidate;
                    HashTableElement *hashTableElement;

                    findCandidate(genomeLocationOfThisHit,true,&candidate,&hashTableElement);
                    if (NULL != hashTableElement) {
                        incrementWeight(hashTableElement);
                    } else if (lowestPossibleRCScoreOfAnyUnseenLocation <= scoreLimit) {
                        allocateNewCandidate(genomeLocationOfThisHit, true, lowestPossibleRCScoreOfAnyUnseenLocation, &candidate, &hashTableElement);
                    }
                }
                nRCSeedsApplied++;
                appliedEitherSeed = true;
            }
        }

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
                        rcRead,
                        &finalResult,
                        finalScore,
                        genomeLocation,
                        hitIsRC,
                        nSeedsApplied,
                        nRCSeedsApplied,
                        mostSeedsContainingAnyParticularBase,
                        mostRCSeedsContainingAnyParticularBase,
                        lowestPossibleScoreOfAnyUnseenLocation,
                        lowestPossibleRCScoreOfAnyUnseenLocation,
                        candidates,
                        bestScore,
                        bestScoreGenomeLocation,
                        secondBestScore,
                        secondBestScoreGenomeOffset,
                        secondBestScoreIsRC,
                        scoreLimit,
                        lvScores,
                        lvScoresAfterBestFound,
                        maxHitsToGet)) {
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
                              multiHitLocations, multiHitRCs, multiHitScores);
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
            rcRead,
            &finalResult,
            finalScore,
            genomeLocation,
            hitIsRC,
            nSeedsApplied,
            nRCSeedsApplied,
            mostSeedsContainingAnyParticularBase,
            mostRCSeedsContainingAnyParticularBase,
            lowestPossibleScoreOfAnyUnseenLocation,
            lowestPossibleRCScoreOfAnyUnseenLocation,
            candidates,
            bestScore,
            bestScoreGenomeLocation,
            secondBestScore,
            secondBestScoreGenomeOffset,
            secondBestScoreIsRC,
            scoreLimit,
            lvScores,
            lvScoresAfterBestFound,
            maxHitsToGet);
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
                  multiHitLocations, multiHitRCs, multiHitScores);
    return finalResult;
}

    bool
BaseAligner::score(
    bool             forceResult,
    Read            *read,
    Read            *rcRead,
    AlignmentResult *result,
    int             *finalScore,
    unsigned        *singleHitGenomeLocation,
    bool            *hitIsRC,
    unsigned         nSeedsApplied,
    unsigned         nRCSeedsApplied,
    unsigned         mostSeedsContainingAnyParticularBase,
    unsigned         mostRCSeedsContainingAnyParticularBase,
    unsigned        &lowestPossibleScoreOfAnyUnseenLocation,
    unsigned        &lowestPossibleRCScoreOfAnyUnseenLocation,
    Candidate       *candidates,
    unsigned        &bestScore,
    unsigned        &bestScoreGenomeLocation,
    unsigned        &secondBestScore,
    unsigned        &secondBestScoreGenomeLocation,
    bool            &secondBestScoreIsRC,
    unsigned        &scoreLimit,
    unsigned        &lvScores,
    unsigned        &lvScoresAfterBestFound,
    unsigned         maxHitsToGet)
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
    read                                    - the read we're aligning
    rcRead                                  - the reverse complement of read
    result                                  - returns the result if we reach one
    singleHitGenomeLocation                 - returns the location in the genome if we return a single hit
    hitIsRC                                 - if we return a single hit, indicates whether it's on the reverse complement
    nSeedsApplied                           - how may seeds have we looked up and not ignored in this alignment
    nRCSeedsApplied                         - how may reverse complement seeds have we looked up and not ignored in this alignment
    mostSeedsContainingAnyParticularBase    - what's the largest number of seeds we've used that contain any single base from the read?
    mostRCSeedsContainingAnyParticularBase  - same thing for the reverse complement read
    lowestPossibleScoreOfAnyUnseenLocation  - in/out the best score that any genome location that we haven't hit could possinly have
    lowestPossibleRCScoreOfAnyUnseenLocation- same thing for the reverse complement read
    candidates                              - in/out the array of candidates that have hit and possibly been scored
    bestScore                               - in/out the best score we've seen so far (recall that score is string distance, so best is least).
    bestScoreGenomeLocation                 - in/out where in the genome was the best score
    secondBestScore                         - in/out what's the second best score we've seen
    scoreLimit                              - in/out maximum edit distance we are need to consider for further scores (bigger ones won't help)
    lvScores                                - in/out how many calls to LandauVishkin have we made (instrumentation)
    lvScoresAfterBestFound                  - in/out how many calls to LandauVishkin have we made after finding the best score (instrumentation)
    popularSeedsSkipped                     - how many overly popular seeds we skipped during this alignment
    maxHitsToGet                            - whether the user wants to get multiple non-confident hits in the MultiHit case

Return Value:

    true iff we've reached a result.  When called with forceResult, we'll always return true.

--*/
{
#ifdef TRACE_ALIGNER
    printf("score() called with force=%d nsa=%d nrcsa=%d best=%u bestloc=%u 2nd=%u\n",
        forceResult, nSeedsApplied, nRCSeedsApplied, bestScore, bestScoreGenomeLocation, secondBestScore);
    printf("Read:\n  ");
    for (int i = 0; i < read->getDataLength(); i++) printf("%c", read->getData()[i]);
    printf("\n");
    //printf("Candidates:\n");
    //for (int i = 0; i < nCandidates; i++) {
    //    Candidate* c = candidates + i;
    //    printf("  loc=%u rc=%d weight=%u minps=%u scored=%d score=%u r=%u-%u\n",
    //        c->genomeLocation, c->isRC, c->weight, c->minPossibleScore, c->scored,
    //        c->score, c->minRange, c->maxRange);
    //}
    //printf("\n\n");
#endif

    if (0 == mostSeedsContainingAnyParticularBase && 0 == mostRCSeedsContainingAnyParticularBase) {
        //
        // The only way we can get here is if we've tried all of the seeds that we're willing
        // to try and every one of them generated too many hits to process.  Declare
        // a multi hit and give up.
        //
        _ASSERT(forceResult);
        *result = MultipleHits;
        return true;
    }

    //
    // Recompute lowestPossibleScore.
    //
    if (0 != mostSeedsContainingAnyParticularBase) {
        lowestPossibleScoreOfAnyUnseenLocation = 
            __max(lowestPossibleScoreOfAnyUnseenLocation, 
                  nSeedsApplied / mostSeedsContainingAnyParticularBase);
    }
    if (0 != mostRCSeedsContainingAnyParticularBase) {
        lowestPossibleRCScoreOfAnyUnseenLocation = 
            __max(lowestPossibleRCScoreOfAnyUnseenLocation, 
                  nRCSeedsApplied / mostRCSeedsContainingAnyParticularBase);
    }
#ifdef TRACE_ALIGNER
    printf("Lowest possible scores for unseen locations: %d (fwd), %d (RC)\n",
        lowestPossibleScoreOfAnyUnseenLocation,
        lowestPossibleRCScoreOfAnyUnseenLocation);
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

        if (__min(lowestPossibleScoreOfAnyUnseenLocation,lowestPossibleRCScoreOfAnyUnseenLocation) > scoreLimit || forceResult) {
            if (weightListToCheck == 0) {
                //
                // We've scored all live candidates and excluded all non-candidates.  We have our
                // answer.
                //
                int realConfDiff = confDiff + (popularSeedsSkipped >= adaptiveConfDiffThreshold ? 1 : 0);
                if (NULL != finalScore) {
                    *finalScore = bestScore;
                }
                if (bestScore + realConfDiff <= secondBestScore && bestScore <= maxK) {
                    if (popularSeedsSkipped == 0 && !forceResult) {
                        *result = CertainHit;
                    } else {
                        *result = SingleHit;
                    }
                    *singleHitGenomeLocation = bestScoreGenomeLocation;
                    return true;
                } else if (bestScore > maxK) {
                    // If none of our seeds was below the popularity threshold, report this as MultipleHits; otherwise,
                    // report it as NotFound
                    *result = (nSeedsApplied == 0 && nRCSeedsApplied == 0) ? MultipleHits : NotFound;
                    return true;
                } else {
                    *result = MultipleHits;
                    return true;
                }
            }
            //
            // Nothing that we haven't already looked up can possibly be the answer.  Score what we've got and exit.
            //
            forceResult = true;
        } else if (weightListToCheck == 0) {
            if (forceResult) {
                if (NULL != finalScore) {
                    *finalScore = bestScore;
                }
                int realConfDiff = confDiff + (popularSeedsSkipped >= adaptiveConfDiffThreshold ? 1 : 0);
                if (bestScore + realConfDiff <= secondBestScore && bestScore <= maxK) {
                    *result = SingleHit;
                    *singleHitGenomeLocation = bestScoreGenomeLocation;
                    return true;
                } else if (bestScore > maxK) {
                    // If none of our seeds was below the popularity threshold, report this as MultipleHits; otherwise,
                    // report it as NotFound
                    *result = (nSeedsApplied == 0 && nRCSeedsApplied == 0) ? MultipleHits : NotFound;
                    return true;
                } else {
                    *result = MultipleHits;
                    return true;
                }
            } else {
                //
                // No candidates, look for more.
                //
                return false;
            }
        }

        elementToScore = weightLists[weightListToCheck].weightNext;
        _ASSERT(!elementToScore->allExtantCandidatesScored);
        _ASSERT(elementToScore->candidatesUsed != 0);
        _ASSERT(elementToScore != &weightLists[weightListToCheck]);

        if (elementToScore->lowestPossibleScore <= scoreLimit) {
            unsigned long candidateIndexToScore;
            _uint64 candidatesMask = elementToScore->candidatesUsed;
            while (_BitScanForward64(&candidateIndexToScore,candidatesMask)) {
                candidatesMask &= ~(_uint64)1 << candidateIndexToScore;
                Candidate *candidateToScore = &elementToScore->candidates[candidateIndexToScore];
                if (candidateToScore->scoredInEpoch == hashTableEpoch) {
                    //
                    // This can happen if we find a candidate, score it, and then hit on an indel
                    // near it.
                    //
                    continue;
                }
    
                unsigned genomeLocation = elementToScore->baseGenomeLocation + candidateIndexToScore;
    
                unsigned score = -1;
                if (elementToScore->isRC) {
                    const char *data = genome->getSubstring(genomeLocation, rcRead->getDataLength() + MAX_K);
                    if (data != NULL) {
#ifdef USE_NEW_DISTANCE
                        score = bsd->compute(data, rcRead->getData(), rcRead->getDataLength(), scoreLimit);
#else
                        _uint64 cacheKey = 0;
                        if (readId != -1) {
                            cacheKey = ((_uint64) readId) << 33 | ((_uint64) elementToScore->isRC) << 32 | genomeLocation;
                        }
                        score = landauVishkin->computeEditDistance(
                            data, rcRead->getDataLength() + MAX_K,
                            rcRead->getData(), rcRead->getDataLength(),
                            scoreLimit, cacheKey);
#endif
                    }
#ifdef TRACE_ALIGNER
                    printf("Computing distance at %u (RC) with limit %d: %d\n", genomeLocation, scoreLimit, score);
#endif
                } else {
                    const char *data = genome->getSubstring(genomeLocation, read->getDataLength() + MAX_K);
                    if (data != NULL) {
#ifdef USE_NEW_DISTANCE
                        score = bsd->compute(data, read->getData(), read->getDataLength(), scoreLimit);
#else
                        _uint64 cacheKey = 0;
                        if (readId != -1) {
                            cacheKey = ((_uint64) readId) << 33 | ((_uint64) elementToScore->isRC) << 32 | genomeLocation;
                        }
                        score = landauVishkin->computeEditDistance(
                            data, read->getDataLength() + MAX_K,
                            read->getData(), read->getDataLength(),
                            scoreLimit, cacheKey);
#endif
                    }
#ifdef TRACE_ALIGNER
                    printf("Computing distance at %u (fwd) with limit %d: %d\n", genomeLocation, scoreLimit, score);
#endif
                }
                
                if (maxHitsToGet > 0 && score != -1 && hitCount[score] < maxHitsToGet) {
                    // Remember the location of this hit because we don't have enough at this distance
                    hitLocations[score][hitCount[score]] = genomeLocation;
                    hitRCs[score][hitCount[score]] = elementToScore->isRC;
                    hitCount[score]++;
                }
                
                candidateToScore->score = score;
                candidateToScore->scoredInEpoch = hashTableEpoch;
                elementToScore->bestScore = __min(elementToScore->bestScore, score);
    
                nLocationsScored++;
                lvScores++;
                lvScoresAfterBestFound++;
        
// off until we fix the insert/delete problem            _ASSERT(candidates[candidateToScore].score >= candidates[candidateToScore].minPossibleScore); // Else we messed up minPossibleScore (or LV or something)
    
                if (bestScore > score) {
		            if ((secondBestScore == UnusedScoreValue || !(secondBestScoreGenomeLocation + maxK > genomeLocation && secondBestScoreGenomeLocation < genomeLocation + maxK)) &&
		                (bestScore == UnusedScoreValue || !(bestScoreGenomeLocation + maxK > genomeLocation && bestScoreGenomeLocation < genomeLocation + maxK))) 
                    {
                            secondBestScore = bestScore;
                            secondBestScoreGenomeLocation = bestScoreGenomeLocation;
                            secondBestScoreIsRC = *hitIsRC;
                    }
    
                    //
                    // We have a new best score.
                    //
                    bestScore = score;
                    bestScoreGenomeLocation = genomeLocation;
                    *singleHitGenomeLocation = bestScoreGenomeLocation;
                    if (NULL != finalScore) {
                        *finalScore = bestScore;
                    }
                    *hitIsRC = elementToScore->isRC;
    
                    lvScoresAfterBestFound = 0;
    
                } else if (secondBestScore > score && 
                          (secondBestScore == UnusedScoreValue || !(secondBestScoreGenomeLocation + maxK > genomeLocation && secondBestScoreGenomeLocation < genomeLocation + maxK)) &&
                          !(bestScoreGenomeLocation + maxK > genomeLocation && bestScoreGenomeLocation < genomeLocation + maxK)) {
                    //
                    // A new second best.
                    //
                    secondBestScore = score;
                    secondBestScoreGenomeLocation = genomeLocation;
                    secondBestScoreIsRC = elementToScore->isRC;
                }

                if (stopOnFirstHit && bestScore <= maxK) {
                    // The user just wanted to find reads that match the database within some distance, but doesn't
                    // care about the best alignment. Stop now but mark the result as MultipleHits because we're not
                    // confident that it's the best one.
                    *result = MultipleHits;
                    return true;
                }

                // Update scoreLimit since we may have improved bestScore or secondBestScore
                if (maxHitsToGet > 0) {
                    scoreLimit = min(bestScore, maxK) + 3;
                } else {
                    scoreLimit = min(bestScore, maxK) + confDiff - 1 + (popularSeedsSkipped >= adaptiveConfDiffThreshold ? 1 : 0);
                    if (bestScore != UnusedScoreValue && secondBestScore < bestScore + confDiff) {
                        // Since secondBestScore already means that our best location so far won't be a SingleHit,
                        // we really just care about finding better locations than that. However, still check for
                        // scores up to bestScore - 1 since those might become a new secondBestScore.
                        // _ASSERT(bestScore >= 1); // Otherwise we would've returned MultipleHits
                        scoreLimit = bestScore - 1;
                    }
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
        
        if (lvScores >= lvCutoff) {
            //
            // Don't score any more, just go with what we've got now.
            //
            if (bestScore > maxK) {
                // If none of our seeds was below the popularity threshold, report this as MultipleHits; otherwise,
                // report it as NotFound
                *result = (nSeedsApplied == 0 && nRCSeedsApplied == 0) ? MultipleHits : NotFound;
                return true;
            }

            int realConfDiff = confDiff + (popularSeedsSkipped >= adaptiveConfDiffThreshold ? 1 : 0);
            if (bestScore + realConfDiff <= secondBestScore && bestScore <= maxK) {
                *result = SingleHit;
                *singleHitGenomeLocation = bestScoreGenomeLocation;
                if (NULL != finalScore) {
                    *finalScore = bestScore;
                }
                return true;
            }

            *result = MultipleHits;
            return true;
        }
    } while (forceResult);

    return false;
}

    void 
BaseAligner::fillHitsFound(
    unsigned    maxHitsToGet, 
    int        *multiHitsFound, 
    unsigned   *multiHitLocations,
    bool       *multiHitRCs,
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
                multiHitRCs[*multiHitsFound] = hitRCs[dist][i];
                multiHitScores[*multiHitsFound] = dist;
                *multiHitsFound += 1;
                if (*multiHitsFound == maxHitsToGet) {
                    return;
                }
            }
        }
    }
}

    void
BaseAligner::findCandidate(
    unsigned         genomeLocation, 
    bool             isRC,
    Candidate        **candidate,
    HashTableElement **hashTableElement)  
/*++

Routine Description:

    Find a candidate in the hash table.  If it doesn't exist, create, initialize and insert it.

Arguments:

    genomeLocation - the location of the candidate we'd like to look up
    candidate - The candidate that was found or created
    hashTableElement - the hashTableElement for the candidate that was found or created.

--*/
{
    HashTableAnchor *hashTable = candidateHashTable[isRC];

    unsigned lowOrderGenomeLocation = genomeLocation % (2 * maxMergeDist);
    unsigned highOrderGenomeLocation = genomeLocation - lowOrderGenomeLocation;

    unsigned hashTableIndex = hash(highOrderGenomeLocation) % candidateHashTablesSize;
    HashTableAnchor *anchor = &hashTable[hashTableIndex];
    if (anchor->epoch != hashTableEpoch) {
        //
        // It's empty.
        //
        *hashTableElement = NULL;
        *candidate = NULL;
        return;
    }

    HashTableElement *lookedUpElement = anchor->element;
    while (NULL != lookedUpElement && lookedUpElement->baseGenomeLocation != highOrderGenomeLocation) {
        lookedUpElement = lookedUpElement->next;
    }

    if (NULL != lookedUpElement) {
        *hashTableElement = lookedUpElement;
        *candidate = &lookedUpElement->candidates[lowOrderGenomeLocation];

        _uint64 bitForThisCandidate = (_uint64)1 << lowOrderGenomeLocation;
        lookedUpElement->allExtantCandidatesScored = lookedUpElement->allExtantCandidatesScored && (lookedUpElement->candidatesUsed & bitForThisCandidate);
        lookedUpElement->candidatesUsed |= bitForThisCandidate;
    } else {
        *hashTableElement = NULL;
        *candidate = NULL;
    }
}

bool doAlignerPrefetch = true;

    void
BaseAligner::allocateNewCandidate(unsigned genomeLocation, bool isRC, unsigned lowestPossibleScore, Candidate **candidate, HashTableElement **hashTableElement)
/*++

Routine Description:

Arguments:

Return Value:

--*/
{
    HashTableAnchor *hashTable = candidateHashTable[isRC];

    unsigned lowOrderGenomeLocation = genomeLocation % (2 * maxMergeDist);
    unsigned highOrderGenomeLocation = genomeLocation - lowOrderGenomeLocation;

    unsigned hashTableIndex = hash(highOrderGenomeLocation) % candidateHashTablesSize;
    HashTableAnchor *anchor = &hashTable[hashTableIndex];
    HashTableElement *element;

    if (doAlignerPrefetch) {
        genomeIndex->prefetchGenomeData(genomeLocation);
    }

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
    element->lowestPossibleScore = lowestPossibleScore;
    element->isRC = isRC;
    element->candidatesUsed = (_uint64)1 << lowOrderGenomeLocation;
    element->weight = 1;
    element->baseGenomeLocation = highOrderGenomeLocation;
    element->bestScore = UnusedScoreValue;
    element->allExtantCandidatesScored = false;

    //
    // And insert it at the end of weight list 1.
    //
    element->weightNext = &weightLists[1];
    element->weightPrev = weightLists[1].weightPrev;
    element->weightNext->weightPrev = element;
    element->weightPrev->weightNext = element;

    *candidate = &element->candidates[lowOrderGenomeLocation];
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

    BigDealloc(seedUsedAsAllocated);
    
    BigDealloc(candidates);
    BigDealloc(rcReadData);

    BigDealloc(hashTableElementPool);
    BigDealloc(candidateHashTable[0]);
    BigDealloc(candidateHashTable[1]);
    BigDealloc(weightLists);

    delete bsd;

    if (ownLandauVishkin) {
        delete landauVishkin;
    }
}

    void
BaseAligner::ComputeHitDistribution(
        Read        *read,
        unsigned     correctGenomeLocation,
        bool         correctHitIsRC,
        unsigned    *hitCountBySeed,
        unsigned    *rcHitCountBySeed,
        unsigned    &nSeedsApplied,
        unsigned    &nRCSeedsApplied,
        unsigned    *hitsContainingCorrectLocation)
{
    nSeedsApplied = 0;
    nRCSeedsApplied = 0;

    for (unsigned i = 0; i < maxSeedsToUse; i++) {
        hitCountBySeed[i] = 0;
        rcHitCountBySeed[i] = 0;
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

    int readLen = (int)read->getDataLength();
    const char *readData = read->getData();
    unsigned countOfNs = 0;
    for (int i = 0; i < readLen; i++) {
        char baseByte = readData[i];
        rcReadData[readLen - i - 1] = rcTranslationTable[baseByte];
        countOfNs += nTable[baseByte];
    }

    if (countOfNs > maxK) {
        return;
    }

    //
    // Block off any seeds that would contain an N.
    //
    if (countOfNs > 0) {
        int minSeedToConsiderNing = 0; // In English, any word can be verbed. Including, apparently, "N."
        for (int i = 0; i < readLen; i++) {
            if (readData[i] == 'N') {
                int limit = __max(i + (int)seedLen - 1, readLen-1);
                for (int j = __max(minSeedToConsiderNing, i - (int)seedLen + 1); j <= limit; j++) {
                    SetSeedUsed(j);
                }
                minSeedToConsiderNing = limit+1;
            }
        }
    }

    Read rcRead[1];
    rcRead->init(NULL,0,rcReadData,NULL,read->getDataLength());

    unsigned nPossibleSeeds = read->getDataLength() - seedLen + 1;

    unsigned nextSeedToTest = 0;
    unsigned wrapCount = 0;

    while (nSeedsApplied + nRCSeedsApplied < maxSeedsToUse) {
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

        const unsigned  *hits;         // The actual hits (of size nHits)
        const unsigned  *rcHits;       // The actual reverse complement hits
    
        genomeIndex->lookupSeed(seed, &hitCountBySeed[nSeedsApplied], &hits, &rcHitCountBySeed[nRCSeedsApplied], &rcHits);

        SetSeedUsed(nextSeedToTest);

        if (hitCountBySeed[nSeedsApplied] <= maxHitsToConsider) {

            if (!correctHitIsRC) {
                //
                // See if the correct answer is in this hit.
                //
                _int64 min = 0;
                _int64 max = ((_int64)hitCountBySeed[nSeedsApplied]) - 1;

                while (min <= max) {
                    _int64 probe = (min + max) / 2;
                    _ASSERT(probe >= min && probe <= max);

                    if (hits[probe] == correctGenomeLocation + nextSeedToTest) {
                        hitsContainingCorrectLocation[nSeedsApplied] = hitCountBySeed[nSeedsApplied];
                        break;
                    }

                    //
                    // Slice off half of the region.  Recall that the hits in the table are sorted backwards, with
                    // hits[0] > hits[1], etc.
                    //
                    if (hits[probe] > correctGenomeLocation) {
                        min = probe+1;
                    } else {
                        max = probe-1;
                    }
                }
            }
            nSeedsApplied++;
        }

        if (rcHitCountBySeed[nRCSeedsApplied] <= maxHitsToConsider) {
            if (correctHitIsRC) {
                //
                // See if the correct answer is in this hit.
                //
                _int64 min = 0;
                _int64 max = ((_int64)rcHitCountBySeed[nRCSeedsApplied]) - 1;

                while (min <= max) {
                    _int64 probe = (min + max) / 2;
                    _ASSERT(probe >= min && probe <= max);

                    if (rcHits[probe] - (read->getDataLength() - seedLen - nextSeedToTest) == correctGenomeLocation) {
                        hitsContainingCorrectLocation[nRCSeedsApplied] = rcHitCountBySeed[nRCSeedsApplied];
                        break;
                    }

                    //
                    // Slice off half of the region.  Recall that the hits in the table are sorted backwards, with
                    // hits[0] > hits[1], etc.
                    //
                    if (rcHits[probe] > correctGenomeLocation) {
                        min = probe+1;
                    } else {
                        max = probe-1;
                    }
                }
            }
            nRCSeedsApplied++;
        }

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
    isRC = false;
    allExtantCandidatesScored = false;
}

    void
BaseAligner::Candidate::init()
{
    score = UnusedScoreValue;
    scoredInEpoch = 0;
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
