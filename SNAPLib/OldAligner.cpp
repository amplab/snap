/*++

Module Name:

    Aligner.cpp

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
#include "LandauVishkin.h"
#include "OldAligner.h"

_int64 nCandidatesCreated = 0;
_int64 nCandidatesEliminated = 0;
_int64 nScoredCandidatesEliminated = 0;
_int64 nScoredCandidatesSoftEliminated = 0;


OldAligner::OldAligner(
    GenomeIndex *i_genomeIndex, 
    unsigned     i_confDiff, 
    unsigned     i_maxHitsToConsider, 
    unsigned     i_maxK,
    unsigned     i_maxReadSize,
    unsigned     i_maxSeedsToUse,
    unsigned     i_lvCutoff) : 
        genomeIndex(i_genomeIndex), confDiff(i_confDiff), maxHitsToConsider(i_maxHitsToConsider), maxK(i_maxK), 
        maxReadSize(i_maxReadSize), maxSeedsToUse(i_maxSeedsToUse), lvCutoff(i_lvCutoff)
/*++

Routine Description:

    Constructor for the OldAligner class.  Aligners align reads against an indexed genome.

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

#if     USE_CAREFUL_BASE_STATUS
    //
    // Allocate the seed status and candidate arrays.  Logically, they are new for each call to AlignRead, but by doing
    // the allocations here we avoid a new/delete cycle for each read.
    //
    bases = new BaseStatus[maxReadSize];
    if (NULL == bases) {
        fprintf(stderr,"OldAligner: unable to allocate bases.\n");
        exit(1);
    }
#endif  // USE_CAREFUL_BASE_STATUS

    candidates = new Candidate[(maxHitsToConsider * (maxReadSize - seedLen + 1)) * 2];  // *2 is for reverse compliment
    if (NULL == candidates) {
        fprintf(stderr,"OldAligner: unable to allocate candidates.\n");
        exit(1);
    }

    rcReadData = new char[maxReadSize];
    if (NULL == rcReadData) {
        fprintf(stderr,"OldAligner: unable to allocate rc read data.");
        exit(1);
    }

    rcTranslationTable['A'] = 'T';
    rcTranslationTable['G'] = 'C';
    rcTranslationTable['C'] = 'G';
    rcTranslationTable['T'] = 'A';
    rcTranslationTable['N'] = 'N';

    for (unsigned i = 0; i < 256; i++) {
        nTable[i] = 0;
    }

    nTable['N'] = 1;

    seedUsed = new BYTE[(maxReadSize + 7 + 128) / 8];    // +128 to make sure it extends at both before and after at least an _int64
    if (NULL == seedUsed) {
        fprintf(stderr,"OldAligner: unable to allocate seedUsed array.\n");
        exit(1);
    }
    seedUsedAsAllocated = seedUsed; // Save the pointer for the delete.
    seedUsed += 8;  // This moves the pointer up an _int64, so we now have the appropriate before buffer.
}


    AlignmentResult
OldAligner::AlignRead(Read *read, unsigned *genomeLocation, bool *hitIsRC, int *finalScore)
/*++

Routine Description:

    Align a particular read.

Arguments:

    read                - the read to align
    genomeLocation      - if this aligned to a SingleHit, the 0-based offset in the genome that this hit.  The aligner treats the entire
                          genome as a single string, even though it's really a set of chrosomes.  It just makes the code simpler.
                          The caller can convert to chromosome+offset by looking up the piece boudaries using the Genome object.
    hitIsRC             - the aligner tries to align both the given read and its reverse compliment.  If we found a SIngleHit this
                          bool is set to indicate whether that hit was on the reverse compliment.
    finalScore          - if a single or confident hit, this is the score of the hit (i.e., the LV distance from the read)

Return Value:

    ConfidentHit, SingleHit, MultiHit or NotFound depending on how the alignment went.

--*/
{
    unsigned lookupsThisRun = 0;
    unsigned lvScores = 0;
    unsigned lvScoresAfterBestFound = 0;
    bool anyHitsSkippedBecauseOfTooHighPopularityThisRun = false;

    AlignmentResult finalResult;

    //
    // A bitvector for used seeds, indexed on the starting location of the seed within the read.
    //
    if (read->getDataLength() > maxReadSize) {
        fprintf(stderr,"OldAligner:: got too big read (%d > %d)", read->getDataLength(),maxReadSize);
        exit(1);
    }

    if (read->getDataLength() < seedLen) {
        //
        // Too short to have any seeds, it's hopeless.
        //
        return NotFound;
    }

    //
    // Clear out the seed used array.
    //
    memset(seedUsed,0,(read->getDataLength() + 7) / 8);

    int rdl = (int)read->getDataLength();
    const char *readData = read->getData();
    unsigned countOfNs = 0;
    for (int i = 0; i < rdl; i++) {
        char baseByte = readData[i];
        rcReadData[rdl - i - 1] = rcTranslationTable[baseByte];
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
        for (int i = 0; i < rdl; i++) {
            if (readData[i] == 'N') {
                int limit = __max(i + (int)seedLen - 1, rdl-1);
                for (int j = __max(minSeedToConsiderNing, i - (int)seedLen + 1); j <= limit; j++) {
                    SetSeedUsed(j);
                }
                minSeedToConsiderNing = limit+1;
            }
        }
    }

    Read rcRead[1];
    rcRead->init(NULL,0,rcReadData,NULL,read->getDataLength());

    clearCandidates();

    //
    // Initialize the bases table, which represents which bases we've checked.
    // We have readSize - seeds size + 1 possible seeds.
    //
    unsigned nPossibleSeeds = read->getDataLength() - seedLen + 1;
#if     USE_CAREFUL_BASE_STATUS
    memset(bases,0,sizeof(BaseStatus) * read->getDataLength());
#endif  // USE_CAREFUL_BASE_STATUS

    unsigned nextSeedToTest = 0;
    unsigned wrapCount = 0;
    unsigned lowestPossibleScoreOfAnyUnseenLocation = 0;
    unsigned lowestPossibleRCScoreOfAnyUnseenLocation = 0;
#if     USE_CAREFUL_BASE_STATUS
    unsigned mostSeedsContainingAnyParticularBase = 0;
    unsigned mostRCSeedsContainingAnyParticularBase = 0;
#else   // USE_CAREFUL_BASE_STATUS
    unsigned mostSeedsContainingAnyParticularBase = 1;  // Instead of tracking this for real, we're just conservative and use wrapCount+1.  It's faster.
    unsigned mostRCSeedsContainingAnyParticularBase = 1;// ditto
#endif  // USE_CAREFUL_BASE_STATUS
    unsigned bestScore = UnusedScoreValue;
    unsigned bestScoreGenomeLocation;
    unsigned secondBestScore = UnusedScoreValue;
    unsigned nSeedsApplied = 0;
    unsigned nRCSeedsApplied = 0;

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
                    lvScores,
                    lvScoresAfterBestFound,
                    anyHitsSkippedBecauseOfTooHighPopularityThisRun);
#if     MAINTAIN_HISTOGRAMS
                lvHistogram->addToCount(lvScores,lvScores);
                lookupHistogram->addToCount(lookupsThisRun,lookupsThisRun);
                if (MultipleHits == finalResult) {
                    lvHistogramForMulti->addToCount(lvScores,lvScores);
                } else if (SingleHit == finalResult) {
                    lvCountWhenBestFound->addToCount(lvScores-lvScoresAfterBestFound);
                }
#endif  // MAINTAIN_HISTOGRAMS

                return finalResult;
            }
            nextSeedToTest = getWrappedNextSeedToTest(wrapCount);

#if     !USE_CAREFUL_BASE_STATUS
            mostSeedsContainingAnyParticularBase = wrapCount + 1;
            mostRCSeedsContainingAnyParticularBase = wrapCount + 1;
#endif  // !USE_CAREFUL_BASE_STATUS
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
        unsigned        nRCHits;       // Number of hits for the seed's reverse compliment
        const unsigned  *rcHits;       // The actual reverse compliment hits
    
        genomeIndex->lookupSeed(seed, &nHits, &hits, &nRCHits, &rcHits);
        nHashTableLookups++;
        lookupsThisRun++;

        SetSeedUsed(nextSeedToTest);

        bool appliedEitherSeed = false;

        if (nHits > maxHitsToConsider) {
            //
            // This seed is matching too many places.  Just pretend we never looked and keep going.
            //
            nHitsIgnoredBecauseOfTooHighPopularity++;
            anyHitsSkippedBecauseOfTooHighPopularityThisRun = true;
        } else {
            //
            // Apply this seed.  First update the table of seeds per base.
            //
#if     USE_CAREFUL_BASE_STATUS
            for (unsigned i = nextSeedToTest; i < nextSeedToTest + seedLen; i++) {
                bases[i].nSeedsIncludingThisBase++;
                mostSeedsContainingAnyParticularBase = __max(mostSeedsContainingAnyParticularBase,bases[i].nSeedsIncludingThisBase);
            }
#endif  // USE_CAREFUL_BASE_STATUS
    
            //
            // And now update the candidates list with any hits from this seed.  If lowest possible score of any unseen location is
            // more than best_score + confDiff then we know that if this location is newly seen then its location won't ever be a
            // winner, and we can ignore it.
            //
            for (unsigned i = 0 ; i < nHits; i++) {
                //
                // Find the genome location where the beginning of the read would hit, given a match on this seed.
                //
                unsigned genomeLocationOfThisHit = hits[i] - nextSeedToTest;
    
                Candidate *candidate = findCandidate(genomeLocationOfThisHit,false);
                if (NULL != candidate) {
                    candidate->weight++;
                    if (candidate->minRange > genomeLocationOfThisHit) {
                        candidate->minRange = genomeLocationOfThisHit;
                    } else if (candidate->maxRange < genomeLocationOfThisHit) {
                        candidate->maxRange = genomeLocationOfThisHit;
                    }
                } else if (bestScore + confDiff > lowestPossibleScoreOfAnyUnseenLocation && 
                           secondBestScore > lowestPossibleScoreOfAnyUnseenLocation) {
                    candidate = allocateNewCandidate(genomeLocationOfThisHit, false);
                    candidate->weight = 1;
                    candidate->minPossibleScore = lowestPossibleScoreOfAnyUnseenLocation;
                    candidate->scored = false;
                }
            }
            nSeedsApplied++;
            appliedEitherSeed = true;
        }

        if (nRCHits > maxHitsToConsider) {
            //
            // This seed is matching too many places.  Just pretend we never looked and keep going.
            //
            nHitsIgnoredBecauseOfTooHighPopularity++;
            anyHitsSkippedBecauseOfTooHighPopularityThisRun = true;
        } else {
            //
            // Apply this seed.  First update the table of seeds per base.
            // The RC seed is at offset ReadSize - SeedSize - seed offset in the RC seed.
            //
            // To see why, imagine that you had a read that looked like 0123456 (where the digits
            // represented some particular bases, and digit' is the base's compliment). Then the
            // RC of that read is 6'5'4'3'2'1'.  So, when we look up the hits for the seed at
            // offset 0 in the forward read (i.e. 012 assuming a seed size of 3) then the index
            // will also return the results for the seed's reverse compliment, i.e., 3'2'1'.
            // This happens as the last seed in the RC read.
            //
            unsigned rcOffset = read->getDataLength() - seedLen - nextSeedToTest;
#if     USE_CAREFUL_BASE_STATUS
            for (unsigned i = rcOffset; i < rcOffset + seedLen; i++) {
                bases[i].nRCSeedsIncludingThisBase++;
                mostRCSeedsContainingAnyParticularBase = __max(mostRCSeedsContainingAnyParticularBase,bases[i].nRCSeedsIncludingThisBase);
            }
#endif  // USE_CAREFUL_BASE_STATUS
    
            //
            // And now update the candidates list with any hits from this seed.  If lowest possible score of any unseen RC location is
            // more than best_score + confDiff then we know that if this location is newly seen then its location won't ever be a
            // winner, and we can ignore it.
            //
            for (unsigned i = 0 ; i < nRCHits; i++) {
                //
                // Find the genome location where the beginning of the read would hit, given a match on this seed.
                //
                unsigned genomeLocationOfThisHit = rcHits[i] - rcOffset;
    
                Candidate *candidate = findCandidate(genomeLocationOfThisHit,true);
                if (NULL != candidate) {
                    candidate->weight++;
                    if (candidate->minRange > genomeLocationOfThisHit) {
                        candidate->minRange = genomeLocationOfThisHit;
                    } else if (candidate->maxRange < genomeLocationOfThisHit) {
                        candidate->maxRange = genomeLocationOfThisHit;
                    }
                } else if (bestScore + confDiff > lowestPossibleRCScoreOfAnyUnseenLocation &&
                           secondBestScore > lowestPossibleRCScoreOfAnyUnseenLocation) {
                    candidate = allocateNewCandidate(genomeLocationOfThisHit, true);
                    candidate->weight = 1;
                    candidate->minPossibleScore = lowestPossibleRCScoreOfAnyUnseenLocation;
                    candidate->scored = false;
                }
            }
            nRCSeedsApplied++;
            appliedEitherSeed = true;
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
                        lvScores,
                        lvScoresAfterBestFound,
                        anyHitsSkippedBecauseOfTooHighPopularityThisRun)) {
#if     MAINTAIN_HISTOGRAMS
                lvHistogram->addToCount(lvScores,lvScores);
                lookupHistogram->addToCount(lookupsThisRun,lookupsThisRun);
                if (MultipleHits == finalResult) {
                    lvHistogramForMulti->addToCount(lvScores,lvScores);
                } else if (SingleHit == finalResult) {
                    lvCountWhenBestFound->addToCount(lvScores-lvScoresAfterBestFound);
                }
#endif  // MAINTAIN_HISTOGRAMS

                return finalResult;
            }
        }
    }

    //
    // Do the best with what we've got.
    //
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
            lvScores,
            lvScoresAfterBestFound,
            anyHitsSkippedBecauseOfTooHighPopularityThisRun);
#if     MAINTAIN_HISTOGRAMS
        lvHistogram->addToCount(lvScores,lvScores);
        lookupHistogram->addToCount(lookupsThisRun,lookupsThisRun);
        if (MultipleHits == finalResult) {
            lvHistogramForMulti->addToCount(lvScores,lvScores);
        } else if (SingleHit == finalResult) {
            lvCountWhenBestFound->addToCount(lvScores-lvScoresAfterBestFound);
        }
#endif  // MAINTAIN_HISTOGRAMS

    return finalResult;
}

    bool
OldAligner::score(
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
    unsigned        &lvScores,
    unsigned        &lvScoresAfterBestFound,
    bool             anyHitsSkippedBecauseOfTooHighPopularityThisRun)
/*++

Routine Description:

    Make progress in scoring a possibly partial alignment.  This is a private method of the OldAligner class that's used
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
    rcRead                                  - the reverse compliment of read
    result                                  - returns the result if we reach one
    singleHitGenomeLocation                 - returns the location in the genome if we return a single hit
    hitIsRC                                 - if we return a single hit, indicates whether it's on the reverse compliment
    nSeedsApplied                           - how may seeds have we looked up and not ignored in this alignment
    nRCSeedsApplied                         - how may reverse compliment seeds have we looked up and not ignored in this alignment
    mostSeedsContainingAnyParticularBase    - what's the largest number of seeds we've used that contain any single base from the read?
    mostRCSeedsContainingAnyParticularBase  - same thing for the reverse compliment read
    lowestPossibleScoreOfAnyUnseenLocation  - in/out the best score that any genome location that we haven't hit could possinly have
    lowestPossibleRCScoreOfAnyUnseenLocation- same thing for the revrse compliment read
    candidates                              - in/out the array of candidates that have hit and possibly been scored
    bestScore                               - in/out the best score we've seen so far (recall that score is string distance, so best is least).
    bestScoreGenomeLocation                 - in/out where in the genome was the best score
    secondBestScore                         - in/out what's the second best score we've seen
    lvScores                                - in/out how many calls to LandauVishkin have we made (instrumentation)
    lvScoresAfterBestFound                  - in/out how many calls to LandauVishkin have we made after finding the best score (instrumentation)

Return Value:

    true iff we've reached a result.  When called with forceResult, we'll always return true.

--*/
{
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

#if     1
    //
    // If we haven't seen enough seeds to exclude unseen locations (i.e., places
    // in the genome where we didn't get a hit), just say no and wait for more
    // lookups.
    //

    if (lowestPossibleScoreOfAnyUnseenLocation < confDiff && lowestPossibleRCScoreOfAnyUnseenLocation < confDiff) {
        if (forceResult) {
            //
            // We can't exclude things that haven't hit, but just blow through and do the best we can.
            //
        } else {
            return false;
        }
    }
#endif  // 0

    do {
        //
        // Run through the candidates doing two things:
        //  1) Excluding those that have a score (or minPossibleScore) >= bestScore + confDiff,
        //  2) Finding the highest weight candidate that's not been scored (if there is one).
        //
        unsigned biggestWeight = 0;
        unsigned minPossibleScoreOfBiggestWeight = 1000;
        Candidate *candidateToScore;
        int i = 0;
        unsigned nUnscoredCandidates = 0;
        int indelPartner1 = nCandidates;
        int indelPartner2 = nCandidates;

        while (i < nCandidates) {   // Use a while loop because sometimes we increase i, sometimes we decrease nCandidates
            unsigned minScoreForCandidate;
            if (candidates[i].scored) {
                minScoreForCandidate = candidates[i].score;
                _ASSERT(minScoreForCandidate >= bestScore);
            } else {
                if (candidates[i].isRC) {
                    if (lowestPossibleRCScoreOfAnyUnseenLocation >= candidates[i].weight) {
                        candidates[i].minPossibleScore = 
                            __max(candidates[i].minPossibleScore,lowestPossibleRCScoreOfAnyUnseenLocation - candidates[i].weight);
                    } 
                } else {
                    if (lowestPossibleScoreOfAnyUnseenLocation >= candidates[i].weight) {
                        candidates[i].minPossibleScore = 
                            __max(candidates[i].minPossibleScore,lowestPossibleScoreOfAnyUnseenLocation - candidates[i].weight);
                    } 
                }

                minScoreForCandidate = candidates[i].minPossibleScore;
            }
    
            if (minScoreForCandidate > bestScore + confDiff || minScoreForCandidate > secondBestScore) {
                //
                // This can't possibly win, exclude it.
                //
nCandidatesEliminated++;
if (candidates[i].scored) {
    nScoredCandidatesEliminated++;
    if (minScoreForCandidate < bestScore + confDiff) nScoredCandidatesSoftEliminated++;
}

                removeCandidateAtIndex(i);
                if (indelPartner2 == i) {
                    indelPartner2 = nCandidates;
                }
                continue;
            }

            if (!candidates[i].scored) {
                nUnscoredCandidates++;
                if (candidates[i].weight > biggestWeight || 
                    candidates[i].weight == biggestWeight && candidates[i].minPossibleScore < minPossibleScoreOfBiggestWeight) {

                    biggestWeight = candidates[i].weight;
                    minPossibleScoreOfBiggestWeight = candidates[i].minPossibleScore;
                    candidateToScore = &candidates[i];   // This is OK, because we've already done whatever excluding we need to <= i

                    if (i > 0 && 
                        candidates[i-1].scored && 
                        candidates[i-1].genomeLocation + confDiff >= candidates[i].genomeLocation) {

                        indelPartner1 = i-1;
                    } else {
                        indelPartner1 = nCandidates;
                    }

                    if (i < nCandidates-1 &&
                        candidates[i+1].scored &&
                        candidates[i+1].genomeLocation - confDiff <= candidates[i].genomeLocation) {

                        indelPartner2 = i+1;
                    } else {
                        indelPartner2 = nCandidates;
                    }
                }
            }

            i++;
        }

        if (biggestWeight > 0) {
            //
            // We have an unscored candidate.  Score it.
            //
            unsigned score;
            if (candidateToScore->isRC) {
                score = candidateToScore->score = landauVishkin.computeEditDistance(
                    genome->getSubstring(candidateToScore->genomeLocation,rcRead->getDataLength()),
                    rcRead->getDataLength(),
                    rcRead->getData(),
                    rcRead->getDataLength(),
                    __min(maxK,bestScore)+confDiff-1);
            } else {
                score = candidateToScore->score = landauVishkin.computeEditDistance(
                    genome->getSubstring(candidateToScore->genomeLocation,read->getDataLength()),
                    read->getDataLength(),
                    read->getData(),
                    read->getDataLength(),
                    __min(maxK,bestScore)+confDiff-1);
            }
            candidateToScore->scored = true;
            nLocationsScored++;
            lvScores++;
            lvScoresAfterBestFound++;
    
            //
            // See if there's an indel partner for this candidate.  If so, delete the one
            // that's got the lower score.  We carefully do indelPartner2 before 1 because
            // if we delete 1 (or candidateToScore) then the index for indelPartner2
            // won't be right anymore.
            //
            bool deletedCandidateToScore = false;
            bool bestScoreWasIndelPartner = false;
            if (indelPartner2 < nCandidates) {
                _ASSERT(candidates[indelPartner2].scored);
                _ASSERT(candidates[indelPartner2].genomeLocation <= candidateToScore->genomeLocation + confDiff);
                _ASSERT(&candidates[indelPartner2-1] == candidateToScore);
                if (candidates[indelPartner2].score <= score) {
                    //
                    // The newly scored candidate is the loser.  Delete it.
                    //
                    _ASSERT(score >= bestScore);
                    removeCandidateAtIndex(indelPartner2-1);
                    deletedCandidateToScore = true;
                } else {
                    bestScoreWasIndelPartner = (bestScore == candidates[indelPartner2].score);
                    removeCandidateAtIndex(indelPartner2);
                }
                nIndelsMerged++;
            }

            if (indelPartner1 < nCandidates) {
                _ASSERT(candidates[indelPartner1].scored);
                _ASSERT(deletedCandidateToScore || candidates[indelPartner1].genomeLocation + confDiff >= candidateToScore->genomeLocation);
                _ASSERT(deletedCandidateToScore || &candidates[indelPartner1+1] == candidateToScore);

                if (candidates[indelPartner1].score <= score && !deletedCandidateToScore) {
                    _ASSERT(score >= bestScore);
                    removeCandidateAtIndex(indelPartner1+1); // i.e., candidateToScore
                    deletedCandidateToScore = true;
                    nIndelsMerged++;
                } else if (!deletedCandidateToScore) {
                    bestScoreWasIndelPartner = (bestScore == candidates[indelPartner1].score);
                    removeCandidateAtIndex(indelPartner1);
                    _ASSERT(indelPartner1 >= 0 && indelPartner1 < nCandidates);
                    candidateToScore = &candidates[indelPartner1];  // We just slid the candidate down into this position.
                    nIndelsMerged++;
                }
            }

            _ASSERT(nUnscoredCandidates > 0);
            nUnscoredCandidates--;

// off until we fix the insert/delete problem            _ASSERT(candidates[candidateToScore].score >= candidates[candidateToScore].minPossibleScore); // Else we messed up minPossibleScore (or LV or something)

            if (!deletedCandidateToScore) {
                if (bestScore > candidateToScore->score) {
                    //
                    // We have a new best score.
                    //
                    if (!bestScoreWasIndelPartner) {
                        secondBestScore = bestScore;
                    }
        
                    bestScore = candidateToScore->score;
                    bestScoreGenomeLocation = candidateToScore->genomeLocation;
                    *hitIsRC = candidateToScore->isRC;

                    lvScoresAfterBestFound = 0;

                } else if (secondBestScore > candidateToScore->score) {
                    //
                    // A new second best.
                    //
                    secondBestScore = candidateToScore->score;
                }

                if (secondBestScore < confDiff) {
                    //
                    // We've scored two different places that are good enough that this must be an ambiguous read.
                    //
                    *result = MultipleHits;
                    return true;
                }
            }
        }
        

        if (bestScore + confDiff <= __min(lowestPossibleScoreOfAnyUnseenLocation,lowestPossibleRCScoreOfAnyUnseenLocation) || forceResult) {
            if (0 == nUnscoredCandidates) {
                //
                // We've scored all live candidates and excluded all non-candidates.  We have our
                // answer.
                //
                if (nCandidates > 0 && bestScore + confDiff <= secondBestScore && bestScore <= maxK) {
                    if (!anyHitsSkippedBecauseOfTooHighPopularityThisRun && !forceResult) {
                        *result = CertainHit;
                    } else {
                        *result = SingleHit;
                    }
                    *singleHitGenomeLocation = bestScoreGenomeLocation;
                    if (NULL != finalScore) {
                        *finalScore = bestScore;
                    }
                    return true;
                } else if (nCandidates == 0 || bestScore > maxK) {
                    *result = NotFound;
                    return true;
                } else {
                    *result = MultipleHits;
                    return true;
                }
            }
            //
            // We don't have an unambiguous winner, but we do know that there's no need to
            // look elsewhere in the genome.  We can simply score all the remaining candidates
            // and go with that, or else we can keep looking for more hits on the theory that it
            // will eliminate candidates without scoring them.
            //
            if (nUnscoredCandidates < 10) {
                //
                // Let's go with what we've got.
                //
                forceResult = true;
            }
        }

        if (lvScores >= lvCutoff) {
            //
            // Don't score any more, just go with what we've got now.
            //
            if (UnusedScoreValue == bestScore) {
                *result = NotFound;
                return true;
            }

            if (bestScore + confDiff <= secondBestScore) {
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

    OldAligner::Candidate *
OldAligner::findCandidate(
    unsigned         genomeLocation, 
    bool             isRC)
/*++

Routine Description:

    Find an existing candidate with the given genome location within the existing set
    of candidates.

    This is a simple binary search.

Arguments:

    genomeLocation - the location of the candidate we'd like to look up
    isRC - Are we searching for the ordinary or reverse compliment candidate (they're treated separately)

Return Value:

    A pointer to the candidate if found, NULL otherwise.

--*/
{
    AssertCandidatesAreSorted();

    Candidate candidateToFind;
    candidateToFind.genomeLocation = genomeLocation;
    candidateToFind.isRC = isRC;

    int min = 0;
    int max = nCandidates -1;

    //
    // Binary search.
    //

    while (min <= max) {
        int probe = (min + max) / 2;
        _ASSERT(probe >= 0 && probe <= nCandidates);
        if (candidates[probe] == candidateToFind) {
            return &candidates[probe];
        } else if (candidates[probe] < candidateToFind) {
            min = probe+1;
        } else {
            max = probe-1;
        }
    }

    return NULL;
}

    OldAligner::Candidate *
OldAligner::allocateNewCandidate(unsigned genomeLocation, bool isRC)
/*++

Routine Description:

Arguments:

Return Value:

--*/
{
    _ASSERT(NULL == findCandidate(genomeLocation, isRC));   // Caller must ensure this

nCandidatesCreated ++;

    Candidate newCandidate;
    newCandidate.genomeLocation = genomeLocation;
    newCandidate.isRC = isRC;
    
    int locationToInsertBefore = -1;

    if (0 == nCandidates) {
        locationToInsertBefore = 0;
    } else {
        int min = 0;
        int max = nCandidates - 1;
    
    
        //
        // Binary search.
        //
        while (min <= max) {
            int probe = (min + max) / 2;
            _ASSERT(probe >= 0 && probe < nCandidates);
            _ASSERT(candidates[probe] != newCandidate);
    
            if (candidates[probe] < newCandidate) {
                if (probe+1 == nCandidates || candidates[probe+1] > newCandidate) {
                    locationToInsertBefore = probe+1;
                    break;
                }
                min = probe+1;
            } else {
                _ASSERT(candidates[probe] > newCandidate);
                if (probe == 0 || candidates[probe-1] < newCandidate) {
                    locationToInsertBefore = probe;
                    break;
                }
                max = probe-1;
            }
        }
    }
    _ASSERT(locationToInsertBefore >= 0);   // i.e., we should always find a spot
    _ASSERT(locationToInsertBefore <= nCandidates + 1);

    //
    // Slide everything down one location.
    //
    memmove(candidates + locationToInsertBefore + 1, candidates + locationToInsertBefore,(nCandidates - locationToInsertBefore) * sizeof(Candidate));

    (nCandidates)++;

    _ASSERT(locationToInsertBefore >= 0 && locationToInsertBefore < nCandidates);
    candidates[locationToInsertBefore].genomeLocation = genomeLocation;
    candidates[locationToInsertBefore].isRC = isRC;

    AssertCandidatesAreSorted();

    return &candidates[locationToInsertBefore];
}

    void
OldAligner::removeCandidateAtIndex(int index)
/*++

Routine Description:

Arguments:

Return Value:

--*/
{
    _ASSERT(index < nCandidates);
    memmove(candidates + index, candidates + index + 1, (nCandidates - index -1) * sizeof(Candidate));
    nCandidates--;
}

#if     DBG
    bool
OldAligner::AssertCandidatesAreSorted()
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

OldAligner::~OldAligner()
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

    delete [] seedUsedAsAllocated;
    
#if      USE_CAREFUL_BASE_STATUS
    delete [] bases;
#endif  // USE_CAREFUL_BASE_STATUS

    delete [] candidates;
    delete [] rcReadData;
}

    void
OldAligner::ComputeHitDistribution(
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
        fprintf(stderr,"OldAligner:: got too big read (%d > %d)", read->getDataLength(),maxReadSize);
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

    int rdl = (int)read->getDataLength();
    const char *readData = read->getData();
    unsigned countOfNs = 0;
    for (int i = 0; i < rdl; i++) {
        char baseByte = readData[i];
        rcReadData[rdl - i - 1] = rcTranslationTable[baseByte];
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
        for (int i = 0; i < rdl; i++) {
            if (readData[i] == 'N') {
                int limit = __max(i + (int)seedLen - 1, rdl-1);
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
        const unsigned  *rcHits;       // The actual reverse compliment hits
    
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
