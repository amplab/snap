/*++

Module Name:

    IntersectingAligner.cpp

Abstract:

    SNAP genome aligner

Authors:

    Bill Bolosky, September, 2011

Environment:

    User mode service.

    This class is NOT thread safe.  It's the caller's responsibility to ensure that
    at most one thread uses an instance at any time.

Revision History:

    Adapted from Matei Zaharia's Scala implementation.

--*/

#include "stdafx.h"
#include "Aligner.h"
#include "IntersectingAligner.h"
#include "LandauVishkin.h"

using std::sort;


IntersectingAligner::IntersectingAligner(
    GenomeIndex *i_genomeIndex, 
    unsigned     i_confDiff, 
    unsigned     i_maxHitsToConsider, 
    unsigned     i_maxK,
    unsigned     i_maxReadSize,
    unsigned     i_maxSeedsToUse,
    unsigned     i_lvCutoff) : 
        genomeIndex(i_genomeIndex), confDiff((int)i_confDiff), maxHitsToConsider(i_maxHitsToConsider), maxK((int)i_maxK), 
        maxReadSize((int)i_maxReadSize), maxSeedsToUse((int)i_maxSeedsToUse), lvCutoff(i_lvCutoff)
/*++

Routine Description:

    Constructor for the Aligner class.  Aligners align reads against an indexed genome.

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

    lookupResults = new LookupResult[maxSeedsToUse];  
    rcLookupResults = new LookupResult[maxSeedsToUse];
    intersectionResultsSpaceSize = maxHitsToConsider * maxSeedsToUse * 2;  // *2 is for RC
    intersectionResultSpace = new unsigned[intersectionResultsSpaceSize];
    intersectionResults = new IntersectionResult[(__max(maxK,maxSeedsToUse)+1) * 2];    // +1 is for 0 misses
    rcIntersectionResults = new IntersectionResult[maxSeedsToUse * 2];    // again, *2 is for RC
    intersectionOffsets = new int[maxSeedsToUse];
    alreadyIntersectedOffsets = new int[maxSeedsToUse];

    if (NULL == lookupResults || NULL == intersectionResultSpace || NULL == intersectionResults ||
        NULL == rcLookupResults || NULL == rcIntersectionResults || NULL == intersectionOffsets) {

        fprintf(stderr,"IntersectingAligner: allocation failed in constructor\n");
        exit(1);
    }

    //
    // Initialize the isRC fields of the lookup results, since they never change.
    //
    for (int i = 0; i < maxSeedsToUse; i++) {
        lookupResults[i].isRC = false;
        rcLookupResults[i].isRC = true;
    }

    genome = genomeIndex->getGenome();
    seedLen = genomeIndex->getSeedLength();

    rcReadData = new char[maxReadSize];
    if (NULL == rcReadData) {
        fprintf(stderr,"IntersectingAligner: unable to allocate rc read data.");
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
        fprintf(stderr,"IntersectingAligner: unable to allocate seedUsed array.\n");
        exit(1);
    }
    seedUsedAsAllocated = seedUsed; // Save the pointer for the delete.
    seedUsed += 8;  // This moves the pointer up an _int64, so we now have the appropriate before buffer.
}


    AlignmentResult
IntersectingAligner::AlignRead(Read *read, unsigned *genomeLocation, bool *hitIsRC, int *finalScore)
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

Return Value:

    SingleHit, MultiHit or NotFound depending on how the alignment went.

--*/
{
    unsigned lookupsThisRun = 0;
    unsigned lvScores = 0;
    unsigned lvScoresAfterBestFound = 0;
    bool     truncatedIntersection = false;

    //
    // A bitvector for used seeds, indexed on the starting location of the seed within the read.
    //
    if ((int)read->getDataLength() > maxReadSize) {
        fprintf(stderr,"IntersectingAligner:: got too big read (%d > %d)", read->getDataLength(),maxReadSize);
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
    int countOfNs = 0;
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

    nLookupResults = 0;
    nRCLookupResults = 0;
    usedEntriesInIntersectionResultSpace = 0;
    nIntersectionResults = 0;
    nRCIntersectionResults = 0;

    //
    // Initialize the bases table, which represents which bases we've checked.
    // We have readSize - seeds size + 1 possible seeds.
    //
    unsigned nPossibleSeeds = read->getDataLength() - seedLen + 1;

    unsigned nextSeedToTest = 0;
    unsigned wrapCount = 0;
    unsigned mostSeedsContainingAnyParticularBase = 1;  // Instead of tracking this for real, we're just conservative and use wrapCount+1.  It's faster.
    unsigned mostRCSeedsContainingAnyParticularBase = 1;// ditto
    int nSeedsApplied = 0;
    int nRCSeedsApplied = 0;

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
                break;
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

        if (wrapCount > 0) {
            //
            // Just finish, since all we're doing is looking at the set size.
            //
            break;
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

        if (nHits < maxHitsToConsider) {
            _ASSERT(nLookupResults < maxSeedsToUse);
            LookupResult *result = &lookupResults[nLookupResults];

            lookupResults[nLookupResults].nHits = nHits;
            lookupResults[nLookupResults].hits = hits;
            lookupResults[nLookupResults].offset = nextSeedToTest;
            appliedEitherSeed = true;

            nLookupResults++;
        } 

        if (nRCHits < maxHitsToConsider) {
            _ASSERT(nRCLookupResults < maxSeedsToUse);
            LookupResult *result = &rcLookupResults[nRCLookupResults];

            rcLookupResults[nRCLookupResults].nHits = nRCHits;
            rcLookupResults[nRCLookupResults].hits = rcHits;
            rcLookupResults[nRCLookupResults].offset = read->getDataLength() - seedLen - nextSeedToTest;
            appliedEitherSeed = true;

            nRCLookupResults++;
        }

        //
        // Move us along.
        //
        nextSeedToTest += seedLen;
    }

    unsigned nIntersections = 0;
    unsigned nRCIntersections = 0;
    if (nLookupResults > 0) {
        nIntersections = __min(nLookupResults-1,(int)maxIntersections);
        computeIntersection(NULL,0,0,nIntersections,lookupResults,nLookupResults,intersectionResults,&truncatedIntersection);
    }
    if (nRCLookupResults > 0) {
        nRCIntersections = __min(nRCLookupResults-1,(int)maxIntersections);
        computeIntersection(NULL,0,0,nRCIntersections,rcLookupResults,nRCLookupResults,rcIntersectionResults,&truncatedIntersection);
    }

    //
    // Now score what we've got to find the best match.
    //
    int bestScore = maxK+1;
    int secondBestScore = maxK+1;
    bool bestScoreIsRC;
    unsigned bestScoreGenomeLocation;
    for (unsigned i = 0; i < __max(nIntersections,nRCIntersections); i++) {
        //
        // See if we can eliminate the need to score the next set of matches.  We know that their
        // best possible score is nMisses, so if we've got a score that's less than nMisses + confDiff
        // we can quit now.
        //
        if (bestScore + confDiff <= intersectionResults[i].nMisses) {
            if (bestScore + confDiff < secondBestScore) {
                *genomeLocation = bestScoreGenomeLocation;
                *hitIsRC = bestScoreIsRC;
                if (NULL != finalScore) {
                    *finalScore = bestScore;
                }

                if (truncatedIntersection) {
                    return SingleHit;
                } else {
                    return CertainHit;
                }
            }
            return MultipleHits;
        }
        
        for (unsigned j = 0; j < intersectionResults[i].nHits; j++) {
            int score = landauVishkin.computeEditDistance(
                            genome->getSubstring(intersectionResults[i].hits[j],read->getDataLength()),
                            read->getDataLength(),
                            read->getData(),
                            read->getDataLength(),
                            __min(maxK,bestScore+confDiff));

            nLocationsScored++;
            lvScores++;
            lvScoresAfterBestFound++;            

            if (score > -1) {
                if (score < bestScore) {
                    secondBestScore = bestScore;
                    bestScore = score;
                    bestScoreIsRC = false;
                    bestScoreGenomeLocation = intersectionResults[i].hits[j];

                    if (secondBestScore <= confDiff) {
                        return MultipleHits;
                    }
                    lvScoresAfterBestFound = 0;
                } else if (score < secondBestScore) {
                    secondBestScore = score;
                    if (secondBestScore <= confDiff) {
                        return MultipleHits;
                    }
                }
            }
        } // For each hit

        for (unsigned j = 0; j < rcIntersectionResults[i].nHits; j++) {
            int score = landauVishkin.computeEditDistance(
                            genome->getSubstring(rcIntersectionResults[i].hits[j],rcRead->getDataLength()),
                            rcRead->getDataLength(),
                            rcRead->getData(),
                            rcRead->getDataLength(),
                            __min(maxK,bestScore+confDiff));

            nLocationsScored++;
            lvScores++;
            lvScoresAfterBestFound++;

            if (score > -1) {
                if (score < bestScore) {
                    secondBestScore = bestScore;
                    bestScore = score;
                    bestScoreIsRC = true;
                    bestScoreGenomeLocation = rcIntersectionResults[i].hits[j];

                    if (secondBestScore <= confDiff) {
                        return MultipleHits;
                    }
                    lvScoresAfterBestFound = 0;
                } else if (score < secondBestScore) {
                    secondBestScore = score;
                    if (secondBestScore <= confDiff) {
                        return MultipleHits;
                    }
                }
            }
        } // For each RC hit
    }

    if (bestScore >= maxK) {
        return NotFound;
    } else if (bestScore + confDiff <= secondBestScore) {
        *genomeLocation = bestScoreGenomeLocation;
        *hitIsRC = bestScoreIsRC;
        if (NULL != finalScore) {
            *finalScore = bestScore;
        }
        return SingleHit;
    }
    return MultipleHits;        
}

IntersectingAligner::~IntersectingAligner()
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
    
    delete [] rcReadData;
}

    bool
LookupResultComparitor(const IntersectingAligner::LookupResult& result1, const IntersectingAligner::LookupResult& result2)
{
    return result1.nHits < result2.nHits;
}


    void
IntersectingAligner::computeIntersection(
    IntersectionResult      *alreadyIntersected,
    int                      nAlreadyIntersected,
    int                      totalSetsRepresentedByAlreadyIntersected,
    int                      maxMisses,
    LookupResult            *lookups,
    int                      nLookups,
    IntersectionResult      *results,
    bool                    *truncatedIntersection
    )
{
    _ASSERT(nLookups >= 1);
    _ASSERT(maxMisses < nLookups);

    //
    // Start by sorting the input (lookups) array.
    //
    sort(lookups, lookups+nLookups, LookupResultComparitor);

    //
    // Initialize the results arrays.
    // 
    _ASSERT(maxMisses <= __max(maxK,maxSeedsToUse));

    for (int i = 0; i <= maxMisses; i++) {     // Note: <= maxMisses results in initializing maxMisses+1 results (0 .. maxMisses inclusive)
        results[i].nMisses = i;
        results[i].nHits = 0;
        results[i].isRC = lookups[0].isRC;
        results[i].hits = intersectionResultSpace + usedEntriesInIntersectionResultSpace;
        usedEntriesInIntersectionResultSpace += maxHitsToConsider; // This is overkill, but sufficient
        _ASSERT(usedEntriesInIntersectionResultSpace <= intersectionResultsSpaceSize);
    }

    //
    // Initialize the working offsets into the sets to be intersected.
    //
    memset(intersectionOffsets,0,sizeof(*intersectionOffsets) * nLookups);
    memset(alreadyIntersectedOffsets,0,sizeof(*alreadyIntersectedOffsets) * nAlreadyIntersected);

    //
    // Run the intersection.  At each step, find the largest of the entries among the
    // maxMisses - totalSetsRepresentedByAlreadyIntersected smallest sets and the already intersected sets, 
    // and then look for that in all of the sets.  
    // (Recall that the sets are all sorted from largest to smallest, so we find things backwards.)
    //
    for (;;) {
        unsigned valueToFind = 0;
        int setContainingValueToFind = maxMisses+1;
        int alreadyIntersectedContainingValueToFind = nAlreadyIntersected+1;
        //
        // Check the maxMisses 
        for (int i = 0; i <= maxMisses - totalSetsRepresentedByAlreadyIntersected; i++) {
            //
            // Assert that they're in (descending) order.
            //
            _ASSERT(intersectionOffsets[i]+1 >= (int)lookups[i].nHits || lookups[i].hits[intersectionOffsets[i]] > lookups[i].hits[intersectionOffsets[i]+1]);

            if (intersectionOffsets[i] < (int)lookups[i].nHits && lookups[i].hits[intersectionOffsets[i]] - lookups[i].offset >= valueToFind) {
                valueToFind = lookups[i].hits[intersectionOffsets[i]] - lookups[i].offset;
                setContainingValueToFind = i;
            }
        }

        for (int i = 0; i < nAlreadyIntersected; i++) {
            //
            // Assert that they're in (descending) order.
            //
            _ASSERT(alreadyIntersectedOffsets[i]+1 >= (int)alreadyIntersected[i].nHits || 
                    alreadyIntersected[i].hits[alreadyIntersectedOffsets[i]] > alreadyIntersected[i].hits[alreadyIntersectedOffsets[i]+1]);

            //
            // Assert that the alreadyIntersected sets are disjoint (or at least that valueToFind is not in both this one and
            // some other one).
            //
            _ASSERT(alreadyIntersectedOffsets[i] >= (int)alreadyIntersected[i].nHits || 
                    alreadyIntersectedContainingValueToFind == maxMisses+1 ||
                    valueToFind != alreadyIntersected[i].hits[alreadyIntersectedOffsets[i]]);

            if (alreadyIntersectedOffsets[i] < (int)alreadyIntersected[i].nHits && 
                alreadyIntersected[i].hits[alreadyIntersectedOffsets[i]] >= valueToFind) {

                if (valueToFind != alreadyIntersected[i].hits[alreadyIntersectedOffsets[i]]) {
                    setContainingValueToFind = maxMisses+1;
                }
                valueToFind = alreadyIntersected[i].hits[alreadyIntersectedOffsets[i]];
                alreadyIntersectedContainingValueToFind = i;
            }
        }

        if (setContainingValueToFind == maxMisses+1 && alreadyIntersectedContainingValueToFind == nAlreadyIntersected+1) {
            //
            // We've intersected them all, we're done.
            //
            return;
        }

        //
        // Count the misses.  Do this in two stages, first for the maxMisses smallest sets (in which if there's
        // a hit it must be the next element to check because we selected the greatest of these to be 
        // valueToFind).  Then, run binary search on the remainder.  We don't need to look through the
        // alreadyIntersected, because there's at most one of them containing valueToFind and we already
        // know which one it is.
        //
        int nMisses;
        if (alreadyIntersectedContainingValueToFind == nAlreadyIntersected+1) {
            //
            // It's not in alreadyIntersected, which means that it missed on all of the sets represented by alreadyIntersected.
            //
            nMisses = totalSetsRepresentedByAlreadyIntersected;
        } else {
            nMisses = alreadyIntersected[alreadyIntersectedContainingValueToFind].nMisses;
        }

        for (int i = 0; i <= maxMisses; i++) {
            if (intersectionOffsets[i] >= (int)lookups[i].nHits) {
                //
                // We've already checked this entire set.  It's a miss.
                //
                nMisses++;
                continue;
            }
            _ASSERT(lookups[i].hits[intersectionOffsets[i]] - lookups[i].offset <= valueToFind);    // Because we selected the largest one

            if (lookups[i].hits[intersectionOffsets[i]] - lookups[i].offset == valueToFind) {
                intersectionOffsets[i]++;
            } else {
                nMisses++;
                if (nMisses > maxMisses) {
                    break;
                }
            }
        }
        _ASSERT(nMisses <= maxMisses);   // It had to be in at least one set, that's why we selected it.

        for (int i = maxMisses+1; i < nLookups; i++) {
            LookupResult *result = &lookups[i];
            unsigned offsetValueToFind = valueToFind + result->offset;
            int min = intersectionOffsets[i];
            int max = result->nHits-1;
            int probe = min;    // Need to init this because of the case where we've already read the whole thing and fall through the while loop

            while (min <= max) {
                probe = (min + max) / 2;

                if (result->hits[probe] > offsetValueToFind) {
                    min = probe+1;
                } else if (result->hits[probe] < offsetValueToFind) {
                    max = probe-1;
                } else {
                    //
                    // We have a hit.
                    //
                    break;
                }
            }
            if (min > max) {
                nMisses++;
                if (nMisses > maxMisses) {
                    break;
                } 
            }
            intersectionOffsets[i] = probe+1;

        } // checking the larger sets

        if (nMisses <= maxMisses) {
            //
            // Add the value to the appropriate intersection result set.
            //
            IntersectionResult *result = &results[nMisses];
            if (result->nHits >= maxHitsToConsider) {
                //
                // Too many hits in this one result.  Truncate.
                //
                *truncatedIntersection = true;
            } else {
                result->hits[result->nHits] = valueToFind;
                _ASSERT(result->nHits == 0 || result->hits[result->nHits-1] > valueToFind);  // Inserting in order
                result->nHits++;
            }
        }
    }   // for ever
}
