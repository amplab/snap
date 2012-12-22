/*++

Module Name:

    mapq.h

Abstract:

    Support functions for mapping quality

Authors:

    Bill Bolosky, December, 2012

Environment:

    User mode service.

Revision History:

 
--*/

#pragma once

void initializeMapqTables();

double mapqToProbability(int mapq);

inline int computeMAPQ(
    double      probabilityOfAllCandidates, 
    double      probabilityOfBestCandidate, 
    int         score, 
    int         firstPassSeedsNotSkipped, 
    int         firstPassRCSeedsNotSkipped, 
    unsigned    smallestSkippedSeed, 
    unsigned    smallestSkippedRCSeed
    )
{
    _ASSERT(probabilityOfAllCandidates >= probabilityOfBestCandidate);
    _ASSERT(probabilityOfBestCandidate >= 0.0);

    //
    // Skipped seeds are ones that the aligner didn't look at because they had too many hits during the first
    // pass throygh the read (i.e., they're disjoint).  Any genome location that was ignored because of
    // maxHits could have at least score of this number (because it would have to have a difference in each
    // of them).  Assume that there are as many reads as the smallest of the sets at this edit distance
    // away from the read (i.e., assume the worst case).  Use a probability of .001 for migrating an edit
    // distance (this is, of course, just a guess since it really depends on the base call qualities, etc.,
    // but since we didn't look at the genome locations at all, this will have to do).
    //

    double probabilityOfSkippedLocations = 0.0;
    if (0xffffffff != smallestSkippedSeed) {
        probabilityOfSkippedLocations = pow(.001, firstPassSeedsNotSkipped) * smallestSkippedSeed;
    }
    if (0xffffffff != smallestSkippedRCSeed) {
        probabilityOfSkippedLocations += pow(.001, firstPassRCSeedsNotSkipped) * smallestSkippedRCSeed;
    }

    double correctnessProbability = probabilityOfBestCandidate / (probabilityOfAllCandidates + probabilityOfSkippedLocations);
    int baseMAPQ;
    if (correctnessProbability >= 1) {
        baseMAPQ =  70;
    } else {
        baseMAPQ = __min(70,(int)(-10 * log10(1 - correctnessProbability)));
    }

    //
    // If the score is higher than the number of non-skipped seeds in the first pass, then there's a chance
    // that we didn't find genome locations with a score at least as low as what we found because of just
    // not having any matching seeds.  Decrement mapq when this happens.
    //
    if (score > __max(firstPassSeedsNotSkipped, firstPassRCSeedsNotSkipped)) {
        baseMAPQ = __max(0,baseMAPQ - 10 * (score - __max(firstPassSeedsNotSkipped, firstPassRCSeedsNotSkipped)));
    }

    return baseMAPQ;

}
