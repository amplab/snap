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

static const double GenericBaseDifferenceProbability = .001;	// The assumed probability of a base differing from the reference, given we don't know it's phred score.

void initializeMapqTables();

double mapqToProbability(int mapq);
const int maxSeedsForUnseenProbability = 10;	// The chance is effectively 0 by here.
extern double firstOrderUnseenLocationProbability[];	// Chance that there are exactly nSeeds differences, one in each seed

double factorial(int n);		// Shouldn't this be in the C library?


inline int computeMAPQ(
    double      probabilityOfAllCandidates, 
    double      probabilityOfBestCandidate, 
    int         score, 
    int         bestPassSeedsNotSkipped, 
    int         bestPassRCSeedsNotSkipped, 
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
        probabilityOfSkippedLocations = pow(.001, bestPassSeedsNotSkipped) * smallestSkippedSeed;
    }
    if (0xffffffff != smallestSkippedRCSeed) {
        probabilityOfSkippedLocations += pow(.001, bestPassRCSeedsNotSkipped) * smallestSkippedRCSeed;
    }

	//
	// Because of the basic SNAP algorithm, there's always a chance that we never saw the correct location
	// because there's a difference in each seed.  In order for this to happen, there would have to be
	// exactly nSeeds differences spread perfectly, which has chance 
	// (difference probability)^nSeeds / ((nSeeds!) / nSeeds^nSeeds).  The term in the denominator comes
	// from the fact that the first difference can be anywhere, the second one can be in (nSeeds-1)/nSeeds places,
	// the second in (nSeeds - 2) / nSeeds, etc.  We've precomputed this.
	// because 
	double probabilityOfMissedLocations = (firstOrderUnseenLocationProbability[bestPassSeedsNotSkipped] + firstOrderUnseenLocationProbability[bestPassRCSeedsNotSkipped]) / 2;

    double correctnessProbability = probabilityOfBestCandidate / (probabilityOfAllCandidates + probabilityOfSkippedLocations + probabilityOfMissedLocations);
    int baseMAPQ;
    if (correctnessProbability >= 1) {
        baseMAPQ =  70;
    } else {
        baseMAPQ = __min(70,(int)(-10 * log10(1 - correctnessProbability)));
    }

    //
    // If the score is higher than the number of non-skipped seeds in the first pass, then there's a chance
    // that we didn't find genome locations with a score at least as low as what we found because of just
    // not having any matching seeds.
    //
    if (score > __max(bestPassSeedsNotSkipped, bestPassRCSeedsNotSkipped)) {
        baseMAPQ = __max(0,baseMAPQ - 10 * (score - __max(bestPassSeedsNotSkipped, bestPassRCSeedsNotSkipped)));
    }

    return baseMAPQ;

}
