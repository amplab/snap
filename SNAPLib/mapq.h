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

#include "directions.h"

void initializeMapqTables();

double mapqToProbability(int mapq); // The probability of a match for the given MAPQ

inline int computeMAPQ(
    double probabilityOfAllCandidates,
    double probabilityOfBestCandidate,
    int score,
    int popularSeedsSkipped)
{
    probabilityOfAllCandidates = __max(probabilityOfAllCandidates, probabilityOfBestCandidate); // You'd think this wouldn't be necessary, but floating point limited precision causes it to be.
    _ASSERT(probabilityOfBestCandidate >= 0.0);
    // Special case for MAPQ 70, which we generate only if there is no evidence of a mismatch at all.

    // cheese is off, so no special casing MAPQ 70.  If you want to turn is back on, return these three lines and then change the two instance of 70 to 69 in the 
    // next set below (baseMAPQ =, twice).
    //
//    if (probabilityOfAllCandidates == probabilityOfBestCandidate && popularSeedsSkipped == 0 && score < 5) {
//        return 70;
//    }

    double correctnessProbability = probabilityOfBestCandidate / probabilityOfAllCandidates;
    int baseMAPQ;
    if (correctnessProbability >= 1) {
        baseMAPQ =  70;
    } else {
        baseMAPQ = __min(70, (int)(-10 * log10(1 - correctnessProbability)));
    }

    //
    // Apply a penalty based on the number of overly popular seeds in the read
    //
    baseMAPQ = __max(0, baseMAPQ - __max(0, popularSeedsSkipped-10) / 2);

#ifdef TRACE_ALIGNER
    printf("computeMAPQ called at %u: score %d, pThis %g, pAll %g, result %d\n",
            location, score, probabilityOfBestCandidate, probabilityOfAllCandidates, baseMAPQ);
#endif

    return baseMAPQ;
}
