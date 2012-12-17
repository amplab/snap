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

inline int computeMAPQ(double probabilityOfAllCandidates, double probabilityOfBestCandidate, int score)
{
    _ASSERT(probabilityOfAllCandidates >= probabilityOfBestCandidate);
    _ASSERT(probabilityOfBestCandidate >= 0.0);
    double correctnessProbability = probabilityOfBestCandidate / probabilityOfAllCandidates;
    int baseMAPQ;
    if (correctnessProbability >= 1) {
        baseMAPQ =  70;
    } else {
        baseMAPQ = __min(70,(int)(-10 * log10(1 - correctnessProbability)));
    }

    return baseMAPQ;
    //
    // Apply a penalty based on the absolute difference between the read and the place it matched, as expressed
    // by its score.
    //
    if (3 * score >= baseMAPQ) {
        return 0;
    } else {
        return baseMAPQ - 3 * score;
    }
}
