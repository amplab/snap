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
#include "Compat.h"

void initializeMapqTables();

double mapqToProbability(int mapq);

inline int computeMAPQ(double probabilityOfAllCandidates, double probabilityOfBestCandidate)
{
    _ASSERT(probabilityOfAllCandidates >= probabilityOfBestCandidate);
    _ASSERT(probabilityOfBestCandidate >= 0.0);
    double correctnessProbability = probabilityOfBestCandidate / probabilityOfAllCandidates;
    if (correctnessProbability >= 1) {
        return 70;
    }
    return __min(70,(int)(-10 * log10(1 - correctnessProbability)));
}
