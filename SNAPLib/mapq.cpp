/*++

Module Name:

    mapq.cpp

Abstract:

    Support functions for mapping quality

Authors:

    Bill Bolosky, December, 2012

Environment:

    User mode service.

Revision History:

 
--*/

#include "stdafx.h"
#include "Compat.h"
#include "mapq.h"

const int maxMAPQ = 70;

static double mapqToProbabilityTable[maxMAPQ+1];

void initializeMapqTables()
{
    mapqToProbabilityTable[0] = .1;  // This should technically be 0, but in practice it's a little better than that, so leave some chance here.
    for (int i = 1; i <= maxMAPQ; i++) {
        mapqToProbabilityTable[i] = 1- pow(10.0,((double)i) / -10.0);
    }

}

double mapqToProbability(int mapq)
{
    _ASSERT(mapq >= 0 && mapq <= maxMAPQ);
    return mapqToProbabilityTable[mapq];
}
