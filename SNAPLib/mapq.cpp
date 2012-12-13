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

#include "stdafx.h"
#include "mapq.h"

const int maxMAPQ = 70;

static double mapqToProbabilityTable[maxMAPQ+1];

void initializeMapqTables()
{
    for (int i = 0; i <= maxMAPQ; i++) {
        mapqToProbabilityTable[i] = pow(10.0,((double)i) / -10.0);
    }
}

double mapqToProbability(int mapq)
{
    _ASSERT(mapq >= 0 && mapq <= maxMAPQ);
    return mapqToProbabilityTable[mapq];
}
