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
#include "Compat.h"
#include "mapq.h"
#include "Compat.h"

const int maxMAPQ = 70;

static double mapqToProbabilityTable[maxMAPQ+1];
double firstOrderUnseenLocationProbability[maxSeedsForUnseenProbability+1];	// Chance that there are exactly nSeeds differences, one in each seed

void initializeMapqTables()
{
    for (int i = 0; i <= maxMAPQ; i++) {
        mapqToProbabilityTable[i] = pow(10.0,((double)i) / -10.0);
    }

	firstOrderUnseenLocationProbability[0] = 1;	// Not that we should ever look here, since this is the case where we used no seeds
	for (int i = 1; i < maxSeedsForUnseenProbability; i++) {
		firstOrderUnseenLocationProbability[i] = pow(GenericBaseDifferenceProbability,i) / (pow((double)i, i) / factorial(i));
	}
}

double mapqToProbability(int mapq)
{
    _ASSERT(mapq >= 0 && mapq <= maxMAPQ);
    return mapqToProbabilityTable[mapq];
}

double factorial(int n)
{
	double retVal = 1.0;
	for (int i = 2; i <= n; i++) {
		retVal *= n;
	}
	return retVal;
}
