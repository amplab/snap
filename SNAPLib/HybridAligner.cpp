/*++

Module Name:

    HybridAligner.h

Abstract:

    Header for SNAP genome aligner

Authors:

    Bill Bolosky, October, 2011

Environment:

    User mode service.

    This class is NOT thread safe.  It's the caller's responsibility to ensure that
    at most one thread uses an instance at any time.

Revision History:

    Adapted from Matei Zaharia's Scala implementation.

--*/

#include "stdafx.h"
#include "HybridAligner.h"

HybridAligner::HybridAligner(
    GenomeIndex *genomeIndex, 
    unsigned     confDiff, 
    unsigned     maxHitsToConsider, 
    unsigned     maxK,
    unsigned     maxReadSize,
    unsigned     maxSeedsToUse,
    unsigned     lvCutoff)
{
    baseAligner = new BaseAligner(genomeIndex,confDiff,maxHitsToConsider,maxK,maxReadSize,maxSeedsToUse,lvCutoff,100000);
    intersectingAligner = new IntersectingAligner(genomeIndex,confDiff,280000,5,maxReadSize,20,lvCutoff);

    if (NULL == baseAligner || NULL == intersectingAligner) {
        fprintf(stderr,"HybridAligner: Unable to allocate underlying aligner\n");
        exit(1);
    }
}

    AlignmentResult
HybridAligner::AlignRead(
    Read        *read,
    unsigned    *genomeLocation,
    bool        *hitIsRC,
    int         *score)
{
    unsigned    genomeLocation1;
    bool        hitIsRC1;
    int         score1;
    AlignmentResult result1 = baseAligner->AlignRead(read,&genomeLocation1,&hitIsRC1,&score1);

    if (CertainHit == result1 /*|| SingleHit == result1 && score1 <= 4*/) {
        *genomeLocation = genomeLocation1;
        *hitIsRC = hitIsRC1;
        if (NULL != score) {
            *score = score1;
        }
        return result1;
    }

    unsigned    genomeLocation2;
    bool        hitIsRC2;
    int         score2;
    AlignmentResult result2 = intersectingAligner->AlignRead(read,&genomeLocation2,&hitIsRC2,&score2);

    if (SingleHit == result1 && SingleHit == result2 && genomeLocation1 != genomeLocation2 &&
        (score1 <= score2 && score1 + baseAligner->getConfDiff() > score2 ||
         score2 <= score1 && score2 + baseAligner->getConfDiff() > score1)) {
        //
        // Both aligners had single hits at different locations with scores within
        // confDiff of one another.  Declare a MultiHit.
        //
        return MultipleHits;
    }

    if (CertainHit == result2 || 
        SingleHit == result2 && (SingleHit == result1 && score2 < score1 || SingleHit != result1)) {
        *genomeLocation = genomeLocation2;
        *hitIsRC = hitIsRC2;
        if (NULL != score) {
            *score = score2;
        }
        return result2;
    }

    *genomeLocation = genomeLocation1;
    *hitIsRC = hitIsRC1;
    if (NULL != score) {
        *score = score1;
    }

    return result1;
}
    



