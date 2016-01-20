/*++

Module Name:

    AlignerStats.cpp

Abstract:

    Common statistics for running single & paired alignment.

Authors:

    Ravi Pandya, May, 2012

Environment:
`
    User mode service.

Revision History:

    Integrated from SingleAligner.cpp & PairedAligner.cpp

--*/

#include "stdafx.h"
#include "options.h"
#include "AlignerStats.h"

AbstractStats::~AbstractStats()
{}

AlignerStats::AlignerStats(AbstractStats* i_extra)
:
    totalReads(0),
    uselessReads(0),
    singleHits(0), 
    multiHits(0),
    notFound(0),
    alignedAsPairs(0),
    extra(i_extra),
    lvCalls(0),
    filtered(0),
    extraAlignments(0)
{
    for (int i = 0; i <= AlignerStats::maxMapq; i++) {
        mapqHistogram[i] = 0;
     }

    for (int i = 0; i < maxMaxHits; i++) {
        countOfBestHitsByWeightDepth[i] = 0;
        countOfAllHitsByWeightDepth[i] = 0;
        probabilityMassByWeightDepth[i] = 0;
    }

#if     TIME_HISTOGRAM
    for (unsigned i = 0; i < 31; i++) {
        countByTimeBucket[i] = nanosByTimeBucket[i] = 0;
    }
#endif  // TIME_HISTOGRAM

}

AlignerStats::~AlignerStats()
{
    if (extra != NULL) {
        delete extra;
    }
}

    void
AlignerStats::printHistograms(
    FILE* out)
{
    // nothing
    if (extra != NULL) {
        extra->printHistograms(out);
    }
}

    void
AlignerStats::add(
    const AbstractStats* i_other)
{
    AlignerStats* other = (AlignerStats*) i_other;
    totalReads += other->totalReads;
    uselessReads += other->uselessReads;
    singleHits += other->singleHits;
    multiHits += other->multiHits;
    notFound += other->notFound;
    alignedAsPairs += other->alignedAsPairs;
    lvCalls += other->lvCalls;
    filtered += other->filtered;
    extraAlignments += other->extraAlignments;

    if (extra != NULL && other->extra != NULL) {
        extra->add(other->extra);
    }

    for (int i = 0; i <= AlignerStats::maxMapq; i++) {
        mapqHistogram[i] += other->mapqHistogram[i];
    }

    for (int i = 0; i < maxMaxHits; i++) {
        countOfBestHitsByWeightDepth[i] += other->countOfBestHitsByWeightDepth[i];
        countOfAllHitsByWeightDepth[i] += other->countOfAllHitsByWeightDepth[i];
        probabilityMassByWeightDepth[i] = other->probabilityMassByWeightDepth[i];
    }

#if TIME_HISTOGRAM
    for (unsigned i = 0; i < 31; i++) {
        countByTimeBucket[i] += other->countByTimeBucket[i];
        nanosByTimeBucket[i] += other->nanosByTimeBucket[i];
    }
#endif // TIME_HISTOGRAM
}
