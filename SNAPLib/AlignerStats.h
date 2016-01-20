/*++

Module Name:

    AlignerStats.h

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

#define     TIME_HISTOGRAM  0

#pragma once
#include "stdafx.h"
#include "Compat.h"

struct AbstractStats
{
    virtual ~AbstractStats();

    virtual void add(const AbstractStats* other) = 0;

    virtual void printHistograms(FILE* out) = 0;
};

//#define TIME_STRING_DISTANCE    1

struct AlignerStats : public AbstractStats
{
    AlignerStats(AbstractStats* i_extra = NULL);

    // TODO: This should also count both-aligned vs one-aligned etc.
    _int64 totalReads;
    _int64 uselessReads;    // Too short or too many Ns, so unalignable
    _int64 singleHits;
    _int64 multiHits;
    _int64 notFound;
    _int64 alignedAsPairs;
    _int64 lvCalls;
    _int64 filtered;
    _int64 extraAlignments;
    static const unsigned maxMapq = 70;
    unsigned mapqHistogram[maxMapq+1];

#if TIME_HISTOGRAM
    //
    // Histogram of alignment times.  Time buckets are divided by powers-of-two nanoseconds, so time bucket 0 is 
    // <= 1 ns, time bucket 10 is <= 1.024 us, etc.  Time bucket 30 is > 1s.
    //
    _int64 countByTimeBucket[31];
    _int64 nanosByTimeBucket[31];
#endif // TIME_HISTOGRAM


    static const unsigned maxMaxHits = 50;
    unsigned countOfBestHitsByWeightDepth[maxMaxHits];
    unsigned countOfAllHitsByWeightDepth[maxMaxHits];
    double probabilityMassByWeightDepth[maxMaxHits];

    AbstractStats* extra;

    virtual ~AlignerStats();

    virtual void add(const AbstractStats* other);

    virtual void printHistograms(FILE* out);
};

