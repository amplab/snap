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
    _int64 usefulReads;
    _int64 singleHits;
    _int64 multiHits;
    _int64 notFound;
    _int64 errors;
    static const unsigned maxMapq = 70;
    unsigned mapqHistogram[maxMapq+1];
    unsigned mapqErrors[maxMapq+1];

    static const unsigned maxMaxHits = 50;
    unsigned countOfBestHitsByWeightDepth[maxMaxHits];
    unsigned countOfAllHitsByWeightDepth[maxMaxHits];
    double probabilityMassByWeightDepth[maxMaxHits];

#ifdef  TIME_STRING_DISTANCE
    _int64 nanosTimeInBSD[2][2]; // [nearby/notNearby][aligned/notAligned]
    _int64 BSDCounts[2][2]; // same
    _int64 hammingCount;    // Number of times used Hamming distance rather than BSD/LV
    _int64 hammingNanos;
#endif  // TIME_STRING_DISTANCE

    AbstractStats* extra;

    virtual ~AlignerStats();

    virtual void add(const AbstractStats* other);

    virtual void printHistograms(FILE* out);
};

