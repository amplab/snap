/*++

Module Name:

    histogram.cpp

Abstract:

    Cheezy histogram class

Authors:

    Bill Bolosky, September, 2011

Environment:

    User mode service.

    This class is NOT thread safe.  It's the caller's responsibility to ensure that
    at most one thread uses an instance at any time.

Revision History:

--*/

#include "stdafx.h"
#include "Compat.h"
#include "Histogram.h"

Histogram::Histogram(unsigned i_nBuckets, bool i_isExponential) :
    nBuckets(i_nBuckets), isExponential(i_isExponential)
{
    _ASSERT(nBuckets > 0);

    buckets = new Bucket[nBuckets];
    if (NULL == buckets) {
        fprintf(stderr,"Histogram: unable to allocate buckets.\n");
        exit(1);
    }

    buckets[0].maxValue = 1;
    buckets[0].count = 0;

    for (unsigned i = 1; i < nBuckets; i++) {
        if (isExponential) {
            buckets[i].maxValue = buckets[i-1].maxValue * 2;
        } else {
            buckets[i].maxValue = buckets[i-1].maxValue +1;
        }
        buckets[i].count = 0;
    }
}

    unsigned
Histogram::getBucketCount(unsigned whichBucket) const
{
    _ASSERT(whichBucket < nBuckets);
    return buckets[whichBucket].count;
}

    unsigned
Histogram::getBucketMax(unsigned whichBucket) const
{
    _ASSERT(whichBucket < nBuckets);
    return buckets[whichBucket].maxValue;
}

    unsigned
Histogram::getBucketMin(unsigned whichBucket) const
{
    _ASSERT(whichBucket < nBuckets);
    if (0 == whichBucket) {
        return 0;
    }
    return buckets[whichBucket-1].maxValue + 1;
}

    void
Histogram::addToCount(unsigned value, unsigned amountToAdd)
{
    for (unsigned i = 0 ; i < nBuckets; i++) {  // It is called "cheezy" after all
        if (value <= buckets[i].maxValue) {
            buckets[i].count += amountToAdd;
            return;
        }
    }
    //
    // Overflow.  Just drop it.
    //
}

    void
Histogram::print() const
{
    printf("MaxValue    Count\n");
    printf("-------- --------\n");
    for (unsigned i = 0; i < nBuckets; i++) {
        printf("%8d %8d\n",buckets[i].maxValue, buckets[i].count);
    }
}

Histogram::~Histogram()
{
    delete [] buckets;
    buckets = NULL;
}
