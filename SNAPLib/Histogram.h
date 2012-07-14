/*++

Module Name:

    histogram.h

Abstract:

    Header for cheezy histogram class

Authors:

    Bill Bolosky, September, 2011

Environment:

    User mode service.

    This class is NOT thread safe.  It's the caller's responsibility to ensure that
    at most one thread uses an instance at any time.

Revision History:

--*/

#pragma once

class Histogram {
public:
    Histogram(unsigned i_nBuckets, bool i_isExponential);

    ~Histogram();

    unsigned getNBuckets() const {return nBuckets;}
    bool getIsExponential() const {return isExponential;}
    unsigned getBucketCount(unsigned whichBucket) const;
    unsigned getBucketMax(unsigned whichBucket) const;
    unsigned getBucketMin(unsigned whichBucket) const;

    void print() const;

    void addToCount(unsigned value, unsigned amountToAdd = 1);

private:

    unsigned nBuckets;
    bool isExponential;

    struct Bucket {
        unsigned    maxValue;
        unsigned    count;
    };

    Bucket *buckets;
};
