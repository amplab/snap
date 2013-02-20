/*++

Module Name:

    BloomFilter.h

Abstract:

    Headers for a simple Bloom Filter.

Authors:

    Bill Bolosky, Feburary, 2013

Environment:

    User mode service.

Revision History:

--*/

#pragma once

#include "BigAlloc.h"
#include "Compat.h"

class BloomFilter {
public:
    BloomFilter(unsigned maxBitsEver_, unsigned baseDivisor_);
    ~BloomFilter();

    void init(unsigned nFunctions_, unsigned maxBits_);

    void addToSet(unsigned value);

    bool mightValueBeInSet(unsigned value);
private:

    unsigned        maxBitsEver;
    unsigned        maxBits;
    unsigned        nFunctions;
    unsigned        baseDivisor;

    unsigned char   *bits;

    inline unsigned arrayOffsetForBit(unsigned whichBit) {return whichBit / 8;}
    inline unsigned char bitOffsetForBit(unsigned whichBit) {return whichBit % 8;}
};