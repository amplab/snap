/*++

Module Name:

    BloomFilter.cpp

Abstract:

    A simple Bloom Filter.

Authors:

    Bill Bolosky, Feburary, 2013

Environment:

    User mode service.

Revision History:

--*/

#include "stdafx.h"
#include "BloomFilter.h"

BloomFilter::BloomFilter(unsigned maxBitsEver_, unsigned baseDivisor_) : maxBitsEver(maxBitsEver_), baseDivisor(baseDivisor_), nFunctions(0)
{
    bits = (unsigned char *) BigAlloc((maxBitsEver + 63) / 8);    // Round up to a 64 bit boundary
}
    
BloomFilter::~BloomFilter()
{
    BigDealloc(bits);
}

    void 
BloomFilter::init(unsigned nFunctions_, unsigned maxBits_)
{
    nFunctions = nFunctions_;
    maxBits = maxBits_;
    _ASSERT(maxBits <= maxBitsEver);

    //
    // Only initialize enough bits for this run.
    //
    memset(bits,0,(maxBits + 63) / 8);
}

    void 
BloomFilter::addToSet(unsigned value)
{
    unsigned effectiveValue = value / baseDivisor;

    for (unsigned i = 0 ; i < nFunctions; i++) {
        unsigned bitToSet = bit_rotate_right(effectiveValue, i) % maxBits;
        bits[arrayOffsetForBit(bitToSet)] |= 1 << bitOffsetForBit(bitToSet);
    }
}

    bool 
BloomFilter::mightValueBeInSet(unsigned value)
{
    unsigned effectiveValue = value / baseDivisor;

    bool anyBitsMissed = false;
    for (unsigned i = 0; i < nFunctions; i++) {
        unsigned bitToCheck = bit_rotate_right(effectiveValue, i);
        anyBitsMissed &= (bits[arrayOffsetForBit(bitToCheck)] & (1 << bitOffsetForBit(effectiveValue))) == 0;
    }

    return !anyBitsMissed;
}
