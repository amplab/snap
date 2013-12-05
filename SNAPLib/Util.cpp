/*++

Module Name:

    Util.cpp

Abstract:

    Generic support routines that don't seem to belong elsewhere.

Authors:

    Bill Bolosky, March, 2013

Environment:

    User mode service.

Revision History:

    Factored from other places

--*/

#include "stdafx.h"
#include "Util.h"


int FirstPowerOf2GreaterThanOrEqualTo(int value)
{
    int highestBitSet;
    for (highestBitSet = 0; highestBitSet <= 30; highestBitSet++) { // Only go to 31, since this is signed
        if (!(value & ~((1 << highestBitSet) - 1))) {
            highestBitSet -= 1;
            break;
        }
    }
    if (1 << highestBitSet == value) return value;
    return 1 << (highestBitSet + 1);
}

    void
util::memrevcpy(
    void* dst,
    const void* src,
    size_t bytes)
{
    size_t dwords = bytes >> 3;
    _uint64* p = (_uint64*) dst;
    const _uint64* q = (const _uint64*) ((const char*)src + bytes - 8);
    for (size_t i = 0; i < dwords; i++) {
        *p++ = ByteSwapUI64(*q--);
    }
    int left = (int) (bytes & 7);
    if (left > 0) {
        char* p2 = (char*) p;
        const char* q2 = (left - 1) + (const char*) src;
        for (int i = 0; i < left; i++) {
            *p2++ = *q2--;
        }
    }
}

