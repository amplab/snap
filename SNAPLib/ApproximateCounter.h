#pragma once

#include "Compat.h"

// Counts the number of distinct items in a stream approximately using Flajolet-Martin.
class ApproximateCounter
{
public:
    ApproximateCounter();

    void add(_uint64 value);

    unsigned getCount();

private:
    static const int SHIFT = 9;
    static const int BUCKETS = 1 << SHIFT;

    std::vector<_uint64> buckets;

    // MurmurHash3 finalization step from http://sites.google.com/site/murmurhash
    inline _uint64 hash(_uint64 value) {
        value ^= (value >> 33);
        value *= 0xff51afd7ed558ccdLL;
        value ^= (value >> 33);
        value *= 0xc4ceb9fe1a85ec53LL;
        value ^= (value >> 33);
        return value;
    }
};
