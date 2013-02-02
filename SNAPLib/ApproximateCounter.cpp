#include "stdafx.h"
#include <math.h>
#include "ApproximateCounter.h"

using namespace std;


ApproximateCounter::ApproximateCounter()
{
    buckets.resize(BUCKETS);
}

void ApproximateCounter::add(_uint64 value)
{
    _uint64 h = hash(value);
    unsigned bucket = (unsigned) h % BUCKETS;
    unsigned rest = (unsigned)(h >> SHIFT);
    unsigned long firstZero;
    if (rest == 0) {
        firstZero = 64 - SHIFT;
    } else {
        CountTrailingZeroes(rest, firstZero);
    }

    buckets[bucket] |= (1LL << firstZero);
}


unsigned ApproximateCounter::getCount()
{
    double s = 0;
    for (int i = 0; i < BUCKETS; i++) {
        int r = 0;
        while (r < 64 && (buckets[i] & (1LL << r)) != 0) {
            r++;
        }
        s += r;
    }
    return (unsigned) (BUCKETS / 0.77351 * pow(2, s / BUCKETS));
}
