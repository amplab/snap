/*++

Module Name:

    Hamming.cpp

Abstract:

    Hamming distance computation

Authors:

    Bill Bolosky, January, 2013

Environment:

    User mode service.

Revision History:

--*/

#include "stdafx.h"
#include "Hamming.h"
#include "BigAlloc.h"
#include "Compat.h"

static double *phred33ToProbability = NULL;
static double *perfectMatchProbability = NULL;

int ComputeHammingDistance(const char *text, const char *pattern, size_t len, const char *quality, int scoreLimit, double *mapProbability)
{
    const char *p = pattern;
    const char *t = text;
    const char *tend = text + len;
    int score = 0;
    *mapProbability = 1;

    while (t < tend) {
        unsigned _int64 diff = (*(const unsigned _int64 *)p) ^ (*(const unsigned _int64 *)t);
        if (diff != 0) {
            int offset = 0;
            size_t offsetForEndOfString = tend - t;
            do {
                //
                // Skip ahead to the next difference in the strings.
                //
                unsigned long zeroes;
                CountTrailingZeroes(diff, zeroes);
                unsigned matchingBytes = zeroes >> 3;
                offset += matchingBytes;
                if (offset >= offsetForEndOfString) {
                    break;
                }
                _ASSERT(offset < 8);
                diff = diff >> (matchingBytes << 3);    // Get rid of the matching bytes.
                diff = diff >> 8;                       // And get rid of one byte for the one we're scoring now.  DON'T COMBINE THIS WITH PREVIOUS LINE
                // The reason we need two sifts is because diff >> 64 is undefined (and in fact in VC returns diff, not 0 as you'd expect).

                _ASSERT(p[offset] != t[offset]);

                score++;
                if (score > scoreLimit) {
                    *mapProbability = 0;
                    return -1;
                }
                *mapProbability *= phred33ToProbability[quality[t-text + offset]];
                offset++;   // For the character we just compared.
            } while (diff != 0);
        }
        t += 8;
        p += 8;
    }
    if (0 == score) {
        *mapProbability = perfectMatchProbability[len];
    }
    return score;
}

void InitializeHammingDistance(double snpProbability)
{
    if (NULL == perfectMatchProbability) {
        perfectMatchProbability = (double *)BigAlloc(20000 * sizeof(double));
    }
    perfectMatchProbability[0] = 1;
    for (int i = 1;  i < 20000; i++) {
        perfectMatchProbability[i] = perfectMatchProbability[i-1] * (1.0 - snpProbability);
    }

    if (NULL != phred33ToProbability) {
        return;
    }

    phred33ToProbability = (double *)BigAlloc(sizeof(double) * 256);

    for (int i = 0; i < 33; i++) {
        phred33ToProbability[i] = snpProbability;  // This isn't a sensible Phred score
    }
    for (int i = 33; i <= 93 + 33; i++) {
         phred33ToProbability[i] = 1.0-(1.0 - pow(10.0,-1.0 * (i - 33.0) / 10.0)) * (1.0 - snpProbability);
    }
    for (int i = 93 + 33 + 1; i < 256; i++) {
        phred33ToProbability[i] = snpProbability;   // This isn't a sensible Phred score
    }
}