#include "stdafx.h"
#include "ProbabilityDistance.h"
#include "Compat.h"


#ifdef TRACE_PROBABILITY_DISTANCE
#define TRACE printf
#else
#define TRACE(...) {}
#endif


namespace {
    inline double max3(double d1, double d2, double d3) {
        if (d1 > d2) {
            return (d1 > d3) ? d1 : d3;
        } else {
            return (d2 > d3) ? d2 : d3;
        }
    }
}


ProbabilityDistance::ProbabilityDistance(double snpProb, double gapOpenProb, double gapExtensionProb)
{
    snpLogProb = log(snpProb);
    gapOpenLogProb = log(gapOpenProb);
    gapExtensionLogProb = log(gapExtensionProb);

    // Fill in the matchLogProb and mismatchLogProb tables for base quality values; assumes Phred+33 encoding
    for (int q = 0; q < 256; q++) {
        if (q < 33) {
            matchLogProb[q] = NO_PROB;
            mismatchLogProb[q] = NO_PROB;
        } else {
            // Probability that we misread this base
            double errorProb = pow(10.0, -(q - 33) / 10.0);
            // A match occurs if we didn't misread the base and it wasn't a SNP (technically it could also
            // be that we misread it and it *was* a SNP, but that's pretty unlikely)
            double matchProb = (1.0 - errorProb) * (1.0 - snpProb);
            double mismatchProb = 1.0 - matchProb;
            matchLogProb[q] = log(matchProb);
            mismatchLogProb[q] = log(mismatchProb);
            if (q >= 33 && q <= 'J') {
                TRACE("q=%c: match %.04f, mismatch %.04f\n", q, matchProb, mismatchProb);
            }
            // TODO: Might be good to check the numerical stability of this (1 - small) stuff
        }
    }
}


int ProbabilityDistance::compute(
        const char *reference,
        const char *read,
        const char *quality,
        int readLen,
        int maxStartShift,
        int maxShift,               // Maximum overall shift to consider
        double *matchProbability)
{
    _ASSERT(maxStartShift < MAX_SHIFT);
    _ASSERT(maxShift < MAX_SHIFT);
    _ASSERT(maxStartShift <= maxShift);

    // Fill in the readPos = 0 row to allow us to start only at -maxStartShift..+maxStartShift
    for (int s = -maxShift-1; s <= maxShift+1; s++) {
        d[0][MAX_SHIFT+s][READ_GAP] = NO_PROB;
        d[0][MAX_SHIFT+s][REF_GAP] = NO_PROB;
        if (s < -maxStartShift || s > maxStartShift) {
            d[0][MAX_SHIFT+s][NO_GAP] = NO_PROB;
        } else {
            d[0][MAX_SHIFT+s][NO_GAP] = log(1.0);
        }
    }

    // Now go through each readPos from 1 to readLen and compute how to best get there
    for (int r = 1; r <= readLen; r++) {
        // Add sentinels at the end of the array
        d[r][MAX_SHIFT-maxShift-1][READ_GAP] = NO_PROB;
        d[r][MAX_SHIFT+maxShift+1][READ_GAP] = NO_PROB;
        d[r][MAX_SHIFT-maxShift-1][REF_GAP] = NO_PROB;
        d[r][MAX_SHIFT+maxShift+1][REF_GAP] = NO_PROB;
        d[r][MAX_SHIFT-maxShift-1][NO_GAP] = NO_PROB;
        d[r][MAX_SHIFT+maxShift+1][NO_GAP] = NO_PROB;

        // Fill in the rest of the values using dynamic program recurrence
        for (int s = -maxShift; s <= maxShift; s++) {
            // The NO_GAP case; we get here either from a previous NO_GAP or by closing a gap from the
            // previous readPos, and in either case, we need to match the current base
            double thisBaseProb = (read[r-1] == reference[r-1+s]) ? matchLogProb[quality[r-1]] : mismatchLogProb[quality[r-1]];
            d[r][MAX_SHIFT+s][NO_GAP] = max3(d[r-1][MAX_SHIFT+s][NO_GAP] + thisBaseProb,
                                             d[r-1][MAX_SHIFT+s][REF_GAP] + thisBaseProb,
                                             d[r-1][MAX_SHIFT+s][READ_GAP] + thisBaseProb);

            // The READ_GAP case; we can either open a new gap from the previous NO_GAP or REF_GAP cases, or
            // extend a gap computed in the previous READ_GAP case
            d[r][MAX_SHIFT+s][READ_GAP] = max3(d[r-1][MAX_SHIFT+s+1][NO_GAP] + gapOpenLogProb,
                                               d[r-1][MAX_SHIFT+s+1][REF_GAP] + gapOpenLogProb,
                                               d[r-1][MAX_SHIFT+s+1][READ_GAP] + gapExtensionLogProb);

            // The REF_GAP case; we can either open a new gap from NO_GAP/READ_GAP, or extend one
            d[r][MAX_SHIFT+s][REF_GAP] = max3(d[r][MAX_SHIFT+s-1][NO_GAP] + gapOpenLogProb,
                                              d[r][MAX_SHIFT+s-1][REF_GAP] + gapExtensionLogProb,
                                              d[r][MAX_SHIFT+s-1][READ_GAP] + gapOpenLogProb);
        }
    }

#ifdef TRACE_PROBABILITY_DISTANCE
    printf("Here is the final matrix:\n");
    for (int r = 0; r <= readLen; r++) {
        printf("%d: ", r);
        for (int g = 0; g < 3; g++) {
            for (int s = -maxShift; s <= maxShift; s++) {
                printf("%7.2g ", d[r][MAX_SHIFT+s][g]);
            }
            if (g < 2) {
                printf("| ");
            }
        }
        printf("\n");
    }
#endif

    // Return the best probability, and a somewhat arbitrary score for it (TODO: need to actually compute # of edits)
    double best = NO_PROB;
    for (int s = -maxShift; s <= maxShift; s++) {
        for (int g = 0; g < 3; g++) {
            best = __max(best, d[readLen][MAX_SHIFT+s][g]);
        }
    }
    *matchProbability = exp(best);
    TRACE("Best match probability: %g (log: %.2g)\n", exp(best), best);
    return 5;
}
