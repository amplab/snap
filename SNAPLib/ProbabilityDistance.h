#pragma once

#include "Read.h"

//
// Similar to BoundedStringDistance and LandauVishkin, but computes the probability of a read
// string being generated from a reference sequence given an error model, mutation model and
// the quality scores of the bases in the read. For now, we assume that the statistical model
// for errors is one with a "gap open probability" and a fixed "gap extension probability"
// per base. For substitions, we could have different probabilities for each transition, but
// we assume that they all have the same probability right now.
//
class ProbabilityDistance {
public:
    static const int MAX_READ = MAX_READ_LENGTH;
    static const int MAX_SHIFT = 20;

    ProbabilityDistance(double snpProb, double gapOpenProb, double gapExtensionProb);

    int compute(
            const char *reference,
            const char *read,
            const char *quality,
            int readLen,
            int maxStartShift,
            int maxTotalShift,
            double *matchProbability);

private:
    double snpLogProb;
    double gapOpenLogProb;
    double gapExtensionLogProb;

    double matchLogProb[256];      // [baseQuality]
    double mismatchLogProb[256];   // [baseQuality]

#define NO_PROB  -1000000.0;  // A really negative log probability -- basically zero.  VC compiler won't allow static const double in a class.

    enum GapStatus { NO_GAP, READ_GAP, REF_GAP };

    // d[readPos][shift][gapStatus] is the best possible log probability for aligning the
    // substring read[0..readPos] to reference[?..readPos + shift]. The "?" in reference is
    // because we allow starting an alignment from reference[-maxStartShift..maxStartShift]
    // instead of just reference[0], to deal with indels toward the start of the read.
    double d[MAX_READ][2*MAX_SHIFT+1][3];   // [readPos][shift][gapStatus]

    // A state in the D array, used for backtracking pointers
    struct State {
        int readPos;
        int shift;
        int gapStatus;
    };

    // Previous state in our dynamic program, for backtracking to print CIGAR strings
    State prev[MAX_READ][2*MAX_SHIFT+1][3];   // [readPos][shift][gapStatus]
};
