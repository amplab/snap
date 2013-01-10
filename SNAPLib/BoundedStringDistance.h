#pragma once

#include "Compat.h"
#include "BigAlloc.h"

/**
 * A generalization of Ukkonen's bounded edit distance algorithm for affine gap penalties.
 * To keep the code simple, we write one version of the class that supports computing or
 * not computing the edit string through a template parameter, and encodes other constants
 * as template parameters too. Hopefully the compiler will optimize this properly.
 *
 * This algorithm uses an affine gap penalty with a cost of gapOpenPenalty for the first
 * base in an insertion/deletion and an additional cost of 1 for each extension. It would
 * be possible to make the mismatch and gap extension penalties configurable too but it
 * doesn't seem common to make them too different from each other.
 *
 * Parameters:
 * - MAX_DISTANCE: maximum distance score supported.
 * - MAX_SHIFT: maximum shift supported. This can be less than MAX_DISTANCE if one
 *   is not expecting huge indels, in order to save computation time.
 *
 * TODO: Allow making the L array "cyclic" through another template parameter (i.e.
 * so that only the last few rows will be kept) to save cache space for big arrays.
 */
template<bool ALLOW_START_SHIFTS=false, int MAX_DISTANCE=31, int MAX_SHIFT=MAX_DISTANCE>
class BoundedStringDistance
{
private:
    enum GapStatus { NO_GAP, TEXT_GAP, PATTERN_GAP };

    static const int MAX_READ = 10000;           // Arbitrary, could easily be bigger
    static const int SHIFT_OFF = MAX_SHIFT + 1;  // For indexing into arrays by shift

    // Edit distance penalties for substitutions and gap open; we assume that gap extend is 1
    int substitutionPenalty;
    int gapOpenPenalty;       // We assume that the gap extend penalty is 1

    // MAPQ probabilities for substitutions and indels from the reference genome. We assume
    // that indels have geometrically distributed lengths once they are open (i.e. the same
    // probability of extension at each step).
    double snpProbability;
    double indelOpenProbability;
    double indelExtendProbability;

    // L[distance][shift][gapStatus] is the largest value x such that the distance
    // between pattern[0..x] and text[0..x+shift] is "distance" and the gap status
    // at the end is "gapStatus". To avoid branches, we keep L as a pointer into
    // a bigger array, padded, with enough space for negative distance indices.
    int padded[2*MAX_DISTANCE+1][2*MAX_SHIFT+3][3];
    int (*L)[2*MAX_SHIFT+3][3];
   
    // Previous gap status we had before getting to each particular state; this is enough
    // to infer the location we came from and the action we took to get here
    int prevGapStatus[MAX_DISTANCE+1][2*MAX_SHIFT+3][3];

    // Mismatch probability as a function of base quality
    double mismatchProbability[256];

    // Match probability as a function of number of bases matching (ignores base quality for now)
    double matchProbability[MAX_READ+1];

public:
    BoundedStringDistance(
            int substitutionPenalty_,
            int gapOpenPenalty_,
            double snpProbability_ = 0.001,
            double indelOpenProbability_ = 0.001,
            double indelExtendProbability_ = 0.2)
        : substitutionPenalty(substitutionPenalty_), gapOpenPenalty(gapOpenPenalty_),
          snpProbability(snpProbability_), indelOpenProbability(indelOpenProbability_),
          indelExtendProbability(indelExtendProbability_)
    {
        _ASSERT(substitutionPenalty <= MAX_DISTANCE);
        _ASSERT(gapOpenPenalty <= MAX_DISTANCE);
        L = padded + MAX_DISTANCE;

        // Clear out the L array
        for (int d = -MAX_DISTANCE; d <= MAX_DISTANCE; d++) {
            for (int s = 0; s < 2 * MAX_SHIFT + 3; s++) {
                for (int g = 0; g < 3; g++) {
                    // Use a very negative value to avoid having to check for when a parent state
                    // has L = 0 in compute, avoiding expensive branches
                    L[d][s][g] = -2 * MAX_DISTANCE;
                }
            }
        }
        
        // Add sentinel values for prevGapStatus
        for (int s = 0; s < 2 * MAX_SHIFT + 3; s++) {
            for (int g = 0; g < 3; g++) {
                prevGapStatus[0][s][g] = -1;
            }
        }
        
        // Fill in probability tables for fast computation of MAPQ; assumes Phred+33 qualities
        for (int q = 0; q < 256; q++) {
            if (q < 33) {
                mismatchProbability[q] = 1.0;
            } else {
                // A correct match will happen only if we don't read the base incorrectly *and* there's
                // no SNP, so invert that to get the mismatch probability. (Ignores the possibility that
                // a SNP actually got misread as the correct base for now.)
                double baseProb = pow(10.0, -0.1 * (q - 33));
                mismatchProbability[q] = 1.0 - (1.0 - baseProb) * (1.0 - snpProbability);
            }
        }
        matchProbability[0] = 1.0;
        for (int len = 1; len <= MAX_READ; len++) {
            matchProbability[len] = matchProbability[len-1] * (1.0 - snpProbability);
        }
    }

    /**
     * Compute the string distance between text[0...] and pattern[0..patternLen], if
     * it is less than limit, or return -1 if it is greater. Optionally, write the
     * edit string (in CIGAR format) into editStringBuf, write the map probability in
     * mapProbability, and write the shift at the end of the alignment into endShift
     * if these pointers are not null. For map probability, qualityString must also
     * not be null. For edit strings, one can specify whether to use 'M' characters.
     *
     * Assumptions:
     * - text is indexable up to patternLen + MAX_SHIFT - 1 (so it's a pointer
     *   into a longer string that we don't need to bounds-check).
     * - pattern is indexable up to 7 characters past its length.
     * - editStringBuf is long enough to hold the edit string, if set. We check that
     *   it's at least 2 * patternLen in case we must make one edit per character.
     */
    int compute(const char *text, const char *pattern, int patternLen, int maxStartShift, int limit, int *endShift,
            double *mapProbability, const char *qualityString, char *editStringBuf, int editStringBufLen, bool useM)
    {
        if (limit < 0 || limit > MAX_DISTANCE) {
            fprintf(stderr, "Invalid distance limit: %d\n", limit);
            exit(1);
        }
        if (mapProbability != NULL && qualityString == NULL) {
            fprintf(stderr, "Non-null mapProbability given but qualityString is NULL\n");
            exit(1);
        }
        if (!ALLOW_START_SHIFTS && maxStartShift != 0) {
            fprintf(stderr, "Non-zero maxStartShift given, but ALLOW_START_SHIFTS is false\n");
            exit(1);
        }

        const char *patternEnd = pattern + patternLen;

        // For distances less than gapOpenPenalty, just search for mismatches with shift = 0
        for (int s = SHIFT_OFF-maxStartShift; s <= SHIFT_OFF+maxStartShift; s++) {
            int pos = 0;
            for (int d = 0; d <= __min(limit, gapOpenPenalty - 1); d += substitutionPenalty) {
                pos += longestPrefix(text+pos+s-SHIFT_OFF, pattern+pos, patternEnd);
                L[d][s][NO_GAP] = pos;
                if (d > 0) {
                    prevGapStatus[d][s][NO_GAP] = NO_GAP;
                }
                if (pos == patternLen) {
                    if (!writeEditString(d, 0, editStringBuf, editStringBufLen, useM)) {
                        return -2;
                    }
                    computeMapProbability(d, 0, mapProbability, qualityString);
                    if (endShift != NULL) {
                        *endShift = s - SHIFT_OFF;
                    }
                    return d;
                } else {
                    pos += 1; // Skip the mismatch so we can go to the next d
                }
            }
        }

        if (ALLOW_START_SHIFTS) {
            // Add back sentinel values in case we were called with a bigger maxStartShift earlier
            for (int d = 0; d <= __min(limit, gapOpenPenalty - 1); d++) {
                for (int g = 0; g < 3; g++) {
                    L[d][SHIFT_OFF-maxStartShift-1][g] = -2 * MAX_DISTANCE;
                    L[d][SHIFT_OFF-maxStartShift-2][g] = -2 * MAX_DISTANCE;
                    L[d][SHIFT_OFF-maxStartShift-3][g] = -2 * MAX_DISTANCE;
                    L[d][SHIFT_OFF+maxStartShift+1][g] = -2 * MAX_DISTANCE;
                    L[d][SHIFT_OFF+maxStartShift+2][g] = -2 * MAX_DISTANCE;
                    L[d][SHIFT_OFF+maxStartShift+3][g] = -2 * MAX_DISTANCE;
                }
            }
        }

        // For distances >= gapOpenPenalty, also allow creating gaps in the text or pattern
        for (int d = gapOpenPenalty; d <= limit; d++) {
            int maxShift = __min(maxStartShift + d - gapOpenPenalty + 1, MAX_SHIFT);
            if (ALLOW_START_SHIFTS) {
                for (int g = 0; g < 3; g++) {
                    L[d][SHIFT_OFF-maxShift-1][g] = -2 * MAX_DISTANCE;
                    L[d][SHIFT_OFF-maxShift-2][g] = -2 * MAX_DISTANCE;
                    L[d][SHIFT_OFF-maxShift-3][g] = -2 * MAX_DISTANCE;
                    L[d][SHIFT_OFF+maxShift+1][g] = -2 * MAX_DISTANCE;
                    L[d][SHIFT_OFF+maxShift+2][g] = -2 * MAX_DISTANCE;
                    L[d][SHIFT_OFF+maxShift+3][g] = -2 * MAX_DISTANCE;
                }
            }
            for (int s = SHIFT_OFF-maxShift; s <= SHIFT_OFF+maxShift; s++) {
                // See whether we can start / extend the "text gap" case
                fillBest(d, s, TEXT_GAP,
                         d-gapOpenPenalty, s-1, NO_GAP, 0,  // Open a gap from NO_GAP
                         d-1, s-1, TEXT_GAP, 0);            // Continue a TEXT_GAP
                // Note that continuing from a PATTERN_GAP is covered because the NO_GAP
                // case at d-gapOpenPenalty considered closing the PATTERN_GAP there

                // See whether we can start / extend the "pattern gap" case
                fillBest(d, s, PATTERN_GAP,
                         d-gapOpenPenalty, s+1, NO_GAP, 1,  // Open a new gap
                         d-1, s+1, PATTERN_GAP, 1);         // Continue a PATTERN_GAP
                // Note that continuing from a TEXT_GAP is covered because the NO_GAP
                // case at d-gapOpenPenalty considered closing the TEXT_GAP there

                // Find the best way to continue the "neither in gap" case; this could be either
                // adding a substitution from the previous "neither in gap" state, or closing a
                // text gap, or a pattern gap. Note that since we compute the best ways to end up
                // in TEXT_GAP and PATTERN_GAP above, this automatically includes using those as
                // starting points and just closing the gap that got us to this distance.
                fillBest(d, s, NO_GAP,
                         d-substitutionPenalty, s, NO_GAP, 1,  // Add subsitution after NO_GAP
                         d, s, TEXT_GAP, 0,                    // Close a TEXT_GAP and go from there
                         d, s, PATTERN_GAP, 0);                // Close a PATTERN_GAP
                // For NO_GAP, also add any string of matching characters after where we end up
                int location = L[d][s][NO_GAP];
                if (location >= 0) {
                    location += longestPrefix(text+location+s-SHIFT_OFF, pattern+location, patternEnd);
                    L[d][s][NO_GAP] = location;
                }

                if (location == patternLen) {
                    // Reached end of pattern; we're done!
                    if (!writeEditString(d, s - SHIFT_OFF, editStringBuf, editStringBufLen, useM)) {
                      return -2;
                    }
                    computeMapProbability(d, s - SHIFT_OFF, mapProbability, qualityString);
                    if (endShift != NULL) {
                        *endShift = s - SHIFT_OFF;
                    }
                    return d;
                }
            }
        }

        return -1;
    }

    // Just compute distance, without MAPQ or edit string
    int compute(const char *text, const char *pattern, int patternLen, int maxStartShift, int limit) {
        return compute(text, pattern, patternLen, maxStartShift, limit, NULL, NULL, NULL, NULL, 0, false);
    }

    // Compute distance and MAPQ but not edit string
    int compute(const char *text, const char *pattern, int patternLen, int maxStartShift, int limit,
            double *matchProb, const char *qualityString)
    {
        return compute(text, pattern, patternLen, maxStartShift, limit, NULL, matchProb, qualityString, NULL, 0, false);
    }

    // Compute distance and MAPQ but not edit string, and return the final shift
    int compute(const char *text, const char *pattern, int patternLen, int maxStartShift, int limit,
            double *matchProb, const char *qualityString, int *endShift)
    {
        return compute(text, pattern, patternLen, maxStartShift, limit, endShift, matchProb, qualityString, NULL, 0, false);
    }

    // Compute distance and edit string but not MAPQ
    int compute(const char *text, const char *pattern, int patternLen, int maxStartShift, int limit,
                char *editStrBuf, int editStrBufLen, bool useM = false)
    {
        return compute(text, pattern, patternLen, maxStartShift, limit, NULL, NULL, NULL, editStrBuf, editStrBufLen, useM);
    }

    // Use BigAlloc when allocating this object in case our arrays are big
    void *operator new(size_t size) {return BigAlloc(size);}
    void operator delete(void *ptr) {BigDealloc(ptr);}

private:
    /**
     * Fill in L[destD][destS][destG] and prev[destD][destS][destG] with the largest of two
     * choices of states we could've come from. For each state, we get a previous d, previous s,
     * and previous gapStatus, as well as a delta to add to the length stored in its L entry.
     * This method saves some repetition that would otherwise happen in compute().
     */
    inline void fillBest(int destD, int destS, int destG,
                         int d1, int s1, int g1, int delta1,
                         int d2, int s2, int g2, int delta2)
    {
        int val1 = L[d1][s1][g1] + delta1;
        int val2 = L[d2][s2][g2] + delta2;
        if (val1 >= val2) {
            L[destD][destS][destG] = val1;
            prevGapStatus[destD][destS][destG] = g1;
        } else {
            L[destD][destS][destG] = val2;
            prevGapStatus[destD][destS][destG] = g2;
        }
    }

    /**
     * Version of fillBest() for three possible parent states.
     */
    inline void fillBest(int destD, int destS, int destG,
                         int d1, int s1, int g1, int delta1,
                         int d2, int s2, int g2, int delta2,
                         int d3, int s3, int g3, int delta3)
    {
        int val1 = L[d1][s1][g1] + delta1;
        int val2 = L[d2][s2][g2] + delta2;
        int val3 = L[d3][s3][g3] + delta3;
        if (val1 >= val2) {
            if (val1 >= val3) {
                L[destD][destS][destG] = val1;
                prevGapStatus[destD][destS][destG] = g1;
            } else {
                L[destD][destS][destG] = val3;
                prevGapStatus[destD][destS][destG] = g3;
            }
        } else {
            if (val2 >= val3) {
                L[destD][destS][destG] = val2;
                prevGapStatus[destD][destS][destG] = g2;
            } else {
                L[destD][destS][destG] = val3;
                prevGapStatus[destD][destS][destG] = g3;
            }
        }
    }

    /**
     * Find the length of the longest common prefix of two strings; compares 8 bytes at a time
     */
    inline int longestPrefix(const char *text, const char *pattern, const char *patternEnd) {
        const char* p = pattern;
        const char* t = text;
        if (*p != *t) {
            return 0;
        }

        while (p < patternEnd) {
            _uint64 x = *((_uint64*) p) ^ *((_uint64*) t);
            if (x) {
                unsigned long zeroes;
                CountTrailingZeroes(x, zeroes);
                zeroes >>= 3;
                return (int)(__min((p - pattern) + zeroes, patternEnd - pattern));
            }
            p += 8;
            t += 8;
        }
        return (int)(patternEnd - pattern);
    }

    /**
     * Fill in the CIGAR string for a successfully matched pair of strings, given the distance and
     * shift we ended up at when we matched the end of the pattern. Returns true if successful.
     *
     * If useM is set, this will report both 'X' and '=' characters as 'M' for older SAM versions.
     */
    bool writeEditString(int endDistance, int endShift, char *editStringBuf, int editStringBufLen, bool useM) {
        if (editStringBuf == NULL) {
            return true;
        }
        if (editStringBufLen == 0) {
            return false;
        }

        // First, we'll backtrack through the L array and find out the actions we used to get to
        // endDistance / endShift. We'll fill these in reverse order using the prev* arrays.
        int numActions = 0;
        char actions[2*MAX_DISTANCE + 2];
        int matchCounts[2*MAX_DISTANCE + 2];

        int d = endDistance;
        int s = SHIFT_OFF + endShift;
        int g = NO_GAP;

        while (g != -1) {
            // Figure out which action we took to get here, and from which state
            int prevD, prevS, prevG;
            char action;
            getPreviousState(d, s, g, &prevD, &prevS, &prevG, &action);
            // Figure out how this affects our CIGAR string
            int matchCount = L[d][s][g] - (prevG == -1 ? 0 : L[prevD][prevS][prevG]); // Exact matches after this
            if (action == 'X' || action == 'I') {
                matchCount -= 1;
            }
            if (matchCount > 0) {
                // Add the '=' action for the exact matches after this edit
                actions[numActions] = (useM ? 'M' : '=');
                matchCounts[numActions] = matchCount;
                numActions++;
            }
            if (action != '=') {
                char stringEntry = (useM && action == 'X') ? 'M' : action;
                actions[numActions] = stringEntry;
                matchCounts[numActions] = 1;
                numActions++;
            }
            // Move to the previous state
            d = prevD; s = prevS; g = prevG;
        }

        // Now let's write the CIGAR string
        char *out = editStringBuf;
        int spaceLeft = editStringBufLen;
        for (int p = numActions - 1; p >= 0; p--) {
            int count = matchCounts[p];
            while (p >= 1 && actions[p-1] == actions[p]) {
                // Multiple actions of the same type occur near each other; coalesce them
                count += matchCounts[p-1];
                p--;
            }
            int written = snprintf(out, spaceLeft, "%d%c", count, actions[p]);
            if (written < 0 || written == spaceLeft) {
                return false;    // snprintf failed to write enough bytes
            } else {
                out += written;
                spaceLeft -= written;
            }
        }
        if (numActions == 0) {
            // The CIGAR is an empty string, but add a \0 to mark the end; note that we checked
            // that editStringBufLen > 0 at the top
            *out = '\0';
        }
        
        return true;
    }

    /**
     * Fill in the match probability for a successfully matched pair of strings, given the distance
     * and shift we ended up at when we matched the end of the pattern. Does nothing if the output
     * pointer is null.
     */
    void computeMapProbability(int endDistance, int endShift, double *mapProbability, const char *qualityString) {
        if (mapProbability == NULL) {
            return;
        }

        *mapProbability = 1.0;

        // Backtrack through the L and prev arrays, multiplying the map probability as we go along
        int d = endDistance;
        int s = SHIFT_OFF + endShift;
        int g = NO_GAP;
        char oldAction = '=';   // Used to tell whether we're expanding or starting a gap

        while (g != -1) {
            // Figure out which action we took to get here, and from which state
            int prevD, prevS, prevG;
            char action;
            getPreviousState(d, s, g, &prevD, &prevS, &prevG, &action);
            // Figure out how this affects the map probability
            int matchCount = L[d][s][g] - (prevG == -1 ? 0 : L[prevD][prevS][prevG]); // Exact matches after this
            if (action == 'X' || action == 'I') {
                matchCount -= 1;
            }
            if (matchCount > 0) {
                // Multiply in the probability of this many matches after the edit
                *mapProbability *= matchProbability[matchCount];
                // Reset oldAction since there's no gap after the edit
                oldAction = '=';
            }
            switch (action) {
                case 'I':
                    *mapProbability *= (oldAction == 'I') ? indelExtendProbability : indelOpenProbability;
                    break;
                case 'D':
                    *mapProbability *= (oldAction == 'D') ? indelExtendProbability : indelOpenProbability;
                    break;
                case 'X':
                    *mapProbability *= mismatchProbability[qualityString[L[d][s][g]-1]];
                    break;
                default:
                    // We can also get the '=' action in some cases, but we its probability above
                    break;
            }
            // Move to the previous state
            d = prevD; s = prevS; g = prevG;
            oldAction = action;
        }
    }

    // Find the distance / shift / gapStatus state that got us to a current state, as well as the
    // edit action we took to get there. Uses the prevGapStatus to figure this out.
    inline void getPreviousState(int curD, int curS, int curG, int *prevD, int *prevS, int *prevG, char *action) {
        *prevG = prevGapStatus[curD][curS][curG];
        switch (curG) {
            case NO_GAP:
                if (*prevG == NO_GAP) {
                    *action = 'X';
                    *prevD = curD - substitutionPenalty;
                    *prevS = curS;
                } else { // Came from TEXT_GAP or PATTERN_GAP at the same position
                    *action = '=';
                    *prevD = curD;
                    *prevS = curS;
                }
                break;
            case TEXT_GAP:
                if (*prevG == TEXT_GAP) {
                    *action = 'D';
                    *prevD = curD - 1;
                    *prevS = curS - 1;
                } else {  // Came from NO_GAP
                    *action = 'D';
                    *prevD = curD - gapOpenPenalty;
                    *prevS = curS - 1;
                }
                break;
            case PATTERN_GAP:
                if (*prevG == PATTERN_GAP) {
                    *action = 'I';
                    *prevD = curD - 1;
                    *prevS = curS + 1;
                } else {  // Came from NO_GAP
                    *action = 'I';
                    *prevD = curD - gapOpenPenalty;
                    *prevS = curS + 1;
                }
                break;
        }
    }
};
