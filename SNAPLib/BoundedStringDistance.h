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
 * - COMPUTE_EDIT_STRING: whether to support computing edit strings (in CIGAR format).
 * - MAX_DISTANCE: maximum distance score supported.
 * - MAX_SHIFT: maximum shift supported. This can be less than MAX_DISTANCE if one
 *   is not expecting huge indels, in order to save computation time.
 *
 * TODO: Allow making the L array "cyclic" through another template parameter (i.e.
 * so that only the last few rows will be kept) to save cache space for big arrays.
 */
template<bool COMPUTE_EDIT_STRING=false, int MAX_DISTANCE=50, int MAX_SHIFT=MAX_DISTANCE>
class BoundedStringDistance
{
public:
    BoundedStringDistance(int gapOpenPenalty_) : gapOpenPenalty(gapOpenPenalty_) {
        // Clear out the L array.
        for (int d = 0; d < MAX_DISTANCE + 1; d++) {
            for (int s = 0; s < 2 * MAX_SHIFT + 3; s++) {
                for (int g = 0; g < 3; g++) {
                    L[d][s][g] = -2;
                }
            }
        }
    }

    /**
     * Compute the string distance between text[0...] and pattern[0..patternLen], if
     * it is less than limit, or return -1 if it is greater. Optionally, write the
     * edit string (in CIGAR format) into editStringBuf if it's not null and if
     * COMPUTE_EDIT_STRING was set.
     *
     * Assumptions:
     * - text is indexable up to patternLen + MAX_SHIFT - 1 (so it's a pointer
     *   into a longer string that we don't need to bounds-check).
     * - pattern is indexable up to 7 characters past its length.
     * - editStringBuf is long enough to hold the edit string, if set. We check that
     *   it's at least 2 * patternLen in case we must make one edit per character.
     */
    int compute(const char *text, const char *pattern, int patternLen, int limit,
                char *editStringBuf, int editStringBufLen)
    {
        if (limit < 0 || limit > MAX_DISTANCE) {
            fprintf(stderr, "Invalid distance limit: %d\n", limit);
            exit(1);
        }

        if (editStringBuf != NULL && !COMPUTE_EDIT_STRING) {
            fprintf(stderr, "Non-null editStringBuf given but COMPUTE_EDIT_STRING is false\n");
            exit(1);
        }

        const char *patternEnd = pattern + patternLen;

        // For distances less than gapOpenPenalty, just search for mismatches with shift = 0
        int pos = 0;
        for (int d = 0; d <= __min(limit, gapOpenPenalty - 1); d++) {
            pos += longestPrefix(text+pos, pattern+pos, patternEnd);
            L[d][SHIFT_OFF][NO_GAP] = pos;
            if (pos == patternLen) {
                return d;
            }
            pos += 1; // Skip the mismatch so we can go to the next d
        }

        // For distances >= gapOpenPenalty, also allow creating gaps in the text or pattern
        for (int d = gapOpenPenalty; d <= limit; d++) {
            int maxShift = __min(d - gapOpenPenalty + 1, MAX_SHIFT);
            for (int s = -maxShift; s <= maxShift; s++) {
                // See whether we can start / extend the "text gap" case
                int bestTextGap = __max(L[d-gapOpenPenalty][SHIFT_OFF+s-1][NO_GAP],
                                        L[d-1][SHIFT_OFF+s-1][TEXT_GAP]);
                L[d][SHIFT_OFF+s][TEXT_GAP] = bestTextGap;

                // See whether we can start / extend the "pattern gap" case
                int bestPatternGap = __max(L[d-gapOpenPenalty][SHIFT_OFF+s+1][NO_GAP] + 1,
                                           L[d-1][SHIFT_OFF+s+1][PATTERN_GAP] + 1);
                if (bestPatternGap < 0)
                    bestPatternGap = -2; // Just in case we did a -2 + 1
                L[d][SHIFT_OFF+s][PATTERN_GAP] = bestPatternGap;

                // Find the best way to continue the "neither in gap" case
                int bestNoGap = L[d-1][SHIFT_OFF+s][NO_GAP] + 1;
                int closeTextGap = L[d][SHIFT_OFF+s][TEXT_GAP];
                if (closeTextGap > bestNoGap)
                    bestNoGap = closeTextGap;
                int closePatternGap = L[d][SHIFT_OFF+s][PATTERN_GAP];
                if (closePatternGap > bestNoGap)
                    bestNoGap = closePatternGap;
                bestNoGap += longestPrefix(text+bestNoGap+s, pattern+bestNoGap, patternEnd);
                L[d][SHIFT_OFF+s][NO_GAP] = bestNoGap;
                if (bestNoGap == patternLen) {
                    return d;
                }
            }
        }

        return -1;
    }

    int compute(const char *text, const char *pattern, int patternLen, int limit) {
        return compute(text, pattern, patternLen, limit, NULL, 0);
    }

    void *operator new(size_t size) {return BigAlloc(size);}
    void operator delete(void *ptr) {BigDealloc(ptr);}

private:
    enum GapStatus { NO_GAP, TEXT_GAP, PATTERN_GAP };
    static const int SHIFT_OFF = MAX_SHIFT + 1;

    int gapOpenPenalty; // We assume that the gap extend penalty is 1.

    // L[distance][shift][gapStatus] is the largest value x such that the score
    // between pattern[0..x] and text[0..x+shift] is distance and the gap status
    // at the end is gapStatus.
    int L[MAX_DISTANCE+1][2*MAX_SHIFT+3][3];

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
                return __min((p - pattern) + zeroes, patternEnd - pattern);
            }
            p += 8;
            t += 8;
        }
        return patternEnd - pattern;
    }
};
