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
template<bool COMPUTE_EDIT_STRING=false, bool COMPUTE_MAPQ=false, int MAX_DISTANCE=50, int MAX_SHIFT=MAX_DISTANCE>
class BoundedStringDistance
{
private:
    enum GapStatus { NO_GAP, TEXT_GAP, PATTERN_GAP };
    static const int SHIFT_OFF = MAX_SHIFT + 1;   // For indexing into arrays by shift
    static const bool COMPUTE_PREV = COMPUTE_EDIT_STRING || COMPUTE_MAPQ;  // Whether we need to backtrack

    int substitutionPenalty;
    int gapOpenPenalty;       // We assume that the gap extend penalty is 1

    // L[distance][shift][gapStatus] is the largest value x such that the distance
    // between pattern[0..x] and text[0..x+shift] is "distance" and the gap status
    // at the end is "gapStatus".
    int L[MAX_DISTANCE+1][2*MAX_SHIFT+3][3];

    // A state in our search through the strings, i.e. a set of coordinates in the L array,
    // as well as the edit action (e.g. 'X', 'D') we took to move away from them
    struct State { 
        int distance;
        int shift;
        int gapStatus;
        char action;

        State(): distance(-1), shift(-1), gapStatus(-1), action('=') {}

        State(int distance_, int shift_, int gapStatus_, char action_)
            : distance(distance_), shift(shift_), gapStatus(gapStatus_), action(action_) {}
    };

    // Previous state we were in before getting to each particular state
    State prev[MAX_DISTANCE+1][2*MAX_SHIFT][3];

public:
    BoundedStringDistance(int substitutionPenalty_, int gapOpenPenalty_)
            : substitutionPenalty(substitutionPenalty_), gapOpenPenalty(gapOpenPenalty_)
    {
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
                char *editStringBuf, int editStringBufLen, bool useM = false)
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
        for (int d = 0; d <= __min(limit, gapOpenPenalty - 1); d += substitutionPenalty) {
            pos += longestPrefix(text+pos, pattern+pos, patternEnd);
            L[d][SHIFT_OFF][NO_GAP] = pos;
            if (COMPUTE_EDIT_STRING && d != 0) {
                prev[d][SHIFT_OFF][NO_GAP] = State(d-substitutionPenalty, SHIFT_OFF, NO_GAP, 'X');
            }
            if (pos == patternLen) {
                if (!writeEditString(d, 0, editStringBuf, editStringBufLen, useM)) {
                    return -2;
                }
                return d;
            } else {
                pos += 1; // Skip the mismatch so we can go to the next d
            }
        }

        // For distances >= gapOpenPenalty, also allow creating gaps in the text or pattern
        for (int d = gapOpenPenalty; d <= limit; d++) {
            int maxShift = __min(d - gapOpenPenalty + 1, MAX_SHIFT);
            for (int s = SHIFT_OFF - maxShift; s <= SHIFT_OFF + maxShift; s++) {
                // See whether we can start / extend the "text gap" case
                fillBest(d, s, TEXT_GAP,
                         d-gapOpenPenalty, s-1, NO_GAP, 0, 'D',  // Open a gap from NO_GAP
                         d-1, s-1, TEXT_GAP, 0, 'D');            // Continue a TEXT_GAP
                // Note that continuing from a PATTERN_GAP is covered because the NO_GAP
                // case at d-gapOpenPenalty considered closing the PATTERN_GAP there

                // See whether we can start / extend the "pattern gap" case
                fillBest(d, s, PATTERN_GAP,
                         d-gapOpenPenalty, s+1, NO_GAP, 1, 'I',  // Open a new gap
                         d-1, s+1, PATTERN_GAP, 1, 'I');         // Continue a PATTERN_GAP
                // Note that continuing from a TEXT_GAP is covered because the NO_GAP
                // case at d-gapOpenPenalty considered closing the TEXT_GAP there

                // Find the best way to continue the "neither in gap" case; this could be either
                // adding a substitution from the previous "neither in gap" state, or closing a
                // text gap, or a pattern gap. Note that since we compute the best ways to end up
                // in TEXT_GAP and PATTERN_GAP above, this automatically includes using those as
                // starting points and just closing the gap that got us to this distance.
                fillBest(d, s, NO_GAP,
                         d-substitutionPenalty, s, NO_GAP, 1, 'X',  // Add subsitution after NO_GAP
                         d, s, TEXT_GAP, 0, '=',                    // Close a TEXT_GAP and go from there
                         d, s, PATTERN_GAP, 0, '=');                // Close a PATTERN_GAP
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
                    return d;
                }
            }
        }

        return -1;
    }

    int compute(const char *text, const char *pattern, int patternLen, int limit) {
        return compute(text, pattern, patternLen, limit, NULL, 0);
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
    inline void fillBest(
            int destD, int destS, int destG,
            int d1, int s1, int g1, int delta1, char action1,
            int d2, int s2, int g2, int delta2, char action2)
    {
        int val1 = (d1 >= 0) ? L[d1][s1][g1] + delta1 : -2;
        int val2 = (d2 >= 0) ? L[d2][s2][g2] + delta2 : -2;
        if (val1 >= val2) {
            L[destD][destS][destG] = (val1 >= 0) ? val1 : -2;
            if (COMPUTE_PREV) {
                prev[destD][destS][destG] = State(d1, s1, g1, action1);
            }
        } else {
            L[destD][destS][destG] = (val2 >= 0) ? val2 : -2;
            if (COMPUTE_PREV) {
                prev[destD][destS][destG] = State(d2, s2, g2, action2);
            }
        }
    }

    /**
     * Version of fillBest() for three possible parent states.
     */
    inline void fillBest(
            int destD, int destS, int destG,
            int d1, int s1, int g1, int delta1, char action1,
            int d2, int s2, int g2, int delta2, char action2,
            int d3, int s3, int g3, int delta3, char action3)
    {
        int val1 = (d1 >= 0) ? L[d1][s1][g1] + delta1 : -2;
        int val2 = (d2 >= 0) ? L[d2][s2][g2] + delta2 : -2;
        int val3 = (d3 >= 0) ? L[d3][s3][g3] + delta3 : -2;
        if (val1 >= val2) {
            if (val1 >= val3) {
                L[destD][destS][destG] = (val1 >= 0) ? val1 : -2;
                if (COMPUTE_PREV) {
                    prev[destD][destS][destG] = State(d1, s1, g1, action1);
                }
            } else {
                L[destD][destS][destG] = (val3 >= 0) ? val3 : -2;
                if (COMPUTE_PREV) {
                    prev[destD][destS][destG] = State(d3, s3, g3, action3);
                }
            }
        } else {
            if (val2 >= val3) {
                L[destD][destS][destG] = (val2 >= 0) ? val2 : -2;
                if (COMPUTE_PREV) {
                    prev[destD][destS][destG] = State(d2, s2, g2, action2);
                }
            } else {
                L[destD][destS][destG] = (val3 >= 0) ? val3 : -2;
                if (COMPUTE_PREV) {
                    prev[destD][destS][destG] = State(d3, s3, g3, action3);
                }
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
        if (!COMPUTE_EDIT_STRING || editStringBuf == NULL) {
            return true;
        }

        // First, we'll backtrack through the L array and find out the actions we used to get to
        // endDistance / endShift. We'll fill these in reverse order using the prev* arrays.
        int numActions = 0;
        char actions[2*MAX_DISTANCE + 2];
        int matchCounts[2*MAX_DISTANCE + 2];

        int d = endDistance;
        int s = SHIFT_OFF + endShift;
        int g = NO_GAP;

        while (d != -1) {
            int prevD = prev[d][s][g].distance;
            int prevS = prev[d][s][g].shift;
            int prevG = prev[d][s][g].gapStatus;
            char action = prev[d][s][g].action;
            int matchCount = L[d][s][g] - (prevD == -1 ? 0 : L[prevD][prevS][prevG]); // Exact matches after this
            if (action == 'X' || action == 'I') {
                matchCount -= 1;
            }
            if (matchCount > 0) {
                // Add the '=' action for the exact matches after this edit
                actions[numActions] = (useM ? 'M' : '=');
                matchCounts[numActions] = matchCount;
                numActions++;
            }
            if (prev[d][s][g].action != '=') {
                char stringEntry = (useM && action == 'X') ? 'M' : action;
                actions[numActions] = stringEntry;
                matchCounts[numActions] = 1;
                numActions++;
            }
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
        
        return true;
    }
};
