#pragma once
#include "Compat.h"
#include "FixedSizeMap.h"
#include "BigAlloc.h"

const int MAX_K = 31;

//
// These are global so there are only one for both senses of the template
//

extern double *lv_indelProbabilities;  // Maps indels by length to probability of occurance.
extern double *lv_phredToProbability;  // Maps ASCII phred character to probability of error, including 
extern double *lv_perfectMatchProbability; // Probability that a read of this length has no mutations


// Computes the edit distance between two strings without returning the edits themselves.
// Set TEXT_DIRECTION to -1 to run backwards through the text.
template<int TEXT_DIRECTION = 1> class LandauVishkin {
public:
    LandauVishkin(int cacheSize = 0)
{
    if (TEXT_DIRECTION != 1 && TEXT_DIRECTION != -1) {
        printf("You can't possibly be serious.\n");
        exit(1);
    }

    for (int i = 0; i < MAX_K+1; i++) {
        for (int j = 0; j < 2*MAX_K+1; j++) {
        L[i][j] = -2;
        }
    }
    if (cacheSize > 0) {
        cache = new FixedSizeMap<_uint64, LVResult>();
        cache->reserve(cacheSize);
    } else {
        cache = NULL;
    }
}

    static size_t getBigAllocatorReservation() {return sizeof(LandauVishkin<TEXT_DIRECTION>);} // maybe we should worry about allocating the cache with a BigAllocator, but not for now.

    ~LandauVishkin()
{
    if (cache != NULL) {
        delete cache;
    }
}

    // Compute the edit distance between two strings, if it is <= k, or return -1 otherwise.
    // For LandauVishkin instances with a cache, the cacheKey should be a unique identifier for
    // the text and pattern combination (e.g. (readID << 33) | direction << 32 | genomeLocation).
    int computeEditDistance(
            const char* text,
            int textLen, 
            const char* pattern,
            const char *qualityString,
            int patternLen,
            int k,
            double *matchProbability,
            _uint64 cacheKey = 0)
{
    _ASSERT(k < MAX_K);
    k = __min(MAX_K - 1, k); // enforce limit even in non-debug builds
    if (NULL == text) {
        // This happens when we're trying to read past the end of the genome.
        if (NULL != matchProbability) {
            *matchProbability = 0.0;
        }
        return -1;
    }
    if (cache != NULL && cacheKey != 0) {
        LVResult old = cache->get(cacheKey);
        if (old.isValid() && (old.result != -1 || old.k >= k)) {
            if (NULL != matchProbability) {
                *matchProbability = old.matchProbability;
            }
            return old.result;
        }
    }
    if (NULL != matchProbability) {
        //
        // Start with perfect match probability and work our way down.
        //
        *matchProbability = 1.0;    
    }
    if (TEXT_DIRECTION == -1) {
        text--; // so now it points at the "first" character of t, not after it.
    }
    const char* p = pattern;
    const char* t = text;
    int end = __min(patternLen, textLen);
    const char* pend = pattern + end;
    while (p < pend) {
        _uint64 x;
        if (TEXT_DIRECTION == 1) {
            x = *((_uint64*) p) ^ *((_uint64*) t);
        } else {
            _uint64 T = *(_uint64 *)(t - 7);
            _uint64 tSwap = ByteSwapUI64(T);
            x = *((_uint64*) p) ^ tSwap;
        }

        if (x) {
            unsigned long zeroes;
            CountTrailingZeroes(x, zeroes);
            zeroes >>= 3;
            L[0][MAX_K] = __min((int)(p - pattern) + (int)zeroes, end);
            goto done1;
        }
        p += 8;
        t += 8 * TEXT_DIRECTION;
    }
    L[0][MAX_K] = end;
    done1:
    if (L[0][MAX_K] == end) {
        int result = (patternLen > end ? patternLen - end : 0); // Could need some deletions at the end
        if (NULL != matchProbability) {
            *matchProbability = lv_perfectMatchProbability[patternLen];    // Becuase the chance of a perfect match is < 1
            if (cache != NULL && cacheKey != 0) {
                cache->put(cacheKey, LVResult(k, result, lv_indelProbabilities[result]));
            }
        }
        return result;
    }

    for (int e = 1; e <= k; e++) {
        // Search d's in the order 0, 1, -1, 2, -2, etc to find an alignment with as few indels as possible.
        for (int d = 0; d != e+1; d = (d > 0 ? -d : -d+1)) {
            int best = L[e-1][MAX_K+d] + 1; // up
            A[e][MAX_K+d] = 'X';
            int left = L[e-1][MAX_K+d-1];
            if (left > best) {
                best = left;
                A[e][MAX_K+d] = 'D';
            }
            int right = L[e-1][MAX_K+d+1] + 1;
            if (right > best) {
                best = right;
                A[e][MAX_K+d] = 'I';
            }

            const char* p = pattern + best;
            const char* t = (text + d * TEXT_DIRECTION) + best * TEXT_DIRECTION;
            if (*p == *t) {
                int end = __min(patternLen, textLen - d);
                const char* pend = pattern + end;

                while (true) {
                    _uint64 x;
                    if (TEXT_DIRECTION == 1) {
                        x = *((_uint64*) p) ^ *((_uint64*) t);
                    } else {
                        _uint64 T = *(_uint64 *)(t - 7);
                        _uint64 tSwap = ByteSwapUI64(T);
                        x = *((_uint64*) p) ^ tSwap;
                    }
                    if (x) {
                        unsigned long zeroes;
                        CountTrailingZeroes(x, zeroes);
                        zeroes >>= 3;
                        best = __min((int)(p - pattern) + (int)zeroes, end);
                        break;
                    }
                    p += 8;
                    if (p >= pend) {
                        best = end;
                        break;
                    }
                    t += 8 * TEXT_DIRECTION;
                }
            }

            if (best == patternLen) {
                if (NULL != matchProbability) {
                    _ASSERT(*matchProbability == 1.0);
                    //
                    // We're done.  Compute the match probability.
                    //
                    int straightMismatches = 0;
                    for (int i = 0; i < end && straightMismatches <= e; i++) { // do this 8 at a time, like in the other loops!
                        if (pattern[i] != text[i * TEXT_DIRECTION]) {
                            straightMismatches++;
                            *matchProbability *= lv_phredToProbability[qualityString[i]];
                        }
                    }

                    if (straightMismatches != e) {
                        //
                        // There are indels.  
                        // Trace backward to build up the CIGAR string.  We do this by filling in the backtraceAction,
                        // backtraceMatched and backtraceD arrays, then going through them in the forward direction to
                        // figure out our string.
                        *matchProbability = 1.0;
                        int curD = d;
                        for (int curE = e; curE >= 1; curE--) {
                            backtraceAction[curE] = A[curE][MAX_K+curD];
                            if (backtraceAction[curE] == 'I') {
                                backtraceD[curE] = curD + 1;
                                backtraceMatched[curE] = L[curE][MAX_K+curD] - L[curE-1][MAX_K+curD+1] - 1;
                            } else if (backtraceAction[curE] == 'D') {
                                backtraceD[curE] = curD - 1;
                                backtraceMatched[curE] = L[curE][MAX_K+curD] - L[curE-1][MAX_K+curD-1];
                            } else { // backtraceAction[curE] == 'X'
                                backtraceD[curE] = curD;
                                backtraceMatched[curE] = L[curE][MAX_K+curD] - L[curE-1][MAX_K+curD] - 1;
                            }
                            curD = backtraceD[curE];
        #ifdef TRACE_LV
                            printf("%d %d: %d %c %d %d\n", curE, curD, L[curE][MAX_K+curD], 
                                backtraceAction[curE], backtraceD[curE], backtraceMatched[curE]);
        #endif
                        }

                        int curE = 1;
                        int offset = L[0][MAX_K+0];
                        while (curE <= e) {
                            // First write the action, possibly with a repeat if it occurred multiple times with no exact matches
                            char action = backtraceAction[curE];
                            int actionCount = 1;
                            while (curE+1 <= e && backtraceMatched[curE] == 0 && backtraceAction[curE+1] == action) {
                                actionCount++;
                                curE++;
                            }
                            if (action == 'I') {
                                *matchProbability *= lv_indelProbabilities[actionCount];
                                offset += actionCount; 
                            } else if (action =='D') {
                                *matchProbability *= lv_indelProbabilities[actionCount];
                                offset -= actionCount;
                            }  else {
                                _ASSERT(action == 'X');
                                for (int i = 0; i < actionCount; i++) {
                                    *matchProbability *= lv_phredToProbability[qualityString[/*BUGBUG - think about what to do here*/__min(patternLen-1,__max(offset,0))]];
                                    offset++;
                                }
                            }

                            offset += backtraceMatched[curE];   // Skip over the matching bases.
                            curE++;
                        }
                    } // if straightMismatches != e (i.e., the indel case)
                    if (cache != NULL && cacheKey != 0) {
                        cache->put(cacheKey, LVResult(k, e, *matchProbability));
                    } 
                    *matchProbability *= lv_perfectMatchProbability[patternLen-e]; // Accounting for the < 1.0 chance of no changes for matching bases
                } else {
                    //
                    // Not tracking match probability.
                    //
                    if (cache != NULL && cacheKey != 0) {
                        cache->put(cacheKey, LVResult(k, e, -1.0));
                    }
                }
                return e;
            }

            L[e][MAX_K+d] = best;
        }
    }

    if (cache != NULL && cacheKey != 0) {
        cache->put(cacheKey, LVResult(k, -1, 0.0));
    }
    return -1;
}


    // Version that does not requre match probability and quality string
    inline int computeEditDistance(
            const char* text,
            int textLen,
            const char* pattern,
            int patternLen,
            int k)
    {
        return computeEditDistance(text, textLen, pattern, NULL, patternLen, k, NULL);
    }

    // Clear the cache of distances computed.
    void clearCache()
    {
        if (cache != NULL) {
            cache->clear();
        }
    }

    void *operator new(size_t size) {return BigAlloc(size);}
    void operator delete(void *ptr) {BigDealloc(ptr);}

    void *operator new(size_t size, BigAllocator *allocator) {_ASSERT(size == sizeof(LandauVishkin<TEXT_DIRECTION>)); return allocator->allocate(size);}
    void operator delete(void *ptr, BigAllocator *allocator) {/*Do nothing.  The memory is freed when the allocator is deleted.*/}
 
private:
    // TODO: For long reads, we should include a version that only has L be 2 x (2*MAX_K+1) cells
    int L[MAX_K+1][2 * MAX_K + 1];
    
    // Action we did to get to each position: 'D' = deletion, 'I' = insertion, 'X' = substitution.  This is needed to compute match probability.
    char A[MAX_K+1][2 * MAX_K + 1];

    // Arrays for backtracing the actions required to match two strings
    char backtraceAction[MAX_K+1];
    int  backtraceMatched[MAX_K+1];
    int  backtraceD[MAX_K+1];

    struct LVResult {
        int k;
        int result;
        double matchProbability;

        LVResult() { k = -1; result = -1; }

        LVResult(int k_, int result_, double matchProbability_) { k = k_; result = result_; matchProbability = matchProbability_;}

        inline bool isValid() { return k != -1; }
    };

    FixedSizeMap<_uint64, LVResult> *cache;
};

void setLVProbabilities(double *i_indelProbabilities, double *i_phredToProbability, double mutationProbability);
void initializeLVProbabilitiesToPhredPlus33();


// Computes the edit distance between two strings and returns a CIGAR string for the edits.

enum CigarFormat
{
    COMPACT_CIGAR_STRING = 0,
    EXPANDED_CIGAR_STRING = 1,
    COMPACT_CIGAR_BINARY = 2,
};

class LandauVishkinWithCigar {
public:
    LandauVishkinWithCigar();

    // Compute the edit distance between two strings and write the CIGAR string in cigarBuf.
    // Returns -1 if the edit distance exceeds k or -2 if we run out of space in cigarBuf.
    int computeEditDistance(const char* text, int textLen, const char* pattern, int patternLen, int k,
                            char* cigarBuf, int cigarBufLen, bool useM, CigarFormat format = COMPACT_CIGAR_STRING);
private:
    int L[MAX_K+1][2 * MAX_K + 1];
    
    // Action we did to get to each position: 'D' = deletion, 'I' = insertion, 'X' = substitution.
    char A[MAX_K+1][2 * MAX_K + 1];

    // Arrays for backtracing the actions required to match two strings
    char backtraceAction[MAX_K+1];
    int backtraceMatched[MAX_K+1];
    int backtraceD[MAX_K+1];
};
