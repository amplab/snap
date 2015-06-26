#pragma once
#include "Compat.h"
#include "FixedSizeMap.h"
#include "BigAlloc.h"
#include "exit.h"
#include "Genome.h"

const int MAX_K = 63;

//
// These are global so there are only one for both senses of the template
//

extern double *lv_indelProbabilities;  // Maps indels by length to probability of occurance.
extern double *lv_phredToProbability;  // Maps ASCII phred character to probability of error, including 
extern double *lv_perfectMatchProbability; // Probability that a read of this length has no mutations

struct LVResult {
    short k;
    short result;
    short netIndel;
    double matchProbability;

    LVResult() { k = -1; result = -1; netIndel = 0;}

    LVResult(short k_, short result_, short netIndel_, double matchProbability_) { 
        k = k_; 
        result = result_; 
        netIndel = netIndel_; 
        matchProbability = matchProbability_;
    }

    inline bool isValid() { return k != -1; }
};



static inline void memsetint(int* p, int value, int count)
{
// this is required to get around a GCC optimization bug
#ifndef _MSC_VER
  volatile
#endif
  int * q = p;
  for (int i = 0; i < count; i++) {
    q[i] = value;
  }
}

// Computes the edit distance between two strings without returning the edits themselves.
// Set TEXT_DIRECTION to -1 to run backwards through the text.
template<int TEXT_DIRECTION = 1> class LandauVishkin {

//
// Macros to make arrays with negative indices seem "natural" in the code.
//
#define L(e,d)			L_zero			[(e) * (2 * MAX_K + 1) + (d)]
#define A(e,d)			A_zero			[(e) * (2 * MAX_K + 1) + (d)]

public:
    LandauVishkin()
{
    if (TEXT_DIRECTION != 1 && TEXT_DIRECTION != -1) {
        fprintf(stderr, "You can't possibly be serious.\n");
        soft_exit(1);
    }

    memsetint(L_space, -2, (MAX_K + 1) * (2 * MAX_K + 1));

    L_zero = L_space + MAX_K;   // The address of L(0,0)
    A_zero = A_space + MAX_K;   // The address of A(0,0)

    //
    // Initialize dTable, which is used to avoid a branch misprediction in our inner loop.
    // The d values are 0, -1, 1, -2, 2, etc.
    //
    for (int i = 0, d = 0; i < 2 * (MAX_K + 1) + 1; i++, d = (d > 0 ? -d : -d+1)) {
        dTable[i] = d;
    }
}
/*
    void pushBackCacheStats()
    {
        if (NULL != cache) {
            cache->pushBackCacheStats();
        }
    }
*/
    static size_t getBigAllocatorReservation() {return sizeof(LandauVishkin<TEXT_DIRECTION>);} // maybe we should worry about allocating the cache with a BigAllocator, but not for now.

    ~LandauVishkin()
{
}

    // Compute the edit distance between two strings, if it is <= k, or return -1 otherwise.

	//
	// The essential method is to build up the L array row by row.  L[e][d] is the farthest that you can get
	// through the pattern (read data) with e changes (single base substition, insert or delete) and a net indel of d.
	// Once you get to the end of the read, you've computed the best edit distance (e).  L[e][d] can be computed
	// by looking at L[e-1][d-1 .. d+1], depending on whether the next change is a deletion, insertion or 
	// substitution.
	//
	// Because d can be negative, the L array doesn't really use L[e][d].  Instead, it uses L[e][MAX_K+d], because MAX_K is
	// the largest possible edit distance, and hence d can never be less than -MAX_K, so MAX_K + d >= 0.  However, the L and A macros
    // conceal this internally.
	//
	// Also, because of the way the alignment algorithms work, sometimes SNAP wants to run the edit distance
	// backward.  This is built as a template with TEXT_DIRECTION either 1 for forward or -1 for backward, just to make
	// it extra confusing.
	//

    int computeEditDistance(
                const char* text,
                int textLen, 
                const char* pattern,
                const char *qualityString,
                int patternLen,
                int k,
                double *matchProbability,
                int *o_netIndel = NULL)   // the net of insertions and deletions in the alignment.  Negative for insertions, positive for deleteions (and 0 if there are non in net).  Filled in only if matchProbability is non-NULL
{
    int localNetIndel;
	int d;
    if (NULL == o_netIndel) {
        //
        // If the user doesn't want netIndel, just use a stack local to avoid
        // having to check it all the time.
        //
        o_netIndel = &localNetIndel;
    }
    _ASSERT(k < MAX_K);

    *o_netIndel = 0;

    k = __min(MAX_K - 1, k); // enforce limit even in non-debug builds
    if (NULL == text) {
        // This happens when we're trying to read past the end of the genome.
        if (NULL != matchProbability) {
            *matchProbability = 0.0;
        }
        return -1;
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

    L(0, 0) = countPerfectMatch(p, t, end);

    if (L(0, 0) == end) {
        int result = (patternLen > end ? patternLen - end : 0); // Could need some deletions at the end
        if (NULL != matchProbability) {
            *matchProbability = lv_perfectMatchProbability[patternLen];    // Becuase the chance of a perfect match is < 1
        }
        if (result > k) {
            //
            // The deletions at the end pushed us oevr the score limit.
            //
            return -1;
        }
        return result;
    }

	int lastBestD = MAX_K + 1;
	int e;

    for (e = 1; e <= k; e++) {
        // Search d's in the order 0, 1, -1, 2, -2, etc to find an alignment with as few indels as possible.
        // dTable is just precomputed d = (d > 0 ? -d : -d+1) to save the branch misprediction from (d > 0)
        int i =0;
        for (d = 0; d != e+1 ; i++, d = dTable[i]) {
            int best = L(e-1, d) + 1; // up
            A(e, d) = 'X';

            const char* p = pattern + best;
            const char* t = (text + d * TEXT_DIRECTION) + best * TEXT_DIRECTION;
            if (*p == *t && best >= 0) {
                int end = __min(patternLen, textLen - d);
                const char* pend = pattern + end;

                best += countPerfectMatch(p, t, (int)(end - (p - pattern)));
            }


            int left = L(e-1, d-1);
            p = pattern + left;
            t = (text + d * TEXT_DIRECTION) + left * TEXT_DIRECTION;
            if (*p == *t && left >= 0) {
                int end = __min(patternLen, textLen - d);
                const char* pend = pattern + end;

                left += countPerfectMatch(p, t, (int)(end - (p - pattern)));
            }

            if (left > best) {
                best = left;
                A(e, d) = 'D';
            }

            int right = L(e-1, d+1) + 1;
            p = pattern + right;
            t = (text + d * TEXT_DIRECTION) + right * TEXT_DIRECTION;
            if (*p == *t && right >= 0) {
                int end = __min(patternLen, textLen - d);
                const char* pend = pattern + end;

                right += countPerfectMatch(p, t, (int)(end - (p - pattern)));
            }

            if (right > best) {
                best = right;
                A(e, d) = 'I';
            }

			if (best == patternLen) {
				//
				// We're through on this iteration.
				//

				if ('X' == A(e, d)) {
					//
					// The last step wasn't an indel, so we're sure it's the right one.
					//
					lastBestD = d;
					goto got_answer;
				} else {
					//
					// We're done on this round, but maybe there's a better answer, so keep looking.
					//
					if (abs(d) < abs(lastBestD)) {
						lastBestD = d;
					}
				}
			} // if best==patternLen

            L(e, d) = best;
        } // for d

		if (MAX_K + 1 != lastBestD) {
			break;
		}
    } // for e

	if (MAX_K + 1 == lastBestD) {
		return -1;
	}

got_answer:

	_ASSERT(abs(lastBestD) < MAX_K + 1);

	if (NULL != matchProbability) {

		_ASSERT(*matchProbability == 1.0);
		//
		// We're done.  Compute the match probability.
		//

		//
		// Trace backward to build up the CIGAR string.  We do this by filling in the backtraceAction,
		// backtraceMatched and backtraceD arrays, then going through them in the forward direction to
		// figure out our string.
		int curD = lastBestD;
		for (int curE = e; curE >= 1; curE--) {
			backtraceAction[curE] = A(curE, curD);
			if (backtraceAction[curE] == 'I') {
				backtraceD[curE] = curD + 1;
				backtraceMatched[curE] = L(curE, curD) - L(curE - 1, curD + 1) - 1;
			} else if (backtraceAction[curE] == 'D') {
				backtraceD[curE] = curD - 1;
				backtraceMatched[curE] = L(curE, curD) - L(curE - 1, curD - 1);
			} else { // backtraceAction[curE] == 'X'
				backtraceD[curE] = curD;
				backtraceMatched[curE] = L(curE, curD) - L(curE - 1, curD) - 1;
			}
			curD = backtraceD[curE];
#ifdef TRACE_LV
			printf("%d %d: %d %c %d %d\n", curE, curD, L(curE, curD),
				backtraceAction[curE], backtraceD[curE], backtraceMatched[curE]);
#endif
		}

		int curE = 1;
		int offset = L(0, 0);
		_ASSERT(*o_netIndel == 0);
		while (curE <= e) {
			// First write the action, possibly with a repeat if it occurred multiple times with no exact matches
			char action = backtraceAction[curE];
			int actionCount = 1;
			while (curE + 1 <= e && backtraceMatched[curE] == 0 && backtraceAction[curE + 1] == action) {
				actionCount++;
				curE++;
			}
			if (action == 'I') {
				*matchProbability *= lv_indelProbabilities[actionCount];
				offset += actionCount;
				*o_netIndel += actionCount;
			}
			else if (action == 'D') {
				*matchProbability *= lv_indelProbabilities[actionCount];
				offset -= actionCount;
				*o_netIndel -= actionCount;
			}
			else {
				_ASSERT(action == 'X');
				for (int i = 0; i < actionCount; i++) {
					*matchProbability *= lv_phredToProbability[qualityString[/*BUGBUG - think about what to do here*/__min(patternLen - 1, __max(offset, 0))]];
					offset++;
				}
			}

			offset += backtraceMatched[curE];   // Skip over the matching bases.
			curE++;
		}

		*matchProbability *= lv_perfectMatchProbability[patternLen - e]; // Accounting for the < 1.0 chance of no changes for matching bases
	} else {
		//
		// Not tracking match probability.
		//
	}

	_ASSERT(e <= k);
	return e;
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

    void *operator new(size_t size) {return BigAlloc(size);}
    void operator delete(void *ptr) {BigDealloc(ptr);}

    void *operator new(size_t size, BigAllocator *allocator) {_ASSERT(size == sizeof(LandauVishkin<TEXT_DIRECTION>)); return allocator->allocate(size);}
    void operator delete(void *ptr, BigAllocator *allocator) {/*Do nothing.  The memory is freed when the allocator is deleted.*/}
 
private:
    //
    // Count characters of a perfect match until a mismatch or the end of one or the other string, the
    // minimum length of which is represented by the end parameter.  Advances p & t to the first mismatch
    // or first character beyond the end.
    //
    inline int countPerfectMatch(const char *& p, const char *& t, int availBytes)      // This is essentially duplicated in LandauVishkinWithCigar
    {
	    const char *pBase = p;
	    const char *pend = p + availBytes;
	    while (true) {
		    _uint64 x;
		    if (TEXT_DIRECTION == 1) {
			    x = *((_uint64*)p) ^ *((_uint64*)t);
		    } else {
			    _uint64 T = *(_uint64 *)(t - 7);
			    _uint64 tSwap = ByteSwapUI64(T);
			    x = *((_uint64*)p) ^ tSwap;
		    }

		    if (x) {
			    unsigned long zeroes;
			    CountTrailingZeroes(x, zeroes);
			    zeroes >>= 3;
			    return __min((int)(p - pBase) + (int)zeroes, availBytes);
		    } // if (x)

		    p += 8;
		    if (p >= pend) {
			    return availBytes;
		    }

		    t += 8 * TEXT_DIRECTION;
	    } // while true

	    return 0;
    }


    //
    // Table of d values for the inner loop in computeEditDistance.  This allows us to avoid the line d = (d > 0 ? -d : -d+1), which causes
    // a branch misprediction every time.
    //
    int dTable[2 * (MAX_K + 1) + 1];
    
	//
	// Note on state arrays:
	// 
	// We have several arrays that need to be indexed on edit distance and net indels.  Because net indels is signed, we want them
	// to have their second coordinate (d) run from [-MAX_K .. MAX_K].  When computing the opening or closing of an indel, we add more than
	// one edit distance, which means we compute LInsert[e][x] based on L[e-OpenPenalty][x+1].  When we're at edit distance < the gap open
	// penalty, of course we can't fill in LInsert; however, rather than just checking it each time (and incurring a branch prediction miss), 
	// we just let it happily index into negative space, which is initialized with -2 and so will never be used.  So, we want our arrays to
	// run from [-MAX_GAP..MAX_K][-MAX_K .. MAX_K].  To do this, we just allocate the space and compute a pointer that would be at [0][0].
	// We use macros to do the indexing, because it's tricky to convince C++ to do this kind of thing statically.
	//
	// Also, conceptually these arrays are local to each computation.  They're here to save memory allocation and initialization overhead.
	// Note that the important parts for the initialization is never overwritten, though the rest is.
	//


    // TODO: For long reads, we should include a version that only has L be 2 x (2*MAX_K+1) cells
    int L_space[(MAX_K + 1) * (2 * MAX_K + 1)];
	int *L_zero;

    // Action we did to get to each position: 'D' = deletion, 'I' = insertion, 'X' = substitution.  This is needed to compute match probability.
	char A_space[(MAX_K + 1) * (2 * MAX_K + 1)];
	char *A_zero;

    // Arrays for backtracing the actions required to match two strings
    char backtraceAction[MAX_K+1];
    int  backtraceMatched[MAX_K+1];
    int  backtraceD[MAX_K+1];

#undef  L
#undef  A
};

void setLVProbabilities(double *i_indelProbabilities, double *i_phredToProbability, double mutationProbability);
void initializeLVProbabilitiesToPhredPlus33();


// Computes the edit distance between two strings and returns a CIGAR string for the edits.

enum CigarFormat
{
    COMPACT_CIGAR_STRING = 0,
    EXPANDED_CIGAR_STRING = 1,
    COMPACT_CIGAR_BINARY = 2,
    BAM_CIGAR_OPS = 3,
};

// express cigar as 2 byte per reference base summarizing the changes
// at that location; may lose information for longer inserts
enum LinearCigarFlags
{
    CigarInsertFlags = 0xfff8,  // 13 bits for insert
    CigarInsertCShift= 3,       // shift to get count
    CigarInsertCount = 0x7,     // after shifting, mask to get # insertions
    CigarInsertBShift= 6,       // shift to get bases
    CigarInsertBases = 0xffc0,   // up to 5 inserted bases, low-order bits first
    CigarInsertNBases= 5,       

    CigarOpcode      = 0x07,    // opcode for ref base
    CigarNoop        = 0x00,    // no change
    CigarReplace     = 0x01,    // base is n-1
    CigarDelete      = 0x05,    // delete
};

class LandauVishkinWithCigar {
public:
    LandauVishkinWithCigar();

    // Compute the edit distance between two strings and write the CIGAR string in cigarBuf.
    // Returns -1 if the edit distance exceeds k or -2 if we run out of space in cigarBuf.
    int computeEditDistance(const char* text, int textLen, const char* pattern, int patternLen, int k,
                            char* cigarBuf, int cigarBufLen, bool useM,
                            CigarFormat format = COMPACT_CIGAR_STRING,
                            int* o_cigarBufUsed = NULL,
                            int* o_textUsed = NULL,
                            int *o_netIndel = NULL);

    // same, but places indels as early as possible, following BWA & VCF conventions
    int computeEditDistanceNormalized(const char* text, int textLen, const char* pattern, int patternLen, int k,
                            char* cigarBuf, int cigarBufLen, bool useM,
                            CigarFormat format = COMPACT_CIGAR_STRING,
                            int* o_cigarBufUsed = NULL,
                            int* o_textUsed = NULL,
                            int *o_netIndel = NULL);

    // take a compact cigar binary format and turn it into one byte per reference base
    // describing the difference from the reference at that location
    // might lose information for large inserts
    // returns number of bytes in result
    static int linearizeCompactBinary(_uint16* o_linear, int referenceSize,
        char* cigar, int cigarSize, char* sample, int sampleSize);

    static void printLinear(char* buffer, int bufferSize, unsigned variant);
private:
    int L[MAX_K+1][2 * MAX_K + 1];
    
    // Action we did to get to each position: 'D' = deletion, 'I' = insertion, 'X' = substitution.
    char A[MAX_K+1][2 * MAX_K + 1];

    //
    // Total (not net) indels at this point.  Parallel to L and A arrays.  Used to select the least-indel path
    // consistent with the lowest edit distance.
    //
    int totalIndels[MAX_K + 1][2 * MAX_K + 1];

    // Arrays for backtracing the actions required to match two strings
    char backtraceAction[MAX_K+1];
    int backtraceMatched[MAX_K+1];
    int backtraceD[MAX_K+1];
};
