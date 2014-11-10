#pragma once
#include "Compat.h"
#include "FixedSizeMap.h"
#include "BigAlloc.h"
#include "exit.h"
#include "Genome.h"
#include "GTFReader.h"
#include <vector>

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
public:
    LandauVishkin()
{
    if (TEXT_DIRECTION != 1 && TEXT_DIRECTION != -1) {
        fprintf(stderr, "You can't possibly be serious.\n");
        soft_exit(1);
    }

    memsetint(L[0], -2, (MAX_K+1)*(2*MAX_K+1));

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
	// the largest possible edit distance, and hence d can never be less than -MAX_K, so MAX_K + d >= 0.
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
                int *netIndel = NULL)   // the net of insertions and deletions in the alignment.  Negative for insertions, positive for deleteions (and 0 if there are non in net).  Filled in only if matchProbability is non-NULL
{
    int localNetIndel;
	int d;
    if (NULL == netIndel) {
        //
        // If the user doesn't want netIndel, just use a stack local to avoid
        // having to check it all the time.
        //
        netIndel = &localNetIndel;
    }
    _ASSERT(k < MAX_K);

    *netIndel = 0;

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
            L[0][MAX_K] = __min((int)(p - pattern + zeroes), end);
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
				//
				// We're through on this iteration.
				//

				if ('X' == A[e][MAX_K + d]) {
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

            L[e][MAX_K+d] = best;
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
			backtraceAction[curE] = A[curE][MAX_K + curD];
			if (backtraceAction[curE] == 'I') {
				backtraceD[curE] = curD + 1;
				backtraceMatched[curE] = L[curE][MAX_K + curD] - L[curE - 1][MAX_K + curD + 1] - 1;
			} else if (backtraceAction[curE] == 'D') {
				backtraceD[curE] = curD - 1;
				backtraceMatched[curE] = L[curE][MAX_K + curD] - L[curE - 1][MAX_K + curD - 1];
			} else { // backtraceAction[curE] == 'X'
				backtraceD[curE] = curD;
				backtraceMatched[curE] = L[curE][MAX_K + curD] - L[curE - 1][MAX_K + curD] - 1;
			}
			curD = backtraceD[curE];
#ifdef TRACE_LV
			printf("%d %d: %d %c %d %d\n", curE, curD, L[curE][MAX_K + curD],
				backtraceAction[curE], backtraceD[curE], backtraceMatched[curE]);
#endif
		}

		int curE = 1;
		int offset = L[0][MAX_K + 0];
		_ASSERT(*netIndel == 0);
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
				*netIndel += actionCount;
			}
			else if (action == 'D') {
				*matchProbability *= lv_indelProbabilities[actionCount];
				offset -= actionCount;
				*netIndel -= actionCount;
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
    // TODO: For long reads, we should include a version that only has L be 2 x (2*MAX_K+1) cells
    int L[MAX_K+1][2 * MAX_K + 1];

    //
    // Table of d values for the inner loop in computeEditDistance.  This allows us to avoid the line d = (d > 0 ? -d : -d+1), which causes
    // a branch misprediction every time.
    //
    int dTable[2 * (MAX_K + 1) + 1];
    
    // Action we did to get to each position: 'D' = deletion, 'I' = insertion, 'X' = substitution.  This is needed to compute match probability.
    char A[MAX_K+1][2 * MAX_K + 1];

    // Arrays for backtracing the actions required to match two strings
    char backtraceAction[MAX_K+1];
    int  backtraceMatched[MAX_K+1];
    int  backtraceD[MAX_K+1];
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
                            char* cigarBuf, int cigarBufLen, bool useM, std::vector<unsigned> &tokens,
                            CigarFormat format = COMPACT_CIGAR_STRING, int* o_cigarBufUsed = NULL, int* o_textUsed = NULL);
     
    // same, but places indels as early as possible, following BWA & VCF conventions
    int computeEditDistanceNormalized(const char* text, int textLen, const char* pattern, int patternLen, int k,
                            char* cigarBuf, int cigarBufLen, bool useM, 
                            std::vector<unsigned> &tokens,
                            CigarFormat format = COMPACT_CIGAR_STRING, 
                            int* cigarBufUsed = NULL, 
                            int* o_addFrontClipping = NULL);
                       
    int insertSpliceJunctions(const GTFReader *gtf, std::vector<unsigned> tokens, std::string transcript_id, unsigned pos, char *cigarNew, int cigarNewLen, CigarFormat format = COMPACT_CIGAR_STRING);

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

    // Arrays for backtracing the actions required to match two strings
    char backtraceAction[MAX_K+1];
    int backtraceMatched[MAX_K+1];
    int backtraceD[MAX_K+1];
};
