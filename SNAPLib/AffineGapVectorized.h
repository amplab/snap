#pragma once

#include "Compat.h"
#include "FixedSizeMap.h"
#include "BigAlloc.h"
#include "exit.h"
#include "Genome.h"
#include "Read.h"
#include "Tables.h"
#include "LandauVishkin.h"
#include "AffineGap.h"


const int VEC_SIZE = 8; // Number of 16-bit elements in vector

const int MAX_VEC_SEGMENTS = (MAX_READ_LENGTH + VEC_SIZE - 1) / (VEC_SIZE);

//
// These are global so there are only one for both senses of the template
//

extern double *lv_indelProbabilities;  // Maps indels by length to probability of occurance.
extern double *lv_phredToProbability;  // Maps ASCII phred character to probability of error, including 
extern double *lv_perfectMatchProbability; // Probability that a read of this length has no mutations

// #define TRACE_AG

enum BacktraceActionType {
    M = 0, 
    D = 1,
    I = 2
};

// 
// Computes the affine gap score between two strings based on Farrar's algorithm
//
template<int TEXT_DIRECTION = 1> class AffineGapVectorized {

public:
    AffineGapVectorized()
    {
        init(1, 4, 6, 1);
    }

    AffineGapVectorized(
        int i_matchReward,
        int i_subPenalty,
        int i_gapOpenPenalty,
        int i_gapExtendPenalty)
    {
        init(i_matchReward, i_subPenalty, i_gapOpenPenalty, i_gapExtendPenalty);
    }

    static size_t getBigAllocatorReservation() { return sizeof(AffineGapVectorized<TEXT_DIRECTION>); }

    ~AffineGapVectorized()
    {
    }

    void init(
        int i_matchReward,
        int i_subPenalty,
        int i_gapOpenPenalty,
        int i_gapExtendPenalty)
    {
        matchReward = i_matchReward;
        subPenalty = -i_subPenalty;
        gapOpenPenalty = i_gapOpenPenalty + i_gapExtendPenalty; // First gap costs (gapOpen + gapExtend)
        gapExtendPenalty = i_gapExtendPenalty;

        //
        // Initialize nucleotide <-> nucleotide transition matrix
        //
        int i, j, k = 0;
        for (i = 0; i < (MAX_ALPHABET_SIZE - 1); i++) {
            for (j = 0; j < (MAX_ALPHABET_SIZE - 1); j++) {
                ntTransitionMatrix[k++] = (i == j) ? matchReward : subPenalty;
            }
            ntTransitionMatrix[k++] = -1; // FIXME: What penalty to use for N ?
        }
        for (i = 0; i < MAX_ALPHABET_SIZE; i++) {
            ntTransitionMatrix[k++] = -1;
        }
    }

    int computeScore(
        const char* text,
        int textLen,
        const char* pattern,
        const char *qualityString,
        int patternLen,
        int w,
        int scoreInit,
        int *o_textOffset = NULL,
        int *o_patternOffset = NULL,
        int *o_nEdits = NULL,
        double *matchProbability = NULL)
    {

#ifdef TRACE_AG
        printf("\n");
#endif

        _ASSERT(w < MAX_K);
        _ASSERT(textLen <= MAX_READ_LENGTH + MAX_K);

        int localTextOffset, localPatternOffset;
        if (NULL == o_textOffset) {
            //
            // If the user doesn't want textOffset, just use a stack local to avoid
            // having to check it all the time.
            //
            o_textOffset = &localTextOffset;
        }

        if (NULL == o_patternOffset) {
            o_patternOffset = &localPatternOffset;
        }

        w = __min(MAX_K - 1, w); // enforce limit even in non-debug builds

        if (NULL == text) {
            // This happens when we're trying to read past the end of the genome.
            if (NULL != matchProbability) {
                *matchProbability = 0.0;
            }
            if (NULL != o_nEdits) {
                *o_nEdits = -1;
            }
            return -1;
        }

        // FIXME: Find better way to do this based on usage
        double localMatchProbability;
        if (NULL == matchProbability) {
            //
            // If the user doesn't want matchProbability, just use a stack local to avoid
            // having to check it all the time.
            //
            matchProbability = &localMatchProbability;
        }

        int localnEdits;
        if (NULL == o_nEdits) {
            //
            // If the user doesn't want matchProbability, just use a stack local to avoid
            // having to check it all the time.
            //
            o_nEdits = &localnEdits;
        }

        //
        // Start with perfect match probability and work our way down.
        //
        *matchProbability = 1.0;

        if (TEXT_DIRECTION == -1) {
            text--; // so now it points at the "first" character of t, not after it.
        }
        
        //
        //       PATTERN (p)
        //  VEC0            VEC1            
        //  p[0]p[1]..      p[8]p[9]..      
        // 
        //
        //   STRIPED_PATTERN (p)
        //  VEC0            VEC1
        //  p[0]p[8]..      p[1]p[9]..
        //  
        int numVec = (patternLen + VEC_SIZE - 1) / VEC_SIZE; // Number of vector segments
        int paddedPatternLen = numVec * VEC_SIZE;
        int patternIdx = 0;

        // 
        // Generate query profile
        //
        int16_t* queryResult = (int16_t*)qProfile;
        for (int i = 0; i < MAX_ALPHABET_SIZE; i++) {
            for (int j = 0; j < numVec; j++) {
                for (int k = j; k < paddedPatternLen; k += numVec) {
                    if (k < patternLen) {
                        uint8_t bp = BASE_VALUE[pattern[k]];
                        queryResult[patternIdx] = ntTransitionMatrix[i * MAX_ALPHABET_SIZE + bp];
                    }
                    else {
                        queryResult[patternIdx] = INT16_MIN;
                    }
                    patternIdx++;
                }
            }
        }

        // 
        // Define constants in their vector form
        //
        __m128i v_zero = _mm_setzero_si128();
        __m128i v_one = _mm_set1_epi16(1);
        __m128i v_two = _mm_set1_epi16(2);
        __m128i v_four = _mm_set1_epi16(4);
        __m128i v_thirtytwo = _mm_set1_epi16(32);
        __m128i v_gapOpen = _mm_set1_epi16(gapOpenPenalty);
        __m128i v_gapExtend = _mm_set1_epi16(gapExtendPenalty);
        __m128i v_scoreInit = _mm_set1_epi16(scoreInit);

        //
        // Initialize scores of first row
        //
        uint16_t scoreFirstRow[VEC_SIZE] = {};
        for (int vecIdx = 0; vecIdx < numVec; vecIdx++) {
            for (int elemIdx = 0; elemIdx < VEC_SIZE; elemIdx++) {
                int patternIdx = elemIdx * numVec + vecIdx;
                if (patternIdx < patternLen) {
                    scoreFirstRow[elemIdx] = __max(0, scoreInit - gapOpenPenalty - patternIdx * gapExtendPenalty);
                }
            }
            _ASSERT(VEC_SIZE == 8);  // FIXME: Initialization below works only when VEC_SIZE = 8
            _mm_store_si128(H + vecIdx, _mm_setr_epi16(scoreFirstRow[0], scoreFirstRow[1], scoreFirstRow[2], scoreFirstRow[3],
                                                        scoreFirstRow[4], scoreFirstRow[5], scoreFirstRow[6], scoreFirstRow[7]));
            _mm_store_si128(E + vecIdx, v_zero);
        }

        int score = -1; // Final alignment score to be returned. 

        // We keep track of the best score and text offset for global and local alignment separately.
        // These are used to choose between global and local alignment for the {text, pattern} pair
        int bestGlobalAlignmentScore = -1;
        int bestGlobalAlignmentTextOffset = -1; 
        int bestLocalAlignmentScore = -1;
        int bestLocalAlignmentTextOffset = -1;
        int bestLocalAlignmentPatternOffset = -1;

        *o_textOffset = -1; // # Characters of text used for aligning against the pattern
        *o_patternOffset = -1; // # Characters of pattern used for obtaining the maximum score. Ideally we would like to use the full pattern
        *o_nEdits = -1;
        
        __m128i* Hptr = H;
        __m128i* Hminus1ptr = Hminus1;

        // Iterate over all rows of text
        for (int i = 0; i < textLen; i++) {

            const char* t = (text + i * TEXT_DIRECTION);

            // Get the query profile for the row
            __m128i* qRowProfile = qProfile + BASE_VALUE[*t] * numVec;

            // Registers to hold intermediate scores for each row
            __m128i m, h, temp, e = v_zero, f = v_zero, max = v_zero;

            int maxScoreRow = 0;

            // Load h from the previous row
            h = _mm_load_si128(Hptr + numVec - 1);

            // Shift left h and blend in initial values
            h = _mm_slli_si128(h, 2); // shift h left by 2 * 8 bits
            int hInit = scoreInit;
            if (i > 0) {
                hInit = __max(0, scoreInit - gapOpenPenalty - (i - 1) * gapExtendPenalty);
            }
            __m128i v_hInit = _mm_set1_epi16(hInit);
            h = _mm_blend_epi16(h, v_hInit, 1);

            for (int j = 0; j < numVec; j++) {

                __m128i backtraceActionVec;

                __m128i mask = _mm_cmpgt_epi16(h, v_zero); // h > 0

                // Below implements: m = m > 0 : m + qRowProfile[j] : 0 
                m = _mm_adds_epi16(h, *qRowProfile++);
                m = _mm_and_si128(m, mask);

                e = _mm_load_si128(E + j);

                // h = max{m, e, f}
                backtraceActionVec = _mm_and_si128(_mm_cmpgt_epi16(e, m), v_one); // action = e > m ? 1 : 0
                h = _mm_max_epi16(m, e);
                __m128i tmpResult = _mm_and_si128(_mm_cmpgt_epi16(f, h), v_two); 
                backtraceActionVec = _mm_or_si128(tmpResult, _mm_andnot_si128(tmpResult, backtraceActionVec)); // action = f > h ? 2 : action 
                h = _mm_max_epi16(h, f);

                max = _mm_max_epi16(max, h);

                // Store h for the next row
                _mm_store_si128(Hminus1ptr + j, h);

                // e = max{m - gapOpen, e - gapExtend}
                e = _mm_subs_epi16(e, v_gapExtend);
                temp = _mm_subs_epi16(m, v_gapOpen);
                temp = _mm_max_epi16(temp, v_zero);
                tmpResult = _mm_and_si128(_mm_cmpgt_epi16(e, temp), v_four);
                backtraceActionVec = _mm_or_si128(backtraceActionVec, tmpResult);
                e = _mm_max_epi16(e, temp);
                _mm_store_si128(E + j, e);

                // f = max{m - gapOpen, f- gapExtend}
                f = _mm_subs_epi16(f, v_gapExtend);
                tmpResult = _mm_and_si128(_mm_cmpgt_epi16(f, temp), v_thirtytwo);
                backtraceActionVec = _mm_or_si128(backtraceActionVec, tmpResult);
                f = _mm_max_epi16(f, temp);

                // Store traceback information
                _mm_store_si128(backtraceAction + (i * numVec + j), backtraceActionVec);

                // Load the next score vector
                h = _mm_load_si128(Hptr + j);

            } // end pattern 

            //
            // Farrar's algorithm does lazy f evaluation. Since f rarely influences final score h, the algorithm speculates f = zero initially 
            // Re-evaluate if f could influence h after we have made a first pass through the row
            //
            for (int k = 0; k < VEC_SIZE; k++) {

                // Extract the last f vector
                f = _mm_slli_si128(f, 2);

                for (int j = 0; j < numVec; j++) {
                    
                    h = _mm_load_si128(Hminus1ptr + j);

                    // action = f > h ? 2 : action
                    __m128i tmpResult = _mm_and_si128(_mm_cmpgt_epi16(f, h), v_two);
                    __m128i backtraceActionVec = _mm_load_si128(backtraceAction + (i * numVec + j));
                    __m128i tmpResult2 = _mm_andnot_si128(tmpResult, backtraceActionVec);
                    backtraceActionVec = _mm_or_si128(tmpResult, tmpResult2); 

                    h = _mm_max_epi16(h, f);
                    _mm_store_si128(Hminus1ptr + j, h);

                    max = _mm_max_epi16(max, h);

                    temp = _mm_subs_epu16(h, v_gapOpen);
                    f = _mm_subs_epu16(f, v_gapExtend);

                    tmpResult = _mm_and_si128(_mm_cmpgt_epi16(f, temp), v_thirtytwo);
                    backtraceActionVec = _mm_or_si128(backtraceActionVec, tmpResult);
                    _mm_store_si128(backtraceAction + (i * numVec + j), backtraceActionVec);

                    // Converged if no element of f can influence h
                    bool converged = !(_mm_movemask_epi8(_mm_cmpgt_epi16(f, temp))); // !! converts int to bool

                    if (converged) goto got_answer;
                }
            }

        got_answer:
            maxScoreRow = getMax(max);
#ifdef TRACE_AG
            __m128i rowScore;
            for (int i = 0; i < numVec; i++) {
                rowScore = _mm_load_si128(Hminus1ptr + i);
                printElem(rowScore);
            }
            printf("\n");
#endif
            
            // Global alignment score (i.e., score when aligning to the end of the pattern)
            __m128i v_globalAlignmentScore = _mm_load_si128(Hminus1ptr + ((patternLen - 1) % numVec));
            int globalAlignmentScore = getElem((patternLen - 1) / numVec, v_globalAlignmentScore);
            if (globalAlignmentScore >= bestGlobalAlignmentScore) {
                bestGlobalAlignmentScore = globalAlignmentScore;
                bestGlobalAlignmentTextOffset = i;
            }

            if (maxScoreRow == 0) break;

            if (maxScoreRow > bestLocalAlignmentScore) { // If we obtained a better score this round  
                // Get index in pattern where maximum score was obtained
                for (int j = 0; j < numVec; ++j) {
                    __m128i h = _mm_load_si128(Hminus1ptr + j); // Load a vector of score values 
                    __m128i mask = _mm_cmpeq_epi16(h, _mm_set1_epi16(maxScoreRow)); // Create a mask of 1's in cells where score = maxScoreRow
                    int result = _mm_movemask_epi8(_mm_packs_epi16(mask, v_zero)); // Convert vector result to 8-bit
                    unsigned long zeroes;
                    if (result) {
                        CountLeadingZeroes(result, zeroes);
                        int patternOffset = zeroes * numVec + j;
                        bestLocalAlignmentPatternOffset = (bestLocalAlignmentPatternOffset > patternOffset) ? bestLocalAlignmentPatternOffset : patternOffset;
                        _ASSERT(bestLocalAlignmentPatternOffset < patternLen); // Scores outside the pattern boundaries should never be the maximum in the row
                    }
                }
                bestLocalAlignmentScore = maxScoreRow;
                bestLocalAlignmentTextOffset = i;
            }
            
            // Swap roles of H and Hminus1 for the next row
            __m128i* hTemp = Hminus1ptr; 
            Hminus1ptr = Hptr; 
            Hptr = hTemp;

        } // end text

        score = bestGlobalAlignmentScore;
        *o_textOffset = bestGlobalAlignmentTextOffset;

        // Choose between local and global alignment for patternOffset
        if ((bestLocalAlignmentScore != bestGlobalAlignmentScore) && (bestLocalAlignmentScore >= bestGlobalAlignmentScore + 5)) { // FIXME: Change 5 to a clipping penalty
            // Local alignment preferred
            *o_patternOffset = bestLocalAlignmentPatternOffset;
        }
        else {
            // Global alignment preferred
            *o_patternOffset = patternLen - 1;
        }

        if (score > scoreInit) {
            // Perform traceback and compute nEdits assuming the entire pattern is aligned
            int rowIdx = *o_textOffset, colIdx = patternLen - 1, matrixIdx;
            BacktraceActionType action = M, prevAction = M;
            int actionCount = 1;
            int nMatches = 0, nMismatches = 0, nGaps = 0;

            // Start traceback from the cell (i,j) with the maximum score
            while (rowIdx >= 0 && colIdx >= 0) {
                uint16_t* backtracePointersRow = (uint16_t*)(backtraceAction + (rowIdx * numVec));
                // 
                // The traceback matrix (H, E, or F) we need to look at depends on the current action.
                // We index the corresponding matrix using the current action as described below:
                // If current action is M, bits[1:0] are used 
                // If current action is D, bits[3:2] are used
                // If current action is I, bits[5:4] are used
                //
                matrixIdx = action << 1; 
                
                // 
                // Two bits are used to encode the backtrace action type
                //
                int stripedColIdx = (colIdx % numVec) * VEC_SIZE + (colIdx / numVec);
                action = (BacktraceActionType) ((backtracePointersRow[stripedColIdx] >> matrixIdx) & 3);
                if (action == M) {
                    if (pattern[colIdx] != text[rowIdx * TEXT_DIRECTION]) {
                        // Compute probabilties of mismatches
                        *matchProbability *= lv_phredToProbability[qualityString[colIdx]];
                        nMismatches++;
                    }
                    else {
                        nMatches++;
                    }
                    rowIdx--; colIdx--;
                }
                else if (action == D) {
                    rowIdx--;
                }
                else { // action == I
                    colIdx--;
                    action = I; // FIXME: we choose I, whenever action == 2 or 3. Make action == 2, since matrixIdx depends on it
                }
                // Ignore transitions from H. Compute number of indels in sequence for matchProbability
                if (prevAction != M) {
                    if (prevAction == action) {
                        actionCount++;
                    }
                    else {
                        // Update match probabilties with indel probabilties here
                        nGaps += actionCount;
                        *matchProbability *= lv_indelProbabilities[actionCount];
                        actionCount = 1;
                    }
                }
                prevAction = action;
            }
            if (rowIdx >= 0) { // deletion just after the seed
                actionCount = rowIdx + 1;
                nGaps += actionCount;
                *matchProbability *= lv_indelProbabilities[actionCount];
            }
            if (colIdx >= 0) { // insertion just after the seed
                actionCount = colIdx + 1;
                nGaps += actionCount;
                *matchProbability *= lv_indelProbabilities[actionCount];
            }

            *o_nEdits = nMismatches + nGaps;

            *o_nEdits = (*o_nEdits <= w) ? *o_nEdits : -1; // return -1 if we have more edits than threshold w

            *matchProbability *= lv_perfectMatchProbability[nMatches];

            *o_textOffset += 1;
            *o_patternOffset += 1;

            *o_textOffset = patternLen - *o_textOffset;
            *o_patternOffset = patternLen - *o_patternOffset;
            
            return score;
        }
        else {
            return -1; 
        }
    }

    void *operator new(size_t size) { return BigAlloc(size); }
    void operator delete(void *ptr) { BigDealloc(ptr); }

    void *operator new(size_t size, BigAllocator *allocator) { _ASSERT(size == sizeof(AffineGapVectorized<TEXT_DIRECTION>)); return allocator->allocate(size); }
    void operator delete(void *ptr, BigAllocator *allocator) {/*Do nothing.  The memory is freed when the allocator is deleted.*/ }

private:

    //
    // Scores for each nucleotide <-> nucleotide transition
    //
    int ntTransitionMatrix[MAX_ALPHABET_SIZE * MAX_ALPHABET_SIZE];

    //
    // Affine gap scoring parameters
    //
    int matchReward;
    int subPenalty;
    int gapOpenPenalty;
    int gapExtendPenalty;

    __m128i qProfile[MAX_ALPHABET_SIZE * MAX_VEC_SEGMENTS];
    __m128i H[MAX_VEC_SEGMENTS];
    __m128i Hminus1[MAX_VEC_SEGMENTS];
    __m128i E[MAX_VEC_SEGMENTS];

    //
    // Pointers to traceback alignment, one for each of the three affine-gap matrices, F, E and H
    // Traceback actions are encoded as follows: Bit 5 to Bit 0 - F[5:4], E[3:2], H[1:0]
    // Although we need only 1 bit to encode paths into E and F matrix, we use 2 bits to simplify decoding
    //
    __m128i backtraceAction[(MAX_READ_LENGTH + MAX_K) * MAX_VEC_SEGMENTS];

    //
    // Helper functions for operations on SSE registers
    //
    inline int getMax(__m128i x) {
        x = _mm_max_epi16(x, _mm_srli_si128((x), 8));
        x = _mm_max_epi16(x, _mm_srli_si128((x), 4));
        x = _mm_max_epi16(x, _mm_srli_si128((x), 2));
        return _mm_extract_epi16(x, 0);
    }

    inline int getElem(int index, __m128i x) {
        int elem;
        switch (index) {
            case 0: elem = _mm_extract_epi16(x, 0); break;
            case 1: elem = _mm_extract_epi16(x, 1); break;
            case 2: elem = _mm_extract_epi16(x, 2); break;
            case 3: elem = _mm_extract_epi16(x, 3); break;
            case 4: elem = _mm_extract_epi16(x, 4); break;
            case 5: elem = _mm_extract_epi16(x, 5); break;
            case 6: elem = _mm_extract_epi16(x, 6); break;
            case 7: elem = _mm_extract_epi16(x, 7); break;
            default: elem = -1; break;
        }
        return elem;
    }

    inline void printElem(__m128i x) {
        printf("%u,%u,%u,%u,%u,%u,%u,%u,", _mm_extract_epi16(x, 0), _mm_extract_epi16(x, 1), _mm_extract_epi16(x, 2), _mm_extract_epi16(x, 3),
            _mm_extract_epi16(x, 4), _mm_extract_epi16(x, 5), _mm_extract_epi16(x, 6), _mm_extract_epi16(x, 7));
    }
};
