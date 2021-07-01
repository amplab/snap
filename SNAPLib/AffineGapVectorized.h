//
// DP formulation based on https://figshare.com/articles/preprint/Notes_on_pairwise_alignment_with_dynamic_programming/5223973
//

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

//
// Helper functions for operations on SSE registers
//
static inline int getMax(__m128i x) {
    x = _mm_max_epi16(x, _mm_srli_si128((x), 8));
    x = _mm_max_epi16(x, _mm_srli_si128((x), 4));
    x = _mm_max_epi16(x, _mm_srli_si128((x), 2));
    return _mm_extract_epi16(x, 0);
}

static inline int16_t getElem(int index, __m128i x) {
    int16_t elem;
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

static inline void printElem(__m128i x) {
    printf("%hi,%hi,%hi,%hi,%hi,%hi,%hi,%hi,", _mm_extract_epi16(x, 0), _mm_extract_epi16(x, 1), _mm_extract_epi16(x, 2), _mm_extract_epi16(x, 3),
        _mm_extract_epi16(x, 4), _mm_extract_epi16(x, 5), _mm_extract_epi16(x, 6), _mm_extract_epi16(x, 7));
}

// return (~mask)*v1 | (mask)*v2
static inline __m128i blend_sse(__m128i v1, __m128i v2, __m128i mask) {
    v1 = _mm_andnot_si128(mask, v1);
    v1 = _mm_or_si128(v1, _mm_and_si128(mask, v2));
    return v1;
}

enum BacktraceActionType {
    M = 0, 
    D = 1,
    I = 2,
    X = 3
};

// 
// Computes the affine gap score between two strings based on Farrar's algorithm
//
template<int TEXT_DIRECTION = 1> class AffineGapVectorized {

public:
    AffineGapVectorized()
    {
        init(1, 4, 6, 1, 10, 5);
    }

    AffineGapVectorized(
        int i_matchReward,
        int i_subPenalty,
        int i_gapOpenPenalty,
        int i_gapExtendPenalty,
        int i_fivePrimeEndBonus,
        int i_threePrimeEndBonus)
    {
        init(i_matchReward, i_subPenalty, i_gapOpenPenalty, i_gapExtendPenalty, i_fivePrimeEndBonus, i_threePrimeEndBonus);
    }

    static size_t getBigAllocatorReservation() { return sizeof(AffineGapVectorized<TEXT_DIRECTION>); }

    ~AffineGapVectorized()
    {
    }

    void init(
        int i_matchReward,
        int i_subPenalty,
        int i_gapOpenPenalty,
        int i_gapExtendPenalty,
        int i_fivePrimeEndBonus,
        int i_threePrimeEndBonus)
    {
        matchReward = i_matchReward;
        subPenalty = -i_subPenalty;
        gapOpenPenalty = i_gapOpenPenalty + i_gapExtendPenalty; // First gap costs (gapOpen + gapExtend)
        gapExtendPenalty = i_gapExtendPenalty;
        fivePrimeEndBonus = i_fivePrimeEndBonus;
        threePrimeEndBonus = i_threePrimeEndBonus;

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

    //
    // Compute Hamming distance between text and pattern. Poorly matching sections of pattern are clipped out
    // TODO: Make faster version
    // 
    int computeGaplessScore(
        const char* text,
        int textLen,
        const char* pattern,
        const char* qualityString,
        int patternLen,
        int scoreInit,
        int scoreLimit,
        int* o_nEdits = NULL,
        int* o_textOffset = NULL,
        int* o_patternOffset = NULL,
        double* matchProbability = NULL,
        int* o_nEditsGapless = NULL)
    {
        _ASSERT(textLen <= MAX_READ_LENGTH + MAX_K);
        _ASSERT(patternLen <= textLen);

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

        int localnEditsGapless;
        if (NULL == o_nEditsGapless) {
            //
            // If the user doesn't want matchProbability, just use a stack local to avoid
            // having to check it all the time.
            //
            o_nEditsGapless = &localnEditsGapless;
        }

        if (scoreLimit < 0 || NULL == text) {
            *o_nEdits = ScoreAboveLimit;
            *o_nEditsGapless = ScoreAboveLimit;
            return ScoreAboveLimit;
        }

        //
        // Start with perfect match probability and work our way down.
        //
        *matchProbability = 1.0;

        if (TEXT_DIRECTION == -1) {
            text--; // so now it points at the "first" character of t, not after it.
        }

        int gapLessScore = scoreInit;
        int maxScore = scoreInit;
        for (int i = 0; i < patternLen; i++) {
            const char* p = pattern + i;
            const char* t = text + i * TEXT_DIRECTION;
            gapLessScore += (*p == *t) ? matchReward : subPenalty;

            // Keep track of maximum score and corresponding length of pattern matched
            if (gapLessScore > maxScore) {
                maxScore = gapLessScore;
                *o_patternOffset = i;
            }

        } // patternLen

        // We were able to extend the seed. Now compute matchProbability for the portion that was extended
        if (maxScore > scoreInit) {
            *o_nEdits = 0;
            int nMatches = 0;
            for (int i = 0; i <= *o_patternOffset; i++) {
                const char* p = pattern + i;
                const char* t = text + i * TEXT_DIRECTION;
                if (*p != *t) {
                    *o_nEdits += 1;
                    *matchProbability *= lv_phredToProbability[qualityString[i]];
                } else {
                    nMatches++;
                }
            }
            *matchProbability *= lv_perfectMatchProbability[nMatches];
            *o_patternOffset += 1;
            *o_patternOffset = patternLen - *o_patternOffset;
            *o_textOffset = *o_patternOffset;
            *o_nEditsGapless = (*o_nEdits <= scoreLimit) ? *o_nEdits : -1;
            *o_nEdits += *o_patternOffset; // Add gaps at the end
            *matchProbability *= lv_indelProbabilities[*o_patternOffset];
            return maxScore;
        }
        else { // We could not find a longer gapless alignment beyond seedLen
            *o_nEdits = ScoreAboveLimit;
            *o_nEditsGapless = ScoreAboveLimit;
            return ScoreAboveLimit;
        }
    }
    
    int computeScoreBanded(
        const char* text,
        int textLen,
        const char* pattern,
        const char *qualityString,
        int patternLen,
        int w,
        int scoreInit,
        bool isRC,
        int *o_textOffset = NULL,
        int *o_patternOffset = NULL,
        int *o_nEdits = NULL,
        double *matchProbability = NULL,
        bool* o_preferIndel = NULL)
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

        if (w < 0) {
            *o_nEdits = ScoreAboveLimit;
            return -1;
        }

        bool localPreferIndel;
        if (NULL == o_preferIndel) {
            o_preferIndel = &localPreferIndel;
        }

        //
        // Start with perfect match probability and work our way down.
        //
        *matchProbability = 1.0;

        if (TEXT_DIRECTION == -1) {
            text--; // so now it points at the "first" character of t, not after it.
        }
         
        int bandWidth = (2 * w + 1) < patternLen ? (2 * w + 1) : patternLen;
        int numVec = (bandWidth + VEC_SIZE - 1) / VEC_SIZE; // Number of vectors per Farrar segment
        int segLen = numVec * VEC_SIZE;
        int numSeg = (patternLen + segLen - 1) / segLen;
        int patternIdx = 0;

        // 
        // Generate query profile
        //
        int16_t* queryResult = (int16_t*)qProfile;
        for (int a = 0; a < MAX_ALPHABET_SIZE; a++) {
            for (int i = 0; i < numSeg; i++) {
                for (int j = 0; j < numVec; j++) {
                    for (int k = j; k < segLen; k += numVec) {
                        int idx = i * segLen + k;
                        if (idx < patternLen) {
                            uint8_t bp = BASE_VALUE[pattern[idx]];
                            queryResult[patternIdx] = ntTransitionMatrix[a * MAX_ALPHABET_SIZE + bp];
                        }
                        else {
                            queryResult[patternIdx] = INT16_MIN;
                        }
                        patternIdx++;
                    }
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
        __m128i v_mask = _mm_cmpgt_epi16(_mm_set_epi16(0, 0, 0, 0, 0, 0, 0, 1), v_zero);

        int endBonus;
        if (!isRC) {
            if (TEXT_DIRECTION == -1) {
                endBonus = fivePrimeEndBonus;
            }
            else {
                endBonus = threePrimeEndBonus;
            }
        }
        else {
            if (TEXT_DIRECTION == -1) {
                endBonus = threePrimeEndBonus;
            }
            else {
                endBonus = fivePrimeEndBonus;
            }
        }

        //
        // Initialize scores of first row
        //
        uint16_t scoreFirstRow[VEC_SIZE] = {};
        for (int segIdx = 0; segIdx < numSeg; segIdx++) {
            for (int vecIdx = 0; vecIdx < numVec; vecIdx++) {
                for (int elemIdx = 0; elemIdx < VEC_SIZE; elemIdx++) {
                    int patternIdx = (segIdx * segLen) + (elemIdx * numVec) + vecIdx;
                    if (patternIdx < patternLen) {
                        scoreFirstRow[elemIdx] = __max(0, scoreInit - gapOpenPenalty - patternIdx * gapExtendPenalty);
                    }
                }
                _ASSERT(VEC_SIZE == 8);  // FIXME: Initialization below works only when VEC_SIZE = 8
                _mm_store_si128(H + (segIdx * numVec) + vecIdx, _mm_setr_epi16(scoreFirstRow[0], scoreFirstRow[1], scoreFirstRow[2], scoreFirstRow[3],
                                                                scoreFirstRow[4], scoreFirstRow[5], scoreFirstRow[6], scoreFirstRow[7]));
                _mm_store_si128(Hminus1 + (segIdx * numVec) + vecIdx, v_zero);
                _mm_store_si128(E + (segIdx * numVec) + vecIdx, v_zero);
            }
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
            __m128i* qRowProfile = qProfile + BASE_VALUE[*t] * numSeg * numVec;

            // Registers to hold intermediate scores for each row
            __m128i m, h, temp, e = v_zero, f = v_zero, max = v_zero, X = v_zero;

            int maxScoreRow = 0;
            int localAlignmentPatternOffset = -1;

            int bandBeg = __max(i - w, 0);
            int bandEnd = __min(i + w, patternLen -1); 
            int segBeg = bandBeg / segLen;
            int segEnd = bandEnd / segLen;

            for (int j = segBeg; j <= segEnd; j++) {

                // Load h from the previous row
                h = _mm_load_si128(Hptr + (j * numVec) + numVec - 1);
                
                // Shift left h and blend in initial values
                h = _mm_slli_si128(h, 2); // shift h left by 2 * 8 bits

                __m128i v_hInit;
                if (j == 0) {
                    int hInit = scoreInit;
                    if (i > 0) {
                        hInit = __max(0, scoreInit - gapOpenPenalty - (i - 1) * gapExtendPenalty);
                    }
                    v_hInit = _mm_set1_epi16(hInit);
                }
                else {
                    if (bandBeg > j * segLen) {
                        v_hInit = v_zero;
                    }
                    else {
                        v_hInit = _mm_load_si128(Hptr + j * numVec - 1);
                        v_hInit = _mm_srli_si128(v_hInit, 14);
                    }
                }
                
                // h = _mm_blend_epi16(h, v_hInit, 1); // Does not work with SSE !
                h = blend_sse(h, v_hInit, v_mask);

                __m128i* backtraceActionSeg = backtraceAction + (i * numVec * numSeg) + (j * numVec);

                for (int k = 0; (k < numVec) && (j * segLen + k) <= bandEnd; k++) {

                    __m128i backtraceActionVec;

                    __m128i mask = _mm_cmpgt_epi16(h, v_zero); // h > 0

                    __m128i qProfileVec = _mm_load_si128(qRowProfile + j * numVec + k);
                    m = _mm_adds_epi16(h, qProfileVec);
                    m = _mm_and_si128(m, mask);

                    e = _mm_load_si128(E + j * numVec + k);

                    // h = max{m, e, f}
                    backtraceActionVec = _mm_and_si128(_mm_cmpgt_epi16(e, m), v_one); // action = e > m ? 1 : 0
                    h = _mm_max_epi16(m, e);
                    __m128i tmpResult = _mm_and_si128(_mm_cmpgt_epi16(f, h), v_two); 
                    backtraceActionVec = _mm_or_si128(tmpResult, _mm_andnot_si128(tmpResult, backtraceActionVec)); // action = f > h ? 2 : action 
                    h = _mm_max_epi16(h, f);

                    max = _mm_max_epi16(max, h);

                    // Store h for the next row
                    _mm_store_si128(Hminus1ptr + j * numVec + k, h);

                    // e = max{m - gapOpen, e - gapExtend}
                    e = _mm_subs_epi16(e, v_gapExtend);
                    temp = _mm_subs_epi16(m, v_gapOpen);
                    temp = _mm_max_epi16(temp, v_zero);
                    tmpResult = _mm_and_si128(_mm_cmpgt_epi16(e, temp), v_four);
                    backtraceActionVec = _mm_or_si128(backtraceActionVec, tmpResult);
                    e = _mm_max_epi16(e, temp);
                    _mm_store_si128(E + j * numVec + k, e);

                    // f = max{m - gapOpen, f- gapExtend}
                    f = _mm_subs_epi16(f, v_gapExtend);
                    tmpResult = _mm_and_si128(_mm_cmpgt_epi16(f, temp), v_thirtytwo);
                    backtraceActionVec = _mm_or_si128(backtraceActionVec, tmpResult);
                    f = _mm_max_epi16(f, temp);

                    // Store traceback information
                    _mm_store_si128(backtraceActionSeg + k, backtraceActionVec);

                    // Load the next score vector
                    h = _mm_load_si128(Hptr + j * numVec + k);

                } // end segment

                //
                // Farrar's algorithm does lazy f evaluation. Since f rarely influences final score h, the algorithm speculates f = zero initially 
                // Re-evaluate if f could influence h after we have made a first pass through the row
                //
                for (int k = 0; k < (VEC_SIZE - 1); k++) {
                    
                    // Save f value at boundary for next segment
                    X = _mm_max_epi16(X, _mm_srli_si128(f, 14));
                    
                    // Extract the last f vector
                    f = _mm_slli_si128(f, 2);

                    for (int v = 0; (v < numVec) && (j * segLen + v) <= bandEnd; v++) {

                         h = _mm_load_si128(Hminus1ptr + j * numVec + v);

                        // action = f > h ? 2 : action
                        __m128i tmpResult = _mm_and_si128(_mm_cmpgt_epi16(f, h), v_two);
                        __m128i backtraceActionVec = _mm_load_si128(backtraceActionSeg + v);
                        __m128i tmpResult2 = _mm_andnot_si128(tmpResult, backtraceActionVec);
                        backtraceActionVec = _mm_or_si128(tmpResult, tmpResult2); 

                        h = _mm_max_epi16(h, f);
                        _mm_store_si128(Hminus1ptr + j * numVec + v, h);

                        max = _mm_max_epi16(max, h);

                        temp = _mm_subs_epu16(h, v_gapOpen);
                        f = _mm_subs_epu16(f, v_gapExtend);

                        tmpResult = _mm_and_si128(_mm_cmpgt_epi16(f, temp), v_thirtytwo);
                        backtraceActionVec = _mm_or_si128(backtraceActionVec, tmpResult);
                        _mm_store_si128(backtraceActionSeg + v, backtraceActionVec);

                        // Converged if no element of f can influence h
                        bool converged = !(_mm_movemask_epi8(_mm_cmpgt_epi16(f, temp))); // !! converts int to bool

                        if (converged) goto got_answer;
                    }
                }
got_answer:
                // Pass on f for next segment
                f = X;

            } // end pattern

            maxScoreRow = getMax(max);
#ifdef TRACE_AG
            __m128i rowScore;
            for (int j = 0; j < segBeg; j++) {
                for (int k = 0; k < numVec; k++) {
                    printElem(v_zero);
                }
            }
            for (int j = segBeg; j <= segEnd; j++) {
                for (int k = 0; k < numVec && (j * segLen + k) <= bandEnd; k++) {
                    rowScore = _mm_load_si128(Hminus1ptr + j * numVec + k);
                    printElem(rowScore);   
                }
            }
            printf("\n");
#endif
            
            if (bandEnd == patternLen - 1) {
                //
                // Global alignment score (i.e., score when aligning to the end of the pattern)
                // 
                int vecIdx = ((bandEnd / segLen) * numVec) + ((bandEnd % segLen) % numVec);
                int elemIdx = (bandEnd % segLen) / numVec;
                __m128i v_globalAlignmentScore = _mm_load_si128(Hminus1ptr + vecIdx);
                int globalAlignmentScore = getElem(elemIdx, v_globalAlignmentScore);
                if (globalAlignmentScore >= bestGlobalAlignmentScore) {
                    bestGlobalAlignmentScore = globalAlignmentScore;
                    bestGlobalAlignmentTextOffset = i;
                }
                
            }

            if (maxScoreRow == 0) break;

            if (maxScoreRow > bestLocalAlignmentScore) { // If we obtained a better score this round  
                // 
                // Get index in pattern where maximum score was obtained
                // 
                for (int j = segBeg; j <= segEnd; j++) {
                    for (int k = 0; (k < numVec) && (j * segLen + k) <= bandEnd; k++) {
                        __m128i h = _mm_load_si128(Hminus1ptr + j * numVec + k); // Load a vector of score values 
                        __m128i mask = _mm_cmpeq_epi16(h, _mm_set1_epi16(maxScoreRow)); // Create a mask of 1's in cells where score = maxScoreRow
                        int result = _mm_movemask_epi8(_mm_packs_epi16(mask, v_zero)); // Convert vector result to 8-bit
                        unsigned long zeroes;
                        if (result) {
                            CountLeadingZeroes(result, zeroes);
                            int patternOffset = j * segLen + zeroes * numVec + k;
                            localAlignmentPatternOffset = (localAlignmentPatternOffset > patternOffset) ? localAlignmentPatternOffset : patternOffset;
                            _ASSERT(localAlignmentPatternOffset < patternLen); // Scores outside the pattern boundaries should never be the maximum in the row
                        }
                    }
                }

                bestLocalAlignmentScore = maxScoreRow;
                bestLocalAlignmentTextOffset = i;
                bestLocalAlignmentPatternOffset = localAlignmentPatternOffset;
                
            }
            
            // Swap roles of H and Hminus1 for the next row
            __m128i* hTemp = Hminus1ptr; 
            Hminus1ptr = Hptr; 
            Hptr = hTemp;

        } // end text

        // Choose between local and global alignment for patternOffset
        if ((bestLocalAlignmentScore != bestGlobalAlignmentScore) && (bestLocalAlignmentScore >= bestGlobalAlignmentScore + endBonus)) { // FIXME: Change 5 to a clipping penalty
            // Local alignment preferred
            *o_patternOffset = bestLocalAlignmentPatternOffset;
            *o_textOffset = bestLocalAlignmentTextOffset;
            score = bestLocalAlignmentScore;

            //
            // Check if we can match more bases near the end of the soft-clipped read by deleting a character from the reference.
            // This helps reduce some INDEL false negatives introduced by using a high gap open penalty
            //
            int patternOffsetAdj = *o_patternOffset - 1;
            int textOffsetAdj = *o_textOffset;
            int countEndMatches = 0;

            while ((patternOffsetAdj + 1 != patternLen) && pattern[patternOffsetAdj + 1] == *(text + (textOffsetAdj + 1) * TEXT_DIRECTION)) {
                countEndMatches++;
                patternOffsetAdj++;
                textOffsetAdj++;
            }

            if (countEndMatches >= 3) { // matched too few bases, don't use the alignment with deletion. FIXME: check if 3 is enough!
                *o_patternOffset = patternOffsetAdj;
                *o_textOffset = textOffsetAdj;
            }
            else {
                //
                // Check if we can match more bases near the end of the soft-clipped read by inserting a character to the pattern.
                // This helps reduce some INDEL false negatives introduced by using a high gap open penalty
                //
                patternOffsetAdj = *o_patternOffset + 1;
                textOffsetAdj = *o_textOffset;
                countEndMatches = 0;

                while ((patternOffsetAdj < patternLen) && pattern[patternOffsetAdj] == *(text + (textOffsetAdj)*TEXT_DIRECTION)) {
                    countEndMatches++;
                    patternOffsetAdj++;
                    textOffsetAdj++;
                }

                if (countEndMatches >= 3) { // matched too few bases, don't use the alignment with insertion. FIXME: check if 3 is enough!
                    *o_patternOffset = patternOffsetAdj - 1;
                    *o_textOffset = textOffsetAdj - 1;
                }
                else if (countEndMatches >= 1) {
                    *o_preferIndel = true;
                    *o_patternOffset = patternOffsetAdj - 1;
                    *o_textOffset = textOffsetAdj - 1;
                }
            }

            if (*o_patternOffset == bestLocalAlignmentPatternOffset && *o_textOffset == bestLocalAlignmentTextOffset) {
                patternOffsetAdj = *o_patternOffset;
                textOffsetAdj = *o_textOffset;
                //
                // Try not to clip high quality bases (>= 65) from the read. These will be reported as insertions in the final alignment
                //
                while (patternOffsetAdj != patternLen - 1 && qualityString[patternOffsetAdj] >= 65 && qualityString[patternOffsetAdj + 1] >= 65) {
                    patternOffsetAdj += 1;
                }

                if (patternOffsetAdj == patternLen - 1) {
                    *o_patternOffset = patternOffsetAdj;
                }
                else if (patternOffsetAdj >= *o_patternOffset + 2) { // a 1bp sub is preferred over a 1bp ins with the current scoring scheme
                    int tmpOffset = patternOffsetAdj + 1;
                    int countRemHighQualityBases = 0;
                    int remPatternLen = patternLen - tmpOffset;

                    while (tmpOffset != patternLen - 1) {
                        if (qualityString[tmpOffset] >= 65) {
                            countRemHighQualityBases++;
                        }
                        tmpOffset++;
                    }

                    _ASSERT(remPatternLen != 0);
                    if (((float)countRemHighQualityBases) / remPatternLen < 0.1) {
                        *o_patternOffset = patternOffsetAdj;
                    }
                }
            }
        }
        else {
            // Global alignment preferred
            *o_patternOffset = patternLen - 1;
            *o_textOffset = bestGlobalAlignmentTextOffset;
            score = bestGlobalAlignmentScore;
        }

        if (score > scoreInit) { 
            // Perform traceback and compute nEdits
            int rowIdx = *o_textOffset, colIdx = *o_patternOffset, matrixIdx;
            BacktraceActionType action = M, prevAction = M;
            int actionCount = 1;
            int nMatches = 0, nMismatches = 0, nGaps = 0;

            // Start traceback from the cell (i,j) with the maximum score
            while (rowIdx >= 0 && colIdx >= 0) {
                int vecIdx = (colIdx / segLen) * numVec + ((colIdx % segLen) % numVec);
                int elemIdx = (colIdx % segLen) / numVec;
                uint16_t* backtracePointersVec = (uint16_t*)(backtraceAction + (rowIdx * numVec * numSeg) + vecIdx);
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
                action = (BacktraceActionType) ((backtracePointersVec[elemIdx] >> matrixIdx) & 3);
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

            // *o_nEdits = (*o_nEdits <= w) ? *o_nEdits : -1; // return -1 if we have more edits than threshold w

            *matchProbability *= lv_perfectMatchProbability[nMatches];

            *o_textOffset += 1;
            *o_patternOffset += 1;

            *o_textOffset = patternLen - *o_textOffset;
            *o_patternOffset = patternLen - *o_patternOffset;

            *matchProbability *= lv_indelProbabilities[*o_patternOffset];
            
            return score;
        }
        else {
            return -1; 
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
        bool isRC,
        int *o_textOffset = NULL,
        int *o_patternOffset = NULL,
        int *o_nEdits = NULL,
        double *matchProbability = NULL,
        bool *o_preferIndel = NULL)
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

        bool localPreferIndel;
        if (NULL == o_preferIndel) {
            o_preferIndel = &localPreferIndel;
        }

        if (w < 0) {
            *o_nEdits = ScoreAboveLimit;
            return -1;
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
        __m128i v_mask = _mm_cmpgt_epi16(_mm_set_epi16(0, 0, 0, 0, 0, 0, 0, 1), v_zero);

        int endBonus;
        if (!isRC) {
            if (TEXT_DIRECTION == -1) {
                endBonus = fivePrimeEndBonus;
            }
            else {
                endBonus = threePrimeEndBonus;
            }
        }
        else {
            if (TEXT_DIRECTION == -1) {
                endBonus = threePrimeEndBonus;
            }
            else {
                endBonus = fivePrimeEndBonus;
            }
        }

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
            int localAlignmentPatternOffset = -1;

            // Load h from the previous row
            h = _mm_load_si128(Hptr + numVec - 1);

            // Shift left h and blend in initial values
            h = _mm_slli_si128(h, 2); // shift h left by 2 * 8 bits
            int hInit = scoreInit;
            if (i > 0) {
                hInit = __max(0, scoreInit - gapOpenPenalty - (i - 1) * gapExtendPenalty);
            }
            __m128i v_hInit = _mm_set1_epi16(hInit);
            // h = _mm_blend_epi16(h, v_hInit, 1); // Does not work with SSE !
            h = blend_sse(h, v_hInit, v_mask);

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
                        localAlignmentPatternOffset = (localAlignmentPatternOffset > patternOffset) ? localAlignmentPatternOffset : patternOffset;
                        _ASSERT(localAlignmentPatternOffset < patternLen); // Scores outside the pattern boundaries should never be the maximum in the row
                    }
                }
                bestLocalAlignmentScore = maxScoreRow;
                bestLocalAlignmentTextOffset = i;
                bestLocalAlignmentPatternOffset = localAlignmentPatternOffset;
            }
            
            // Swap roles of H and Hminus1 for the next row
            __m128i* hTemp = Hminus1ptr; 
            Hminus1ptr = Hptr; 
            Hptr = hTemp;

        } // end text

        // Choose between local and global alignment for patternOffset
        if ((bestLocalAlignmentScore != bestGlobalAlignmentScore) && (bestLocalAlignmentScore >= bestGlobalAlignmentScore + endBonus)) { // FIXME: Change 5 to a clipping penalty
            // Local alignment preferred
            *o_patternOffset = bestLocalAlignmentPatternOffset;
            *o_textOffset = bestLocalAlignmentTextOffset;
            score = bestLocalAlignmentScore;

            //
            // Check if we can match more bases near the end of the soft-clipped read by deleting a character from the reference.
            // This helps reduce some INDEL false negatives introduced by using a high gap open penalty
            //
            int patternOffsetAdj = *o_patternOffset - 1;
            int textOffsetAdj = *o_textOffset;
            int countEndMatches = 0;

            while ((patternOffsetAdj + 1 != patternLen) && pattern[patternOffsetAdj + 1] == *(text + (textOffsetAdj + 1) * TEXT_DIRECTION)) {
                countEndMatches++;
                patternOffsetAdj++;
                textOffsetAdj++;
            }

            if (countEndMatches >= 3) { // matched too few bases, don't use the alignment with deletion. FIXME: check if 3 is enough!
                *o_patternOffset = patternOffsetAdj;
                *o_textOffset = textOffsetAdj;
            }
            else {
                //
                // Check if we can match more bases near the end of the soft-clipped read by inserting a character to the pattern.
                // This helps reduce some INDEL false negatives introduced by using a high gap open penalty
                //
                patternOffsetAdj = *o_patternOffset + 1;
                textOffsetAdj = *o_textOffset;
                countEndMatches = 0;

                while ((patternOffsetAdj < patternLen) && pattern[patternOffsetAdj] == *(text + (textOffsetAdj)*TEXT_DIRECTION)) {
                    countEndMatches++;
                    patternOffsetAdj++;
                    textOffsetAdj++;
                }

                if (countEndMatches >= 3) { // matched too few bases, don't use the alignment with insertion. FIXME: check if 3 is enough!
                    *o_patternOffset = patternOffsetAdj - 1;
                    *o_textOffset = textOffsetAdj - 1;
                }
                else if (countEndMatches >= 1) {
                    *o_preferIndel = true;
                    *o_patternOffset = patternOffsetAdj - 1;
                    *o_textOffset = textOffsetAdj - 1;
                }
            }

            if (*o_patternOffset == bestLocalAlignmentPatternOffset && *o_textOffset == bestLocalAlignmentTextOffset) {
                patternOffsetAdj = *o_patternOffset;
                textOffsetAdj = *o_textOffset;
                //
                // Try not to clip high quality bases (>= 65) from the read. These will be reported as insertions in the final alignment
                //
                while (patternOffsetAdj != patternLen - 1 && qualityString[patternOffsetAdj] >= 65 && qualityString[patternOffsetAdj + 1] >= 65) {
                    patternOffsetAdj += 1;
                }

                if (patternOffsetAdj == patternLen - 1) {
                    *o_patternOffset = patternOffsetAdj;
                }
                else if (patternOffsetAdj >= *o_patternOffset + 2) {
                    int tmpOffset = patternOffsetAdj + 1;
                    int countRemHighQualityBases = 0;
                    int remPatternLen = patternLen - tmpOffset;

                    while (tmpOffset != patternLen - 1) {
                        if (qualityString[tmpOffset] >= 65) {
                            countRemHighQualityBases++;
                        }
                        tmpOffset++;
                    }

                    _ASSERT(remPatternLen != 0);
                    if (((float)countRemHighQualityBases) / remPatternLen < 0.1) {
                        *o_patternOffset = patternOffsetAdj;
                    }
                }
            }
        }
        else {
            // Global alignment preferred
            *o_patternOffset = patternLen - 1;
            *o_textOffset = bestGlobalAlignmentTextOffset;
            score = bestGlobalAlignmentScore;
        }

        if (score > scoreInit) {
            // Perform traceback and compute nEdits
            int rowIdx = *o_textOffset, colIdx = *o_patternOffset, matrixIdx;
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

            // *o_nEdits = (*o_nEdits <= w) ? *o_nEdits : -1; // return -1 if we have more edits than threshold w

            *matchProbability *= lv_perfectMatchProbability[nMatches];

            *o_textOffset += 1;
            *o_patternOffset += 1;

            *o_textOffset = patternLen - *o_textOffset;
            *o_patternOffset = patternLen - *o_patternOffset;

            *matchProbability *= lv_indelProbabilities[*o_patternOffset];
            
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
    int fivePrimeEndBonus;
    int threePrimeEndBonus;

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
};

class AffineGapVectorizedWithCigar {
public:
    AffineGapVectorizedWithCigar(); // FIXME: Pass scoring parameters

    AffineGapVectorizedWithCigar(int i_matchReward, int i_subPenalty, int i_gapOpenPenalty, int i_gapExtendPenalty);

    // Compute the affine gap score between two strings and write the CIGAR string in cigarBuf.
    // Returns -1 if the edit distance exceeds k or -2 if we run out of space in cigarBuf.
    int computeGlobalScore(const char* text, int textLen, const char* pattern, const char* quality, int patternLen, int w,
        char* cigarBuf, int cigarBufLen, bool useM,
        CigarFormat format = COMPACT_CIGAR_STRING,
        int* o_cigarBufUsed = NULL, int *o_netDel = NULL, int *o_tailIns = NULL);

    int computeGlobalScoreBanded(const char* text, int textLen, const char* pattern, const char* quality, int patternLen, int w, int scoreInit,
        char* cigarBuf, int cigarBufLen, bool useM,
        CigarFormat format = COMPACT_CIGAR_STRING,
        int* o_cigarBufUsed = NULL, int* o_netDel = NULL, int* o_tailIns = NULL);

    int computeGlobalScoreNormalized(const char* text, int textLen,
        const char* pattern, const char* quality, int patternLen,
        int k,
        char *cigarBuf, int cigarBufLen, bool useM,
        CigarFormat format, int* o_cigarBufUsed,
        int* o_addFrontClipping, int *o_netDel = NULL, int* o_tailIns = NULL);

    bool writeCigar(char** o_buf, int* o_buflen, int count, char code, CigarFormat format);

    void setGapOpenPenalty(int gapOpenPenalty_) {
        gapOpenPenalty = gapOpenPenalty_;
    }

    void setGapExtendPenalty(int gapExtendPenalty_) {
        gapExtendPenalty = gapExtendPenalty_;
    }

    int getGapOpenPenalty() {
        return gapOpenPenalty;
    }

    int getGapExtendPenalty() {
        return gapExtendPenalty;
    }

    int getSubPenalty() {
        return subPenalty;
    }

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
    int defaultGapOpenPenalty;
    int gapOpenPenalty;
    int gapExtendPenalty;
    int minScoreParam;
    int maxScoreParam;

    //
    // Precompute query profile which is a table containing the result of
    // matching each letter of the alphabet with each character of the pattern
    // and filling in the appropriate transition score from the ntTransMatrix.
    // This precomputation saves looking up the ntTransition matrix for each
    // (ref_i, query_j) in the inner dynamic programming loop.
    //
    __m128i qProfile[MAX_ALPHABET_SIZE * MAX_VEC_SEGMENTS];

    // 
    // H and E arrays used for storing scores in every row
    //
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
    // Structure used for storing (action, count) pairs from backtracking
    //
    typedef struct {
        BacktraceActionType action[MAX_READ_LENGTH];
        int count[MAX_READ_LENGTH];
    } LocalCigarResult;

    LocalCigarResult res;

};
