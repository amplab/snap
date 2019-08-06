#pragma once

#include "Compat.h"
#include "FixedSizeMap.h"
#include "BigAlloc.h"
#include "exit.h"
#include "Genome.h"
#include "Read.h"
#include "Tables.h"
#include "LandauVishkin.h"

const int MAX_ALPHABET_SIZE = 5;

//
// These are global so there are only one for both senses of the template
//

extern double *lv_indelProbabilities;  // Maps indels by length to probability of occurance.
extern double *lv_phredToProbability;  // Maps ASCII phred character to probability of error, including 
extern double *lv_perfectMatchProbability; // Probability that a read of this length has no mutations

// #define PRINT_SCORES

// 
// Computes the affine gap score between two strings
//
template<int TEXT_DIRECTION = 1> class AffineGap {

public:
    AffineGap()
{
        init(1, 4, 6, 1);
}

    AffineGap(
        int i_matchReward,
        int i_subPenalty,
        int i_gapOpenPenalty,
        int i_gapExtendPenalty) 
{
        init(i_matchReward, i_subPenalty, i_gapOpenPenalty, i_gapExtendPenalty);
}

    static size_t getBigAllocatorReservation() { return sizeof(AffineGap<TEXT_DIRECTION>); }

    ~AffineGap() 
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

    // Initialize score arrays
    memsetint(H, 0, MAX_READ_LENGTH + 1);
    memsetint(E, 0, MAX_READ_LENGTH + 1);

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

#ifdef PRINT_SCORES
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
        // Generate query profile
        //
        for (int i = 0, k = 0; i < MAX_ALPHABET_SIZE; i++) {
            for (int j = 0; j < patternLen; j++) {
                _uint8 bp = BASE_VALUE[pattern[j]]; // 2b encoded base pair
                qProfile[k++] = ntTransitionMatrix[i * MAX_ALPHABET_SIZE + bp];
            }
        }

        // Initialize score arrays
        memsetint(H, 0, patternLen + 1);
        memsetint(E, 0, patternLen + 1);

        // 
        // Initialize H array
        //
        H[0] = scoreInit;
        H[1] = scoreInit > gapOpenPenalty ? scoreInit - gapOpenPenalty : 0;
        for (int j = 2; j < patternLen && H[j - 1] > gapExtendPenalty; ++j) {
            H[j] = H[j - 1] - gapExtendPenalty;
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

        int beg = 0, end = patternLen;
        // Iterate over all rows of text
        for (int i = 0; i < textLen; i++) {
            // Compute only a 2w+1 band
            if (beg < i - w) beg = i - w;
            if (end > i + w + 1) end = i + w + 1;
            if (end > patternLen) end = patternLen;

            const char* t = (text + i * TEXT_DIRECTION);

            // Get the query profile for the row
            _int8* qRowProfile = &qProfile[BASE_VALUE[*t] * patternLen];

            _mm_prefetch((const char*)&qProfile[BASE_VALUE[*(t + TEXT_DIRECTION)] * patternLen], _MM_HINT_T0);

            int hLeft = 0, fCurr = 0;

            // Initialize first column
            if (beg == 0) {
                hLeft = scoreInit - (gapOpenPenalty + i * gapExtendPenalty);
            }

            //  H    -1      0       j-1     j       .. (pattern)
            // -1     
            //  0        
            //  1
            //  i-1                hDiag
            //  i                  hLeft  hCurr
            //  ...
            //
            //  (text)

            int maxScoreRow = 0;

#ifdef PRINT_SCORES
            for (int j = 0; j < beg; ++j) {
                printf("0,");
            }
#endif

            // Iterate over all columns of pattern (within the band)
            for (int j = beg; j < end; ++j) {

                // 
                // Compute H(i,j).
                // H(i,j) = max(H(i-1,j-1) + qRowProfile[j], E(i,j), F(i,j))
                //
                int hDiag = H[j];
                hDiag = hDiag > 0 ? hDiag + qRowProfile[j] : 0;
                char action = (hDiag >= E[j]) ? 'M' : 'D';
                int hCurr = __max(hDiag, E[j]);
                action = (hCurr >= fCurr) ? action : 'I';
                hCurr = __max(hCurr, fCurr);
                backtraceAction[i][j][0] = action;

                bestLocalAlignmentPatternOffset = (hCurr >= maxScoreRow) ? j : bestLocalAlignmentPatternOffset;
                maxScoreRow = __max(maxScoreRow, hCurr);                

#ifdef PRINT_SCORES
                printf("%d,", hCurr);
#endif
                // hLeft will be used as hDiag for the next row i+1
                H[j] = hLeft;

                // hCurr will be used as hLeft for the next column j+1
                hLeft = hCurr;

                //
                // Compute E(i+1,j) and F(i,j+1)
                // E(i+1,j) = max(H(i,j) - gapOpen, E(i,j) - gapExtend)
                // F(i,j+1) = max(H(i,j) - gapOpen, F(i,j) - gapExtend)
                // Note hDiag is used instead of hCurr below for H(i,j) to disallow insertions followed by deletions (e.g., 1I1D) in CIGAR
                //
                int hTemp = __max(hDiag - gapOpenPenalty, 0);
                E[j] -= gapExtendPenalty;
                action = (E[j] > hTemp) ? 'D' : 'M';
                E[j] = __max(hTemp, E[j]);
                backtraceAction[i][j][1] = action;

                fCurr -= gapExtendPenalty;
                action = (fCurr > hTemp) ? 'I' : 'M';
                fCurr = __max(hTemp, fCurr);
                backtraceAction[i][j][2] = action;

            } // end pattern
            H[end] = hLeft;
            E[end] = 0;

#ifdef PRINT_SCORES
            printf("\n");
#endif

            if (end == patternLen) { // we got a semi-global alignment
                int globalAlignmentScore = H[end];
                if (globalAlignmentScore >= bestGlobalAlignmentScore) {
                    bestGlobalAlignmentScore = globalAlignmentScore;
                    bestGlobalAlignmentTextOffset = i;
                }
            }

            // All scores in row are zero. Break out and report local alignment
            if (maxScoreRow == 0) {
                break;
            }

            if (maxScoreRow > bestLocalAlignmentScore) { // If we obtained a better score this round  
                bestLocalAlignmentScore = maxScoreRow;
                bestLocalAlignmentTextOffset = i;
            }

            // BWA-MEM's band narrowing heuristic
            int j;
            for (j = beg; j < end && H[j] == 0 && E[j] == 0; ++j);
            beg = j;
            for (j = end; j >= beg && H[j] == 0 && E[j] == 0; --j);
            end = j + 2; 

        } // end text

        // Choose between local and global alignment for patternOffset
        if ((bestLocalAlignmentScore != bestGlobalAlignmentScore) && (bestLocalAlignmentScore >= bestGlobalAlignmentScore + 5)) { // FIXME: Change 5 to a clipping penalty
            // Local alignment preferred
            *o_patternOffset = bestLocalAlignmentPatternOffset;
            *o_textOffset = bestLocalAlignmentTextOffset;
            score = bestLocalAlignmentScore;
        }
        else {
            // Global alignment preferred
            *o_patternOffset = patternLen - 1;
            *o_textOffset = bestGlobalAlignmentTextOffset;
            score = bestGlobalAlignmentScore;
        }

        if (bestGlobalAlignmentScore > scoreInit) {
            int rowIdx = bestGlobalAlignmentTextOffset, colIdx = patternLen - 1, matrixIdx = 0;
            if (abs(rowIdx - colIdx) > w) return -1;
            char action = 'M', prevAction = 'M';
            int actionCount = 1; 
            int nMatches = 0, nMismatches = 0, nGaps = 0;
 
            // Start traceback from the cell (i,j) with the maximum score
            while (rowIdx >= 0 && colIdx >= 0) {
                action = backtraceAction[rowIdx][colIdx][matrixIdx];
                _ASSERT(action == 'M' || action == 'I' || action == 'D');
                if (action == 'M') {
                    if (pattern[colIdx] != text[rowIdx * TEXT_DIRECTION]) {
                        // Compute probabilties of mismatches
                        *matchProbability *= lv_phredToProbability[qualityString[colIdx]];
                        nMismatches++;
                    }
                    else {
                        nMatches++;
                    }
                    rowIdx--; colIdx--; matrixIdx = 0;
                } else if (action == 'D') {
                    rowIdx--; matrixIdx = 1;
                } else if (action == 'I') {
                    colIdx--; matrixIdx = 2;
                }
                // Ignore transitions from H. Compute number of indels in sequence for matchProbability
                if (prevAction != 'M') {
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
                /*
                if (nMismatches + nGaps > w) {
                    return -1;
                }
                */
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

void *operator new(size_t size, BigAllocator *allocator) { _ASSERT(size == sizeof(AffineGap<TEXT_DIRECTION>)); return allocator->allocate(size); }
void operator delete(void *ptr, BigAllocator *allocator) {/*Do nothing.  The memory is freed when the allocator is deleted.*/ }

private:

    //
    // Scores for each nucleotide <-> nucleotide transition
    //
    int ntTransitionMatrix[MAX_ALPHABET_SIZE * MAX_ALPHABET_SIZE];

    //
    // Precompute query profile which is a table containing the result of
    // matching each letter of the alphabet with each character of the pattern
    // and filling in the appropriate transition score from the ntTransMatrix.
    // This precomputation saves looking up the ntTransition matrix for each
    // (ref_i, query_j) in the inner dynamic programming loop.
    //
    _int8 qProfile[MAX_ALPHABET_SIZE * MAX_READ_LENGTH];

    // 
    // H and E arrays used for storing scores in every row
    //
    int H[MAX_READ_LENGTH + 1];
    int E[MAX_READ_LENGTH + 1];

    //
    // Pointers to traceback alignment, one for each of the three affine-gap matrices, H, E and F
    //
    char backtraceAction[(MAX_READ_LENGTH + MAX_K)][MAX_READ_LENGTH][3];

    //
    // Affine gap scoring parameters
    //
    int matchReward;
    int subPenalty;
    int gapOpenPenalty;
    int gapExtendPenalty;

};

class AffineGapWithCigar {
public:
    AffineGapWithCigar(); // FIXME: Pass scoring parameters

    AffineGapWithCigar(int i_matchReward, int i_subPenalty, int i_gapOpenPenalty, int i_gapExtendPenalty);

    // Compute the affine gap score between two strings and write the CIGAR string in cigarBuf.
    // Returns -1 if the edit distance exceeds k or -2 if we run out of space in cigarBuf.
    int computeGlobalScore(const char* text, int textLen, const char* pattern, int patternLen, int w,
        char* cigarBuf, int cigarBufLen, bool useM,
        CigarFormat format = COMPACT_CIGAR_STRING,
        int* o_cigarBufUsed = NULL, int *o_netDel = NULL, int *o_tailIns = NULL);

    int computeGlobalScoreNormalized(const char* text, int textLen,
        const char* pattern, int patternLen,
        int k,
        char *cigarBuf, int cigarBufLen, bool useM,
        CigarFormat format, int* o_cigarBufUsed,
        int* o_addFrontClipping, int *o_netDel = NULL, int* o_tailIns = NULL);

    bool writeCigar(char** o_buf, int* o_buflen, int count, char code, CigarFormat format);

private:

    //
    // Scores for each nucleotide <-> nucleotide transition
    //
    int ntTransitionMatrix[MAX_ALPHABET_SIZE * MAX_ALPHABET_SIZE];

    //
    // Precompute query profile which is a table containing the result of
    // matching each letter of the alphabet with each character of the pattern
    // and filling in the appropriate transition score from the ntTransMatrix.
    // This precomputation saves looking up the ntTransition matrix for each
    // (ref_i, query_j) in the inner dynamic programming loop.
    //
    _int8 qProfile[MAX_ALPHABET_SIZE * MAX_READ_LENGTH];

    // 
    // H and E arrays used for storing scores in every row
    //
    int H[MAX_READ_LENGTH + 1];
    int E[MAX_READ_LENGTH + 1];

    //
    // Pointers to traceback alignment, one for each of the three affine-gap matrices, H, E and F
    //S
    char backtraceAction[(MAX_READ_LENGTH + MAX_K)][MAX_READ_LENGTH][3];

    //
    // Structure used for storing (action, count) pairs from backtracking
    //
    typedef struct {
        char action[MAX_READ_LENGTH];
        int count[MAX_READ_LENGTH];
    } LocalCigarResult;

    //
    // Affine gap scoring parameters
    //
    int matchReward;
    int subPenalty;
    int gapOpenPenalty;
    int gapExtendPenalty;

};