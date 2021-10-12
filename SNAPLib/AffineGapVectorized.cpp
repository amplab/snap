#include "stdafx.h"
#include "Compat.h"
#include "AffineGapVectorized.h"
#include "LandauVishkin.h"
#include "mapq.h"
#include "Read.h"
#include "BaseAligner.h"
#include "Bam.h"
#include "exit.h"
#include "Error.h"

// #define PRINT_SCORES 1
// #define TRACE_AG 1

AffineGapVectorizedWithCigar::AffineGapVectorizedWithCigar(
    int i_matchReward,
    int i_subPenalty,
    int i_gapOpenPenalty,
    int i_gapExtendPenalty) :
    matchReward(i_matchReward),
    subPenalty(-i_subPenalty),
    gapOpenPenalty(i_gapOpenPenalty + i_gapExtendPenalty),
    gapExtendPenalty(i_gapExtendPenalty)
{
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
    minScoreParam = subPenalty;
    maxScoreParam = matchReward;
}

AffineGapVectorizedWithCigar::AffineGapVectorizedWithCigar() :
    matchReward(1),
    subPenalty(-4),
    gapOpenPenalty(7),
    gapExtendPenalty(1)
{
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

/*++
    Write cigar to buffer, return true if it fits
    null-terminates buffer if it returns false (i.e. fills up buffer)
--*/
bool
AffineGapVectorizedWithCigar::writeCigar(char** o_buf, int* o_buflen, int count, char code, CigarFormat format)
{
    _ASSERT(count >= 0);
    if (count <= 0) {
        return true;
    }
    switch (format) {
    case EXPANDED_CIGAR_STRING: {
        int n = min(*o_buflen, count);
        for (int i = 0; i < n; i++) {
            *(*o_buf)++ = code;
        }
        *o_buflen -= n;
        if (*o_buflen == 0) {
            *(*o_buf - 1) = '\0';
        }
        return *o_buflen > 0;
    }
    case COMPACT_CIGAR_STRING: {
        if (*o_buflen == 0) {
            *(*o_buf - 1) = '\0';
            return false;
        }
        int written = snprintf(*o_buf, *o_buflen, "%d%c", count, code);
        if (written > *o_buflen - 1) {
            *o_buf = '\0';
            return false;
        }
        else {
            *o_buf += written;
            *o_buflen -= written;
            return true;
        }
    }
    case COMPACT_CIGAR_BINARY:
        // binary format with non-zero count byte followed by char (easier to examine programmatically)
        while (true) {
            if (*o_buflen < 3) {
                *(*o_buf) = '\0';
                return false;
            }
            *(*o_buf)++ = min(count, 255);
            *(*o_buf)++ = code;
            *o_buflen -= 2;
            if (count <= 255) {
                return true;
            }
            count -= 255;
        }
    case BAM_CIGAR_OPS:
        if (*o_buflen < 4 || count >= (1 << 28)) {
            return false;
        }
        *(_uint32*)*o_buf = (count << 4) | BAMAlignment::CigarToCode[code];
        *o_buf += 4;
        *o_buflen -= 4;
        return true;
    default:
        WriteErrorMessage("invalid cigar format %d\n", format);
        soft_exit(1);
        return false;        // Not reached.  This is just here to suppress a compiler warning.
    } // switch
}

int AffineGapVectorizedWithCigar::computeGlobalScore(const char* text, int textLen, const char* pattern, const char* quality, int patternLen, int w,
    char* cigarBuf, int cigarBufLen, bool useM,
    CigarFormat format,
    int* o_cigarBufUsed, int *o_netDel, int *o_tailIns)
{
    _ASSERT(w < MAX_K);
    _ASSERT(textLen <= MAX_READ_LENGTH + MAX_K);

    w = __min(MAX_K - 1, w); // enforce limit even in non-debug builds

    if (NULL == text) {
        return -1;
    }

    int localNetDel = 0;

    if (o_netDel == NULL) {
        o_netDel = &localNetDel;
    }

    int localTailIns = 0;
    if (o_tailIns == NULL) {
        o_tailIns = &localTailIns;
    }

    *o_netDel = *o_tailIns = 0;

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
    __m128i v_intmin = _mm_set1_epi16(INT16_MIN);
    __m128i v_one = _mm_set1_epi16(1);
    __m128i v_two = _mm_set1_epi16(2);
    __m128i v_four = _mm_set1_epi16(4);
    __m128i v_thirtytwo = _mm_set1_epi16(32);
    __m128i v_gapOpen = _mm_set1_epi16(gapOpenPenalty);
    __m128i v_gapExtend = _mm_set1_epi16(gapExtendPenalty);
    __m128i v_mask = _mm_cmpgt_epi16(_mm_set_epi16(0, 0, 0, 0, 0, 0, 0, 1), v_zero);

    //
    // Initialize scores of first row
    //
    int16_t scoreFirstRow[VEC_SIZE] = {};
    for (int vecIdx = 0; vecIdx < numVec; vecIdx++) {
        for (int elemIdx = 0; elemIdx < VEC_SIZE; elemIdx++) {
            int patternIdx = elemIdx * numVec + vecIdx;
            if (patternIdx < patternLen) {
                scoreFirstRow[elemIdx] = -(gapOpenPenalty + patternIdx * gapExtendPenalty);
            }
            else {
                scoreFirstRow[elemIdx] = INT16_MIN;
            }
        }
        _ASSERT(VEC_SIZE == 8);  // FIXME: Initialization below works only when VEC_SIZE = 8
        _mm_store_si128(H + vecIdx, _mm_setr_epi16(scoreFirstRow[0], scoreFirstRow[1], scoreFirstRow[2], scoreFirstRow[3],
            scoreFirstRow[4], scoreFirstRow[5], scoreFirstRow[6], scoreFirstRow[7]));
        _mm_store_si128(E + vecIdx, v_intmin);
    }

    int16_t score = INT16_MIN, nEdits = -1; // Final alignment score and edit distance to be returned
    int textUsed = -1;

    __m128i* Hptr = H;
    __m128i* Hminus1ptr = Hminus1;

    // Iterate over all rows of text
    for (int i = 0; i < textLen; i++) {

        const char* t = (text + i);

        // Get the query profile for the row
        __m128i* qRowProfile = qProfile + BASE_VALUE[*t] * numVec;

        // Registers to hold intermediate scores for each row
        __m128i m, h, temp, e, f = v_intmin;

        // Load h from the previous row
        h = _mm_load_si128(Hptr + numVec - 1);

        // Shift left h and blend in initial values
        h = _mm_slli_si128(h, 2); // shift h left by 2 * 8 bits

        int16_t hInit = 0;
        if (i > 0) {
            hInit = -(gapOpenPenalty + (i - 1) * gapExtendPenalty);
        }
        __m128i v_hInit = _mm_set1_epi16(hInit);
        // h = _mm_blend_epi16(h, v_hInit, 1); // Does not work with SSE !
        h = blend_sse(h, v_hInit, v_mask);

        // Iterate over all columns of pattern (within the band)
        for (int j = 0; j < numVec; ++j) {
            __m128i backtraceActionVec;

            m = _mm_adds_epi16(h, *qRowProfile++);

            e = _mm_load_si128(E + j);

            // h = max{m, e, f}
            backtraceActionVec = _mm_and_si128(_mm_cmpgt_epi16(e, m), v_one); // action = e > m ? 1 : 0
            h = _mm_max_epi16(m, e);
            __m128i tmpResult = _mm_and_si128(_mm_cmpgt_epi16(f, h), v_two);
            backtraceActionVec = _mm_or_si128(tmpResult, _mm_andnot_si128(tmpResult, backtraceActionVec)); // action = f > h ? 2 : action 
            h = _mm_max_epi16(h, f);

            // Store h for the next row
            _mm_store_si128(Hminus1ptr + j, h);

            // e = max{m - gapOpen, e - gapExtend}
            e = _mm_subs_epi16(e, v_gapExtend);
            temp = _mm_subs_epi16(m, v_gapOpen);
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

        } // end of pattern

        //
        // Farrar's algorithm does lazy f evaluation. Since f rarely influences final score h, the algorithm speculates f = zero initially 
        // Re-evaluate if f could influence h after we have made a first pass through the row
        //
        for (int k = 0; k < (VEC_SIZE - 1); k++) {

            // Extract the last f vector
            f = _mm_slli_si128(f, 2);
            f = blend_sse(f, v_intmin, v_mask);

            for (int j = 0; j < numVec; j++) {

                h = _mm_load_si128(Hminus1ptr + j);

                // action = f > h ? 2 : action
                __m128i tmpResult = _mm_and_si128(_mm_cmpgt_epi16(f, h), v_two);
                __m128i backtraceActionVec = _mm_load_si128(backtraceAction + (i * numVec + j));
                __m128i tmpResult2 = _mm_andnot_si128(tmpResult, backtraceActionVec);
                backtraceActionVec = _mm_or_si128(tmpResult, tmpResult2);

                h = _mm_max_epi16(h, f);
                _mm_store_si128(Hminus1ptr + j, h);

                temp = _mm_subs_epi16(h, v_gapOpen);
                f = _mm_subs_epi16(f, v_gapExtend);

                tmpResult = _mm_and_si128(_mm_cmpgt_epi16(f, temp), v_thirtytwo);
                backtraceActionVec = _mm_or_si128(backtraceActionVec, tmpResult);
                _mm_store_si128(backtraceAction + (i * numVec + j), backtraceActionVec);

                // Converged if no element of f can influence h
                bool converged = !(_mm_movemask_epi8(_mm_cmpgt_epi16(f, temp)));

                if (converged) goto got_answer;
            }
        }

    got_answer:
        // Global alignment score (i.e., score when aligning to the end of the pattern)
        __m128i v_globalAlignmentScore = _mm_load_si128(Hminus1ptr + ((patternLen - 1) % numVec));
        int16_t globalAlignmentScore = getElem((patternLen - 1) / numVec, v_globalAlignmentScore);
        if (globalAlignmentScore >= score) {
            score = globalAlignmentScore;
            textUsed = i;
        }

#ifdef PRINT_SCORES
        __m128i rowScore;
        for (int j = 0; j < numVec; j++) {
            rowScore = _mm_load_si128(Hminus1ptr + j);
            printElem(rowScore);
        }
        printf("\n");
#endif
    
        // Swap roles of H and Hminus1 for the next row
        __m128i* hTemp = Hminus1ptr;
        Hminus1ptr = Hptr;
        Hptr = hTemp;

    } // end of text

    int n_res = 0; // Number of (action, count) pairs

    if (score > INT16_MIN) {
        int rowIdx = textUsed;
        int colIdx = patternLen - 1, matrixIdx = 0;
        BacktraceActionType action = M, prevAction = X;
        int actionCount = 1;

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
            action = (BacktraceActionType)((backtracePointersRow[stripedColIdx] >> matrixIdx) & 3);

            if (action == M) {
                rowIdx--;
                colIdx--;
                matrixIdx = 0;
            }
            else if (action == D) {
                rowIdx--;
                matrixIdx = 1;
            }
            else {
                colIdx--;
                action = I;
                matrixIdx = 2;
            }
            if (prevAction == action) {
                actionCount++;
            }
            else {
                if (prevAction != X) {
                    res.action[n_res] = prevAction;
                    res.count[n_res] = actionCount;
                    n_res++;
                    actionCount = 1;
                }
            }
            prevAction = action;
        }
        // We went past either pattern or text. Dump out any action counts that we have been tracking
        if (prevAction == action) {
            res.action[n_res] = prevAction;
            res.count[n_res] = actionCount;
            n_res++;
        }
        if (rowIdx >= 0) {
            actionCount = rowIdx + 1;
            res.action[n_res] = D;
            res.count[n_res] = actionCount;
            n_res++;
        }
        if (colIdx >= 0) {
            actionCount = colIdx + 1;
            res.action[n_res] = I;
            res.count[n_res] = actionCount;
            *o_tailIns = actionCount;
            n_res++;
        }

        char* cigarBufStart = cigarBuf;
        rowIdx = 0; colIdx = 0;  nEdits = 0;

        // Special handling of tail insertions which will be soft-clipped
        _ASSERT(n_res > 0);
        int min_i = 0;
        if (res.action[0] == I) {
            min_i = 1;
            *o_tailIns = res.count[0];
        }

        // Flip order of insertions followed by substitutions in CIGAR string
        for (int i = n_res - 1; i >= min_i; --i) {
            if (res.action[i] == M) {
                rowIdx += res.count[i];
                colIdx += res.count[i];
            }
            else if (res.action[i] == D) {
                rowIdx += res.count[i];
            }
            else {
                if (i > 0 && rowIdx < textUsed && colIdx < patternLen - 1) {
                    if ((pattern[colIdx + 1] == pattern[colIdx]) && (pattern[colIdx + 1] != text[rowIdx]) && (quality[colIdx] < 65)) { // only insert high quality bases if possible
                        if ((i + 1 <= n_res - 1) && res.action[i + 1] == M && res.count[i - 1] > 1) {
                            res.count[i + 1] += 1;
                            rowIdx++;
                            colIdx++;
                        }
                        if (res.action[i - 1] == M && res.count[i - 1] > 1) {
                            res.count[i - 1] -= 1;
                        }
                    }
                }
                colIdx += res.count[i];
            }
        }

        rowIdx = 0; colIdx = 0;
        // Flip order of insertions and substitution with match in between in CIGAR string
        for (int i = n_res - 1; i >= min_i; --i) {
            if (res.action[i] == M) {
                rowIdx += res.count[i];
                colIdx += res.count[i];
            }
            else if (res.action[i] == D) {
                rowIdx += res.count[i];
            }
            else {
                if (i > 0 && rowIdx + 1 < textUsed && colIdx + res.count[i] < patternLen - 1) {
                    if ((pattern[colIdx + res.count[i]] == pattern[colIdx]) && (pattern[colIdx + res.count[i] + 1] != text[rowIdx + 1]) && (quality[colIdx] < 65)) { // only insert high quality bases (>= 65) if possible
                        if ((i + 1 <= n_res - 1) && res.action[i + 1] == M && res.count[i - 1] > 2) {
                            res.count[i + 1] += 2;
                            rowIdx += 2;
                            colIdx += 2;
                        }
                        if (res.action[i - 1] == M && res.count[i - 1] > 2) {
                            res.count[i - 1] -= 2;
                        }
                    }
                }
                colIdx += res.count[i];
            }
        }

        rowIdx = 0; colIdx = 0;
        for (int i = n_res - 1; i >= min_i; --i) { // Reverse local result to obtain CIGAR in correct order
            if (useM) {
                // Compute edit distance NM:i flag
                if (res.action[i] == M) {
                    int j = 0, countM = 0;
                    for (int j = 0; j < res.count[i]; ++j) {
                        if (text[rowIdx + j] != pattern[colIdx + j]) {
                            nEdits++;
                        }
                    }
                    rowIdx += res.count[i];
                    colIdx += res.count[i];
                    if (!writeCigar(&cigarBuf, &cigarBufLen, res.count[i], 'M', format)) {
                        return -2;
                    }
                }
                else {
                    if (res.action[i] == D) {
                        rowIdx += res.count[i];
                        *o_netDel += res.count[i];
                        if (!writeCigar(&cigarBuf, &cigarBufLen, res.count[i], 'D', format)) {
                            return -2;
                        }
                    }
                    else if (res.action[i] == I) {
                        colIdx += res.count[i];
                        if (!writeCigar(&cigarBuf, &cigarBufLen, res.count[i], 'I', format)) {
                            return -2;
                        }
                    }
                    nEdits += res.count[i];
                }
            }
            else {
                if (res.action[i] == M) {
                    int j = 0, countM = 0;
                    for (int j = 0; j < res.count[i]; ++j) {
                        if (text[rowIdx + j] != pattern[colIdx + j]) {
                            // Write the matches '='
                            if (countM > 0) {
                                if (!writeCigar(&cigarBuf, &cigarBufLen, countM, '=', format)) {
                                    return -2;
                                }
                            }
                            // Write the mismatching 'X'
                            if (!writeCigar(&cigarBuf, &cigarBufLen, 1, 'X', format)) {
                                return -2;
                            }
                            countM = 0;
                        }
                        else {
                            countM++;
                        }
                    }
                    rowIdx += res.count[i];
                    colIdx += res.count[i];
                }
                else {
                    if (!writeCigar(&cigarBuf, &cigarBufLen, res.count[i], res.action[i], format)) {
                        return -2;
                    }
                    if (res.action[i] == D) {
                        rowIdx += res.count[i];
                    }
                    else if (res.action[i] == I) {
                        colIdx += res.count[i];
                    }
                }
            }
        }

        if (format != BAM_CIGAR_OPS) {
            *(cigarBuf - (cigarBufLen == 0 ? 1 : 0)) = '\0'; // terminate string
        }
        if (o_cigarBufUsed != NULL) {
            *o_cigarBufUsed = (int)(cigarBuf - cigarBufStart);
        }
        // nEdits = (nEdits <= w) ? nEdits : -1; // return -1 if we have more edits than threshold w
        return nEdits;

    } // score > 0
    else {
        // Could not align strings with at most K edits
        *(cigarBuf - (cigarBufLen == 0 ? 1 : 0)) = '\0'; // terminate string
        return -1;
    }
}

int AffineGapVectorizedWithCigar::computeGlobalScoreBanded(
    const char* text,
    int textLen,
    const char* pattern,
    const char* quality,
    int patternLen,
    int w,
    int scoreInit,
    char* cigarBuf,
    int cigarBufLen,
    bool useM,
    CigarFormat format,
    int* o_cigarBufUsed,
    int* o_netDel,
    int* o_tailIns)
{

#ifdef TRACE_AG
    printf("\n");
#endif

    _ASSERT(w < MAX_K);
    _ASSERT(textLen <= MAX_READ_LENGTH + MAX_K);

    w = __min(MAX_K - 1, w); // enforce limit even in non-debug builds

    if (NULL == text) {
        return -1;
    }

    int localNetDel = 0;

    if (o_netDel == NULL) {
        o_netDel = &localNetDel;
    }

    int localTailIns = 0;
    if (o_tailIns == NULL) {
        o_tailIns = &localTailIns;
    }

    *o_netDel = *o_tailIns = 0;

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
    __m128i v_intmin = _mm_set1_epi16(INT16_MIN);
    __m128i v_one = _mm_set1_epi16(1);
    __m128i v_two = _mm_set1_epi16(2);
    __m128i v_four = _mm_set1_epi16(4);
    __m128i v_thirtytwo = _mm_set1_epi16(32);
    __m128i v_gapOpen = _mm_set1_epi16(gapOpenPenalty);
    __m128i v_gapExtend = _mm_set1_epi16(gapExtendPenalty);
    __m128i v_mask = _mm_cmpgt_epi16(_mm_set_epi16(0, 0, 0, 0, 0, 0, 0, 1), v_zero);

    //
    // Initialize scores of first row
    //
    int16_t scoreFirstRow[VEC_SIZE] = { 0 };
    for (int segIdx = 0; segIdx < numSeg; segIdx++) {
        for (int vecIdx = 0; vecIdx < numVec; vecIdx++) {
            for (int elemIdx = 0; elemIdx < VEC_SIZE; elemIdx++) {
                int patternIdx = (segIdx * segLen) + (elemIdx * numVec) + vecIdx;
                if (patternIdx < patternLen) {
                    scoreFirstRow[elemIdx] =  __max(0, scoreInit - (gapOpenPenalty + patternIdx * gapExtendPenalty));
                }
            }
            _ASSERT(VEC_SIZE == 8);  // FIXME: Initialization below works only when VEC_SIZE = 8
            _mm_store_si128(H + (segIdx * numVec) + vecIdx, _mm_setr_epi16(scoreFirstRow[0], scoreFirstRow[1], scoreFirstRow[2], scoreFirstRow[3],
                scoreFirstRow[4], scoreFirstRow[5], scoreFirstRow[6], scoreFirstRow[7]));
            _mm_store_si128(Hminus1 + (segIdx * numVec) + vecIdx, v_zero);
            _mm_store_si128(E + (segIdx * numVec) + vecIdx, v_zero);
        }
    }

    int score = scoreInit, nEdits = -1; // Final alignment score and edit distance to be returned. 
    int textUsed = -1;

    __m128i* Hptr = H;
    __m128i* Hminus1ptr = Hminus1;

    // Iterate over all rows of text
    for (int i = 0; i < textLen; i++) {

        const char* t = text + i;

        // Get the query profile for the row
        __m128i* qRowProfile = qProfile + BASE_VALUE[*t] * numSeg * numVec;

        // Registers to hold intermediate scores for each row
        __m128i m, h, temp, e = v_zero, f = v_zero, max = v_zero, X = v_zero;

        int maxScoreRow = 0;
        int localAlignmentPatternOffset = -1;

        int bandBeg = __max(i - w, 0);
        int bandEnd = __min(i + w, patternLen - 1);
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
                    hInit = scoreInit - (gapOpenPenalty + (i - 1) * gapExtendPenalty);
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

                __m128i qProfileVec = _mm_load_si128(qRowProfile + j * numVec + k);
                m = _mm_adds_epi16(h, qProfileVec);

                e = _mm_load_si128(E + j * numVec + k);

                // h = max{m, e, f}
                backtraceActionVec = _mm_and_si128(_mm_cmpgt_epi16(e, m), v_one); // action = e > m ? 1 : 0
                h = _mm_max_epi16(m, e);
                __m128i tmpResult = _mm_and_si128(_mm_cmpgt_epi16(f, h), v_two);
                backtraceActionVec = _mm_or_si128(tmpResult, _mm_andnot_si128(tmpResult, backtraceActionVec)); // action = f > h ? 2 : action 
                h = _mm_max_epi16(h, f);

                // Store h for the next row
                _mm_store_si128(Hminus1ptr + j * numVec + k, h);

                // e = max{m - gapOpen, e - gapExtend}
                e = _mm_subs_epi16(e, v_gapExtend);
                temp = _mm_subs_epi16(m, v_gapOpen);

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

                    temp = _mm_subs_epi16(h, v_gapOpen);
                    f = _mm_subs_epi16(f, v_gapExtend);

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
            if (globalAlignmentScore > score) {
                score = globalAlignmentScore;
                textUsed = i;
            }
        }

        // Swap roles of H and Hminus1 for the next row
        __m128i* hTemp = Hminus1ptr;
        Hminus1ptr = Hptr;
        Hptr = hTemp;

    } // end text

    int n_res = 0; // Number of (action, count) pairs

    if (score > 0) {
        int rowIdx = textUsed;
        int colIdx = patternLen - 1, matrixIdx = 0;
        BacktraceActionType action = M, prevAction = X;
        int actionCount = 1;

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
            action = (BacktraceActionType)((backtracePointersVec[elemIdx] >> matrixIdx) & 3);

            if (action == M) {
                rowIdx--;
                colIdx--;
                matrixIdx = 0;
            }
            else if (action == D) {
                rowIdx--;
                matrixIdx = 1;
            }
            else {
                colIdx--;
                action = I;
                matrixIdx = 2;
            }
            if (prevAction == action) {
                actionCount++;
            }
            else {
                if (prevAction != X) {
                    res.action[n_res] = prevAction;
                    res.count[n_res] = actionCount;
                    n_res++;
                    actionCount = 1;
                }
            }
            prevAction = action;
        }
        // We went past either pattern or text. Dump out any action counts that we have been tracking
        if (prevAction == action) {
            res.action[n_res] = prevAction;
            res.count[n_res] = actionCount;
            n_res++;
        }
        if (rowIdx >= 0) {
            actionCount = rowIdx + 1;
            res.action[n_res] = D;
            res.count[n_res] = actionCount;
            n_res++;
        }
        if (colIdx >= 0) {
            actionCount = colIdx + 1;
            res.action[n_res] = I;
            res.count[n_res] = actionCount;
            *o_tailIns = actionCount;
            n_res++;
        }

        char* cigarBufStart = cigarBuf;
        rowIdx = 0; colIdx = 0;  nEdits = 0;

        // Special handling of tail insertions which will be soft-clipped
        _ASSERT(n_res > 0);
        int min_i = 0;
        if (res.action[0] == I) {
            min_i = 1;
            *o_tailIns = res.count[0];
        }

        // Flip order of insertions followed by substitutions in CIGAR string
        for (int i = n_res - 1; i >= min_i; --i) {
            if (res.action[i] == M) {
                rowIdx += res.count[i];
                colIdx += res.count[i];
            }
            else if (res.action[i] == D) {
                rowIdx += res.count[i];
            }
            else {
                if (i > 0 && rowIdx < textUsed && colIdx < patternLen - 1) {
                    if ((pattern[colIdx + 1] == pattern[colIdx]) && (pattern[colIdx + 1] != text[rowIdx]) && (quality[colIdx] < 65)) { // only insert high quality bases (>= 65) if possible
                        if ((i + 1 <= n_res - 1) && res.action[i + 1] == M && res.count[i - 1] > 1) {
                            res.count[i + 1] += 1;
                            rowIdx++;
                            colIdx++;
                        }
                        if (res.action[i - 1] == M && res.count[i - 1] > 1) {
                            res.count[i - 1] -= 1;
                        }
                    }
                }
                colIdx += res.count[i];
            }
        }

        rowIdx = 0; colIdx = 0;
        // Flip order of insertions and substitution with match in between in CIGAR string
        for (int i = n_res - 1; i >= min_i; --i) {
            if (res.action[i] == M) {
                rowIdx += res.count[i];
                colIdx += res.count[i];
            }
            else if (res.action[i] == D) {
                rowIdx += res.count[i];
            }
            else {
                if (i > 0 && rowIdx + 1 < textUsed && colIdx + res.count[i] < patternLen - 1) {
                    if ((pattern[colIdx + res.count[i]] == pattern[colIdx]) && (pattern[colIdx + res.count[i] + 1] != text[rowIdx + 1]) && (quality[colIdx] < 65)) { // only insert high quality bases (>= 65) if possible
                        if ((i + 1 <= n_res - 1) && res.action[i + 1] == M && res.count[i - 1] > 2) {
                            res.count[i + 1] += 2;
                            rowIdx += 2;
                            colIdx += 2;
                        }
                        if (res.action[i - 1] == M && res.count[i - 1] > 2) {
                            res.count[i - 1] -= 2;
                        }
                    }
                }
                colIdx += res.count[i];
            }
        }

        rowIdx = 0; colIdx = 0;
        for (int i = n_res - 1; i >= min_i; --i) { // Reverse local result to obtain CIGAR in correct order
            if (useM) {
                // Compute edit distance NM:i flag
                if (res.action[i] == M) {
                    int j = 0, countM = 0;
                    for (int j = 0; j < res.count[i]; ++j) {
                        if (text[rowIdx + j] != pattern[colIdx + j]) {
                            nEdits++;
                        }
                    }
                    rowIdx += res.count[i];
                    colIdx += res.count[i];
                    if (!writeCigar(&cigarBuf, &cigarBufLen, res.count[i], 'M', format)) {
                        return -2;
                    }
                }
                else {
                    if (res.action[i] == D) {
                        rowIdx += res.count[i];
                        *o_netDel += res.count[i];
                        if (!writeCigar(&cigarBuf, &cigarBufLen, res.count[i], 'D', format)) {
                            return -2;
                        }
                    }
                    else if (res.action[i] == I) {
                        colIdx += res.count[i];
                        if (!writeCigar(&cigarBuf, &cigarBufLen, res.count[i], 'I', format)) {
                            return -2;
                        }
                    }
                    nEdits += res.count[i];
                }
            }
            else {
                if (res.action[i] == M) {
                    int j = 0, countM = 0;
                    for (int j = 0; j < res.count[i]; ++j) {
                        if (text[rowIdx + j] != pattern[colIdx + j]) {
                            // Write the matches '='
                            if (countM > 0) {
                                if (!writeCigar(&cigarBuf, &cigarBufLen, countM, '=', format)) {
                                    return -2;
                                }
                            }
                            // Write the mismatching 'X'
                            if (!writeCigar(&cigarBuf, &cigarBufLen, 1, 'X', format)) {
                                return -2;
                            }
                            countM = 0;
                        }
                        else {
                            countM++;
                        }
                    }
                    rowIdx += res.count[i];
                    colIdx += res.count[i];
                }
                else {
                    if (!writeCigar(&cigarBuf, &cigarBufLen, res.count[i], res.action[i], format)) {
                        return -2;
                    }
                    if (res.action[i] == D) {
                        rowIdx += res.count[i];
                    }
                    else if (res.action[i] == I) {
                        colIdx += res.count[i];
                    }
                }
            }
        }

        if (format != BAM_CIGAR_OPS) {
            *(cigarBuf - (cigarBufLen == 0 ? 1 : 0)) = '\0'; // terminate string
        }
        if (o_cigarBufUsed != NULL) {
            *o_cigarBufUsed = (int)(cigarBuf - cigarBufStart);
        }
        // nEdits = (nEdits <= w) ? nEdits : -1; // return -1 if we have more edits than threshold w
        return nEdits;

    } // score > scoreInit
    else {
        // Could not align strings with at most K edits
        *(cigarBuf - (cigarBufLen == 0 ? 1 : 0)) = '\0'; // terminate string
        return -1;
    }
}

int AffineGapVectorizedWithCigar::computeGlobalScoreNormalized(const char* text, int textLen,
    const char* pattern, const char* quality, int patternLen,
    int k,
    char *cigarBuf, int cigarBufLen, bool useM,
    CigarFormat format, int* o_cigarBufUsed,
    int* o_addFrontClipping, int *o_netDel, int *o_tailIns)
{
    if (format != BAM_CIGAR_OPS && format != COMPACT_CIGAR_STRING) {
        WriteErrorMessage("AffineGapWithCigar::computeGlobalScoreNormalized invalid parameter\n");
        soft_exit(1);
    }
    int bamBufLen = (format == BAM_CIGAR_OPS ? 1 : 2) * cigarBufLen; // should be enough
    char* bamBuf = (char*)alloca(bamBufLen);
    int bamBufUsed;
    int score;
    if (patternLen >= (3 * (2 * k + 1))) {
        score = computeGlobalScoreBanded(text, (int)textLen, pattern, quality, (int)patternLen, k, MAX_READ_LENGTH, bamBuf, bamBufLen,
            useM, BAM_CIGAR_OPS, &bamBufUsed, o_netDel, o_tailIns);
        if (score < 0 || score > k) {
            score = computeGlobalScore(text, (int)textLen, pattern, quality, (int)patternLen, k, bamBuf, bamBufLen,
                useM, BAM_CIGAR_OPS, &bamBufUsed, o_netDel, o_tailIns);
        }
    } else {
        score = computeGlobalScore(text, (int)textLen, pattern, quality, (int)patternLen, k, bamBuf, bamBufLen,
            useM, BAM_CIGAR_OPS, &bamBufUsed, o_netDel, o_tailIns);
    }

    if (score < 0) {
        return score;
    }

    _uint32* bamOps = (_uint32*)bamBuf;
    int bamOpCount = bamBufUsed / sizeof(_uint32);

    // FIXME: Need to check if this is still valid with affine gap
    // _ASSERT('I' != BAMAlignment::CodeToCigar[BAMAlignment::GetCigarOpCode(bamOps[bamOpCount - 1])]);

    if (o_addFrontClipping != NULL) {
        char firstCode = BAMAlignment::CodeToCigar[BAMAlignment::GetCigarOpCode(bamOps[0])];
        if (firstCode == 'D') {
            *o_addFrontClipping = BAMAlignment::GetCigarOpCount(bamOps[0]);
            if (*o_addFrontClipping != 0) {
                return 0; // can fail, will be rerun with new clipping
            }
        }
        else if (firstCode == 'I') {
            *o_addFrontClipping = -1 * BAMAlignment::GetCigarOpCount(bamOps[0]);
        }
        else {
            *o_addFrontClipping = 0;
        }
    }

    // FIXME: Need to check if this is still valid with affine gap
    // _ASSERT(bamOpCount <= 1 || BAMAlignment::CodeToCigar[BAMAlignment::GetCigarOpCode(bamOps[bamOpCount - 1])] != 'I');	// We should have cleared all of these out
    _ASSERT(bamOpCount <= 1 || BAMAlignment::CodeToCigar[BAMAlignment::GetCigarOpCode(bamOps[bamOpCount - 1])] != 'D');	// And none of these should happen, either.
    // Seems to happen; TODO: fix this	_ASSERT(bamOpCount <= 1 || BAMAlignment::CodeToCigar[BAMAlignment::GetCigarOpCode(bamOps[0])] != 'D');
    // Seems to happen; TODO: fix this	_ASSERT(bamOpCount <= 1 || BAMAlignment::CodeToCigar[BAMAlignment::GetCigarOpCode(bamOps[0])] != 'I');

    // copy out cigar info
    if (format == BAM_CIGAR_OPS) {
        memcpy(cigarBuf, bamOps, bamBufUsed);
        if (o_cigarBufUsed != NULL) {
            *o_cigarBufUsed = bamBufUsed;
        }
    }
    else {
        bool ok = BAMAlignment::decodeCigar(cigarBuf, cigarBufLen, bamOps, bamOpCount);
        if (!ok) {
            return -1;
        }
        if (o_cigarBufUsed != NULL) {
            *o_cigarBufUsed = (int)strlen(cigarBuf) + 1;
        }
    }
    return score;
}