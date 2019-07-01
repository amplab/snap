#include "stdafx.h"
#include "Compat.h"
#include "AffineGap.h"
#include "LandauVishkin.h"
#include "mapq.h"
#include "Read.h"
#include "BaseAligner.h"
#include "Bam.h"
#include "exit.h"
#include "Error.h"

// #define PRINT_SCORES 1

AffineGapWithCigar::AffineGapWithCigar(
    int i_matchReward,
    int i_subPenalty,
    int i_gapOpenPenalty,
    int i_gapExtendPenalty) : 
        matchReward(i_matchReward),
        subPenalty(-i_subPenalty),
        gapOpenPenalty(i_gapOpenPenalty + i_gapExtendPenalty),
        gapExtendPenalty(i_gapExtendPenalty)
{
    // Initialize score arrays
    memsetint(H, INT16_MIN, MAX_READ_LENGTH + 1);
    memsetint(E, INT16_MIN, MAX_READ_LENGTH + 1);

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

AffineGapWithCigar::AffineGapWithCigar() :
    matchReward(1),
    subPenalty(-4),
    gapOpenPenalty(7),
    gapExtendPenalty(1)
{
    // Initialize score arrays
    memsetint(H, INT16_MIN, MAX_READ_LENGTH + 1);
    memsetint(E, INT16_MIN, MAX_READ_LENGTH + 1);

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
AffineGapWithCigar::writeCigar(char** o_buf, int* o_buflen, int count, char code, CigarFormat format)
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

int AffineGapWithCigar::computeGlobalScore(const char* text, int textLen, const char* pattern, int patternLen, int w,
    char* cigarBuf, int cigarBufLen, bool useM,
    CigarFormat format,
    int* o_cigarBufUsed, int *o_netDel, int *o_tailIns) 
{
#ifdef PRINT_SCORES
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
    memsetint(H, INT16_MIN, MAX_READ_LENGTH + 1);
    memsetint(E, INT16_MIN, MAX_READ_LENGTH + 1);

    // 
    // Initialize H array
    //
    H[0] = 0;
    for (int j = 1; j <= w; ++j) {
        H[j] = -(gapOpenPenalty + (j - 1) * gapExtendPenalty);
    }

    int score = INT16_MIN, nEdits = -1; // Final alignment score and edit distance to be returned
    int textUsed = -1;

    // Iterate over all rows of text
    for (int i = 0; i < textLen; i++) {
        int beg = __max(i - w, 0);
        int end = __min(i + w + 1, patternLen);

        const char* t = (text + i);

        // Get the query profile for the row
        _int8* qRowProfile = &qProfile[BASE_VALUE[*t] * patternLen];

        int hLeft = INT16_MIN, fCurr = INT16_MIN;

        // Initialize first column
        if (beg == 0) {
            hLeft = -(gapOpenPenalty + i * gapExtendPenalty);
        }

#ifdef PRINT_SCORES
        for (int j = 0; j < beg; ++j) {
            printf(" ,");
        }
#endif

        // Iterate over all columns of pattern (within the band)
        for (int j = beg; j < end; ++j) {

            // 
            // Compute H(i,j).
            // H(i,j) = max(H(i-1,j-1) + qRowProfile[j], E(i,j), F(i,j))
            //
            int hDiag = H[j];
            hDiag = hDiag + qRowProfile[j];
            char action = (hDiag >= E[j]) ? 'M' : 'D';
            int hCurr = __max(hDiag, E[j]);
            action = (hCurr >= fCurr) ? action : 'I';
            hCurr = __max(hCurr, fCurr);
            backtraceAction[i][j - beg][0] = action;

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
            int hTemp = hDiag - gapOpenPenalty;
            E[j] -= gapExtendPenalty;
            action = (E[j] > hTemp) ? 'D' : 'M';
            E[j] = __max(hTemp, E[j]);
            backtraceAction[i][j - beg][1] = action;

            fCurr -= gapExtendPenalty;
            action = (fCurr > hTemp) ? 'I' : 'M';
            fCurr = __max(hTemp, fCurr);
            backtraceAction[i][j - beg][2] = action;

        } // end of pattern
        H[end] = hLeft;
        E[end] = INT16_MIN;

        if (end == patternLen) {
            if (score <= H[end]) {
                score = H[end];
                textUsed = i;
            }
        }

#ifdef PRINT_SCORES
        printf("\n");
#endif
    } // end of text
   
    LocalCigarResult res; // Track each (action, count) pair for generating CIGAR
    int n_res = 0; // Number of (action, count) pairs

    if (score > INT16_MIN) {
        int rowIdx = textUsed;
        int colIdx = min(rowIdx + w + 1, patternLen) - 1, matrixIdx = 0;
        char action = 'X', prevAction = 'X'; // initial
        int actionCount = 1;
        // Start traceback from the cell (i,j) with the maximum score
        while (rowIdx >= 0 && colIdx >= 0) {
            int actionIdx = (rowIdx > w) ? colIdx - (rowIdx - w) : colIdx; // ColIdx can span patternLen, while backtrack actions are stored only within band w.
            //if (actionIdx < 0) {
            //    WriteErrorMessage("Invalid traceback action: %d. (r,c,w) : (%d,%d,%d). Text:%.*s. Pattern:%.*s\n", actionIdx, rowIdx, colIdx, w, textLen, text, patternLen, pattern);
            //}
            _ASSERT(actionIdx >= 0);
            action = backtraceAction[rowIdx][actionIdx][matrixIdx];
            if (action == 'M') {
                rowIdx--;
                colIdx--;
                matrixIdx = 0;
            }
            else if (action == 'D') {
                rowIdx--;
                matrixIdx = 1;
            }
            else if (action == 'I') {
                colIdx--;
                matrixIdx = 2;
            }
            if (prevAction == action) {
                actionCount++;
            }
            else {
                if (prevAction != 'X') {
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
            res.action[n_res] = 'D';
            res.count[n_res] = actionCount;
            n_res++;
        }
        if (colIdx >= 0) {
            actionCount = colIdx + 1;
            res.action[n_res] = 'I';
            res.count[n_res] = actionCount;
            *o_tailIns = actionCount;
            n_res++;
        }

        char* cigarBufStart = cigarBuf;
        rowIdx = 0; colIdx = 0;  nEdits = 0;

        // Special handling of tail insertions which will be soft-clipped
        _ASSERT(n_res > 0);
        int min_i = 0;
        if (res.action[0] == 'I') {
            min_i = 1;
            *o_tailIns = res.count[0];
        }

        for (int i = n_res - 1; i >= min_i; --i) { // Reverse local result to obtain CIGAR in correct order
            if (useM) {
                if (!writeCigar(&cigarBuf, &cigarBufLen, res.count[i], res.action[i], format)) {
                    return -2;
                }
                // Compute edit distance NM:i flag
                if (res.action[i] == 'M') {
                    int j = 0, countM = 0;
                    for (int j = 0; j < res.count[i]; ++j) {
                        if (text[rowIdx + j] != pattern[colIdx + j]) {
                            nEdits++;
                        }
                    }
                    rowIdx += res.count[i];
                    colIdx += res.count[i];
                }
                else {
                    if (res.action[i] == 'D') {
                        rowIdx += res.count[i];
                        *o_netDel += res.count[i];
                    }
                    else if (res.action[i] == 'I') {
                        colIdx += res.count[i];
                    }
                    nEdits += res.count[i];
                }
            }
            else {
                if (res.action[i] == 'M') {
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
                    if (res.action[i] == 'D') {
                        rowIdx += res.count[i];
                    }
                    else if (res.action[i] == 'I') {
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

int AffineGapWithCigar::computeGlobalScoreNormalized(const char* text, int textLen,
    const char* pattern, int patternLen,
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
    int score = computeGlobalScore(text, (int)textLen, pattern, (int)patternLen, k, bamBuf, bamBufLen,
        useM, BAM_CIGAR_OPS, &bamBufUsed, o_netDel, o_tailIns);
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