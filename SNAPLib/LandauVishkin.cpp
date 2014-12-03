#include "stdafx.h"
#include "Compat.h"
#include "LandauVishkin.h"
#include "mapq.h"
#include "Read.h"
#include "BaseAligner.h"
#include "Bam.h"
#include "exit.h"
#include "Error.h"

using std::make_pair;
using std::min;

 
LandauVishkinWithCigar::LandauVishkinWithCigar()
{
    for (int i = 0; i < MAX_K+1; i++) {
        for (int j = 0; j < 2*MAX_K+1; j++) {
            L[i][j] = -2;
        }
    }
    totalIndels[0][MAX_K] = 0;
}

/*++
    Write cigar to buffer, return true if it fits
    null-terminates buffer if it returns false (i.e. fills up buffer)
--*/
bool writeCigar(char** o_buf, int* o_buflen, int count, char code, CigarFormat format)
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
        } else {
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
        WriteErrorMessage( "invalid cigar format %d\n", format);
        soft_exit(1);
        return false;        // Not reached.  This is just here to suppress a compiler warning.
    } // switch
}

#if 0
inline void validateAction(char& last, char current)
{
    _ASSERT(last != current);
    last = current;
}
#else
inline void validateAction(char& last, char current) {}
#endif


    void
LandauVishkinWithCigar::printLinear(
    char* buffer,
    int bufferSize,
    unsigned variant)
{
    _ASSERT(bufferSize >= 12);
    int inserts = (variant >> CigarInsertCShift) & CigarInsertCount;
    if (inserts > 0) {
        *buffer++ = '0' + inserts;
        *buffer++ = 'I';
        for (int i = 0; i < inserts; i++) {
            *buffer++ = VALUE_BASE[(variant >> (CigarInsertBShift + 2 * i)) & 3];
        }
    }
    unsigned op = variant & CigarOpcode;
    if (op >= CigarReplace && op < CigarDelete) {
        *buffer++ = 'X';
        *buffer++ = VALUE_BASE[op - CigarReplace];
    } else if (op == CigarDelete) {
        *buffer++ = 'D';
    }
    *buffer++ = 0;
}

#if 0
static const int PrevDelta[3][3] =  // Version that minimizes NET indels (ie., |#ins - #del|
    {{+1, 0, -1},    // d < 0
    {0, -1, +1},    // d == 0
    {-1, 0, +1}};   // d > 0
#else    // 0
    static const int PrevDelta[3][3] = // Version that minimizes absolute indels (ie., #ins + #del)
    { { 0, +1, -1},      // d < 0
      { 0, +1, -1 },     // d == 0
      { 0, -1, +1 } };   // d > 0
#endif // 0

int LandauVishkinWithCigar::computeEditDistance(
    const char* text, int textLen,
    const char* pattern, int patternLen,
    int k,
    char *cigarBuf, int cigarBufLen, bool useM, 
    CigarFormat format, int* o_cigarBufUsed, int* o_textUsed,
    int *o_netIndel)
{
    int localNetIndel;
    if (NULL == o_netIndel) {
        //
        // If the user doesn't want netIndel, just use a stack local to avoid
        // having to check it all the time.
        //
        o_netIndel = &localNetIndel;
    }
    _ASSERT(k < MAX_K);

    *o_netIndel = 0;
    
    _ASSERT(patternLen >= 0 && textLen >= 0);
    _ASSERT(k < MAX_K);
    const char* p = pattern;
    const char* t = text;
    char* cigarBufStart = cigarBuf;
    if (NULL == text) {
        return -1;            // This happens when we're trying to read past the end of the genome.
    }

    int end = min(patternLen, textLen);
    const char* pend = pattern + end;
    while (p < pend) {
        _uint64 x = *((_uint64*) p) ^ *((_uint64*) t);
        if (x) {
            unsigned long zeroes;
            CountTrailingZeroes(x, zeroes);
            zeroes >>= 3;
            L[0][MAX_K] = min((int)(p - pattern) + (int)zeroes, end);
            goto done1;
        }
        p += 8;
        t += 8;
    }
    L[0][MAX_K] = end;
done1:
    if (L[0][MAX_K] == end) {
        // We matched the text exactly; fill the CIGAR string with all ='s (or M's)
		if (useM) {
			if (! writeCigar(&cigarBuf, &cigarBufLen, patternLen, 'M', format)) {
				return -2;
			}
            // todo: should this also write X's like '=' case? or is 'M' special?
		} else {
			if (! writeCigar(&cigarBuf, &cigarBufLen, end, '=', format)) {
				return -2;
			}
			if (patternLen > end) {
				// Also need to write a bunch of X's past the end of the text
				if (! writeCigar(&cigarBuf, &cigarBufLen, patternLen - end, 'X', format)) {
					return -2;
				}
			}
		}
        // todo: should this null-terminate?
        if (o_cigarBufUsed != NULL) {
            *o_cigarBufUsed = (int)(cigarBuf - cigarBufStart);
        }
        if (o_textUsed != NULL) {
            *o_textUsed = end;
        }
        return 0;
    }

    char lastAction = '*';

	int e;
	int lastBestIndels = MAX_K + 1;
    int lastBestD = MAX_K + 1;
	int lastBestBest;

    for (e = 1; e <= k; e++) {
        // Go through the offsets, d, in the order 0, -1, 1, -2, 2, etc, in order to find CIGAR strings
        // with few indels first if possible.
        for (int d = 0; d != -(e+1); d = (d >= 0 ? -(d+1) : -d)) {
            int bestdelta = 0;
            int bestbest = -1;
            int bestBestIndels = MAX_K + 1;
            //  extend previous solutions as far as possible, pick best, minimizing indels
            int dy = (d >= 0) + (d > 0);
            for (int dx = 0; dx < 3; dx++) {
                int delta = PrevDelta[dy][dx];
                int best = L[e-1][MAX_K+d + delta] + (delta >= 0);
                int bestIndels = totalIndels[e - 1][MAX_K + d + delta] + (delta != 0);  // Our parent, plus one if this is an indel
                if (best < 0) {
                    continue;
                }
                const char* p = pattern + best;
                const char* t = (text + d) + best;
                if (*p == *t) {
                    int end = min(patternLen, textLen - d);
                    const char* pend = pattern + end;

                    while (true) {
                        _uint64 x = *((_uint64*) p) ^ *((_uint64*) t);
                        if (x) {
                            unsigned long zeroes;
                            CountTrailingZeroes(x, zeroes);
                            zeroes >>= 3;
                            best = min((int)(p - pattern) + (int)zeroes, end);
                            break;
                        }
                        p += 8;
                        if (p >= pend) {
                            best = end;
                            break;
                        }
                        t += 8;
                    }
                }
                if (best > bestbest || best == bestbest && bestIndels < bestBestIndels) {
                    bestbest = best;
                    bestdelta = delta;
                    bestBestIndels = bestIndels;
                }
            }
            int best = bestbest;
            A[e][MAX_K+d] = "DXI"[bestdelta + 1];

            L[e][MAX_K+d] = best;
            totalIndels[e][MAX_K + d] = bestBestIndels;

			if (best == patternLen) {

				if (bestBestIndels == 0) {
					lastBestIndels = bestBestIndels;
                    lastBestD = d;
					lastBestBest = best;
					goto got_answer;
				}

				if (abs(lastBestIndels) > bestBestIndels) {
                    lastBestIndels = bestBestIndels;
                    lastBestD = d;
					lastBestBest = best;
				}
            } // if best == patternlen
        } // for d

		if (lastBestD != MAX_K + 1) {
			goto got_answer;
		}
    } // for e

    // Could not align strings with at most K edits
    *(cigarBuf - (cigarBufLen == 0 ? 1 : 0)) = '\0'; // terminate string
    return -1;

got_answer:
	// We're done. First, let's see whether we can reach e errors with no indels. Otherwise, we'll
	// trace back through the dynamic programming array to build up the CIGAR string.

	int straightMismatches = 0;
	for (int i = 0; i < end; i++) {
		if (pattern[i] != text[i]) {
			straightMismatches++;
		}
	}
	straightMismatches += patternLen - end;
	if (straightMismatches == e) {
		// We can match with no indels; let's do that
		if (useM) {
			//
			// No inserts or deletes, and with useM equal and SNP look the same, so just
			// emit a simple string.
			//
			validateAction(lastAction, 'M');
			if (!writeCigar(&cigarBuf, &cigarBufLen, patternLen, 'M', format)) {
				return -2;
			}
		}
		else {
			int streakStart = 0;
			bool matching = (pattern[0] == text[0]);
			for (int i = 0; i < end; i++) {
				bool newMatching = (pattern[i] == text[i]);
				if (newMatching != matching) {
					validateAction(lastAction, matching ? '=' : 'X');
					if (!writeCigar(&cigarBuf, &cigarBufLen, i - streakStart, (matching ? '=' : 'X'), format)) {
						return -2;
					}
					matching = newMatching;
					streakStart = i;
				}
			}

			// Write the last '=' or 'X' streak
			if (patternLen > streakStart) {
				if (!matching) {
					// Write out X's all the way to patternLen
					validateAction(lastAction, 'X');
					if (!writeCigar(&cigarBuf, &cigarBufLen, patternLen - streakStart, 'X', format)) {
						return -2;
					}
				}
				else {
					// Write out some ='s and then possibly X's if pattern is longer than text
					validateAction(lastAction, '=');
					if (!writeCigar(&cigarBuf, &cigarBufLen, end - streakStart, '=', format)) {
						return -2;
					}
					if (patternLen > end) {
						validateAction(lastAction, 'X');
						if (!writeCigar(&cigarBuf, &cigarBufLen, patternLen - end, 'X', format)) {
							return -2;
						}
					}
				}
			}
		}
		*(cigarBuf - (cigarBufLen == 0 ? 1 : 0)) = '\0'; // terminate string
		if (o_cigarBufUsed != NULL) {
			*o_cigarBufUsed = (int)(cigarBuf - cigarBufStart);
		}
		if (o_textUsed != NULL) {
			*o_textUsed = end;
		}
		return e;
	}

#ifdef TRACE_LV
	// Dump the contents of the various arrays
	printf("Done with e=%d, d=%d\n", e, d);
	for (int ee = 0; ee <= e; ee++) {
		for (int dd = -e; dd <= e; dd++) {
			if (dd >= -ee && dd <= ee)
				printf("%3d ", L[ee][MAX_K + dd]);
			else
				printf("    ");
		}
		printf("\n");
	}
	for (int ee = 0; ee <= e; ee++) {
		for (int dd = -e; dd <= e; dd++) {
			if (dd >= -ee && dd <= ee)
				printf("%3c ", A[ee][MAX_K + dd]);
			else
				printf("    ");
		}
		printf("\n");
	}
#endif

	// Trace backward to build up the CIGAR string.  We do this by filling in the backtraceAction,
	// backtraceMatched and backtraceD arrays, then going through them in the forward direction to
	// figure out our string.
	int curD = lastBestD;
	for (int curE = e; curE >= 1; curE--) {
		backtraceAction[curE] = A[curE][MAX_K + curD];
		if (backtraceAction[curE] == 'I') {
			backtraceD[curE] = curD + 1;
			backtraceMatched[curE] = L[curE][MAX_K + curD] - L[curE - 1][MAX_K + curD + 1] - 1;
		}
		else if (backtraceAction[curE] == 'D') {
			backtraceD[curE] = curD - 1;
			backtraceMatched[curE] = L[curE][MAX_K + curD] - L[curE - 1][MAX_K + curD - 1];
		}
		else { // backtraceAction[curE] == 'X'
			backtraceD[curE] = curD;
			backtraceMatched[curE] = L[curE][MAX_K + curD] - L[curE - 1][MAX_K + curD] - 1;
		}
		curD = backtraceD[curE];
#ifdef TRACE_LV
		printf("%d %d: %d %c %d %d\n", curE, curD, L[curE][MAX_K + curD],
			backtraceAction[curE], backtraceD[curE], backtraceMatched[curE]);
#endif
	}

	int accumulatedMs;	// Count of Ms that we need to emit before an I or D (or ending).
	if (useM) {
		accumulatedMs = L[0][MAX_K + 0];
	}
	else {
		// Write out ='s for the first patch of exact matches that brought us to L[0][0]
		if (L[0][MAX_K + 0] > 0) {
			validateAction(lastAction, '=');
			if (!writeCigar(&cigarBuf, &cigarBufLen, L[0][MAX_K + 0], '=', format)) {
				return -2;
			}
		}
	}

	int curE = 1;
	while (curE <= e) {
		// First write the action, possibly with a repeat if it occurred multiple times with no exact matches
		char action = backtraceAction[curE];
		int actionCount = 1;
		while (curE + 1 <= e && backtraceMatched[curE] == 0 && backtraceAction[curE + 1] == action) {
			actionCount++;
			curE++;
		}

        if (action == 'I') {
            *o_netIndel -= actionCount;
        } else if (action == 'D') {
            *o_netIndel += actionCount;
        }

		if (useM) {
			if (action == '=' || action == 'X') {
				accumulatedMs += actionCount;
			}
			else {
				if (accumulatedMs != 0) {
					validateAction(lastAction, 'M');
					if (!writeCigar(&cigarBuf, &cigarBufLen, accumulatedMs, 'M', format)) {
						return -2;
					}
					accumulatedMs = 0;
				}
				validateAction(lastAction, action);
				if (!writeCigar(&cigarBuf, &cigarBufLen, actionCount, action, format)) {
					return -2;
				}
			}
		}
		else {
			validateAction(lastAction, action);
			if (!writeCigar(&cigarBuf, &cigarBufLen, actionCount, action, format)) {
				return -2;
			}
		}
		// Next, write out ='s for the exact match
		if (backtraceMatched[curE] > 0) {
			if (useM) {
				accumulatedMs += backtraceMatched[curE];
			}
			else {
				validateAction(lastAction, '=');
				if (!writeCigar(&cigarBuf, &cigarBufLen, backtraceMatched[curE], '=', format)) {
					return -2;
				}
			}
		}
		curE++;
	}
	if (useM && accumulatedMs != 0) {
		//
		// Write out the trailing Ms.
		//
		validateAction(lastAction, 'M');
		if (!writeCigar(&cigarBuf, &cigarBufLen, accumulatedMs, 'M', format)) {
			return -2;
		}
	}
	if (format != BAM_CIGAR_OPS) {
		*(cigarBuf - (cigarBufLen == 0 ? 1 : 0)) = '\0'; // terminate string
	}
	if (o_cigarBufUsed != NULL) {
		*o_cigarBufUsed = (int)(cigarBuf - cigarBufStart);
	}
	if (o_textUsed != NULL) {
		*o_textUsed = min(textLen, lastBestBest + lastBestD);
	}
	return e; 
}

int LandauVishkinWithCigar::computeEditDistanceNormalized(
    const char* text, int textLen,
    const char* pattern, int patternLen,
    int k,
    char *cigarBuf, int cigarBufLen, bool useM,
    CigarFormat format, int* o_cigarBufUsed,
    int* o_addFrontClipping,
    int *o_netIndel)
{
    if (format != BAM_CIGAR_OPS && format != COMPACT_CIGAR_STRING) {
        WriteErrorMessage("LandauVishkinWithCigar::computeEditDistanceNormalized invalid parameter\n");
        soft_exit(1);
    }
    int bamBufLen = (format == BAM_CIGAR_OPS ? 1 : 2) * cigarBufLen; // should be enough
    char* bamBuf = (char*)alloca(bamBufLen);
    int bamBufUsed, textUsed;
    int score = computeEditDistance(text, (int)textLen, pattern, (int)patternLen, k, bamBuf, bamBufLen,
        useM, BAM_CIGAR_OPS, &bamBufUsed, &textUsed, o_netIndel);
    if (score < 0) {
        return score;
    }

    _uint32* bamOps = (_uint32*)bamBuf;
    int bamOpCount = bamBufUsed / sizeof(_uint32);

#if  0 // Not sure this is necessary, and it seems to cause problems with the new LV that won't put indels at the end
	bool hasIndels = false;
    for (int i = 0; i < bamOpCount; i++) {
        char c = BAMAlignment::CodeToCigar[BAMAlignment::GetCigarOpCode(bamOps[i])];
        if (c == 'I' || c == 'D') {
            hasIndels = true;
            break;
        }
    }
	

    if (hasIndels) {
        // run it again in reverse so it pushes indels towards the beginning
        char* text2 = (char*)alloca(textLen + 1);
        _ASSERT(textUsed <= textLen);
        util::memrevcpy(text2, text, textUsed);
        char* pattern2 = (char*)alloca(patternLen + 1);
        util::memrevcpy(pattern2, pattern, patternLen);
        char* bamBuf2 = (char*)alloca(bamBufLen);
        int bamBufUsed2, textUsed2;
        int score2 = computeEditDistance(text2, textUsed, pattern2, patternLen, k, bamBuf2, bamBufLen,
            useM, BAM_CIGAR_OPS, &bamBufUsed2, &textUsed2);
        if (score == score2 /* && bamBufUsed2 == bamBufUsed && textUsed2 == textUsed*/) {
            bamBuf = bamBuf2;
            bamBufUsed = bamBufUsed2;
            bamOpCount = bamBufUsed2 / sizeof(_uint32);
            textUsed = textUsed2;
            // reverse the operations
            for (int i = 0; i < bamOpCount; i++) {
                bamOps[i] = ((_uint32*)bamBuf2)[bamOpCount - 1 - i];
            }
        } else if (false) { // debugging
            text2[textUsed2] = 0;
            pattern2[patternLen] = 0;
            WriteErrorMessage("inconsistent forward/reverse comparison\nreverse score %d, textUsed %d, bamUsed %d, text/pattern:\n%s\n%s\n",
                score2, textUsed2, bamBufUsed2, text2, pattern2);
            memcpy(text2, text, textLen);
            text2[textLen] = 0;
            memcpy(pattern2, pattern, patternLen);
            pattern2[patternLen] = 0;
            WriteErrorMessage("forward score %d, textUsed %d, bamUsed %d, text/pattern:\n%s\n%s\n",
                score, textUsed, bamBufUsed, text2, pattern2);
        }
    }
#endif //  0 // Not sure this is necessary, and it seems to cause problems with the new LV that won't put indels at the end

#if 0 // This shouldn't happen anymore, the basic computeEditDistance doesn't allow it

    //
    // Trim out any trailing insertions, which can just be changed or merge in to X (or M as the case may be).
    //
    _ASSERT(bamOpCount > 0);
    char lastCode = BAMAlignment::CodeToCigar[BAMAlignment::GetCigarOpCode(bamOps[bamOpCount - 1])];
    if ('I' == lastCode) {
        //
        // See if it merges into the previous cigar code (which it will if it's M or X, but not D or =).
        //
        if (bamOpCount != 1) {
            char previousOp = BAMAlignment::CodeToCigar[BAMAlignment::GetCigarOpCode(bamOps[bamOpCount - 2])];
            if ('X' == previousOp || 'M' == previousOp) {
                int newCount = BAMAlignment::GetCigarOpCount(bamOps[bamOpCount - 1]) + BAMAlignment::GetCigarOpCount(bamOps[bamOpCount - 2]);
                bamOps[bamOpCount - 2] = (newCount << 4) | BAMAlignment::CigarToCode[previousOp];
                bamOpCount--;
                bamBufUsed -= sizeof(_uint32);
            } else if ('=' == previousOp) {
                //
                // The previous op was =, which this obviously doesn't.  Convert the final code to X or M.
                //
                bamOps[bamOpCount - 1] = (BAMAlignment::GetCigarOpCount(bamOps[bamOpCount - 1]) << 4) | BAMAlignment::CigarToCode[useM ? 'M' : 'X'];
            }
        }
    }

#endif // 0 // This shouldn't happen anymore, the basic computeEditDistance doesn't allow it.  Just assert it
	_ASSERT('I' != BAMAlignment::CodeToCigar[BAMAlignment::GetCigarOpCode(bamOps[bamOpCount - 1])]);

    //
    // Turn leading 'D' into soft clipping, and 'I' into an alignment change followed by an X.
    //
    if (o_addFrontClipping != NULL) {
        char firstCode = BAMAlignment::CodeToCigar[BAMAlignment::GetCigarOpCode(bamOps[0])];
		if (firstCode == 'D') {
			*o_addFrontClipping = BAMAlignment::GetCigarOpCount(bamOps[0]);
			if (*o_addFrontClipping != 0) {
				return 0; // can fail, will be rerun with new clipping
			}
		} else if (firstCode == 'I') {
			*o_addFrontClipping = -1 * BAMAlignment::GetCigarOpCount(bamOps[0]);
        } else {
            *o_addFrontClipping = 0;
        }
    }


	_ASSERT(bamOpCount <= 1 || BAMAlignment::CodeToCigar[BAMAlignment::GetCigarOpCode(bamOps[bamOpCount - 1])] != 'I');	// We should have cleared all of these out
	_ASSERT(bamOpCount <= 1 || BAMAlignment::CodeToCigar[BAMAlignment::GetCigarOpCode(bamOps[bamOpCount - 1])] != 'D');	// And none of these should happen, either.
// Seems to happen; TODO: fix this	_ASSERT(bamOpCount <= 1 || BAMAlignment::CodeToCigar[BAMAlignment::GetCigarOpCode(bamOps[0])] != 'D');
// Seems to happen; TODO: fix this	_ASSERT(bamOpCount <= 1 || BAMAlignment::CodeToCigar[BAMAlignment::GetCigarOpCode(bamOps[0])] != 'I');
	
	// copy out cigar info
    if (format == BAM_CIGAR_OPS) {
        memcpy(cigarBuf, bamOps, bamBufUsed);
        if (o_cigarBufUsed != NULL) {
            *o_cigarBufUsed = bamBufUsed;
        }
    } else {
        bool ok = BAMAlignment::decodeCigar(cigarBuf, cigarBufLen, bamOps, bamOpCount);
        if (! ok) {
            return -1;
        }
        if (o_cigarBufUsed != NULL) {
            *o_cigarBufUsed = (int)strlen(cigarBuf) + 1;
        }


	}
    return score;
}

    int
LandauVishkinWithCigar::linearizeCompactBinary(
    _uint16* o_linear,
    int referenceSize,
    char* cigar,
    int cigarSize,
    char* sample,
    int sampleSize)
{
    memset(o_linear, 0, referenceSize * 2); // zero-init
    int ic = 0, ir = 0, is = 0; // index into cigar, linear/reference, and sample
    while (ic < cigarSize) {
        int n = (unsigned char) cigar[ic++];
        char code = cigar[ic++];
        int ii, base;
        for (int i = 0; i < n; i++) {
            if ((code != 'I' && ir >= referenceSize) || (code != 'D' && is >= sampleSize)) {
                return ir;
            }
            if (code != 'D') {
                base = sample[is] != 'N' ? BASE_VALUE[sample[is]] : 0;
                is++;
            }
            switch (code) {
            case '=':
                ir++;
                break;
            case 'X':
                o_linear[ir++] |= CigarReplace + base;
                break;
            case 'D':
                o_linear[ir++] |= CigarDelete;
                break;
            case 'I':
                ii = (o_linear[ir] >> CigarInsertCShift) & CigarInsertCount;
                if (ii < 4) {
                    o_linear[ir] = (base << (2 * ii + CigarInsertBShift)) | ((ii + 1) << CigarInsertCShift);
                } else if (ii < 7) {
                    o_linear[ir] = (o_linear[ir] & CigarInsertBases) | ((ii + 1) << CigarInsertCShift);
                }
                break;
            default:
                _ASSERT(false);
            }
        }
    }
    return ir;
}

    void 
setLVProbabilities(double *i_indelProbabilities, double *i_phredToProbability, double mutationProbability)
{
    lv_indelProbabilities = i_indelProbabilities;

    //
    // Compute the phred table to incorporate the mutation probability, assuming that machine errors and mutations
    // are independent (which there's no reason not to think is the case).  If P(A) and P(B) are independent, then
    // P(A or B) = P(not (not-A and not-B)) = 1-(1-P(A))(1-P(B)).
    //
    for (unsigned i = 0; i < 255; i++) {
        lv_phredToProbability[i] = 1.0-(1.0 - i_phredToProbability[i]) * (1.0 - mutationProbability);
    }
}

    void
initializeLVProbabilitiesToPhredPlus33()
{
    static bool alreadyInitialized = false;
    if (alreadyInitialized) {
        return;
    }
    alreadyInitialized = true;

    //
    // indel probability is .0001 for any indel (10% of a SNP real difference), and then 10% worse for each longer base.
    //
    _ASSERT(NULL == lv_phredToProbability);
    lv_phredToProbability = (double *)BigAlloc(sizeof(double) * 256);

    static const int maxIndels = 10000; // Way more than we'll see, and in practice enough to result in p=0.0;
    _ASSERT(NULL == lv_indelProbabilities);
    lv_indelProbabilities = (double *)BigAlloc(sizeof(double) * maxIndels);
 
    const double mutationRate = SNP_PROB;
    lv_indelProbabilities = new double[maxIndels+1];
    lv_indelProbabilities[0] = 1.0;
    lv_indelProbabilities[1] = GAP_OPEN_PROB;
    for (int i = 2; i <= maxIndels; i++) {
        lv_indelProbabilities[i] = lv_indelProbabilities[i-1] * GAP_EXTEND_PROB;
    }

    //
    // Use 0.001 as the probability of a real SNP, then or it with the Phred+33 probability.
    //
    for (int i = 0; i < 33; i++) {
        lv_phredToProbability[i] = mutationRate;  // This isn't a sensible Phred score
    }
    for (int i = 33; i <= 93 + 33; i++) {
         lv_phredToProbability[i] = 1.0-(1.0 - pow(10.0,-1.0 * (i - 33.0) / 10.0)) * (1.0 - mutationRate);
    }
    for (int i = 93 + 33 + 1; i < 256; i++) {
        lv_phredToProbability[i] = mutationRate;   // This isn't a sensible Phred score
    }

    _ASSERT(NULL == lv_perfectMatchProbability);
    lv_perfectMatchProbability = new double[MaxReadLength+1];
    lv_perfectMatchProbability[0] = 1.0;
    for (unsigned i = 1; i <= MaxReadLength; i++) {
        lv_perfectMatchProbability[i] = lv_perfectMatchProbability[i - 1] * (1 - SNP_PROB);
    }

    initializeMapqTables();
}

double *lv_phredToProbability = NULL;
double *lv_indelProbabilities = NULL;
double *lv_perfectMatchProbability = NULL;
