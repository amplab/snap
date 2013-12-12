#include "stdafx.h"
#include "Compat.h"
#include "LandauVishkin.h"
#include "mapq.h"
#include "Read.h"
#include "BaseAligner.h"
#include "Bam.h"
#include "exit.h"
#include <vector>

using std::make_pair;
using std::min;

LandauVishkinWithCigar::LandauVishkinWithCigar()
{
    for (int i = 0; i < MAX_K+1; i++) {
        for (int j = 0; j < 2*MAX_K+1; j++) {
            L[i][j] = -2;
        }
    }
}

/*++
    Write cigar to buffer, return true if it fits
    null-terminates buffer if it returns false (i.e. fills up buffer)
--*/
bool writeCigar(char** o_buf, int* o_buflen, int count, char code, CigarFormat format, std::vector<unsigned> &tokens)
{

    tokens.push_back(count);
    tokens.push_back(code);

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
        printf("invalid cigar format %d\n", format);
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

int LandauVishkinWithCigar::insertSpliceJunctions(
    const GTFReader *gtf,
    std::vector<unsigned> tokens, 
    std::string transcript_id, 
    unsigned pos,
    char *cigarNew, 
    int cigarNewLen, 
    CigarFormat format)
{

    //Create array to store new tokens
    std::vector<unsigned> final;

    //Keep pointer to beginning
    char* cigarBufStart = cigarNew;
    
    //If the length of tokens is odd, something is wrong
    if (tokens.size() % 2 != 0) {
        fprintf(stderr, "Error: Size of CIGAR tokens is odd (%d)\n", tokens.size());
        exit(1);
    }

    //printf("%d %s %u\n", gtf->Size(), transcript_id.c_str(), pos);  
    unsigned prev = pos;
    unsigned current = pos;
    std::string new_cigar;
    int numCigarOps = 0;

    //Loop over all the operators
    for (std::vector<unsigned>::iterator it = tokens.begin(); it != tokens.end(); it+=2) {
       
        //Get the size of the operator and the operator character
        unsigned length = *it;
        char op = *(it+1);
   
        //printf("%u %c\n", length, op);
   
        if (op == 'I') {
            if (! writeCigar(&cigarNew, &cigarNewLen, length, op, format, final)) {
                return -2;
            }
            numCigarOps++;
        }
        else if (op == 'S') {
            if (! writeCigar(&cigarNew, &cigarNewLen, length, op, format, final)) {
                return -2;
            }
            numCigarOps++;
        }
        else {
      
            //Get the end of this operator
            current += length-1;
            
            //Get the junctions
            //printf("Querying with pos: %d length: %d\n", prev, length);
            std::vector<junction> junctions;
            gtf->GetTranscript(transcript_id).Junctions(prev, length, junctions);

            //If this operator crosses a splice junction, we must add it into the CIGAR operator
            if (junctions.size() > 0) {
            
                //printf("Junctions: %d Remainder: %d\n", junctions.size(), length);
            
                // Add any junctions to the CIGAR operator
                unsigned remainder = length;
                          
                //Iterate over the junctions
                for (std::vector<junction>::iterator jit = junctions.begin(); jit != junctions.end(); ++jit) {
                                    
                    // Special case 1: the read begins right on the splice
                    // junction: in this case we do not insert the junction
                    if (jit->first == pos) {
                        continue; 
                    }
                                            
                    // Junction[0] is the position in transcript space where the
                    // splice junction begins (end of exon + 1)
                    int step = jit->first - prev;
                    remainder = remainder - step;
                    
                    //printf("%u %u %u %d %u\n", jit->first, jit->second, prev, step, remainder);
                              
                    // Add in portion of operator
                    if (step > 0) {
                        if (! writeCigar(&cigarNew, &cigarNewLen, step, op, format, final)) {
                            return -2;
                        }
                        numCigarOps++;
                    }   
                                                
                    // Add in splice junction spacer only if this is not the first operator
                    if (! writeCigar(&cigarNew, &cigarNewLen, jit->second->Length(), 'N', format, final)) {
                        return -2;
                    }
                    numCigarOps++;
                                               
                    // Increment the amount we move through the transcript
                    prev = prev + step;  
                    
                }                      
                        
                // Add in remainder of current operator
                if (remainder > 0) {
                    if (! writeCigar(&cigarNew, &cigarNewLen, remainder, op, format, final)) {
                        return -2;
                    }
                    numCigarOps++;
                }

        
            // If there are no splice junctions, just print the operator as-is
            } else {            
                if (! writeCigar(&cigarNew, &cigarNewLen, length, op, format, final)) {
                    return -2;
                }
                numCigarOps++;
            }
            
            // We begin at the next position for the next operator
            current = current + 1;           
            prev = current;
        }
    }

    if (format != BAM_CIGAR_OPS) {
        *(cigarNew - (cigarNewLen == 0 ? 1 : 0)) = '\0'; // terminate string
    }

    return numCigarOps;

}

static const int PrevDelta[3][3] = 
    {{+1, 0, -1},    // d < 0
    {0, -1, +1},    // d == 0
    {-1, 0, +1}};   // d > 0

int LandauVishkinWithCigar::computeEditDistance(
    const char* text, int textLen,
    const char* pattern, int patternLen,
    int k,
    char *cigarBuf, int cigarBufLen, bool useM,
    std::vector<unsigned> &tokens,
    CigarFormat format, int* o_cigarBufUsed, int* o_textUsed)
{
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
			if (! writeCigar(&cigarBuf, &cigarBufLen, patternLen, 'M', format, tokens)) {
				return -2;
			}
            // todo: should this also write X's like '=' case? or is 'M' special?
		} else {
			if (! writeCigar(&cigarBuf, &cigarBufLen, end, '=', format, tokens)) {
				return -2;
			}
			if (patternLen > end) {
				// Also need to write a bunch of X's past the end of the text
				if (! writeCigar(&cigarBuf, &cigarBufLen, patternLen - end, 'X', format, tokens)) {
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

    for (int e = 1; e <= k; e++) {
        // Go through the offsets, d, in the order 0, -1, 1, -2, 2, etc, in order to find CIGAR strings
        // with few indels first if possible.
        for (int d = 0; d != -(e+1); d = (d >= 0 ? -(d+1) : -d)) {
            int bestdelta = 0;
            int bestbest = -1;
            //  extend previous solutions as far as possible, pick best, minimizing indels
            int dy = (d >= 0) + (d > 0);
            for (int dx = 0; dx < 3; dx++) {
                int delta = PrevDelta[dy][dx];
                int best = L[e-1][MAX_K+d + delta] + (delta >= 0);
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
                if (best > bestbest) {
                    bestbest = best;
                    bestdelta = delta;
                }
            }
            int best = bestbest;
            A[e][MAX_K+d] = "DXI"[bestdelta + 1];

            L[e][MAX_K+d] = best;

            if (best == patternLen) {
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
						if (!writeCigar(&cigarBuf, &cigarBufLen, patternLen, 'M', format, tokens)) {
							return -2;
						}
					} else {
						int streakStart = 0;
						bool matching = (pattern[0] == text[0]);
						for (int i = 0; i < end; i++) {
							bool newMatching = (pattern[i] == text[i]);
							if (newMatching != matching) {
                                validateAction(lastAction, matching ? '=' : 'X');
								if (!writeCigar(&cigarBuf, &cigarBufLen, i - streakStart, (matching ? '=' : 'X'), format, tokens)) {
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
								if (!writeCigar(&cigarBuf, &cigarBufLen, patternLen - streakStart, 'X', format, tokens)) {
									return -2;
								}
							} else {
								// Write out some ='s and then possibly X's if pattern is longer than text
                                validateAction(lastAction, '=');
								if (!writeCigar(&cigarBuf, &cigarBufLen, end - streakStart, '=', format, tokens)) {
									return -2;
								}
								if (patternLen > end) {
                                    validateAction(lastAction, 'X');
									if (!writeCigar(&cigarBuf, &cigarBufLen, patternLen - end, 'X', format, tokens)) {
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
                            printf("%3d ", L[ee][MAX_K+dd]);
                        else
                            printf("    ");
                    }
                    printf("\n");
                }
                for (int ee = 0; ee <= e; ee++) {
                    for (int dd = -e; dd <= e; dd++) {
                        if (dd >= -ee && dd <= ee)
                            printf("%3c ", A[ee][MAX_K+dd]);
                        else
                            printf("    ");
                    }
                    printf("\n");
                }
#endif

                // Trace backward to build up the CIGAR string.  We do this by filling in the backtraceAction,
                // backtraceMatched and backtraceD arrays, then going through them in the forward direction to
                // figure out our string.
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

				int accumulatedMs;	// Count of Ms that we need to emit before an I or D (or ending).
				if (useM) {
					accumulatedMs = L[0][MAX_K+0];
				} else {
					// Write out ='s for the first patch of exact matches that brought us to L[0][0]
					if (L[0][MAX_K+0] > 0) {
                        validateAction(lastAction, '=');
						if (! writeCigar(&cigarBuf, &cigarBufLen, L[0][MAX_K+0], '=', format, tokens)) {
							return -2;
						}
					}
				}

                int curE = 1;
                while (curE <= e) {
                    // First write the action, possibly with a repeat if it occurred multiple times with no exact matches
                    char action = backtraceAction[curE];
                    int actionCount = 1;
                    while (curE+1 <= e && backtraceMatched[curE] == 0 && backtraceAction[curE+1] == action) {
                        actionCount++;
                        curE++;
                    }
					if (useM) {
						if (action == '=' || action == 'X') {
							accumulatedMs += actionCount;
						} else {
							if (accumulatedMs != 0) {
                                validateAction(lastAction, 'M');
								if (!writeCigar(&cigarBuf, &cigarBufLen, accumulatedMs, 'M', format, tokens)) {
									return -2;
								}
								accumulatedMs = 0;
							}
                            validateAction(lastAction, action);
						if (!writeCigar(&cigarBuf, &cigarBufLen, actionCount, action, format, tokens)) {
								return -2;
							}
						}
					} else {
                        validateAction(lastAction, action);
					if (! writeCigar(&cigarBuf, &cigarBufLen, actionCount, action, format, tokens)) {
							return -2;
						}
					}
                    // Next, write out ='s for the exact match
                    if (backtraceMatched[curE] > 0) {
						if (useM) {
							accumulatedMs += backtraceMatched[curE];
						} else {
                            validateAction(lastAction, '=');
			    if (! writeCigar(&cigarBuf, &cigarBufLen, backtraceMatched[curE], '=', format, tokens)) {
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
					if (!writeCigar(&cigarBuf, &cigarBufLen, accumulatedMs, 'M', format, tokens)) {
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
                    *o_textUsed = min(textLen, best + d);
                }
                return e;
            }
        }
    }

    // Could not align strings with at most K edits
    *(cigarBuf - (cigarBufLen == 0 ? 1 : 0)) = '\0'; // terminate string
    return -1;
}

int LandauVishkinWithCigar::computeEditDistanceNormalized(
    const char* text, int textLen,
    const char* pattern, int patternLen,
    int k,
    char *cigarBuf, int cigarBufLen, bool useM, 
    std::vector<unsigned> &tokens,
    CigarFormat format, int* cigarBufUsed)
{
    if (format != BAM_CIGAR_OPS && format != COMPACT_CIGAR_STRING) {
        fprintf(stderr, "LandauVishkinWithCigar::computeEditDistanceNormalized invalid parameter\n");
        soft_exit(1);
    }
    int bamBufLen = (format == BAM_CIGAR_OPS ? 1 : 2) * cigarBufLen; // should be enough
    char* bamBuf = (char*) alloca(bamBufLen);
    int bamBufUsed, textUsed;
    int score = computeEditDistance(text, textLen, pattern, patternLen, k, bamBuf, bamBufLen,
        useM, tokens, BAM_CIGAR_OPS, &bamBufUsed, &textUsed);
    if (score < 0) {
        return score;
    }
    _uint32* bamOps = (_uint32*) bamBuf;
    int bamOpCount = bamBufUsed / sizeof(_uint32);
    bool hasIndels = false;
    for (int i = 0; i < bamOpCount && ! hasIndels; i++) {
        char c = BAMAlignment::CodeToCigar[BAMAlignment::GetCigarOpCode(bamOps[i])];
        hasIndels = c == 'I' || c == 'D';
    }
    if (hasIndels) {
        // run it again in reverse so it pushes indels towards the beginning
        char* text2 = (char*) alloca(textLen + 1);
        _ASSERT(textUsed <= textLen);
        util::memrevcpy(text2, text, textUsed);
        char* pattern2 = (char*) alloca(patternLen + 1);
        util::memrevcpy(pattern2, pattern, patternLen);
        char* bamBuf2 = (char*) alloca(bamBufLen);
        int bamBufUsed2, textUsed2;
        int score2 = computeEditDistance(text2, textUsed, pattern2, patternLen, k, bamBuf2, bamBufLen,
            useM, tokens, BAM_CIGAR_OPS, &bamBufUsed2, &textUsed2);
        if (score == score2 /* && bamBufUsed2 == bamBufUsed && textUsed2 == textUsed*/) {
            bamBuf = bamBuf2;
            bamBufUsed = bamBufUsed2;
            bamOpCount = bamBufUsed2 / sizeof(_uint32);
            textUsed = textUsed2;
            // reverse the operations
            for (int i = 0; i < bamOpCount; i++) {
                bamOps[i] = ((_uint32*) bamBuf2)[bamOpCount - 1 - i];
            }
        } else if (false) { // debugging
            text2[textUsed2] = 0;
            pattern2[patternLen] = 0;
            printf("inconsistent forward/reverse comparison\nreverse score %d, textUsed %d, bamUsed %d, text/pattern:\n%s\n%s\n",
                score2, textUsed2, bamBufUsed2, text2, pattern2);
            memcpy(text2, text, textLen);
            text2[textLen] = 0;
            memcpy(pattern2, pattern, patternLen);
            pattern2[patternLen] = 0;
            printf("forward score %d, textUsed %d, bamUsed %d, text/pattern:\n%s\n%s\n",
                score, textUsed, bamBufUsed, text2, pattern2);
        }
    }
    // copy out cigar info
    if (format == BAM_CIGAR_OPS) {
        memcpy(cigarBuf, bamOps, bamBufUsed);
        if (cigarBufUsed != NULL) {
            *cigarBufUsed = bamBufUsed;
        }
    } else {
        bool ok = BAMAlignment::decodeCigar(cigarBuf, cigarBufLen, bamOps, bamOpCount);
        if (! ok) {
            return -1;
        }
        if (cigarBufUsed != NULL) {
            *cigarBufUsed = (int) strlen(cigarBuf) + 1;
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
