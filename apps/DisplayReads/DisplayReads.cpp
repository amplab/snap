#include "stdafx.h"
#include "BigAlloc.h"
#include "exit.h"
#include "Compat.h"
#include "Genome.h"
#include "Read.h"
#include "Bam.h"
#include "SAM.h"

void
usage()
{
    fprintf(stderr, "usage: DisplayReads index inputSamFileName contig/chromosome startAddress endAddress\n");
    fprintf(stderr, "Start and end address refer to the aligned address of the read, so the display will extend farther.\n");
    soft_exit(1);
}

bool
BreakLineIntoPiecesBasedOnTabs(char *inputBuffer, char **piece, int nPieces) {
    //
    // Parse the line into tab separated pieces.
    //

    piece[0] = inputBuffer;

    for (int whichPiece = 1; whichPiece < nPieces; whichPiece++) {
        piece[whichPiece] = strchr(piece[whichPiece - 1], '\t');
        if (NULL == piece[whichPiece]) {
            return false;
        }
        *piece[whichPiece] = '\0';  // Turn the tab into a null, so we have separate strings for each piece
        piece[whichPiece]++; // And point at the real beginning of the string.
    }

    return true;
}
enum CigarOpType { Equal, Unequal, Match, Insert, Delete, SoftClip, HardClip, N };

struct CigarOp {
    CigarOpType opType;
    int         count;
    int         pos;    // Of the beginning of the region described by this cigar op
};

const size_t maxCigarOps = 50;

struct ReadListEntry {
    char *bases;
    char *cigar;
    int pos;
    CigarOp cigarOps[maxCigarOps];
    size_t nCigarOps;

    ReadListEntry *next;
};


int compareReadListEntriesByPos(const void *r1, const void *r2)
{
    ReadListEntry *entry1 = *(ReadListEntry **)r1;
    ReadListEntry *entry2 = *(ReadListEntry **)r2;

    if (entry1->pos > entry2->pos) {
        return 1;
    }
    else if (entry1->pos == entry2->pos) {
        return 0;
    }
    else {
        return -1;
    }
}



bool
parseCigarString(const char *cigarString, size_t maxOps, CigarOp *ops, size_t *nOps, int alignedPos, int *finalPos)
{
    const char *curPtr = cigarString;
    int currentPos = alignedPos;
    for (size_t i = 0; i < maxOps; i++) {
        //
        // Cigar ops are of the form #OpType where # is, well, a number, and opType is a single character.
        //
        if (*curPtr == '\0') {
            *nOps = i;
            *finalPos = currentPos;
            return (i != 0);    // We're OK iff we're not the empty string
        }
        const size_t maxDigits = 10;
        char digits[maxDigits + 1];
        size_t whichDigit;
        for (whichDigit = 0; whichDigit < maxDigits && *curPtr >= '0' && *curPtr <= '9'; whichDigit++) {
            digits[whichDigit] = *curPtr;
            curPtr++;
        }

        if (whichDigit == maxDigits) return false;
        digits[whichDigit] = '\0';

        ops[i].count = atoi(digits);
        if (ops[i].count <= 0) return false;
        ops[i].pos = currentPos;

        switch (*curPtr) {
            case '=': ops[i].opType = Equal;        currentPos += ops[i].count; break;
            case 'X': ops[i].opType = Unequal;      currentPos += ops[i].count; break;
            case 'M': ops[i].opType = Match;        currentPos += ops[i].count; break;
            case 'I': ops[i].opType = Insert;                                   break;
            case 'D': ops[i].opType = Delete;       currentPos += ops[i].count; break;
            case 'S': ops[i].opType = SoftClip;                                 break;
            case 'H': ops[i].opType = HardClip;                                 break;
            case 'N': ops[i].opType = N;            currentPos += ops[i].count; break;

            default: return false;
        }
        curPtr++;
    }

    return false;
}

int main(int argc, char * argv[])
{
    BigAllocUseHugePages = false;

    if (argc != 6) usage();

    static const char *genomeSuffix = "Genome";
    size_t filenameLen = strlen(argv[1]) + 1 + strlen(genomeSuffix) + 1;
    char *fileName = new char[strlen(argv[1]) + 1 + strlen(genomeSuffix) + 1];
    snprintf(fileName, filenameLen, "%s%c%s", argv[1], PATH_SEP, genomeSuffix);
    const Genome *genome = Genome::loadFromFile(fileName, 0);
    if (NULL == genome) {
        fprintf(stderr, "Unable to load genome from file '%s'\n", fileName);
        return -1;
    }
    delete[] fileName;
    fileName = NULL;

    //
    // Now translate the start and end addresses into GenomeLocations.
    //
    GenomeLocation contigBase;
    if (!genome->getLocationOfContig(argv[3], &contigBase)) {
        WriteErrorMessage("Unable to find contig '%s'\n", argv[3]);
        soft_exit(1);
    }

    int startPos = atoi(argv[4]);
    int endPos = atoi(argv[5]);
    if (startPos >= endPos) {
        WriteErrorMessage("Start location must be less than end location\n");
        soft_exit(1);
    }

    const char *inputFileName = argv[2];
    FILE *inputFile = fopen(inputFileName, "r");
    if (NULL == inputFile) {
        fprintf(stderr, "Unable to open input file %s\n", inputFileName);
        soft_exit(1);
    }

    //
    // Run through all of the reads and save those that fit within our range.
    //

    int nSavedReads = 0;
    ReadListEntry *savedReads = NULL;

    const size_t maxCigarOps = 50;
    const size_t inputBufferSize = 40000;
    char inputBuffer[inputBufferSize];
    int maxPosForAnyBase = 0;
    int minPosForAnyBase = 0x7fffffff;


    while (fgets(inputBuffer, inputBufferSize - 1, inputFile)) {
        if ('@' == inputBuffer[0]) {
            //
            // Ignore the header lines.
            //
            continue;
        }

        const int Flag = 1;
        const int RName = 2;
        const int Pos = 3;
        const int Cigar = 5;
        const int Seq = 9;
        const int nPieces = 10;

        char *piece[nPieces];
        if (!BreakLineIntoPiecesBasedOnTabs(inputBuffer, piece, nPieces)) {
            continue;   // Badly formatted
        }

        int flag = atoi(piece[Flag]);
        if (flag & SAM_UNMAPPED) {
            continue;
        }

        if (strcmp(piece[RName], argv[3])) {
            continue;   // Wrong chromosome
        }

        ReadListEntry *entry = new ReadListEntry;
        int finalPos;
        int pos = atoi(piece[Pos]);
        if (pos < startPos || pos > endPos || !parseCigarString(piece[Cigar], maxCigarOps, entry->cigarOps, &entry->nCigarOps, pos, &finalPos)) {
            delete entry;
            continue;
        }

        maxPosForAnyBase = __max(maxPosForAnyBase, finalPos);
        minPosForAnyBase = __min(minPosForAnyBase, pos);

        entry->bases = new char[strlen(piece[Seq]) + 1];
        strcpy(entry->bases, piece[Seq]);
        entry->cigar = new char[strlen(piece[Cigar]) + 1];
        strcpy(entry->cigar, piece[Cigar]);
        entry->pos = pos;

        entry->next = savedReads;
        savedReads = entry;
        nSavedReads++;
    }

    if (0 == nSavedReads) {
        printf("No reads aligned to the specified range\n");
        return 0;
    }

    ReadListEntry **sortedReads = new ReadListEntry *[nSavedReads];

    for (int i = 0; i < nSavedReads; i++) {
        sortedReads[i] = savedReads;
        savedReads = savedReads->next;  
    }

    qsort(sortedReads, nSavedReads, sizeof(ReadListEntry *), compareReadListEntriesByPos);

    //
    // In order to display the reads properly, we need to add columns where there are insertions.  So, we build an array that holds
    // the maximum insertion count after any base.
    //

    unsigned *maxInsertionsAfterBase = new unsigned[maxPosForAnyBase - minPosForAnyBase + 1];
    for (GenomeDistance i = 0; i <maxPosForAnyBase - minPosForAnyBase + 1; i++) {
        maxInsertionsAfterBase[i] = 0;
    }

    for (int i = 0; i < nSavedReads; i++) {
        ReadListEntry *entry = sortedReads[i];
        for (size_t op = 0; op < entry->nCigarOps; op++) {
            if (entry->cigarOps[op].opType == Insert) {
                maxInsertionsAfterBase[entry->cigarOps[op].pos - minPosForAnyBase] = __max(maxInsertionsAfterBase[entry->cigarOps[op].pos], entry->cigarOps[op].count);
            }
        }
    }

    //
    // Print out the header with locations (and insertions).  We do three rows: hundreds, tens and units.
    //

    for (int radix = 100; radix > 0; radix /= 10) {
        printf("%d", (minPosForAnyBase / radix) % 10);    // Always print in the first column.
        for (int pos = minPosForAnyBase + 1; pos <= maxPosForAnyBase; pos++) {
            if (pos % radix == 0) {
                printf("%d", (pos / radix) % 10);
            }
            else {
                printf(" ");
            }
            for (unsigned i = 0; i < maxInsertionsAfterBase[pos - minPosForAnyBase]; i++) {
                if (1 == radix) {
                    printf("I");
                }
                else {
                    printf(" ");
                }
            }
        }
        printf("\n");
    }

    printf("\n");

    //
    // Finally....print out the reads.
    //
    for (int i = 0; i < nSavedReads; i++) {
        ReadListEntry *entry = sortedReads[i];

        //
        // Print spaces to get to to the first of its bases.
        //
        for (int pos = minPosForAnyBase; pos < entry->pos; pos++) {
            for (unsigned j = 0; j <= maxInsertionsAfterBase[pos - minPosForAnyBase]; j++) { // The <= makes it print one space for the normal (non-insertion) base
                printf(" ");
            }
        }

        int offsetInSeq = 0;
        int pos = entry->cigarOps[0].pos;
        int nBasesInsertedAtPriorPos = 0;

        for (size_t op = 0; op < entry->nCigarOps; op++) {
            switch (entry->cigarOps[op].opType) {
                case Match:
                case Equal:
                case Unequal:
                    if (pos > minPosForAnyBase) {
                        for (unsigned j = nBasesInsertedAtPriorPos; j < maxInsertionsAfterBase[pos - 1 - minPosForAnyBase]; j++) {
                            printf(" ");
                        }
                    }

                    for (int j = 0; j < entry->cigarOps[op].count; j ++ ) {
                        const char *refBase = genome->getSubstring(contigBase + pos - 1, 1);    // contigBase + pos - 1 has -1 because pos is 1-based.
                        if (NULL == refBase) {
                            printf("X");
                        }
                        else {
                            printf("%c", (*refBase == entry->bases[offsetInSeq]) ? '.' : entry->bases[offsetInSeq]);
                        }
                        pos++;
                        offsetInSeq++;
                    }

                    nBasesInsertedAtPriorPos = 0;
                    break;

                case SoftClip:
                    //
                    // Just skip the bases in the read.
                    //
                    offsetInSeq += entry->cigarOps[op].count;
                    break;
                    
                case HardClip:
                    //
                    // The bases aren't even in the read.  Just ignore this entirely.
                    //
                    break;
                    
                case Insert:
                    //
                    // Just print out the bases.
                    //
                    for (int j = 0; j < entry->cigarOps[op].count; j++) {
                        printf("%c", entry->bases[offsetInSeq]);
                        offsetInSeq++;
                    }
                    nBasesInsertedAtPriorPos = entry->cigarOps[op].count;
                    break;

                case Delete:
                    //
                    // Just skip over the pos.
                    //
                    for (int j = 0; j < entry->cigarOps[op].count; j++) {
                        printf(" ");
                     }
                    pos += entry->cigarOps[op].count;
                    nBasesInsertedAtPriorPos = 0;
                    break;

                case N:
                    //
                    // This probably needs to be taken into account in the basic formatting.  Just ignore it, which will doubtless produce garbage.
                    //
                    break;

            } // switch
        } // for each cigar Op
        printf("\n");
    } // for each read



    return 0;
}

