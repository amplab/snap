#include "stdafx.h"

#include "GenomeIndex.h"
#include "exit.h"
#include "SAM.h"
#include "FixedSizeMap.h"


const Genome *ReferenceGenome = NULL;

void usage()
{
    fprintf(stderr, "usage: mutant-expression index inputFile {-h}\n");
    fprintf(stderr, "inputFile is tab separated text like the TCGA LAML spreadsheet, with the last column (AM) specifying a \n");
    fprintf(stderr, "line number that can be converted into a filename containing all relevant reads\n");
    fprintf(stderr, "-h means to print the header\n");
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

const char *translateAlternateChromName(const char *inputName)
{
    if (!strcmp(inputName, "GL000191.1")) return "chr1_gl000191_random";
    if (!strcmp(inputName, "GL000192.1")) return "chr1_gl000192_random";
    if (!strcmp(inputName, "GL000193.1")) return "chr4_gl000193_random";
    if (!strcmp(inputName, "GL000194.1")) return "chr4_gl000194_random";
    if (!strcmp(inputName, "GL000195.1")) return "chr7_gl000195_random";
    if (!strcmp(inputName, "GL000205.1")) return "chr17_gl000205_random";
    if (!strcmp(inputName, "GL000209.1")) return "chr19_gl000209_random";
    if (!strcmp(inputName, "GL000211.1")) return "chrUn_gl000211";
    if (!strcmp(inputName, "GL000212.1")) return "chrUn_gl000212";
    if (!strcmp(inputName, "GL000213.1")) return "chrUn_gl000213";
    if (!strcmp(inputName, "GL000214.1")) return "chrUn_gl000214";
    if (!strcmp(inputName, "GL000218.1")) return "chrUn_gl000218";
    if (!strcmp(inputName, "GL000219.1")) return "chrUn_gl000219";
    if (!strcmp(inputName, "GL000220.1")) return "chrUn_gl000220";
    if (!strcmp(inputName, "GL000229.1")) return "chrUn_gl000229";
    if (!strcmp(inputName, "GL000237.1")) return "chrUn_gl000237";
    if (!strcmp(inputName, "GL000241.1")) return "chrUn_gl000241";


    return NULL;
}

template<typename T>
class StringHash
{
public:
    inline _uint64 operator() (T string) {
        _uint64 value = 0;

        size_t len = strlen(string);
        for (int i = 0; i < len; i++) {
            value = value * 131 + string[i];
        }

        return value;
    }
};

struct RefSeqList {
    char *refSeqId;
    RefSeqList *next;
};

FixedSizeMap<const char *, RefSeqList *, StringHash<const char *> > *HugoToRefSeq = NULL;

int main(int argc, char* argv[])
{
    if (3 != argc && 4 != argc) usage();
    int lineNumber = 0;

    bool printHeader = false;
    if (4 == argc) {
        if (strcmp(argv[3], "-h")) {
            usage();
        }
        printHeader = true;
    }

#if 0
    HugoToRefSeq = new FixedSizeMap<const char *, RefSeqList *, StringHash<const char *> >(131072);

    const char *hugoFileName = "f:\\sequence\\gene_data\\protein-coding_gene.txt";
    FILE *hugoFile = fopen(hugoFileName, "r");
    if (NULL == hugoFile) {
        fprintf(stderr, "Unable to open hugo file '%s'\n", hugoFileName);
        soft_exit(1);
    }
    const size_t hugoBufferSize = 2000;
    char hugoBuffer[hugoBufferSize];
    fscanf(hugoFile, "%[^\n]\n", hugoBuffer);    // Skip the header line.  And yes, it's a buffer overflow.  This is only for friendly use.


    while (1 == fscanf(hugoFile, "%[^\n]\n", hugoBuffer)) {
        lineNumber++;
        const unsigned nPieces = 10;
        char *pieces[nPieces];
        if (!BreakLineIntoPiecesBasedOnTabs(hugoBuffer, pieces, nPieces)) {
            fprintf(stderr, "Unable to parse hugo file line number %d, '%s'\n", lineNumber, hugoBuffer);
            continue;
        }

        char *hugoSymbol = new char[strlen(pieces[1]) + 1];
        strcpy(hugoSymbol, pieces[1]);

        //
        // Break up the list of RefSeq IDs.
        //
        const unsigned maxRefSeqIds = 20;
        char *refSeqIds[maxRefSeqIds];
        unsigned nRefSeqIds = 1;
        refSeqIds[0] = pieces[9];
        char *comma;
        while (NULL != (comma = strchr(refSeqIds[nRefSeqIds - 1], ','))) {
            refSeqIds[nRefSeqIds] = comma + 1; // Yes, another buffer overflow.  Sue me.
            while (refSeqIds[nRefSeqIds][0] == ',' || refSeqIds[nRefSeqIds][0] == ' ') {
                refSeqIds[nRefSeqIds]++;
            }
            nRefSeqIds++;

            *comma = '\0';
        }

        RefSeqList *prevList = NULL;
        RefSeqList *firstList = NULL;
        for (unsigned i = 0; i < nRefSeqIds; i++) {
            RefSeqList *list = new RefSeqList;
            if (NULL == firstList) {
                firstList = list;
            }
            list->refSeqId = new char[strlen(refSeqIds[i]) + 1];
            strcpy(list->refSeqId, refSeqIds[i]);
            list->next = NULL;
            if (NULL != prevList) {
                prevList->next = list;
            }
            prevList = list;
        }
        

        HugoToRefSeq->put(hugoSymbol, firstList);
    }

#endif // 0

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


    if (printHeader) {
        printf("Cancer_Type\tLocal_file_name\tHugo_Symbol\tEntrez_Gene_Id\tCenter\tNCBI_Build\tChromosome\tStart_Position\tEnd_Position\tStrand\tVariant_Classification\tVariant_Type\tReference_Allele\tTumor_Seq_Allele_1\tTumor_Seq_Allele_2\tdbSNP_RS\tdbSNP_Val_Status\t"
            "Tumor_Sample_Barcode\tMatched_Norm_Sample_Barcode\tMatch_Norm_Seq_Allele1\tMatch_Norm_Seq_Allele2\tTumor_Validation_Allele1\tTumor_Validation_Allele2\tMatch_Norm_Validation_Allele1\tMatch_Norm_Validation_Allele2\tVerification_Status\t"
            "Validation_Status\tMutation_Status\tSequencing_Phase\tSequence_Source\tValidation_Method\tScore\tBAM_File\tSequencer\tTumor_Sample_UUID\tMatched_Norm_Sample_UUID\tFile_Name\tArchive_Name\tLine_Number\t"
            "n_Matching_Reference\tn_Matching_Tumor\tn_Matching_Neither\tn_Matching_Both\n");
    }

    ReferenceGenome = genome;

    FILE *inputFile = fopen(argv[2], "r");

    if (NULL == inputFile) {
        fprintf(stderr, "Unable to open input file '%s'\n", argv[1]);
        soft_exit(1);
    }

    const size_t inputLineSize = 20000; // We're among friends, don't worry about buffer overflows

    char inputBuffer[inputLineSize];
    char rawInputBuffer[inputLineSize];
    //fgets(inputBuffer, inputLineSize - 1, inputFile);   // Skip the header

    const int nInputLinePieces = 38;
    char *piece[nInputLinePieces];
    FILE *samFile = NULL;

    int mutationNumber = 1; // Because excel starts at 1

    int noSAMFile = 0;  
    int nMT = 0;
    int nWrongClass = 0;
    int nValidationProblems = 0;
    int nSAMProblems = 0;
    lineNumber = 0;

    while (fgets(inputBuffer, inputLineSize - 1, inputFile)) {
        int nMatchingReference = 0;
        int nMatchingMutation = 0;
        int nMatchingNeither = 0;
        int nMatchingBoth = 0;

        memcpy(rawInputBuffer, inputBuffer, inputLineSize);

        if (!BreakLineIntoPiecesBasedOnTabs(inputBuffer, piece, nInputLinePieces)) {
            fprintf(stderr, "Malformed line '%s', number %d ignoring.\n", inputBuffer, mutationNumber);
            goto doneWithMutation;
        }

        const int RNAFileName = 0;
        const int GeneName = 1;
        const int Chrom = 5;
        const int StartPosition = 6;
        const int EndPosition = 7;
        const int VariantClass = 9;
        const int VariantType = 10;
        const int Reference_Allele = 11;
        const int Tumor_Seq_Allele_1 = 12;
        const int Tumor_Seq_Allele_2 = 13;
        const int LineNumber = 37;

        //
        // Now process the mutation.
        //
//        int lineNumber = atoi(piece[LineNumber]);   // This is the mutation line number from the spreadsheet, used to match it with its SAM file
//        const size_t SAMFileNameSize = 1000;
        lineNumber++;
        char *SAMFileName = piece[RNAFileName];

        samFile = fopen(SAMFileName, "r");
        if (NULL == samFile) {
            //
            //
            fprintf(stderr, "Unable to open SAM file %s\n", SAMFileName);
            noSAMFile++;
            goto doneWithMutation;
        }

        unsigned mutationStartPos = atoi(piece[StartPosition]);
        if (0 == mutationStartPos) {
            fprintf(stderr, "Got 0 mutation start pos from field '%s' on mutation line number %d, ignoting mutation\n", piece[StartPosition], mutationNumber);
            goto doneWithMutation;
        }

        if (!strcmp("M", piece[Chrom])) {
            //
            // Mitochondrial variations aren't useful with this method.  Most of them were stripped earlier, but I only
            // checked for "MT" as the chromosome name there, and some of them are just "M".  Drop them now.
            //
            nMT++;
            goto doneWithMutation;
        }

#if 0

        if (!strcmp("RNA", piece[VariantClass]) || !strcmp("Silent", piece[VariantClass])) {
            //
            // I'm not quite sure what an RNA variant is, but oddly, it doesn't seem to show up in RNA.
            // Silent variants aren't interesting, because they don't affect the protein and so probably
            // aren't real variants at all, and so you wouldn't expect to see loss of the other copy.
            //
            nWrongClass++;
            goto doneWithMutation;
        }

#endif // 0

        GenomeLocation mutationContigStartLocation;
        if (!ReferenceGenome->getLocationOfContig(piece[Chrom], &mutationContigStartLocation)) {
            //
            // try prefixing it with "chr"
            //
            char chromosomeName[1000];
            sprintf(chromosomeName, "chr%s", piece[Chrom]);
            if (!ReferenceGenome->getLocationOfContig(chromosomeName, &mutationContigStartLocation)) {
                //
                // Some of the alternate alleles have different names in the MAF files from our build.  Translate them.
                //
                const char *alternateChromName = translateAlternateChromName(piece[Chrom]);
                if (NULL == alternateChromName || !ReferenceGenome->getLocationOfContig(alternateChromName, &mutationContigStartLocation)) {
                    fprintf(stderr, "Couldn't find contig '%s' for mutation number %d\n", piece[Chrom], mutationNumber);
                    goto doneWithMutation;
                }
            }
        }

        GenomeLocation mutationLocation = mutationContigStartLocation + mutationStartPos - 1;   // -1 is because of 1-based reference

        const char *refAtMutation = ReferenceGenome->getSubstring(mutationLocation, 1);
        if (NULL == refAtMutation) {
            fprintf(stderr, "Couldn't get reference data for contig '%s', location %d, mutation number %d\n", piece[Chrom], mutationStartPos, mutationNumber);
            goto doneWithMutation;
        }

        const size_t matchAreaSize = 30;
        char referenceMatch[matchAreaSize];
        char mutantMatch[matchAreaSize];
        memset(mutantMatch, 0, matchAreaSize);

        unsigned matchAreaStartPos;
        matchAreaStartPos = (unsigned)__max(1, ((int)mutationStartPos) - 10);
        const char *refAtMatchStart;
        refAtMatchStart = ReferenceGenome->getSubstring(mutationLocation + matchAreaStartPos - mutationStartPos, 1);

        memcpy(referenceMatch, refAtMatchStart, matchAreaSize);

        bool mutationIsInsertion = false;
        bool mutationIsDeletion = false;
        bool mutationIsSNV = false;

        if (!strcmp(piece[VariantType], "DEL")) {
            mutationIsDeletion = true;
            if (strcmp(piece[Tumor_Seq_Allele_2], "-")) {
                fprintf(stderr, "Mutation is a deletion, but has a non-empty tumor allele '%s', line number %d\n", piece[Tumor_Seq_Allele_2], lineNumber);
                nValidationProblems++;
                goto doneWithMutation;
            }
            size_t deletionLength = strlen(piece[Reference_Allele]);
            memcpy(mutantMatch, referenceMatch, mutationStartPos - matchAreaStartPos);
            memcpy(mutantMatch + mutationStartPos - matchAreaStartPos, refAtMutation + deletionLength, matchAreaSize - (mutationStartPos - matchAreaStartPos));
        } else if (!strcmp(piece[VariantType], "INS")) {
            mutationIsInsertion = true;
            if (strcmp(piece[Reference_Allele], "-")) {
                fprintf(stderr, "Mutation is insertion but has a reference allele, mutation number %d\n", mutationNumber);
                nValidationProblems++;
                goto doneWithMutation;
            }
            size_t insertionLength = strlen(piece[Tumor_Seq_Allele_2]);
            size_t preMutantBaseCount = (size_t)(mutationStartPos - matchAreaStartPos+1);
            memcpy(mutantMatch, referenceMatch, preMutantBaseCount);
            if (insertionLength + preMutantBaseCount >= matchAreaSize) {
                //
                // The insertion fills (or overflows) the rest of the match area
                //
                memcpy(mutantMatch + preMutantBaseCount, piece[Tumor_Seq_Allele_2], matchAreaSize - preMutantBaseCount);
            } else {
                //
                // The insertion foes not fill the rest of the match area
                //
                memcpy(mutantMatch + preMutantBaseCount, piece[Tumor_Seq_Allele_2], insertionLength);
                memcpy(mutantMatch + preMutantBaseCount + insertionLength, refAtMutation + 1, matchAreaSize - preMutantBaseCount - insertionLength);
            }
        } else if (!strcmp(piece[VariantType], "SNP") || !strcmp(piece[VariantType], "DNP") || !strcmp(piece[VariantType], "TNP")) {
            // Single, double or triple nucleotide polymorphism.

            int polymorphismLength;
            if (!strcmp(piece[VariantType], "SNP")) {
                polymorphismLength = 1;
            } else if (!strcmp(piece[VariantType], "DNP")) {
                polymorphismLength = 2;
            } else {
                polymorphismLength = 3;
            }
            mutationIsSNV = true;
            if (strlen(piece[Reference_Allele]) != polymorphismLength && strlen(piece[Tumor_Seq_Allele_2]) != polymorphismLength) {
                fprintf(stderr, "S/D/TNP with ref and/or tumor wrong length: type %s ref %s, tumor %s, mutation number %d\n", piece[VariantType], piece[Reference_Allele], piece[Tumor_Seq_Allele_2], mutationNumber);
                nValidationProblems++;
                goto doneWithMutation;
            }
            if (memcmp(refAtMutation, piece[Reference_Allele], polymorphismLength)) {
                fprintf(stderr, "Reference doesn't match expected reference allele (type %s), %c != %c, mutation number %d, chromosome %s, location %d\n", 
                    piece[VariantType], *refAtMutation, piece[Reference_Allele][0], mutationNumber, piece[Chrom], mutationStartPos);
                nValidationProblems++;
                goto doneWithMutation;
            }
            if (!memcmp(piece[Tumor_Seq_Allele_2], refAtMutation, polymorphismLength)) {
                fprintf(stderr, "Mutation is %s, but matches reference, mutation number %d\n", piece[VariantType], mutationNumber);
                nValidationProblems++;
                goto doneWithMutation;
            }
            memcpy(mutantMatch, referenceMatch, matchAreaSize);
            mutantMatch[mutationStartPos - matchAreaStartPos] = piece[Tumor_Seq_Allele_2][0];
        } else {
            fprintf(stderr, "Unknown mutation type: '%s', mutation number %d\n", piece[VariantType], mutationNumber);
            nValidationProblems++;
            goto doneWithMutation;
        }



        char samLine[inputLineSize];
        int samFileLineNumber = 0;
        while (fgets(samLine, inputLineSize - 1, samFile)) {
            const int Flag = 1;
            const int RName = 2;
            const int Pos = 3;
            const int Cigar = 5;
            const int Seq = 9;
            bool matchesReference = true;
            bool matchesMutant = true;
            bool checkedMutationLocation = false;

            samFileLineNumber++;    // Increment first for 1-based

            char *samPieces[Seq + 2];
            if (!BreakLineIntoPiecesBasedOnTabs(samLine, samPieces, Seq + 2)) {
                fprintf(stderr, "Unparsable SAM line in file '%s', line number %d, line '%s'\n", SAMFileName, samFileLineNumber, samLine);
                nSAMProblems++;
                goto doneWithMutation;
            }

            unsigned samFlags = atoi(samPieces[Flag]);
            if (samFlags & SAM_UNMAPPED) {
                //
                // This read isn't mapped, it's probably the mate pair of a read that was mapped. Ignore it.
                //
                continue;
            }

            unsigned pos = atoi(samPieces[Pos]);
            if (pos > mutationStartPos) {
                //
                // This read is mapped after the end of the mutation.  Ignore it.
                //
                continue;
            }

            //
            // Walk through the read and CIGAR string and find the base that's mapped at the mutation position.
            //

            int offsetInRead = 0;
            unsigned currentReferencePos = atoi(samPieces[Pos]);
            unsigned currentReadPos = currentReferencePos;
            if (currentReferencePos == 0) {
                fprintf(stderr, "Got illegal POS ('%s') in mapped sam line, file '%s', line %d\n", samPieces[Pos], SAMFileName, samFileLineNumber);
                nSAMProblems++;
                goto doneWithMutation;
            }

            size_t seqLength = strlen(samPieces[Seq]);

            //
            // Walk through the read and CIGAR string
            //
            size_t offsetInCigarString = 0;
            const char *cigarString = samPieces[Cigar];

            while (cigarString[offsetInCigarString] != '\0') {
                _ASSERT(currentReadPos <= mutationStartPos);

                unsigned cigarCount = atoi(cigarString + offsetInCigarString);
                if (0 == cigarCount && cigarString[offsetInCigarString] != '0') { // Yes, there's actually a CIGAR string that's 27M0N48M
                    fprintf(stderr, "Unable to properly parse CIGAR string '%s', sam file '%s', line %d\n", cigarString, SAMFileName, samFileLineNumber);
                    nSAMProblems++;
                    goto doneWithMutation;
                }

                while (cigarString[offsetInCigarString] >= '0' && cigarString[offsetInCigarString] <= '9') {
                    offsetInCigarString++;
                }

                char op = cigarString[offsetInCigarString];
                offsetInCigarString++;

                switch (op) {

                case 'S':
                    offsetInRead += cigarCount;
                    break;

#if 0
                case 'M':
                    //
                    // See if this includes the mutant base.
                    //
                    if (currentReadPos + cigarCount > mutationStartPos) {
                        //
                        // The mutation is in this chunk of read.
                        //
                        if (samPieces[Seq][offsetInRead + mutationStartPos - currentReadPos] == *refAtMutation) {
                            //
                            // It matches the reference.
                            //
                            nMatchingReference++;
                        }
                        else if (mutationIsSNV && samPieces[Seq][offsetInRead + mutationStartPos - currentReadPos] == piece[Tumor_Seq_Allele_2][0]) {
                            nMatchingMutation++;
                        }
                        else {
                            nMatchingNeither++;
                        }

                    }

                    //
                    // Move along equal amounts in read, reference and read string.
                    //
                    currentReadPos += cigarCount;
                    offsetInRead += cigarCount;
                    currentReferencePos += cigarCount;
                    break;

                case 'I':
                    if (currentReadPos == mutationStartPos) {
                        if (mutationIsInsertion) {
                            //
                            // See if the bases in the read match the bases in the mutation.
                            //
                            if (seqLength - offsetInRead >= strlen(piece[Tumor_Seq_Allele_2]) && !memcmp(samPieces[Seq] + offsetInRead, piece[Tumor_Seq_Allele_2], strlen(piece[Tumor_Seq_Allele_2]))) {
                                nMatchingMutation++;
                            }
                            else {
                                nMatchingNeither++; // Insertions never match the reference.
                            }
                        }
                        else {
                            nMatchingNeither++;
                        }
                        currentReadPos++;   // Just to get us out of processing this read
                    }

                    offsetInRead += cigarCount;
                    break;


                case 'D':
                    if (currentReadPos + cigarCount > mutationStartPos) {
                        if (mutationIsDeletion) {
                            nMatchingMutation++;
                        }
                        else {
                            nMatchingNeither++;
                        }
                    }

                    currentReadPos += cigarCount;
                    offsetInRead += cigarCount;
                    break;


#else // 0
                case 'M':
                    if (currentReadPos + cigarCount > matchAreaStartPos) {
                        //
                        // This covers at least some of the match area
                        //
                        unsigned firstPos = __max(currentReadPos, matchAreaStartPos);
                        unsigned lastPos = __min(currentReadPos + cigarCount, matchAreaStartPos + matchAreaSize);

                        for (unsigned pos = firstPos; pos < lastPos; pos++) {
                            matchesReference &= (samPieces[Seq][offsetInRead + pos - currentReadPos] == referenceMatch[pos - matchAreaStartPos]);
                            matchesMutant    &= (samPieces[Seq][offsetInRead + pos - currentReadPos] == mutantMatch   [pos - matchAreaStartPos]);
                        }

                        if (currentReadPos <= mutationStartPos && currentReadPos + cigarCount > mutationStartPos) {
                            checkedMutationLocation = true;
                        }
                    }

                    //
                    // Move along equal amounts in read, reference and read string.
                    //
                    currentReadPos += cigarCount;
                    offsetInRead += cigarCount;
                    currentReferencePos += cigarCount;
                    break;

                case 'D':
                    
                    currentReferencePos += cigarCount;

                    break;


                case 'I':

                    if (currentReadPos > matchAreaStartPos) {
                        //
                        // This covers at least some of the match area
                        //
                        unsigned firstPos = __max(currentReadPos, matchAreaStartPos);
                        unsigned lastPos = __min(currentReadPos + cigarCount, matchAreaStartPos + matchAreaSize);

                        for (unsigned pos = firstPos; pos < lastPos; pos++) {
                            matchesReference &= (samPieces[Seq][offsetInRead + pos - currentReadPos] == referenceMatch[pos - matchAreaStartPos]);
                            matchesMutant &= (samPieces[Seq][offsetInRead + pos - currentReadPos] == mutantMatch[pos - matchAreaStartPos]);
                        }

                        if (currentReadPos == mutationStartPos) {
                            checkedMutationLocation = true;
                        }
                    }

                    offsetInRead += cigarCount;
                    break;

#endif // 0
                case 'N':
                    currentReadPos += cigarCount;
                    break;

                default:
                    fprintf(stderr, "Invalid cigar operation %c in cigar string '%s' on line %d of SAM file '%s'\n", op, cigarString, samFileLineNumber, SAMFileName);
                    nSAMProblems++;
                    goto doneWithMutation;
                } // switch


                if (currentReadPos > mutationStartPos) {
                    break;
                }
            } // While still processing CIGAR string

            if (checkedMutationLocation) {
                if (matchesReference) {
                    if (matchesMutant) {
                        nMatchingBoth++;
                    } else {
                        nMatchingReference++;
                    }
                } else if (matchesMutant) {
                    nMatchingMutation++;
                } else {
                    nMatchingNeither++;
                }
            }
        } // For each line in the SAM file

        //
        // Write out the input line, but with the counts appended.  First strip the newline out of the raw input buffer.
        //

        char *newline = strchr(rawInputBuffer, '\n');
        if (NULL != newline) {
            *newline = '\0';
        }

        //
        // Now figure out the tumor type by extracting the first component of the pathname.
        //

        char *pathSep = strchr(piece[RNAFileName], '\\');
        char cancerType[100];
        if (NULL == pathSep) {
            strcpy(cancerType, "Unknown");
        } else {
            *pathSep = '\0';
            strcpy(cancerType, piece[RNAFileName]);
        }


        printf("%s\t%s\t%d\t%d\t%d\t%d\n", cancerType, rawInputBuffer, nMatchingReference, nMatchingMutation, nMatchingNeither, nMatchingBoth);

    doneWithMutation:
        if (NULL != samFile) {
            fclose(samFile);
            samFile = NULL;
        }
        mutationNumber++;
    }

    fprintf(stderr, "%d mutations, %d mitochondrial, %d no SAM file, %d wrong class, %d validation problems, %d SAM file problems\n", mutationNumber, nMT, noSAMFile, nWrongClass, nValidationProblems, nSAMProblems);
}

