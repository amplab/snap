/*++

Module Name:

    ComputeROC.cpp

Abstract:

   Take a SAM file with simulated reads and compute a ROC curve from it.

Authors:

    Bill Bolosky, December, 2012

Environment:
`
    User mode service.

Revision History:

    
--*/

#include "stdafx.h"
#include "SAM.h"
#include "Genome.h"
#include "Compat.h"
#include "Read.h"
#include "RangeSplitter.h"
#include "BigAlloc.h"

void usage()
{
    fprintf(stderr,"usage: ComputeROC genomeDirectory inputFile {-b}\n");
    fprintf(stderr,"       -b means to accept reads that match either end of the range regardless of RC\n");
  	exit(1);
}

RangeSplitter *rangeSplitter;
volatile _int64 nRunningThreads;
SingleWaiterObject allThreadsDone;
const char *inputFileName;
const Genome *genome;
bool matchBothWays = false;

static const int MaxMAPQ = 70;
const unsigned MaxEditDistance = 100;
struct ThreadContext {
    unsigned    whichThread;

    _int64 countOfReads[MaxMAPQ+1];
    _int64 countOfMisalignments[MaxMAPQ+1];
    _int64 nUnaligned;

    _int64 countOfReadsByEditDistance[MaxMAPQ+1][MaxEditDistance+1];
    _int64 countOfMisalignmentsByEditDistance[MaxMAPQ+1][MaxEditDistance+1];

    ThreadContext() {
        nUnaligned = 0;
        for (int i = 0; i <= MaxMAPQ; i++) {
            countOfReads[i] = countOfMisalignments[i] = 0;
            for (int j = 0; j <= MaxEditDistance; j++) {
                countOfReadsByEditDistance[i][j] = 0;
                countOfMisalignmentsByEditDistance[i][j] = 0;
            }
        }
    }

};

inline bool isWithin(unsigned a, unsigned b, unsigned distance)
{
    return a <= b && a+distance >= b || a >= b && a <= b + distance;
}

bool inline isADigit(char x) {
    return x >= '0' && x <= '9';
}

void
WorkerThreadMain(void *param)
{
    ThreadContext *context = (ThreadContext *)param;

    _int64 rangeStart, rangeLength;

    SAMReader *samReader = NULL;
    while (rangeSplitter->getNextRange(&rangeStart, &rangeLength)) {
        if (NULL == samReader) {
            samReader = SAMReader::create(DataSupplier::Default[false], inputFileName, genome, rangeStart, rangeLength);
        } else {
            ((ReadReader *)samReader)->reinit(rangeStart, rangeLength);
        }

        AlignmentResult alignmentResult;
        unsigned genomeLocation;
        Direction isRC;
        unsigned mapQ;
        unsigned flag;
        const char *cigar;
        unsigned nextFileToWrite = 0;
        Read read;
        LandauVishkinWithCigar lv;
        while (samReader->getNextRead(&read, &alignmentResult, &genomeLocation, &isRC, &mapQ, &flag, &cigar)) {

            if (mapQ < 0 || mapQ > MaxMAPQ) {
                fprintf(stderr,"Invalid MAPQ: %d\n",mapQ);
                exit(1);
            }

            if (0xffffffff == genomeLocation) {
                context->nUnaligned++;
            } else {
                const Genome::Piece *piece = genome->getPieceAtLocation(genomeLocation);
                if (NULL == piece) {
                    fprintf(stderr,"couldn't find genome piece for offset %u\n",genomeLocation);
                    exit(1);
                }
                unsigned offsetA, offsetB;
                bool matched;

                const unsigned cigarBufLen = 1000;
                char cigarForAligned[cigarBufLen];
                const char *alignedGenomeData = genome->getSubstring(genomeLocation, 1); 
                int editDistance = lv.computeEditDistance(alignedGenomeData, read.getDataLength() + 20, read.getData(), read.getDataLength(), 30, cigarForAligned, cigarBufLen, false);

                if (editDistance == -1 || editDistance > MaxEditDistance) {
                    editDistance = MaxEditDistance;
                }

                //
                // Parse the read ID.  The format is ChrName_OffsetA_OffsetB_?:<more stuff>.  This would be simple to parse, except that
                // ChrName can include "_".  So, we parse it by looking for the first : and then working backward.
                //
                char idBuffer[10000];   // Hopefully big enough.  I'm not worried about malicious input data here.

                memcpy(idBuffer,read.getId(),read.getIdLength());
                idBuffer[read.getIdLength()] = 0;
                    
                const char *firstColon = strchr(idBuffer,':');
                bool badParse = true;
                size_t chrNameLen;
                const char *beginningOfSecondNumber;
                const char *beginningOfFirstNumber; int stage = 0;
                unsigned offsetOfCorrectChromosome;
 
                if (NULL != firstColon && firstColon - 3 > idBuffer && (*(firstColon-1) == '?' || isADigit(*(firstColon - 1)))) {
                    //
                    // We've parsed backwards to see that we have at least #: or ?: where '#' is a digit and ? is literal.  If it's
                    // a digit, then scan backwards through that number.
                    //
                    const char *underscoreBeforeFirstColon = firstColon - 2;
                    while (underscoreBeforeFirstColon > idBuffer && isADigit(*underscoreBeforeFirstColon)) {
                        underscoreBeforeFirstColon--;
                    }

                    if (*underscoreBeforeFirstColon == '_' && (isADigit(*(underscoreBeforeFirstColon - 1)) || *(underscoreBeforeFirstColon - 1) == '_')) {
                        stage = 1;
                        if (isADigit(*(underscoreBeforeFirstColon - 1))) {
                            beginningOfSecondNumber = firstColon - 3;
                            while (beginningOfSecondNumber > idBuffer && isADigit(*beginningOfSecondNumber)) {
                                beginningOfSecondNumber--;
                            }
                            beginningOfSecondNumber++; // That loop actually moved us back one char before the beginning;
                        } else {
                            //
                            // There's only one number,  we have two consecutive underscores.
                            //
                            beginningOfSecondNumber = underscoreBeforeFirstColon;
                        }
                        if (beginningOfSecondNumber - 2 > idBuffer && *(beginningOfSecondNumber - 1) == '_' && isADigit(*(beginningOfSecondNumber - 2))) {
                            stage = 2;
                            beginningOfFirstNumber = beginningOfSecondNumber - 2;
                            while (beginningOfFirstNumber > idBuffer && isADigit(*beginningOfFirstNumber)) {
                                beginningOfFirstNumber--;
                            }
                            beginningOfFirstNumber++; // Again, we went one too far.

                            offsetB = -1;

                            if (*(beginningOfFirstNumber - 1) == '_' && 1 == sscanf(beginningOfFirstNumber,"%u",&offsetA) &&
                                ('_' == *beginningOfSecondNumber || 1 == sscanf(beginningOfSecondNumber,"%u", &offsetB))) {
                                    stage = 3;

                                chrNameLen = (beginningOfFirstNumber - 1) - idBuffer;
                                char correctChromosomeName[1000];
                                memcpy(correctChromosomeName, idBuffer, chrNameLen);
                                correctChromosomeName[chrNameLen] = '\0';

                                if (!genome->getOffsetOfPiece(correctChromosomeName, &offsetOfCorrectChromosome)) {
                                    fprintf(stderr, "Couldn't parse chromosome name '%s' from read id\n", correctChromosomeName);
                                } else {
                                    badParse = false;
                                }
                            }
                        }
                    }

 

                    if (badParse) {
                        fprintf(stderr,"Unable to parse read ID '%s', perhaps this isn't simulated data.  piecelen = %d, pieceName = '%s', piece offset = %u, genome offset = %u\n", idBuffer, strlen(piece->name), piece->name, piece->beginningOffset, genomeLocation);
                        exit(1);
                    }

                    offsetB = -1;
                    bool match0 = false;
                    bool match1 = false;
                    if (strncmp(piece->name, idBuffer, __min(read.getIdLength(), chrNameLen))) {
                        matched = false;
                    } else if (!(('_' == *beginningOfSecondNumber && 1 == sscanf(idBuffer + strlen(piece->name), "_%u", &offsetA) || 
                                (2 == sscanf(idBuffer + strlen(piece->name),"_%u_%u",&offsetA,&offsetB))))) {
                        matched = false;
                    } else {
#if 1
                        if (isWithin(offsetA, genomeLocation - piece->beginningOffset, 50)) {
                            matched = true;
                            match0 = true;
                        } else if (isWithin(offsetB, genomeLocation - piece->beginningOffset, 50)) {
                            matched = true;
                            match1 = true;
                        } else {
                            matched = false;
                            if (flag & SAM_FIRST_SEGMENT) {
                                match0 = true;
                            } else {
                                match1 = true;
                            }
                        }
#else   // 1

                        if (matchBothWays || flag & SAM_NEXT_UNMAPPED) {
                            match0 = true;
                            match1 = true;
                        } else if (flag & SAM_FIRST_SEGMENT) {
                            match0 = true;
                        } else if (flag & SAM_LAST_SEGMENT) {
                            match1 = true;
                        } else {
                            //
                            // Neither first nor last segement.  
                            match0 = true;
                            match1 = true;
                        }

                       matched = match0 && isWithin(offsetA, genomeLocation - piece->beginningOffset, 50) || match1 && isWithin(offsetB, genomeLocation - piece->beginningOffset, 50);

#endif // 1


                    }

                    context->countOfReads[mapQ]++;
                    context->countOfReadsByEditDistance[mapQ][editDistance]++;

                    if (!matched) {
                        context->countOfMisalignments[mapQ]++;
                        context->countOfMisalignmentsByEditDistance[mapQ][editDistance]++;

                        if (70 == mapQ || 69 == mapQ) {
                           if (flag & SAM_REVERSE_COMPLEMENT) {
                                read.becomeRC();
                            }
                            
                            unsigned correctLocation = (match0 ? offsetA : offsetB) + offsetOfCorrectChromosome;
 
                            char cigarForCorrect[cigarBufLen];


                            const char *correctGenomeData = genome->getSubstring(correctLocation, 1);
                           
                            bool foo = isWithin(offsetA, genomeLocation - piece->beginningOffset, 50) || match1 && isWithin(offsetB, genomeLocation - piece->beginningOffset, 50);

                            lv.computeEditDistance(correctGenomeData, read.getDataLength() + 20, read.getData(), read.getDataLength(), 30, cigarForCorrect, cigarBufLen, false);

 
                            printf("%s\t%d\t%s\t%u\t%d\t%s\t*\t*\t100\t%.*s\t%.*s\tAlignedGenomeLocation:%u\tCorrectGenomeLocation: %u\tCigarForCorrect: %s\tCorrectData: %.*s\tAlignedData: %.*s\n", 
                                idBuffer, flag, piece->name, genomeLocation - piece->beginningOffset, mapQ, cigarForAligned, read.getDataLength(), read.getData(), 
                                read.getDataLength(), read.getQuality(),  genomeLocation, correctLocation, cigarForCorrect, read.getDataLength(),
                                correctGenomeData, read.getDataLength(), alignedGenomeData);
                        }
                    }
                }
            } // if it was mapped
        } // for each read from the sam reader
    }

     if (0 == InterlockedAdd64AndReturnNewValue(&nRunningThreads, -1)) {
        SignalSingleWaiterObject(&allThreadsDone);
    }
}


int main(int argc, char * argv[])
{
    BigAllocUseHugePages = false;

    if (3 != argc && 4 != argc) {
		usage();
	}

    if (4 == argc) {
        if (!strcmp(argv[3],"-b")) {
            matchBothWays = true;
        } else {
            usage();
        }
    }

    static const char *genomeSuffix = "Genome";
	size_t filenameLen = strlen(argv[1]) + 1 + strlen(genomeSuffix) + 1;
	char *fileName = new char[strlen(argv[1]) + 1 + strlen(genomeSuffix) + 1];
	snprintf(fileName,filenameLen,"%s%c%s",argv[1],PATH_SEP,genomeSuffix);
	genome = Genome::loadFromFile(fileName); 
	if (NULL == genome) {
		fprintf(stderr,"Unable to load genome from file '%s'\n",fileName);
		return -1;
	}
	delete [] fileName;
	fileName = NULL;

    inputFileName = argv[2];

    nRunningThreads = GetNumberOfProcessors();

    _int64 fileSize = QueryFileSize(argv[2]);
    rangeSplitter = new RangeSplitter(fileSize, GetNumberOfProcessors());
    
    CreateSingleWaiterObject(&allThreadsDone);
    ThreadContext *contexts = new ThreadContext[GetNumberOfProcessors()];

    for (unsigned i = 0; i < GetNumberOfProcessors(); i++) {
        contexts[i].whichThread = i;

        StartNewThread(WorkerThreadMain, &contexts[i]);
    }

    WaitForSingleWaiterObject(&allThreadsDone);

    _int64 nUnaligned = 0;
    for (unsigned i = 0; i < GetNumberOfProcessors(); i++) {
        nUnaligned += contexts[i].nUnaligned;
    }
    printf("%lld total unaligned\nMAPQ\tnReads\tnMisaligned\n",nUnaligned);

    for (int i = 0; i <= MaxMAPQ; i++) {
        _int64 nReads = 0;
        _int64 nMisaligned = 0;
        for (unsigned j = 0; j < GetNumberOfProcessors(); j++) {
            nReads += contexts[j].countOfReads[i];
            nMisaligned += contexts[j].countOfMisalignments[i];
        }
        printf("%d\t%lld\t%lld\n", i, nReads, nMisaligned);
    }

    int maxEditDistanceSeen = 0;
    for (int i = 0; i < GetNumberOfProcessors(); i++) {
    }

	return 0;
}

