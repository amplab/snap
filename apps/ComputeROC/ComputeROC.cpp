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
#include "Bam.h"
#include "Genome.h"
#include "Compat.h"
#include "Read.h"
#include "RangeSplitter.h"
#include "BigAlloc.h"

void usage()
{
    fprintf(stderr,"usage: ComputeROC genomeDirectory inputFile {-b}\n");
    fprintf(stderr,"       -b means to accept reads that match either end of the range regardless of RC\n");
    fprintf(stderr,"       -c means to just count the number of reads that are aligned, not to worry about correctness\n");
    fprintf(stderr,"       -v means to correct for the error in generating the wgsim coordinates in the Venter data\n");
    fprintf(stderr,"       -e means to print out misaligned reads where the aligned location has a lower edit distance than the 'correct' one.\n");
    fprintf(stderr,"       -70 means to print out any misaligned reads with MAPQ 70.\n");
    fprintf(stderr,"You can specify only one of -b or -c\n");
  	exit(1);
}

ReadSupplierGenerator *readSupplierGenerator;
volatile _int64 nRunningThreads;
SingleWaiterObject allThreadsDone;
const char *inputFileName;
const Genome *genome;
bool matchBothWays = false;
bool justCount = false;
bool venter = false;
unsigned slackAmount = 151;
bool printBetterErrors = false;
bool printErrorsAtMAPQ70 = false;

static const int MaxMAPQ = 70;
const unsigned MaxEditDistance = 100;
struct ThreadContext {
    unsigned    whichThread;

    _int64 countOfReads[MaxMAPQ+1];
    _int64 countOfMisalignments[MaxMAPQ+1];
    _int64 countOfMisalignetsWithBetterEditDistance[MaxMAPQ+1];
    _int64 nUnaligned;
    _int64 totalReads;

    _int64 countOfReadsByEditDistance[MaxMAPQ+1][MaxEditDistance+1];
    _int64 countOfMisalignmentsByEditDistance[MaxMAPQ+1][MaxEditDistance+1];
    _int64 countOfMisalignetsWithBetterEditDistanceByEditDistance[MaxMAPQ+1][MaxEditDistance+1];
 

    ThreadContext() {
        nUnaligned = 0;
        totalReads = 0;
        for (int i = 0; i <= MaxMAPQ; i++) {
            countOfReads[i] = countOfMisalignments[i] = countOfMisalignetsWithBetterEditDistance[i] = 0;
            for (int j = 0; j <= MaxEditDistance; j++) {
                countOfReadsByEditDistance[i][j] = 0;
                countOfMisalignmentsByEditDistance[i][j] = 0;
                countOfMisalignetsWithBetterEditDistanceByEditDistance[i][j] = 0;
            }
        }
    }

};

bool inline isADigit(char x) {
    return x >= '0' && x <= '9';
}

void
WorkerThreadMain(void *param)
{
    ThreadContext *context = (ThreadContext *)param;

    ReadSupplier *readSupplier = readSupplierGenerator->generateNewReadSupplier();

    Read *read;
    LandauVishkinWithCigar lv;
    while (NULL != (read = readSupplier->getNextRead())) {
        unsigned mapQ = read->getOriginalMAPQ();
        unsigned genomeLocation = read->getOriginalAlignedLocation();
        unsigned flag = read->getOriginalSAMFlags();

        if (flag & SAM_UNMAPPED) {
            genomeLocation = 0xffffffff;
        }

        if (mapQ < 0 || mapQ > MaxMAPQ) {
            fprintf(stderr,"Invalid MAPQ: %d\n",mapQ);
            exit(1);
        }

        context->totalReads++;

        if (0xffffffff == genomeLocation) {
            context->nUnaligned++;
        } else if (justCount) {
            context->countOfReads[mapQ]++;
	    } else if (!justCount) {
            if (flag & SAM_REVERSE_COMPLEMENT) {
                read->becomeRC();
            }
                            
            const Genome::Contig *contig = genome->getContigAtLocation(genomeLocation);
            if (NULL == contig) {
                fprintf(stderr,"couldn't find genome contig for offset %u\n",genomeLocation);
                exit(1);
            }
            unsigned offsetA, offsetB;
            bool matched;

            const unsigned cigarBufLen = 1000;
            char cigarForAligned[cigarBufLen];
            const char *alignedGenomeData = genome->getSubstring(genomeLocation, 1); 
            int editDistance = lv.computeEditDistance(alignedGenomeData, read->getDataLength() + 20, read->getData(), read->getDataLength(), 30, cigarForAligned, cigarBufLen, false);

            if (editDistance == -1 || editDistance > MaxEditDistance) {
                editDistance = MaxEditDistance;
            }

            //
            // Parse the read ID.  The format is ChrName_OffsetA_OffsetB_?:<more stuff>.  This would be simple to parse, except that
            // ChrName can include "_".  So, we parse it by looking for the first : and then working backward.
            //
            char idBuffer[10000];   // Hopefully big enough.  I'm not worried about malicious input data here.

            memcpy(idBuffer,read->getId(),read->getIdLength());
            idBuffer[read->getIdLength()] = 0;
                    
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

                        offsetA = -1;
                        offsetB = -1;

                        if (*(beginningOfFirstNumber - 1) == '_' && 1 == sscanf(beginningOfFirstNumber,"%u",&offsetA) &&
                            ('_' == *beginningOfSecondNumber || 1 == sscanf(beginningOfSecondNumber,"%u", &offsetB))) {
                                stage = 3;

                            chrNameLen = (beginningOfFirstNumber - 1) - idBuffer;
                            char correctChromosomeName[1000];
                            memcpy(correctChromosomeName, idBuffer, chrNameLen);
                            correctChromosomeName[chrNameLen] = '\0';

                            if (venter && offsetB >= read->getDataLength()) {
                                offsetB -= read->getDataLength();
                            }

                            if (!genome->getOffsetOfContig(correctChromosomeName, &offsetOfCorrectChromosome)) {
                                fprintf(stderr, "Couldn't parse chromosome name '%s' from read id\n", correctChromosomeName);
                            } else {
                                badParse = false;
                            }
                        }
                    }
                }

                if (badParse) {
                    fprintf(stderr,"Unable to parse read ID '%s', perhaps this isn't simulated data.  contiglen = %d, contigName = '%s', contig offset = %u, genome offset = %u\n", idBuffer, strlen(contig->name), contig->name, contig->beginningOffset, genomeLocation);
                    exit(1);
                }

 
                bool match0 = false;
                bool match1 = false;
                if (-1 == offsetA || -1 == offsetB) {
                    matched = false;
                }  else if(strncmp(contig->name, idBuffer, __min(read->getIdLength(), chrNameLen))) {
                    matched = false;
                } else {
                    if (isWithin(offsetA, genomeLocation - contig->beginningOffset, slackAmount)) {
                        matched = true;
                        match0 = true;
                    } else if (isWithin(offsetB, genomeLocation - contig->beginningOffset, slackAmount)) {
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
                }

                context->countOfReads[mapQ]++;
                context->countOfReadsByEditDistance[mapQ][editDistance]++;

                if (!matched) {
                    context->countOfMisalignments[mapQ]++;
                    context->countOfMisalignmentsByEditDistance[mapQ][editDistance]++;

                    if ((70 == mapQ && printErrorsAtMAPQ70) || printBetterErrors) {

                        //
                        // We don't know which offset is correct, because neither one matched.  Just take the one with the lower edit distance.
                        //
                        unsigned correctLocationA = offsetOfCorrectChromosome + offsetA;
                        unsigned correctLocationB = offsetOfCorrectChromosome + offsetB;

                        unsigned correctLocation = 0;
                        const char *correctData = NULL;

                        const char *dataA = genome->getSubstring(correctLocationA, 1);
                        const char *dataB = genome->getSubstring(correctLocationB, 1);
                        int distanceA, distanceB;
                        char cigarA[cigarBufLen];
                        char cigarB[cigarBufLen];

                        cigarA[0] = '*'; cigarA[1] = '\0';
                        cigarB[0] = '*'; cigarB[1] = '\0';

                        if (dataA == NULL) {
                            distanceA = -1;
                        } else {
                            distanceA = lv.computeEditDistance(dataA, read->getDataLength() + 20, read->getData(), read->getDataLength(), 30, cigarA, cigarBufLen, false);
                        }

                        if (dataB == NULL) {
                            distanceB = -1;
                        } else {
                            distanceB = lv.computeEditDistance(dataB, read->getDataLength() + 20, read->getData(), read->getDataLength(), 30, cigarB, cigarBufLen, false);
                        }

                        const char *correctGenomeData;
                        char *cigarForCorrect;

                        if (distanceA != -1 && distanceA <= distanceB || distanceB == -1) {
                            correctGenomeData = dataA;
                            correctLocation = correctLocationA;
                            cigarForCorrect = cigarA;
                        } else {
                            correctGenomeData = dataB;
                            correctLocation = correctLocationB;
                            cigarForCorrect = cigarB;
                        }

                        bool betterEditDistance = ((distanceA > editDistance && distanceB > editDistance) || (-1 == distanceA && -1 == distanceB));
                        if (betterEditDistance) {
                            context->countOfMisalignetsWithBetterEditDistanceByEditDistance[mapQ][editDistance]++;
                            context->countOfMisalignetsWithBetterEditDistance[mapQ]++;
                        }

                        // if (!printBetterErrors || (printBetterErrors && betterEditDistance)) {
                           
                        //     printf("%s\t%d\t%s\t%u\t%d\t%s\t*\t*\t100\t%.*s\t%.*s\tAlignedGenomeLocation:%u\tCorrectGenomeLocation: %u\tCigarForCorrect: %s\tCorrectData: %.*s\tAlignedData: %.*s\n", 
                        //         idBuffer, flag, contig->name, genomeLocation - contig->beginningOffset, mapQ, cigarForAligned, read.getDataLength(), read.getData(), 
                        //         read.getDataLength(), read.getQuality(),  genomeLocation, correctLocation, cigarForCorrect, read.getDataLength(),
                        //         correctGenomeData, read.getDataLength(), alignedGenomeData);
                        //}
                    }
                }
            }
        } // if it was mapped
    } // for each read from the sam reader

     if (0 == InterlockedAdd64AndReturnNewValue(&nRunningThreads, -1)) {
        SignalSingleWaiterObject(&allThreadsDone);
    }
}


int main(int argc, char * argv[])
{
    BigAllocUseHugePages = false;

    if (argc < 3) usage();

    for (int i = 3; i < argc; i++) {
        if (!strcmp(argv[i], "-b")) {
            matchBothWays = true;
        } else if (!strcmp(argv[i], "-c")) {
            justCount = true;
        } else if (!strcmp(argv[i], "-v")) {
            venter = true;        
        } else if (!strcmp(argv[i], "-e")) {
            printBetterErrors = true;        
        } else if (!strcmp(argv[i], "-70")) {
            printErrorsAtMAPQ70 = true;
        } else {
            usage();
        }
    }

    static const char *genomeSuffix = "Genome";
	size_t filenameLen = strlen(argv[1]) + 1 + strlen(genomeSuffix) + 1;
	char *fileName = new char[strlen(argv[1]) + 1 + strlen(genomeSuffix) + 1];
	snprintf(fileName,filenameLen,"%s%c%s",argv[1],PATH_SEP,genomeSuffix);
	genome = Genome::loadFromFile(fileName, 0);
	if (NULL == genome) {
		fprintf(stderr,"Unable to load genome from file '%s'\n",fileName);
		return -1;
	}
	delete [] fileName;
	fileName = NULL;

    inputFileName = argv[2];

    unsigned nThreads;
#ifdef _DEBUG
    nThreads = 1;
#else   // _DEBUG
    nThreads = GetNumberOfProcessors();
#endif // _DEBUG


    DataSupplier::ThreadCount = nThreads;
    nRunningThreads = nThreads;

    ReaderContext readerContext;
    readerContext.clipping = NoClipping;
    readerContext.defaultReadGroup = "";
    readerContext.genome = genome;
    readerContext.ignoreSecondaryAlignments = true;
    readerContext.ignoreSupplementaryAlignments = true;
	readerContext.header = NULL;
	readerContext.headerLength = 0;
	readerContext.headerBytes = 0;

    if (NULL != strrchr(inputFileName, '.') && !_stricmp(strrchr(inputFileName, '.'), ".bam")) {
        readSupplierGenerator = BAMReader::createReadSupplierGenerator(inputFileName, nThreads, readerContext);
    } else {
        readSupplierGenerator = SAMReader::createReadSupplierGenerator(inputFileName, nThreads, readerContext);
    }

    CreateSingleWaiterObject(&allThreadsDone);
    ThreadContext *contexts = new ThreadContext[nThreads];

    for (unsigned i = 0; i < nThreads; i++) {
        contexts[i].whichThread = i;

        StartNewThread(WorkerThreadMain, &contexts[i]);
    }

    WaitForSingleWaiterObject(&allThreadsDone);

    _int64 nUnaligned = 0;
    _int64 totalReads = 0;
    for (unsigned i = 0; i < nThreads; i++) {
        nUnaligned += contexts[i].nUnaligned;
        totalReads += contexts[i].totalReads;
    }
    printf("%lld reads, %lld unaligned (%0.2f%%)\n", (long long)totalReads, (long long)nUnaligned, 100. * (double)nUnaligned / (double)totalReads);

    printf("MAPQ\tnReads\tnMisaligned");
    if (printBetterErrors) {
        printf("\tBetterMisaligned");
    }
    printf("\n");
    for (int i = 0; i <= MaxMAPQ; i++) {
        _int64 nReads = 0;
        _int64 nMisaligned = 0;
        _int64 betterMisaligned = 0;
        for (unsigned j = 0; j < nThreads; j++) {
            nReads += contexts[j].countOfReads[i];
            nMisaligned += contexts[j].countOfMisalignments[i];
            betterMisaligned += contexts[j].countOfMisalignetsWithBetterEditDistance[i];
        }
        printf("%d\t%lld\t%lld", i, (long long)nReads,  (long long)nMisaligned);
        if (printBetterErrors) {
            printf("\t%lld",  (long long)betterMisaligned);
        }
        printf("\n");
    }

    int maxEditDistanceSeen = 0;
    for (unsigned i = 0; i < nThreads; i++) {
    }

	return 0;
}

