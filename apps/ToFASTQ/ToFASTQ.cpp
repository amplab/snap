/*++

Module Name:

    ToFASTQ.cpp

Abstract:

   Take a set of reads in SAM or BAM format and convert them to FASTQ

Authors:

    Bill Bolosky, January, 2014

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
#include "FASTQ.h"

void usage()
{
    fprintf(stderr,"usage: ToFASTQ genomeIndex inputFile outputFile {outputFile2}\n");
    fprintf(stderr,"       Specifying two output files means that the input is paired.  If you specify only one output file, then\n");
    fprintf(stderr,"       ToFASTQ will generate a single-ended FASTQ even for a paired input.\n");
    fprintf(stderr,"       The genomeIndex must contain the same set of contigs used to align the input file.\n");
    fprintf(stderr,"       To produce interleaved paired-end FASTQ, specify outputFile2 as '-i'.\n");
  	soft_exit(1);
}

ReadSupplierGenerator *readSupplierGenerator = NULL;

volatile _int64 nRunningThreads;
SingleWaiterObject allThreadsDone;
const char *inputFileName;
const Genome *genome;
FASTQWriter *fastqWriter[2] = {NULL, NULL};

struct ThreadContext {
    unsigned    whichThread;

    _int64 totalReads;

    ThreadContext() {
        totalReads = 0;
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
     while (NULL != (read = readSupplier->getNextRead())) {
        context->totalReads++;
        fastqWriter[0]->writeRead(read);
     } // for each read from the reader

     if (0 == InterlockedAdd64AndReturnNewValue(&nRunningThreads, -1)) {
        SignalSingleWaiterObject(&allThreadsDone);
    }
}

_int64
ProcessPairedInput(PairedReadSupplierGenerator *pairedReadSupplierGenerator)
{
    //
    // This runs single threaded, because it's the only easy way to assure that the two output files match.
    //
    PairedReadSupplier *readSupplier = pairedReadSupplierGenerator->generateNewPairedReadSupplier();
    
    Read *read[NUM_READS_PER_PAIR];

    const size_t idBufferSize = 10000;
    char idBuffer[idBufferSize];

    _int64 totalReads = 0;

    while (readSupplier->getNextReadPair(&read[0], &read[1])) {
        for (int i = 0; i < NUM_READS_PER_PAIR; i++) {
            //
            // Create a new read with an ID that includes /1 or /2.  Don't bother to fill in the stuff that the
            // FASTQ writer doesn't care about anyway.
            //
            snprintf(idBuffer, idBufferSize-1,"%.*s/%d", read[i]->getIdLength(), read[i]->getId(), i+1);
            Read local;
            local.init(idBuffer, (unsigned)strlen(idBuffer),read[i]->getUnclippedData(), read[i]->getUnclippedQuality(), read[i]->getUnclippedLength(),
                0,0,0,0,0,0,0,NULL,0,0);
            fastqWriter[i]->writeRead(&local);
        }
        totalReads += 2;
    }

    return totalReads;
}


int main(int argc, char * argv[])
{
    BigAllocUseHugePages = false;

    if (4 != argc && 5 != argc) usage();

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

    fastqWriter[0] = FASTQWriter::Factory(argv[3]);
    if (NULL == fastqWriter[0]) {
        fprintf(stderr,"Unable to open FASTQ writer for output file '%s'\n", argv[3]);
        soft_exit(1);
    }

    _int64 totalReads = 0;
    unsigned nThreads;

#ifdef _DEBUG
    nThreads = 1;
#else   // _DEBUG
    if (5 == argc) {
        nThreads = 1;
    } else {
        nThreads = GetNumberOfProcessors();
    }
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

    if (5 == argc) {
        if (!strcmp(argv[4], "-i")) {
            fastqWriter[1] = fastqWriter[0];    // Interleaved
        } else {
            fastqWriter[1] = FASTQWriter::Factory(argv[4]);
            if (NULL == fastqWriter[1]) {
                fprintf(stderr,"Unable to open FASTQ writer for output file '%s'\n", argv[4]);
                soft_exit(1);
            }
        }

        PairedReadSupplierGenerator *pairedReadSupplierGenerator;
        if (NULL != strrchr(inputFileName, '.') && !_stricmp(strrchr(inputFileName, '.'), ".bam")) {
            pairedReadSupplierGenerator = BAMReader::createPairedReadSupplierGenerator(inputFileName, nThreads, true, readerContext);
        } else {
            pairedReadSupplierGenerator = SAMReader::createPairedReadSupplierGenerator(inputFileName, nThreads, true, readerContext);
        }


        totalReads = ProcessPairedInput(pairedReadSupplierGenerator);
        if (fastqWriter[1] != fastqWriter[0]) {
            delete fastqWriter[1];
        } else {
            fastqWriter[1] = NULL;
        }
    } else {
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

        for (unsigned i = 0; i < nThreads; i++) {
            totalReads += contexts[i].totalReads;
        }
    }
    printf("%lld reads\n", totalReads);
    delete fastqWriter[0];

	return 0;
}

