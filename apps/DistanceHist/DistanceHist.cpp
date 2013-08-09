/*++

Module Name:

    DistanceHist.cpp

Abstract:

    Compute a histogram of the edit distances between simulated reads and their correct 
    alignments.

Authors:

    Bill Bolosky, May 2013

Environment:
`
    User mode service.

Revision History:


--*/


#include "stdafx.h"
#include "Compat.h"
#include "Genome.h"
#include "exit.h"
#include "SAM.h"
#include "FASTQ.h"
#include "WGsim.h"
#include "LandauVishkin.h"
#include "Tables.h"

const Genome *genome = NULL;

struct DistHistogram {
    static const unsigned MaxDistance = 100;
    unsigned counts[MaxDistance+2];

    DistHistogram() {
        for (unsigned i = 0 ; i < MaxDistance+2; i++) {
            counts[i] = 0;
        }
    }

    void addIn(const DistHistogram &peer) {
        for (unsigned i = 0; i < MaxDistance+2; i++) {
            counts[i] += peer.counts[i];
        }
    }
};

ReadSupplierGenerator *readSupplierGenerator = NULL;
volatile int nRunningThreads;
SingleWaiterObject threadsDone;

void
usage()
{
    fprintf(stderr,"usage: DistanceHist index inputFile\n");
    soft_exit(1);
}

void workerThreadMain(void *context)
{
    DistHistogram histogram;    // Don't use the context one until the end to avoid false sharing

    LandauVishkinWithCigar lv;

    ReadSupplier *readSupplier = readSupplierGenerator->generateNewReadSupplier();

    const unsigned maxReadLen = MAX_READ_LENGTH;
    char *rcBuffer = new char[maxReadLen];

    Read *read;
    while (NULL != (read = readSupplier->getNextRead())) {
        unsigned readLen = read->getDataLength();
        unsigned highOffset, lowOffset;
        const char *readData = read->getData();
        const char *quality = read->getQuality();

        bool lowQual = false;
        for (unsigned i = 0 ; i < readLen; i++) {
            if (quality[i] < '?') {
                lowQual = true;
                break;
            }
        }
        if (lowQual) {
            continue;
        }

        char cigar[4][100];
        const char *genomeData[4];
        int edit[4];

        //
        // We don't care if it's misaligned (or in fact if it's aligned at all).  Just get the
        // offsets from the wgsim name.
        //
        wgsimReadMisaligned(read, 0, genome, 0, &lowOffset, &highOffset);

        unsigned bestAt = 0;
        cigar[0][0] = '\0';
        int bestDistance = read->getDataLength();
        genomeData[0] = genome->getSubstring(lowOffset, readLen + 20);
        int dist = edit[0] = lv.computeEditDistance(genome->getSubstring(lowOffset, readLen + 20), readLen + 20, readData, readLen, MAX_K - 1, cigar[0], 100, false);
        if (dist >= 0) {
            bestDistance = dist;
        }

        genomeData[1] = genome->getSubstring(highOffset, readLen + 20);
        edit[1] = dist = lv.computeEditDistance(genome->getSubstring(highOffset, readLen + 20), readLen + 20, readData, readLen, MAX_K - 1, cigar[1], 100, false);

        if (dist >= 0 && dist < bestDistance) {
            bestDistance = dist;
            bestAt = 1;
        }

        for (unsigned i = 0; i < readLen; i++) {
            rcBuffer[readLen - i - 1] = COMPLEMENT[readData[i]];
        }

        genomeData[2] = genome->getSubstring(lowOffset, readLen + 20);
        edit[2] = dist = lv.computeEditDistance(genome->getSubstring(lowOffset, readLen + 20), readLen + 20, rcBuffer, readLen, MAX_K - 1, cigar[2], 100, false);
        if (dist >= 0 && dist < bestDistance) {
            bestDistance =dist;
            bestAt = 2;
        }

        genomeData[3] = genome->getSubstring(highOffset, readLen + 20);
        edit[3] = dist = lv.computeEditDistance(genome->getSubstring(highOffset, readLen + 20), readLen + 20, rcBuffer, readLen, MAX_K - 1, cigar[3], 100, false);

        if (dist >= 0 && dist < bestDistance) {
            bestDistance = dist;
            bestAt = 3;
        }

        bool containsIndels = false;
        for (size_t i = 0; i < strlen(cigar[bestAt]); i++) {
            if (cigar[bestAt][i] == 'I' || cigar[bestAt][i] == 'D') {
                containsIndels = true;
                break;
            }
        }

        if (containsIndels) {
            continue;
        }

        if (bestDistance <0 || bestDistance > DistHistogram::MaxDistance) {
            histogram.counts[DistHistogram::MaxDistance]++;
        } else {
            histogram.counts[bestDistance]++;
        }
    }

    ((DistHistogram *)context)->addIn(histogram);
    if (0 == InterlockedDecrementAndReturnNewValue(&nRunningThreads)) {
        SignalSingleWaiterObject(&threadsDone);
    }

    delete rcBuffer;    // Cause I'm just that kinda guy.
}

void main(int argc, char * argv[])
{
	if (3 != argc) usage();

    BigAllocUseHugePages = false;

    const char *genomeFileName = "Genome";
    char *pathname = new char[strlen(argv[1]) + 1 /* for directory separator */ + strlen(genomeFileName) + 1 /* for null */];
    sprintf(pathname, "%s%c%s", argv[1], PATH_SEP, genomeFileName);

    _int64 start = timeInMillis();
    printf("Loading genome...");
    genome = Genome::loadFromFile(pathname, 0);
    if (NULL == genome) {
        fprintf(stderr,"Unable to load genome from file '%s'\n",pathname);
        soft_exit(1);
    }
    printf("%llds.\n", (timeInMillis() + 500 - start) / 1000);

    unsigned threadCount = GetNumberOfProcessors();
#ifdef _DEBUG
    threadCount = 1; // BJB
#endif // _DEBUG

    const char *lastDot = strchr(argv[2], '.');
    if (NULL != lastDot && !_stricmp(lastDot,".sam")) {
        readSupplierGenerator = SAMReader::createReadSupplierGenerator(argv[2], threadCount, genome);
    } else {
        readSupplierGenerator = FASTQReader::createReadSupplierGenerator(argv[2], threadCount);
    }

    if (NULL == readSupplierGenerator) {
        fprintf(stderr,"Unable to open file '%s' to get reads\n", argv[2]);
        soft_exit(1);
    }
    
    nRunningThreads = threadCount;
    DistHistogram *histograms = new DistHistogram[threadCount];
    CreateSingleWaiterObject(&threadsDone);

    for (unsigned i = 0; i < threadCount; i++) {
        StartNewThread(workerThreadMain, &histograms[i]);
    }

    WaitForSingleWaiterObject(&threadsDone);

    for (unsigned i = 1; i < threadCount; i++) {
        histograms[0].addIn(histograms[i]);
    }


    unsigned totalReads = 0;
    for (unsigned i = 0; i < DistHistogram::MaxDistance+1; i++) {
        printf("%d\t%d\n",i,histograms[0].counts[i]);
        totalReads += histograms[0].counts[i];
    }

    if (histograms[0].counts[DistHistogram::MaxDistance+1] != 0) {
        printf("More\t%d\n", histograms[0].counts[DistHistogram::MaxDistance+1]);
        totalReads += histograms[0].counts[DistHistogram::MaxDistance+1];
    }

    _int64 stop = timeInMillis();
    printf("\nProcessed %d reads in %llds, %lld reads/s\n", totalReads, (stop + 500 - start) / 1000, ((_int64) totalReads) * 1000 / (stop - start));

}

