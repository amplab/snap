// CheckChr.cpp : Defines the entry point for the console application.
//

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
    fprintf(stderr, "usage: CheckChr inputListFile outputFile\n");
    fprintf(stderr, "Look at each of the bam files in the input list to see whether they express chromosomes as 'chrN' or just 'N', and\n");
    fprintf(stderr, "write the result to outputFile\n");

    soft_exit(1);
}

int main(int argc, char* argv[])
{
    if (3 != argc) usage();

    FILE *inputListFile = fopen(argv[1], "r");
    if (NULL == inputListFile) {
        fprintf(stderr, "Unable to open input list file '%s'\n", argv[1]);
        soft_exit(1);
    }

    FILE *outputFile = fopen(argv[2], "w");
    if (NULL == outputFile) {
        fprintf(stderr, "Unable to open output file '%s'\n", argv[2]);
        soft_exit(1);
    }

    const size_t inputBufferSize = 4096;
    char inputBuffer[inputBufferSize];

    int nProcessed = 0;

    while (NULL != fgets(inputBuffer, inputBufferSize - 1, inputListFile)) {
        inputBuffer[inputBufferSize - 1] = '\0';    // Just to be safe

        //
        // Trim the \n that fgets leaves.
        //
        char *newline = strchr(inputBuffer, '\n');
        if (NULL != newline) {
            *newline = '\0';
        }

        DataSupplier::ThreadCount = 1;

        ReaderContext readerContext;
        readerContext.clipping = NoClipping;
        readerContext.defaultReadGroup = "";
        readerContext.genome = NULL;
        readerContext.ignoreSecondaryAlignments = true;
        readerContext.ignoreSupplementaryAlignments = true;
        readerContext.header = NULL;
        readerContext.headerLength = 0;
        readerContext.headerBytes = 0;

        BAMReader *reader = BAMReader::create(inputBuffer, 1, 0, 0, readerContext);

        bool sawAnyChr = false;
        for (int i = 0; i < reader->GetNRef(); i++) {
            const char *refName = reader->getRefName(i);
            if (strlen(refName) < 4) {
                continue;
            }

            if ((refName[0] == 'c' || refName[0] == 'C') && (refName[1] == 'h' || refName[1] == 'H') && (refName[2] == 'r' || refName[2] == 'R')) {
                sawAnyChr = true;
                break;
            }
        }

        //
        // Extract the analysis ID from the pathname.  It's the last directory in the path.
        //

        char *lastBackslash = strrchr(inputBuffer, '\\');
        if (NULL == lastBackslash) {
            fprintf(stderr, "Must use fully qualified pathnames for input files: '%s'\n", inputBuffer);
            soft_exit(1);
        }
        *lastBackslash = '\0';

        char *nextToLastBackslash = strrchr(inputBuffer, '\\');
        if (NULL == nextToLastBackslash) {
            fprintf(stderr, "Couldn't find analysis ID in pathname (truncated) '%s'\n", inputBuffer);
            soft_exit(1);
        }

        char *analysisID = nextToLastBackslash + 1;
        if (strlen(analysisID) != 36) {
            fprintf(stderr, "Directory name doesn't look like analysis ID (truncated) '%s'\n", inputBuffer);
            soft_exit(1);
        }


        fprintf(outputFile, "%s\t%s\n", analysisID, sawAnyChr ? "chr" : "noChr");

        delete reader;

        nProcessed++;
        if (nProcessed % 100 == 0) {
            printf(".");
            fflush(stdout);
        }
    }
}

