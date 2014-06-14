/*++

Module Name:

    flt3ItdCaller.cpp

Abstract:

    Hacky code to try to call FLT3 ITDs from a set of reads already filtered by the FLT3-ITD SNAP version.

Authors:

    Bill Bolosky, June, 2014

Environment:

    User mode service.

Revision History:

   

--*/

#include "stdafx.h"
#include "options.h"
#include "FASTA.h"
#include "GenomeIndex.h"
#include "SingleAligner.h"
#include "PairedAligner.h"
#include "exit.h"
#include "SeedSequencer.h"
#include "AlignerOptions.h"
#include "SAM.h"


using namespace std;

static void usage()
{
    fprintf(stderr,
            "Usage: flt3ItdCaller genomeIndex inputFile\n");
    soft_exit(1);
}

unsigned hamming(const char *str1, const char *str2, size_t length)
{
    unsigned distance = 0;

    for (int i = 0; i < length; i++) {
        if (str1[i] != str2[i]) {
            distance++;
        }
    }

    return distance;
}

unsigned flt3itdLowerBound, flt3itdUpperBound;
int main(int argc, const char **argv)
{
    if (3 != argc) usage();

    size_t genomeFileNameLength = strlen(argv[1]) + 8;
    char *genomeFileName = new char[genomeFileNameLength];
    strcpy(genomeFileName, argv[1]);
    strcat(genomeFileName, "\\Genome");

    const Genome *genome = Genome::loadFromFile(genomeFileName, 500 /* told you it was hacky*/);

    if (NULL == genome) {
        fprintf(stderr,"Unable to open genome file '%s', aborting\n", genomeFileName);
        soft_exit(1);
    }

    //
    // Fill in the globals that tell us where to look for FLT3 ITDs.
    //
    unsigned chr13;
    if (!genome->getOffsetOfContig("chr13", &chr13)) {
        WriteErrorMessage("Unable to find chr13.  Is this hg19?\n");
        soft_exit(1);
    }
    flt3itdLowerBound = chr13 + 28608200;
    flt3itdUpperBound = chr13 + 28608500;

    const int maxReads = 100000;

    Read **reads = new Read*[maxReads];

    ReaderContext readerContext;
    readerContext.clipping = NoClipping;
    readerContext.defaultReadGroup = "";
    readerContext.genome = genome;
    readerContext.ignoreSecondaryAlignments = true;
    readerContext.ignoreSupplementaryAlignments = true;
	readerContext.header = NULL;
	readerContext.headerLength = 0;
	readerContext.headerBytes = 0;

    ReadSupplierGenerator *readSupplierGenerator = SAMReader::createReadSupplierGenerator(argv[2], 1, readerContext);
    if (NULL == readSupplierGenerator) {
        fprintf(stderr,"unable to open input file '%s'\n", argv[2]);
        soft_exit(1);
    }

    ReadSupplier *readSupplier = readSupplierGenerator->generateNewReadSupplier();

    Read *read;
    unsigned nReads = 0;
    while (NULL != (read = readSupplier->getNextRead())) {
        if (nReads == maxReads) {
            fprintf(stderr,"Too many reads in input file\n");
            soft_exit(1);
        }

        if (read->getOriginalAlignedLocation() == 1) {
            //
            // Ignore these, they don't map at all.
            //
            continue;
        }

        reads[nReads] = new Read;
        reads[nReads]->copyFromOtherRead(*read);

        nReads++;
    }

    unsigned *bestAlignments = new unsigned[maxReads];  // The lowest hamming distance alignment for each read against the 

}
