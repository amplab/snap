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

GenomeLocation flt3itdLowerBound, flt3itdUpperBound;
int main(int argc, const char **argv)
{
    if (3 != argc) usage();

	extern bool BigAllocUseHugePages;
	BigAllocUseHugePages = false;

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
    GenomeLocation chr13;
    if (!genome->getLocationOfContig("chr13", &chr13)) {
        WriteErrorMessage("Unable to find chr13.  Is this hg19?\n");
        soft_exit(1);
    }
    flt3itdLowerBound = chr13 + 28608000;
    flt3itdUpperBound = chr13 + 28609000;

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

	struct ReferenceLocation {
		ReferenceLocation() : mappedCount(0), mismatchCount(0), mismatchAs(0), mismatchTs(0), mismatchCs(0), mismatchGs(0) {}
		unsigned			mappedCount;
		unsigned			mismatchCount;
		unsigned			mismatchAs, mismatchTs, mismatchCs, mismatchGs;
	};

	ReferenceLocation *referenceLocations = new ReferenceLocation[flt3itdUpperBound - flt3itdLowerBound + 1];

    Read *read;
    unsigned nReads = 0;
    while (NULL != (read = readSupplier->getNextRead())) {
        if (nReads == maxReads) {
            fprintf(stderr,"Too many reads in input file\n");
            soft_exit(1);
        }

		if (read->getOriginalSAMFlags() & SAM_REVERSE_COMPLEMENT) {
			read->becomeRC();
		}

		for (unsigned readOffset = 0; readOffset < read->getDataLength(); readOffset++) {
			GenomeLocation location = read->getOriginalAlignedLocation() + readOffset;
			if (location >= flt3itdLowerBound && location <= flt3itdUpperBound) {
				referenceLocations[location - flt3itdLowerBound].mappedCount++;
				char readBase = read->getData()[readOffset];
				char referenceBase = *genome->getSubstring(location, 1);
				if (readBase != 'N' && referenceBase != 'N' && readBase != referenceBase) {
					referenceLocations[location - flt3itdLowerBound].mismatchCount++;
					switch (readBase) {
						case 'A': referenceLocations[location - flt3itdLowerBound].mismatchAs++; break;
						case 'T': referenceLocations[location - flt3itdLowerBound].mismatchTs++; break;
						case 'C': referenceLocations[location - flt3itdLowerBound].mismatchCs++; break;
						case 'G': referenceLocations[location - flt3itdLowerBound].mismatchGs++; break;
					}
				}
			}
		}

        reads[nReads] = new Read;
        reads[nReads]->copyFromOtherRead(*read);

        nReads++;
    }


	const _int64 lineWidth = 40;
	for (GenomeLocation lineStart = flt3itdLowerBound; lineStart <= flt3itdUpperBound; lineStart += lineWidth) {
		printf("%lld - %lld:\n", lineStart - chr13, __min(lineStart - chr13 + lineWidth, flt3itdUpperBound - chr13));
		for (GenomeLocation location = lineStart; location < lineStart + lineWidth && location <= flt3itdUpperBound; location++) {
			printf("   %c ", *genome->getSubstring(location, 1));
		}
		printf("\n");
		for (GenomeLocation location = lineStart; location < lineStart + lineWidth && location <= flt3itdUpperBound; location++) {
				printf("%4d ", referenceLocations[location - flt3itdLowerBound].mappedCount);
		}
		printf("\n");
		for (GenomeLocation location = lineStart; location < lineStart + lineWidth && location <= flt3itdUpperBound; location++) {
			printf(" %2d%% ", referenceLocations[location - flt3itdLowerBound].mismatchCount * 100 / referenceLocations[location - flt3itdLowerBound].mappedCount);
		}
		printf("\n");
		for (GenomeLocation location = lineStart; location < lineStart + lineWidth && location <= flt3itdUpperBound; location++) {
			char charToPrint;
			if (referenceLocations[location - flt3itdLowerBound].mismatchCount == 0) {
				charToPrint = ' ';
			}
			else {
				unsigned biggestCountSoFar = referenceLocations[location - flt3itdLowerBound].mismatchAs;
				charToPrint = 'A';
				if (referenceLocations[location - flt3itdLowerBound].mismatchTs > biggestCountSoFar) {
					charToPrint = 'T';
					biggestCountSoFar = referenceLocations[location - flt3itdLowerBound].mismatchTs;
				}
				if (referenceLocations[location - flt3itdLowerBound].mismatchCs > biggestCountSoFar) {
					charToPrint = 'C';
					biggestCountSoFar = referenceLocations[location - flt3itdLowerBound].mismatchCs;
				}
				if (referenceLocations[location - flt3itdLowerBound].mismatchGs > biggestCountSoFar) {
					charToPrint = 'G';
					biggestCountSoFar = referenceLocations[location - flt3itdLowerBound].mismatchGs;
				}

				if (biggestCountSoFar * 10 < referenceLocations[location - flt3itdLowerBound].mismatchCount * 8) {
					charToPrint += 'a' - 'A';
				}
			}
			printf("   %c ", charToPrint);
		}
		printf("\n\n");
	}
	
}
