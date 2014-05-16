/*++

Module Name:

    snap.cpp

Abstract:

    Main entry point for the snap binary. Calls into the indexer, single-end aligner or
    paired-end aligner functions defined in other source files in this directory.

Authors:

    Matei Zaharia & Bill Bolosky, February, 2012

Environment:

    User mode service.

Revision History:

    Adapted from cSNAP, which was in turn adapted from the scala prototype

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


using namespace std;

const char *SNAP_VERSION = "1.0dev.23"; 

static void usage()
{
    fprintf(stderr,
            "Usage: snapr <command> [<options>]\n"
            "Commands:\n"
            "   index         build a genome index\n"
            "   transcriptome build a transcriptome index\n"
            "   single        align single-end reads\n"
            "   paired        align paired-end reads\n"
            "Type a command without arguments to see its help.\n");
    soft_exit(1);
}

int main(int argc, const char **argv)
{
    printf("Welcome to SNAPR version %s.\n\n", SNAP_VERSION);

    InitializeSeedSequencers();

    if (argc < 2) {
        usage();
    } else if (strcmp(argv[1], "index") == 0) {
        GenomeIndex::runIndexer(argc - 2, argv + 2);
    } else if (strcmp(argv[1], "transcriptome") == 0) {
        GenomeIndex::runTranscriptomeIndexer(argc - 2, argv + 2);
    } else if (strcmp(argv[1], "single") == 0 || strcmp(argv[1], "paired") == 0) {
        for (int i = 1; i < argc; /* i is increased below */) {
            unsigned nArgsConsumed;
            if (strcmp(argv[i], "single") == 0) {
                SingleAlignerContext single;
                single.runAlignment(argc - (i + 1), argv + i + 1, SNAP_VERSION, &nArgsConsumed);
            } else if (strcmp(argv[i], "paired") == 0) {
                PairedAlignerContext paired;
                paired.runAlignment(argc - (i + 1), argv + i + 1, SNAP_VERSION, &nArgsConsumed);
            } else {
                fprintf(stderr, "Invalid command: %s\n\n", argv[i]);
                usage();
            }
            _ASSERT(nArgsConsumed > 0);
            i += nArgsConsumed + 1;  // +1 for single or paired
        }
    } else {
        fprintf(stderr, "Invalid command: %s\n\n", argv[1]);
        usage();
    }
}
