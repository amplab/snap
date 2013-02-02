/*++

Module Name:

    snap.cpp

Abstract:

    Main entry point for the snap binary. Calls into the indexer, single-end aligner or
    paired-end aligner functions defined in other source files in this directory.

Authors:

    Matei Zaharia, February, 2012

Environment:
`
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

using namespace std;

const char *SNAP_VERSION = "0.15.4";

static void usage()
{
    fprintf(stderr,
            "Usage: snap <command> [<options>]\n"
            "Commands:\n"
            "   index    build a genome index\n"
            "   single   align single-end reads\n"
            "   paired   align paired-end reads\n"
            "Type a command without arguments to see its help.\n");
    exit(1);
}


int main(int argc, const char **argv)
{
    printf("Welcome to SNAP version %s.\n\n", SNAP_VERSION);

    if (argc < 2) {
        usage();
    } else if (strcmp(argv[1], "index") == 0) {
        GenomeIndex::runIndexer(argc - 2, argv + 2);
    } else if (strcmp(argv[1], "single") == 0) {
        SingleAlignerContext single;
        single.runAlignment(argc - 2, argv + 2, SNAP_VERSION);
    } else if (strcmp(argv[1], "paired") == 0) {
        PairedAlignerContext paired;
        paired.runAlignment(argc - 2, argv + 2, SNAP_VERSION);
    } else {
        fprintf(stderr, "Invalid command: %s\n\n", argv[1]);
        usage();
    }
}
