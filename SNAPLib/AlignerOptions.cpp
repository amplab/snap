/*++

Module Name:

    AlignerOptions.cpp

Abstract:

    Common parameters for running single & paired alignment.

Authors:

    Ravi Pandya, May, 2012

Environment:
`
    User mode service.

Revision History:

    Integrated from SingleAligner.cpp & PairedAligner.cpp

--*/

#include "stdafx.h"
#include "options.h"
#include "AlignerOptions.h"

AlignerOptions::AlignerOptions(
    const char* i_commandLine,
    bool forPairedEnd)
    :
    commandLine(i_commandLine),
    indexDir(NULL),
    numThreads(1),
    computeError(false),
    bindToProcessors(false),
    ignoreMismatchedIDs(false),
    selectivity(1),
    samFileTemplate(NULL),
    doAlignerPrefetch(false),
    inputFilename(NULL),
    inputFileIsFASTQ(true),
    clipping(ClipBack),
    sortOutput(false),
    sortMemory(0),
    filterFlags(0),
    explorePopularSeeds(false),
    stopOnFirstHit(false),
	useM(false),
    gapPenalty(0),
	extra(NULL),
    rgLineContents(NULL)
{
    if (forPairedEnd) {
        maxDist             = 15;
        numSeeds            = 25;
        maxHits             = 250;
        confDiff            = 1;
        adaptiveConfDiff    = 7;
    } else {
        maxDist             = 8;
        numSeeds            = 25;
        maxHits             = 250;
        confDiff            = 2;
        adaptiveConfDiff    = 4;
    }
}

    void
AlignerOptions::usage()
{
    usageMessage();
    exit(1);
}

    void
AlignerOptions::usageMessage()
{
    fprintf(stderr,
        "Usage: %s\n"
        "Options:\n"
        "  -o   output alignments to a given SAM file\n"
#ifndef _MSC_VER
        "       (note: will create one output file per thread)\n"  // Linux guys: you should fix up your SAM writer!
#endif  // _MSC_VER 
        "  -d   maximum edit distance allowed per read or pair (default: %d)\n"
        "  -n   number of seeds to use per read (default: %d)\n"
        "  -h   maximum hits to consider per seed (default: %d)\n"
        "  -c   confidence threshold (default: %d)\n"
        "  -a   confidence adaptation threshold (default: %d)\n"
        "  -t   number of threads\n"
        "  -b   bind each thread to its processor (off by default)\n"
        "  -e   compute error rate assuming wgsim-generated reads\n"
        "  -P   disables cache prefetching in the genome; may be helpful for machines\n"
        "       with small caches or lots of cores/cache\n"
        "  -so  sort output file by alignment location\n"
        "  -sm  memory to use for sorting in Gb\n"
        "  -x   explore some hits of overly popular seeds (useful for filtering)\n"
        "  -f   stop on first match within edit distance limit (filtering mode)\n"
        "  -F   filter output (a=aligned only, s=single hit only, u=unaligned only)\n"
#if     USE_DEVTEAM_OPTIONS
        "  -I   ignore IDs that don't match in the paired-end aligner\n"
        "  -S   selectivity; randomly choose 1/selectivity of the reads to score\n"
#endif  // USE_DEVTEAM_OPTIONS
        "  -Cxx must be followed by two + or - symbols saying whether to clip low-quality\n"
        "       bases from front and back of read respectively; default: back only (-C-+)\n"
		"  -M   indicates that CIGAR strings in the generated SAM file should use M (alignment\n"
		"       match) rather than = and X (sequence (mis-)match)\n"
        "  -G   specify a gap penalty to use when generating CIGAR strings\n"
// not written yet        "  -r   Specify the content of the @RG line in the SAM header.\n"
            ,
            commandLine,
            maxDist.start,
            numSeeds.start,
            maxHits.start,
            confDiff.start,
            adaptiveConfDiff.start);
    if (extra != NULL) {
        extra->usageMessage();
    }
}

    bool
AlignerOptions::parse(
    const char** argv,
    int argc,
    int& n)
{
    if (strcmp(argv[n], "-d") == 0) {
        if (n + 1 < argc) {
            maxDist = Range::parse(argv[n+1]);
            n++;
            return true;
        }
    } else if (strcmp(argv[n], "-n") == 0) {
        if (n + 1 < argc) {
            numSeeds = Range::parse(argv[n+1]);
            n++;
            return true;
        }
    } else if (strcmp(argv[n], "-h") == 0) {
        if (n + 1 < argc) {
            maxHits = Range::parse(argv[n+1]);
            n++;
            return true;
        }
    } else if (strcmp(argv[n], "-c") == 0) {
        if (n + 1 < argc) {
            confDiff = Range::parse(argv[n+1]);
            n++;
            return true;
        }
    } else if (strcmp(argv[n], "-a") == 0) {
        if (n + 1 < argc) {
            adaptiveConfDiff = Range::parse(argv[n+1]);
            n++;
            return true;
        }
    } else if (strcmp(argv[n], "-t") == 0) {
        if (n + 1 < argc) {
            numThreads = atoi(argv[n+1]);
            n++;
            return true;
        }
    } else if (strcmp(argv[n], "-o") == 0) {
        if (n + 1 < argc) {
            samFileTemplate = argv[n+1];
            n++;
            return true;
        }
    } else if (strcmp(argv[n], "-e") == 0) {
        computeError = true;
        return true;
    } else if (strcmp(argv[n], "-P") == 0) {
        doAlignerPrefetch = false;
        return true;
    } else if (strcmp(argv[n], "-b") == 0) {
        bindToProcessors = true;
        return true;
    } else if (strcmp(argv[n], "-so") == 0) {
        sortOutput = true;
        return true;
    } else if (strcmp(argv[n], "-sm") == 0) {
        if (n + 1 < argc && argv[n+1][0] >= '0' && argv[n+1][0] <= '9') {
            sortMemory = atoi(argv[n+1]);
            n++;
            return true;
        }
    } else if (strcmp(argv[n], "-F") == 0) {
        if (n + 1 < argc) {
            n++;
            if (strcmp(argv[n], "a") == 0) {
                filterFlags = FilterSingleHit | FilterMultipleHits;
            } else if (strcmp(argv[n], "s") == 0) {
                filterFlags = FilterSingleHit;
            } else if (strcmp(argv[n], "u") == 0) {
                filterFlags = FilterUnaligned;
            } else {
                return false;
            }
            return true;
        }
    } else if (strcmp(argv[n], "-x") == 0) {
        explorePopularSeeds = true;
        return true;
    } else if (strcmp(argv[n], "-f") == 0) {
        stopOnFirstHit = true;
        return true;
#if     USE_DEVTEAM_OPTIONS
    } else if (strcmp(argv[n], "-I") == 0) {
        ignoreMismatchedIDs = true;
        return true;
    } else if (strcmp(argv[n], "-S") == 0) {
        if (n + 1 < argc) {
            selectivity = atoi(argv[n+1]);
            if (selectivity < 2) {
                fprintf(stderr,"Selectivity must be at least 2.\n");
                exit(1);
            }
            n++;
            return true;
        } else {
            fprintf(stderr,"Must have the selectivity value after -S\n");
        }
#endif  // USE_DEVTEAM_OPTIONS
	} else if (strcmp(argv[n], "-M") == 0) {
		useM = true;
		return true;
	} else if (strcmp(argv[n], "-G") == 0) {
        if (n + 1 < argc) {
            gapPenalty = atoi(argv[n+1]);
            if (gapPenalty < 1) {
                fprintf(stderr,"Gap penalty must be at least 1.\n");
                exit(1);
            }
            n++;
            return true;
        } else {
            fprintf(stderr,"Must have the gap penalty value after -G\n");
        }
    } else if (strcmp(argv[n], "-r") == 0) {
#if 0   // This isn't ready yet.
        if (n + 1 < argc) {
            //
            // Check the line for sanity.  It must consist either of @RG\t<fields> or just <fields> (in which
            // case we add the @RG part).  It must contain a field called ID. Fields are separated by tabs.
            // We don't require that the fields be things that are listed in the SAM spec, however, because
            // new ones might be added.
            //
            const char *parsePointer = argv[n+1];
            if (*parsePointer == '@') {
            }

            rgLineContents = argv[n+1];
            n++;
            return true;
        }
#endif  // 0
    } else if (strlen(argv[n]) >= 2 && '-' == argv[n][0] && 'C' == argv[n][1]) {
        if (strlen(argv[n]) != 4 || '-' != argv[n][2] && '+' != argv[n][2] ||
            '-' != argv[n][3] && '+' != argv[n][3]) {

            fprintf(stderr,"Invalid -C argument.\n\n");
            return false;
        }

        if ('-' == argv[n][2]) {
            if ('-' == argv[n][3]) {
                clipping = NoClipping;
            } else {
                clipping = ClipBack;
            }
        } else {
            if ('-' == argv[n][3]) {
                clipping = ClipFront;
            } else {
                clipping = ClipFrontAndBack;
            }
        }
        return true;
    } else if (extra != NULL) {
        return extra->parse(argv, argc, n);
    }
    return false;
}

    bool
AlignerOptions::passFilter(
    Read* read,
    AlignmentResult result)
{
    if (filterFlags == 0) {
        return true;
    }
    switch (result) {
    case NotFound:
    case UnknownAlignment:
        return (filterFlags & FilterUnaligned) != 0;
    case SingleHit:
    case CertainHit:
        return (filterFlags & FilterSingleHit) != 0;
    case MultipleHits:
        return (filterFlags & FilterMultipleHits) != 0;
    default:
        return false; // shouldn't happen!
    }
}
