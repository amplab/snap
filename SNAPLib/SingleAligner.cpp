/*++

Module Name:

    SingleAligner.cpp

Abstract:

    Functions for running the single end aligner sub-program.

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
#include "BaseAligner.h"
#include "Compat.h"
#include "RangeSplitter.h"
#include "GenomeIndex.h"
#include "Range.h"
#include "SAM.h"
#include "Tables.h"
#include "WGsim.h"
#include "GoodRandom.h"
#include "AlignerContext.h"
#include "AlignerOptions.h"
#include "FASTQ.h"
#include "SingleAligner.h"

using namespace std;

// Check whether a string str ends with a given pattern
static bool stringEndsWith(const char* str, const char* pattern) {
    if (strlen(str) < strlen(pattern)) {
        return false;
    } else {
        return strcmp(str + (strlen(str) - strlen(pattern)), pattern) == 0;
    }
}

SingleAlignerContext::SingleAlignerContext(AlignerExtension* i_extension)
    : AlignerContext(0, NULL, NULL, i_extension)
{
}

    AlignerOptions*
SingleAlignerContext::parseOptions(
    int i_argc,
    const char **i_argv,
    const char *i_version)
{
    argc = i_argc;
    argv = i_argv;
    version = i_version;

    AlignerOptions* options = new AlignerOptions(
        "snap single <index-dir> <inputFile> [-o output.sam] [<options>]",false);
    options->extra = extension->extraOptions();
    if (argc < 2) {
        options->usage();
    }

    options->indexDir = argv[0];
    options->inputFilename = argv[1];
    options->inputFileIsFASTQ = !stringEndsWith(argv[1], ".sam");

    for (int n = 2; n < argc; n++) {
        if (! options->parse(argv, argc, n)) {
            options->usage();
        }
    }
    
    return options;
}

    AlignerStats*
SingleAlignerContext::newStats()
{
    return new AlignerStats();
}

    void
SingleAlignerContext::runTask()
{
    ParallelTask<SingleAlignerContext> task(this);
    task.run();
}
    
    void
SingleAlignerContext::runIterationThread()
{
    int maxReadSize = 10000;
    int lvLimit = 1000000;
    BaseAligner *aligner = new BaseAligner(
            index,
            confDiff,
            maxHits,
            maxDist,
            maxReadSize,
            numSeeds,
            lvLimit,
            adaptiveConfDiff);
    if (aligner == NULL) {
        fprintf(stderr, "Failed to create aligner!\n");
        exit(1);
    }
    aligner->setExplorePopularSeeds(options->explorePopularSeeds);
    aligner->setStopOnFirstHit(options->stopOnFirstHit);

    // Keep grabbing ranges of the file and processing them.
    ReadReader *reader = NULL;
    _int64 rangeStart, rangeLength;
    while (fileSplitter->getNextRange(&rangeStart, &rangeLength)) {
        if (NULL == reader) {
            if (inputFileIsFASTQ) {
                reader = FASTQReader::create(inputFilename, rangeStart, rangeLength, clipping);
            } else {
                reader = SAMReader::create(inputFilename, index->getGenome(), rangeStart, rangeLength, clipping);
            }
            if (NULL == reader) {
                fprintf(stderr, "Failed to create input file reader for %s.\n", inputFilename);
                exit(1);
            }
        } else {
            reader->reinit(rangeStart, rangeLength);
        }

        // Align the reads.
        Read read(reader);
        while (reader->getNextRead(&read)) {
            if (1 != selectivity && GoodFastRandom(selectivity-1) != 0) {
                //
                // Skip this read.
                //
                continue;
            }
            stats->totalReads++;

            // Skip the read if it has too many Ns or trailing 2 quality scores.
            if (read.getDataLength() < 50 || read.countOfNs() > maxDist) {
                if (samWriter != NULL && options->passFilter(&read, NotFound)) {
                    samWriter->write(&read, NotFound, 0xFFFFFFFF, false);
                }
                continue;
            } else {
                stats->usefulReads++;
            }

            unsigned location = 0xFFFFFFFF;
            bool isRC;
            int score;
            AlignmentResult result = aligner->AlignRead(&read, &location, &isRC, &score);

            writeRead(&read, result, location, isRC, score);

            updateStats(stats, &read, result, location, score);
        }
    }

    delete aligner;
    if (reader != NULL) {
        delete reader;
    }
}
    
    void
SingleAlignerContext::writeRead(
    Read* read,
    AlignmentResult result,
    unsigned location,
    bool isRC,
    int score)
{
    if (samWriter != NULL && options->passFilter(read, result)) {
        samWriter->write(read, result, location, isRC);
    }
}

    void
SingleAlignerContext::updateStats(
    AlignerStats* stats,
    Read* read,
    AlignmentResult result,
    unsigned location, 
    int score)
{
    if (isOneLocation(result)) {
        stats->singleHits++;
        if (computeError) {
            if (wgsimReadMisaligned(read, location, index, maxDist)) {
                stats->errors++;
            }
        }
    } else if (result == MultipleHits) {
        stats->multiHits++;
    } else {
        _ASSERT(result == NotFound);
        stats->notFound++;
    }
}
