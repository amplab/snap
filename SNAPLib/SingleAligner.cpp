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
    stats->threadEntry->threadId = GetThreadId();
    stats->threadEntry->threadNumber = threadNum;

    int maxReadSize = 10000;
 
    BigAllocator *allocator = new BigAllocator(BaseAligner::getBigAllocatorReservation(true, maxHits, maxReadSize, index->getSeedLength(),numSeeds));
   
    BaseAligner *aligner = new (allocator) BaseAligner(
            index,
            confDiff,
            maxHits,
            maxDist,
            maxReadSize,
            numSeeds,
            adaptiveConfDiff,
            NULL,               // LV (no need to cache in the single aligner)
            similarityMap,
            stats,
            allocator);

    allocator->assertAllMemoryUsed();

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
            Direction direction;
            int score;
            int mapq;

            AlignmentResult result = aligner->AlignRead(&read, &location, &direction, &score, &mapq);

            bool wasError = false;
            if (result != NotFound && computeError) {
                wasError = wgsimReadMisaligned(&read, location, index, options->misalignThreshold);
            }

            //if (wasError) {
                writeRead(&read, result, location, direction, score, mapq);
            //}

            updateStats(stats, &read, result, location, score, mapq, wasError);
        }
    }
    stats->threadEntry->nReads = stats->usefulReads;        // This preserves the total across add()
    stats->threadEntry->lvCalls = aligner->getLocationsScored();

    aligner->~BaseAligner(); // This calls the destructor without calling operator delete, allocator owns the memory.
    if (reader != NULL) {
        delete reader;
    }

    delete allocator;   // This is what actually frees the memory.
}
    
    void
SingleAlignerContext::writeRead(
    Read* read,
    AlignmentResult result,
    unsigned location,
    Direction direction,
    int score,
    int mapq)
{
    if (samWriter != NULL && options->passFilter(read, result)) {
        samWriter->write(read, result, location, direction, mapq);
    }
}

    void
SingleAlignerContext::updateStats(
    AlignerStats* stats,
    Read* read,
    AlignmentResult result,
    unsigned location, 
    int score,
    int mapq,
    bool wasError)
{
    if (isOneLocation(result)) {
        stats->singleHits++;
        if (computeError) {
            stats->errors += wasError ? 1 : 0;
        }
    } else if (result == MultipleHits) {
        stats->multiHits++;
    } else {
        _ASSERT(result == NotFound);
        stats->notFound++;
    }

    if (result != NotFound) {
        _ASSERT(mapq >= 0 && mapq <= AlignerStats::maxMapq);
        stats->mapqHistogram[mapq]++;
        stats->mapqErrors[mapq] += wasError ? 1 : 0;
    }
}
