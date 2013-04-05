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
#include "MultiInputReadSupplier.h"

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
    const char *i_version,
    unsigned *argsConsumed)
{
    argc = i_argc;
    argv = i_argv;
    version = i_version;

    AlignerOptions* options = new AlignerOptions(
        "snap single <index-dir> <inputFile(s)> [<options>]"
		"   where <input file(s)> is a list of files to process.\n",false);
    options->extra = extension->extraOptions();
    if (argc < 2) {
        options->usage();
    }

    options->indexDir = argv[0];

	int nInputs = 0;
	for (int i = 1; i < argc; i++) {
		if (argv[i][0] == '-' || argv[i][0] == ',' && argv[i][1] == '\0') {
			break;
		}
		nInputs++;
	}

	if (0 == nInputs) {
		options->usage();
	}

	options->nInputs = nInputs;
	options->inputs = new SNAPInput[nInputs];

	for (int i = 1; i < argc; i++) {
		if (argv[i][0] == '-' || argv[i][0] == ',' && argv[i][1] == '\0') {
			break;
		}

		options->inputs[i-1].fileName = argv[i];
		options->inputs[i-1].fileType =
            stringEndsWith(argv[i],".sam") ? SAMFile :
            stringEndsWith(argv[i],".bam") ? BAMFile :
            stringEndsWith(argv[i], ".fastq.gz") || stringEndsWith(argv[i], ".fq.gz") ? GZipFASTQFile :
            FASTQFile;
	}

    int n;
    for (n = 1 + nInputs; n < argc; n++) {
        bool done;
        if (! options->parse(argv, argc, n, &done)) {
            options->usage();
        }

        if (done) {
            n++;    // for the ',' arg
            break;
        }
    }
    
    *argsConsumed = n;
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

    ReadSupplier *supplier = readSupplierGenerator->generateNewReadSupplier();
    if (NULL == supplier) {
        //
        // No work for this thread to do.
        //
        return;
    }

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
            NULL,               // reverse LV
            NULL,
            stats,
            allocator);

    allocator->assertAllMemoryUsed();
    allocator->checkCanaries();

    aligner->setExplorePopularSeeds(options->explorePopularSeeds);
    aligner->setStopOnFirstHit(options->stopOnFirstHit);

#ifdef  _MSC_VER
    if (options->useTimingBarrier) {
        if (0 == InterlockedDecrementAndReturnNewValue(nThreadsAllocatingMemory)) {
            AllowEventWaitersToProceed(memoryAllocationCompleteBarrier);
        } else {
            WaitForEvent(memoryAllocationCompleteBarrier);
        }
    }
#endif  // _MSC_VER

    // Align the reads.
    Read *read;
    while (NULL != (read = supplier->getNextRead())) {
        if (1 != selectivity && GoodFastRandom(selectivity-1) != 0) {
            //
            // Skip this read.
            //
            continue;
        }
        stats->totalReads++;

        // Skip the read if it has too many Ns or trailing 2 quality scores.
        if (read->getDataLength() < 50 || read->countOfNs() > maxDist) {
            if (readWriter != NULL && options->passFilter(read, NotFound)) {
                readWriter->writeRead(read, NotFound, 0, 0xFFFFFFFF, false);
            }
            continue;
        } else {
            stats->usefulReads++;
        }

        unsigned location = 0xFFFFFFFF;
        Direction direction;
        int score;
        int mapq;

        AlignmentResult result = aligner->AlignRead(read, &location, &direction, &score, &mapq);

        allocator->checkCanaries();

        bool wasError = false;
        if (result != NotFound && computeError) {
            wasError = wgsimReadMisaligned(read, location, index, options->misalignThreshold);
        }

        writeRead(read, result, location, direction, score, mapq);
        
        updateStats(stats, read, result, location, score, mapq, wasError);
    }
    stats->threadEntry->nReads = stats->usefulReads;        // This preserves the total across add()
    stats->threadEntry->lvCalls = aligner->getLocationsScored();

    aligner->~BaseAligner(); // This calls the destructor without calling operator delete, allocator owns the memory.
 
    if (supplier != NULL) {
        delete supplier;
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
    if (readWriter != NULL && options->passFilter(read, result)) {
        readWriter->writeRead(read, result, mapq, location, direction);
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

    void 
SingleAlignerContext::typeSpecificBeginIteration()
{
    if (1 == options->nInputs) {
        //
        // We've only got one input, so just connect it directly to the consumer.
        //
        readSupplierGenerator = options->inputs[0].createReadSupplierGenerator(options->numThreads, index->getGenome(), options->clipping);
    } else {
        //
        // We've got multiple inputs, so use a MultiInputReadSupplier to combine the individual inputs.
        //
        ReadSupplierGenerator **generators = new ReadSupplierGenerator *[options->nInputs];
        for (int i = 0; i < options->nInputs; i++) {
            generators[i] = options->inputs[i].createReadSupplierGenerator(options->numThreads, index->getGenome(), options->clipping);
        }
        readSupplierGenerator = new MultiInputReadSupplierGenerator(options->nInputs,generators);
    }
}
    void 
SingleAlignerContext::typeSpecificNextIteration()
{
     delete readSupplierGenerator;
     readSupplierGenerator = NULL;
}
