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
    const char *i_version)
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
		if (argv[i][0] == '-') {
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
		if (argv[i][0] == '-') {
			break;
		}

		options->inputs[i-1].fileName = argv[i];
		options->inputs[i-1].fileType =
            stringEndsWith(argv[i],".sam") ? SAMFile :
            stringEndsWith(argv[i],".bam") ? BAMFile :
            stringEndsWith(argv[i], ".fastq.gz") || stringEndsWith(argv[i], ".fq.gz") ? GZipFASTQFile :
            FASTQFile;
	}

    for (int n = 1 + nInputs; n < argc; n++) {
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
    ReadSupplier *supplier = readSupplierGenerator->generateNewReadSupplier();
    if (NULL == supplier) {
        //
        // No work for this thread to do.
        //
        return;
    }
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
            if (samWriter != NULL && options->passFilter(read, NotFound)) {
                samWriter->write(read, NotFound, 0xFFFFFFFF, false);
            }
            continue;
        } else {
            stats->usefulReads++;
        }

        unsigned location = 0xFFFFFFFF;
        bool isRC;
        int score;
        AlignmentResult result = aligner->AlignRead(read, &location, &isRC, &score);

        writeRead(read, result, location, isRC, score);

        updateStats(stats, read, result, location, score);
    }
  

    delete aligner;
    if (supplier != NULL) {
        delete supplier;
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
    //delete readSupplierGenerator;
    //readSupplierGenerator = NULL;
}
