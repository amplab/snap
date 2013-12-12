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
#include "AlignerContext.h"
#include "AlignerOptions.h"
#include "FASTQ.h"
#include "Util.h"
#include "SingleAligner.h"
#include "MultiInputReadSupplier.h"
#include "AlignmentFilter.h"

using namespace std;
using util::stringEndsWith;

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
        "snapr single <genome-dir> <transcriptome-dir> <annotation> <inputFile(s)> [<options>]"
		"   where <input file(s)> is a list of files to process.\n",false);
    options->extra = extension->extraOptions();
    if (argc < 4) {
        fprintf(stderr,"Too few parameters\n");
        options->usage();
    }

    options->indexDir = argv[0];
    options->transcriptomeDir = argv[1];
    options->annotation = argv[2];

	int nInputs = 0;
	for (int i = 3; i < argc; i++) {
		if (argv[i][0] == '-' || argv[i][0] == ',' && argv[i][1] == '\0') {
			break;
		}
		nInputs++;
	}

	if (0 == nInputs) {
        fprintf(stderr,"Didn't see any input files\n");
		options->usage();
	}

	options->nInputs = nInputs;
	options->inputs = new SNAPInput[nInputs];

	for (int i = 3; i < argc; i++) {
		if (argv[i][0] == '-' || argv[i][0] == ',' && argv[i][1] == '\0') {
			break;
		}

		options->inputs[i-3].fileName = argv[i];
		options->inputs[i-3].fileType =
            stringEndsWith(argv[i],".sam") ? SAMFile :
            stringEndsWith(argv[i],".bam") ? BAMFile :
            stringEndsWith(argv[i], ".fastq.gz") || stringEndsWith(argv[i], ".fq.gz") ? GZipFASTQFile :
            FASTQFile;
	}

    int n;
    for (n = 3 + nInputs; n < argc; n++) {
        bool done;
        int oldN = n;
        if (!options->parse(argv, argc, n, &done)) {
            fprintf(stderr,"Didn't understand options starting at %s\n", argv[oldN]);
            options->usage();
        }

        if (done) {
            n++;    // for the ',' arg
            break;
        }
    }
    
    if (options->maxDist.end + options->extraSearchDepth >= MAX_K) {
        fprintf(stderr,"You specified too large of a maximum edit distance combined with extra search depth.  The must add up to less than %d.\n", MAX_K);
        fprintf(stderr,"Either reduce their sum, or change MAX_K in LandauVishkin.h and recompile.\n");
        soft_exit(1);
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
    ReadSupplier *supplier = readSupplierGenerator->generateNewReadSupplier();
    if (NULL == supplier) {
        //
        // No work for this thread to do.
        //
        return;
    }
	if (extension->runIterationThread(supplier, this)) {
		delete supplier;
		return;
	}
    if (index == NULL) {
        // no alignment, just input/output
        Read *read;
        while (NULL != (read = supplier->getNextRead())) {
            stats->totalReads++;
            writeRead(read, NotFound, InvalidGenomeLocation, FORWARD, false, 0, 0, 0);
        }
        delete supplier;
        return;
    }

    int maxReadSize = MAX_READ_LENGTH;
 
    BigAllocator *g_allocator = new BigAllocator(BaseAligner::getBigAllocatorReservation(true, maxHits, maxReadSize, index->getSeedLength(), numSeedsFromCommandLine, seedCoverage));

    BaseAligner *g_aligner = new (g_allocator) BaseAligner(
            index,
            maxHits,
            maxDist,
            maxReadSize,
            numSeedsFromCommandLine,
            seedCoverage,
            extraSearchDepth,
            NULL,               // LV (no need to cache in the single aligner)
            NULL,               // reverse LV
            stats,
            g_allocator);

    g_allocator->checkCanaries();
    g_aligner->setExplorePopularSeeds(options->explorePopularSeeds);
    g_aligner->setStopOnFirstHit(options->stopOnFirstHit);
    
    BigAllocator *t_allocator = new BigAllocator(BaseAligner::getBigAllocatorReservation(true, maxHits, maxReadSize, transcriptome->getSeedLength(), numSeedsFromCommandLine, seedCoverage));    
    BaseAligner *t_aligner = new (t_allocator) BaseAligner(
            transcriptome,
            maxHits,
            maxDist,
            maxReadSize,
            numSeedsFromCommandLine,
            seedCoverage,
            extraSearchDepth,
            NULL,               // LV (no need to cache in the single aligner)
            NULL,               // reverse LV
            stats,
            t_allocator);

    t_allocator->checkCanaries();
    t_aligner->setExplorePopularSeeds(options->explorePopularSeeds);
    t_aligner->setStopOnFirstHit(options->stopOnFirstHit);

    BigAllocator *c_allocator = NULL;
    BaseAligner *c_aligner = NULL;

    if (contamination != NULL) {

        c_allocator = new BigAllocator(BaseAligner::getBigAllocatorReservation(true, maxHits, maxReadSize, contamination->getSeedLength(), numSeedsFromCommandLine, seedCoverage));
        c_aligner = new (c_allocator) BaseAligner(
                contamination,
                maxHits,
                maxDist,
                maxReadSize,
                numSeedsFromCommandLine,
                seedCoverage,
                extraSearchDepth,
                NULL,               // LV (no need to cache in the single aligner)
                NULL,               // reverse LV
                stats,
                c_allocator);

        c_allocator->checkCanaries();
        c_aligner->setExplorePopularSeeds(options->explorePopularSeeds);
        c_aligner->setStopOnFirstHit(options->stopOnFirstHit);

    }

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
        stats->totalReads++;
        
        //Quality filtering
        bool quality = read->qualityFilter(options->minPercentAbovePhred, options->minPhred, options->phredOffset);

        // Skip the read if it has too many Ns or trailing 2 quality scores.
        if (read->getDataLength() < 50 || read->countOfNs() > maxDist || !quality) {
            if (readWriter != NULL && options->passFilter(read, NotFound)) {
                readWriter->writeRead(read, NotFound, 0, InvalidGenomeLocation, false, false, 0);
            }
            continue;
        } else {
            stats->usefulReads++;
        }

        unsigned location = InvalidGenomeLocation;
        Direction direction;
        int score;
        int mapq;

        //Set transcriptome and flag here
        bool isTranscriptome = false;
        unsigned tlocation = 0;

        AlignmentFilter filter(NULL, read, index->getGenome(), transcriptome->getGenome(), gtf, 0, 0, options->confDiff, options->maxDist.start, index->getSeedLength(), g_aligner);

        AlignmentResult t_result = t_aligner->AlignRead(read, &location, &direction, &score, &mapq);
        t_allocator->checkCanaries();
        filter.AddAlignment(location, direction, score, mapq, true, true);

        AlignmentResult g_result = g_aligner->AlignRead(read, &location, &direction, &score, &mapq);
        g_allocator->checkCanaries();
        filter.AddAlignment(location, direction, score, mapq, false, true);
        
        //Filter the results
        AlignmentResult result = filter.FilterSingle(&location, &direction, &score, &mapq, &isTranscriptome, &tlocation);

        //If the read is still unaligned
        if (result == NotFound) {

            //If the contamination database is present
            if (c_aligner != NULL) {

                AlignmentResult c_result = c_aligner->AlignRead(read, &location, &direction, &score, &mapq);
                c_allocator->checkCanaries();

                if (c_result != NotFound) {
                    c_filter->AddAlignment(location, string(read->getId(), read->getIdLength()), string(read->getData(), read->getDataLength()));

                }
            }
        }

        bool wasError = false;
        if (result != NotFound && computeError) {
            wasError = wgsimReadMisaligned(read, location, index, options->misalignThreshold);
        }

        writeRead(read, result, location, direction, isTranscriptome, tlocation, score, mapq);
        
        updateStats(stats, read, result, location, score, mapq, wasError);
    }

    g_aligner->~BaseAligner(); // This calls the destructor without calling operator delete, allocator owns the memory.
    t_aligner->~BaseAligner(); 
 
    if (supplier != NULL) {
        delete supplier;
    }

    delete g_allocator;   // This is what actually frees the memory.
    delete t_allocator;
}
    
    void
SingleAlignerContext::writeRead(
    Read* read,
    AlignmentResult result,
    unsigned location,
    Direction direction,
    bool isTranscriptome,
    unsigned tlocation,
    int score,
    int mapq)
{
    if (readWriter != NULL && options->passFilter(read, result)) {
        readWriter->writeRead(read, result, mapq, location, direction, isTranscriptome, tlocation);
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
        options->inputs[0].readHeader(readerContext);
        readSupplierGenerator = options->inputs[0].createReadSupplierGenerator(options->numThreads, readerContext);
    } else {
        //
        // We've got multiple inputs, so use a MultiInputReadSupplier to combine the individual inputs.
        //
        ReadSupplierGenerator **generators = new ReadSupplierGenerator *[options->nInputs];
        // use separate context for each supplier, initialized from common
        for (int i = 0; i < options->nInputs; i++) {
            ReaderContext context(readerContext);
            options->inputs[i].readHeader(context);
            generators[i] = options->inputs[i].createReadSupplierGenerator(options->numThreads, context);
        }
        readSupplierGenerator = new MultiInputReadSupplierGenerator(options->nInputs,generators);
    }
}
    void 
SingleAlignerContext::typeSpecificNextIteration()
    {
    if (readerContext.header != NULL) {
        delete [] readerContext.header;
        readerContext.header = NULL;
        readerContext.headerLength = readerContext.headerBytes = 0;
    }
    delete readSupplierGenerator;
    readSupplierGenerator = NULL;
}
