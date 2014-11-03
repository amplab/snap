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
#include "SAM.h"
#include "Tables.h"
#include "AlignerContext.h"
#include "AlignerOptions.h"
#include "FASTQ.h"
#include "Util.h"
#include "SingleAligner.h"
#include "MultiInputReadSupplier.h"

using namespace std;
using util::stringEndsWith;

SingleAlignerContext::SingleAlignerContext(AlignerExtension* i_extension)
    : AlignerContext(0, NULL, NULL, i_extension)
{
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
	PreventMachineHibernationWhileThisThreadIsAlive();

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
            SingleAlignmentResult result;
            result.status = NotFound;
            result.direction = FORWARD;
            result.mapq = 0;
            result.score = 0;
            result.location = InvalidGenomeLocation;
            writeRead(read, result, false);
        }
        delete supplier;
        return;
    }

    int maxReadSize = MAX_READ_LENGTH;

    SingleAlignmentResult *secondaryAlignments = NULL;
    unsigned secondaryAlignmentBufferCount;
    if (maxSecondaryAlignmentAdditionalEditDistance < 0) {
        secondaryAlignmentBufferCount = 0;
    } else {
        secondaryAlignmentBufferCount = BaseAligner::getMaxSecondaryResults(numSeedsFromCommandLine, seedCoverage, maxReadSize, maxHits, index->getSeedLength());
    }
    size_t secondaryAlignmentBufferSize = sizeof(*secondaryAlignments) * secondaryAlignmentBufferCount;
 
    BigAllocator *allocator = new BigAllocator(BaseAligner::getBigAllocatorReservation(true, maxHits, maxReadSize, index->getSeedLength(), numSeedsFromCommandLine, seedCoverage) + secondaryAlignmentBufferSize);
   
    BaseAligner *aligner = new (allocator) BaseAligner(
            index,
            maxHits,
            maxDist,
            maxReadSize,
            numSeedsFromCommandLine,
            seedCoverage,
	    minWeightToCheck,
            extraSearchDepth,
            noUkkonen,
            noOrderedEvaluation,
			noTruncation,
            NULL,               // LV (no need to cache in the single aligner)
            NULL,               // reverse LV
            stats,
            allocator);

    if (maxSecondaryAlignmentAdditionalEditDistance >= 0) {
        secondaryAlignments = (SingleAlignmentResult *)allocator->allocate(secondaryAlignmentBufferSize);
    }

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
    _uint64 lastReportTime = timeInMillis();
    _uint64 readsWhenLastReported = 0;

    while (NULL != (read = supplier->getNextRead())) {
        stats->totalReads++;

        if (AlignerOptions::useHadoopErrorMessages && stats->totalReads % 10000 == 0 && timeInMillis() - lastReportTime > 10000) {
            fprintf(stderr,"reporter:counter:SNAP,readsAligned,%lu\n",stats->totalReads - readsWhenLastReported);
            readsWhenLastReported = stats->totalReads;
            lastReportTime = timeInMillis();
        }

        // Skip the read if it has too many Ns or trailing 2 quality scores.
        if (read->getDataLength() < minReadLength || read->countOfNs() > maxDist) {
            if (readWriter != NULL && options->passFilter(read, NotFound)) {
                readWriter->writeRead(readerContext, read, NotFound, 0, InvalidGenomeLocation, FORWARD, false);
            }
            continue;
        } else {
            stats->usefulReads++;
        }

#if     TIME_HISTOGRAM
        _int64 startTime = timeInNanos();
#endif // TIME_HISTOGRAM

        SingleAlignmentResult result;
        int nSecondaryResults = 0;

#ifdef LONG_READS
        int oldMaxK = aligner->getMaxK();
        if (options->maxDistFraction > 0.0) {
            aligner->setMaxK(min(MAX_K, (int)(read->getDataLength() * options->maxDistFraction)));
        }
#endif

		aligner->AlignRead(read, &result, maxSecondaryAlignmentAdditionalEditDistance, secondaryAlignmentBufferCount, &nSecondaryResults, secondaryAlignments);
#ifdef LONG_READS
        aligner->setMaxK(oldMaxK);
#endif

#if     TIME_HISTOGRAM
        _int64 runTime = timeInNanos() - startTime;
        int timeBucket = min(30, cheezyLogBase2(runTime));
        stats->countByTimeBucket[timeBucket]++;
        stats->nanosByTimeBucket[timeBucket] += runTime;
#endif // TIME_HISTOGRAM

        allocator->checkCanaries();


		writeRead(read, result, false);

		for (int i = 0; i < __min(nSecondaryResults, maxSecondaryAlignments); i++) {
            writeRead(read, secondaryAlignments[i], true);
        }
        
        updateStats(stats, read, result.status, result.score, result.mapq);
    }

    aligner->~BaseAligner(); // This calls the destructor without calling operator delete, allocator owns the memory.
 
    if (supplier != NULL) {
        delete supplier;
    }

    delete allocator;   // This is what actually frees the memory.
}
    
    void
SingleAlignerContext::writeRead(
    Read* read,
    const SingleAlignmentResult &result,
    bool secondaryAlignment
    )
{
    if (readWriter != NULL && options->passFilter(read, result.status)) {
        readWriter->writeRead(readerContext, read, result.status, result.mapq, result.location, result.direction, secondaryAlignment);
    }
}

    void
SingleAlignerContext::updateStats(
    AlignerStats* stats,
    Read* read,
    AlignmentResult result,
    int score,
    int mapq)
{
    if (isOneLocation(result)) {
        stats->singleHits++;
    } else if (result == MultipleHits) {
        stats->multiHits++;
    } else {
        _ASSERT(result == NotFound);
        stats->notFound++;
    }

    if (result != NotFound) {
        _ASSERT(mapq >= 0 && mapq <= AlignerStats::maxMapq);
        stats->mapqHistogram[mapq]++;
    }
}

    void 
SingleAlignerContext::typeSpecificBeginIteration()
{
    if (1 == options->nInputs) {
        //
        // We've only got one input, so just connect it directly to the consumer.
        //
        readSupplierGenerator = options->inputs[0].createReadSupplierGenerator(options->numThreads, readerContext);
    } else {
        //
        // We've got multiple inputs, so use a MultiInputReadSupplier to combine the individual inputs.
        //
        ReadSupplierGenerator **generators = new ReadSupplierGenerator *[options->nInputs];
        // use separate context for each supplier, initialized from common
        for (int i = 0; i < options->nInputs; i++) {
            ReaderContext context(readerContext);
            generators[i] = options->inputs[i].createReadSupplierGenerator(options->numThreads, context);
        }
        readSupplierGenerator = new MultiInputReadSupplierGenerator(options->nInputs,generators);
    }
    ReaderContext* context = readSupplierGenerator->getContext();
    readerContext.header = context->header;
    readerContext.headerBytes = context->headerBytes;
    readerContext.headerLength = context->headerLength;
    readerContext.headerMatchesIndex = context->headerMatchesIndex;
}
    void 
SingleAlignerContext::typeSpecificNextIteration()
    {
    if (readerContext.header != NULL) {
        delete [] readerContext.header;
        readerContext.header = NULL;
        readerContext.headerLength = readerContext.headerBytes = 0;
        readerContext.headerMatchesIndex = false;
    }
    delete readSupplierGenerator;
    readSupplierGenerator = NULL;
}

 
