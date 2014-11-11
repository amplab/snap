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
#include "AlignmentFilter.h"

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

    SingleAlignmentResult *g_secondaryAlignments = NULL;
    unsigned g_secondaryAlignmentBufferCount;
    if (maxSecondaryAlignmentAdditionalEditDistance < 0) {
        g_secondaryAlignmentBufferCount = 0;
    } else {
        g_secondaryAlignmentBufferCount = BaseAligner::getMaxSecondaryResults(numSeedsFromCommandLine, seedCoverage, maxReadSize, maxHits, index->getSeedLength());
    }
    size_t g_secondaryAlignmentBufferSize = sizeof(*g_secondaryAlignments) * g_secondaryAlignmentBufferCount;
 
    BigAllocator *g_allocator = new BigAllocator(BaseAligner::getBigAllocatorReservation(true, maxHits, maxReadSize, index->getSeedLength(), numSeedsFromCommandLine, seedCoverage) + g_secondaryAlignmentBufferSize);
   
    BaseAligner *g_aligner = new (g_allocator) BaseAligner(
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
            g_allocator);


    if (maxSecondaryAlignmentAdditionalEditDistance >= 0) {
        g_secondaryAlignments = (SingleAlignmentResult *)g_allocator->allocate(g_secondaryAlignmentBufferSize);
    }

    g_allocator->checkCanaries();

    g_aligner->setExplorePopularSeeds(options->explorePopularSeeds);
    g_aligner->setStopOnFirstHit(options->stopOnFirstHit);
    
    SingleAlignmentResult *t_secondaryAlignments = NULL;
    unsigned t_secondaryAlignmentBufferCount;
    if (maxSecondaryAlignmentAdditionalEditDistance < 0) {
        t_secondaryAlignmentBufferCount = 0;
    } else {
        t_secondaryAlignmentBufferCount = BaseAligner::getMaxSecondaryResults(numSeedsFromCommandLine, seedCoverage, maxReadSize, maxHits, transcriptome->getSeedLength());
    }
    size_t t_secondaryAlignmentBufferSize = sizeof(*t_secondaryAlignments) * t_secondaryAlignmentBufferCount;

    BigAllocator *t_allocator = new BigAllocator(BaseAligner::getBigAllocatorReservation(true, maxHits, maxReadSize, transcriptome->getSeedLength(), numSeedsFromCommandLine, seedCoverage) + t_secondaryAlignmentBufferSize);

    BaseAligner *t_aligner = new (t_allocator) BaseAligner(
            transcriptome,
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
            t_allocator);


    if (maxSecondaryAlignmentAdditionalEditDistance >= 0) {
        t_secondaryAlignments = (SingleAlignmentResult *)t_allocator->allocate(t_secondaryAlignmentBufferSize);
    }

    t_allocator->checkCanaries();

    t_aligner->setExplorePopularSeeds(options->explorePopularSeeds);
    t_aligner->setStopOnFirstHit(options->stopOnFirstHit);

    BigAllocator *c_allocator = NULL;
    BaseAligner *c_aligner = NULL;
    SingleAlignmentResult *c_secondaryAlignments = NULL;
    unsigned c_secondaryAlignmentBufferCount;
    size_t c_secondaryAlignmentBufferSize;
    
    if (contamination != NULL) {

      if (maxSecondaryAlignmentAdditionalEditDistance < 0) {
          c_secondaryAlignmentBufferCount = 0;
      } else {
          c_secondaryAlignmentBufferCount = BaseAligner::getMaxSecondaryResults(numSeedsFromCommandLine, seedCoverage, maxReadSize, maxHits, contamination->getSeedLength());
      }  
      c_secondaryAlignmentBufferSize = sizeof(*c_secondaryAlignments) * c_secondaryAlignmentBufferCount;

      c_allocator = new BigAllocator(BaseAligner::getBigAllocatorReservation(true, maxHits, maxReadSize, contamination->getSeedLength(), numSeedsFromCommandLine, seedCoverage) + c_secondaryAlignmentBufferSize);

      c_aligner = new (c_allocator) BaseAligner(
              contamination,
              maxHits,
              maxDist,
              maxReadSize,
              numSeedsFromCommandLine,
              seedCoverage,
              extraSearchDepth,
              minWeightToCheck,
              noUkkonen,
              noOrderedEvaluation,
              noTruncation,
              NULL,               // LV (no need to cache in the single aligner)
              NULL,               // reverse LV
              stats,
              c_allocator);


      if (maxSecondaryAlignmentAdditionalEditDistance >= 0) {
          c_secondaryAlignments = (SingleAlignmentResult *)c_allocator->allocate(c_secondaryAlignmentBufferSize);
      }

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
    _uint64 lastReportTime = timeInMillis();
    _uint64 readsWhenLastReported = 0;

    while (NULL != (read = supplier->getNextRead())) {
        stats->totalReads++;
        
        //Quality filtering
        bool quality = read->qualityFilter(options->minPercentAbovePhred, options->minPhred, options->phredOffset);

        if (AlignerOptions::useHadoopErrorMessages && stats->totalReads % 10000 == 0 && timeInMillis() - lastReportTime > 10000) {
            fprintf(stderr,"reporter:counter:SNAP,readsAligned,%lu\n",stats->totalReads - readsWhenLastReported);
            readsWhenLastReported = stats->totalReads;
            lastReportTime = timeInMillis();
        }

        // Skip the read if it has too many Ns or trailing 2 quality scores.
        if (read->getDataLength() < minReadLength || read->countOfNs() > maxDist) {
            if (readWriter != NULL && options->passFilter(read, NotFound)) {
                readWriter->writeRead(readerContext, read, NotFound, 0, InvalidGenomeLocation, FORWARD, false, false, 0);
            }
            continue;
        } else {
            stats->usefulReads++;
        }

#if     TIME_HISTOGRAM
        _int64 startTime = timeInNanos();
#endif // TIME_HISTOGRAM

        SingleAlignmentResult result, contaminantResult;
        result.isTranscriptome = false;
        result.tlocation = 0;
        contaminantResult.isTranscriptome = false;
        contaminantResult.tlocation = 0;
        
        int nSecondaryResults = 0;

#ifdef LONG_READS
        int oldMaxK = aligner->getMaxK();
        if (options->maxDistFraction > 0.0) {
            t_aligner->setMaxK(min(MAX_K, (int)(read->getDataLength() * options->maxDistFraction)));
            g_aligner->setMaxK(min(MAX_K, (int)(read->getDataLength() * options->maxDistFraction)));
            if (c_aligner != NULL) {
              c_aligner->setMaxK(min(MAX_K, (int)(read->getDataLength() * options->maxDistFraction)));
            }
        }
#endif

        AlignmentFilter filter(NULL, read, index->getGenome(), transcriptome->getGenome(), gtf, 0, 0, options->confDiff, options->maxDist, index->getSeedLength(), g_aligner);

        t_aligner->AlignRead(read, &result, maxSecondaryAlignmentAdditionalEditDistance, t_secondaryAlignmentBufferCount, &nSecondaryResults, t_secondaryAlignments);

        t_allocator->checkCanaries();

        filter.AddAlignment(result.location, result.direction, result.score, result.mapq, true, true);
        for (int i = 0; i < nSecondaryResults; i++) {
          filter.AddAlignment(t_secondaryAlignments[i].location, t_secondaryAlignments[i].direction, t_secondaryAlignments[i].score, t_secondaryAlignments[i].mapq, true, true);
        }

        g_aligner->AlignRead(read, &result, maxSecondaryAlignmentAdditionalEditDistance, g_secondaryAlignmentBufferCount, &nSecondaryResults, g_secondaryAlignments);

        g_allocator->checkCanaries();

        filter.AddAlignment(result.location, result.direction, result.score, result.mapq, false, true);
        for (int i = 0; i < nSecondaryResults; i++) {
          filter.AddAlignment(g_secondaryAlignments[i].location, g_secondaryAlignments[i].direction, g_secondaryAlignments[i].score, g_secondaryAlignments[i].mapq, true, true);
        }

        //Filter the results
        AlignmentResult status = filter.FilterSingle(&result.location, &result.direction, &result.score, &result.mapq, &result.isTranscriptome, &result.tlocation);

        //If the read is still unaligned
        if (status == NotFound) {

            //If the contamination database is present
            if (c_aligner != NULL) {

                c_aligner->AlignRead(read, &contaminantResult, maxSecondaryAlignmentAdditionalEditDistance, c_secondaryAlignmentBufferCount, &nSecondaryResults, c_secondaryAlignments);
                c_allocator->checkCanaries();

                if (contaminantResult.status != NotFound) {
                    c_filter->AddAlignment(contaminantResult.location, string(read->getId(), read->getIdLength()), string(read->getData(), read->getDataLength()));

                }
            }
        }

#ifdef LONG_READS
        t_aligner->setMaxK(oldMaxK);
        g_aligner->setMaxK(oldMaxK);
        if (c_aligner != NULL) {
          c_aligner->setMaxK(oldMaxK);
        }
#endif

#if     TIME_HISTOGRAM
        _int64 runTime = timeInNanos() - startTime;
        int timeBucket = min(30, cheezyLogBase2(runTime));
        stats->countByTimeBucket[timeBucket]++;
        stats->nanosByTimeBucket[timeBucket] += runTime;
#endif // TIME_HISTOGRAM

        writeRead(read, result, false);

        /*
	for (int i = 0; i < __min(nSecondaryResults, maxSecondaryAlignments); i++) {
            writeRead(read, secondaryAlignments[i], true);
        }
        */   
    }

    g_aligner->~BaseAligner(); // This calls the destructor without calling operator delete, allocator owns the memory.
    t_aligner->~BaseAligner(); 
 
    if (supplier != NULL) {
        delete supplier;
    }

    delete g_allocator;   // This is what actually frees the memory.
    delete t_allocator;

    if (c_allocator != NULL) {
        c_aligner->~BaseAligner();
        delete c_allocator;
    }
}
    
    void
SingleAlignerContext::writeRead(
    Read* read,
    const SingleAlignmentResult &result,
    bool secondaryAlignment)
{
    if (readWriter != NULL && options->passFilter(read, result.status)) {
        readWriter->writeRead(readerContext, read, result.status, result.mapq, result.location, result.direction, secondaryAlignment, result.isTranscriptome, result.tlocation);
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

 
