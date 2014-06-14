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

extern unsigned flt3itdLowerBound, flt3itdUpperBound;

    bool
isFLT3ITDOneDirection(Read *read, GenomeIndex *genomeIndex, LandauVishkin<1> *lv, unsigned *alignedLocation)
{
    *alignedLocation = -1;
    const unsigned minMatchLength = 50;
    const int maxDifferences = 3;

    int bestDifference = maxDifferences + 1;

    unsigned readLen = read->getDataLength();

    if (readLen < minMatchLength) {
        return false;
    }

    unsigned seedLen = genomeIndex->getSeedLength();
    for (unsigned nextSeedToTest = 0; nextSeedToTest < readLen - seedLen; nextSeedToTest++) {
        if (!Seed::DoesTextRepresentASeed(read->getData() + nextSeedToTest, seedLen)) {
            continue;
        }

        Seed seed(read->getData() + nextSeedToTest, seedLen);

        unsigned        nHits[NUM_DIRECTIONS];      // Number of times this seed hits in the genome
        const unsigned  *hits[NUM_DIRECTIONS];      // The actual hits (of size nHits)

        genomeIndex->lookupSeed(seed, flt3itdLowerBound, flt3itdUpperBound, &nHits[0], &hits[0], &nHits[1], &hits[1]);

        //
        // Since we're only going one direction, ignore the RC versions.
        //
        for (unsigned whichHit = 0; whichHit < nHits[FORWARD]; whichHit++) {
            int editDistance;
            //
            // Or at the end.
            //
            editDistance = lv->computeEditDistance(genomeIndex->getGenome()->getSubstring(hits[FORWARD][whichHit] + readLen - nextSeedToTest - minMatchLength, minMatchLength + maxDifferences), minMatchLength + maxDifferences, read->getData() + readLen - minMatchLength, 
                __min(readLen, minMatchLength), maxDifferences);

            if (editDistance < bestDifference && editDistance >= 0) {
                bestDifference = editDistance;
                *alignedLocation = hits[FORWARD][whichHit] - nextSeedToTest;
            }

            //
            // See if it matches at the beginning of the read.
            //
            editDistance = lv->computeEditDistance(genomeIndex->getGenome()->getSubstring(hits[FORWARD][whichHit] - nextSeedToTest, minMatchLength + maxDifferences), minMatchLength + maxDifferences, read->getData(), __min(readLen, minMatchLength), maxDifferences);

            if (editDistance < bestDifference && editDistance >= 0) {
                *alignedLocation = hits[FORWARD][whichHit] - nextSeedToTest;
                bestDifference = editDistance;
            }

        }
    }

    if (bestDifference <= maxDifferences) return true;

    //
    // Check for the pattern you'd expect from the end of 2812's ITD
    //
    const char *itdEnd = "AGTACTCATTATCTGAGGAG";
    const char *itdEnd2 = "ATTCTCTGAAATCAACGTAG";
    size_t itdEndLen = strlen(itdEnd);
    size_t itdEnd2Len = strlen(itdEnd2);
    for (int i = 0; i < readLen - itdEndLen; i++) {
        if (!memcmp(read->getData() + i, itdEnd, itdEndLen) || !memcmp(read->getData() + i, itdEnd2, itdEnd2Len)) {
            *alignedLocation = 1;
            return true;
        }
    }

    return bestDifference <= maxDifferences;
}

    bool
isFLT3ITD(Read *read, GenomeIndex *index, LandauVishkin<1> *lv, SingleAlignmentResult *result)
{
    bool forward = isFLT3ITDOneDirection(read, index, lv, &result->location);
    if (forward) {
        result->status = SingleHit;
        result->direction = FORWARD;
        if (result->location != 1) {
            return true;
        }
    }

    unsigned rcLocation;
    read->becomeRC();
    bool is = isFLT3ITDOneDirection(read, index, lv, &rcLocation);
    if (is) {
        result->status = SingleHit;
        result->location = rcLocation;
        result->direction = RC;
    } else if (!forward) {
        result->status = NotFound;
    }
    read->becomeRC();

    return is || forward;
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
    if (maxSecondaryAligmmentAdditionalEditDistance < 0) {
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
            extraSearchDepth,
            noUkkonen,
            noOrderedEvaluation,
            NULL,               // LV (no need to cache in the single aligner)
            NULL,               // reverse LV
            stats,
            allocator);

    if (maxSecondaryAligmmentAdditionalEditDistance >= 0) {
        secondaryAlignments = (SingleAlignmentResult *)allocator->allocate(secondaryAlignmentBufferSize);
    }

    LandauVishkin<1> lv;

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
        if (read->getDataLength() < 50 || read->countOfNs() > maxDist) {
            if (readWriter != NULL && options->passFilter(read, NotFound)) {
                //readWriter->writeRead(read, NotFound, 0, InvalidGenomeLocation, FORWARD, false);
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

        aligner->AlignRead(read, &result, maxSecondaryAligmmentAdditionalEditDistance, secondaryAlignmentBufferCount, &nSecondaryResults, secondaryAlignments);

#if     TIME_HISTOGRAM
        _int64 runTime = timeInNanos() - startTime;
        int timeBucket = min(30, cheezyLogBase2(runTime));
        stats->countByTimeBucket[timeBucket]++;
        stats->nanosByTimeBucket[timeBucket] += runTime;
#endif // TIME_HISTOGRAM

        allocator->checkCanaries();

        bool wasError = false;
        if (result.status != NotFound && computeError) {
            wasError = wgsimReadMisaligned(read, result.location, index, options->misalignThreshold);
        }

        if (NotFound == result.status && isFLT3ITD(read, index, &lv, &result) || result.location >= flt3itdLowerBound && result.location <= flt3itdUpperBound) {
            writeRead(read, result, false);
            result.status = SingleHit;
        } else {
            result.status = NotFound;
            result.location = -1;
        }

        updateStats(stats, read, result.status, result.location, result.score, result.mapq, wasError);
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
        readWriter->writeRead(read, result.status, result.mapq, result.location, result.direction, secondaryAlignment);
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

 
