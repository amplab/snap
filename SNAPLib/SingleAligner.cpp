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

extern GenomeLocation flt3itdLowerBound, flt3itdUpperBound;
const unsigned minMatchLength = 40;
const int maxDifferences = 3;

    bool
isFLT3ITDOneDirection(Read *read, GenomeIndex *genomeIndex, LandauVishkin<1> *lv, GenomeLocation *alignedLocation)
{
    *alignedLocation = -1;

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

        _int64					nHits[NUM_DIRECTIONS];      // Number of times this seed hits in the genome
        const GenomeLocation	*hits[NUM_DIRECTIONS];      // The actual hits (of size nHits)
		GenomeLocation			singletonHits[NUM_DIRECTIONS];
		const unsigned			*hits32[NUM_DIRECTIONS];

		if (genomeIndex->doesGenomeIndexHave64BitLocations()) {
			genomeIndex->lookupSeed(seed, flt3itdLowerBound, flt3itdUpperBound, &nHits[FORWARD], &hits[FORWARD], &nHits[RC], &hits[RC], &singletonHits[FORWARD], &singletonHits[RC]);
		} else {
			genomeIndex->lookupSeed32(seed, flt3itdLowerBound, flt3itdUpperBound, &nHits[FORWARD], &hits32[FORWARD], &nHits[RC], &hits32[RC]);
		}

        //
        // Since we're only going one direction, ignore the RC versions.
        //
        for (unsigned whichHit = 0; whichHit < nHits[FORWARD]; whichHit++) {
            int editDistance;
            //
            // Or at the end.
            //
			GenomeLocation hitLocation;
			if (genomeIndex->doesGenomeIndexHave64BitLocations()) {
				hitLocation = hits[FORWARD][whichHit];
			} else {
				hitLocation = hits32[FORWARD][whichHit];
			}
			editDistance = lv->computeEditDistance(genomeIndex->getGenome()->getSubstring(hitLocation + readLen - nextSeedToTest - minMatchLength, minMatchLength + maxDifferences), minMatchLength + maxDifferences, read->getData() + readLen - minMatchLength,
                __min(readLen, minMatchLength), maxDifferences);

            if (editDistance < bestDifference && editDistance >= 0) {
                bestDifference = editDistance;
				if (genomeIndex->doesGenomeIndexHave64BitLocations()) {
					*alignedLocation = hits[FORWARD][whichHit] - nextSeedToTest;
				} else {
					*alignedLocation = hits32[FORWARD][whichHit] - nextSeedToTest;
				}
            }

            //
            // See if it matches at the beginning of the read.
            //
			editDistance = lv->computeEditDistance(genomeIndex->getGenome()->getSubstring(hitLocation - nextSeedToTest, minMatchLength + maxDifferences), minMatchLength + maxDifferences, read->getData(), __min(readLen, minMatchLength), maxDifferences);

            if (editDistance < bestDifference && editDistance >= 0) {
				if (genomeIndex->doesGenomeIndexHave64BitLocations()) {
					*alignedLocation = hits[FORWARD][whichHit] - nextSeedToTest;
				} else {
					*alignedLocation = hits32[FORWARD][whichHit] - nextSeedToTest;
				}
				bestDifference = editDistance;
            }
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

    GenomeLocation rcLocation;
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
		result.correctAlignmentForSoftClipping(read, index->getGenome());

		for (int i = 0; i < nSecondaryResults; i++) {
			secondaryAlignments[i].correctAlignmentForSoftClipping(read, index->getGenome());
		}

#if     TIME_HISTOGRAM
        _int64 runTime = timeInNanos() - startTime;
        int timeBucket = min(30, cheezyLogBase2(runTime));
        stats->countByTimeBucket[timeBucket]++;
        stats->nanosByTimeBucket[timeBucket] += runTime;
#endif // TIME_HISTOGRAM

        allocator->checkCanaries();

		if (NotFound == result.status && isFLT3ITD(read, index, &lv, &result) || result.location >= flt3itdLowerBound && result.location <= flt3itdUpperBound) {
			extern ExclusiveLock BJBLock;
			AcquireExclusiveLock(&BJBLock);
#if 1
			if (result.direction == RC) {
				read->becomeRC();
			}

			unsigned editDistances[MAX_READ_LENGTH];
			for (int prefixSize = minMatchLength - maxDifferences; prefixSize < read->getDataLength(); prefixSize++) {
				editDistances[prefixSize] = lv.computeEditDistance(index->getGenome()->getSubstring(result.location, prefixSize + maxDifferences), prefixSize + maxDifferences, read->getData(), prefixSize, 62);
			}

			unsigned editDistanceDeltas[MAX_READ_LENGTH];
			for (int prefixSize = minMatchLength - maxDifferences + 1; prefixSize < read->getDataLength(); prefixSize++) {
				if (editDistances[prefixSize] == -1) {
					editDistanceDeltas[prefixSize] = 1;
				} else {
					editDistanceDeltas[prefixSize] = editDistances[prefixSize] - editDistances[prefixSize - 1];
				}
			}

			double biggestDifference = 0;
			GenomeLocation breakpoint = 0;
			for (_int64 prefixSize = minMatchLength - maxDifferences + 2; prefixSize < read->getDataLength() - 1; prefixSize++) {
				double prefixTotal = 0;
				double prefixCount = 0;
				double postfixTotal = 0;
				double postfixCount = 0;

				for (int i = minMatchLength - maxDifferences + 1; i < read->getDataLength(); i++) {
					if (i < prefixSize) {
						prefixCount++;
						prefixTotal += editDistanceDeltas[i];
					} else {
						postfixCount++;
						postfixTotal += editDistanceDeltas[i];
					}
				}

				double prefixAverage = prefixTotal / prefixCount;
				double postfixAverage = postfixTotal / postfixCount;

				if (postfixAverage - prefixAverage > biggestDifference) {
					biggestDifference = postfixAverage - prefixAverage;
					breakpoint = prefixSize + GenomeLocationAsInt64(result.location);
				}
			}

#if 0
			if (0 == GenomeLocationAsInt64(breakpoint)) {
				printf("Read aligned ITD with no breakpoint, possibly mapped to ITD rather than normal.  aligned to %lld, data %.*s\n", result.location, read->getDataLength(), read->getData());
			} else {
				printf("Possible breakpoint at %lld, first 5 bases %.*s, whole read %.*s\n", GenomeLocationAsInt64(breakpoint), 5, read->getData() + GenomeLocationAsInt64(breakpoint) - GenomeLocationAsInt64(result.location), read->getDataLength(), read->getData());
			}
#endif // 0

			if (result.direction == RC) {
				read->becomeRC();
			}
#endif	// 1

            writeRead(read, result, false);
            result.status = SingleHit;
			ReleaseExclusiveLock(&BJBLock);
        } else {
            result.status = NotFound;
            result.location = -1;
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
        readWriter->writeRead(read, result.status, result.mapq, result.location, result.direction, secondaryAlignment);
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

 
