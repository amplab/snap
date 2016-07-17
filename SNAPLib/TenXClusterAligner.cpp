/*++

Module Name:

	TenXClusterAligner.cpp

Abstract:

	A paired-end aligner calls into a different paired-end aligner, and if
	it fails to find an alignment, aligns each of the reads singly.  This handles
	chimeric reads that would otherwise be unalignable.

Authors:

	Bill Bolosky, June, 2013

Environment:

	User mode service.

Revision History:

--*/


#include "stdafx.h"
#include "TenXClusterAligner.h"
#include "mapq.h"
#include "directions.h"
#include "BigAlloc.h"
#include "Util.h"

using namespace std;

#ifdef TRACE_PAIRED_ALIGNER
#define TRACE printf
#else
#define TRACE(...) {}
#endif

TenXClusterAligner::TenXClusterAligner(
	GenomeIndex         *index_,
	unsigned            maxReadSize,
	unsigned            maxHits,
	unsigned            maxK,
	unsigned            maxSeedsFromCommandLine,
	double              seedCoverage,
	unsigned            minWeightToCheck,
	bool                forceSpacing_,
	unsigned            extraSearchDepth,
	bool                noUkkonen,
	bool                noOrderedEvaluation,
	bool				noTruncation,
	bool                ignoreAlignmentAdjustmentsForOm,
	TenXSingleAligner    **underlyingTenXSingleAligner_,
	unsigned			maxBarcodeSize_,
	unsigned            minReadLength_,
	int                 maxSecondaryAlignmentsPerContig,
	BigAllocator        *allocator)
	: underlyingTenXSingleAligner(underlyingTenXSingleAligner_), maxBarcodeSize(maxBarcodeSize_), forceSpacing(forceSpacing_), index(index_), minReadLength(minReadLength_)
{
	// Create single-end aligners.
	singleAligner = new (allocator) BaseAligner(index, maxHits, maxK, maxReadSize,
		maxSeedsFromCommandLine, seedCoverage, minWeightToCheck, extraSearchDepth, noUkkonen, noOrderedEvaluation, noTruncation, ignoreAlignmentAdjustmentsForOm, maxSecondaryAlignmentsPerContig, &lv, &reverseLV, NULL, allocator);
	for (unsigned i = 0; i < maxBarcodeSize; i++)
		underlyingTenXSingleAligner[i]->setLandauVishkin(&lv, &reverseLV);

	singleSecondary[0] = singleSecondary[1] = NULL;
}

size_t
TenXClusterAligner::getBigAllocatorReservation(
	GenomeIndex *   index,
	unsigned        maxReadSize,
	unsigned        maxHits,
	unsigned        seedLen,
	unsigned        maxSeedsFromCommandLine,
	double          seedCoverage,
	unsigned        maxEditDistanceToConsider,
	unsigned        maxExtraSearchDepth,
	unsigned        maxCandidatePoolSize,
	int             maxSecondaryAlignmentsPerContig)
{
	return BaseAligner::getBigAllocatorReservation(index, false, maxHits, maxReadSize, seedLen, maxSeedsFromCommandLine, seedCoverage, maxSecondaryAlignmentsPerContig, maxExtraSearchDepth) + sizeof(TenXClusterAligner) + sizeof(_uint64);
}


TenXClusterAligner::~TenXClusterAligner()
{
	singleAligner->~BaseAligner();
}

#ifdef _DEBUG
extern bool _DumpAlignments;
#endif // _DEBUG


bool TenXClusterAligner::align(
	Read					**pairedReads,
	unsigned				clusterSize,
	PairedAlignmentResult	**result,
	int						maxEditDistanceForSecondaryResults,
	_int64					*secondaryResultBufferSize,
	_int64					*nSecondaryResults,
	_int64					*singleSecondaryBufferSize,
	_int64					maxSecondaryAlignmentsToReturn,
	_int64					*nSingleEndSecondaryResults,
	SingleAlignmentResult	**singleEndSecondaryResults,    // Single-end secondary alignments for when the paired-end alignment didn't work properly
	bool					*notFinished
)
{
	bool barcodeFinished = true;
	for (int readIdx = 0; readIdx < clusterSize; readIdx++) {
		if (notFinished[readIdx]) {
			result[readIdx]->status[0] = result[readIdx]->status[1] = NotFound;
			nSecondaryResults[readIdx] = 0;

			_int64 *nSingleEndSecondaryResultsForFirstRead = nSingleEndSecondaryResults + NUM_READS_PER_PAIR * readIdx;
			_int64 *nSingleEndSecondaryResultsForSecondRead = nSingleEndSecondaryResults + NUM_READS_PER_PAIR * readIdx + 1;

			*nSingleEndSecondaryResultsForFirstRead = 0;
			*nSingleEndSecondaryResultsForSecondRead = 0;

			Read *read0 = pairedReads[readIdx * NUM_READS_PER_PAIR];
			Read *read1 = pairedReads[readIdx * NUM_READS_PER_PAIR + 1];

			if (read0->getDataLength() < minReadLength && read1->getDataLength() < minReadLength) {
				TRACE("Reads are both too short -- returning");
				for (int whichRead = 0; whichRead < NUM_READS_PER_PAIR; whichRead++) {
					result[readIdx]->location[whichRead] = 0;
					result[readIdx]->mapq[whichRead] = 0;
					result[readIdx]->score[whichRead] = 0;
					result[readIdx]->status[whichRead] = NotFound;
				}
				result[readIdx]->alignedAsPair = false;
				result[readIdx]->fromAlignTogether = false;
				result[readIdx]->nanosInAlignTogether = 0;
				result[readIdx]->nLVCalls = 0;
				result[readIdx]->nSmallHits = 0;
				notFinished[readIdx] = false;
				continue;//return true;
			}

			_int64 start = timeInNanos();
			if (read0->getDataLength() >= minReadLength && read1->getDataLength() >= minReadLength) {
				//
				// Let the LVs use the cache that we built up.
				//
				bool fitInSecondaryBuffer =
					underlyingTenXSingleAligner[readIdx]->align(read0, read1, result[readIdx], maxEditDistanceForSecondaryResults, secondaryResultBufferSize[readIdx], &nSecondaryResults[readIdx], result[readIdx] + 1,
						maxSecondaryAlignmentsToReturn);

				if (!fitInSecondaryBuffer) {
					*nSingleEndSecondaryResultsForFirstRead = *nSingleEndSecondaryResultsForSecondRead = 0;
					nSecondaryResults[readIdx] = secondaryResultBufferSize[readIdx] + 1; // So the caller knows it's the paired secondary buffer that overflowed
					barcodeFinished = false;
					continue;//return false;
				}

				_int64 end = timeInNanos();

				result[readIdx]->nanosInAlignTogether = end - start;
				result[readIdx]->fromAlignTogether = true;
				result[readIdx]->alignedAsPair = true;

				if (forceSpacing) {
					if (result[readIdx]->status[0] == NotFound) {
						result[readIdx]->fromAlignTogether = false;
					}
					else {
						_ASSERT(result[readIdx]->status[1] != NotFound); // If one's not found, so is the other
					}
					notFinished[readIdx] = false;
					continue;//return true;
				}

				if (result[readIdx]->status[0] != NotFound && result[readIdx]->status[1] != NotFound) {
					//
					// Not a chimeric read.
					//
					notFinished[readIdx] = false;
					continue;//return true;
				}
			}

			/*
			//
			// If the intersecting aligner didn't find an alignment for these reads, then they may be
			// chimeric and so we should just align them with the single end aligner and apply a MAPQ penalty.
			//
			Read *read[NUM_READS_PER_PAIR] = { read0, read1 };
			_int64 *resultCount[2] = { nSingleEndSecondaryResultsForFirstRead, nSingleEndSecondaryResultsForSecondRead };

			for (int r = 0; r < NUM_READS_PER_PAIR; r++) {
				SingleAlignmentResult singleResult;
				_int64 singleEndSecondaryResultsThisTime = 0;

				if (read[r]->getDataLength() < minReadLength) {
					result->status[r] = NotFound;
					result->mapq[r] = 0;
					result->direction[r] = FORWARD;
					result->location[r] = 0;
					result->score[r] = 0;
				}
				else {
					// We're using *nSingleEndSecondaryResultsForFirstRead because it's either 0 or what all we've seen (i.e., we know NUM_READS_PER_PAIR is 2)
					bool fitInSecondaryBuffer = //true;
						singleAligner->AlignRead(read[r], &singleResult, maxEditDistanceForSecondaryResults,
							singleSecondaryBufferSize - *nSingleEndSecondaryResultsForFirstRead, &singleEndSecondaryResultsThisTime,
							maxSecondaryAlignmentsToReturn, singleEndSecondaryResults + *nSingleEndSecondaryResultsForFirstRead);

					if (!fitInSecondaryBuffer) {
						*nSecondaryResults = 0;
						*nSingleEndSecondaryResultsForFirstRead = singleSecondaryBufferSize + 1;
						*nSingleEndSecondaryResultsForSecondRead = 0;
						return false;
					}

					*(resultCount[r]) = singleEndSecondaryResultsThisTime;

					result->status[r] = singleResult.status;
					result->mapq[r] = singleResult.mapq / 3;   // Heavy quality penalty for chimeric reads
					result->direction[r] = singleResult.direction;
					result->location[r] = singleResult.location;
					result->score[r] = singleResult.score;
					result->scorePriorToClipping[r] = singleResult.scorePriorToClipping;
				}
			}

			result->fromAlignTogether = false;
			result->alignedAsPair = false;
			*/

#ifdef _DEBUG
			if (_DumpAlignments) {
				printf("TenXClusterAligner: (%u, %u) score (%d, %d), MAPQ (%d, %d)\n\n\n", result[readIdx]->location[0].location, result[readIdx]->location[1].location,
					result[readIdx]->score[0], result[readIdx]->score[1], result[readIdx]->mapq[0], result[readIdx]->mapq[1]);
			}
#endif // _DEBUG

			// This pair is done processing.
			notFinished[readIdx] = false;
		}
	}
	return barcodeFinished;
}
