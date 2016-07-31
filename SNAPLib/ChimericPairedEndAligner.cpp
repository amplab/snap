/*++

Module Name:

    ChimericPairedEndAligner.cpp

Abstract:

    A paired-end aligner calls into a different paired-end aligner, and if
    it fails to find an alignment, aligns each of the reads singly.  This handles
    chimeric reads that would otherwise be unalignable.

Authors:

    Hongyi Xin and Bill Bolosky, June, 2016

Environment:

    User mode service.

Revision History:

--*/


#include "stdafx.h"
#include "ChimericPairedEndAligner.h"
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

ChimericPairedEndAligner::ChimericPairedEndAligner(
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
        PairedEndAligner    *underlyingPairedEndAligner_,
	    unsigned            minReadLength_,
        int                 maxSecondaryAlignmentsPerContig,
	    bool                noSingle_,
        BigAllocator        *allocator)
		: noSingle(noSingle_), underlyingPairedEndAligner(underlyingPairedEndAligner_), forceSpacing(forceSpacing_), index(index_), minReadLength(minReadLength_)
{
    // Create single-end aligners.
    singleAligner = new (allocator) BaseAligner(index, maxHits, maxK, maxReadSize,
        maxSeedsFromCommandLine, seedCoverage, minWeightToCheck, extraSearchDepth, noUkkonen, noOrderedEvaluation, noTruncation, ignoreAlignmentAdjustmentsForOm, maxSecondaryAlignmentsPerContig, &lv, &reverseLV, NULL, allocator);
    
    underlyingPairedEndAligner->setLandauVishkin(&lv, &reverseLV);

    singleSecondary[0] = singleSecondary[1] = NULL;
}

    size_t 
ChimericPairedEndAligner::getBigAllocatorReservation(
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
    return BaseAligner::getBigAllocatorReservation(index, false, maxHits, maxReadSize, seedLen, maxSeedsFromCommandLine, seedCoverage, maxSecondaryAlignmentsPerContig, maxExtraSearchDepth) + sizeof(ChimericPairedEndAligner)+sizeof(_uint64);
}


ChimericPairedEndAligner::~ChimericPairedEndAligner()
{
    singleAligner->~BaseAligner();
}

#ifdef _DEBUG
extern bool _DumpAlignments;
#endif // _DEBUG


bool ChimericPairedEndAligner::align(
        Read                  *read0,
        Read                  *read1,
        PairedAlignmentResult *result,
        int                    maxEditDistanceForSecondaryResults,
        _int64                 secondaryResultBufferSize,
        _int64                *nSecondaryResults,
        PairedAlignmentResult *secondaryResults,             // The caller passes in a buffer of secondaryResultBufferSize and it's filled in by AlignRead()
        _int64                 singleSecondaryBufferSize,
        _int64                 maxSecondaryAlignmentsToReturn,
        _int64                *nSingleEndSecondaryResultsForFirstRead,
        _int64                *nSingleEndSecondaryResultsForSecondRead,
        SingleAlignmentResult *singleEndSecondaryResults     // Single-end secondary alignments for when the paired-end alignment didn't work properly        
        )
{
	result->status[0] = result->status[1] = NotFound;
    *nSecondaryResults = 0;
    *nSingleEndSecondaryResultsForFirstRead = 0;
    *nSingleEndSecondaryResultsForSecondRead = 0;

	if (read0->getDataLength() < minReadLength && read1->getDataLength() < minReadLength) {
        TRACE("Reads are both too short -- returning");
		for (int whichRead = 0; whichRead < NUM_READS_PER_PAIR; whichRead++) {
			result->location[whichRead] = 0;
			result->mapq[whichRead] = 0;
			result->score[whichRead] = 0;
			result->status[whichRead] = NotFound;
		}
		result->alignedAsPair = false;
		result->fromAlignTogether = false;
		result->nanosInAlignTogether = 0;
		result->nLVCalls = 0;
		result->nSmallHits = 0;
		return true;
    }

    _int64 start = timeInNanos();
	if (read0->getDataLength() >= minReadLength && read1->getDataLength() >= minReadLength) {
		//
		// Let the LVs use the cache that we built up.
		//
        bool fitInSecondaryBuffer = 
		underlyingPairedEndAligner->align(read0, read1, result, maxEditDistanceForSecondaryResults, secondaryResultBufferSize, nSecondaryResults, secondaryResults,
            singleSecondaryBufferSize, maxSecondaryAlignmentsToReturn, nSingleEndSecondaryResultsForFirstRead, nSingleEndSecondaryResultsForSecondRead, 
            singleEndSecondaryResults);

        if (!fitInSecondaryBuffer) {
            *nSingleEndSecondaryResultsForFirstRead = *nSingleEndSecondaryResultsForSecondRead = 0;
            *nSecondaryResults = secondaryResultBufferSize + 1; // So the caller knows it's the paired secondary buffer that overflowed
            return false;
        }

		_int64 end = timeInNanos();

		result->nanosInAlignTogether = end - start;
		result->fromAlignTogether = true;
		result->alignedAsPair = true;

		if (forceSpacing) {
			if (result->status[0] == NotFound) {
				result->fromAlignTogether = false;
			}
			else {
				_ASSERT(result->status[1] != NotFound); // If one's not found, so is the other
			}
			return true;
		}

		if (result->status[0] != NotFound && result->status[1] != NotFound) {
			//
			// Not a chimeric read.
			//
			return true;
		}
	}

    //
    // If the intersecting aligner didn't find an alignment for these reads, then they may be
    // chimeric and so we should just align them with the single end aligner and apply a MAPQ penalty.
    //
	if (!noSingle) {
		Read *read[NUM_READS_PER_PAIR] = {read0, read1};
		_int64 *resultCount[2] = {nSingleEndSecondaryResultsForFirstRead, nSingleEndSecondaryResultsForSecondRead};

		for (int r = 0; r < NUM_READS_PER_PAIR; r++) {
			SingleAlignmentResult singleResult;
			_int64 singleEndSecondaryResultsThisTime = 0;

			if (read[r]->getDataLength() < minReadLength) {
				result->status[r] = NotFound;
				result->mapq[r] = 0;
				result->direction[r] = FORWARD;
				result->location[r] = 0;
				result->score[r] = 0;
			} else {
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
	}

#ifdef _DEBUG
    if (_DumpAlignments) {
        printf("ChimericPairedEndAligner: (%u, %u) score (%d, %d), MAPQ (%d, %d)\n\n\n",result->location[0].location, result->location[1].location,
            result->score[0], result->score[1], result->mapq[0], result->mapq[1]);
    }
#endif // _DEBUG

    return true;                    
}
