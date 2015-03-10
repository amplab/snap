/*++

Module Name:

    ChimericPairedEndAligner.cpp

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
       PairedEndAligner    *underlyingPairedEndAligner_,
	   unsigned				minReadLength_,
       int                  maxSecondaryAlignmentsPerContig,
        BigAllocator        *allocator)
		: underlyingPairedEndAligner(underlyingPairedEndAligner_), forceSpacing(forceSpacing_), index(index_), minReadLength(minReadLength_)
{
    // Create single-end aligners.
    singleAligner = new (allocator) BaseAligner(index, maxHits, maxK, maxReadSize,
        maxSeedsFromCommandLine, seedCoverage, minWeightToCheck, extraSearchDepth, noUkkonen, noOrderedEvaluation, noTruncation, maxSecondaryAlignmentsPerContig, &lv, &reverseLV, NULL, allocator);
    
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
    return BaseAligner::getBigAllocatorReservation(index, false, maxHits, maxReadSize, seedLen, maxSeedsFromCommandLine, seedCoverage, maxSecondaryAlignmentsPerContig) + sizeof(ChimericPairedEndAligner)+sizeof(_uint64);
}


ChimericPairedEndAligner::~ChimericPairedEndAligner()
{
    singleAligner->~BaseAligner();
}

#ifdef _DEBUG
extern bool _DumpAlignments;
#endif // _DEBUG


void ChimericPairedEndAligner::align(
        Read                  *read0,
        Read                  *read1,
        PairedAlignmentResult *result,
        int                    maxEditDistanceForSecondaryResults,
        int                    secondaryResultBufferSize,
        int                   *nSecondaryResults,
        PairedAlignmentResult *secondaryResults,             // The caller passes in a buffer of secondaryResultBufferSize and it's filled in by AlignRead()
        int                    singleSecondaryBufferSize,
        int                    maxSecondaryAlignmentsToReturn,
        int                   *nSingleEndSecondaryResultsForFirstRead,
        int                   *nSingleEndSecondaryResultsForSecondRead,
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
		return;
    }

    _int64 start = timeInNanos();
	if (read0->getDataLength() >= minReadLength && read1->getDataLength() >= minReadLength) {
		//
		// Let the LVs use the cache that we built up.
		//
		underlyingPairedEndAligner->align(read0, read1, result, maxEditDistanceForSecondaryResults, secondaryResultBufferSize, nSecondaryResults, secondaryResults,
            singleSecondaryBufferSize, maxSecondaryAlignmentsToReturn, nSingleEndSecondaryResultsForFirstRead, nSingleEndSecondaryResultsForSecondRead, 
            singleEndSecondaryResults);

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
			return;
		}

		if (result->status[0] != NotFound && result->status[1] != NotFound) {
			//
			// Not a chimeric read.
			//
			return;
		}
	}

    //
    // If the intersecting aligner didn't find an alignment for these reads, then they may be
    // chimeric and so we should just align them with the single end aligner and apply a MAPQ penalty.
    //
    Read *read[NUM_READS_PER_PAIR] = {read0, read1};
    int *resultCount[2] = {nSingleEndSecondaryResultsForFirstRead, nSingleEndSecondaryResultsForSecondRead};

    for (int r = 0; r < NUM_READS_PER_PAIR; r++) {
        SingleAlignmentResult singleResult;
        int singleEndSecondaryResultsThisTime = 0;

		if (read[r]->getDataLength() < minReadLength) {
			result->status[r] = NotFound;
			result->mapq[r] = 0;
			result->direction[r] = FORWARD;
			result->location[r] = 0;
			result->score[r] = 0;
		} else {
			// We're using *nSingleEndSecondaryResultsForFirstRead because it's either 0 or what all we've seen (i.e., we know NUM_READS_PER_PAIR is 2)
			singleAligner->AlignRead(read[r], &singleResult, maxEditDistanceForSecondaryResults,
				singleSecondaryBufferSize - *nSingleEndSecondaryResultsForFirstRead, &singleEndSecondaryResultsThisTime,
                maxSecondaryAlignmentsToReturn, singleEndSecondaryResults + *nSingleEndSecondaryResultsForFirstRead);

			*(resultCount[r]) = singleEndSecondaryResultsThisTime;

			result->status[r] = singleResult.status;
			result->mapq[r] = singleResult.mapq / 3;   // Heavy quality penalty for chimeric reads
			result->direction[r] = singleResult.direction;
			result->location[r] = singleResult.location;
			result->score[r] = singleResult.score;
		}
    }

    result->fromAlignTogether = false;
    result->alignedAsPair = false;

#ifdef _DEBUG
    if (_DumpAlignments) {
        printf("ChimericPairedEndAligner: (%u, %u) score (%d, %d), MAPQ (%d, %d)\n\n\n",result->location[0], result->location[1],
            result->score[0], result->score[1], result->mapq[0], result->mapq[1]);
    }
#endif // _DEBUG
                    
}
