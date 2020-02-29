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
        bool                useAffineGap,
        bool                ignoreAlignmentAdjustmentsForOm,
		bool				altAwareness,
        PairedEndAligner    *underlyingPairedEndAligner_,
	    unsigned            minReadLength_,
        int                 maxSecondaryAlignmentsPerContig,
        int                 maxScoreGapToPreferNonAltAlignment,
        unsigned            matchReward,
        unsigned            subPenalty,
        unsigned            gapOpenPenalty,
        unsigned            gapExtendPenalty,
        unsigned            minAGScore,
        BigAllocator        *allocator)
		: underlyingPairedEndAligner(underlyingPairedEndAligner_), forceSpacing(forceSpacing_), index(index_), minReadLength(minReadLength_)
{
    // Create single-end aligners.
    singleAligner = new (allocator) BaseAligner(index, maxHits, maxK, maxReadSize,
                                                maxSeedsFromCommandLine, seedCoverage, minWeightToCheck, extraSearchDepth, 
                                                noUkkonen, noOrderedEvaluation, noTruncation, useAffineGap, 
                                                ignoreAlignmentAdjustmentsForOm, altAwareness, maxSecondaryAlignmentsPerContig, maxScoreGapToPreferNonAltAlignment , &lv, &reverseLV,
                                                matchReward, subPenalty, gapOpenPenalty, gapExtendPenalty, minAGScore,
                                                NULL, allocator);
    
    underlyingPairedEndAligner->setLandauVishkin(&lv, &reverseLV);

    ag.init(matchReward, subPenalty, gapOpenPenalty, gapExtendPenalty);
    reverseAG.init(matchReward, subPenalty, gapOpenPenalty, gapExtendPenalty);

    underlyingPairedEndAligner->setAffineGap(&ag, &reverseAG);

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
        PairedAlignmentResult *firstALTResult,
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
    result->usedAffineGapScoring[0] = result->usedAffineGapScoring[1] = false;
    result->basesClippedBefore[0] = result->basesClippedBefore[1] = 0;
    result->basesClippedAfter[0] = result->basesClippedAfter[1] = 0;
    result->agScore[0] = result->agScore[1] = 0;

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
    int pairAGScore = 0;
    bool compareWithSingleEndAlignment = false;
	if (read0->getDataLength() >= minReadLength && read1->getDataLength() >= minReadLength) {
		//
		// Let the LVs use the cache that we built up.
		//
        bool fitInSecondaryBuffer = 
		underlyingPairedEndAligner->align(read0, read1, result, firstALTResult, maxEditDistanceForSecondaryResults, secondaryResultBufferSize, nSecondaryResults, secondaryResults,
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
        
        int maxScore = __max(result->score[0], result->score[1]);

        if ((result->usedAffineGapScoring[0] || result->usedAffineGapScoring[1]) && maxScore >= 3) { // FIXME: Replace 3 with parameter
            compareWithSingleEndAlignment = true;
        }
        
		if (result->status[0] != NotFound && result->status[1] != NotFound && (!compareWithSingleEndAlignment)) {
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
    Read *read[NUM_READS_PER_PAIR] = {read0, read1};
    _int64 *resultCount[2] = {nSingleEndSecondaryResultsForFirstRead, nSingleEndSecondaryResultsForSecondRead};

    SingleAlignmentResult singleResult[NUM_READS_PER_PAIR];
    SingleAlignmentResult firstSingleALTResult[NUM_READS_PER_PAIR];
    int singleEndAGScore = 0;
    for (int r = 0; r < NUM_READS_PER_PAIR; r++) {
        _int64 singleEndSecondaryResultsThisTime = 0;

        if (compareWithSingleEndAlignment) {
            pairAGScore += result->agScore[r];
        }

        if (read[r]->getDataLength() < minReadLength) {
            result->status[r] = NotFound;
            result->mapq[r] = 0;
            result->direction[r] = FORWARD;
            result->location[r] = 0;
            result->score[r] = 0;
            result->usedAffineGapScoring[r] = false;
            result->basesClippedBefore[r] = 0;
            result->basesClippedAfter[r] = 0;
            result->agScore[r] = 0;
            result->fromAlignTogether = false;
            result->alignedAsPair = false;
        } else {
            // We're using *nSingleEndSecondaryResultsForFirstRead because it's either 0 or what all we've seen (i.e., we know NUM_READS_PER_PAIR is 2)
            bool fitInSecondaryBuffer =
                singleAligner->AlignRead(read[r], &singleResult[r], &firstSingleALTResult[r], maxEditDistanceForSecondaryResults,
                    singleSecondaryBufferSize - *nSingleEndSecondaryResultsForFirstRead, &singleEndSecondaryResultsThisTime,
                    maxSecondaryAlignmentsToReturn, singleEndSecondaryResults + *nSingleEndSecondaryResultsForFirstRead);

            if (!fitInSecondaryBuffer) {
                *nSecondaryResults = 0;
                *nSingleEndSecondaryResultsForFirstRead = singleSecondaryBufferSize + 1;
                *nSingleEndSecondaryResultsForSecondRead = 0;
                return false;
            }

            *(resultCount[r]) = singleEndSecondaryResultsThisTime;

            if (compareWithSingleEndAlignment) {
                singleEndAGScore += singleResult[r].agScore;
            }
        } // Not too short
    } // For each read in the pair

    if (!compareWithSingleEndAlignment || (singleEndAGScore >= pairAGScore + 12)) { // FIXME: Add threshold for choosing single-end alignments
        for (int r = 0; r < NUM_READS_PER_PAIR; r++) {
            if (read[r]->getDataLength() < minReadLength) {
                result->status[r] = NotFound;
                result->mapq[r] = 0;
                result->direction[r] = FORWARD;
                result->location[r] = 0;
                result->score[r] = 0;
                result->usedAffineGapScoring[r] = false;
                result->basesClippedBefore[r] = 0;
                result->basesClippedAfter[r] = 0;
                result->agScore[r] = 0;
            }
            else {
                result->status[r] = singleResult[r].status;
                result->mapq[r] = singleResult[r].mapq / 3;   // Heavy quality penalty for chimeric reads
                result->direction[r] = singleResult[r].direction;
                result->location[r] = singleResult[r].location;
                result->score[r] = singleResult[r].score;
                result->scorePriorToClipping[r] = singleResult[r].scorePriorToClipping;
                result->usedAffineGapScoring[r] = singleResult[r].usedAffineGapScoring;
                result->basesClippedBefore[r] = singleResult[r].basesClippedBefore;
                result->basesClippedAfter[r] = singleResult[r].basesClippedAfter;
                result->agScore[r] = singleResult[r].agScore;
                _ASSERT(result->basesClippedAfter[r] >= 0);
                _ASSERT(result->basesClippedBefore[r] >= 0);
            }
        }
        result->fromAlignTogether = false;
        result->alignedAsPair = false;
    }

#ifdef _DEBUG
    if (_DumpAlignments) {
        printf("ChimericPairedEndAligner: (%llu, %llu) score (%d, %d), MAPQ (%d, %d)\n\n\n",result->location[0].location, result->location[1].location,
            result->score[0], result->score[1], result->mapq[0], result->mapq[1]);
    }
#endif // _DEBUG

    return true;                    
}
