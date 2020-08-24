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
        bool                emitALTAlignments_,
        PairedEndAligner    *underlyingPairedEndAligner_,
	    unsigned            minReadLength_,
        int                 maxSecondaryAlignmentsPerContig,
        int                 maxScoreGapToPreferNonAltAlignment,
        unsigned            matchReward,
        unsigned            subPenalty,
        unsigned            gapOpenPenalty,
        unsigned            gapExtendPenalty,
        int                 minScoreRealignment_,
        int                 minScoreGapRealignmentALT_,
        int                 minAGScoreImprovement_,
        BigAllocator        *allocator)
		: underlyingPairedEndAligner(underlyingPairedEndAligner_), forceSpacing(forceSpacing_), index(index_), minReadLength(minReadLength_), emitALTAlignments(emitALTAlignments_), 
           maxKSingleEnd(maxK / 2),
           minScoreRealignment(minScoreRealignment_), minScoreGapRealignmentALT(minScoreGapRealignmentALT_), minAGScoreImprovement(minAGScoreImprovement_)
{
    // Create single-end aligners.
    singleAligner = new (allocator) BaseAligner(index, maxHits, maxK / 2  /* allocate half to each end instead of letting it float like when they're aligned together */, maxReadSize,
                                                maxSeedsFromCommandLine, seedCoverage, minWeightToCheck, extraSearchDepth, 
                                                noUkkonen, noOrderedEvaluation, noTruncation, useAffineGap, 
                                                ignoreAlignmentAdjustmentsForOm, altAwareness, emitALTAlignments, maxScoreGapToPreferNonAltAlignment,
                                                maxSecondaryAlignmentsPerContig, &lv, &reverseLV,
                                                matchReward, subPenalty, gapOpenPenalty, gapExtendPenalty,
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
        SingleAlignmentResult *singleEndSecondaryResults,     // Single-end secondary alignments for when the paired-end alignment didn't work properly
        _int64				  maxLVCandidatesForAffineGapBufferSize,
        _int64				  *nLVCandidatesForAffineGap,
        PairedAlignmentResult *lvCandidatesForAffineGap
	)
{
	result->status[0] = result->status[1] = NotFound;
    *nSecondaryResults = 0;
    *nSingleEndSecondaryResultsForFirstRead = 0;
    *nSingleEndSecondaryResultsForSecondRead = 0;
    result->usedAffineGapScoring[0] = result->usedAffineGapScoring[1] = false;
    result->basesClippedBefore[0] = result->basesClippedBefore[1] = 0;
    result->basesClippedAfter[0] = result->basesClippedAfter[1] = 0;
    result->clippingForReadAdjustment[0] = result->clippingForReadAdjustment[1] = 0;
    result->agScore[0] = result->agScore[1] = 0;

    firstALTResult->status[0] = firstALTResult->status[1] = NotFound;

    *nLVCandidatesForAffineGap = 0;

	if (read0->getDataLength() < minReadLength && read1->getDataLength() < minReadLength) {
        TRACE("Reads are both too short -- returning");
		for (int whichRead = 0; whichRead < NUM_READS_PER_PAIR; whichRead++) {
			result->location[whichRead] = InvalidGenomeLocation;
			result->mapq[whichRead] = 0;
			result->score[whichRead] = 0;
			result->status[whichRead] = NotFound;
		}
		result->alignedAsPair = false;
		result->nanosInAlignTogether = 0;
		result->nLVCalls = 0;
		result->nSmallHits = 0;
		return true;
    }

    _int64 start = timeInNanos();
    int pairAGScore = 0, sumPairScore = 0, sumPairScoreAlt = 0;
    bool compareWithSingleEndAlignment = false;
	if (read0->getDataLength() >= minReadLength && read1->getDataLength() >= minReadLength) {
		//
		// Let the LVs use the cache that we built up.
		//
        bool fitInSecondaryBuffer = 
		underlyingPairedEndAligner->align(read0, read1, result, firstALTResult, maxEditDistanceForSecondaryResults, secondaryResultBufferSize, nSecondaryResults, secondaryResults,
            singleSecondaryBufferSize, maxSecondaryAlignmentsToReturn, nSingleEndSecondaryResultsForFirstRead, nSingleEndSecondaryResultsForSecondRead, 
            singleEndSecondaryResults, maxLVCandidatesForAffineGapBufferSize, nLVCandidatesForAffineGap, lvCandidatesForAffineGap);

        if (*nLVCandidatesForAffineGap > maxLVCandidatesForAffineGapBufferSize) {
            *nSecondaryResults = *nSingleEndSecondaryResultsForFirstRead = *nSingleEndSecondaryResultsForSecondRead = 0;
            *nLVCandidatesForAffineGap = maxLVCandidatesForAffineGapBufferSize + 1; // So the caller knows it's the paired LV candidate buffer that overflowed
            return false;
        }

        if (!fitInSecondaryBuffer) {
            *nSingleEndSecondaryResultsForFirstRead = *nSingleEndSecondaryResultsForSecondRead = 0;
            *nLVCandidatesForAffineGap = 0;
            *nSecondaryResults = secondaryResultBufferSize + 1; // So the caller knows it's the paired secondary buffer that overflowed
            return false;
        }

		_int64 end = timeInNanos();

		result->nanosInAlignTogether = end - start;
		result->alignedAsPair = true;

		if (forceSpacing) {
			if (result->status[0] == NotFound) {
				result->alignedAsPair = false;
			} else {
				_ASSERT(result->status[1] != NotFound); // If one's not found, so is the other
			}
			return true;
		}
        
        int maxScore = __max(result->score[0], result->score[1]);
        sumPairScore = result->score[0] + result->score[1];
        sumPairScoreAlt = firstALTResult->score[0] + firstALTResult->score[1];

        // If we have already seen a good ALT pair, don't separate them with by running the single end aligner
        bool seenBetterAltResult = (firstALTResult->status[0] != NotFound)
                                   && (firstALTResult->status[1] != NotFound)
                                   && (sumPairScoreAlt <= sumPairScore - minScoreGapRealignmentALT);

        if ((result->usedAffineGapScoring[0] || result->usedAffineGapScoring[1]) && maxScore >= minScoreRealignment && !seenBetterAltResult) {
            compareWithSingleEndAlignment = true;
        }
        
		if (result->status[0] != NotFound && result->status[1] != NotFound && (!compareWithSingleEndAlignment)) {
			//
			// Not a chimeric read.
			//
			return true;
		}
	}

    int scoreLimitLeft = maxKSingleEnd;
    if (compareWithSingleEndAlignment) {
        scoreLimitLeft = sumPairScore - 3;
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
    bool chooseSingleEndMapq = true;
    for (int r = 0; r < NUM_READS_PER_PAIR; r++) {
        _int64 singleEndSecondaryResultsThisTime = 0;

        if (compareWithSingleEndAlignment) {
            pairAGScore += result->agScore[r];
        }

        // Reset max edit distance for single end aligner
        singleAligner->setMaxK(maxKSingleEnd);

        if (read[r]->getDataLength() < minReadLength) {
            result->status[r] = NotFound;
            result->mapq[r] = 0;
            result->direction[r] = FORWARD;
            result->location[r] = InvalidGenomeLocation;
            result->score[r] = 0;
            result->usedAffineGapScoring[r] = false;
            result->basesClippedBefore[r] = 0;
            result->basesClippedAfter[r] = 0;
            result->agScore[r] = 0;
            result->alignedAsPair = false;
            result->clippingForReadAdjustment[r] = 0;

            firstALTResult->status[r] = NotFound;   // Don't need to fill in the rest, this suppresses writing it
        } else {
            if (compareWithSingleEndAlignment) {
                // Single-end alignments are not good enough to be considered
                if (scoreLimitLeft < 0) {
                    break;
                }
                int maxKForRead = __min(maxKSingleEnd, __min(result->score[r], scoreLimitLeft));
                singleAligner->setMaxK(maxKForRead);
            }

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
                if (singleResult[r].score != ScoreAboveLimit) {
                    scoreLimitLeft -= singleResult[r].score;
                }
                else {
                    scoreLimitLeft = ScoreAboveLimit;
                }
                singleEndAGScore += singleResult[r].agScore;
                if (result->agScore[r] >= singleResult[r].agScore) {
                    chooseSingleEndMapq = false;
                }
            }
        } // Not too short
    } // For each read in the pair

    if (chooseSingleEndMapq) {
        for (int r = 0; r < NUM_READS_PER_PAIR; r++) {
            //
            // If the single-end aligner returns a lower MAPQ choose that for the result
            //
            result->mapq[r] = __min(result->mapq[r], singleResult[r].mapq);
        }
    }


    if (!compareWithSingleEndAlignment || (singleEndAGScore >= pairAGScore + minAGScoreImprovement)) {
        for (int r = 0; r < NUM_READS_PER_PAIR; r++) {
            if (read[r]->getDataLength() < minReadLength) {
                result->status[r] = NotFound;
                result->mapq[r] = 0;
                result->direction[r] = FORWARD;
                result->location[r] = InvalidGenomeLocation;
                result->score[r] = 0;
                result->usedAffineGapScoring[r] = false;
                result->basesClippedBefore[r] = 0;
                result->basesClippedAfter[r] = 0;
                result->agScore[r] = 0;
                result->clippingForReadAdjustment[r] = 0;

                firstALTResult->status[r] = NotFound;   // Don't need to fill in the rest, this suppresses writing it
            } else {
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
        result->alignedAsPair = false;
    }

#ifdef _DEBUG
    if (_DumpAlignments) {
        printf("ChimericPairedEndAligner: (%s:%llu, %s:%llu) score (%d, %d), MAPQ (%d, %d)\n\n\n",
            index->getGenome()->getContigAtLocation(result->location[0])->name, result->location[0] - index->getGenome()->getContigAtLocation(result->location[0])->beginningLocation,
            index->getGenome()->getContigAtLocation(result->location[1])->name, result->location[1] - index->getGenome()->getContigAtLocation(result->location[1])->beginningLocation,
            result->score[0], result->score[1], result->mapq[0], result->mapq[1]);
    }
#endif // _DEBUG

    return true;                    
}
