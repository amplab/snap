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
#include "AlignerOptions.h"

using namespace std;

#ifdef TRACE_PAIRED_ALIGNER
#define TRACE printf
#else
#define TRACE(...) {}
#endif

ChimericPairedEndAligner::ChimericPairedEndAligner(
        GenomeIndex             *index_,
        unsigned                 maxReadSize,
        unsigned                 maxHits,
        unsigned                 maxK,
        unsigned                 maxSeedsFromCommandLine,
        double                   seedCoverage,
		unsigned                 minWeightToCheck,
        bool                     forceSpacing_,
        unsigned                 extraSearchDepth_,
        DisabledOptimizations    disabledOptimizations,  
        bool                     useAffineGap_,
        bool                     ignoreAlignmentAdjustmentsForOm,
		bool				     altAwareness,
        bool                     emitALTAlignments_,
        PairedEndAligner        *underlyingPairedEndAligner_,
	    unsigned                 minReadLength_,
        int                      maxSecondaryAlignmentsPerContig,
        int                      maxScoreGapToPreferNonAltAlignment,
        int                      flattenMAPQAtOrBelow_,
        bool                     useSoftClipping_,
        unsigned                 matchReward,
        unsigned                 subPenalty,
        unsigned                 gapOpenPenalty,
        unsigned                 gapExtendPenalty,
        unsigned                 fivePrimeEndBonus,
        unsigned                 threePrimeEndBonus,
        int                      minScoreRealignment_,
        int                      minScoreGapRealignmentALT_,
        int                      minAGScoreImprovement_,
        bool                     enableHammingScoringBaseAligner_,
        BigAllocator            *allocator)
		: underlyingPairedEndAligner(underlyingPairedEndAligner_), forceSpacing(forceSpacing_), index(index_), minReadLength(minReadLength_), emitALTAlignments(emitALTAlignments_), 
           maxKSingleEnd(maxK / 2), maxKPairedEnd(maxK), extraSearchDepth(extraSearchDepth_),
           minScoreRealignment(minScoreRealignment_), minScoreGapRealignmentALT(minScoreGapRealignmentALT_), minAGScoreImprovement(minAGScoreImprovement_), useSoftClipping(useSoftClipping_),
           flattenMAPQAtOrBelow(flattenMAPQAtOrBelow_), enableHammingScoringBaseAligner(enableHammingScoringBaseAligner_), useAffineGap(useAffineGap_)
{
    // Create single-end aligners.
    singleAligner = new (allocator) BaseAligner(index, maxHits, maxK / 2  /* allocate half to each end instead of letting it float like when they're aligned together */, maxReadSize,
                                                maxSeedsFromCommandLine, seedCoverage, minWeightToCheck, extraSearchDepth, 
                                                disabledOptimizations, useAffineGap, 
                                                ignoreAlignmentAdjustmentsForOm, altAwareness, emitALTAlignments, maxScoreGapToPreferNonAltAlignment,
                                                maxSecondaryAlignmentsPerContig, &lv, &reverseLV,
                                                matchReward, subPenalty, gapOpenPenalty, gapExtendPenalty, fivePrimeEndBonus, threePrimeEndBonus,
                                                NULL, allocator);
    
    underlyingPairedEndAligner->setLandauVishkin(&lv, &reverseLV);

    ag.init(matchReward, subPenalty, gapOpenPenalty, gapExtendPenalty, fivePrimeEndBonus, threePrimeEndBonus);
    reverseAG.init(matchReward, subPenalty, gapOpenPenalty, gapExtendPenalty, fivePrimeEndBonus, threePrimeEndBonus);

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
extern volatile bool _DumpAlignments;
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
        _int64                 maxPairedCandidatesForAffineGapBufferSize,
        _int64                *nPairedCandidatesForAffineGap,
        PairedAlignmentResult *pairedCandidatesForAffineGap,
        _int64                 maxSingleCandidatesForAffineGapBufferSize,
        _int64                *nSingleCandidatesForAffineGapFirstRead,
        _int64                *nSingleCandidatesForAffineGapSecondRead,
        SingleAlignmentResult *singleCandidatesForAffineGap,
        int                    maxK
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
    result->liftover[0] = result->liftover[1] = false;
    result->agForcedSingleAlignerCall = false;

    firstALTResult->status[0] = firstALTResult->status[1] = NotFound;

    *nPairedCandidatesForAffineGap = 0;
    *nSingleCandidatesForAffineGapFirstRead = 0;
    *nSingleCandidatesForAffineGapSecondRead = 0;

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
            singleEndSecondaryResults, maxPairedCandidatesForAffineGapBufferSize, nPairedCandidatesForAffineGap, pairedCandidatesForAffineGap, maxSingleCandidatesForAffineGapBufferSize,
            nSingleCandidatesForAffineGapFirstRead, nSingleCandidatesForAffineGapSecondRead, singleCandidatesForAffineGap, maxKPairedEnd);

        if (*nPairedCandidatesForAffineGap > maxPairedCandidatesForAffineGapBufferSize) {
            *nSecondaryResults = *nSingleEndSecondaryResultsForFirstRead = *nSingleEndSecondaryResultsForSecondRead = 0;
            *nPairedCandidatesForAffineGap = maxPairedCandidatesForAffineGapBufferSize + 1; // So the caller knows it's the paired LV candidate buffer that overflowed
            return false;
        }

        if (!fitInSecondaryBuffer) {
            *nSingleEndSecondaryResultsForFirstRead = *nSingleEndSecondaryResultsForSecondRead = 0;
            *nPairedCandidatesForAffineGap = 0;
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

        //
        // If the command line requests zeroing out low mapqs, do it now.  There's one other place we need to check below.
        //
        result->mapq[0] = result->mapq[0] <= flattenMAPQAtOrBelow ? 0 : result->mapq[0];
        result->mapq[1] = result->mapq[1] <= flattenMAPQAtOrBelow ? 0 : result->mapq[1];

        // If we have already seen a good ALT pair, don't separate them by running the single end aligner
        bool seenBetterAltResult = (firstALTResult->status[0] != NotFound)
                                   && (firstALTResult->status[1] != NotFound)
                                   && (sumPairScoreAlt <= sumPairScore - minScoreGapRealignmentALT);

        bool useAltLiftover = result->liftover[0] && result->liftover[1];

        if ((result->usedAffineGapScoring[0] || result->usedAffineGapScoring[1]) && maxScore >= minScoreRealignment && !seenBetterAltResult && !useAltLiftover) {
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
        scoreLimitLeft = sumPairScore;
        if (result->status[0] != NotFound && result->status[1] != NotFound) {
            result->agForcedSingleAlignerCall = true;   // Only set this if we wouldn't have done it anyway.
        }
    }

    //
    // If the intersecting aligner didn't find an alignment for these reads (or we're double checking because of affine gap alignments), then they may be
    // chimeric and so we should just align them with the single end aligner and apply a MAPQ penalty.
    //
    Read *read[NUM_READS_PER_PAIR] = {read0, read1};
    _int64 *resultCount[2] = {nSingleEndSecondaryResultsForFirstRead, nSingleEndSecondaryResultsForSecondRead};
    _int64 *affineCandidates[2] = { nSingleCandidatesForAffineGapFirstRead, nSingleCandidatesForAffineGapSecondRead };

    SingleAlignmentResult singleResult[NUM_READS_PER_PAIR];
    SingleAlignmentResult firstSingleALTResult[NUM_READS_PER_PAIR];
    int singleEndAGScore = 0;
    bool chooseSingleEndMapq = true;
    for (int r = 0; r < NUM_READS_PER_PAIR; r++) {
        _int64 singleEndSecondaryResultsThisTime = 0;
        _int64 singleEndAffineCandidatesThisTime = 0;

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
            chooseSingleEndMapq = false;
        } else {
            if (compareWithSingleEndAlignment) {
                // Single-end alignments are not good enough to be considered
                if (scoreLimitLeft < 0) {
                    break;
                }
                int maxKForRead = __min((int)maxKSingleEnd, __min(result->score[r], scoreLimitLeft));
                singleAligner->setMaxK(maxKForRead);
            }

            // We're using *nSingleEndSecondaryResultsForFirstRead because it's either 0 or what all we've seen (i.e., we know NUM_READS_PER_PAIR is 2)
            bool fitInSecondaryBuffer =
                singleAligner->AlignRead(read[r], &singleResult[r], &firstSingleALTResult[r], maxEditDistanceForSecondaryResults,
                    singleSecondaryBufferSize - *nSingleEndSecondaryResultsForFirstRead, &singleEndSecondaryResultsThisTime,
                    maxSecondaryAlignmentsToReturn, singleEndSecondaryResults + *nSingleEndSecondaryResultsForFirstRead,
                    maxSingleCandidatesForAffineGapBufferSize - *nSingleCandidatesForAffineGapFirstRead, &singleEndAffineCandidatesThisTime,
                    singleCandidatesForAffineGap + *nSingleCandidatesForAffineGapFirstRead);

            if (singleEndAffineCandidatesThisTime > (maxSingleCandidatesForAffineGapBufferSize - *nSingleCandidatesForAffineGapFirstRead)) {
                *nSecondaryResults = *nSingleEndSecondaryResultsForFirstRead = *nSingleEndSecondaryResultsForSecondRead = 0;
                *nSingleCandidatesForAffineGapFirstRead = maxSingleCandidatesForAffineGapBufferSize + 1; // So the caller knows it's the single end candidate buffer that overflowed
                return false;
            }

            if (!fitInSecondaryBuffer) {
                *nSecondaryResults = 0;
                *nSingleEndSecondaryResultsForFirstRead = singleSecondaryBufferSize + 1;
                *nSingleEndSecondaryResultsForSecondRead = 0;
                *nSingleCandidatesForAffineGapFirstRead = 0;
                return false;
            }

            bool usedHammingScoreingBaseAligner = false;
            if (useSoftClipping && enableHammingScoringBaseAligner) {
                //
                // Try Hamming distance scoring to align any unmapped reads
                //
                if (singleResult[r].status == NotFound && result->status[r] == NotFound) {
                    
                    usedHammingScoreingBaseAligner = true;

                    singleAligner->AlignRead(read[r], &singleResult[r], &firstSingleALTResult[r], maxEditDistanceForSecondaryResults,
                        singleSecondaryBufferSize - *nSingleEndSecondaryResultsForFirstRead, &singleEndSecondaryResultsThisTime,
                        maxSecondaryAlignmentsToReturn, singleEndSecondaryResults + *nSingleEndSecondaryResultsForFirstRead,
                        maxSingleCandidatesForAffineGapBufferSize - *nSingleCandidatesForAffineGapFirstRead, &singleEndAffineCandidatesThisTime,
                        singleCandidatesForAffineGap + *nSingleCandidatesForAffineGapFirstRead, true);

                    if (singleEndAffineCandidatesThisTime > (maxSingleCandidatesForAffineGapBufferSize - *nSingleCandidatesForAffineGapFirstRead)) {
                        *nSecondaryResults = *nSingleEndSecondaryResultsForFirstRead = *nSingleEndSecondaryResultsForSecondRead = 0;
                        *nSingleCandidatesForAffineGapFirstRead = maxSingleCandidatesForAffineGapBufferSize + 1; // So the caller knows it's the single end candidate buffer that overflowed
                        return false;
                    }

                    if (!fitInSecondaryBuffer) {
                        *nSecondaryResults = 0;
                        *nSingleEndSecondaryResultsForFirstRead = singleSecondaryBufferSize + 1;
                        *nSingleEndSecondaryResultsForSecondRead = 0;
                        *nSingleCandidatesForAffineGapFirstRead = 0;
                        return false;
                    }

                    _ASSERT(useAffineGap);
                    singleAligner->alignAffineGap(read[r], &singleResult[r], &firstSingleALTResult[r], 
                        singleEndAffineCandidatesThisTime, singleCandidatesForAffineGap + *nSingleCandidatesForAffineGapFirstRead);
                }
            }

            *(resultCount[r]) = singleEndSecondaryResultsThisTime;
            *(affineCandidates[r]) = singleEndAffineCandidatesThisTime;

            if (compareWithSingleEndAlignment) {
                //
                // Do not lower scoreLimit for mate if we used Hamming scoring for the base aligner
                //
                if (!usedHammingScoreingBaseAligner) {
                    if (singleResult[r].score != ScoreAboveLimit && singleResult[r].score != BaseAligner::UnusedScoreValue) {
                        scoreLimitLeft -= singleResult[r].score;
                    } else {
                        scoreLimitLeft = ScoreAboveLimit;
                    }
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

            //
            // If the command line requests zeroing out low mapqs, do it now.  There's one other place we need to check below.
            //  
            if (result->mapq[r] <= flattenMAPQAtOrBelow) {
                result->mapq[r] = 0;
            } // if we flatten
        } // for each read
    } // choose single end mapq


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
                result->mapq[r] = result->mapq[r] <= 3 ? 0 : result->mapq[r];
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
