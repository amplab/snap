/*++

Module Name:

    TenXClusterAligner.h

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

#ifndef __TENX_CLUSTER_ALIGNER__
#define __TENX_CLUSTER_ALIGNER__

#pragma once

#include "TenXSingleAligner.h"
#include "BaseAligner.h"
#include "BigAlloc.h"

// declare trackers first
struct TenXMultiTracker
{
    // Debugging Field
    //std::string* id[NUM_DIRECTIONS];
    // Debugging Field

    bool                  pairNotDone;
    bool                  singleNotDone;
    TenXSingleAligner     *aligner;
    GenomeLocation        nextLocus; // Keep it here so that hopefully nextLocus be in cache.
    GenomeLocation        currentLocus; // Keep it here so that hopefully currentLocus be in cache.
    TenXMultiTracker      *nextTracker; // linked list next link

    Read                  *pairedReads[NUM_READS_PER_PAIR];
    PairedAlignmentResult *results;
    SingleAlignmentResult *singleEndSecondaryResults;
    _int64                secondaryResultBufferSize;
    _int64                singleSecondaryBufferSize;
    _int64                nSecondaryResults;
    _int64                nSingleEndSecondaryResults[NUM_READS_PER_PAIR];
    unsigned              popularSeedsSkipped[NUM_READS_PER_PAIR];
    bool                  useful[NUM_READS_PER_PAIR];
    int                   bestPairScore;
    GenomeLocation        bestResultGenomeLocation[NUM_READS_PER_PAIR];
    Direction             bestResultDirection[NUM_READS_PER_PAIR];
    double                probabilityOfAllPairs;
    unsigned              bestResultScore[NUM_READS_PER_PAIR];
    double                probabilityOfBestPair;
    int                   bestClusterIdx;
    
    _uint8                *clusterCounterAry;
    bool                  *clusterToggle;
    
    bool                  updatedBestScore;

    static int compare(const void *a_raw, const void *b_raw)
    {
        TenXMultiTracker* a = (TenXMultiTracker*)a_raw;
        TenXMultiTracker* b = (TenXMultiTracker*)b_raw;
        if (!a->pairNotDone && !b->pairNotDone)
            return 0;
        else if (!a->pairNotDone)
            return 1;
        else if (!b->pairNotDone)
            return -1;
        else if (a->nextLocus == b->nextLocus)
            return 0;
        else if (a->nextLocus < b->nextLocus)
            return 1;
        else
            return -1;
    }

};

struct TenXAnchorTracker {
    Read                  *pairedReads[NUM_READS_PER_PAIR];
	PairedAlignmentResult result;

    static int compare(const void *first_, const void *second_)
	{
		const TenXAnchorTracker *first = (TenXAnchorTracker *)first_;
		const TenXAnchorTracker *second = (TenXAnchorTracker *)second_;

		if (first->result.location[0] < second->result.location[0]) {
			return 1;
		} else if (first->result.location[0] > second->result.location[0]) {
			return -1;
		} else {
			return 0;
		}
	}

};


class TenXClusterAligner : public PairedEndAligner {
public:
    TenXClusterAligner(
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
        bool                noTruncation,
        bool                ignoreALignmentAdjustmentsForOm,
		TenXAnchorTracker   *anchorTracker_,
        TenXMultiTracker    *multiTracker_,
        unsigned            minPairsPerCluster_,
        _uint64             minClusterSpan_,
        double              unclusteredPenalty_,
        unsigned            clusterEDCompensation_,
        unsigned            minReadLength_,
        int                 maxSecondaryAlignmentsPerContig,
        unsigned            printStatsMapQLimit,
        int                 maxEditDistanceForSecondaryResults,
        _int64              maxSecondaryAlignmentsToReturn,
        BigAllocator        *allocator);

    virtual ~TenXClusterAligner();

    static size_t getBigAllocatorReservation(GenomeIndex * index, unsigned maxReadSize, unsigned maxHits, unsigned seedLen, unsigned maxSeedsFromCommandLine,
        double seedCoverage, unsigned maxEditDistanceToConsider, unsigned maxExtraSearchDepth, unsigned maxCandidatePoolSize,
        int maxSecondaryAlignmentsPerContig);

    void *operator new(size_t size, BigAllocator *allocator) { _ASSERT(size == sizeof(TenXClusterAligner)); return allocator->allocate(size); }
    void operator delete(void *ptr, BigAllocator *allocator) {/* do nothing.  Memory gets cleaned up when the allocator is deleted.*/ }

    // First stage will call underlyingAligner->phase1 and phase2. First stage should be only called once
    // Return true is no single read is worthy of further examination
    bool align_first_stage(
        unsigned                anchorNum_,
        unsigned                multiPairNum_
    );

    // Second stage will call underlyingAligner->phase3 but will not generate results.
    void align_second_stage_clustering (
    );

    // Second stage will call underlyingAligner->phase3 but not generate results. First stage should be only called once
    // Return ture if no pair requires memory reallocation.
    bool align_second_stage_check_reallocate (
    );
    
    // Third stage will generate results and call phase4.
    // Return ture if no pair requires memory reallocation.
    void align_second_stage_generate_results (
    );

    // Third stage will clean up cluster mapping results.
    void align_third_stage();

    // Forth stage will handle single mappings.
    bool align_forth_stage();

    // TenX function. Move all trackers between start and end pointers beyond clusterBoundary, while registering each location with clusterID: clusterIdx.
    // All temporary changed will be held in updateHolder and will need to be merged with mergeUpdate().
    void registerClusterForReads(TenXMultiTracker *preStart, TenXMultiTracker *start, TenXMultiTracker *end, GenomeLocation clusterBoundary, int clusterIdx);
    // TenX function. Sort all progress trackers and link them and initialize root.
    void sortAndLink();
    // For debug purpose, quickly forward all pairs beyond the forward loc.
    void fastForward(GenomeLocation fowardLoc);

    // return true if all pairs within the cluster have been processed. No secondary results overflow.
    bool align(
        Read                    **pairedReads,
        unsigned                multiPairNum,
        PairedAlignmentResult   **result,
        int                     maxEditDistanceForSecondaryResults,
        _int64                  *secondaryResultBufferSize,
        _int64                  *nSecondaryResults,
        _int64                  *singleSecondaryBufferSize,
        _int64                  maxSecondaryAlignmentsToReturn,
        _int64                  *nSingleEndSecondaryResults,
        SingleAlignmentResult   **singleEndSecondaryResults,    // Single-end secondary alignments for when the paired-end alignment didn't work properly
        unsigned                *popularSeedsSkipped
    );

    /*
     * this align is just a place holder, to comply with the pure virtual function in PairedEndAligner
     */
    virtual bool align(
        Read                    *read0,
        Read                    *read1,
        PairedAlignmentResult   *result,
        int                     maxEditDistanceForSecondaryResults,
        _int64                  secondaryResultBufferSize,
        _int64                  *nSecondaryResults,
        PairedAlignmentResult   *secondaryResults,                // The caller passes in a buffer of secondaryResultBufferSize and it's filled in by AlignRead()
        _int64                  singleSecondaryBufferSize,
        _int64                  maxSecondaryAlignmentsToReturn,
        _int64                  *nSingleEndSecondaryResultsForFirstRead,
        _int64                  *nSingleEndSecondaryResultsForSecondRead,
        SingleAlignmentResult   *singleEndSecondaryResults        // Single-end secondary alignments for when the paired-end alignment didn't work properly
    ) {
        return true;
    };

    void *operator new(size_t size) { return BigAlloc(size); }
    void operator delete(void *ptr) { BigDealloc(ptr); }

    virtual _int64 getLocationsScored() const {
        _int64 locationsScored = 0;
        for (int readIdx = 0; readIdx < multiPairNum; readIdx++) {
            locationsScored += multiTracker[readIdx].aligner->getLocationsScored();
        }
        return locationsScored + singleAligner->getLocationsScored();
    }

private:
    // Push newUpdate to the updateList at the correct ordered position. 
    void pushToUpdate(TenXMultiTracker *newUpdate);
    // Merge trackers from rootList and updateList in a sorted fashion.
    void mergeUpdate();
    // Find the anchor next that is alone, out of the cluster range
	unsigned endOfAnchorChain (unsigned anchorIdx);
	unsigned moveAnchorPassLocus (unsigned anchorIdx, const GenomeLocation& targetLocus);

    bool                forceSpacing;
    BaseAligner         *singleAligner;
    TenXAnchorTracker   *anchorTracker;
    TenXMultiTracker    *multiTracker;
    int                 maxEditDistanceForSecondaryResults;
    _int64              maxSecondaryAlignmentsToReturn;

    // avoid allocation in aligner calls
    IdPairVector*       singleSecondary[2];

    LandauVishkin<1>    lv;
    LandauVishkin<-1>   reverseLV;

    GenomeIndex         *index;
    unsigned            minReadLength;
    
    // 10x data
    unsigned            anchorNum;
    unsigned            multiPairNum;
    _uint8              minPairsPerCluster;
    _uint64             coverageScanRange;
    _uint64             magnetRange;
    double              unclusteredPenalty;
    unsigned            clusterEDCompensation;

    // 10x update pointers
    TenXMultiTracker *trackerRoot;
    TenXMultiTracker *updateHolder;

    // 10x cluster validation tracker
    int                 maxClusterIdx;

};

#endif //#define __TENX_CLUSTER_ALIGNER__