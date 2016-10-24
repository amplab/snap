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

#pragma once

#include "TenXSingleAligner.h"
#include "BaseAligner.h"
#include "BigAlloc.h"

class TenXClusterAligner : public PairedEndAligner {
public:
    TenXClusterAligner(
        GenomeIndex            *index_,
        unsigned            maxReadSize,
        unsigned            maxHits,
        unsigned            maxK,
        unsigned            maxSeedsFromCommandLine,
        double                seedCoverage,
        unsigned            minWeightToCheck,
        bool                forceSpacing_,
        unsigned            extraSearchDepth,
        bool                noUkkonen,
        bool                noOrderedEvaluation,
        bool                noTruncation,
        bool                ignoreALignmentAdjustmentsForOm,
        TenXProgressTracker    *progressTracker_,
        bool                *clusterIsValid_,
        unsigned            *mappedPairsPerCluster_,
        unsigned            maxBarcodeSize_,
        unsigned            minPairsPerCluster_,
        _uint64                minClusterSpan_,
        double                unclusteredPenalty_,
        unsigned            clusterEDCompensation_,
        unsigned            minReadLength_,
        int                    maxSecondaryAlignmentsPerContig,
        unsigned            printStatsMapQLimit,
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
        unsigned                barcodeSize_
    );

    // Second stage will call underlyingAligner->phase3 but will not generate results.
    void align_second_stage_clustering (
    );

    // Second stage will call underlyingAligner->phase3 but not generate results. First stage should be only called once
    // Return ture if no pair requires memory reallocation.
    bool align_second_stage_check_reallocate (
    int                        maxEditDistanceForSecondaryResults
    );
    
    // Third stage will generate results and call phase4.
    // Return ture if no pair requires memory reallocation.
    void align_second_stage_generate_results (
    int                        maxEditDistanceForSecondaryResults,
    _int64                    maxSecondaryAlignmentsToReturn
    );


    // Return true if no more cluster change from true to false.
    bool checkClusterStabilized();
    void clusterResultCleanUp();

    // Third stage will clean up cluster mapping results.
    void align_thrid_stage(
    int maxEditDistanceForSecondaryResults,
    _int64 maxSecondaryAlignmentsToReturn
    );

    // Forth stage will handle single mappings.
    bool align_forth_stage(
    int maxEditDistanceForSecondaryResults,
    _int64 maxSecondaryAlignmentsToReturn
    );

    // TenX function. Move all trackers between start and end pointers beyond clusterBoundary, while registering each location with clusterID: clusterIdx.
    // All temporary changed will be held in updateHolder and will need to be merged with mergeUpdate().
    void registerClusterForReads(TenXProgressTracker *preStart, TenXProgressTracker *start, TenXProgressTracker *end, GenomeLocation clusterBoundary, int clusterIdx);
    // TenX function. Sort all progress trackers and link them and initialize root.
    void sortAndLink();
    // For debug purpose, quickly forward all pairs beyond the forward loc.
    void fastForward(GenomeLocation fowardLoc);

    // return true if all pairs within the cluster have been processed. No secondary results overflow.
    bool align(
        Read                    **pairedReads,
        unsigned                barcodeSize,
        PairedAlignmentResult    **result,
        int                        maxEditDistanceForSecondaryResults,
        _int64                    *secondaryResultBufferSize,
        _int64                    *nSecondaryResults,
        _int64                    *singleSecondaryBufferSize,
        _int64                    maxSecondaryAlignmentsToReturn,
        _int64                    *nSingleEndSecondaryResults,
        SingleAlignmentResult    **singleEndSecondaryResults,    // Single-end secondary alignments for when the paired-end alignment didn't work properly
        unsigned                *popularSeedsSkipped
    );

    /*
     * this align is just a place holder, to comply with the pure virtual function in PairedEndAligner
     */
    virtual bool align(
        Read                    *read0,
        Read                    *read1,
        PairedAlignmentResult    *result,
        int                        maxEditDistanceForSecondaryResults,
        _int64                    secondaryResultBufferSize,
        _int64                    *nSecondaryResults,
        PairedAlignmentResult    *secondaryResults,                // The caller passes in a buffer of secondaryResultBufferSize and it's filled in by AlignRead()
        _int64                    singleSecondaryBufferSize,
        _int64                    maxSecondaryAlignmentsToReturn,
        _int64                    *nSingleEndSecondaryResultsForFirstRead,
        _int64                    *nSingleEndSecondaryResultsForSecondRead,
        SingleAlignmentResult    *singleEndSecondaryResults        // Single-end secondary alignments for when the paired-end alignment didn't work properly
    ) {
        return true;
    };

    void *operator new(size_t size) { return BigAlloc(size); }
    void operator delete(void *ptr) { BigDealloc(ptr); }

    virtual _int64 getLocationsScored() const {
        _int64 locationsScored = 0;
        for (int readIdx = 0; readIdx < barcodeSize; readIdx++) {
            locationsScored += progressTracker[readIdx].aligner->getLocationsScored();
        }
        return locationsScored + singleAligner->getLocationsScored();
    }

private:
    // Push newUpdate to the updateList at the correct ordered position. 
    void pushToUpdate(TenXProgressTracker *newUpdate);
    // Merge trackers from rootList and updateList in a sorted fashion.
    void mergeUpdate();

    bool                forceSpacing;
    BaseAligner         *singleAligner;
    unsigned            maxBarcodeSize;
    TenXProgressTracker *progressTracker;

    // avoid allocation in aligner calls
    IdPairVector*       singleSecondary[2];

    LandauVishkin<1>    lv;
    LandauVishkin<-1>   reverseLV;

    GenomeIndex         *index;
    unsigned            minReadLength;
    
    // 10x data
    unsigned            barcodeSize;
    _uint8              minPairsPerCluster;
    _uint64             minClusterSpan;
    double              unclusteredPenalty;
    unsigned            clusterEDCompensation;

    // 10x update pointers
    TenXProgressTracker *trackerRoot;
    TenXProgressTracker *updateHolder;

    // 10x cluster validation tracker
    bool                *clusterIsValid;
    unsigned            *mappedPairsPerCluster;
    int                 maxClusterIdx;

};
