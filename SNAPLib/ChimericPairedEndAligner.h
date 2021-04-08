/*++

Module Name:

    ChimericPairedEndAligner.h

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

#pragma once

#include "PairedEndAligner.h"
#include "BaseAligner.h"
#include "BigAlloc.h"

class ChimericPairedEndAligner : public PairedEndAligner {
public:
    ChimericPairedEndAligner(
        GenomeIndex         *index_,
        unsigned            maxReadSize,
        unsigned            maxHits,
        unsigned            maxK,
        unsigned            maxSeedsFromCommandLine,
        double              seedCoverage,
	    unsigned            minWeightToCheck,
        bool                forceSpacing_,
        unsigned            extraSearchDepth_,
        bool                noUkkonen,
        bool                noOrderedEvaluation,
		bool				noTruncation,
        bool                useAffineGap,
        bool                ignoreAlignmentAdjustmentsForOm,
		bool                altAwareness,
        bool                emitALTAlignments,
        PairedEndAligner    *underlyingPairedEndAligner_,
		unsigned			minReadLength_,
        int                 maxSecondaryAlignmentsPerContig,
        int                 maxScoreGapToPreferNonAltAlignment,
        int                 flattenMAPQAtOrBelow_,
        bool                useSoftClipping_,
        unsigned            matchReward = 1,
        unsigned            subPenalty = 4,
        unsigned            gapOpenPenalty = 6,
        unsigned            gapExtendPenalty = 1,
        int                 minScoreRealignment_ = 3,
        int                 minScoreGapRealignmentALT_ = 3,
        int                 minAGScoreImprovement_ = 24,
        bool                enableHammingScoringBaseAligner = false,
        BigAllocator        *allocator = NULL);
    
    virtual ~ChimericPairedEndAligner();
    
    static size_t getBigAllocatorReservation(GenomeIndex * index, unsigned maxReadSize, unsigned maxHits, unsigned seedLen, unsigned maxSeedsFromCommandLine, 
                                             double seedCoverage, unsigned maxEditDistanceToConsider, unsigned maxExtraSearchDepth, unsigned maxCandidatePoolSize,
                                             int maxSecondaryAlignmentsPerContig);

    void *operator new(size_t size, BigAllocator *allocator) {_ASSERT(size == sizeof(ChimericPairedEndAligner)); return allocator->allocate(size);}
    void operator delete(void *ptr, BigAllocator *allocator) {/* do nothing.  Memory gets cleaned up when the allocator is deleted.*/}

    virtual bool align(
        Read                  *read0,
        Read                  *read1,
        PairedAlignmentResult *result,
        PairedAlignmentResult *firatALTResult,
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
	);

    void *operator new(size_t size) {return BigAlloc(size);}
    void operator delete(void *ptr) {BigDealloc(ptr);}

    virtual _int64 getLocationsScored() const {
        return underlyingPairedEndAligner->getLocationsScored() + singleAligner->getLocationsScored();
    }

private:
   
    bool        forceSpacing;
    BaseAligner *singleAligner;
    PairedEndAligner *underlyingPairedEndAligner;

    // avoid allocation in aligner calls
    IdPairVector* singleSecondary[2];

    LandauVishkin<1> lv;
    LandauVishkin<-1> reverseLV;

    // AffineGap<1> ag;
    // AffineGap<-1> reverseAG;
    AffineGapVectorized<1> ag;
    AffineGapVectorized<-1> reverseAG;

	GenomeIndex *index;
	unsigned	minReadLength;

    bool emitALTAlignments;
    bool useSoftClipping;

    unsigned maxKSingleEnd;
    unsigned maxKPairedEnd;

    int extraSearchDepth;
    int minScoreRealignment;
    int minScoreGapRealignmentALT;
    int minAGScoreImprovement;
    int flattenMAPQAtOrBelow;

    bool enableHammingScoringBaseAligner;
    bool useAffineGap;

};
