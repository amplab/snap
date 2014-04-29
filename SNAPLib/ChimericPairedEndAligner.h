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

  friend class PairedAlignerContext;

public:
    ChimericPairedEndAligner(
        GenomeIndex         *index,
        unsigned            maxReadSize,
        unsigned            maxHits,
        unsigned            maxK,
        unsigned            maxSeedsFromCommandLine,
        double              seedCoverage,
        bool                forceSpacing_,
        unsigned            extraSearchDepth,
        bool                noUkkonen,
        bool                noOrderedEvaluation,
        PairedEndAligner    *underlyingPairedEndAligner_,
        BigAllocator        *allocator);
    
    virtual ~ChimericPairedEndAligner();
    
    static unsigned getMaxSingleEndSecondaryResults(unsigned maxSeedsToUse, double maxSeedCoverage, unsigned maxReadSize, unsigned maxHits, unsigned seedLength)
    {
        return BaseAligner::getMaxSecondaryResults(maxSeedsToUse, maxSeedCoverage, maxReadSize, maxHits, seedLength) * 2;   // *2 for each end
    }

    static size_t getBigAllocatorReservation(GenomeIndex * index, unsigned maxReadSize, unsigned maxHits, unsigned seedLen, unsigned maxSeedsFromCommandLine, 
                                             double seedCoverage, unsigned maxEditDistanceToConsider, unsigned maxExtraSearchDepth, unsigned maxCandidatePoolSize);

    void *operator new(size_t size, BigAllocator *allocator) {_ASSERT(size == sizeof(ChimericPairedEndAligner)); return allocator->allocate(size);}
    void operator delete(void *ptr, BigAllocator *allocator) {/* do nothing.  Memory gets cleaned up when the allocator is deleted.*/}

    virtual void align(
        Read                  *read0,
        Read                  *read1,
        PairedAlignmentResult *result,
        int                    maxEditDistanceForSecondaryResults,
        int                    secondaryResultBufferSize,
        int                   *nSecondaryResults,
        PairedAlignmentResult *secondaryResults,             // The caller passes in a buffer of secondaryResultBufferSize and it's filled in by AlignRead()
        int                    singleSecondaryBufferSize,
        int                   *nSingleEndSecondaryResultsForFirstRead,
        int                   *nSingleEndSecondaryResultsForSecondRead,
        SingleAlignmentResult *singleEndSecondaryResults     // Single-end secondary alignments for when the paired-end alignment didn't work properly
        );

    void *operator new(size_t size) {return BigAlloc(size);}
    void operator delete(void *ptr) {BigDealloc(ptr);}

    virtual _int64 getLocationsScored() const {
        return underlyingPairedEndAligner->getLocationsScored() + singleAligner->getLocationsScored();
    }

protected:
   
    bool        forceSpacing;
    BaseAligner *singleAligner;
    PairedEndAligner *underlyingPairedEndAligner;

    // avoid allocation in aligner calls
    IdPairVector* singleSecondary[2];

    LandauVishkin<1> lv;
    LandauVishkin<-1> reverseLV;
};
