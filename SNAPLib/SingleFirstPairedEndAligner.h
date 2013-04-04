/*++

Module Name:

    SingleFirstPairedEndAligner.h

Abstract:

    A paired-end aligner that uses a single-end aligner to try to find the best match
    for a pair.  If that doesn't work, it falls back to a different paired-end aligner.

Authors:

    Bill Bolosky, March, 2013

Environment:

    User mode service.

Revision History:

    Factored from Matei's SmarterPariedEndAligner

--*/

#pragma once

#include "PairedEndAligner.h"
#include "BaseAligner.h"

class SingleFirstPairedEndAligner {
public:
    SingleFirstPairedEndAligner(
        GenomeIndex         *index,
        unsigned            maxReadSize_,
        unsigned            confDiff_,
        unsigned            maxHits,
        unsigned            maxK_,
        unsigned            maxSeeds,
        unsigned            minSpacing_,                // Minimum distance to allow between the two ends.
        unsigned            maxSpacing_,                // Maximum distance to allow between the two ends.
        bool                forceSpacing_,
        unsigned            maxReadSize,
        unsigned            adaptiveConfDiffThreshold,  // Increase confDiff if this many seeds in the read have multiple hits.
        bool                skipAlignToegther_,
        PairedEndAligner    *underlyingPairedEndAligner_);
    
    virtual ~SingleFirstPairedEndAligner();
    
    virtual void align(
        Read                  *read0,
        Read                  *read1,
        PairedAlignmentResult *result);

    void *operator new(size_t size) {return BigAlloc(size);}
    void operator delete(void *ptr) {BigDealloc(ptr);}

    static size_t getBigAllocatorReservation(unsigned maxHitsToConsider, unsigned maxReadSize, unsigned seedLen, unsigned maxSeedsToUse);
    void *operator new(size_t size, BigAllocator *allocator) {_ASSERT(size == sizeof(SingleFirstPairedEndAligner)); return allocator->allocate(size);}
    void operator delete(void *ptr, BigAllocator *allocator) {/* do nothing, the owner of the allocator is responsible for cleaning up the memory */}

private:
    static const int INFINITE_SCORE = 0x7FFF;

    unsigned confDiff;
    unsigned maxK;
    unsigned minSpacing;
    unsigned maxSpacing;
    bool forceSpacing;
    bool skipAlignTogether;
    
    BaseAligner *singleAligner;
    BaseAligner *mateAligner;
    PairedEndAligner *underlyingPairedEndAligner;

    LandauVishkin<1> lv;
    LandauVishkin<-1> reverseLV;
};
