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
        GenomeIndex         *index,
        unsigned            maxReadSize,
        unsigned            maxHits,
        unsigned            maxK,
        unsigned            maxSeedsFromCommandLine,
        double              seedCoverage,
        bool                forceSpacing_,
        unsigned            extraSearchDepth,
        PairedEndAligner    *underlyingPairedEndAligner_,
        BigAllocator        *allocator)
 :  underlyingPairedEndAligner(underlyingPairedEndAligner_), forceSpacing(forceSpacing_), lv(0), reverseLV(0)
{
    // Create single-end aligners.
    singleAligner = new (allocator) BaseAligner(index, maxHits, maxK, maxReadSize,
                                    maxSeedsFromCommandLine,  seedCoverage, extraSearchDepth, &lv, &reverseLV, NULL, allocator);
    
    underlyingPairedEndAligner->setLandauVishkin(&lv, &reverseLV);
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
        unsigned        maxCandidatePoolSize)
{
    return BaseAligner::getBigAllocatorReservation(false, maxHits, maxReadSize, seedLen, maxSeedsFromCommandLine, seedCoverage) + sizeof(ChimericPairedEndAligner) + sizeof(_uint64);
}


ChimericPairedEndAligner::~ChimericPairedEndAligner()
{
    singleAligner->~BaseAligner();
}

#ifdef _DEBUG
extern bool _DumpAlignments;
#endif // _DEBUG


void ChimericPairedEndAligner::align(Read *read0, Read *read1, PairedAlignmentResult *result)
{
	result->status[0] = result->status[1] = NotFound;
     if (read0->getDataLength() < 50 && read1->getDataLength() < 50) {
        TRACE("Reads are both too short -- returning");
        return;
    }

    _int64 start = timeInNanos();
    //
    // Let the LVs use the cache that we built up.
    //
    underlyingPairedEndAligner->align(read0, read1, result); 
    _int64 end = timeInNanos();

    result->nanosInAlignTogether = end - start;
    result->fromAlignTogether = true;
    result->alignedAsPair = true;

    if (forceSpacing) {
        if (result->status[0] == NotFound) {
            result->fromAlignTogether = false;
        } else {
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

    //
    // If the intersecting aligner didn't find an alignment for these reads, then they may be
    // chimeric and so we should just align them with the single end aligner and apply a MAPQ penalty.
    //
    Read *read[NUM_READS_PER_PAIR] = {read0, read1};
    for (int r = 0; r < NUM_READS_PER_PAIR; r++) {
        result->status[r] = singleAligner->AlignRead(read[r], &result->location[r], &result->direction[r], &result->score[r], &result->mapq[r]);
        result->mapq[r] /= 3;   // Heavy quality penalty for chimeric reads
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
