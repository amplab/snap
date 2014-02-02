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


void ChimericPairedEndAligner::align(Read *read0, Read *read1, PairedAlignmentResult *result, IdPairVector* secondary)
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
    underlyingPairedEndAligner->align(read0, read1, result, secondary); 
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
    IdPairVector* singleSecondary[2] = {NULL, NULL};
    if (secondary != NULL) {
        singleSecondary[0] = new IdPairVector;
        singleSecondary[1] = new IdPairVector;
    }
    for (int r = 0; r < NUM_READS_PER_PAIR; r++) {
        result->status[r] = singleAligner->AlignRead(read[r], &result->location[r], &result->direction[r], &result->score[r], &result->mapq[r], singleSecondary[r]);
        result->mapq[r] /= 3;   // Heavy quality penalty for chimeric reads
    }
    if (secondary != NULL && singleSecondary[0]->size() + singleSecondary[1]->size() > 0) {
        // loop through all combinations of secondary alignments
        secondary->clear();
        for (int i = -1; i < singleSecondary[0]->size(); i++) {
            for (int j = -1; j < singleSecondary[1]->size(); j++) {
                if (i > -1 || j > -1) {
                    if (i == -1) {
                        secondary->push_back(IdPair(result->location[0], result->direction[0]));
                    } else {
                        secondary->push_back((*singleSecondary[0])[i]);
                    }
                    if (j == -1) {
                        secondary->push_back(IdPair(result->location[1], result->direction[1]));
                    } else {
                        secondary->push_back((*singleSecondary[1])[j]);
                    }
                }
            }
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
