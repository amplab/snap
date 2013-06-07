/*++

Module Name:

    SingleFirstPairedEndAligner.cpp

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

#include "stdafx.h"
#include "SingleFirstPairedEndAligner.h"
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

SingleFirstPairedEndAligner::SingleFirstPairedEndAligner(
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
        unsigned            extraSearchDepth,
        bool                skipAlignTogether_,
        PairedEndAligner    *underlyingPairedEndAligner_)
 :   minSpacing(minSpacing_), maxSpacing(maxSpacing_), forceSpacing(forceSpacing_), 
     skipAlignTogether(skipAlignTogether_), 
     underlyingPairedEndAligner(underlyingPairedEndAligner_), lv(0/*2 * FirstPowerOf2GreaterThanOrEqualTo(maxSeeds * maxHits * 4)*/),
     reverseLV(0/*2 * FirstPowerOf2GreaterThanOrEqualTo(maxSeeds * maxHits * 4)*/),
     confDiff(confDiff_), maxK(maxK_)
{
    // Create single-end aligners.
    singleAligner = new BaseAligner(index, confDiff + 1, maxHits, maxK, maxReadSize,
                                    maxSeeds, adaptiveConfDiffThreshold, extraSearchDepth, &lv, &reverseLV);

    mateAligner = new BaseAligner(index, confDiff, maxHits, maxK, maxReadSize,
                                  maxSeeds, adaptiveConfDiffThreshold, extraSearchDepth, &lv, &reverseLV);

    underlyingPairedEndAligner->setLandauVishkin(&lv, &reverseLV);
}


SingleFirstPairedEndAligner::~SingleFirstPairedEndAligner()
{
    delete singleAligner;
    delete mateAligner;
}

#ifdef _DEBUG
extern bool _DumpAlignments;
#endif // _DEBUG


void SingleFirstPairedEndAligner::align(Read *read0, Read *read1, PairedAlignmentResult *result)
{
    Read *reads[NUM_READS_PER_PAIR] = {read0, read1};
    int numNotFound = 0;
    int numCertainlyNotFound = 0;
    int numIgnoredMulti = 0;  // MultiHits where all seeds returned too many hits
    int numSingleWithNotFound = 0;

    result->status[0] = NotFound;
    result->status[1] = NotFound;

    result->location[0] = InvalidGenomeLocation;
    result->location[1] = InvalidGenomeLocation;

    result->fromAlignTogether = false;
    result->alignedAsPair = false;
    result->nanosInAlignTogether = 0;

    unsigned bestScore[NUM_READS_PER_PAIR] = {INFINITE_SCORE, INFINITE_SCORE};
    unsigned mateScore[NUM_READS_PER_PAIR] = {INFINITE_SCORE, INFINITE_SCORE};
    unsigned bestMapq[NUM_READS_PER_PAIR] = {0, 0};
    unsigned mateMapq[NUM_READS_PER_PAIR] = {0, 0};
    unsigned singleLoc[NUM_READS_PER_PAIR] = {InvalidGenomeLocation, InvalidGenomeLocation};
    unsigned mateLoc[NUM_READS_PER_PAIR] = {InvalidGenomeLocation, InvalidGenomeLocation};
    Direction singleDirection[NUM_READS_PER_PAIR] = {FORWARD, FORWARD};
    Direction mateDirection[NUM_READS_PER_PAIR] = {FORWARD, FORWARD};
    bool certainlyNotFound[NUM_READS_PER_PAIR] = {false, false};
    bool bestScoreCertain[NUM_READS_PER_PAIR] = {false, false};
    AlignmentResult singleStatus[NUM_READS_PER_PAIR] = {NotFound, NotFound};
    AlignmentResult mateStatus[NUM_READS_PER_PAIR] = {NotFound, NotFound};
    
    TRACE("Aligning read pair %.*s:\n%.*s\n%.*s\n",
            reads[0]->getIdLength(), reads[0]->getId(),
            reads[0]->getDataLength(), reads[0]->getData(),
            reads[1]->getDataLength(), reads[1]->getData());

    if (read0->getDataLength() < 50 && read1->getDataLength() < 50) {
        TRACE("Reads are both too short -- returning");
        return;
    }
   
    lv.clearCache();
    reverseLV.clearCache();
#if 0

    // First try aligning each end alone and looking for its mate nearby if it was confident.
    // TODO: Guess which end will be faster to search for first based on it seeds?
    // TODO: We'll want to search for chimeras here when the second read is not found.
    for (int r = 0; r < NUM_READS_PER_PAIR; r++) {
        unsigned loc0, loc1;
        Direction direction0, direction1;
        int score0, score1;
        int mapq0, mapq1;
        AlignmentResult status0, status1;

        singleAligner->setReadId(r);
        status0 = singleAligner->AlignRead(reads[r], &loc0, &direction0, &score0, &mapq0);

        if (score0 <= (int) (maxK + confDiff)) {
            bestScore[r] = score0;
            bestMapq[r] = mapq0;
            singleLoc[r] = loc0;
            singleDirection[r] = direction0;
            singleStatus[r] = status0;
        } else {
            bestScore[r] = INFINITE_SCORE;
            bestMapq[r] = 0;
        }

        bestScoreCertain[r] = singleAligner->checkedAllSeeds();
        TRACE("Read %d returned %s (%d) at loc %u-%d\n", r, AlignmentResultToString(status0), score0, loc0, rc0);

        if (isOneLocation(status0)) {
            // Let's just search for the other read nearby.
            if (reads[1-r]->getDataLength() < 50) {
                TRACE("Read %d was found but %d is too short, so returning SingleHit/NotFound\n", r, 1-r);
                result->status[r] = SingleHit;
                result->location[r] = loc0;
                result->direction[r] = direction0;
                result->score[r] = score0;
                result->mapq[r] = mapq0;
                result->status[1-r] = NotFound;
                return;
            }

            mateAligner->setMaxK(maxK - score0 + 1);
            mateAligner->setReadId(1-r);

            status1 = mateAligner->AlignRead(reads[1-r], &loc1, &direction1, &score1, &mapq1, maxSpacing, loc0, OppositeDirection(direction0));

            TRACE("Mate %d returned %s at loc %u-%d\n", 1-r, AlignmentResultToString(status1), loc1, rc1);
            if (/*status1 != MultipleHits &&*/ score0 + score1 <= (int)maxK) {
                if (mapq0 >= 70 && mapq1 >= 70) {
                    result->status[r] = SingleHit;
                    result->location[r] = loc0;
                    result->direction[r] = direction0;
                    result->score[r] = score0;
                    result->mapq[r] = mapq0;
                    result->status[1-r] = isOneLocation(status1) ? SingleHit : status1;
                    result->location[1-r] = loc1;
                    result->direction[1-r] = direction1;
                    result->score[1-r] = score1;
                    result->mapq[1-r] = min(mapq0,mapq1);
                    result->alignedAsPair = true;
#ifdef _DEBUG
                    if (_DumpAlignments) {
                        printf("SingleFirstPairedEndAligner: too good to pass up (%u, %u) score (%d, %d), MAPQ (%d, %d)\n\n\n",result->location[0], result->location[1],
                            result->score[0], result->score[1], result->mapq[0], result->mapq[1]);
                    }
#endif // _DEBUG
                    return;
                }

                mateStatus[1-r] = status1;
                mateLoc[1-r] = loc1;
                mateMapq[1-r] = min(mapq0,mapq1);
                mateDirection[1-r] = direction1;
                mateScore[1-r] = score1;
           } else if(status1 == NotFound) {
                // We found read r at one location and didn't find the mate nearby. Let's remember because
                // if the mate is not found anywhere else, we can just return a location for read r.

                numSingleWithNotFound++;
            }
        } else if (status0 == NotFound) {
            numNotFound++;
            certainlyNotFound[r] = singleAligner->checkedAllSeeds();
            numCertainlyNotFound += certainlyNotFound[r] ? 1 : 0;
        } else if (status0 == MultipleHits && loc0 == InvalidGenomeLocation) {
            numIgnoredMulti++;
            result->status[r] = MultipleHits;
            result->mapq[r] = 0;
        }
    }


    
    if (numNotFound + numIgnoredMulti == 2) {
        // We couldn't find either read even as MultipleHits, so there's no way we'll do the pair.
        return;
    }

    if (numCertainlyNotFound == 1) {
        // One read was NotFound and we checked all its seeds, so there's no way we will find it as
        // either a bucket element or bucket mate in alignTogether. Let's just return the result for
        // the other one.
        TRACE("One read is certainly not found: %d %d %d, %d %d %d\n",
                certainlyNotFound[0], singleLoc[0], bestScore[0],
                certainlyNotFound[1], singleLoc[1], bestScore[1]);
        for (int r = 0; r < NUM_READS_PER_PAIR; r++) {
            if (certainlyNotFound[r]) {
                result->status[r] = NotFound;
            } else {
                result->status[r] = (singleLoc[r] != InvalidGenomeLocation && bestScore[r] <= 0.6 * maxK) ? SingleHit : NotFound;
                result->location[r] = singleLoc[r];
                result->direction[r] = singleDirection[r];
                result->score[r] = bestScore[r];
                result->mapq[r] = bestMapq[r];
            }
        }

        return;
    }

    if (numSingleWithNotFound == 1 && numNotFound == 1) {
        // One read was Single and the other was NotFound, neither near its mate nor elsewhere.
        TRACE("One read was single and its mate was NotFound either nearby or by itself\n");
        for (int r = 0; r < NUM_READS_PER_PAIR; r++) {
            if (singleLoc[r] == InvalidGenomeLocation) {
                result->status[r] = NotFound;
            } else {
                result->status[r] = (bestScore[r] <= 0.6 * maxK ? SingleHit : NotFound);
                result->location[r] = (bestScore[r] <= 0.6 * maxK ? singleLoc[r] : InvalidGenomeLocation);
                result->direction[r] = singleDirection[r];
                result->score[r] = bestScore[r];
                result->mapq[r] = bestMapq[r];
            }
        }

        return;
    }

    if (read0->getDataLength() < 50 || read1->getDataLength() < 50) {
        TRACE("Can't go further because one read is too short");
        return;
    }

 
    if (isOneLocation(singleStatus[0]) && bestScore[0] + mateScore[1] <= (int)maxK || isOneLocation(singleStatus[1]) && bestScore[1] + mateScore[0] <= (int)maxK) {
        //
        // Good enough to quit.  Take the best one.
        //
        int r;
        if (bestMapq[0] + mateMapq[1] >= bestMapq[1] + mateMapq[0]) {
            r = 0;
        } else {
            r = 1;
        }

        result->status[r] = SingleHit;
        result->location[r] = singleLoc[r];
        result->direction[r] = singleDirection[r];
        result->score[r] = bestScore[r];
        result->mapq[r] = bestMapq[r];
        result->status[1-r] = isOneLocation(mateStatus[1-r]) ? SingleHit : mateStatus[1-r];
        result->location[1-r] = mateLoc[1-r];
        result->direction[1-r] = mateDirection[1-r];
        result->score[1-r] = mateScore[1-r];
        result->mapq[1-r] = mateMapq[1-r];
        result->alignedAsPair = true;

#ifdef _DEBUG
                    if (_DumpAlignments) {
                        printf("SingleFirstPairedEndAligner: single was good enough (%u, %u) score (%d, %d), MAPQ (%d, %d)\n\n\n",result->location[0], result->location[1],
                            result->score[0], result->score[1], result->mapq[0], result->mapq[1]);
                    }
#endif // _DEBUG

        return;
    }

    if (skipAlignTogether) {
        return;
    }

#endif // 0

    // At least one read was MultipleHits; let's look for them simultaneously.
    _int64 start = timeInNanos();
    //
    // Let the LVs use the cache that we built up.
    //
    lv.enterCacheLookupPhase();
    reverseLV.enterCacheLookupPhase();
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

    //
    // If the intersecting aligner didn't find an alignment for these reads, maybe we've got
    // something strange like reads around a structural variation of some kind.  Just use the
    // results of the single aligner, but with MAPQ reduced significantly.
    //
    for (int r = 0; r < NUM_READS_PER_PAIR; r++) {
        if (result->status[r] == NotFound) {
            result->fromAlignTogether = false;
            result->location[r] = singleLoc[r];
            result->direction[r] = singleDirection[r];
            result->mapq[r] = bestMapq[r] / 4;
            result->score[r] = bestScore[r];
            result->status[r] = singleStatus[r];
            result->alignedAsPair = false;
        }
    }

#ifdef _DEBUG
    if (_DumpAlignments) {
        printf("SingleFirstPairedEndAligner: from align together/chimeric (%u, %u) score (%d, %d), MAPQ (%d, %d)\n\n\n",result->location[0], result->location[1],
            result->score[0], result->score[1], result->mapq[0], result->mapq[1]);
    }
#endif // _DEBUG
                    
    TRACE("alignTogether took %lld ns and returned %s %s\n", end - start,
        AlignmentResultToString(result->status[0]),
        AlignmentResultToString(result->status[1]));
}
