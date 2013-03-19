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
        bool                skipAlignTogether_,
        unsigned            alignTogetherLVLimit_,
        PairedEndAligner    *underlyingPairedEndAligner_)
 :   minSpacing(minSpacing_), maxSpacing(maxSpacing_), forceSpacing(forceSpacing_), 
     skipAlignTogether(skipAlignTogether_), alignTogetherLVLimit(alignTogetherLVLimit_),
     underlyingPairedEndAligner(underlyingPairedEndAligner_), lv(2 * FirstPowerOf2GreaterThanOrEqualTo(maxSeeds * maxHits * 4)),
     confDiff(confDiff_), maxK(maxK_)
{
    // Create single-end aligners.
    singleAligner = new BaseAligner(index, confDiff + 1, maxHits, maxK, maxReadSize,
                                    maxSeeds, adaptiveConfDiffThreshold);

    mateAligner = new BaseAligner(index, confDiff, maxHits, maxK, maxReadSize,
                                  maxSeeds, adaptiveConfDiffThreshold, &lv);
}


SingleFirstPairedEndAligner::~SingleFirstPairedEndAligner()
{
    delete singleAligner;
    delete mateAligner;
}


void SingleFirstPairedEndAligner::align(Read *read0, Read *read1, PairedAlignmentResult *result)
{
    Read *reads[NUM_READS_PER_PAIR] = {read0, read1};
    int numNotFound = 0;
    int numCertainlyNotFound = 0;
    int numIgnoredMulti = 0;  // MultiHits where all seeds returned too many hits
    int numSingleWithNotFound = 0;

    result->status[0] = NotFound;
    result->status[1] = NotFound;

    result->location[0] = 0xFFFFFFFF;
    result->location[1] = 0xFFFFFFFF;

    result->fromAlignTogether = false;
    result->nanosInAlignTogether = 0;

    unsigned bestScore[NUM_READS_PER_PAIR] = {INFINITE_SCORE, INFINITE_SCORE};
    unsigned bestMapq[NUM_READS_PER_PAIR] = {0, 0};
    unsigned singleLoc[NUM_READS_PER_PAIR] = {0xFFFFFFFF, 0xFFFFFFFF};
    Direction singleDirection[NUM_READS_PER_PAIR] = {FORWARD, FORWARD};
    bool certainlyNotFound[NUM_READS_PER_PAIR] = {false, false};
    bool bestScoreCertain[NUM_READS_PER_PAIR] = {false, false};
    AlignmentResult singleStatus[NUM_READS_PER_PAIR] = {NotFound, NotFound};
    
    TRACE("Aligning read pair %.*s:\n%.*s\n%.*s\n",
            reads[0]->getIdLength(), reads[0]->getId(),
            reads[0]->getDataLength(), reads[0]->getData(),
            reads[1]->getDataLength(), reads[1]->getData());

    if (read0->getDataLength() < 50 && read1->getDataLength() < 50) {
        TRACE("Reads are both too short -- returning");
        return;
    }

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

            status1 = mateAligner->AlignRead(reads[1-r], &loc1, &direction1, &score1, &mapq1, maxSpacing, loc0, OppositeDirection(direction0), NULL, NULL);

            TRACE("Mate %d returned %s at loc %u-%d\n", 1-r, AlignmentResultToString(status1), loc1, rc1);
            if (/*status1 != MultipleHits &&*/ score0 + score1 <= (int)maxK) {
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
                return;
            } else if(status1 == NotFound) {
                // We found read r at one location and didn't find the mate nearby. Let's remember because
                // if the mate is not found anywhere else, we can just return a location for read r.

                numSingleWithNotFound++;
            }
        } else if (status0 == NotFound) {
            numNotFound++;
            certainlyNotFound[r] = singleAligner->checkedAllSeeds();
            numCertainlyNotFound += certainlyNotFound[r] ? 1 : 0;
        } else if (status0 == MultipleHits && loc0 == 0xFFFFFFFF) {
            numIgnoredMulti++;
            result->status[r] = MultipleHits;
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
                result->status[r] = (singleLoc[r] != 0xFFFFFFFF && bestScore[r] <= 0.6 * maxK) ? SingleHit : NotFound;
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
            if (singleLoc[r] == 0xFFFFFFFF) {
                result->status[r] = NotFound;
            } else {
                result->status[r] = (bestScore[r] <= 0.6 * maxK ? SingleHit : NotFound);
                result->location[r] = (bestScore[r] <= 0.6 * maxK ? singleLoc[r] : 0xFFFFFFFF);
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

    if (bestScore[0] != INFINITE_SCORE && bestScore[1] != INFINITE_SCORE &&
            bestScore[0] + bestScore[1] > 1.5 * maxK) {
        TRACE("Best scores in each direction add up to more than 1.5 * maxK; returning");
        return;
    }

    int lowerBound[NUM_READS_PER_PAIR] = {0, 0};
    if (bestScoreCertain[0]) {
        lowerBound[0] = bestScore[0] < confDiff + 1 ? 0 : bestScore[0];
    }
    if (bestScoreCertain[1]) {
        lowerBound[1] = bestScore[1] < confDiff + 1 ? 0 : bestScore[1];
    }

    if (lowerBound[0] + lowerBound[1] > (int)maxK) {
        TRACE("Best scores in each direction certainly add up to more than maxK; returning");
        return;
    }

    if (skipAlignTogether) {
        return;
    } 
    // At least one read was MultipleHits; let's look for them simultaneously.
    _int64 start = timeInNanos();
    underlyingPairedEndAligner->align(read0, read1, result); 
    _int64 end = timeInNanos();

    result->nanosInAlignTogether = end - start;
    result->fromAlignTogether = true;

    if (forceSpacing) {
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
        }
    }

    TRACE("alignTogether took %lld ns and returned %s %s\n", end - start,
            AlignmentResultToString(result->status[0]),
            AlignmentResultToString(result->status[1]));
}
