/*++

Module Name:

    SimplePairedEndAligner.cpp

Abstract:

    First cut at a paired-end aligner.

Authors:

    Matei Zaharia, December, 2011

Environment:

    User mode service.

Revision History:

--*/

#include "stdafx.h"
#include <math.h>
#include "SimplePairedEndAligner.h"
#include "LandauVishkin.h"

using namespace std;


namespace {

struct Hit {
    unsigned location;
    bool isRC;
    int score;
    
    Hit() {}
    
    Hit(unsigned _location, bool _isRC, int _score)
        : location(_location), isRC(_isRC), score(_score) {}
};

bool compareHits(const Hit& hit1, const Hit& hit2) {
    return hit1.location < hit2.location;
}

}

SimplePairedEndAligner::SimplePairedEndAligner(
         GenomeIndex  *_index,
         unsigned      _maxReadSize,
         unsigned      _confDiff, 
         unsigned      _maxHits, 
         unsigned      _maxK,
         unsigned      _maxSeeds,
         unsigned      _maxSpacing,
         unsigned      _adaptiveConfDiffThreshold,
         unsigned      _multiHitConfDiff,
         unsigned      _multiHitConfDiffThreshold)
 :   index(_index), maxReadSize(_maxReadSize), confDiff(_confDiff), maxHits(_maxHits), maxK(_maxK),
     maxSeeds(_maxSeeds), maxSpacing(_maxSpacing), adaptiveConfDiffThreshold(_adaptiveConfDiffThreshold),
     multiHitConfDiff(_multiHitConfDiff), multiHitConfDiffThreshold(_multiHitConfDiffThreshold)
{
    singleAligner = new BaseAligner(index, confDiff, maxHits, maxK, maxReadSize, maxSeeds, 
                                    1000000 /* lvLimit */, adaptiveConfDiffThreshold);
}

SimplePairedEndAligner::~SimplePairedEndAligner()
{
    
}

void SimplePairedEndAligner::align(Read *read0, Read *read1, PairedAlignmentResult *result)
{   
    // First, try to align both reads single-endedly (?).
    AlignmentResult singleStatus[2];
    bool singleIsRC[2];
    unsigned singleLoc[2];
    
    singleStatus[0] = singleAligner->AlignRead(read0, &singleLoc[0], &singleIsRC[0], NULL, 0, 0, false,
                        MAX_HITS_PER_READ, &numMultiHits[0], multiHitLocs[0], multiHitRCs[0], multiHitScores[0]);
    singleStatus[1] = singleAligner->AlignRead(read1, &singleLoc[1], &singleIsRC[1], NULL, 0, 0, false,
                        MAX_HITS_PER_READ, &numMultiHits[1], multiHitLocs[1], multiHitRCs[1], multiHitScores[1]);

#ifdef TRACE_PAIRED_ALIGNER
    printf("\nLooking at read %.*s and its mate:\n%.*s\n%.*s\n",
        read0->getIdLength(), read0->getId(), read0->getDataLength(), read0->getData(),
        read1->getDataLength(), read1->getData());
    printf("Single-end alignment returned %s %s (%u %u)\n",
        AlignmentResultToString(singleStatus[0]),
        AlignmentResultToString(singleStatus[1]),
        isOneLocation(singleStatus[0]) ? singleLoc[0] : 0,
        isOneLocation(singleStatus[1]) ? singleLoc[1] : 0);
#endif
    
    if (isOneLocation(singleStatus[0]) && isOneLocation(singleStatus[1]))
    {
        //
        // Both reads aligned to single locations. If they're close enough to each other and in opposing
        // orientations, as is to be expected for paired-end reads, go with those. Otherwise, it's quite
        // likely that one read is wrong, or we might have a translocation event (very rare). Check for
        // good alignments of each read near the other read too.
        //        
        result->status[0] = SingleHit;
        result->location[0] = singleLoc[0];
        result->isRC[0] = singleIsRC[0];
        result->status[1] = SingleHit;
        result->location[1] = singleLoc[1];
        result->isRC[1] = singleIsRC[1];
        unsigned distance = max(singleLoc[0], singleLoc[1]) - min(singleLoc[0], singleLoc[1]);
        if (distance > maxSpacing || singleIsRC[0] == singleIsRC[1]) {
            // Find where each read aligns around the other read.
            unsigned newLoc[2];
            bool newIsRC[2];
            AlignmentResult newStatus[2];
            newStatus[0] = singleAligner->AlignRead(
                read0, &newLoc[0], &newIsRC[0], NULL, maxSpacing, singleLoc[1], !singleIsRC[1]);
            newStatus[1] = singleAligner->AlignRead(
                read1, &newLoc[1], &newIsRC[1], NULL, maxSpacing, singleLoc[0], !singleIsRC[0]);
            if (isOneLocation(newStatus[1]) && !isOneLocation(newStatus[0])) {
                // Move read 1 to the location found around read 0
                result->status[1] = SingleHit;
                result->location[1] = newLoc[1];
                result->isRC[1] = newIsRC[1];
            } else if (isOneLocation(newStatus[0]) && !isOneLocation(newStatus[1])) {
                // Move read 1 to the location found around read 0
                result->status[0] = SingleHit;
                result->location[0] = newLoc[0];
                result->isRC[0] = newIsRC[0];
            } else if (isOneLocation(newStatus[0]) && isOneLocation(newStatus[1])) {
                // Well, this is weird; the reads aligned confidently alone, but they also
                // align well to two different sets of locations next to each other. We should
                // return MultipleHits to be safe.
                result->status[0] = MultipleHits;
                result->status[1] = MultipleHits;
                //fprintf(stderr, "WARNING: Read %.*s and its mate match two regions in a skewed way\n",
                //    read0->getIdLength(), read0->getId());
            } else {
                // Neither read aligns well near the other one. Maybe this really is a translocation.
                // Go with the original positions found by single-end alignment.
                fprintf(stderr, "WARNING: Read %.*s and its mate match overly distant locations\n",
                    read0->getIdLength(), read0->getId());
            }
        }
    } else if (isOneLocation(singleStatus[0]) && singleStatus[1] == MultipleHits) {
        //
        // Read 0 aligned in one place but read 1 didn't; look for read 1 in a band around read 0.
        //
        // TODO: Also consider searching with a higher maxSeeds if the status of read 1 is NotFound,
        // in case we missed it due to including many errors; in general we may want to use LV against
        // a larger interval around read 1 here.
        //
        unsigned newLoc;
        bool newIsRC;
#ifdef TRACE_PAIRED_ALIGNER
        //const char *data = index->getGenome()->getSubstring(singleLoc[0] - maxSpacing, 2 * maxSpacing);
        //if (data != NULL) {
        //    printf("Searching the following region for mate:\n%.*s\n", 2 * maxSpacing, data);
        //} else {
        //    printf("Searching the following region for mate:\n<NULL>\n");
        //}
#endif
        AlignmentResult newStatus = singleAligner->AlignRead(
            read1, &newLoc, &newIsRC, NULL, maxSpacing, singleLoc[0], !singleIsRC[0]);
        if (isOneLocation(newStatus)) {
            // Hurray! We've narrowed it down to a single location.
            result->status[1] = SingleHit;
            result->location[1] = newLoc;
            result->isRC[1] = !singleIsRC[0];
        } else {
            // Guess we'll have to return MultipleHits for the other read.
            result->status[1] = MultipleHits;
        }
        // In either case, report our one location for read 0.
        result->status[0] = SingleHit;
        result->location[0] = singleLoc[0];
        result->isRC[0] = singleIsRC[0];
    } else if (isOneLocation(singleStatus[1]) && singleStatus[0] == MultipleHits) {
        //
        // Read 1 aligned in one place but read 0 didn't; look for read 0 in a band around read 1.
        //
        // TODO: Also consider searching with a higher maxSeeds if the status of read 1 is NotFound,
        // in case we missed it due to including many errors; in general we may want to use LV against
        // a larger interval around read 1 here.
        //
        unsigned newLoc;
        bool newIsRC;
        AlignmentResult newStatus = singleAligner->AlignRead(
            read0, &newLoc, &newIsRC, NULL, maxSpacing, singleLoc[1], !singleIsRC[1]);
        if (isOneLocation(newStatus)) {
            // Hurray! We've narrowed it down to a single location.
            result->status[0] = SingleHit;
            result->location[0] = newLoc;
            result->isRC[0] = !singleIsRC[1];
        } else {
            // Guess we'll have to return MultipleHits for the other read.
            result->status[0] = MultipleHits;
        }
        // In either case, report our one location for read 1.
        result->status[1] = SingleHit;
        result->location[1] = singleLoc[1];
        result->isRC[1] = singleIsRC[1];
    } else if (singleStatus[0] == MultipleHits && singleStatus[1] == MultipleHits) {
        //
        // Both reads gave multiple hits. Let's try to find a pair of locations out of these hits where
        // they both work well, but if there are too many such pairs of locations, return MultipleHits.
        //
        // TODO: We really should be more careful in dealing with very close hits due to indels here.
        //       Skipping this right now for simplicity.
        //

        handleMultiHits2(read0, read1, result);
        //result->status[0] = MultipleHits;
        //result->status[1] = MultipleHits;
    } else {
        //
        // At least one read gave NotFound; we can't do anything for the other but treat it as single-ended.
        //
        result->status[0] = isOneLocation(singleStatus[0]) ? SingleHit : singleStatus[0];
        result->location[0] = singleLoc[0];
        result->isRC[0] = singleIsRC[0];
        result->status[1] = isOneLocation(singleStatus[1]) ? SingleHit : singleStatus[1];
        result->location[1] = singleLoc[1];
        result->isRC[1] = singleIsRC[1];
    }
}


void SimplePairedEndAligner::handleMultiHits1(Read *read0, Read *read1, PairedAlignmentResult *result) {
    // First put the (location, RC) pairs found into arrays to sort them simultaneously
    Hit hits[2][MAX_HITS_PER_READ];
    for (int i = 0; i < numMultiHits[0]; i++) {
        hits[0][i].location = multiHitLocs[0][i];
        hits[0][i].isRC = multiHitRCs[0][i];
        hits[0][i].score = multiHitScores[0][i];
    }
    for (int i = 0; i < numMultiHits[1]; i++) {
        hits[1][i].location = multiHitLocs[1][i];
        hits[1][i].isRC = multiHitRCs[1][i];
        hits[1][i].score = multiHitScores[1][i];
    }
    
    sort(hits[0], hits[0] + numMultiHits[0], compareHits);
    sort(hits[1], hits[1] + numMultiHits[1], compareHits);
    
#ifdef TRACE_PAIRED_ALIGNER
    printf("Sorted hit locations (%d / %d):\nRead0:", numMultiHits[0], numMultiHits[1]);
    for (int i = 0; i < numMultiHits[0]; i++) {
        printf(" %u-%d", hits[0][i].location, hits[0][i].isRC);
    }
    printf("\nRead1:");
    for (int i = 0; i < numMultiHits[1]; i++) {
        printf(" %u-%d", hits[1][i].location, hits[1][i].isRC);
    }
    printf("\n");
#endif
    
    int bestScore = 100000;
    int secondBestScore = 100000;
    int goodPairs = 0;

    int maxScore = maxK;
    
    int start1 = 0; // Where in array 1 to start searching for matches (to avoid quadratic-time search).
    for (int pos0 = 0; pos0 < numMultiHits[0]; pos0++) {
        unsigned loc0 = hits[0][pos0].location;
        bool rc0 = hits[0][pos0].isRC;
        // Advance start1 until we are in the right distance of the hit at pos0
        while (start1 < numMultiHits[1] && hits[1][start1].location + maxSpacing < loc0) {
            start1++;
        }
        for (int pos1 = start1; pos1 < numMultiHits[1]; pos1++) {
            unsigned loc1 = hits[1][pos1].location;
            bool rc1 = hits[1][pos1].isRC;
            unsigned distance = max(loc0, loc1) - min(loc0, loc1);
            if (distance <= maxSpacing && rc0 != rc1) {
                // We found a good pair of hits. Let's score it.
                int score = hits[0][pos0].score + hits[1][pos1].score;
#ifdef TRACE_PAIRED_ALIGNER
                printf("Scoring pair %u %u [%u], %d %d => %d\n",
                    loc0, loc1, distance, hits[0][pos0].score, hits[1][pos1].score, score);
#endif              
                if (score <= bestScore) {
                    secondBestScore = bestScore;
                    bestScore = score;
                    result->location[0] = loc0;
                    result->isRC[0] = rc0;
                    result->location[1] = loc1;
                    result->isRC[1] = rc1;
                } else if (score < secondBestScore) {
                    secondBestScore = score;
                }
                if (score <= maxScore + (int)multiHitConfDiff + 1) {
                    goodPairs++;
                }
            } else if (loc1 > loc0 + maxSpacing) {
                // We've gone too far into array 1 relative to our position in array 0, so break.
                break;
            }
        }
    }
    
    int realConfDiff = multiHitConfDiff + (goodPairs > (int)multiHitConfDiffThreshold ? 1 : 0);
    if (bestScore <= maxScore && bestScore + realConfDiff <= secondBestScore) {
        result->status[0] = SingleHit;
        result->status[1] = SingleHit;
#ifdef TRACE_PAIRED_ALIGNER
        printf("Returning SingleHits %u %u\n", result->location[0], result->location[1]);
#endif
    } else if (bestScore > maxScore) {
        result->status[0] = NotFound;
        result->status[1] = NotFound;
#ifdef TRACE_PAIRED_ALIGNER
        printf("Could not find pair of matching hits for read %.*s\n", read0->getIdLength(), read0->getId());
#endif
    } else {
        result->status[0] = MultipleHits;
        result->status[1] = MultipleHits;
#ifdef TRACE_PAIRED_ALIGNER
        printf("Returning MultipleHits because scores were too close\n");
#endif
    }
}


void SimplePairedEndAligner::handleMultiHits2(Read *read0, Read *read1, PairedAlignmentResult *result) {
    // First put the (location, RC) pairs found into arrays
    Hit hits[2][MAX_HITS_PER_READ];
    for (int i = 0; i < numMultiHits[0]; i++) {
        hits[0][i].location = multiHitLocs[0][i];
        hits[0][i].isRC = multiHitRCs[0][i];
        hits[0][i].score = multiHitScores[0][i];
    }
    for (int i = 0; i < numMultiHits[1]; i++) {
        hits[1][i].location = multiHitLocs[1][i];
        hits[1][i].isRC = multiHitRCs[1][i];
        hits[1][i].score = multiHitScores[1][i];
    }
    
#ifdef TRACE_PAIRED_ALIGNER
    printf("Hit locations (%d / %d):\nRead0:", numMultiHits[0], numMultiHits[1]);
    for (int i = 0; i < numMultiHits[0]; i++) {
        printf(" %u-%d-%d", hits[0][i].location, hits[0][i].isRC, hits[0][i].score);
    }
    printf("\nRead1:");
    for (int i = 0; i < numMultiHits[1]; i++) {
        printf(" %u-%d-%d", hits[1][i].location, hits[1][i].isRC, hits[1][i].score);
    }
    printf("\n");
#endif
    
    int bestScore = 100000;
    int secondBestScore = 100000;
    unsigned goodPairs = 0;

    int maxScore = maxK;

    // Look for read 1 around each location of read 0
    if (numMultiHits[0] < 20 && numMultiHits[1] < 20) {
        for (int i = 0; i < numMultiHits[0]; i++) {
            unsigned loc0 = hits[0][i].location;
            bool rc0 = hits[0][i].isRC;
            unsigned otherLoc;
            bool otherIsRC;
            int otherScore;
            AlignmentResult otherStatus = singleAligner->AlignRead(
                read1, &otherLoc, &otherIsRC, &otherScore, maxSpacing, loc0, !rc0);
            if (isOneLocation(otherStatus)) {
                unsigned distance = max(loc0, otherLoc) - min(loc0, otherLoc);
                if (distance <= maxSpacing && rc0 != otherIsRC) {
                    // We found a good pair of hits. Let's score it.
                    int score = hits[0][i].score + otherScore;
#ifdef TRACE_PAIRED_ALIGNER
                    printf("Scoring pair %u %u [%u], %d %d => %d\n",
                        loc0, otherLoc, distance, hits[0][i].score, otherScore, score);
#endif              
                    if (score <= bestScore) {
                        secondBestScore = bestScore;
                        bestScore = score;
                        result->location[0] = loc0;
                        result->isRC[0] = rc0;
                        result->location[1] = otherLoc;
                        result->isRC[1] = otherIsRC;
                    } else if (score < secondBestScore) {
                        secondBestScore = score;
                    }
                    if (score <= maxScore + (int)multiHitConfDiff + 1) {
                        goodPairs++;
                    }
                }
            }
        }
    }
   
    // Look for read 0 around each location of read 1
    if (numMultiHits[1] < 20 && numMultiHits[0] < 20) {
        for (int i = 0; i < numMultiHits[1]; i++) {
            unsigned loc1 = hits[1][i].location;
            bool rc1 = hits[1][i].isRC;
            unsigned otherLoc;
            bool otherIsRC;
            int otherScore;
            AlignmentResult otherStatus = singleAligner->AlignRead(
                read0, &otherLoc, &otherIsRC, &otherScore, maxSpacing, loc1, !rc1);
            if (isOneLocation(otherStatus)) {
                unsigned distance = max(loc1, otherLoc) - min(loc1, otherLoc);
                if (distance <= maxSpacing && rc1 != otherIsRC) {
                    // We found a good pair of hits. Let's score it.
                    int score = hits[1][i].score + otherScore;
#ifdef TRACE_PAIRED_ALIGNER
                    printf("Scoring pair %u %u [%u], %d %d => %d\n",
                        otherLoc, loc1, distance, otherScore, hits[1][i].score, score);
#endif              
                    if (otherLoc == result->location[0] && loc1 == result->location[1]) {
                        continue;
                    } else if (score <= bestScore) {
                        secondBestScore = bestScore;
                        bestScore = score;
                        result->location[0] = otherLoc;
                        result->isRC[0] = otherIsRC;
                        result->location[1] = loc1;
                        result->isRC[1] = rc1;
                    } else if (score < secondBestScore) {
                        secondBestScore = score;
                    }
                    if (score <= maxScore + (int)multiHitConfDiff + 1) {
                        goodPairs++;
                    }
                }
            }
        }
    }
    
    int realConfDiff = multiHitConfDiff + (goodPairs > multiHitConfDiffThreshold ? 1 : 0);
    if (bestScore <= maxScore && bestScore + realConfDiff <= secondBestScore) {
        result->status[0] = SingleHit;
        result->status[1] = SingleHit;
#ifdef TRACE_PAIRED_ALIGNER
        printf("Returning SingleHits %u %u\n", result->location[0], result->location[1]);
#endif
    } else if (bestScore > maxScore) {
        result->status[0] = NotFound;
        result->status[1] = NotFound;
#ifdef TRACE_PAIRED_ALIGNER
        printf("Could not find pair of matching hits for read %.*s\n", read0->getIdLength(), read0->getId());
#endif
    } else {
        result->status[0] = MultipleHits;
        result->status[1] = MultipleHits;
#ifdef TRACE_PAIRED_ALIGNER
        printf("Returning MultipleHits because scores were too close\n");
#endif
    }
}
