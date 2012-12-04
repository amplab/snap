/*++

Module Name:

    SmarterPairedEndAligner.cpp

Abstract:

    A more sophisticated paired-end aligner.

Authors:

    Matei Zaharia, February, 2012

Environment:

    User mode service.

Revision History:

--*/

#include "stdafx.h"
#include <math.h>
#include "SmarterPairedEndAligner.h"
#include "LandauVishkin.h"

using namespace std;


#ifdef TRACE_PAIRED_ALIGNER
#define TRACE printf
#else
#define TRACE(...) {}
#endif


SmarterPairedEndAligner::SmarterPairedEndAligner(
         GenomeIndex  *index_,
         unsigned      maxReadSize_,
         unsigned      confDiff_,
         unsigned      maxHits_,
         unsigned      maxK_,
         unsigned      maxSeeds_,
         unsigned      minSpacing_,
         unsigned      maxSpacing_,
         unsigned      adaptiveConfDiffThreshold_)
 :   index(index_), maxReadSize(maxReadSize_), confDiff(confDiff_), maxHits(maxHits_),
     maxK(maxK_), maxSeeds(maxSeeds_), minSpacing(minSpacing_), maxSpacing(maxSpacing_),
     adaptiveConfDiffThreshold(adaptiveConfDiffThreshold_), seedLen(index_->getSeedLength()),
     lv(2 * MAX_BUCKETS)
{
    // Initialize the bucket data structures.
    buckets = new Bucket[MAX_BUCKETS];
    for (int r = 0; r < 2; r++) {
        for (int rc = 0; rc < 2; rc++) {
            bucketLocations[r][rc].reserve(MAX_BUCKETS);
            bucketTable[r][rc].reserve(2 * MAX_BUCKETS);
        }
    }
    candidates.reserve(MAX_BUCKETS);
    
    // Initialize the complements array.
    memset(complement, 0, sizeof(complement));
    complement['A'] = 'T';
    complement['C'] = 'G';
    complement['G'] = 'C';
    complement['T'] = 'A';
    complement['N'] = 'N';
    complement['n'] = 'n';
    
    // Initialize the seed offsets used when we wrap around in testing seeds.
    vector<bool> used(seedLen + 1, false);
    used[0] = true;
    used[seedLen] = true;      // Sentinel to make the calculation work nicely
    wrapOffset[0] = 0;
    for (int i = 1; i < seedLen; i++) {
        // For the next wrap offset, pick the offset that is most distant from used offsets.
        int bestDist = 0;
        int bestP = 0;
        for (int p = 1; p < seedLen; p++) {
            int dist = 0;
            while (!used[p - dist] && !used[p + dist]) {
                dist++;
            }
            if (dist > bestDist) {
                bestDist = dist;
                bestP = p;
            }
        }
        wrapOffset[i] = bestP;
        used[bestP] = 1;
    }
#ifdef TRACE_PAIRED_ALIGNER
    printf("Wrap offsets:");
    for (int i = 0; i < seedLen; i++) {
        printf(" %d", wrapOffset[i]);
    }
    printf("\n");
#endif

    // Create single-end aligners.
    singleAligner = new BaseAligner(index, confDiff + 1, maxHits, maxK, maxReadSize,
                                    maxSeeds, 1000000 /* lvLimit */, adaptiveConfDiffThreshold);
    mateAligner = new BaseAligner(index, confDiff, maxHits, maxK, maxReadSize,
                                  maxSeeds, 1000000 /* lvLimit */, adaptiveConfDiffThreshold, &lv);
}


SmarterPairedEndAligner::~SmarterPairedEndAligner()
{
    delete[] buckets;
}


void SmarterPairedEndAligner::align(Read *read0, Read *read1, PairedAlignmentResult *result)
{
    Read *reads[2] = {read0, read1};
    int numNotFound = 0;
    int numCertainlyNotFound = 0;
    int numIgnoredMulti = 0;  // MultiHits where all seeds returned too many hits
    int numSingleWithNotFound = 0;

    clearState();
    
    result->status[0] = NotFound;
    result->status[1] = NotFound;

    result->location[0] = 0xFFFFFFFF;
    result->location[1] = 0xFFFFFFFF;

    unsigned bestScore[2] = {INFINITE_SCORE, INFINITE_SCORE};
    unsigned singleLoc[2] = {0xFFFFFFFF, 0xFFFFFFFF};
    bool singleIsRC[2] = {false, false};
    bool certainlyNotFound[2] = {false, false};
    bool bestScoreCertain[2] = {false, false};
    
    TRACE("Aligning read pair %.*s:\n%.*s\n%.*s\n",
            reads[0]->getIdLength(), reads[0]->getId(),
            reads[0]->getDataLength(), reads[0]->getData(),
            reads[1]->getDataLength(), reads[1]->getData());

    if (read0->getDataLength() < 50 && read1->getDataLength() < 50) {
        TRACE("Reads are both too short -- returning");
        return;
    }

#define SINGLE_FIRST
#ifdef SINGLE_FIRST
    // First try aligning each end alone and looking for its mate nearby if it was confident.
    // TODO: Guess which end will be faster to search for first based on it seeds?
    // TODO: We'll want to search for chimeras here when the second read is not found.
    for (int r = 0; r < 2; r++) {
        unsigned loc0, loc1;
        bool rc0, rc1;
        int score0, score1;
        AlignmentResult status0, status1;
        singleAligner->setReadId(r);
        status0 = singleAligner->AlignRead(reads[r], &loc0, &rc0, &score0);
        bestScore[r] = (score0 <= (int) (maxK + confDiff ? score0 : INFINITE_SCORE));
        bestScoreCertain[r] = singleAligner->checkedAllSeeds();
        TRACE("Read %d returned %s (%d) at loc %u-%d\n", r, AlignmentResultToString(status0), score0, loc0, rc0);
        if (isOneLocation(status0)) {
            // Let's just search for the other read nearby.
            if (reads[1-r]->getDataLength() < 50) {
                TRACE("Read %d was found but %d is too short, so returning SingleHit/NotFound\n", r, 1-r);
                result->status[r] = SingleHit;
                result->location[r] = loc0;
                result->isRC[r] = rc0;
                result->score[r] = score0;
                result->status[1-r] = NotFound;
                return;
            }
            mateAligner->setMaxK(maxK - score0 + 1);
            mateAligner->setReadId(1-r);
            status1 = mateAligner->AlignRead(reads[1-r], &loc1, &rc1, &score1, maxSpacing, loc0, !rc0);
            TRACE("Mate %d returned %s at loc %u-%d\n", 1-r, AlignmentResultToString(status1), loc1, rc1);
            if (/*status1 != MultipleHits &&*/ score0 + score1 <= (int)maxK) {
                result->status[r] = SingleHit;
                result->location[r] = loc0;
                result->isRC[r] = rc0;
                result->score[r] = score0;
                result->status[1-r] = isOneLocation(status1) ? SingleHit : status1;
                result->location[1-r] = loc1;
                result->isRC[1-r] = rc1;
                result->score[1-r] = score1;
                return;
            } else if(status1 == NotFound) {
                // We found read r at one location and didn't find the mate nearby. Let's remember because
                // if the mate is not found anywhere else, we can just return a location for read r.
                singleLoc[r] = loc0;
                singleIsRC[r] = rc0;
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
        for (int r = 0; r < 2; r++) {
            if (certainlyNotFound[r]) {
                result->status[r] = NotFound;
            } else {
                result->status[r] = (singleLoc[r] != 0xFFFFFFFF && bestScore[r] <= 0.6 * maxK) ? SingleHit : NotFound;
                result->location[r] = singleLoc[r];
                result->isRC[r] = singleIsRC[r];
                result->score[r] = bestScore[r];
            }
        }
        return;
    }

    if (numSingleWithNotFound == 1 && numNotFound == 1) {
        // One read was Single and the other was NotFound, neither near its mate nor elsewhere.
        TRACE("One read was single and its mate was NotFound either nearby or by itself\n");
        for (int r = 0; r < 2; r++) {
            if (singleLoc[r] == 0xFFFFFFFF) {
                result->status[r] = NotFound;
            } else {
                result->status[r] = (bestScore[r] <= 0.6 * maxK ? SingleHit : NotFound);
                result->location[r] = (bestScore[r] <= 0.6 * maxK ? singleLoc[r] : 0xFFFFFFFF);
                result->isRC[r] = singleIsRC[r];
                result->score[r] = bestScore[r];
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

    int lowerBound[2] = {0, 0};
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

#endif // SINGLE_FIRST

    // At least one read was MultipleHits; let's look for them simultaneously.
    _int64 start = timeInNanos();
    alignTogether(reads, result, lowerBound);
    _int64 end = timeInNanos();
    TRACE("alignTogether took %lld ns and returned %s %s\n", end - start,
            AlignmentResultToString(result->status[0]),
            AlignmentResultToString(result->status[1]));
}


void SmarterPairedEndAligner::alignTogether(Read *reads[2], PairedAlignmentResult *result, int lowerBound[2])
{
    int readLen[2];
    const char *readData[2][2];     // [read][isRC]
    char rcData[2][MAX_READ_SIZE];
    
    for (int r = 0; r < 2; r++) {
        readLen[r] = reads[r]->getDataLength();
        readData[r][0] = reads[r]->getData();
        computeRC(reads[r], &rcData[r][0]);
        readData[r][1] = &rcData[r][0];
    }
    
    int bestScore = INFINITE_SCORE;
    int secondBestScore = INFINITE_SCORE;
    
    // Locations of best scoring pair
    int bestLoc[2] = {0xFFFFFFFF, 0xFFFFFFFF};
    int bestRC[2];
    int bestIsMulti[2] = {0, 0};
    int bestScores[2] = {INFINITE_SCORE, INFINITE_SCORE};
    
    int disjointSeedsUsed[2][2] = {{0, 0}, {0, 0}};   // [read][isRC]
    int popularSeeds[2][2] = {{0, 0}, {0, 0}};        // [read][isRC]
    int seedHits[2][2] = {{0, 0}, {0, 0}};            // [read][isRC]
    int totalPopularSeeds = 0;
    
    int nextSeed[2] = {0, 0};        // Which seed to try next on each read
    int wrapCount[2] = {0, 0};       // How many times we've wraped around on each read

    int scoreLimit = maxK + confDiff + 1;

    unsigned candidatesTried = 0;

    unsigned seedsTried = 0;
    bool shouldStop = false;
    while (seedsTried < maxSeeds && !shouldStop) {
        // Try one seed from each read just to keep things simple.
        for (int r = 0; r < 2; r++) {
            // TODO: We could stop searching for one of the reads early if we rule it out via confDiff?
            if (wrapCount[r] >= seedLen) {
                continue;           // We've tried all possible seeds for this read
            }
            
            // Look up the current seed.
            int seedPos = nextSeed[r];
            TRACE("Trying read %d seed %d: %.*s\n", r, seedPos, seedLen, readData[r][0] + seedPos);
            if (Seed::DoesTextRepresentASeed(readData[r][0] + seedPos, seedLen)) {
                Seed seed(readData[r][0] + seedPos, seedLen);
                unsigned nHits[2];          // Number of hits in forward and RC directions
                const unsigned *hits[2];    // Actual hit locations in each direction
                index->lookupSeed(seed, &nHits[0], &hits[0], &nHits[1], &hits[1]);
                
#ifdef TRACE_PAIRED_ALIGNER
                printf("  %u fwd hits:", nHits[0]);
                //for (int i = 0; i < nHits[0]; i++) { printf(" %u", hits[0][i]); }
                printf("\n");
                printf("  %u RC hits:", nHits[1]);
                //for (int i = 0; i < nHits[1]; i++) { printf(" %u", hits[1][i]); }
                printf("\n");
#endif
            
                // Add in all the hits.
                for (int rc = 0; rc < 2; rc++) {
                    seedHits[r][rc] += nHits[rc];
                    if (nHits[rc] <= maxHits) {
                        for (unsigned i = 0; i < nHits[rc]; i++) {
                            // Get the true genome location of the read and mark the hit in its bucket.
                            unsigned location = hits[rc][i] - (rc ? (readLen[r] - seedLen - seedPos) : seedPos);
                            Bucket *bucket = getBucket(r, rc, location);
                            unsigned offset = location % BUCKET_SIZE;
                            bucket->found |= (1 << offset);
                            bucket->seedHits++;
                            if (wrapCount[r] == 0) {
                                bucket->disjointSeedHits++;
                            }
                            //TRACE("Adding seed hit at %d %d %u\n", r, rc, location, bucket);
                        }
                        if (wrapCount[r] == 0) {
                            disjointSeedsUsed[r][rc]++;
                        }
                    } else {
                        popularSeeds[r][rc]++;
                        totalPopularSeeds++;
                    }
                }
            } else {             
                TRACE("Seed contains N, so skipping it\n");
            }
            
            // Advance to the next seed.
            nextSeed[r] += seedLen;
            while (nextSeed[r] + seedLen >= readLen[r] && wrapCount[r] < seedLen) {
                wrapCount[r]++;
                nextSeed[r] = wrapOffset[wrapCount[r]];
            }
        }
        seedsTried += 1;
        
        // A good lower bound on best possible score is the min number of disjoint seeds used in
        // each direction, because any unseen location in that direction must have this many errors.
        int bestScoreOfAnyUnseenPair = min(
                min(lowerBound[1] + disjointSeedsUsed[0][0], lowerBound[0] + disjointSeedsUsed[1][1]),
                min(lowerBound[1] + disjointSeedsUsed[0][1], lowerBound[0] + disjointSeedsUsed[1][0]));
        
        bool finishNow = (seedsTried >= maxSeeds || bestScoreOfAnyUnseenPair > scoreLimit);
        
        // If we've tried a bunch more seeds or we should finish up, score some or all viable pairs.
        if (seedsTried == 10 || finishNow) {
            TRACE("Scoring time! seedsTried=%d, finishNow=%d, bestScoreOfAnyUnseenPair=%d\n",
                    seedsTried, finishNow, bestScoreOfAnyUnseenPair);
            
            // First, sort the bucket vectors to let us find viable pairs faster.
            for (int r = 0; r < 2; r++) {
                for (int rc = 0; rc < 2; rc++) {
                    sort(bucketLocations[r][rc].begin(), bucketLocations[r][rc].end());
#ifdef TRACE_PAIRED_ALIGNER
                    printf("Bucket locations %d %d:", r, rc);
                    for (int i = 0; i < bucketLocations[r][rc].size(); i++) {
                        printf(" %u", bucketLocations[r][rc][i]);
                    }
                    printf("\n");
#endif
                }
            }
            
            // Build a list of candidates in both orientations. Each candidate is a bucket along with a
            // number of seed hits (for it + its mate). We then sort the buckets by seed hits and sort
            // them, depending on the current scoreLimit.
            // Note that there are two orientations to try: read0 fwd + read1 RC or read1 fwd + read0 RC.
            // This is because the two ends of a read should normally match in opposite orientations.
            candidates.clear();
            unsigned pairsConsidered = 0;
            int bestPossibleScore = bestScoreOfAnyUnseenPair;
            for (int r0 = 0; r0 < 2; r0++) {
                for (int rc0 = 0; rc0 < 2; rc0++) {
                    int r1  = 1 - r0;
                    int rc1 = 1 - rc0;
                    FixedSizeVector<unsigned>& locs0 = bucketLocations[r0][rc0];
                    FixedSizeVector<unsigned>& locs1 = bucketLocations[r1][rc1];
                    int maxDist = maxSpacing + 2 * BUCKET_SIZE;
                    int start1 = 0;     // Position in locs1 to start searching for pairs at.
                    for (int pos0 = 0; pos0 < locs0.size(); pos0++) {
                        unsigned loc0 = locs0[pos0];
                        Bucket *b0 = getBucket(r0, rc0, loc0);
                        if (b0->allScored() && b0->mateStatus != UnknownAlignment)
                            continue;
                        b0->minPairScore = max(
                                (int) b0->minPairScore, 
                                max(lowerBound[r0], disjointSeedsUsed[r0][rc0] - b0->disjointSeedHits));
                        if (b0->minPairScore > scoreLimit)
                            continue;
                        while (start1 < locs1.size() && locs1[start1] + maxDist < loc0) {
                            start1++;
                        }
                        int bestSeedHits = b0->seedHits;
                        int bestMateDisjointSeedHits = 0;
                        for (int pos1 = start1; pos1 < locs1.size() && locs1[pos1] <= loc0 + maxDist; pos1++) {
                            unsigned loc1 = locs1[pos1];
                            Bucket *b1 = getBucket(r1, rc1, loc1);
                            bestSeedHits = max(bestSeedHits, b0->seedHits + b1->seedHits);
                            bestMateDisjointSeedHits = max(bestMateDisjointSeedHits, (int) b1->disjointSeedHits);
                            pairsConsidered++;
                        }
                        int bestMateScore = max(lowerBound[r1], disjointSeedsUsed[r1][rc1] - bestMateDisjointSeedHits);
                        b0->minPairScore = max(
                                (int) b0->minPairScore,
                                max(lowerBound[r0], disjointSeedsUsed[r0][rc0] - b0->disjointSeedHits) + bestMateScore);
                        bestPossibleScore = min(bestPossibleScore, (int) b0->minPairScore);
                        if (b0->minPairScore > scoreLimit)
                            continue;
                        candidates.push_back(Candidate(r0, rc0, loc0, b0, bestSeedHits));
                    }
                }
            }
            TRACE("Bucket pairs considered: %u\n", pairsConsidered);
            TRACE("Candidates found: %u\n", candidates.size());
            sort(candidates.begin(), candidates.end(), compareCandidates);

            for (int i = 0; i < candidates.size(); i++) {
                int r0 = candidates[i].read;
                int rc0 = candidates[i].isRC;
                int r1  = 1 - r0;
                int rc1 = 1 - rc0;
                unsigned loc0 = candidates[i].bucketLoc;
                Bucket *b0 = candidates[i].bucket;
                int maxDist = maxSpacing + 2 * BUCKET_SIZE;
            
                TRACE("Looking at bucket %d %d %u\n", r0, rc0, loc0);
                
                // Figure out which confDiff to use for this orientation.
                int realConfDiff = getConfDiff(seedsTried, popularSeeds, seedHits);

                candidatesTried++;
                if (candidatesTried > 60 && bestScore > (int)maxK) {
                    // We've tried a lot of candidates without finding anything remotely good.
                    // It's unlikely that trying more seeds will help.
                    TRACE("Breaking early because we tried a lot of candidates\n");
                    shouldStop = true;
                    break;
                }
                if (candidatesTried > 250 && bestScore <= (int)maxK && secondBestScore < bestScore + realConfDiff) {
                    // We've tried a lot of candidates without finding anything remotely good.
                    // It's unlikely that trying more seeds will help.
                    TRACE("Breaking early because we tried a lot of candidates\n");
                    shouldStop = true;
                    break;
                }
                if (candidatesTried > 700) {
                    shouldStop = true;
                    break;
                }
                
                if (b0->minPairScore > scoreLimit)
                    continue;

                int bestMateScore = b0->minPairScore - (disjointSeedsUsed[r0][rc0] - b0->disjointSeedHits);
                if (b0->mateStatus == SingleHit)
                    bestMateScore = b0->mateScore;
                scoreBucket(b0, r0, rc0 != 0, loc0, readData[r0][rc0], readLen[r0], scoreLimit - bestMateScore + 1);
                if ((int)b0->score > scoreLimit) {
                    continue;
                }
                
                scoreBucketMate(b0, r0, rc0 != 0, loc0, reads[r1], scoreLimit);
                if (b0->mateStatus == NotFound /*|| b0->mateStatus == MultipleHits*/)
                    continue;
                
                // If we got here, the bucket has a good score and its mate is also OK.
                loc0 += b0->bestOffset;
                unsigned loc1 = b0->mateLocation;
                int pairScore = b0->score + b0->mateScore;
                int multi1 = b0->mateStatus == MultipleHits;
                TRACE("Looking at pair %u-%d %u-%d [multi: %d]\n", loc0, rc0, loc1, rc1, multi1);
                TRACE("  Score: %d + %d = %d\n", b0->score, b0->mateScore, pairScore);
                
                // Check that this isn't the same pair with slightly different offsets.
                unsigned dist0 = distance(bestLoc[r0], loc0);
                unsigned dist1 = distance(bestLoc[r1], loc1);
                if (bestRC[r0] == rc0 && bestRC[r1] == rc1 &&
                    (dist0 <= maxK /*|| bestIsMulti[r0] && dist0 <= maxSpacing*/) &&
                    (dist1 <= maxK || ((bestIsMulti[r1] || multi1) && dist1 <= (unsigned)maxDist)))
                {
                    if (pairScore < bestScore) {
                        TRACE("  Updating score of best pair to %d\n", pairScore);
                        bestScore = pairScore;
                        bestScores[0] = b0->score;
                        bestScores[1] = b0->mateScore;
                        bestLoc[r0] = loc0; bestRC[r0] = rc0; bestIsMulti[r0] = 0;
                        bestLoc[r1] = loc1; bestRC[r1] = rc1; bestIsMulti[r1] = multi1;
                    }
                } else {
                    if (pairScore < bestScore) {
                        secondBestScore = bestScore;
                        bestScore = pairScore;
                        bestScores[0] = b0->score;
                        bestScores[1] = b0->mateScore;
                        bestLoc[r0] = loc0; bestRC[r0] = rc0; bestIsMulti[r0] = 0;
                        bestLoc[r1] = loc1; bestRC[r1] = rc1; bestIsMulti[r1] = multi1;
                    } else if (pairScore < secondBestScore) {
                        secondBestScore = pairScore;
                    }
                }
                
                // Return MultipleHits right away if bestScore and secondBestScore are both too small.
                if (bestScore < bestPossibleScore + realConfDiff && secondBestScore < bestScore + realConfDiff) {
                    TRACE("Returning 2 * MultipleHits because best and secondBest are both too small\n");
                    for (int i = 0; i < 2; i++) {
                        result->status[i] = MultipleHits;
                        result->location[i] = bestLoc[i];
                        result->isRC[i] = bestRC[i] != 0;
                        result->score[i] = bestScores[i];
                    }
                    return;
                }
                
                // Update scoreLimit.
                scoreLimit = min(bestScore, (int) maxK) + confDiff + 1;
                if (bestScore != INFINITE_SCORE && secondBestScore < bestScore + realConfDiff) {
                    // Since secondBestScore already means that our best location so far won't be a SingleHit,
                    // we really just care about finding better locations than that. However, still check for
                    // scores up to bestScore - 1 since those might become a new secondBestScore.
                    _ASSERT(bestScore >= 1); // Otherwise we would've returned MultipleHits
                    scoreLimit = bestScore - 1;
                }
                
                TRACE("  New best: %d, secondBest: %d, realConfDiff: %d, scoreLimit: %d\n", 
                    bestScore, secondBestScore, realConfDiff, scoreLimit);

                if (bestPossibleScore > scoreLimit) {
                    shouldStop = true;
                    break;
                }
            }
        } // End scoring section
    
        // Recompute bestScoreOfAnyUnseenPair. Now that we've scored all buckets found that have a
        // chance, the only way we can have an unseen location is if it hasn't been seen in *both*
        // of its orientations. Thus we add up the disjoint seeds for the two orientation combinations.
        bestScoreOfAnyUnseenPair = min(
                max(lowerBound[0], disjointSeedsUsed[0][0]) + max(lowerBound[1], disjointSeedsUsed[1][1]),
                max(lowerBound[0], disjointSeedsUsed[0][1]) + max(lowerBound[1], disjointSeedsUsed[1][0]));
        if (bestScoreOfAnyUnseenPair > scoreLimit) {
            break;
        }
        
    } // End loop over seeds
    
    int realConfDiff = getConfDiff(seedsTried, popularSeeds, seedHits);
    
    TRACE("End state: best=%d, secondBest=%d, realConfDiff=%d\n", bestScore, secondBestScore, realConfDiff);
    
    // Return based on the best and second-best scores so far.
    if ((unsigned)bestScore <= maxK) {
        if (bestScore + realConfDiff <= secondBestScore) {
            for (int i = 0; i < 2; i++) {
                result->status[i] = bestIsMulti[i] ? MultipleHits : SingleHit;
                result->location[i] = bestLoc[i];
                result->isRC[i] = bestRC[i] != 0;
                result->score[i] = bestScores[i];
            }
        } else {
            // There are multiple hits but we should still report the positions of a good one.
            for (int i = 0; i < 2; i++) {
                result->status[i] = MultipleHits;
                result->location[i] = bestLoc[i];
                result->isRC[i] = bestRC[i] != 0;
                result->score[i] = bestScores[i];
            }
        }
    } else {
        // We didn't find any viable pair, but let's see whether we have an alignment for at least one end.
        for (int r = 0; r < 2; r++) {
            unsigned best = INFINITE_SCORE;
            unsigned second = INFINITE_SCORE;
            unsigned bestLoc = 0xFFFFFFFF;
            bool bestRC;
            for (int rc = 0; rc < 2; rc++) {
                for (int i = 0; i < bucketLocations[r][rc].size(); i++) {
                    unsigned loc = bucketLocations[r][rc][i];
                    Bucket *bucket = getBucket(r, rc, loc);
                    if (bucket->allScored()) {
                        if (bucket->score < best) {
                            second = best;
                            best = bucket->score;
                            bestLoc = loc + bucket->bestOffset;
                            bestRC = rc != 0;
                        } else if (bucket->score < second) {
                            second = bucket->score;
                        }
                    }
                }
            }
            if (best <= 0.6 * maxK && best + realConfDiff + 1 <= second) {
                result->status[r] = SingleHit;
                result->location[r] = bestLoc;
                result->isRC[r] = bestRC;
                result->score[r] = best;
            } else {
                result->status[r] = NotFound;   
            }
        }
    }
    
    // TODO: return both multi if all seeds were too popular
}


void SmarterPairedEndAligner::scoreBucket(
    Bucket *bucket, int readId, bool isRC, unsigned location, const char *readData, int readLen, int scoreLimit)
{
    //_ASSERT(scoreLimit >= 0);
    if (scoreLimit < 0) {
        bucket->scored = bucket->found;
        bucket->score = INFINITE_SCORE;
        return;
    }
    //scoreLimit = min(scoreLimit, int(2 * maxK / 3));
    if (bucket->found != bucket->scored) {
        _uint64 unscored = bucket->found ^ bucket->scored;  // Bit vector
        unsigned long offset;
        // Keep looking for the first non-zero bit in unscored, and score that offset.
        while (_BitScanForward64(&offset, unscored)) {
            unscored ^= ((_int64)1 << offset);
            unsigned score = INFINITE_SCORE;
            const char *refData = index->getGenome()->getSubstring(location + offset, readLen + MAX_K);
            if (refData != NULL) {
                TRACE("  Genome: %.*s\n  Read:   %.*s\n", readLen, refData, readLen, readData);
                _uint64 cacheKey = ((_uint64) readId) << 33 | ((_uint64) isRC) << 32 | (location + offset);
                score = lv.computeEditDistance(refData, readLen + MAX_K, readData, readLen, scoreLimit, cacheKey);
                TRACE("  Called LV at %lu with limit %d: %d\n", location + offset, scoreLimit, score);
                if (score < 0) {
                    score = INFINITE_SCORE;
                }
            }
            if (score < bucket->score) {
                bucket->score = score;
                bucket->bestOffset = (unsigned short)offset;
            }
        }
        bucket->scored = bucket->found;
        TRACE("Updated score of bucket %u to %d (at %u)\n",
            location, bucket->score, location + bucket->bestOffset);
    }
}


void SmarterPairedEndAligner::scoreBucketMate(
    Bucket *bucket, int readId, bool isRC, unsigned location, Read *mate, int scoreLimit)
{
    _ASSERT(scoreLimit >= 0);
    if (bucket->mateStatus == UnknownAlignment) {
        //mateAligner->setMaxK(max(min(int(scoreLimit - bucket->score), int(2 * maxK / 3)), 0));
        mateAligner->setMaxK(max(int(scoreLimit - bucket->score), 0));
        mateAligner->setReadId(1 - readId);
        bool mateRC;
        bucket->mateStatus = mateAligner->AlignRead(
            mate, &bucket->mateLocation, &mateRC, &bucket->mateScore, maxSpacing + BUCKET_SIZE, location, !isRC);
        
/*
        // Search for the mate independently in two intervals, at minSpacing/maxSpacing before
        // and after the bucket itself.
        // TODO: If minSpacing is too small we should just search a single interval.
        int dist = (maxSpacing - minSpacing) / 2 + 2 * BUCKET_SIZE;
        AlignmentResult res1, res2;
        unsigned loc1, loc2;
        int score1, score2;
        bool rc1, rc2;
        mateAligner->setMaxK(max((int) (scoreLimit - bucket->score), 0));
        mateAligner->setReadId(1 - readId);
        res1 = mateAligner->AlignRead(mate, &loc1, &rc1, &score1, dist, location + BUCKET_SIZE - minSpacing - dist, !isRC);
        res2 = mateAligner->AlignRead(mate, &loc2, &rc2, &score2, dist, location - BUCKET_SIZE + minSpacing + dist, !isRC);
        TRACE("Mate statuses of bucket %u: %s %d at %u / %s %d at %u\n",
            location, AlignmentResultToString(res1), score1, loc1, AlignmentResultToString(res2), score2, loc2);
        if (res1 != NotFound && (res2 == NotFound || score2 >= score1 + (int)confDiff)) {
            bucket->mateStatus = res1;
            bucket->mateLocation = loc1;
            bucket->mateScore = score1;
        } else if (res2 != NotFound && (res1 == NotFound || score1 >= score2 + (int)confDiff)) {
            bucket->mateStatus = SingleHit;
            bucket->mateLocation = loc2;
            bucket->mateScore = score2;
        } else if (res1 != NotFound && res2 != NotFound) {
            bucket->mateStatus = MultipleHits;
            bucket->mateLocation = (score1 < score2 ? loc1 : loc2);
            bucket->mateScore = (score1 < score2 ? score1 : score2);
        } else {
            bucket->mateStatus = NotFound;
        }
*/
    }
}


void SmarterPairedEndAligner::computeRC(Read *read, char *outputBuf)
{
    const char *data = read->getData();
    int len = read->getDataLength();
    for (int i = 0; i < len; i++) {
        outputBuf[i] = complement[data[len-1-i]];
    }
}


void SmarterPairedEndAligner::clearState()
{
    bucketsUsed = 0;
    for (int r = 0; r < 2; r++) {
        for (int rc = 0; rc < 2; rc++) {
            bucketLocations[r][rc].clear();
            bucketTable[r][rc].clear();
        }
    }
    candidates.clear();
    lv.clearCache();
}


SmarterPairedEndAligner::Bucket *SmarterPairedEndAligner::newBucket()
{
    _ASSERT(bucketsUsed < MAX_BUCKETS);
    Bucket *b = &buckets[bucketsUsed++];
    b->found = 0;
    b->scored = 0;
    b->seedHits = 0;
    b->disjointSeedHits = 0;
    b->minPairScore = 0;
    b->score = INFINITE_SCORE;
    b->mateStatus = UnknownAlignment;
    return b;
}


SmarterPairedEndAligner::Bucket *SmarterPairedEndAligner::getBucket(int read, int rc, unsigned location)
{
    int key = location / BUCKET_SIZE;  // Important so that our hashes don't get confused
    Bucket *old = bucketTable[read][rc].get(key);
    if (old != NULL) {
        return old;
    } else {
        Bucket *b = newBucket();
        bucketTable[read][rc].put(key, b);
        bucketLocations[read][rc].push_back(key * BUCKET_SIZE);
        return b;
    }
}


int SmarterPairedEndAligner::getConfDiff(int seedsTried, int popularSeeds[2][2], int seedHits[2][2])
{
    unsigned total = max(popularSeeds[0][0] + popularSeeds[1][1], popularSeeds[1][0] + popularSeeds[0][1]);
    if (total > 2 * adaptiveConfDiffThreshold) {
        return confDiff + 2;
    } else if (total > adaptiveConfDiffThreshold) {
        return confDiff + 1;
    } else {
        return confDiff;
    }
}
