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


SmarterPairedEndAligner::SmarterPairedEndAligner(
         GenomeIndex  *index_,
         unsigned      maxReadSize_,
         unsigned      confDiff_,
         unsigned      maxHits_,
         unsigned      maxK_,
         unsigned      maxSeeds_,
         unsigned      minSpacing_,
         unsigned      maxSpacing_,
         unsigned      adaptiveConfDiffThreshold_,
         bool          skipAlignTogether_)
 :   index(index_), maxReadSize(maxReadSize_), confDiff(confDiff_), maxHits(maxHits_),
     maxK(maxK_), maxSeeds(maxSeeds_), minSpacing(minSpacing_), maxSpacing(maxSpacing_),
     adaptiveConfDiffThreshold(adaptiveConfDiffThreshold_), seedLen(index_->getSeedLength()),
     maxBuckets(FirstPowerOf2GreaterThanOrEqualTo(maxSeeds_ * maxHits_ * 4)), lv(2 * maxBuckets),
     skipAlignTogether(skipAlignTogether_)
{
    // Initialize the bucket data structures.
    buckets = new Bucket[maxBuckets];
    for (int r = 0; r < 2; r++) {
        for (int rc = 0; rc < 2; rc++) {
            bucketLocations[r][rc].reserve(maxBuckets);
            bucketTable[r][rc].reserve(2 * maxBuckets);
        }
    }
    candidates.reserve(maxBuckets);

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

    mateAligner = new BaseAligner(index, confDiff, maxHits, maxK, maxReadSize, maxSeeds, adaptiveConfDiffThreshold, &lv);

    //
    // Allocate these all in one big allocation.
    //
    reversedRead[FORWARD][0] = (char *)BigAlloc(maxReadSize * 6);   // 4 for each read and forward/rc + 2 for reverse quality (RC quality == reverse quality)
    reversedRead[RC][0] = reversedRead[FORWARD][0] + maxReadSize;
    reversedRead[FORWARD][1] = reversedRead[FORWARD][0] + 2 * maxReadSize;
    reversedRead[RC][1] = reversedRead[FORWARD][0] + 3 * maxReadSize;
    reversedQuality[0] = reversedRead[FORWARD][0] + 4 * maxReadSize;
    reversedQuality[1] = reversedRead[FORWARD][0] + 5 * maxReadSize;
}


SmarterPairedEndAligner::~SmarterPairedEndAligner()
{
    delete[] buckets;
    BigDealloc(reversedRead[FORWARD][0]);
}


void SmarterPairedEndAligner::align(Read *read0, Read *read1, PairedAlignmentResult *result)
{
    Read *reads[NUM_READS_PER_PAIR] = {read0, read1};
 
    result->status[0] = NotFound;
    result->status[1] = NotFound;

    result->location[0] = 0xFFFFFFFF;
    result->location[1] = 0xFFFFFFFF;

    result->fromAlignTogether = true;
    result->nanosInAlignTogether = 0;

    TRACE("Aligning read pair %.*s:\n%.*s\n%.*s\n",
            reads[0]->getIdLength(), reads[0]->getId(),
            reads[0]->getDataLength(), reads[0]->getData(),
            reads[1]->getDataLength(), reads[1]->getData());

    if (read0->getDataLength() < 50 && read1->getDataLength() < 50) {
        TRACE("Reads are both too short -- returning");
        return;
    }

    int lowerBound[2] = {maxK, maxK};
    alignTogether(reads, result, lowerBound);
}

void SmarterPairedEndAligner::alignTogether(Read *reads[NUM_READS_PER_PAIR], PairedAlignmentResult *result, int lowerBound[NUM_READS_PER_PAIR])
{
    int readLen[2];
    const char *readData[2][2];     // [read][direction]
    const char *quality[2];         // [read]
    char rcData[2][MAX_READ_SIZE];

    clearState();

    for (int r = 0; r < 2; r++) {
        readLen[r] = reads[r]->getDataLength();
        readData[r][0] = reads[r]->getData();
        quality[r] = reads[r]->getQuality();
        computeRC(reads[r], &rcData[r][0]);
        readData[r][1] = &rcData[r][0];
    }
    
    int bestScore = INFINITE_SCORE;
    double probabilityOfBestPair = 0;
    double probabilityOfAllPairs = 0;
    double probabilityOfAllSingles[NUM_READS_PER_PAIR] = {0, 0};
    int secondBestScore = INFINITE_SCORE;

    // Locations of best scoring pair
    int bestLoc[NUM_READS_PER_PAIR] = {0xFFFFFFFF, 0xFFFFFFFF};
    Direction bestDirection[NUM_READS_PER_PAIR];
    int bestIsMulti[NUM_READS_PER_PAIR] = {0, 0};
    int bestScores[NUM_READS_PER_PAIR] = {INFINITE_SCORE, INFINITE_SCORE};
    int bestMapqs[NUM_READS_PER_PAIR] = {0, 0};
  
    //
    // Probabilities for everything measured
    //
    double allProbabilities[NUM_READS_PER_PAIR] = {0.0, 0.0};
    
    int disjointSeedsUsed[NUM_READS_PER_PAIR][NUM_DIRECTIONS] = {{0, 0}, {0, 0}};   // [read][direction]
    int popularSeeds[NUM_READS_PER_PAIR][NUM_DIRECTIONS] = {{0, 0}, {0, 0}};        // [read][direction]
    int seedHits[NUM_READS_PER_PAIR][NUM_DIRECTIONS] = {{0, 0}, {0, 0}};            // [read][direction]
    int totalPopularSeeds = 0;
    
    int nextSeed[NUM_READS_PER_PAIR] = {0, 0};        // Which seed to try next on each read
    int wrapCount[NUM_READS_PER_PAIR] = {0, 0};       // How many times we've wraped around on each read
    int firstPassNotSkippedSeeds[NUM_READS_PER_PAIR] = {0, 0};
    unsigned smallestSkippedSeed[NUM_READS_PER_PAIR] = {0xffffffff, 0xffffffff};

    int scoreLimit = maxK + confDiff + 1;

    unsigned candidatesTried = 0;

    //
    // Compute the reversed data.
    //
    for (int r = 0; r < NUM_READS_PER_PAIR; r++) {
        for (Direction dir = FORWARD; dir < NUM_DIRECTIONS; dir++) {
            for (int i = 0; i < readLen[r]; i++) {
                reversedRead[r][dir][i] = readData[r][dir][readLen[r] - i - 1];
            }
        }
        for (int i = 0; i < readLen[r]; i++) {
            reversedQuality[r][i] = quality[r][readLen[r] - i - 1];
        }
    }

    unsigned seedsTried = 0;
    bool shouldStop = false;
    while (seedsTried < maxSeeds && !shouldStop) {
        // Try one seed from each read just to keep things simple.
        for (int r = 0; r < NUM_READS_PER_PAIR; r++) {
            // TODO: We could stop searching for one of the reads early if we rule it out via confDiff?
            if (wrapCount[r] >= seedLen) {
                continue;           // We've tried all possible seeds for this read
            }
            
            // Look up the current seed.
            int seedPos = nextSeed[r];
            TRACE("Trying read %d seed %d: %.*s\n", r, seedPos, seedLen, readData[r][0] + seedPos);
            if (Seed::DoesTextRepresentASeed(readData[r][FORWARD] + seedPos, seedLen)) {
                Seed seed(readData[r][FORWARD] + seedPos, seedLen);
                unsigned nHits[NUM_DIRECTIONS];          // Number of hits in forward and RC directions
                const unsigned *hits[NUM_DIRECTIONS];    // Actual hit locations in each direction

                //
                // Find hits for this seed (in both Directions) in the big hash table.
                //
                index->lookupSeed(seed, &nHits[FORWARD], &hits[FORWARD], &nHits[RC], &hits[RC]);
                
#ifdef TRACE_PAIRED_ALIGNER
                printf("  %u fwd hits:", nHits[0]);
                //for (int i = 0; i < nHits[0]; i++) { printf(" %u", hits[0][i]); }
                printf("\n");
                printf("  %u RC hits:", nHits[1]);
                //for (int i = 0; i < nHits[1]; i++) { printf(" %u", hits[1][i]); }
                printf("\n");
#endif
            
                // Add in all the hits.
                for (Direction direction = 0; direction < NUM_DIRECTIONS; direction++) {
                    seedHits[r][direction] += nHits[direction];
                    if (nHits[direction] <= maxHits) {
                        if (wrapCount[direction] == 0) {
                            firstPassNotSkippedSeeds[direction]++;
                        }
                        for (unsigned i = 0; i < nHits[direction]; i++) {
                            // Get the true genome location of the read and mark the hit in its bucket.
                            unsigned location = hits[direction][i] - (direction == RC ? (readLen[r] - seedLen - seedPos) : seedPos);
                            Bucket *bucket = getBucket(r, direction, location);
                            unsigned offset = location % BUCKET_SIZE;
                            bucket->found |= (1 << offset);
                            bucket->seedHits++;
                            if (wrapCount[r] == 0) {
                                bucket->disjointSeedHits++;
                            }
                            //TRACE("Adding seed hit at %d %d %u\n", r, rc, location, bucket);
                        }
                        if (wrapCount[r] == 0) {
                            disjointSeedsUsed[r][direction]++;
                        }
                    } else {
                        smallestSkippedSeed[direction] = __min(smallestSkippedSeed[direction], nHits[direction]);
                        popularSeeds[r][direction]++;
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
        } // for each read
        seedsTried += 1;
        
        // A good lower bound on best possible score is the min number of disjoint seeds used in
        // each direction, because any unseen location in that direction must have this many errors.
        int bestScoreOfAnyUnseenPair = min(
                min(lowerBound[1] + disjointSeedsUsed[0][FORWARD], lowerBound[0] + disjointSeedsUsed[1][RC]),
                min(lowerBound[1] + disjointSeedsUsed[0][RC], lowerBound[0] + disjointSeedsUsed[1][FORWARD]));
        
        bool finishNow = (seedsTried >= maxSeeds || bestScoreOfAnyUnseenPair > scoreLimit);
        
        // If we've tried a bunch more seeds or we should finish up, score some or all viable pairs.
        if (seedsTried == 10 || finishNow) {
            TRACE("Scoring time! seedsTried=%d, finishNow=%d, bestScoreOfAnyUnseenPair=%d\n",
                    seedsTried, finishNow, bestScoreOfAnyUnseenPair);
            
            // First, sort the bucket vectors to let us find viable pairs faster.
            for (int r = 0; r < NUM_READS_PER_PAIR; r++) {
                for (Direction direction = 0; direction < 2; direction++) {
                    sort(bucketLocations[r][direction].begin(), bucketLocations[r][direction].end());
#ifdef TRACE_PAIRED_ALIGNER
                    printf("Bucket locations %d %d:", r, rc);
                    for (int i = 0; i < bucketLocations[r][direction].size(); i++) {
                        printf(" %u", bucketLocations[r][direction][i]);
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
            for (int r0 = 0; r0 < NUM_READS_PER_PAIR; r0++) {
                for (Direction direction0 = 0; direction0 < 2; direction0++) {
                    int r1  = 1 - r0;
                    Direction direction1 = OppositeDirection(direction0);
                    FixedSizeVector<unsigned>& locs0 = bucketLocations[r0][direction0];
                    FixedSizeVector<unsigned>& locs1 = bucketLocations[r1][direction1];
                    int maxDist = maxSpacing + 2 * BUCKET_SIZE;
                    int start1 = 0;     // Position in locs1 to start searching for pairs at.
                    for (int pos0 = 0; pos0 < locs0.size(); pos0++) {
                        unsigned loc0 = locs0[pos0];
                        Bucket *b0 = getBucket(r0, direction0, loc0);
                        if (b0->allScored() && b0->mateStatus != UnknownAlignment)
                            continue;
                        b0->minPairScore = max(
                                (int) b0->minPairScore, 
                                max(lowerBound[r0], disjointSeedsUsed[r0][direction0] - b0->disjointSeedHits));
                        if (b0->minPairScore > scoreLimit)
                            continue;
                        while (start1 < locs1.size() && locs1[start1] + maxDist < loc0) {
                            start1++;
                        }
                        int bestSeedHits = b0->seedHits;
                        int bestMateDisjointSeedHits = 0;
                        for (int pos1 = start1; pos1 < locs1.size() && locs1[pos1] <= loc0 + maxDist; pos1++) {
                            unsigned loc1 = locs1[pos1];
                            Bucket *b1 = getBucket(r1, direction1, loc1);
                            bestSeedHits = max(bestSeedHits, b0->seedHits + b1->seedHits);
                            bestMateDisjointSeedHits = max(bestMateDisjointSeedHits, (int) b1->disjointSeedHits);
                            pairsConsidered++;
                        }
                        int bestMateScore = max(lowerBound[r1], disjointSeedsUsed[r1][direction1] - bestMateDisjointSeedHits);
                        b0->minPairScore = max(
                                (int) b0->minPairScore,
                                max(lowerBound[r0], disjointSeedsUsed[r0][direction0] - b0->disjointSeedHits) + bestMateScore);
                        bestPossibleScore = min(bestPossibleScore, (int) b0->minPairScore);
                        if (b0->minPairScore > scoreLimit)
                            continue;
                        candidates.push_back(Candidate(r0, direction0, loc0, b0, bestSeedHits));
                    }
                }
            }
            TRACE("Bucket pairs considered: %u\n", pairsConsidered);
            TRACE("Candidates found: %u\n", candidates.size());
            sort(candidates.begin(), candidates.end(), compareCandidates);

            for (int i = 0; i < candidates.size(); i++) {
                int r0 = candidates[i].read;
                Direction direction0 = candidates[i].direction;
                int r1  = 1 - r0;
                Direction direction1 = OppositeDirection(direction0);
                unsigned loc0 = candidates[i].bucketLoc;
                Bucket *b0 = candidates[i].bucket;
                int maxDist = maxSpacing + 2 * BUCKET_SIZE;
            
                TRACE("Looking at bucket %d %d %u\n", r0, direction0, loc0);
                
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

                int bestMateScore = b0->minPairScore - (disjointSeedsUsed[r0][direction0] - b0->disjointSeedHits);
                if (b0->mateStatus == SingleHit)
                    bestMateScore = b0->mateScore;
                double matchProbability = 0.0;

                scoreBucket(b0, r0, direction0, loc0, 0, 0, readData[r0][direction0], reversedRead[r0][direction0], quality[r0], reversedQuality[r0], readLen[r0], scoreLimit - bestMateScore + 1, &matchProbability);
                _ASSERT(matchProbability >= 0);
                probabilityOfAllSingles[r0] += matchProbability;

                if ((int)b0->score > scoreLimit) {
                    continue;
                }
                
                int mateMapq;
                scoreBucketMate(b0, r0, direction0, loc0, reads[r1], scoreLimit, &mateMapq);
                if (b0->mateStatus == NotFound /*|| b0->mateStatus == MultipleHits*/)
                    continue;
                
                //
                // The match probability is the product of the probabilities of each end, which we get for the
                // mate by using its mapq (which was computed only within the region).
                //
                double pairMatchProbability = matchProbability * b0->mateProbability;
                probabilityOfAllPairs += pairMatchProbability;  

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
                if (bestDirection[r0] == direction0 && bestDirection[r1] == direction1 &&
                    (dist0 <= maxK /*|| bestIsMulti[r0] && dist0 <= maxSpacing*/) &&
                    (dist1 <= maxK || ((bestIsMulti[r1] || multi1) && dist1 <= (unsigned)maxDist)))
                {
                    if (pairScore < bestScore || pairScore == bestScore && pairMatchProbability > probabilityOfBestPair) {
                        TRACE("  Updating score of best pair to %d\n", pairScore);
                        bestScore = pairScore;
                        probabilityOfBestPair = pairMatchProbability;
                        bestScores[0] = b0->score;
                        bestScores[1] = b0->mateScore;
                        bestLoc[r0] = loc0; bestDirection[r0] = direction0; bestIsMulti[r0] = 0;
                        bestLoc[r1] = loc1; bestDirection[r1] = direction1; bestIsMulti[r1] = multi1;
                    }
                } else {
                    if (pairScore < bestScore || pairScore == bestScore && pairMatchProbability > probabilityOfBestPair) {
                        secondBestScore = bestScore;
                        bestScore = pairScore;
                        probabilityOfBestPair = pairMatchProbability;
                        bestScores[0] = b0->score;
                        bestScores[1] = b0->mateScore;
                        bestLoc[r0] = loc0; bestDirection[r0] = direction0; bestIsMulti[r0] = 0;
                        bestLoc[r1] = loc1; bestDirection[r1] = direction1; bestIsMulti[r1] = multi1;
                    } else if (pairScore < secondBestScore) {
                        secondBestScore = pairScore;
                    }
                }
                
                // Update scoreLimit.
                scoreLimit = min(bestScore + 5, (int) maxK) + confDiff + 1;
#if     0 // Off because we want to compute probabilities for plausible pairs so that we can get a decent mapq estimate
                if (bestScore != INFINITE_SCORE && secondBestScore < bestScore + realConfDiff) {
                    // Since secondBestScore already means that our best location so far won't be a SingleHit,
                    // we really just care about finding better locations than that. However, still check for
                    // scores up to bestScore - 1 since those might become a new secondBestScore.
                    _ASSERT(bestScore >= 1); // Otherwise we would've returned MultipleHits
                    scoreLimit = bestScore - 1;
                }
#endif  // 0
                
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
                max(lowerBound[0], disjointSeedsUsed[0][FORWARD]) + max(lowerBound[1], disjointSeedsUsed[1][RC]),
                max(lowerBound[0], disjointSeedsUsed[0][RC]) + max(lowerBound[1], disjointSeedsUsed[1][FORWARD]));
        if (bestScoreOfAnyUnseenPair > scoreLimit) {
            break;
        }
        
    } // End loop over seeds
    
    int realConfDiff = getConfDiff(seedsTried, popularSeeds, seedHits);
    
    TRACE("End state: best=%d, secondBest=%d, realConfDiff=%d\n", bestScore, secondBestScore, realConfDiff);
    
    // Return based on the best and second-best scores so far.
    if ((unsigned)bestScore <= maxK) {
        result->mapq[0] = result->mapq[1] = computeMAPQ(probabilityOfAllPairs, probabilityOfBestPair, bestScore, totalPopularSeeds, true);

        if (bestScore + realConfDiff <= secondBestScore) {
            for (int i = 0; i < NUM_READS_PER_PAIR; i++) {
                result->status[i] = bestIsMulti[i] ? MultipleHits : SingleHit;
                result->location[i] = bestLoc[i];
                result->direction[i] = bestDirection[i];
                result->score[i] = bestScores[i];
            }
        } else {
            // There are multiple hits but we should still report the positions of a good one.
            for (int i = 0; i < NUM_READS_PER_PAIR; i++) {
                result->status[i] = MultipleHits;
                result->location[i] = bestLoc[i];
                result->direction[i] = bestDirection[i];
                result->score[i] = bestScores[i];
            }
        }
    } else {
        // We didn't find any viable pair, but let's see whether we have an alignment for at least one end.
        for (int r = 0; r < NUM_READS_PER_PAIR; r++) {
            unsigned best = INFINITE_SCORE;
            unsigned second = INFINITE_SCORE;
            unsigned bestLoc = 0xFFFFFFFF;
            Direction bestDirection;
            int bestMapq = 0;
            for (Direction direction = 0; direction < 2; direction++) {
                for (int i = 0; i < bucketLocations[r][direction].size(); i++) {
                    unsigned loc = bucketLocations[r][direction][i];
                    Bucket *bucket = getBucket(r, direction, loc);
                    if (bucket->allScored()) {
                        if (bucket->score < best || bucket->score == best && bucket->matchProbability > mapqToProbability(bestMapq)) {
                            second = best;
                            best = bucket->score;
                            bestLoc = loc + bucket->bestOffset;
                            bestDirection = direction;
                            bestMapq = computeMAPQ(probabilityOfAllSingles[r], bucket->matchProbability, bestScore, totalPopularSeeds, true);
                        } else if (bucket->score < second) {
                            second = bucket->score;
                        }
                    }
                }
            }
            if (best <= 0.6 * maxK && best + realConfDiff + 1 <= second) {
                result->status[r] = (best + realConfDiff + 1 <= second) ? SingleHit : MultipleHits;
                result->location[r] = bestLoc;
                result->direction[r] = bestDirection;
                result->score[r] = best;
                result->mapq[r] = bestMapq;
            } else {
                result->status[r] = NotFound;   
            }
        }
    }
    
    // TODO: return both multi if all seeds were too popular
}


void SmarterPairedEndAligner::scoreBucket(
        Bucket *            bucket, 
        int                 readId, 
        Direction           direction, 
        unsigned            location, 
        int                 seedOffset,
        int                 seedLength,
        const char *        readData, 
        const char *        reverseReadData,
        const char *        qualityString, 
        const char *        reverseQualityString,
        int                 readLen, 
        int                 scoreLimit, 
        double *            matchProbability)
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
            int genomeDataLength = readLen + MAX_K; // Leave extra space in case the read has deletions 
            const char *refData = index->getGenome()->getSubstring(location + offset, readLen + MAX_K);
            double localMatchProbability = 0.0;

            if (NULL == refData) {
                //
                // We're up against the end of a chromosome.  Reduce the extra space enough that it isn't too
                // long.  We're willing to reduce it to less than the length of a read, because the read could
                // but up against the end of the chromosome and have insertions in it.
                //
                const Genome::Piece *piece = index->getGenome()->getPieceAtLocation(location + offset);
                const Genome::Piece *nextPiece = index->getGenome()->getPieceAtLocation(location + offset + readLen + MAX_K);
                _ASSERT(NULL != piece && piece->beginningOffset <= location + offset && piece != nextPiece);
                unsigned endOffset;
                if (NULL != nextPiece) {
                    endOffset = nextPiece->beginningOffset;
                } else {
                    endOffset = index->getGenome()->getCountOfBases();
                }
                genomeDataLength = endOffset - (location + offset) - 1;
                if (genomeDataLength >= readLen - MAX_K) {
                    refData = index->getGenome()->getSubstring(location + offset, genomeDataLength);
                    _ASSERT(NULL != refData);
                }
            }

            if (refData != NULL) {
                TRACE("  Genome: %.*s\n  Read:   %.*s\n", readLen, refData, readLen, readData);
                _uint64 cacheKey = ((_uint64) readId) << 33 | ((_uint64) direction) << 32 | (location + offset);

                TRACE("  Called LV at %lu with limit %d: %d\n", location + offset, scoreLimit, score);
                // Compute the distance separately in the forward and backward directions from the seed, to allow
                // arbitrary offsets at both the start and end but not have to pay the cost of exploring all start
                // shifts in BoundedStringDistance
                double matchProb1, matchProb2;
                int score1, score2;
                // First, do the forward direction from where the seed aligns to past of it

                int tailStart = seedOffset + seedLen;

                score1 = lv.computeEditDistance(refData + tailStart, genomeDataLength - tailStart, readData + tailStart, qualityString + tailStart, readLen - tailStart,
                    scoreLimit, &matchProb1);
                if (score1 == -1) {
                    score = -1;
                } else {
                    // The tail of the read matched; now let's reverse the reference genome data and match the head
                    int limitLeft = scoreLimit - score1;
                    score2 = reverseLV.computeEditDistance(refData + seedOffset, seedOffset + MAX_K, reverseReadData + readLen - seedOffset,
                                                                                reverseQualityString + readLen - seedOffset, seedOffset, limitLeft, &matchProb2);
                    if (score2 == -1) {
                        score = -1;
                    } else {
                        score = score1 + score2;
                        // Map probabilities for substrings can be multiplied, but make sure to count seed too
                        localMatchProbability = matchProb1 * matchProb2 * pow(1 - SNP_PROB, seedLen);
                    }
                }

                if (score == -1) {
                    localMatchProbability = 0;
                    score = INFINITE_SCORE;
                }
            }

            if (score < bucket->score || score == bucket->score && localMatchProbability > bucket->matchProbability) {
                bucket->score = score;
                bucket->matchProbability = localMatchProbability;
                bucket->bestOffset = (unsigned short)offset;
            }
        }
        bucket->scored = bucket->found;
        TRACE("Updated score of bucket %u to %d (at %u)\n",
            location, bucket->score, location + bucket->bestOffset);
        *matchProbability = bucket->matchProbability;
    }
}


void SmarterPairedEndAligner::scoreBucketMate(
    Bucket *bucket, int readId, Direction direction, unsigned location, Read *mate, int scoreLimit, int *mateMapq)
{
    _ASSERT(scoreLimit >= 0);
    if (bucket->mateStatus == UnknownAlignment) {
        //mateAligner->setMaxK(max(min(int(scoreLimit - bucket->score), int(2 * maxK / 3)), 0));
        mateAligner->setMaxK(max(int(scoreLimit - bucket->score), 0));
        mateAligner->setReadId(1 - readId);
        Direction mateDirection;
        bucket->mateStatus = mateAligner->AlignRead(
            mate, &bucket->mateLocation, &mateDirection, &bucket->mateScore, mateMapq, maxSpacing + BUCKET_SIZE, location, OppositeDirection(direction), &bucket->mateProbability, NULL);
        
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
        for (Direction direction = 0; direction < 2; direction++) {
            bucketLocations[r][direction].clear();
            bucketTable[r][direction].clear();
        }
    }
    candidates.clear();
    lv.clearCache();
}


SmarterPairedEndAligner::Bucket *SmarterPairedEndAligner::newBucket()
{
    _ASSERT(bucketsUsed < maxBuckets);
    Bucket *b = &buckets[bucketsUsed++];
    b->found = 0;
    b->scored = 0;
    b->seedHits = 0;
    b->disjointSeedHits = 0;
    b->minPairScore = 0;
    b->score = INFINITE_SCORE;
    b->mateStatus = UnknownAlignment;
    b->matchProbability = 0;
    return b;
}


SmarterPairedEndAligner::Bucket *SmarterPairedEndAligner::getBucket(int read, Direction direction, unsigned location)
{
    int key = location / BUCKET_SIZE;  // Important so that our hashes don't get confused
    Bucket *old = bucketTable[read][direction].get(key);
    if (old != NULL) {
        return old;
    } else {
        Bucket *b = newBucket();
        bucketTable[read][direction].put(key, b);
        bucketLocations[read][direction].push_back(key * BUCKET_SIZE);
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
