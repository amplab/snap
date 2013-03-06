/*++

Module Name:

    IntersectingPairedEndAligner.cpp

Abstract:

    A paired-end aligner based on set intersections to narrow down possible candidate locations.

Authors:

    Bill Bolosky, February, 2013

Environment:

    User mode service.

Revision History:

--*/

#include "stdafx.h"
#include "IntersectingPairedEndAligner.h"
#include "SeedSequencer.h"
#include "mapq.h"
#include "exit.h"

#ifdef  COMPILE_INTERSECTING

IntersectingPairedEndAligner::IntersectingPairedEndAligner(
        GenomeIndex  *index_,
        unsigned      maxReadSize_,
        unsigned      maxHits_,
        unsigned      maxK_,
        unsigned      maxSeeds_,
        unsigned      minSpacing_,                 // Minimum distance to allow between the two ends.
        unsigned      maxSpacing_,                 // Maximum distance to allow between the two ends.
        BigAllocator  *allocator) :
    index(index_), maxReadSize(maxReadSize_), maxHits(maxHits_), maxK(maxK_), maxSeeds(maxSeeds_), minSpacing(minSpacing_), maxSpacing(maxSpacing_),
    landauVishkin(NULL), reverseLandauVishkin(NULL), maxBigHits(500000), extraScoreLimit(5) /* should be a parameter */, maxMergeDistance(31) /*also should be a parameter*/
{
    allocateDynamicMemory(allocator, maxReadSize, maxHits, maxSeeds);

    rcTranslationTable['A'] = 'T';
    rcTranslationTable['G'] = 'C';
    rcTranslationTable['C'] = 'G';
    rcTranslationTable['T'] = 'A';
    rcTranslationTable['N'] = 'N';

    for (unsigned i = 0; i < 256; i++) {
        nTable[i] = 0;
    }

    nTable['N'] = 1;

    seedLen = index->getSeedLength();

    genome = index->getGenome();
    genomeSize = genome->getCountOfBases();

    distanceToSearchBeyondBestScore = 2;    // If probability goes down by about 1000x per mismatch, then going more than two steps away make an effect of 1 part per billion, which is too small to care about
}
    
IntersectingPairedEndAligner::~IntersectingPairedEndAligner()
{
    if (landauVishkin) {
        landauVishkin->~LandauVishkin<>();
    }

    if (reverseLandauVishkin) {
        reverseLandauVishkin->~LandauVishkin<-1>();
    }

    if (NULL != baseAligner) {
        baseAligner->~BaseAligner();
    }
}
    
    size_t 
IntersectingPairedEndAligner::getBigAllocatorReservation(GenomeIndex * index, unsigned maxHitsToConsider, unsigned maxReadSize, unsigned seedLen, unsigned maxSeedsToUse)
{
    CountingBigAllocator countingAllocator;
    {
        IntersectingPairedEndAligner aligner; // This has to be in a nested scope so it's destructor is called before that of the countingAllocator
        aligner.index = index;

        aligner.allocateDynamicMemory(&countingAllocator, maxReadSize, maxHitsToConsider, maxSeedsToUse);
        return sizeof(aligner) + countingAllocator.getMemoryUsed();
    }
}

    void
IntersectingPairedEndAligner::allocateDynamicMemory(BigAllocator *allocator, unsigned maxReadSize, unsigned maxHitsToConsider, unsigned maxSeedsToUse)
{
    seedUsed = (BYTE *) allocator->allocate(index->getSeedLength() + 7 / 8);

    landauVishkin = new(allocator) LandauVishkin<>();
    reverseLandauVishkin = new(allocator) LandauVishkin<-1>();

    for (unsigned whichRead = 0; whichRead < NUM_READS_PER_PAIR; whichRead++) {
        rcReadData[whichRead] = (char *)allocator->allocate(maxReadSize);
        rcReadQuality[whichRead] = (char *)allocator->allocate(maxReadSize);

        for (Direction dir = 0; dir < NUM_DIRECTIONS; dir++) {
            reversedRead[whichRead][dir] = (char *)allocator->allocate(maxReadSize);
            hashTableHitSets[whichRead][dir] = new HashTableHitSet(maxSeeds);
        }
    }

    for (int i = 0; i < NUM_SET_PAIRS; i++) {
        hitLocations[i] = new HitLocationRingBuffer(maxMergeDistance * 2 + 2);      // *2 is for before & after, +2 is just to make sure the ring buffer tail doesn't overrun the head
        mateHitLocations[i] = new HitLocationRingBuffer(2 * (maxSpacing + 1) + 2);  // Likewise.
    }

    baseAligner = new(allocator) BaseAligner(index, 1, maxHitsToConsider, maxK/2, maxReadSize, maxSeedsToUse, 4, landauVishkin, NULL, NULL, allocator);
}

    unsigned BJBMinLoc = 0xffffffff, BJBMaxLoc = 0;

    void 
IntersectingPairedEndAligner::align(
        Read                  *read0,
        Read                  *read1,
        PairedAlignmentResult *result)
{
    Read rcReads[NUM_READS_PER_PAIR];

    unsigned bestResultGenomeLocation[NUM_READS_PER_PAIR];
    Direction bestResultDirection[NUM_READS_PER_PAIR];
    unsigned bestResultScore[NUM_READS_PER_PAIR];
    unsigned popularSeedsSkipped[NUM_READS_PER_PAIR];

    reads[0][FORWARD] = read0;
    reads[1][FORWARD] = read1;

    //
    // Don't bother if one or both reads are too short.
    //
    if (read0->getDataLength() < 50 || read1->getDataLength() < 50) {
        alignWithBaseAligner(read0, read1, result, 70);
        return;
    }

    //
    // Build the RC reads.
    //
    unsigned countOfNs = 0;

    for (unsigned whichRead = 0; whichRead < NUM_READS_PER_PAIR; whichRead++) {
        Read *read = reads[whichRead][FORWARD];
        readLen[whichRead] = read->getDataLength();
        popularSeedsSkipped[whichRead] = 0;
        countOfHashTableLookups[whichRead] = 0;
        hitLocations[whichRead]->clear();
        mateHitLocations[whichRead]->clear();


        for (Direction dir = FORWARD; dir < NUM_DIRECTIONS; dir++) {
            totalHashTableHits[whichRead][dir] = 0;
            largestHashTableHit[whichRead][dir] = 0;
            hashTableHitSets[whichRead][dir]->init();
        }

        if (readLen[whichRead] > maxReadSize) {
            fprintf(stderr,"IntersectingPairedEndAligner:: got too big read (%d > %d)", readLen[whichRead], maxReadSize);
            soft_exit(1);
        }

        for (unsigned i = 0; i < reads[whichRead][FORWARD]->getDataLength(); i++) {
            rcReadData[whichRead][i] = rcTranslationTable[read->getData()[readLen[whichRead] - i - 1]];
            rcReadQuality[whichRead][i] = read->getQuality()[readLen[whichRead] - i - 1];
            countOfNs += nTable[read->getData()[i]];
        }
        reads[whichRead][RC] = &rcReads[whichRead];
        reads[whichRead][RC]->init(read->getId(), read->getIdLength(), rcReadData[whichRead], rcReadQuality[whichRead], read->getDataLength());
    }

    if (countOfNs > maxK) {
        alignWithBaseAligner(read0, read1, result, 70);
        return;
    }

    //
    // Build the reverse data for both reads in both directions for the backwards LV to use.
    //
    for (unsigned whichRead = 0; whichRead < NUM_READS_PER_PAIR; whichRead++) {
        for (Direction dir = 0; dir < NUM_DIRECTIONS; dir++) {
            Read *read = reads[whichRead][dir];

            for (unsigned i = 0; i < read->getDataLength(); i++) {
                reversedRead[whichRead][dir][i] = read->getData()[read->getDataLength() - i - 1];
            }
        }
    }

    unsigned thisPassSeedsNotSkipped[NUM_READS_PER_PAIR][NUM_DIRECTIONS] = {{0,0}, {0,0}}; 

    //
    // Initialize the member variables that are effectively stack locals, but are in the object
    // to aviod having to pass them to score.
    //
    double probabilityOfBestPair = 0;
    localBestPairProbability[0] = 0;
    localBestPairProbability[1] = 0;
    double probabilityOfAllPairs = 0;
    unsigned bestPairScore = 65536;
    unsigned scoreLimit = maxK + distanceToSearchBeyondBestScore;

    //
    // Phase 1: do the hash table lookups for each of the seeds for each of the reads and add them to the hit sets.
    //
    for (unsigned whichRead = 0; whichRead < NUM_READS_PER_PAIR; whichRead++) {
        unsigned nextSeedToTest = 0;
        unsigned wrapCount = 0;
        unsigned nPossibleSeeds = readLen[whichRead] - seedLen + 1;
        memset(seedUsed, 0, (__max(readLen[0], readLen[1]) + 7) / 8);
    
        while (countOfHashTableLookups[whichRead] < nPossibleSeeds && countOfHashTableLookups[whichRead] < maxSeeds) {
            if (nextSeedToTest >= nPossibleSeeds) {
                wrapCount++;
                if (wrapCount >= seedLen) {
                    //
                    // There aren't enough valid seeds in this read to reach our target.
                    //
                    break;
                }
                nextSeedToTest = GetWrappedNextSeedToTest(seedLen, wrapCount);
            }


            while (nextSeedToTest < nPossibleSeeds && IsSeedUsed(nextSeedToTest)) {
                //
                // This seed is already used.  Try the next one.
                //
                nextSeedToTest++;
            }

            if (nextSeedToTest >= nPossibleSeeds) {
                //
                // Unusable seeds have pushed us past the end of the read.  Go back around the outer loop so we wrap properly.
                //
                continue;
            }

            SetSeedUsed(nextSeedToTest);

            if (!Seed::DoesTextRepresentASeed(reads[whichRead][FORWARD]->getData() + nextSeedToTest, seedLen)) {
                //
                // It's got Ns in it, so just skip it.
                //
                nextSeedToTest++;
                continue;
            }

            Seed seed(reads[whichRead][FORWARD]->getData() + nextSeedToTest, seedLen);
            //
            // Find all instances of this seed in the genome.
            //
            unsigned nHits[NUM_DIRECTIONS];
            const unsigned *hits[NUM_DIRECTIONS];

            index->lookupSeed(seed, 0, genomeSize, &nHits[FORWARD], &hits[FORWARD], &nHits[RC], &hits[RC]);

            countOfHashTableLookups[whichRead]++;
            for (Direction dir = FORWARD; dir < NUM_DIRECTIONS; dir++) {
                unsigned offset;
                if (dir == FORWARD) {
                    offset = nextSeedToTest;
                } else {
                    offset = readLen[whichRead] - seedLen - nextSeedToTest;
                }
                if (nHits[dir] < maxBigHits) {
                    totalHashTableHits[whichRead][dir] += nHits[dir];
                    hashTableHitSets[whichRead][dir]->recordLookup(offset, nHits[dir], hits[dir]);
                } else {
                    popularSeedsSkipped[whichRead]++;
                }
            }

            nextSeedToTest += seedLen;
        } // while we need to lookup seeds for this read
    } // for each read

    readWithMoreHits = totalHashTableHits[0][FORWARD] + totalHashTableHits[0][RC] > totalHashTableHits[1][FORWARD] + totalHashTableHits[1][RC] ? 0 : 1;
    readWithFewerHits = 1 - readWithMoreHits;

    //
    // Phase 2: intersect the sets to find pairs, which we score.  The basic strategy is to look for seed hits in the read with fewer hits that also have
    // hits on the other read in the appropriate range.  When we find one, we score it and then if it has a low enough score then score the other
    // possibilities for the pair.  We interleave steps between read0 FORWARD/read1 RC and vice versa as a way of pushing score limit down more quickly.
    //

    scoreLimit = maxK + extraScoreLimit;
    struct IntersectionState {
        unsigned            lastSeedOffsetForReadWithFewerHits;
        unsigned            lastGenomeLocationForReadWithFewerHits;
        unsigned            lastGenomeLocationForReadWithMoreHits;
    } intersectionState[NUM_SET_PAIRS]; // One for each set pair
 
    bool setPairDone[NUM_SET_PAIRS] = {false, false};
    unsigned whichSetPairToCheck = 0;
    HashTableHitSet *setPair[NUM_SET_PAIRS][NUM_READS_PER_PAIR] = {{hashTableHitSets[0][FORWARD], hashTableHitSets[1][RC]}, {hashTableHitSets[0][RC], hashTableHitSets[1][FORWARD]}};
    Direction setPairDirection[NUM_SET_PAIRS][NUM_READS_PER_PAIR] = {{FORWARD, RC}, {RC, FORWARD}};

    //
    // Seed the intersection state by doing a first lookup for each pair.
    //
    for (unsigned i = 0; i < 2; i++) {
        setPairDone[i] = !setPair[i][readWithFewerHits]->getFirstHit(&intersectionState[i].lastGenomeLocationForReadWithFewerHits,
            &intersectionState[i].lastSeedOffsetForReadWithFewerHits);
        intersectionState[i].lastGenomeLocationForReadWithMoreHits = 0xffffffff;
    }

//    setPairDone[1] = true;  // Just do one while debugging, alternating's too confusing.
 
    if (setPairDone[0]) {
        whichSetPairToCheck = 1;
    }

    while (!(setPairDone[0] && setPairDone[1])) {
        //
        // Each iteration of this loop considers a single hit possibility on the read with fewer hits.  We've already looked it up,
        // but have not yet inserted it in the ring buffer.
        //
        unsigned smallReadHitGenomeLocation = intersectionState[whichSetPairToCheck].lastGenomeLocationForReadWithFewerHits;
if (smallReadHitGenomeLocation >= BJBMinLoc && smallReadHitGenomeLocation <= BJBMaxLoc) {
    printf("Here!\n");
}

        //
        // We get here after doing a lookup in the smaller read without checking the larger read for potential mate pairs.
        // Do that now.  We may already have one.  If not, then look.
        //
        if (intersectionState[whichSetPairToCheck].lastGenomeLocationForReadWithMoreHits > smallReadHitGenomeLocation + maxSpacing) {
            //
            // The last lookup for the big read is too high in the genome to be a mate pair for the new small read lookup.  Move it to the next hit that's in range.
            // We can dump the scored mate pairs wholesale here, too, since they're clearly all out of range now.
            //
            _ASSERT(mateHitLocations[whichSetPairToCheck]->isEmpty() || mateHitLocations[whichSetPairToCheck]->getHead()->genomeLocation > 
                         intersectionState[whichSetPairToCheck].lastGenomeLocationForReadWithFewerHits + maxSpacing);
            
            mateHitLocations[whichSetPairToCheck]->clear();

            unsigned seedOffset;
            if (!setPair[whichSetPairToCheck][readWithMoreHits]->getNextHitLessThanOrEqualTo(smallReadHitGenomeLocation + maxSpacing,
                                                                    &intersectionState[whichSetPairToCheck].lastGenomeLocationForReadWithMoreHits, &seedOffset)) {
                //
                // There's nothing left for the mate, and we know that there's no mates already found.  We're done with this set pair.
                //
                setPairDone[whichSetPairToCheck] = true;
                whichSetPairToCheck = 1 - whichSetPairToCheck;
                continue;
            }

            //
            // Insert this mate location into the scoring list (even if we don't need it for this small read hit, it might be necessary for a later one, and the
            // invariant is that the looked up mate locations are in the ring buffer).
            //
            mateHitLocations[whichSetPairToCheck]->insertHead(intersectionState[whichSetPairToCheck].lastGenomeLocationForReadWithMoreHits, seedOffset);
        } else {
            //
            // Get rid of any remembered mate hits that are too high to be possible mate pairs.
            //
            mateHitLocations[whichSetPairToCheck]->trimAboveLocation(smallReadHitGenomeLocation + maxSpacing);
        }

        unsigned minLocationToCheck;
        if (smallReadHitGenomeLocation >= maxSpacing) {
            minLocationToCheck = smallReadHitGenomeLocation - maxSpacing;
        } else {
            minLocationToCheck = 0;
        }

        // BUGBUG: don't let mate pairs cross chromosome boundaries.

        if (mateHitLocations[whichSetPairToCheck]->isEmpty() || mateHitLocations[whichSetPairToCheck]->getTail()->genomeLocation < minLocationToCheck) {
            //
            // There's no possible mate pair here.  Look for the next possible hit in the read with fewer hits, and then loop around.
            //
            unsigned moreHitsLocation;
            if (mateHitLocations[whichSetPairToCheck]->isEmpty()) {
                moreHitsLocation = intersectionState[whichSetPairToCheck].lastGenomeLocationForReadWithMoreHits;
            } else {
                moreHitsLocation = mateHitLocations[whichSetPairToCheck]->getTail()->genomeLocation;
            }
            if (!setPair[whichSetPairToCheck][readWithFewerHits]->getNextHitLessThanOrEqualTo(moreHitsLocation + maxSpacing,
                                                                   &intersectionState[whichSetPairToCheck].lastGenomeLocationForReadWithFewerHits,
                                                                   &intersectionState[whichSetPairToCheck].lastSeedOffsetForReadWithFewerHits)) {
                //
                // Out of candidates on the read with fewer hits.  Stop looking at this set pair.
                //
                setPairDone[whichSetPairToCheck] = true;
                whichSetPairToCheck = 1 - whichSetPairToCheck;
            } else {
                if (!setPairDone[1 - whichSetPairToCheck]) {
                    whichSetPairToCheck = 1 - whichSetPairToCheck;
                }
            }
            continue;
        }

        //
        // We've got two hits close enough to be mate pairs.  Score the hit on the read with fewer hits.
        //
        unsigned fewerHitScore;
        double fewerHitProbability;

        scoreLocation(readWithFewerHits, setPairDirection[whichSetPairToCheck][readWithFewerHits], smallReadHitGenomeLocation, intersectionState[whichSetPairToCheck].lastSeedOffsetForReadWithFewerHits, 
            scoreLimit, &fewerHitScore, &fewerHitProbability);

        if (-1 == fewerHitScore) {
            //
            // This location was too far off to be useful.  Skip to the next lower location, which might match the current set pair.
            //
            if (!setPair[whichSetPairToCheck][readWithFewerHits]->getNextLowerHit(&intersectionState[whichSetPairToCheck].lastGenomeLocationForReadWithFewerHits,
                &intersectionState[whichSetPairToCheck].lastSeedOffsetForReadWithFewerHits)) {
                //
                // Out of candidates on the read with fewer hits.  Stop looking at this set pair.
                //
                setPairDone[whichSetPairToCheck] = true;
                whichSetPairToCheck = 1 - whichSetPairToCheck;
            } else {
                if (!setPairDone[1 - whichSetPairToCheck]) {
                    whichSetPairToCheck = 1 - whichSetPairToCheck;
                }
            }
            continue;
        } // If the location had a -1 score

        hitLocations[whichSetPairToCheck]->trimAboveLocation(smallReadHitGenomeLocation + maxMergeDistance);
        hitLocations[whichSetPairToCheck]->insertHead(smallReadHitGenomeLocation, intersectionState[whichSetPairToCheck].lastSeedOffsetForReadWithFewerHits, fewerHitScore, fewerHitProbability);

        //
        // Add potential mate pairs to the ring buffer.
        //
        while (intersectionState[whichSetPairToCheck].lastGenomeLocationForReadWithMoreHits >= minLocationToCheck) {
            unsigned seedOffset;
            if (!setPair[whichSetPairToCheck][readWithMoreHits]->getNextLowerHit(&intersectionState[whichSetPairToCheck].lastGenomeLocationForReadWithMoreHits, &seedOffset)) {
                //
                // We're out of hits for the mate pair candidate, but we still need to score what we've got.
                //
                break;
            }

            mateHitLocations[whichSetPairToCheck]->insertHead(intersectionState[whichSetPairToCheck].lastGenomeLocationForReadWithMoreHits, seedOffset);
        } 

        //
        // Run through the saved mate hit locations.  For any that are in the appropriate range to be mate pairs, score them if they haven't already been scored, and
        // consider them as match candidates.
        //
        unsigned index;
        for (HitLocation *mateHitLocation = mateHitLocations[whichSetPairToCheck]->getTail(&index); NULL != mateHitLocation && mateHitLocation->genomeLocation >= minLocationToCheck; mateHitLocation = mateHitLocations[whichSetPairToCheck]->getNext(&index)) {
            unsigned largerGenomeOffset = __max(mateHitLocation->genomeLocation, smallReadHitGenomeLocation);
            unsigned smallerGenomeOffset = __min(mateHitLocation->genomeLocation, smallReadHitGenomeLocation);

            unsigned delta = largerGenomeOffset - smallerGenomeOffset;
            if (delta <= maxSpacing && delta >= minSpacing) {
                //
                // This is in a location that makes it a potential mate candidate.
                //
                if (!mateHitLocation->isScored) {
                    //
                    // It's not already scored, so score it.
                    //
                    scoreLocation(readWithMoreHits, setPairDirection[whichSetPairToCheck][readWithMoreHits], mateHitLocation->genomeLocation, mateHitLocation->seedOffset, scoreLimit, &mateHitLocation->score, &mateHitLocation->matchProbability);
                    mateHitLocation->isScored = true;
                }

                if (mateHitLocation->score != -1) {
                    double pairProbability = mateHitLocation->matchProbability * fewerHitProbability;
                    unsigned pairScore = mateHitLocation->score + fewerHitScore;
                    if (pairScore <= maxK && (pairScore < bestPairScore || pairScore == bestPairScore && pairProbability > probabilityOfBestPair)) {
                        //
                        // A new best hit.
                        //
                        bestPairScore = pairScore;
                        probabilityOfBestPair = pairProbability;
                        bestResultGenomeLocation[readWithFewerHits] = smallReadHitGenomeLocation;
                        bestResultGenomeLocation[readWithMoreHits] = mateHitLocation->genomeLocation;
                        bestResultScore[readWithFewerHits] = fewerHitScore;
                        bestResultScore[readWithMoreHits] = mateHitLocation->score;
                        bestResultDirection[readWithFewerHits] = setPairDirection[whichSetPairToCheck][readWithFewerHits];
                        bestResultDirection[readWithMoreHits] = setPairDirection[whichSetPairToCheck][readWithMoreHits];

                        scoreLimit = bestPairScore + distanceToSearchBeyondBestScore;
                    }

                    probabilityOfAllPairs += pairProbability;

                    if (probabilityOfAllPairs >= 4.9) {
                        //
                        // Nothing will rescue us from a 0 MAPQ, so just stop looking.
                        //
                        goto doneScoring;
                    }
                } // If the mate candidate had a non -1 score
            }
        } // for each potential mate candidate

        //
        // Advance the lookup in the smaller read, and then loop around checking the other set pair.
        //
        if (!setPair[whichSetPairToCheck][readWithFewerHits]->getNextLowerHit(&intersectionState[whichSetPairToCheck].lastGenomeLocationForReadWithFewerHits,
            &intersectionState[whichSetPairToCheck].lastSeedOffsetForReadWithFewerHits)) {
            //
            // Done with this set pair.
            //
            setPairDone[whichSetPairToCheck] = true;
        }

        if (!setPairDone[1 - whichSetPairToCheck]) {
            whichSetPairToCheck = 1 - whichSetPairToCheck;
        }
    } // loop alternating set pairs over 
     
doneScoring:

    if (bestPairScore == 65536) {
        //
        // Found nothing.
        //
        for (unsigned whichRead = 0; whichRead < NUM_READS_PER_PAIR; whichRead++) {
            result->location[whichRead] = 0xffffffff;
            result->mapq[whichRead] = 0;
            result->score[whichRead] = -1;
            result->status[whichRead] = NotFound;
        }
    } else {
        for (unsigned whichRead = 0; whichRead < NUM_READS_PER_PAIR; whichRead++) {
            result->location[whichRead] = bestResultGenomeLocation[whichRead];
            result->direction[whichRead] = bestResultDirection[whichRead];
            result->mapq[whichRead] = computeMAPQ(probabilityOfAllPairs, probabilityOfBestPair, bestResultScore[whichRead], popularSeedsSkipped[whichRead], true);
            result->status[whichRead] = result->mapq[whichRead] > 10 ? SingleHit : MultipleHits;
        }
    }
}

    void
IntersectingPairedEndAligner::alignWithBaseAligner(Read *read0, Read *read1, PairedAlignmentResult *result, int maxMapq)
{
    //
    // For whatever reason we can't align these reads singly.  Align them individually with the base aligner.
    //
    baseAligner->AlignRead(read0 ,&result->location[0], &result->direction[0], &result->score[0], &result->mapq[0]);
    baseAligner->AlignRead(read1 ,&result->location[1], &result->direction[1], &result->score[1], &result->mapq[1]);

    for (unsigned whichRead = 0; whichRead < NUM_READS_PER_PAIR; whichRead++) {
        result->mapq[whichRead] = __min(maxMapq, result->mapq[whichRead]);
    }
}

    void 
IntersectingPairedEndAligner::scoreLocation(
    unsigned             whichRead,
    Direction            direction,
    unsigned             genomeLocation,
    unsigned             seedOffset,
    unsigned             scoreLimit,
    unsigned            *score,
    double              *matchProbability)
{
    Read *readToScore = reads[whichRead][direction];
    unsigned readDataLength = readToScore->getDataLength();
    unsigned genomeDataLength = readDataLength + MAX_K; // Leave extra space in case the read has deletions
    const char *data = genome->getSubstring(genomeLocation, genomeDataLength);
    if (NULL == data) {
        //
        // We're up against the end of a chromosome.  Reduce the extra space enough that it isn't too
        // long.  We're willing to reduce it to less than the length of a read, because the read could
        // but up against the end of the chromosome and have insertions in it.
        //
        const Genome::Piece *piece = genome->getPieceAtLocation(genomeLocation);
                    
        unsigned endOffset;
        if (genomeLocation + readDataLength + MAX_K >= genome->getCountOfBases()) {
            endOffset = genome->getCountOfBases();
        } else {
            const Genome::Piece *nextPiece = genome->getPieceAtLocation(genomeLocation + readDataLength + MAX_K);
            _ASSERT(NULL != piece && piece->beginningOffset <= genomeLocation && piece != nextPiece);

            endOffset = nextPiece->beginningOffset;
        }
        genomeDataLength = endOffset - genomeLocation - 1;
        if (genomeDataLength >= readDataLength - MAX_K) {
            data = genome->getSubstring(genomeLocation, genomeDataLength);
            _ASSERT(NULL != data);
        }
    }

    if (NULL == data) {
        *score = -1;
        *matchProbability = 0;
        return;
    }
                
                
    // Compute the distance separately in the forward and backward directions from the seed, to allow
    // arbitrary offsets at both the start and end but not have to pay the cost of exploring all start
    // shifts in BoundedStringDistance
    double matchProb1, matchProb2;
    int score1, score2;
    // First, do the forward direction from where the seed aligns to past of it
    int readLen = readToScore->getDataLength();
    int seedLen = index->getSeedLength();
    int tailStart = seedOffset + seedLen;

    score1 = landauVishkin->computeEditDistance(data + tailStart, genomeDataLength - tailStart, readToScore->getData() + tailStart, readToScore->getQuality() + tailStart, readLen - tailStart,
        scoreLimit, &matchProb1);
    if (score1 == -1) {
        *score = -1;
    } else {
        // The tail of the read matched; now let's reverse the reference genome data and match the head
        int limitLeft = scoreLimit - score1;
        score2 = reverseLandauVishkin->computeEditDistance(data + seedOffset, seedOffset + MAX_K, reversedRead[whichRead][direction] + readLen - seedOffset,
                                                                    reads[whichRead][OppositeDirection(direction)]->getQuality() + readLen - seedOffset, seedOffset, limitLeft, &matchProb2);

        if (score2 == -1) {
            *score = -1;
        } else {
            *score = score1 + score2;
            // Map probabilities for substrings can be multiplied, but make sure to count seed too
            *matchProbability = matchProb1 * matchProb2 * pow(1 - SNP_PROB, seedLen);
        }
    }

    if (*score == -1) {
        *matchProbability = 0;
    }
}


 IntersectingPairedEndAligner::HashTableHitSet::HashTableHitSet(unsigned maxSeeds_) : maxSeeds(maxSeeds_), nLookupsUsed(0)
 {
      lookups = new HashTableLookup[maxSeeds];
 }
    void 
IntersectingPairedEndAligner::HashTableHitSet::init()
{
    nLookupsUsed = 0;
}
    void 
IntersectingPairedEndAligner::HashTableHitSet::recordLookup(unsigned seedOffset, unsigned nHits, const unsigned *hits)
{
    _ASSERT(nLookupsUsed < maxSeeds);
    if (nHits == 0) {
        //
        // Empty sets don't add anything, since we're really considering the union of all of them.
        //
        return;
    }
    lookups[nLookupsUsed].currentHitForIntersection = 0;
    lookups[nLookupsUsed].hits = hits;
    lookups[nLookupsUsed].nHits = nHits;
    lookups[nLookupsUsed].seedOffset = seedOffset;

    nLookupsUsed++;
}

    bool    
IntersectingPairedEndAligner::HashTableHitSet::getNextHitLessThanOrEqualTo(unsigned maxGenomeOffsetToFind, unsigned *actualGenomeOffsetFound, unsigned *seedOffsetFound)
{
    //
    // Iterate through all of the lookups and search for the best candidate.
    //
    bool anyFound = false;
    unsigned bestOffsetFound = 0;
    for (unsigned i = 0; i < nLookupsUsed; i++) {
        //
        // Binary search from the current starting offset to either the right place or the end.
        //
        int limit[2] = {(int)lookups[i].currentHitForIntersection, (int)lookups[i].nHits - 1};
        unsigned maxGenomeOffsetToFindThisSeed = maxGenomeOffsetToFind - lookups[i].seedOffset;
        while (limit[0] <= limit[1]) {
            unsigned probe = (limit[0] + limit[1]) / 2;
            //
            // Recall that the hit sets are sorted from largest to smallest, so the strange looking logic is actually right.
            //
            if (lookups[i].hits[probe] <= maxGenomeOffsetToFindThisSeed && (probe == 0 || lookups[i].hits[probe-1] > maxGenomeOffsetToFindThisSeed)) {
                anyFound = true;
                mostRecentLocationReturned = *actualGenomeOffsetFound = bestOffsetFound = __max(lookups[i].hits[probe] - lookups[i].seedOffset, bestOffsetFound);
                *seedOffsetFound = lookups[i].seedOffset;
                lookups[i].currentHitForIntersection = probe;
                break;
            }

            if (lookups[i].hits[probe] > maxGenomeOffsetToFindThisSeed) {   // Recode this without the if to avoid the hard-to-predict branch.
                limit[0] = probe + 1;
            } else {
                limit[1] = probe - 1;
            }
        } // While we're looking

        if (limit[0] > limit[1]) {
            lookups[i].currentHitForIntersection = lookups[i].nHits;    // We're done with this lookup.
        }
    } // For each lookup

    return anyFound;
}

    bool    
IntersectingPairedEndAligner::HashTableHitSet::getFirstHit(unsigned *genomeLocation, unsigned *seedOffsetFound)
{
    bool anyFound = false;
    *genomeLocation = 0;
    for (unsigned i = 0; i < nLookupsUsed; i++) {
        if (lookups[i].nHits > 0 && lookups[i].hits[0] - lookups[i].seedOffset > *genomeLocation) {
            mostRecentLocationReturned = *genomeLocation = lookups[i].hits[0] - lookups[i].seedOffset;
            *seedOffsetFound = lookups[i].seedOffset;
            anyFound = true;
        }
    }

    return anyFound;
}

    bool    
IntersectingPairedEndAligner::HashTableHitSet::getNextLowerHit(unsigned *genomeLocation, unsigned *seedOffsetFound)
{
    //
    // Look through all of the lookups and find the one with the highest location smaller than the
    //
    unsigned foundLocation = 0;
    bool anyFound = false;
 
    //
    // Run through the lookups pushing up any that are at the most recently returned 
    //
    for (unsigned i = 0; i < nLookupsUsed; i++) {
        _ASSERT(lookups[i].currentHitForIntersection == lookups[i].nHits || lookups[i].hits[lookups[i].currentHitForIntersection] - lookups[i].seedOffset <= mostRecentLocationReturned);
        if (lookups[i].currentHitForIntersection != lookups[i].nHits && lookups[i].hits[lookups[i].currentHitForIntersection] - lookups[i].seedOffset == mostRecentLocationReturned) {
            lookups[i].currentHitForIntersection++;
        }

        if (lookups[i].currentHitForIntersection != lookups[i].nHits) {
            *genomeLocation = foundLocation = __max(foundLocation, lookups[i].hits[lookups[i].currentHitForIntersection] - lookups[i].seedOffset);
            *seedOffsetFound = lookups[i].seedOffset;
            anyFound = true;
        }
    }
 
    if (anyFound) {
        mostRecentLocationReturned = foundLocation;
    }
    return anyFound;
}

#endif  // COMPILE_INTERSECTING
