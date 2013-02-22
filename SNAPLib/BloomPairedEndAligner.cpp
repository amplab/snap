/*++

Module Name:

    BloomPairedEndAligner.cpp

Abstract:

    A Bloom filter based paired-end aligner.

Authors:

    Bill Bolosky, February, 2013

Environment:

    User mode service.

Revision History:

--*/

#include "stdafx.h"
#include "BloomPairedEndAligner.h"
#include "SeedSequencer.h"
#include "mapq.h"
#include "BloomFilter.h"

#ifdef  COMPILE_BLOOM

volatile _int64          totalBloomFilterAdds = 0;
volatile _int64          totalBloomFilterLookups = 0;
volatile _int64          totalBloomFilterPasses = 0;


BloomPairedEndAligner::BloomPairedEndAligner(
        GenomeIndex  *index_,
        unsigned      maxReadSize_,
        unsigned      maxHits_,
        unsigned      maxK_,
        unsigned      maxSeeds_,
        unsigned      minSpacing_,                 // Minimum distance to allow between the two ends.
        unsigned      maxSpacing_,                 // Maximum distance to allow between the two ends.
        BigAllocator  *allocator) :
    index(index_), maxReadSize(maxReadSize_), maxHits(maxHits_), maxK(maxK_), maxSeeds(maxSeeds_), minSpacing(minSpacing_), maxSpacing(maxSpacing_),
    landauVishkin(NULL), reverseLandauVishkin(NULL), maxBigHits(50000)
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

    hashTableEpochNumber = 0;
    distanceToSearchBeyondBestScore = 2;    // If probability goes down by about 1000x per mismatch, then going more than two steps away make an effect of 1 part per billion, which is too small to care about
}
    
BloomPairedEndAligner::~BloomPairedEndAligner()
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
BloomPairedEndAligner::getBigAllocatorReservation(GenomeIndex * index, unsigned maxHitsToConsider, unsigned maxReadSize, unsigned seedLen, unsigned maxSeedsToUse)
{
    CountingBigAllocator countingAllocator;
    {
        BloomPairedEndAligner aligner; // This has to be in a nested scope so it's destructor is called before that of the countingAllocator
        aligner.index = index;

        aligner.allocateDynamicMemory(&countingAllocator, maxReadSize, maxHitsToConsider, maxSeedsToUse);
        return sizeof(aligner) + countingAllocator.getMemoryUsed();
    }
}

    void
BloomPairedEndAligner::allocateDynamicMemory(BigAllocator *allocator, unsigned maxReadSize, unsigned maxHitsToConsider, unsigned maxSeedsToUse)
{
    seedUsed = (BYTE *) allocator->allocate(index->getSeedLength() + 7 / 8);

    landauVishkin = new(allocator) LandauVishkin<>();
    reverseLandauVishkin = new(allocator) LandauVishkin<-1>();

    weightLists = (CandidateGroup *)allocator->allocate(sizeof(CandidateGroup) * (maxSeedsToUse + 1));
 
    candidateGroupHashTableSize = maxHitsToConsider * maxSeedsToUse;
    superGroupHashTableSize = maxHitsToConsider * maxSeedsToUse / 2;    // Probably some of these will be close together

    candidateGroupPoolSize = maxHitsToConsider * maxSeedsToUse * NUM_READS_PER_PAIR * NUM_DIRECTIONS;
    superGroupPoolSize = maxHitsToConsider * maxSeedsToUse * NUM_READS_PER_PAIR * NUM_DIRECTIONS;

    for (unsigned whichRead = 0; whichRead < NUM_READS_PER_PAIR; whichRead++) {
        rcReadData[whichRead] = (char *)allocator->allocate(maxReadSize);
        rcReadQuality[whichRead] = (char *)allocator->allocate(maxReadSize);
        bloomFilter[whichRead] = new BloomFilter(maxBigHits * maxSeeds, maxSpacing - minSpacing);

        for (Direction dir = 0; dir < NUM_DIRECTIONS; dir++) {
            reversedRead[whichRead][dir] = (char *)allocator->allocate(maxReadSize);
            superGroupHashTable[whichRead][dir] = (SuperGroupHashTableAnchor *)allocator->allocate(sizeof(SuperGroupHashTableAnchor) * superGroupHashTableSize);
            candidateGroupHashTable[whichRead][dir] = (CandidateGroupHashTableAnchor *)allocator->allocate(sizeof(CandidateGroupHashTableAnchor) * candidateGroupHashTableSize);
 
        }
    }

    superGroupPool = (SuperGroup *) allocator->allocate(sizeof(SuperGroup) * superGroupPoolSize);
    candidateGroupPool = (CandidateGroup *)allocator->allocate(sizeof(CandidateGroup) * candidateGroupPoolSize);

    baseAligner = new(allocator) BaseAligner(index, 1, maxHitsToConsider, maxK/2, maxReadSize, maxSeedsToUse, 4, landauVishkin, NULL, NULL, allocator);
}

    void 
BloomPairedEndAligner::align(
        Read                  *read0,
        Read                  *read1,
        PairedAlignmentResult *result)
{
    Read rcReads[NUM_READS_PER_PAIR];

    unsigned nSeedsApplied[NUM_READS_PER_PAIR][NUM_DIRECTIONS];

    _int64 bloomFilterAdds = 0;
    _int64 bloomFilterLookups = 0;
    _int64 bloomFilterPasses = 0;    

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

        for (Direction dir = FORWARD; dir < NUM_DIRECTIONS; dir++) {
            totalHashTableHits[whichRead][dir] = 0;
            largestHashTableHit[whichRead][dir] = 0;
        }

        if (readLen[whichRead] > maxReadSize) {
            fprintf(stderr,"BloomPairedEndAligner:: got too big read (%d > %d)", readLen[whichRead], maxReadSize);
            exit(1);
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
            nSeedsApplied[whichRead][dir] = 0;
            bestPassSeedsNotSkipped[whichRead][dir] = 0;
            Read *read = reads[whichRead][dir];

            for (unsigned i = 0; i < read->getDataLength(); i++) {
                reversedRead[whichRead][dir][i] = read->getData()[read->getDataLength() - i - 1];
            }
        }
    }

    clearCandidates();



    unsigned nPossibleSeeds = __max(readLen[0], readLen[1]) - seedLen + 1;
    unsigned thisPassSeedsNotSkipped[NUM_READS_PER_PAIR][NUM_DIRECTIONS] = {{0,0}, {0,0}}; 

    //
    // Initialize the member variables that are effectively stack locals, but are in the object
    // to aviod having to pass them to score.
    //
    probabilityOfBestPair = 0;
    localBestPairProbability[0] = 0;
    localBestPairProbability[1] = 0;
    probabilityOfAllPairs = 0;
    bestPairScore = -1;
    scoreLimit = maxK + distanceToSearchBeyondBestScore;



    //
    // Phase 1: do the hash table lookups for each of the seeds for each of the reads, and count how many there are.  This will allow us to determine which
    // read in each direction is larger and so needs to go into the Bloom filter.
    //
    for (unsigned whichRead = 0; whichRead < NUM_READS_PER_PAIR; whichRead++) {
        unsigned nextSeedToTest = 0;
        unsigned wrapCount = 0;
        memset(seedUsed, 0, (__max(readLen[0], readLen[1]) + 7) / 8);
    
        while (countOfHashTableLookups[whichRead] < nPossibleSeeds && countOfHashTableLookups[whichRead] < maxSeeds) {
            if (nextSeedToTest > nPossibleSeeds) {
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
                continue;
            }

            Seed seed(reads[whichRead][FORWARD]->getData() + nextSeedToTest, seedLen);
            //
            // Find all instances of this seed in the genome.
            //
            HashTableHit *hitRecord = &hashTableHits[whichRead][countOfHashTableLookups[whichRead]];
            hitRecord->seedOffset = nextSeedToTest;
            index->lookupSeed(seed, 0, genomeSize, &hitRecord->nHits[FORWARD], &hitRecord->hits[FORWARD], &hitRecord->nHits[RC], &hitRecord->hits[RC]);
            countOfHashTableLookups[whichRead]++;
            for (Direction dir = FORWARD; dir < NUM_DIRECTIONS; dir++) {
                totalHashTableHits[whichRead][dir] += hitRecord->nHits[dir];
                largestHashTableHit[whichRead][dir] = __max(largestHashTableHit[whichRead][dir], hitRecord->nHits[dir]);
            }

            nextSeedToTest++;
        } // while we need to lookup seeds for this read
    } // for each read

    //
    // Now select which read has the most hits.  This is the read that we'll put into the Bloom filter.
    //
    readWithMostHits = totalHashTableHits[0][FORWARD] + totalHashTableHits[0][RC] > totalHashTableHits[1][FORWARD] + totalHashTableHits[1][RC] ? 0 : 1;
    readWithFewestHits = 1 - readWithMostHits;

    //
    // Phase 2: build the Bloom fiters.
    //
    for (Direction dir = FORWARD; dir < NUM_DIRECTIONS; dir++) {
        //
        // Set up the Bloom filter with 8 bits per hash table hit.  If all of the hits are unique,
        // this will result in about a 1/4 of the bits being set and so a false positive rate of
        // 1/16.  However, in practice many of the hits will be for the same location and so the false
        // positive rate will be much lower.
        //
        bloomFilter[dir]->init(2, __max(256, 8 * totalHashTableHits[readWithMostHits][dir]));

        for (unsigned whichHit = 0; whichHit < countOfHashTableLookups[readWithMostHits]; whichHit++) {
            unsigned offset;
            HashTableHit *hitRecord = &hashTableHits[readWithMostHits][whichHit];
            if (dir == FORWARD) {
                offset = hitRecord->seedOffset;
            } else {
                offset = readLen[readWithMostHits] - seedLen - hitRecord->seedOffset;
            }
            for (unsigned whichSeedHit = 0; whichSeedHit < hitRecord->nHits[dir]; whichSeedHit++) {
                bloomFilter[dir]->addToSet(hitRecord->hits[dir][readWithMostHits] - offset);
            }
            bloomFilterAdds += hitRecord->nHits[dir];
        }
    }

    //
    // Phase 3: build up candidate lists for reads, excluding ones that certainly do not have mates.
    //
    for (Direction dir = FORWARD; dir < NUM_DIRECTIONS; dir++) {
        BloomFilter *ourBloomFilter = bloomFilter[OppositeDirection(dir)];
        for (unsigned whichHit = 0; whichHit < countOfHashTableLookups[readWithFewestHits]; whichHit++) {
            unsigned offset;
            HashTableHit *hitRecord = &hashTableHits[readWithFewestHits][whichHit];
            if (dir == FORWARD) {
                offset = hitRecord->seedOffset;
            } else {
                offset = readLen[readWithFewestHits] - seedLen - hitRecord->seedOffset;
            }

            bloomFilterLookups += hitRecord->nHits[dir];
            for (unsigned whichSeedHit = 0; whichSeedHit < hitRecord->nHits[dir]; whichSeedHit++) {
                unsigned hitLocation = hitRecord->hits[dir][whichSeedHit] - offset;
                unsigned minRange[2];
                unsigned maxRange[2];
                if (hitLocation < minSpacing) {
                    minRange[0] = maxRange[0] = 0;
                } else if (hitLocation < maxSpacing) {
                    minRange[0] = 0;
                    maxRange[0] = hitLocation - minSpacing;
                }

                //
                // Don't worry about the "above" range wrapping, because there must be some room past the size of the genome for
                //
                minRange[1] = hitLocation + minSpacing;
                maxRange[1] = hitLocation + maxSpacing;

                //
                // TODO: prevent cross chromosome alignments.
                //
                if (ourBloomFilter->mightRangeBeInSet(minRange[0], maxRange[0]) || ourBloomFilter->mightRangeBeInSet(minRange[1], maxRange[1])) {
                    //
                    // This might be good.  Enter it in the candidate set.
                    //
                    bloomFilterPasses++;
                }
            }
        }
    }

    InterlockedAdd64AndReturnNewValue(&totalBloomFilterAdds, bloomFilterAdds);
    InterlockedAdd64AndReturnNewValue(&totalBloomFilterLookups, bloomFilterLookups);
    InterlockedAdd64AndReturnNewValue(&totalBloomFilterPasses, bloomFilterPasses);

}

    void
BloomPairedEndAligner::alignWithBaseAligner(Read *read0, Read *read1, PairedAlignmentResult *result, int maxMapq)
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
BloomPairedEndAligner::clearCandidates()
{
    if (0 == hashTableEpochNumber) {
        //
        // Clear out the hash tables the slow way, because we're new or we've wrapped epochs.  Wraps only happen every 4B reads, so it's not
        // much overhead.
        //
        for (unsigned whichRead = 0; whichRead < NUM_READS_PER_PAIR; whichRead++) {
            for (Direction dir = 0; dir < NUM_DIRECTIONS; dir++) {
                for (unsigned i = 0; i < candidateGroupHashTableSize; i++) {
                    candidateGroupHashTable[whichRead][dir][i].candidateGroups = NULL;
                }
                for (unsigned i = 0; i < superGroupHashTableSize; i++) {
                    superGroupHashTable[whichRead][dir][i].supergroups = NULL;
                }
            }
        }
    } // if we've clearing the tables
    hashTableEpochNumber++;
    maxWeightListUsed = 0;
}

    bool 
BloomPairedEndAligner::isThereAMateCandidate(unsigned genomeLocation, unsigned whichRead, Direction dir)
{
    //
    // We need to look in two ranges: below and above the genome location.
    //
    unsigned minLoc0, maxLoc0, minLoc1, maxLoc1;
    if (genomeLocation < maxSpacing) {
        minLoc0 = 0;
    } else {
        minLoc0 = genomeLocation - maxSpacing;
    }
   
    if (genomeLocation < minSpacing) {
        maxLoc0 = 0;
    } else {
        maxLoc0 = genomeLocation - minSpacing;
    }

   if (maxSpacing > genomeSize - genomeLocation) {
       maxLoc1 = genomeSize;
   } else {
       maxLoc1 = genomeLocation + maxSpacing;
   }

   if (minSpacing > genomeSize - genomeLocation) {
       minLoc1 = genomeSize;
   } else {
       minLoc1 = genomeLocation + minSpacing;
   }

   return isThereACandidateInRange(minLoc0, maxLoc0, whichRead, dir) || isThereACandidateInRange(minLoc1, maxLoc1, whichRead, dir);
}

    bool
BloomPairedEndAligner::isThereACandidateInRange(unsigned minLoc, unsigned maxLoc, unsigned whichRead, Direction dir)
{
    for (unsigned baseOffset = minLoc - minLoc % SuperGroupSpan; baseOffset < ((maxLoc + SuperGroupSpan - 1) / SuperGroupSpan * SuperGroupSpan); baseOffset += SuperGroupSpan) {
        SuperGroupHashTableAnchor *anchor = &superGroupHashTable[whichRead][dir][hash(baseOffset / SuperGroupSpan, superGroupHashTableSize)];
        if (anchor->hashTableEpoch != hashTableEpochNumber) {
            anchor->supergroups = NULL;
        }

        SuperGroup *supergroup = anchor->supergroups;

        while (NULL != supergroup && supergroup->baseOffset != baseOffset) {
            // Walk down the hash chain looking for the right place.
            supergroup = supergroup->next;
        }

        if (NULL != supergroup) {
            //
            // Figure out the sub-part of this entry we're interested in.  First get rid of any high order bits
            // that are unneeded because we're not in that range.
            _uint64 mask = getSupergroupBitMaskForRange(baseOffset, minLoc, maxLoc);

            if (mask & supergroup->usedGroups) {
                return true;
            }
        }
    } // for each super group

    return false;
}
    void 
BloomPairedEndAligner::findCandidateAndCreateIfNotExtant(unsigned whichRead, Direction dir, unsigned genomeLocation, CandidateGroup **group, Candidate **candidate, unsigned bestPossibleScoreIfNew)
{
    unsigned baseGenomeLocation = genomeLocation - genomeLocation % GroupSpan;
    unsigned candidateIndex = genomeLocation - baseGenomeLocation;
    CandidateGroupHashTableAnchor *anchor = &candidateGroupHashTable[whichRead][dir][hash(genomeLocation / GroupSpan, candidateGroupHashTableSize)];

    if (anchor->hashTableEpoch != hashTableEpochNumber || anchor->candidateGroups == NULL) {
        //
        // This hash table bucket is stale.  Clear it out.
        //
        anchor->candidateGroups = NULL;
        anchor->hashTableEpoch = hashTableEpochNumber;
    }

    for (*group = anchor->candidateGroups; *group != NULL && (*group)->baseGenomeLocation != baseGenomeLocation; *group = (*group)->next) {
        // This loop body intentionally left blank
    }

    if (NULL == *group) {
        //
        // Allocate a new group.
        //
        if (nUsedCandidateGroups >= candidateGroupPoolSize) {
            fprintf(stderr, "BloomPairedEndAligner: out of candidate groups (%d).  This shouldn't happen.\n", candidateGroupPoolSize);
            exit(1);
        }

        *group = &candidateGroupPool[nUsedCandidateGroups];
        nUsedCandidateGroups++;

        (*group)->init(baseGenomeLocation, whichRead, dir);
        (*group)->next = anchor->candidateGroups;
        (*group)->bestPossibleScore = bestPossibleScoreIfNew;
        anchor->candidateGroups = *group;

        //
        // Mark it allocated in its supergroup.
        //
        SuperGroup *supergroup = findSuperGroupAndCreateIfNotExtant(whichRead, dir, genomeLocation);
        supergroup->usedGroups |= ((_uint64)1) << ((genomeLocation % SuperGroupSpan) / (SuperGroupSpan/GroupSpan));
    }

    *candidate = &(*group)->candidates[genomeLocation - baseGenomeLocation];
    (*group)->usedCandidates |= (((_uint64)1) << candidateIndex);
}

    BloomPairedEndAligner::SuperGroup *
BloomPairedEndAligner::findSuperGroupAndCreateIfNotExtant(unsigned whichRead, Direction dir, unsigned genomeLocation)
{
    unsigned baseGenomeLocation = genomeLocation - genomeLocation % SuperGroupSpan;
    SuperGroupHashTableAnchor *anchor = &superGroupHashTable[whichRead][dir][hash(genomeLocation/SuperGroupSpan, superGroupHashTableSize)];

    if (anchor->hashTableEpoch != hashTableEpochNumber || anchor->supergroups == NULL) {
        anchor->supergroups = NULL;
        anchor->hashTableEpoch = hashTableEpochNumber;
    }

    SuperGroup *supergroup;
    for (supergroup = anchor->supergroups; NULL != supergroup && supergroup->baseOffset != baseGenomeLocation; supergroup = supergroup->next) {
        // This loop body intentionally left blank.
    }

    if (NULL == supergroup) {
        if (nUsedSuperGroups >= superGroupPoolSize) {
            fprintf(stderr, "BloomPairedEndAligner: out of super groups (%d).  This shouldn't happen.\n", superGroupPoolSize);
            exit(1);
        }
        supergroup = &superGroupPool[nUsedSuperGroups];
        nUsedSuperGroups++;

        supergroup->init();
        supergroup->next = anchor->supergroups;
        supergroup->baseOffset = baseGenomeLocation;
        anchor->supergroups = supergroup;
    }

    return supergroup;
}

    bool 
BloomPairedEndAligner::score(
    bool                     forceResult,
    Read                    *read[NUM_READS_PER_PAIR][NUM_DIRECTIONS],
    PairedAlignmentResult   *result)
{
    do {
        CandidateGroup *groupToScore;

        do { // Loop around looking for a candidate group to score that's not already scored

            while (maxWeightListUsed > 0 && weightLists[maxWeightListUsed].weightListNext == &weightLists[maxWeightListUsed]) {
                maxWeightListUsed--;
            }

            if (0 == maxWeightListUsed) {
                if (forceResult || __min(bestPassSeedsNotSkipped[0][FORWARD] + bestPassSeedsNotSkipped[1][RC], bestPassSeedsNotSkipped[0][RC] + bestPassSeedsNotSkipped[1][FORWARD]) > scoreLimit) {
                    //
                    // It won't help us to find any more candidates from the hash table.
                    //
                    if (0 == maxWeightListUsed) {
                        //
                        // And we've scored everything that we've got.
                        //
                        if (bestPairScore != -1) {

                            for (unsigned whichRead = 0; whichRead < NUM_READS_PER_PAIR; whichRead++) {
                                unsigned mapq = 
                                result->direction[whichRead] = bestResultDirection[whichRead];
                                result->location[whichRead] = bestResultGenomeLocation[whichRead];
                                result->score[whichRead] = bestResultScore[whichRead];
                                result->mapq[whichRead] = computeMAPQ(probabilityOfAllPairs, probabilityOfBestPair, bestPairScore, popularSeedsSkipped[whichRead]);
                                result->status[whichRead] = result->mapq[whichRead] >= 10 ? SingleHit : MultipleHits;

                            }
                        } else {
                        }
                        return true;
                    }
                    forceResult = true;
                } else {
                    //
                    // Nothing to score, but not done generating candidates.
                    //
                    return false;
                }
            }

            groupToScore = weightLists[maxWeightListUsed].weightListNext;
            groupToScore->removeFromWeightList();
            if (groupToScore->bestPossibleScore > scoreLimit || groupToScore->allExtantCandidatesScored) {
                //
                // This group doens't work for us, because it's guaranteed to have too high of a score to be
                // interesting.  Don't bother scoring it, just loop around looking for another.
                //
                groupToScore = NULL;
            }
        } while (NULL == groupToScore);


        //
        // Build the list of possible mates for this candidate.
        //
        CandidateGroup *mates = NULL;
        unsigned bestPossibleScoreOfAnyMate = -1;
        addCandidateGroupsInRangeToScoringList(&mates, groupToScore->baseGenomeLocation, false, 1-groupToScore->whichRead, OppositeDirection(groupToScore->direction), &bestPossibleScoreOfAnyMate);
        addCandidateGroupsInRangeToScoringList(&mates, groupToScore->baseGenomeLocation + readLen[groupToScore->whichRead], true, 1-groupToScore->whichRead, OppositeDirection(groupToScore->direction), &bestPossibleScoreOfAnyMate);

        if (bestPossibleScoreOfAnyMate == -1 || bestPossibleScoreOfAnyMate + groupToScore->bestPossibleScore < scoreLimit) {

            scoreGroup(groupToScore, groupToScore->whichRead, groupToScore->direction);

            if (groupToScore->bestScore <= scoreLimit) {
                //
                // Run through all of the mates and score them, looking for good matches.
                //
                CandidateGroup *bestMate = NULL;
                unsigned bestMateScore = -1;
                double bestMateMatchProbability = 0;
                while (mates != NULL) {
                    if (mates->bestPossibleScore <= scoreLimit - groupToScore->bestScore) {

                        scoreGroup(mates, mates->whichRead, mates->direction);

                        if (mates->bestScore < bestMateScore || mates->bestScore == bestMateScore && mates->matchProbabilityForBestScore >= bestMateMatchProbability) {
                            bestMate = mates;
                            bestMateScore = mates->bestScore;
                            bestMateMatchProbability = mates->matchProbabilityForBestScore;
                        }
                        mates = mates->scoringNext;
                    }
                } // for each mate candidate

                if (bestMateScore != -1) {
                    double matchProbability = groupToScore->matchProbabilityForBestScore * bestMate->matchProbabilityForBestScore;
                    unsigned overallScore = bestMateScore + groupToScore->bestScore;
                    if (overallScore < bestPairScore || overallScore == bestPairScore && matchProbability > probabilityOfBestPair) {
                        //
                        // A new best score.  Remember it.
                        //
                        probabilityOfBestPair = matchProbability;
                        bestPairScore = overallScore;
                        _ASSERT(groupToScore->whichRead != bestMate->whichRead);
                        _ASSERT(groupToScore->direction != bestMate->direction);
                        bestResultGenomeLocation[groupToScore->whichRead] = groupToScore->bestScoreGenomeLocation;
                        bestResultGenomeLocation[bestMate->whichRead] = bestMate->bestScoreGenomeLocation;
                        bestResultDirection[groupToScore->whichRead] = groupToScore->direction;
                        bestResultDirection[bestMate->whichRead] = bestMate->direction;
                        bestResultScore[groupToScore->whichRead] = groupToScore->bestScore;
                        bestResultScore[bestMate->whichRead] = bestMate->bestScore;

                        scoreLimit = bestPairScore + 2;
                    }
                    probabilityOfAllPairs += matchProbability;
                } // If the mate had a low enough score
            } // If the group was scored low enough to matter
        } // If there was any hope for this group
    } while (forceResult);

    return false;
}

    void 
BloomPairedEndAligner::addCandidateGroupsInRangeToScoringList(
    CandidateGroup **       listHead, 
    unsigned                baseGenomeLocation,
    bool                    lookBelow, 
    unsigned                whichRead, 
    Direction               dir, 
    unsigned *              bestPossibleScoreOfAnyMate)
{
    unsigned rangeMin;
    unsigned rangeMax;

    *bestPossibleScoreOfAnyMate = -1;

    if (lookBelow) {
        const Genome::Piece *piece = genome->getPieceAtLocation(baseGenomeLocation);
        if (maxSpacing > baseGenomeLocation) {
            rangeMin = piece->beginningOffset;
        } else {
            rangeMin = __max(piece->beginningOffset, baseGenomeLocation - maxSpacing);
        }
        if (minSpacing > baseGenomeLocation) {
            //
            // The entire range falls off the beginning of the genome.  Just quit.
            //
            return;
        }
        rangeMax = __max(piece->beginningOffset, baseGenomeLocation - minSpacing);
    } else {
        const Genome::Piece *piece = genome->getNextPieceAfterLocation(baseGenomeLocation);
        rangeMin = __min(piece->beginningOffset - 1, baseGenomeLocation + minSpacing);
        rangeMax = __min(piece->beginningOffset - 1, baseGenomeLocation + maxSpacing);
    }

    //
    // Now run throught the supergroups in the range
    //
    for (unsigned offsetToConsider = rangeMin; offsetToConsider < rangeMax; offsetToConsider += SuperGroupSpan) {
        unsigned supergroupBase = offsetToConsider - offsetToConsider % SuperGroupSpan;
        SuperGroupHashTableAnchor *anchor = &superGroupHashTable[whichRead][dir][hash(supergroupBase, superGroupHashTableSize)];
        SuperGroup *supergroup;
        if (anchor->hashTableEpoch == hashTableEpochNumber) {
            for (supergroup = anchor->supergroups; NULL != supergroup && supergroup->baseOffset != supergroupBase; supergroup = supergroup->next) {
                // This loop body intentionally left blank.
            }

            if (NULL != supergroup) {
                //
                // See if there are any groups in the supergroup that meet the range criteria.
                //
                _uint64 candidateGroupBits = supergroup->usedGroups & getSupergroupBitMaskForRange(supergroupBase, rangeMin, rangeMax);
                unsigned long candidateToScore;
                while (_BitScanForward64(&candidateToScore, candidateGroupBits)) {
                    CandidateGroup *group = lookupCandidateGroup(supergroupBase + GroupSpan * candidateToScore, whichRead, dir);
                    _ASSERT(NULL != group); // Else its bit shouldn't have been set
                    if (!group->allExtantCandidatesScored) {
                        group->scoringNext = *listHead;
                        *listHead = group->scoringNext;
                    }

                    *bestPossibleScoreOfAnyMate = __min(*bestPossibleScoreOfAnyMate, group->bestPossibleScore);
                    candidateGroupBits ^= (((_uint64)1) << candidateToScore);
                } 
            } // if we have a supergroup
        } // if we have an anchor
    } // for each supergroup offset
}

    _uint64 
BloomPairedEndAligner::getSupergroupBitMaskForRange(unsigned supergroupBaseOffset, unsigned minLoc, unsigned maxLoc)
{
    _ASSERT(minLoc <= maxLoc);
    _ASSERT(supergroupBaseOffset % SuperGroupSpan == 0);

    //
    // Truncate minLoc and maxLoc to fit within this supergroup.
    //
    minLoc = (minLoc < supergroupBaseOffset) ? supergroupBaseOffset :  minLoc;
    maxLoc = (maxLoc >= supergroupBaseOffset + SuperGroupSpan) ? supergroupBaseOffset + SuperGroupSpan - 1 : maxLoc;

    _uint64 mask;
    //
    // First, get the right number of set bits in mask.  We special case where it's the entire range, because
    // we want to avoid doing a right shift of a value with the high bit set.
    //
    _ASSERT(maxLoc >= minLoc);
    _ASSERT(maxLoc - minLoc <= SuperGroupSpan - 1);
    if (maxLoc - minLoc >= SuperGroupSpan - GroupSpan) {
        return -1;
    }

    mask = 0x7fffffffffffffff >> (63 - (maxLoc - minLoc + GroupSpan - 1) / GroupSpan);

    //
    // Now shift it to the right place.
    //
    mask <<= (minLoc - supergroupBaseOffset) / GroupSpan;

    return mask;
}

    BloomPairedEndAligner::CandidateGroup *
BloomPairedEndAligner::lookupCandidateGroup(unsigned genomeOffset, unsigned whichRead, Direction dir)
{
    unsigned groupOffset = genomeOffset - genomeOffset % GroupSpan;

    CandidateGroupHashTableAnchor *anchor = &candidateGroupHashTable[whichRead][dir][hash(groupOffset, candidateGroupHashTableSize)];
    if (anchor->hashTableEpoch != hashTableEpochNumber) {
        return NULL;
    }

    for (CandidateGroup *group = anchor->candidateGroups; NULL != group; group = group->next) {
        if (group->baseGenomeLocation == groupOffset) {
            return group;
        }
    }
    return NULL;
}

    void 
BloomPairedEndAligner::scoreGroup(CandidateGroup *group, unsigned whichRead, Direction direction)
{
    Read *readToScore = reads[whichRead][direction];
    unsigned long candidateIndexToScore;
    while (_BitScanForward64(&candidateIndexToScore, group->usedCandidates & ~group->scoredCandidates)) {
        _uint64 candidateBit = ((_uint64)1) << candidateIndexToScore;
        group->scoredCandidates |= candidateBit;
        Candidate *candidateToScore = &group->candidates[candidateIndexToScore];

        unsigned genomeLocation = group->baseGenomeLocation + candidateIndexToScore;
        unsigned readDataLength = readToScore->getDataLength();
        unsigned genomeDataLength = readDataLength + MAX_K; // Leave extra space in case the read has deletions
        const char *data = genome->getSubstring(genomeLocation, genomeDataLength);
        
        _ASSERT(candidateToScore->seedOffset + seedLen <= readToScore->getDataLength());

        int maxStartShift = __min(scoreLimit, 7);
        for (int shift = 1; shift <= maxStartShift; shift++) {
            group->scoredCandidates |= (candidateBit << shift);
            group->scoredCandidates |= (candidateBit >> shift);
        }
        
        // Compute the distance separately in the forward and backward directions from the seed, to allow
        // arbitrary offsets at both the start and end but not have to pay the cost of exploring all start
        // shifts in BoundedStringDistance
        double matchProb1, matchProb2;
        int score1, score2;
        // First, do the forward direction from where the seed aligns to past of it
        int readLen = readToScore->getDataLength();
        int seedLen = index->getSeedLength();
        int seedOffset = candidateToScore->seedOffset; // Since the data is reversed
        int tailStart = seedOffset + seedLen;
        double matchProbability;

        score1 = landauVishkin->computeEditDistance(data + tailStart, genomeDataLength - tailStart, readToScore->getData() + tailStart, readToScore->getQuality() + tailStart, readLen - tailStart,
            scoreLimit, &matchProb1);
        if (score1 == -1) {
            candidateToScore->score = -1;
            matchProbability = 0;
        } else {
            // The tail of the read matched; now let's reverse the reference genome data and match the head
            int limitLeft = scoreLimit - score1;
            score2 = reverseLandauVishkin->computeEditDistance(data + seedOffset, seedOffset + MAX_K, reversedRead[whichRead][direction] + readLen - seedOffset,
                                                                        reads[whichRead][OppositeDirection(direction)]->getQuality() + readLen - seedOffset, seedOffset, limitLeft, &matchProb2);

            if (score2 == -1) {
                candidateToScore->score = -1;
            } else {
                candidateToScore->score = score1 + score2;
                // Map probabilities for substrings can be multiplied, but make sure to count seed too
                matchProbability = matchProb1 * matchProb2 * pow(1 - SNP_PROB, seedLen);
            }
        }

        if (candidateToScore->score < group->bestScore || candidateToScore->score == group->bestScore && matchProbability > group->matchProbabilityForBestScore) {
            group->bestScore = candidateToScore->score;
            group->matchProbabilityForBestScore = matchProbability;
            group->bestScoreGenomeLocation = genomeLocation;
        }
    }// for each unscored candidate bit
    group->allExtantCandidatesScored = true;
}

#endif  // COMPILE_BLOOM
