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
#include <xmmintrin.h>

#ifdef  _DEBUG
extern bool _DumpAlignments;    // From BaseAligner.cpp
#endif  // _DEBUG

IntersectingPairedEndAligner::IntersectingPairedEndAligner(
        GenomeIndex  *index_,
        unsigned      maxReadSize_,
        unsigned      maxHits_,
        unsigned      maxK_,
        unsigned      numSeedsFromCommandLine_,
        double        seedCoverage_,
        unsigned      minSpacing_,                 // Minimum distance to allow between the two ends.
        unsigned      maxSpacing_,                 // Maximum distance to allow between the two ends.
        unsigned      maxBigHits_,
        unsigned      extraSearchDepth_,
        BigAllocator  *allocator) :
    index(index_), maxReadSize(maxReadSize_), maxHits(maxHits_), maxK(maxK_), numSeedsFromCommandLine(__min(MAX_MAX_SEEDS,numSeedsFromCommandLine_)), minSpacing(minSpacing_), maxSpacing(maxSpacing_),
    landauVishkin(NULL), reverseLandauVishkin(NULL), maxBigHits(maxBigHits_), maxMergeDistance(31), seedCoverage(seedCoverage_) /*also should be a parameter*/,
    extraSearchDepth(extraSearchDepth_), nLocationsScored(0)
{
    unsigned maxSeedsToUse;
    if (0 != numSeedsFromCommandLine) {
        maxSeedsToUse = numSeedsFromCommandLine;
    } else {
        maxSeedsToUse = (unsigned)(maxReadSize * seedCoverage / index->getSeedLength());
    }
    allocateDynamicMemory(allocator, maxReadSize, maxBigHits, maxSeedsToUse, maxK, extraSearchDepth);

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
}
    
IntersectingPairedEndAligner::~IntersectingPairedEndAligner()
{
}
    
    size_t 
IntersectingPairedEndAligner::getBigAllocatorReservation(GenomeIndex * index, unsigned maxBigHitsToConsider, unsigned maxReadSize, unsigned seedLen, unsigned numSeedsFromCommandLine, 
                                                         double seedCoverage, unsigned maxEditDistanceToConsider, unsigned maxExtraSearchDepth)
{
    unsigned maxSeedsToUse;
    if (0 != numSeedsFromCommandLine) {
        maxSeedsToUse = numSeedsFromCommandLine;
    } else {
        maxSeedsToUse = (unsigned)(maxReadSize * seedCoverage / index->getSeedLength());
    }
    CountingBigAllocator countingAllocator;
    {
        IntersectingPairedEndAligner aligner; // This has to be in a nested scope so it's destructor is called before that of the countingAllocator
        aligner.index = index;

        aligner.allocateDynamicMemory(&countingAllocator, maxReadSize, maxBigHitsToConsider, maxSeedsToUse, maxEditDistanceToConsider, maxExtraSearchDepth);
        return sizeof(aligner) + countingAllocator.getMemoryUsed();
    }
}

    void
IntersectingPairedEndAligner::allocateDynamicMemory(BigAllocator *allocator, unsigned maxReadSize, unsigned maxBigHitsToConsider, unsigned maxSeedsToUse, 
                                                    unsigned maxEditDistanceToConsider, unsigned maxExtraSearchDepth)
{
    seedUsed = (BYTE *) allocator->allocate(100 + (maxReadSize + 7) / 8);

    for (unsigned whichRead = 0; whichRead < NUM_READS_PER_PAIR; whichRead++) {
        rcReadData[whichRead] = (char *)allocator->allocate(maxReadSize);
        rcReadQuality[whichRead] = (char *)allocator->allocate(maxReadSize);

        for (Direction dir = 0; dir < NUM_DIRECTIONS; dir++) {
            reversedRead[whichRead][dir] = (char *)allocator->allocate(maxReadSize);
            hashTableHitSets[whichRead][dir] =(HashTableHitSet *)allocator->allocate(sizeof(HashTableHitSet));/* new HashTableHitSet()*/;
            hashTableHitSets[whichRead][dir]->firstInit(maxSeedsToUse, maxEditDistanceToConsider);
        }
    }

#if 0


    for (int i = 0; i < NUM_SET_PAIRS; i++) {
        hitLocations[i] = new(allocator) HitLocationRingBuffer(maxEditDistanceToConsider * 2 + 2);      // *2 is for before & after, +2 is just to make sure the ring buffer tail doesn't overrun the head
        mateHitLocations[i] = new(allocator) HitLocationRingBuffer(2 * (maxSpacing + 1) + 2);  // Likewise.
    }

#endif // 0

    scoringCandidatePoolSize = maxBigHitsToConsider * maxSeedsToUse * NUM_READS_PER_PAIR;

    scoringCandidates = (ScoringCandidate **) allocator->allocate(sizeof(ScoringCandidate *) * (maxEditDistanceToConsider + maxExtraSearchDepth + 1));  //+1 is for 0.
    scoringCandidatePool = (ScoringCandidate *)allocator->allocate(sizeof(ScoringCandidate) * scoringCandidatePoolSize);

    for (unsigned i = 0; i < NUM_READS_PER_PAIR; i++) {
        scoringMateCandidates[i] = (ScoringMateCandidate *) allocator->allocate(sizeof(ScoringMateCandidate) * scoringCandidatePoolSize / NUM_READS_PER_PAIR);
    }

    mergeAnchorPoolSize = scoringCandidatePoolSize;
    mergeAnchorPool = (MergeAnchor *)allocator->allocate(sizeof(MergeAnchor) * mergeAnchorPoolSize);
}

    void 
IntersectingPairedEndAligner::align(
        Read                  *read0,
        Read                  *read1,
        PairedAlignmentResult *result)
{
    result->nLVCalls = 0;
    result->nSmallHits = 0;

    unsigned maxSeeds;
    if (numSeedsFromCommandLine != 0) {
        maxSeeds = numSeedsFromCommandLine;
    } else {
        maxSeeds = (unsigned)(max(read0->getDataLength(), read1->getDataLength()) * seedCoverage / index->getSeedLength());
    }

#ifdef  _DEBUG
    if (_DumpAlignments) {
        printf("\nIntersectingAligner aligning reads '%*.s' and '%.*s' with data '%.*s' and '%.*s'\n", read0->getIdLength(), read0->getId(), read1->getIdLength(), read1->getId(), read0->getDataLength(), read0->getData(), read1->getDataLength(), read1->getData());
    }
#endif  // _DEBUG

    lowestFreeScoringCandidatePoolEntry = 0;
    for (unsigned k = 0; k <= maxK + extraSearchDepth; k++) {
        scoringCandidates[k] = NULL;
    }

    for (unsigned i = 0; i < NUM_SET_PAIRS; i++) {
        lowestFreeScoringMateCandidate[i] = 0;
    }
    firstFreeMergeAnchor = 0;

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
#if 0
        hitLocations[whichRead]->clear();
        mateHitLocations[whichRead]->clear();
#endif // 0

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
    // to avoid having to pass them to score.
    //
    double probabilityOfBestPair = 0;
    localBestPairProbability[0] = 0;
    localBestPairProbability[1] = 0;
    double probabilityOfAllPairs = 0;
    unsigned bestPairScore = 65536;
    unsigned scoreLimit = maxK + extraSearchDepth;

    //
    // Phase 1: do the hash table lookups for each of the seeds for each of the reads and add them to the hit sets.
    //
    for (unsigned whichRead = 0; whichRead < NUM_READS_PER_PAIR; whichRead++) {
        unsigned nextSeedToTest = 0;
        unsigned wrapCount = 0;
        unsigned nPossibleSeeds = readLen[whichRead] - seedLen + 1;
        memset(seedUsed, 0, (__max(readLen[0], readLen[1]) + 7) / 8);
		bool beginsDisjointHitSet = whichRead == 0;
    
        while (countOfHashTableLookups[whichRead] < nPossibleSeeds && countOfHashTableLookups[whichRead] < maxSeeds) {
            if (nextSeedToTest >= nPossibleSeeds) {
                wrapCount++;
				beginsDisjointHitSet = true;
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

            index->lookupSeed(seed, 0, InvalidGenomeLocation, &nHits[FORWARD], &hits[FORWARD], &nHits[RC], &hits[RC]);

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
                    hashTableHitSets[whichRead][dir]->recordLookup(offset, nHits[dir], hits[dir], beginsDisjointHitSet);
                } else {
                    popularSeedsSkipped[whichRead]++;
                }
            }

            nextSeedToTest += seedLen;
        } // while we need to lookup seeds for this read
    } // for each read

    readWithMoreHits = totalHashTableHits[0][FORWARD] + totalHashTableHits[0][RC] > totalHashTableHits[1][FORWARD] + totalHashTableHits[1][RC] ? 0 : 1;
    readWithFewerHits = 1 - readWithMoreHits;

#ifdef  _DEBUG
    if (_DumpAlignments) {
        printf("Read 0 has %d hits, read 1 has %d hits\n", totalHashTableHits[0][FORWARD] + totalHashTableHits[0][RC], totalHashTableHits[1][FORWARD] + totalHashTableHits[1][RC]);
    }
#endif  // _DEBUG
        
    Direction setPairDirection[NUM_SET_PAIRS][NUM_READS_PER_PAIR] = {{FORWARD, RC}, {RC, FORWARD}};


    //
    // Phase 2: find all possible candidates and add them to candidate lists (for the reads with fewer and more hits).
    //
    unsigned maxUsedBestPossibleScoreList = 0;

    for (unsigned whichSetPair = 0; whichSetPair < NUM_SET_PAIRS; whichSetPair++) {
        HashTableHitSet *setPair[NUM_READS_PER_PAIR];
        if (whichSetPair == 0) {
            setPair[0] = hashTableHitSets[0][FORWARD];
            setPair[1] = hashTableHitSets[1][RC];
        } else {
            setPair[0] = hashTableHitSets[0][RC];
            setPair[1] = hashTableHitSets[1][FORWARD];
        }

        unsigned            lastSeedOffsetForReadWithFewerHits;
        unsigned            lastGenomeLocationForReadWithFewerHits;
        unsigned            lastGenomeLocationForReadWithMoreHits;
		unsigned			bestPossibleScoreForReadWithFewerHits;
        unsigned            lastSeedOffsetForReadWithMoreHits;
        unsigned            bestPossibleScoreForReadWithMoreHits;

        bool                outOfMoreHitsLocations = false;

        //
        // Seed the intersection state by doing a first lookup.
        //
        if (!setPair[readWithFewerHits]->getFirstHit(&lastGenomeLocationForReadWithFewerHits, &lastSeedOffsetForReadWithFewerHits, &bestPossibleScoreForReadWithFewerHits)) {
            //
            // No hits in this direction.
            //
            continue;   // The outer loop over set pairs.
        }

        lastGenomeLocationForReadWithMoreHits = InvalidGenomeLocation;

        //
        // Loop over the candidates in for the read with more hits.  At the top of the loop, we have a candidate but don't know if it has
        // a mate.  Each pass through the loop considers a single hit on the read with fewer hits.
        //
        for (;;) {

            //
            // Loop invariant: lastGenomeLocationForReadWithFewerHits is the highest genome offset that has not been considered.
            // lastGenomeLocationForReadWithMoreHits is also the highest genome offset on that side that has not been
            // considered (or is InvalidGenomeLocation), but higher ones within the appropriate range might already be in scoringMateCandidates.
            // We go once through this loop for each 
            //

            if (lastGenomeLocationForReadWithMoreHits > lastGenomeLocationForReadWithFewerHits + maxSpacing) {
                //
                // The more hits side is too high to be a mate candidate for the fewer hits side.  Move it down to the largest
                // location that's not too high.
                //
                if (!setPair[readWithMoreHits]->getNextHitLessThanOrEqualTo(lastGenomeLocationForReadWithFewerHits + maxSpacing, 
                                                                             &lastGenomeLocationForReadWithMoreHits, &lastSeedOffsetForReadWithMoreHits,
                                                                             &bestPossibleScoreForReadWithMoreHits)) {
                    break;  // End of all of the mates.  We're done with this set pair.
                }
            } 
            
            if (lastGenomeLocationForReadWithMoreHits + maxSpacing < lastGenomeLocationForReadWithFewerHits && 
                (0 == lowestFreeScoringMateCandidate[whichSetPair] || 
                !isWithin(scoringMateCandidates[whichSetPair][lowestFreeScoringMateCandidate[whichSetPair]-1].readWithMoreHitsGenomeLocation, lastGenomeLocationForReadWithFewerHits, maxSpacing))) {
                //
                // No mates for the hit on the read with fewer hits.  Skip to the next candidate.
                //
                if (!setPair[readWithFewerHits]->getNextHitLessThanOrEqualTo(lastGenomeLocationForReadWithMoreHits + maxSpacing, &lastGenomeLocationForReadWithFewerHits, &lastSeedOffsetForReadWithFewerHits,
                                                                             &bestPossibleScoreForReadWithFewerHits)) {
                    //
                    // No more candidates on the read with fewer hits side.  We're done with this set pair.
                    //
                    break;
                }
                continue;
            }

            //
            // Add all of the mate candidates for this fewer side hit.
            //

            unsigned previousMoreHitsLocation = lastGenomeLocationForReadWithMoreHits;
            while (lastGenomeLocationForReadWithMoreHits + maxSpacing >= lastGenomeLocationForReadWithFewerHits && !outOfMoreHitsLocations) {
                _ASSERT(lowestFreeScoringMateCandidate[whichSetPair] < scoringCandidatePoolSize / NUM_SET_PAIRS);   // Because we allocated an upper bound number of them

                scoringMateCandidates[whichSetPair][lowestFreeScoringMateCandidate[whichSetPair]].init(
                                lastGenomeLocationForReadWithMoreHits, bestPossibleScoreForReadWithMoreHits, lastSeedOffsetForReadWithMoreHits);

#ifdef _DEBUG
                if (_DumpAlignments) {
                    printf("SetPair %d, added more hits candidate %d at genome location %u, bestPossibleScore %d, seedOffset %d\n",
                            whichSetPair, lowestFreeScoringMateCandidate[whichSetPair], lastGenomeLocationForReadWithMoreHits, 
                            bestPossibleScoreForReadWithMoreHits,
                            lastSeedOffsetForReadWithMoreHits);
                }
#endif // _DEBUG
                    
                lowestFreeScoringMateCandidate[whichSetPair]++;

                previousMoreHitsLocation = lastGenomeLocationForReadWithMoreHits;

                if (!setPair[readWithMoreHits]->getNextLowerHit(&lastGenomeLocationForReadWithMoreHits, &lastSeedOffsetForReadWithMoreHits, &bestPossibleScoreForReadWithMoreHits)) {
                    lastGenomeLocationForReadWithMoreHits = 0;
                    outOfMoreHitsLocations = true;
                    break; // out of the loop looking for candidates on the more hits side.
                }
            }

            //
            // And finally add the hit from the fewer hit side.  To compute its best possible score, we need to look at all of the mates; we couldn't do it in the
            // loop immediately above because some of them might have already been in the mate list from a different, nearby fewer hit location.
            //

            unsigned lowestBestPossibleScoreOfAnyPossibleMate = maxK + extraSearchDepth;
            for (int i = lowestFreeScoringMateCandidate[whichSetPair] - 1; i >= 0; i--) {
                if (scoringMateCandidates[whichSetPair][i].readWithMoreHitsGenomeLocation > lastGenomeLocationForReadWithFewerHits + maxSpacing) {
                    break;
                }
                lowestBestPossibleScoreOfAnyPossibleMate = __min(lowestBestPossibleScoreOfAnyPossibleMate, scoringMateCandidates[whichSetPair][i].bestPossibleScore);
            }

            if (lowestBestPossibleScoreOfAnyPossibleMate + bestPossibleScoreForReadWithFewerHits <= maxK + extraSearchDepth) {
                //
                // There's a set of ends that we can't prove doesn't have too large of a score.  Allocate a fewer hit candidate and stick it in the
                // correct weight list.
                //
                _ASSERT(lowestFreeScoringCandidatePoolEntry < scoringCandidatePoolSize); // because we allocated an upper bound size of them.
                scoringCandidatePool[lowestFreeScoringCandidatePoolEntry].init(lastGenomeLocationForReadWithFewerHits, whichSetPair, lowestFreeScoringMateCandidate[whichSetPair] - 1,
                                                                                lastSeedOffsetForReadWithFewerHits, bestPossibleScoreForReadWithFewerHits,
                                                                                scoringCandidates[lowestBestPossibleScoreOfAnyPossibleMate + bestPossibleScoreForReadWithFewerHits]);

                scoringCandidates[lowestBestPossibleScoreOfAnyPossibleMate + bestPossibleScoreForReadWithFewerHits] = &scoringCandidatePool[lowestFreeScoringCandidatePoolEntry];
 
#ifdef _DEBUG
                if (_DumpAlignments) {
                    printf("SetPair %d, added fewer hits candidate %d at genome location %u, bestPossibleScore %d, seedOffset %d\n",
                            whichSetPair, lowestFreeScoringCandidatePoolEntry, lastGenomeLocationForReadWithFewerHits, 
                            lowestBestPossibleScoreOfAnyPossibleMate + bestPossibleScoreForReadWithFewerHits,
                            lastSeedOffsetForReadWithFewerHits);
                }
#endif // _DEBUG

                lowestFreeScoringCandidatePoolEntry++;
                maxUsedBestPossibleScoreList = max(maxUsedBestPossibleScoreList, lowestBestPossibleScoreOfAnyPossibleMate + bestPossibleScoreForReadWithFewerHits);

                if (!setPair[readWithFewerHits]->getNextLowerHit(&lastGenomeLocationForReadWithFewerHits, &lastSeedOffsetForReadWithFewerHits, &bestPossibleScoreForReadWithFewerHits)) {
                    break;
                }
             }
        }
    } // For each set pair

    //
    // Phase 3: score and merge the candidates we've found.
    //
    unsigned currentBestPossibleScoreList = 0;
    scoreLimit = maxK + extraSearchDepth;
    //
    // Loop until we've scored all of the candidates, or proven that what's left must have too high of a score to be interesting.
    //
    while (currentBestPossibleScoreList <= maxUsedBestPossibleScoreList && currentBestPossibleScoreList <= scoreLimit) {
        if (scoringCandidates[currentBestPossibleScoreList] == NULL) {
            //
            // No more candidates on this list.  Skip to the next one.
            //
            currentBestPossibleScoreList++;
            continue;
        }

        //
        // Grab the first candidate on the highest list and score it.
        //
        ScoringCandidate *candidate = scoringCandidates[currentBestPossibleScoreList];
 
        unsigned fewerEndScore;
        double fewerEndMatchProbability;
        int fewerEndGenomeLocationOffset;

        scoreLocation(readWithFewerHits, setPairDirection[candidate->whichSetPair][readWithFewerHits], candidate->readWithFewerHitsGenomeLocation, 
            candidate->seedOffset, scoreLimit, &fewerEndScore, &fewerEndMatchProbability, &fewerEndGenomeLocationOffset);

        _ASSERT(-1 == fewerEndScore || fewerEndScore >= candidate->bestPossibleScore);

#ifdef _DEBUG
        if (_DumpAlignments) {
            printf("Scored fewer end candidate %d, set pair %d, read %d, location %u, seed offset %d, score limit %d, score %d, offset %d\n", (int)(candidate - scoringCandidatePool),
                candidate->whichSetPair, readWithFewerHits, candidate->readWithFewerHitsGenomeLocation, candidate->seedOffset,
                scoreLimit, fewerEndScore, fewerEndGenomeLocationOffset);
        }
#endif // DEBUG

        if (fewerEndScore != -1) {
            //
            // Find and score mates.  The index in scoringMateCandidateIndex is the lowest mate (i.e., the highest index number).
            //
            unsigned mateIndex = candidate->scoringMateCandidateIndex;

            for (;;) {

                ScoringMateCandidate *mate = &scoringMateCandidates[candidate->whichSetPair][mateIndex];
                _ASSERT(isWithin(mate->readWithMoreHitsGenomeLocation, candidate->readWithFewerHitsGenomeLocation, maxSpacing));
                if (!isWithin(mate->readWithMoreHitsGenomeLocation, candidate->readWithFewerHitsGenomeLocation, minSpacing) && mate->bestPossibleScore <= scoreLimit - fewerEndScore) {
                    //
                    // It's within the range and not necessarily too poor of a match.  Consider it.
                    //

                    //
                    // If we haven't yet scored this mate, or we've scored it and not gotten an answer, but had a higher score limit than we'd
                    // use now, score it.
                    //
                    if (mate->score == -2 || mate->score == -1 && mate->scoreLimit < scoreLimit - fewerEndScore) {
                        scoreLocation(readWithMoreHits, setPairDirection[candidate->whichSetPair][readWithMoreHits], mate->readWithMoreHitsGenomeLocation,
                            mate->seedOffset, scoreLimit - fewerEndScore, &mate->score, &mate->matchProbability, 
                            &mate->genomeOffset);
#ifdef _DEBUG
                        if (_DumpAlignments) {
                            printf("Scored mate candidate %d, set pair %d, read %d, location %u, seed offset %d, score limit %d, score %d, offset %d\n",
                                (int)(mate - scoringMateCandidates[candidate->whichSetPair]), candidate->whichSetPair, readWithMoreHits, mate->readWithMoreHitsGenomeLocation,
                                mate->seedOffset, scoreLimit - fewerEndScore, mate->score, mate->genomeOffset);
                        }
#endif // _DEBUG

                        _ASSERT(-1 == mate->score || mate->score >= mate->bestPossibleScore);

                        mate->scoreLimit = scoreLimit - fewerEndScore;
                    }

                    if (mate->score != -1) {
                        double pairProbability = mate->matchProbability * fewerEndMatchProbability;
                        unsigned pairScore = mate->score + fewerEndScore;
                        //
                        // See if this should be ignored as a merge, or if we need to back out a previously scored location
                        // because it's a worse version of this location.
                        //
                        MergeAnchor *mergeAnchor = candidate->mergeAnchor;

                        if (NULL == mergeAnchor) {
                            //
                            // Look up and down the array of candidates to see if we have possible merge candidates.
                            //
                            for (ScoringCandidate *mergeCandidate = candidate - 1; 
                                        mergeCandidate >= scoringCandidatePool && 
                                        isWithin(mergeCandidate->readWithFewerHitsGenomeLocation, candidate->readWithFewerHitsGenomeLocation + fewerEndGenomeLocationOffset, 50) &&
                                        mergeCandidate->whichSetPair == candidate->whichSetPair;
                                        mergeCandidate--) {

                                if (mergeCandidate->mergeAnchor != NULL) {
                                    candidate->mergeAnchor = mergeAnchor = mergeCandidate->mergeAnchor;
                                    break;
                                }
                            }

                            if (NULL == mergeAnchor) {
                                for (ScoringCandidate *mergeCandidate = candidate + 1; 
                                            mergeCandidate < scoringCandidatePool + lowestFreeScoringCandidatePoolEntry && 
                                            isWithin(mergeCandidate->readWithFewerHitsGenomeLocation, candidate->readWithFewerHitsGenomeLocation + fewerEndGenomeLocationOffset, 50) &&
                                            mergeCandidate->whichSetPair == candidate->whichSetPair;
                                            mergeCandidate--) {

                                    if (mergeCandidate->mergeAnchor != NULL) {
                                        candidate->mergeAnchor = mergeAnchor = mergeCandidate->mergeAnchor;
                                        break;
                                    }
                                }
                            }
                        }

                        bool merged;

                        double oldPairProbability;

                        if (NULL == mergeAnchor) {
                            _ASSERT(firstFreeMergeAnchor < mergeAnchorPoolSize);
                            mergeAnchor = &mergeAnchorPool[firstFreeMergeAnchor];
                            firstFreeMergeAnchor++;

                            mergeAnchor->init(mate->readWithMoreHitsGenomeLocation + mate->genomeOffset, candidate->readWithFewerHitsGenomeLocation + fewerEndGenomeLocationOffset, 
                                pairProbability, pairScore);

                            merged = false;
                            oldPairProbability = 0;
                            candidate->mergeAnchor = mergeAnchor;
                        } else {
                            merged = mergeAnchor->checkMerge(mate->readWithMoreHitsGenomeLocation + mate->genomeOffset, candidate->readWithFewerHitsGenomeLocation + fewerEndGenomeLocationOffset, 
                                pairProbability, pairScore, &oldPairProbability);
                        }

                        if (!merged) {
                            //
                            // Back out the probability of the old match that we're merged with, if any.  The max
                            // is necessary because a + b - b is not necessarily a in floating point.  If there
                            // was no merge, the oldPairProbability is 0.
                            //
                            probabilityOfAllPairs = __max(0, probabilityOfAllPairs - oldPairProbability);

                            bool isBestHit = false;

                            if (pairScore <= maxK && (pairScore < bestPairScore || pairScore == bestPairScore && pairProbability > probabilityOfBestPair)) {
                                //
                                // A new best hit.
                                //
                                bestPairScore = pairScore;
                                probabilityOfBestPair = pairProbability;
                                bestResultGenomeLocation[readWithFewerHits] = candidate->readWithFewerHitsGenomeLocation + fewerEndGenomeLocationOffset;
                                bestResultGenomeLocation[readWithMoreHits] = mate->readWithMoreHitsGenomeLocation + mate->genomeOffset;
                                bestResultScore[readWithFewerHits] = fewerEndScore;
                                bestResultScore[readWithMoreHits] = mate->score;
                                bestResultDirection[readWithFewerHits] = setPairDirection[candidate->whichSetPair][readWithFewerHits];
                                bestResultDirection[readWithMoreHits] = setPairDirection[candidate->whichSetPair][readWithMoreHits];

                                scoreLimit = bestPairScore + extraSearchDepth;

                                isBestHit = true;
                            }

                            probabilityOfAllPairs += pairProbability;
    #ifdef  _DEBUG
                            if (_DumpAlignments) {
                                printf("Added %e (= %e * %e) @ (%u, %u), giving new probability of all pairs %e, score %d = %d + %d%s\n", 
                                    pairProbability, mate->matchProbability , fewerEndMatchProbability, 
                                    candidate->readWithFewerHitsGenomeLocation + fewerEndGenomeLocationOffset, mate->readWithMoreHitsGenomeLocation + mate->genomeOffset,
                                    probabilityOfAllPairs,
                                    pairScore, fewerEndScore, mate->score, isBestHit ? " New best hit" : "");
                            }
    #endif  // _DEBUG

                            if (probabilityOfAllPairs >= 4.9) {
                                //
                                // Nothing will rescue us from a 0 MAPQ, so just stop looking.
                                //
                                goto doneScoring;
                            }
                        }
                    }// if the mate has a non -1 score
                }

                if (mateIndex == 0 || !isWithin(scoringMateCandidates[candidate->whichSetPair][mateIndex-1].readWithMoreHitsGenomeLocation, candidate->readWithFewerHitsGenomeLocation, maxSpacing)) {
                    //
                    // Out of mate candidates.
                    //
                    break;
                }

                mateIndex--;
            } 
        }

        //
        // Remove us from the head of the list and proceed to the next candidate to score.
        //
        scoringCandidates[currentBestPossibleScoreList] = candidate->scoreListNext;
     }

doneScoring:

    if (bestPairScore == 65536) {
        //
        // Found nothing.
        //
        for (unsigned whichRead = 0; whichRead < NUM_READS_PER_PAIR; whichRead++) {
            result->location[whichRead] = InvalidGenomeLocation;
            result->mapq[whichRead] = 0;
            result->score[whichRead] = -1;
            result->status[whichRead] = NotFound;
#ifdef  _DEBUG
            if (_DumpAlignments) {
                printf("No sufficiently good pairs found.\n");
            }
#endif  // DEBUG
        }
    } else {
        for (unsigned whichRead = 0; whichRead < NUM_READS_PER_PAIR; whichRead++) {
            result->location[whichRead] = bestResultGenomeLocation[whichRead];
            result->direction[whichRead] = bestResultDirection[whichRead];
            result->mapq[whichRead] = computeMAPQ(probabilityOfAllPairs, probabilityOfBestPair, bestResultScore[whichRead], popularSeedsSkipped[0] + popularSeedsSkipped[1]);
            result->status[whichRead] = result->mapq[whichRead] > 10 ? SingleHit : MultipleHits;
            result->score[whichRead] = bestResultScore[whichRead];
        }
#ifdef  _DEBUG
            if (_DumpAlignments) {
                printf("Returned %u %s %u %s with MAPQ %d and %d, probability of all pairs %e, probability of best pair %e\n", 
                    result->location[0], result->direction[0] == RC ? "RC" : "", result->location[1], result->direction[1] == RC ? "RC" : "", result->mapq[0], result->mapq[1],
                    probabilityOfAllPairs, probabilityOfBestPair);
            }
#endif  // DEBUG
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
    double              *matchProbability,
    int                 *genomeLocationOffset)
{
    nLocationsScored++;

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

    _ASSERT(!memcmp(data+seedOffset, readToScore->getData() + seedOffset, seedLen));    // that the seed actually matches

    // NB: This cacheKey computation MUST match the one in BaseAligner or all hell will break loose.
    _uint64 cacheKey = (genomeLocation + tailStart) | (((_uint64) direction) << 32) | (((_uint64) whichRead) << 33) | (((_uint64)tailStart) << 34);

    score1 = landauVishkin->computeEditDistance(data + tailStart, genomeDataLength - tailStart, readToScore->getData() + tailStart, readToScore->getQuality() + tailStart, readLen - tailStart,
        scoreLimit, &matchProb1, cacheKey);
    if (score1 == -1) {
        *score = -1;
    } else {
        // The tail of the read matched; now let's reverse the reference genome data and match the head
        int limitLeft = scoreLimit - score1;
        score2 = reverseLandauVishkin->computeEditDistance(data + seedOffset, seedOffset + MAX_K, reversedRead[whichRead][direction] + readLen - seedOffset,
                                                                    reads[whichRead][OppositeDirection(direction)]->getQuality() + readLen - seedOffset, seedOffset, limitLeft, &matchProb2, cacheKey, genomeLocationOffset);

        if (score2 == -1) {
            *score = -1;
        } else {
            *score = score1 + score2;
            _ASSERT(*score <= scoreLimit);
            // Map probabilities for substrings can be multiplied, but make sure to count seed too
            *matchProbability = matchProb1 * matchProb2 * pow(1 - SNP_PROB, seedLen);
        }
    }

    if (*score == -1) {
        *matchProbability = 0;
    }
}

    void
 IntersectingPairedEndAligner::HashTableHitSet::firstInit(unsigned maxSeeds_, unsigned maxMergeDistance_) 
 {
    maxSeeds = maxSeeds_;
    maxMergeDistance = maxMergeDistance_;
    nLookupsUsed = 0;
    lookups = new HashTableLookup[maxSeeds];
 }
    void 
IntersectingPairedEndAligner::HashTableHitSet::init()
{
    nLookupsUsed = 0;
}
    void 
IntersectingPairedEndAligner::HashTableHitSet::recordLookup(unsigned seedOffset, unsigned nHits, const unsigned *hits, bool beginsDisjointHitSet)
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
	lookups[nLookupsUsed].beginsDisjointHitSet = beginsDisjointHitSet;

    //
    // Trim off any hits that are smaller than seedOffset, since they are clearly meaningless.
    //
    while (lookups[nLookupsUsed].nHits > 0 && lookups[nLookupsUsed].hits[lookups[nLookupsUsed].nHits - 1] < lookups[nLookupsUsed].seedOffset) {
        lookups[nLookupsUsed].nHits--;
    }

    nLookupsUsed++;
}

	unsigned 
IntersectingPairedEndAligner::HashTableHitSet::computeBestPossibleScoreForCurrentHit()
{
 	//
	// Now compute the best possible score for the hit.  This is the largest number of misses in any disjoint hit set.
	//
	unsigned bestPossibleScoreSoFar = 0;
	unsigned bestPossibleScoreForThisHitSet = 0;
	for (unsigned i = 0; i < nLookupsUsed; i++) {
		if (lookups[i].beginsDisjointHitSet) {
			bestPossibleScoreSoFar = max(bestPossibleScoreSoFar, bestPossibleScoreForThisHitSet);
			bestPossibleScoreForThisHitSet = 0;
		}
		if (!(lookups[i].currentHitForIntersection != lookups[i].nHits && 
				isWithin(lookups[i].hits[lookups[i].currentHitForIntersection], mostRecentLocationReturned + lookups[i].seedOffset,  maxMergeDistance) ||	
			lookups[i].currentHitForIntersection != 0 &&
				isWithin(lookups[i].hits[lookups[i].currentHitForIntersection-1], mostRecentLocationReturned + lookups[i].seedOffset,  maxMergeDistance))) {
			//
			// This one was not close enough.
			//
			bestPossibleScoreForThisHitSet++;
		}
	}

	return bestPossibleScoreSoFar;
}

	bool    
IntersectingPairedEndAligner::HashTableHitSet::getNextHitLessThanOrEqualTo(unsigned maxGenomeOffsetToFind, unsigned *actualGenomeOffsetFound, unsigned *seedOffsetFound, unsigned *bestPossibleScore)
{
#if     0
    //
    // Iterate through all of the lookups and search for the best candidate.  Instead of doing the obvious thing and binary searching each hit set in turn,
    // we interleave them.  This is so that we can do cache prefetches for the index locations and then do computation on the other hit sets while
    // the prefetches execute.
    //
    int limit[MAX_MAX_SEEDS][2];
    bool done[MAX_MAX_SEEDS];
    unsigned maxGenomeOffsetToFindThisSeed[MAX_MAX_SEEDS];
    unsigned nRemaining = nLookupsUsed;
    unsigned liveHitSets[MAX_MAX_SEEDS];

    for (unsigned i = 0; i < nLookupsUsed; i++) {
        limit[i][0] = (int)lookups[i].currentHitForIntersection;
        limit[i][1] = (int)lookups[i].nHits - 1;
        maxGenomeOffsetToFindThisSeed[i] = maxGenomeOffsetToFind - lookups[i].seedOffset;
        liveHitSets[i] = i;

        if (doAlignerPrefetch && 0) {
            //
            // Prefetch the first location we'll look at now.
            //
            _mm_prefetch((const char *)&lookups[i].hits[(limit[i][0] + limit[i][1])/2], _MM_HINT_T2);
        }
        done[i] = false;
    }

    bool anyFound = false;
    unsigned bestOffsetFound = 0;
    unsigned liveHitSetIndex = 0;
    while (nRemaining > 0) {    // Actually, we always exit out of the middle of the loop...  This could be for(;;), but the while conveys the idea more clearly.
        unsigned whichHitSet = liveHitSets[liveHitSetIndex];
        _ASSERT(!done[whichHitSet]);

        int probe = (limit[whichHitSet][0] + limit[whichHitSet][1]) / 2;
        //
        // Recall that the hit sets are sorted from largest to smallest, so the strange looking logic is actually right.
        // We're evaluating the expression "lookups[i].hits[probe] <= maxGenomeOffsetToFindThisSeed && (probe == 0 || lookups[i].hits[probe-1] > maxGenomeOffsetToFindThisSeed)"
        // It's written in this strange way just so the profile tool will show us where the time's going.
        //
        unsigned clause1 = lookups[whichHitSet].hits[probe] <= maxGenomeOffsetToFindThisSeed[whichHitSet];
        unsigned clause2 = probe == 0;

        if (clause1 && (clause2 || lookups[whichHitSet].hits[probe-1] > maxGenomeOffsetToFindThisSeed[whichHitSet])) {
            anyFound = true;
            mostRecentLocationReturned = *actualGenomeOffsetFound = bestOffsetFound = __max(lookups[whichHitSet].hits[probe] - lookups[whichHitSet].seedOffset, bestOffsetFound);
            *seedOffsetFound = lookups[whichHitSet].seedOffset;
            lookups[whichHitSet].currentHitForIntersection = probe;
            done[whichHitSet] = true;
            nRemaining--;
            if (0 == nRemaining) {
                break;
            }
            //
            // Swap us out of the live hit sets.
            //
            liveHitSets[liveHitSetIndex] = liveHitSets[nRemaining]; // Recall that we've already decremented nRemaining.
            liveHitSetIndex = (liveHitSetIndex + 1)%nRemaining;
            continue;
        }

        if (lookups[whichHitSet].hits[probe] > maxGenomeOffsetToFindThisSeed[whichHitSet]) {   // Recode this without the if to avoid the hard-to-predict branch.
            limit[whichHitSet][0] = probe + 1;
        } else {
            limit[whichHitSet][1] = probe - 1;
        }

        if (limit[whichHitSet][0] > limit[whichHitSet][1]) {
            lookups[whichHitSet].currentHitForIntersection = lookups[whichHitSet].nHits;    // We're done with this lookup.
            done[whichHitSet] = true;
            nRemaining--;
            if (0 == nRemaining) {
                break;
            }
            //
            // Swap us out of the live hit sets.
            //
            liveHitSets[liveHitSetIndex] = liveHitSets[nRemaining]; // Recall that we've already decremented nRemaining.
            liveHitSetIndex = (liveHitSetIndex + 1)%nRemaining;
            continue;
        }
        //
        // Launch a prefetch for the next hit we'll look at in this hit set.
        //
        if (doAlignerPrefetch && 0) {
            _mm_prefetch((const char *)&lookups[whichHitSet].hits[(limit[whichHitSet][0] + limit[whichHitSet][1])/2], _MM_HINT_T2);
        }

        //
        // And proceed on to the next hit set.
        //
       liveHitSetIndex = (liveHitSetIndex + 1) % nRemaining;
    }

#else     // The old way
    bool anyFound = false;
    unsigned bestOffsetFound = 0;
    for (unsigned i = 0; i < nLookupsUsed; i++) {
        //
        // Binary search from the current starting offset to either the right place or the end.
        //
        int limit[2] = {(int)lookups[i].currentHitForIntersection, (int)lookups[i].nHits - 1};
        unsigned maxGenomeOffsetToFindThisSeed = maxGenomeOffsetToFind + lookups[i].seedOffset;
        while (limit[0] <= limit[1]) {
            unsigned probe = (limit[0] + limit[1]) / 2;
            //
            // Recall that the hit sets are sorted from largest to smallest, so the strange looking logic is actually right.
            // We're evaluating the expression "lookups[i].hits[probe] <= maxGenomeOffsetToFindThisSeed && (probe == 0 || lookups[i].hits[probe-1] > maxGenomeOffsetToFindThisSeed)"
            // It's written in this strange way just so the profile tool will show us where the time's going.
            //
            unsigned clause1 = lookups[i].hits[probe] <= maxGenomeOffsetToFindThisSeed;
            unsigned clause2 = probe == 0;

            if (clause1 && (clause2 || lookups[i].hits[probe-1] > maxGenomeOffsetToFindThisSeed)) {
                anyFound = true;
                if (lookups[i].hits[probe] - lookups[i].seedOffset >  bestOffsetFound) {
                    mostRecentLocationReturned = *actualGenomeOffsetFound = bestOffsetFound = lookups[i].hits[probe] - lookups[i].seedOffset;
                    *seedOffsetFound = lookups[i].seedOffset;
                }
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
#endif  // 0

	if (anyFound) {
		*bestPossibleScore = computeBestPossibleScoreForCurrentHit();
        _ASSERT(*actualGenomeOffsetFound <= maxGenomeOffsetToFind);
     }

    return anyFound;
}


    bool    
IntersectingPairedEndAligner::HashTableHitSet::getFirstHit(unsigned *genomeLocation, unsigned *seedOffsetFound, unsigned *bestPossibleScore)
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

	if (anyFound) {
		*bestPossibleScore = computeBestPossibleScoreForCurrentHit();
	}

	return anyFound;
}

    bool    
IntersectingPairedEndAligner::HashTableHitSet::getNextLowerHit(unsigned *genomeLocation, unsigned *seedOffsetFound, unsigned *bestPossibleScore)
{
    //
    // Look through all of the lookups and find the one with the highest location smaller than the current one.
    //
    unsigned foundLocation = 0;
    bool anyFound = false;
 
    //
    // Run through the lookups pushing up any that are at the most recently returned 
    //
    for (unsigned i = 0; i < nLookupsUsed; i++) {
        _ASSERT(lookups[i].currentHitForIntersection == lookups[i].nHits || lookups[i].hits[lookups[i].currentHitForIntersection] - lookups[i].seedOffset <= mostRecentLocationReturned || 
            lookups[i].hits[lookups[i].currentHitForIntersection] < lookups[i].seedOffset);

        if (lookups[i].currentHitForIntersection != lookups[i].nHits && lookups[i].hits[lookups[i].currentHitForIntersection] - lookups[i].seedOffset == mostRecentLocationReturned) {
            lookups[i].currentHitForIntersection++;
        }

        if (lookups[i].currentHitForIntersection != lookups[i].nHits) {
            if (foundLocation < lookups[i].hits[lookups[i].currentHitForIntersection] - lookups[i].seedOffset && // found location is OK
                lookups[i].hits[lookups[i].currentHitForIntersection] >= lookups[i].seedOffset) // found location isn't too small to push us before the beginning of the genome
            { 
                *genomeLocation = foundLocation = lookups[i].hits[lookups[i].currentHitForIntersection] - lookups[i].seedOffset;
                *seedOffsetFound = lookups[i].seedOffset;
                anyFound = true;
            }
        }
    }
 
    if (anyFound) {
        mostRecentLocationReturned = foundLocation;
		*bestPossibleScore = computeBestPossibleScoreForCurrentHit();
    }

    return anyFound;
}

            bool 
IntersectingPairedEndAligner::MergeAnchor::checkMerge(unsigned newMoreHitLocation, unsigned newFewerHitLocation, double newMatchProbability, int newPairScore, 
                        double *oldMatchProbability) 
{
    if (locationForReadWithMoreHits == InvalidGenomeLocation || !doesRangeMatch(newMoreHitLocation, newFewerHitLocation)) {
        //
        // No merge.  Remember the new one.
        //
        locationForReadWithMoreHits = newMoreHitLocation;
        locationForReadWithFewerHits = newFewerHitLocation;
        matchProbability = newMatchProbability;
        pairScore = newPairScore;
        *oldMatchProbability = 0.0;
        return false;
    }  else {
        //
        // Within merge distance.  Keep the better score (or if they're tied the better match probability).
        //
        if (newPairScore < pairScore || newPairScore == pairScore && newMatchProbability > matchProbability) {
#ifdef _DEBUG
            if (_DumpAlignments) {
                printf("Merge replacement at anchor (%u, %u), loc (%u, %u), old match prob %e, new match prob %e, old pair score %d, new pair score %d\n",
                    locationForReadWithMoreHits, locationForReadWithFewerHits, newMoreHitLocation, newFewerHitLocation,
                    matchProbability, newMatchProbability, pairScore, newPairScore);
            }
#endif // DEBUG

            *oldMatchProbability = matchProbability;
            matchProbability = newMatchProbability;
            pairScore = newPairScore;
            return false;
        } else {
            //
            // The new one should just be ignored.
            //
#ifdef _DEBUG
            if (_DumpAlignments) {
                printf("Merged at anchor (%u, %u), loc (%u, %u), old match prob %e, new match prob %e, old pair score %d, new pair score %d\n",
                    locationForReadWithMoreHits, locationForReadWithFewerHits, newMoreHitLocation, newFewerHitLocation,
                    matchProbability, newMatchProbability, pairScore, newPairScore);
            }
#endif // DEBUG
            return true;
        }
    }

    _ASSERT(!"NOTREACHED");
}
