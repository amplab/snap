/*++

Module Name:

	TenXSingleAligner.cpp

Abstract:

	A paired-end aligner based on set intersections to narrow down possible candidate locations.

Authors:

	Hongyi Xin and Bill Bolosky, June, 2016

Environment:

	User mode service.

Revision History:

--*/

#include "stdafx.h"
#include "TenXSingleAligner.h"
#include "SeedSequencer.h"
#include "mapq.h"
#include "exit.h"
#include "Error.h"
#include "BigAlloc.h"
#include "AlignerOptions.h"
#include <iostream>

#ifdef  _DEBUG
extern bool _DumpAlignments;    // From BaseAligner.cpp
#endif  // _DEBUG

TenXSingleAligner::TenXSingleAligner(
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
	unsigned      maxCandidatePoolSize,
	int           maxSecondaryAlignmentsPerContig_,
	BigAllocator  *allocator,
	bool          noUkkonen_,
	bool          noOrderedEvaluation_,
	bool          noTruncation_,
	bool          ignoreAlignmentAdjustmentsForOm_,
	unsigned      printStatsMapQLimit_) :
	index(index_), maxReadSize(maxReadSize_), maxHits(maxHits_), maxK(maxK_), numSeedsFromCommandLine(__min(MAX_MAX_SEEDS, numSeedsFromCommandLine_)), minSpacing(minSpacing_), maxSpacing(maxSpacing_),
	landauVishkin(NULL), reverseLandauVishkin(NULL), maxBigHits(maxBigHits_), seedCoverage(seedCoverage_),
	extraSearchDepth(extraSearchDepth_), nLocationsScored(0), noUkkonen(noUkkonen_), noOrderedEvaluation(noOrderedEvaluation_), noTruncation(noTruncation_),
	maxSecondaryAlignmentsPerContig(maxSecondaryAlignmentsPerContig_), alignmentAdjuster(index->getGenome()), ignoreAlignmentAdjustmentsForOm(ignoreAlignmentAdjustmentsForOm_), printStatsMapQLimit(printStatsMapQLimit_)
{
	doesGenomeIndexHave64BitLocations = index->doesGenomeIndexHave64BitLocations();

	unsigned maxSeedsToUse;
	if (0 != numSeedsFromCommandLine) {
		maxSeedsToUse = numSeedsFromCommandLine;
	}
	else {
		maxSeedsToUse = (unsigned)(maxReadSize * seedCoverage / index->getSeedLength());
	}
	allocateDynamicMemory(allocator, maxReadSize, maxBigHits, maxSeedsToUse, maxK, extraSearchDepth, maxCandidatePoolSize, maxSecondaryAlignmentsPerContig);

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

	setPair[0][0] = hashTableHitSets[0][FORWARD];
	setPair[0][1] = hashTableHitSets[1][RC];
	setPair[1][0] = hashTableHitSets[0][RC];
	setPair[1][1] = hashTableHitSets[1][FORWARD];
}

TenXSingleAligner::~TenXSingleAligner()
{
}

size_t
TenXSingleAligner::getBigAllocatorReservation(GenomeIndex * index, unsigned maxBigHitsToConsider, unsigned maxReadSize, unsigned seedLen, unsigned numSeedsFromCommandLine,
	double seedCoverage, unsigned maxEditDistanceToConsider, unsigned maxExtraSearchDepth, unsigned maxCandidatePoolSize,
	int maxSecondaryAlignmentsPerContig)
{
	unsigned maxSeedsToUse;
	if (0 != numSeedsFromCommandLine) {
		maxSeedsToUse = numSeedsFromCommandLine;
	}
	else {
		maxSeedsToUse = (unsigned)(maxReadSize * seedCoverage / index->getSeedLength());
	}
	CountingBigAllocator countingAllocator;
	{
		TenXSingleAligner aligner; // This has to be in a nested scope so its destructor is called before that of the countingAllocator
		aligner.index = index;
		aligner.doesGenomeIndexHave64BitLocations = index->doesGenomeIndexHave64BitLocations();

		aligner.allocateDynamicMemory(&countingAllocator, maxReadSize, maxBigHitsToConsider, maxSeedsToUse, maxEditDistanceToConsider, maxExtraSearchDepth, maxCandidatePoolSize,
			maxSecondaryAlignmentsPerContig);

		fprintf(stderr, "****sizeof(aligner): %lld  getMemoryUsed(): %lld\n", sizeof(aligner), countingAllocator.getMemoryUsed());

		return sizeof(aligner) + countingAllocator.getMemoryUsed();
	}
}

void
TenXSingleAligner::allocateDynamicMemory(BigAllocator *allocator, unsigned maxReadSize, unsigned maxBigHitsToConsider, unsigned maxSeedsToUse,
	unsigned maxEditDistanceToConsider, unsigned maxExtraSearchDepth, unsigned maxCandidatePoolSize,
	int maxSecondaryAlignmentsPerContig)
{
	CountingBigAllocator *allocatorCast = (CountingBigAllocator*) allocator;
	seedUsed = (BYTE *)allocator->allocate(100 + (maxReadSize + 7) / 8);

	for (unsigned whichRead = 0; whichRead < NUM_READS_PER_PAIR; whichRead++) {
		rcReadData[whichRead] = (char *)allocator->allocate(maxReadSize);
		rcReadQuality[whichRead] = (char *)allocator->allocate(maxReadSize);

		for (Direction dir = 0; dir < NUM_DIRECTIONS; dir++) {
			reversedRead[whichRead][dir] = (char *)allocator->allocate(maxReadSize);
			hashTableHitSets[whichRead][dir] = (HashTableHitSet *)allocator->allocate(sizeof(HashTableHitSet)); /*new HashTableHitSet();*/
			hashTableHitSets[whichRead][dir]->firstInit(maxSeedsToUse, maxMergeDistance, allocator, doesGenomeIndexHave64BitLocations);
		}
	}

	//fprintf(stderr, "**stage 1 getMemoryUsed(): %lld\n", allocatorCast->getMemoryUsed());

	scoringCandidatePoolSize = min(maxCandidatePoolSize, maxBigHitsToConsider * maxSeedsToUse * NUM_READS_PER_PAIR);

	scoringCandidates = (ScoringCandidate **)allocator->allocate(sizeof(ScoringCandidate *) * (maxEditDistanceToConsider + maxExtraSearchDepth + 1));  //+1 is for 0.
    probabilityForED = (double*) allocator->allocate(sizeof(double) * (maxEditDistanceToConsider + maxExtraSearchDepth + 1));  //+1 is for 0.
	
	//fprintf(stderr, "**stage 1.1 getMemoryUsed(): %lld\n", allocatorCast->getMemoryUsed());
	
	scoringCandidatePool = (ScoringCandidate *)allocator->allocate(sizeof(ScoringCandidate) * scoringCandidatePoolSize);
	
	//fprintf(stderr, "**stage 1.2 getMemoryUsed(): %lld\n", allocatorCast->getMemoryUsed());

	for (unsigned i = 0; i < NUM_READS_PER_PAIR; i++) {
		scoringMateCandidates[i] = (ScoringMateCandidate *)allocator->allocate(sizeof(ScoringMateCandidate) * scoringCandidatePoolSize / NUM_READS_PER_PAIR);
	}

	//fprintf(stderr, "**stage 2 getMemoryUsed(): %lld\n", allocatorCast->getMemoryUsed());

	mergeAnchorPoolSize = scoringCandidatePoolSize;
	mergeAnchorPool = (MergeAnchor *)allocator->allocate(sizeof(MergeAnchor) * mergeAnchorPoolSize);
	
	//fprintf(stderr, "**stage 2.1 getMemoryUsed(): %lld\n", allocatorCast->getMemoryUsed());

	if (maxSecondaryAlignmentsPerContig > 0) {
		size_t size = sizeof(*hitsPerContigCounts) * index->getGenome()->getNumContigs();
		hitsPerContigCounts = (HitsPerContigCounts *)allocator->allocate(size);
		memset(hitsPerContigCounts, 0, size);
		contigCountEpoch = 0;
	}
	else {
		hitsPerContigCounts = NULL;
	}

	//fprintf(stderr, "**stage 3 getMemoryUsed(): %lld\n", allocatorCast->getMemoryUsed());

}

bool
TenXSingleAligner::align_phase_1(Read* read0, Read* read1, unsigned *popularSeedsSkipped)
{
	int maxSeeds;
	if (numSeedsFromCommandLine != 0) {
		maxSeeds = (int)numSeedsFromCommandLine;
	}
	else {
		maxSeeds = (int)(max(read0->getDataLength(), read1->getDataLength()) * seedCoverage / index->getSeedLength());
	}

#ifdef  _DEBUG
	if (_DumpAlignments) {
		printf("\nIntersectingAligner aligning reads '%*.s' and '%.*s' with data '%.*s' and '%.*s'\n", read0->getIdLength(), read0->getId(), read1->getIdLength(), read1->getId(), read0->getDataLength(), read0->getData(), read1->getDataLength(), read1->getData());
	}
#endif  // _DEBUG

	lowestFreeScoringCandidatePoolEntry = 0;
	for (unsigned k = 0; k <= maxK + extraSearchDepth; k++) {
		scoringCandidates[k] = NULL;
        probabilityForED[k] = 0;
	}

	for (unsigned i = 0; i < NUM_SET_PAIRS; i++) {
		lowestFreeScoringMateCandidate[i] = 0;
	}
	firstFreeMergeAnchor = 0;

	reads[0][FORWARD] = read0;
	reads[1][FORWARD] = read1;

	//
	// Don't bother if one or both reads are too short.  The minimum read length here is the seed length, but usually there's a longer
	// minimum enforced by our called
	//
	if (read0->getDataLength() < seedLen || read1->getDataLength() < seedLen) {
		return true;
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
			WriteErrorMessage("TenXSingleAligner:: got too big read (%d > %d)\n"
				"Change MAX_READ_LENTH at the beginning of Read.h and recompile.\n", readLen[whichRead], maxReadSize);
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
		return true;
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

	unsigned thisPassSeedsNotSkipped[NUM_READS_PER_PAIR][NUM_DIRECTIONS] = { {0,0}, {0,0} };

	//
	// Phase 1: do the hash table lookups for each of the seeds for each of the reads and add them to the hit sets.
	//
	for (unsigned whichRead = 0; whichRead < NUM_READS_PER_PAIR; whichRead++) {
		int nextSeedToTest = 0;
		unsigned wrapCount = 0;
		int nPossibleSeeds = (int)readLen[whichRead] - seedLen + 1;
		memset(seedUsed, 0, (__max(readLen[0], readLen[1]) + 7) / 8);
		bool beginsDisjointHitSet[NUM_DIRECTIONS] = { true, true };

		while (countOfHashTableLookups[whichRead] < nPossibleSeeds && countOfHashTableLookups[whichRead] < maxSeeds) {
			if (nextSeedToTest >= nPossibleSeeds) {
				wrapCount++;
				beginsDisjointHitSet[FORWARD] = beginsDisjointHitSet[RC] = true;
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
			_int64 nHits[NUM_DIRECTIONS];
			const GenomeLocation *hits[NUM_DIRECTIONS];
			const unsigned *hits32[NUM_DIRECTIONS];

			if (doesGenomeIndexHave64BitLocations) {
				index->lookupSeed(seed, &nHits[FORWARD], &hits[FORWARD], &nHits[RC], &hits[RC],
					hashTableHitSets[whichRead][FORWARD]->getNextSingletonLocation(), hashTableHitSets[whichRead][RC]->getNextSingletonLocation());
			}
			else {
				index->lookupSeed32(seed, &nHits[FORWARD], &hits32[FORWARD], &nHits[RC], &hits32[RC]);
			}

			countOfHashTableLookups[whichRead]++;
			for (Direction dir = FORWARD; dir < NUM_DIRECTIONS; dir++) {
				int offset;
				if (dir == FORWARD) {
					offset = nextSeedToTest;
				}
				else {
					offset = readLen[whichRead] - seedLen - nextSeedToTest;
				}
				if (nHits[dir] < maxBigHits) {
					totalHashTableHits[whichRead][dir] += nHits[dir];
					if (doesGenomeIndexHave64BitLocations) {
						hashTableHitSets[whichRead][dir]->recordLookup(offset, nHits[dir], hits[dir], beginsDisjointHitSet[dir]);
					}
					else {
						hashTableHitSets[whichRead][dir]->recordLookup(offset, nHits[dir], hits32[dir], beginsDisjointHitSet[dir]);
					}
					beginsDisjointHitSet[dir] = false;
				}
				else {
					popularSeedsSkipped[whichRead]++;
				}
			}

			//
			// If we don't have enough seeds left to reach the end of the read, space out the seeds more-or-less evenly.
			//
			if ((maxSeeds - countOfHashTableLookups[whichRead] + 1) * (int)seedLen + nextSeedToTest < nPossibleSeeds) {
				_ASSERT((nPossibleSeeds - nextSeedToTest - 1) / (maxSeeds - countOfHashTableLookups[whichRead] + 1) >= (int)seedLen);
				nextSeedToTest += (nPossibleSeeds - nextSeedToTest - 1) / (maxSeeds - countOfHashTableLookups[whichRead] + 1);
				_ASSERT(nextSeedToTest < nPossibleSeeds);   // We haven't run off the end of the read.
			}
			else {
				nextSeedToTest += seedLen;
			}
		} // while we need to lookup seeds for this read
	} // for each read

	readWithMoreHits = totalHashTableHits[0][FORWARD] + totalHashTableHits[0][RC] > totalHashTableHits[1][FORWARD] + totalHashTableHits[1][RC] ? 0 : 1;
	readWithFewerHits = 1 - readWithMoreHits;

#ifdef  _DEBUG
	if (_DumpAlignments) {
		printf("Read 0 has %d hits, read 1 has %d hits\n", totalHashTableHits[0][FORWARD] + totalHashTableHits[0][RC], totalHashTableHits[1][FORWARD] + totalHashTableHits[1][RC]);
	}
#endif  // _DEBUg

	// default, not exiting early
	return false;
}


int
TenXSingleAligner::align_phase_2_move_locus(unsigned whichSetPair)
{
	//
		// Loop invariant: lastGenomeLocationForReadWithFewerHits is the highest genome offset that has not been considered.
		// lastGenomeLocationForReadWithMoreHits is also the highest genome offset on that side that has not been
		// considered (or is InvalidGenomeLocation), but higher ones within the appropriate range might already be in scoringMateCandidates.
		// We go once through this loop for each
		//
		// default, not exiting early

	if (lastGenomeLocationForReadWithMoreHits[whichSetPair] > lastGenomeLocationForReadWithFewerHits[whichSetPair] + maxSpacing) {
		//
		// The more hits side is too high to be a mate candidate for the fewer hits side.  Move it down to the largest
		// location that's not too high.
		//
		if (!setPair[whichSetPair][readWithMoreHits]->getNextHitLessThanOrEqualTo(lastGenomeLocationForReadWithFewerHits[whichSetPair] + maxSpacing,
			&lastGenomeLocationForReadWithMoreHits[whichSetPair], &lastSeedOffsetForReadWithMoreHits[whichSetPair])) {
			return 1;  // End of all of the mates.  We're done with this set pair.
		}
	}

	//
	// Even though we are out of more hit locations, we might still backtrack!
	//
	if ((lastGenomeLocationForReadWithMoreHits[whichSetPair] + maxSpacing < lastGenomeLocationForReadWithFewerHits[whichSetPair] || outOfMoreHitsLocations[whichSetPair]) &&
		(0 == lowestFreeScoringMateCandidate[whichSetPair] ||
			!genomeLocationIsWithin(scoringMateCandidates[whichSetPair][lowestFreeScoringMateCandidate[whichSetPair] - 1].readWithMoreHitsGenomeLocation, lastGenomeLocationForReadWithFewerHits[whichSetPair], maxSpacing))) {
		//
		// No mates for the hit on the read with fewer hits.  Skip to the next candidate.
		//
		if (outOfMoreHitsLocations[whichSetPair]) {
			//
			// Nothing left on the more hits side, we're done with this set pair.
			//
			return 1;
		}

		if (!setPair[whichSetPair][readWithFewerHits]->getNextHitLessThanOrEqualTo(lastGenomeLocationForReadWithMoreHits[whichSetPair] + maxSpacing, &lastGenomeLocationForReadWithFewerHits[whichSetPair],
			&lastSeedOffsetForReadWithFewerHits[whichSetPair])) {
			//
			// No more candidates on the read with fewer hits side.  We're done with this set pair.
			//
			return 1;
		}
		return -1;
	}
	return 0;
}


bool
TenXSingleAligner::align_phase_2_single_step_add_candidate(unsigned whichSetPair, int clusterIdx, unsigned clusterEDCompensation)
{
	//
	// Add all of the mate candidates for this fewer side hit.
	//
	while (lastGenomeLocationForReadWithMoreHits[whichSetPair] + maxSpacing >= lastGenomeLocationForReadWithFewerHits[whichSetPair] && !outOfMoreHitsLocations[whichSetPair]) {
		unsigned bestPossibleScoreForReadWithMoreHits;
		if (noTruncation) {
			bestPossibleScoreForReadWithMoreHits = 0;
		}
		else {
			bestPossibleScoreForReadWithMoreHits = setPair[whichSetPair][readWithMoreHits]->computeBestPossibleScoreForCurrentHit();
		}

		if (lowestFreeScoringMateCandidate[whichSetPair] >= scoringCandidatePoolSize / NUM_READS_PER_PAIR) {
			WriteErrorMessage("Ran out of scoring candidate pool entries.  Perhaps trying with a larger value of -mcp will help.\n");
			soft_exit(1);
		}
		scoringMateCandidates[whichSetPair][lowestFreeScoringMateCandidate[whichSetPair]].init(
			lastGenomeLocationForReadWithMoreHits[whichSetPair], bestPossibleScoreForReadWithMoreHits, lastSeedOffsetForReadWithMoreHits[whichSetPair]);

#ifdef _DEBUG
		if (_DumpAlignments) {
			printf("SetPair %d, added more hits candidate %d at genome location %u, bestPossibleScore %d, seedOffset %d\n",
				whichSetPair, lowestFreeScoringMateCandidate[whichSetPair], lastGenomeLocationForReadWithMoreHits[whichSetPair].location,
				bestPossibleScoreForReadWithMoreHits,
				lastSeedOffsetForReadWithMoreHits[whichSetPair]);
		}
#endif // _DEBUG

		lowestFreeScoringMateCandidate[whichSetPair]++;

		if (!setPair[whichSetPair][readWithMoreHits]->getNextLowerHit(&lastGenomeLocationForReadWithMoreHits[whichSetPair], &lastSeedOffsetForReadWithMoreHits[whichSetPair])) {
			lastGenomeLocationForReadWithMoreHits[whichSetPair] = 0;
			outOfMoreHitsLocations[whichSetPair] = true;
			break; // out of the loop looking for candidates on the more hits side.
		}
	}

	//
	// And finally add the hit from the fewer hit side.  To compute its best possible score, we need to look at all of the mates; we couldn't do it in the
	// loop immediately above because some of them might have already been in the mate list from a different, nearby fewer hit location.
	//
	unsigned bestPossibleScoreForReadWithFewerHits;

	if (noTruncation) {
		bestPossibleScoreForReadWithFewerHits = 0;
	}
	else {
		bestPossibleScoreForReadWithFewerHits = setPair[whichSetPair][readWithFewerHits]->computeBestPossibleScoreForCurrentHit();
	}

	unsigned lowestBestPossibleScoreOfAnyPossibleMate = maxK + extraSearchDepth;
	for (int i = lowestFreeScoringMateCandidate[whichSetPair] - 1; i >= 0; i--) {
		if (scoringMateCandidates[whichSetPair][i].readWithMoreHitsGenomeLocation > lastGenomeLocationForReadWithFewerHits[whichSetPair] + maxSpacing) {
			break;
		}
		lowestBestPossibleScoreOfAnyPossibleMate = __min(lowestBestPossibleScoreOfAnyPossibleMate, scoringMateCandidates[whichSetPair][i].bestPossibleScore);
	}

	if (lowestBestPossibleScoreOfAnyPossibleMate + bestPossibleScoreForReadWithFewerHits <= maxK + extraSearchDepth) {
		//
		// There's a set of ends that we can't prove doesn't have too large of a score.  Allocate a fewer hit candidate and stick it in the
		// correct weight list.
		//
		if (lowestFreeScoringCandidatePoolEntry >= scoringCandidatePoolSize) {
			WriteErrorMessage("Ran out of scoring candidate pool entries.  Perhaps rerunning with a larger value of -mcp will help.\n");
			soft_exit(1);
		}

		// Add 10X cluster penalty
		unsigned clusterScorePenalty = clusterIdx == -1 ? clusterEDCompensation : 0;

		//
		// If we have noOrderedEvaluation set, just stick everything on list 0, regardless of what it really is.  This will cause us to
		// evaluate the candidates in more-or-less inverse genome order.
		//
		unsigned bestPossibleScore = noOrderedEvaluation ? 0 : lowestBestPossibleScoreOfAnyPossibleMate + bestPossibleScoreForReadWithFewerHits + clusterScorePenalty;

		scoringCandidatePool[lowestFreeScoringCandidatePoolEntry].init(lastGenomeLocationForReadWithFewerHits[whichSetPair], whichSetPair, lowestFreeScoringMateCandidate[whichSetPair] - 1,
			lastSeedOffsetForReadWithFewerHits[whichSetPair], bestPossibleScoreForReadWithFewerHits,
			scoringCandidates[bestPossibleScore], clusterIdx);


		scoringCandidates[bestPossibleScore] = &scoringCandidatePool[lowestFreeScoringCandidatePoolEntry];

#ifdef _DEBUG
		if (_DumpAlignments) {
			printf("SetPair %d, added fewer hits candidate %d at genome location %u, bestPossibleScore %d, seedOffset %d\n",
				whichSetPair, lowestFreeScoringCandidatePoolEntry, lastGenomeLocationForReadWithFewerHits[whichSetPair].location,
				lowestBestPossibleScoreOfAnyPossibleMate + bestPossibleScoreForReadWithFewerHits,
				lastSeedOffsetForReadWithFewerHits[whichSetPair]);
		}
#endif // _DEBUG

		lowestFreeScoringCandidatePoolEntry++;
		maxUsedBestPossibleScoreList = max(maxUsedBestPossibleScoreList, bestPossibleScore);
	}

	if (!setPair[whichSetPair][readWithFewerHits]->getNextLowerHit(&lastGenomeLocationForReadWithFewerHits[whichSetPair], &lastSeedOffsetForReadWithFewerHits[whichSetPair])) {
		return true;
	}

#ifdef  _DEBUG
	if (_DumpAlignments) {
		printf("Stepping function is working alright\n");
	}
#endif

	return false;
}


/*
bool
TenXSingleAligner::align_phase_2_single_step(unsigned whichSetPair)
{
	int check_range_result = align_phase_2_move_locus(whichSetPair);
	if (check_range_result == 1)
		return true;
	else if (check_range_result == -1)
		return false;
	else 
		return align_phase_2_single_step_add_candidate(whichSetPair, -1, 0);
	//return align_phase_2_move_locus(setPair, whichSetPair, outOfMoreHitsLocations, lastSeedOffsetForReadWithFewerHits, lastGenomeLocationForReadWithFewerHits, lastSeedOffsetForReadWithMoreHits, lastGenomeLocationForReadWithMoreHits, maxUsedBestPossibleScoreList, NULL);
}
*/


bool
TenXSingleAligner::align_phase_2_to_target_loc(const GenomeLocation &clusterTargetLoc, int clusterIdx, unsigned clusterEDCompensation)
{
	bool keepGoing = true;
	bool targetNotMet = false;
	bool targetNotMetSingleSet;

	for (int whichSetPair = 0; whichSetPair < NUM_DIRECTIONS; whichSetPair++)
	{
		//fprintf(stderr, "beginning: targetLoc: %lld, ReadLoc: %lld\n", clusterTargetLoc.location, lastGenomeLocationForReadWithFewerHits[whichSetPair].location);
		if (!noMoreLoci[whichSetPair]) {
			targetNotMetSingleSet = lastGenomeLocationForReadWithFewerHits[whichSetPair] > clusterTargetLoc;
			targetNotMet = targetNotMet || targetNotMetSingleSet;
		}
	}

//) {//
	while (keepGoing && targetNotMet) {
		keepGoing = false;
		for (int whichSetPair = 0; whichSetPair < NUM_DIRECTIONS; whichSetPair++) {
			//std::cout << "setPair[" << whichSetPair << "]" << setPair[whichSetPair] << std::endl;
			//printf("setPair[%d]: %p\n", whichSetPair, setPair[whichSetPair]);
			if (!noMoreLoci[whichSetPair]) {
				int check_range_result = align_phase_2_move_locus(whichSetPair);
				if (check_range_result == 1) {
					noMoreLoci[whichSetPair] = true;
				}
				else {
					noMoreLoci[whichSetPair] = false;
					//
					// keep going until we have a good loc pair
					//
					if (check_range_result == -1) {
						keepGoing = keepGoing || !noMoreLoci[whichSetPair];
						continue;
					}
				}

				if (!noMoreLoci[whichSetPair]) {
#ifdef _DEBUG
					if (_DumpAlignments) {
						printf("Pair: %d  beginning: targetLoc: %lld, ReadLoc: %lld\n", whichSetPair, clusterTargetLoc.location, lastGenomeLocationForReadWithFewerHits[whichSetPair].location);
					}
#endif
					targetNotMetSingleSet = lastGenomeLocationForReadWithFewerHits[whichSetPair] > clusterTargetLoc;
					targetNotMet = targetNotMet || targetNotMetSingleSet;
					if (targetNotMetSingleSet) {
#ifdef _DEBUG
						if (_DumpAlignments) {
							printf("Pair: %d  targetNotMetSingleSet: %s\n", whichSetPair, (targetNotMetSingleSet ? "true" : "false"));
						}
#endif //_DEBUG
						noMoreLoci[whichSetPair] = align_phase_2_single_step_add_candidate(whichSetPair, clusterIdx, clusterEDCompensation);
						//
						// We keep working on the loop as long as one set is still not stopped
						//
						keepGoing = keepGoing || !noMoreLoci[whichSetPair];
					}
				} // add_candidate
			} // check_range and add_candidate
		} // whichSetPair
	} // while

	return keepGoing;
}


GenomeLocation* TenXSingleAligner::align_phase_2_get_loci()
{
	GenomeLocation* nextLoci = NULL;
	for (unsigned direction = 0; direction < NUM_DIRECTIONS; direction++) {
		if (!noMoreLoci[direction] && (nextLoci == NULL || *nextLoci < lastGenomeLocationForReadWithFewerHits[direction]) )
			nextLoci = &lastGenomeLocationForReadWithFewerHits[direction];
	}
	return nextLoci;
}


bool TenXSingleAligner::align_phase_2_init()
{
	bool keepGoing = false;
	maxUsedBestPossibleScoreList = 0;

	//
	// Initialize variables
	//
	for (int whichSetPair = 0; whichSetPair < NUM_DIRECTIONS; whichSetPair++)
	{
		lastGenomeLocationForReadWithMoreHits[whichSetPair] = InvalidGenomeLocation;
		outOfMoreHitsLocations[whichSetPair] = false;

		if (setPair[whichSetPair][readWithFewerHits]->getFirstHit(&lastGenomeLocationForReadWithFewerHits[whichSetPair], &lastSeedOffsetForReadWithFewerHits[whichSetPair]))
			//
			// No hits in this direction.
			//
			noMoreLoci[whichSetPair] = true;   // The outer loop over set pairs.
		else
			noMoreLoci[whichSetPair] = false;

		keepGoing = keepGoing || !noMoreLoci[whichSetPair];
	}
	return keepGoing;
}


void
TenXSingleAligner::align_phase_2()
{
	//
	// Phase 2: find all possible candidates and add them to candidate lists (for the reads with fewer and more hits).
	//

	//
	// Loop over the candidates in for the read with more hits.  At the top of the loop, we have a candidate but don't know if it has
	// a mate.  Each pass through the loop considers a single hit on the read with fewer hits.
	//
	if(align_phase_2_init() ) {
		GenomeLocation clusterTargetLoc = GenomeLocation(0000000000);
		align_phase_2_to_target_loc(clusterTargetLoc, -1, 0);
	}
}


bool
TenXSingleAligner::align_phase_3(int maxEditDistanceForSecondaryResults, _int64 secondaryResultBufferSize, _int64* nSecondaryResults, PairedAlignmentResult* secondaryResults, _int64 maxSecondaryResultsToReturn,
	unsigned &bestPairScore, GenomeLocation *bestResultGenomeLocation, Direction *bestResultDirection, double &probabilityOfAllPairs, unsigned *bestResultScore, unsigned *popularSeedsSkipped, double &probabilityOfBestPair, int &bestClusterIdx, double unclusteredPenalty, unsigned clusterEDCompensation) //This is a hack. Pass in result by reference
{
	//
	// Phase 3: score and merge the candidates we've found.
	//
	//Initialize results
	*nSecondaryResults = 0;
	//
	// Initialize the member variables that are effectively stack locals, but are in the object
	// to avoid having to pass them to score.
	//
	localBestPairProbability[0] = 0;
	localBestPairProbability[1] = 0;
	Direction setPairDirection[NUM_SET_PAIRS][NUM_READS_PER_PAIR] = { {FORWARD, RC}, {RC, FORWARD} };

	unsigned currentBestPossibleScoreList = 0;
	unsigned scoreLimit = maxK + extraSearchDepth;
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

		//**** 10X surrogates
		unsigned	compensatedScoreLimit;
		double		probabilityCorrectionFactor;

		if (candidate->clusterIdx == -1)
			probabilityCorrectionFactor = unclusteredPenalty;
		else
			probabilityCorrectionFactor = 1;
		//**** 10X surrogates

		unsigned fewerEndScore;
		double fewerEndMatchProbability;
		int fewerEndGenomeLocationOffset;

		scoreLocation(readWithFewerHits, setPairDirection[candidate->whichSetPair][readWithFewerHits], candidate->readWithFewerHitsGenomeLocation,
			candidate->seedOffset, scoreLimit, &fewerEndScore, &fewerEndMatchProbability, &fewerEndGenomeLocationOffset);

		_ASSERT(-1 == fewerEndScore || fewerEndScore >= candidate->bestPossibleScore);

#ifdef _DEBUG
		if (_DumpAlignments) {
			printf("Scored fewer end candidate %d, set pair %d, read %d, location %u, seed offset %d, score limit %d, score %d, offset %d\n", (int)(candidate - scoringCandidatePool),
				candidate->whichSetPair, readWithFewerHits, candidate->readWithFewerHitsGenomeLocation.location, candidate->seedOffset,
				compensatedScoreLimit, fewerEndScore, fewerEndGenomeLocationOffset);
		}
#endif // DEBUG

		if (fewerEndScore != -1) {
			//
			// Find and score mates.  The index in scoringMateCandidateIndex is the lowest mate (i.e., the highest index number).
			//
			unsigned mateIndex = candidate->scoringMateCandidateIndex;

			for (;;) {

				//**** 10X surrogates
				unsigned EDCompensation = 0;

				if (candidate->clusterIdx != -1)
					EDCompensation = clusterEDCompensation;

				compensatedScoreLimit = scoreLimit + EDCompensation;
				//**** 10X surrogates

				ScoringMateCandidate *mate = &scoringMateCandidates[candidate->whichSetPair][mateIndex];
				_ASSERT(genomeLocationIsWithin(mate->readWithMoreHitsGenomeLocation, candidate->readWithFewerHitsGenomeLocation, maxSpacing));
				if (!genomeLocationIsWithin(mate->readWithMoreHitsGenomeLocation, candidate->readWithFewerHitsGenomeLocation, minSpacing) && mate->bestPossibleScore <= scoreLimit - fewerEndScore) {
					//
					// It's within the range and not necessarily too poor of a match.  Consider it.
					//

					//
					// If we haven't yet scored this mate, or we've scored it and not gotten an answer, but had a higher score limit than we'd
					// use now, score it.
					//
					if (mate->score == -2 || mate->score == -1 && mate->scoreLimit < compensatedScoreLimit - fewerEndScore) {
						scoreLocation(readWithMoreHits, setPairDirection[candidate->whichSetPair][readWithMoreHits], GenomeLocationAsInt64(mate->readWithMoreHitsGenomeLocation),
							mate->seedOffset, compensatedScoreLimit - fewerEndScore, &mate->score, &mate->matchProbability,
							&mate->genomeOffset);
#ifdef _DEBUG
						if (_DumpAlignments) {
							printf("Scored mate candidate %d, set pair %d, read %d, location %u, seed offset %d, score limit %d, score %d, offset %d\n",
								(int)(mate - scoringMateCandidates[candidate->whichSetPair]), candidate->whichSetPair, readWithMoreHits, mate->readWithMoreHitsGenomeLocation.location,
								mate->seedOffset, compensatedScoreLimit - fewerEndScore, mate->score, mate->genomeOffset);
						}
#endif // _DEBUG

						_ASSERT(-1 == mate->score || mate->score >= mate->bestPossibleScore);

						mate->scoreLimit = compensatedScoreLimit - fewerEndScore;
					}

					if (mate->score != -1) {
						double pairProbability = mate->matchProbability * fewerEndMatchProbability;
		
						//**** 10X apply probability correction
						pairProbability *= probabilityCorrectionFactor;
						//**** 10X apply probability correction

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
								(mergeCandidate->whichSetPair != candidate->whichSetPair ||
								genomeLocationIsWithin(mergeCandidate->readWithFewerHitsGenomeLocation, candidate->readWithFewerHitsGenomeLocation + fewerEndGenomeLocationOffset, 50) );
								mergeCandidate--) {

								// Because now the candidates from both directions are mixed up, we no longer terminates the loop when we see opposite direction
								if (mergeCandidate->whichSetPair == candidate->whichSetPair && mergeCandidate->mergeAnchor != NULL) {
									candidate->mergeAnchor = mergeAnchor = mergeCandidate->mergeAnchor;
									break;
								}
							}

							if (NULL == mergeAnchor) {
								for (ScoringCandidate *mergeCandidate = candidate + 1;
									mergeCandidate < scoringCandidatePool + lowestFreeScoringCandidatePoolEntry &&
									(mergeCandidate->whichSetPair != candidate->whichSetPair ||
									genomeLocationIsWithin(mergeCandidate->readWithFewerHitsGenomeLocation, candidate->readWithFewerHitsGenomeLocation + fewerEndGenomeLocationOffset, 50) );
									mergeCandidate++) {

									// Because now the candidates from both directions are mixed up, we no longer terminates the loop when we see opposite direction
									if (mergeCandidate->whichSetPair == candidate->whichSetPair && mergeCandidate->mergeAnchor != NULL) {
										candidate->mergeAnchor = mergeAnchor = mergeCandidate->mergeAnchor;
										break;
									}
								}
							}
						}

						bool merged;
						unsigned compensatedEDIdx = pairScore + clusterEDCompensation - EDCompensation;
						unsigned oldCompensatedEDIdx = 0;
						double oldPairProbability;

						if (NULL == mergeAnchor) {
							if (firstFreeMergeAnchor >= mergeAnchorPoolSize) {
								WriteErrorMessage("Ran out of merge anchor pool entries.  Perhaps rerunning with a larger value of -mcp will help\n");
								soft_exit(1);
							}

							mergeAnchor = &mergeAnchorPool[firstFreeMergeAnchor];

							firstFreeMergeAnchor++;

							mergeAnchor->init(mate->readWithMoreHitsGenomeLocation + mate->genomeOffset, candidate->readWithFewerHitsGenomeLocation + fewerEndGenomeLocationOffset,
								pairProbability, pairScore);

							merged = false;
							oldPairProbability = 0;
							candidate->mergeAnchor = mergeAnchor;
						}
						else {
							//return true if this mapping should be ignored
							merged = mergeAnchor->checkMerge(mate->readWithMoreHitsGenomeLocation + mate->genomeOffset, candidate->readWithFewerHitsGenomeLocation + fewerEndGenomeLocationOffset,
								pairProbability, compensatedEDIdx, &oldCompensatedEDIdx, &oldPairProbability);
						}

						if (!merged) {
							//
							// Back out the probability of the old match that we're merged with, if any.  The max
							// is necessary because a + b - b is not necessarily a in floating point.  If there
							// was no merge, the oldPairProbability is 0.
							//

							probabilityForED[oldCompensatedEDIdx] = __max(0, probabilityForED[oldCompensatedEDIdx] - oldPairProbability);

							// New ED lower bound. Need to do some cleaning up.
							if (pairScore + extraSearchDepth < worstED) {
								for (unsigned scoreIdx = pairScore + extraSearchDepth + 1; scoreIdx <= worstED; scoreIdx++)
									probabilityOfAllPairs = __max(0, probabilityOfAllPairs - probabilityForED[scoreIdx]);

								// update the new worstED.
								worstED = pairScore + extraSearchDepth;
							}

							bool isBestHit = false;
						
							//****10X debug switch
							//if (pairScore <= maxK && (pairScore < bestPairScore || (pairScore == bestPairScore && pairProbability > probabilityOfBestPair))) {
							if (pairScore <= maxK && pairProbability > probabilityOfBestPair) {
							//****10X debug switch
								//
								// A new best hit.
								//
								if (maxEditDistanceForSecondaryResults != -1 && (unsigned)maxEditDistanceForSecondaryResults >= bestPairScore - pairScore) {
									//
									// Move the old best to be a secondary alignment.  This won't happen on the first time we get a valid alignment,
									// because bestPairScore is initialized to be very large.
									//
									//
									if (*nSecondaryResults >= secondaryResultBufferSize) {
										*nSecondaryResults = secondaryResultBufferSize + 1;
										return true;
									}

									PairedAlignmentResult *secondaryResult = &secondaryResults[*nSecondaryResults];
									secondaryResult->alignedAsPair = true;
									secondaryResult->fromAlignTogether = true;

									for (int r = 0; r < NUM_READS_PER_PAIR; r++) {
										secondaryResult->direction[r] = bestResultDirection[r];
										secondaryResult->location[r] = bestResultGenomeLocation[r];
										secondaryResult->mapq[r] = 0;
										secondaryResult->score[r] = bestResultScore[r];
										secondaryResult->status[r] = MultipleHits;
										secondaryResult->clusterIdx = bestClusterIdx;
									}

									(*nSecondaryResults)++;

								}
								bestPairScore = pairScore;
								probabilityOfBestPair = pairProbability;
								bestResultGenomeLocation[readWithFewerHits] = candidate->readWithFewerHitsGenomeLocation + fewerEndGenomeLocationOffset;
								bestResultGenomeLocation[readWithMoreHits] = mate->readWithMoreHitsGenomeLocation + mate->genomeOffset;
								bestResultScore[readWithFewerHits] = fewerEndScore;
								bestResultScore[readWithMoreHits] = mate->score;
								bestResultDirection[readWithFewerHits] = setPairDirection[candidate->whichSetPair][readWithFewerHits];
								bestResultDirection[readWithMoreHits] = setPairDirection[candidate->whichSetPair][readWithMoreHits];
								bestClusterIdx = candidate->clusterIdx;

								if (!noUkkonen) {
									scoreLimit = bestPairScore + extraSearchDepth;
								}

								isBestHit = true;
							}
							else {
								if (maxEditDistanceForSecondaryResults != -1 && pairScore <= maxK && (unsigned)maxEditDistanceForSecondaryResults >= pairScore - bestPairScore) {
									//
									// A secondary result to save.
									//
									if (*nSecondaryResults >= secondaryResultBufferSize) {
										*nSecondaryResults = secondaryResultBufferSize + 1;
										return true;
									}

									PairedAlignmentResult *secondaryResult = &secondaryResults[*nSecondaryResults];
									secondaryResult->alignedAsPair = true;
									secondaryResult->direction[readWithMoreHits] = setPairDirection[candidate->whichSetPair][readWithMoreHits];
									secondaryResult->direction[readWithFewerHits] = setPairDirection[candidate->whichSetPair][readWithFewerHits];
									secondaryResult->fromAlignTogether = true;
									secondaryResult->location[readWithMoreHits] = mate->readWithMoreHitsGenomeLocation + mate->genomeOffset;
									secondaryResult->location[readWithFewerHits] = candidate->readWithFewerHitsGenomeLocation + fewerEndGenomeLocationOffset;
									secondaryResult->mapq[0] = secondaryResult->mapq[1] = 0;
									secondaryResult->score[readWithMoreHits] = mate->score;
									secondaryResult->score[readWithFewerHits] = fewerEndScore;
									secondaryResult->status[readWithFewerHits] = secondaryResult->status[readWithMoreHits] = MultipleHits;
									secondaryResult->clusterIdx = candidate->clusterIdx;
									secondaryResult->probability = pairProbability;

									(*nSecondaryResults)++;
								}
							}

							if (pairScore <= maxK) {
								probabilityOfAllPairs += pairProbability; // This is only a vague estimate of probability.
								probabilityForED[pairScore] += pairProbability;
							}
#ifdef  _DEBUG
							if (_DumpAlignments) {
								printf("Added %e (= %e * %e) @ (%u, %u), giving new probability of all pairs %e, score %d = %d + %d%s\n",
									pairProbability, mate->matchProbability, fewerEndMatchProbability,
									candidate->readWithFewerHitsGenomeLocation.location + fewerEndGenomeLocationOffset, mate->readWithMoreHitsGenomeLocation.location + mate->genomeOffset,
									probabilityOfAllPairs,
									pairScore, fewerEndScore, mate->score, isBestHit ? " New best hit" : "");
							}
#endif  // _DEBUG

							if (probabilityOfAllPairs >= 4.9 && -1 == maxEditDistanceForSecondaryResults) {
								//
								// Nothing will rescue us from a 0 MAPQ, so just stop looking. Perceed to next stage.
								// default, not exiting early
								return false;
							}
						}
					}// if the mate has a non -1 score
				}

				if (mateIndex == 0 || !genomeLocationIsWithin(scoringMateCandidates[candidate->whichSetPair][mateIndex - 1].readWithMoreHitsGenomeLocation, candidate->readWithFewerHitsGenomeLocation, maxSpacing)) {
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
	// default, not exiting early
	return false;
}


void TenXSingleAligner::align_phase_4(
	Read *read0,
	Read *read1,
	PairedAlignmentResult *result,
	int maxEditDistanceForSecondaryResults,
	_int64 *nSecondaryResults,
	PairedAlignmentResult *secondaryResults,             // The caller passes in a buffer of secondaryResultBufferSize and it's filled in by align()
	_int64 maxSecondaryResultsToReturn,
	unsigned *popularSeedsSkipped,
	unsigned bestPairScore,
	GenomeLocation *bestResultGenomeLocation,
	Direction *bestResultDirection,
	double probabilityOfAllPairs,
	unsigned *bestResultScore,
	double probabilityOfBestPair,
	int bestClusterIdx
)
{
	if (bestPairScore == 65536) {
		//
		// Found nothing.
		//
		for (unsigned whichRead = 0; whichRead < NUM_READS_PER_PAIR; whichRead++) {
			result->location[whichRead] = InvalidGenomeLocation;
			result->mapq[whichRead] = 0;
			result->score[whichRead] = -1;
			result->status[whichRead] = NotFound;
			result->clippingForReadAdjustment[whichRead] = 0;
			result->clusterIdx = -1;
#ifdef  _DEBUG
			if (_DumpAlignments) {
				printf("No sufficiently good pairs found.\n");
			}
#endif  // DEBUG
		}
	}
	else {
		for (unsigned whichRead = 0; whichRead < NUM_READS_PER_PAIR; whichRead++) {
			result->location[whichRead] = bestResultGenomeLocation[whichRead];
			result->direction[whichRead] = bestResultDirection[whichRead];
			result->mapq[whichRead] = computeMAPQ(probabilityOfAllPairs, probabilityOfBestPair, bestResultScore[whichRead], popularSeedsSkipped[0] + popularSeedsSkipped[1]);
			result->status[whichRead] = result->mapq[whichRead] > printStatsMapQLimit ? SingleHit : MultipleHits;
			result->score[whichRead] = bestResultScore[whichRead];
			result->clippingForReadAdjustment[whichRead] = 0;
			result->clusterIdx = bestClusterIdx;
		}
#ifdef  _DEBUG
		if (_DumpAlignments) {
			printf("Returned %u %s %u %s with MAPQ %d and %d, probability of all pairs %e, probability of best pair %e\n",
				result->location[0].location, result->direction[0] == RC ? "RC" : "", result->location[1].location, result->direction[1] == RC ? "RC" : "", result->mapq[0], result->mapq[1],
				probabilityOfAllPairs, probabilityOfBestPair);
		}
#endif  // DEBUG
	}

	//
	// Get rid of any secondary results that are too far away from the best score.  (NB: the rest of the code in align() is very similar to BaseAligner::finalizeSecondaryResults.  Sorry)
	//


	Read *inputReads[2] = { read0, read1 };
	for (int whichRead = 0; whichRead < NUM_READS_PER_PAIR; whichRead++) {
		result->scorePriorToClipping[whichRead] = result->score[whichRead];
	}

	if (!ignoreAlignmentAdjustmentsForOm) {
		//
		// Start by adjusting the alignments.
		//
		alignmentAdjuster.AdjustAlignments(inputReads, result);
		if (result->status[0] != NotFound && result->status[1] != NotFound && !ignoreAlignmentAdjustmentsForOm) {
			bestPairScore = result->score[0] + result->score[1];
		}

		for (int i = 0; i < *nSecondaryResults; i++) {
			for (int whichRead = 0; whichRead < NUM_READS_PER_PAIR; whichRead++) {
				secondaryResults[i].scorePriorToClipping[whichRead] = secondaryResults[i].score[whichRead];
			}
			alignmentAdjuster.AdjustAlignments(inputReads, &secondaryResults[i]);
			if (secondaryResults[i].status[0] != NotFound && secondaryResults[i].status[1] != NotFound && !ignoreAlignmentAdjustmentsForOm) {
				bestPairScore = __min(bestPairScore, (unsigned)(secondaryResults[i].score[0] + secondaryResults[i].score[1]));
			}
		}
	}
	else {
		for (int i = 0; i < *nSecondaryResults; i++) {
			for (int whichRead = 0; whichRead < NUM_READS_PER_PAIR; whichRead++) {
				secondaryResults[i].scorePriorToClipping[whichRead] = secondaryResults[i].score[whichRead];
			}
		}
	}

	//****10X unfinished yet! Needs to take cluster info into consideration!
	int i = 0;
	while (i < *nSecondaryResults) {
		if ((int)(secondaryResults[i].score[0] + secondaryResults[i].score[1]) > (int)bestPairScore + maxEditDistanceForSecondaryResults ||
			secondaryResults[i].status[0] == NotFound || secondaryResults[i].status[1] == NotFound) {

			secondaryResults[i] = secondaryResults[(*nSecondaryResults) - 1];
			(*nSecondaryResults)--;
		}
		else {
			i++;
		}
	}

	//
	// Now check to see if there are too many for any particular contig.
	//
	if (maxSecondaryAlignmentsPerContig > 0 && result->status[0] != NotFound) {
		//
		// Run through the results and count the number of results per contig, to see if any of them are too big.
		// First, record the primary result.
		//

		bool anyContigHasTooManyResults = false;
		contigCountEpoch++;

		int primaryContigNum = genome->getContigNumAtLocation(result->location[0]);
		hitsPerContigCounts[primaryContigNum].hits = 1;
		hitsPerContigCounts[primaryContigNum].epoch = contigCountEpoch;


		for (i = 0; i < *nSecondaryResults; i++) {
			int contigNum = genome->getContigNumAtLocation(secondaryResults[i].location[0]);    // We know they're on the same contig, so either will do
			if (hitsPerContigCounts[contigNum].epoch != contigCountEpoch) {
				hitsPerContigCounts[contigNum].epoch = contigCountEpoch;
				hitsPerContigCounts[contigNum].hits = 0;
			}

			hitsPerContigCounts[contigNum].hits++;
			if (hitsPerContigCounts[contigNum].hits > maxSecondaryAlignmentsPerContig) {
				anyContigHasTooManyResults = true;
				break;
			}
		}

		if (anyContigHasTooManyResults) {
			//
			// Just sort them all, in order of contig then hit depth.
			//
			qsort(secondaryResults, *nSecondaryResults, sizeof(*secondaryResults), PairedAlignmentResult::compareByContigAndScore);

			//
			// Now run through and eliminate any contigs with too many hits.  We can't use the same trick at the first loop above, because the
			// counting here relies on the results being sorted.  So, instead, we just copy them as we go.
			//
			int currentContigNum = -1;
			int currentContigCount = 0;
			int destResult = 0;

			for (int sourceResult = 0; sourceResult < *nSecondaryResults; sourceResult++) {
				int contigNum = genome->getContigNumAtLocation(secondaryResults[sourceResult].location[0]);
				if (contigNum != currentContigNum) {
					currentContigNum = contigNum;
					currentContigCount = (contigNum == primaryContigNum) ? 1 : 0;
				}

				currentContigCount++;

				if (currentContigCount <= maxSecondaryAlignmentsPerContig) {
					//
					// Keep it.  If we don't get here, then we don't copy the result and
					// don't increment destResult.  And yes, this will sometimes copy a
					// result over itself.  That's harmless.
					//
					secondaryResults[destResult] = secondaryResults[sourceResult];
					destResult++;
				}
			} // for each source result
			*nSecondaryResults = destResult;
		}
	} // if we're limiting by contig


	if (*nSecondaryResults > maxSecondaryResultsToReturn) {
		qsort(secondaryResults, *nSecondaryResults, sizeof(*secondaryResults), PairedAlignmentResult::compareByScore);
		*nSecondaryResults = maxSecondaryResultsToReturn;   // Just truncate it
	}


}


bool
TenXSingleAligner::align(
	Read                  *read0,
	Read                  *read1,
	PairedAlignmentResult *result,
	int                    maxEditDistanceForSecondaryResults,
	_int64                 secondaryResultBufferSize,
	_int64                *nSecondaryResults,
	PairedAlignmentResult *secondaryResults,             // The caller passes in a buffer of secondaryResultBufferSize and it's filled in by align()
	_int64                 maxSecondaryResultsToReturn
)
{
	// Initialize data before phase 1
	result->nLVCalls = 0;
	result->nSmallHits = 0;
	result->clippingForReadAdjustment[0] = result->clippingForReadAdjustment[1] = 0;
	unsigned popularSeedsSkipped[NUM_READS_PER_PAIR];

	//**** Phase 1
	if (align_phase_1(read0, read1, popularSeedsSkipped))
		return true;
	
	//**** Phase 2
	align_phase_2();

	// Declare and initialize data before phase 3
	unsigned bestPairScore = 65536;
	GenomeLocation bestResultGenomeLocation[NUM_READS_PER_PAIR];
	Direction bestResultDirection[NUM_READS_PER_PAIR];
	double probabilityOfAllPairs = 0;
	unsigned bestResultScore[NUM_READS_PER_PAIR];
	double probabilityOfBestPair = 0;
	int bestClusterIdx = -1;

	//**** Phase 3
	if (align_phase_3(maxEditDistanceForSecondaryResults, secondaryResultBufferSize, nSecondaryResults, secondaryResults, maxSecondaryResultsToReturn,
		bestPairScore, bestResultGenomeLocation, bestResultDirection, probabilityOfAllPairs, bestResultScore, popularSeedsSkipped, probabilityOfBestPair, bestClusterIdx, 1, 0))
		return false; // Not enough space for secondary alignment. Flag is raised

	//**** Phase 4
	align_phase_4(read0, read1, result, maxEditDistanceForSecondaryResults, nSecondaryResults, secondaryResults, maxSecondaryResultsToReturn, popularSeedsSkipped, bestPairScore, bestResultGenomeLocation, bestResultDirection, probabilityOfAllPairs, bestResultScore, probabilityOfBestPair, bestClusterIdx);

	return true;
}

void
TenXSingleAligner::scoreLocation(
	unsigned             whichRead,
	Direction            direction,
	GenomeLocation       genomeLocation,
	unsigned             seedOffset,
	unsigned             scoreLimit,
	unsigned            *score,
	double              *matchProbability,
	int                 *genomeLocationOffset)
{
	nLocationsScored++;

	Read *readToScore = reads[whichRead][direction];
	unsigned readDataLength = readToScore->getDataLength();
	GenomeDistance genomeDataLength = readDataLength + MAX_K; // Leave extra space in case the read has deletions
	const char *data = genome->getSubstring(genomeLocation, genomeDataLength);

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

	_ASSERT(!memcmp(data + seedOffset, readToScore->getData() + seedOffset, seedLen));    // that the seed actually matches

	int textLen;
	if (genomeDataLength - tailStart > INT32_MAX) {
		textLen = INT32_MAX;
	}
	else {
		textLen = (int)(genomeDataLength - tailStart);
	}
	score1 = landauVishkin->computeEditDistance(data + tailStart, textLen, readToScore->getData() + tailStart, readToScore->getQuality() + tailStart, readLen - tailStart,
		scoreLimit, &matchProb1);
	if (score1 == -1) {
		*score = -1;
	}
	else {
		// The tail of the read matched; now let's reverse the reference genome data and match the head
		int limitLeft = scoreLimit - score1;
		score2 = reverseLandauVishkin->computeEditDistance(data + seedOffset, seedOffset + MAX_K, reversedRead[whichRead][direction] + readLen - seedOffset,
			reads[whichRead][OppositeDirection(direction)]->getQuality() + readLen - seedOffset, seedOffset, limitLeft, &matchProb2, genomeLocationOffset);

		if (score2 == -1) {
			*score = -1;
		}
		else {
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
TenXSingleAligner::HashTableHitSet::firstInit(unsigned maxSeeds_, unsigned maxMergeDistance_, BigAllocator *allocator, bool doesGenomeIndexHave64BitLocations_)
{
	maxSeeds = maxSeeds_;
	maxMergeDistance = maxMergeDistance_;
	doesGenomeIndexHave64BitLocations = doesGenomeIndexHave64BitLocations_;
	nLookupsUsed = 0;
	if (doesGenomeIndexHave64BitLocations) {
		lookups64 = (HashTableLookup<GenomeLocation> *)allocator->allocate(sizeof(HashTableLookup<GenomeLocation>) * maxSeeds);
		lookups32 = NULL;
	}
	else {
		lookups32 = (HashTableLookup<unsigned> *)allocator->allocate(sizeof(HashTableLookup<unsigned>) * maxSeeds);
		lookups64 = NULL;
	}
	disjointHitSets = (DisjointHitSet *)allocator->allocate(sizeof(DisjointHitSet) * maxSeeds);
}
void
TenXSingleAligner::HashTableHitSet::init()
{
	nLookupsUsed = 0;
	currentDisjointHitSet = -1;
	if (doesGenomeIndexHave64BitLocations) {
		lookupListHead64->nextLookupWithRemainingMembers = lookupListHead64->prevLookupWithRemainingMembers = lookupListHead64;
		lookupListHead32->nextLookupWithRemainingMembers = lookupListHead32->prevLookupWithRemainingMembers = NULL;
	}
	else {
		lookupListHead32->nextLookupWithRemainingMembers = lookupListHead32->prevLookupWithRemainingMembers = lookupListHead32;
		lookupListHead64->nextLookupWithRemainingMembers = lookupListHead64->prevLookupWithRemainingMembers = NULL;
	}
}


//
// I apologize for this, but I had to do two versions of recordLookup, one for the 32 bit and one for the 64 bit version.  The options were
// copying the code or doing a macro with the types as parameters.  I chose macro, so you get ugly but unlikely to accidentally diverge.
// At least it's just isolated to the HashTableHitSet class.
//

#define RL(lookups, glType, lookupListHead)                                                                                                                 \
    void                                                                                                                                                    \
TenXSingleAligner::HashTableHitSet::recordLookup(unsigned seedOffset, _int64 nHits, const glType *hits, bool beginsDisjointHitSet)               \
{                                                                                                                                                           \
    _ASSERT(nLookupsUsed < maxSeeds);                                                                                                                       \
    if (beginsDisjointHitSet) {                                                                                                                             \
        currentDisjointHitSet++;                                                                                                                            \
        _ASSERT(currentDisjointHitSet < (int)maxSeeds);                                                                                                     \
        disjointHitSets[currentDisjointHitSet].countOfExhaustedHits = 0;                                                                                    \
    }                                                                                                                                                       \
                                                                                                                                                            \
    if (0 == nHits) {                                                                                                                                       \
        disjointHitSets[currentDisjointHitSet].countOfExhaustedHits++;                                                                                      \
    } else {                                                                                                                                                \
        _ASSERT(currentDisjointHitSet != -1);    /* Essentially that beginsDisjointHitSet is set for the first recordLookup call */                         \
        lookups[nLookupsUsed].currentHitForIntersection = 0;                                                                                                \
        lookups[nLookupsUsed].hits = hits;                                                                                                                  \
        lookups[nLookupsUsed].nHits = nHits;                                                                                                                \
        lookups[nLookupsUsed].seedOffset = seedOffset;                                                                                                      \
        lookups[nLookupsUsed].whichDisjointHitSet = currentDisjointHitSet;                                                                                  \
                                                                                                                                                            \
        /* Trim off any hits that are smaller than seedOffset, since they are clearly meaningless. */                                                       \
                                                                                                                                                            \
        while (lookups[nLookupsUsed].nHits > 0 && lookups[nLookupsUsed].hits[lookups[nLookupsUsed].nHits - 1] < lookups[nLookupsUsed].seedOffset) {         \
            lookups[nLookupsUsed].nHits--;                                                                                                                  \
        }                                                                                                                                                   \
                                                                                                                                                            \
        /* Add this lookup into the non-empty lookup list. */                                                                                               \
                                                                                                                                                            \
        lookups[nLookupsUsed].prevLookupWithRemainingMembers = lookupListHead;                                                                              \
        lookups[nLookupsUsed].nextLookupWithRemainingMembers = lookupListHead->nextLookupWithRemainingMembers;                                              \
        lookups[nLookupsUsed].prevLookupWithRemainingMembers->nextLookupWithRemainingMembers =                                                              \
            lookups[nLookupsUsed].nextLookupWithRemainingMembers->prevLookupWithRemainingMembers = &lookups[nLookupsUsed];                                  \
                                                                                                                                                            \
        if (doAlignerPrefetch) {                                                                                                                            \
            _mm_prefetch((const char *)&lookups[nLookupsUsed].hits[lookups[nLookupsUsed].nHits / 2], _MM_HINT_T2);                                          \
        }                                                                                                                                                   \
                                                                                                                                                            \
        nLookupsUsed++;                                                                                                                                     \
    }                                                                                                                                                       \
}

RL(lookups32, unsigned, lookupListHead32)
RL(lookups64, GenomeLocation, lookupListHead64)

#undef RL


unsigned
TenXSingleAligner::HashTableHitSet::computeBestPossibleScoreForCurrentHit()
{
	//
	// Now compute the best possible score for the hit.  This is the largest number of misses in any disjoint hit set.
	//
	for (int i = 0; i <= currentDisjointHitSet; i++) {
		disjointHitSets[i].missCount = disjointHitSets[i].countOfExhaustedHits;
	}

	//
	// Another macro.  Sorry again.
	//
#define loop(glType, lookupListHead)                                                                                                                                \
    for (HashTableLookup<glType> *lookup = lookupListHead->nextLookupWithRemainingMembers; lookup != lookupListHead;                                                \
         lookup = lookup->nextLookupWithRemainingMembers) {                                                                                                         \
                                                                                                                                                                    \
        if (!(lookup->currentHitForIntersection != lookup->nHits &&                                                                                                 \
                genomeLocationIsWithin(lookup->hits[lookup->currentHitForIntersection], mostRecentLocationReturned + lookup->seedOffset,  maxMergeDistance) ||      \
            lookup->currentHitForIntersection != 0 &&                                                                                                               \
                genomeLocationIsWithin(lookup->hits[lookup->currentHitForIntersection-1], mostRecentLocationReturned + lookup->seedOffset,  maxMergeDistance))) {   \
                                                                                                                                                                    \
            /* This one was not close enough. */                                                                                                                    \
                                                                                                                                                                    \
            disjointHitSets[lookup->whichDisjointHitSet].missCount++;                                                                                               \
        }                                                                                                                                                           \
    }

	if (doesGenomeIndexHave64BitLocations) {
		loop(GenomeLocation, lookupListHead64);
	}
	else {
		loop(unsigned, lookupListHead32);
	}
#undef loop

	unsigned bestPossibleScoreSoFar = 0;
	for (int i = 0; i <= currentDisjointHitSet; i++) {
		bestPossibleScoreSoFar = max(bestPossibleScoreSoFar, disjointHitSets[i].missCount);
	}

	return bestPossibleScoreSoFar;
}

bool
TenXSingleAligner::HashTableHitSet::getNextHitLessThanOrEqualTo(GenomeLocation maxGenomeLocationToFind, GenomeLocation *actualGenomeLocationFound, unsigned *seedOffsetFound)
{

	bool anyFound = false;
	GenomeLocation bestLocationFound = 0;
	for (unsigned i = 0; i < nLookupsUsed; i++) {
		//
		// Binary search from the current starting offset to either the right place or the end.
		//
		_int64 limit[2];
		GenomeLocation maxGenomeLocationToFindThisSeed;

		if (doesGenomeIndexHave64BitLocations) {
			limit[0] = (_int64)lookups64[i].currentHitForIntersection;
			limit[1] = (_int64)lookups64[i].nHits - 1;
			maxGenomeLocationToFindThisSeed = maxGenomeLocationToFind + lookups64[i].seedOffset;
		}
		else {
			limit[0] = (_int64)lookups32[i].currentHitForIntersection;
			limit[1] = (_int64)lookups32[i].nHits - 1;
			maxGenomeLocationToFindThisSeed = maxGenomeLocationToFind + lookups32[i].seedOffset;
		}

		while (limit[0] <= limit[1]) {
			_int64 probe = (limit[0] + limit[1]) / 2;
			if (doAlignerPrefetch) { // not clear this helps.  We're probably not far enough ahead.
				if (doesGenomeIndexHave64BitLocations) {
					_mm_prefetch((const char *)&lookups64[i].hits[(limit[0] + probe) / 2 - 1], _MM_HINT_T2);
					_mm_prefetch((const char *)&lookups64[i].hits[(limit[1] + probe) / 2 + 1], _MM_HINT_T2);
				}
				else {
					_mm_prefetch((const char *)&lookups32[i].hits[(limit[0] + probe) / 2 - 1], _MM_HINT_T2);
					_mm_prefetch((const char *)&lookups32[i].hits[(limit[1] + probe) / 2 + 1], _MM_HINT_T2);
				}
			}
			//
			// Recall that the hit sets are sorted from largest to smallest, so the strange looking logic is actually right.
			// We're evaluating the expression "lookups[i].hits[probe] <= maxGenomeOffsetToFindThisSeed && (probe == 0 || lookups[i].hits[probe-1] > maxGenomeOffsetToFindThisSeed)"
			// It's written in this strange way just so the profile tool will show us where the time's going.
			//
			GenomeLocation probeHit;
			GenomeLocation probeMinusOneHit;
			unsigned seedOffset;
			if (doesGenomeIndexHave64BitLocations) {
				probeHit = lookups64[i].hits[probe];
				probeMinusOneHit = lookups64[i].hits[probe - 1];
				seedOffset = lookups64[i].seedOffset;
			}
			else {
				probeHit = lookups32[i].hits[probe];
				probeMinusOneHit = lookups32[i].hits[probe - 1];
				seedOffset = lookups32[i].seedOffset;
			}
			unsigned clause1 = probeHit <= maxGenomeLocationToFindThisSeed;
			unsigned clause2 = probe == 0;

			if (clause1 && (clause2 || probeMinusOneHit > maxGenomeLocationToFindThisSeed)) {
				if (probeHit - seedOffset > bestLocationFound) {
					anyFound = true;
					mostRecentLocationReturned = *actualGenomeLocationFound = bestLocationFound = probeHit - seedOffset;
					*seedOffsetFound = seedOffset;
				}

				if (doesGenomeIndexHave64BitLocations) {
					lookups64[i].currentHitForIntersection = probe;
				}
				else {
					lookups32[i].currentHitForIntersection = probe;
				}
				break;
			}

			if (probeHit > maxGenomeLocationToFindThisSeed) {   // Recode this without the if to avoid the hard-to-predict branch.
				limit[0] = probe + 1;
			}
			else {
				limit[1] = probe - 1;
			}
		} // While we're looking

		if (limit[0] > limit[1]) {
			// We're done with this lookup.
			if (doesGenomeIndexHave64BitLocations) {
				lookups64[i].currentHitForIntersection = lookups64[i].nHits;
			}
			else {
				lookups32[i].currentHitForIntersection = lookups32[i].nHits;
			}
		}
	} // For each lookup

	_ASSERT(!anyFound || *actualGenomeLocationFound <= maxGenomeLocationToFind);

	return anyFound;
}


bool
TenXSingleAligner::HashTableHitSet::getFirstHit(GenomeLocation *genomeLocation, unsigned *seedOffsetFound)
{
	bool anyFound = false;
	*genomeLocation = 0;

	//
	// Yet another macro.  This makes me want to write in a better language sometimes.  But then it would be too slow.  :-(
	//

#define LOOP(lookups)                                                                                                                       \
    for (unsigned i = 0; i < nLookupsUsed; i++) {                                                                                           \
        if (lookups[i].nHits > 0 && lookups[i].hits[0] - lookups[i].seedOffset > GenomeLocationAsInt64(*genomeLocation)) {                  \
            mostRecentLocationReturned = *genomeLocation = lookups[i].hits[0] - lookups[i].seedOffset;                                      \
            *seedOffsetFound = lookups[i].seedOffset;                                                                                       \
            anyFound = true;                                                                                                                \
        }                                                                                                                                   \
    }

	if (doesGenomeIndexHave64BitLocations) {
		LOOP(lookups64);
	}
	else {
		LOOP(lookups32);
	}

#undef LOOP

	return !anyFound;
}

bool
TenXSingleAligner::HashTableHitSet::getNextLowerHit(GenomeLocation *genomeLocation, unsigned *seedOffsetFound)
{
	//
	// Look through all of the lookups and find the one with the highest location smaller than the current one.
	//
	GenomeLocation foundLocation = 0;
	bool anyFound = false;

	//
	// Run through the lookups pushing up any that are at the most recently returned
	//

	for (unsigned i = 0; i < nLookupsUsed; i++) {
		_int64 *currentHitForIntersection;
		_int64 nHits;
		GenomeLocation hitLocation;
		unsigned seedOffset;

		//
		// A macro to initialize stuff that we need to avoid a bigger macro later.
		//
#define initVars(lookups)                                                                                               \
        currentHitForIntersection = &lookups[i].currentHitForIntersection;                                              \
        nHits = lookups[i].nHits;                                                                                       \
        seedOffset = lookups[i].seedOffset;                                                                             \
        if (nHits != *currentHitForIntersection) {                                                                      \
            hitLocation = lookups[i].hits[*currentHitForIntersection];                                                  \
        }


		if (doesGenomeIndexHave64BitLocations) {
			initVars(lookups64);
		}
		else {
			initVars(lookups32);
		}
#undef  initVars

		_ASSERT(*currentHitForIntersection == nHits || hitLocation - seedOffset <= mostRecentLocationReturned || hitLocation < seedOffset);

		if (*currentHitForIntersection != nHits && hitLocation - seedOffset == mostRecentLocationReturned) {
			(*currentHitForIntersection)++;
			if (*currentHitForIntersection == nHits) {
				continue;
			}
			if (doesGenomeIndexHave64BitLocations) {
				hitLocation = lookups64[i].hits[*currentHitForIntersection];
			}
			else {
				hitLocation = lookups32[i].hits[*currentHitForIntersection];
			}
		}

		if (*currentHitForIntersection != nHits) {
			if (foundLocation < hitLocation - seedOffset && // found location is OK
				hitLocation >= seedOffset) // found location isn't too small to push us before the beginning of the genome
			{
				*genomeLocation = foundLocation = hitLocation - seedOffset;
				*seedOffsetFound = seedOffset;
				anyFound = true;
			}
		}
	}

	if (anyFound) {
		mostRecentLocationReturned = foundLocation;
	}

	return anyFound;
}

bool
TenXSingleAligner::MergeAnchor::checkMerge(GenomeLocation newMoreHitLocation, GenomeLocation newFewerHitLocation, double newMatchProbability, int newPairScore,
	unsigned *oldMatchED, double *oldMatchProbability)
{
	if (locationForReadWithMoreHits == InvalidGenomeLocation || !doesRangeMatch(newMoreHitLocation, newFewerHitLocation)) {
		//
		// No merge.  Remember the new one.
		//
		locationForReadWithMoreHits = newMoreHitLocation;
		locationForReadWithFewerHits = newFewerHitLocation;
		matchProbability = newMatchProbability;
		pairScore = newPairScore;
		*oldMatchED = 0;
		*oldMatchProbability = 0.0;
		return false;
	}
	else {
		//
		// Within merge distance.  Keep the better score (or if they're tied the better match probability).
		//
		//****10X debug switch
		//if (newPairScore < pairScore || newPairScore == pairScore && newMatchProbability > matchProbability) {
		if (newPairScore < pairScore || newMatchProbability > matchProbability) {
		//****10X debug switch
#ifdef _DEBUG
			if (_DumpAlignments) {
				printf("Merge replacement at anchor (%u, %u), loc (%u, %u), old match prob %e, new match prob %e, old pair score %d, new pair score %d\n",
					locationForReadWithMoreHits.location, locationForReadWithFewerHits.location, newMoreHitLocation.location, newFewerHitLocation.location,
					matchProbability, newMatchProbability, pairScore, newPairScore);
			}
#endif // DEBUG
			*oldMatchED = pairScore;
			*oldMatchProbability = matchProbability;
			matchProbability = newMatchProbability;
			pairScore = newPairScore;
			return false;
		}
		else {
			//
			// The new one should just be ignored.
			//
#ifdef _DEBUG
			if (_DumpAlignments) {
				printf("Merged at anchor (%u, %u), loc (%u, %u), old match prob %e, new match prob %e, old pair score %d, new pair score %d\n",
					locationForReadWithMoreHits.location, locationForReadWithFewerHits.location, newMoreHitLocation.location, newFewerHitLocation.location,
					matchProbability, newMatchProbability, pairScore, newPairScore);
			}
#endif // DEBUG
			return true;
		}
	}

	_ASSERT(!"NOTREACHED");
}

const unsigned TenXSingleAligner::maxMergeDistance = 31;
