/*++

Module Name:

	TenXClusterAligner.cpp

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
#include "TenXClusterAligner.h"
#include "mapq.h"
#include "directions.h"
#include "BigAlloc.h"
#include <cstdlib>
#include "Util.h"

using namespace std;

#ifdef TRACE_PAIRED_ALIGNER
#define TRACE printf
#else
#define TRACE(...) {}
#endif

TenXClusterAligner::TenXClusterAligner(
	GenomeIndex			*index_,
	unsigned			maxReadSize,
	unsigned			maxHits,
	unsigned			maxK,
	unsigned			maxSeedsFromCommandLine,
	double				seedCoverage,
	unsigned			minWeightToCheck,
	bool				forceSpacing_,
	unsigned			extraSearchDepth,
	bool				noUkkonen,
	bool				noOrderedEvaluation,
	bool				noTruncation,
	bool				ignoreAlignmentAdjustmentsForOm,
	TenXProgressTracker	*progressTracker_,
	unsigned			maxBarcodeSize_,
	unsigned			minPairsPerCluster_,
	_uint64				minClusterSpan_,
	double				unclusteredPenalty_,
	unsigned			clusterEDCompensation_,
	unsigned			minReadLength_,
	int					maxSecondaryAlignmentsPerContig,
	unsigned			printStatsMapQLimit,
	BigAllocator		*allocator)
	: progressTracker(progressTracker_), unclusteredPenalty(unclusteredPenalty_), clusterEDCompensation(clusterEDCompensation_), maxBarcodeSize(maxBarcodeSize_), minPairsPerCluster(minPairsPerCluster_), minClusterSpan(minClusterSpan_), forceSpacing(forceSpacing_), index(index_), minReadLength(minReadLength_)
{
	// Create single-end aligners.
	singleAligner = new (allocator) BaseAligner(index, maxHits, maxK, maxReadSize,
		maxSeedsFromCommandLine, seedCoverage, minWeightToCheck, extraSearchDepth, noUkkonen, noOrderedEvaluation, noTruncation, ignoreAlignmentAdjustmentsForOm, maxSecondaryAlignmentsPerContig, printStatsMapQLimit, &lv, &reverseLV, NULL, allocator);
	for (unsigned i = 0; i < maxBarcodeSize; i++)
		progressTracker[i].aligner->setLandauVishkin(&lv, &reverseLV);

	singleSecondary[0] = singleSecondary[1] = NULL;
}

size_t
TenXClusterAligner::getBigAllocatorReservation(
	GenomeIndex		*index,
	unsigned		maxReadSize,
	unsigned		maxHits,
	unsigned		seedLen,
	unsigned		maxSeedsFromCommandLine,
	double			seedCoverage,
	unsigned		maxEditDistanceToConsider,
	unsigned		maxExtraSearchDepth,
	unsigned		maxCandidatePoolSize,
	int				maxSecondaryAlignmentsPerContig)
{
	return BaseAligner::getBigAllocatorReservation(index, false, maxHits, maxReadSize, seedLen, maxSeedsFromCommandLine, seedCoverage, maxSecondaryAlignmentsPerContig, maxExtraSearchDepth) + sizeof(TenXClusterAligner) + sizeof(_uint64);
}


TenXClusterAligner::~TenXClusterAligner()
{
	singleAligner->~BaseAligner();
}

#ifdef _DEBUG
extern bool _DumpAlignments;
#endif // _DEBUG


void TenXClusterAligner::sortAndLink()
{
	// TenX code. Initialize and sort all trackers
	qsort(progressTracker, maxBarcodeSize, sizeof(TenXProgressTracker), TenXProgressTracker::compare);

	// Link all tracker entries based on their GenomeLocation.
	unsigned nActiveTrackers = 0;
	for (unsigned pairIdx = 1; pairIdx < maxBarcodeSize; pairIdx++) {

		if (progressTracker[pairIdx - 1].pairNotDone && progressTracker[pairIdx - 1].nextLoci >= 0)
			nActiveTrackers++;

		progressTracker[pairIdx - 1].nextTracker = &progressTracker[pairIdx];
	}
	
	if (progressTracker[maxBarcodeSize - 1].pairNotDone && progressTracker[maxBarcodeSize - 1].nextLoci >= 0)
		nActiveTrackers++;

	trackerRoot = &progressTracker[0];
	updateHolder = NULL;
	//fprintf(stderr, "activeTrackers: %d\n", nActiveTrackers);
}

// Moves the cursor (cursor is modified!) to the first tracker that has a loci that's greater than the target.
unsigned trackersToMeetTargetLoci(
	TenXProgressTracker *&cursor,
	GenomeLocation clusterBoundary)
{
	unsigned cursorCounter = 0;

	while (cursor->nextLoci >= 0 && cursor->nextLoci >= clusterBoundary) {
		cursorCounter++;
		cursor = cursor->nextTracker;
	}
	return cursorCounter;
}


//reference 
GenomeLocation resolveLociPtr(GenomeLocation *lociPtr) {
	if (lociPtr == NULL)
		return -1;
	else
		return *lociPtr;
}


TenXProgressTracker* traverseProgressPtr(TenXProgressTracker *cursor, unsigned steps) 
{
	while (steps > 0 && cursor->pairNotDone && cursor->nextLoci >= 0) {
		cursor = cursor->nextTracker;
		steps--;
	}
	return cursor;
}


void TenXClusterAligner::pushToUpdate(TenXProgressTracker *newUpdate)
{
	if (updateHolder == NULL || updateHolder->nextLoci <= newUpdate->nextLoci) {
		newUpdate->nextTracker = updateHolder;
		updateHolder = newUpdate;
		return;
	}

	// At least updateHolder has a bigger loci
	TenXProgressTracker* cursor = updateHolder;

	// Invariant: the loop stops when cursor->nextLoci is still > newUpdate->nextLoci but cursor->nextTracker is not!!
	while (cursor->nextTracker != NULL && cursor->nextTracker->nextLoci > newUpdate->nextLoci)
		cursor = cursor->nextTracker;

	// Insert the newUpdate
	newUpdate->nextTracker = cursor->nextTracker;
	cursor->nextTracker = newUpdate;
}


void TenXClusterAligner::mergeUpdate()
{
	if (updateHolder == NULL)
		return;

	if (trackerRoot == NULL) {
		trackerRoot = updateHolder;
		updateHolder = NULL;
		return;
	}
		
	TenXProgressTracker *parentRootCursor = trackerRoot;
	TenXProgressTracker *rootCursor = trackerRoot->nextTracker;
	TenXProgressTracker *parentUpdateCursor = NULL;
	TenXProgressTracker *updateCursor = updateHolder;

	// First update the new root, if necessary
	if (updateHolder->nextLoci >= trackerRoot->nextLoci)
		trackerRoot = updateHolder;

	// Move updateCursor to the first merge agent
	while (updateCursor != NULL && updateCursor->nextLoci >= parentRootCursor->nextLoci) {
		parentUpdateCursor = updateCursor;
		updateCursor = updateCursor->nextTracker;
	}

	if (parentUpdateCursor != NULL)
		parentUpdateCursor->nextTracker = parentRootCursor;

	// ****10X debug
	int updateListSize = 1;
	TenXProgressTracker *updateCursorTemp = updateHolder;
	while (updateCursorTemp != NULL) {
		updateListSize++;
		updateCursorTemp = updateCursorTemp->nextTracker;
	}
	int rootListSize = 1;
	TenXProgressTracker *rootCursorTemp = rootCursor;
	while (rootCursorTemp != NULL) {
		rootListSize++;
		rootCursorTemp = rootCursorTemp->nextTracker;
	}
	fprintf(stderr, "updateList size: %d\n", updateListSize);
	fflush(stderr);
	// ****10X debug
	
	while (rootCursor != NULL && updateCursor != NULL) {
		if (rootCursor->nextLoci < updateCursor->nextLoci) {
			parentRootCursor->nextTracker = updateCursor;
			// Remember that updateList is sorted!
			updateCursor = updateCursor->nextTracker;
			parentRootCursor->nextTracker->nextTracker = rootCursor;
			parentRootCursor = parentRootCursor->nextTracker;
			updateListSize--;
		}
		else {
			parentRootCursor = rootCursor;
			rootCursor = rootCursor->nextTracker;
			rootListSize--;
		}
	}

	if (rootCursor == NULL)
		parentRootCursor->nextTracker = updateCursor;

	updateHolder = NULL;
}



// Progress each single aligner to move pass targetLoc, while registering the candidate with clusterIdx. The process stops BEFORE processing end. [start, end)
// It also terminates when cursor->nextLoci >= targetLoci, whichever comes first.
void TenXClusterAligner::registerClusterForReads(struct TenXProgressTracker* preStart, struct TenXProgressTracker* start, struct TenXProgressTracker* end, class GenomeLocation clusterBoundary, int clusterIdx)
{
	if (preStart == NULL)
		_ASSERT(start == trackerRoot);

	TenXProgressTracker *cursor = start;
	while (cursor != end && cursor->nextLoci >= 0 && cursor->nextLoci > clusterBoundary) {
		cursor->aligner->align_phase_2_to_target_loc(clusterBoundary, clusterIdx);
		cursor->nextLoci = resolveLociPtr(cursor->aligner->align_phase_2_get_loci() );

		// A temporary holder, because cursor will be modified in pushToUpdate
		TenXProgressTracker *nextCursor = cursor->nextTracker;
		// Remove from root tracker list
		if (preStart != NULL)
			preStart->nextTracker = nextCursor;
		else
			trackerRoot = nextCursor;


		// Push to update list
		pushToUpdate(cursor);

		// Move to next tracker
		cursor = nextCursor;
	}
}



void TenXClusterAligner::fastForward(GenomeLocation forwardLoc) {
	for (unsigned pairIdx = 0; pairIdx < barcodeSize; pairIdx++) {
		if (progressTracker[pairIdx].pairNotDone) {
			//fprintf(stderr, "lastLoci: %lld\n", progressTracker[pairIdx].nextLoci.location);
			progressTracker[pairIdx].aligner->align_phase_2_to_target_loc(forwardLoc, -1);
			progressTracker[pairIdx].nextLoci = resolveLociPtr(progressTracker[pairIdx].aligner->align_phase_2_get_loci() );
		}
		//fprintf(stderr, "forwardLoci: %lld  afterForward: %lld\n", forwardLoc.location, resolveLociPtr(progressTracker[pairIdx].aligner->align_phase_2_get_loci() ).location);
	}
}

bool TenXClusterAligner::align_first_stage(
	unsigned		barcodeSize_)
{
	barcodeSize = barcodeSize_;
	bool barcodeFinished = true;
	for (unsigned pairIdx = 0; pairIdx < barcodeSize; pairIdx++) {
		if (progressTracker[pairIdx].pairNotDone) {
			progressTracker[pairIdx].results[0].status[0] = progressTracker[pairIdx].results[0].status[1] = NotFound;

			Read *read0 = progressTracker[pairIdx].pairedReads[0];
			Read *read1 = progressTracker[pairIdx].pairedReads[1];

			if (read0->getDataLength() < minReadLength && read1->getDataLength() < minReadLength) {
				TRACE("Reads are both too short -- returning");
				// Update primary result
				for (int whichRead = 0; whichRead < NUM_READS_PER_PAIR; whichRead++) {
					progressTracker[pairIdx].results[0].location[whichRead] = 0;
					progressTracker[pairIdx].results[0].mapq[whichRead] = 0;
					progressTracker[pairIdx].results[0].score[whichRead] = 0;
					progressTracker[pairIdx].results[0].status[whichRead] = NotFound;
				}
				progressTracker[pairIdx].results[0].alignedAsPair = false;
				progressTracker[pairIdx].results[0].fromAlignTogether = false;
				progressTracker[pairIdx].results[0].nanosInAlignTogether = 0;
				progressTracker[pairIdx].results[0].nLVCalls = 0;
				progressTracker[pairIdx].results[0].nSmallHits = 0;

				progressTracker[pairIdx].pairNotDone = false;
				progressTracker[pairIdx].singleNotDone = false;
				continue;//return true;
			}

			//At least one read of the pair is worthy of further examination
			barcodeFinished = false;

			if (read0->getDataLength() >= minReadLength && read1->getDataLength() >= minReadLength) {
				//
				// Let the LVs use the cache that we built up.
				//
				progressTracker[pairIdx].pairNotDone = !progressTracker[pairIdx].aligner->align_phase_1(read0, read1, progressTracker[pairIdx].popularSeedsSkipped);
				
				// Initialize for phase_2 if the alginer if not stopped
				if (progressTracker[pairIdx].pairNotDone) {
					progressTracker[pairIdx].pairNotDone = progressTracker[pairIdx].aligner->align_phase_2_init();
					progressTracker[pairIdx].nextLoci = resolveLociPtr(progressTracker[pairIdx].aligner->align_phase_2_get_loci() );
					progressTracker[pairIdx].nextTracker = NULL;
				}
			}
		}
	}

	// Sort all trackers and link them based on location order.
	sortAndLink();
	// debuging
	//fastForward(150000000);
	//sortAndLink();

	//fastForward(0);

	bool registeringCluster = false;
	int globalClusterId = 0;
	int clusterId;
	GenomeLocation clusterBoundary;

	// Intitialize boundary
	TenXProgressTracker *cursor; // pointer that keeps track of the progress of walking down progress trackers

	while (trackerRoot->pairNotDone && trackerRoot->nextLoci != -1) {
		// Initialization.
		cursor = trackerRoot;
		//expirationCursor = trackerRoot;

		unsigned nPotentialPairs = trackersToMeetTargetLoci(cursor, trackerRoot->nextLoci - minClusterSpan);

		if (nPotentialPairs > minPairsPerCluster || cursor == NULL || cursor->nextLoci == -1) { //this is a clustered pair, obsorb it! Note that if nextLoci is -1 we will obsorb all pairs into the last cluster.
			registeringCluster = true;
			clusterId = globalClusterId;

			if (cursor != NULL && cursor->nextLoci != -1)
				clusterBoundary = cursor->nextLoci + minClusterSpan;
			else
				clusterBoundary = -1;

			registerClusterForReads(NULL, trackerRoot, cursor, clusterBoundary, clusterId); //register the pairs and update the loci pointers.
		}
		else { //the cluster ends here.
			if (registeringCluster) { //when we were half way of adding a cluster, we need to finish it with the old targetLoc.
				fprintf(stderr, "clusterBoundary: %lld	globalClusterId: %d\n", clusterBoundary.location, globalClusterId);
				fflush(stderr);
				registeringCluster = false;
				registerClusterForReads(NULL, trackerRoot, cursor, clusterBoundary, clusterId); //use the previous id.
				globalClusterId++;
			}
			else { //we were not adding a cluster, just tag these locus as not clustered (-1) and move the locus pointer over the new targetLoc.
				cursor = traverseProgressPtr(cursor, minPairsPerCluster - nPotentialPairs);
				clusterBoundary = cursor->nextLoci + minClusterSpan;
				registerClusterForReads(NULL, trackerRoot, cursor, clusterBoundary, -1); //use the previous id.
			}
		}
		//fprintf(stderr, "clusterBoundary: %lld\n", clusterBoundary.location);
		//sortAndLink(); //fix the order. ****If we use this too often, it would be a potential performance problem. Might need to fix it later. (That's why I kept the linked list pointer!).
		mergeUpdate();
		//trackerRoot = &progressTracker[0];
	}

	return barcodeFinished;

}


bool TenXClusterAligner::align_second_stage(
	int						maxEditDistanceForSecondaryResults,
	_int64					maxSecondaryAlignmentsToReturn
)
{
	bool barcodeFinished = true;
	for (unsigned pairIdx = 0; pairIdx < barcodeSize; pairIdx++) {
		if (progressTracker[pairIdx].pairNotDone) {// && progressTracker[pairIdx].notDone) {
			Read *read0 = progressTracker[pairIdx].pairedReads[0];
			Read *read1 = progressTracker[pairIdx].pairedReads[1];

			progressTracker[pairIdx].nSingleEndSecondaryResults[0] = 0;
			progressTracker[pairIdx].nSingleEndSecondaryResults[1] = 0;

			unsigned bestPairScore = 65536;
			GenomeLocation bestResultGenomeLocation[NUM_READS_PER_PAIR];
			Direction bestResultDirection[NUM_READS_PER_PAIR];
			double probabilityOfAllPairs = 0;
			unsigned bestResultScore[NUM_READS_PER_PAIR];
			double probabilityOfBestPair = 0;

			bool secondaryBufferOverflow = progressTracker[pairIdx].aligner->align_phase_3(maxEditDistanceForSecondaryResults, progressTracker[pairIdx].secondaryResultBufferSize, &progressTracker[pairIdx].nSecondaryResults, &progressTracker[pairIdx].results[1], maxSecondaryAlignmentsToReturn, bestPairScore, bestResultGenomeLocation, bestResultDirection, probabilityOfAllPairs, bestResultScore, progressTracker[pairIdx].popularSeedsSkipped, probabilityOfBestPair, unclusteredPenalty, clusterEDCompensation);

			if (secondaryBufferOverflow) {
				progressTracker[pairIdx].nSingleEndSecondaryResults[0] = progressTracker[pairIdx].nSingleEndSecondaryResults[1] = 0;
				progressTracker[pairIdx].nSecondaryResults = progressTracker[pairIdx].secondaryResultBufferSize + 1; // So the caller knows it's the paired secondary buffer that overflowed
				barcodeFinished = false;
				continue;//return false;
			}
			else {
				progressTracker[pairIdx].aligner->align_phase_4(read0, read1, &progressTracker[pairIdx].results[0], maxEditDistanceForSecondaryResults, &progressTracker[pairIdx].nSecondaryResults, &progressTracker[pairIdx].results[1], maxSecondaryAlignmentsToReturn, progressTracker[pairIdx].popularSeedsSkipped, bestPairScore, bestResultGenomeLocation, bestResultDirection, probabilityOfAllPairs, bestResultScore, probabilityOfBestPair);
			}

			/* timing no longer makes sence
			_int64 start = timeInNanos();
			_int64 end = timeInNanos();
			result[pairIdx]->nanosInAlignTogether = end - start;
			*/

			progressTracker[pairIdx].results[0].nanosInAlignTogether = 0; //timing no longer makes sence
			progressTracker[pairIdx].results[0].fromAlignTogether = true;
			progressTracker[pairIdx].results[0].alignedAsPair = true;

			if (forceSpacing) {
				if (progressTracker[pairIdx].results[0].status[0] == NotFound) {
					progressTracker[pairIdx].results[0].fromAlignTogether = false;
				}
				else {
					_ASSERT(progressTracker[pairIdx].results[0].status[1] != NotFound); // If one's not found, so is the other
				}
				progressTracker[pairIdx].pairNotDone = false;
				progressTracker[pairIdx].singleNotDone = false;
				continue;//return true;
			}

			if (progressTracker[pairIdx].results[0].status[0] != NotFound && progressTracker[pairIdx].results[0].status[1] != NotFound) {
				//
				// Not a chimeric read.
				//
				progressTracker[pairIdx].pairNotDone = false;
				progressTracker[pairIdx].singleNotDone = false;
				continue;//return true;
			}
			//paired analysis is done anyways
			progressTracker[pairIdx].pairNotDone = false;
		}
	}
	return barcodeFinished;
}


bool TenXClusterAligner::align_third_stage(
	int maxEditDistanceForSecondaryResults,
	_int64 maxSecondaryAlignmentsToReturn
)
{
	bool barcodeFinished = true;
	for (unsigned pairIdx = 0; pairIdx < barcodeSize; pairIdx++) {
		if (progressTracker[pairIdx].singleNotDone) {// && progressTracker[pairIdx].notDone) {
			Read *read0 = progressTracker[pairIdx].pairedReads[0];
			Read *read1 = progressTracker[pairIdx].pairedReads[1];

			Read *read[NUM_READS_PER_PAIR] = { read0, read1 };
			_int64 *resultCount[2] = { &progressTracker[pairIdx].nSingleEndSecondaryResults[0], &progressTracker[pairIdx].nSingleEndSecondaryResults[1] };

			bool noOverflow = true;
			for (int r = 0; r < NUM_READS_PER_PAIR; r++) {
				SingleAlignmentResult singleResult;
				_int64 singleEndSecondaryResultsThisTime = 0;

				if (read[r]->getDataLength() < minReadLength) {
					progressTracker[pairIdx].results[0].status[r] = NotFound;
					progressTracker[pairIdx].results[0].mapq[r] = 0;
					progressTracker[pairIdx].results[0].direction[r] = FORWARD;
					progressTracker[pairIdx].results[0].location[r] = 0;
					progressTracker[pairIdx].results[0].score[r] = 0;
				}
				else {
					// We're using *nSingleEndSecondaryResultsForFirstRead because it's either 0 or what all we've seen (i.e., we know NUM_READS_PER_PAIR is 2)
					bool fitInSecondaryBuffer = //true;
						singleAligner->AlignRead(read[r], &singleResult, maxEditDistanceForSecondaryResults,
							progressTracker[pairIdx].singleSecondaryBufferSize - progressTracker[pairIdx].nSingleEndSecondaryResults[0], &singleEndSecondaryResultsThisTime,
							maxSecondaryAlignmentsToReturn, &progressTracker[pairIdx].singleEndSecondaryResults[progressTracker[pairIdx].nSingleEndSecondaryResults[0]]);

					if (!fitInSecondaryBuffer) {
						progressTracker[pairIdx].nSecondaryResults = 0;
						progressTracker[pairIdx].nSingleEndSecondaryResults[0] = progressTracker[pairIdx].singleSecondaryBufferSize + 1;
						progressTracker[pairIdx].nSingleEndSecondaryResults[1] = 0;
						barcodeFinished = false;
						noOverflow = false;
						break;//return false;
					}

					*(resultCount[r]) = singleEndSecondaryResultsThisTime;

					progressTracker[pairIdx].results[0].status[r] = singleResult.status;
					progressTracker[pairIdx].results[0].mapq[r] = singleResult.mapq / 3;   // Heavy quality penalty for chimeric reads
					progressTracker[pairIdx].results[0].direction[r] = singleResult.direction;
					progressTracker[pairIdx].results[0].location[r] = singleResult.location;
					progressTracker[pairIdx].results[0].score[r] = singleResult.score;
					progressTracker[pairIdx].results[0].scorePriorToClipping[r] = singleResult.scorePriorToClipping;
				}
			}

			// This pair is done processing, only if both directions has no overflow.
			if (noOverflow) {
				progressTracker[pairIdx].singleNotDone = false;
				progressTracker[pairIdx].results[0].fromAlignTogether = false;
				progressTracker[pairIdx].results[0].alignedAsPair = false;
			}

#ifdef _DEBUG
			if (_DumpAlignments) {
				printf("TenXClusterAligner: (%u, %u) score (%d, %d), MAPQ (%d, %d)\n\n\n", progressTracker[pairIdx].results[0].location[0].location, progressTracker[pairIdx].results[0].location[1].location,
					progressTracker[pairIdx].results[0].score[0], progressTracker[pairIdx].results[0].score[1], progressTracker[pairIdx].results[0].mapq[0], progressTracker[pairIdx].results[0].mapq[1]);
			}
#endif // _DEBUG

		}
	}
	return barcodeFinished;
}



bool TenXClusterAligner::align(
	Read					**pairedReads,
	unsigned				barcodeSize,
	PairedAlignmentResult	**result,
	int						maxEditDistanceForSecondaryResults,
	_int64					*secondaryResultBufferSize,
	_int64					*nSecondaryResults,
	_int64					*singleSecondaryBufferSize,
	_int64					maxSecondaryAlignmentsToReturn,
	_int64					*nSingleEndSecondaryResults,
	SingleAlignmentResult	**singleEndSecondaryResults,    // Single-end secondary alignments for when the paired-end alignment didn't work properly
	unsigned				*popularSeedsSkipped
)
{
	if (align_first_stage(barcodeSize))
		return true;
	if (!align_second_stage(maxEditDistanceForSecondaryResults, maxSecondaryAlignmentsToReturn))
		return false;
	if (align_third_stage(maxEditDistanceForSecondaryResults, maxSecondaryAlignmentsToReturn))
		return true;
	return false;
}
