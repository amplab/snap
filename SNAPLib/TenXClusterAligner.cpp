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
    GenomeIndex         *index_,
    unsigned            maxReadSize,
    unsigned            maxHits,
    unsigned            maxK,
    unsigned            maxSeedsFromCommandLine,
    double              seedCoverage,
    unsigned            minWeightToCheck,
    bool                forceSpacing_,
    unsigned            extraSearchDepth,
    bool                noUkkonen,
    bool                noOrderedEvaluation,
    bool                noTruncation,
    bool                ignoreAlignmentAdjustmentsForOm,
    TenXProgressTracker *progressTracker_,
    bool                *clusterIsValid_,
    unsigned            *mappedPairsPerCluster_,
    unsigned            maxBarcodeSize_,
    unsigned            minPairsPerCluster_,
    _uint64             minClusterSpan_,
    double              unclusteredPenalty_,
    unsigned            clusterEDCompensation_,
    unsigned            minReadLength_,
    int                 maxSecondaryAlignmentsPerContig,
    unsigned            printStatsMapQLimit,
    int                 maxEditDistanceForSecondaryResults_,
    _int64              maxSecondaryAlignmentsToReturn_,
    BigAllocator        *allocator)
    : progressTracker(progressTracker_), clusterIsValid(clusterIsValid_), mappedPairsPerCluster(mappedPairsPerCluster_), unclusteredPenalty(unclusteredPenalty_), clusterEDCompensation(clusterEDCompensation_), maxBarcodeSize(maxBarcodeSize_), minPairsPerCluster(minPairsPerCluster_), minClusterSpan(minClusterSpan_), forceSpacing(forceSpacing_), index(index_), minReadLength(minReadLength_), maxEditDistanceForSecondaryResults(maxEditDistanceForSecondaryResults_), maxSecondaryAlignmentsToReturn(maxSecondaryAlignmentsToReturn)
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
    GenomeIndex        *index,
    unsigned        maxReadSize,
    unsigned        maxHits,
    unsigned        seedLen,
    unsigned        maxSeedsFromCommandLine,
    double            seedCoverage,
    unsigned        maxEditDistanceToConsider,
    unsigned        maxExtraSearchDepth,
    unsigned        maxCandidatePoolSize,
    int                maxSecondaryAlignmentsPerContig)
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

        if (progressTracker[pairIdx - 1].pairNotDone && progressTracker[pairIdx - 1].nextLocus >= 0)
            nActiveTrackers++;

        progressTracker[pairIdx - 1].nextTracker = &progressTracker[pairIdx];
    }
    
    if (progressTracker[maxBarcodeSize - 1].pairNotDone && progressTracker[maxBarcodeSize - 1].nextLocus >= 0)
        nActiveTrackers++;

    trackerRoot = &progressTracker[0];
    updateHolder = NULL;
    //fprintf(stderr, "activeTrackers: %d\n", nActiveTrackers);
}

// Moves the cursor (cursor is modified!) to the first tracker that has a locus that's greater than the target.
unsigned trackersToMeetTargetLocus(
    TenXProgressTracker *&cursor,
    GenomeLocation clusterBoundary)
{
    unsigned cursorCounter = 0;

    while (cursor->nextLocus >= 0 && cursor->nextLocus >= clusterBoundary) {
        cursorCounter++;
        cursor = cursor->nextTracker;
    }
    return cursorCounter;
}


//reference 
GenomeLocation resolveLocusPtr(GenomeLocation *locusPtr) {
    if (locusPtr == NULL)
        return -1;
    else
        return *locusPtr;
}


TenXProgressTracker* traverseProgressPtr(TenXProgressTracker *cursor, unsigned steps) 
{
    while (steps > 0 && cursor->pairNotDone && cursor->nextLocus >= 0) {
        cursor = cursor->nextTracker;
        steps--;
    }
    return cursor;
}


void TenXClusterAligner::pushToUpdate(TenXProgressTracker *newUpdate)
{
    if (updateHolder == NULL || updateHolder->nextLocus <= newUpdate->nextLocus) {
        newUpdate->nextTracker = updateHolder;
        updateHolder = newUpdate;
        return;
    }

    // At least updateHolder has a bigger locus
    TenXProgressTracker* cursor = updateHolder;

    // Invariant: the loop stops when cursor->nextLocus is still > newUpdate->nextLocus but cursor->nextTracker is not!!
    while (cursor->nextTracker != NULL && cursor->nextTracker->nextLocus > newUpdate->nextLocus)
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

    /*
    */
    if (trackerRoot->nextLocus == -1) {
        TenXProgressTracker *updateCursor = updateHolder;

        while (updateCursor->nextTracker != NULL)
            updateCursor = updateCursor->nextTracker;

        updateCursor->nextTracker = trackerRoot;
        trackerRoot = updateHolder;
        updateHolder = NULL;
        return;
    }
        
    TenXProgressTracker *parentRootCursor = trackerRoot;
    TenXProgressTracker *rootCursor = trackerRoot->nextTracker;
    TenXProgressTracker *parentUpdateCursor = NULL;
    TenXProgressTracker *updateCursor = updateHolder;

    // First update the new root, if necessary
    if (updateHolder->nextLocus >= trackerRoot->nextLocus)
        trackerRoot = updateHolder;

    // Move updateCursor to the first merge agent
    while (updateCursor != NULL && updateCursor->nextLocus >= parentRootCursor->nextLocus) {
        parentUpdateCursor = updateCursor;
        updateCursor = updateCursor->nextTracker;
    }
    // Fix the upstart link.
    if (parentUpdateCursor != NULL)
        parentUpdateCursor->nextTracker = parentRootCursor;

    /*
    //****10X debug
    int updateListSize = 1;
    TenXProgressTracker *updateCursorTemp = updateCursor;
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
    //****10X debug
    */

    // Merge updateList with rootList based on nextLocus in sorted order
    while (updateCursor != NULL && updateCursor->nextLocus != -1 && rootCursor != NULL && rootCursor->nextLocus != -1) {
    //while (updateCursor != NULL && rootCursor != NULL) {
        if (rootCursor->nextLocus < updateCursor->nextLocus) {
            parentRootCursor->nextTracker = updateCursor;
            // Remember that updateList is sorted!
            updateCursor = updateCursor->nextTracker;
            parentRootCursor->nextTracker->nextTracker = rootCursor;
            parentRootCursor = parentRootCursor->nextTracker;
            //updateListSize--; //****10X debug counter
        }
        else {
            parentRootCursor = rootCursor;
            rootCursor = rootCursor->nextTracker;
            //rootListSize--; //****10X debug counter
        }
    }

    // Repair the chain, in case it is terminated because of nextLocus == -1
    if (rootCursor == NULL) {
        parentRootCursor->nextTracker = updateCursor;
    }
    /* debug */
    else if (rootCursor->nextLocus == -1) {
        parentRootCursor->nextTracker = updateCursor;
        
        while (updateCursor->nextTracker != NULL)
            updateCursor = updateCursor->nextTracker;
        
        updateCursor->nextTracker = rootCursor;
    }
    else if (updateCursor != NULL && updateCursor->nextLocus == -1) {
        while (rootCursor->nextTracker != NULL && rootCursor->nextTracker->nextLocus != -1)
            rootCursor = rootCursor->nextTracker;

        TenXProgressTracker *nextRootCursor = rootCursor->nextTracker;
        
        rootCursor->nextTracker = updateCursor;
        
        while (updateCursor->nextTracker != NULL)
            updateCursor = updateCursor->nextTracker;

        updateCursor->nextTracker = nextRootCursor;
    }
    updateHolder = NULL;
}



// Progress each single aligner to move its locus pointer pass targetLoc, while registering the candidate with clusterIdx. The process stops BEFORE processing end. [start, end)
// It also terminates when cursor->nextLocus >= targetLocus, whichever comes first.
// PreStart is the node before start. We will need to repair the breaking link of prestart unless start is the root.
void TenXClusterAligner::registerClusterForReads(struct TenXProgressTracker* preStart, struct TenXProgressTracker* start, struct TenXProgressTracker* end, GenomeLocation clusterBoundary, int clusterIdx)
{
    if (preStart == NULL)
        _ASSERT(start == trackerRoot);

    TenXProgressTracker *cursor = start;
    while (cursor != end && cursor->nextLocus >= 0 && cursor->nextLocus > clusterBoundary) {
        cursor->aligner->align_phase_2_to_target_loc(clusterBoundary, clusterIdx);
        cursor->nextLocus = resolveLocusPtr(cursor->aligner->align_phase_2_get_locus() );

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
            //fprintf(stderr, "lastLocus: %lld\n", progressTracker[pairIdx].nextLocus.location);
            progressTracker[pairIdx].aligner->align_phase_2_to_target_loc(forwardLoc, -1);
            progressTracker[pairIdx].nextLocus = resolveLocusPtr(progressTracker[pairIdx].aligner->align_phase_2_get_locus() );
        }
        //fprintf(stderr, "forwardLocus: %lld  afterForward: %lld\n", forwardLoc.location, resolveLocusPtr(progressTracker[pairIdx].aligner->align_phase_2_get_locus() ).location);
    }
}


bool TenXClusterAligner::align_first_stage(
    unsigned        barcodeSize_)
{
    barcodeSize = barcodeSize_;
    bool barcodeFinished = true;
    maxClusterIdx = -1;

    for (unsigned pairIdx = 0; pairIdx < barcodeSize; pairIdx++) {
        clusterIsValid[pairIdx] = false;
        mappedPairsPerCluster[pairIdx] = 0;

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
                    progressTracker[pairIdx].nextLocus = resolveLocusPtr(progressTracker[pairIdx].aligner->align_phase_2_get_locus() );
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

    while (trackerRoot->pairNotDone && trackerRoot->nextLocus != -1) {
        // Initialization.
        cursor = trackerRoot;
        //expirationCursor = trackerRoot;

        unsigned nPotentialPairs = trackersToMeetTargetLocus(cursor, trackerRoot->nextLocus - minClusterSpan);

        if (nPotentialPairs > minPairsPerCluster || cursor == NULL || cursor->nextLocus == -1) { // this is a clustered pair,
        // tag it with clusterId! Note that if we are at the end (nextLocus is -1) we will add remaining pairs into the last
        // cluster (even though nPotentialPairs is not big enough).
            registeringCluster = true;
            clusterId = globalClusterId;

            if (cursor != NULL && cursor->nextLocus != -1)
                clusterBoundary = cursor->nextLocus + minClusterSpan;
            else
                clusterBoundary = -1;

            registerClusterForReads(NULL, trackerRoot, cursor, clusterBoundary, clusterId); //register the pairs and update the locus pointers.
        }
        else { //the cluster ends here.
            if (registeringCluster) { //when we were half way of adding a cluster, we need to finish it with the old targetLoc.
                //fprintf(stderr, "clusterBoundary: %lld    globalClusterId: %d\n", clusterBoundary.location, globalClusterId);
                //fflush(stderr);
                registerClusterForReads(NULL, trackerRoot, cursor, clusterBoundary, clusterId); //use the previous id.
                globalClusterId++;
                registeringCluster = false;

            }
            else { //we were not adding a cluster, just tag these locus as not clustered (-1) and move the locus pointer over the new targetLoc.
                cursor = traverseProgressPtr(cursor, minPairsPerCluster - nPotentialPairs);
                clusterBoundary = cursor->nextLocus + minClusterSpan;
                registerClusterForReads(NULL, trackerRoot, cursor, clusterBoundary, -1); //use the previous id.
            }
        }
        //fprintf(stderr, "clusterBoundary: %lld\n", clusterBoundary.location);
        //sortAndLink(); //fix the order. ****If we use this too often, it would be a potential performance problem. Might need to fix it later. (That's why I kept the linked list pointer!).
        //trackerRoot = &progressTracker[0];
        mergeUpdate();
    }

    return barcodeFinished;

}


void TenXClusterAligner::align_second_stage_clustering()
{
    bool barcodeFinished = true;
    for (unsigned pairIdx = 0; pairIdx < barcodeSize; pairIdx++) {
        if (progressTracker[pairIdx].pairNotDone) {
            Read *read0 = progressTracker[pairIdx].pairedReads[0];
            Read *read1 = progressTracker[pairIdx].pairedReads[1];

            progressTracker[pairIdx].nSingleEndSecondaryResults[0] = 0;
            progressTracker[pairIdx].nSingleEndSecondaryResults[1] = 0;

            progressTracker[pairIdx].bestPairScore = 65536;
            progressTracker[pairIdx].probabilityOfAllPairs = 0;
            progressTracker[pairIdx].probabilityOfBestPair = 0;
            progressTracker[pairIdx].bestClusterIdx = 0;

            progressTracker[pairIdx].aligner->align_phase_3_score(progressTracker[pairIdx].bestPairScore, false);
            progressTracker[pairIdx].aligner->align_phase_3_increment_cluster(progressTracker[pairIdx].bestPairScore);
        }
    }

    for (unsigned pairIdx = 0; pairIdx < barcodeSize; pairIdx++) {
        if (progressTracker[pairIdx].pairNotDone) {
            progressTracker[pairIdx].updatedBestScore = progressTracker[pairIdx].aligner->align_phase_3_correct_best_score(progressTracker[pairIdx].bestPairScore, minPairsPerCluster);
        }
    }

    // If the bestScore is updated, we have to add more potential candidates
    for (unsigned pairIdx = 0; pairIdx < barcodeSize; pairIdx++) {
        if (progressTracker[pairIdx].pairNotDone && progressTracker[pairIdx].updatedBestScore) {
            progressTracker[pairIdx].aligner->align_phase_3_score(progressTracker[pairIdx].bestPairScore, true);
            progressTracker[pairIdx].aligner->align_phase_3_increment_cluster(progressTracker[pairIdx].bestPairScore);
        }
    }
    
    // This time, best result can only change to be better
    for (unsigned pairIdx = 0; pairIdx < barcodeSize; pairIdx++) {
        if (progressTracker[pairIdx].pairNotDone && progressTracker[pairIdx].updatedBestScore) {
            progressTracker[pairIdx].aligner->align_phase_3_correct_best_score(progressTracker[pairIdx].bestPairScore, minPairsPerCluster);
        }
    }
}


bool TenXClusterAligner::align_second_stage_check_reallocate() 
{
    bool noOverflow = true;
    
    for (unsigned pairIdx = 0; pairIdx < barcodeSize; pairIdx++) {
        if (progressTracker[pairIdx].pairNotDone) {
            bool overflow = progressTracker[pairIdx].aligner->align_phase_3_count_results(maxEditDistanceForSecondaryResults,
            progressTracker[pairIdx].bestPairScore, minPairsPerCluster, &progressTracker[pairIdx].nSecondaryResults,
            progressTracker[pairIdx].secondaryResultBufferSize, progressTracker[pairIdx].probabilityOfAllPairs);
            
            if (overflow)
                noOverflow = false;
        }
    }
    
    return noOverflow;
}


void TenXClusterAligner::align_second_stage_generate_results()
{
    bool barcodeFinished = true;
    for (unsigned pairIdx = 0; pairIdx < barcodeSize; pairIdx++) {
        if (progressTracker[pairIdx].pairNotDone) {
            progressTracker[pairIdx].aligner->align_phase_3_generate_results(
                minPairsPerCluster,
                maxEditDistanceForSecondaryResults,
                progressTracker[pairIdx].bestPairScore,
                &progressTracker[pairIdx].nSecondaryResults,
                &progressTracker[pairIdx].results[1],
                &progressTracker[pairIdx].results[0]
            );
        }
    }
}


void TenXClusterAligner::clusterResultCleanUp()
{
    for (unsigned pairIdx = 0; pairIdx < barcodeSize; pairIdx++) {
        if (progressTracker[pairIdx].pairNotDone) {

            // Only examine the results if it has any
            if (progressTracker[pairIdx].nSecondaryResults > 0) {
                // Sort all results (best result might be corrected later)
                qsort(&progressTracker[pairIdx].results[1], progressTracker[pairIdx].nSecondaryResults, sizeof(PairedAlignmentResult), PairedAlignmentResult::compareByClusterIdx);

                int prevClusterIdx = -1;
                unsigned resultIdx;

                for (resultIdx = 1; resultIdx < progressTracker[pairIdx].nSecondaryResults + 1; resultIdx++) {
                    // We have a valid cluster and it is different from the bestClusterIdx (because bestClusterIdx is already registered)
                    if (progressTracker[pairIdx].results[resultIdx].clusterIdx != -1 && progressTracker[pairIdx].results[resultIdx].clusterIdx != progressTracker[pairIdx].bestClusterIdx) {
                        if (progressTracker[pairIdx].results[resultIdx].clusterIdx != prevClusterIdx)
                            mappedPairsPerCluster[progressTracker[pairIdx].results[resultIdx].clusterIdx]++;

                        prevClusterIdx = progressTracker[pairIdx].results[resultIdx].clusterIdx;
                    }

                    // Update the maxClusterIdx if necessary
                    if (progressTracker[pairIdx].results[resultIdx].clusterIdx > maxClusterIdx)
                        maxClusterIdx = progressTracker[pairIdx].results[resultIdx].clusterIdx;
                }
            }
        }
    }

    // Validate clusters
    for (unsigned clusterIdx = 0; clusterIdx < maxClusterIdx; clusterIdx++)
    {
        // Check > 0 because we don't care if a cluster has 0 reads (should not happen but might happen later on)
        if (mappedPairsPerCluster[clusterIdx] >= minPairsPerCluster) {
            clusterIsValid[clusterIdx] = true;
        }
    }
}


void TenXClusterAligner::align_third_stage()
{
    //First count the number of reads per cluster.
    for (unsigned pairIdx = 0; pairIdx < barcodeSize; pairIdx++) {
        Read *read0 = progressTracker[pairIdx].pairedReads[0];
        Read *read1 = progressTracker[pairIdx].pairedReads[1];

        progressTracker[pairIdx].aligner->align_phase_4(
            read0,
            read1,
            maxEditDistanceForSecondaryResults, 
            maxSecondaryAlignmentsToReturn,
            progressTracker[pairIdx].popularSeedsSkipped,
            progressTracker[pairIdx].bestPairScore,
            progressTracker[pairIdx].probabilityOfAllPairs,
            &progressTracker[pairIdx].nSecondaryResults,
            &progressTracker[pairIdx].results[1],
            &progressTracker[pairIdx].results[0]
        );

        /* timing no longer makes sence
        _int64 start = timeInNanos();
        _int64 end = timeInNanos();
        result[pairIdx]->nanosInAlignTogether = end - start;
        */

        progressTracker[pairIdx].results[0].nanosInAlignTogether = 0; //timing for individual read no longer makes sence
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


bool TenXClusterAligner::align_forth_stage(
)
{
    bool barcodeFinished = true;
    for (unsigned pairIdx = 0; pairIdx < barcodeSize; pairIdx++) {
        if (progressTracker[pairIdx].singleNotDone) {// && progressTracker[pairIdx].notDone) {
            Read *read0 = progressTracker[pairIdx].pairedReads[0];
            Read *read1 = progressTracker[pairIdx].pairedReads[1];

            Read *read[NUM_READS_PER_PAIR] = { read0, read1 };
            _int64 *resultCount[2] = { &progressTracker[pairIdx].nSingleEndSecondaryResults[0], &progressTracker[pairIdx].nSingleEndSecondaryResults[1] };

            *resultCount[0] = *resultCount[1] = 0;

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
    Read                    **pairedReads,
    unsigned                barcodeSize_,
    PairedAlignmentResult   **result,
    int                     maxEditDistanceForSecondaryResults,
    _int64                  *secondaryResultBufferSize,
    _int64                  *nSecondaryResults,
    _int64                  *singleSecondaryBufferSize,
    _int64                  maxSecondaryAlignmentsToReturn,
    _int64                  *nSingleEndSecondaryResults,
    SingleAlignmentResult   **singleEndSecondaryResults,    // Single-end secondary alignments for when the paired-end alignment didn't work properly
    unsigned                *popularSeedsSkipped
)
{
    /* this should be useless now
    align_init_stage(barcodeSize_);
    if (align_first_stage())
        return true;
    if (!align_second_stage(maxEditDistanceForSecondaryResults, maxSecondaryAlignmentsToReturn))
        return false;

    align_second_stage(maxEditDistanceForSecondaryResults, maxSecondaryAlignmentsToReturn);

    if (align_forth_stage(maxEditDistanceForSecondaryResults, maxSecondaryAlignmentsToReturn))
        return true;
        */
    return false;
}
