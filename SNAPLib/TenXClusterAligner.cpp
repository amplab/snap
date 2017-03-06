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
	TenXAnchorTracker   *anchorTracker_,
    TenXMultiTracker    *multiTracker_,
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
    : multiTracker(multiTracker_), anchorTracker(anchorTracker_), unclusteredPenalty(unclusteredPenalty_), clusterEDCompensation(clusterEDCompensation_), minPairsPerCluster(minPairsPerCluster_), minClusterSpan(minClusterSpan_), forceSpacing(forceSpacing_), index(index_), minReadLength(minReadLength_), maxEditDistanceForSecondaryResults(maxEditDistanceForSecondaryResults_), maxSecondaryAlignmentsToReturn(maxSecondaryAlignmentsToReturn_)
{
    // Create single-end aligners.
    singleAligner = new (allocator) BaseAligner(index, maxHits, maxK, maxReadSize,
        maxSeedsFromCommandLine, seedCoverage, minWeightToCheck, extraSearchDepth, noUkkonen, noOrderedEvaluation, noTruncation, ignoreAlignmentAdjustmentsForOm, maxSecondaryAlignmentsPerContig, printStatsMapQLimit, &lv, &reverseLV, NULL, allocator);
    for (unsigned i = 0; i < multiPairNum; i++)
        multiTracker[i].aligner->setLandauVishkin(&lv, &reverseLV);

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
    qsort(multiTracker, multiPairNum, sizeof(TenXMultiTracker), TenXMultiTracker::compare);

    // Link all tracker entries based on their GenomeLocation.
    unsigned nActiveTrackers = 0;
    for (unsigned pairIdx = 1; pairIdx < multiPairNum; pairIdx++) {

        if (multiTracker[pairIdx - 1].pairNotDone && multiTracker[pairIdx - 1].nextLocus >= 0)
            nActiveTrackers++;

        multiTracker[pairIdx - 1].nextTracker = &multiTracker[pairIdx];
    }
    
    if (multiTracker[multiPairNum - 1].pairNotDone && multiTracker[multiPairNum - 1].nextLocus >= 0)
        nActiveTrackers++;

    trackerRoot = &multiTracker[0];
    updateHolder = NULL;
    //fprintf(stderr, "activeTrackers: %d\n", nActiveTrackers);
}

// Moves the cursor (cursor is modified!) to the first tracker that has a locus that's greater than the target.
unsigned trackersToMeetTargetLocus(
    TenXMultiTracker *&cursor,
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


TenXMultiTracker* traverseProgressPtr(TenXMultiTracker *cursor, unsigned steps) 
{
    while (steps > 0 && cursor->pairNotDone && cursor->nextLocus >= 0) {
        cursor = cursor->nextTracker;
        steps--;
    }
    return cursor;
}


void TenXClusterAligner::pushToUpdate(TenXMultiTracker *newUpdate)
{
    if (updateHolder == NULL || updateHolder->nextLocus <= newUpdate->nextLocus) {
        newUpdate->nextTracker = updateHolder;
        updateHolder = newUpdate;
        return;
    }

    // At least updateHolder has a bigger locus
    TenXMultiTracker* cursor = updateHolder;

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
        TenXMultiTracker *updateCursor = updateHolder;

        while (updateCursor->nextTracker != NULL)
            updateCursor = updateCursor->nextTracker;

        updateCursor->nextTracker = trackerRoot;
        trackerRoot = updateHolder;
        updateHolder = NULL;
        return;
    }
        
    TenXMultiTracker *parentRootCursor = trackerRoot;
    TenXMultiTracker *rootCursor = trackerRoot->nextTracker;
    TenXMultiTracker *parentUpdateCursor = NULL;
    TenXMultiTracker *updateCursor = updateHolder;

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
    TenXMultiTracker *updateCursorTemp = updateCursor;
    while (updateCursorTemp != NULL) {
        updateListSize++;
        updateCursorTemp = updateCursorTemp->nextTracker;
    }
    int rootListSize = 1;
    TenXMultiTracker *rootCursorTemp = rootCursor;
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

        TenXMultiTracker *nextRootCursor = rootCursor->nextTracker;
        
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
void TenXClusterAligner::registerClusterForReads(struct TenXMultiTracker* preStart, struct TenXMultiTracker* start, struct TenXMultiTracker* end, GenomeLocation clusterBoundary, int clusterIdx)
{
    if (preStart == NULL)
        _ASSERT(start == trackerRoot);

    TenXMultiTracker *cursor = start;
    while (cursor != end && cursor->nextLocus >= 0 && cursor->nextLocus > clusterBoundary) {

        /*/ debug
        string **id = cursor->id;
        if (clusterIdx == 0) {
            printf("clusterIdx is 0, read id is %s\n", id[0]->c_str());
            fflush(stdout);
        }

        if (*id[0] == "HISEQ-002:213:HCGHVADXX:2:1207:11698:6700" || *id[1] == "HISEQ-002:213:HCGHVADXX:2:1207:11698:6700")
            printf("found the debugging cursor. ClusterId: %d\n", clusterIdx);
        */// debug

        cursor->aligner->align_phase_2_to_target_loc(clusterBoundary, clusterIdx);
        cursor->nextLocus = resolveLocusPtr(cursor->aligner->align_phase_2_get_locus() );

        // A temporary holder, because cursor will be modified in pushToUpdate
        TenXMultiTracker *nextCursor = cursor->nextTracker;
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
    for (unsigned pairIdx = 0; pairIdx < multiPairNum; pairIdx++) {
        if (multiTracker[pairIdx].pairNotDone) {
            //fprintf(stderr, "lastLocus: %lld\n", multiTracker[pairIdx].nextLocus.location);
            multiTracker[pairIdx].aligner->align_phase_2_to_target_loc(forwardLoc, -1);
            multiTracker[pairIdx].nextLocus = resolveLocusPtr(multiTracker[pairIdx].aligner->align_phase_2_get_locus() );
        }
        //fprintf(stderr, "forwardLocus: %lld  afterForward: %lld\n", forwardLoc.location, resolveLocusPtr(multiTracker[pairIdx].aligner->align_phase_2_get_locus() ).location);
    }
}


bool TenXClusterAligner::align_first_stage(
    unsigned        anchorNum_,
    unsigned        multiPairNum_
	)
{
    multiPairNum = anchorNum_,
    multiPairNum = multiPairNum_;
    bool barcodeFinished = true;
    maxClusterIdx = -1;

    for (unsigned pairIdx = 0; pairIdx < multiPairNum; pairIdx++) {

        if (multiTracker[pairIdx].pairNotDone) {
            multiTracker[pairIdx].results[0].status[0] = multiTracker[pairIdx].results[0].status[1] = NotFound;

            Read *read0 = multiTracker[pairIdx].pairedReads[0];
            Read *read1 = multiTracker[pairIdx].pairedReads[1];

            if (read0->getDataLength() < minReadLength && read1->getDataLength() < minReadLength) {
                TRACE("Reads are both too short -- returning");
                // Update primary result
                for (int whichRead = 0; whichRead < NUM_READS_PER_PAIR; whichRead++) {
                    multiTracker[pairIdx].results[0].location[whichRead] = 0;
                    multiTracker[pairIdx].results[0].mapq[whichRead] = 0;
                    multiTracker[pairIdx].results[0].score[whichRead] = 0;
                    multiTracker[pairIdx].results[0].status[whichRead] = NotFound;
                }
                multiTracker[pairIdx].results[0].alignedAsPair = false;
                multiTracker[pairIdx].results[0].fromAlignTogether = false;
                multiTracker[pairIdx].results[0].nanosInAlignTogether = 0;
                multiTracker[pairIdx].results[0].nLVCalls = 0;
                multiTracker[pairIdx].results[0].nSmallHits = 0;

                multiTracker[pairIdx].pairNotDone = false;
                multiTracker[pairIdx].singleNotDone = false;
                continue;//return true;
            }

            //At least one read of the pair is worthy of further examination
            barcodeFinished = false;

            if (read0->getDataLength() >= minReadLength && read1->getDataLength() >= minReadLength) {
                //
                // Let the LVs use the cache that we built up.
                //
                multiTracker[pairIdx].pairNotDone = !multiTracker[pairIdx].aligner->align_phase_1(read0, read1, multiTracker[pairIdx].popularSeedsSkipped);
                
                // Initialize for phase_2 if the alginer if not stopped
                if (multiTracker[pairIdx].pairNotDone) {
                    multiTracker[pairIdx].pairNotDone = multiTracker[pairIdx].aligner->align_phase_2_init();
                    multiTracker[pairIdx].nextLocus = resolveLocusPtr(multiTracker[pairIdx].aligner->align_phase_2_get_locus() );
                    multiTracker[pairIdx].nextTracker = NULL;
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
    TenXMultiTracker *cursor; // pointer that keeps track of the progress of walking down progress trackers

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
                printf("globalClusterId++: %d\n", globalClusterId);
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
        //trackerRoot = &multiTracker[0];
        mergeUpdate();
    }

    return barcodeFinished;

}


void TenXClusterAligner::align_second_stage_clustering()
{
    bool barcodeFinished = true;
    for (unsigned pairIdx = 0; pairIdx < multiPairNum; pairIdx++) {
        if (multiTracker[pairIdx].pairNotDone) {
            Read *read0 = multiTracker[pairIdx].pairedReads[0];
            Read *read1 = multiTracker[pairIdx].pairedReads[1];

            multiTracker[pairIdx].nSingleEndSecondaryResults[0] = 0;
            multiTracker[pairIdx].nSingleEndSecondaryResults[1] = 0;

            multiTracker[pairIdx].bestPairScore = 65536;
            multiTracker[pairIdx].probabilityOfAllPairs = 0;
            multiTracker[pairIdx].probabilityOfBestPair = 0;
            multiTracker[pairIdx].bestClusterIdx = -1;

            multiTracker[pairIdx].aligner->align_phase_3_score(multiTracker[pairIdx].bestPairScore, false);
            multiTracker[pairIdx].aligner->align_phase_3_increment_cluster(multiTracker[pairIdx].bestPairScore);
        }
    }

    for (unsigned pairIdx = 0; pairIdx < multiPairNum; pairIdx++) {
        if (multiTracker[pairIdx].pairNotDone) {
            multiTracker[pairIdx].updatedBestScore = multiTracker[pairIdx].aligner->align_phase_3_correct_best_score(multiTracker[pairIdx].bestPairScore, minPairsPerCluster);
        }
    }

    // If the bestScore is updated, we have to add more potential candidates
    for (unsigned pairIdx = 0; pairIdx < multiPairNum; pairIdx++) {
        if (multiTracker[pairIdx].pairNotDone && multiTracker[pairIdx].updatedBestScore) {
            multiTracker[pairIdx].aligner->align_phase_3_score(multiTracker[pairIdx].bestPairScore, true);
            multiTracker[pairIdx].aligner->align_phase_3_increment_cluster(multiTracker[pairIdx].bestPairScore);
        }
    }
    
    // This time, best result can only change to be better
    for (unsigned pairIdx = 0; pairIdx < multiPairNum; pairIdx++) {
        if (multiTracker[pairIdx].pairNotDone && multiTracker[pairIdx].updatedBestScore) {
            multiTracker[pairIdx].aligner->align_phase_3_correct_best_score(multiTracker[pairIdx].bestPairScore, minPairsPerCluster);
        }
    }
}


bool TenXClusterAligner::align_second_stage_check_reallocate() 
{
    bool noOverflow = true;
    
    for (unsigned pairIdx = 0; pairIdx < multiPairNum; pairIdx++) {
        if (multiTracker[pairIdx].pairNotDone) {
            bool overflow = multiTracker[pairIdx].aligner->align_phase_3_count_results(maxEditDistanceForSecondaryResults,
            multiTracker[pairIdx].bestPairScore, minPairsPerCluster, &multiTracker[pairIdx].nSecondaryResults,
            multiTracker[pairIdx].secondaryResultBufferSize, multiTracker[pairIdx].probabilityOfAllPairs);
            
            if (overflow)
                noOverflow = false;
        }
    }
    
    return noOverflow;
}


void TenXClusterAligner::align_second_stage_generate_results()
{
    bool barcodeFinished = true;
    for (unsigned pairIdx = 0; pairIdx < multiPairNum; pairIdx++) {
        if (multiTracker[pairIdx].pairNotDone) {
            multiTracker[pairIdx].aligner->align_phase_3_generate_results(
                minPairsPerCluster,
                maxEditDistanceForSecondaryResults,
                multiTracker[pairIdx].bestPairScore,
                &multiTracker[pairIdx].nSecondaryResults,
                &multiTracker[pairIdx].results[1],
                &multiTracker[pairIdx].results[0]
            );
        }
    }
}


void TenXClusterAligner::align_third_stage()
{
    //First count the number of reads per cluster.
    for (unsigned pairIdx = 0; pairIdx < multiPairNum; pairIdx++) {
        Read *read0 = multiTracker[pairIdx].pairedReads[0];
        Read *read1 = multiTracker[pairIdx].pairedReads[1];

        multiTracker[pairIdx].aligner->align_phase_4(
            read0,
            read1,
            maxEditDistanceForSecondaryResults, 
            maxSecondaryAlignmentsToReturn,
            multiTracker[pairIdx].popularSeedsSkipped,
            multiTracker[pairIdx].bestPairScore,
            multiTracker[pairIdx].probabilityOfAllPairs,
            &multiTracker[pairIdx].nSecondaryResults,
            &multiTracker[pairIdx].results[1],
            &multiTracker[pairIdx].results[0]
        );

        /* timing no longer makes sence
        _int64 start = timeInNanos();
        _int64 end = timeInNanos();
        result[pairIdx]->nanosInAlignTogether = end - start;
        */

        multiTracker[pairIdx].results[0].nanosInAlignTogether = 0; //timing for individual read no longer makes sence
        multiTracker[pairIdx].results[0].fromAlignTogether = true;
        multiTracker[pairIdx].results[0].alignedAsPair = true;

        if (forceSpacing) {
            if (multiTracker[pairIdx].results[0].status[0] == NotFound) {
                multiTracker[pairIdx].results[0].fromAlignTogether = false;
            }
            else {
                _ASSERT(multiTracker[pairIdx].results[0].status[1] != NotFound); // If one's not found, so is the other
            }
            multiTracker[pairIdx].pairNotDone = false;
            multiTracker[pairIdx].singleNotDone = false;
            continue;//return true;
        }

        if (multiTracker[pairIdx].results[0].status[0] != NotFound && multiTracker[pairIdx].results[0].status[1] != NotFound) {
            //
            // Not a chimeric read.
            //
            multiTracker[pairIdx].pairNotDone = false;
            multiTracker[pairIdx].singleNotDone = false;
            continue;//return true;
        }
        //paired analysis is done anyways
        multiTracker[pairIdx].pairNotDone = false;
    }
}


bool TenXClusterAligner::align_forth_stage(
)
{
    bool barcodeFinished = true;
    for (unsigned pairIdx = 0; pairIdx < multiPairNum; pairIdx++) {
        if (multiTracker[pairIdx].singleNotDone) {// && multiTracker[pairIdx].notDone) {
            // Tag single mapping's clusterIdx to -1
            multiTracker[pairIdx].results[0].clusterIdx = -2;

            Read *read0 = multiTracker[pairIdx].pairedReads[0];
            Read *read1 = multiTracker[pairIdx].pairedReads[1];

            Read *read[NUM_READS_PER_PAIR] = { read0, read1 };
            _int64 *resultCount[2] = { &multiTracker[pairIdx].nSingleEndSecondaryResults[0], &multiTracker[pairIdx].nSingleEndSecondaryResults[1] };

            *resultCount[0] = *resultCount[1] = 0;

            bool noOverflow = true;
            for (int r = 0; r < NUM_READS_PER_PAIR; r++) {
                SingleAlignmentResult singleResult;
                _int64 singleEndSecondaryResultsThisTime = 0;

                if (read[r]->getDataLength() < minReadLength) {
                    multiTracker[pairIdx].results[0].status[r] = NotFound;
                    multiTracker[pairIdx].results[0].mapq[r] = 0;
                    multiTracker[pairIdx].results[0].direction[r] = FORWARD;
                    multiTracker[pairIdx].results[0].location[r] = 0;
                    multiTracker[pairIdx].results[0].score[r] = 0;
                }
                else {
                    // We're using *nSingleEndSecondaryResultsForFirstRead because it's either 0 or what all we've seen (i.e., we know NUM_READS_PER_PAIR is 2)
                    bool fitInSecondaryBuffer = //true;
                        singleAligner->AlignRead(read[r], &singleResult, maxEditDistanceForSecondaryResults,
                            multiTracker[pairIdx].singleSecondaryBufferSize - multiTracker[pairIdx].nSingleEndSecondaryResults[0], &singleEndSecondaryResultsThisTime,
                            maxSecondaryAlignmentsToReturn, &multiTracker[pairIdx].singleEndSecondaryResults[multiTracker[pairIdx].nSingleEndSecondaryResults[0]]);

                    if (!fitInSecondaryBuffer) {
                        multiTracker[pairIdx].nSecondaryResults = 0;
                        multiTracker[pairIdx].nSingleEndSecondaryResults[0] = multiTracker[pairIdx].singleSecondaryBufferSize + 1;
                        multiTracker[pairIdx].nSingleEndSecondaryResults[1] = 0;
                        barcodeFinished = false;
                        noOverflow = false;
                        break;//return false;
                    }

                    *(resultCount[r]) = singleEndSecondaryResultsThisTime;

                    multiTracker[pairIdx].results[0].status[r] = singleResult.status;
                    multiTracker[pairIdx].results[0].mapq[r] = singleResult.mapq / 3;   // Heavy quality penalty for chimeric reads
                    multiTracker[pairIdx].results[0].direction[r] = singleResult.direction;
                    multiTracker[pairIdx].results[0].location[r] = singleResult.location;
                    multiTracker[pairIdx].results[0].score[r] = singleResult.score;
                    multiTracker[pairIdx].results[0].scorePriorToClipping[r] = singleResult.scorePriorToClipping;
    
                }
            }

            // This pair is done processing, only if both directions has no overflow.
            if (noOverflow) {
                multiTracker[pairIdx].singleNotDone = false;
                multiTracker[pairIdx].results[0].fromAlignTogether = false;
                multiTracker[pairIdx].results[0].alignedAsPair = false;
            }

#ifdef _DEBUG
            if (_DumpAlignments) {
                printf("TenXClusterAligner: (%u, %u) score (%d, %d), MAPQ (%d, %d)\n\n\n", multiTracker[pairIdx].results[0].location[0].location, multiTracker[pairIdx].results[0].location[1].location,
                    multiTracker[pairIdx].results[0].score[0], multiTracker[pairIdx].results[0].score[1], multiTracker[pairIdx].results[0].mapq[0], multiTracker[pairIdx].results[0].mapq[1]);
            }
#endif // _DEBUG

        }
    }
    return barcodeFinished;
}



bool TenXClusterAligner::align(
    Read                    **pairedReads,
    unsigned                multiPairNum_,
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
    align_init_stage(multiPairNum_);
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
