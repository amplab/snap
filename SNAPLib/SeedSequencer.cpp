/*++

Module Name:

    SeedSequencer.cpp

Abstract:

    Code for determining the order of seeds in a read.

Authors:

    Bill Bolosky, August, 2013

Environment:

    User mode service.

Revision History:

    
--*/

#include "stdafx.h"
#include "SeedSequencer.h"

static SeedSequencer *Sequencers[LargestSeedSize + 1];

void InitializeSeedSequencers()
{
    for (unsigned i = 1; i <= LargestSeedSize; i++) {
        Sequencers[i] = new SeedSequencer(i);
    }
}

SeedSequencer::SeedSequencer(unsigned i_seedSize) : seedSize(i_seedSize)
{
    offsets = new unsigned[seedSize];
    for (unsigned i = 0; i < seedSize; i++) {
        offsets[i] = 0;
    }

    if (seedSize == 1) return;  // Not that seed size 1 makes any sense or is in any way supported.  But it's in the array, so we generate it.

    struct WorkItem {
        unsigned lowerBound;
        unsigned upperBound;
        WorkItem    *next;
    };

    unsigned nFilledOffsets = 1;
    WorkItem *workItems = new WorkItem;
    WorkItem *tail = NULL;
    workItems->next = NULL;
    workItems->lowerBound = 1;
    workItems->upperBound = seedSize - 1;

    while (NULL != workItems) {
        WorkItem *itemToProcess = workItems;
        workItems = itemToProcess->next;
        if (NULL == workItems) {
            tail = NULL;
        }

        unsigned selectedLocation = (itemToProcess->lowerBound + itemToProcess->upperBound) / 2;
        _ASSERT(offsets[selectedLocation] == 0);
        offsets[selectedLocation] = nFilledOffsets;
        nFilledOffsets++;

        //
        // Add the upper half to the tail.  It will be the larger of the two if they're not equal, since / 2 rounds down.
        //
        if (itemToProcess->upperBound > selectedLocation) {
            WorkItem *upperItem = new WorkItem;
            upperItem->next = NULL;
            upperItem->upperBound = itemToProcess->upperBound;
            upperItem->lowerBound = selectedLocation + 1;

            //
            // Add it to the tail.
            //
            if (NULL != tail) {
                tail->next = upperItem;
                tail = upperItem;
            } else {
                _ASSERT(workItems == NULL);
                workItems = tail = upperItem;
            }
        }

        if (itemToProcess->lowerBound < selectedLocation) {
            _ASSERT(workItems != NULL); // We must have already added an upper half.
            itemToProcess->upperBound = selectedLocation - 1;
            itemToProcess->next = NULL;
            tail->next = itemToProcess;
            tail = itemToProcess;
        } else {
            delete itemToProcess;
        }
    }

    _ASSERT(nFilledOffsets == seedSize);
}

unsigned GetWrappedNextSeedToTest(unsigned seedLen, unsigned wrapCount) 
{
    _ASSERT(seedLen <= LargestSeedSize);
    return Sequencers[seedLen]->GetWrappedNextSeedToTest(wrapCount);
}
