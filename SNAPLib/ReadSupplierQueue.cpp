/*++

Module Name:

    ReadSupplierQueue.cpp

Abstract:

    Code for parallel queue of reads

Authors:

    Bill Bolosky, November, 2012

Environment:

    User mode service.

Revision History:


--*/

#include "stdafx.h"
#include "Read.h"
#include "Compat.h"
#include "ReadSupplierQueue.h"

 ReadSupplierQueue::ReadSupplierQueue(int i_nReaders, ReadReader **readers)
 {
     commonInit(i_nReaders);

     for (int i = 0; i < nReaders; i++) {
         readerGroups[i].singleReader[0] = readers[i];
     }
 }

ReadSupplierQueue::ReadSupplierQueue(int i_nReaders, ReadReader **firstHalfReaders, ReadReader **secondHalfReaders)
{
     commonInit(i_nReaders);

     for (int i = 0; i < nReaders; i++) {
         readerGroups[i].singleReader[0] = firstHalfReaders[i];
         readerGroups[i].singleReader[1] = secondHalfReaders[i];
     }
     nReadersRunning = 2 * nReaders;
}

ReadSupplierQueue::ReadSupplierQueue(int i_nReaders, PairedReadReader **pairedReaders)
{
     commonInit(i_nReaders);

     for (int i = 0; i < nReaders; i++) {
         readerGroups[i].pairedReader = pairedReaders[i];
     }
}

void
ReadSupplierQueue::commonInit(int i_nReaders)
{
    nReaders = i_nReaders;
    readerGroups = new ReaderGroup[nReaders];
    nReadersRunning = nReaders;
    nSuppliersRunning = 0;
    allReadsQueued = false;

    emptyQueue->next = emptyQueue->prev = emptyQueue;
    readerGroupsWithReadyReads->next = readerGroupsWithReadyReads->prev = readerGroupsWithReadyReads;

    InitializeExclusiveLock(&lock);
    CreateSingleWaiterObject(&readsReady);
    CreateSingleWaiterObject(&emptyBuffersAvailable);
    CreateSingleWaiterObject(&allReadsConsumed);

    //
    // Create 2 buffers per reader.  We'll add more buffers as we add suppliers.
    //
    for (int i = 0 ; i < nReaders * 2; i++) {
        ReadQueueElement *element = new ReadQueueElement;
        element->addToTail(emptyQueue);
    }

    SignalSingleWaiterObject(&emptyBuffersAvailable);
}

ReadSupplierQueue::~ReadSupplierQueue()
{
    for (int i = 0; i < nReaders; i++) {
        delete readerGroups[i].singleReader[0];
        delete readerGroups[i].singleReader[1];
        delete readerGroups[i].pairedReader;
    }
    delete [] readerGroups;
}


    bool 
ReadSupplierQueue::startReaders()
{
    bool worked = true;
    for (int i = 0; i < nReaders; i++) {
        ReaderThreadParams *readerParams = new ReaderThreadParams;
        readerParams->group = &readerGroups[i];
        readerParams->isSecondReader = false;
        readerParams->queue = this;
        worked &= StartNewThread(ReaderThreadMain, readerParams);

        if (readerGroups[i].singleReader[1] != NULL) {
            readerParams = new ReaderThreadParams;
            readerParams->group = &readerGroups[i];
            readerParams->isSecondReader = true;
            readerParams->queue = this;
            worked &= StartNewThread(ReaderThreadMain, readerParams);            
        }
    }

    return worked;
}

    void 
ReadSupplierQueue::waitUntilFinished()
{
    WaitForSingleWaiterObject(&allReadsConsumed);
}


    ReadSupplierFromQueue *
ReadSupplierQueue::createSupplier()
{
    AcquireExclusiveLock(&lock);
    nSuppliersRunning++;
    //
    // Add two more queue elements for this supplier.
    //
    for (int i = 0; i < 2; i++) {
        ReadQueueElement *element = new ReadQueueElement;
        element->addToTail(emptyQueue);
    }

    SignalSingleWaiterObject(&emptyBuffersAvailable);
    ReleaseExclusiveLock(&lock);
   
    return new ReadSupplierFromQueue(this);
}

        PairedReadSupplierFromQueue *
ReadSupplierQueue::createPairedSupplier()
{
    AcquireExclusiveLock(&lock);
    nSuppliersRunning++;
    //
    // Add two more queue elements (four for paired-end, double file).
    //
    for (int i = 0; i < ((readerGroups[0].singleReader[1] == NULL) ? 2 : 4); i++) {
        ReadQueueElement *element = new ReadQueueElement;
        element->addToTail(emptyQueue);
    }

    SignalSingleWaiterObject(&emptyBuffersAvailable);
    ReleaseExclusiveLock(&lock);
   
    return new PairedReadSupplierFromQueue(this, readerGroups[0].singleReader[1] != NULL);
}

    ReadQueueElement *
ReadSupplierQueue::getElement()
{
    AcquireExclusiveLock(&lock);
    while (readerGroupsWithReadyReads->next == readerGroupsWithReadyReads) {
        ReleaseExclusiveLock(&lock);
        if (allReadsQueued) {
            //
            // Everything's queued and the queue is empty.  No more work.
            //
            return NULL;
        }
        WaitForSingleWaiterObject(&readsReady);
        AcquireExclusiveLock(&lock);
    }

    ReaderGroup *group = readerGroupsWithReadyReads->next;
    _ASSERT(group->readyQueue->next != group->readyQueue);  // There are reads ready.
    ReadQueueElement *element = group->readyQueue->next;
    element->removeFromQueue();
    if (group->readyQueue->next == group->readyQueue) {
        //
        // No reads left in this group.
        //
        group->removeFromQueue();
 
        if (readerGroupsWithReadyReads->next == readerGroupsWithReadyReads) {
            //
            // No groups left with reads ready.
            //
            ResetSingleWaiterObject(&readsReady);
        }
    } else {
        //
        // Move the group to the end of the queue, so we don't starve readers.
        //
        group->removeFromQueue();
        group->addToQueue(readerGroupsWithReadyReads);
    }
 
    ReleaseExclusiveLock(&lock);

    return element;
}

        bool 
ReadSupplierQueue::getElements(ReadQueueElement **element1, ReadQueueElement **element2)
{
    AcquireExclusiveLock(&lock);
    while (readerGroupsWithReadyReads->next == readerGroupsWithReadyReads) {
        ReleaseExclusiveLock(&lock);
        if (allReadsQueued) {
            //
            // Everything's queued and the queue is empty.  No more work.
            //
            return NULL;
        }
        WaitForSingleWaiterObject(&readsReady);
        AcquireExclusiveLock(&lock);
    }

    ReaderGroup *group = readerGroupsWithReadyReads->next;
    _ASSERT(group->readyQueue->next != group->readyQueue);  // There are reads ready.
    _ASSERT(group->readyQueue[1].next != &group->readyQueue[1]); // From both readers
    *element1 = group->readyQueue->next;
    *element2 = group->readyQueue[1].next;
    (*element1)->removeFromQueue();
    (*element2)->removeFromQueue();
    if (group->readyQueue->next == group->readyQueue || group->readyQueue[1].next == &group->readyQueue[1]) {
        //
        // No reads left in this group.
        //
        group->removeFromQueue();
 
        if (readerGroupsWithReadyReads->next == readerGroupsWithReadyReads) {
            //
            // No groups left with reads ready.
            //
            ResetSingleWaiterObject(&readsReady);
        }
    } else {
        //
        // Move the group to the end of the queue, so we don't starve readers.
        //
        group->removeFromQueue();
        group->addToQueue(readerGroupsWithReadyReads);
    }
 
    ReleaseExclusiveLock(&lock);
    return true;
}

    void 
ReadSupplierQueue::doneWithElement(ReadQueueElement *element)
{
    AcquireExclusiveLock(&lock);
    element->addToTail(emptyQueue);
    SignalSingleWaiterObject(&emptyBuffersAvailable);
    ReleaseExclusiveLock(&lock);

}
    void 
ReadSupplierQueue::supplierFinished()
{
    AcquireExclusiveLock(&lock);
    _ASSERT(allReadsQueued);
    _ASSERT(nSuppliersRunning > 0);
    nSuppliersRunning--;
    if (0 == nSuppliersRunning) {
        SignalSingleWaiterObject(&allReadsConsumed);
    }
    ReleaseExclusiveLock(&lock);
}

    void
ReadSupplierQueue::ReaderThreadMain(void *param)
{
    ReaderThreadParams *params = (ReaderThreadParams *)param;
    params->queue->ReaderThread(params);
    delete params;
}

    void
ReadSupplierQueue::ReaderThread(ReaderThreadParams *params)
{
    AcquireExclusiveLock(&lock);
    bool done = false;
    ReadReader *reader;
    if (params->isSecondReader) { // In the pairedReader case, this will just be NULL
        reader = params->group->singleReader[1];
    } else {
        reader = params->group->singleReader[0];
    }
    int increment = (NULL == reader) ? 2 : 1;

    while (!done) {
        while (emptyQueue->next == emptyQueue) {
            // Wait for a buffer.
            ReleaseExclusiveLock(&lock);
            WaitForSingleWaiterObject(&emptyBuffersAvailable);
            AcquireExclusiveLock(&lock);
        }

        ReadQueueElement *element = emptyQueue->next;
        element->removeFromQueue();
        if (emptyQueue->next == emptyQueue) {
            ResetSingleWaiterObject(&emptyBuffersAvailable);
        }
        ReleaseExclusiveLock(&lock);

        //
        // Now fill in the reads from the reader into the element until it's
        // full or the reader finishes.
        //
        for (element->totalReads = 0; element->totalReads < element->nReads; element->totalReads += increment) {
            if (NULL != reader) {
                if (!reader->getNextRead(&element->reads[element->totalReads])) {
                   done = true;
                   break;
                } 
            } else if (NULL != params->group->pairedReader) {
                if (!params->group->pairedReader->getNextReadPair(&element->reads[element->totalReads], &element->reads[element->totalReads+1])) {
                    done = true;
                    break;
                }
           }
        }
        
        AcquireExclusiveLock(&lock);

        if (element->totalReads > 0) {
            element->addToTail(&params->group->readyQueue[params->isSecondReader ? 1 : 0]);
            if (params->group->next == NULL && 
                (NULL == params->group->singleReader[1] || params->group->readyQueue[params->isSecondReader ? 0 : 1].next != &params->group->readyQueue[params->isSecondReader ? 0 : 1])) {
 
                //
                // Our group is not already on the ready queue.  Furthermore, we now have data ready, either because we
                // are a single file reader or because the queue for the other file is non-empty.  Add us to the
                // ready queue and signal that there's data ready.
                //
                params->group->addToQueue(readerGroupsWithReadyReads);
                SignalSingleWaiterObject(&readsReady);
            }
        }
    } // While ! done
    _ASSERT(nReadersRunning > 0);
    nReadersRunning--;
    if (0 == nReadersRunning) {
        allReadsQueued = true;
    }
    ReleaseExclusiveLock(&lock);
}

ReadSupplierFromQueue::ReadSupplierFromQueue(ReadSupplierQueue *i_queue) : queue(i_queue), outOfReads(false), currentElement(NULL), nextReadIndex(0), done(false)
{
}

Read *
ReadSupplierFromQueue::getNextRead()
{
    if (done) {
        return false;
    }

    if (NULL != currentElement && nextReadIndex >= currentElement->totalReads) {
        queue->doneWithElement(currentElement);
        currentElement = NULL;
    }

    if (NULL == currentElement) {
        currentElement = queue->getElement();
        if (NULL == currentElement) {
            done = true;
            queue->supplierFinished();
            return NULL;
        }
        nextReadIndex = 0;
    }

    return &currentElement->reads[nextReadIndex++]; // Note the post increment.
}

PairedReadSupplierFromQueue::PairedReadSupplierFromQueue(ReadSupplierQueue *i_queue, bool i_twoFiles) : queue(i_queue), twoFiles(i_twoFiles), done(false), 
    currentElement(NULL), currentSecondElement(NULL), nextReadIndex(0) {}

    bool
PairedReadSupplierFromQueue::getNextReadPair(Read **read0, Read **read1)
{
    if (done) {
        *read0 = NULL;
        *read1 = NULL;
        return false;
    }

    if (NULL != currentElement && nextReadIndex >= currentElement->totalReads) {
        queue->doneWithElement(currentElement);
        currentElement = NULL;
        if (twoFiles) {
            queue->doneWithElement(currentSecondElement);
        }
        currentSecondElement = NULL;
    }

    if (NULL == currentElement) {
        if ((twoFiles && !queue->getElements(&currentElement, &currentSecondElement)) || 
            (!twoFiles && NULL == (currentElement = queue->getElement()))) {

            done = true;
            queue->supplierFinished();
            *read0 = NULL;
            *read1 = NULL;
            return false;
        }
        if (twoFiles) {
            // Assert that both elements match.
            _ASSERT(currentSecondElement->totalReads == currentElement->totalReads);
        } else {
            //
            // Assert that there are an even number of reads (since they're in pairs)
            //
            _ASSERT(currentElement->totalReads % 2 == 0);
        }
        nextReadIndex = 0;
    }

    if (twoFiles) {
        *read0 = &currentElement->reads[nextReadIndex];
        *read1 = &currentSecondElement->reads[nextReadIndex];
        nextReadIndex++;
    } else {
        *read0 = &currentElement->reads[nextReadIndex];
        *read1 = &currentElement->reads[nextReadIndex+1];
        nextReadIndex += 2;
    }

    return true;
}
