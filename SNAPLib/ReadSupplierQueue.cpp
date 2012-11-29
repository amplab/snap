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

 ReadSupplierQueue::ReadSupplierQueue(ReadReader *reader)
 {
     commonInit();

     singleReader[0] = reader;
 }

ReadSupplierQueue::ReadSupplierQueue(ReadReader *firstHalfReader, ReadReader *secondHalfReader)
{
     commonInit();

     singleReader[0] = firstHalfReader;
     singleReader[1] = secondHalfReader;
}

ReadSupplierQueue::ReadSupplierQueue(PairedReadReader *i_pairedReader)
{
     commonInit();
    pairedReader = i_pairedReader;
}

void
ReadSupplierQueue::commonInit()
{
    nReadersRunning = 0;
    nSuppliersRunning = 0;
    allReadsQueued = false;

    emptyQueue->next = emptyQueue->prev = emptyQueue;
    readyQueue[0].next = readyQueue[0].prev = &readyQueue[0];
    readyQueue[1].next = readyQueue[1].prev = &readyQueue[1];

    InitializeExclusiveLock(&lock);
    CreateEventObject(&readsReady);
    CreateEventObject(&emptyBuffersAvailable);
    CreateEventObject(&allReadsConsumed);

    //
    // Create 2 buffers for the reader.  We'll add more buffers as we add suppliers.
    //
    for (int i = 0 ; i < 2; i++) {
        ReadQueueElement *element = new ReadQueueElement;
        element->addToTail(emptyQueue);
    }

    AllowEventWaitersToProceed(&emptyBuffersAvailable);

    for (int i = 0; i < 2; i++) {
        CreateEventObject(&throttle[i]);
        AllowEventWaitersToProceed(&throttle[i]);
        singleReader[i] = NULL;
    }
    pairedReader = NULL;
}

ReadSupplierQueue::~ReadSupplierQueue()
{
    delete singleReader[0];
    delete singleReader[1];
    delete pairedReader;

    DestroyEventObject(&throttle[0]);
    DestroyEventObject(&throttle[1]);
}


    bool 
ReadSupplierQueue::startReaders()
{
    bool worked = true;

    if (singleReader[1] == NULL) {
        nReadersRunning = 1;
    } else {
        nReadersRunning = 2;
    }

    ReaderThreadParams *readerParams = new ReaderThreadParams;
    readerParams->isSecondReader = false;
    readerParams->queue = this;
    if (!StartNewThread(ReaderThreadMain, readerParams)) {
        return false;
    }

    if (singleReader[1] == NULL) {
        return true;
    }

    readerParams = new ReaderThreadParams;
    readerParams->isSecondReader = true;
    readerParams->queue = this;
    return (StartNewThread(ReaderThreadMain, readerParams));
}

    void 
ReadSupplierQueue::waitUntilFinished()
{
    WaitForEvent(&allReadsConsumed);
}


    ReadSupplier *
ReadSupplierQueue::generateNewReadSupplier()
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

    AllowEventWaitersToProceed(&emptyBuffersAvailable);
    ReleaseExclusiveLock(&lock);
   
    return new ReadSupplierFromQueue(this);
}

        PairedReadSupplier *
ReadSupplierQueue::generateNewPairedReadSupplier()
{
    ReadQueueElement * newElements[4];
    for (int i = 0 ; i < ((singleReader[1] == NULL) ? 2 : 4); i++) {
        newElements[i] = new ReadQueueElement;
    }

    AcquireExclusiveLock(&lock);
    nSuppliersRunning++;
    //
    // Add two more queue elements (four for paired-end, double file).
    //
    for (int i = 0; i < ((singleReader[1] == NULL) ? 2 : 4); i++) {
        ReadQueueElement *element = newElements[i];
        element->addToTail(emptyQueue);
    }

    AllowEventWaitersToProceed(&emptyBuffersAvailable);
    ReleaseExclusiveLock(&lock);
   
    return new PairedReadSupplierFromQueue(this, singleReader[1] != NULL);
}

    ReadQueueElement *
ReadSupplierQueue::getElement()
{
    _ASSERT(singleReader[1] == NULL);   // i.e., we're doing file (but possibly single or paired end) reads
    AcquireExclusiveLock(&lock);
    while (!areAnyReadsReady()) {
        ReleaseExclusiveLock(&lock);
        if (allReadsQueued) {
            //
            // Everything's queued and the queue is empty.  No more work.
            //
            return NULL;
        }
        WaitForEvent(&readsReady);
        AcquireExclusiveLock(&lock);
    }

    ReadQueueElement *element = readyQueue[0].next;
    _ASSERT(element != &readyQueue[0]);
    element->removeFromQueue();

    if (!areAnyReadsReady()) {
        PreventEventWaitersFromProceeding(&readsReady);
    }
 
    ReleaseExclusiveLock(&lock);

    return element;
}

        bool 
ReadSupplierQueue::getElements(ReadQueueElement **element1, ReadQueueElement **element2)
{
   _ASSERT(singleReader[1] != NULL);   // i.e., we're doing paired file reads

   AcquireExclusiveLock(&lock);
    while (!areAnyReadsReady()) {
        ReleaseExclusiveLock(&lock);
        if (allReadsQueued) {
            //
            // Everything's queued and the queue is empty.  No more work.
            //
            return NULL;
        }
        WaitForEvent(&readsReady);
        AcquireExclusiveLock(&lock);
    }

    *element1 = readyQueue[0].next;
    *element2 = readyQueue[1].next;
    (*element1)->removeFromQueue();
    (*element2)->removeFromQueue();

    if (!areAnyReadsReady()) {
        PreventEventWaitersFromProceeding(&readsReady);
    }
 
    ReleaseExclusiveLock(&lock);
    return true;
}

    bool
ReadSupplierQueue::areAnyReadsReady() // must hold the lock to call this.
{
    if (readyQueue[0].next == &readyQueue[0]) {
        return false;
    }

    if (singleReader[1] == NULL) {
        return true;
    }

    return readyQueue[1].next != &readyQueue[1];
}

    void 
ReadSupplierQueue::doneWithElement(ReadQueueElement *element)
{
    AcquireExclusiveLock(&lock);
    element->addToTail(emptyQueue);
    AllowEventWaitersToProceed(&emptyBuffersAvailable);
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
        AllowEventWaitersToProceed(&allReadsConsumed);
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
    if (params->isSecondReader) { 
        reader = singleReader[1];
    } else {
        reader = singleReader[0]; // In the pairedReader case, this will just be NULL
    }
    int increment = (NULL == reader) ? 2 : 1;
    int balanceIncrement = params->isSecondReader ? -1 : 1;
    int firstOrSecond = params->isSecondReader ? 0 : 1;
    bool isSingleReader = (NULL == singleReader[1]);

    while (!done) {
        if ((!isSingleReader) && balance * balanceIncrement > MaxImbalance) {
            //
            // We're over full.  Wait to get back in balance.
            //
            ReleaseExclusiveLock(&lock);
            WaitForEvent(&throttle[firstOrSecond]);
            AcquireExclusiveLock(&lock);
            _ASSERT(balance * balanceIncrement <= MaxImbalance);
        }

        while (emptyQueue->next == emptyQueue) {
            // Wait for a buffer.
            ReleaseExclusiveLock(&lock);
            WaitForEvent(&emptyBuffersAvailable);
            AcquireExclusiveLock(&lock);
        }

        ReadQueueElement *element = emptyQueue->next;
        element->removeFromQueue();
        if (emptyQueue->next == emptyQueue) {
            PreventEventWaitersFromProceeding(&emptyBuffersAvailable);
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
            } else if (NULL != pairedReader) {
                if (!pairedReader->getNextReadPair(&element->reads[element->totalReads], &element->reads[element->totalReads+1])) {
                    done = true;
                    break;
                }
           }
        }
        
        AcquireExclusiveLock(&lock);

        if (element->totalReads > 0) {
            element->addToTail(&readyQueue[firstOrSecond]);
            if (isSingleReader || &readyQueue[1-firstOrSecond] != readyQueue[1-firstOrSecond].next) {
                //
                // Signal that an element is ready.
                //
                AllowEventWaitersToProceed(&readsReady);
            }

            if (!isSingleReader) {
                balance += balanceIncrement;
                if (balance * balanceIncrement > MaxImbalance) {
                    _ASSERT(balance * balanceIncrement == MaxImbalance + 1);  // We can get at most one past the limit
                    //
                    // We're too far ahead.  Close our throttle.
                    //
                    PreventEventWaitersFromProceeding(&throttle[firstOrSecond]);
                } else if (balance * -1 * balanceIncrement == MaxImbalance) {
                    //
                    // We just pushed it back into balance (barely) for the other guy.  Allow him to
                    // proceed.
                    //
                    AllowEventWaitersToProceed(&throttle[1-firstOrSecond]);
                }
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

PairedReadSupplierFromQueue::~PairedReadSupplierFromQueue()
{}


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
