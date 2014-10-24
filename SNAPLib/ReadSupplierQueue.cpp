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
#include "exit.h"
#include "SAM.h"

//#define PAIR_MATCH_DEBUG

 ReadSupplierQueue::ReadSupplierQueue(ReadReader *reader)
     : tracker(64)
{
    commonInit();

    singleReader[0] = reader;
}

ReadSupplierQueue::ReadSupplierQueue(ReadReader *firstHalfReader, ReadReader *secondHalfReader)
     : tracker(64)
{
    commonInit();

    singleReader[0] = firstHalfReader;
    singleReader[1] = secondHalfReader;
}

ReadSupplierQueue::ReadSupplierQueue(PairedReadReader *i_pairedReader)
     : tracker(128)
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

    balance = 0;

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
    elementSize = ReadQueueElement::MaxReadsPerElement;
}

ReadSupplierQueue::~ReadSupplierQueue()
{
    delete singleReader[0];
    delete singleReader[1];
    delete pairedReader;

    DestroyEventObject(&throttle[0]);
    DestroyEventObject(&throttle[1]);
    DestroyExclusiveLock(&lock);
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
    // Add more queue elements for this supplier.
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
    const int addElements = (singleReader[1] == NULL) ? 2 : 4 + MaxImbalance;
    ReadQueueElement * newElements[MaxImbalance + 4];
    for (int i = 0 ; i < addElements; i++) {
        newElements[i] = new ReadQueueElement;
    }

    AcquireExclusiveLock(&lock);
    nSuppliersRunning++;
    //
    // Add two more queue elements (4+MaxImbalance for paired-end, double file).
    //
    for (int i = 0; i < addElements; i++) {
        ReadQueueElement *element = newElements[i];
        element->addToTail(emptyQueue);
    }

    AllowEventWaitersToProceed(&emptyBuffersAvailable);
    ReleaseExclusiveLock(&lock);
   
    return new PairedReadSupplierFromQueue(this, singleReader[1] != NULL);
}

    ReaderContext*
ReadSupplierQueue::getContext()
{
    return singleReader[0] != NULL ? singleReader[0]->getContext() : pairedReader->getContext();
}

    ReadQueueElement *
ReadSupplierQueue::getElement()
{
    _ASSERT(singleReader[1] == NULL);   // i.e., we're doing file (but possibly single or paired end) reads
    //WriteErrorMessage("Thread %u: getElement wait acquire lock\n", GetThreadId());
    AcquireExclusiveLock(&lock);
    //WriteErrorMessage("Thread %u: getElement acquired lock\n", GetThreadId());
    while (!areAnyReadsReady()) {
        //WriteErrorMessage("Thread %u: getElement loop releasing lock\n", GetThreadId());
        ReleaseExclusiveLock(&lock);
        //WriteErrorMessage("Thread %u: getElement loop released lock\n", GetThreadId());
        if (allReadsQueued) {
            //
            // Everything's queued and the queue is empty.  No more work.
            //
            //WriteErrorMessage("Thread %u: getElement loop exit allReadsQueued\n", GetThreadId());
            return NULL;
        }
        //WriteErrorMessage("Thread %u: getElement loop wait readsReady\n", GetThreadId());
        WaitForEvent(&readsReady);
        //WriteErrorMessage("Thread %u: getElement loop wait acquire lock\n", GetThreadId());
        AcquireExclusiveLock(&lock);
        //WriteErrorMessage("Thread %u: getElement loop acquired lock\n", GetThreadId());
    }

    ReadQueueElement *element = readyQueue[0].next;
    _ASSERT(element != &readyQueue[0]);
    element->removeFromQueue();

    if (!areAnyReadsReady() && !allReadsQueued) {
        //WriteErrorMessage("Thread %u: getElement block readsReady\n", GetThreadId());
        PreventEventWaitersFromProceeding(&readsReady);
    }
    ReleaseExclusiveLock(&lock);
    //WriteErrorMessage("Thread %u: getElement released lock\n", GetThreadId());

    return element;
}

        bool 
ReadSupplierQueue::getElements(ReadQueueElement **element1, ReadQueueElement **element2)
{
   _ASSERT(singleReader[1] != NULL);   // i.e., we're doing paired file reads

   //WriteErrorMessage("Thread %u: getElements wait acquire lock\n", GetThreadId());
    AcquireExclusiveLock(&lock);
    //WriteErrorMessage("Thread %u: getElements acquired lock\n", GetThreadId());
    while (!areAnyReadsReady()) {
        //WriteErrorMessage("Thread %u: getElements loop releasing lock\n", GetThreadId());
        ReleaseExclusiveLock(&lock);
        //WriteErrorMessage("Thread %u: getElements loop released lock\n", GetThreadId());
        if (allReadsQueued) {
            //
            // Everything's queued and the queue is empty.  No more work.
            //
            //WriteErrorMessage("Thread %u: getElement loop exit allReadsQueued\n", GetThreadId());
            return NULL;
        }
        //WriteErrorMessage("Thread %u: getElements loop wait readsReady\n", GetThreadId());
        WaitForEvent(&readsReady);
        //WriteErrorMessage("Thread %u: getElements loop wait acquire lock\n", GetThreadId());
        AcquireExclusiveLock(&lock);
        //WriteErrorMessage("Thread %u: getElements loop acquired lock\n", GetThreadId());
    }

    ReadQueueElement* copyOut = NULL; // for adjusting sizes
    *element1 = readyQueue[0].next;
    *element2 = readyQueue[1].next;
    if ((*element1)->totalReads != (*element2)->totalReads) {
        copyOut = getEmptyElement(); // might release/acquire lock so state could change!
        *element1 = readyQueue[0].next;
        *element2 = readyQueue[1].next;
        if ((*element1)->totalReads == (*element2)->totalReads) {
            doneWithElement(copyOut);
        }
    }

    if ((*element1)->totalReads == (*element2)->totalReads) {
        (*element1)->removeFromQueue();
        (*element2)->removeFromQueue();
    } else {
        //fprintf(stderr,"getElements different sizes %d %d\n", (*element1)->totalReads, (*element2)->totalReads);
        // need to balance out reads between the two
        // make a copy of the min# of reads from larger element
        // shrink the larger element and leave it there
        ReadQueueElement* elements[2] = {*element1, *element2};
        int sizes[2] = {elements[0]->totalReads, elements[1]->totalReads};
        int largerOne = elements[1]->totalReads > elements[0]->totalReads;
        int minReads = elements[1-largerOne]->totalReads;
        memcpy(copyOut->reads, elements[largerOne]->reads, minReads * sizeof(Read));
        _ASSERT(elements[0]->totalReads == sizes[0] && elements[1]->totalReads == sizes[1] && elements[largerOne]->totalReads > elements[1-largerOne]->totalReads);
        copyOut->totalReads = minReads;
        _ASSERT(elements[0]->totalReads == sizes[0] && elements[1]->totalReads == sizes[1] && elements[largerOne]->totalReads > elements[1-largerOne]->totalReads);
        memmove(elements[largerOne]->reads, &elements[largerOne]->reads[minReads],
            (elements[largerOne]->totalReads - minReads) * sizeof(Read));
        elements[largerOne]->totalReads -= minReads;
        copyOut->batches.append(&elements[largerOne]->batches);
        for (BatchVector::iterator i = copyOut->batches.begin(); i != copyOut->batches.end(); i++) {
            holdBatch(*i);
        }
        if (largerOne == 0) {
            *element1 = copyOut;
            (*element2)->removeFromQueue();
        } else {
            (*element1)->removeFromQueue();
            *element2 = copyOut;
        }
        //WriteErrorMessage("Thread %u: balanced sizes %d %d\n", GetThreadId(), sizes[0], sizes[1]);
    }
    //fprintf(stderr,"getElements %x/%x with %d/%d reads\n", (int) (*element1), (int) (*element2), (*element1)->totalReads, (*element2)->totalReads);

    if (!areAnyReadsReady() && !allReadsQueued) {
        //WriteErrorMessage("Thread %u: getElements block readsReady\n", GetThreadId());
        PreventEventWaitersFromProceeding(&readsReady);
    }
 
    //WriteErrorMessage("Thread %u: getElements releasing lock\n", GetThreadId());
    ReleaseExclusiveLock(&lock);
    //WriteErrorMessage("Thread %u: getElements released lock", GetThreadId());
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
    //WriteErrorMessage("Thread %u: doneWithElement wait acquire lock\n", GetThreadId());
    AcquireExclusiveLock(&lock);
    //WriteErrorMessage("Thread %u: doneWithElement acquired lock\n", GetThreadId());
    _ASSERT(element->totalReads > 0);
    VariableSizeVector<DataBatch> batches = element->batches;
    element->batches.clear();
    element->addToTail(emptyQueue);
    AllowEventWaitersToProceed(&emptyBuffersAvailable);
    //WriteErrorMessage("Thread %u: doneWithElement releasing lock\n", GetThreadId());
    ReleaseExclusiveLock(&lock);
    //WriteErrorMessage("Thread %u: doneWithElement released lock\n", GetThreadId());
    for (VariableSizeVector<DataBatch>::iterator b = batches.begin(); b != batches.end(); b++) {
        releaseBatch(*b);
    }
}

    void 
ReadSupplierQueue::supplierFinished()
{
    //WriteErrorMessage("Thread %u: supplierFinished wait acquire lock\n", GetThreadId());
    AcquireExclusiveLock(&lock);
    //WriteErrorMessage("Thread %u: supplierFinished acquired lock\n", GetThreadId());
    _ASSERT(allReadsQueued);
    _ASSERT(nSuppliersRunning > 0);
    nSuppliersRunning--;
    if (0 == nSuppliersRunning) {
        AllowEventWaitersToProceed(&allReadsConsumed);
    }
    //WriteErrorMessage("Thread %u: supplierFinished releasing lock\n", GetThreadId());
    ReleaseExclusiveLock(&lock);
    //WriteErrorMessage("Thread %u: supplierFinished released lock\n", GetThreadId());
}
    
    void
ReadSupplierQueue::holdBatch(
    DataBatch batch)
{
    if (pairedReader != NULL) {
        pairedReader->holdBatch(batch);
    } else if (singleReader[1] == NULL) {
        singleReader[0]->holdBatch(batch);
    } else {
        singleReader[batch.fileID % 2]->holdBatch(DataBatch(batch.batchID, batch.fileID / 2));
    }
}

    bool
ReadSupplierQueue::releaseBatch(
    DataBatch batch)
{
    if (pairedReader != NULL) {
        return pairedReader->releaseBatch(batch);
    } else if (singleReader[1] == NULL) {
        return singleReader[0]->releaseBatch(batch);
    } else {
        return singleReader[batch.fileID % 2]->releaseBatch(DataBatch(batch.batchID, batch.fileID / 2));
    }
}

    void
ReadSupplierQueue::ReaderThreadMain(void *param)
{
    ReaderThreadParams *params = (ReaderThreadParams *)param;
    params->queue->ReaderThread(params);
    delete params;
}

    ReadQueueElement*
ReadSupplierQueue::getEmptyElement()
{
    while (emptyQueue->next == emptyQueue) {
        // Wait for a buffer.
        //WriteErrorMessage("Thread %u: getEmptyElement releasing lock\n", GetThreadId());
        ReleaseExclusiveLock(&lock);
        //WriteErrorMessage("Thread %u: getEmptyElement released lock\n", GetThreadId());
        WaitForEvent(&emptyBuffersAvailable);
        //WriteErrorMessage("Thread %u: getEmptyElement acquiring lock\n", GetThreadId());
        AcquireExclusiveLock(&lock);
        //WriteErrorMessage("Thread %u: getEmptyElement acquired lock\n", GetThreadId());
    }

    ReadQueueElement *element = emptyQueue->next;
    element->removeFromQueue();
    if (emptyQueue->next == emptyQueue) {
        PreventEventWaitersFromProceeding(&emptyBuffersAvailable);
    }

    return element;
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
    int firstOrSecond = params->isSecondReader ? 1 : 0;
    bool isSingleReader = (NULL == singleReader[1]);

    _int64 balanceTime = 0;
    _int64 bufferWaitTime = 0;
    _int64 processingTime = 0;
    _int64 startTime = timeInNanos();

    // may pass reads forward from element loops to maintain batching or element size
    Read firstReadForNextElement[2];
    bool hasFirstReadForNextElement = false;

    while (!done) {
        if ((!isSingleReader) && balance * balanceIncrement > MaxImbalance) {
            //
            // We're over full.  Wait to get back in balance.
            //
            ReleaseExclusiveLock(&lock);

            _int64 now = timeInNanos();
            processingTime += now - startTime;
            startTime = now;

            WaitForEvent(&throttle[firstOrSecond]);

            now = timeInNanos();
            balanceTime += now - startTime;
            startTime = now;

            AcquireExclusiveLock(&lock);
            _ASSERT(balance * balanceIncrement <= MaxImbalance);
        }

        // pull an empty element from the queue
        _int64 now = timeInNanos();
        processingTime += now - startTime;
        startTime = now;
        ReadQueueElement* element = getEmptyElement();
        now = timeInNanos();
        bufferWaitTime += now - startTime;
        startTime = now;

        //
        // Now fill in the reads from the reader into the element until it's
        // full or the reader finishes or it exceeds batch count
        //
        ReleaseExclusiveLock(&lock);
        element->totalReads = 0;
        for (; element->totalReads <= (int) elementSize - increment; element->totalReads += increment) {
            
            if (NULL != reader) {
                Read* read = &element->reads[element->totalReads];
                if (hasFirstReadForNextElement) {
                    _ASSERT(element->totalReads == 0);
                    *read = firstReadForNextElement[0];
                    // already called holdBatch when firstReadForNextElement set
                    element->batches.push_back(read->getBatch());
                    hasFirstReadForNextElement = false;
                } else {
                    done = ! reader->getNextRead(read);
                    if (done) {
                        break;
                    }
                    if (! isSingleReader) {
                        read->setBatch(DataBatch(read->getBatch().batchID, read->getBatch().fileID * 2 + firstOrSecond));
                    }
                    bool newBatch = element->totalReads == 0 ||
                        (read->getBatch() != read[-1].getBatch() && element->batches.search(read->getBatch()) == element->batches.end());
                    if (element->batches.size() + newBatch <= BatchesPerElement) {
                        // won't exceed limit for this element
                        if (newBatch && element->batches.add(read->getBatch())) {
                            holdBatch(read->getBatch());
                        }
                    } else {
                        // too many batches, hold for next queue element
                        firstReadForNextElement[0] = *read;
                        hasFirstReadForNextElement = true;
                        holdBatch(read->getBatch());
                        break;
                    }
                }
            } else if (NULL != pairedReader) {
                Read* read = &element->reads[element->totalReads];
                if (hasFirstReadForNextElement) {
                    _ASSERT(element->totalReads == 0);
                    read[0] = firstReadForNextElement[0];
                    read[1] = firstReadForNextElement[1];
                    // already called holdBatch when firstReadForNextElement set
                    element->batches.push_back(read[0].getBatch());
                    if (read[1].getBatch() != read[0].getBatch()) {
                        element->batches.push_back(read[1].getBatch());
                    }
                    hasFirstReadForNextElement = false;
                } else {
                    done = !pairedReader->getNextReadPair(&read[0], &read[1]);
                    if (done) {
                        break;
                    }
                    DataBatch b[2] = {read[0].getBatch(), read[1].getBatch()};
                    bool newBatch[2] =
                        {(element->totalReads == 0 || read[-2].getBatch() != b[0]) &&
                            element->batches.search(b[0]) == element->batches.end(),
                        b[0] != b[1] && (element->totalReads == 0 || read[-1].getBatch() != b[1]) &&
                            element->batches.search(b[1]) == element->batches.end()};
                    if (element->batches.size() + newBatch[0] + newBatch[1] <= BatchesPerElement) {
                        if (newBatch[0] && element->batches.add(b[0])) {
                            holdBatch(b[0]);
                        }
                        if (newBatch[1] && element->batches.add(b[1])) {
                            holdBatch(b[1]);
                        }
                    } else {
                        firstReadForNextElement[0] = read[0];
                        firstReadForNextElement[1] = read[1];
                        holdBatch(b[0]);
                        if (b[1] != b[0]) {
                            holdBatch(b[1]);
                        }
                        hasFirstReadForNextElement = true;
                        break;
                    }
                    //WriteErrorMessage("ReadSupplierQueue::ReaderThread element %x batches %d:%d, %d:%d\n", (int) element, b[0].fileID, b[0].batchID, b[1].fileID, b[1].batchID);
                }
           }
        }

        //WriteErrorMessage("ReadSupplierQueue element[%d] %x with %d reads %d batches\n", firstOrSecond, (int) element, element->totalReads, element->batches.size());
        
        AcquireExclusiveLock(&lock);
        
        // do this before AllowEventWaitersToProceed to avoid race condition
        if (done && 1 == nReadersRunning) {
            //WriteErrorMessage("Thread %u: set allReadsQueued (%d) in ReaderThread...\n", GetThreadId(), element->totalReads);
            allReadsQueued = true;
            AllowEventWaitersToProceed(&readsReady);    // Even if we have nothing to queue, allow the consumers to wake up so they can exit
        }

        if (element->totalReads > 0) {
            element->addToTail(&readyQueue[firstOrSecond]);
            if (isSingleReader || &readyQueue[1-firstOrSecond] != readyQueue[1-firstOrSecond].next) {
                //
                // Signal that an element is ready.
                //
                //WriteErrorMessage("Thread %u: signal readsReady in ReaderThread...\n", GetThreadId());
                AllowEventWaitersToProceed(&readsReady);
            }

            if (!isSingleReader) {
                //WriteErrorMessage("Thread %u: balance %d %+d = %d...\n", GetThreadId(), balance, balanceIncrement, balance + balanceIncrement);
                balance += balanceIncrement;
                if (balance * balanceIncrement > MaxImbalance) {
                    _ASSERT(balance * balanceIncrement == MaxImbalance + 1);  // We can get at most one past the limit
                    //
                    // We're too far ahead.  Close our throttle.
                    //
		    //WriteErrorMessage("Thread %u: close throttle %d in ReaderThread...\n", GetThreadId(), firstOrSecond);
                    PreventEventWaitersFromProceeding(&throttle[firstOrSecond]);
                } else if (balance * -1 * balanceIncrement == MaxImbalance) {
                    //
                    // We just pushed it back into balance (barely) for the other guy.  Allow him to
                    // proceed.
                    //
		    //WriteErrorMessage("Thread %u: release throttle %d in ReaderThread...\n", GetThreadId(), 1-firstOrSecond);
                    AllowEventWaitersToProceed(&throttle[1-firstOrSecond]);
                }
            }
        }
    } // While ! done

    processingTime += timeInNanos() - startTime;

    //WriteErrorMessage("ReadSupplier: %llds processing, %llds waiting for balance, %llds waiting for buffer\n", processingTime / 1000000000, balanceTime / 1000000000, bufferWaitTime / 1000000000);

    _ASSERT(nReadersRunning > 0);
    nReadersRunning--;
    ReleaseExclusiveLock(&lock);
}

ReadSupplierFromQueue::ReadSupplierFromQueue(
    ReadSupplierQueue *i_queue)
    :
    queue(i_queue),
    outOfReads(false),
    currentElement(NULL),
    nextReadIndex(0),
    done(false)
{
}

    Read *
ReadSupplierFromQueue::getNextRead()
{
    if (done) {
        return NULL;
    }

    ReadQueueElement* doneElement = NULL;
    if (NULL != currentElement && nextReadIndex >= currentElement->totalReads) {
        doneElement = currentElement;
        currentElement = NULL;
    }

    if (NULL == currentElement) {
        currentElement = queue->getElement();
        if (doneElement != NULL) {
            queue->doneWithElement(doneElement);
        }
        if (NULL == currentElement) {
            done = true;
            queue->supplierFinished();
            return NULL;
        }
        nextReadIndex = 0;
    }

    return &currentElement->reads[nextReadIndex++]; // Note the post increment.
}

PairedReadSupplierFromQueue::PairedReadSupplierFromQueue(ReadSupplierQueue *i_queue, bool i_twoFiles) :
    queue(i_queue), twoFiles(i_twoFiles), done(false), 
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
        //fprintf(stderr,"PairedReadSupplierFromQueue finished element %x with %d reads %d batches:", (int) currentElement, currentElement->totalReads, currentElement->batches.size()); for (BatchVector::iterator i = currentElement->batches.begin(); i != currentElement->batches.end(); i++) { fprintf(stderr," %d:%d", i->fileID, i->batchID); } fprintf(stderr,"\n");
        queue->doneWithElement(currentElement);
        currentElement = NULL;
        if (twoFiles) {
            //fprintf(stderr,"PairedReadSupplierFromQueue finished 2nd element %x with %d reads %d batches:", (int) currentSecondElement, currentSecondElement->totalReads, currentSecondElement->batches.size()); for (BatchVector::iterator i = currentSecondElement->batches.begin(); i != currentSecondElement->batches.end(); i++) { printf(" %d:%d", i->fileID, i->batchID); } printf("\n");
            queue->doneWithElement(currentSecondElement);
            currentSecondElement = NULL;
        }
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
		nextReadIndex = 0;
    }
    if (twoFiles) {
        // Assert that both elements match.
        _ASSERT(currentSecondElement->totalReads == currentElement->totalReads);
#ifdef PAIR_MATCH_DEBUG
        for (int i = 0; i < currentElement->totalReads; i++) {
            Read::checkIdMatch(&currentElement->reads[i], &currentSecondElement->reads[i]);
        }
#endif
    } else {
        //
        // Assert that there are an even number of reads (since they're in pairs)
        //
        _ASSERT(currentElement->totalReads % 2 == 0);
#ifdef PAIR_MATCH_DEBUG
        for (int i = 0; i < currentElement->totalReads; i += 2) {
            Read::checkIdMatch(&currentElement->reads[i], &currentElement->reads[i+1]);
        }
#endif
    }

    if (twoFiles) {
        *read0 = &currentElement->reads[nextReadIndex];
        *read1 = &currentSecondElement->reads[nextReadIndex];
#ifdef PAIR_MATCH_DEBUG
		Read::checkIdMatch(*read0, *read1);
#endif
        nextReadIndex++;
    } else {
        *read0 = &currentElement->reads[nextReadIndex];
        *read1 = &currentElement->reads[nextReadIndex+1];
#ifdef PAIR_MATCH_DEBUG
		Read::checkIdMatch(*read0, *read1);
#endif
        nextReadIndex += 2;
    }

    return true;
}
    
