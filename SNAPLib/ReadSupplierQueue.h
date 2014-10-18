/*++

Module Name:

    ReadSupplierQueue.h

Abstract:

    Headers for parallel queue of reads

Authors:

    Bill Bolosky, November, 2012

Environment:

    User mode service.

Revision History:


--*/

#pragma once
#include "Read.h"
#include "Compat.h"
#include "VariableSizeVector.h"
#include "VariableSizeMap.h"

using std::pair;

class ReadSupplierFromQueue;
class PairedReadSupplierFromQueue;

typedef VariableSizeVector<DataBatch> BatchVector;

struct ReadQueueElement {
    ReadQueueElement()
        : next(NULL), prev(NULL)
    {
        reads = (Read*) BigAlloc(MaxReadsPerElement * sizeof(Read));
    }

    ~ReadQueueElement()
    {
        BigDealloc(reads);
        reads = NULL;
    }

    // note this should be about read buffer size for input reads
#ifdef LONG_READS
    static const int    MaxReadsPerElement = 400; 
#else
    static const int    MaxReadsPerElement = 5000; 
#endif
    ReadQueueElement    *next;
    ReadQueueElement    *prev;
    int                 totalReads;
    Read*               reads;
    BatchVector         batches;

    void addToTail(ReadQueueElement *queueHead) {
        next = queueHead;
        prev = queueHead->prev;
        prev->next = this;
        next->prev = this;
    }

    void removeFromQueue() {
        prev->next = next;
        next->prev = prev;
        prev = next = NULL;
    }
};
    
class ReadSupplierQueue: public ReadSupplierGenerator, public PairedReadSupplierGenerator {
public:
    //
    // This queue can handle several different kinds of inputs and outputs.  It will do either single
    // ended or paired reads.  In both cases, it can accept multiple independent readers (typically
    // one per (pair of) input file(s).  For paired reads that come from pairs of input files (think
    // FASTQ) it will run them independently and then combine the results as they're extracted.  For
    // paired reads that come from single files (SAM/BAM/CRAM, etc.) it still uses two queues internally,
    // but they're both written by a single PairedReadReader.
    //

    //
    // The version for single ended reads.  This is useful for formats that can't be divided by the
    // RangeSplitter, like BAM (though that's theoretically possible, so maybe..)  It takes a set 
    // of readers (presumably for different files), each of which runs independently and in parallel.
    // 
    ReadSupplierQueue(ReadReader *i_reader);

    //
    // The version for paired reads for which each end comes from a different Reader (and presumably
    // file, think FASTQ).  This is mostly useful for cases where the RangeSplitter can't handle
    // the files, probably because they FASTQ files with unequal length reads).
    //
    ReadSupplierQueue(ReadReader *i_firstHalfReader, ReadReader *i_secondHalfReader);

    //
    // The version for paired reads that come from a single file but for which RangeSplitter won't
    // work (BAM, CRAM, compressed FASTQ, maybe SRA).
    //
    ReadSupplierQueue(PairedReadReader *pairedReader);
    
    virtual ~ReadSupplierQueue();

    bool startReaders();
    void waitUntilFinished();
    ReadSupplier *generateNewReadSupplier();
    PairedReadSupplier *generateNewPairedReadSupplier();
    ReaderContext* getContext();

    ReadQueueElement *getElement();     // Called from the supplier threads
    bool getElements(ReadQueueElement **element1, ReadQueueElement **element2);   // Called from supplier threads
    void doneWithElement(ReadQueueElement *element);
    void supplierFinished();

    void holdBatch(DataBatch batch);
    bool releaseBatch(DataBatch batch);

    static int BufferCount(int numThreads)
    { return (__max(numThreads,2) + 1) * BatchesPerElement; }

private:

    static const int BatchesPerElement = 4;

    void commonInit();

    ReadReader          *singleReader[2];   // Only [0] is filled in for single ended reads
    PairedReadReader    *pairedReader;      // This is filled in iff there are no single readers

    ReadQueueElement    readyQueue[2];      // Queue [1] is used only when there are two single end readers

    BatchTracker        tracker;            // track batches used in queues, use refcount per element (not per read)

    EventObject         throttle[2];        // Two throttles, one for each of the readers.  At least one must be open at all times.
    int balance;                            // The size of readyQueue[0] - the size of readyQueue[1].  This is used to throttle.
    static const int MaxImbalance = 5;      // Engage the throttle when |balance| > MaxImbalance

    volatile unsigned   elementSize;        // reads per element, used to ensure paired single readers use same size that is ~ buffer size
 
    int                 nReadersRunning;
    int                 nSuppliersRunning;
    volatile bool       allReadsQueued;

    ReadQueueElement* getEmptyElement(); // must hold the lock to call this

    bool areAnyReadsReady(); // must hold the lock to call this.

    //
    // Empty buffers waiting for the readers.
    //
    ReadQueueElement    emptyQueue[1];
  
    //
    // Just one lock for all of the shared objects (the queues and Waiter objects, and counts of
    // readers and suppliers running, as well as allReadsQueued).
    //
    ExclusiveLock       lock;
    EventObject         readsReady;
    EventObject         emptyBuffersAvailable;

    EventObject         allReadsConsumed;

    struct ReaderThreadParams {
        ReadSupplierQueue       *queue;
        bool                     isSecondReader;
    };

    static void ReaderThreadMain(void *);
    void ReaderThread(ReaderThreadParams *params);
};

//
// A read supplier that takes its data from a ReadSupplierQueue.
//
class ReadSupplierFromQueue: public ReadSupplier {
public:
    ReadSupplierFromQueue(ReadSupplierQueue *i_queue);
    ~ReadSupplierFromQueue() {}

    Read *getNextRead();
    
    virtual void holdBatch(DataBatch batch)
    { queue->holdBatch(batch); }

    virtual bool releaseBatch(DataBatch batch)
    { return queue->releaseBatch(batch); }

private:
    bool                done;
    ReadSupplierQueue   *queue;
    bool                outOfReads;
    ReadQueueElement    *currentElement;
    int                 nextReadIndex;          
};

class PairedReadSupplierFromQueue: public PairedReadSupplier {
public:
    PairedReadSupplierFromQueue(ReadSupplierQueue *i_queue, bool i_twoFiles);
    ~PairedReadSupplierFromQueue();

    bool getNextReadPair(Read **read0, Read **read1);

    virtual void holdBatch(DataBatch batch)
    { queue->holdBatch(batch); }

    virtual bool releaseBatch(DataBatch batch)
    { return queue->releaseBatch(batch); }

private:
    ReadSupplierQueue   *queue;
    bool                done;
    bool                twoFiles;
    ReadQueueElement    *currentElement;
    ReadQueueElement    *currentSecondElement;
    int                 nextReadIndex;          
};
