/*++

Module Name:

    SortedDataWriter.cpp

Abstract:

    File writer that sorts records using a temporary file.

Environment:

    User mode service.

    Not thread safe.

--*/

#include "stdafx.h"
#include "BigAlloc.h"
#include "Compat.h"
#include "Util.h"
#include "DataWriter.h"
#include "BufferedAsync.h"
#include "VariableSizeVector.h"
#include "FileFormat.h"
#include "PriorityQueue.h"
#include "exit.h"
#include "Bam.h"
#include "Error.h"

//#define VALIDATE_SORT 1

using std::max;

#pragma pack(push, 4)
struct SortEntry
{
    SortEntry() : offset(0), length(0) {}
    SortEntry(size_t i_offset, GenomeDistance i_length,ContigAndPos i_contigAndPos)
        : offset(i_offset), length(i_length), contigAndPos(i_contigAndPos) {}
    size_t                      offset; // offset in file
    _int64                      length; // number of bytes
    //GenomeLocation              location; // location in genome
    ContigAndPos                contigAndPos;
    static bool comparator(const SortEntry& e1, const SortEntry& e2)
    {
        return e1.contigAndPos < e2.contigAndPos;
    }
};
#pragma pack(pop)

typedef VariableSizeVector<SortEntry,150,true> SortVector;

struct SortBlock
{
#ifdef VALIDATE_SORT
    SortBlock() : start(0), bytes(0), location(0), length(0), reader(NULL), minLocation(0), maxLocation(0) {}
#else
    SortBlock() : start(0), bytes(0), length(0), reader(NULL), dataReaderIsBuffer(false), data(NULL) {}
    SortBlock(DataReader* bufferDataReader) : start(0), bytes(0), length(0), reader(bufferDataReader), dataReaderIsBuffer(bufferDataReader != NULL), data(NULL) {}
#endif
	SortBlock(const SortBlock& other) { *this = other; }
    void operator=(const SortBlock& other);

    size_t      start;
    size_t      bytes;
#ifdef VALIDATE_SORT
	GenomeLocation	minLocation, maxLocation;
#endif
    // for mergesort phase
    DataReader* reader;
    ContigAndPos    contigAndPos;   // the sort key
    char*       data; // read data in read buffer
    GenomeDistance    length; // length in bytes
    bool dataReaderIsBuffer;
};

    void
SortBlock::operator=(
    const SortBlock& other)
{
    start = other.start;
    bytes = other.bytes;
    length = other.length;
    reader = other.reader;
    dataReaderIsBuffer = other.dataReaderIsBuffer;
#ifdef VALIDATE_SORT
	minLocation = other.minLocation;
	maxLocation = other.maxLocation;
#endif
}

typedef VariableSizeVector<SortBlock> SortBlockVector;
    
class SortedDataFilterSupplier;

class SortedDataFilter : public DataWriter::Filter
{
public:
    SortedDataFilter(SortedDataFilterSupplier* i_parent)
        : Filter(DataWriter::CopyFilter), parent(i_parent), locations(10000000), seenLastBatch(false)
    {}

    virtual ~SortedDataFilter() 
    {
        _ASSERT(seenLastBatch);
    }

    virtual void onAdvance(DataWriter* writer, size_t batchOffset, char* data, GenomeDistance bytes, GenomeLocation location);

    virtual size_t onNextBatch(DataWriter* writer, size_t offset, size_t bytes, bool lastBatch = false, bool* needMoreBuffer = NULL, size_t* fromBufferUsed = NULL);

private:
    SortedDataFilterSupplier*   parent;
    SortVector                  locations;
    bool                        seenLastBatch;
};

class SortedDataFilterSupplier : public DataWriter::FilterSupplier
{
public:

    SortedDataFilterSupplier(
        const FileFormat* i_fileFormat,
        const Genome* i_genome,
        const char* i_tempFileName,
        const char* i_sortedFileName,
        DataWriter::FilterSupplier* i_sortedFilterSupplier,
        size_t i_bufferSize,
        size_t i_bufferSpace,
        bool i_emitInternalScore,
        char *i_internalScoreTag,
        FileEncoder* i_encoder,
        int i_numThreads)
        :
        format(i_fileFormat),
        genome(i_genome),
        FilterSupplier(DataWriter::CopyFilter),
        encoder(i_encoder),
        tempFileName(i_tempFileName),
        sortedFileName(i_sortedFileName),
        sortedFilterSupplier(i_sortedFilterSupplier),
        bufferSize(i_bufferSize),
        bufferSpace(i_bufferSpace),
        blocks(),
        emitInternalScore(i_emitInternalScore),
        totalReadsSorted(0),
        numThreads(i_numThreads)
    {
        if (emitInternalScore) {
            if (strlen(i_internalScoreTag) != 2) {  // This should never happen, since the command line parser should catch it first.  Still, since we're about to strcpy into a fixed-length buffer, safety first.
                WriteErrorMessage("SortedDataFilterSupplier: improper internal score tag\n");
                soft_exit(1);
            }
            strcpy(internalScoreTag, i_internalScoreTag);
        } else {
            internalScoreTag[0] = '\0';
        }
        InitializeExclusiveLock(&lock);
    }

    virtual ~SortedDataFilterSupplier()
    {
        DestroyExclusiveLock(&lock);
        delete encoder;
    }

    virtual DataWriter::Filter* getFilter();

    virtual void onClosing(DataWriterSupplier* supplier) {}
    virtual void onClosed(DataWriterSupplier* supplier);

    void setHeaderSize(size_t bytes)
    { headerSize = bytes; }

    inline const Genome* getGenome() {
        return genome;
    }

#ifndef VALIDATE_SORT
	void addBlock(size_t start, size_t bytes, DataReader *reader = NULL);
#else
    void addBlock(size_t start, size_t bytes, GenomeLocationOrderedByOriginalContigs minLocation, GenomeLocationOrderedByOriginalContigs maxLocation);
#endif

private:
    bool mergeSort();
    void mergeSortThread(SortBlockVector* blocksForThisThread, DataWriter *writer);
    static void MergeSortThreadMain(void* threadParameter);

    void mergeSortNode(DataReader* readers, int nReaders, DataWriter* writer);  // Merges the data from the readers, and writes it to the writers.

    const Genome*                   genome;
    const FileFormat*               format;
    const char*                     tempFileName;
    const char*                     sortedFileName;
    DataWriter::FilterSupplier*     sortedFilterSupplier;
    FileEncoder*                    encoder;
    size_t                          headerSize;
    ExclusiveLock                   lock; // for adding blocks
    SortBlockVector                 blocks;
    size_t                          bufferSize;
    size_t                          bufferSpace;
    bool                            emitInternalScore;
    char                            internalScoreTag[3];
    _int64                          totalReadsSorted;
    int                             numThreads;

	friend class SortedDataFilter;
};

class ParallelQueue
{
public:
    ParallelQueue(int maxQueueDepth_);
    ~ParallelQueue();

    void enqueue(void* data);
    void* dequeue();

    //
    // It's born with one writer.  When all the writers are gone and all the elements
    // have been dequeued, dequeue will return NULL.
    //
    void addWriter();
    void releaseWriter();

private:
 
    int maxQueueDepth;

    struct QueueElement
    {
        void* data;
        QueueElement* next;
        QueueElement* prev;

        QueueElement()
        {
            data = NULL;
            next = prev = NULL;
        }

        void add(QueueElement* head) 
        {
            next = head;
            prev = head->prev;
            next->prev = this;
            prev->next = this;
        }

        void remove()
        {
            next->prev = prev;
            prev->next = next;

            next = prev = NULL;
        }

        //
        // These methods are only to be called on queue heads.
        //
        void initQueueHead()
        {
            next = prev = this;
        }

        bool isEmpty()
        {
            return next == this;
        }

        QueueElement* removeFirst()
        {
            _ASSERT(!isEmpty());

            QueueElement* element = next;
            element->remove();
            return element;
        }
    };

    QueueElement readyQueueHead;
    QueueElement freeQueueElementHead;
    QueueElement* queueElements;    // maxQueueDepth total elements

    EventObject queueNotFullEvent;
    EventObject queueHasDataEvent;

    int nWriters;

    ExclusiveLock queueLock;
};

ParallelQueue::ParallelQueue(int maxQueueDepth_)
{
    maxQueueDepth = maxQueueDepth_;
    queueElements = new QueueElement[maxQueueDepth];

    readyQueueHead.initQueueHead();
    freeQueueElementHead.initQueueHead();

    for (int i = 0; i < maxQueueDepth; i++) 
    {
        queueElements[i].add(&freeQueueElementHead);
    }

    InitializeExclusiveLock(&queueLock);
    CreateEventObject(&queueNotFullEvent);
    CreateEventObject(&queueHasDataEvent);

    AllowEventWaitersToProceed(&queueNotFullEvent);
    PreventEventWaitersFromProceeding(&queueHasDataEvent);

    nWriters = 1;
} // ParallelQueue::ParallelQueue

ParallelQueue::~ParallelQueue()
{
    DestroyExclusiveLock(&queueLock);
    DestroyEventObject(&queueNotFullEvent);
    DestroyEventObject(&queueHasDataEvent);
    delete[] queueElements;
}

    void
ParallelQueue::enqueue(void* data) 
{
        AcquireExclusiveLock(&queueLock);
        while (freeQueueElementHead.isEmpty()) 
        {
            ReleaseExclusiveLock(&queueLock);
            WaitForEvent(&queueNotFullEvent);
            AcquireExclusiveLock(&queueLock);
        }

        ParallelQueue::QueueElement* element = freeQueueElementHead.removeFirst();
        
        if (freeQueueElementHead.isEmpty()) 
        {
            PreventEventWaitersFromProceeding(&queueNotFullEvent);
        }

        element->data = data;
        element->add(&readyQueueHead);
        AllowEventWaitersToProceed(&queueHasDataEvent);

        ReleaseExclusiveLock(&queueLock);

} // ParallelQueue::enqueue


    void *
ParallelQueue::dequeue()
{
    AcquireExclusiveLock(&queueLock);
    while (readyQueueHead.isEmpty() && nWriters > 0)
    {
        ReleaseExclusiveLock(&queueLock);
        WaitForEvent(&queueHasDataEvent);
        AcquireExclusiveLock(&queueLock);
    }

    if (readyQueueHead.isEmpty()) 
    {
        ReleaseExclusiveLock(&queueLock);
        return NULL;
    }

    QueueElement* element = readyQueueHead.removeFirst();
    if (readyQueueHead.isEmpty()) 
    {
        PreventEventWaitersFromProceeding(&queueHasDataEvent);
    }
    
    void* data = element->data;
    element->add(&freeQueueElementHead);
    AllowEventWaitersToProceed(&queueNotFullEvent);

    ReleaseExclusiveLock(&queueLock);
    return data;
}

void
ParallelQueue::addWriter()
{
    AcquireExclusiveLock(&queueLock);
    _ASSERT(nWriters > 0);
    nWriters++;
    ReleaseExclusiveLock(&queueLock);
} // ParallelQueue::addWriter()

void
ParallelQueue::releaseWriter()
{
    AcquireExclusiveLock(&queueLock);
    _ASSERT(nWriters > 0);
    nWriters--;
    if (0 == nWriters)
    {
        AllowEventWaitersToProceed(&queueHasDataEvent);
    }
    ReleaseExclusiveLock(&queueLock);
} // ParallelQueue::releaseWriter()

//
// A fixed-capacity multi-thread queue where the input and output are implemented as DataWriters and DataReaders.
//
class DataQueueReader;
class DataQueueWriter;
class DataQueue 
{
public:
    DataQueue(int nBuffers_, size_t bufferSize_);
    ~DataQueue();

    DataQueueReader* getReader();
    DataQueueWriter* getWriter();

private:

    ParallelQueue* freeBufferQueue;
    ParallelQueue* readyBufferQueue;

    //
    // These are what're on the queues.
    //
    struct Buffer {
        char* buffer;
        size_t usedBufferSpace;
    };

    int nBuffers;
    size_t bufferSize;
    Buffer* buffers;
 
    ExclusiveLock queueLock;
    int nWriters;
    int nReaders;

    void releaseReader();
    void releaseWriter();

    friend class DataQueueWriter;
    friend class DataQueueReader;
};

    void
DataQueue::releaseReader()
{
    AcquireExclusiveLock(&queueLock);
    _ASSERT(nReaders > 0);
    nReaders--;

    if (nReaders == 0) 
    {
        freeBufferQueue->releaseWriter();
    }
    ReleaseExclusiveLock(&queueLock);
}

    void
DataQueue::releaseWriter()
{
    AcquireExclusiveLock(&queueLock);
    _ASSERT(nWriters > 0);
    nWriters--;

    if (nWriters == 0)
    {
        readyBufferQueue->releaseWriter();
    }
    ReleaseExclusiveLock(&queueLock);
}


DataQueue::DataQueue(int nBuffers_, size_t bufferSize_)
{
    _ASSERT(bufferSize > 0);
    nBuffers = nBuffers_;
    bufferSize = bufferSize_;
    nWriters = 0;
    nReaders = 0;

    freeBufferQueue = new ParallelQueue(nBuffers);
    readyBufferQueue = new ParallelQueue(nBuffers);

    buffers = new Buffer[nBuffers];
    for (int i = 0; i < nBuffers; i++)
    {
        buffers[i].buffer = (char*)BigAlloc(bufferSize);
        buffers[i].usedBufferSpace = 0; // Though this is meaningless for free buffers anyway.

        freeBufferQueue->enqueue(&buffers[i]);
    }

    InitializeExclusiveLock(&queueLock);
}

DataQueue::~DataQueue()
{
    _ASSERT(nReaders == 0 && nWriters == 0);

    DestroyExclusiveLock(&queueLock);
    delete freeBufferQueue;
    delete readyBufferQueue;

    for (int i = 0; i < nBuffers; i++) 
    {
        if (buffers[i].buffer != NULL) {
            BigDealloc(buffers[i].buffer);
            buffers[i].buffer = NULL;
        }
    }

    delete[] buffers;
}

class DataQueueWriter : public DataWriter 
{
public:
    DataQueueWriter(DataQueue *queue_) : DataWriter(NULL), currentBuffer(NULL)
    {
        queue = queue_;
        nextBatch();
    }

    // get remaining space in current buffer for writing
    virtual bool getBuffer(char** o_buffer, size_t* o_size)
    {
        if (currentBuffer == NULL || currentBuffer->usedBufferSpace >= queue->bufferSize) 
        {
            *o_buffer = NULL;
            *o_size = 0;
            return false;
        }

        *o_buffer = currentBuffer->buffer + currentBuffer->usedBufferSpace;
        *o_size = queue->bufferSize - currentBuffer->usedBufferSpace;

        return true;
    }

    // advance within current buffer, reducing available space
    // should be called on each read, with the location
    virtual void advance(_int64 bytes, GenomeLocation location = 0)
    {
        _ASSERT(currentBuffer != NULL && currentBuffer->usedBufferSpace + bytes <= queue->bufferSize);
        currentBuffer->usedBufferSpace += bytes;
    }

    // get complete data buffer in batch, relative==0 is current, relative==-1 is previous, etc.
    // if negative gets old data written, else waits for write to complete so you can write into it
    // o_offset gets physical offset (e.g. compressed), o_logical gets data offset (e.g. uncompressed)
    virtual bool getBatch(int relative, char** o_buffer, size_t* o_size = NULL, size_t* o_used = NULL, size_t* o_offset = NULL, size_t* o_logicalUsed = 0, size_t* o_logicalOffset = NULL)
    {
        WriteErrorMessage("DataQueueWriter: getBatch not implemented\n");
        soft_exit(1);
        return false;
    }

    // advance to next buffer
    virtual bool nextBatch(bool lastBatch = false)
    {
        if (currentBuffer != NULL) 
        {
            queue->readyBufferQueue->enqueue(currentBuffer);
            currentBuffer = NULL;
        }

        currentBuffer = (DataQueue::Buffer *)queue->freeBufferQueue->dequeue();
        currentBuffer->usedBufferSpace = 0;

        return true;
    }

    // this thread is complete
    virtual void close()
    {
        if (currentBuffer != NULL && currentBuffer->usedBufferSpace != 0) 
        {   
            //
            // Write out what we've got.
            //
            nextBatch(true);
        }
        queue->releaseWriter();
    }
private:

    DataQueue* queue;
    DataQueue::Buffer* currentBuffer;
};

class DataQueueReader : public DataReader
{
public:
    DataQueueReader(DataQueue* queue_) : queue(queue_), currentBuffer(NULL), readOffsetInCurrentBuffer(0)
    {
    }

    ~DataQueueReader()
    {
        if (NULL != queue)
        {
            queue->releaseReader();
        }
    }

    bool init(const char* fileName) { return true;  }
    char* readHeader(_int64* io_headerSize)
    {
        if (io_headerSize != NULL)
        {
            *io_headerSize = 0;
        }

        return NULL;
    }

    // seek to a particular range in the file
    virtual void reinit(_int64 startingOffset, _int64 amountOfFileToProcess)
    {
        WriteErrorMessage("DataQueueReader: reinit not supported\n");
        soft_exit(1);
    }

    // get all remaining data in current batch
    // return false if no more data in current batch
    // startBytes is data "owned" by this block in which reads may start
    // validBytes may also include overflow bytes to handle records spanning batches
    // if you advance() past startBytes, nextBatch() will start offset at that point
    virtual bool getData(char** o_buffer, _int64* o_validBytes, _int64* o_startBytes = NULL);

    // advance through data in current batch, reducing results from next getData call
    virtual void advance(_int64 bytes)
    {
        _ASSERT(readOffsetInCurrentBuffer + bytes <= currentBuffer->usedBufferSpace);
        readOffsetInCurrentBuffer += bytes;
    }

    // advance to next batch
    // by default automatically releases previous batch
    virtual void nextBatch();

    // whether current batch is last in file
    virtual bool isEOF()
    {
        WriteErrorMessage("DataQueueReader: isEOF not supported\n");
        soft_exit(1);
        return false;
    }

    // get current batch identifier
    virtual DataBatch getBatch()
    {
        return DataBatch(); // No real identifier here
    }

    // hold buffers associated with this batch for reuse, increments refcount
    // NOTE: this may be called from another thread,
    // so anything it touches must be thread-safe!
    virtual void holdBatch(DataBatch batch)
    {
        WriteErrorMessage("DataQueueReader: holdbatch not supported\n");
        soft_exit(1);
    }

    // release buffers associated with this batch for reuse
    // decrements refcount, returns true if last release
    // NOTE: this may be called from another thread,
    // so anything it touches must be thread-safe!
    virtual bool releaseBatch(DataBatch batch)
    {
        WriteErrorMessage("DataQueueReader: releaseBatch not supported\n");
        soft_exit(1);
        return true;
    }
    // get current offset into file
    virtual _int64 getFileOffset()
    {
        WriteErrorMessage("DataQueueReader: getFileOffset not supported\n");
        soft_exit(1);
        return -1;
    }

    // get pointer to extra data area for current batch
    // todo: allow this to grow dynamically while keeping stable pointers to previous data
    virtual void getExtra(char** o_extra, _int64* o_length)
    {
        WriteErrorMessage("DataQueueReader: getExtra not supported\n");
        soft_exit(1);
    }

    // get filename for debugging / error printing
    virtual const char* getFilename()
    {
        return "DataQueue";
    }
private:

    DataQueue* queue;
    DataQueue::Buffer* currentBuffer;
    size_t readOffsetInCurrentBuffer;
};


DataQueueReader*
DataQueue::getReader()
{
    AcquireExclusiveLock(&queueLock);
    nReaders++;
    ReleaseExclusiveLock(&queueLock);

    return new DataQueueReader(this);
}

DataQueueWriter*
DataQueue::getWriter()
{
    AcquireExclusiveLock(&queueLock);
    nWriters++;
    ReleaseExclusiveLock(&queueLock);

    return new DataQueueWriter(this);
}

    bool 
DataQueueReader::getData(char** o_buffer, _int64* o_validBytes, _int64* o_startBytes)
{
    _ASSERT(o_startBytes == NULL);  // this isn't used, so we don't fill it in

    if (currentBuffer == NULL || readOffsetInCurrentBuffer >= currentBuffer->usedBufferSpace) 
    {
        *o_buffer = NULL;
        *o_validBytes = 0;
        return false;
    }

    *o_buffer = currentBuffer->buffer + readOffsetInCurrentBuffer;
    *o_validBytes = currentBuffer->usedBufferSpace - readOffsetInCurrentBuffer;

    return true;
} // DataQueueReader::getData


    void
DataQueueReader::nextBatch()
{
    if (queue == NULL) 
    {
        _ASSERT(currentBuffer == NULL);
        return;
    }

    if (currentBuffer != NULL) 
    {
        queue->freeBufferQueue->enqueue(currentBuffer);
    }

    currentBuffer = (DataQueue::Buffer *)queue->readyBufferQueue->dequeue();

    if (currentBuffer == NULL) 
    {
        queue->releaseReader();
        queue = NULL;
    }
    
    readOffsetInCurrentBuffer = 0;
} // DataQueueReader::nextBatch()

class BufferDataReader : public DataReader
{
public:
    BufferDataReader(size_t dataSize_) : dataSize(dataSize_), readOffsetInBuffer(0)
    {
        buffer = (char*)BigAlloc(dataSize + 4096);   // Allow a little empty space at the end
        memset(buffer, 0, dataSize + 4096);         // Write it sequentially, because random causes a lot of system work.
        //
        // It's up to the caller to write the data into the buffer.
        //
    }

    char* getBuffer() 
    {
        return buffer;
    }

    ~BufferDataReader()
    {
        if (buffer != NULL) {
            BigDealloc(buffer);
            buffer = NULL;
        }
    }

    bool init(const char* fileName) { return true; }
    char* readHeader(_int64* io_headerSize)
    {
        if (io_headerSize != NULL)
        {
            *io_headerSize = 0;
        }

        return NULL;
    }

    // seek to a particular range in the file
    virtual void reinit(_int64 startingOffset, _int64 amountOfFileToProcess)
    {
        WriteErrorMessage("BufferDataReader: reinit() called.\\n");
        soft_exit(1);
    }

    // get all remaining data in current batch
    // return false if no more data in current batch
    // startBytes is data "owned" by this block in which reads may start
    // validBytes may also include overflow bytes to handle records spanning batches
    // if you advance() past startBytes, nextBatch() will start offset at that point
    virtual bool getData(char** o_buffer, _int64* o_validBytes, _int64* o_startBytes = NULL)
    {
        *o_buffer = buffer + readOffsetInBuffer;
        *o_validBytes = dataSize - readOffsetInBuffer;

        return *o_validBytes > 0;
    }

    // advance through data in current batch, reducing results from next getData call
    virtual void advance(_int64 bytes)
    {
        _ASSERT(readOffsetInBuffer + bytes <= dataSize);
        readOffsetInBuffer += bytes;
    }

    // advance to next batch
    // by default automatically releases previous batch
    // We never have more than our initial data.
    virtual void nextBatch()
    {
        readOffsetInBuffer = dataSize;
    }

    // whether current batch is last in file
    virtual bool isEOF()
    {
        return readOffsetInBuffer >= dataSize;
    }

    // get current batch identifier
    virtual DataBatch getBatch()
    {
        return DataBatch(); // No real identifier here
    }

    // hold buffers associated with this batch for reuse, increments refcount
    // NOTE: this may be called from another thread,
    // so anything it touches must be thread-safe!
    virtual void holdBatch(DataBatch batch)
    {
        WriteErrorMessage("BufferDataReader: holdbatch not supported\n");
        soft_exit(1);
    }

    // release buffers associated with this batch for reuse
    // decrements refcount, returns true if last release
    // NOTE: this may be called from another thread,
    // so anything it touches must be thread-safe!
    virtual bool releaseBatch(DataBatch batch)
    {
        WriteErrorMessage("BufferDataReader: releaseBatch not supported\n");
        soft_exit(1);
        return true;
    }
    // get current offset into file
    virtual _int64 getFileOffset()
    {
        WriteErrorMessage("BufferDataReader: getFileOffset not supported\n");
        soft_exit(1);
        return -1;
    }

    // get pointer to extra data area for current batch
    // todo: allow this to grow dynamically while keeping stable pointers to previous data
    virtual void getExtra(char** o_extra, _int64* o_length)
    {
        WriteErrorMessage("BufferDataReader: getExtra not supported\n");
        soft_exit(1);
    }

    // get filename for debugging / error printing
    virtual const char* getFilename()
    {
        return "BufferDataReader";
    }
private:

    char* buffer;
    size_t dataSize;
    size_t readOffsetInBuffer;
}; // BufferDataReader


    void
SortedDataFilter::onAdvance(
    DataWriter* writer,
    size_t batchOffset,
    char* data,
    GenomeDistance bytes,
    GenomeLocation location)
{

    OriginalContigNum originalContigNum;
    int pos;

    if (location == 0) {
        originalContigNum = OriginalContigNum(0);
        pos = 0;
    } else if (location == InvalidGenomeLocation) {
        originalContigNum = OriginalContigNum(-1);
        pos = 0;
    } else {
        const Genome::Contig* contig = parent->getGenome()->getContigAtLocation(location);
        originalContigNum = contig->originalContigNumber;
        pos = (int)(location - contig->beginningLocation + 1);
    }

    SortEntry entry(batchOffset, bytes, ContigAndPos(originalContigNum, pos));

#ifdef VALIDATE_SORT
		if (memcmp(data, "BAM", 3) != 0 && memcmp(data, "@HD", 3) != 0) { // skip header block
            GenomeLocation loc;
			GenomeDistance len;
			parent->format->getSortInfo(parent->genome, data, bytes, &loc, &len);
			_ASSERT(loc == location);
		}
#endif
    locations.push_back(entry);
}

    size_t
SortedDataFilter::onNextBatch(
    DataWriter* writer,
    size_t offset,
    size_t bytes,
    bool lastBatch,
    bool* needMoreBuffer,
    size_t* fromBufferUsed)
{
    _ASSERT(!seenLastBatch);
    seenLastBatch |= lastBatch;

    // sort buffered reads by location for later merge sort
    std::stable_sort(locations.begin(), locations.end(), SortEntry::comparator);
    
    // copy from previous buffer into current in sorted order
    char* fromBuffer;
    size_t fromSize, fromUsed;
    char* toBuffer;
    size_t toSize, toUsed;
    BufferDataReader* reader;

    if (!writer->getBatch(-1, &fromBuffer, &fromSize, &fromUsed)) 
    {
        WriteErrorMessage("SortedDataFilter::onNextBatch getBatch of old buffer failed\n");
    }

    if (!lastBatch || offset == 0 || bytes == 0 || *needMoreBuffer) {    // Don't do the last batch optimization at offset 0, because we have special handling for the header.
        if (!writer->getBatch(0, &toBuffer, &toSize, &toUsed))
        {
            WriteErrorMessage("SortedDataFilter::onNextBatch getBatch of new buffer failed\n");
            soft_exit(1);
        }
        reader = NULL;
    } else {
        //
        // For the last batch, we just copy the data into memory instead of writing it to disk and use a BufferDataReader.
        // Get the data reader here, which allocates the buffer that we'll then copy the data into in sorted order.
        //
        reader = new BufferDataReader(bytes);
        toSize = fromSize;
        toBuffer = reader->getBuffer();
        toUsed = 0;
    }

    size_t target = 0;
	GenomeLocation previous = 0;
    for (VariableSizeVector<SortEntry>::iterator i = locations.begin(); i != locations.end(); i++) {
#ifdef VALIDATE_SORT
		if (locations.size() > 1) { // skip header block
            GenomeLocation loc;
            GenomeDistance len;
			parent->format->getSortInfo(parent->genome, fromBuffer + i->offset, i->length, &loc, &len);
			_ASSERT(loc == i->location && loc >= previous);
			previous = loc;
		}
#endif
        memcpy(toBuffer + target, fromBuffer + i->offset, i->length);
        target += i->length;
    }
    
    // remember block extent for later merge sort

    // handle header specially
    size_t header = offset > 0 ? 0 : locations[0].length;
    if (header > 0) {
        parent->setHeaderSize(header);
    }
	int first = offset == 0;
#ifdef VALIDATE_SORT
    GenomeLocationOrderedByOriginalContigs minLocation = locations.size() > first ? locations[first].location : 0;
    GenomeLocationOrderedByOriginalContigs maxLocation = GenomeLocationOrderedByOriginalContigs(locations.size() > first ? locations[locations.size() - 1].location : UINT32_MAX, genome);
    parent->addBlock(offset + header, bytes - header, minLocation, maxLocation);
#else
    parent->addBlock(offset + header, bytes - header, reader);
#endif
    locations.clear();

    return reader == NULL ? target : UINT64_MAX;
}
    
    DataWriter::Filter*
SortedDataFilterSupplier::getFilter()
{
    return new SortedDataFilter(this);
}

    void
SortedDataFilterSupplier::onClosed(
    DataWriterSupplier* supplier)
{
    if (blocks.size() == 1 && sortedFilterSupplier == NULL && false /* this doens't work anymore with the keep-the-last-block-in-memory code*/) {
        // just rename/move temp file to real file, we're done
        DeleteSingleFile(sortedFileName); // if it exists
        if (! MoveSingleFile(tempFileName, sortedFileName)) {
            WriteErrorMessage( "unable to move temp file %s to final sorted file %s\n", tempFileName, sortedFileName);
            soft_exit(1);
        }
        return;
    }
    // merge sort into final file
    if (! mergeSort()) {
        WriteErrorMessage( "merge sort failed\n");
        soft_exit(1);
    }
}

    void
SortedDataFilterSupplier::addBlock(
    size_t start,
    size_t bytes
#ifdef VALIDATE_SORT
	, GenomeLocationOrderedByOriginalContigs minLocation
	, GenomeLocationOrderedByOriginalContigs maxLocation
#endif
    , DataReader *reader
	)
{
    if (bytes > 0) {
        AcquireExclusiveLock(&lock);
#if VALIDATE_SORT
		for (SortBlockVector::iterator i = blocks.begin(); i != blocks.end(); i++) {
			_ASSERT(i->start + i->length <= start || start + bytes <= i->start);
		}
#endif
        SortBlock block(reader);
        block.start = start;
        block.bytes = bytes;
#if VALIDATE_SORT
		block.minLocation = minLocation;
		block.maxLocation = maxLocation;
#endif
        blocks.push_back(block);
        ReleaseExclusiveLock(&lock);
    }
}

struct MergeSortThreadState
{
    SortedDataFilterSupplier* filterSupplier;
    bool deleteSortBlockVector;
    SortBlockVector* blocksForThisThread;
    DataWriter* writer;
};

_int64 mergeSortStartTime;

    void
SortedDataFilterSupplier::MergeSortThreadMain(void* threadParameter)
{
    MergeSortThreadState* state = (MergeSortThreadState*)threadParameter;
 //   fprintf(stderr, "%lld: MergeSortThread %d, deleteSortBlockVector %d\n", timeInMillis() - mergeSortStartTime, GetCurrentThreadId(), state->deleteSortBlockVector);
    state->filterSupplier->mergeSortThread(state->blocksForThisThread, state->writer);

    if (state->deleteSortBlockVector)
    {
        delete state->blocksForThisThread;
    }
    delete state;
}
    //
    // Merge a set of reads coming from readers (either queue or file) into a writer (also either a queue or a file).
    //
    void
SortedDataFilterSupplier::mergeSortThread(SortBlockVector* blocksForThisThread, DataWriter* writer) 
{
    _int64 readWaitTime = 0;
    _int64 writeWaitTime = 0;

    // merge temp blocks into output
    _int64 total = 0;
    // get initial merge sort data
    typedef PriorityQueue<ContigAndPos, _int64> BlockQueue;
    BlockQueue queue;
    for (SortBlockVector::iterator b = blocksForThisThread->begin(); b != blocksForThisThread->end(); b++) {
        _int64 bytes;
        if (!b->reader->getData(&b->data, &bytes))
        {
            _int64 start = timeInMillis();
            b->reader->nextBatch();
            if (!b->reader->getData(&b->data, &bytes))
            {
                WriteErrorMessage("mergeSortThread: unable to get initial data from reader\n");
                soft_exit(1);
            }
            readWaitTime += timeInMillis() - start;
        }

        OriginalContigNum originalContigNum;
        int pos;
        format->getSortInfo(genome, b->data, bytes, NULL, &b->length, &originalContigNum, &pos);
        b->contigAndPos = ContigAndPos(originalContigNum, pos);

        queue.add((_uint32)(b - blocksForThisThread->begin()), b->contigAndPos);
    }

    ContigAndPos current = ContigAndPos(0,0); // current location for validation.  Starts at minimum value.

    int lastRefID = -1, lastPos = 0;
    while (queue.size() > 0) {
#if VALIDATE_SORT
        GenomeLocation check;
        queue.peek(&GenomeLocationOrderedByOriginalContigs(check, genome));
        _ASSERT(check >= current);
#endif
        ContigAndPos secondLocation;
        _int64 smallestIndex = queue.pop();
        _int64 secondIndex = queue.size() > 0 ? queue.peek(&secondLocation) : -1;
        ContigAndPos limit = secondIndex != -1 ? secondLocation : ContigAndPos(-1,0);
        SortBlock* b = &((*blocksForThisThread)[smallestIndex]);
        char* writeBuffer;
        size_t writeBytes;
        writer->getBuffer(&writeBuffer, &writeBytes);
        const int NBLOCKS = 20;
        SortBlock oldBlocks[NBLOCKS];
        int oldBlockIndex = 0;
        while (b->contigAndPos <= limit) {
#if VALIDATE_SORT
            _ASSERT(b->location >= b->minLocation && b->location <= b->maxLocation);
#endif
            if (writeBytes < (size_t)b->length) {
                _int64 start = timeInMillis();
                writer->nextBatch();
                writeWaitTime += timeInMillis() - start;
                writer->getBuffer(&writeBuffer, &writeBytes);
                if (writeBytes < (size_t)b->length) {
                    WriteErrorMessage("mergeSort: buffer size too small\n");
                    soft_exit(1);
                }
            }
            memcpy(writeBuffer, b->data, b->length);
#ifdef VALIDATE_BAM
            if (format == FileFormat::BAM[0] || format == FileFormat::BAM[1]) {
                ((BAMAlignment*)b->data)->validate();
            }
#endif
#if VALIDATE_SORT
            int refID, pos;
            format->getSortInfo(genome, b->data, b->length, NULL, NULL, &refID, &pos);
            _ASSERT(refID == -1 || refID > lastRefID || (refID == lastRefID && pos >= lastPos));
            if (refID != -1) {
                lastRefID = refID;
                lastPos = pos;
            }
#endif
            total++;
            writer->advance(b->length);
            writeBytes -= b->length;
            writeBuffer += b->length;
            oldBlocks[oldBlockIndex] = *b;
            oldBlockIndex = (oldBlockIndex + 1) % NBLOCKS;
            b->reader->advance(b->length);
            _ASSERT(b->contigAndPos >= current);
            current = b->contigAndPos;
            _int64 readBytes;
            if (!b->reader->getData(&b->data, &readBytes)) {
                _int64 start = timeInMillis();
                b->reader->nextBatch();
                readWaitTime += timeInMillis() - start;
                if (!b->reader->getData(&b->data, &readBytes)) {
                    // This isn't supported in DataQueueReader, so this assert is off_ASSERT(b->reader->isEOF());
                    delete b->reader;
                    b->reader = NULL;
                    break;
                }
            }
            ContigAndPos previous = b->contigAndPos;
            OriginalContigNum originalContigNum;
            int pos;

            format->getSortInfo(genome, b->data, readBytes, NULL, &b->length, &originalContigNum, &pos);
            b->contigAndPos = ContigAndPos(originalContigNum, pos);
            _ASSERT(b->length <= readBytes && b->contigAndPos >= previous);
        }
        if (b->reader != NULL) {
            queue.add(smallestIndex, b->contigAndPos);
        }
    }

    InterlockedAdd64AndReturnNewValue(&totalReadsSorted, total);

    // writer->nextBatch();
    // close everything
    _int64 start = timeInMillis();
    writer->close();
    writeWaitTime += timeInMillis() - start;
    delete writer;

    //fprintf(stderr, "%lld: Thread %d read %lldms, write %lldms\n", timeInMillis() - mergeSortStartTime, GetCurrentThreadId(), readWaitTime, writeWaitTime);
}

    bool
SortedDataFilterSupplier::mergeSort()
{
    mergeSortStartTime = timeInMillis();
    // merge sort from temp file into sorted file

    WriteStatusMessage("sorting...");
    _int64 start = timeInMillis();
    _int64 startReadWaitTime = DataReader::ReadWaitTime;
    _int64 startReleaseWaitTime = DataReader::ReleaseWaitTime;
    _int64 startWriteWaitTime = DataWriter::WaitTime;
    _int64 startWriteFilterTime = DataWriter::FilterTime;

    // set up buffered output
    DataWriterSupplier* writerSupplier = DataWriterSupplier::create(sortedFileName, bufferSize, emitInternalScore, internalScoreTag ,sortedFilterSupplier,
        encoder, encoder != NULL ? 6 : 4); // use more buffers to let encoder run async
    DataWriter* writer = writerSupplier->getWriter();
    if (writer == NULL) {
        WriteErrorMessage( "open sorted file for write failed\n");
        return false;
    }
    DataSupplier* readerSupplier = DataSupplier::Default; // autorelease
    // setup - open all files, read first block, begin read for second
    if (blocks.size() > 5000) {
        WriteErrorMessage("warning: merging %d blocks could be slow, try increasing sort memory with -sm option\n", blocks.size());
    }
    for (SortBlockVector::iterator i = blocks.begin(); i != blocks.end(); i++) {
        if (i->reader == NULL) // Otheriwse, it's the last block that's in memory
        {
            i->reader = readerSupplier->getDataReader(1, MAX_READ_LENGTH * 8, 0.0,
                __min(1UL << 23, __max(1UL << 17, bufferSpace / blocks.size()))); // 128kB to 8MB buffer space per block
            if (!i->reader->init(tempFileName)) {
                WriteErrorMessage("SortedDataFilterSupplier::mergeSort: reader->init(%s) failed\n", tempFileName);
                soft_exit(1);
            }
            i->reader->reinit(i->start, i->bytes);
        }
    }

    // write out header
    if (headerSize > 0xffffffff) {
        WriteErrorMessage("SortedDataFilterSupplier: headerSize too big\n");
        soft_exit(1);
    }
    if (headerSize > 0) {
        DataReader* headerReader;
        if (blocks[0].dataReaderIsBuffer) 
        {
            headerReader = readerSupplier->getDataReader(1, MAX_READ_LENGTH * 8, 0.0, headerSize + 4096);
            if (!headerReader->init(tempFileName)) {
                WriteErrorMessage("SortedDataFilterSupplier::mergeSort: reader->init(%s) failed for headerReader\n", tempFileName);
                soft_exit(1);
            }
        } else {
            headerReader = blocks[0].reader;
        }
        headerReader->reinit(0, headerSize);
		writer->inHeader(true);
        char* rbuffer;
        _int64 rbytes;
        char* wbuffer;
        size_t wbytes;
		for (size_t left = headerSize; left > 0; ) {
			if ((!headerReader->getData(&rbuffer, &rbytes)) || rbytes == 0) {
                headerReader->nextBatch();
				if (!headerReader->getData(&rbuffer, &rbytes)) {
					WriteErrorMessage( "read header failed, left %lld, headerSize %lld\n", left, headerSize);
                    headerReader->dumpState();
					soft_exit(1);
				}
			}
			if ((! writer->getBuffer(&wbuffer, &wbytes)) || wbytes == 0) {
				writer->nextBatch();
				if (! writer->getBuffer(&wbuffer, &wbytes)) {
					WriteErrorMessage( "write header failed\n");
					soft_exit(1);
				}
			}
			size_t xfer = min(left, min((size_t) rbytes, wbytes));
			_ASSERT(xfer > 0 && xfer <= UINT32_MAX);
			memcpy(wbuffer, rbuffer, xfer);
            headerReader->advance(xfer);
			writer->advance((unsigned) xfer);
			left -= xfer;
		}

        if (blocks[0].dataReaderIsBuffer) 
        {
            delete headerReader;
        } else {
            blocks[0].reader->reinit(blocks[0].start, blocks[0].bytes);
        }
        headerReader = NULL;

		writer->nextBatch();
		writer->inHeader(false);
    }

    //
    // Set up the merge sort threads.  We do a two level tree: the leaves split the input files evenly and the root merges them all together.
    // Unless there are very few input files, in which case we just run one thread.
    //

    if (blocks.size() < 1 || numThreads == 1 || true /* Multi-thread doesn't seem to help, so just stick with single thread */)  
    {
        //
        // The single thread case.
        //
        MergeSortThreadState* rootThreadState = new MergeSortThreadState();
        rootThreadState->blocksForThisThread = &blocks;
        rootThreadState->writer = writer;
        rootThreadState->deleteSortBlockVector = false;
        rootThreadState->filterSupplier = this;

        MergeSortThreadMain(rootThreadState);
    }
    else 
    {
        //
        // Multiple threads.
        //

        int minPerThread = 2;   // BJB - low for testing.
        int nLeafThreads = numThreads - 1;
        int nBlocksAssigned = 0;

        MergeSortThreadState* rootThreadState = new MergeSortThreadState();
        rootThreadState->blocksForThisThread = new SortBlockVector(nLeafThreads);
        rootThreadState->filterSupplier = this;
        rootThreadState->writer = writer;   // The writer that actually writes to the output file.
        rootThreadState->deleteSortBlockVector = false;

        DataQueue** dataQueues = new DataQueue * [nLeafThreads];

        for (int i = 0; i < nLeafThreads; i++)
        {
            int nBlocksThisThread = (int)(__min(blocks.size() - nBlocksAssigned, __max(minPerThread, (blocks.size() - nBlocksAssigned + (nLeafThreads - i - 1)) / (nLeafThreads - i))));
            if (nBlocksThisThread == 0) 
            {
                continue;
                dataQueues[i] = NULL;
            }

            dataQueues[i] = new DataQueue(5, (size_t)16 * 1024 * 1024); // Buffer count/size is kinda arbitrary

            MergeSortThreadState* leafThreadState = new MergeSortThreadState();
            leafThreadState->deleteSortBlockVector = true;
            leafThreadState->filterSupplier = this;

            size_t totalSize = 0;
            leafThreadState->blocksForThisThread = new SortBlockVector(nBlocksThisThread);
            for (int blockIndex = nBlocksAssigned; blockIndex < nBlocksAssigned + nBlocksThisThread; blockIndex++) 
            {
                leafThreadState->blocksForThisThread->push_back(blocks[blockIndex]);
                totalSize += blocks[blockIndex].bytes;
            }

            leafThreadState->writer = dataQueues[i]->getWriter();
            SortBlock outputBlock;
            outputBlock.bytes = totalSize;
            outputBlock.reader = dataQueues[i]->getReader();

            rootThreadState->blocksForThisThread->push_back(outputBlock);

            nBlocksAssigned += nBlocksThisThread;

            if (!StartNewThread(MergeSortThreadMain, leafThreadState)) 
            {
                WriteErrorMessage("merge sort: StartNewThread failed.\n");
                soft_exit(1);
            }
        } // for each worker thread

        //
        // Just run the root on this thread.
        //
        MergeSortThreadMain(rootThreadState);
    } // The multi-thread case

    writerSupplier->close();
    delete writerSupplier;
    if (! DeleteSingleFile(tempFileName)) {
        WriteErrorMessage( "warning: failure deleting temp file %s\n", tempFileName);
    }

    WriteStatusMessage("sorted %lld reads in %u blocks, %lld s\n"
        /*"read wait align %.3f s + merge %.3f s, read release align %.3f s + merge %.3f s\n"
        "write wait %.3f s align + %.3f s merge, write filter %.3f s align + %.3f s merge\n"*/,
        totalReadsSorted, blocks.size(), (timeInMillis() - start)/1000 /*,
        startReadWaitTime * 1e-9, (DataReader::ReadWaitTime - startReadWaitTime) * 1e-9,
        startReleaseWaitTime * 1e-9, (DataReader::ReleaseWaitTime - startReleaseWaitTime) * 1e-9,
        startWriteWaitTime * 1e-9, (DataWriter::WaitTime - startWriteWaitTime) * 1e-9,
        startWriteFilterTime * 1e-9, (DataWriter::FilterTime - startWriteFilterTime) * 1e-9*/);
    return true;
}

    DataWriterSupplier*
DataWriterSupplier::sorted(
    const FileFormat* format,
    const Genome* genome,
    const char* tempFileName,
    size_t tempBufferMemory,
    int numThreads,
    const char* sortedFileName,
    DataWriter::FilterSupplier* sortedFilterSuppler,
    size_t maxBufferSize,
    bool emitInternalScore,
    char *internalScoreTag,
    FileEncoder* encoder)
{
    const int bufferCount = 3;
    const size_t bufferSpace = tempBufferMemory > 0 ? tempBufferMemory : (numThreads * (size_t)1 << 30);
    const size_t bufferSize = bufferSpace / ((size_t)bufferCount * numThreads);
    DataWriter::FilterSupplier* filterSupplier =
        new SortedDataFilterSupplier(format, genome, tempFileName, sortedFileName, sortedFilterSuppler, bufferSize, bufferSpace, emitInternalScore, internalScoreTag, encoder, numThreads);
    return DataWriterSupplier::create(tempFileName, bufferSize, emitInternalScore, internalScoreTag, filterSupplier, NULL, bufferCount);
}
