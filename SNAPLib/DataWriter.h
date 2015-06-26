/*++

Module Name:

    DataWriter.h

Abstract:

    Headers for the DataWriter & related classes for the SNAP sequencer

Authors:

    Ravi Pandya, Feb 2013

Environment:

    User mode service.

Revision History:

--*/

#pragma once

#include "Compat.h"
#include "Read.h"
#include "ParallelTask.h"
#include "Genome.h"

class DataWriterSupplier;

// per-thread writer for data into a single destination
class DataWriter
{
public:
    
    enum FilterType
    {
        ReadFilter, // reads data but does not modify it
        ModifyFilter, // modifies data in place
        CopyFilter, // copies data into new buffer, same size
        TransformFilter, // copies data into new buffer, possibly different size
        ResizeFilter, // rewrites data in same buffer, possibly different size
    };
    // single filter instance per thread
    // points to filterSupplier for common data
    class Filter
    {
    public:
        Filter (FilterType i_filterType) : filterType(i_filterType) {}

        const FilterType filterType;

        virtual ~Filter() {}

		// called to set whether we're writing a header vs. individual reads
		virtual void inHeader(bool flag) {} // default do nothing

        // called when a chunk of data (i.e. a single read) has been written into the file
        virtual void onAdvance(DataWriter* writer, size_t batchOffset, char* data, GenomeDistance bytes, GenomeLocation location) = 0;

        // called when a batch has been completed, after advancing to the next
        // e.g. so use getBatch(-1, ...) to get the one that was just completed
        // TransformFilters return #byte of transformed data in current buffer, so we need to advance again
        // TransformFilters should call getBatch(0) to ensure current buffer has been written before they write into it
        virtual size_t onNextBatch(DataWriter* writer, size_t offset, size_t bytes) = 0;
    };
    
    // factory for per-thread filters
    class FilterSupplier
    {
    public:
        FilterSupplier (FilterType i_filterType) : filterType(i_filterType) {}
        
        const FilterType filterType;

        virtual ~FilterSupplier() {}

        FilterSupplier* compose(FilterSupplier* other);

        virtual Filter* getFilter() = 0;

        // called when entire file is done; onClosing before file is closed, onClosed after
        virtual void onClosing(DataWriterSupplier* supplier) = 0;
        virtual void onClosed(DataWriterSupplier* supplier) = 0;
    };

    DataWriter(Filter* i_filter) : filter(i_filter) {}

    virtual ~DataWriter() {}

	void inHeader(bool flag)
	{ if (filter != NULL) { filter->inHeader(flag); } }

    // get remaining space in current buffer for writing
    virtual bool getBuffer(char** o_buffer, size_t* o_size) = 0;

    // advance within current buffer, reducing available space
    // should be called on each read, with the location
    virtual void advance(GenomeDistance bytes, GenomeLocation location = 0) = 0;

    // get complete data buffer in batch, relative==0 is current, relative==-1 is previous, etc.
    // if negative gets old data written, else waits for write to complete so you can write into it
    // o_offset gets physical offset (e.g. compressed), o_logical gets data offset (e.g. uncompressed)
    virtual bool getBatch(int relative, char** o_buffer, size_t* o_size = NULL, size_t* o_used = NULL, size_t* o_offset = NULL, size_t* o_logicalUsed = 0, size_t* o_logicalOffset = NULL) = 0;

    // advance to next buffer
    virtual bool nextBatch() = 0;

    // this thread is complete
    virtual void close() = 0;

    // nanosecond timers
    static volatile _int64 FilterTime;
    static volatile _int64 WaitTime;

protected:
    Filter* filter;
};

class FileFormat;
class Genome;
class GzipWriterFilterSupplier;
class FileEncoder;

// creates writers for multiple threads
class DataWriterSupplier
{
public:
    virtual ~DataWriterSupplier() {}

    virtual DataWriter* getWriter() = 0;

    // call when all threads are done, all filters destroyed
    virtual void close() = 0;
    
    static DataWriterSupplier* create(
        const char* filename,
        size_t bufferSize,
        DataWriter::FilterSupplier* filterSupplier = NULL,
        FileEncoder* encoder = NULL,
        int count = 4);
    
    static DataWriterSupplier* sorted(
        const FileFormat* format,
        const Genome* genome,
        const char* tempFileName,
        size_t tempBufferMemory,
        int numThreads,
        const char* sortedFileName,
        DataWriter::FilterSupplier* sortedFilterSupplier,
        size_t maxBufferSize,
        FileEncoder* encoder = NULL);

    // defaults follow BAM output spec
    static GzipWriterFilterSupplier* gzip(bool bamFormat, size_t chunkSize, int numThreads, bool bindToProcessors, bool multiThreaded);

    static DataWriter::FilterSupplier* markDuplicates(const Genome* genome);

    static DataWriter::FilterSupplier* bamIndex(const char* indexFileName, const Genome* genome, GzipWriterFilterSupplier* gzipSupplier);
};

class AsyncDataWriter;

class FileEncoder
{
public:
    FileEncoder(int numThreads, bool bindToProcessors, ParallelWorkerManager* i_supplier);

    ~FileEncoder()
    {
        if (coworker != NULL) {
            _ASSERT(! encoderRunning); coworker->stop(); delete coworker;
        }
    }

    static FileEncoder* gzip(GzipWriterFilterSupplier* filterSupplier, int numThreads, bool bindToProcessor, size_t chunkSize = 65536, bool bam = true);

    // post-construction initialization
    void initialize(AsyncDataWriter* i_writer);

    // called by writer when there is data to encode; threadsafe
    void inputReady();

    void close();

    // called by supplier to get/set information about current batch

    void setupEncode(int relative);
    
    void getEncodeBatch(char** o_batch, size_t* o_batchSize, size_t* o_batchUsed);

    void getOffsets(size_t* o_logicalOffset, size_t* o_physicalOffset);

    void setEncodedBatchSize(size_t newSize);

private:
    // static callback for encoder; threadsafe
    static void outputReadyCallback(void *p);

    // called by encoder when a block of data has been encoded; threadsafe
    void outputReady();

    // scans writer and kicks off encoder if there is something ready; must hold lock
    void checkForInput();

    AsyncDataWriter* writer;
    ParallelCoworker* coworker;
    ExclusiveLock* lock;
    bool encoderRunning;
    int encoderBatch;

    friend class AsyncDataWriter;
};

class StdoutAsyncFile : public AsyncFile
{
public:
    StdoutAsyncFile();
    virtual ~StdoutAsyncFile();

    bool close();

    static StdoutAsyncFile *open(const char *filename, bool write);

    AsyncFile::Writer* getWriter();
    AsyncFile::Reader* getReader();

    void beginWrite(void *buffer, size_t length, size_t offset, size_t *o_bytesWritten);
    void waitForCompletion(size_t offset);

private:
    ExclusiveLock   lock;

    struct WriteElement {
        void                *buffer;
        size_t               length;
        size_t               offset;
        size_t              *o_bytesWritten;

        WriteElement        *next;
        WriteElement        *prev;

        void enqueue(WriteElement *previous);
        void dequeue();
    };

    bool isQueueEmpty() {
        return writeElementQueue->next == writeElementQueue;
    }

    size_t          highestOffsetCompleted;

    //
    // The queue is kept in order, and the writer writes without gaps, so if you put on blocks 10 and 12, the writer will write
    // 10, and then leave 12 on the queue and wait for 11 to be added and written before processing 12.
    //
    WriteElement writeElementQueue[1];

    EventObject  unexaminedElementsOnQueue;     // This gets set when a writer puts a block on the queue, and cleared when the consumer has seen it.
    EventObject  elementsCompleted;             // This gets set when any element is completed by the consumer, and reset when a waiter starts

    SingleWaiterObject  consumerThreadDone;

    bool                closing;

    static void ConsumerThreadMain(void *param);
    void runConsumer();

    static bool anyCreated;              // Because there's no way to multiplex stdout, you only get one per run of SNAP
};
