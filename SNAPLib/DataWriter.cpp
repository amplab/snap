/*++

Module Name:

    DataWriter.cpp

Abstract:

    General file writer.

Environment:

    User mode service.

    Not thread safe.

--*/

#include "stdafx.h"
#include "BigAlloc.h"
#include "Compat.h"
#include "DataWriter.h"
#include "ParallelTask.h"
#include "exit.h"
#include "Bam.h"

using std::min;
using std::max;

class AsyncDataWriterSupplier : public DataWriterSupplier
{
public:
    AsyncDataWriterSupplier(const char* i_filename, DataWriter::FilterSupplier* i_filterSupplier,
        FileEncoder* i_encoder, int i_bufferCount, size_t i_bufferSize);

    virtual DataWriter* getWriter();

    virtual void close();

private:
    friend class AsyncDataWriter;
    friend class FileEncoder;
    void advance(size_t physical, size_t logical, size_t* o_physical, size_t* o_logical);

    const char* filename;
    AsyncFile* file;
    DataWriter::FilterSupplier* filterSupplier;
    FileEncoder* encoder;
    const int bufferCount;
    const size_t bufferSize;
    ExclusiveLock lock;
    size_t sharedOffset;
    size_t sharedLogical;
    bool closing;
};

class AsyncDataWriter : public DataWriter
{
public:

    AsyncDataWriter(AsyncFile* i_file, AsyncDataWriterSupplier* i_supplier,
        int i_count, size_t i_bufferSize, Filter* i_filter, FileEncoder* i_encoder);

    virtual ~AsyncDataWriter()
    {
        for (int i = 0; i < count; i++) {
            delete batches[i].file;
        }
        BigDealloc(batches[0].buffer); // all in one big block
        delete [] batches;
        if (encoder != NULL) {
            delete encoder;
        }
        DestroyExclusiveLock(&lock);
    }

    virtual bool getBuffer(char** o_buffer, size_t* o_size);

    virtual void advance(unsigned bytes, unsigned location = 0);

    virtual bool getBatch(int relative, char** o_buffer, size_t* o_size, size_t* o_used, size_t* o_offset, size_t* o_logicalUsed = 0, size_t* o_logicalOffset = NULL);

    virtual bool nextBatch();
    
    virtual void close();

private:

    void acquireLock()
    { if (encoder != NULL) { AcquireExclusiveLock(&lock); } }

    void releaseLock()
    { if (encoder != NULL) { ReleaseExclusiveLock(&lock); } }

    struct Batch
    {
        char* buffer;
        AsyncFile::Writer* file;
        size_t used;
        size_t fileOffset;
        size_t logicalUsed;
        size_t logicalOffset;
        EventObject encoded;
    };
    Batch* batches;
    const int count;
    const size_t bufferSize;
    AsyncDataWriterSupplier* supplier;
    int current;
    FileEncoder* encoder;
    ExclusiveLock lock;

    friend class FileEncoder;
};

FileEncoder::FileEncoder(
    int numThreads,
    bool bindToProcessors,
    ParallelWorkerManager* i_manager)
    :
    encoderRunning(false),
    coworker(numThreads == 0 ? NULL
        : new ParallelCoworker(numThreads, bindToProcessors, i_manager, FileEncoder::outputReadyCallback, this))
{}

    void
FileEncoder::initialize(
    AsyncDataWriter* i_writer)
{
    writer = i_writer;
    lock = &writer->lock;
    encoderBatch = writer->count - 1;
    if (coworker != NULL) {
        coworker->getManager()->initialize(this);
        coworker->start();
    }
}

    void
FileEncoder::inputReady()
{
    AcquireExclusiveLock(lock);
    if (! encoderRunning) {
        checkForInput();
    }
    ReleaseExclusiveLock(lock);
}

    void
FileEncoder::close()
{
    // wait for pending encodes
    AcquireExclusiveLock(lock);
    int start = encoderBatch;
    int pending = (writer->current + writer->count - start) % writer->count;
    ReleaseExclusiveLock(lock);
    for (int i = 0; i < pending; i++) {
        WaitForEvent(&writer->batches[(start + i) % writer->count].encoded);
    }
    coworker->stop();
}

    void
FileEncoder::outputReadyCallback(
    void *p)
{
    ((FileEncoder*) p)->outputReady();
}

    void
FileEncoder::outputReady()
{
    AcquireExclusiveLock(lock);

    encoderRunning = false;

    // begin writing the buffer to disk
    AsyncDataWriter::Batch* write = &writer->batches[encoderBatch];
    writer->supplier->advance(write->used, 0, &write->fileOffset, &write->logicalOffset);
    //printf("outputReady write batch %d @%lld:%lld\n", encoderBatch, write->fileOffset, write->used);
    if (! write->file->beginWrite(write->buffer, write->used, write->fileOffset, NULL)) {
        fprintf(stderr, "error: file write %lld bytes at offset %lld failed\n", write->used, write->fileOffset);
        soft_exit(1);
    }
    AllowEventWaitersToProceed(&write->encoded);

    // check for more work
    checkForInput();

    ReleaseExclusiveLock(lock);
}

    void
FileEncoder::checkForInput()
{
    // look for another block ready to encode
    while (true) {
        int nextBatch = (encoderBatch + 1) % writer->count;
        if (nextBatch == writer->current) {
            break;
        }
        encoderBatch = nextBatch;
        AsyncDataWriter::Batch* encode = &writer->batches[encoderBatch];
        if (encode->used > 0) {
            encoderRunning = true;
            coworker->step();
            break;
        }
    }
}

    void
FileEncoder::setupEncode(
    int relative)
{
    encoderBatch = (writer->current + relative + writer->count) % writer->count;
}

    void
FileEncoder::getEncodeBatch(
    char** o_batch,
    size_t* o_batchSize,
    size_t* o_batchUsed)
{
    AsyncDataWriter::Batch* batch = &writer->batches[encoderBatch];
    *o_batch = batch->buffer;
    *o_batchSize = writer->bufferSize;
    *o_batchUsed = batch->used;
    //printf("getEncodeBatch #%d: %lld/%lld\n", encoderBatch, batch->used, writer->bufferSize);
}

    void
FileEncoder::getOffsets(
    size_t* o_logicalOffset,
    size_t* o_physicalOffset)
{
    // logical has already been set correctly in batch
    *o_logicalOffset = writer->batches[encoderBatch].logicalOffset;
    // physical is not yet updated, use shared
    *o_physicalOffset = writer->supplier->sharedOffset;
}

    void
FileEncoder::setEncodedBatchSize(
    size_t newSize)
{
    size_t old = writer->batches[encoderBatch].used;
//printf("setEncodedBatchSize #%d %lld -> %lld\n", encoderBatch, old, newSize);
    if (newSize != old) {
        AcquireExclusiveLock(lock);
        AsyncDataWriter::Batch* batch = &writer->batches[encoderBatch];
        batch->logicalUsed = batch->used;
        batch->used = newSize;
        ReleaseExclusiveLock(lock);
    }
}

AsyncDataWriter::AsyncDataWriter(
    AsyncFile* i_file,
    AsyncDataWriterSupplier* i_supplier, 
    int i_count,
    size_t i_bufferSize,
    Filter* i_filter,
    FileEncoder* i_encoder)
    :
    DataWriter(i_filter),
    encoder(i_encoder),
    supplier(i_supplier),
    count(i_count),
    bufferSize(i_bufferSize),
    current(0)
{
    _ASSERT(count >= 2);
    char* block = (char*) BigAlloc(count * bufferSize);
    if (block == NULL) {
        fprintf(stderr, "Unable to allocate %lld bytes for write buffers\n", count * bufferSize);
        soft_exit(1);
    }
    batches = new Batch[count];
    for (int i = 0; i < count; i++) {
        batches[i].buffer = block + i * bufferSize;
        batches[i].file = i_file->getWriter();
        batches[i].used = 0;
        batches[i].fileOffset = 0;
        batches[i].logicalUsed = 0;
        batches[i].logicalOffset = 0;
        if (encoder != NULL) {
            CreateEventObject(&batches[i].encoded);
            AllowEventWaitersToProceed(&batches[i].encoded); // initialize so empty bufs are available
        }
    }

    InitializeExclusiveLock(&lock);
    if (encoder != NULL) {
        encoder->initialize(this);
    }
}
    
    bool
AsyncDataWriter::getBuffer(
    char** o_buffer,
    size_t* o_size)
{
    *o_buffer = batches[current].buffer + batches[current].used;
    *o_size = bufferSize - batches[current].used;
    return true;
}

    void
AsyncDataWriter::advance(
    unsigned bytes,
    unsigned location)
{
    _ASSERT(bytes <= bufferSize - batches[current].used);
    char* data = batches[current].buffer + batches[current].used;
    size_t batchOffset = batches[current].used;
    batches[current].used = min(bufferSize, batchOffset + bytes);
    if (filter != NULL) {
        //_int64 start = timeInNanos();
        filter->onAdvance(this, batchOffset, data, bytes, location);
        //InterlockedAdd64AndReturnNewValue(&FilterTime, timeInNanos() - start);
    }
}

    bool
AsyncDataWriter::getBatch(
    int relative,
    char** o_buffer,
    size_t* o_size,
    size_t* o_used,
    size_t* o_offset,
    size_t* o_logicalUsed,
    size_t* o_logicalOffset)
{
    if (relative < 1 - count || relative > count - 1) {
        return false;
    }
    if (encoder != NULL && relative <= ((encoder->encoderBatch - current + count) % count) - count) {
        return false;
    }
    int index = (current + relative + count) % count; // ensure non-negative
    Batch* batch = &batches[index];
    *o_buffer = batch->buffer;
    if (o_size != NULL) {
        *o_size = bufferSize;
    }
    if (o_used != NULL) {
        *o_used = relative <= 0 ? batch->used : 0;
    }
    if (o_offset != NULL) {
        *o_offset = relative <= 0 ? batch->fileOffset : 0;
    }
    if (o_logicalUsed != NULL) {
        *o_logicalUsed = relative <=0 ? batch->logicalUsed: 0;
    }
    if (o_logicalOffset != NULL) {
        *o_logicalOffset = relative <=0 ? batch->logicalOffset : 0;
    }
    if (relative >= 0) {
        if (encoder != NULL) {
            WaitForEvent(&batch->encoded);
        }
        batch->file->waitForCompletion();
    }
    return true;
}

    bool
AsyncDataWriter::nextBatch()
{
    _int64 start = timeInNanos();
    if (encoder != NULL) {
        WaitForEvent(&batches[(current + 1) % count].encoded);
    }
    acquireLock();
    int written = current;
    Batch* write = &batches[written];
    write->logicalUsed = write->used;
    current = (current + 1) % count;
    //printf("nextBatch reset %d used=0\n", current);
    batches[current].used = 0;
    bool newBuffer = filter != NULL && (filter->filterType == CopyFilter || filter->filterType == TransformFilter);
    bool newSize = filter != NULL && (filter->filterType == TransformFilter || filter->filterType == ResizeFilter);
    if (newSize) {
        // advisory only
        write->fileOffset = supplier->sharedOffset;
        write->logicalOffset = supplier->sharedLogical;
    } else {
        supplier->advance(encoder == NULL ? write->used : 0, write->logicalUsed, &write->fileOffset, &write->logicalOffset);
    }
    if (filter != NULL) {
        size_t n = filter->onNextBatch(this, write->fileOffset, write->used);
	    if (newSize) {
	        write->used = n;
            supplier->advance(encoder == NULL ? write->used : 0, write->logicalUsed, &write->fileOffset, &write->logicalOffset);
	    }
        if (newBuffer) {
            // current has used>0, written has logicalUsed>0, for compressed & uncompressed data respectively
            batches[current].used = write->used;
            batches[current].fileOffset = write->fileOffset;
            batches[current].logicalUsed = 0;
            batches[current].logicalOffset = write->logicalOffset;
            write->used = 0;
            written = current;
            write = &batches[written];
            current = (current + 1) % count;
            batches[current].used = 0;
            batches[current].logicalUsed = 0;
        }
    }
    _int64 start2 = timeInNanos();
    releaseLock();

    InterlockedAdd64AndReturnNewValue(&FilterTime, start2 - start);
    if (encoder == NULL) {
        //printf("nextBatch beginWrite #%d @%lld: %lld bytes\n", write-batches, write->fileOffset, write->used);
        //_ASSERT(BgzfHeader::validate(write->buffer, write->used)); //!! remove before checkin
        if (! write->file->beginWrite(write->buffer, write->used, write->fileOffset, NULL)) {
            fprintf(stderr, "error: file write %lld bytes at offset %lld failed\n", write->used, write->fileOffset);
            soft_exit(1);
        }
    } else {
        PreventEventWaitersFromProceeding(&write->encoded);
        encoder->inputReady();
    }
    if (! batches[current].file->waitForCompletion()) {
        fprintf(stderr, "error: file write failed\n");
        soft_exit(1);
    }
    InterlockedAdd64AndReturnNewValue(&WaitTime, timeInNanos() - start2);
    return true;
}

    void
AsyncDataWriter::close()
{
    nextBatch(); // ensure last buffer gets written
    if (encoder != NULL) {
        encoder->close();
        for (int i = 0; i < count; i++) {
            DestroyEventObject(&batches[i].encoded);
        }
    }
    for (int i = 0; i < count; i++) {
        batches[i].file->close();
    }
}

AsyncDataWriterSupplier::AsyncDataWriterSupplier(
    const char* i_filename,
    DataWriter::FilterSupplier* i_filterSupplier,
    FileEncoder* i_encoder,
    int i_bufferCount,
    size_t i_bufferSize)
    :
    filename(i_filename),
    filterSupplier(i_filterSupplier),
    encoder(i_encoder),
    bufferCount(i_bufferCount),
    bufferSize(i_bufferSize),
    sharedOffset(0),
    sharedLogical(0),
    closing(false)
{
    file = AsyncFile::open(filename, true);
    if (file == NULL) {
        fprintf(stderr, "failed to open %s for write\n", filename);
        soft_exit(1);
    }
    InitializeExclusiveLock(&lock);
}

    DataWriter*
AsyncDataWriterSupplier::getWriter()
{
    return new AsyncDataWriter(file, this, bufferCount, bufferSize,
        filterSupplier && ! closing ? filterSupplier->getFilter() : NULL,
        closing ? NULL : encoder);
}

    void
AsyncDataWriterSupplier::close()
{
    closing = true;
    if (filterSupplier != NULL) {
        filterSupplier->onClosing(this);
    }
    file->close();
    if (filterSupplier != NULL) {
        filterSupplier->onClosed(this);
    }
    DestroyExclusiveLock(&lock);
}
    void
AsyncDataWriterSupplier::advance(
    size_t physical,
    size_t logical,
    size_t* o_physical,
    size_t* o_logical)
{
    AcquireExclusiveLock(&lock);
    *o_physical = sharedOffset;
    sharedOffset += physical;
    *o_logical = sharedLogical;
    sharedLogical += logical;
    //printf("advance %lld + %lld = %lld, logical %lld + %lld = %lld\n", *o_physical, physical, sharedOffset, *o_logical, logical, sharedLogical);
    ReleaseExclusiveLock(&lock);
}

    DataWriterSupplier*
DataWriterSupplier::create(
    const char* filename,
    DataWriter::FilterSupplier* filterSupplier,
    FileEncoder* encoder,
    int count,
    size_t bufferSize)
{
    return new AsyncDataWriterSupplier(filename, filterSupplier, encoder, count, bufferSize);
}

class ComposeFilter : public DataWriter::Filter
{
public:
    ComposeFilter(DataWriter::Filter* i_a, DataWriter::Filter* i_b) :
        Filter(max(i_a->filterType, i_b->filterType)), a(i_a), b(i_b) {}

    virtual ~ComposeFilter()
    { delete a; delete b; }
    
	virtual void inHeader(bool flag)
	{
		a->inHeader(flag);
		b->inHeader(flag);
	}

    virtual void onAdvance(DataWriter* writer, size_t batchOffset, char* data, unsigned bytes, unsigned location)
    {
        a->onAdvance(writer, batchOffset, data, bytes, location);
        b->onAdvance(writer, batchOffset, data, bytes, location);
    }

    virtual size_t onNextBatch(DataWriter* writer, size_t offset, size_t bytes)
    {
        size_t sa = a->onNextBatch(writer, offset, bytes);
        size_t sb = b->onNextBatch(writer, offset, sa);
        return sb;
    }

private:
    DataWriter::Filter* a;
    DataWriter::Filter* b;
};

class ComposeFilterSupplier : public DataWriter::FilterSupplier
{
public:
    ComposeFilterSupplier(DataWriter::FilterSupplier* i_a, DataWriter::FilterSupplier* i_b) :
        FilterSupplier(max(i_a->filterType, i_b->filterType)), a(i_a), b(i_b) {}

    virtual ~ComposeFilterSupplier()
    { delete a; delete b; }
    
    virtual DataWriter::Filter* getFilter()
    { return new ComposeFilter(a->getFilter(), b->getFilter()); }

    virtual void onClosing(DataWriterSupplier* supplier)
    {
        a->onClosing(supplier);
        b->onClosing(supplier);
    }
    
    virtual void onClosed(DataWriterSupplier* supplier)
    {
        a->onClosed(supplier);
        b->onClosed(supplier);
    }

private:
    DataWriter::FilterSupplier* a;
    DataWriter::FilterSupplier* b;
};

    DataWriter::FilterSupplier*
DataWriter::FilterSupplier::compose(
    DataWriter::FilterSupplier* other)
{
    return new ComposeFilterSupplier(this, other);
}

volatile _int64 DataWriter::WaitTime = 0;
volatile _int64 DataWriter::FilterTime = 0;
