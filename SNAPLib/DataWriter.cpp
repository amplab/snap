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
#include "exit.h"

using std::min;
using std::max;

class AsyncDataWriterSupplier : public DataWriterSupplier
{
public:
    AsyncDataWriterSupplier(const char* i_filename, DataWriter::FilterSupplier* i_filterSupplier,
        int i_bufferCount, size_t i_bufferSize);

    virtual DataWriter* getWriter();

    virtual void close();

private:
    friend class AsyncDataWriter;
    void advance(size_t physical, size_t logical, size_t* o_physical, size_t* o_logical);

    const char* filename;
    AsyncFile* file;
    DataWriter::FilterSupplier* filterSupplier;
    const int bufferCount;
    const size_t bufferSize;
    ExclusiveLock lock;
    size_t sharedOffset;
    size_t sharedLogical;
};

class AsyncDataWriter : public DataWriter
{
public:

    AsyncDataWriter(AsyncFile* i_file, AsyncDataWriterSupplier* i_supplier, int i_count, size_t i_bufferSize, Filter* i_filter);

    virtual ~AsyncDataWriter()
    {
        for (int i = 0; i < count; i++) {
            delete batches[i].file;
        }
        BigDealloc(batches[0].buffer); // all in one big block
        delete [] batches;
    }

    virtual bool getBuffer(char** o_buffer, size_t* o_size);

    virtual void advance(unsigned bytes, unsigned location = 0);

    virtual bool getBatch(int relative, char** o_buffer, size_t* o_size, size_t* o_used, size_t* o_offset, size_t* o_logicalUsed = 0, size_t* o_logicalOffset = NULL);

    virtual bool nextBatch();

    virtual void close();

private:
    struct Batch
    {
        char* buffer;
        AsyncFile::Writer* file;
        size_t used;
        size_t fileOffset;
        size_t logicalUsed;
        size_t logicalOffset;
    };
    Batch* batches;
    const int count;
    const size_t bufferSize;
    AsyncDataWriterSupplier* supplier;
    int current;
};

AsyncDataWriter::AsyncDataWriter(
    AsyncFile* i_file,
    AsyncDataWriterSupplier* i_supplier, 
    int i_count,
    size_t i_bufferSize,
    Filter* i_filter)
    :
    DataWriter(i_filter),
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
        batch->file->waitForCompletion();
    }
    return true;
}

    bool
AsyncDataWriter::nextBatch()
{
    _int64 start = timeInNanos();
    int written = current;
    Batch* write = &batches[written];
    write->logicalUsed = write->used;
    current = (current + 1) % count;
    batches[current].used = 0;
    bool newBuffer = filter != NULL && (filter->filterType == CopyFilter || filter->filterType == TransformFilter);
    bool newSize = newBuffer && filter->filterType == TransformFilter;
    if (newSize) {
        // advisory only
        write->fileOffset = supplier->sharedOffset;
        write->logicalOffset = supplier->sharedLogical;
    } else {
        supplier->advance(write->used, write->logicalUsed, &write->fileOffset, &write->logicalOffset);
    }
    if (filter != NULL) {
        size_t n = filter->onNextBatch(this, write->fileOffset, write->used);
        if (newSize) {
            write->used = n;
            supplier->advance(n, write->logicalUsed, &write->fileOffset, &write->logicalOffset);
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
    InterlockedAdd64AndReturnNewValue(&FilterTime, start2 - start);
    if (! write->file->beginWrite(write->buffer, write->used, write->fileOffset, NULL)) {
        fprintf(stderr, "error: file write %lld bytes at offset %lld failed\n", write->used, write->fileOffset);
        soft_exit(1);
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
    for (int i = 0; i < count; i++) {
        batches[i].file->close();
    }
}

AsyncDataWriterSupplier::AsyncDataWriterSupplier(
    const char* i_filename,
    DataWriter::FilterSupplier* i_filterSupplier,
    int i_bufferCount,
    size_t i_bufferSize)
    :
    filename(i_filename),
    filterSupplier(i_filterSupplier),
    bufferCount(i_bufferCount),
    bufferSize(i_bufferSize),
    sharedOffset(0),
    sharedLogical(0)
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
        filterSupplier ? filterSupplier->getFilter() : NULL);
}

    void
AsyncDataWriterSupplier::close()
{
    if (filterSupplier != NULL && filterSupplier->filterType == DataWriter::TransformFilter) {
        filterSupplier->onClose(this);
    }
    file->close();
    if (filterSupplier != NULL && filterSupplier->filterType != DataWriter::TransformFilter) {
        filterSupplier->onClose(this);
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
    ReleaseExclusiveLock(&lock);
}

    DataWriterSupplier*
DataWriterSupplier::create(
    const char* filename,
    DataWriter::FilterSupplier* filterSupplier,
    int count,
    size_t bufferSize)
{
    return new AsyncDataWriterSupplier(filename, filterSupplier, count, bufferSize);
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

    virtual void onClose(DataWriterSupplier* supplier)
    {
        a->onClose(supplier);
        b->onClose(supplier);
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
