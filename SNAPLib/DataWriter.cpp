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

using std::min;

class AsyncDataWriter : public DataWriter
{
public:

    AsyncDataWriter(AsyncFile* i_file, volatile size_t* i_sharedOffset, int i_count, size_t i_bufferSize, Filter* i_filter);

    virtual ~AsyncDataWriter()
    {
        for (int i = 0; i < count; i++) {
            delete batches[i].file;
        }
        BigDealloc(batches[0].buffer); // all in one big block
        delete [] batches;
    }

    virtual bool getBuffer(char** o_buffer, size_t* o_size);

    virtual void advance(size_t bytes, unsigned location = 0);

    virtual bool getBatch(int relative, char** o_buffer, size_t* o_size, size_t* o_used);

    virtual bool nextBatch();

    virtual void close();

private:
    struct Batch
    {
        char* buffer;
        AsyncFile::Writer* file;
        size_t offset;
    };
    Batch* batches;
    const int count;
    const size_t bufferSize;
    volatile size_t* sharedOffset;
    int current;
};

AsyncDataWriter::AsyncDataWriter(
    AsyncFile* i_file,
    volatile size_t* i_sharedOffset, 
    int i_count,
    size_t i_bufferSize,
    Filter* i_filter)
    :
    DataWriter(i_filter),
    sharedOffset(i_sharedOffset),
    count(i_count),
    bufferSize(i_bufferSize),
    current(0)
{
    _ASSERT(count >= 2);
    char* block = (char*) BigAlloc(count * bufferSize);
    if (block == NULL) {
        fprintf(stderr, "Unable to allocate %lld bytes for write buffers\n", count * bufferSize);
        exit(1);
    }
    batches = new Batch[count];
    for (int i = 0; i < count; i++) {
        batches[i].buffer = block + i * bufferSize;
        batches[i].file = i_file->getWriter();
        batches[i].offset = 0;
    }
}
    
    bool
AsyncDataWriter::getBuffer(
    char** o_buffer,
    size_t* o_size)
{
    *o_buffer = batches[current].buffer + batches[current].offset;
    *o_size = bufferSize - batches[current].offset;
    return true;
}

    void
AsyncDataWriter::advance(
    size_t bytes,
    unsigned location)
{
    _ASSERT(bytes <= bufferSize - batches[current].offset);
    char* data = batches[current].buffer + batches[current].offset;
    size_t batchOffset = batches[current].offset;
    batches[current].offset = min(bufferSize, batchOffset + bytes);
    if (filter != NULL) {
        filter->onAdvance(this, batchOffset, data, bytes, location);
    }
}

    bool
AsyncDataWriter::getBatch(
    int relative,
    char** o_buffer,
    size_t* o_size,
    size_t* o_used)
{
    if (relative <= -count || relative >= count) {
        return false;
    }
    int index = (current + relative + count) % count; // ensure non-negative
    *o_buffer = batches[index].buffer;
    *o_size = bufferSize;
    *o_used = batches[index].offset;
    if (relative >= 0) {
        batches[index].file->waitForCompletion();
    }
    return true;
}

    bool
AsyncDataWriter::nextBatch()
{
    size_t bytes = batches[current].offset;
    current = (current + 1) % count;
    batches[current].offset = 0;
    size_t offset = 0;
    bool newSize = filter != NULL && filter->filterType == TransformFilter;
    bool newBuffer = filter != NULL && 
        (filter->filterType == CopyFilter || filter->filterType == TransformFilter);
    if (! newSize) {
        offset = (size_t) InterlockedAdd64AndReturnNewValue((_int64*)sharedOffset, bytes) - bytes;
    } else if (filter != NULL) {
        offset = *sharedOffset; // advisory only
    }
    if (filter != NULL) {
        size_t n = filter->onNextBatch(this, offset, bytes);
        if (newSize) {
            bytes = n;
            offset = (size_t) InterlockedAdd64AndReturnNewValue((_int64*)sharedOffset, bytes) - bytes;
        }
    }
    if (newBuffer) {
        current = (current + 1) % count;
    }
    int writeIndex = (current + count - 1) % count;
    if (! batches[writeIndex].file->beginWrite(batches[writeIndex].buffer, bytes, offset, NULL)) {
        return false;
    }
    if (! batches[current].file->waitForCompletion()) {
        return false;
    }
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

class AsyncDataWriterSupplier : public DataWriterSupplier
{
public:
    AsyncDataWriterSupplier(const char* i_filename, DataWriter::FilterSupplier* i_filterSupplier,
        int i_bufferCount, size_t i_bufferSize);

    virtual DataWriter* getWriter();

    virtual void close();

private:
    const char* filename;
    AsyncFile* file;
    DataWriter::FilterSupplier* filterSupplier;
    const int bufferCount;
    const size_t bufferSize;
    volatile size_t sharedOffset;
};

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
    sharedOffset(0)
{
    file = AsyncFile::open(filename, true);
    if (file == NULL) {
        fprintf(stderr, "failed to open %s for write\n", filename);
        exit(1);
    }
}

    DataWriter*
AsyncDataWriterSupplier::getWriter()
{
    return new AsyncDataWriter(file, &sharedOffset,
        bufferCount, bufferSize,
        filterSupplier ? filterSupplier->getFilter() : NULL);
}

    void
AsyncDataWriterSupplier::close()
{
    if (filterSupplier != NULL && filterSupplier->filterType == DataWriter::TransformFilter) {
        DataWriter* writer = new AsyncDataWriter(file, &sharedOffset, bufferCount, bufferSize, NULL);
        filterSupplier->onClose(this, writer);
        writer->close();
        delete writer;
    }
    file->close();
    if (filterSupplier != NULL && filterSupplier->filterType != DataWriter::TransformFilter) {
        filterSupplier->onClose(this, NULL);
    }
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
