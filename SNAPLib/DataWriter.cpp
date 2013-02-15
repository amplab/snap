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

    AsyncDataWriter(AsyncFile* i_file, volatile size_t* i_sharedOffset, int i_count, size_t i_bufferSize, Watcher* i_watcher);

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
    Watcher* i_watcher)
    :
    DataWriter(i_watcher),
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
    batches[current].offset = min(bufferSize, batches[current].offset + bytes);
    if (watcher != NULL) {
        watcher->onAdvance(this, data, bytes, location);
    }
}

    bool
AsyncDataWriter::getBatch(
    int relative,
    char** o_buffer,
    size_t* o_size,
    size_t* o_used)
{
    int index = (((current + relative) % count) + count) % count; // ensure non-negative
    *o_buffer = batches[index].buffer;
    *o_size = bufferSize;
    *o_used = batches[index].offset;
    return true;
}

    bool
AsyncDataWriter::nextBatch()
{
    size_t bytes = batches[current].offset;
    size_t offset = (size_t) InterlockedAdd64AndReturnNewValue((_int64*)sharedOffset, bytes) - bytes;
    if (bytes > 0) {
        if (! batches[current].file->beginWrite(batches[current].buffer, bytes, offset, NULL)) {
            return false;
        }
    }
    current = (current + 1) % count;
    if (! batches[current].file->waitForCompletion()) {
        return false;
    }
    if (watcher != NULL) {
        watcher->onNextBatch(this, offset, bytes);
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
    AsyncDataWriterSupplier(const char* i_filename, DataWriter::WatcherSupplier* i_watcherSupplier,
        int i_bufferCount, size_t i_bufferSize);

    virtual DataWriter* getWriter();

    virtual void close();

private:
    const char* filename;
    AsyncFile* file;
    DataWriter::WatcherSupplier* watcherSupplier;
    const int bufferCount;
    const size_t bufferSize;
    volatile size_t sharedOffset;
};

AsyncDataWriterSupplier::AsyncDataWriterSupplier(
    const char* i_filename,
    DataWriter::WatcherSupplier* i_watcherSupplier,
    int i_bufferCount,
    size_t i_bufferSize)
    :
    filename(i_filename),
    watcherSupplier(i_watcherSupplier),
    bufferCount(i_bufferCount),
    bufferSize(i_bufferSize),
    sharedOffset(0)
{
    file = AsyncFile::open(filename, true);
}

    DataWriter*
AsyncDataWriterSupplier::getWriter()
{
    return new AsyncDataWriter(file, &sharedOffset,
        bufferCount, bufferSize,
        watcherSupplier ? watcherSupplier->getWatcher() : NULL);
}

    void
AsyncDataWriterSupplier::close()
{
    file->close();
    if (watcherSupplier != NULL) {
        watcherSupplier->onClose(this);
    }
}

    DataWriterSupplier*
DataWriterSupplier::create(
    const char* filename,
    DataWriter::WatcherSupplier* watcherSupplier,
    int count,
    size_t bufferSize)
{
    return new AsyncDataWriterSupplier(filename, watcherSupplier, count, bufferSize);
}