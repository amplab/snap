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

class DataWriterSupplier;

// per-thread writer for data into a single destination
class DataWriter
{
public:
    
    class Watcher
    {
    public:
        virtual ~Watcher() {}

        // called when a chunk of data has been written into the file
        virtual void onAdvance(DataWriter* writer, char* data, size_t bytes, unsigned location) = 0;

        // called when a batch has been completed, after advancing to the next
        // e.g. so use getBatch(-1, ...) to get the one that was just written
        virtual void onNextBatch(DataWriter* writer, size_t offset, size_t bytes) = 0;
    };
    
    class WatcherSupplier
    {
    public:
        virtual ~WatcherSupplier();

        virtual Watcher* getWatcher() = 0;

        // called when entire file is done, all Watchers destroyed
        virtual void onClose(DataWriterSupplier* supplier) = 0;
    };

    DataWriter(Watcher* i_watcher) : watcher(i_watcher) {}

    virtual ~DataWriter() {}

    // get remaining space in current buffer for writing
    virtual bool getBuffer(char** o_buffer, size_t* o_size) = 0;

    // advance within current buffer, reducing available space
    // should be called on each read, with the location
    virtual void advance(size_t bytes, unsigned location = 0) = 0;

    // get complete data buffer in batch, relative==0 is current, relative==-1 is previous, etc.
    virtual bool getBatch(int relative, char** o_buffer, size_t* o_size, size_t* o_used) = 0;

    // advance to next buffer
    virtual bool nextBatch() = 0;

    // this thread is complete
    virtual void close() = 0;

protected:
    Watcher* watcher;
};

// creates writers for multiple threads
class DataWriterSupplier
{
public:
    virtual ~DataWriterSupplier() {}

    virtual DataWriter* getWriter() = 0;

    // call when all threads are done, all watchers destroyed
    virtual void close() = 0;
    
    static DataWriterSupplier* create(const char* filename,
        DataWriter::WatcherSupplier* watcherSupplier = NULL,
        int count = 2, size_t bufferSize = 16 * 1024 * 1024);
};
