/*++

Module Name:

    BufferedAsync.h

Abstract:

    Double-buffered asynchronous file I/O

Environment:

    User mode service.

    This class is NOT thread safe.  It's the caller's responsibility to ensure that
    at most one thread uses an instance at any time.

--*/

#pragma once

#include "stdafx.h"
#include "Compat.h"
 
class BufferedAsyncReader
{
public:
    bool                open(AsyncFile* file, size_t offset, size_t length, size_t bufferSize, bool async = false, void* buffer0 = NULL, void* buffer1 = NULL); 
    void                endOpen();
    bool                atEnd();
    bool                read(void* data, size_t bytes);
    bool                close();
    _int64              getWaitTimeInMillis() { return waitTime / 1000000; }

private:
    int                 reading; // 0 or 1
    size_t              readOffset;
    size_t              bufferSize;
    size_t              nextFileOffset;
    size_t              fileSize;
    AsyncFile::Reader*  reader[2];
    char*               buffer[2];
    size_t              length[2];
    bool                ownBuffer;
    _int64              waitTime; // in nanos
};

class BufferedAsyncWriter
{
public:
    bool                open(AsyncFile* file, size_t bufferSize, volatile _int64* sharedOffset = NULL); 
    bool                write(void* data, size_t bytes);
    void*               forWrite(size_t bytes);
    bool                close();
    _int64              getWaitTimeInMillis() { return waitTime / 1000000; }

private:
    int                 writing; // 0 or 1
    size_t              writeOffset; // within buffer
    size_t              bufferSize;
    volatile _int64     privateFileOffset; // used by default
    volatile _int64*    nextFileOffset; // for current buffer to write
    AsyncFile::Writer*  writer[2];
    char*               buffer[2];
    _int64              waitTime; // in nanos
};
