/*++

Module Name:

    BufferedAsync.cpp

Abstract:

    Double-buffered asynchronous file I/O

Environment:

    User mode service.

    This class is NOT thread safe.  It's the caller's responsibility to ensure that
    at most one thread uses an instance at any time.

--*/

#include "stdafx.h"
#include "Compat.h"
#include "BigAlloc.h"
#include "BufferedAsync.h"

using std::min;


    bool
BufferedAsyncReader::open(
    AsyncFile* file,
    size_t offset,
    size_t bytes,
    size_t i_bufferSize,
    bool async,
    void* buffer0,
    void* buffer1)
{
    _ASSERT((buffer0 == NULL) == (buffer1 == NULL)); // both null or both non-null
    waitTime = 0;
    reader[0] = file->getReader();
    reader[1] = file->getReader();
    bufferSize = i_bufferSize;
    fileSize = offset + bytes;
    ownBuffer = buffer0 == NULL;
    buffer[0] = (char*) (buffer0 != NULL ? buffer0 : BigAlloc(bufferSize));
    buffer[1] = (char*) (buffer1 != NULL ? buffer1 : BigAlloc(bufferSize));
    if (reader[0] == NULL || reader[1] == NULL || buffer[0] == NULL || buffer[1] == NULL) {
        fprintf(stderr, "unable to setup temp file reader\n");
        return false;
    }
    length[0] = min(bytes, bufferSize);
    _int64 start = timeInNanos();
    reader[0]->beginRead(buffer[0],length[0], offset, NULL);
    if (bytes > bufferSize) {
        length[1] = min(bytes - bufferSize, bufferSize);
        reader[1]->beginRead(buffer[1], length[1], offset + bufferSize, NULL);
    } else {
        length[1] = 0;
    }
    waitTime += timeInNanos() - start;
    reading = 0;
    readOffset = 0;
    nextFileOffset = offset + min(bytes, 2 * bufferSize);
    if (! async) {
        endOpen();
    }
    return true;
}

    void
BufferedAsyncReader::endOpen()
{
    _int64 start = timeInNanos();
    reader[reading]->waitForCompletion();
    waitTime += timeInNanos() - start;
}

    bool
BufferedAsyncReader::atEnd()
{
    return length[reading] == 0;
}

    bool
BufferedAsyncReader::read(
    void* data,
    size_t bytes)
{
    if (bytes == 0) {
        return true;
    }
    if (atEnd()) {
        return false;
    }
    size_t first = min(bytes, length[reading] - readOffset);
    memcpy(data, buffer[reading] + readOffset, first);
    readOffset += bytes;
    if (readOffset >= length[reading]) {
        // switch buffers
        readOffset -= length[reading];
        reading = 1 - reading;
        _int64 start = timeInNanos();
        reader[reading]->waitForCompletion();
        waitTime += timeInNanos() - start;
        if (readOffset > 0) {
            // copy second part of read segment
            // todo: allow for longer reads
            if (readOffset > length[reading]) {
                fprintf(stderr, "read length too big\n");
                return false;
            }
            memcpy((char*) data + first, buffer[reading], readOffset);
        }
        // begin read of next block
        if (nextFileOffset < fileSize) {
            length[1 - reading] = min(fileSize - nextFileOffset, bufferSize);
            start = timeInNanos();
            reader[1 - reading]->beginRead(buffer[1 - reading], length[1 - reading], nextFileOffset, NULL);
            waitTime += timeInNanos() - start;
            nextFileOffset += length[1 - reading];
        } else {
            length[1 - reading] = 0;
        }
    }
    return true;
}

    bool
BufferedAsyncReader::close()
{
    delete reader[0];
    delete reader[1];
    if (ownBuffer) {
        BigDealloc(buffer[0]);
        BigDealloc(buffer[1]);
    }
    return true;
}
   
    bool
BufferedAsyncWriter::open(
    AsyncFile* file,
    size_t i_bufferSize,
    volatile _int64* sharedOffset)
{
    waitTime = 0;
    writer[0] = file->getWriter();
    writer[1] = file->getWriter();
    bufferSize = i_bufferSize;
    buffer[0] = (char*) BigAlloc(bufferSize);
    buffer[1] = (char*) BigAlloc(bufferSize);
    writeOffset = 0;
    writing = 0;
    privateFileOffset = 0;
    nextFileOffset = sharedOffset != NULL ? sharedOffset : &privateFileOffset;
    return writer[0] != NULL && writer[1] != NULL && buffer[0] != NULL && buffer[1] != NULL;
}

    bool
BufferedAsyncWriter::write(
    void* data,
    size_t bytes)
{
    void* p = forWrite(bytes);
    if (p == NULL) {
        // todo: allow bytes > bufferSize using synchronous writes
        return false;
    }
    memcpy(p, data, bytes);
    return true;
}
    
    void*
BufferedAsyncWriter::forWrite(
    size_t bytes)
{
    bool ok = true;
    if (writeOffset + bytes <= bufferSize) {
        writeOffset += bytes;
        return buffer[writing] + writeOffset - bytes;
    } else {
        _int64 fileOffset = InterlockedAdd64AndReturnNewValue(nextFileOffset, writeOffset) - writeOffset;
        _int64 start = timeInNanos();
        ok = writer[writing]->beginWrite(buffer[writing], writeOffset, fileOffset, NULL);
        ok &= writer[1 - writing]->waitForCompletion();
        waitTime += timeInNanos() - start;
        writing = 1 - writing;
        if (! ok) {
            fprintf(stderr, "BufferedAsyncWriter write failed\n");
            return NULL;
        }
        if (bytes <= bufferSize) {
            writeOffset = bytes;
            return buffer[writing];
        } else {
            // too big to be async
            fprintf(stderr, "BufferedAsyncWriter write too large\n");
            return NULL;
        }
    }
}

    bool
BufferedAsyncWriter::close()
{
    _int64 start = timeInNanos();
    bool ok = writer[1 - writing]->waitForCompletion();
    waitTime += timeInNanos() - start;
    if (writeOffset > 0) {
        _int64 fileOffset = InterlockedAdd64AndReturnNewValue(nextFileOffset, writeOffset) - writeOffset;
        _int64 start = timeInNanos();
        ok &= writer[writing]->beginWrite(buffer[writing], writeOffset, fileOffset, NULL);
        ok &= writer[writing]->waitForCompletion();
        waitTime += timeInNanos() - start;
    }
    delete writer[0];
    delete writer[1];
    BigDealloc(buffer[0]);
    BigDealloc(buffer[1]);
    return ok;
}
    