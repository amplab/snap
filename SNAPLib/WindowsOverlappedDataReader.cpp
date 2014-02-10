/*++

Module Name:

    WindowsOverlappedDataReader.cpp

Abstract:

        Concrete implementation classes for WindowsOverlappedDataReader.

        These are completely opaque, and are only exposed through static supplier objects
        defined in DataReader.h

Environment:

    User mode service.

--*/

#include "stdafx.h"
#include "BigAlloc.h"
#include "Compat.h"
#include "RangeSplitter.h"
#include "ParallelTask.h"
#include "DataReader.h"
#include "Bam.h"
#include "zlib.h"
#include "exit.h"

#if 0   // This shoudln't really be here at all, it's just kept around as a reference as I refactor.  --bb


using std::max;
using std::min;
using std::map;
using std::string;

#ifdef _MSC_VER

//
// WindowsOverlapped
//

class WindowsOverlappedDataReader : public DataReader
{
public:

    WindowsOverlappedDataReader(unsigned i_nBuffers, _int64 i_overflowBytes, double extraFactor, bool autoRelease);

    virtual ~WindowsOverlappedDataReader();
    
    virtual bool init(const char* fileName);

    virtual char* readHeader(_int64* io_headerSize);

    virtual void reinit(_int64 startingOffset, _int64 amountOfFileToProcess);

    virtual bool getData(char** o_buffer, _int64* o_validBytes, _int64* o_startBytes = NULL);

    virtual void advance(_int64 bytes);

    virtual void nextBatch();

    virtual bool isEOF();

    virtual DataBatch getBatch();

    virtual void releaseBatch(DataBatch batch);

    virtual _int64 getFileOffset();

    virtual void getExtra(char** o_extra, _int64* o_length);

private:
    
    // must hold the lock to call
    void startIo();

    // must hold the lock to call
    void waitForBuffer(unsigned bufferNumber);

    // must hold the lock to call
    void addBuffer();

    const char*         fileName;
    LARGE_INTEGER       fileSize;
    HANDLE              hFile;
  
    static const unsigned bufferSize = 4 * 1024 * 1024 - 4096;

    enum BufferState {Empty, Reading, Full, InUse};

    struct BufferInfo
    {
        char            *buffer;
        BufferState     state;
        DWORD           validBytes;
        DWORD           nBytesThatMayBeginARead;
        bool            isEOF;
        unsigned        offset;     // How far has the consumer gotten in current buffer
        _int64          fileOffset;
        _uint32         batchID;
        OVERLAPPED      lap;
        char*           extra;
        int             next, previous; // index of next/previous in free/ready list, -1 if end
    };

    unsigned            nBuffers;
    const unsigned      maxBuffers;
    _int64              extraBytes;
    _int64              overflowBytes;
    BufferInfo*         bufferInfo;
    LARGE_INTEGER       readOffset;
    _int64              endingOffset;
    _uint32             nextBatchID;
    int                 nextBufferForReader; // list head (singly linked), -1 if empty
    int                 nextBufferForConsumer; // list head (doubly linked), -1 if empty
    int                 lastBufferForConsumer; // list tail, -1 if empty
    HANDLE              releaseEvent;
    DWORD               releaseWait;
    ExclusiveLock       lock;
};

WindowsOverlappedDataReader::WindowsOverlappedDataReader(
    unsigned i_nBuffers,
    _int64 i_overflowBytes,
    double extraFactor,
    bool autoRelease)
    : DataReader(autoRelease), nBuffers(i_nBuffers), overflowBytes(i_overflowBytes), maxBuffers(4 * i_nBuffers)
{
    //
    // Initialize the buffer info struct.
    //
    
    // allocate all the data in one big block
    // NOTE: buffers are not null-terminated (since memmap version can't do it)
    _ASSERT(extraFactor >= 0 && i_nBuffers > 0);
    bufferInfo = new BufferInfo[maxBuffers];
    extraBytes = max(0LL, (_int64) ((bufferSize + overflowBytes) * extraFactor));
    char* allocated = (char*) BigReserve(maxBuffers * (bufferSize + extraBytes + overflowBytes));
    if (NULL == allocated) {
        fprintf(stderr,"WindowsOverlappedDataReader: unable to allocate IO buffer\n");
        soft_exit(1);
    }
    if (! BigCommit(allocated, nBuffers * (bufferSize + extraBytes + overflowBytes))) {
        fprintf(stderr, "WindowsOverlappedDataReader: ubable to commit IO  buffer\n");
        soft_exit(1);
    }
    for (unsigned i = 0 ; i < nBuffers; i++) {
        bufferInfo[i].buffer = allocated;
        allocated += bufferSize + overflowBytes;
        bufferInfo[i].extra = extraBytes > 0 ? allocated : NULL;
        allocated += extraBytes;

        bufferInfo[i].lap.hEvent = CreateEvent(NULL,TRUE,FALSE,NULL);
        if (NULL == bufferInfo[i].lap.hEvent) {
            fprintf(stderr,"WindowsOverlappedDataReader: Unable to create event\n");
            soft_exit(1);
        }

        bufferInfo[i].state = Empty;
        bufferInfo[i].isEOF = false;
        bufferInfo[i].offset = 0;
        bufferInfo[i].next = i < nBuffers - 1 ? i + 1 : -1;
        bufferInfo[i].previous = i > 0 ? i - 1 : -1;
    }
    nextBatchID = 1;
    hFile = INVALID_HANDLE_VALUE;
    nextBufferForConsumer = -1;
    lastBufferForConsumer = -1;
    nextBufferForReader = 0;
    releaseEvent = CreateEvent(NULL, TRUE, TRUE, NULL);
    releaseWait = 5; // wait up to 5 ms before allocating a new buffer
    InitializeExclusiveLock(&lock);
}

WindowsOverlappedDataReader::~WindowsOverlappedDataReader()
{
    BigDealloc(bufferInfo[0].buffer);
    for (unsigned i = 0; i < nBuffers; i++) {
        bufferInfo[i].buffer = bufferInfo[i].extra = NULL;
        CloseHandle(bufferInfo[i].lap.hEvent);
    }
    CloseHandle(hFile);
    hFile = INVALID_HANDLE_VALUE;
    DestroyExclusiveLock(&lock);
}
    
    bool
WindowsOverlappedDataReader::init(
    const char* i_fileName)
{
    fileName = i_fileName;
    hFile = CreateFile(fileName,GENERIC_READ,FILE_SHARE_READ,NULL,OPEN_EXISTING,FILE_FLAG_OVERLAPPED,NULL);
    if (INVALID_HANDLE_VALUE == hFile) {
        return false;
    }

    if (!GetFileSizeEx(hFile,&fileSize)) {
        fprintf(stderr,"WindowsOverlappedDataReader: unable to get file size of '%s', %d\n",fileName,GetLastError());
        return false;
    }
    return true;
}

    char*
WindowsOverlappedDataReader::readHeader(
    _int64* io_headerSize)
{
    BufferInfo *info = &bufferInfo[0];
    info->fileOffset = 0;
    info->offset = 0;
    info->lap.Offset = 0;
    info->lap.OffsetHigh = 0;
    _ASSERT(nextBufferForReader == 0 && nextBufferForConsumer == -1 && lastBufferForConsumer == -1 && info->next == 1 && info->previous == -1);
    nextBufferForReader = 1;
    nextBufferForConsumer = lastBufferForConsumer = 0;
    info->next = info->previous = -1;

    if (*io_headerSize > 0xffffffff) {
        fprintf(stderr,"WindowsOverlappedDataReader: trying to read too many bytes at once: %lld\n", *io_headerSize);
        soft_exit(1);
    }

    if (!ReadFile(hFile,info->buffer,(DWORD)*io_headerSize,&info->validBytes,&info->lap)) {
        if (GetLastError() != ERROR_IO_PENDING) {
            fprintf(stderr,"WindowsOverlappedSAMReader::init: unable to read header of '%s', %d\n",fileName,GetLastError());
            return false;
        }
    }

    if (!GetOverlappedResult(hFile,&info->lap,&info->validBytes,TRUE)) {
        fprintf(stderr,"WindowsOverlappedSAMReader::init: error reading header of '%s', %d\n",fileName,GetLastError());
        return false;
    }

    *io_headerSize = info->validBytes;
    return info->buffer;
}

    void
WindowsOverlappedDataReader::reinit(
    _int64 i_startingOffset,
    _int64 amountOfFileToProcess)
{
    _ASSERT(INVALID_HANDLE_VALUE != hFile);  // Must call init() before reinit()

    AcquireExclusiveLock(&lock);

    //
    // First let any pending IO complete.
    //
    for (unsigned i = 0; i < nBuffers; i++) {
        if (bufferInfo[i].state == Reading) {
            waitForBuffer(i);
        }
        bufferInfo[i].state = Empty;
        bufferInfo[i].isEOF= false;
        bufferInfo[i].offset = 0;
        bufferInfo[i].next = i < nBuffers - 1 ? i + 1 : -1;
        bufferInfo[i].previous = i > 0 ? i - 1 : -1;
    }

    nextBufferForConsumer = -1; 
    lastBufferForConsumer = -1;
    nextBufferForReader = 0;

    readOffset.QuadPart = i_startingOffset;
    if (amountOfFileToProcess == 0) {
        //
        // This means just read the whole file.
        //
        endingOffset = fileSize.QuadPart;
    } else {
        endingOffset = min(fileSize.QuadPart,i_startingOffset + amountOfFileToProcess);
    }

    //
    // Kick off IO, wait for the first buffer to be read
    //
    startIo();
    waitForBuffer(nextBufferForConsumer);

    ReleaseExclusiveLock(&lock);
}

    bool
WindowsOverlappedDataReader::getData(
    char** o_buffer,
    _int64* o_validBytes,
    _int64* o_startBytes)
{
    _ASSERT(nextBufferForConsumer >= 0 && nextBufferForConsumer < (int) nBuffers);
    BufferInfo *info = &bufferInfo[nextBufferForConsumer];
    if (info->isEOF && info->offset >= info->validBytes) {
        //
        // EOF.
        //
        return false;
    }

    if (info->offset >= info->nBytesThatMayBeginARead) {
        //
        // Past the end of our section.
        //
        return false;
    }

    if (info->state != Full) {
        _ASSERT(info->state != InUse);
        AcquireExclusiveLock(&lock);
        waitForBuffer(nextBufferForConsumer);
        ReleaseExclusiveLock(&lock);
    }
    
    *o_buffer = info->buffer + info->offset;
    *o_validBytes = info->validBytes - info->offset;
    if (o_startBytes != NULL) {
        *o_startBytes = info->nBytesThatMayBeginARead - info->offset;
    }
    return true;
}

    void
WindowsOverlappedDataReader::advance(
    _int64 bytes)
{
    BufferInfo* info = &bufferInfo[nextBufferForConsumer];
    _ASSERT(info->validBytes >= info->offset && bytes >= 0 && bytes <= info->validBytes - info->offset);
    info->offset += min(info->validBytes - info->offset, (unsigned)max(0, bytes));
}

    void
WindowsOverlappedDataReader::nextBatch()
{
    AcquireExclusiveLock(&lock);
    _ASSERT(nextBufferForConsumer >= 0 && nextBufferForConsumer < (int) nBuffers);
    BufferInfo* info = &bufferInfo[nextBufferForConsumer];
    if (info->isEOF) {
        ReleaseExclusiveLock(&lock);
        if (autoRelease) {
            releaseBatch(DataBatch(info->batchID));
        }
        return;
    }
    DataBatch priorBatch = DataBatch(info->batchID);

    info->state = InUse;
    _uint32 overflow = max((DWORD) info->offset, info->nBytesThatMayBeginARead) - info->nBytesThatMayBeginARead;
    _int64 nextStart = info->fileOffset + info->nBytesThatMayBeginARead; // for validation

    nextBufferForConsumer = info->next;

    bool first = true;
    while (nextBufferForConsumer == -1) {
        nextStart = 0; // can no longer count on getting sequential buffers from file
        ReleaseExclusiveLock(&lock);
        if (! first) {
            //printf("WindowsOverlappedDataReader::nextBatch thread %d wait for release\n", GetCurrentThreadId());
            _int64 start = timeInNanos();
            DWORD result = WaitForSingleObject(releaseEvent, releaseWait);
            InterlockedAdd64AndReturnNewValue(&ReleaseWaitTime, timeInNanos() - start);
            //printf("WindowsOverlappedDataReader::nextBatch thread %d released\n", GetCurrentThreadId());
            if (result == WAIT_TIMEOUT) {
                AcquireExclusiveLock(&lock);
                addBuffer();
                ReleaseExclusiveLock(&lock);
            }
        }
        first = false;
        startIo();
        AcquireExclusiveLock(&lock);
    }
    if (bufferInfo[nextBufferForConsumer].state != Full) {
        waitForBuffer(nextBufferForConsumer);
    }
    bufferInfo[nextBufferForConsumer].offset = overflow;
    _ASSERT(nextStart == 0 || nextStart == bufferInfo[nextBufferForConsumer].fileOffset);

    ReleaseExclusiveLock(&lock);

    if (autoRelease) {
        releaseBatch(priorBatch);
    }
}

    bool
WindowsOverlappedDataReader::isEOF()
{
    return bufferInfo[nextBufferForConsumer].isEOF;
}
    
    DataBatch
WindowsOverlappedDataReader::getBatch()
{
    return DataBatch(bufferInfo[nextBufferForConsumer].batchID);
}

    void
WindowsOverlappedDataReader::releaseBatch(
    DataBatch batch)
{
    AcquireExclusiveLock(&lock);

    bool released = false;
    for (unsigned i = 0; i < nBuffers; i++) {
        BufferInfo* info = &bufferInfo[i];
        if (info->batchID == batch.batchID) {
            switch (info->state) {
            case Empty:
                // nothing
                break;

            case Reading:
                // todo: cancel read operation?
                _ASSERT(false);
                break;

            case InUse:
                released = true;
                // fall through
            case Full:
                //printf("releaseBatch batch %d, releasing %s buffer %d\n", batch.batchID, info->state == InUse ? "InUse" : "Full", i);
                info->state = Empty;
                // remove from ready list
                if (i == nextBufferForConsumer) {
                    nextBufferForConsumer = info->next;
                }
                if (i == lastBufferForConsumer) {
                    lastBufferForConsumer = info->previous;
                }
                if (info->next != -1) {
                    bufferInfo[info->next].previous = info->previous;
                }
                if (info->previous != -1) {
                    bufferInfo[info->previous].next = info->next;
                }
                // add to head of free list
                info->next = nextBufferForReader;
                nextBufferForReader = i;
                break;

            default:
                fprintf(stderr, "invalid enum\n");
                soft_exit(1);
            }
        }
    }

    startIo();

    if (released) {
        //printf("releaseBatch set releaseEvent\n");
        SetEvent(releaseEvent);
    }

    ReleaseExclusiveLock(&lock);
}

    _int64
WindowsOverlappedDataReader::getFileOffset()
{
    return bufferInfo[nextBufferForConsumer].fileOffset + bufferInfo[nextBufferForConsumer].offset;
}

    void
WindowsOverlappedDataReader::getExtra(
    char** o_extra,
    _int64* o_length)
{
    // hack: return valid buffer even when no consumer buffers - this may happen when reading header
    *o_extra = bufferInfo[max(0,nextBufferForConsumer)].extra;
    *o_length = extraBytes;
}
    
    void
WindowsOverlappedDataReader::startIo()
{
    //
    // Launch reads on whatever buffers are ready.
    //
    AcquireExclusiveLock(&lock);
    while (nextBufferForReader != -1) {
        // remove from free list
        BufferInfo* info = &bufferInfo[nextBufferForReader];
        _ASSERT(info->state == Empty);
        int index = nextBufferForReader;
        nextBufferForReader = info->next;
        info->batchID = nextBatchID++;
        // add to end of consumer list
        if (lastBufferForConsumer != -1) {
            _ASSERT(bufferInfo[lastBufferForConsumer].next == -1);
            bufferInfo[lastBufferForConsumer].next = index;
        }
        info->next = -1;
        info->previous = lastBufferForConsumer;
        lastBufferForConsumer = index;
		if (nextBufferForConsumer == -1) {
			nextBufferForConsumer = index;
		}

        if (readOffset.QuadPart >= fileSize.QuadPart || readOffset.QuadPart >= endingOffset) {
            info->validBytes = 0;
            info->nBytesThatMayBeginARead = 0;
            info->isEOF = true;
            info->state = Full;
            SetEvent(info->lap.hEvent);
            ReleaseExclusiveLock(&lock);
            return;
        }

        unsigned amountToRead;
        _int64 finalOffset = min(fileSize.QuadPart, endingOffset + overflowBytes);
        _int64 finalStartOffset = min(fileSize.QuadPart, endingOffset);
        amountToRead = (unsigned)min(finalOffset - readOffset.QuadPart, (_int64) bufferSize);   // Cast OK because can't be longer than unsigned bufferSize
        info->isEOF = readOffset.QuadPart + amountToRead == finalOffset;
        info->nBytesThatMayBeginARead = (unsigned)min(bufferSize - overflowBytes, finalStartOffset - readOffset.QuadPart);

        _ASSERT(amountToRead >= info->nBytesThatMayBeginARead && (!info->isEOF || finalOffset == readOffset.QuadPart + amountToRead));
        ResetEvent(info->lap.hEvent);
        info->lap.Offset = readOffset.LowPart;
        info->lap.OffsetHigh = readOffset.HighPart;
        info->fileOffset = readOffset.QuadPart;

        readOffset.QuadPart += info->nBytesThatMayBeginARead;
        info->state = Reading;
        info->offset = 0;
         
        //printf("startIo on %d at %lld for %uB\n", index, readOffset, amountToRead);
        ReleaseExclusiveLock(&lock);
        if (!ReadFile(
                hFile,
                info->buffer,
                amountToRead,
                &info->validBytes,
                &info->lap)) {

            if (GetLastError() != ERROR_IO_PENDING) {
                fprintf(stderr,"FASTQReader::startIo(): readFile failed, %d\n",GetLastError());
                soft_exit(1);
            }
        }
        AcquireExclusiveLock(&lock);
    }
    if (nextBufferForConsumer == -1) {
        if (nextBufferForConsumer == -1) {
            //printf("startIo thread %x reset releaseEvent\n", GetCurrentThreadId());
            ResetEvent(releaseEvent);
        }
    }
    ReleaseExclusiveLock(&lock);
}

    void
WindowsOverlappedDataReader::waitForBuffer(
    unsigned bufferNumber)
{
    _ASSERT(bufferNumber >= 0 && bufferNumber < nBuffers);
    BufferInfo *info = &bufferInfo[bufferNumber];

    while (info->state == InUse) {
        //printf("WindowsOverlappedDataReader::waitForBuffer %d InUse...\n", bufferNumber);
        // must already have lock to call, release & wait & reacquire
        ReleaseExclusiveLock(&lock);
        _int64 start = timeInNanos();
        DWORD result = WaitForSingleObject(releaseEvent, releaseWait);
        InterlockedAdd64AndReturnNewValue(&ReleaseWaitTime, timeInNanos() - start);
        AcquireExclusiveLock(&lock);
        if (result == WAIT_TIMEOUT) {
            // this isn't going to directly make this buffer available, but will reduce pressure
            addBuffer();
        }
    }

    if (info->state == Full) {
        return;
    }

    if (info->state != Reading) {
        startIo();
    }

    _int64 start = timeInNanos();
    if (!GetOverlappedResult(hFile,&info->lap,&info->validBytes,TRUE)) {
        fprintf(stderr,"Error reading FASTQ file, %d\n",GetLastError());
        soft_exit(1);
    }
    InterlockedAdd64AndReturnNewValue(&ReadWaitTime, timeInNanos() - start);

    info->state = Full;
    info->buffer[info->validBytes] = 0;
    ResetEvent(info->lap.hEvent);
}

    void
WindowsOverlappedDataReader::addBuffer()
{
    if (nBuffers == maxBuffers) {
        //printf("WindowsOverlappedDataReader: addBuffer at limit\n");
        return;
    }
    _ASSERT(nBuffers < maxBuffers);
    //printf("WindowsOverlappedDataReader: addBuffer %d of %d\n", nBuffers, maxBuffers);
    size_t bytes = bufferSize + extraBytes + overflowBytes;
    bufferInfo[nBuffers].buffer = bufferInfo[nBuffers-1].buffer + bytes;
    if (! BigCommit(bufferInfo[nBuffers].buffer, bytes)) {
        fprintf(stderr, "WindowsOverlappedDataReader: unable to commit IO buffer\n");
        soft_exit(1);
    }
    bufferInfo[nBuffers].extra = extraBytes > 0 ? bufferInfo[nBuffers].buffer + bytes - extraBytes : NULL;

    bufferInfo[nBuffers].lap.hEvent = CreateEvent(NULL,TRUE,FALSE,NULL);
    if (NULL == bufferInfo[nBuffers].lap.hEvent) {
        fprintf(stderr,"WindowsOverlappedDataReader: Unable to create event\n");
        soft_exit(1);
    }

    bufferInfo[nBuffers].state = Empty;
    bufferInfo[nBuffers].isEOF= false;
    bufferInfo[nBuffers].offset = 0;
    bufferInfo[nBuffers].next = nextBufferForReader;
    bufferInfo[nBuffers].previous = -1;
    nextBufferForReader = nBuffers;
    nBuffers++;
    _ASSERT(nBuffers <= maxBuffers);
    if (nBuffers == maxBuffers) {
        releaseWait = INFINITE;
    }
}
    
class WindowsOverlappedDataSupplier : public DataSupplier
{
public:
    WindowsOverlappedDataSupplier(bool autoRelease) : DataSupplier(autoRelease) {}
    virtual DataReader* getDataReader(_int64 overflowBytes, double extraFactor = 0.0)
    {
        int buffers = autoRelease ? 2 : (ThreadCount + max(ThreadCount * 3 / 4, 3));
        return new WindowsOverlappedDataReader(buffers, overflowBytes, extraFactor, autoRelease);
    }
};

DataSupplier* DataSupplier::WindowsOverlapped[2] =
{ new WindowsOverlappedDataSupplier(false), new WindowsOverlappedDataSupplier(true) };

#endif // _MSC_VER


#endif // 0

