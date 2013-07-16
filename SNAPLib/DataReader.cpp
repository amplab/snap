/*++

Module Name:

    DataReader.cpp

Abstract:

    Concrete implementation classes for DataReader and DataSupplier.

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
#include "zlib.h"
#include "exit.h"

using std::max;
using std::min;

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
    virtual DataReader* getDataReader(_int64 overflowBytes, double extraFactor = 0.0) const
    {
        int buffers = autoRelease ? 2 : (ThreadCount + max(ThreadCount * 3 / 4, 3));
        return new WindowsOverlappedDataReader(buffers, overflowBytes, extraFactor, autoRelease);
    }
};

const DataSupplier* DataSupplier::WindowsOverlapped[2] =
{ new WindowsOverlappedDataSupplier(false), new WindowsOverlappedDataSupplier(true) };

#endif // _MSC_VER

//
// Gzip
//

static const int windowBits = 15;
static const int ENABLE_ZLIB_GZIP = 32;

static const double MIN_FACTOR = 1.2;
static const double MAX_FACTOR = 10.0;

class GzipDataReader : public DataReader
{
public:

    GzipDataReader(_int64 i_overflowBytes, _int64 i_extraBytes, DataReader* i_inner, bool autoRelease);

    virtual ~GzipDataReader();
    
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

    enum DecompressMode { SingleBlock, ContinueMultiBlock, StartMultiBlock };
    bool decompress(char* input, _int64 inputBytes, _int64* o_inputRead,
        char* output, _int64 outputBytes, _int64* o_outputWritten,
        DecompressMode mode);

    void decompressBatch();
    
    const _int64    extraBytes; // extra bytes I have to use
    const _int64    overflowBytes; // overflow bytes I may need to allow between chunks
    DataReader*     inner; // inner reader to get data chunks
    z_stream        zstream; // stateful zlib stream
    _int64          startBytes; // start bytes in extra data
    _int64          validBytes; // valid bytes in extra data
    _int64          offset; // current offset in extra data
    _int64          priorBytes; // bytes copied from end of prior buffer
    bool            gotBatchData; // whether data has been read & decompressed for current batch
    bool            continueBlock; // whether to continue a block or start a new one
    ThreadHeap      heap;
};
 
GzipDataReader::GzipDataReader(_int64 i_overflowBytes, _int64 i_extraBytes, DataReader* i_inner, bool autoRelease)
    : DataReader(autoRelease), overflowBytes(i_overflowBytes),
    extraBytes(i_extraBytes),
    inner(i_inner),
    zstream(),
    heap(65536)
{
}

GzipDataReader::~GzipDataReader()
{
    delete inner;
}

    bool
GzipDataReader::init(
    const char* i_fileName)
{
    if (! inner->init(i_fileName)) {
        return false;
    }
    zstream.zalloc = zalloc;
    zstream.zfree = zfree;
    zstream.opaque = &heap;
    return true;
}

    char*
GzipDataReader::readHeader(
    _int64* io_headerSize)
{
    _int64 compressedBytes = (_int64)(*io_headerSize / MIN_FACTOR);
    char* compressed = inner->readHeader(&compressedBytes);
    char* header;
    _int64 extraBytes;
    inner->getExtra(&header, &extraBytes);
    _int64 compressedHeaderSize;
    decompress(compressed, compressedBytes, &compressedHeaderSize, header, *io_headerSize, io_headerSize, SingleBlock);
    inner->advance(compressedHeaderSize);
    return header;
}

    void
GzipDataReader::reinit(
    _int64 i_startingOffset,
    _int64 amountOfFileToProcess)
{
    // todo: transform start/amount to add for compression? I don't think so...
    inner->reinit(i_startingOffset, amountOfFileToProcess);
    gotBatchData = false;
    continueBlock = false;
    offset = 0;
    priorBytes = 0;
}

    bool
GzipDataReader::getData(
    char** o_buffer,
    _int64* o_validBytes,
    _int64* o_startBytes)
{
    if (! gotBatchData) {
        decompressBatch();
        gotBatchData = true;
    }
    if (offset >= startBytes) {
        return false;
    }
    inner->getExtra(o_buffer, o_validBytes);
    *o_buffer += offset;
    *o_validBytes = validBytes - offset;
    if (o_startBytes) {
        *o_startBytes = startBytes - offset;
    }
    return true;
}

    void
GzipDataReader::advance(
    _int64 bytes)
{
    offset = min(offset + max((_int64) 0, bytes), validBytes);
}

    void
GzipDataReader::nextBatch()
{
    // copy trailing data off the end of prior buffer
    priorBytes = validBytes - max(offset, startBytes);
    char* priorData; _int64 n;
    inner->getExtra(&priorData, &n);
    priorData += validBytes - priorBytes;
    DataBatch innerBatch = inner->getBatch();
    inner->nextBatch();
    char* currentData;
    inner->getExtra(&currentData, &n);
    memcpy(currentData, priorData, priorBytes);
    if (autoRelease) {
        inner->releaseBatch(innerBatch);
    }
    offset = 0;
    gotBatchData = false;
    continueBlock = true;
}

    bool
GzipDataReader::isEOF()
{
    return inner->isEOF();
}
    
    DataBatch
GzipDataReader::getBatch()
{
    return inner->getBatch();
}

    void
GzipDataReader::releaseBatch(
    DataBatch batch)
{
    inner->releaseBatch(batch);
}

    _int64
GzipDataReader::getFileOffset()
{
    return inner->getFileOffset();
}

    void
GzipDataReader::getExtra(
    char** o_extra,
    _int64* o_length)
{
    inner->getExtra(o_extra, o_length);
    *o_extra += extraBytes;
    *o_length -= extraBytes;
}
    
    bool
GzipDataReader::decompress(
    char* input,
    _int64 inputBytes,
    _int64* o_inputRead,
    char* output,
    _int64 outputBytes,
    _int64* o_outputWritten,
    DecompressMode mode)
{
    if (inputBytes > 0xffffffff || outputBytes > 0xffffffff) {
        fprintf(stderr,"GzipDataReader: inputBytes or outputBytes > max unsigned int\n");
        soft_exit(1);
    }
    zstream.next_in = (Bytef*) input;
    zstream.avail_in = (uInt)inputBytes;
    zstream.next_out = (Bytef*) output;
    zstream.avail_out = (uInt)outputBytes;
    uInt oldAvailOut, oldAvailIn;
    int block = 0;
    int status;
    do {
	    if (mode != ContinueMultiBlock || block != 0) {
            heap.reset();
            status = inflateInit2(&zstream, windowBits | ENABLE_ZLIB_GZIP);
            if (status < 0) {
                fprintf(stderr, "GzipDataReader: inflateInit2 failed with %d\n", status);
                return false;
            }
        }
        oldAvailOut = zstream.avail_out;
        oldAvailIn = zstream.avail_out;
        status = inflate(&zstream, mode == SingleBlock ? Z_NO_FLUSH : Z_FINISH);
        //printf("decompress block #%d %lld -> %lld = %d\n", block, zstream.next_in - lastIn, zstream.next_out - lastOut, status);
        block++;
        if (status < 0 && status != Z_BUF_ERROR) {
            fprintf(stderr, "GzipDataReader: inflate failed with %d\n", status);
            soft_exit(1);
        }
        if (status < 0 && zstream.avail_out == 0 && zstream.avail_in > 0) {
            fprintf(stderr, "GzipDataReader: insufficient decompression buffer space\n");
            soft_exit(1);
        }
    } while (zstream.avail_in != 0 && (zstream.avail_out != oldAvailOut || zstream.avail_in != oldAvailIn) && mode != SingleBlock);
    //printf("end decompress status=%d, avail_in=%lld, last block=%lld->%lld, avail_out=%lld\n", status, zstream.avail_in, zstream.next_in - lastIn, zstream.next_out - lastOut, zstream.avail_out);
    if (o_inputRead) {
        *o_inputRead = inputBytes - zstream.avail_in;
    }
    if (o_outputWritten) {
        *o_outputWritten = outputBytes - zstream.avail_out;
    }
    return zstream.avail_in == 0;
}
    
    void
GzipDataReader::decompressBatch()
{
    size_t fileOffset = inner->getFileOffset();
    char* compressed;
    _int64 compressedBytes;
    if (! inner->getData(&compressed, &compressedBytes)) {
	    inner->nextBatch();
	    if (! inner->getData(&compressed, &compressedBytes)) {
            if (inner->isEOF()) {
                offset = 0;
                validBytes = startBytes = 0;
                return;
            }
            fprintf(stderr, "GzipDataReader:decompressBatch failed getData at %lld\n", fileOffset);
            soft_exit(1);
	    }
    }

    char* uncompressed;
    _int64 uncompressedBytes;
    inner->getExtra(&uncompressed, &uncompressedBytes);
    uncompressedBytes = extraBytes - priorBytes; // limit to just mine, subtracting prior data already copied in
    bool all = decompress(compressed, compressedBytes, NULL, uncompressed + priorBytes, uncompressedBytes, &validBytes,
        continueBlock ? ContinueMultiBlock : StartMultiBlock);
    if (! all) {
        // todo: handle this situation!!
        fprintf(stderr, "GzipDataReader:decompressBatch too big at %lld\n", fileOffset);
        soft_exit(1);
    }
    validBytes += priorBytes; // add back offset
    //printf("file offset %lld decompress %lld -> %lld bytes, carry over %lld\n", fileOffset, compressedBytes, validBytes - priorBytes, priorBytes);
    startBytes = inner->isEOF() ? validBytes : validBytes - overflowBytes ;
    inner->advance(compressedBytes);
    offset = 0;
}

class GzipDataSupplier : public DataSupplier
{
public:
    GzipDataSupplier(const DataSupplier* i_inner, bool autoRelease)
        : DataSupplier(autoRelease), inner(i_inner)
    {}

    virtual DataReader* getDataReader(_int64 overflowBytes, double extraFactor = 0.0) const
    {
        // adjust extra factor for compression ratio
        double totalFactor = MAX_FACTOR * (1.0 + extraFactor);
        // get inner reader with no overflow since zlib can't deal with it
        DataReader* data = inner->getDataReader(0, totalFactor);
        // compute how many extra bytes are owned by this layer
        char* p;
        _int64 mine;
        data->getExtra(&p, &mine);
        mine = (_int64)(mine * MAX_FACTOR / totalFactor);
        // create new reader, telling it how many bytes it owns
        // it will subtract overflow off the end of each batch
        return new GzipDataReader(overflowBytes, mine, data, autoRelease);
    }

private:
    const DataSupplier* inner;
};

    DataSupplier*
DataSupplier::Gzip(
    const DataSupplier* inner,
    bool autoRelease)
{
    return new GzipDataSupplier(inner, autoRelease);
}

//
// MemMap
//

class MemMapDataReader : public DataReader
{
public:

    MemMapDataReader(int i_batchCount, _int64 i_batchSize, _int64 i_overflowBytes, _int64 i_batchExtra, bool autoRelease);

    virtual ~MemMapDataReader();
    
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
    
    void acquireLock()
    {
        if (batchCount != 1) {
            AcquireExclusiveLock(&lock);
        }
    }

    void releaseLock()
    {
        if (batchCount != 1) {
            ReleaseExclusiveLock(&lock);
        }
    }

    const int       batchCount; // number of batches
    const _int64    batchSizeParam; // bytes per batch, 0 for entire file
    _int64          batchSize; // actual batch size for this file
    const _int64    overflowBytes;
    const _int64    batchExtra; // extra bytes per batch
    const char*     fileName; // current file name for diagnostics
    _int64          fileSize; // total size of current file
    char*           currentMap; // currently mapped block if non-NULL
    _int64          currentMapOffset; // current file offset of mapped region
    _int64          currentMapStartSize; // start size of mapped region (not incl overflow)
    _int64          currentMapSize; // total valid size of mapped region (incl overflow)
    char*           extra; // extra data buffer
    int             extraUsed; // number of extra data buffers in use
    DataBatch*      extraBatches; // non-zero for each extra buffer that is in use
    int             currentExtraIndex; // index of extra block for current batch
    _int64          offset; // into current batch
    _uint32         currentBatch; // current batch number starting at 1
    _int64          startBytes; // in current batch
    _int64          validBytes; // in current batch
    FileMapper      mapper;
    SingleWaiterObject waiter; // flow control
    ExclusiveLock   lock; // lock around flow control members (currentBatch, extraUsed, etc.)
};
 

MemMapDataReader::MemMapDataReader(int i_batchCount, _int64 i_batchSize, _int64 i_overflowBytes, _int64 i_batchExtra, bool autoRelease)
    : DataReader(autoRelease),
    batchCount(i_batchCount),
        batchSizeParam(i_batchSize),
        overflowBytes(i_overflowBytes),
        batchExtra(i_batchExtra),
        currentBatch(1),
        extraUsed(0),
        currentMap(NULL),
        currentMapOffset(0),
        currentMapSize(0),
        currentExtraIndex(0),
        mapper()
{
    _ASSERT(batchCount > 0 && batchSizeParam >= 0 && batchExtra >= 0);
    if (batchExtra > 0) {
        extra = (char*) BigAlloc(batchCount * batchExtra);
        extraBatches = new DataBatch[batchCount];
        memset(extraBatches, 0, sizeof(DataBatch));
    } else {
        extra = NULL;
        extraBatches = NULL;
    }
    if (! (CreateSingleWaiterObject(&waiter) && InitializeExclusiveLock(&lock))) {
        fprintf(stderr, "MemMapDataReader: CreateSingleWaiterObject failed\n");
        soft_exit(1);
    }
}

MemMapDataReader::~MemMapDataReader()
{
    if (extra != NULL) {
        BigDealloc(extra);
        extra = NULL;
    }
    if (extraBatches != NULL) {
        delete [] extraBatches;
    }
    DestroyExclusiveLock(&lock);
    DestroySingleWaiterObject(&waiter);
}

    bool
MemMapDataReader::init(
    const char* i_fileName)
{
    if (! mapper.init(i_fileName)) {
        return false;
    }
    fileName = i_fileName;
    fileSize = mapper.getFileSize();
    batchSize = batchSizeParam == 0 ? fileSize : batchSizeParam;
    return true;
}

    char*
MemMapDataReader::readHeader(
    _int64* io_headerSize)
{
    *io_headerSize = min(*io_headerSize, fileSize);
    reinit(0, *io_headerSize);
    return currentMap;
}

    void
MemMapDataReader::reinit(
    _int64 i_startingOffset,
    _int64 amountOfFileToProcess)
{
    _ASSERT(i_startingOffset >= 0 && amountOfFileToProcess >= 0);
    if (currentMap != NULL) {
        mapper.unmap();
    }
    _int64 oldAmount = amountOfFileToProcess;
    _int64 startSize = amountOfFileToProcess == 0 ? fileSize - i_startingOffset
	: max((_int64) 0, min(fileSize - i_startingOffset, amountOfFileToProcess));
    amountOfFileToProcess = max((_int64)0, min(startSize + overflowBytes, fileSize - i_startingOffset));
    currentMap = mapper.createMapping(i_startingOffset, amountOfFileToProcess);
    if (currentMap == NULL) {
        fprintf(stderr, "MemMapDataReader: fail to map %s at %lld,%lld\n", fileName, i_startingOffset, amountOfFileToProcess);
        soft_exit(1);
    }
    acquireLock();
    currentMapOffset = i_startingOffset;
    currentMapStartSize = startSize;
    currentMapSize = amountOfFileToProcess;
    offset = 0;
    startBytes = min(batchSize, currentMapStartSize - (currentBatch - 1) * batchSize);
    validBytes = min(batchSize + overflowBytes, currentMapSize - (currentBatch - 1) * batchSize);
    currentBatch = 1;
    extraUsed = 1;
    currentExtraIndex = 0;
    if (extraBatches != NULL) {
        memset(extraBatches, 0, sizeof(DataBatch) * batchCount);
        extraBatches[currentExtraIndex] = currentBatch;
    }
    releaseLock();
    if (batchCount != 1) {
        SignalSingleWaiterObject(&waiter);
    }
}

    bool
MemMapDataReader::getData(
    char** o_buffer,
    _int64* o_validBytes,
    _int64* o_startBytes)
{
    if (offset >= startBytes) {
        return false;
    }
    *o_buffer = currentMap + (currentBatch - 1) * batchSize + offset;
    *o_validBytes = validBytes - offset;
    if (o_startBytes) {
        *o_startBytes = max((_int64)0, startBytes - offset);
    }
    return *o_validBytes > 0;
}

    void
MemMapDataReader::advance(
    _int64 bytes)
{
    _ASSERT(bytes >= 0);
    offset = min(offset + max((_int64)0, bytes), validBytes);
}

    void
MemMapDataReader::nextBatch()
{
    if (isEOF()) {
        return;
    }
    while (true) {
        acquireLock();
        if (extraBatches == NULL || extraUsed < batchCount) {
            currentBatch++;
            if (extraBatches != NULL) {
                bool found = false;
                for (int i = 0; i < batchCount; i++) {
                    if (extraBatches[i].batchID == 0) {
                        extraBatches[i].batchID = currentBatch;
                        currentExtraIndex = i;
                        found = true;
                        break;
                    }
                }
                _ASSERT(found);
                extraUsed++;
                //printf("MemMap nextBatch %d:%d = index %d used %d of %d\n", 0, currentBatch, currentExtraIndex, extraUsed, batchCount); 
                if (extraUsed == batchCount) {
                    ResetSingleWaiterObject(&waiter);
                }
            }
            releaseLock();
	    offset = max(offset, startBytes) - startBytes;
            startBytes = min(batchSize, currentMapStartSize - (currentBatch - 1) * batchSize);
            validBytes = min(batchSize + overflowBytes, currentMapSize - (currentBatch - 1) * batchSize);
            _ASSERT(validBytes >= 0);
            return;
        }
        releaseLock();
        WaitForSingleWaiterObject(&waiter);
    }
}

    bool
MemMapDataReader::isEOF()
{
    return currentBatch * batchSize  >= currentMapSize;
}
    
    DataBatch
MemMapDataReader::getBatch()
{
    return DataBatch(currentBatch);
}

    void
MemMapDataReader::releaseBatch(
    DataBatch batch)
{
    if (extraBatches == NULL) {
        return;
    }
    acquireLock();
    for (int i = 0; i < batchCount; i++) {
        if (extraBatches[i] == batch) {
            extraBatches[i].batchID = 0;
            _ASSERT(extraUsed > 0);
            extraUsed--;
	    //printf("MemMap: releaseBatch %d:%d = index %d now using %d of %d\n", batch.fileID, batch.batchID, i, extraUsed, batchCount);
            if (extraUsed == batchCount - 1) {
                SignalSingleWaiterObject(&waiter);
            }
            break;
        }
    }
    releaseLock();
}

    _int64
MemMapDataReader::getFileOffset()
{
    return currentMapOffset + (currentBatch - 1) * batchSize + offset;
}

    void
MemMapDataReader::getExtra(
    char** o_extra,
    _int64* o_length)
{
    if (extra == NULL) {
        *o_extra = NULL;
        *o_length = 0;
    } else {
        *o_extra = extra + currentExtraIndex * batchExtra;
        *o_length = batchExtra;
    }
}

class MemMapDataSupplier : public DataSupplier
{
public:
    MemMapDataSupplier(bool autoRelease) : DataSupplier(autoRelease) {}
    virtual DataReader* getDataReader(_int64 overflowBytes, double extraFactor = 0.0) const
    {
        _ASSERT(extraFactor >= 0 && overflowBytes >= 0);
        if (extraFactor == 0) {
            // no per-batch expansion factor, so can read entire file as a batch
            return new MemMapDataReader(1, 0, overflowBytes, 0, autoRelease);
        } else {
            // break up into 4Mb batches
            _int64 batch = 4 * 1024 * 1024;
            _int64 extra = (_int64)(batch * extraFactor);
            return new MemMapDataReader(autoRelease ? 2 : ThreadCount + min(ThreadCount * 3 / 4, 3), batch, overflowBytes, extra, autoRelease);
        }
    }
};

//
// BatchTracker
//

BatchTracker::BatchTracker(int i_capacity)
    : pending(i_capacity)
{
}

    void
BatchTracker::addRead(
    DataBatch batch)
{
    DataBatch::Key key = batch.asKey();
    unsigned* p = pending.tryFind(key);
    int n = 1;
    if (p != NULL) {
        n = ++(*p);
    } else {
        pending.put(key, 1);
    }
    //_ASSERT(pending.tryFind(DataBatch(batch.batchID, 1^batch.fileID).asKey) != p);
    //unsigned* q = pending.tryFind(key); _ASSERT(q && (p == NULL || p == q) && *q == n);
    //printf("thread %d tracker %lx addRead %u:%u = %d\n", GetCurrentThreadId(), this, batch.fileID, batch.batchID, n);
}

    bool
BatchTracker::removeRead(
    DataBatch removed)
{
    DataBatch::Key key = removed.asKey();
    unsigned* p = pending.tryFind(key);
    //printf("thread %d tracker %lx removeRead %u:%u = %d\n", GetCurrentThreadId(), this, removed.fileID, removed.batchID, p ? *p - 1 : -1);
    _ASSERT(p != NULL && *p > 0);
    if (p != NULL) {
        if (*p > 1) {
            pending.put(key, *p - 1);
            //unsigned* q = pending.tryFind(key); _ASSERT(q == p);
            //_ASSERT(pending.tryFind(DataBatch(removed.batchID, 1^removed.fileID).asKey) != p);
            return false;
        }
        pending.erase(key);
    }
    return true;
}

//
// public static suppliers
//

const DataSupplier* DataSupplier::MemMap[] =
{ new MemMapDataSupplier(false), new MemMapDataSupplier(true) };

#ifdef _MSC_VER
const DataSupplier* DataSupplier::Default[2] =
{ DataSupplier::WindowsOverlapped[false], DataSupplier::WindowsOverlapped[true] };
#else
const DataSupplier* DataSupplier::Default[2] = 
{ DataSupplier::MemMap[false], DataSupplier::MemMap[true] };
#endif

const DataSupplier* DataSupplier::GzipDefault[2] =
{ DataSupplier::Gzip(DataSupplier::Default[false], false), DataSupplier::Gzip(DataSupplier::Default[true], true) };

int DataSupplier::ThreadCount = 1;

volatile _int64 DataReader::ReadWaitTime = 0;
volatile _int64 DataReader::ReleaseWaitTime = 0;
