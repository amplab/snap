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
#include "Bam.h"
#include "zlib.h"
#include "exit.h"

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

//
// Decompress
//

static const int windowBits = 15;
static const int ENABLE_ZLIB_GZIP = 32;

static const double MIN_FACTOR = 1.2;
static const double MAX_FACTOR = 10.0;

class DecompressDataReader : public DataReader
{
public:

    DecompressDataReader(bool autoRelease, DataReader* i_inner, int i_count, _int64 totalExtra, _int64 i_extraBytes, _int64 i_overflowBytes, int i_chunkSize = BAM_BLOCK);

    virtual ~DecompressDataReader();

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

    enum DecompressMode { SingleBlock, ContinueMultiBlock, StartMultiBlock };

    static bool decompress(z_stream* zstream, ThreadHeap* heap, char* input, _int64 inputSize, _int64* o_inputUsed,
        char* output, _int64 outputSize, _int64* o_outputUsed, DecompressMode mode);

    // debugging
    char* findPointer(void* p);

private:


    static void decompressThread(void *context);

    static void decompressThreadContinuous(void *context);

    friend class DecompressManager;
    friend class DecompressWorker;

    enum EntryState
    {
        EntryReady, // for reading by client, on first list
        EntryHeld, // finished reading but not released, not on a list
        EntryAvailable, // released by client, on available list
        EntryReading // reading or decompressing, not on a list
    };

    struct Entry
    {
        Entry* next; // next entry on first/available list
        EntryState state;
        DataBatch batch;
        char* compressed;
        _int64 compressedStart; // limit to start a new zip block
        _int64 compressedValid; // total available data
        char* decompressed;
        _int64 decompressedStart;
        _int64 decompressedValid;
        bool allocated; // if decompressed has been allocated specially, not from inner extra data
    };

    // use only these routines to manipulate the linked  lists
    Entry* peekReady(); // from first, block if none
    void popReady(); // from first
    void enqueueReady(Entry* entry); // as last
    Entry* dequeueAvailable(); // as available, block if none
    void enqueueAvailable(Entry* entry); // from available

    DataReader* inner; // inner reader for compressed data
    const _int64 extraBytes; // number of bytes of extra that I get to use
    const _int64 overflowBytes; // overflow between batches
    const _int64 totalExtra; // total extra data
    const int chunkSize; // max size of decompressed data
    _int64 offset; // into current entry
    bool threadStarted; // whether thread has been started
    bool eof; // true when we've read to eof of previous
    volatile bool stopping; // set to stop everything
    EventObject decompressThreadDone; // signalled by background thread on exit

    // entry lists
    Entry* entries; // ring buffer of batches from inner reader
    int count; // # of entries
    Entry* first; // first ready buffer, NULL if none, currently being read by client
    Entry* last; // last ready buffer, NULL if none
    EventObject readyEvent; // signalled by bg thread when first goes NULL->non-NULL
    Entry* available; // first non-ready buffer (head of freelist), NULL if none
    EventObject availableEvent; // signalled by main thread when available goes NULL->non-NULL
    ExclusiveLock lock; // lock on linked list pointers in this object and in Entry
};


DecompressDataReader::DecompressDataReader(
    bool autoRelease,
    DataReader* i_inner,
    int i_count,
    _int64 i_totalExtra,
    _int64 i_extraBytes,
    _int64 i_overflowBytes,
    int i_chunkSize)
    : DataReader(autoRelease), inner(i_inner), count(i_count), offset(i_overflowBytes),
    totalExtra(i_totalExtra), extraBytes(i_extraBytes), overflowBytes(i_overflowBytes),
    chunkSize(i_chunkSize), threadStarted(false), eof(false), stopping(false)
{
    entries = new Entry[count];
    for (int i = 0; i < count; i++) {
        Entry* entry = &entries[i];
        entry->state = EntryAvailable;
        entry->next = i < count - 1 ? &entries[i + 1] : NULL;
        entry->decompressed = NULL;
        entry->allocated = false;
    }
    available = entries;
    first = last = NULL;
    CreateEventObject(&readyEvent);
    PreventEventWaitersFromProceeding(&readyEvent);
    CreateEventObject(&availableEvent);
    AllowEventWaitersToProceed(&availableEvent);
    CreateEventObject(&decompressThreadDone);
    PreventEventWaitersFromProceeding(&decompressThreadDone);
    InitializeExclusiveLock(&lock);
}

DecompressDataReader::~DecompressDataReader()
{
    if (threadStarted) {
        stopping = true;
        AllowEventWaitersToProceed(&availableEvent);
        WaitForEvent(&decompressThreadDone);
    }
    for (int i = 0;  i < count; i++) {
        if (entries[i].allocated) {
            BigDealloc(entries[i].decompressed);
        }
    }
    DestroyExclusiveLock(&lock);
    delete inner;
}

    bool
DecompressDataReader::init(
    const char* fileName)
{
    return inner->init(fileName);
}

    char*
DecompressDataReader::readHeader(
    _int64* io_headerSize)
{
    z_stream zstream;
    ThreadHeap heap(max(chunkSize,1000));
    _int64 compressedBytes = (_int64)(*io_headerSize / MIN_FACTOR);
    char* compressed = inner->readHeader(&compressedBytes);
    char* header;
    _int64 total;
    inner->getExtra(&header, &total);
    _ASSERT(total == totalExtra);
    _int64 headerSize = 0;
    while (headerSize < *io_headerSize && compressedBytes > 0) {
        _int64 compressedBlockSize, decompressedBlockSize;
        decompress(&zstream, chunkSize != 0 ? &heap : NULL,
            compressed, compressedBytes, &compressedBlockSize,
            header + headerSize, totalExtra - headerSize, &decompressedBlockSize,
            StartMultiBlock);
        inner->advance(compressedBlockSize);
        compressed += compressedBlockSize;
        compressedBytes -= compressedBlockSize;
        headerSize += decompressedBlockSize;
    }
    *io_headerSize = headerSize;
    return header;
}

    void
DecompressDataReader::reinit(
    _int64 startingOffset,
    _int64 amountOfFileToProcess)
{
    if (threadStarted) {
        fprintf(stderr, "DecompressDataReader reinit called twice\n");
        soft_exit(1);
    }
    // todo: transform start/amount to add for compression? I don't think so...
    inner->reinit(startingOffset, amountOfFileToProcess);
    threadStarted = true;
    if (! StartNewThread(chunkSize > 0 ? decompressThread : decompressThreadContinuous, this)) {
        fprintf(stderr, "failed to start decompressThread\n");
        soft_exit(1);
    }
}

    bool
DecompressDataReader::getData(
    char** o_buffer,
    _int64* o_validBytes,
    _int64* o_startBytes)
{
    if (eof) {
        return false;
    }
    Entry* entry = peekReady();
    if (offset >= entry->decompressedStart) {
        return false;
    }
    *o_buffer = entry->decompressed + offset;
    *o_validBytes = entry->decompressedValid - offset;
    if (o_startBytes != NULL) {
        *o_startBytes = entry->decompressedStart- offset;
    }
    return true;
}

    void
DecompressDataReader::advance(
    _int64 bytes)

{
    offset = min(offset + max(bytes, (_int64) 0), peekReady()->decompressedValid);
}

    void
DecompressDataReader::nextBatch()
{
    if (eof) {
        return;
    }
    Entry* old = peekReady();
    popReady();
    if (old->decompressedValid == overflowBytes) {
        eof = true;
        return;
    }
    Entry* next = peekReady();
    _ASSERT(next->state == EntryReady);
    for (int i = 0; i < count; i++) {
        _ASSERT(next == &entries[i] || entries[i].state == EntryAvailable || entries[i].decompressed != next->decompressed);
    }
    _int64 copy = old->decompressedValid - max(offset, old->decompressedStart);
    memcpy(next->decompressed + overflowBytes - copy, old->decompressed + old->decompressedValid - copy, copy);
    offset = overflowBytes - copy;
    //printf("DecompressDataReader nextBatch %d:%d #%d -> %d:%d #%d copy %lld + %lld/%lld\n", old->batch.fileID, old->batch.batchID, old-entries, next->batch.fileID, next->batch.batchID, next-entries, copy, next->decompressedStart, next->decompressedValid);
    if (autoRelease) {
        releaseBatch(old->batch);
    }
    if (offset == next->decompressedValid) {
        eof = true;
        _ASSERT(inner->isEOF());
    }
}

    bool
DecompressDataReader::isEOF()
{
    return eof;
}

    DataBatch
DecompressDataReader::getBatch()
{
    return peekReady()->batch;
}

    void
DecompressDataReader::releaseBatch(DataBatch batch)
{
    // loop backward from previous since that is most likely to be released
    for (int i = 0; i < count; i++) {
        Entry* entry = &entries[i];
        if (entry->batch == batch) {
            //printf("DecompressDataReader releaseBatch %d:%d #%d\n", batch.fileID, batch.batchID, i);
            enqueueAvailable(entry);
            break;
        }
    }
    inner->releaseBatch(batch);
}

    _int64
DecompressDataReader::getFileOffset()
{
    return inner->getFileOffset();
}

    void
DecompressDataReader::getExtra(
    char** o_extra,
    _int64* o_length)
{
    *o_extra = peekReady()->decompressed + extraBytes;
    *o_length = totalExtra - extraBytes;
}
    
    bool
DecompressDataReader::decompress(
    z_stream* zstream,
    ThreadHeap* heap,
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
    zstream->next_in = (Bytef*) input;
    zstream->avail_in = (uInt)inputBytes;
    zstream->next_out = (Bytef*) output;
    zstream->avail_out = (uInt)outputBytes;
    if (heap != NULL) {
        zstream->zalloc = zalloc;
        zstream->zfree = zfree;
        zstream->opaque = heap;
    } else {
        zstream->zalloc = NULL;
        zstream->zfree = NULL;
    }
    uInt oldAvailOut, oldAvailIn;
    int block = 0;
    bool multiBlock = true;
    int status;
    do {
            if (mode != ContinueMultiBlock || block != 0) {
            if (heap != NULL) {
                heap->reset();
            }
            status = inflateInit2(zstream, windowBits | ENABLE_ZLIB_GZIP);
            if (status < 0) {
                fprintf(stderr, "GzipDataReader: inflateInit2 failed with %d\n", status);
                return false;
            }
        }
        oldAvailOut = zstream->avail_out;
        oldAvailIn = zstream->avail_out;
        status = inflate(zstream, mode == SingleBlock ? Z_NO_FLUSH : Z_FINISH);
        // printf("decompress block #%d %lld -> %lld = %d\n", block, zstream.next_in - lastIn, zstream.next_out - lastOut, status);
        block++;
        if (status < 0 && status != Z_BUF_ERROR) {
            fprintf(stderr, "GzipDataReader: inflate failed with %d\n", status);
            soft_exit(1);
        }
        if (status < 0 && zstream->avail_out == 0 && zstream->avail_in > 0) {
            fprintf(stderr, "insufficient decompression buffer space - increase expansion factor, currently -xf %.1f\n", DataSupplier::ExpansionFactor);
            soft_exit(1);
        }
    } while (zstream->avail_in != 0 && (zstream->avail_out != oldAvailOut || zstream->avail_in != oldAvailIn) && mode != SingleBlock);
    // printf("end decompress status=%d, avail_in=%lld, last block=%lld->%lld, avail_out=%lld\n", status, zstream.avail_in, zstream.next_in - lastIn, zstream.next_out - lastOut, zstream.avail_out);
    if (o_inputRead) {
        *o_inputRead = inputBytes - zstream->avail_in;
    }
    if (o_outputWritten) {
        *o_outputWritten = outputBytes - zstream->avail_out;
    }
    return zstream->avail_in == 0;
}
  
    char*
DecompressDataReader::findPointer(
    void* p)
{
    static char result[100];
    sprintf(result, "not found");
    for (int i = 0; i < count; i++) {
        Entry* e = &entries[i];
        if (e->compressed <= p && p < e->compressed + e->compressedValid) {
            sprintf(result, "compressed #%d @ %lld", i, (char*)p - e->compressed);
            break;
        }
        if (e->decompressed <= p && p < e->decompressed + extraBytes) {
            sprintf(result, "decompressed #%d %lld", i, (char*) p - e->decompressed);
            break;
        }
        if (e->decompressed + extraBytes <= p && p < e->decompressed + totalExtra) {
            sprintf(result, "extra #%d %lld", i, (char*) p - e->decompressed - extraBytes);
            break;
        }
    }
    return result;
}

typedef VariableSizeVector<_int64> OffsetVector;

class DecompressWorker : public ParallelWorker
{
public:
    DecompressWorker();

    virtual void step();

private:
    z_stream zstream;
    ThreadHeap heap;
};
    
class DecompressManager: public ParallelWorkerManager
{
public:
    DecompressManager(OffsetVector* i_inputs, OffsetVector* i_outputs)
        : inputs(i_inputs), outputs(i_outputs)
    {}

    virtual ParallelWorker* createWorker()
    { return new DecompressWorker(); }

    OffsetVector* inputs;
    OffsetVector* outputs;
    DecompressDataReader::Entry* entry;

    friend class DecompressWorker;
};

DecompressWorker::DecompressWorker()
    : heap(BAM_BLOCK)
{
    zstream.zalloc = zalloc;
    zstream.zfree = zfree;
    zstream.opaque = &heap;
}

    void
DecompressWorker::step()
{
    DecompressManager* manager = (DecompressManager*) getManager();
    for (int i = getThreadNum(); i < manager->inputs->size() - 1; i += getNumThreads()) {
        _int64 inputUsed, outputUsed;
        DecompressDataReader::decompress(&zstream,
            &heap,
            manager->entry->compressed + (*manager->inputs)[i],
            (*manager->inputs)[i + 1] - (*manager->inputs)[i],
            &inputUsed,
            manager->entry->decompressed + (*manager->outputs)[i],
            (*manager->outputs)[i + 1] - (*manager->outputs)[i],
            &outputUsed,
            DecompressDataReader::SingleBlock);
        _ASSERT(inputUsed == (*manager->inputs)[i + 1] - (*manager->inputs)[i] &&
            outputUsed == (*manager->outputs)[i + 1] - (*manager->outputs)[i]);
    }
}

    void
DecompressDataReader::decompressThread(
    void* context)
{
    DecompressDataReader* reader = (DecompressDataReader*) context;
    OffsetVector inputs, outputs;
    DecompressManager manager(&inputs, &outputs);
    ParallelCoworker coworker(min(4, DataSupplier::ThreadCount), false, &manager);
    coworker.start();
    // keep reading & decompressing entries until stopped
    bool stop = false;
    while (! stop) {
        Entry* entry = reader->dequeueAvailable();
        if (reader->stopping) {
            break;
        }
        // always starts with a fresh batch - advances after reading it all
        bool ok = reader->inner->getData(&entry->compressed, &entry->compressedValid, &entry->compressedStart);
        int index = (int) (entry - reader->entries);
        if (! ok) {
            //printf("decompressThread #%d %d:%d eof\n", index, reader->inner->getBatch().fileID, reader->inner->getBatch().batchID);
            if (! reader->inner->isEOF()) {
                fprintf(stderr, "error reading file at offset %lld\n", reader->getFileOffset());
                soft_exit(1);
            }
            // mark as eof - no data
            entry->decompressedValid = entry->decompressedStart = reader->overflowBytes;
            DataBatch b = reader->inner->getBatch();
            entry->batch = DataBatch(b.batchID + 1, b.fileID);
            // decompressed buffer is same as next-to-last batch, need to allocate own buffer
            entry->decompressed = (char*) BigAlloc(reader->totalExtra);
            entry->allocated = true;
            stop = true;
        } else {
            _int64 ignore;
            reader->inner->getExtra(&entry->decompressed, &ignore);
            _ASSERT(ignore >= reader->extraBytes && ignore >= reader->overflowBytes);
            // figure out offsets and advance inner data
            inputs.clear();
            outputs.clear();
            _int64 input = 0;
            _int64 output = reader->overflowBytes;
            do {
                inputs.push_back(input);
                outputs.push_back(output);
                BgzfHeader* zip = (BgzfHeader*) (entry->compressed + input);
                input += zip->BSIZE() + 1;
                output += zip->ISIZE();
                _ASSERT(output <= reader->extraBytes && input <= entry->compressedValid);
                _ASSERT(zip->BSIZE() < BAM_BLOCK && zip->ISIZE() <= BAM_BLOCK);
            } while (input < entry->compressedStart);
            // append final offsets
            inputs.push_back(input);
            outputs.push_back(output);
            //printf("decompressThread read #%d %lld->%lld\n", index, input, output);
            reader->inner->advance(input);
            entry->decompressedValid = output;
            entry->decompressedStart = output - reader->overflowBytes;
            entry->batch = reader->inner->getBatch();
            reader->inner->nextBatch(); // start reading next batch
            // decompress all chunks synchronously on multiple threads
            manager.entry = entry;
            coworker.step();
        }
        // make buffer available for clients & go on to next
        //printf("decompressThread #%d %d:%d ready\n", index, entry->batch.fileID, entry->batch.batchID);
        reader->enqueueReady(entry);
    }
    coworker.stop();
    AllowEventWaitersToProceed(&reader->decompressThreadDone);
}

    void
DecompressDataReader::decompressThreadContinuous(
    void* context)
{
    DecompressDataReader* reader = (DecompressDataReader*) context;
    z_stream zstream;
    bool first = true;
    bool stop = false;
    while (! stop) {
        Entry* entry = reader->dequeueAvailable();
        if (reader->stopping) {
            break;
        }
        // always starts with a fresh batch - advances after reading it all
        bool ok = reader->inner->getData(&entry->compressed, &entry->compressedValid, &entry->compressedStart);
        int index = (int) (entry - reader->entries);
        if (! ok) {
            //printf("decompressThread #%d %d:%d eof\n", index, reader->inner->getBatch().fileID, reader->inner->getBatch().batchID);
            if (! reader->inner->isEOF()) {
                fprintf(stderr, "error reading file at offset %lld\n", reader->getFileOffset());
                soft_exit(1);
            }
            // mark as eof - no data
            entry->decompressedValid = entry->decompressedStart = reader->overflowBytes;
            DataBatch b = reader->inner->getBatch();
            entry->batch = DataBatch(b.batchID + 1, b.fileID);
            entry->decompressed = (char*) BigAlloc(reader->totalExtra);
            entry->allocated = true;
            stop = true;
        } else {
            // figure out offsets and advance inner data
            _int64 ignore;
            reader->inner->getExtra(&entry->decompressed, &ignore);
            _ASSERT(ignore >= reader->extraBytes && ignore >= reader->overflowBytes);
            _int64 compressedRead, decompressedWritten;
            entry->batch = reader->inner->getBatch();
            reader->inner->advance(entry->compressedValid);
            reader->inner->nextBatch(); // start reading next batch
            decompress(&zstream, NULL,
                entry->compressed, entry->compressedValid, &compressedRead,
                entry->decompressed + reader->overflowBytes, reader->extraBytes - reader->overflowBytes, &decompressedWritten,
                first ? StartMultiBlock : ContinueMultiBlock);
            _ASSERT(compressedRead == entry->compressedValid && decompressedWritten <= reader->extraBytes - reader->overflowBytes);
            entry->decompressedValid = reader->overflowBytes + decompressedWritten;
            entry->decompressedStart = decompressedWritten;
            first = false;
        }
        // make buffer available for clients & go on to next
        //printf("decompressThread #%d %d:%d ready\n", index, entry->batch.fileID, entry->batch.batchID);
        reader->enqueueReady(entry);
    }
    AllowEventWaitersToProceed(&reader->decompressThreadDone);
}

    DecompressDataReader::Entry*
DecompressDataReader::peekReady()
{
    // not thread-safe relative to popReady!
    if (first == NULL) {
        WaitForEvent(&readyEvent);
    }
    _ASSERT(first->state == EntryReady);
    return first;
}

    void
DecompressDataReader::popReady()
{
    while (true) {
        AcquireExclusiveLock(&lock);
        if (first != NULL) {
            _ASSERT(first->state == EntryReady);
            //printf("popReady %d:%d #%d -> held\n", first->batch.fileID, first->batch.batchID, first - entries);
            first->state = EntryHeld;
            if (first->next == NULL) {
                _ASSERT(last == first);
                last = NULL;
                PreventEventWaitersFromProceeding(&readyEvent);
            }
            first = first->next;
            _ASSERT(first == NULL || first->state == EntryReady);
            ReleaseExclusiveLock(&lock);
            return;
        }
        ReleaseExclusiveLock(&lock);
        WaitForEvent(&readyEvent);
    }
}

    void
DecompressDataReader::enqueueReady(Entry* entry)
{
    AcquireExclusiveLock(&lock);
    _ASSERT(entry->state == EntryReading);
    entry->next = NULL;
    entry->state = EntryReady;
    if (last == NULL) {
        first = last = entry;
        AllowEventWaitersToProceed(&readyEvent);
    } else {
        last->next = entry;
        last = entry;
    }
    ReleaseExclusiveLock(&lock);
}

    DecompressDataReader::Entry*
DecompressDataReader::dequeueAvailable()
{
    while (true) {
        AcquireExclusiveLock(&lock);
        //printf("dequeueAvailable #%d\n", available == NULL ? -1 : available - entries);
        if (available!= NULL) {
            _ASSERT(available->state == EntryAvailable);
            available->state = EntryReading;
            Entry* result = available;
            available = available->next;
            if (available == NULL) {
                PreventEventWaitersFromProceeding(&availableEvent);
            }
            ReleaseExclusiveLock(&lock);
            return result;
        }
        ReleaseExclusiveLock(&lock);
        WaitForEvent(&availableEvent);
        if (stopping) {
            return NULL;
        }
    }
}

    void
DecompressDataReader::enqueueAvailable(Entry* entry)
{
    AcquireExclusiveLock(&lock);
    _ASSERT(entry->state == EntryHeld);
    entry->state = EntryAvailable;
    entry->next = available;
    available = entry;
    if (entry->next == NULL) {
        AllowEventWaitersToProceed(&availableEvent);
    }
    ReleaseExclusiveLock(&lock);
}

class DecompressDataReaderSupplier : public DataSupplier
{
public:
    DecompressDataReaderSupplier(DataSupplier* i_inner, bool autoRelease, int i_blockSize = BAM_BLOCK)
        : DataSupplier(autoRelease), inner(i_inner), blockSize(i_blockSize)
    {}

    virtual DataReader* getDataReader(_int64 overflowBytes = 0, double extraFactor = 0.0);

private:
    DataSupplier* inner;
    const int blockSize;
};

    DataReader*
DecompressDataReaderSupplier::getDataReader(
    _int64 overflowBytes,
    double extraFactor)
{
    // adjust extra factor for compression ratio
    double expand = MAX_FACTOR * DataSupplier::ExpansionFactor;
    double totalFactor = expand * (1.0 + extraFactor);
    // get inner reader with no overflow since zlib can't deal with it
    DataReader* data = inner->getDataReader(blockSize, totalFactor);
    // compute how many extra bytes are owned by this layer
    char* p;
    _int64 totalExtra;
    data->getExtra(&p, &totalExtra);
    _int64 mine = (_int64)(totalExtra * expand / totalFactor);
    // create new reader, telling it how many bytes it owns
    // it will subtract overflow off the end of each batch
    // batch count here is cheap, so make it high enough that it never hits the limit
    return new DecompressDataReader(autoRelease, data, max(8, DataSupplier::ThreadCount*4), totalExtra, mine, overflowBytes, blockSize);
}
    
    DataSupplier*
DataSupplier::GzipBam(
    DataSupplier* inner,
    bool autoRelease)
{
    return new DecompressDataReaderSupplier(inner, autoRelease, BAM_BLOCK);
}
    
    DataSupplier*
DataSupplier::Gzip(
    DataSupplier* inner,
    bool autoRelease)
{
    return new DecompressDataReaderSupplier(inner, autoRelease, 0);
}

//
// MemMap
//

class MemMapDataSupplier : public DataSupplier
{
public:
    MemMapDataSupplier(bool autoRelease);

    virtual ~MemMapDataSupplier();

    virtual DataReader* getDataReader(_int64 overflowBytes, double extraFactor = 0.0);

    FileMapper* getMapper(const char* fileName);
    void releaseMapper(const char* fileName);

private:
    ExclusiveLock lock;
    map<string,FileMapper*> mappers;
    map<string,int> refcounts;
};

class MemMapDataReader : public DataReader
{
public:

    MemMapDataReader(MemMapDataSupplier* i_supplier, int i_batchCount, _int64 i_batchSize, _int64 i_overflowBytes, _int64 i_batchExtra, bool autoRelease);

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
    void*           currentMappedBase; // mapped base for unmap
    char*           extra; // extra data buffer
    int             extraUsed; // number of extra data buffers in use
    DataBatch*      extraBatches; // non-zero for each extra buffer that is in use
    int             currentExtraIndex; // index of extra block for current batch
    _int64          offset; // into current batch
    _uint32         currentBatch; // current batch number starting at 1
    _int64          startBytes; // in current batch
    _int64          validBytes; // in current batch
    MemMapDataSupplier* supplier;
    FileMapper*     mapper;
    SingleWaiterObject waiter; // flow control
    ExclusiveLock   lock; // lock around flow control members (currentBatch, extraUsed, etc.)
};
 

MemMapDataReader::MemMapDataReader(MemMapDataSupplier* i_supplier, int i_batchCount, _int64 i_batchSize, _int64 i_overflowBytes, _int64 i_batchExtra, bool autoRelease)
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
        supplier(i_supplier),
        mapper(NULL)
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
    if (mapper != NULL) {
        if (currentMap != NULL) {
            mapper->unmap(currentMappedBase);
        }
        supplier->releaseMapper(fileName);
    }
}

    bool
MemMapDataReader::init(
    const char* i_fileName)
{
    if (mapper != NULL) {
        supplier->releaseMapper(fileName);
    }
    mapper = supplier->getMapper(i_fileName);
    if (mapper == NULL) {
        return false;
    }
    fileName = i_fileName;
    fileSize = mapper->getFileSize();
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
        mapper->unmap(currentMappedBase);
    }
    _int64 oldAmount = amountOfFileToProcess;
    _int64 startSize = amountOfFileToProcess == 0 ? fileSize - i_startingOffset
        : max((_int64) 0, min(fileSize - i_startingOffset, amountOfFileToProcess));
    amountOfFileToProcess = max((_int64)0, min(startSize + overflowBytes, fileSize - i_startingOffset));
    currentMap = mapper->createMapping(i_startingOffset, amountOfFileToProcess, &currentMappedBase);
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

MemMapDataSupplier::MemMapDataSupplier(bool autoRelease)
    : DataSupplier(autoRelease)
{
    InitializeExclusiveLock(&lock);
}

MemMapDataSupplier::~MemMapDataSupplier()
{
    DestroyExclusiveLock(&lock);
}

    DataReader* 
MemMapDataSupplier::getDataReader(
    _int64 overflowBytes,
    double extraFactor)
{
    _ASSERT(extraFactor >= 0 && overflowBytes >= 0);
    if (extraFactor == 0) {
        // no per-batch expansion factor, so can read entire file as a batch
        return new MemMapDataReader(this, 1, 0, overflowBytes, 0, autoRelease);
    } else {
        // break up into 4Mb batches
        _int64 batch = 4 * 1024 * 1024;
        _int64 extra = (_int64)(batch * extraFactor);
        return new MemMapDataReader(this, autoRelease ? 2 : ThreadCount + min(ThreadCount * 3 / 4, 3), batch, overflowBytes, extra, autoRelease);
    }
}

    FileMapper*
MemMapDataSupplier::getMapper(
    const char* fileName)
{
    AcquireExclusiveLock(&lock);
    FileMapper* result = mappers[fileName];
    if (result == NULL) {
        result = new FileMapper();
        result->init(fileName);
        mappers[fileName] = result;
    }
    ++refcounts[fileName];
    ReleaseExclusiveLock(&lock);
    return result;
}

    void
MemMapDataSupplier::releaseMapper(
    const char* fileName)
{
    AcquireExclusiveLock(&lock);
    int n = refcounts[fileName];
    if (n > 0 && 0 == --refcounts[fileName]) {
        delete mappers[fileName];
        mappers[fileName] = NULL;
    }
    ReleaseExclusiveLock(&lock);
}

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

DataSupplier* DataSupplier::MemMap[] =
{ new MemMapDataSupplier(false), new MemMapDataSupplier(true) };

#ifdef _MSC_VER
DataSupplier* DataSupplier::Default[2] =
{ DataSupplier::WindowsOverlapped[false], DataSupplier::WindowsOverlapped[true] };
#else
DataSupplier* DataSupplier::Default[2] = 
{ DataSupplier::MemMap[false], DataSupplier::MemMap[true] };
#endif

// inner supplier must NOT be autorelease since wrapper keeps multiple buffers alive
DataSupplier* DataSupplier::GzipDefault[2] =
{ DataSupplier::Gzip(DataSupplier::Default[false], false), DataSupplier::Gzip(DataSupplier::Default[false], true) };

DataSupplier* DataSupplier::GzipBamDefault[2] =
{ DataSupplier::GzipBam(DataSupplier::Default[false], false), DataSupplier::GzipBam(DataSupplier::Default[false], true) };

int DataSupplier::ThreadCount = 1;

double DataSupplier::ExpansionFactor = 1.0;

volatile _int64 DataReader::ReadWaitTime = 0;
volatile _int64 DataReader::ReleaseWaitTime = 0;
