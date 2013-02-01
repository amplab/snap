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

using std::max;
using std::min;

//
// WindowsOverlapped
//

class WindowsOverlappedDataReader : public DataReader
{
public:

    WindowsOverlappedDataReader(double extraFactor, _int64 extraBytes);

    virtual ~WindowsOverlappedDataReader();
    
    virtual bool init(const char* fileName);

    virtual char* readHeader(_int64* io_headerSize);

    virtual void reinit(_int64 startingOffset, _int64 amountOfFileToProcess);

    virtual bool getData(char** o_buffer, _int64* o_validBytes);

    virtual void advance(_int64 bytes);

    virtual void nextBatch(bool dontRelease = false);

    virtual bool isEOF();

    virtual DataBatch getBatch();

    virtual void release(DataBatch batch);

    virtual _int64 getFileOffset();

    virtual void getExtra(char** o_extra, _int64* o_length);

private:
    
    void startIo();

    void waitForBuffer(unsigned bufferNumber);

    const char*         fileName;
    LARGE_INTEGER       fileSize;
    HANDLE              hFile;
  
    static const unsigned nBuffers = 3;
    static const unsigned bufferSize = 32 * 1024 * 1024 - 4096;
    static const unsigned maxReadSizeInBytes = 1024;

    enum BufferState {Empty, Reading, Full};

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
    };

    _int64              extraBytes;
    BufferInfo          bufferInfo[nBuffers];
    LARGE_INTEGER       readOffset;
    _int64              endingOffset;
    _uint32             nextBatchID;
    unsigned            nextBufferForReader;
    unsigned            nextBufferForConsumer;
};

WindowsOverlappedDataReader::WindowsOverlappedDataReader(double extraFactor, _int64 extraBytes)
{
    //
    // Initialize the buffer info struct.
    //
    
    // allocate all the data in one big block
    _ASSERT(extraFactor >= 0 && extraBytes >= 0);
    extraBytes = max(0LL, (_int64) (bufferSize * extraFactor + extraBytes));
    char* allocated = (char*) BigAlloc(nBuffers * (bufferSize + 1 + extraBytes));
    if (NULL == allocated) {
        fprintf(stderr,"WindowsOverlappedDataReader: unable to allocate IO buffer\n");
        exit(1);
    }
    for (unsigned i = 0 ; i < nBuffers; i++) {
        bufferInfo[i].buffer = allocated;
        allocated += bufferSize + 1; // +1 gives us a place to put a terminating null
        bufferInfo[i].extra = extraBytes > 0 ? allocated : NULL;
        allocated += extraBytes;

        bufferInfo[i].buffer[bufferSize] = 0;       // The terminating null.
        
        bufferInfo[i].lap.hEvent = CreateEvent(NULL,TRUE,FALSE,NULL);
        if (NULL == bufferInfo[i].lap.hEvent) {
            fprintf(stderr,"WindowsOverlappedDataReader: Unable to create event\n");
            exit(1);
        }

        bufferInfo[i].state = Empty;
        bufferInfo[i].isEOF = false;
        bufferInfo[i].offset = 0;
    }
    nextBatchID = 1;
    hFile = INVALID_HANDLE_VALUE;
    nextBufferForConsumer = nextBufferForReader = 0;
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

    if (!ReadFile(hFile,info->buffer,*io_headerSize,&info->validBytes,&info->lap)) {
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
    }

    nextBufferForReader = 0;
    nextBufferForConsumer = 0;    

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
    // Kick off IO, wait for the first buffer to be read and then skip until hitting the first newline.
    //
    startIo();
    waitForBuffer(nextBufferForConsumer);
}

    bool
WindowsOverlappedDataReader::getData(
    char** o_buffer,
    _int64* o_validBytes)
{
    BufferInfo *info = &bufferInfo[nextBufferForConsumer];
    if (info->isEOF && info->offset >= info->validBytes) {
        //
        // EOF.
        //
        return false;
    }

    if (info->offset > info->nBytesThatMayBeginARead /* todo: && !ignoreEndOfRange*/) {
        //
        // Past the end of our section.
        //
        return false;
    }

    if (info->state != Full) {
        waitForBuffer(nextBufferForConsumer);
    }
    
    *o_buffer = info->buffer + info->offset;
    *o_validBytes = info->validBytes - info->offset;
    return true;
}

    void
WindowsOverlappedDataReader::advance(
    _int64 bytes)
{
    BufferInfo* info = &bufferInfo[nextBufferForConsumer];
    _ASSERT(info->validBytes >= info->offset && bytes >= 0 && bytes <= info->validBytes - info->offset);
    info->offset += min((_int64) info->validBytes - info->offset, max((_int64) 0, bytes));
}

    void
WindowsOverlappedDataReader::nextBatch(bool dontRelease)
{
    _uint32 released = bufferInfo[nextBufferForConsumer].batchID;
    const unsigned advance = (nextBufferForConsumer + 1) % nBuffers;
    printf("nextBatch %d state %d\n", advance, bufferInfo[advance].state);
    if (bufferInfo[advance].state != Full) {
        waitForBuffer(advance);
    }

    _ASSERT(bufferInfo[nextBufferForConsumer].fileOffset + bufferInfo[nextBufferForConsumer].validBytes == 
                bufferInfo[advance].fileOffset);
        
    nextBufferForConsumer = advance;
    startIo();
    if (! dontRelease) {
        release(DataBatch(released));
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
WindowsOverlappedDataReader::release(
    DataBatch batch)
{
    for (int i = 0; i < nBuffers; i++) {
        BufferInfo* info = &bufferInfo[i];
        if (info->batchID <= batch.batchID) {
            switch (info->state) {
            case Empty:
                // nothing
                break;
            case Reading:
                // todo: cancel read operation?
                _ASSERT(false);
            case Full:
                info->state = Empty;
            }
        }
    }
    startIo();
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
    *o_extra = bufferInfo[nextBufferForConsumer].extra;
    *o_length = extraBytes;
}
    
    void
WindowsOverlappedDataReader::startIo()
{
    //
    // Launch reads on whatever buffers are ready.
    //
    while (bufferInfo[nextBufferForReader].state == Empty) {
        BufferInfo *info = &bufferInfo[nextBufferForReader];
        info->batchID = nextBatchID++;

        if (readOffset.QuadPart >= fileSize.QuadPart || readOffset.QuadPart >= endingOffset + maxReadSizeInBytes) {
            info->validBytes = 0;
            info->nBytesThatMayBeginARead = 0;
            info->isEOF = readOffset.QuadPart >= fileSize.QuadPart;
            info->state = Full;
            SetEvent(info->lap.hEvent);
            return;
        }

        unsigned amountToRead;
        if (fileSize.QuadPart - readOffset.QuadPart > bufferSize && endingOffset + maxReadSizeInBytes - readOffset.QuadPart > bufferSize) {
            amountToRead = bufferSize;

            if (readOffset.QuadPart + amountToRead > endingOffset) {
                info->nBytesThatMayBeginARead = (unsigned)(endingOffset - readOffset.QuadPart);
            } else {
                info->nBytesThatMayBeginARead = amountToRead;
            }
            info->isEOF = false;
        } else {
            amountToRead = (unsigned)__min(fileSize.QuadPart - readOffset.QuadPart,endingOffset+maxReadSizeInBytes - readOffset.QuadPart);
            if (endingOffset <= readOffset.QuadPart) {
                //
                // We're only reading this for overflow buffer.
                //
                info->nBytesThatMayBeginARead = 0;
            } else {
                info->nBytesThatMayBeginARead = __min(amountToRead,(unsigned)(endingOffset - readOffset.QuadPart));    // Don't begin a read past endingOffset
            }
            info->isEOF = readOffset.QuadPart + amountToRead >= fileSize.QuadPart;
        }

        _ASSERT(amountToRead >= info->nBytesThatMayBeginARead || !info->isEOF || fileSize.QuadPart == readOffset.QuadPart + amountToRead);
        ResetEvent(info->lap.hEvent);
        info->lap.Offset = readOffset.LowPart;
        info->lap.OffsetHigh = readOffset.HighPart;
        info->fileOffset = readOffset.QuadPart;
         
        printf("startIo on %d at %lld for %uB\n", nextBufferForReader, readOffset, amountToRead);
        if (!ReadFile(
                hFile,
                info->buffer,
                amountToRead,
                &info->validBytes,
                &info->lap)) {

            if (GetLastError() != ERROR_IO_PENDING) {
                fprintf(stderr,"FASTQReader::startIo(): readFile failed, %d\n",GetLastError());
                exit(1);
            }
        }

        readOffset.QuadPart += amountToRead;
        info->state = Reading;
        info->offset = 0;

        nextBufferForReader = (nextBufferForReader + 1) % nBuffers;
    }
}

    void
WindowsOverlappedDataReader::waitForBuffer(
    unsigned bufferNumber)
{
    BufferInfo *info = &bufferInfo[bufferNumber];

    if (info->state == Full) {
        return;
    }

    if (info->state != Reading) {
        startIo();
    }

    if (!GetOverlappedResult(hFile,&info->lap,&info->validBytes,TRUE)) {
        fprintf(stderr,"Error reading FASTQ file, %d\n",GetLastError());
        exit(1);
    }

    info->state = Full;
    info->buffer[info->validBytes] = 0;
    ResetEvent(info->lap.hEvent);
}

    
class WindowsOverlappedDataSupplier : public DataSupplier
{
public:
    virtual DataReader* getDataReader(double extraFactor, _int64 extraBytes) const
    {
        return new WindowsOverlappedDataReader(extraFactor, extraBytes);
    }
};

const DataSupplier* DataSupplier::WindowsOverlapped = new WindowsOverlappedDataSupplier();

//
// Gzip
//

static const int windowBits = 15;
static const int ENABLE_ZLIB_GZIP = 32;

static const double MIN_FACTOR = 1.2;
static const double MAX_FACTOR = 4.0;

class GzipDataReader : public DataReader
{
public:

    GzipDataReader(DataReader* i_inner, _int64 i_extraBytes);

    virtual ~GzipDataReader();
    
    virtual bool init(const char* fileName);

    virtual char* readHeader(_int64* io_headerSize);

    virtual void reinit(_int64 startingOffset, _int64 amountOfFileToProcess);

    virtual bool getData(char** o_buffer, _int64* o_validBytes);

    virtual void advance(_int64 bytes);

    virtual void nextBatch(bool dontRelease = false);

    virtual bool isEOF();

    virtual DataBatch getBatch();

    virtual void release(DataBatch batch);

    virtual _int64 getFileOffset();

    virtual void getExtra(char** o_extra, _int64* o_length);

private:

    enum DecompressMode { SingleBlock, ContinueMultiBlock, StartMultiBlock };
    bool decompress(char* input, _int64 inputBytes, _int64* o_inputRead,
        char* output, _int64 outputBytes, _int64* o_outputWritten,
        DecompressMode mode);

    void decompressBatch();
    
    const _int64    extraBytes; // extra bytes I have to use
    DataReader*     inner; // inner reader to get data chunks
    z_stream        zstream; // stateful zlib stream
    _int64          validBytes; // valid bytes in extra data
    _int64          offset; // current offset in extra data
    bool            gotBatchData; // whether data has been read & decompressed for current batch
    bool            continueBlock; // whether to continue a block or start a new one
};
 
GzipDataReader::GzipDataReader(DataReader* i_inner, _int64 i_extraBytes)
    : extraBytes(i_extraBytes),
    inner(i_inner),
    zstream()
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
    zstream.zalloc = Z_NULL;
    zstream.zfree = Z_NULL;
    zstream.opaque = Z_NULL;
    return true;
}

    char*
GzipDataReader::readHeader(
    _int64* io_headerSize)
{
    _int64 compressedBytes = *io_headerSize / MIN_FACTOR;
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
}

    bool
GzipDataReader::getData(
    char** o_buffer,
    _int64* o_validBytes)
{
    if (! gotBatchData) {
        decompressBatch();
        gotBatchData = true;
    }
    inner->getExtra(o_buffer, o_validBytes);
    *o_buffer += offset;
    *o_validBytes = validBytes - offset;
    return offset < validBytes;
}

    void
GzipDataReader::advance(
    _int64 bytes)
{
    offset = min(offset + max(0LL, bytes), validBytes);
}

    void
GzipDataReader::nextBatch(bool dontRelease)
{
    inner->nextBatch(dontRelease);
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
GzipDataReader::release(
    DataBatch batch)
{
    inner->release(batch);
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
    zstream.next_in = (Bytef*) input;
    zstream.avail_in = inputBytes;
    zstream.next_out = (Bytef*) output;
    zstream.avail_out = outputBytes;
    uInt oldAvail;
    bool first = true;
    do {
        if (mode != ContinueMultiBlock || ! first) {
            int status = inflateInit2(&zstream, windowBits | ENABLE_ZLIB_GZIP);
            if (status < 0) {
                fprintf(stderr, "GzipDataReader: inflateInit2 failed with %d\n", status);
                return false;
            }
        }
        oldAvail = zstream.avail_out;
        int status = inflate(&zstream, mode == SingleBlock ? Z_NO_FLUSH : Z_FINISH);
        if (status < 0 && status != Z_BUF_ERROR) {
            fprintf(stderr, "GzipDataReader: inflate failed with %d\n", status);
            exit(1);
        }
        first = false;
    } while (zstream.avail_in != 0 && zstream.avail_out != oldAvail && mode != SingleBlock);
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
    char* compressed;
    _int64 compressedBytes;
    if (! inner->getData(&compressed, &compressedBytes)) {
        fprintf(stderr, "GzipDataReader:decompressBatch failed getData at %lld\n", inner->getFileOffset());
        exit(1);
    }

    char* uncompressed;
    _int64 uncompressedBytes;
    inner->getExtra(&uncompressed, &uncompressedBytes);
    uncompressedBytes = extraBytes; // limit to just mine

    bool all = decompress(compressed, compressedBytes, NULL, uncompressed, uncompressedBytes, &validBytes,
        continueBlock ? ContinueMultiBlock : StartMultiBlock);
    if (! all) {
        // todo: handle this situation!!
        fprintf(stderr, "GzipDataReader:decompressBatch too big at %lld\n", inner->getFileOffset());
        exit(1);
    }
    inner->advance(compressedBytes);
    offset = 0;
}

class GzipDataSupplier : public DataSupplier
{
public:
    GzipDataSupplier(const DataSupplier* i_inner)
        : inner(i_inner)
    {}

    virtual DataReader* getDataReader(double extraFactor, _int64 extraBytes) const
    {
        double totalFactor = MAX_FACTOR * (1.0 + extraFactor);
        DataReader* data = inner->getDataReader(totalFactor, extraBytes);
        char* p;
        _int64 mine;
        data->getExtra(&p, &mine);
        mine = (mine - extraBytes) * MAX_FACTOR / totalFactor;
        return new GzipDataReader(data, mine);
    }

private:
    const DataSupplier* inner;
};

    DataSupplier*
DataSupplier::Gzip(
    const DataSupplier* inner)
{
    return new GzipDataSupplier(inner);
}

//
// MemMap
//

class MemMapDataReader : public DataReader
{
public:

    MemMapDataReader(int i_batchCount, _int64 i_batchSize, _int64 i_batchExtra);

    virtual ~MemMapDataReader();
    
    virtual bool init(const char* fileName);

    virtual char* readHeader(_int64* io_headerSize);

    virtual void reinit(_int64 startingOffset, _int64 amountOfFileToProcess);

    virtual bool getData(char** o_buffer, _int64* o_validBytes);

    virtual void advance(_int64 bytes);

    virtual void nextBatch(bool dontRelease = false);

    virtual bool isEOF();

    virtual DataBatch getBatch();

    virtual void release(DataBatch batch);

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
    const _int64    batchExtra; // extra bytes per batch
    const char*     fileName; // current file name for diagnostics
    _int64          fileSize; // total size of current file
    char*           currentMap; // currently mapped block if non-NULL
    _int64          currentMapOffset; // current file offset of mapped region
    _int64          currentMapSize; // current size of mapped region
    char*           extra; // extra data buffer
    _int64          offset; // into current batch
    _uint32         currentBatch; // current batch number starting at 1
    _int64          validBytes; // in current batch
    _uint32         releasedBatch; // largest batch number released
    FileMapper      mapper;
    SingleWaiterObject waiter; // flow control
    ExclusiveLock   lock; // lock around flow control members (currentBatch, releasedBatch)
};
 

MemMapDataReader::MemMapDataReader(int i_batchCount, _int64 i_batchSize, _int64 i_batchExtra)
    : batchCount(i_batchCount),
        batchSizeParam(i_batchSize),
        batchExtra(i_batchExtra),
        currentBatch(1),
        releasedBatch(0),
        currentMap(NULL),
        currentMapOffset(0),
        currentMapSize(0),
        mapper()
{
    _ASSERT(batchCount > 0 && batchSizeParam >= 0 && batchExtra >= 0);
    extra = batchExtra > 0 ? (char*) BigAlloc(batchCount * batchExtra) : NULL;
    if (batchCount != 1) {
        if (! (CreateSingleWaiterObject(&waiter) && InitializeExclusiveLock(&lock))) {
            fprintf(stderr, "MemMapDataReader: CreateSingleWaiterObject failed\n");
            exit(1);
        }
    }
}

MemMapDataReader::~MemMapDataReader()
{
    if (extra != NULL) {
        BigDealloc(extra);
        extra = NULL;
    }
    if (batchCount != 1) {
        DestroyExclusiveLock(&lock);
        DestroySingleWaiterObject(&waiter);
    }
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
    amountOfFileToProcess = max(0LL, min(amountOfFileToProcess == 0 ? fileSize : amountOfFileToProcess, fileSize - i_startingOffset));
    currentMap = mapper.createMapping(i_startingOffset, amountOfFileToProcess);
    if (currentMap == NULL) {
        fprintf(stderr, "MemMapDataReader: fail to map %s at %lld,%lld\n", fileName, i_startingOffset, amountOfFileToProcess);
        exit(1);
    }
    acquireLock();
    currentMapOffset = i_startingOffset;
    currentMapSize = amountOfFileToProcess;
    offset = 0;
    validBytes = min(currentMapSize, batchSize);
    currentBatch = 1;
    releasedBatch = 0;
    releaseLock();
    if (batchCount != 1) {
        SignalSingleWaiterObject(&waiter);
    }
}

    bool
MemMapDataReader::getData(
    char** o_buffer,
    _int64* o_validBytes)
{
    if (batchCount != 1 && currentBatch - releasedBatch > batchCount) {
        WaitForSingleWaiterObject(&waiter);
    }
    *o_buffer = currentMap + (currentBatch - 1) * batchSize + offset;
    *o_validBytes = validBytes - offset;
    return *o_validBytes > 0;
}

    void
MemMapDataReader::advance(
    _int64 bytes)
{
    _ASSERT(bytes >= 0);
    offset = min(offset + max(0LL, bytes), validBytes);
}

    void
MemMapDataReader::nextBatch(bool dontRelease)
{
    if (isEOF()) {
        return;
    }
    acquireLock();
    currentBatch++;
    if (dontRelease) {
        _ASSERT(batchCount != 1);
        if (batchCount != 1 && currentBatch - releasedBatch > batchCount) {
            ResetSingleWaiterObject(&waiter);
        }
    } else {
        release(DataBatch(currentBatch - 1));
    }
    validBytes = min(batchSize, currentMapSize - (currentBatch - 1) * batchSize);
    releaseLock();
    _ASSERT(validBytes >= 0);
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
MemMapDataReader::release(
    DataBatch batch)
{
    if (batch.batchID > releasedBatch) {
        acquireLock();
        if (batchCount != 1 && currentBatch - releasedBatch > batchCount && currentBatch - batch.batchID <= batchCount) {
            SignalSingleWaiterObject(&waiter);
        }
        releasedBatch = min(currentBatch - 1, batch.batchID);
        releaseLock();
    }
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
        int index = (currentBatch - 1) % batchCount;
        *o_extra = extra + index * batchExtra;
        *o_length = batchExtra;
    }
}

class MemMapDataSupplier : public DataSupplier
{
public:
    virtual DataReader* getDataReader(double extraFactor, _int64 extraBytes) const
    {
        _ASSERT(extraFactor >= 0 && extraBytes >= 0);
        if (extraFactor == 0) {
            // no per-batch expansion factor, so can read entire file as a batch
            return new MemMapDataReader(1, 0, extraBytes);
        } else {
            // break up into 16Mb batches
            _int64 batch = 16 * 1024 * 1024;
            _int64 extra = batch * extraFactor + extraBytes;
            return new MemMapDataReader(3, batch, extra);
        }
    }
};

const DataSupplier* DataSupplier::MemMap = new MemMapDataSupplier();
