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
#include "Error.h"

using std::max;
using std::min;
using std::map;
using std::string;

//
// Read-Based
//
//
// A data reader that uses a read-type call to get its data (as opposed to memory mapping).
// This class contains the generic implementation, it must be subclassed to implement
// startIo() and waitForBuffer(), which do the actual IO.
//
class ReadBasedDataReader : public DataReader
{
public:

    ReadBasedDataReader(unsigned i_nBuffers, _int64 i_overflowBytes, double extraFactor, size_t bufferSpace = 0);

    virtual ~ReadBasedDataReader();
    
    virtual bool init(const char* fileName) = 0;

    char* readHeader(_int64* io_headerSize);

    virtual void reinit(_int64 startingOffset, _int64 amountOfFileToProcess);

    virtual bool getData(char** o_buffer, _int64* o_validBytes, _int64* o_startBytes = NULL);

    virtual void advance(_int64 bytes);

    virtual void nextBatch();

    virtual bool isEOF();

    virtual DataBatch getBatch();

    virtual void holdBatch(DataBatch batch);

    virtual bool releaseBatch(DataBatch batch);

    virtual _int64 getFileOffset();

    virtual void getExtra(char** o_extra, _int64* o_length);

    virtual const char* getFilename() = 0;

protected:
    
    // must hold the lock to call
    virtual void startIo() = 0;

    // must hold the lock to call
    virtual void waitForBuffer(unsigned bufferNumber) = 0;

    // must hold the lock to call
    virtual void addBuffer();
  
    static const unsigned BUFFER_SIZE = 4 * 1024 * 1024 - 4096;

    enum BufferState {Empty, Reading, Full, InUse};

    struct BufferInfo
    {
        char            *buffer;
        BufferState     state;
        unsigned        validBytes;
        unsigned        nBytesThatMayBeginARead;
        bool            isEOF;
        unsigned        offset;     // How far has the consumer gotten in current buffer
        _int64          fileOffset;
        _uint32         batchID;
        int             holds;
        char*           extra;
        int             next, previous; // index of next/previous in free/ready list, -1 if end
		bool			headerBuffer;	// Set if this is a special buffer that holds the rewound header.  These get read once and deallocated.

		void operator=(BufferInfo &peer) {
			buffer = peer.buffer;
			state = peer.state;
			validBytes = peer.validBytes;
			nBytesThatMayBeginARead = peer.nBytesThatMayBeginARead;
			isEOF = peer.isEOF;
			offset = peer.offset;
			fileOffset = peer.fileOffset;
			batchID = peer.batchID;
			holds = peer.holds;
			extra = peer.extra;
			next = peer.next;
			previous = peer.previous;
			headerBuffer = peer.headerBuffer;
		}
    };

    unsigned            nBuffers;
    const unsigned      maxBuffers;
	int					headerBuffersOutstanding;
	bool				startedReadingHeader;
    _int64              extraBytes;
    _int64              overflowBytes;
    BufferInfo*         bufferInfo;
    _uint32             nextBatchID;
    int                 nextBufferForReader; // list head (singly linked), -1 if empty
    int                 nextBufferForConsumer; // list head (doubly linked), -1 if empty
    int                 lastBufferForConsumer; // list tail, -1 if empty

    EventObject         releaseEvent;
    _int64              releaseWaitInMillis;

	ExclusiveLock       lock;

private:

	virtual bool getDataInternal(char** o_buffer, _int64* o_validBytes, _int64* o_startBytes = NULL);

	//
	// Stuff for handling the header read.  We allow arbitrarily large header reads, and service them by copying data from the underlying
	// data reader into a local buffer.  We might wind up reading more than the actual header, so we serve reads out of the header buffer
	// until it's used up.
	//
	char				*headerBuffer;
	_int64				 headerBufferSize;
	char				*headerExtra;		// Allocated in one go with the headerBuffer
	_int64				 headerExtraSize;
	_int64				 amountAdvancedThroughUnderlyingStoreByUs;
	unsigned			 nHeaderBuffersAllocated;
	bool			     hitEOFReadingHeader;
protected:
    const size_t         bufferSize;
};

ReadBasedDataReader::ReadBasedDataReader(
    unsigned i_nBuffers,
    _int64 i_overflowBytes,
    double extraFactor,
    size_t i_bufferSpace)
    : DataReader(), nBuffers(i_nBuffers), overflowBytes(i_overflowBytes),
    maxBuffers(i_nBuffers * (i_nBuffers == 1 ? 2 : 4)),
    bufferSize(i_bufferSpace > 0 ? i_bufferSpace / (i_nBuffers * 2) : BUFFER_SIZE),
	headerBuffer(NULL), headerBufferSize(0), amountAdvancedThroughUnderlyingStoreByUs(0), 
	headerExtra(NULL), headerExtraSize(0), startedReadingHeader(false), headerBuffersOutstanding(0), nHeaderBuffersAllocated(0),
	hitEOFReadingHeader(false)
{
    //
    // Initialize the buffer info struct.
    //
    
    // allocate all the data in one big block
    // NOTE: buffers are not null-terminated (since memmap version can't do it)
    _ASSERT(extraFactor >= 0 && i_nBuffers > 0);
    bufferInfo = new BufferInfo[maxBuffers];
    extraBytes = max((_int64) 0, (_int64) ((bufferSize + overflowBytes) * extraFactor));
    char* allocated = (char*) BigReserve(maxBuffers * (bufferSize + extraBytes + overflowBytes));
    BigCommit(allocated, nBuffers * (bufferSize + extraBytes + overflowBytes));
    if (NULL == allocated) {
        WriteErrorMessage("ReadBasedDataReader: unable to allocate IO buffer\n");
        soft_exit(1);
    }
    for (unsigned i = 0 ; i < nBuffers; i++) {
        bufferInfo[i].buffer = allocated;
        allocated += bufferSize + overflowBytes;
        bufferInfo[i].extra = extraBytes > 0 ? allocated : NULL;
        allocated += extraBytes;

        bufferInfo[i].state = Empty;
        bufferInfo[i].isEOF = false;
        bufferInfo[i].offset = 0;
        bufferInfo[i].next = i < nBuffers - 1 ? i + 1 : -1;
        bufferInfo[i].previous = i > 0 ? i - 1 : -1;
        bufferInfo[i].holds = 0;
		bufferInfo[i].headerBuffer = false;
    }
    nextBatchID = 1;
 
    nextBufferForConsumer = -1;
    lastBufferForConsumer = -1;
    nextBufferForReader = 0;
    CreateEventObject(&releaseEvent);
    releaseWaitInMillis = 5; // wait up to 5 ms before allocating a new buffer

    InitializeExclusiveLock(&lock);
}

ReadBasedDataReader::~ReadBasedDataReader()
{
    BigDealloc(bufferInfo[0].buffer);
    for (unsigned i = 0; i < nBuffers; i++) {
        bufferInfo[i].buffer = bufferInfo[i].extra = NULL;
    }

	if (NULL != headerBuffer) {
		delete[] headerBuffer;
		headerBuffer = NULL;
	}

	if (NULL != headerExtra) {
		delete[] headerExtra;
		headerExtra = NULL;
	}

    DestroyExclusiveLock(&lock);
    DestroyEventObject(&releaseEvent);
}

	char *
ReadBasedDataReader::readHeader(_int64* io_headerSize)
{
	_ASSERT(!startedReadingHeader);

	_int64 validBytesInHeader;
	if (NULL != headerBuffer) {
		if (*io_headerSize <= headerBufferSize) {
			return headerBuffer;
		}

		//
		// We need more data for the header.  Reallocate the buffer, copy the data from the old buffer into the new, and then get more
		// data from the underlying reader.
		//
		char *newHeaderBuffer = new char[*io_headerSize];
		memcpy(newHeaderBuffer, headerBuffer, headerBufferSize);

		delete[] headerBuffer;

		headerBuffer = newHeaderBuffer;

		validBytesInHeader = headerBufferSize;
		headerBufferSize = *io_headerSize;
	} else {
		reinit(0, 0);
		headerBufferSize = *io_headerSize;
		headerBuffer = new char[headerBufferSize];
		validBytesInHeader = 0;
	}

	//
	// Run through the underlying data provider getting data until we've filled the header buffer or hit EOF.
	//
	_int64 bytesLeftToGet = headerBufferSize - validBytesInHeader;
	_ASSERT(bytesLeftToGet);
	while (bytesLeftToGet != 0) {
		if (amountAdvancedThroughUnderlyingStoreByUs < validBytesInHeader) {
			_int64 amountToAdvance = validBytesInHeader - amountAdvancedThroughUnderlyingStoreByUs - overflowBytes;	// Leave overflowBytes left over, since we 
			if (amountToAdvance <= 0) {
				//
				// We're probably almost at EOF.  Consume the overflow bytes, too.
				//
				amountToAdvance = validBytesInHeader - amountAdvancedThroughUnderlyingStoreByUs;
			}
			advance(amountToAdvance);
			amountAdvancedThroughUnderlyingStoreByUs += amountToAdvance;
		}

		char *dataFromUnderlyingStore;
		_int64 dataSizeFromUnderlyingStore;

		if (!getDataInternal(&dataFromUnderlyingStore, &dataSizeFromUnderlyingStore)) {
			nextBatch();
			if (!getDataInternal(&dataFromUnderlyingStore, &dataSizeFromUnderlyingStore)) {
				//
				// Hit EOF while reading header.
				//
				hitEOFReadingHeader = true;
				headerBufferSize = *io_headerSize = validBytesInHeader;
				return headerBuffer;
			}
		}

		//
		// Adjust for the fact that we don't advance as far as we've read, so that we leave some overlap for
		// subsequent readers who want the data, not the header.
		//
		_ASSERT(amountAdvancedThroughUnderlyingStoreByUs <= validBytesInHeader);	// We haven't advanced over something we need.
		_int64 offsetIntoBuffer = validBytesInHeader - amountAdvancedThroughUnderlyingStoreByUs;
		_ASSERT(dataSizeFromUnderlyingStore >= offsetIntoBuffer);
		dataSizeFromUnderlyingStore -= offsetIntoBuffer;

		_int64 bytesToCopy = __min(dataSizeFromUnderlyingStore, bytesLeftToGet);

		memcpy(headerBuffer + validBytesInHeader, dataFromUnderlyingStore + offsetIntoBuffer, bytesToCopy);
		bytesLeftToGet -= bytesToCopy;
		validBytesInHeader += bytesToCopy;
	}

	return headerBuffer;	// No need to reset *io_headerSize, we read as much as was requested
}


//
// This gets called only for subclasses that can't implement their own.  It's able to put the header reads
// back on the queue, but can't seek anywhere else.
//
     void
ReadBasedDataReader::reinit(
    _int64 startingOffset,
    _int64 amountOfFileToProcess)
{
    AcquireExclusiveLock(&lock);

	if (0 != amountOfFileToProcess) {
		WriteErrorMessage("ReadBasedDataReader:reinit called with non-zero amountOfFileToProcess (%lld, %lld)\n", startingOffset, amountOfFileToProcess);
		soft_exit(1);
	}

	if (startedReadingHeader || 0 != headerBuffersOutstanding) {
		WriteErrorMessage("ReadBasedDataReader:reinit called after reading some data (%lld, %lld)\n", startingOffset, amountOfFileToProcess);
		soft_exit(1);
	}

	//
	// We've already read a bunch of data from the underlying reader during header read.  Create some new virtual buffers that point into the
	// header buffer, and stick them at the head of the "already read" list.
	//
	_ASSERT(!startedReadingHeader && 0 == headerBuffersOutstanding);
	if (0 != headerBufferSize) {
		startedReadingHeader = true;
	}

	//
	// First let any pending IO complete.
	//
	for (unsigned i = 0; i < nBuffers; i++) {
		if (bufferInfo[i].state == Reading) {
			waitForBuffer(i);
		}
	}

	_ASSERT(amountAdvancedThroughUnderlyingStoreByUs <= headerBufferSize);
	if (amountAdvancedThroughUnderlyingStoreByUs == 0) {
		nHeaderBuffersAllocated = 0;
	} else {
		nHeaderBuffersAllocated = (int)((amountAdvancedThroughUnderlyingStoreByUs + bufferSize - 1) / (bufferSize - overflowBytes));	// Round up, in case we read the last buffer.
	}
	int totalBuffersNeeded = (int)(maxBuffers + nHeaderBuffersAllocated);

	//
	// Reallocate the buffers array.
	//
	BufferInfo *newBuffers = new BufferInfo[totalBuffersNeeded];
	for (unsigned i = 0; i < maxBuffers; i++) {
		newBuffers[i] = bufferInfo[i];
	}

	delete[] bufferInfo;
	bufferInfo = newBuffers;
	//
	// Don't increase maxBuffers, so the buffer adder won't use the headerBuffers for anything.
	//

	//
	// Now construct the header buffers.
	//
	headerExtraSize = extraBytes * nHeaderBuffersAllocated;
	headerExtra = new char[headerExtraSize];

	char *headerPointer = headerBuffer;
	char *headerExtraPointer = headerExtra;
	_int64 fileOffset = 0;
	_int64 bytesRemaining = amountAdvancedThroughUnderlyingStoreByUs;

	for (int i = maxBuffers; i < totalBuffersNeeded; i++) {


		bufferInfo[i].state = Full;
		bufferInfo[i].isEOF = false;
		bufferInfo[i].offset = 0;
		bufferInfo[i].next = (i == totalBuffersNeeded - 1) ? nextBufferForConsumer : i + 1;
		bufferInfo[i].previous = (i == maxBuffers) ? -1 : i - 1;
		bufferInfo[i].holds = 0;
		bufferInfo[i].headerBuffer = true;
		bufferInfo[i].validBytes = (int)__min(bytesRemaining, (_int64) bufferSize);
		bufferInfo[i].nBytesThatMayBeginARead = (int)((bytesRemaining <= (_int64) bufferSize) ? bytesRemaining : __max(bufferInfo[i].validBytes - overflowBytes, 0));
		bufferInfo[i].offset = 0;
		bufferInfo[i].fileOffset = fileOffset;
		bufferInfo[i].batchID = nextBatchID++;

		bufferInfo[i].buffer = headerPointer;
		headerPointer += bufferInfo[i].nBytesThatMayBeginARead;	// NB: don't add overflowBytes; these buffers overlap
		fileOffset += bufferInfo[i].nBytesThatMayBeginARead;
		bufferInfo[i].extra = headerExtraPointer;
		headerExtraPointer += extraBytes;

		headerBuffersOutstanding++;
		bytesRemaining -= bufferInfo[i].nBytesThatMayBeginARead;
	}

	if (nHeaderBuffersAllocated > 0) {
		_ASSERT(bufferInfo[nextBufferForConsumer].previous == -1);
		bufferInfo[nextBufferForConsumer].previous = totalBuffersNeeded - 1;
		nextBufferForConsumer = maxBuffers;
		if (hitEOFReadingHeader) {
			bufferInfo[totalBuffersNeeded - 1].isEOF = true;
		}
	}

    //
    // Kick off IO, wait for the first buffer to be read
    //
    startIo();
    waitForBuffer(nextBufferForConsumer);

    ReleaseExclusiveLock(&lock);

	//
	// Now, consume data until we've gotten to startingOffset.
	//
	_int64 bytesToSkip = startingOffset;

	while (bytesToSkip > 0) {
		char *p;
		_int64 valid, start;
		bool ok = getData(&p, &valid, &start);
		if (!ok) {
			WriteErrorMessage("ReadBasedDataReader::init() failure getting data\n");
			soft_exit(1);
		}

		_int64 bytesToSkipThisTime = __min(valid, bytesToSkip);
		advance(bytesToSkipThisTime);
		if (bytesToSkipThisTime > start) {
			nextBatch();
		}
		getData(&p, &valid, &start);

		bytesToSkip -= bytesToSkipThisTime;
	}
}

	 bool
ReadBasedDataReader::getData(
	char** o_buffer,
	_int64* o_validBytes,
	_int64* o_startBytes)
{
	if (NULL != headerBuffer && !startedReadingHeader) {
		delete[] headerBuffer;
		headerBuffer = NULL;
		_ASSERT(NULL == headerExtra);
	}
	return getDataInternal(o_buffer, o_validBytes, o_startBytes);
}


    bool
ReadBasedDataReader::getDataInternal(
    char** o_buffer,
    _int64* o_validBytes,
    _int64* o_startBytes)
{
    _ASSERT(nextBufferForConsumer >= 0);
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
ReadBasedDataReader::advance(
    _int64 bytes)
{
    BufferInfo* info = &bufferInfo[nextBufferForConsumer];
    _ASSERT(info->validBytes >= info->offset && bytes >= 0 && bytes <= info->validBytes - info->offset);
    info->offset += min(info->validBytes - info->offset, (unsigned)max((_int64)0, bytes));
}

    void
ReadBasedDataReader::nextBatch()
{
    AcquireExclusiveLock(&lock);
    _ASSERT(nextBufferForConsumer >= 0);
    BufferInfo* info = &bufferInfo[nextBufferForConsumer];
    if (info->isEOF) {
        ReleaseExclusiveLock(&lock);
        if (info->holds == 0) {
            releaseBatch(DataBatch(info->batchID));
        }
        return;
    }
    DataBatch priorBatch = DataBatch(info->batchID);

    info->state = InUse;
    _uint32 overflow = max((unsigned) info->offset, info->nBytesThatMayBeginARead) - info->nBytesThatMayBeginARead;
    _int64 nextStart = info->fileOffset + info->nBytesThatMayBeginARead; // for validation
    //fprintf(stderr, "ReadBasedDataReader:nextBatch() finished buffer %d, starting buffer %d\n", nextBufferForConsumer, info->next);
    //fprintf(stderr, "ReadBasedDataReader:nextBatch() skipping %u overflow bytes used in previous batch\n", overflow);
    nextBufferForConsumer = info->next;

    bool first = true;
    while (nextBufferForConsumer == -1) {
        nextStart = 0; // can no longer count on getting sequential buffers from file
        ReleaseExclusiveLock(&lock);
        if (! first) {
            //fprintf(stderr, "ReadBasedDataReader::nextBatch thread %d wait for release\n", GetCurrentThreadId());
            _int64 start = timeInNanos();
            bool waitSucceeded = WaitForEventWithTimeout(&releaseEvent, releaseWaitInMillis);
             InterlockedAdd64AndReturnNewValue(&ReleaseWaitTime, timeInNanos() - start);
            //fprintf(stderr, "ReadBasedDataReader::nextBatch thread %d released\n", GetCurrentThreadId());
            if (!waitSucceeded) {
                AcquireExclusiveLock(&lock);
                addBuffer();
                ReleaseExclusiveLock(&lock);
            }
        }
        first = false;
        AcquireExclusiveLock(&lock);
        startIo();
    }

    if (bufferInfo[nextBufferForConsumer].state != Full) {
        waitForBuffer(nextBufferForConsumer);
    }

    bufferInfo[nextBufferForConsumer].offset = overflow;
    bufferInfo[nextBufferForConsumer].holds = 0;
    //fprintf(stderr,"emitting buffer starting at 0x%llx\n", info->fileOffset);
    //if (nextStart != 0) fprintf(stderr, "checking NextStart 0x%llx\n", nextStart);  
    _ASSERT(nextStart == 0 || nextStart == bufferInfo[nextBufferForConsumer].fileOffset || bufferInfo[nextBufferForConsumer].isEOF);

    ReleaseExclusiveLock(&lock);

    if (info->holds == 0) {
        releaseBatch(priorBatch);
    }
}

    bool
ReadBasedDataReader::isEOF()
{
    return bufferInfo[nextBufferForConsumer].isEOF;
}
    
    DataBatch
ReadBasedDataReader::getBatch()
{
    return DataBatch(bufferInfo[nextBufferForConsumer].batchID);
}

    void
ReadBasedDataReader::holdBatch(
	DataBatch batch)
{
	AcquireExclusiveLock(&lock);
	for (unsigned i = 0; i < maxBuffers + nHeaderBuffersAllocated; i = (i == nBuffers - 1) ? maxBuffers : i+1) {	// Goofy loop is because headerBuffers get tacked on beyond maxBuffers
		BufferInfo *info = &bufferInfo[i];
		if (info->batchID == batch.batchID) {
			//fprintf(stderr, "%x holdBatch batch 0x%x, holds on buffer %d now %d\n", (unsigned) this, batch.batchID, i, info->holds);
			info->holds++;
		}
	}
	ReleaseExclusiveLock(&lock);
}


    bool
ReadBasedDataReader::releaseBatch(
    DataBatch batch)
{
    AcquireExclusiveLock(&lock);

    bool released = false;
    bool result = true;
	for (unsigned i = 0; i < maxBuffers + nHeaderBuffersAllocated; i = (i == nBuffers - 1) ? maxBuffers : i + 1) {	// Goofy loop is because headerBuffers get tacked on beyond maxBuffers
        BufferInfo* info = &bufferInfo[i];
        if (info->batchID == batch.batchID) {
            switch (info->state) {
            case Empty:
                // should never happen
                break;

            case Reading:
                // should never happen
                _ASSERT(false);
                break;

            case InUse:
                released = info->holds <= 1;
                // fall through

            case Full:
                if (info->holds > 0) {
                    info->holds--;
                }
                if (info->holds == 0) {
                    //fprintf(stderr,"%x releaseBatch batch %d, releasing %s buffer %d\n", (unsigned) this, batch.batchID, info->state == InUse ? "InUse" : "Full", i);
                    info->state = Empty;
                    // remove from ready list
                    if (i == nextBufferForConsumer) {
                        //fprintf(stderr, "ReadBasedDataReader::releaseBatch change nextBufferForConsumer %d->%d\n", nextBufferForConsumer, info->next);
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

					if (info->headerBuffer) {
						// Header buffers never get reused.  Just get rid of it.
						info->buffer = NULL;
						info->extra = NULL;
						_ASSERT(headerBuffersOutstanding > 0);
                        if (headerBuffersOutstanding > 0) {
                            headerBuffersOutstanding--;
                            if (0 == headerBuffersOutstanding) {
                                delete[] headerBuffer;
                                delete[] headerExtra;
                                headerBuffer = headerExtra = NULL;
                                nHeaderBuffersAllocated = 0;
                            }
                        }
					} else {
						// add to head of free list
						info->next = nextBufferForReader;
						info->batchID = 0;
#ifdef _DEBUG
						memset(info->buffer, 0xde, bufferSize + extraBytes);
#endif
						nextBufferForReader = i;
					}
					result = true;
                } else {
                    //fprintf(stderr,"%x releaseBatch batch %d, holds on buffer %d now %d\n", (unsigned) this, batch.batchID, i, info->holds);
                    result = false;
                }
                break;

            default:
                WriteErrorMessage("ReadBasedDataReader::releaseBatch():invalid enum\n");
                soft_exit(1);
            }
        }
    }

    startIo();

    if (released) {
        //fprintf(stderr, "releaseBatch set releaseEvent\n");
        AllowEventWaitersToProceed(&releaseEvent);
    }

    ReleaseExclusiveLock(&lock);

    return result;
}

    _int64
ReadBasedDataReader::getFileOffset()
{
    return bufferInfo[nextBufferForConsumer].fileOffset + bufferInfo[nextBufferForConsumer].offset;
}

    void
ReadBasedDataReader::getExtra(
    char** o_extra,
    _int64* o_length)
{
	// hack: return valid buffer even when no consumer buffers - this may happen when reading header
	*o_extra = bufferInfo[max(0, nextBufferForConsumer)].extra;
	*o_length = extraBytes;
}
    

    void
ReadBasedDataReader::addBuffer()
{
    if (nBuffers == maxBuffers) {
        //fprintf(stderr, "ReadBasedDataReader: addBuffer at limit\n");
        return;
    }
    _ASSERT(nBuffers < maxBuffers);
    //fprintf(stderr, "ReadBasedDataReader: addBuffer %d of %d\n", nBuffers, maxBuffers);
    size_t bytes = bufferSize + extraBytes + overflowBytes;
    bufferInfo[nBuffers].buffer = bufferInfo[nBuffers-1].buffer + bytes;
    if (! BigCommit(bufferInfo[nBuffers].buffer, bytes)) {
        WriteErrorMessage("ReadBasedDataReader: unable to commit IO buffer\n");
        soft_exit(1);
    }
    bufferInfo[nBuffers].extra = extraBytes > 0 ? bufferInfo[nBuffers].buffer + bytes - extraBytes : NULL;


    bufferInfo[nBuffers].state = Empty;
    bufferInfo[nBuffers].isEOF= false;
    bufferInfo[nBuffers].offset = 0;
    bufferInfo[nBuffers].next = nextBufferForReader;
    bufferInfo[nBuffers].previous = -1;
    bufferInfo[nBuffers].headerBuffer = false;
    nextBufferForReader = nBuffers;
    nBuffers++;
    _ASSERT(nBuffers <= maxBuffers);
    if (nBuffers == maxBuffers) {
        releaseWaitInMillis = 1000 * 3600 * 24 * 7; // A week
    }
}

class StdioDataReader : public ReadBasedDataReader 
{
public:
    StdioDataReader(unsigned i_nBuffers, _int64 i_overflowBytes, double extraFactor);
    ~StdioDataReader();

    virtual bool init(const char* i_fileName);

    virtual const char* getFilename()
    { return "-"; }

 protected:
    
    // must hold the lock to call
    virtual void startIo();

    // must hold the lock to call
    virtual void waitForBuffer(unsigned bufferNumber);

private:
    //
    // Because reads don't necessarily divide evenly into buffers, we have to assure that
    // the buffers that we read can overlap.  In file-IO based readers, we do this by reading
    // a buffer's worth of data each time, but advancing the file pointer only by
    // bufferSize - overflowBytes, so each buffer ovelaps with its predecessor by a little.
    // That doesn't work for stdio, since it can't rewind.  So, instead, we allocate
    // storage on the side to hold a copy of the last overflowBytes
    // and then just copy those bytes into the beginning of the next buffer to read.
    // We also use this buffer to hold the header (the first read), and to allow
    // reading the header plus some extra data, parsing the header, and then seeking
    // backward to the actual end of the header.
    //

    char    *overflowBuffer;
    bool     overflowBufferFilled;   // For the very first read, there may be no overlap buffer data.

    bool    started;
    bool    hitEOF;

    _int64 readOffset;
};

StdioDataReader::StdioDataReader(unsigned i_nBuffers, _int64 i_overflowBytes, double extraFactor) :
    ReadBasedDataReader(i_nBuffers, i_overflowBytes, extraFactor), started(false), hitEOF(false), overflowBufferFilled(false),
    readOffset(0), overflowBuffer(NULL)
{
}

StdioDataReader::~StdioDataReader()
{
    BigDealloc(overflowBuffer);
    overflowBuffer = NULL;
}

bool
StdioDataReader::init(const char * i_fileName)
{
    if (strcmp(i_fileName, "-")) {
        WriteErrorMessage("StdioDataReader: must have filename of '-', got '%s'\n", i_fileName);
        soft_exit(1);
    }

#ifdef _MSC_VER
    int result = _setmode( _fileno( stdin ), _O_BINARY );  // puts stdin in to non-translated mode, so if we're reading compressed data windows' CRLF processing doesn't destroy it.
    if (-1 == result) {
        WriteErrorMessage("StdioDataReader::freopen to change to untranslated mode failed\n");
        soft_exit(1);
    }
#endif // _MSC_VER

    return true;
}

void
StdioDataReader::startIo()
{
	AssertExclusiveLockHeld(&lock);
	
	started = true;

    //
    // Synchronously read data into whatever buffers are ready.
    //
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
            //fprintf(stderr, "StdioDataReader::startIo set nextBufferForConsumder -1 -> %d\n", index);
			nextBufferForConsumer = index;
		}
       
        if (hitEOF) {
            info->validBytes = 0;
            info->buffer[0] = '\0';
            info->nBytesThatMayBeginARead = 0;
            info->isEOF = true;
            info->state = Full;
            return;
        }

        size_t amountToRead;
        size_t bufferOffset;
        if (overflowBufferFilled) {
            //
            // Copy the bytes from the overflow buffer into our buffer.
            //
			memcpy(info->buffer, overflowBuffer, overflowBytes);
			bufferOffset = overflowBytes;
			amountToRead = bufferSize - overflowBytes;
			info->fileOffset = readOffset - overflowBytes;
        } else {
            amountToRead = bufferSize;
            bufferOffset = 0;
            info->fileOffset = readOffset;
        }

        //
        // We have to run this holding the lock, because otherwise there's no way to make the overflow buffer work properly.  
        //
        size_t bytesRead = fread(info->buffer + bufferOffset, 1, amountToRead, stdin);
        //fprintf(stderr,"StdioDataReader:startIO(): Read offset 0x%llx into buffer at 0x%llx, size %d, copied 0x%x overflow bytes, start at 0x%llx, tid %d\n", readOffset, info->buffer, bytesRead, bufferOffset, readOffset - bufferOffset, GetCurrentThreadId());

        readOffset += bytesRead;

        if (bytesRead != amountToRead) {
            if (feof(stdin)) {
                info->isEOF = true;
                hitEOF = true;
            } else {
                WriteErrorMessage("StdinDataReader: Error reading stdin (but not EOF).\n");
                soft_exit(1);
            }
        } else {
            info->isEOF = false;
        }

        info->validBytes = (unsigned)(bytesRead + bufferOffset);

        if (hitEOF) {
            info->nBytesThatMayBeginARead = (unsigned)(bytesRead + bufferOffset);
            overflowBufferFilled = false;
        } else {
            info->nBytesThatMayBeginARead = (unsigned)(bytesRead + bufferOffset - overflowBytes);
            //
            // Fill the overflow buffer with the last bytes from this buffer.
            //
            if (NULL == overflowBuffer) {
                //
                // We can get here if we never called readHeader().  If so, we know we never will and so
                // we can just allocate the overflow buffer to be the size of the header.
                //
                overflowBuffer = (char *)BigAlloc(overflowBytes);
            }
            memcpy(overflowBuffer, info->buffer + bufferOffset + bytesRead - overflowBytes, overflowBytes);
            overflowBufferFilled = true;
        }
        info->state = Full;
    }

    if (nextBufferForConsumer == -1) {
        //fprintf(stderr, "startIo thread %x reset releaseEvent\n", GetCurrentThreadId());
        PreventEventWaitersFromProceeding(&releaseEvent);
    }
}
 
    void
StdioDataReader::waitForBuffer(
    unsigned bufferNumber)
{
    _ASSERT(bufferNumber >= 0 && (bufferNumber < nBuffers || bufferNumber >= maxBuffers && 0 != headerBuffersOutstanding));
    BufferInfo *info = &bufferInfo[bufferNumber];

    while (info->state == InUse) {
        //fprintf(stderr, "StdioDataReader::waitForBuffer %d InUse...\n", bufferNumber);
        // must already have lock to call, release & wait & reacquire
        ReleaseExclusiveLock(&lock);
	// TODO: implement timed wait on Linux
#ifdef _MSC_VER
        _int64 start = timeInNanos();
        _uint32 waitTime;
        if (releaseWaitInMillis > 0xffffffff) {
            waitTime = INFINITE;
        } else {
            waitTime = (_uint32)releaseWaitInMillis;
        }
        _uint32 result = WaitForSingleObject(releaseEvent, waitTime);
        InterlockedAdd64AndReturnNewValue(&ReleaseWaitTime, timeInNanos() - start);
#else
		WaitForEvent(&releaseEvent);
#endif
        AcquireExclusiveLock(&lock);
#ifdef _MSC_VER
        if (result == WAIT_TIMEOUT) {
            // this isn't going to directly make this buffer available, but will reduce pressure
            addBuffer();
        }
#endif
    }

    if (info->state == Full) {
        return;
    }

    _ASSERT(info->state != Reading);    // We're synchronous, we don't use Reading
    startIo();
    
    info->state = Full;
    info->buffer[info->validBytes] = 0;
}

class StdioDataSupplier : public DataSupplier
{
public:
    StdioDataSupplier() : DataSupplier() {}
    virtual DataReader* getDataReader(int bufferCount, _int64 overflowBytes, double extraFactor = 0.0, size_t bufferSpace = 0)
    {
        if (supplied) {
            WriteErrorMessage("You can only use stdin input for one run per execution of SNAP (i.e., if you use ',' to run SNAP more than once without reloading the index, you can only use stdin once)\n");
            soft_exit_no_print(1);
        }

        supplied = true;

        return new StdioDataReader(bufferCount, overflowBytes, extraFactor);
    }
private:

    static bool supplied;
};

bool StdioDataSupplier::supplied = false;


#ifdef _MSC_VER
class WindowsOverlappedDataReader : public ReadBasedDataReader
{
public:

    WindowsOverlappedDataReader(unsigned i_nBuffers, _int64 i_overflowBytes, double extraFactor, size_t bufferSpace);

    virtual ~WindowsOverlappedDataReader();
    
    virtual bool init(const char* i_fileName);

    virtual void reinit(_int64 startingOffset, _int64 amountOfFileToProcess);

//    virtual char* readHeader(_int64* io_headerSize);

    virtual const char* getFilename()
    { return fileName; }

 protected:
    
    // must hold the lock to call
    virtual void startIo();

    // must hold the lock to call
    virtual void waitForBuffer(unsigned bufferNumber);

    // must hold the lock to call
    virtual void addBuffer();

    OVERLAPPED          *bufferLaps;

    const char*         fileName;
    LARGE_INTEGER       fileSize;
    HANDLE              hFile;
  
    LARGE_INTEGER       readOffset;
    _int64              endingOffset;

};

WindowsOverlappedDataReader::WindowsOverlappedDataReader(unsigned i_nBuffers, _int64 i_overflowBytes, double extraFactor, size_t bufferSpace) :
    ReadBasedDataReader(i_nBuffers, i_overflowBytes, extraFactor, bufferSpace), fileName(NULL), hFile(INVALID_HANDLE_VALUE), endingOffset(0)
{
    readOffset.QuadPart = 0;
    bufferLaps = (OVERLAPPED *)malloc(sizeof(OVERLAPPED) * maxBuffers);
    for (unsigned i = 0; i < i_nBuffers; i++) {
        bufferLaps[i].hEvent = CreateEvent(NULL,TRUE,FALSE,NULL);
        if (NULL == bufferLaps[i].hEvent) {
            WriteErrorMessage("WindowsOverlappedDataReader: Unable to create event\n");
            soft_exit(1);
        }
    }
}

WindowsOverlappedDataReader::~WindowsOverlappedDataReader()
{
    for (unsigned i = 0; i < nBuffers; i++) {
       CloseHandle(bufferLaps[i].hEvent);
    }
    free(bufferLaps);
    bufferLaps = NULL;
    CloseHandle(hFile);
}

bool
WindowsOverlappedDataReader::init(const char* i_fileName)
{
    fileName = i_fileName;
    hFile = CreateFile(fileName,GENERIC_READ,FILE_SHARE_READ,NULL,OPEN_EXISTING,FILE_FLAG_OVERLAPPED,NULL);
    if (INVALID_HANDLE_VALUE == hFile) {
        return false;
    }

    if (!GetFileSizeEx(hFile,&fileSize)) {
        WriteErrorMessage("WindowsOverlappedDataReader: unable to get file size of '%s', %d\n",fileName,GetLastError());
        return false;
    }
    return true;
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


        void
WindowsOverlappedDataReader::startIo()
{
    //
    // Launch reads on whatever buffers are ready.
    //
    AssertExclusiveLockHeld(&lock);

    while (nextBufferForReader != -1) {
        // remove from free list
        BufferInfo* info = &bufferInfo[nextBufferForReader];
        OVERLAPPED *bufferLap = &bufferLaps[nextBufferForReader];
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
            SetEvent(bufferLap->hEvent);
            return;
        }

        unsigned amountToRead;
        _int64 finalOffset = min(fileSize.QuadPart, endingOffset + overflowBytes);
        _int64 finalStartOffset = min(fileSize.QuadPart, endingOffset);
        amountToRead = (unsigned)min(finalOffset - readOffset.QuadPart, (_int64) bufferSize);   // Cast OK because can't be longer than unsigned bufferSize
        info->isEOF = readOffset.QuadPart + amountToRead == finalOffset;
        info->nBytesThatMayBeginARead = (unsigned)min((_int64)bufferSize - overflowBytes, finalStartOffset - readOffset.QuadPart);

        _ASSERT(amountToRead >= info->nBytesThatMayBeginARead && (!info->isEOF || finalOffset == readOffset.QuadPart + amountToRead));
        ResetEvent(bufferLap->hEvent);
        bufferLap->Offset = readOffset.LowPart;
        bufferLap->OffsetHigh = readOffset.HighPart;
        info->fileOffset = readOffset.QuadPart;

        readOffset.QuadPart += info->nBytesThatMayBeginARead;
        info->state = Reading;
        info->offset = 0;
         
        //fprintf(stderr, "startIo on %d at %lld for %uB\n", index, readOffset, amountToRead);
        if (!ReadFile(
                hFile,
                info->buffer,
                amountToRead,
                (DWORD *)&info->validBytes,
                bufferLap)) {

            if (GetLastError() != ERROR_IO_PENDING) {
                WriteErrorMessage("WindowsOverlappedDataReader::startIo(): readFile failed, %d\n",GetLastError());
                soft_exit(1);
            }
        }
    }
    if (nextBufferForConsumer == -1) {
        //fprintf(stderr, "startIo thread %x reset releaseEvent\n", GetCurrentThreadId());
        ResetEvent(releaseEvent);
    }
}

    void
WindowsOverlappedDataReader::waitForBuffer(
    unsigned bufferNumber)
{
    _ASSERT(bufferNumber >= 0 && bufferNumber < nBuffers);
    BufferInfo *info = &bufferInfo[bufferNumber];
    OVERLAPPED *bufferLap = &bufferLaps[bufferNumber];

    while (info->state == InUse) {
        //fprintf(stderr, "WindowsOverlappedDataReader::waitForBuffer %d InUse...\n", bufferNumber);
        // must already have lock to call, release & wait & reacquire
        ReleaseExclusiveLock(&lock);
        _int64 start = timeInNanos();
        DWORD waitTime;
        if (releaseWaitInMillis > 0xffffffff) {
            waitTime = INFINITE;
        } else {
            waitTime = (DWORD)releaseWaitInMillis;
        }
        DWORD result = WaitForSingleObject(releaseEvent, waitTime);
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
    if (!GetOverlappedResult(hFile, bufferLap, (DWORD *)&info->validBytes,TRUE)) {
        WriteErrorMessage("Error reading FASTQ file, %d\n",GetLastError());
        soft_exit(1);
    }
    InterlockedAdd64AndReturnNewValue(&ReadWaitTime, timeInNanos() - start);

    info->state = Full;
    info->buffer[info->validBytes] = 0;
    ResetEvent(bufferLap->hEvent);
}

void
WindowsOverlappedDataReader::addBuffer()
{
    if (nBuffers == maxBuffers) {
        WriteErrorMessage("WindowsOverlappedDataReader: addBuffer at limit\n");
        return;
    }
    _ASSERT(nBuffers < maxBuffers);

    bufferLaps[nBuffers].hEvent = CreateEvent(NULL,TRUE,FALSE,NULL);
    if (NULL == bufferLaps[nBuffers].hEvent) {
        WriteErrorMessage("WindowsOverlappedDataReader: Unable to create event\n");
        soft_exit(1);
    }

    ReadBasedDataReader::addBuffer();
}

class WindowsOverlappedDataSupplier : public DataSupplier
{
public:
    WindowsOverlappedDataSupplier() : DataSupplier() {}
    virtual DataReader* getDataReader(int bufferCount, _int64 overflowBytes, double extraFactor, size_t bufferSpace)
    {
        // add some buffers for read-ahead
        return new WindowsOverlappedDataReader(bufferCount + (bufferCount > 1 ? 4 : 0), overflowBytes, extraFactor, bufferSpace);
    }
};

DataSupplier* DataSupplier::WindowsOverlapped = new WindowsOverlappedDataSupplier();

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

    DecompressDataReader(DataReader* i_inner, int i_count, _int64 totalExtra, _int64 i_extraBytes, _int64 i_overflowBytes, int i_chunkSize = BAM_BLOCK);

    virtual ~DecompressDataReader();

    virtual bool init(const char* fileName);

    virtual char* readHeader(_int64* io_headerSize);

    virtual void reinit(_int64 startingOffset, _int64 amountOfFileToProcess);

    virtual bool getData(char** o_buffer, _int64* o_validBytes, _int64* o_startBytes = NULL);

    virtual void advance(_int64 bytes);

    virtual void nextBatch();

    virtual bool isEOF();

    virtual DataBatch getBatch();

    virtual void holdBatch(DataBatch batch);

    virtual bool releaseBatch(DataBatch batch);

    virtual _int64 getFileOffset();

    virtual void getExtra(char** o_extra, _int64* o_length);

    virtual const char* getFilename()
    { return inner->getFilename(); }

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
    DataReader* i_inner,
    int i_count,
    _int64 i_totalExtra,
    _int64 i_extraBytes,
    _int64 i_overflowBytes,
    int i_chunkSize)
    : DataReader(), inner(i_inner), count(i_count), offset(i_overflowBytes),
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
        entry->batch = DataBatch(0, 0);
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
    _ASSERT(total >= totalExtra);
    _int64 headerSize = 0;
    while (headerSize < *io_headerSize && compressedBytes > 0) {
        _int64 compressedBlockSize, decompressedBlockSize;
        //fprintf(stderr,"decompress chunkSize %d compressedBytes %d headerSize %d totalExtra %d\n", chunkSize, compressedBytes, headerSize, totalExtra);
        decompress(&zstream, chunkSize != 0 ? &heap : NULL,
            compressed, compressedBytes, &compressedBlockSize,
            header + headerSize, totalExtra - headerSize, &decompressedBlockSize,
            StartMultiBlock);
        // This just gets reinit()'ed later, and in the interim confuses the non-rewind stdio data reader.  inner->advance(compressedBlockSize);
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
        WriteErrorMessage("DecompressDataReader reinit called twice\n");
        soft_exit(1);
    }
    // todo: transform start/amount to add for compression? I don't think so...
    inner->reinit(startingOffset, amountOfFileToProcess);
    threadStarted = true;
    if (! StartNewThread(chunkSize > 0 ? decompressThread : decompressThreadContinuous, this)) {
        WriteErrorMessage("failed to start decompressThread\n");
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
    _ASSERT(next->state == EntryReady && next->decompressed != NULL);
    _int64 copy = old->decompressedValid - max(offset, old->decompressedStart);
    memcpy(next->decompressed + overflowBytes - copy, old->decompressed + old->decompressedValid - copy, copy);
    offset = overflowBytes - copy;
    //fprintf(stderr,"DecompressDataReader nextBatch %d:%d #%d -> %d:%d #%d copy %lld + %lld/%lld\n", old->batch.fileID, old->batch.batchID, old-entries, next->batch.fileID, next->batch.batchID, next-entries, copy, next->decompressedStart, next->decompressedValid);
    releaseBatch(old->batch); // holdBatch was called in decompress thread, release now if no customers added holds
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
DecompressDataReader::holdBatch(DataBatch batch)
{
    inner->holdBatch(batch);
}

    bool
DecompressDataReader::releaseBatch(DataBatch batch)
{
    if (! inner->releaseBatch(batch)) {
        return false;
    }
    // truly released, find matching entry & put back on available list
    AcquireExclusiveLock(&lock);
    for (int i = 0; i < count; i++) {
        Entry* entry = &entries[i];
        if (entry->batch == batch) {
			//fprintf(stderr,"DecompressDataReader releaseBatch %d:0x%x #%d\n", batch.fileID, batch.batchID, i);
            if (entry->state == EntryHeld) {
                enqueueAvailable(entry);
            } else {
                _ASSERT(entry->state == EntryAvailable);
            }
            break;
        }
    }
    ReleaseExclusiveLock(&lock);
    return true;
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
        WriteErrorMessage("GzipDataReader: inputBytes or outputBytes > max unsigned int\n");
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
                WriteErrorMessage("GzipDataReader: inflateInit2 failed with %d\n", status);
                return false;
            }
        }
        oldAvailOut = zstream->avail_out;
        oldAvailIn = zstream->avail_in;
        status = inflate(zstream, mode == SingleBlock ? Z_NO_FLUSH : Z_FINISH);
        // fprintf(stderr, "decompress block #%d %lld -> %lld = %d\n", block, zstream.next_in - lastIn, zstream.next_out - lastOut, status);
        block++;
        if (status < 0 && status != Z_BUF_ERROR) {
            WriteErrorMessage("GzipDataReader: inflate failed with %d\n", status);
            soft_exit(1);
        }
        if (status < 0 && zstream->avail_out == 0 && zstream->avail_in > 0) {
            WriteErrorMessage("insufficient decompression buffer space - increase expansion factor, currently -xf %.1f\n", DataSupplier::ExpansionFactor);
            soft_exit(1);
        }
    } while (zstream->avail_in != 0 && (zstream->avail_out != oldAvailOut || zstream->avail_in != oldAvailIn) && mode != SingleBlock);
    // fprintf(stderr, "end decompress status=%d, avail_in=%lld, last block=%lld->%lld, avail_out=%lld\n", status, zstream.avail_in, zstream.next_in - lastIn, zstream.next_out - lastOut, zstream.avail_out);
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
    ParallelCoworker coworker(min(8, DataSupplier::ThreadCount), false, &manager);
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
            //fprintf(stderr, "decompressThread #%d %d:%d eof\n", index, reader->inner->getBatch().fileID, reader->inner->getBatch().batchID);
            if (! reader->inner->isEOF()) {
                WriteErrorMessage("error reading file at offset %lld\n", reader->getFileOffset());
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
            _int64 extraSize;
            reader->inner->getExtra(&entry->decompressed, &extraSize);
			_ASSERT(extraSize >= reader->extraBytes && extraSize >= reader->overflowBytes);
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

                if (output > reader->extraBytes) {
                    fprintf(stderr, "insufficient decompression space, increase -xf parameter\n");
                    soft_exit(1);
                }
                if (input > entry->compressedValid || zip->BSIZE() >= BAM_BLOCK || zip->ISIZE() > BAM_BLOCK) {
                    fprintf(stderr, "error reading BAM file at offset %lld\n", reader->getFileOffset());
                    soft_exit(1);
                }
            } while (input < entry->compressedStart);
            // append final offsets
            inputs.push_back(input);
            outputs.push_back(output);
            //fprintf(stderr, "decompressThread read #%d %lld->%lld\n", index, input, output);
            reader->inner->advance(input);
            entry->decompressedValid = output;
            entry->decompressedStart = output - reader->overflowBytes;
            entry->batch = reader->inner->getBatch();
			reader->holdBatch(entry->batch); // hold batch while decompressing
            reader->inner->nextBatch(); // start reading next batch
            // decompress all chunks synchronously on multiple threads
            manager.entry = entry;
            coworker.step();
        }
        // make buffer available for clients & go on to next
        //fprintf(stderr, "decompressThread #%d %d:%d ready\n", index, entry->batch.fileID, entry->batch.batchID);
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
            //fprintf(stderr, "decompressThreadContinuous#%d %d:%d eof\n", index, reader->inner->getBatch().fileID, reader->inner->getBatch().batchID);
            if (! reader->inner->isEOF()) {
                WriteErrorMessage("error reading file at offset %lld\n", reader->getFileOffset());
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
            reader->holdBatch(entry->batch); // hold batch while decompressing
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
        //fprintf(stderr, "decompressThreadContinuous#%d %d:%d ready\n", index, entry->batch.fileID, entry->batch.batchID);
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
            //fprintf(stderr, "popReady %d:%d #%d -> held\n", first->batch.fileID, first->batch.batchID, first - entries);
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
        //fprintf(stderr, "dequeueAvailable #%d\n", available == NULL ? -1 : available - entries);
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
    AssertExclusiveLockHeld(&lock);
    _ASSERT(entry->state == EntryHeld);
    entry->state = EntryAvailable;
    entry->next = available;
    available = entry;
    if (entry->next == NULL) {
        AllowEventWaitersToProceed(&availableEvent);
    }
}

class DecompressDataReaderSupplier : public DataSupplier
{
public:
    DecompressDataReaderSupplier(DataSupplier* i_inner, int i_blockSize = BAM_BLOCK)
        : DataSupplier(), inner(i_inner), blockSize(i_blockSize)
    {}

    virtual DataReader* getDataReader(int bufferCount, _int64 overflowBytes, double extraFactor, size_t bufferSpace);

private:
    DataSupplier* inner;
    const int blockSize;
};

    DataReader*
DecompressDataReaderSupplier::getDataReader(
    int bufferCount,
    _int64 overflowBytes,
    double extraFactor,
    size_t bufferSpace)
{
    // adjust extra factor for compression ratio
    double expand = MAX_FACTOR * DataSupplier::ExpansionFactor;
    double totalFactor = expand * (1.0 + extraFactor);
    // get inner reader with no overflow since zlib can't deal with it
    // add 2 buffers for compression thread
    DataReader* data = inner->getDataReader(bufferCount + 2, blockSize, totalFactor, bufferSpace);
    // compute how many extra bytes are owned by this layer
    char* p;
    _int64 totalExtra;
    data->getExtra(&p, &totalExtra);
    _int64 mine = (_int64)(totalExtra * expand / totalFactor);
    // create new reader, telling it how many bytes it owns
    // it will subtract overflow off the end of each batch
    return new DecompressDataReader(data, bufferCount, totalExtra, mine, overflowBytes, blockSize);
}
    
    DataSupplier*
DataSupplier::GzipBam(
    DataSupplier* inner)
{
    return new DecompressDataReaderSupplier(inner, BAM_BLOCK);
}
    
    DataSupplier*
DataSupplier::Gzip(
    DataSupplier* inner)
{
    return new DecompressDataReaderSupplier(inner, 0);
}
    DataSupplier* 
DataSupplier::StdioSupplier()
{
    return new StdioDataSupplier();
}
//
// MemMap
//

class MemMapDataSupplier : public DataSupplier
{
public:
    MemMapDataSupplier();

    virtual ~MemMapDataSupplier();

    virtual DataReader* getDataReader(int bufferCount, _int64 overflowBytes, double extraFactor, size_t bufferSpace);

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

    MemMapDataReader(MemMapDataSupplier* i_supplier, int i_batchCount, _int64 i_batchSize, _int64 i_overflowBytes, _int64 i_batchExtra);

    virtual ~MemMapDataReader();
    
    virtual bool init(const char* fileName);

    virtual char* readHeader(_int64* io_headerSize);

    virtual void reinit(_int64 startingOffset, _int64 amountOfFileToProcess);

    virtual bool getData(char** o_buffer, _int64* o_validBytes, _int64* o_startBytes = NULL);

    virtual void advance(_int64 bytes);

    virtual void nextBatch();

    virtual bool isEOF();

    virtual DataBatch getBatch();

    virtual void holdBatch(DataBatch batch);

    virtual bool releaseBatch(DataBatch batch);

    virtual _int64 getFileOffset();

    virtual void getExtra(char** o_extra, _int64* o_length);

    virtual const char* getFilename()
    { return fileName; }

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
    int*            extraHolds; // keeps hold count for each extra buffer
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
 

MemMapDataReader::MemMapDataReader(MemMapDataSupplier* i_supplier, int i_batchCount, _int64 i_batchSize, _int64 i_overflowBytes, _int64 i_batchExtra)
    : DataReader(),
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
        memset(extraBatches, 0, batchCount * sizeof(DataBatch));
        extraHolds = new int[batchCount];
        memset(extraHolds, 0, batchCount * sizeof(int));
    } else {
        extra = NULL;
        extraBatches = NULL;
    }
    if (! (CreateSingleWaiterObject(&waiter) && InitializeExclusiveLock(&lock))) {
        WriteErrorMessage("MemMapDataReader: CreateSingleWaiterObject failed\n");
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
        WriteErrorMessage("MemMapDataReader: fail to map %s at %lld,%lld\n", fileName, i_startingOffset, amountOfFileToProcess);
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
        memset(extraHolds, 0, sizeof(int) * batchCount);
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
                //fprintf(stderr, "MemMap nextBatch %d:%d = index %d used %d of %d\n", 0, currentBatch, currentExtraIndex, extraUsed, batchCount); 
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
MemMapDataReader::holdBatch(
    DataBatch batch)
{
    if (extraBatches == NULL) {
        return;
    }
    acquireLock();
    for (int i = 0; i < batchCount; i++) {
        if (extraBatches[i] == batch) {
            extraHolds[i]++;
            break;
        }
    }
    releaseLock();
}

    bool
MemMapDataReader::releaseBatch(
    DataBatch batch)
{
    if (extraBatches == NULL) {
        return true;
    }
    bool result = true;
    acquireLock();
    for (int i = 0; i < batchCount; i++) {
        if (extraBatches[i] == batch) {
            if (extraHolds[i] > 0) {
                extraHolds[i]--;
            }
            if (extraHolds[i] == 0) {
                extraBatches[i].batchID = 0;
                _ASSERT(extraUsed > 0);
                extraUsed--;
                //fprintf(stderr,"MemMap: releaseBatch %d:%d = index %d now using %d of %d\n", batch.fileID, batch.batchID, i, extraUsed, batchCount);
                if (extraUsed == batchCount - 1) {
                    SignalSingleWaiterObject(&waiter);
                }
                releaseLock();
            } else {
                result = false;
            }
            break;
        }
    }
    releaseLock();
    return result;
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

MemMapDataSupplier::MemMapDataSupplier() : DataSupplier()
{
    InitializeExclusiveLock(&lock);
}

MemMapDataSupplier::~MemMapDataSupplier()
{
    DestroyExclusiveLock(&lock);
}

    DataReader* 
MemMapDataSupplier::getDataReader(
    int bufferCount,
    _int64 overflowBytes,
    double extraFactor,
    size_t bufferSpace /*not relevant*/)
{
    _ASSERT(extraFactor >= 0 && overflowBytes >= 0);
    if (extraFactor == 0) {
        // no per-batch expansion factor, so can read entire file as a batch
        return new MemMapDataReader(this, 1, 0, overflowBytes, 0);
    } else {
        // break up into 4Mb batches
        _int64 batch = 4 * 1024 * 1024;
        _int64 extra = (_int64)(batch * extraFactor);
        return new MemMapDataReader(this, bufferCount, batch, overflowBytes, extra);
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

    bool
BatchTracker::holdBatch(
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
    //fprintf(stderr, "thread %d tracker %lx addRead %u:%u = %d\n", GetCurrentThreadId(), this, batch.fileID, batch.batchID, n);
    return p == NULL;
}

    bool
BatchTracker::releaseBatch(
    DataBatch removed)
{
    DataBatch::Key key = removed.asKey();
    unsigned* p = pending.tryFind(key);
    //fprintf(stderr, "thread %d tracker %lx removeRead %u:%u = %d\n", GetCurrentThreadId(), this, removed.fileID, removed.batchID, p ? *p - 1 : -1);
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

DataSupplier* DataSupplier::MemMap = new MemMapDataSupplier();

#ifdef _MSC_VER
DataSupplier* DataSupplier::Default = DataSupplier::WindowsOverlapped;
#else
DataSupplier* DataSupplier::Default = DataSupplier::MemMap;
#endif

DataSupplier* DataSupplier::GzipDefault = DataSupplier::Gzip(DataSupplier::Default);

DataSupplier* DataSupplier::GzipBamDefault = DataSupplier::GzipBam(DataSupplier::Default);

DataSupplier* DataSupplier::Stdio = DataSupplier::StdioSupplier();

DataSupplier* DataSupplier::GzipStdio = DataSupplier::Gzip(DataSupplier::Stdio);

DataSupplier* DataSupplier::GzipBamStdio = DataSupplier::GzipBam(DataSupplier::Stdio);


int DataSupplier::ThreadCount = 1;

double DataSupplier::ExpansionFactor = 1.0;

volatile _int64 DataReader::ReadWaitTime = 0;
volatile _int64 DataReader::ReleaseWaitTime = 0;
