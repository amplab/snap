/*++

Module Name:

    DataReader.h

Abstract:

    Headers for the DataReader & related classes for the SNAP sequencer

Authors:

    Ravi Pandya, Jan 2013

Environment:

    User mode service.

Revision History:

--*/

#pragma once

#include "Compat.h"
#include "VariableSizeMap.h"
//
// This defines a family of composable classes for efficiently reading data with flow control.
//
// DataReader
//      Reads data from one or more files either sequentially, in ranges, or memory-mapped.
//      A DataReader should be accessed by only one thread at a time,
//      except for release() which may be called from any thread.
//      Divides data into sequential batches each of which is identified by a file ID and batch ID.
//      Data in a batch will remain stable until it is released by the consumer.
//      Consumers should release batches as soon as possible to make buffers free for read-ahead.
//      Batches may include extra data for higher layers that also remains stable.
//      Extra data size is defined as a factor of the underlying data size, and/or a fixed number of bytes.
//
// DataSupplier
//      A factory for DataReaders, which may be called from multiple threads.
//

struct DataBatch
{
    _uint32   fileID;
    _uint32   batchID;

    inline DataBatch() : fileID(0), batchID(0) {}

    inline DataBatch(_uint32 i_batchID, _uint32 i_fileID = 0) : fileID(i_fileID), batchID(i_batchID) {}

    inline DataBatch(const DataBatch& o) : fileID(o.fileID), batchID(o.batchID) {}
    
    static bool comparator(const DataBatch& a, const DataBatch& b)
    { return a.fileID < b.fileID || (a.fileID == b.fileID && a.batchID < b.batchID); }

    inline bool operator<=(const DataBatch& b) const
    { return fileID < b.fileID || (fileID == b.fileID && batchID <= b.batchID); }
    
    inline bool operator<(const DataBatch& b) const
    { return fileID < b.fileID || (fileID == b.fileID && batchID < b.batchID); }
    
    inline bool operator==(const DataBatch& b) const
    { return fileID == b.fileID && batchID == b.batchID; }
    
    inline bool operator!=(const DataBatch& b) const
    { return batchID != b.batchID || fileID != b.fileID; }

    inline DataBatch Min(const DataBatch& b) const
    { return *this <= b ? *this : b; }

    inline bool isZero() const
    { return fileID == 0 && batchID == 0; }

    // convert to _int64 for use as a hashtable key

    typedef _int64 Key;

    inline Key asKey()
    { return (((_int64) fileID) << 32) + (_int64) batchID; }

    inline DataBatch(Key key) : fileID((_uint32) (key >> 32)), batchID((_uint32) key) {}
};

// read data from a file or other source
// should all be called from a single thread, except for hold/releaseBatch which are thread-safe
class DataReader
{
public:
    DataReader() {}

    virtual ~DataReader() {}
    
    // initialize to use a specific filename
    virtual bool init(const char* fileName) = 0;

    // read bytes from the beginning of the file for the header
    virtual char* readHeader(_int64* io_headerSize) = 0;

    // seek to a particular range in the file
    virtual void reinit(_int64 startingOffset, _int64 amountOfFileToProcess) = 0;

    // get all remaining data in current batch
    // return false if no more data in current batch
    // startBytes is data "owned" by this block in which reads may start
    // validBytes may also include overflow bytes to handle records spanning batches
    // if you advance() past startBytes, nextBatch() will start offset at that point
    virtual bool getData(char** o_buffer, _int64* o_validBytes, _int64* o_startBytes = NULL) = 0;

    // advance through data in current batch, reducing results from next getData call
    virtual void advance(_int64 bytes) = 0;

    // advance to next batch
    // by default automatically releases previous batch
    virtual void nextBatch() = 0;

    // whether current batch is last in file
    virtual bool isEOF() = 0;

    // get current batch identifier
    virtual DataBatch getBatch() = 0;
    
    // hold buffers associated with this batch for reuse, increments refcount
    // NOTE: this may be called from another thread,
    // so anything it touches must be thread-safe!
    virtual void holdBatch(DataBatch batch) = 0;

    // release buffers associated with this batch for reuse
    // decrements refcount, returns true if last release
    // NOTE: this may be called from another thread,
    // so anything it touches must be thread-safe!
    virtual bool releaseBatch(DataBatch batch) = 0;

    // get current offset into file
    virtual _int64 getFileOffset() = 0;

    // get pointer to extra data area for current batch
    // todo: allow this to grow dynamically while keeping stable pointers to previous data
    virtual void getExtra(char** o_extra, _int64* o_length) = 0;

    // get filename for debugging / error printing
    virtual const char* getFilename() = 0;

    // timing for performance tuning (in nanos)
    static volatile _int64 ReadWaitTime;
    static volatile _int64 ReleaseWaitTime;
};

class DataSupplier
{
public:
    DataSupplier() {}
    
    virtual ~DataSupplier() {}

    virtual DataReader* getDataReader(int bufferCount, _int64 overflowBytes, double extraFactor, size_t bufferSpace) = 0;

    //
    // creating specific factories
    //

    // 
    static DataSupplier* GzipBam(DataSupplier* inner);
    static DataSupplier* Gzip(DataSupplier* inner);
    static DataSupplier* StdioSupplier();

    // memmap works on both platforms (but better on Linux)
    static DataSupplier* MemMap;

#ifdef _MSC_VER
    // overlapped is only on Windows
    static DataSupplier* WindowsOverlapped;
#endif

    // default raw data supplier for platform
    static DataSupplier* Default;
    static DataSupplier* GzipDefault;
    static DataSupplier* GzipBamDefault;

    static DataSupplier* GzipStdio;
    static DataSupplier* Stdio;
    static DataSupplier* GzipBamStdio;

    // hack: must be set to communicate thread count into suppliers
    static int ThreadCount;

    // hack: global for additional expansion factor
    static double ExpansionFactor;
};

// manages lifetime tracking for batches of reads
class BatchTracker
{
public:
    BatchTracker(int i_capacity);

    // reference was added from a batch, increment reference count
    // return true if first hold
    bool holdBatch(DataBatch batch);

    // reference was removed from a batch
    // returns true if the batch has no more references
    bool releaseBatch(DataBatch batch);

private:
    typedef VariableSizeMap<DataBatch::Key,unsigned> BatchMap;
    BatchMap pending;
};

