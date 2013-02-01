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

//
// This defines a family of composable classes for reading data.
//
// DataReader
//      Reads data from one or more files either sequentially, in ranges, or memory-mapped.
//      A DataReader should be accessed by only one thread at a time.
//      Divides data into sequential batches each of which is identified by a file ID and batch ID.
//      Data in a batch will remain stable until it is released by the consumer.
//      Consumers should release batches as soon as possible to make buffers free for read-ahead.
//      Batches may include extra data for higher layers that also remains stable.
//      Extra data size is defined as a factor of the underlying data size.
//
// DataSupplier
//      A factory for DataReaders.
//

struct DataBatch
{
    const _uint32   fileID;
    const _uint32   batchID;

    DataBatch() : fileID(0), batchID(0) {}

    DataBatch(_uint32 i_batchID, _uint32 i_fileID = 0) : fileID(i_fileID), batchID(i_batchID) {}

    DataBatch(const DataBatch& o) : fileID(o.fileID), batchID(o.fileID) {}
    
    static bool comparator(const DataBatch& a, const DataBatch& b)
    { return a.fileID < b.fileID || (a.fileID == b.fileID && a.batchID < b.batchID); }

    bool operator<=(const DataBatch& b)
    { return fileID < b.fileID || (fileID == b.fileID && batchID <= b.batchID); }
    
    bool operator==(const DataBatch& b)
    { return fileID == b.fileID && batchID == b.batchID; }
};

class DataReader
{
public:

    virtual ~DataReader() {}
    
    // initialize to use a specific filename
    // extra data is a factor of buffer size
    virtual bool init(const char* fileName) = 0;

    // read bytes from the beginning of the file for the header
    virtual char* readHeader(_int64* io_headerSize) = 0;

    // seek to a particular location -- todo: keep?
    virtual void reinit(_int64 startingOffset, _int64 amountOfFileToProcess) = 0;

    // get all remaining data in current batch
    // return false if no more data in current batch
    virtual bool getData(char** o_buffer, _int64* o_validBytes) = 0;

    // advance through data in current batch, reducing results from next getData call
    virtual void advance(_int64 bytes) = 0;

    // advance to next batch
    // by default automatically releases previous batch
    virtual void nextBatch(bool dontRelease = false) = 0;

    // whether current batch is last in file
    virtual bool isEOF() = 0;

    // get current batch identifier
    virtual DataBatch getBatch() = 0;

    // release all batches up to and including the given batch
    virtual void release(DataBatch batch) = 0;

    // get current offset into file
    virtual _int64 getFileOffset() = 0;

    // get pointer to extra data area for current batch
    virtual void getExtra(char** o_extra, _int64* o_length) = 0;
};

class DataSupplier
{
public:
    
    virtual ~DataSupplier() {}

    virtual DataReader* getDataReader(double extraFactor = 0.0, const DataSupplier* inner = NULL) const = 0;

    static const DataSupplier* Gzip;

#ifdef _MSC_VER
    static const DataSupplier* WindowsOverlapped;
#endif
};
