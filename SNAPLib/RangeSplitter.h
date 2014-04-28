/*++

Module Name:

    RangeSplitter.h

Abstract:

    Headers for code to split a range into pieces for multiple cores to process.  It's designed
    to handle cores that proceed at varying rates.

Authors:

    Bill Bolosky, 2011

Environment:

    User mode service.

Revision History:

    Pulled out of cSNAP.cpp to make it useful for various versions of the aligner
    Generalized from FileSplitter to ranges, e.g. for scanning the genome /ravip/5/2012/

--*/

#pragma once

#include "Compat.h"
#include "Read.h"
#include "Genome.h"
#include "AlignerOptions.h"

//
// Utility class for letting multiple threads split chunks of a range to process.
// This is used by the parallel versions of the aligners.
//
class RangeSplitter
{
public:
    RangeSplitter(_int64 rangeEnd_, int numThreads_, unsigned divisonSize_ = 5, _int64 rangeBegin_ = 0, unsigned minMillis_ = 200, unsigned minRangeSize_ = 32768);

    // Get the next range for a thread to process, or return false if the whole range is done.
    bool getNextRange(_int64 *rangeStart, _int64 *rangeLength);

private:
    int numThreads;
    _int64 rangeBegin;
    _int64 rangeEnd;
    unsigned divisionSize;
    unsigned minMillis;
    unsigned minRangeSize;
    volatile _int64 position;
    volatile _int64 startTime;
};

class RangeSplittingReadSupplier : public ReadSupplier {
public:
    RangeSplittingReadSupplier(RangeSplitter *i_splitter, ReadReader *i_underlyingReader) : 
      splitter(i_splitter), underlyingReader(i_underlyingReader), read() {}

    virtual ~RangeSplittingReadSupplier();

    Read *getNextRead();
 
    virtual void holdBatch(DataBatch batch)
    { underlyingReader->holdBatch(batch); }
    
    virtual bool releaseBatch(DataBatch batch)
    { return underlyingReader->releaseBatch(batch); }

private:
    RangeSplitter *splitter;
    ReadReader *underlyingReader;
    Read read;
};

class RangeSplittingReadSupplierGenerator: public ReadSupplierGenerator {
public:
    RangeSplittingReadSupplierGenerator(const char *i_fileName, bool i_isSAM, unsigned numThreads, const ReaderContext& context);
    ~RangeSplittingReadSupplierGenerator() {delete splitter; delete [] fileName;}

    ReadSupplier *generateNewReadSupplier();
    ReaderContext* getContext() { return &context; }

private:
    RangeSplitter *splitter;
    char *fileName;
    const bool isSAM;
    const int numThreads;
    ReaderContext context;
};


class RangeSplittingPairedReadSupplier : public PairedReadSupplier {
public:
    RangeSplittingPairedReadSupplier(RangeSplitter *i_splitter, PairedReadReader *i_underlyingReader) : splitter(i_splitter), underlyingReader(i_underlyingReader) {}
    virtual ~RangeSplittingPairedReadSupplier();

    virtual bool getNextReadPair(Read **read1, Read **read2);
       
    virtual void holdBatch(DataBatch batch)
    { underlyingReader->releaseBatch(batch); }

    virtual bool releaseBatch(DataBatch batch)
    { return underlyingReader->releaseBatch(batch); }

 private:
    PairedReadReader *underlyingReader;
    RangeSplitter *splitter;
    Read internalRead1;
    Read internalRead2;
 };

class RangeSplittingPairedReadSupplierGenerator: public PairedReadSupplierGenerator {
public:
    RangeSplittingPairedReadSupplierGenerator(const char *i_fileName1, const char *i_fileName2, enum FileType i_fileType, unsigned numThreads, bool i_quicklyDropUnpairedReads, const ReaderContext& context);
    ~RangeSplittingPairedReadSupplierGenerator();

    PairedReadSupplier *generateNewPairedReadSupplier();
    ReaderContext* getContext() { return &context; }

private:
    RangeSplitter *splitter;
    char *fileName1;
    char *fileName2;
    const int numThreads;
    enum FileType fileType;
    ReaderContext context;
    bool quicklyDropUnpairedReads;
};

