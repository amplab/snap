/*++

Module Name:

    GzipDataWriter.h

Abstract:

    Headers for the GzipDataWriter & related classes for the SNAP sequencer

Authors:

    Ravi Pandya, Mar 2013

Environment:

    User mode service.

Revision History:

--*/

#pragma once

#include "Compat.h"
#include "Read.h"
#include "Compat.h"
#include "DataWriter.h"
#include "VariableSizeVector.h"
#include "stdafx.h"
#include "zlib.h"

using std::pair;

class GzipWriterFilterSupplier;

class GzipWriterFilter : public DataWriter::Filter
{
public:
    GzipWriterFilter(GzipWriterFilterSupplier* i_supplier);
    
    virtual ~GzipWriterFilter() {}

    virtual void onAdvance(DataWriter* writer, size_t batchOffset, char* data, size_t bytes, unsigned location);

    virtual size_t onNextBatch(DataWriter* writer, size_t offset, size_t bytes);

private:
    size_t compressChunk(char* toBuffer, size_t toSize, char* fromBuffer, size_t fromUsed);

    GzipWriterFilterSupplier* supplier;
    const size_t chunkSize;
    const bool bamFormat;
    z_stream zstream;
};

class GzipWriterFilterSupplier : public DataWriter::FilterSupplier
{
public:
    GzipWriterFilterSupplier(bool i_bamFormat, size_t i_chunkSize)
    :
        FilterSupplier(DataWriter::TransformFilter),
        bamFormat(i_bamFormat),
        chunkSize(i_chunkSize)
    {}

    virtual ~GzipWriterFilterSupplier()
    {
    }

    virtual DataWriter::Filter* getFilter()
    {
        return new GzipWriterFilter(this);
    }

    virtual void onClose(DataWriterSupplier* supplier, DataWriter* writer);
    
    bool translate(_uint64 logical, _uint64* o_physical, _uint64* delta);

    // translate to BAM virtual offset format
    _uint64 toVirtualOffset(_uint64 logical)
    {
        _uint64 physical, delta;
        if (logical != UINT64_MAX && translate(logical, &physical, &delta)) {
            if (delta < 65536 && physical < ((_uint64) 1 << 48)) {
                return (physical << 48) | delta;
            }
            fprintf(stderr, "Invalid virtual file offset, logical=%llu, start=%llu, delta=%llu\n", logical, physical, delta);
        }
        return 0;
    }

private:
    friend class GzipWriterFilter;

    void addTranslation(_uint64 logical, _uint64 physical)
    { translation.push_back(pair<_uint64,_uint64>(logical, physical)); }

    static bool translationComparator(const pair<_uint64,_uint64>& a, const pair<_uint64,_uint64>& b);

    const bool bamFormat;
    const size_t chunkSize;
    VariableSizeVector<pair<_uint64,_uint64>> translation;
};