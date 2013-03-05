/*++

Module Name:

    ZipDataWriter.cpp

Abstract:

    File writer that compresses data into zip format.

Environment:

    User mode service.

    Not thread safe.

--*/

#include "stdafx.h"
#include "BigAlloc.h"
#include "Compat.h"
#include "DataWriter.h"
#include "VariableSizeVector.h"
#include "zlib.h"
#include "exit.h"

class GzipWriterFilter : public DataWriter::Filter
{
public:
    GzipWriterFilter(bool i_bamFormat, size_t i_chunkSize) :
        Filter(DataWriter::TransformFilter),
        bamFormat(i_bamFormat),
        chunkSize(i_chunkSize)
    {
        zstream.zalloc = Z_NULL;
        zstream.zfree = Z_NULL;
        zstream.opaque = Z_NULL;
    }
    
    virtual ~GzipWriterFilter()
    {
    }

    virtual void onAdvance(DataWriter* writer, size_t batchOffset, char* data, size_t bytes, unsigned location);

    virtual size_t onNextBatch(DataWriter* writer, size_t offset, size_t bytes);

private:
    size_t compressChunk(char* toBuffer, size_t toSize, char* fromBuffer, size_t fromUsed);

    const size_t chunkSize;
    const bool bamFormat;
    z_stream zstream;
};

    void
GzipWriterFilter::onAdvance(
    DataWriter* writer,
    size_t batchOffset,
    char* data,
    size_t bytes,
    unsigned location)
{
    // nothing
}

    size_t
GzipWriterFilter::onNextBatch(
    DataWriter* writer,
    size_t offset,
    size_t bytes)
{
    char* fromBuffer;
    size_t fromSize, fromUsed;
    writer->getBatch(-1, &fromBuffer, &fromSize, &fromUsed);
    if (fromUsed == 0) {
        return 0;
    }
    char* toBuffer;
    size_t toSize, toUsed;
    writer->getBatch(0, &toBuffer, &toSize, &toUsed);

    for (size_t chunk = 0; chunk < fromUsed; chunk += chunkSize) {
        toUsed += compressChunk(toBuffer + toUsed, toSize - toUsed, fromBuffer + chunk, std::min(fromUsed - chunk, chunkSize));
    }
    return toUsed;
}


    size_t
GzipWriterFilter::compressChunk(
    char* toBuffer,
    size_t toSize,
    char* fromBuffer,
    size_t fromUsed)
{
    if (bamFormat && fromUsed > 0x10000) {
        fprintf(stderr, "exceeded BAM chunk size\n");
        exit(1);
    }
    // set up BAM header structure
    gz_header header;
    _uint8 bamExtraData[6];
    if (bamFormat) {
        header.text = false;
        header.time = 0;
        header.xflags = 0;
        header.os = 0;
        header.extra = bamExtraData;
        header.extra_len = 6;
        header.extra_max = 6;
        header.name = NULL;
        header.name_max = 0;
        header.comment = NULL;
        header.comm_max = 0;
        header.hcrc = false;
        header.done = true;
        bamExtraData[0] = 'B';
        bamExtraData[1] = 'C';
        bamExtraData[2] = 2;
        bamExtraData[3] = 0;
        bamExtraData[4] = 3; // will be filled in later
        bamExtraData[5] = 7; // will be filled in later
    }

    // based on sample code at http://www.lemoda.net/c/zlib-open-write/index.html
    const int windowBits = 15;
    const int GZIP_ENCODING = 16;
    zstream.next_in = (Bytef*) fromBuffer;
    zstream.avail_in = fromUsed;
    zstream.next_out = (Bytef*) toBuffer;
    zstream.avail_out = toSize;
    uInt oldAvail;
    int status;

    status = deflateInit2(&zstream, Z_DEFAULT_COMPRESSION, Z_DEFLATED, windowBits | GZIP_ENCODING, 8, Z_DEFAULT_STRATEGY);
    if (status < 0) {
        fprintf(stderr, "GzipWriterFilter: deflateInit2 failed with %d\n", status);
        soft_exit(1);
    }
    if (bamFormat) {
        status = deflateSetHeader(&zstream, &header);
        if (status != Z_OK) {
            fprintf(stderr, "GzipWriterFilter: defaultSetHeader failed with %d\n", status);
            soft_exit(1);
        }
    }
    oldAvail = zstream.avail_out;
    status = deflate(&zstream, Z_FINISH);
    if (status < 0 && status != Z_BUF_ERROR) {
        fprintf(stderr, "GzipWriterFilter: deflate failed with %d\n", status);
        soft_exit(1);
    }
    
    // make sure it all got written out in a single compressed block
    if (zstream.avail_in != 0) {
        fprintf(stderr, "GzipWriterFilter: default failed to read all input\n");
        soft_exit(1);
    }
    if (zstream.avail_out == oldAvail) {
        fprintf(stderr, "GzipWriterFilter: default failed to write output\n");
        soft_exit(1);
    }
    status = deflateEnd(&zstream);
    if (status < 0) {
        fprintf(stderr, "GzipWriterFilter: deflateEnd failed with %d\n", status);
        soft_exit(1);
    }

    size_t toUsed = toSize - zstream.avail_out;
    if (bamFormat) {
        // backpatch compressed block size into gzip header
        if (toUsed >= 0x10000) {
            fprintf(stderr, "exceeded BAM chunk size\n");
            soft_exit(1);
        }
        * (_uint16*) (toBuffer + 16) = (_uint16) (toUsed - 1);
    }
    return toUsed;
}
    
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
        return new GzipWriterFilter(bamFormat, chunkSize);
    }

    // called when entire file is done, all Filters destroyed
    virtual void onClose(DataWriterSupplier* supplier, DataWriter* writer)
    {
        if (bamFormat) {
            // write empty block as BAM end of file marker
            static _uint8 eof[] = {
                0x1f, 0x8b, 0x08, 0x04, 0x00, 0x00, 0x00, 0x00, 0x00, 0xff, 0x06, 0x00, 0x42, 0x43,
                0x02, 0x00, 0x1b, 0x00, 0x03, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00
            };
            char* buffer;
            size_t bytes;
            if (! (writer->getBuffer(&buffer, &bytes) && bytes >= sizeof(eof))) {
                fprintf(stderr, "no space to write eof marker\n");
                soft_exit(1);
            }
            memcpy(buffer, eof, sizeof(eof));
            writer->advance(sizeof(eof));
        }
    }

private:
    const bool bamFormat;
    const size_t chunkSize;
};

    DataWriter::FilterSupplier*
DataWriterSupplier::gzip(
    bool bamFormat,
    size_t chunkSize)
{
    return new GzipWriterFilterSupplier(bamFormat, chunkSize);
}
