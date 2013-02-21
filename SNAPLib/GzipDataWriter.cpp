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

class GzipWriterFilter : public DataWriter::Filter
{
public:
    GzipWriterFilter(bool i_bamFormat, size_t i_chunkSize, int i_headerChunks, DataWriter::Filter* i_inner)
        :
        Filter(DataWriter::TransformFilter),
        bamFormat(i_bamFormat),
        chunkSize(i_chunkSize),
        headerChunks(i_headerChunks),
        chunkStart(0),
        nextUsed(0),
        inner(i_inner)
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
    void GzipWriterFilter::compressChunk(DataWriter* writer, int relative, size_t ignore);

    int headerChunks;
    const size_t chunkSize;
    const bool bamFormat;
    size_t chunkStart;
    size_t nextUsed;
    z_stream zstream;
    DataWriter::Filter* inner;
};

    void
GzipWriterFilter::onAdvance(
    DataWriter* writer,
    size_t batchOffset,
    char* data,
    size_t bytes,
    unsigned location)
{
    if (headerChunks > 0) {
        headerChunks--;
        compressChunk(writer, 0, 0);
        chunkStart = batchOffset + bytes;
    } else if (batchOffset + bytes - chunkStart > chunkSize) {
        compressChunk(writer, 0, bytes);
        chunkStart = batchOffset;
    }
}

    size_t
GzipWriterFilter::onNextBatch(
    DataWriter* writer,
    size_t offset,
    size_t bytes)
{
    compressChunk(writer, -1, 0);
    size_t result = nextUsed;
    nextUsed = 0;
    chunkStart = 0;
    if (inner != NULL) {
        inner->onNextBatch(NULL, offset, result);
    }
    return result;
}


    void
GzipWriterFilter::compressChunk(
    DataWriter* writer,
    int relative,
    size_t ignore)
{
    char* fromBuffer;
    size_t fromSize, fromUsed;
    writer->getBatch(relative, &fromBuffer, &fromSize, &fromUsed);
    if (fromUsed - chunkStart <= ignore) {
        return;
    }
    int count = (int)(fromUsed - chunkStart - ignore);
    char* toBuffer;
    size_t toSize, toUsed;
    writer->getBatch(relative + 1, &toBuffer, &toSize, &toUsed);

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
        bamExtraData[4] = 0; // will be filled in later
        bamExtraData[5] = 0; // will be filled in later
    }

    // based on sample code at http://www.lemoda.net/c/zlib-open-write/index.html
    const int windowBits = 15;
    const int GZIP_ENCODING = 16;
    zstream.next_in = (Bytef*) (fromBuffer + chunkStart);
    zstream.avail_in = count;
    zstream.next_out = (Bytef*) (toBuffer + nextUsed);
    zstream.avail_out = (uInt)(toSize - nextUsed);
    uInt oldAvail;
    int status;

    status = deflateInit2(&zstream, Z_DEFAULT_COMPRESSION, Z_DEFLATED, windowBits | GZIP_ENCODING, 8, Z_DEFAULT_STRATEGY);
    if (status < 0) {
        fprintf(stderr, "GzipWriterFilter: deflateInit2 failed with %d\n", status);
        exit(1);
    }
    if (bamFormat) {
        status = deflateSetHeader(&zstream, &header);
        if (status != Z_OK) {
            fprintf(stderr, "GzipWriterFilter: defaultSetHeader failed with %d\n", status);
            exit(1);
        }
    }
    oldAvail = zstream.avail_out;
    status = deflate(&zstream, Z_FINISH);
    if (status < 0 && status != Z_BUF_ERROR) {
        fprintf(stderr, "GzipWriterFilter: deflate failed with %d\n", status);
        exit(1);
    }
    
    // make sure it all got written out in a single compressed block
    if (zstream.avail_in != 0) {
        fprintf(stderr, "GzipWriterFilter: default failed to read all input\n");
        exit(1);
    }
    if (zstream.avail_out == oldAvail) {
        fprintf(stderr, "GzipWriterFilter: default failed to write output\n");
        exit(1);
    }
    status = deflateEnd(&zstream);
    if (status < 0) {
        fprintf(stderr, "GzipWriterFilter: deflateEnd failed with %d\n", status);
        exit(1);
    }

    size_t oldNextUsed = nextUsed;
    nextUsed = toSize - zstream.avail_out;
    if (bamFormat) {
        // backpatch compressed block size into gzip header
        if (nextUsed - oldNextUsed >= 0x10000) {
            fprintf(stderr, "exceeded BAM chunk size\n");
            exit(1);
        }
        * (_uint16*) (toBuffer + oldNextUsed + 16) = (_uint16) (nextUsed - oldNextUsed - 1);
    }
    if (inner != NULL) {
        inner->onAdvance(NULL, oldNextUsed, toBuffer + oldNextUsed, nextUsed - oldNextUsed, 0);
    }
}
    
class GzipWriterFilterSupplier : public DataWriter::FilterSupplier
{
public:
    GzipWriterFilterSupplier(bool i_bamFormat, size_t i_chunkSize,
        int i_headerChunks, DataWriter::FilterSupplier* i_inner)
    :
        FilterSupplier(DataWriter::TransformFilter),
        bamFormat(i_bamFormat),
        chunkSize(i_chunkSize),
        headerChunks(i_headerChunks),
        inner(i_inner)
    {}

    virtual ~GzipWriterFilterSupplier()
    {
    }

    virtual DataWriter::Filter* getFilter()
    {
        int n = headerChunks;
        headerChunks = 0;
        return new GzipWriterFilter(bamFormat, chunkSize, n, inner != NULL ? inner->getFilter() : NULL);
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
                exit(1);
            }
            memcpy(buffer, eof, sizeof(eof));
            writer->advance(sizeof(eof));
        }
        if (inner != NULL) {
            inner->onClose(supplier, writer);
        }
    }

private:
    const bool bamFormat;
    const size_t chunkSize;
    int headerChunks;
    DataWriter::FilterSupplier* inner;
};

    DataWriter::FilterSupplier*
DataWriterSupplier::gzip(
    bool bamFormat,
    size_t chunkSize,
    int headerChunks,
    DataWriter::FilterSupplier* inner)
{
    return new GzipWriterFilterSupplier(bamFormat, chunkSize, headerChunks, inner);
}
