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
#include "GzipDataWriter.h"
#include "BigAlloc.h"
#include "VariableSizeVector.h"
#include "ParallelTask.h"
#include "RangeSplitter.h"
#include "zlib.h"
#include "exit.h"

using std::min;
using std::max;
using std::pair;

class GzipSharedContext
{
public:
    GzipSharedContext(bool i_bam, size_t i_chunkSize, int i_numThreads, bool i_bindToProcessors)
        : stop(false), bam(i_bam), numThreads(i_numThreads),
        nChunks(0), chunkSize(i_chunkSize), bindToProcessors(i_bindToProcessors)
    {
        workReady = new EventObject[numThreads];
        workDone = new EventObject[numThreads];
        for (int i = 0; i < numThreads; i++) {
            CreateEventObject(&workReady[i]);
            PreventEventWaitersFromProceeding(&workReady[i]);
            CreateEventObject(&workDone[i]);
            PreventEventWaitersFromProceeding(&workDone[i]);
        }
    }

    ~GzipSharedContext()
    {
        for (int i = 0; i < numThreads; i++) {
            DestroyEventObject(&workReady[i]);
            DestroyEventObject(&workDone[i]);
        }
        delete [] workReady;
        delete [] workDone;
    }

    size_t* run(int i_nChunks, char* i_input, char* i_output)
    {
        nChunks = i_nChunks;
        sizes.clear();
        sizes.extend(nChunks);
 
        input = i_input;
        output = i_output;
        for (int i = 0; i < numThreads; i++) {
            PreventEventWaitersFromProceeding(&workDone[i]);
            AllowEventWaitersToProceed(&workReady[i]);
        }

        //
        // Let the workers run and wait for them all to finish.
        //

        for (int i = 0; i < numThreads; i++) {
            WaitForEvent(&workDone[i]);
        }

        return sizes.begin();
    }

    void done()
    {
        stop = true;
        for (int i = 0; i < numThreads; i++) {
            AllowEventWaitersToProceed(&workReady[i]);
        }
    }


    EventObject *workReady; // One per worker thread
    EventObject *workDone;  // One per worker thread
    volatile bool stop;
    VariableSizeVector<size_t> sizes;
    volatile int nChunks;
    const size_t chunkSize;
    const bool bam;
    const int numThreads;
    const bool bindToProcessors;
    char* input;
    char* output;
};

struct GzipContext : public TaskContextBase
{
    GzipContext()
    {
#ifdef  _MSC_VER
        useTimingBarrier = false;
#endif
    }
    GzipSharedContext* shared;
    void initializeThread() {}
    void runThread();
    void finishThread(GzipContext* common) {}
};

class GzipWriterFilterSupplier;

class GzipWriterFilter : public DataWriter::Filter
{
public:
    GzipWriterFilter(GzipWriterFilterSupplier* i_supplier);
    
    virtual ~GzipWriterFilter()
    {
        shared.done();
        delete task;
    }

    virtual void onAdvance(DataWriter* writer, size_t batchOffset, char* data, unsigned bytes, unsigned location);

    virtual size_t onNextBatch(DataWriter* writer, size_t offset, size_t bytes);

    static size_t compressChunk(z_stream& zstream, bool bamFormat, char* toBuffer, size_t toSize, char* fromBuffer, size_t fromUsed);

private:

    GzipWriterFilterSupplier* supplier;
    const size_t chunkSize;
    const bool bamFormat;
    z_stream zstream;
    ThreadHeap heap;
    GzipSharedContext shared;
    GzipContext context;
    ParallelTask<GzipContext>* task;
};

    void
GzipContext::runThread()
{
    ThreadHeap heap(shared->chunkSize * 5); // seems to use 4x chunk size per run
    z_stream zstream;
    zstream.zalloc = zalloc;
    zstream.zfree = zfree;
    zstream.opaque = &heap;
    while (true) {
        //printf("zip task thread %d waiting to begin\n", GetCurrentThreadId());
        WaitForEvent(&shared->workReady[threadNum]);
        PreventEventWaitersFromProceeding(&shared->workReady[threadNum]);
        if (shared->stop) {
            return;
        }
        //printf("zip task thread %d begin\n", GetCurrentThreadId());
        _int64 start = timeInMillis();
        int n = (shared->nChunks + shared->numThreads - 1) / shared->numThreads;
        int begin = threadNum * n;
        int end = min(begin + n, (int) shared->nChunks);
        for (int i = begin; i < end; i++) {
            shared->sizes[i] = GzipWriterFilter::compressChunk(zstream, shared->bam,
                shared->output + i * shared->chunkSize, shared->chunkSize, 
                shared->input + i * shared->chunkSize, shared->chunkSize);
        }
        //printf("zip task thread %d done %lld ms\n", GetCurrentThreadId(), timeInMillis() - start);
        AllowEventWaitersToProceed(&shared->workDone[threadNum]);
    }
}

GzipWriterFilter::GzipWriterFilter(
    GzipWriterFilterSupplier* i_supplier)
    :
    Filter(DataWriter::TransformFilter),
    supplier(i_supplier),
    bamFormat(i_supplier->bamFormat),
    chunkSize(i_supplier->chunkSize),
    shared(i_supplier->bamFormat, i_supplier->chunkSize, i_supplier->numThreads, i_supplier->bindToProcessors),
    heap(i_supplier->chunkSize * 5)
{
    zstream.zalloc = zalloc;
    zstream.zfree = zfree;
    zstream.opaque = &heap;
    if (i_supplier->multiThreaded) {
        context.totalThreads = i_supplier->numThreads;
        context.bindToProcessors = i_supplier->bindToProcessors;
        context.shared = &shared;
        task = new ParallelTask<GzipContext>(&context);
        task->fork();
    }
}

    void
GzipWriterFilter::onAdvance(
    DataWriter* writer,
    size_t batchOffset,
    char* data,
    unsigned bytes,
    unsigned location)
{
    // nothing
}

    static void
validateZipBlock(
    char* start,
    size_t bytes,
    unsigned uncompressedBytes)
{
    _ASSERT(*(_uint32*)start == 0x04088b1f && *(_uint32*)(start+bytes-4) == uncompressedBytes);
}

    size_t
GzipWriterFilter::onNextBatch(
    DataWriter* writer,
    size_t offset,
    size_t bytes)
{
    char* fromBuffer;
    size_t fromSize, fromUsed, physicalOffset, logicalOffset;
    writer->getBatch(-1, &fromBuffer, &fromSize, &fromUsed, &physicalOffset, NULL, &logicalOffset);
    if (fromUsed == 0) {
        return 0;
    }
    char* toBuffer;
    size_t toSize, toUsed;
    writer->getBatch(0, &toBuffer, &toSize, &toUsed);
    if (toSize - toUsed < fromSize) {
        fprintf(stderr, "GzipWriterFilter: not enough space for compression buffers\n");
        soft_exit(1);
    }
    if (supplier->closing) {
        _ASSERT(toSize - toUsed >= fromUsed);
        memcpy(toBuffer + toUsed, fromBuffer, fromUsed);
        return toUsed + fromUsed;
    }
    int nChunks = (int) (fromUsed / chunkSize);
    size_t extra = fromUsed - chunkSize * nChunks;
    if (supplier->multiThreaded) {
        size_t* sizes = shared.run(nChunks, fromBuffer, toBuffer + toUsed);
        size_t extraUsed = 0;
        if (extra > 0) {
            extraUsed = compressChunk(zstream, supplier->bamFormat,
				      toBuffer + toUsed + nChunks * chunkSize, max(extra, (size_t)1024),
                fromBuffer + nChunks * chunkSize, extra);
        }
        for (int i = 0; i < nChunks; i++) {
            supplier->addTranslation(logicalOffset, physicalOffset + toUsed);
            logicalOffset += chunkSize;
            validateZipBlock(toBuffer + i * chunkSize, sizes[i], (unsigned) chunkSize);
            memmove(toBuffer + toUsed, toBuffer + i * chunkSize, sizes[i]);
            toUsed += sizes[i];
        }
        if (extra > 0) {
            supplier->addTranslation(logicalOffset, physicalOffset + toUsed);
            validateZipBlock(toBuffer + nChunks * chunkSize, extraUsed, (unsigned) extra);
            memmove(toBuffer + toUsed, toBuffer + nChunks * chunkSize, extraUsed);
            toUsed += extraUsed;
        }
    } else {
        for (size_t chunk = 0; chunk < fromUsed; chunk += chunkSize) {
            supplier->addTranslation(logicalOffset + chunk, physicalOffset + toUsed);
            toUsed += compressChunk(zstream, supplier->bamFormat,
                toBuffer + toUsed, toSize - toUsed, 
                fromBuffer + chunk, min(fromUsed - chunk, chunkSize));
        }
    }
    return toUsed;
}

    size_t
GzipWriterFilter::compressChunk(
    z_stream& zstream,
    bool bamFormat,
    char* toBuffer,
    size_t toSize,
    char* fromBuffer,
    size_t fromUsed)
{
    if (bamFormat && fromUsed > 0x10000) {
        fprintf(stderr, "exceeded BAM chunk size\n");
        soft_exit(1);
    }
    if (zstream.opaque != NULL) {
        ((ThreadHeap*)zstream.opaque)->reset();
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

    if (fromUsed > 0xffffffff || toSize > 0xffffffff) {
        fprintf(stderr,"GZipDataWriter: fromUsed or toSize too big\n");
        soft_exit(1);
    }

    // based on sample code at http://www.lemoda.net/c/zlib-open-write/index.html
    const int windowBits = 15;
    const int GZIP_ENCODING = 16;
    zstream.next_in = (Bytef*) fromBuffer;
    zstream.avail_in = (uInt)fromUsed;
    zstream.next_out = (Bytef*) toBuffer;
    zstream.avail_out = (uInt)toSize;
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
    

    GzipWriterFilterSupplier*
DataWriterSupplier::gzip(
    bool bamFormat,
    size_t chunkSize,
    int numThreads,
    bool bindToProcessors,
    bool multiThreaded)
{
    return new GzipWriterFilterSupplier(bamFormat, chunkSize, numThreads, bindToProcessors, multiThreaded);
}
    
    DataWriter::Filter*
GzipWriterFilterSupplier::getFilter()
{
    return new GzipWriterFilter(this);
}

    void
GzipWriterFilterSupplier::onClose(
    DataWriterSupplier* supplier)
{
    if (bamFormat) {
        closing = true;
        DataWriter* writer = supplier->getWriter();
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
        writer->close();
        delete writer;
    }
    
    // sort translations
    std::sort(translation.begin(), translation.end(), translationComparator);
}

    bool
GzipWriterFilterSupplier::translate(
    _uint64 logical,
    _uint64* o_physical,
    _uint64* o_logicalDelta)
{
    pair<_uint64,_uint64> value;
    value.first = logical;
    value.second = 0; //ignored
    pair<_uint64,_uint64>* upper = std::upper_bound(translation.begin(), translation.end(), value, translationComparator);
    if (upper == translation.begin()) {
        return false;
    }
    upper--;
    *o_physical = upper->second;
    *o_logicalDelta = logical - upper->first;
    return true;
}

    bool
GzipWriterFilterSupplier::translationComparator(
    const pair<_uint64,_uint64>& a,
    const pair<_uint64,_uint64>& b)
{
    return a.first < b.first;
}
    
