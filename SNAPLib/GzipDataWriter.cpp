/*++

Module Name:

    GzipDataWriter.cpp

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
#include "Bam.h"
#include "zlib.h"
#include "exit.h"
#include "Error.h"

using std::min;
using std::max;
using std::pair;

class GzipWriterFilterSupplier;

class GzipCompressWorkerManager : public ParallelWorkerManager
{
public:
    GzipCompressWorkerManager(GzipWriterFilterSupplier* i_filterSupplier)
        : filterSupplier(i_filterSupplier), buffer(NULL),
        chunkSize(i_filterSupplier->chunkSize), bam(i_filterSupplier->bamFormat)
    {}

    virtual ~GzipCompressWorkerManager();

    virtual void initialize(void* i_writer);

    virtual ParallelWorker* createWorker();

    virtual void beginStep();

    virtual void finishStep();

private:
    VariableSizeVector<size_t> sizes;
    volatile int nChunks;
    const size_t chunkSize;
    const bool bam;
    FileEncoder* encoder;
    GzipWriterFilterSupplier* filterSupplier;
    char* input;
    size_t inputSize;
    size_t inputUsed;
    char* buffer;
    VariableSizeVector< pair<_uint64,_uint64> > translation;

    friend class GzipCompressWorker;
};

class GzipCompressWorker : public ParallelWorker
{
public:
    GzipCompressWorker() : heap(NULL) {}

    virtual ~GzipCompressWorker() { delete heap; }

    virtual void step();

    static size_t compressChunk(z_stream& zstream, bool bamFormat, char* toBuffer, size_t toSize, char* fromBuffer, size_t fromUsed);

private:
    z_stream zstream;
    ThreadHeap* heap;
};

// used for case where each thread compresses by itself

class GzipWriterFilter : public DataWriter::Filter
{
public:
    GzipWriterFilter(GzipWriterFilterSupplier* i_supplier);

    virtual void onAdvance(DataWriter* writer, size_t batchOffset, char* data, GenomeDistance bytes, GenomeLocation location);

    virtual size_t onNextBatch(DataWriter* writer, size_t offset, size_t bytes);

private:

    GzipWriterFilterSupplier* supplier;
    // if doing inline compression, filled in with minimally initialized objects
    GzipCompressWorkerManager* manager;
    ParallelWorker* worker;
    FileEncoder* encoder;
};

GzipCompressWorkerManager::~GzipCompressWorkerManager()
{
    if (buffer != NULL) {
        delete buffer;
    }
}

    void
GzipCompressWorkerManager::initialize(
    void* i_encoder)
{
    encoder = (FileEncoder*) i_encoder;
}

    ParallelWorker*
GzipCompressWorkerManager::createWorker()
{
    return new GzipCompressWorker();
}

    void
GzipCompressWorkerManager::beginStep()
{
    if (filterSupplier->closing) {
        nChunks = 0;
        return;
    }
    encoder->getEncodeBatch(&input, &inputSize, &inputUsed);
    nChunks = (int) ((inputUsed + chunkSize - 1) / chunkSize);
    sizes.clear();
    sizes.extend(nChunks);

    if (buffer == NULL) {
        buffer = (char*) BigAlloc(inputSize);
    }
}

    void
GzipCompressWorkerManager::finishStep()
{
    if (filterSupplier->closing) {
        return;
    }
    size_t toUsed = 0, logicalOffset, physicalOffset;
    encoder->getOffsets(&logicalOffset, &physicalOffset);
    for (int i = 0; i < nChunks; i++) {
        translation.push_back(pair<_uint64,_uint64>(logicalOffset, physicalOffset + toUsed));
        _ASSERT(i * chunkSize < inputUsed);
        _ASSERT(sizes[i] <= chunkSize);
        size_t logicalChunk = min(chunkSize, inputUsed - i * chunkSize);
        logicalOffset += logicalChunk;
        _ASSERT(((BgzfHeader*)(buffer + i * chunkSize))->validate(sizes[i], logicalChunk));
        memcpy(input + toUsed, buffer + i * chunkSize, sizes[i]);
        toUsed += sizes[i];
    }
    _ASSERT(BgzfHeader::validate(input, toUsed));
    encoder->setEncodedBatchSize(toUsed);
    filterSupplier->addTranslations(&translation);
    translation.clear();
}

    void
GzipCompressWorker::step()
{
    GzipCompressWorkerManager* supplier = (GzipCompressWorkerManager*) getManager();
    if (heap == NULL) {
        heap = new ThreadHeap(supplier->chunkSize * 8); // appears to use 4*chunkSize per run
        zstream.zalloc = zalloc;
        zstream.zfree = zfree;
        zstream.opaque = heap;
    }
    //fprintf(stderr, "zip task thread %d begin\n", GetCurrentThreadId());
    _int64 start = timeInMillis();
    int begin = (getThreadNum() * supplier->nChunks) / getNumThreads();
    int end = ((1 + getThreadNum()) * supplier->nChunks) / getNumThreads();
    for (int i = begin; i < end; i++) {
        size_t bytes = min(supplier->chunkSize, supplier->inputUsed - i * supplier->chunkSize);
        supplier->sizes[i] = compressChunk(zstream, supplier->bam,
            supplier->buffer + i * supplier->chunkSize, supplier->chunkSize,
            supplier->input + i * supplier->chunkSize, bytes);
        _ASSERT(supplier->sizes[i] <= supplier->chunkSize); // can't grow!
    }
}


    size_t
GzipCompressWorker::compressChunk(
    z_stream& zstream,
    bool bamFormat,
    char* toBuffer,
    size_t toSize,
    char* fromBuffer,
    size_t fromUsed)
{
    if (bamFormat && fromUsed > BAM_BLOCK) {
        WriteErrorMessage("exceeded BAM chunk size\n");
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
        WriteErrorMessage("GZipDataWriter: fromUsed or toSize too big\n");
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
        WriteErrorMessage("GzipWriterFilter: deflateInit2 failed with %d\n", status);
        soft_exit(1);
    }
    if (bamFormat) {
        status = deflateSetHeader(&zstream, &header);
        if (status != Z_OK) {
            WriteErrorMessage("GzipWriterFilter: defaultSetHeader failed with %d\n", status);
            soft_exit(1);
        }
    }
    oldAvail = zstream.avail_out;
    status = deflate(&zstream, Z_FINISH);
    if (status < 0 && status != Z_BUF_ERROR) {
        WriteErrorMessage("GzipWriterFilter: deflate failed with %d\n", status);
        soft_exit(1);
    }

    // make sure it all got written out in a single compressed block
    if (zstream.avail_in != 0) {
        WriteErrorMessage("GzipWriterFilter: default failed to read all input\n");
        soft_exit(1);
    }
    if (zstream.avail_out == oldAvail) {
        WriteErrorMessage("GzipWriterFilter: default failed to write output\n");
        soft_exit(1);
    }
    status = deflateEnd(&zstream);
    if (status < 0) {
        WriteErrorMessage("GzipWriterFilter: deflateEnd failed with %d\n", status);
        soft_exit(1);
    }

    size_t toUsed = toSize - zstream.avail_out;
    if (bamFormat) {
        // backpatch compressed block size into gzip header
        if (toUsed >= BAM_BLOCK) {
            WriteErrorMessage("exceeded BAM chunk size\n");
            soft_exit(1);
        }
        * (_uint16*) (toBuffer + 16) = (_uint16) (toUsed - 1);
    }
    return toUsed;
}

GzipWriterFilter::GzipWriterFilter(GzipWriterFilterSupplier* i_supplier)
    : DataWriter::Filter(DataWriter::ResizeFilter), supplier(i_supplier), manager(NULL), worker(NULL)
{}


    void
GzipWriterFilter::onAdvance(
    DataWriter* writer,
    size_t batchOffset,
    char* data,
    GenomeDistance bytes,
    GenomeLocation location)
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
    size_t fromSize, fromUsed, physicalOffset, logicalOffset;
    writer->getBatch(-1, &fromBuffer, &fromSize, &fromUsed, &physicalOffset, NULL, &logicalOffset);
    if (fromUsed == 0 || supplier->multiThreaded || supplier->closing) {
        return fromUsed;
    }
    // do compress buffer synchronously in-place
    if (manager == NULL) {
        manager = new GzipCompressWorkerManager(supplier);
        worker = manager->createWorker();
        encoder = new FileEncoder(0, false, manager);
        encoder->initialize((AsyncDataWriter*) writer);
        manager->initialize(encoder);
        manager->configure(worker, 0, 1);
    }
    encoder->setupEncode(-1);
    manager->beginStep();
    worker->step();
    manager->finishStep();
    writer->getBatch(-1, &fromBuffer, &fromSize, &fromUsed, &physicalOffset, NULL, &logicalOffset);
    return fromUsed;
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
GzipWriterFilterSupplier::onClosing(
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
            WriteErrorMessage("no space to write eof marker\n");
            soft_exit(1);
        }
        memcpy(buffer, eof, sizeof(eof));
        writer->advance(sizeof(eof));

        // add final translation for last empty block
        writer->nextBatch();
        char* ignore;
        pair<_uint64,_uint64> last;
        size_t used;
        writer->getBatch(-1, &ignore, NULL, &used, (size_t*) &last.second, NULL, (size_t*) &last.first);
        last.second += used;
        translation.push_back(last);

        writer->close();
        delete writer;
    }

    // sort translations
    std::sort(translation.begin(), translation.end(), translationComparator);
}

    void
GzipWriterFilterSupplier::addTranslations(
    VariableSizeVector< pair<_uint64,_uint64> >* moreTranslations)
{
    AcquireExclusiveLock(&lock);
    translation.append(moreTranslations);
    ReleaseExclusiveLock(&lock);
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

    FileEncoder*
FileEncoder::gzip(
    GzipWriterFilterSupplier* filterSupplier,
    int numThreads,
    bool bindToProcessor,
    size_t chunkSize,
    bool bam)
{
    return new FileEncoder(numThreads, bindToProcessor, new GzipCompressWorkerManager(filterSupplier));
}
