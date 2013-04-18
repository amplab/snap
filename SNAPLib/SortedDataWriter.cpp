/*++

Module Name:

    SortedDataWriter.cpp

Abstract:

    File writer that sorts records using a temporary file.

Environment:

    User mode service.

    Not thread safe.

--*/

#include "stdafx.h"
#include "BigAlloc.h"
#include "Compat.h"
#include "DataWriter.h"
#include "BufferedAsync.h"
#include "VariableSizeVector.h"
#include "FileFormat.h"
#include "exit.h"

#define USE_DEVTEAM_OPTIONS 1

#pragma pack(push, 4)
struct SortEntry
{
    SortEntry() : offset(0), length(0), location(0) {}
    SortEntry(size_t i_offset, unsigned i_length, unsigned i_location)
        : offset(i_offset), length(i_length), location(i_location) {}
    size_t                      offset; // offset in file
    unsigned                    length; // number of bytes
    unsigned                    location; // location in genome
    static bool comparator(const SortEntry& e1, const SortEntry& e2)
    {
        return e1.location < e2.location;
    }
};
#pragma pack(pop)

typedef VariableSizeVector<SortEntry> SortVector;

struct SortBlock
{
    SortBlock() : start(0), bytes(0), location(0), length(0), reader(NULL) {}
    SortBlock(SortBlock& other) { *this = other; }
    void operator=(SortBlock& other);

    size_t      start;
    size_t      bytes;
    // for mergesort phase
    DataReader* reader;
    unsigned    location; // genome location of current read
    char*       data; // read data in read buffer
    unsigned    length; // length in bytes
};

    void
SortBlock::operator=(
    SortBlock& other)
{
    start = other.start;
    bytes = other.bytes;
    location = other.location;
    length = other.length;
    reader = other.reader;
}

typedef VariableSizeVector<SortBlock> SortBlockVector;
    
class SortedDataFilterSupplier;

class SortedDataFilter : public DataWriter::Filter
{
public:
    SortedDataFilter(SortedDataFilterSupplier* i_parent)
        : Filter(DataWriter::CopyFilter), parent(i_parent), locations()
    {}

    virtual ~SortedDataFilter() {}

    virtual void onAdvance(DataWriter* writer, size_t batchOffset, char* data, unsigned bytes, unsigned location);

    virtual size_t onNextBatch(DataWriter* writer, size_t offset, size_t bytes);

private:
    SortedDataFilterSupplier*   parent;
    SortVector                  locations;
};

class SortedDataFilterSupplier : public DataWriter::FilterSupplier
{
public:

    SortedDataFilterSupplier(
        const Genome* i_genome,
        const char* i_tempFileName,
        const char* i_sortedFileName,
        DataWriter::FilterSupplier* i_sortedFilterSupplier)
        :
        genome(i_genome),
        FilterSupplier(DataWriter::CopyFilter),
        tempFileName(i_tempFileName),
        sortedFileName(i_sortedFileName),
        sortedFilterSupplier(i_sortedFilterSupplier),
        blocks()
    {
        InitializeExclusiveLock(&lock);
    }

    virtual ~SortedDataFilterSupplier()
    {
        DestroyExclusiveLock(&lock);
    }

    virtual DataWriter::Filter* getFilter();

    virtual void onClose(DataWriterSupplier* supplier);

    void setHeaderSize(size_t bytes)
    { headerSize = bytes; }

    void addBlock(size_t start, size_t bytes);

private:
    bool mergeSort();

    const Genome*                   genome;
    const char*                     tempFileName;
    const char*                     sortedFileName;
    DataWriter::FilterSupplier*     sortedFilterSupplier;
    size_t                          headerSize;
    ExclusiveLock                   lock; // for adding blocks
    SortBlockVector                 blocks;
};

    void
SortedDataFilter::onAdvance(
    DataWriter* writer,
    size_t batchOffset,
    char* data,
    unsigned bytes,
    unsigned location)
{
    SortEntry entry(batchOffset, bytes, location);
    locations.push_back(entry);
}

    size_t
SortedDataFilter::onNextBatch(
    DataWriter* writer,
    size_t offset,
    size_t bytes)
{
    // sort buffered reads by location for later merge sort
    std::stable_sort(locations.begin(), locations.end(), SortEntry::comparator);
    
    // copy from previous buffer into current in sorted order
    char* fromBuffer;
    size_t fromSize, fromUsed;
    char* toBuffer;
    size_t toSize, toUsed;
    if (! (writer->getBatch(-1, &fromBuffer, &fromSize, &fromUsed) &&
        writer->getBatch(0, &toBuffer, &toSize, &toUsed)))
    {
        fprintf(stderr, "SortedDataFilter::onNextBatch getBatch failed\n");
        soft_exit(1);
    }
    size_t target = 0;
    for (VariableSizeVector<SortEntry>::iterator i = locations.begin(); i != locations.end(); i++) {
        memcpy(toBuffer + target, fromBuffer + i->offset, i->length);
        i->offset = offset + target; // set correct file offset
        target += i->length;
    }
    
    // remember block extent for later merge sort
    SortBlock block;
    // handle header specially
    size_t header = offset > 0 ? 0 : locations[0].length;
    if (header > 0) {
        parent->setHeaderSize(header);
    }
    parent->addBlock(offset + header, bytes - header);
    locations.clear();

    return target;
}
    
    DataWriter::Filter*
SortedDataFilterSupplier::getFilter()
{
    return new SortedDataFilter(this);
}

    void
SortedDataFilterSupplier::onClose(
    DataWriterSupplier* supplier)
{
    if (blocks.size() == 1 && sortedFilterSupplier == NULL) {
        // just rename/move temp file to real file, we're done
        DeleteSingleFile(sortedFileName); // if it exists
        if (! MoveSingleFile(tempFileName, sortedFileName)) {
            fprintf(stderr, "unable to move temp file %s to final sorted file %s\n", tempFileName, sortedFileName);
            soft_exit(1);
        }
        return;
    }
    // merge sort into final file
    if (! mergeSort()) {
        fprintf(stderr, "merge sort failed\n");
        soft_exit(1);
    }
}

    void
SortedDataFilterSupplier::addBlock(
    size_t start,
    size_t bytes)
{
    if (bytes > 0) {
        AcquireExclusiveLock(&lock);
        SortBlock block;
        block.start = start;
        block.bytes = bytes;
        blocks.push_back(block);
        ReleaseExclusiveLock(&lock);
    }
}

    bool
SortedDataFilterSupplier::mergeSort()
{
    // merge sort from temp file into sorted file
    FileFormat* format;
#if USE_DEVTEAM_OPTIONS
    printf("sorting...");
    _int64 start = timeInMillis();
#endif

    // set up buffered output
    DataWriterSupplier* writerSupplier = DataWriterSupplier::create(sortedFileName, sortedFilterSupplier);
    DataWriter* writer = writerSupplier->getWriter();
    if (writer == NULL) {
        fprintf(stderr, "open sorted file for write failed\n");
        return false;
    }
    const DataSupplier* readerSupplier = DataSupplier::Default[true]; // autorelease
    // setup - open all files, read first block, begin read for second
    for (SortBlockVector::iterator i = blocks.begin(); i != blocks.end(); i++) {
        i->reader = readerSupplier->getDataReader(8192); // todo: parameterize max read len
        i->reader->init(tempFileName);
        i->reader->reinit(i->start, i->length);
    }

    // write out header
    char* buffer;
    size_t bytes;
    if (headerSize > 0) {
        blocks[0].reader->reinit(0, headerSize);
        char* rbuffer;
        _int64 rbytes;
        char* wbuffer;
        size_t wbytes;
        bool ok = blocks[0].reader->getData(&rbuffer, &rbytes) &&
            rbytes >= headerSize &&
            writer->getBuffer(&wbuffer, &wbytes) &&
            wbytes >= headerSize;
        if (! ok ) {
            fprintf(stderr, "read header failed\n");
            return false;
        }
        if (headerSize > 0xffffffff) {
            fprintf(stderr,"SortedDataFilterSupplier: headerSize too big\n");
            soft_exit(1);
        }
        memcpy(wbuffer, rbuffer, headerSize);
        writer->advance((unsigned)headerSize);
        blocks[0].reader->reinit(blocks[0].start, blocks[0].bytes);
    }

    // merge temp blocks into output
    _int64 total = 0;
    // get initial merge sort data
    for (SortBlockVector::iterator b = blocks.begin(); b != blocks.end(); b++) {
        _int64 bytes;
        b->reader->getData(&b->data, &bytes);
        format->getSortInfo(genome, b->data, bytes, &b->location, &b->length);
    }
    while (true) {
        int smallest = -1, second = -1;
        for (int i = 0; i < blocks.size(); i++) {
            SortBlock* b = &blocks[i];
            if (b->reader == NULL) {
                continue;
            }
            if (smallest == -1 || b->location < blocks[smallest].location) {
                second = smallest;
                smallest = i;
            } else if (second == -1 || b->location < blocks[second].location) {
                second = i;
            }
        }
        if (smallest == -1) {
            break;
        }
        unsigned limit = second != -1 ? blocks[second].location : UINT32_MAX;
        SortBlock* b = &blocks[smallest];
        char* writeBuffer;
        size_t writeBytes;
        writer->getBuffer(&writeBuffer, &writeBytes);
        while (b->location <= limit) {
            if (writeBytes < b->length) {
                writer->nextBatch();
                writer->getBuffer(&writeBuffer, &writeBytes);
                if (writeBytes < b->length) {
                    fprintf(stderr, "mergeSort: buffer size too small\n");
                    return false;
                }
            }
            memcpy(writeBuffer, b->data, b->length);
            writer->advance(b->length);
            b->reader->advance(b->length);
            _int64 readBytes;
            if (! b->reader->getData(&b->data, &readBytes)) {
                b->reader->nextBatch();
                if (! b->reader->getData(&b->data, &readBytes)) {
                    _ASSERT(b->reader->isEOF());
                    delete b->reader;
                    b->reader = NULL;
                }
            }
            format->getSortInfo(genome, b->data, readBytes, &b->location, &b->length);
            _ASSERT(b->length <= readBytes);
        }
    }
    
    // close everything
    writer->close();
    delete writer;
    writerSupplier->close();
    delete writerSupplier;
    if (! DeleteSingleFile(tempFileName)) {
        printf("warning: failure deleting temp file %s\n", tempFileName);
    }

#if USE_DEVTEAM_OPTIONS
    printf("sorted %lld reads in %u blocks, %lld s\n",
        total, blocks.size(), (timeInMillis() - start)/1000);
#endif
    return true;
}

    DataWriterSupplier*
DataWriterSupplier::sorted(
    const Genome* genome,
    const char* tempFileName,
    size_t tempBufferMemory,
    int numThreads,
    const char* sortedFileName,
    DataWriter::FilterSupplier* sortedFilterSuppler)
{
    const int bufferCount = 3;
    const size_t bufferSize = tempBufferMemory > 0
        ? tempBufferMemory / (bufferCount * numThreads)
        : max((size_t) 16 * 1024 * 1024, ((size_t) genome->getCountOfBases() / 3) / bufferCount);
    DataWriter::FilterSupplier* filterSupplier =
        new SortedDataFilterSupplier(genome, tempFileName, sortedFileName, sortedFilterSuppler);
    return DataWriterSupplier::create(tempFileName, filterSupplier, bufferCount, bufferSize);
}
