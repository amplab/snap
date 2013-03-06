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

struct SortBlock
{
    SortBlock() : entries(), fileOffset(0), fileBytes(0), index(0), reader() {}
    SortBlock(int capacity) : entries(capacity), fileOffset(0), fileBytes(0), index(0), reader() {}
    SortBlock(SortBlock& other) { *this = other; }
    void operator=(SortBlock& other);

    VariableSizeVector<SortEntry>   entries;
    size_t                          fileOffset;
    size_t                          fileBytes;
    unsigned                        index;
    BufferedAsyncReader             reader; // for reading phase
};

    void
SortBlock::operator=(
    SortBlock& other)
{
    entries = other.entries;
    fileOffset = other.fileOffset;
    fileBytes = other.fileBytes;
    index = other.index;
    reader = other.reader;
}

class SortedDataFilterSupplier;

class SortedDataFilter : public DataWriter::Filter
{
public:
    SortedDataFilter(SortedDataFilterSupplier* i_parent)
        : Filter(DataWriter::CopyFilter), parent(i_parent), largest(0)
    {}

    virtual ~SortedDataFilter() {}

    virtual void onAdvance(DataWriter* writer, size_t batchOffset, char* data, size_t bytes, unsigned location);

    virtual size_t onNextBatch(DataWriter* writer, size_t offset, size_t bytes);

private:
    SortedDataFilterSupplier*  parent;
    int                         largest; // largest location count so far
    SortBlock                   locations;
};

class SortedDataFilterSupplier : public DataWriter::FilterSupplier
{
public:

    SortedDataFilterSupplier(const char* i_tempFileName, const char* i_sortedFileName,
        DataWriter::FilterSupplier* i_sortedFilterSupplier, int i_readBufferCount, size_t i_readBufferSize)
        :
        FilterSupplier(DataWriter::CopyFilter),
        tempFileName(i_tempFileName),
        sortedFileName(i_sortedFileName),
        sortedFilterSupplier(i_sortedFilterSupplier),
        readBufferCount(i_readBufferCount),
        readBufferSize(i_readBufferSize),
        locations(i_readBufferSize/300)
    {
        InitializeExclusiveLock(&lock);
    }

    virtual ~SortedDataFilterSupplier()
    {
        DestroyExclusiveLock(&lock);
    }

    virtual DataWriter::Filter* getFilter();

    virtual void onClose(DataWriterSupplier* supplier, DataWriter* writer);

    void addLocations(SortBlock& block);

private:
    bool mergeSort();

    const char*                     tempFileName;
    const char*                     sortedFileName;
    const int                       readBufferCount;
    const size_t                    readBufferSize;
    DataWriter::FilterSupplier*     sortedFilterSupplier;
    size_t                          headerSize;
    ExclusiveLock                   lock;
    VariableSizeVector<SortBlock>   locations;
};

    void
SortedDataFilter::onAdvance(
    DataWriter* writer,
    size_t batchOffset,
    char* data,
    size_t bytes,
    unsigned location)
{
    SortEntry entry(batchOffset, bytes, location);
    locations.entries.push_back(entry);
}

    size_t
SortedDataFilter::onNextBatch(
    DataWriter* writer,
    size_t offset,
    size_t bytes)
{
    // sort buffered reads by location for later merge sort
    std::stable_sort(locations.entries.begin(), locations.entries.end(), SortEntry::comparator);
    
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
    for (VariableSizeVector<SortEntry>::iterator i = locations.entries.begin(); i != locations.entries.end(); i++) {
        memcpy(toBuffer + target, fromBuffer + i->offset, i->length);
        i->offset = offset + target; // set correct file offset
        target += i->length;
    }
    
    // remember offsets for full-file sort, get a new vector of appropriate size
    locations.fileOffset = offset;
    locations.fileBytes = bytes;
    parent->addLocations(locations);
    if (locations.entries.size() > largest) {
        largest = (locations.entries.size() * 6) / 5; // grow by 20% if it exceeds prior max
    }
    locations.entries.clear();
    locations.entries.reserve(largest);

    return target;
}
    
    DataWriter::Filter*
SortedDataFilterSupplier::getFilter()
{
    return new SortedDataFilter(this);
}

    void
SortedDataFilterSupplier::onClose(
    DataWriterSupplier* supplier,
    DataWriter* writer)
{
    if (locations.size() == 1 && sortedFilterSupplier == NULL) {
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
SortedDataFilterSupplier::addLocations(
    SortBlock& added)
{
    AcquireExclusiveLock(&lock);
    if (added.fileOffset == 0) {
        _ASSERT(added.entries.size() > 0 && added.entries[0].offset == 0);
        headerSize = added.entries[0].length;
        added.entries.erase(0);
    }
    if (added.entries.size() > 0) {
        locations.push_back(added);
    }
    ReleaseExclusiveLock(&lock);
}

    bool
SortedDataFilterSupplier::mergeSort()
{
    // merge sort from temp file into sorted file
#if USE_DEVTEAM_OPTIONS
    printf("sorting...");
    _int64 start = timeInMillis();
#endif

    // first replace offset in entries with block index, and sort them all in one large array
    size_t total = 0;
    for (VariableSizeVector<SortBlock>::iterator i = locations.begin(); i != locations.end(); i++) {
        total += i->entries.size();
    }
    SortEntry* entries = (SortEntry*) BigAlloc(total * sizeof(SortEntry));
    size_t offset = 0;
    for (VariableSizeVector<SortBlock>::iterator i = locations.begin(); i != locations.end(); i++) {
        unsigned n = i->entries.size();
        size_t blockIndex = i - locations.begin();
        memcpy(entries + offset, i->entries.begin(), n * (size_t) sizeof(SortEntry));
        for (unsigned j = 0; j < n; j++) {
            entries[offset + j].offset = blockIndex;
        }
        offset += n;
    }
    std::stable_sort(entries, entries + total, SortEntry::comparator);
#if USE_DEVTEAM_OPTIONS
    printf(" %ld s\nwriting sorted reads...", (timeInMillis() - start) / 1000);
    start = timeInMillis();
#endif

    // setup - open all files, read first block, begin read for second
    AsyncFile* temp = AsyncFile::open(tempFileName, false);
    char* buffers = (char*) BigAlloc(readBufferSize * 2 * locations.size());
    unsigned j = 0;
#if USE_DEVTEAM_OPTIONS
    if (timeInMillis() - start > 1000) {
        printf(" (allocated %lld Mb in %lld s)",
            readBufferSize * 2 * locations.size() / (2 << 20), (timeInMillis() - start) / 1000);
    }
#endif
    for (VariableSizeVector<SortBlock>::iterator i = locations.begin(); i != locations.end(); i++, j++) {
        i->reader.open(temp, i->fileOffset, i->fileBytes, readBufferSize, true,
            buffers + j * 2 * readBufferSize, buffers + (j * 2 + 1) * readBufferSize);
    }
    for (VariableSizeVector<SortBlock>::iterator i = locations.begin(); i != locations.end(); i++) {
        i->reader.endOpen();
    }

    // set up buffered output
    DataWriterSupplier* writerSupplier = DataWriterSupplier::create(sortedFileName, sortedFilterSupplier, readBufferCount, readBufferSize);
    DataWriter* writer = writerSupplier->getWriter();
    if (writer == NULL) {
        fprintf(stderr, "open sorted file for write failed\n");
        return false;
    }

    // write out header
    char* buffer;
    size_t bytes;
    if (headerSize > 0) {
        AsyncFile::Reader* hread = NULL;
        bool ok = writer->getBuffer(&buffer, &bytes) &&
            bytes >= headerSize &&
            (hread = temp->getReader()) != NULL &&
            hread->beginRead(buffer, headerSize, 0, NULL) &&
            hread->waitForCompletion() &&
            hread->close();
        if (hread != NULL) {
            delete hread;
        }
        if (! ok ) {
            fprintf(stderr, "read header failed\n");
            return false;
        }
        writer->advance(headerSize);
    }

    // merge input blocks into output using pre-sorted list
    for (size_t i = 0; i < total; i++) {
        SortEntry* entry = &entries[i];
        writer->getBuffer(&buffer, &bytes);
        if (bytes < entry->length) {
            writer->nextBatch();
            writer->getBuffer(&buffer, &bytes);
            if (bytes < entry->length) {
                fprintf(stderr, "mergeSort: buffer size too small\n");
                return false;
            }
        }
        if (! locations[entry->offset].reader.read(buffer, entry->length)) {
            fprintf(stderr, "merge sort read failed\n");
            return false;
        }
        writer->advance(entry->length, entry->location);
    }

    // close everything
    writer->close();
    delete writer;
    writerSupplier->close();
    delete writerSupplier;
    BigDealloc(entries);
    bool ok = true;
    _int64 readWait = 0;
    for (VariableSizeVector<SortBlock>::iterator i = locations.begin(); i != locations.end(); i++) {
        ok &= i->reader.close();
        readWait += i->reader.getWaitTimeInMillis();
    }
    ok &= temp->close();
    delete temp;
    BigDealloc(buffers);
    if (! ok) {
        printf("files did not close properly\n");
    }
    if (! DeleteSingleFile(tempFileName)) {
        printf("warning: failure deleting temp file %s\n", tempFileName);
    }

#if USE_DEVTEAM_OPTIONS
    printf(" %u reads in %u blocks, %lld s (%lld s read wait)\n",
        total, locations.size(), (timeInMillis() - start)/1000, readWait/1000);
#endif
#ifdef PROFILE_BIGALLOC
    PrintAllocProfile();
#endif
    return true;
}

    DataWriterSupplier*
DataWriterSupplier::sorted(
    const char* tempFileName,
    const char* sortedFileName,
    DataWriter::FilterSupplier* sortedFilterSuppler,
    size_t bufferSize,
    int buffers)
{
    DataWriter::FilterSupplier* filterSupplier =
        new SortedDataFilterSupplier(tempFileName, sortedFileName, sortedFilterSuppler, buffers, bufferSize);
    return create(tempFileName, filterSupplier, buffers, bufferSize);
}
