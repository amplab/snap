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
#include "Util.h"
#include "DataWriter.h"
#include "BufferedAsync.h"
#include "VariableSizeVector.h"
#include "FileFormat.h"
#include "PriorityQueue.h"
#include "exit.h"
#include "Bam.h"

#define USE_DEVTEAM_OPTIONS 1
//#define VALIDATE_SORT 1

using std::max;

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
#ifdef VALIDATE_SORT
    SortBlock() : start(0), bytes(0), location(0), length(0), reader(NULL), minLocation(0), maxLocation(0) {}
#else
    SortBlock() : start(0), bytes(0), location(0), length(0), reader(NULL) {}
#endif
	SortBlock(const SortBlock& other) { *this = other; }
    void operator=(const SortBlock& other);

    size_t      start;
    size_t      bytes;
#ifdef VALIDATE_SORT
	unsigned	minLocation, maxLocation;
#endif
    // for mergesort phase
    DataReader* reader;
    unsigned    location; // genome location of current read
    char*       data; // read data in read buffer
    unsigned    length; // length in bytes
};

    void
SortBlock::operator=(
    const SortBlock& other)
{
    start = other.start;
    bytes = other.bytes;
    location = other.location;
    length = other.length;
    reader = other.reader;
#ifdef VALIDATE_SORT
	minLocation = other.minLocation;
	maxLocation = other.maxLocation;
#endif
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
        const FileFormat* i_fileFormat,
        const Genome* i_genome,
        const char* i_tempFileName,
        const char* i_sortedFileName,
        DataWriter::FilterSupplier* i_sortedFilterSupplier,
        FileEncoder* i_encoder = NULL)
        :
        format(i_fileFormat),
        genome(i_genome),
        FilterSupplier(DataWriter::CopyFilter),
        encoder(i_encoder),
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

    virtual void onClosing(DataWriterSupplier* supplier) {}
    virtual void onClosed(DataWriterSupplier* supplier);

    void setHeaderSize(size_t bytes)
    { headerSize = bytes; }

#ifndef VALIDATE_SORT
	void addBlock(size_t start, size_t bytes);
#else
    void addBlock(size_t start, size_t bytes, unsigned minLocation, unsigned maxLocation);
#endif

private:
    bool mergeSort();

    const Genome*                   genome;
    const FileFormat*               format;
    const char*                     tempFileName;
    const char*                     sortedFileName;
    DataWriter::FilterSupplier*     sortedFilterSupplier;
    FileEncoder*                    encoder;
    size_t                          headerSize;
    ExclusiveLock                   lock; // for adding blocks
    SortBlockVector                 blocks;

	friend class SortedDataFilter;
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
#ifdef VALIDATE_SORT
		if (memcmp(data, "BAM", 3) != 0 && memcmp(data, "@HD", 3) != 0) { // skip header block
			unsigned loc, len;
			parent->format->getSortInfo(parent->genome, data, bytes, &loc, &len);
			_ASSERT(loc == location);
		}
#endif
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
	unsigned previous = 0;
    for (VariableSizeVector<SortEntry>::iterator i = locations.begin(); i != locations.end(); i++) {
#ifdef VALIDATE_SORT
		if (locations.size() > 1) { // skip header block
			unsigned loc, len;
			parent->format->getSortInfo(parent->genome, fromBuffer + i->offset, i->length, &loc, &len);
			_ASSERT(loc >= previous);
			previous = loc;
		}
#endif
        memcpy(toBuffer + target, fromBuffer + i->offset, i->length);
        target += i->length;
    }
    
    // remember block extent for later merge sort
    SortBlock block;
    // handle header specially
    size_t header = offset > 0 ? 0 : locations[0].length;
    if (header > 0) {
        parent->setHeaderSize(header);
    }
	int first = offset == 0;
#ifdef VALIDATE_SORT
	unsigned minLocation = locations.size() > first ? locations[first].location : 0;
	unsigned maxLocation = locations.size() > first ? locations[locations.size()-1].location : UINT32_MAX;
    parent->addBlock(offset + header, bytes - header, minLocation, maxLocation);
#else
    parent->addBlock(offset + header, bytes - header);
#endif
    locations.clear();

    return target;
}
    
    DataWriter::Filter*
SortedDataFilterSupplier::getFilter()
{
    return new SortedDataFilter(this);
}

    void
SortedDataFilterSupplier::onClosed(
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
    size_t bytes
#ifdef VALIDATE_SORT
	, unsigned minLocation
	, unsigned maxLocation
#endif
	)
{
    if (bytes > 0) {
        AcquireExclusiveLock(&lock);
#if VALIDATE_SORT
		for (SortBlockVector::iterator i = blocks.begin(); i != blocks.end(); i++) {
			_ASSERT(i->start + i->length <= start || start + bytes <= i->start);
		}
#endif
        SortBlock block;
        block.start = start;
        block.bytes = bytes;
#if VALIDATE_SORT
		block.minLocation = minLocation;
		block.maxLocation = maxLocation;
#endif
        blocks.push_back(block);
        ReleaseExclusiveLock(&lock);
    }
}

    bool
SortedDataFilterSupplier::mergeSort()
{
    // merge sort from temp file into sorted file
#if USE_DEVTEAM_OPTIONS
    printf("sorting...");
    _int64 start = timeInMillis();
    _int64 startReadWaitTime = DataReader::ReadWaitTime;
    _int64 startReleaseWaitTime = DataReader::ReleaseWaitTime;
    _int64 startWriteWaitTime = DataWriter::WaitTime;
    _int64 startWriteFilterTime = DataWriter::FilterTime;
#endif

    // set up buffered output
    DataWriterSupplier* writerSupplier = DataWriterSupplier::create(sortedFileName, sortedFilterSupplier,
        encoder, encoder != NULL ? 6 : 4); // use more buffers to let encoder run async
    DataWriter* writer = writerSupplier->getWriter();
    if (writer == NULL) {
        fprintf(stderr, "open sorted file for write failed\n");
        return false;
    }
    DataSupplier* readerSupplier = DataSupplier::Default[true]; // autorelease
    // setup - open all files, read first block, begin read for second
    for (SortBlockVector::iterator i = blocks.begin(); i != blocks.end(); i++) {
        i->reader = readerSupplier->getDataReader(MAX_READ_LENGTH * 8); // todo: standardize max length
        i->reader->init(tempFileName);
        i->reader->reinit(i->start, i->bytes);
    }

    // write out header
    if (headerSize > 0xffffffff) {
        fprintf(stderr,"SortedDataFilterSupplier: headerSize too big\n");
        soft_exit(1);
    }
    if (headerSize > 0) {
        blocks[0].reader->reinit(0, headerSize);
		writer->inHeader(true);
        char* rbuffer;
        _int64 rbytes;
        char* wbuffer;
        size_t wbytes;
		for (size_t left = headerSize; left > 0; ) {
			if ((! blocks[0].reader->getData(&rbuffer, &rbytes)) || rbytes == 0) {
				blocks[0].reader->nextBatch();
				if (! blocks[0].reader->getData(&rbuffer, &rbytes)) {
					fprintf(stderr, "read header failed\n");
					soft_exit(1);
				}
			}
			if ((! writer->getBuffer(&wbuffer, &wbytes)) || wbytes == 0) {
				writer->nextBatch();
				if (! writer->getBuffer(&wbuffer, &wbytes)) {
					fprintf(stderr, "write header failed\n");
					soft_exit(1);
				}
			}
			size_t xfer = min(left, min((size_t) rbytes, wbytes));
			_ASSERT(xfer > 0 && xfer <= UINT32_MAX);
			memcpy(wbuffer, rbuffer, xfer);
			blocks[0].reader->advance(xfer);
			writer->advance((unsigned) xfer);
			left -= xfer;
		}
        blocks[0].reader->reinit(blocks[0].start, blocks[0].bytes);
		writer->nextBatch();
		writer->inHeader(false);
    }

    // merge temp blocks into output
    _int64 total = 0;
    // get initial merge sort data
    // queue using complement of location since priority queue is largest first
    typedef PriorityQueue<unsigned,int,-3> BlockQueue;
    BlockQueue queue;
    for (SortBlockVector::iterator b = blocks.begin(); b != blocks.end(); b++) {
        _int64 bytes;
        b->reader->getData(&b->data, &bytes);
        format->getSortInfo(genome, b->data, bytes, &b->location, &b->length);
        queue.put((_uint32) (b - blocks.begin()), ~b->location); 
    }
    unsigned current = 0; // current location for validation
	int lastRefID = -1, lastPos = 0;
    while (queue.size() > 0) {
#if VALIDATE_SORT
		unsigned check;
		queue.peek(&check);
		_ASSERT(~check >= current);
#endif
        unsigned secondLocation;
        int smallestIndex = queue.pop();
        int secondIndex = queue.size() > 0 ? queue.peek(&secondLocation) : -1;
        unsigned limit = secondIndex != -1 ? ~secondLocation : UINT32_MAX;
        SortBlock* b = &blocks[smallestIndex];
        char* writeBuffer;
        size_t writeBytes;
        writer->getBuffer(&writeBuffer, &writeBytes);
        while (b->location <= limit) {
#if VALIDATE_SORT
			_ASSERT(b->location >= b->minLocation && b->location <= b->maxLocation);
#endif
            if (writeBytes < b->length) {
                writer->nextBatch();
                writer->getBuffer(&writeBuffer, &writeBytes);
                if (writeBytes < b->length) {
                    fprintf(stderr, "mergeSort: buffer size too small\n");
                    return false;
                }
            }
            memcpy(writeBuffer, b->data, b->length);
#ifdef VALIDATE_BAM
            if (format == FileFormat::BAM[0] || format == FileFormat::BAM[1]) {
                ((BAMAlignment*)b->data)->validate();
            }
#endif
#if VALIDATE_SORT
			int refID, pos;
			format->getSortInfo(genome, b->data, b->length, NULL, NULL, &refID, &pos);
			_ASSERT(refID == -1 || refID > lastRefID || (refID == lastRefID && pos >= lastPos));
			if (refID != -1) {
				lastRefID = refID;
				lastPos = pos;
			}
#endif
			total++;
            writer->advance(b->length);
            writeBytes -= b->length;
            writeBuffer += b->length;
            b->reader->advance(b->length);
            _ASSERT(b->location >= current);
            current = b->location;
            _int64 readBytes;
            if (! b->reader->getData(&b->data, &readBytes)) {
                b->reader->nextBatch();
                if (! b->reader->getData(&b->data, &readBytes)) {
                    _ASSERT(b->reader->isEOF());
                    delete b->reader;
                    b->reader = NULL;
                    break;
                }
            }
            unsigned previous = b->location;
            format->getSortInfo(genome, b->data, readBytes, &b->location, &b->length);
            _ASSERT(b->length <= readBytes && b->location >= previous);
        }
        if (b->reader != NULL) {
            queue.put(smallestIndex, ~b->location);
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
    printf("sorted %lld reads in %u blocks, %lld s\n"
        "read wait align %.3f s + merge %.3f s, read release align %.3f s + merge %.3f s\n"
        "write wait %.3f s align + %.3f s merge, write filter %.3f s align + %.3f s merge\n",
        total, blocks.size(), (timeInMillis() - start)/1000,
        startReadWaitTime * 1e-9, (DataReader::ReadWaitTime - startReadWaitTime) * 1e-9,
        startReleaseWaitTime * 1e-9, (DataReader::ReleaseWaitTime - startReleaseWaitTime) * 1e-9,
        startWriteWaitTime * 1e-9, (DataWriter::WaitTime - startWriteWaitTime) * 1e-9,
        startWriteFilterTime * 1e-9, (DataWriter::FilterTime - startWriteFilterTime) * 1e-9);
#endif
    return true;
}

    DataWriterSupplier*
DataWriterSupplier::sorted(
    const FileFormat* format,
    const Genome* genome,
    const char* tempFileName,
    size_t tempBufferMemory,
    int numThreads,
    const char* sortedFileName,
    DataWriter::FilterSupplier* sortedFilterSuppler,
    FileEncoder* encoder)
{
    const int bufferCount = 3;
    const size_t bufferSize = tempBufferMemory > 0
        ? tempBufferMemory / (bufferCount * numThreads)
        : max((size_t) 16 * 1024 * 1024, ((size_t) (genome ? genome->getCountOfBases() : 0) / 3) / bufferCount);
    DataWriter::FilterSupplier* filterSupplier =
        new SortedDataFilterSupplier(format, genome, tempFileName, sortedFileName, sortedFilterSuppler, encoder);
    return DataWriterSupplier::create(tempFileName, filterSupplier, NULL, bufferCount, bufferSize);
}
