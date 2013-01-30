/*++

Module Name:

    Bam.cpp

Abstract:

    Binary Alignment Map (BAM) file writer and reader.

Environment:

    User mode service.

    BamWriter and BamReader (and their subclasses) aren't thread safe.

--*/

#include "stdafx.h"
#include "BigAlloc.h"
#include "Compat.h"
#include "Read.h"
#include "Bam.h"
#include "Tables.h"
#include "RangeSplitter.h"
#include "ParallelTask.h"
#include "ReadSupplierQueue.h"

using std::max;
using std::min;

BAMReader::BAMReader(
    ReadClippingType i_clipping,
    const Genome* i_genome,
    bool i_paired)
    : clipping(i_clipping), genome(i_genome), paired(i_paired)
{
}

BAMReader::~BAMReader()
{
}

    bool
BAMReader::getNextReadPair(
    Read *read1,
    Read *read2,
    PairedAlignmentResult *alignmentResult,
    unsigned *mapQ,
    const char **cigar)
{
    return false;
}


    void
BAMReader::readDoneWithBuffer(unsigned *referenceCount)
{
}

    bool
BAMReader::init(
    const char *fileName,
    _int64 startingOffset,
    _int64 amountOfFileToProcess)
{
    if (! buffers.init(fileName)) {
        return false;
    }
    size_t headerSize = 1024 * 1024; // 1M header max
    char* buffer = buffers.readHeader(&headerSize);
    BAMHeader* header = (BAMHeader*) buffer;
    if (header->magic != BAMHeader::BAM_MAGIC) {
        fprintf(stderr, "BAMReader: Not a valid BAM file\n");
        return false;
    }
    size_t textHeaderSize = header->l_text;
    if (!SAMReader::parseHeader(fileName, header->text(), header->text() + textHeaderSize, genome, &textHeaderSize)) {
        fprintf(stderr,"BAMReader: failed to parse header on '%s'\n",fileName);
        return false;
    }
    n_ref = header->n_ref();
    refOffset = new unsigned[n_ref];
    BAMHeaderRefSeq* refSeq = header->firstRefSeq();
    for (int i = 0; i < n_ref; i++, refSeq = refSeq->next()) {
        if (! genome->getOffsetOfPiece(refSeq->name(), &refOffset[i])) {
            // fprintf(stderr, "BAMReader: unknown ref seq name %s\n", refSeq->name());
            refOffset[i] = UINT32_MAX;
            // exit(1); ??
        }
    }
    _int64 skip = (char*) refSeq - buffer;
    buffers.reinit(startingOffset + skip, max(0LL, amountOfFileToProcess - skip));
    return true;
}

    BAMReader*
BAMReader::create(
    const char *fileName,
    const Genome *genome,
    _int64 startingOffset,
    _int64 amountOfFileToProcess,
    ReadClippingType clipping,
    bool paired)
{
    BAMReader* reader = new BAMReader(clipping, genome, paired);
    if (! reader->init(fileName, startingOffset, amountOfFileToProcess)) {
        delete reader;
        return NULL;
    }
    return reader;
}
    
    void
BAMReader::reinit(
    _int64 startingOffset,
    _int64 amountOfFileToProcess)
{
    buffers.reinit(startingOffset, amountOfFileToProcess);
}

    ReadSupplierGenerator *
BAMReader::createReadSupplierGenerator(
    const char *fileName,
    int numThreads,
    const Genome *genome,
    ReadClippingType clipping)
{
    BAMReader* reader = create(fileName, genome, 0, 0, clipping);
    ReadSupplierQueue* queue = new ReadSupplierQueue((ReadReader*)reader);
    queue->startReaders();
    return queue;
}

    PairedReadSupplierGenerator *
BAMReader::createPairedReadSupplierGenerator(
    const char *fileName,
    int numThreads,
    const Genome *genome,
    ReadClippingType clipping)
{
    BAMReader* reader = create(fileName, genome, 0, 0, clipping);
    ReadSupplierQueue* queue = new ReadSupplierQueue((PairedReadReader*)reader);
    return queue;
}
    
    void
BAMReader::expandSeq(
    char* o_sequence,
    _uint8* nibbles,
    int bases)
{
    for (int i = 0; i < bases; i++) {
        int n = (i & 1) ? *nibbles & 0xf : *nibbles >> 4;
        nibbles += i & 1;
        *o_sequence++ = "=ACMGRSVTWYHKDBN"[n];
    }
}
    
    void
BAMReader::expandQual(
    char* o_qual,
    char* quality,
    int bases)
{
    for (int i = 0; i < bases; i++) {
        char q = quality[i];
        o_qual[i] = q < 0 || q >= 64 ? '!' : q + '!';
    }
}

    void
BAMReader::expandCigar(
    char* o_sequence,
    _uint32* cigar,
    int ops)
{
    int i = 0;
    while (ops > 0 && i < MAX_SEQ_LENGTH - 20) {
        i += sprintf(o_sequence + i, "%u", *cigar >> 4);
        _ASSERT(*cigar & 0xf <= 8);
        o_sequence[i++] = "MIDNSHP=X"[*cigar & 0xf];
        ops--;
        cigar++;
    }
    o_sequence[i++] = 0;
}

    bool
BAMReader::getNextRead(
    Read *read,
    AlignmentResult *alignmentResult, 
    unsigned *genomeLocation,
    bool *isRC,
    unsigned *mapQ,
    unsigned *flag,
    bool ignoreEndOfRange,
    const char **cigar)
{
    char* buffer;
    DWORD bytes;
    if (! buffers.getData(&buffer, &bytes)) {
        return false;
    }
    BAMAlignment* bam = (BAMAlignment*) buffer;
    if (bytes < sizeof(bam->block_size) || bytes < bam->size()) {
        if (buffers.isEOF()) {
            fprintf(stderr, "Unexpected end of BAM file at %lld\n", buffers.getFileOffset());
            exit(1);
        }
        char* overflow = buffers.useOverflowBuffer();
        DWORD copied = bytes;
        memcpy(overflow, buffer, bytes);
        buffers.nextBuffer();
        buffers.getData(&buffer, &bytes, ignoreEndOfRange);
        _ASSERT(buffer != NULL && bytes > 0);
        if (copied < sizeof(bam->block_size)) {
            memcpy(overflow + copied, buffer, sizeof(bam->block_size) - copied);
            _ASSERT(copied < bam->size());
        }
        _ASSERT(bam->size() < MAX_RECORD_LENGTH && bam->size() <= copied + bytes);
        memcpy(overflow + copied, buffer, bam->size() - copied);
        buffers.addOffset(bam->size() - copied);
        buffer = overflow;
        bytes += bam->size();
    } else {
        buffers.addOffset(bam->size());
    }
    size_t lineLength;
    getReadFromLine(genome, buffer, buffer + bytes, read, alignmentResult, genomeLocation,
        isRC, mapQ, &lineLength, flag, cigar, clipping);
    return true;
}

    void
BAMReader::getReadFromLine(
    const Genome *genome,
    char *line,
    char *endOfBuffer,
    Read *read,
    AlignmentResult *alignmentResult,
    unsigned *genomeLocation,
    bool *isRC,
    unsigned *mapQ, 
    size_t *lineLength,
    unsigned *flag,
    const char **cigar,
    ReadClippingType clipping)
{
    _ASSERT(endOfBuffer - line >= sizeof(BAMHeader));
    BAMAlignment* bam = (BAMAlignment*) line;
    _ASSERT(endOfBuffer - line >= bam->size());
    
    if (NULL != genomeLocation) {
        _ASSERT(-1 <= bam->refID && bam->refID < n_ref);
        *genomeLocation = bam->refID != -1 && refOffset[bam->refID] != UINT32_MAX
            ? refOffset[bam->refID] + bam->pos : 0xffffffff;
    }

    int mateIndex = (bam->FLAG & SAM_MULTI_SEGMENT) && (bam->FLAG && SAM_LAST_SEGMENT) ? 1 : 0;
    if (NULL != read) {
        _ASSERT(bam->l_seq < MAX_SEQ_LENGTH);
        expandSeq(seqBuffer[mateIndex], bam->seq(), bam->l_seq);
        expandQual(qualBuffer[mateIndex], bam->qual(), bam->l_seq);
        read->init(bam->read_name(), bam->l_read_name, seqBuffer[mateIndex], qualBuffer[mateIndex], bam->l_seq);
        if (bam->FLAG & SAM_REVERSE_COMPLEMENT) {
            read->becomeRC();
        }
        read->clip(clipping);
    }

    if (NULL != alignmentResult) {
        _ASSERT(bam->FLAG & SAM_UNMAPPED || bam->refID >= 0);
        *alignmentResult = bam->FLAG & SAM_UNMAPPED ? NotFound : SingleHit; // todo: look at MAPQ?
    }

    if (NULL != isRC) {
        *isRC = (bam->FLAG & SAM_REVERSE_COMPLEMENT) != 0;
    }

    if (NULL != mapQ) {
        *mapQ = bam->MAPQ;
    }

    if (NULL != flag) {
        *flag = bam->FLAG;
    }

    if (NULL != cigar) {
        expandCigar(cigarBuffer[mateIndex], bam->cigar(), bam->n_cigar_op);
        *cigar = cigarBuffer[mateIndex];
    }
}
