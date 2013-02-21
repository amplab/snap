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
#include "Util.h"
#include "FileFormat.h"
#include "AlignerOptions.h"

using std::max;
using std::min;
using util::strnchr;

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

    bool
BAMReader::init(
    const char *fileName,
    _int64 startingOffset,
    _int64 amountOfFileToProcess)
{
    // todo: integrate supplier models
    // might need up to 2x extra for expanded sequence + quality + cigar data
    data = DataSupplier::GzipDefault->getDataReader(MAX_RECORD_LENGTH, 2.5);
    if (! data->init(fileName)) {
        return false;
    }
    _int64 headerSize = 1024 * 1024; // 1M header max
    char* buffer = data->readHeader(&headerSize);
    BAMHeader* header = (BAMHeader*) buffer;
    if (header->magic != BAMHeader::BAM_MAGIC) {
        fprintf(stderr, "BAMReader: Not a valid BAM file\n");
        return false;
    }
    _int64 textHeaderSize = header->l_text;
    if (!SAMReader::parseHeader(fileName, header->text(), header->text() + textHeaderSize, genome, &textHeaderSize)) {
        fprintf(stderr,"BAMReader: failed to parse header on '%s'\n",fileName);
        return false;
    }
    n_ref = header->n_ref();
    refOffset = new unsigned[n_ref];
    BAMHeaderRefSeq* refSeq = header->firstRefSeq();
    for (unsigned i = 0; i < n_ref; i++, refSeq = refSeq->next()) {
        if (! genome->getOffsetOfPiece(refSeq->name(), &refOffset[i])) {
            // fprintf(stderr, "BAMReader: unknown ref seq name %s\n", refSeq->name());
            refOffset[i] = UINT32_MAX;
            // exit(1); ??
        }
    }
    data->reinit(data->getFileOffset(), amountOfFileToProcess == 0 ? 0 : amountOfFileToProcess - data->getFileOffset());
    extraOffset = 0;
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
    data->reinit(startingOffset, amountOfFileToProcess);
    extraOffset = 0;
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
    ReadClippingType clipping,
    int matchBufferSize)
{
    BAMReader* reader = create(fileName, genome, 0, 0, clipping);
    PairedReadReader* matcher = PairedReadReader::PairMatcher(5000, reader);
    ReadSupplierQueue* queue = new ReadSupplierQueue(matcher);
    return queue;
}
    
const char* BAMAlignment::CodeToSeq = "=ACMGRSVTWYHKDBN";
_uint8 BAMAlignment::SeqToCode[256];
const char* BAMAlignment::CodeToCigar = "MIDNSHP=X";
_uint8 BAMAlignment::CigarToCode[256];

BAMAlignment::_init BAMAlignment::_init_;

    void
BAMAlignment::decodeSeq(
    char* o_sequence,
    _uint8* nibbles,
    int bases)
{
    for (int i = 0; i < bases; i++) {
        int n = (i & 1) ? *nibbles & 0xf : *nibbles >> 4;
        nibbles += i & 1;
        *o_sequence++ = BAMAlignment::CodeToSeq[n];
    }
}
    
    void
BAMAlignment::decodeQual(
    char* o_qual,
    char* quality,
    int bases)
{
    for (int i = 0; i < bases; i++) {
        char q = quality[i];
        o_qual[i] = q < 0 || q >= 64 ? '!' : q + '!';
    }
}

    bool
BAMAlignment::decodeCigar(
    char* o_cigar,
    int cigarSize,
    _uint32* cigar,
    int ops)
{
    int i = 0;
    while (ops > 0 && i < cigarSize - 11) { // 9 decimal digits (28 bits) + 1 cigar char + null terminator
        i += sprintf(o_cigar + i, "%u", *cigar >> 4);
        _ASSERT((*cigar & 0xf) <= 8);
        o_cigar[i++] = BAMAlignment::CodeToCigar[*cigar & 0xf];
        ops--;
        cigar++;
    }
    o_cigar[i++] = 0;
    return ops == 0;
}

    void
BAMAlignment::encodeSeq(
    _uint8* encoded,
    char* ascii,
    int length)
{
    _uint8* p = encoded;
    for (int i = 0; i + 1 < length; i += 2) {
        *p++ = (BAMAlignment::SeqToCode[ascii[i]] << 4) | BAMAlignment::SeqToCode[ascii[i+1]];
    }
    if (length % 2) {
        *p = BAMAlignment::SeqToCode[ascii[length - 1]] << 4;
    }
}

// static initializer
BAMAlignment::_init::_init()
{
    memset(SeqToCode, 0, 256);
    for (int i = 1; i < 16; i++) {
        SeqToCode[CodeToSeq[i]] = i;
    }
    memset(CigarToCode, 0, 256);
    for (int i = 1; i < 9; i++) {
        CigarToCode[CodeToCigar[i]] = i;
    }
}

    int
BAMAlignment::reg2bin(
    int beg,
    int end)
{
    --end;
    if (beg>>14 == end>>14) return ((1<<15)-1)/7 + (beg>>14);
    if (beg>>17 == end>>17) return ((1<<12)-1)/7 + (beg>>17);
    if (beg>>20 == end>>20) return ((1<<9)-1)/7 + (beg>>20);
    if (beg>>23 == end>>23) return ((1<<6)-1)/7 + (beg>>23);
    if (beg>>26 == end>>26) return ((1<<3)-1)/7 + (beg>>26);
    return 0;
}

    int
BAMAlignment::reg2bins(
    int beg,
    int end,
    _uint16* list)
{
    int i = 0, k;
    --end;
    list[i++] = 0;
    for (k = 1 + (beg>>26); k <= 1 + (end>>26); ++k) list[i++] = k;
    for (k = 9 + (beg>>23); k <= 9 + (end>>23); ++k) list[i++] = k;
    for (k = 73 + (beg>>20); k <= 73 + (end>>20); ++k) list[i++] = k;
    for (k = 585 + (beg>>17); k <= 585 + (end>>17); ++k) list[i++] = k;
    for (k = 4681 + (beg>>14); k <= 4681 + (end>>14); ++k) list[i++] = k;
    return i;
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
    _int64 bytes;
    if (! data->getData(&buffer, &bytes)) {
        data->nextBatch();
        if (! data->getData(&buffer, &bytes)) {
            return false;
        }
    }
    BAMAlignment* bam = (BAMAlignment*) buffer;
    if ((unsigned _int64)bytes < sizeof(bam->block_size) || (unsigned _int64)bytes < bam->size()) {
        fprintf(stderr, "Unexpected end of BAM file at %lld\n", data->getFileOffset());
        exit(1);
    }
    static size_t lastSize = bam->size();//!!
    data->advance(bam->size());
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

    if (NULL != read) {
        _ASSERT(bam->l_seq < MAX_SEQ_LENGTH);
        char* seqBuffer = getExtra(bam->l_seq);
        char* qualBuffer = getExtra(bam->l_seq);
        BAMAlignment::decodeSeq(seqBuffer, bam->seq(), bam->l_seq);
        BAMAlignment::decodeQual(qualBuffer, bam->qual(), bam->l_seq);
        read->init(bam->read_name(), bam->l_read_name - 1, seqBuffer, qualBuffer, bam->l_seq);
        read->setBatch(data->getBatch());
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
        char* cigarBuffer = getExtra(MAX_SEQ_LENGTH);
        if (! BAMAlignment::decodeCigar(cigarBuffer, MAX_SEQ_LENGTH, bam->cigar(), bam->n_cigar_op)) {
            *cigar = ""; // todo: fail?
        }
        *cigar = cigarBuffer;
    }
}

    char*
BAMReader::getExtra(
    _int64 bytes)
{
    char* extra;
    _int64 limit;
    data->getExtra(&extra, &limit);
    _ASSERT(extra != NULL && bytes >= 0 && limit - extraOffset >= 2 * bytes);
    char* result = extra + extraOffset;
    extraOffset += max((_int64) 0, bytes);
    return result;
}

    
class BAMFormat : public FileFormat
{
public:
    BAMFormat(bool i_useM) : useM(i_useM) {}

    virtual bool isFormatOf(const char* filename) const;
    
    virtual ReadWriterSupplier* getWriterSupplier(AlignerOptions* options, const Genome* genome) const;

    virtual bool writeHeader(
        const Genome *genome, char *header, size_t headerBufferSize, size_t *headerActualSize,
        bool sorted, int argc, const char **argv, const char *version, const char *rgLine) const;

    virtual bool writeRead(
        const Genome * genome, LandauVishkinWithCigar * lv, char * buffer, size_t bufferSpace, 
        size_t * spaceUsed, size_t qnameLen, Read * read, AlignmentResult result, 
        unsigned genomeLocation, Direction direction,
        bool hasMate = false, bool firstInPair = false, Read * mate = NULL, 
        AlignmentResult mateResult = NotFound, unsigned mateLocation = 0, Direction mateDirection = FORWARD) const; 

private:

    static int computeCigarOps(const Genome * genome, LandauVishkinWithCigar * lv,
        char * cigarBuf, int cigarBufLen,
        const char * data, unsigned dataLength, unsigned basesClippedBefore, unsigned basesClippedAfter,
        unsigned genomeLocation, bool isRC, bool useM, int * editDistance);

    const bool useM;
};

const FileFormat* FileFormat::BAM[] = { new BAMFormat(false), new BAMFormat(true) };

    bool
BAMFormat::isFormatOf(
    const char* filename) const
{
    return util::stringEndsWith(filename, ".bam");
}
     
    ReadWriterSupplier*
BAMFormat::getWriterSupplier(
    AlignerOptions* options,
    const Genome* genome) const
{
    DataWriterSupplier* dataSupplier;
    if (options->sortOutput) {
        int len = strlen(options->outputFileTemplate);
        // todo: this is going to leak, but there's no easy way to free it, and it's small...
        char* tempFileName = (char*) malloc(5 + len);
        strcpy(tempFileName, options->outputFileTemplate);
        strcpy(tempFileName + len, ".tmp");
        dataSupplier = DataWriterSupplier::sorted(tempFileName, options->outputFileTemplate, DataWriterSupplier::gzip());
    } else {
        dataSupplier = DataWriterSupplier::create(options->outputFileTemplate, DataWriterSupplier::gzip(), 3);
    }
    return ReadWriterSupplier::create(this, dataSupplier, genome);
}

    bool
BAMFormat::writeHeader(
    const Genome *genome,
    char *header,
    size_t headerBufferSize,
    size_t *headerActualSize,
    bool sorted,
    int argc,
    const char **argv,
    const char *version,
    const char *rgLine) const
{
    if (headerBufferSize < BAMHeader::size(0)) {
        return false;
    }
    size_t cursor = 0;
    BAMHeader* bamHeader = (BAMHeader*) header;
    bamHeader->magic = BAMHeader::BAM_MAGIC;
    size_t samHeaderSize;
    bool ok = FileFormat::SAM[0]->writeHeader(genome, bamHeader->text(), headerBufferSize - BAMHeader::size(0), &samHeaderSize,
        sorted, argc, argv, version, rgLine);
    if (! ok) {
        return false;
    }
    bamHeader->l_text = samHeaderSize;
    cursor = BAMHeader::size(samHeaderSize);

    // Write a RefSeq record for each chromosome / piece in the genome
    const Genome::Piece *pieces = genome->getPieces();
    int numPieces = genome->getNumPieces();
    bamHeader->n_ref() = numPieces;
    BAMHeaderRefSeq* refseq = bamHeader->firstRefSeq();
    unsigned genomeLen = genome->getCountOfBases();
    for (int i = 0; i < numPieces; i++) {
        int len = strlen(pieces[i].name) + 1;
        cursor += BAMHeaderRefSeq::size(len);
        if (cursor > headerBufferSize) {
            return false;
        }
        refseq->l_name = len;
        memcpy(refseq->name(), pieces[i].name, len);
        unsigned start = pieces[i].beginningOffset;
        unsigned end = (i + 1 < numPieces) ? pieces[i+1].beginningOffset : genomeLen;
        refseq->l_ref() = end - start;
        refseq = refseq->next();
        _ASSERT((char*) refseq - header == cursor);
    }
    *headerActualSize = cursor;
    return true;
}

    bool
BAMFormat::writeRead(
    const Genome * genome,
    LandauVishkinWithCigar * lv,
    char * buffer,
    size_t bufferSpace, 
    size_t * spaceUsed,
    size_t qnameLen,
    Read * read,
    AlignmentResult result, 
    unsigned genomeLocation,
    Direction direction,
    bool hasMate,
    bool firstInPair,
    Read * mate, 
    AlignmentResult mateResult,
    unsigned mateLocation,
    Direction mateDirection) const
{
    const int MAX_READ = 10000;
    const int cigarBufSize = MAX_READ * 2;
    char cigarBuf[cigarBufSize];

    int flags = 0;
    const char *pieceName = "*";
    int pieceIndex = -1;
    unsigned positionInPiece = 0;
    int mapQuality = 0;
    int cigarOps;
    const char *matePieceName = "*";
    unsigned matePositionInPiece = 0;
    _int64 templateLength = 0;

    char data[MAX_READ];
    char quality[MAX_READ];

    const char* clippedData;
    unsigned fullLength;
    unsigned clippedLength;
    unsigned basesClippedBefore;
    unsigned basesClippedAfter;
    int editDistance;

    if (! getSAMData(genome, lv, data, quality, MAX_READ, pieceName, pieceIndex, 
        flags, positionInPiece, mapQuality, matePieceName, matePositionInPiece, templateLength,
        fullLength, clippedData, clippedLength, basesClippedBefore, basesClippedAfter,
        qnameLen, read, result, genomeLocation, direction, useM,
        hasMate, firstInPair, mate, mateResult, mateLocation, mateDirection))
    {
        return false;
    }
    if (genomeLocation != 0xFFFFFFFF) {
        cigarOps = computeCigarOps(genome, lv, cigarBuf, cigarBufSize,
                                   clippedData, clippedLength, basesClippedBefore, basesClippedAfter,
                                   genomeLocation, direction, useM, &editDistance);
    }

    // Write the BAM entry
    size_t bamSize = BAMAlignment::size(qnameLen + 1, cigarOps, fullLength);
    if (bamSize > bufferSpace) {
        return false;
    }
    BAMAlignment* bam = (BAMAlignment*) buffer;
    bam->block_size = bamSize - 4;
    bam->refID = pieceIndex;
    bam->pos = positionInPiece - 1;
    bam->l_read_name = qnameLen + 1;
    bam->MAPQ = mapQuality;
    // todo: what is bin for unmapped reads?
    bam->bin = genomeLocation != 0xFFFFFFFF ? BAMAlignment::reg2bin(genomeLocation, genomeLocation + fullLength) : 0;
    bam->n_cigar_op = cigarOps;
    bam->FLAG = flags;
    bam->l_seq = fullLength;
    if (mateLocation != 0xFFFFFFFF && mateResult != NotFound) {
        const Genome::Piece* matePiece = genome->getPieceAtLocation(mateLocation);
        _ASSERT(matePiece != NULL);
        bam->next_refID = matePiece - genome->getPieces();
        bam->next_pos = mateLocation - matePiece->beginningOffset;
    } else {
        bam->next_refID = -1;
        bam->next_pos = -1;
    }
    bam->tlen = templateLength;
    memcpy(bam->read_name(), read->getId(), qnameLen);
    bam->read_name()[qnameLen] = 0;
    memcpy(bam->cigar(), cigarBuf, cigarOps * 4);
    BAMAlignment::encodeSeq(bam->seq(), data, fullLength);
    for (int i = 0; i < fullLength; i++) {
        quality[i] -= '!';
    }
    memcpy(bam->qual(), quality, fullLength);

    if (NULL != spaceUsed) {
        *spaceUsed = bamSize;
    }
    return true;
}

// Compute the CIGAR edit sequence operations in BAM format for a read against a given genome location
// Returns number of operations (or 0 if there was a problem)
    int
BAMFormat::computeCigarOps(
    const Genome *              genome,
    LandauVishkinWithCigar *    lv,
    char *                      cigarBuf,
    int                         cigarBufLen,
    const char *                data,
    unsigned                    dataLength,
    unsigned                    basesClippedBefore,
    unsigned                    basesClippedAfter,
    unsigned                    genomeLocation,
    bool                        isRC,
	bool						useM,
    int *                       editDistance
)
{
    const char *reference = genome->getSubstring(genomeLocation, dataLength);
    int used;
    if (NULL != reference) {
        *editDistance = lv->computeEditDistance(
                            genome->getSubstring(genomeLocation, dataLength),
                            dataLength,
                            data,
                            dataLength,
                            MAX_K - 1,
                            cigarBuf + 4 * (basesClippedBefore > 0),
                            cigarBufLen - 4 * ((basesClippedBefore > 0) + (basesClippedAfter > 0)),
						    useM, BAM_CIGAR_OPS, &used);
    } else {
        //
        // Fell off the end of the chromosome.
        //
        return 0;
    }

    if (*editDistance == -2) {
        fprintf(stderr, "WARNING: computeEditDistance returned -2; cigarBuf may be too small\n");
        return 0;
    } else if (*editDistance == -1) {
        static bool warningPrinted = false;
        if (!warningPrinted) {
            fprintf(stderr, "WARNING: computeEditDistance returned -1; this shouldn't happen\n");
            warningPrinted = true;
        }
        return 0;
    } else {
        // Add some CIGAR instructions for soft-clipping if we've ignored some bases in the read.
        if (basesClippedBefore > 0) {
            *(_uint32*)cigarBuf = (basesClippedBefore << 4) | BAMAlignment::CigarToCode['S'];
            used += 4;
        }
        if (basesClippedAfter > 0) {
            *(_uint32*)(cigarBuf + used) = (basesClippedAfter << 4) | BAMAlignment::CigarToCode['S'];
            used += 4;
        }
        return used / 4;
    }
}
