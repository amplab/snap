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
#include "SAM.h"
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
#include "exit.h"
#include "VariableSizeMap.h"
#include "PairedAligner.h"
#include "GzipDataWriter.h"
#include "Error.h"

using std::max;
using std::min;
using util::strnchr;

GenomeLocation 
BAMAlignment::getLocation(const BAMReader * bamReader) const
{
    return bamReader == NULL || pos < 0 || refID < 0 || (FLAG & SAM_UNMAPPED)
        ? InvalidGenomeLocation : (bamReader->getLocationForRefAndOffset(refID, pos));
}

GenomeLocation 
BAMAlignment::getNextLocation(const BAMReader * bamReader) const
{
    return next_pos < 0 || next_refID < 0 || (FLAG & SAM_NEXT_UNMAPPED)
        ? InvalidGenomeLocation : (bamReader->getLocationForRefAndOffset(next_refID, next_pos));
}


BAMReader::BAMReader(const ReaderContext& i_context) : ReadReader(i_context), refLocation(0), n_ref(0), data(NULL),
    numRGLines(0), rgLines(NULL), rgLineOffsets(NULL)
{
}

BAMReader::~BAMReader()
{
    BigDealloc(refLocation);
    refLocation = NULL;

    if (NULL != data) {
        delete data;
    }
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
BAMReader::init(
    const char *fileName,
    int bufferCount,
    _int64 startingOffset,
    _int64 amountOfFileToProcess)
{
    // todo: integrate supplier models
    // might need up to 3x extra for expanded sequence + quality + cigar data
    if (!strcmp("-", fileName)) {
        data = DataSupplier::GzipBamStdio->getDataReader(bufferCount, MAX_RECORD_LENGTH, 3.0 * DataSupplier::ExpansionFactor, 0);
    } else {
        data = DataSupplier::GzipBamDefault->getDataReader(bufferCount, MAX_RECORD_LENGTH, 3.0 * DataSupplier::ExpansionFactor, 0);
    }

    if (! data->init(fileName)) {
        WriteErrorMessage("Unable to read file %s\n", fileName);
        soft_exit(1);
    }

    if (startingOffset == 0) {
        readHeader(fileName);
    }

    _ASSERT(context.headerBytes > 0);
    reinit(startingOffset, amountOfFileToProcess);
    if ((size_t) startingOffset < context.headerBytes) {
		_int64 bytesToSkip = context.headerBytes - startingOffset;

		while (bytesToSkip > 0) {
			char* p;
			_int64 valid, start;
			bool ok = data->getData(&p, &valid, &start);
			if (!ok) {
				WriteErrorMessage("failure reading file %s\n", fileName);
				soft_exit(1);
			}

			_int64 bytesToSkipThisTime = __min(valid, bytesToSkip);
			data->advance(bytesToSkipThisTime);
			if (bytesToSkipThisTime > start) {
				data->nextBatch();
			}
			data->getData(&p, &valid, &start);

			bytesToSkip -= bytesToSkipThisTime;
		}
    }
}

    void
BAMReader::readHeader(
    const char* fileName)
{
    _ASSERT(context.header == NULL);
    _int64 headerSize = 1024 * 1024; // 1M header initially 

	bool sawWholeHeader;
	BAMHeader* header;
	_int64 textHeaderSize;
	char* buffer;

	buffer = data->readHeader(&headerSize);

	if (headerSize < sizeof(BAMHeader)) {
		WriteErrorMessage("Malformed BAM file '%s', too small to conatain even a header.\n", fileName);
		soft_exit(1);
	}

	header = (BAMHeader*)buffer;
	if (header->magic != BAMHeader::BAM_MAGIC) {
		WriteErrorMessage("BAMReader: Not a valid BAM file\n");
		soft_exit(1);
	}
	textHeaderSize = header->l_text;

	if (textHeaderSize + (_int64)sizeof(BAMHeader) > headerSize) {
		headerSize = textHeaderSize + (_int64)sizeof(BAMHeader);
		buffer = data->readHeader(&headerSize);
		if (textHeaderSize + (_int64)sizeof(BAMHeader) >  headerSize) {
			WriteErrorMessage("Unable to read entire header of BAM file '%s', it may be malformed.\n", fileName);
			soft_exit(1);
		}

		header = (BAMHeader*)buffer;
		if (header->magic != BAMHeader::BAM_MAGIC) {
			WriteErrorMessage("BAMReader: Not a valid BAM file\n");
			soft_exit(1);
		}

		_ASSERT(textHeaderSize == header->l_text);	// We got the same thing this time
	}

	if (!SAMReader::parseHeader(fileName, header->text(), header->text() + headerSize - sizeof(BAMHeader), context.genome, &textHeaderSize, &context.headerMatchesIndex, &sawWholeHeader, &n_ref, &refLocation, 
        &numRGLines, &rgLines, &rgLineOffsets)) {
		WriteErrorMessage("BAMReader: failed to parse header on '%s'\n", fileName);
		soft_exit(1);
	}

	if (!sawWholeHeader) {
		WriteErrorMessage("We had the entire header loaded for file '%s', but it didn't parse correctly\n", fileName);
		soft_exit(1);
	}

    if (header->n_ref() != n_ref) { // We got the same value from the SAM header parser as is in the BAM header
        WriteErrorMessage("Truncated or corrupt BAM file near offset %lld\n", data->getFileOffset());
        soft_exit(1);
    }
	BAMHeaderRefSeq* refSeq = header->firstRefSeq();
	for (int i = 0; i < n_ref; i++, refSeq = refSeq->next()) {
		// just advance
	}
	
	char* p = new char[textHeaderSize + 1];
    memcpy(p, header->text(), textHeaderSize);
    p[textHeaderSize] = 0;
    context.header = p;
    context.headerLength = textHeaderSize;
    context.headerBytes = (char*) refSeq - buffer;
    context.numRGLines = numRGLines;
    context.rgLines = rgLines;
    context.rgLineOffsets = rgLineOffsets;
}

    BAMReader*
BAMReader::create(
    const char *fileName,
    int bufferCount,
    _int64 startingOffset,
    _int64 amountOfFileToProcess,
    const ReaderContext& context)
{
    BAMReader* reader = new BAMReader(context);
    reader->init(fileName, bufferCount, startingOffset, amountOfFileToProcess);
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
    const ReaderContext& context)
{
    BAMReader* reader = create(fileName, ReadSupplierQueue::BufferCount(numThreads), 0, 0, context);
    ReadSupplierQueue* queue = new ReadSupplierQueue((ReadReader*)reader);
    queue->startReaders();
    return queue;
}

    PairedReadSupplierGenerator *
BAMReader::createPairedReadSupplierGenerator(
    const char *fileName,
    int numThreads,
    bool quicklyDropUnmatchedReads,
    const ReaderContext& context,
    int matchBufferSize)
{
    BAMReader* reader = create(fileName, 
        ReadSupplierQueue::BufferCount(numThreads) + PairedReadReader::MatchBuffers, 0, 0, context);
    PairedReadReader* matcher = PairedReadReader::PairMatcher(reader, quicklyDropUnmatchedReads);
    ReadSupplierQueue* queue = new ReadSupplierQueue(matcher);
    queue->startReaders();
    return queue;
}

const char* BAMAlignment::CodeToSeq =   "=ACMGRSVTWYHKDBN";
const char *BAMAlignment::CodeToSeqRC = "NTGKCYWBASRDMHVN"; // Bill's best guess for things other than ATCG, not that it matters for SNAP
_uint16 BAMAlignment::CodeToSeqPair[256];
_uint16 BAMAlignment::CodeToSeqPairRC[256];
_uint8 BAMAlignment::SeqToCode[256];
const char* BAMAlignment::CodeToCigar = "MIDNSHP=X";
_uint8 BAMAlignment::CigarToCode[256];
_uint8 BAMAlignment::CigarCodeToRefBase[9] = {1, 0, 1, 1, 0, 0, 1, 1, 1};

const _uint8 BAM_CIGAR_M = 0;
const _uint8 BAM_CIGAR_I = 1;
const _uint8 BAM_CIGAR_D = 2;
const _uint8 BAM_CIGAR_N = 3;
const _uint8 BAM_CIGAR_S = 4;
const _uint8 BAM_CIGAR_H = 5;
const _uint8 BAM_CIGAR_P = 6;
const _uint8 BAM_CIGAR_EQUAL = 7;
const _uint8 BAM_CIGAR_X = 8;

BAMAlignment::_init BAMAlignment::_init_;

    void
BAMAlignment::decodeSeq(
    char* o_sequence,
    const _uint8* nibbles,
    int bases)
{


    _uint16 *o_sequence_pairs = (_uint16 *)o_sequence;
    int pairs = bases / 2;
    for (int i = 0; i < pairs; i++) {
        o_sequence_pairs[i] = CodeToSeqPair[nibbles[i]];
    }

    if (bases % 2 == 1) {
        o_sequence[bases - 1] = CodeToSeq[nibbles[bases / 2] >> 4];
    }

#ifdef _DEBUG   // Make sure the new one does the same thing as the old.
    for (int i = 0; i < bases; i++) {
        int bit = 1 ^ (i & 1);
        int n = (*nibbles >> (bit << 2)) & 0xf; // extract nibble without branches
        nibbles += 1 ^ bit;
        _ASSERT(o_sequence[i] == BAMAlignment::CodeToSeq[n]);
    }
#endif // _DEBUG
}
    void
BAMAlignment::decodeSeqRC(
char* o_sequence,
const _uint8* nibbles,
int bases)
{
    _uint16 *o_sequence_pairs = (_uint16 *)&o_sequence[bases % 2];
    int pairs = bases / 2;
    for (int i = 0; i < pairs; i++) {
        o_sequence_pairs[pairs-i-1] = CodeToSeqPairRC[nibbles[i]];
    }

    if (bases % 2 == 1) {
        o_sequence[0] = CodeToSeqRC[nibbles[bases / 2] >> 4];
    }
}
    void
BAMAlignment::decodeQual(
    char* o_qual,
    char* quality,
    int bases)
{
    for (int i = 0; i < bases; i++) {
        o_qual[i] = CIGAR_QUAL_TO_SAM[((_uint8*)quality)[i]];
    }
}

    void
BAMAlignment::decodeQualRC(
    char* o_qual,
    char* quality,
    int bases)
{
    for (int i = 0; i < bases; i++) {
        o_qual[bases-i-1] = CIGAR_QUAL_TO_SAM[((_uint8*)quality)[i]];
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
    _uint32 lastOp = 99999;
    while (ops > 0 && i < cigarSize - 11) { // 9 decimal digits (28 bits) + 1 cigar char + null terminator
        i += sprintf(o_cigar + i, "%u", *cigar >> 4);
        _ASSERT((*cigar & 0xf) <= 8);
        _uint32 op = *cigar & 0xf;
        o_cigar[i++] = BAMAlignment::CodeToCigar[op];
        _ASSERT(op != lastOp);
        lastOp = op;
        ops--;
        cigar++;
    }
    o_cigar[i++] = 0;
    return ops == 0;
}

    void 
BAMAlignment::getClippingFromCigar(
    _uint32 *cigar, 
    int ops, 
    unsigned *o_frontClipping, 
    unsigned *o_backClipping, 
    unsigned *o_frontHardClipping, 
    unsigned *o_backHardClipping)
{
    *o_frontHardClipping = 0;   // Gets overwritten if we have any
    *o_frontClipping = 0;
    *o_backHardClipping = 0;
    *o_backClipping = 0;

    if (0 == ops) return;

    if ((*cigar & 0xf) == BAM_CIGAR_H) {
        *o_frontHardClipping = *cigar >> 4;
        cigar++;
        ops--;
        if (0 == ops) {
            return;   // What a strange cigar string, all hard clip!
        }
    }

    if ((*cigar & 0xf) == BAM_CIGAR_S) {
        *o_frontClipping = *cigar >> 4;
        cigar++;
        ops--;
        if (0 == ops) {
            return;
        }
    }

    if ((cigar[ops - 1] & 0xf) == BAM_CIGAR_H) {
        *o_backHardClipping = cigar[ops - 1] >> 4;
        ops--;
        if (0 == ops) {
            return;
        }
    }

    if ((cigar[ops - 1] & 0xf) == BAM_CIGAR_S) {
        *o_backClipping = cigar[ops - 1] >> 4;
    }
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

    int
BAMAlignment::l_ref()
{
    if (FLAG & SAM_UNMAPPED) {
        return 0;
    }
    if (n_cigar_op == 0) {
        return l_seq;
    }
    _uint32* p = cigar();
    int len = 0;
    static const int op_ref[16] = {1, 0, 1, 1, 1, 1, 1, 1, 1, 0, 0, 0, 0, 0, 0, 0};
    for (int i = 0; i < n_cigar_op; i++) {
        _uint32 op = *p++;
        len += op_ref[(op & 15)] * (op >> 4);
    }
    return len;
}

    GenomeLocation
BAMAlignment::getUnclippedStart(GenomeLocation loc)
{
    if (FLAG & SAM_UNMAPPED) {
        return loc;
    }
    if (n_cigar_op == 0) {
        return loc;
    }
    _uint32* p = cigar();
    _uint32 op = *p & 0xf;
    int len = 0;
    if (op == BAM_CIGAR_H || op == BAM_CIGAR_S) {
        len += (*p >> 4);
    }
    return loc - len;
}

    GenomeLocation
BAMAlignment::getUnclippedEnd(GenomeLocation loc)
{
    if (FLAG & SAM_UNMAPPED) {
        return loc;
    }
    if (n_cigar_op == 0) {
        return loc;
    }
    _uint32* p = cigar();
    _uint32 op = *p & 0xf;
    int len = 0;
    if (op != BAM_CIGAR_H && op != BAM_CIGAR_S) {
        len += (*p >> 4);
    }
    p++;
    static const int op_ref[16] = {1, 0, 1, 1, 1, 1, 1, 1, 1, 0, 0, 0, 0, 0, 0, 0};
    for (int i = 1; i < n_cigar_op; i++) {
        _uint32 op = *p++;
        len += op_ref[(op & 15)] * (op >> 4);
    }
    return loc + len;
}

// static initializer
BAMAlignment::_init::_init()
{
    memset(SeqToCode, 0, 256);
    for (int i = 1; i < 16; i++) {
        SeqToCode[CodeToSeq[i]] = i;
    }

    for (int i = 0; i < 256; i++) {
        CodeToSeqPair[i] = CodeToSeq[i >> 4] | (CodeToSeq[i & 0xf] << 8); // If this looks backwards, recall that the machines are little-endian
        CodeToSeqPairRC[i] = (CodeToSeqRC[i >> 4] << 8) | CodeToSeqRC[i & 0xf]; // Doubled backwards == forward
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

#ifdef VALIDATE_BAM
    void
BAMAlignment::validate()
{
    if (block_size >= 0x100000) { // sanity check, should be <1MB!
        WriteErrorMessage("Read name: %.*s. Block size: %d greater than 1 MB", l_read_name, read_name(), block_size);
        soft_exit(1);
    }
    if (size(l_read_name, n_cigar_op, l_seq, 0) > block_size + sizeof(block_size)) {
        WriteErrorMessage("Read name: %.*s. Size of BAM record %lld larger than allocated %lld\n", l_read_name, read_name(), size(l_read_name, n_cigar_op, l_seq, 0), block_size + sizeof(block_size));
        soft_exit(1);
    }
    if (refID < -1 || refID >(int)0x100000) {
        WriteErrorMessage("Read name: %.*s. refID: %d out of range\n", l_read_name, read_name(), refID);
        soft_exit(1);
    }
    // todo: validate bin, requires more info
    if (MAPQ > 80 && MAPQ != 255) {
        WriteErrorMessage("Read name: %.*s. MAPQ %u incorrect\n", l_read_name, read_name(), MAPQ);
        soft_exit(1);
    }
    if (FLAG > 0x7fff) {
        WriteErrorMessage("Read name: %.*s. FLAG %u incorrect\n", l_read_name, read_name(), FLAG);
        soft_exit(1);
    }
    if (next_refID < -1 || next_refID >(int)0x100000) {
        WriteErrorMessage("Read name: %.*s. next_refID: %d out of range\n", l_read_name, read_name(), next_refID);
        soft_exit(1);
    }
    for (char* p = read_name(); p < read_name() + l_read_name - 1; p++) {
        if (*p < ' ' || *p > '~') {
            WriteErrorMessage("Read name: %.*s has strange character %c\n", l_read_name, read_name(), *p);
            soft_exit(1);
        }
    }
    if (read_name()[l_read_name - 1] != 0) {
        WriteErrorMessage("Read name: %.*s not null terminated.\n", l_read_name, read_name());
        soft_exit(1);
    }
    // can't validate seq, all values are valid (though some are unlikely!)
    char* q = qual();
    for (int i = 0; i < l_seq; i++) {
        if (q[i] < -1 && q[i] > 80) {
            WriteErrorMessage("Read name: %.*s. Base quality %c incorrect\n", l_read_name, read_name(), q[i]);
            soft_exit(1);
        }
    }
    BAMAlignAux* aux = firstAux();
    for (; (char*)aux - (char*)firstAux() < auxLen(); aux = aux->next()) {
        if (aux->tag[0] < ' ' || aux->tag[0] > '~' || aux->tag[1] < ' ' || aux->tag[1] > '~') {
            WriteErrorMessage("Read name: %.*s. Aux tag1 %c tag2 %c incorrect\n", l_read_name, read_name(), aux->tag[0], aux->tag[1]);
            soft_exit(1);
        }
        if (strchr("AcCsSiIfZHB", aux->val_type) == NULL) {
            WriteErrorMessage("Read name: %.*s. Aux val type %c incorrect\n", l_read_name, read_name(), aux->val_type);
            soft_exit(1);
        }
    }
    if ((char*)aux - (char*)firstAux() != auxLen()) {
        WriteErrorMessage("Read name: %.*s. Aux field length mismatch. Found %d. Expected %d\n", l_read_name, read_name(), (char*)aux - (char*)firstAux(), auxLen());
        soft_exit(1);
    }
}
#endif

    bool
BAMReader::getNextRead(
    Read *read,
    AlignmentResult *alignmentResult,
    GenomeLocation *genomeLocation,
    bool *isRC,
    unsigned *mapQ,
    unsigned *flag,
    bool ignoreEndOfRange,
    const char **cigar)
{
    unsigned local_flag;
    if (NULL == flag) {
        flag = &local_flag;
    }

    do {
        char* buffer;
        _int64 bytes;
        if (! data->getData(&buffer, &bytes)) {
            data->nextBatch();
            if (! data->getData(&buffer, &bytes)) {
                return false;
            }
            extraOffset = 0;
        }
        BAMAlignment* bam = (BAMAlignment*) buffer;
        if ((_uint64)bytes < sizeof(bam->block_size) || (_uint64)bytes < bam->size()) {
			WriteErrorMessage("Truncated or corrupt BAM file near offset %lld\n", data->getFileOffset());
            soft_exit(1);
        }
        data->advance(bam->size());
        size_t lineLength;
        getReadFromLine(context.genome, buffer, buffer + bytes, read, alignmentResult, genomeLocation,
            isRC, mapQ, &lineLength, flag, cigar, context.clipping);
        unsigned auxLen = bam->auxLen();
        read->setReadGroup(context.defaultReadGroup);
        if (auxLen > 0) {
            read->setAuxiliaryData((char*) bam->firstAux(), auxLen);
            char* rgFromAux = NULL;
            int rgFromAuxLen = 0;
            for (BAMAlignAux* aux = bam->firstAux(); aux < bam->endAux(); aux = aux->next()) {
                if (aux->val_type == 'Z' && aux->tag[0] == 'R' && aux->tag[1] == 'G') {
                    rgFromAux = (char*)aux->value();
                    rgFromAuxLen = strlen(rgFromAux);
                    read->setReadGroup(READ_GROUP_FROM_AUX);
                    break;
                }
            }

            // LB
            if (context.rgLines != NULL) {
                // get library for read group
                for (int i = 0; i < context.numRGLines; i++) {
                    char* rgStart = context.rgLines + context.rgLineOffsets[i];
                    char* rgEnd = strchr(rgStart, '\t');
                    int rgLen = (int)(rgEnd - rgStart);
                    if (rgFromAuxLen != rgLen) continue;
                    if (!strncmp(rgFromAux, rgStart, rgLen)) {
                        char* lbStart = rgEnd + 1;
                        char* lbEnd = strchr(lbStart, '\t');
                        read->setLibrary(lbStart);
                        read->setLibraryLength((int)(lbEnd - lbStart));
                    }
                }
            }
        }
    } while ((context.ignoreSecondaryAlignments && (*flag & SAM_SECONDARY)) || 
             (context.ignoreSupplementaryAlignments && (*flag & SAM_SUPPLEMENTARY)));
    _ASSERT(read->getData()[0]);
    return true;
}

    void
BAMReader::getReadFromLine(
    const Genome *genome,
    char *line,
    char *endOfBuffer,
    Read *read,
    AlignmentResult *alignmentResult,
    GenomeLocation *out_genomeLocation,
    bool *isRC,
    unsigned *mapQ,
    size_t *lineLength,
    unsigned *flag,
    const char **cigar,
    ReadClippingType clipping)
{
    _ASSERT(endOfBuffer - line >= sizeof(BAMHeader));
    BAMAlignment* bam = (BAMAlignment*) line;
    if ((size_t)(endOfBuffer - line) < bam->size()) {
        WriteErrorMessage("Truncated or corrupt BAM file near offset %lld\n", data->getFileOffset());
        soft_exit(1);
    }
    bam->validate();

    GenomeLocation genomeLocation = bam->getLocation(this);

    if (NULL != out_genomeLocation) {
        _ASSERT(-1 <= bam->refID && bam->refID < (int)genome->getNumContigs());
        *out_genomeLocation = genomeLocation;
    }

    if (NULL != cigar) {
        const char* cigarBuffer;
        {
            char *writableCigarBuffer = getExtra(min(MAX_K * 5, MAX_SEQ_LENGTH));
            if (!BAMAlignment::decodeCigar(writableCigarBuffer, MAX_SEQ_LENGTH, bam->cigar(), bam->n_cigar_op)) {
                cigarBuffer = ""; // todo: fail?
            }
            else {
                cigarBuffer = writableCigarBuffer;
            }
        }
        *cigar = cigarBuffer;
    }

    if (NULL != read) {
        if (bam->l_seq > MAX_SEQ_LENGTH) {
            WriteErrorMessage("Truncated or corrupt BAM file near offset %lld\n", data->getFileOffset());
            soft_exit(1);
        }
		char* seqBuffer = getExtra(bam->l_seq);
        char* qualBuffer = getExtra(bam->l_seq);

        unsigned originalFrontClipping, originalBackClipping, originalFrontHardClipping, originalBackHardClipping;

        if (bam->FLAG & SAM_REVERSE_COMPLEMENT) {
            BAMAlignment::decodeSeqRC(seqBuffer, bam->seq(), bam->l_seq);
            BAMAlignment::decodeQualRC(qualBuffer, bam->qual(), bam->l_seq);

            //
            // Get the clipping, but reverse the outputs front/back because this is an RC read.
            //
            BAMAlignment::getClippingFromCigar(bam->cigar(), bam->n_cigar_op, &originalBackClipping, &originalFrontClipping, &originalBackHardClipping, &originalFrontHardClipping);
        } else {
            BAMAlignment::decodeSeq(seqBuffer, bam->seq(), bam->l_seq);
            BAMAlignment::decodeQual(qualBuffer, bam->qual(), bam->l_seq);

            BAMAlignment::getClippingFromCigar(bam->cigar(), bam->n_cigar_op, &originalFrontClipping, &originalBackClipping, &originalFrontHardClipping, &originalBackHardClipping);
        }

        const char *rnext;
        unsigned rnextLen;
        if (bam->next_refID < 0 || (genome != NULL && bam->next_refID >= genome->getNumContigs())) {
            rnext = "*";
            rnextLen = 1;
        } else {
            rnext = genome->getContigs()[bam->next_refID].name;
            rnextLen = genome->getContigs()[bam->next_refID].nameLength;
        }
        read->init(bam->read_name(), bam->l_read_name - 1, seqBuffer, qualBuffer, bam->l_seq, genomeLocation, bam->MAPQ, bam->FLAG,
            originalFrontClipping, originalBackClipping, originalFrontHardClipping, originalBackHardClipping, rnext, rnextLen, bam->next_pos + 1, true);
        read->setBatch(data->getBatch());
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


}

    char*
BAMReader::getExtra(
    _int64 bytes)
{
    char* extra;
    _int64 limit;
    data->getExtra(&extra, &limit);
    _ASSERT(extra != NULL && bytes >= 0 && limit - extraOffset >= bytes);
    if (limit - extraOffset < bytes) {
        WriteErrorMessage("error: not enough space for expanding BAM file - increase expansion factor, currently -xf %.1f\n", DataSupplier::ExpansionFactor);
        soft_exit(1);
    }
    char* result = extra + extraOffset;
    extraOffset += max((_int64) 0, bytes);
    return result;
}


class BAMFormat : public FileFormat
{
public:
    BAMFormat(bool i_useM) : useM(i_useM) {}

    virtual void getSortInfo(const Genome* genome, char* buffer, _int64 bytes, GenomeLocation* o_location, GenomeDistance* o_readBytes, int* o_refID, int* o_pos) const;

    virtual void setupReaderContext(AlignerOptions* options, ReaderContext* readerContext) const
    { FileFormat::setupReaderContext(options, readerContext, true); }

    virtual ReadWriterSupplier* getWriterSupplier(AlignerOptions* options, const Genome* genome) const;

    virtual bool writeHeader(
        const ReaderContext& context, char *header, size_t headerBufferSize, size_t *headerActualSize,
        bool sorted, int argc, const char **argv, const char *version, const char *rgLine, bool omitSQLines) const;

    virtual bool writePairs(
        const ReaderContext& context, LandauVishkinWithCigar * lv, AffineGapVectorizedWithCigar * ag, 
        bool useAffineGap, char * buffer, size_t bufferSpace,
        size_t * spaceUsed, size_t* qnameLen, Read ** reads, GenomeLocation* locations, PairedAlignmentResult* result,
        bool isSecondary, bool emitInternalScore, char *internalScoreTag, int * writeOrder,
        int* cumulativePositiveAddFrontClipping, bool * secondReadLocationChanged, bool * outOfSpace) const;

    virtual bool writeRead(
        const ReaderContext& context, LandauVishkinWithCigar * lv, char * buffer, size_t bufferSpace,
        size_t * spaceUsed, size_t qnameLen, Read * read, AlignmentResult result,
        int mapQuality, GenomeLocation genomeLocation, Direction direction, bool secondaryAlignment, bool supplementaryAlignment, int * o_addFrontClipping,
        int internalScore, bool emitInternalScore, char *internalScoreTag, int bpClippedBefore = 0, int bpClippedAfter = 0,
        bool hasMate = false, bool firstInPair = false, Read * mate = NULL,
        AlignmentResult mateResult = NotFound, GenomeLocation mateLocation = 0, Direction mateDirection = FORWARD,
        bool alignedAsPair = false, int mateBpClippedBefore = 0, int mateBpClippedAfter = 0) const;

    virtual bool writeRead(
        const ReaderContext& context, AffineGapVectorizedWithCigar * ag, char * buffer, size_t bufferSpace,
        size_t * spaceUsed, size_t qnameLen, Read * read, AlignmentResult result,
        int mapQuality, GenomeLocation genomeLocation, Direction direction, bool secondaryAlignment, bool supplementaryAlignment, int * o_addFrontClipping,
        int internalScore, bool emitInternalScore, char *internalScoreTag, int bpClippedBefore = 0, int bpClippedAfter = 0,
        bool hasMate = false, bool firstInPair = false, Read * mate = NULL,
        AlignmentResult mateResult = NotFound, GenomeLocation mateLocation = 0, Direction mateDirection = FORWARD,
        bool alignedAsPair = false, int mateBpClippedBefore = 0, int mateBpClippedAfter = 0) const;

private:

    static int computeCigarOps(const Genome * genome, LandauVishkinWithCigar * lv,
        char * cigarBuf, int cigarBufLen,
        const char * data, unsigned dataLength, unsigned basesClippedBefore, unsigned extraBasesClippedBefore, unsigned basesClippedAfter,
        unsigned frontHardClipping, unsigned backHardClipping,
        GenomeLocation genomeLocation, bool isRC, bool useM, int * o_editDistance, int * o_addFrontClipping,
        int * o_refSpan);

    static int computeCigarOps(const Genome * genome, AffineGapVectorizedWithCigar * ag,
        char * cigarBuf, int cigarBufLen,
        const char * data, unsigned dataLength, unsigned basesClippedBefore, unsigned extraBasesClippedBefore, unsigned basesClippedAfter,
        unsigned frontHardClipping, unsigned backHardClipping,
        GenomeLocation genomeLocation, bool isRC, bool useM, int * o_editDistance, int * o_addFrontClipping,
        int * o_refSpan);

    const bool useM;
};

const FileFormat* FileFormat::BAM[] = { new BAMFormat(false), new BAMFormat(true) };

    void
BAMFormat::getSortInfo(
    const Genome* genome,
    char* buffer,
    _int64 bytes,
    GenomeLocation* o_location,
	GenomeDistance* o_readBytes,
	int* o_refID,
	int* o_pos) const
{
    BAMAlignment* bam = (BAMAlignment*) buffer;
    _ASSERT((size_t) bytes >= sizeof(BAMAlignment) && bam->size() <= (size_t) bytes && bam->refID < genome->getNumContigs());
	if (o_location != NULL) {
		if (bam->refID < 0 || bam->refID >= genome->getNumContigs() || bam->pos < 0) {
			if (bam->next_refID < 0 || bam->next_refID > genome->getNumContigs() || bam->next_pos < 0) {
				*o_location = InvalidGenomeLocation;
			} else {
				*o_location = genome->getContigs()[bam->next_refID].beginningLocation + bam->next_pos;
			}
		} else {
			*o_location = genome->getContigs()[bam->refID].beginningLocation + bam->pos;
		}
	}
	if (o_readBytes != NULL) {
		*o_readBytes = (unsigned) bam->size();
	}
	if (o_refID != NULL) {
		*o_refID = bam->refID;
	}
	if (o_pos != NULL) {
		*o_pos = bam->pos;
	}
}

    ReadWriterSupplier*
BAMFormat::getWriterSupplier(
    AlignerOptions* options,
    const Genome* genome) const
{
    DataWriterSupplier* dataSupplier;
    GzipWriterFilterSupplier* gzipSupplier =
        DataWriterSupplier::gzip(true, BAM_BLOCK, max(1, options->numThreads - 1), false, options->sortOutput);
        // (leave a thread free for main, and let OS map threads to cores to allow system IO etc.)
    if (options->sortOutput) {
        char *tempFileName = DataWriterSupplier::generateSortIntermediateFilePathName(options);
        // todo: make markDuplicates optional?
        DataWriter::FilterSupplier* filters = gzipSupplier;
        if (!options->noIndex) {
            size_t len = strlen(options->outputFile.fileName);
            char* indexFileName = (char*)malloc(5 + len);
            strcpy(indexFileName, options->outputFile.fileName);
            strcpy(indexFileName + len, ".bai");
            filters = DataWriterSupplier::bamIndex(indexFileName, genome, gzipSupplier)->compose(filters);
        }
        if (! options->noDuplicateMarking) {
            filters = DataWriterSupplier::bamMarkDuplicates(genome)->compose(filters);
        }
        dataSupplier = DataWriterSupplier::sorted(this, genome, tempFileName,
            options->sortMemory * (1ULL << 30),
            options->numThreads, options->outputFile.fileName, filters, options->writeBufferSize,
            options->emitInternalScore, options->internalScoreTag,
            FileEncoder::gzip(gzipSupplier, options->numThreads, options->bindToProcessors));
    } else {
        dataSupplier = DataWriterSupplier::create(options->outputFile.fileName, options->writeBufferSize, options->emitInternalScore, options->internalScoreTag, gzipSupplier);
    }
    return ReadWriterSupplier::create(this, dataSupplier, genome, options->killIfTooSlow, options->emitInternalScore, options->internalScoreTag, options->ignoreAlignmentAdjustmentsForOm);
}

    bool
BAMFormat::writeHeader(
    const ReaderContext& context,
    char *header,
    size_t headerBufferSize,
    size_t *headerActualSize,
    bool sorted,
    int argc,
    const char **argv,
    const char *version,
    const char *rgLine,
	bool omitSQLines) const
{
	_ASSERT(!omitSQLines);	// This is just for SAM files, at least for now.

    if (headerBufferSize < BAMHeader::size(0)) {
        return false;
    }
    size_t cursor = 0;
    BAMHeader* bamHeader = (BAMHeader*) header;
    bamHeader->magic = BAMHeader::BAM_MAGIC;
    size_t samHeaderSize;
    bool ok = FileFormat::SAM[0]->writeHeader(context, bamHeader->text(), headerBufferSize - BAMHeader::size(0), &samHeaderSize,
        sorted, argc, argv, version, rgLine, omitSQLines);
    if (! ok) {
        return false;
    }
    bamHeader->l_text = (int)samHeaderSize;
    cursor = BAMHeader::size((int)samHeaderSize);

    // Write a RefSeq record for each chromosome / contig in the genome
    // todo: handle null genome index case - reparse header & translate into BAM
    bamHeader->n_ref() = 0; // in case of overflow or no genome
	if (context.genome != NULL) {
		const Genome::Contig *contigs = context.genome->getContigs();
		int numContigs = context.genome->getNumContigs();
		bamHeader->n_ref() = numContigs;
		BAMHeaderRefSeq* refseq = bamHeader->firstRefSeq();
		GenomeDistance genomeLen = context.genome->getCountOfBases();
		for (int i = 0; i < numContigs; i++) {
			int len = (int)strlen(contigs[i].name) + 1;
			cursor += BAMHeaderRefSeq::size(len);
			if (cursor > headerBufferSize) {
				return false;
			}
			refseq->l_name = len;
			memcpy(refseq->name(), contigs[i].name, len);
			GenomeLocation start = contigs[i].beginningLocation;
            GenomeLocation end = ((i + 1 < numContigs) ? contigs[i+1].beginningLocation : genomeLen) - context.genome->getChromosomePadding();
            refseq->l_ref() = (int)(end - start);
			refseq = refseq->next();
			_ASSERT((char*) refseq - header == cursor);
		}
	}
    *headerActualSize = cursor;
    return true;
}

    bool 
BAMFormat::writePairs(
    const ReaderContext& context,
    LandauVishkinWithCigar * lv,
    AffineGapVectorizedWithCigar * ag,
    bool useAffineGap,
    char * buffer,
    size_t bufferSpace,
    size_t * spaceUsed,
    size_t * qnameLen,
    Read ** reads,
    GenomeLocation * locations,
    PairedAlignmentResult * result,
    bool isSecondary,
    bool emitInternalScore,
    char *internalScoreTag,
    int * writeOrder,
    int * cumulativePositiveAddFrontClipping,
    bool * secondReadLocationChanged,
    bool * outOfSpace) const
{
    const int MAX_READ = MAX_READ_LENGTH;
    const int cigarBufSize = MAX_READ;
    _uint32 cigarBuf[2][cigarBufSize];
    int cigarOps[2] = {0, 0};

    int flags[2] = {0, 0};
    const char *contigName[2] = {"*", "*"};
    int contigIndex[2] = {-1, -1};
    GenomeDistance positionInContig[2] = {0, 0};
    const char *mateContigName[2] = {"*", "*"};
    int mateContigIndex[2] = {-1, -1};
    GenomeDistance matePositionInContig[2] = {0, 0};
    _int64 templateLength[2] = {0, 0};

    char data[2][MAX_READ];
    char quality[2][MAX_READ];

    const char* clippedData[2];
    unsigned fullLength[2];
    unsigned clippedLength[2];
    unsigned basesClippedBefore[2];
    unsigned basesClippedAfter[2];
    GenomeDistance extraBasesClippedBefore[2];   // Clipping added if we align before the beginning of a chromosome
    int editDistance[2] = {-1, -1};
    int refSpanFromCigar[2] = {0, 0};

    // Create SAM entry and compute CIGAR
    for (int firstOrSecond = 0; firstOrSecond < NUM_READS_PER_PAIR; firstOrSecond++) {
        int whichRead = writeOrder[firstOrSecond];
        Read* read = reads[whichRead];
        bool firstInPair = writeOrder[firstOrSecond] == 0;

        int addFrontClipping;
        do {
            addFrontClipping = 0;
            if (!SAMFormat::createSAMLine(context.genome, data[whichRead], quality[whichRead], MAX_READ, contigName[whichRead], contigIndex[whichRead],
                    flags[whichRead], positionInContig[whichRead], result->mapq[whichRead], contigName[1 - whichRead], contigIndex[1 - whichRead],
                    positionInContig[1 - whichRead], templateLength[whichRead],
                    fullLength[whichRead], clippedData[whichRead], clippedLength[whichRead], basesClippedBefore[whichRead], basesClippedAfter[whichRead],
                    basesClippedBefore[1 - whichRead], basesClippedAfter[1 - whichRead], qnameLen[whichRead], reads[whichRead], 
                    result->status[whichRead], locations[whichRead], result->direction[whichRead], isSecondary, result->supplementary[whichRead], useM,
                    true, firstInPair, result->alignedAsPair, reads[1 - whichRead], result->status[1 - whichRead], locations[1 - whichRead], result->direction[1 - whichRead], 
                    &extraBasesClippedBefore[whichRead], result->basesClippedBefore[whichRead], result->basesClippedAfter[whichRead], 
                    result->basesClippedBefore[1 - whichRead], result->basesClippedAfter[1 - whichRead]))
            {
                return false;
            }

            if (locations[whichRead] != InvalidGenomeLocation) {
                if (useAffineGap && (result->usedAffineGapScoring[whichRead] || result->score[whichRead] > 0)) {
                    cigarOps[whichRead] = computeCigarOps(context.genome, ag, (char*)cigarBuf[whichRead], cigarBufSize * sizeof(_uint32),
                        clippedData[whichRead], clippedLength[whichRead], basesClippedBefore[whichRead], (unsigned)extraBasesClippedBefore[whichRead], basesClippedAfter[whichRead],
                        read->getOriginalFrontHardClipping(), read->getOriginalBackHardClipping(),
                        locations[whichRead], result->direction[whichRead] == RC, useM, &editDistance[whichRead], &addFrontClipping, &refSpanFromCigar[whichRead]);
                    if (addFrontClipping != 0) {
                        *secondReadLocationChanged = firstOrSecond == 1;
                        const Genome::Contig *originalContig = context.genome->getContigAtLocation(locations[whichRead]);
                        const Genome::Contig *newContig = context.genome->getContigAtLocation(locations[whichRead] + addFrontClipping);
                        if (newContig != originalContig || NULL == newContig || locations[whichRead] + addFrontClipping > originalContig->beginningLocation + originalContig->length - context.genome->getChromosomePadding()) {
                            //
                            // Altering this would push us over a contig boundary.  Just give up on the read.
                            //
                            result->status[whichRead] = NotFound;
                            result->location[whichRead] = InvalidGenomeLocation;
                            locations[whichRead] = InvalidGenomeLocation;
                        }
                        else {
                            if (addFrontClipping < 0) { // Insertion (soft-clip)
                                cumulativePositiveAddFrontClipping[firstOrSecond] += addFrontClipping;
                                if (result->direction[whichRead] == FORWARD) {
                                    reads[whichRead]->setAdditionalFrontClipping(-cumulativePositiveAddFrontClipping[firstOrSecond]);
                                }
                                else {
                                    reads[whichRead]->setAdditionalBackClipping(-cumulativePositiveAddFrontClipping[firstOrSecond]);
                                }
                            }
                            else { // Deletion
                                locations[whichRead] += addFrontClipping;
                            }
                        }
                    }
                }
                else {
                    cigarOps[whichRead] = computeCigarOps(context.genome, lv, (char*)cigarBuf[whichRead], cigarBufSize * sizeof(_uint32),
                        clippedData[whichRead], clippedLength[whichRead], basesClippedBefore[whichRead], (unsigned)extraBasesClippedBefore[whichRead], basesClippedAfter[whichRead],
                        read->getOriginalFrontHardClipping(), read->getOriginalBackHardClipping(),
                        locations[whichRead], result->direction[whichRead] == RC, useM, &editDistance[whichRead], &addFrontClipping, &refSpanFromCigar[whichRead]);

                    if (addFrontClipping != 0) {
                        *secondReadLocationChanged = firstOrSecond == 1;
                        const Genome::Contig *originalContig = context.genome->getContigAtLocation(locations[whichRead]);
                        const Genome::Contig *newContig = context.genome->getContigAtLocation(locations[whichRead] + addFrontClipping);
                        if (newContig != originalContig || NULL == newContig || locations[whichRead] + addFrontClipping > originalContig->beginningLocation + originalContig->length - context.genome->getChromosomePadding()) {
                            //
                            // Altering this would push us over a contig boundary.  Just give up on the read.
                            //
                            result->status[whichRead] = NotFound;
                            result->location[whichRead] = InvalidGenomeLocation;
                            locations[whichRead] = InvalidGenomeLocation;
                        }
                        else {
                            if (addFrontClipping > 0) {
                                cumulativePositiveAddFrontClipping[firstOrSecond] += addFrontClipping;
                                reads[whichRead]->setAdditionalFrontClipping(cumulativePositiveAddFrontClipping[firstOrSecond]);
                            }
                            locations[whichRead] += addFrontClipping;
                        }
                    }
                }

                // Uncomment for debug
                if (editDistance[whichRead] == -1) {
                    const char* read_data = read->getUnclippedData();
                    const char* readId = read->getId();
                    for (unsigned i = 0; i < read->getIdLength(); ++i) {
                        printf("%c", readId[i]);
                    }
                    printf(",");
                    for (unsigned i = 0; i < read->getUnclippedLength(); ++i) {
                        printf("%c", read_data[i]);
                    }
                    printf("\n");
                }
            }
		} while (addFrontClipping != 0);
	}

    // Fill mate information
    for (int firstOrSecond = 0; firstOrSecond < NUM_READS_PER_PAIR; firstOrSecond++) {
        int whichRead = writeOrder[firstOrSecond];
        bool firstInPair = writeOrder[firstOrSecond] == 0;
        for (unsigned i = 0; i < fullLength[whichRead]; i++) {
            quality[whichRead][i] -= '!';
        }
        SAMFormat::fillMateInfo(context.genome, flags[whichRead], reads[whichRead], locations[whichRead], result->direction[whichRead], 
            contigName[whichRead], contigIndex[whichRead], positionInContig[whichRead], templateLength[whichRead], basesClippedBefore[whichRead],
            firstInPair, result->alignedAsPair, reads[1 - whichRead], locations[1 - whichRead], result->direction[1 - whichRead],
            mateContigName[whichRead], mateContigIndex[whichRead], matePositionInContig[whichRead], basesClippedBefore[1 - whichRead],
            refSpanFromCigar[whichRead], refSpanFromCigar[1 - whichRead]);
    }
    
    // Write the BAM entry
    for (int firstOrSecond = 0; firstOrSecond < NUM_READS_PER_PAIR; firstOrSecond++) {
        int whichRead = writeOrder[firstOrSecond];
        Read* read = reads[whichRead];
        unsigned auxLen;
        bool auxSAM;
        char* aux = read->getAuxiliaryData(&auxLen, &auxSAM);
        static bool warningPrinted = false;
        bool translateReadGroupFromSAM = false;
        if (aux != NULL && auxSAM) {
            if (!warningPrinted) {
                warningPrinted = true;
                WriteErrorMessage("warning: translating optional data from SAM->BAM is not yet implemented, optional data will not appear in BAM\n");
            }
            if (read->getReadGroup() == READ_GROUP_FROM_AUX) {
                for (char* p = aux; p != NULL && p < aux + auxLen; p = SAMReader::skipToBeyondNextFieldSeparator(p, aux + auxLen)) {
                    if (strncmp(p, "RG:Z:", 5) == 0) {
                        size_t fieldLen;
                        SAMReader::skipToBeyondNextFieldSeparator(p, aux + auxLen, &fieldLen);
                        aux = p;
                        auxLen = (unsigned)fieldLen;
                        translateReadGroupFromSAM = true;
                        break;
                    }
                }
            }
            if (!translateReadGroupFromSAM) {
                aux = NULL;
                auxLen = 0;
            }
        }
        size_t bamSize = BAMAlignment::size((unsigned)qnameLen[whichRead] + 1, cigarOps[whichRead], fullLength[whichRead], !translateReadGroupFromSAM ? auxLen : auxLen - 1);
        if (read->getReadGroup() != NULL && read->getReadGroup() != READ_GROUP_FROM_AUX) {
            if (strcmp(read->getReadGroup(), context.defaultReadGroup) != 0) {
                bamSize += 4 + strlen(read->getReadGroup());
            }
            else {
                bamSize += context.defaultReadGroupAuxLen;
            }
        }

        if (read->getLibrary() != NULL) {
            bamSize += 4 + read->getLibraryLength();
        }

        //
        // Add in the size for the tags now.  We have to do this in this ugly way, because the tags are written directly into the
        // buffer, so we can only write them if there's space.  However, BamAuxAlign::size() depends on the contents of the aux field
        // (obviously), so we can't call it until it's filled in.  Which, of course, we can't do until the space is allocated.  Hence,
        // this plus some asserts below.
        //
        bamSize += 8 + 4 + (emitInternalScore ? 7 : 0); // NM:C PG:Z:SNAP fields and optionally the internal score field (which is 32 bits rather than the 8 used in NM)
        bamSize += 7; // extra space to store mate quality score for duplicate marking
        if (bamSize > bufferSpace) {
            *outOfSpace = true;
            return false;
        }
        BAMAlignment* bam = (BAMAlignment*)buffer;
        bam->block_size = (int)bamSize - 4;
        bam->refID = contigIndex[whichRead];
        if (positionInContig[whichRead] > INT32_MAX || matePositionInContig[whichRead] > INT32_MAX) {
            WriteErrorMessage("Can't write read to BAM file because aligned position (or mate position) within contig > 2^31, which is the limit for the BAM format.\n");
            soft_exit(1);
        }
        bam->pos = (int)(positionInContig[whichRead] - 1);

        if (qnameLen[whichRead] > 254) {
            WriteErrorMessage("BAM format: QNAME field must be less than 254 characters long, instead it's %lld\n", qnameLen[whichRead]);
            soft_exit(1);
        }
        bam->l_read_name = (_uint8)qnameLen[whichRead] + 1;
        bam->MAPQ = result->mapq[whichRead];
        int refLength = cigarOps[whichRead] > 0 ? 0 : fullLength[whichRead];
        for (int i = 0; i < cigarOps[whichRead]; i++) {
            refLength += BAMAlignment::CigarCodeToRefBase[cigarBuf[whichRead][i] & 0xf] * (cigarBuf[whichRead][i] >> 4);
        }
        bam->bin = locations[whichRead] != InvalidGenomeLocation ? BAMAlignment::reg2bin((int)positionInContig[whichRead] - 1, (int)positionInContig[whichRead] - 1 + refLength) :
            // unmapped is at mate's position, length 1
            locations[1 - whichRead] != InvalidGenomeLocation ? BAMAlignment::reg2bin((int)matePositionInContig[whichRead] - 1, (int)matePositionInContig[whichRead]) :
            // otherwise at -1, length 1
            BAMAlignment::reg2bin(-1, 0);
        bam->n_cigar_op = cigarOps[whichRead];
        bam->FLAG = flags[whichRead];
        bam->l_seq = fullLength[whichRead];
        bam->next_refID = mateContigIndex[whichRead];
        bam->next_pos = (int)matePositionInContig[whichRead] - 1;
        bam->tlen = (int)templateLength[whichRead];
        memcpy(bam->read_name(), read->getId(), qnameLen[whichRead]);
        bam->read_name()[qnameLen[whichRead]] = 0;
        memcpy(bam->cigar(), cigarBuf[whichRead], cigarOps[whichRead] * 4);
        BAMAlignment::encodeSeq(bam->seq(), data[whichRead], fullLength[whichRead]);

        memcpy(bam->qual(), quality[whichRead], fullLength[whichRead]);
        if (aux != NULL && auxLen > 0) {
            if (((char*)bam->firstAux()) + auxLen > buffer + bufferSpace) {
                *outOfSpace = true;
                return false;
            }
            if (!translateReadGroupFromSAM) {
                memcpy(bam->firstAux(), aux, auxLen);
            }
            else {
                // hack, build just RG field from SAM opt field
                BAMAlignAux* auxData = bam->firstAux();
                auxData->tag[0] = 'R';
                auxData->tag[1] = 'G';
                auxData->val_type = 'Z';
                memcpy(auxData->value(), aux + 5, auxLen - 5);
                ((char*)auxData->value())[auxLen - 5] = 0;
                auxLen -= 1; // RG:Z:xxx -> RGZxxx\0
            }
        }
        // RG
        if (read->getReadGroup() != NULL && read->getReadGroup() != READ_GROUP_FROM_AUX) {
            if (strcmp(read->getReadGroup(), context.defaultReadGroup) != 0) {
                if ((char*)bam->firstAux() + auxLen + 4 + strlen(read->getReadGroup()) > buffer + bufferSpace) {
                    *outOfSpace = true;
                    return false;
                }
                BAMAlignAux* rg = (BAMAlignAux*)(auxLen + (char*)bam->firstAux());
                rg->tag[0] = 'R'; rg->tag[1] = 'G'; rg->val_type = 'Z';
                strcpy((char*)rg->value(), read->getReadGroup());
                auxLen += (unsigned)rg->size();
            }
            else {
                if ((char*)bam->firstAux() + auxLen + context.defaultReadGroupAuxLen > buffer + bufferSpace) {
                    *outOfSpace = true;
                    return false;
                }
                memcpy((char*)bam->firstAux() + auxLen, context.defaultReadGroupAux, context.defaultReadGroupAuxLen);
                auxLen += context.defaultReadGroupAuxLen;
            }
        }

        // PG
        BAMAlignAux* pg = (BAMAlignAux*)(auxLen + (char*)bam->firstAux());
        pg->tag[0] = 'P'; pg->tag[1] = 'G'; pg->val_type = 'Z';
        strcpy((char*)pg->value(), "SNAP");
        _ASSERT(pg->size() == 8);   // Known above in the bamSize += line
        auxLen += (unsigned)pg->size();
        // NM
        BAMAlignAux* nm = (BAMAlignAux*)(auxLen + (char*)bam->firstAux());
        nm->tag[0] = 'N'; nm->tag[1] = 'M'; nm->val_type = 'C';
        *(_uint8*)nm->value() = (_uint8)editDistance[whichRead];
        _ASSERT(nm->size() == 4);   // Known above in the bamSize += line
        auxLen += (unsigned)nm->size();

        if (emitInternalScore) {
            BAMAlignAux *in = (BAMAlignAux*)(auxLen + (char *)bam->firstAux());
            in->tag[0] = internalScoreTag[0];  in->tag[1] = internalScoreTag[1]; in->val_type = 'i';
            *(_int32*)in->value() = (flags[whichRead] & SAM_UNMAPPED) ? -1 : result->scorePriorToClipping[whichRead];
            _ASSERT(in->size() == 7);   // Known above in the bamSize += line
            auxLen += (unsigned)in->size();
        }

        // QS
        int result = 0;
        _uint8* p = (_uint8*)quality[1 - whichRead];
        for (int i = 0; i < fullLength[1 - whichRead]; i++) {
            int q = *p++;
            // Picard MarkDup uses a score threshold of 15 (default)
            result += (q >= 15) ? (q != 255) * q : 0; // avoid branch?
        }
        BAMAlignAux* mq = (BAMAlignAux*)(auxLen + (char*)bam->firstAux());
        mq->tag[0] = 'Q'; mq->tag[1] = 'S'; mq->val_type = 'i';
        *(_int32*)mq->value() = result;
        _ASSERT(mq->size() == 7);   // Known above in the bamSize += line
        auxLen += (unsigned)mq->size();

        // LB
        if (read->getLibrary() != NULL) {
            if ((char*)bam->firstAux() + auxLen + 4 + read->getLibraryLength() > buffer + bufferSpace) {
                *outOfSpace = true;
                return false;
            }
            BAMAlignAux* lb = (BAMAlignAux*)(auxLen + (char*)bam->firstAux());
            lb->tag[0] = 'L'; lb->tag[1] = 'B'; lb->val_type = 'Z';
            strncpy((char*)lb->value(), read->getLibrary(), read->getLibraryLength());
            ((char*)(lb->value()))[read->getLibraryLength()] = '\0';
            _ASSERT(lb->size() == 4 + read->getLibraryLength());
            auxLen += (unsigned)lb->size();
        }

        if (NULL != spaceUsed) {
            spaceUsed[firstOrSecond] = bamSize;
        }

        buffer += spaceUsed[firstOrSecond];
        bufferSpace -= spaceUsed[firstOrSecond];

        // debugging: _ASSERT(0 == memcmp(bam->firstAux()->tag, "RG", 2) && 0 == memcmp(bam->firstAux()->next()->tag, "PG", 2) && 0 == memcmp(bam->firstAux()->next()->next()->tag, "NM", 2));
        bam->validate();
    }
    return true;
}

    bool
BAMFormat::writeRead(
    const ReaderContext& context,
    LandauVishkinWithCigar * lv,
    char * buffer,
    size_t bufferSpace,
    size_t * spaceUsed,
    size_t qnameLen,
    Read * read,
    AlignmentResult result,
    int mapQuality,
    GenomeLocation genomeLocation,
    Direction direction,
    bool secondaryAlignment,
    bool supplementaryAlignment,
    int *o_addFrontClipping,
    int internalScore, 
    bool emitInternalScore, 
    char *internalScoreTag,
    int bpClippedBefore,
    int bpClippedAfter, 
    bool hasMate,
    bool firstInPair,
    Read * mate,
    AlignmentResult mateResult,
    GenomeLocation mateLocation,
    Direction mateDirection,
    bool alignedAsPair,
    int mateBpClippedBefore,
    int mateBpClippedAfter) const
{
    const int MAX_READ = MAX_READ_LENGTH;
    const int cigarBufSize = MAX_READ;
    _uint32 cigarBuf[cigarBufSize];

    int flags = 0;
    const char *contigName = "*";
    int contigIndex = -1;
    GenomeDistance positionInContig = 0;
    int cigarOps = 0;
    const char *mateContigName = "*";
    int mateContigIndex = -1;
    GenomeDistance matePositionInContig = 0;
    _int64 templateLength = 0;

    char data[MAX_READ];
    char quality[MAX_READ];

    const char* clippedData;
    unsigned fullLength;
    unsigned clippedLength;
    unsigned basesClippedBefore, mateBasesClippedBefore;
    GenomeDistance extraBasesClippedBefore;
    unsigned basesClippedAfter, mateBasesClippedAfter;
    int editDistance = -1;
    int newAddFrontClipping = 0;
    int refSpanFromCigar = 0;

    if (!SAMFormat::createSAMLine(context.genome, 
        // outputs:
        data, quality, MAX_READ, contigName, contigIndex,
        flags, positionInContig, mapQuality, mateContigName, mateContigIndex, matePositionInContig, templateLength,
        fullLength, clippedData, clippedLength, basesClippedBefore, basesClippedAfter, mateBasesClippedBefore, mateBasesClippedAfter,
        // inputs:
        qnameLen, read, result, genomeLocation, direction, secondaryAlignment, supplementaryAlignment, useM,
        hasMate, firstInPair, alignedAsPair, mate, mateResult, mateLocation, mateDirection,
        &extraBasesClippedBefore, bpClippedBefore, bpClippedAfter, mateBpClippedBefore, mateBpClippedAfter))
    {
        return false;
    }
    if (genomeLocation != InvalidGenomeLocation) {
        cigarOps = computeCigarOps(context.genome, lv, (char*)cigarBuf, cigarBufSize * sizeof(_uint32),
                                   clippedData, clippedLength, basesClippedBefore, (unsigned)extraBasesClippedBefore, basesClippedAfter,
                                   read->getOriginalFrontHardClipping(), read->getOriginalBackHardClipping(),
                                   genomeLocation, direction == RC, useM, &editDistance, o_addFrontClipping, &refSpanFromCigar);
        if (*o_addFrontClipping != 0) {
            return false;
        }
    }

    // Write the BAM entry
    unsigned auxLen;
    bool auxSAM;
    char* aux = read->getAuxiliaryData(&auxLen, &auxSAM);
    static bool warningPrinted = false;
    bool translateReadGroupFromSAM = false;
    if (aux != NULL && auxSAM) {
        if (! warningPrinted) {
            warningPrinted = true;
            WriteErrorMessage("warning: translating optional data from SAM->BAM is not yet implemented, optional data will not appear in BAM\n");
        }
        if (read->getReadGroup() == READ_GROUP_FROM_AUX) {
            for (char* p = aux; p != NULL && p < aux + auxLen; p = SAMReader::skipToBeyondNextFieldSeparator(p, aux + auxLen)) {
                if (strncmp(p, "RG:Z:", 5) == 0) {
                    size_t fieldLen;
                    SAMReader::skipToBeyondNextFieldSeparator(p, aux + auxLen, &fieldLen);
                    aux = p;
                    auxLen = (unsigned) fieldLen;
                    translateReadGroupFromSAM = true;
                    break;
                }
            }
        }
        if (! translateReadGroupFromSAM) {
            aux = NULL;
            auxLen = 0;
        }
    }
    size_t bamSize = BAMAlignment::size((unsigned)qnameLen + 1, cigarOps, fullLength, !translateReadGroupFromSAM ? auxLen : auxLen - 1);
    if (read->getReadGroup() != NULL && read->getReadGroup() != READ_GROUP_FROM_AUX) {
        if (strcmp(read->getReadGroup(), context.defaultReadGroup) != 0) {
            bamSize += 4 + strlen(read->getReadGroup());
        } else {
            bamSize += context.defaultReadGroupAuxLen;
        }
    }

    if (read->getLibrary() != NULL) {
        bamSize += 4 + read->getLibraryLength();
    }

    //
    // Add in the size for the tags now.  We have to do this in this ugly way, because the tags are written directly into the
    // buffer, so we can only write them if there's space.  However, BamAuxAlign::size() depends on the contents of the aux field
    // (obviously), so we can't call it until it's filled in.  Which, of course, we can't do until the space is allocated.  Hence,
    // this plus some asserts below.
    //
    bamSize += 8 + 4 + (emitInternalScore ? 7 : 0); // NM:C PG:Z:SNAP fields and optionally the internal score field (which is 32 bits rather than the 8 used in NM)
    if (bamSize > bufferSpace) {
        return false;
    }
    BAMAlignment* bam = (BAMAlignment*) buffer;
    bam->block_size = (int)bamSize - 4;
    bam->refID = contigIndex;
    if (positionInContig > INT32_MAX || matePositionInContig > INT32_MAX) {
        WriteErrorMessage("Can't write read to BAM file because aligned position (or mate position) within contig > 2^31, which is the limit for the BAM format.\n");
        soft_exit(1);
    }
    bam->pos = (int)(positionInContig - 1);

    if (qnameLen > 254) {
        WriteErrorMessage("BAM format: QNAME field must be less than 254 characters long, instead it's %lld\n", qnameLen);
        soft_exit(1);
    }
    bam->l_read_name = (_uint8)qnameLen + 1;
    bam->MAPQ = mapQuality;
    int refLength = cigarOps > 0 ? 0 : fullLength;
    for (int i = 0; i < cigarOps; i++) {
        refLength += BAMAlignment::CigarCodeToRefBase[cigarBuf[i] & 0xf] * (cigarBuf[i] >> 4);
    }
    bam->bin = genomeLocation != InvalidGenomeLocation ? BAMAlignment::reg2bin((int)positionInContig-1, (int)positionInContig-1 + refLength) :
		// unmapped is at mate's position, length 1
		mateLocation != InvalidGenomeLocation ? BAMAlignment::reg2bin((int)matePositionInContig-1, (int)matePositionInContig) :
		// otherwise at -1, length 1
		BAMAlignment::reg2bin(-1, 0);
    bam->n_cigar_op = cigarOps;
    bam->FLAG = flags;
    bam->l_seq = fullLength;
    bam->next_refID = mateContigIndex;
    bam->next_pos = (int)matePositionInContig - 1;
    bam->tlen = (int)templateLength;
    memcpy(bam->read_name(), read->getId(), qnameLen);
    bam->read_name()[qnameLen] = 0;
    memcpy(bam->cigar(), cigarBuf, cigarOps * 4);
    BAMAlignment::encodeSeq(bam->seq(), data, fullLength);
    for (unsigned i = 0; i < fullLength; i++) {
        quality[i] -= '!';
    }
    memcpy(bam->qual(), quality, fullLength);
    if (aux != NULL && auxLen > 0) {
        if (((char*)bam->firstAux()) + auxLen > buffer + bufferSpace) {
            return false;
        }
        if (! translateReadGroupFromSAM) {
            memcpy(bam->firstAux(), aux, auxLen);
        } else {
            // hack, build just RG field from SAM opt field
            BAMAlignAux* auxData = bam->firstAux();
            auxData->tag[0] = 'R';
            auxData->tag[1] = 'G';
            auxData->val_type = 'Z';
            memcpy(auxData->value(), aux + 5, auxLen - 5);
            ((char*)auxData->value())[auxLen-5] = 0;
            auxLen -= 1; // RG:Z:xxx -> RGZxxx\0
        }
    }
    // RG
    if (read->getReadGroup() != NULL && read->getReadGroup() != READ_GROUP_FROM_AUX) {
        if (strcmp(read->getReadGroup(), context.defaultReadGroup) != 0) {
            if ((char*)bam->firstAux() + auxLen + 4 + strlen(read->getReadGroup()) > buffer + bufferSpace) {
                return false;
            }
            BAMAlignAux* rg = (BAMAlignAux*)(auxLen + (char*)bam->firstAux());
            rg->tag[0] = 'R'; rg->tag[1] = 'G'; rg->val_type = 'Z';
            strcpy((char*)rg->value(), read->getReadGroup());
            auxLen += (unsigned)rg->size();
        } else {
            if ((char*)bam->firstAux() + auxLen + context.defaultReadGroupAuxLen > buffer + bufferSpace) {
                return false;
            }
            memcpy((char*)bam->firstAux() + auxLen, context.defaultReadGroupAux, context.defaultReadGroupAuxLen);
            auxLen += context.defaultReadGroupAuxLen;
        }
    }
    // PG
    BAMAlignAux* pg = (BAMAlignAux*) (auxLen + (char*) bam->firstAux());
    pg->tag[0] = 'P'; pg->tag[1] = 'G'; pg->val_type = 'Z';
    strcpy((char*) pg->value(), "SNAP");
    _ASSERT(pg->size() == 8);   // Known above in the bamSize += line
    auxLen += (unsigned) pg->size();
    // NM
    BAMAlignAux* nm = (BAMAlignAux*) (auxLen + (char*) bam->firstAux());
    nm->tag[0] = 'N'; nm->tag[1] = 'M'; nm->val_type = 'C';
    *(_uint8*)nm->value() = (_uint8)editDistance;
    _ASSERT(nm->size() == 4);   // Known above in the bamSize += line
    auxLen += (unsigned)nm->size();

    if (emitInternalScore) {
        BAMAlignAux *in = (BAMAlignAux*)(auxLen + (char *)bam->firstAux());
        in->tag[0] = internalScoreTag[0];  in->tag[1] = internalScoreTag[1]; in->val_type = 'i';
        *(_int32*)in->value() = (flags & SAM_UNMAPPED) ? -1 : internalScore;
        _ASSERT(in->size() == 7);   // Known above in the bamSize += line
        auxLen += (unsigned)in->size();
    }

    // LB
    if (read->getLibrary() != NULL) {
        if ((char*)bam->firstAux() + auxLen + 4 + read->getLibraryLength() > buffer + bufferSpace) {
            return false;
        }
        BAMAlignAux* lb = (BAMAlignAux*)(auxLen + (char*)bam->firstAux());
        lb->tag[0] = 'L'; lb->tag[1] = 'B'; lb->val_type = 'Z';
        strncpy((char*)lb->value(), read->getLibrary(), read->getLibraryLength());
        ((char*)(lb->value()))[read->getLibraryLength()] = '\0';
        _ASSERT(lb->size() == 4 + read->getLibraryLength());
        auxLen += (unsigned)lb->size();
    }

    if (NULL != spaceUsed) {
        *spaceUsed = bamSize;
    }
    // debugging: _ASSERT(0 == memcmp(bam->firstAux()->tag, "RG", 2) && 0 == memcmp(bam->firstAux()->next()->tag, "PG", 2) && 0 == memcmp(bam->firstAux()->next()->next()->tag, "NM", 2));
    bam->validate();
    return true;
}

    bool
BAMFormat::writeRead(
    const ReaderContext& context,
    AffineGapVectorizedWithCigar * ag,
    char * buffer,
    size_t bufferSpace,
    size_t * spaceUsed,
    size_t qnameLen,
    Read * read,
    AlignmentResult result,
    int mapQuality,
    GenomeLocation genomeLocation,
    Direction direction,
    bool secondaryAlignment,
    bool supplementaryAlignment,
    int *o_addFrontClipping,
    int internalScore,
    bool emitInternalScore,
    char *internalScoreTag,
    int bpClippedBefore,
    int bpClippedAfter,
    bool hasMate,
    bool firstInPair,
    Read * mate,
    AlignmentResult mateResult,
    GenomeLocation mateLocation,
    Direction mateDirection,
    bool alignedAsPair,
    int mateBpClippedBefore,
    int mateBpClippedAfter) const
{
    const int MAX_READ = MAX_READ_LENGTH;
    const int cigarBufSize = MAX_READ;
    _uint32 cigarBuf[cigarBufSize];

    int flags = 0;
    const char *contigName = "*";
    int contigIndex = -1;
    GenomeDistance positionInContig = 0;
    int cigarOps = 0;
    const char *mateContigName = "*";
    int mateContigIndex = -1;
    GenomeDistance matePositionInContig = 0;
    _int64 templateLength = 0;

    char data[MAX_READ];
    char quality[MAX_READ];

    const char* clippedData;
    unsigned fullLength;
    unsigned clippedLength;
    unsigned basesClippedBefore, mateBasesClippedBefore;
    GenomeDistance extraBasesClippedBefore;
    unsigned basesClippedAfter, mateBasesClippedAfter;
    int editDistance = -1;
    int newAddFrontClipping = 0;
    int refSpanFromCigar = 0;

    if (!SAMFormat::createSAMLine(context.genome,
        // outputs:
        data, quality, MAX_READ, contigName, contigIndex,
        flags, positionInContig, mapQuality, mateContigName, mateContigIndex, matePositionInContig, templateLength,
        fullLength, clippedData, clippedLength, basesClippedBefore, basesClippedAfter, mateBasesClippedBefore, mateBasesClippedAfter,
        // inputs:
        qnameLen, read, result, genomeLocation, direction, secondaryAlignment, supplementaryAlignment, useM,
        hasMate, firstInPair, alignedAsPair, mate, mateResult, mateLocation, mateDirection,
        &extraBasesClippedBefore, bpClippedBefore, bpClippedAfter, mateBpClippedBefore, mateBpClippedAfter))
    {
        return false;
    }
    if (genomeLocation != InvalidGenomeLocation) {
        cigarOps = computeCigarOps(context.genome, ag, (char*)cigarBuf, cigarBufSize * sizeof(_uint32),
            clippedData, clippedLength, basesClippedBefore, (unsigned)extraBasesClippedBefore, basesClippedAfter,
            read->getOriginalFrontHardClipping(), read->getOriginalBackHardClipping(),
            genomeLocation, direction == RC, useM, &editDistance, o_addFrontClipping, &refSpanFromCigar);
        // Uncomment for debug
        if (editDistance == -1) {
            const char* read_data = read->getUnclippedData();
            const char* readId = read->getId();
            for (unsigned i = 0; i < read->getIdLength(); ++i) {
                printf("%c", readId[i]);
            }
            printf(",");
            for (unsigned i = 0; i < read->getUnclippedLength(); ++i) {
                printf("%c", read_data[i]);
            }
            printf("\n");
        }
        if (*o_addFrontClipping != 0) {
            return false;
        }
    }

    // Write the BAM entry
    unsigned auxLen;
    bool auxSAM;
    char* aux = read->getAuxiliaryData(&auxLen, &auxSAM);
    static bool warningPrinted = false;
    bool translateReadGroupFromSAM = false;
    if (aux != NULL && auxSAM) {
        if (!warningPrinted) {
            warningPrinted = true;
            WriteErrorMessage("warning: translating optional data from SAM->BAM is not yet implemented, optional data will not appear in BAM\n");
        }
        if (read->getReadGroup() == READ_GROUP_FROM_AUX) {
            for (char* p = aux; p != NULL && p < aux + auxLen; p = SAMReader::skipToBeyondNextFieldSeparator(p, aux + auxLen)) {
                if (strncmp(p, "RG:Z:", 5) == 0) {
                    size_t fieldLen;
                    SAMReader::skipToBeyondNextFieldSeparator(p, aux + auxLen, &fieldLen);
                    aux = p;
                    auxLen = (unsigned)fieldLen;
                    translateReadGroupFromSAM = true;
                    break;
                }
            }
        }
        if (!translateReadGroupFromSAM) {
            aux = NULL;
            auxLen = 0;
        }
    }
    size_t bamSize = BAMAlignment::size((unsigned)qnameLen + 1, cigarOps, fullLength, !translateReadGroupFromSAM ? auxLen : auxLen - 1);
    if (read->getReadGroup() != NULL && read->getReadGroup() != READ_GROUP_FROM_AUX) {
        if (strcmp(read->getReadGroup(), context.defaultReadGroup) != 0) {
            bamSize += 4 + strlen(read->getReadGroup());
        }
        else {
            bamSize += context.defaultReadGroupAuxLen;
        }
    }

    if (read->getLibrary() != NULL) {
        bamSize += 4 + read->getLibraryLength();
    }

    //
    // Add in the size for the tags now.  We have to do this in this ugly way, because the tags are written directly into the
    // buffer, so we can only write them if there's space.  However, BamAuxAlign::size() depends on the contents of the aux field
    // (obviously), so we can't call it until it's filled in.  Which, of course, we can't do until the space is allocated.  Hence,
    // this plus some asserts below.
    //
    bamSize += 8 + 4 + (emitInternalScore ? 7 : 0); // NM:C PG:Z:SNAP fields and optionally the internal score field (which is 32 bits rather than the 8 used in NM)
    if (bamSize > bufferSpace) {
        return false;
    }
    BAMAlignment* bam = (BAMAlignment*)buffer;
    bam->block_size = (int)bamSize - 4;
    bam->refID = contigIndex;
    if (positionInContig > INT32_MAX || matePositionInContig > INT32_MAX) {
        WriteErrorMessage("Can't write read to BAM file because aligned position (or mate position) within contig > 2^31, which is the limit for the BAM format.\n");
        soft_exit(1);
    }
    bam->pos = (int)(positionInContig - 1);

    if (qnameLen > 254) {
        WriteErrorMessage("BAM format: QNAME field must be less than 254 characters long, instead it's %lld\n", qnameLen);
        soft_exit(1);
    }
    bam->l_read_name = (_uint8)qnameLen + 1;
    bam->MAPQ = mapQuality;
    int refLength = cigarOps > 0 ? 0 : fullLength;
    for (int i = 0; i < cigarOps; i++) {
        refLength += BAMAlignment::CigarCodeToRefBase[cigarBuf[i] & 0xf] * (cigarBuf[i] >> 4);
    }
    bam->bin = genomeLocation != InvalidGenomeLocation ? BAMAlignment::reg2bin((int)positionInContig - 1, (int)positionInContig - 1 + refLength) :
        // unmapped is at mate's position, length 1
        mateLocation != InvalidGenomeLocation ? BAMAlignment::reg2bin((int)matePositionInContig - 1, (int)matePositionInContig) :
        // otherwise at -1, length 1
        BAMAlignment::reg2bin(-1, 0);
    bam->n_cigar_op = cigarOps;
    bam->FLAG = flags;
    bam->l_seq = fullLength;
    bam->next_refID = mateContigIndex;
    bam->next_pos = (int)matePositionInContig - 1;
    bam->tlen = (int)templateLength;
    memcpy(bam->read_name(), read->getId(), qnameLen);
    bam->read_name()[qnameLen] = 0;
    memcpy(bam->cigar(), cigarBuf, cigarOps * 4);
    BAMAlignment::encodeSeq(bam->seq(), data, fullLength);
    for (unsigned i = 0; i < fullLength; i++) {
        quality[i] -= '!';
    }
    memcpy(bam->qual(), quality, fullLength);
    if (aux != NULL && auxLen > 0) {
        if (((char*)bam->firstAux()) + auxLen > buffer + bufferSpace) {
            return false;
        }
        if (!translateReadGroupFromSAM) {
            memcpy(bam->firstAux(), aux, auxLen);
        }
        else {
            // hack, build just RG field from SAM opt field
            BAMAlignAux* auxData = bam->firstAux();
            auxData->tag[0] = 'R';
            auxData->tag[1] = 'G';
            auxData->val_type = 'Z';
            memcpy(auxData->value(), aux + 5, auxLen - 5);
            ((char*)auxData->value())[auxLen - 5] = 0;
            auxLen -= 1; // RG:Z:xxx -> RGZxxx\0
        }
    }
    // RG
    if (read->getReadGroup() != NULL && read->getReadGroup() != READ_GROUP_FROM_AUX) {
        if (strcmp(read->getReadGroup(), context.defaultReadGroup) != 0) {
            if ((char*)bam->firstAux() + auxLen + 4 + strlen(read->getReadGroup()) > buffer + bufferSpace) {
                return false;
            }
            BAMAlignAux* rg = (BAMAlignAux*)(auxLen + (char*)bam->firstAux());
            rg->tag[0] = 'R'; rg->tag[1] = 'G'; rg->val_type = 'Z';
            strcpy((char*)rg->value(), read->getReadGroup());
            auxLen += (unsigned)rg->size();
        }
        else {
            if ((char*)bam->firstAux() + auxLen + context.defaultReadGroupAuxLen > buffer + bufferSpace) {
                return false;
            }
            memcpy((char*)bam->firstAux() + auxLen, context.defaultReadGroupAux, context.defaultReadGroupAuxLen);
            auxLen += context.defaultReadGroupAuxLen;
        }
    }
    // PG
    BAMAlignAux* pg = (BAMAlignAux*)(auxLen + (char*)bam->firstAux());
    pg->tag[0] = 'P'; pg->tag[1] = 'G'; pg->val_type = 'Z';
    strcpy((char*)pg->value(), "SNAP");
    _ASSERT(pg->size() == 8);   // Known above in the bamSize += line
    auxLen += (unsigned)pg->size();
    // NM
    BAMAlignAux* nm = (BAMAlignAux*)(auxLen + (char*)bam->firstAux());
    nm->tag[0] = 'N'; nm->tag[1] = 'M'; nm->val_type = 'C';
    *(_uint8*)nm->value() = (_uint8)editDistance;
    _ASSERT(nm->size() == 4);   // Known above in the bamSize += line
    auxLen += (unsigned)nm->size();

    if (emitInternalScore) {
        BAMAlignAux *in = (BAMAlignAux*)(auxLen + (char *)bam->firstAux());
        in->tag[0] = internalScoreTag[0];  in->tag[1] = internalScoreTag[1]; in->val_type = 'i';
        *(_int32*)in->value() = (flags & SAM_UNMAPPED) ? -1 : internalScore;
        _ASSERT(in->size() == 7);   // Known above in the bamSize += line
        auxLen += (unsigned)in->size();
    }

    // LB
    if (read->getLibrary() != NULL) {
        if ((char*)bam->firstAux() + auxLen + 4 + read->getLibraryLength() > buffer + bufferSpace) {
            return false;
        }
        BAMAlignAux* lb = (BAMAlignAux*)(auxLen + (char*)bam->firstAux());
        lb->tag[0] = 'L'; lb->tag[1] = 'B'; lb->val_type = 'Z';
        strncpy((char*)lb->value(), read->getLibrary(), read->getLibraryLength());
        ((char*)(lb->value()))[read->getLibraryLength()] = '\0';
        _ASSERT(lb->size() == 4 + read->getLibraryLength());
        auxLen += (unsigned)lb->size();
    }

    if (NULL != spaceUsed) {
        *spaceUsed = bamSize;
    }
    // debugging: _ASSERT(0 == memcmp(bam->firstAux()->tag, "RG", 2) && 0 == memcmp(bam->firstAux()->next()->tag, "PG", 2) && 0 == memcmp(bam->firstAux()->next()->next()->tag, "NM", 2));
    bam->validate();
    return true;
}

// Compute the CIGAR edit sequence operations in BAM format for a read against a given genome location
// Returns number of operations (or 0 if there was a problem)
// if returns with *o_addFrontClipping set non-zero, need to adjust front clipping & rerun
    int
BAMFormat::computeCigarOps(
    const Genome *              genome,
    LandauVishkinWithCigar *    lv,
    char *                      cigarBuf,
    int                         cigarBufLen,
    const char *                data,
    unsigned                    dataLength,
    unsigned                    basesClippedBefore,
    unsigned                    extraBasesClippedBefore,
    unsigned                    basesClippedAfter,
    unsigned                    frontHardClipping,
    unsigned                    backHardClipping,
    GenomeLocation              genomeLocation,
    bool                        isRC,
	bool						useM,
    int *                       o_editDistance,
    int *                       o_addFrontClipping,
    int *                       o_refSpan
)
{
    GenomeDistance extraBasesClippedAfter = 0;
    int used = 0;

    unsigned clippingWordsBefore = ((basesClippedBefore + extraBasesClippedBefore > 0) ? 1 : 0) + ((frontHardClipping > 0) ? 1 : 0);
    unsigned clippingWordsAfter = ((basesClippedAfter + extraBasesClippedAfter > 0) ? 1 : 0) + ((backHardClipping > 0) ? 1 : 0);

    SAMFormat::computeCigar(BAM_CIGAR_OPS, genome, lv, cigarBuf + 4 * clippingWordsBefore, cigarBufLen - 4 * (clippingWordsBefore + clippingWordsAfter), data, dataLength, basesClippedBefore, extraBasesClippedBefore,
        basesClippedAfter, &extraBasesClippedAfter, genomeLocation, useM, o_editDistance, &used,  o_addFrontClipping);

    if (*o_addFrontClipping != 0) {
        return 0;
    }

    if (*o_editDistance == -2) {
        WriteErrorMessage("WARNING: computeEditDistance returned -2; cigarBuf may be too small\n");
        return 0;
    } else if (*o_editDistance == -1) {
        static bool warningPrinted = false;
        if (!warningPrinted) {
            WriteErrorMessage( "WARNING: computeEditDistance returned -1; this shouldn't happen. Read %.*s\n", dataLength, data);
            warningPrinted = true;
        }
        return 0;
    } else {
        //
        // If we have hard clipping, add in the cigar string for it.
        //
        if (frontHardClipping > 0) {
            *(_uint32*)cigarBuf = (frontHardClipping << 4) | BAMAlignment::CigarToCode['H'];
            used += 4;
        }
        // Add some CIGAR instructions for soft-clipping if we've ignored some bases in the read.
        if (basesClippedBefore + extraBasesClippedBefore > 0) {
            *((_uint32*)cigarBuf + ((frontHardClipping > 0) ? 1 : 0))  = ((basesClippedBefore + extraBasesClippedBefore) << 4) | BAMAlignment::CigarToCode['S'];
            used += 4;
        }
        if (basesClippedAfter + extraBasesClippedAfter > 0) {
            *(_uint32*)(cigarBuf + used) = ((int)(basesClippedAfter + extraBasesClippedAfter) << 4) | BAMAlignment::CigarToCode['S'];
            used += 4;
        }

        if (backHardClipping > 0) {
            *(_uint32*)(cigarBuf + used) = (backHardClipping << 4) | BAMAlignment::CigarToCode['H'];
            used += 4;
        }

        // refSpan is used for calculating TLEN (e.g., read start 5' - mate end 5') correctly.
        // This is tricky for reads aligning as RC and needs to be inferred from CIGAR
        *o_refSpan = 0;
        const int cigarOps = used / 4;
        if (cigarOps > 0) {
            _uint32 op = *((_uint32*)cigarBuf) & 0xf;
            if (op != BAM_CIGAR_H && op != BAM_CIGAR_S) {
                *o_refSpan += (*((_uint32*)cigarBuf) >> 4);
            }
            static const int op_ref[16] = {1, 0, 1, 1, 1, 1, 1, 1, 1, 0, 0, 0, 0, 0, 0, 0};
            for (int i = 1; i < cigarOps; i++) {
                _uint32 op = *(_uint32*)(cigarBuf + 4 * i);
                *o_refSpan += op_ref[op & 0xf] * (op >> 4);
            }
        }

        return used / 4;
    }
}

// Compute the CIGAR edit sequence operations in BAM format for a read against a given genome location
// Returns number of operations (or 0 if there was a problem)
// if returns with *o_addFrontClipping set non-zero, need to adjust front clipping & rerun
    int
BAMFormat::computeCigarOps(
    const Genome *              genome,
    AffineGapVectorizedWithCigar *ag,
    char *                      cigarBuf,
    int                         cigarBufLen,
    const char *                data,
    unsigned                    dataLength,
    unsigned                    basesClippedBefore,
    unsigned                    extraBasesClippedBefore,
    unsigned                    basesClippedAfter,
    unsigned                    frontHardClipping,
    unsigned                    backHardClipping,
    GenomeLocation              genomeLocation,
    bool                        isRC,
	bool						useM,
    int *                       o_editDistance,
    int *                       o_addFrontClipping,
    int *                       o_refSpan
)
{
    GenomeDistance extraBasesClippedAfter = 0;
    int used = 0;
    int backClippingMissedByLV = 0;

    unsigned clippingWordsBefore = ((basesClippedBefore + extraBasesClippedBefore > 0) ? 1 : 0) + ((frontHardClipping > 0) ? 1 : 0);
    unsigned clippingWordsAfter = ((basesClippedAfter + extraBasesClippedAfter > 0) ? 1 : 0) + ((backHardClipping > 0) ? 1 : 0);

    SAMFormat::computeCigar(BAM_CIGAR_OPS, genome, ag, cigarBuf + 4 * clippingWordsBefore, cigarBufLen - 4 * (clippingWordsBefore + clippingWordsAfter), data, dataLength, basesClippedBefore, extraBasesClippedBefore,
        basesClippedAfter, &extraBasesClippedAfter, genomeLocation, useM, o_editDistance, &used,  o_addFrontClipping, &backClippingMissedByLV);

    if (*o_addFrontClipping != 0) {
        return 0;
    }

    if (*o_editDistance == -2) {
        WriteErrorMessage("WARNING: computeGlobalScore returned -2; cigarBuf may be too small\n");
        return 0;
    } else if (*o_editDistance == -1) {
        static bool warningPrinted = false;
        if (!warningPrinted) {
            WriteErrorMessage("WARNING: computeGlobalScore returned -1; this shouldn't happen\n");
            warningPrinted = true;
        }
        return 0;
    } else {
        // There may be a better way to do this. For now, whenever we see tail insertions, soft-clip them
        basesClippedAfter += backClippingMissedByLV;
        dataLength -= backClippingMissedByLV;

        //
        // If we have hard clipping, add in the cigar string for it.
        //
        if (frontHardClipping > 0) {
            *(_uint32*)cigarBuf = (frontHardClipping << 4) | BAMAlignment::CigarToCode['H'];
            used += 4;
        }
        // Add some CIGAR instructions for soft-clipping if we've ignored some bases in the read.
        if (basesClippedBefore + extraBasesClippedBefore > 0) {
            *((_uint32*)cigarBuf + ((frontHardClipping > 0) ? 1 : 0))  = ((basesClippedBefore + extraBasesClippedBefore) << 4) | BAMAlignment::CigarToCode['S'];
            used += 4;
        }
        if (basesClippedAfter + extraBasesClippedAfter > 0) {
            *(_uint32*)(cigarBuf + used) = ((int)(basesClippedAfter + extraBasesClippedAfter) << 4) | BAMAlignment::CigarToCode['S'];
            used += 4;
        }

        if (backHardClipping > 0) {
            *(_uint32*)(cigarBuf + used) = (backHardClipping << 4) | BAMAlignment::CigarToCode['H'];
            used += 4;
        }

        // refSpan is used for calculating TLEN (e.g., read start 5' - mate end 5') correctly.
        // This is tricky for reads aligning as RC and needs to be inferred from CIGAR
        *o_refSpan = 0;
        const int cigarOps = used / 4;
        if (cigarOps > 0) {
            _uint32 op = *((_uint32*)cigarBuf) & 0xf;
            if (op != BAM_CIGAR_H && op != BAM_CIGAR_S) {
                *o_refSpan += (*((_uint32*)cigarBuf) >> 4);
            }
            static const int op_ref[16] = {1, 0, 1, 1, 1, 1, 1, 1, 1, 0, 0, 0, 0, 0, 0, 0};
            for (int i = 1; i < cigarOps; i++) {
                _uint32 op = *(_uint32*)(cigarBuf + 4 * i);
                *o_refSpan += op_ref[op & 0xf] * (op >> 4);
            }
        }

        return used / 4;
    }
}

class BAMFilter : public DataWriter::Filter
{
public:
    BAMFilter(DataWriter::FilterType i_type) : Filter(i_type), offsets(1000), header(false) {}

    virtual ~BAMFilter() {
        offsets.clear();
    }

	virtual void inHeader(bool flag)
	{ header = flag; }

    virtual void onAdvance(DataWriter* writer, size_t batchOffset, char* data, GenomeDistance bytes, GenomeLocation location);

    virtual size_t onNextBatch(DataWriter* writer, size_t offset, size_t bytes, bool lastBatch = false, bool* needMoreBuffer = NULL, size_t* fromBufferUsed = NULL);
    
protected:
    virtual void onRead(BAMAlignment* bam, size_t fileOffset, int batchIndex) = 0;

    BAMAlignment* getRead(size_t fileOffset);

    BAMAlignment* getNextRead(BAMAlignment* read, size_t* o_fileOffset = NULL);

    BAMAlignment* tryFindRead(size_t offset, size_t endOffset, const char* id, size_t* o_offset);

    VariableSizeVector<size_t> offsets;
    DataWriter* currentWriter;
    char* currentBuffer;
    size_t currentBufferBytes; // # of valid bytes
    size_t currentOffset; // logical file offset of beginning of current buffer

private:
    bool header;
};

    size_t
BAMFilter::onNextBatch(
    DataWriter* writer,
    size_t offset,
    size_t bytes,
    bool lastBatch,
    bool* needMoreBuffer,
    size_t* fromBufferUsed)
{

    // 
    // Nothing to write
    //
    if (bytes == 0 || *needMoreBuffer) {
        return 0;
    }

    bool ok = writer->getBatch(-1, &currentBuffer, NULL, NULL, NULL, &currentBufferBytes, &currentOffset);
    if (!ok) {
        WriteErrorMessage("Error writing to output file\n");
        soft_exit(1);
    }
    currentWriter = writer;
    int index = 0;

    size_t bytesToRead = min<long long>(bytes, currentBufferBytes);
    VariableSizeVector<size_t>::iterator i;
    for (i = offsets.begin(); i != offsets.end(); i++) {
        if (*i >= bytesToRead) {
            break;
        }
        onRead((BAMAlignment*) (currentBuffer + *i), currentOffset + *i, index++);
    }
    if (i != offsets.end()) {
        VariableSizeVector<size_t> nextBatchOffsets;
        VariableSizeVector<size_t>::iterator j;
        for (j = i; j != offsets.end(); j++) {
            nextBatchOffsets.push_back(*j - *i);
        }
        offsets.clear();
        for (j = nextBatchOffsets.begin(); j != nextBatchOffsets.end(); j++) {
            offsets.push_back(*j);
        }
        nextBatchOffsets.clear();
    }
    else {
        offsets.clear();
    }
    currentWriter = NULL;
    currentBuffer = NULL;
    currentBufferBytes = 0;
    currentOffset = 0;
    if (fromBufferUsed != NULL) {
        *fromBufferUsed = bytes;
    }
    return bytes;
}

    void
BAMFilter::onAdvance(
    DataWriter* writer,
    size_t batchOffset,
    char* data,
    GenomeDistance bytes,
    GenomeLocation location)
{
    if (! header) {
        offsets.push_back(batchOffset);
    }
}

    BAMAlignment*
BAMFilter::getRead(
    size_t offset)
{
    if (offset >= currentOffset && offset < currentOffset + currentBufferBytes) {
        return (BAMAlignment*) (currentBuffer + (offset - currentOffset));
    }
    for (int i = -2; ; i--) {
        char* buffer;
        size_t bufferFileOffset, bufferUsed; // logical
        if (! currentWriter->getBatch(i, &buffer, NULL, NULL, NULL, &bufferUsed, &bufferFileOffset)) {
            break;
        }
        if (offset >= bufferFileOffset && offset < bufferFileOffset + bufferUsed) {
            return (BAMAlignment*) (buffer + (offset - bufferFileOffset));
        }
    }
    return NULL;
}

    BAMAlignment*
BAMFilter::getNextRead(
    BAMAlignment* bam,
    size_t* io_offset)
{
    char* p = (char*) bam;
    size_t size = bam->size();
    size_t oldOffset = *io_offset;
    *io_offset += size;
    if (p >= currentBuffer && p < currentBuffer + currentBufferBytes) {
        p += bam->size();
        if (p >= currentBuffer + currentBufferBytes) {
            return NULL;
        }
        _ASSERT(*io_offset == currentOffset + (p - currentBuffer));
        _ASSERT(((BAMAlignment*)p)->refID >= -1);
        return (BAMAlignment*) p;
    }
    for (int i = -2; ; i--) {
        char* buffer;
        size_t bufferOffset, bufferUsed; // logical
        if (! currentWriter->getBatch(i, &buffer, NULL, NULL, NULL, &bufferUsed, &bufferOffset)) {
            break;
        }
        if (p >= buffer && p < buffer+ bufferUsed) {
            p += size;
            _ASSERT(*io_offset == bufferOffset + (p - buffer));
            return p < buffer + bufferUsed? (BAMAlignment*) p : getRead(*io_offset);
        }
    }
    return NULL;
}

    BAMAlignment*
BAMFilter::tryFindRead(
    size_t offset,
    size_t endOffset,
    const char* id,
    size_t* o_offset)
{
    BAMAlignment* bam = getRead(offset);
    size_t len = strcspn(id, "\0 /");
    while (bam != NULL && offset < endOffset) {
        if (readIdsMatch(bam->read_name(), id, len)) {
            if (o_offset != NULL) {
                *o_offset = offset;
            }
            return bam;
        }
        bam = getNextRead(bam, &offset);
    }
    return NULL;
}

//
// One of these per possible duplicate read (matches on location, next location, RC, next RC).
//
struct DuplicateReadKey
{
    DuplicateReadKey()
    { memset(this, 0, sizeof(DuplicateReadKey)); }

    DuplicateReadKey(BAMAlignment* bam, const Genome *genome, size_t _libraryHash)
    {
        if (bam == NULL) {
            locations[0] = locations[1] = InvalidGenomeLocation;
            isRC[0] = isRC[1] = false;
            libraryHash = 0;
        } else {
            locations[0] = bam->getLocation(genome);
            // locations[1] = bam->getNextLocation(genome);
            isRC[0] = (bam->FLAG & SAM_REVERSE_COMPLEMENT) != 0;
            isRC[1] = (bam->FLAG & SAM_NEXT_REVERSED) != 0;
            locations[0] = isRC[0] ? bam->getUnclippedEnd(locations[0]) : bam->getUnclippedStart(locations[0]);
            locations[1] = locations[0] + bam->tlen;
            libraryHash = _libraryHash;
            if (((((_uint64) GenomeLocationAsInt64(locations[0])) << 1) | (isRC[0] ? 1 : 0)) > ((((_uint64) GenomeLocationAsInt64(locations[1])) << 1) | (isRC[1] ? 1 : 0))) {
                const GenomeLocation t = locations[1];
                locations[1] = locations[0];
                locations[0] = t;
                const bool f = isRC[1];
                isRC[1] = isRC[0];
                isRC[0] = f;
            }
        }
    }

    bool operator==(const DuplicateReadKey& b) const
    {
        return libraryHash == b.libraryHash && locations[0] == b.locations[0] && locations[1] == b.locations[1] &&
            isRC[0] == b.isRC[0] && isRC[1] == b.isRC[1];
    }

    bool operator!=(const DuplicateReadKey& b) const
    {
        return ! ((*this) == b);
    }

    bool operator<(const DuplicateReadKey& b) const
    {
        return libraryHash < b.libraryHash ||
            (libraryHash == b.libraryHash && locations[0] < b.locations[0]) ||
            (locations[0] == b.locations[0] &&
                (locations[1] < b.locations[1] ||
                    (locations[1] == b.locations[1] &&
                        isRC[0] * 2 + isRC[1] <  b.isRC[0] *2 + b.isRC[1])));
    }


    // required for use as a key in VariableSizeMap template
    DuplicateReadKey(int x)
    { locations[0] = locations[1] = x; isRC[0] = isRC[1] = false; }
    bool operator==(int x) const
    { return locations[0] == (_uint32) x && locations[1] == (_uint32) x; }
    bool operator!=(int x) const
    { return locations[0] != (_uint32) x || locations[1] != (_uint32) x; }
    operator _uint64()
    { return ((_uint64) (GenomeLocationAsInt64(locations[1]) ^ (isRC[1] ? 1 : 0))) << 32 | (_uint64) (GenomeLocationAsInt64(locations[0]) ^ (isRC[0] ? 1 : 0)); }

    GenomeLocation locations[NUM_READS_PER_PAIR];
    bool isRC[NUM_READS_PER_PAIR];
    size_t libraryHash; // hash value of library name string LB:Z:xxxx
};

//
// Fragments duplicate marking does not require mate information
//
struct DuplicateFragmentKey
{
    DuplicateFragmentKey()
    {
        memset(this, 0, sizeof(DuplicateFragmentKey));
    }

    DuplicateFragmentKey(BAMAlignment* bam, const Genome* genome, size_t _libraryHash)
    {
        if (bam == NULL) {
            location = InvalidGenomeLocation;
            isRC = false;
            libraryHash = 0;
        }
        else {
            location = bam->getLocation(genome);
            isRC = (bam->FLAG & SAM_REVERSE_COMPLEMENT) != 0;
            location = isRC ? bam->getUnclippedEnd(location) : bam->getUnclippedStart(location);
            libraryHash = _libraryHash;
        }
    }

    bool operator==(const DuplicateFragmentKey& b) const
    {
        return libraryHash == b.libraryHash && location == b.location && isRC == b.isRC;
    }

    bool operator!=(const DuplicateFragmentKey& b) const
    {
        return !((*this) == b);
    }

    bool operator<(const DuplicateFragmentKey& b) const
    {
        return libraryHash < b.libraryHash ||
            (libraryHash == b.libraryHash && location < b.location) ||
            (location == b.location && isRC < b.isRC);
    }

    // required for use as a key in VariableSizeMap template
    DuplicateFragmentKey(int x)
    {
        location = x; isRC = false;
    }
    bool operator==(int x) const
    {
        return location == (_uint32)x;
    }
    bool operator!=(int x) const
    {
        return location != (_uint32)x;
    }
    operator _uint64()
    {
        return (_uint64)(GenomeLocationAsInt64(location) ^ (isRC ? 1 : 0));
    }

    GenomeLocation location;
    bool isRC;
    size_t libraryHash;
};

struct DuplicateMateInfo
{
    DuplicateMateInfo()
    {
        memset(this, 0, sizeof(DuplicateMateInfo));
    }

    ~DuplicateMateInfo()
    {
    }

    bool isMateMapped;
    int bestReadQuality; // total quality of first/both best reads
    char bestReadId[120];
    
    //
    // The below are useful for marking optical duplicates. Helps break ties in Illumina reads
    //
    int tile;
    int x;
    int y;

    void setBestReadId(const char* id) { strncpy(bestReadId, id, sizeof(bestReadId)); }
    const char* getBestReadId() { return bestReadId; }
    void setBestTileXY(int tile_, int x_, int y_) { tile = tile_; x = x_; y = y_; }
    void getBestTileXY(int* tile_, int* x_, int* y_) { *tile_ = tile; *x_ = x; *y_ = y; }

    void checkBestRecord(BAMAlignment* bam, int totalQuality_, int tile_, int x_, int y_) {

        if ((bam->FLAG & SAM_DUPLICATE) != 0) {
            return;
        }
        if (totalQuality_ > bestReadQuality) {
            bestReadQuality = totalQuality_;
            setBestReadId(bam->read_name());
            setBestTileXY(tile_, x_, y_);
        }
        else if (totalQuality_ == bestReadQuality) {
            if (tile_ < tile) {
                bestReadQuality = totalQuality_;
                setBestReadId(bam->read_name());
                setBestTileXY(tile_, x_, y_);
            }
            else if (tile_ == tile) {
                if (x_ < x) {
                    bestReadQuality = totalQuality_;
                    setBestReadId(bam->read_name());
                    setBestTileXY(tile_, x_, y_);
                }
                else if (x_ == x) {
                    if (y_ < y) {
                        bestReadQuality = totalQuality_;
                        setBestReadId(bam->read_name());
                        setBestTileXY(tile_, x_, y_);
                    }
                }
            }
        }
    }
};

struct BamDupMarkEntry
{
    BamDupMarkEntry() : libraryNameHash(0), runOffset(0), mateQual(0), mateInfo(0), info(0) {}

    bool operator<(const BamDupMarkEntry& b) const {
        return (libraryNameHash < b.libraryNameHash) ||
            ((libraryNameHash == b.libraryNameHash) && (info < b.info)) ||
            ((libraryNameHash == b.libraryNameHash && info == b.info && mateInfo < b.mateInfo));
    }

    static size_t hash(char* str) // We could probably use a better hash function, but I doubt it matters
    {
        size_t value = 0x123456789abcdef0;
        char* libraryName = str;

        for (int i = 0; i < strlen(str); i++)
        {
            value = (value * 131) ^ libraryName[i];
        }

        return value;
    } // hash

    size_t libraryNameHash; // hash value of library name string LB:Z:xxxx
    size_t runOffset; // read offset in run
    _int32 mateQual; // mate quality score (only considering bases with phred score >= 15)
    _uint64 mateInfo; // mate information has: matelocation and matedirection
    _uint64 info; // read information has: location and direction
};

class BAMDupMarkFilter : public BAMFilter
{
public:
    BAMDupMarkFilter(const Genome* i_genome) :
        BAMFilter(DataWriter::DupMarkFilter),
        genome(i_genome), runOffset(0), runLocation(InvalidGenomeLocation), prevRunLocation(InvalidGenomeLocation), runCount(0), mates(), fragments()
    {
    }

    ~BAMDupMarkFilter()
    {
#ifdef USE_DEVTEAM_OPTIONS
        if (mates.size() > 0) {
            WriteErrorMessage("duplicate matching ended with %d unmatched reads:\n", mates.size());
            for (MateMap::iterator i = mates.begin(); i != mates.end(); i = mates.next(i)) {
	      WriteErrorMessage("%u%s/%u%s\n", GenomeLocationAsInt64(i->key.locations[0]), i->key.isRC[0] ? "rc" : "", GenomeLocationAsInt64(i->key.locations[1]), i->key.isRC[1] ? "rc" : "");
            }
        }
#endif
        run.clear();
        runFragment.clear();
    }

    static bool isDuplicate(const BAMAlignment* a, const BAMAlignment* b)
    { return a->pos == b->pos && a->refID == b->refID &&
        ((a->FLAG ^ b->FLAG) & (SAM_REVERSE_COMPLEMENT | SAM_NEXT_REVERSED)) == 0; }

    virtual size_t onNextBatch(DataWriter* writer, size_t offset, size_t bytes, bool lastBatch = false, bool* needMoreBuffer = NULL, size_t* fromBufferUsed = NULL);

    void dupMarkBatch(BAMAlignment* lastBam, size_t lastOffset);

protected:
    virtual void onRead(BAMAlignment* bam, size_t fileOffset, int batchIndex);

private:
    static int getTotalQuality(BAMAlignment* bam);

    static void getTileXY(const char* id, int* o_tile, int* o_x, int* o_y);

    const Genome* genome;
    size_t runOffset; // offset in file of first read in run
    GenomeLocation runLocation; // location in genome
    GenomeLocation prevRunLocation; // location in genome
    int runCount; // number of aligned reads

    typedef VariableSizeMap<DuplicateReadKey,DuplicateMateInfo,150,MapNumericHash<DuplicateReadKey>,70,0,-2> MateMap;
    typedef VariableSizeMap<DuplicateFragmentKey, DuplicateMateInfo, 150, MapNumericHash<DuplicateFragmentKey>, 70, 0, -2> FragmentMap;
    typedef VariableSizeVector<BamDupMarkEntry> RunVector;

    RunVector run; // used for paired-end duplicate marking
    RunVector runFragment; // used for single-end duplicate marking 
    MateMap mates;
    FragmentMap fragments;
};

    size_t
BAMDupMarkFilter::onNextBatch(
    DataWriter* writer,
    size_t offset,
    size_t bytes,
    bool lastBatch,
    bool* needMoreBuffer,
    size_t* fromBufferUsed)
{
    // 
    // Nothing to write
    //
    if (bytes == 0 || *needMoreBuffer) {
        return 0;
    }

    bool ok = writer->getBatch(-1, &currentBuffer, NULL, NULL, NULL, &currentBufferBytes, &currentOffset);
    if (!ok) {
        WriteErrorMessage("Error writing to output file\n");
        soft_exit(1);
    }
    currentWriter = writer;

    size_t next_i = 0, unfinishedRunStart = 0;
    bool foundNextBatchStart = false;
    size_t nextBatchStart = 0;
    runCount = 0;
    runOffset = 0;
    runLocation = InvalidGenomeLocation;
    BAMAlignment* lastBam = NULL;
    BAMAlignment* firstBam = NULL;

    for (size_t i = 0; i < offsets.size(); i = next_i) {
        lastBam = (BAMAlignment*) (currentBuffer + offsets[i]);
        GenomeLocation location = lastBam->getLocation(genome);
        GenomeLocation nextLocation = lastBam->getNextLocation(genome);
        GenomeLocation logicalLocation = location != InvalidGenomeLocation ? location : nextLocation;
        logicalLocation = lastBam->getUnclippedStart(logicalLocation);

        //
        // Initialize run
        //
        if (runLocation == InvalidGenomeLocation) {
            runCount = 1;
            runLocation = logicalLocation;
            runOffset = currentOffset + offsets[i];
            unfinishedRunStart = i;
            next_i = i + 1;
            firstBam = (BAMAlignment*)(currentBuffer + offsets[i]);
        }
        else {
            // 
            // Track the read from which we need to start the next run. Next run starts at nextBatchStart
            // 
            if (!foundNextBatchStart && logicalLocation > runLocation + MAX_READ_LENGTH + MAX_K) {
                nextBatchStart = i;
                foundNextBatchStart = true;
                next_i = i + 1;
            }
            //
            // Add read to run
            // 
            if (logicalLocation <= runLocation + (2 * (MAX_READ_LENGTH + MAX_K))) {
                runCount++;
                next_i = i + 1;
            }
            else {
                // 
                // We are done with the run. Begin marking duplicates
                // 
                dupMarkBatch(lastBam, currentOffset + offsets[i]);

                //
                // Reset run parameters in preparation for next run
                // 
                runCount = 0;
                runOffset = 0;
                runLocation = InvalidGenomeLocation;
                next_i = nextBatchStart;
                foundNextBatchStart = false;
            }
        }
    } // end offsets

    size_t bytesRead = min<long long>(bytes, currentBufferBytes);

    //
    // Run could potentially span multiple buffers
    //
    if (runCount > 1) { // runs of size <= 1 cannot have duplicates
        if (lastBatch) {
            //
            // We have reached the lastBatch. Mark all duplicates
            //
            dupMarkBatch(lastBam, currentOffset + offsets[offsets.size() - 1]);
            offsets.clear();
        }
        else if (offsets[unfinishedRunStart] == 0) { // we did not yet mark duplicates for this run
            //
            // If we have a different run from what we have seen before, 
            // simply mark all duplicates in the run we currently have
            //
            if (runLocation != prevRunLocation) {
                dupMarkBatch(lastBam, currentOffset + offsets[offsets.size() - 1]);
                offsets.clear();
            }
            else {
                *needMoreBuffer = true;
                if (fromBufferUsed != NULL) {
                    *fromBufferUsed = 0;
                }
                return 0;
            }
        }
        else {
            //
            // Copy over read offsets for those reads that will be duplicate marked in the next batch
            //
            bytesRead = offsets[unfinishedRunStart];
            VariableSizeVector<size_t> nextBatchOffsets;
            for (size_t i = unfinishedRunStart; i < offsets.size(); i++) {
                nextBatchOffsets.push_back(offsets[i] - offsets[unfinishedRunStart]);
            }
            offsets.clear();
            for (size_t i = 0; i < nextBatchOffsets.size(); i++) {
                offsets.push_back(nextBatchOffsets[i]);
            }
            nextBatchOffsets.clear();
        }
    } // runcount > 1
    else {
        offsets.clear();
    }

    prevRunLocation = runLocation;
    currentWriter = NULL;
    currentBuffer = NULL;
    currentBufferBytes = 0;
    currentOffset = 0;

    if (fromBufferUsed != NULL) {
        *fromBufferUsed = bytesRead;
    }

    return bytesRead; // return bytes consumed in current batch. This may be less than 'bytes', if reads require duplicate marking in subsequent batch
}

    void
BAMDupMarkFilter::dupMarkBatch(BAMAlignment* lastBam, size_t lastOffset) {

    // partition by duplicate key, find best read in each partition
    size_t offset = runOffset;
    int numRecords = 0;
    run.clear();
    runFragment.clear();

    for (BAMAlignment* record = getRead(offset); record != NULL && numRecords < runCount; record = getNextRead(record, &offset)) {
        
        // secondary/supplementary alignments not marked as duplicates in coordinate sorted alignments
        if ((record->FLAG & SAM_SECONDARY) != 0 || (record->FLAG & SAM_SUPPLEMENTARY) != 0) continue;
        
        GenomeLocation loc = record->getLocation(genome);
        bool isRC = (record->FLAG & SAM_REVERSE_COMPLEMENT) != 0;
        bool isMateRC = (record->FLAG & SAM_NEXT_REVERSED) != 0;
        GenomeLocation myLoc = isRC ? record->getUnclippedEnd(loc) : record->getUnclippedStart(loc);
        GenomeLocation mateLoc = myLoc + record->tlen;

        BamDupMarkEntry entry, entryFragment;
        BAMAlignAux* aux = record->firstAux();
        entry.mateQual = entryFragment.mateQual = -1;
        entry.libraryNameHash = entryFragment.libraryNameHash = 0;
        bool foundLibraryTag = false;

        //
        // Extract QS (mate quality score) and LB fields from BAM aux tags
        //
        while (aux != NULL && aux->isValidValType() && aux < record->endAux()) {
            if (aux->tag[0] == 'Q' && aux->tag[1] == 'S' && aux->val_type == 'i') {
                entry.mateQual = entryFragment.mateQual = *(_int32*)aux->value();
            }
            if (!foundLibraryTag && aux->tag[0] == 'L' && aux->tag[1] == 'B' && aux->val_type == 'Z') {
                foundLibraryTag = true;
                // fixme: conflicts from hashing library names to the same value
                entry.libraryNameHash = entryFragment.libraryNameHash = BamDupMarkEntry::hash((char*)aux->value());
            }
            aux = aux->next();
        }

        entry.mateInfo = (((_uint64) GenomeLocationAsInt64(mateLoc)) << 1) | (isMateRC ? 1 : 0);
        entry.info = (((_uint64) GenomeLocationAsInt64(myLoc)) << 1) | (isRC ? 1 : 0);

        entryFragment.mateInfo = 0;
        entryFragment.info = (((_uint64)GenomeLocationAsInt64(myLoc)) << 1) | (isRC ? 1 : 0);

        entry.runOffset = entryFragment.runOffset = offset;
        
        // don't need to look at mate information in single-ended datasets
        if ((record->FLAG & SAM_MULTI_SEGMENT) != 0) {
            run.push_back(entry);
        }

        runFragment.push_back(entryFragment);
        numRecords++;
    }

    if (run.size() == 0 && runFragment.size() == 0) {
        return;
    }
    // ensure that duplicates reads are adjacent
    std::stable_sort(run.begin(), run.end());
    std::stable_sort(runFragment.begin(), runFragment.end());

    bool foundRun = false;
    for (RunVector::iterator i = run.begin(); i != run.end(); i++) {
        offset = i->runOffset;
        BAMAlignment* record = getRead(offset);
        _ASSERT(record->refID >= -1 && record->refID < genome->getNumContigs()); // simple sanity check
        
        // skip unmapped reads and reads with unmapped mates
        if (((record->FLAG & SAM_UNMAPPED) != 0) || (record->FLAG & SAM_NEXT_UNMAPPED) != 0) continue;

        // adjacent entries with different library names
        if ((i == run.begin() || (i->libraryNameHash != ((i - 1)->libraryNameHash))) &&
            (i + 1 == run.end() || (i->libraryNameHash != ((i + 1)->libraryNameHash)))) {
            continue;
        }

        // adjacent entries with different location/orientation
        if ((i == run.begin() || (i->info) != ((i - 1)->info)) &&
            (i + 1 == run.end() || (i->info) != ((i + 1)->info))) {
            continue;
        }

        // check if mate location/orientation matches adjacent entry
        if ((i == run.begin() || (i->mateInfo) != ((i - 1)->mateInfo)) &&
            (i + 1 == run.end() || (i->mateInfo) != ((i + 1)->mateInfo))) {
            continue;
        }
        foundRun = true;
        DuplicateReadKey key(record, genome, i->libraryNameHash);
        MateMap::iterator f = mates.find(key);
        DuplicateMateInfo* info;
        if (f == mates.end()) {
            mates.put(key, DuplicateMateInfo());
            info = &mates[key];
            //fprintf(stderr, "add %u%s/%u%s -> %d\n", key.locations[0], key.isRC[0] ? "rc" : "", key.locations[1], key.isRC[1] ? "rc" : "", mates.size());
            info->isMateMapped = true;
        } else {
            info = &f->value;
        }
        int totalQuality = getTotalQuality(record);
        if ((record->FLAG & SAM_MULTI_SEGMENT) != 0) {
            totalQuality += i->mateQual;
        }
        int tile, x, y;
        getTileXY(record->read_name(), &tile, &x, &y); // parse read name and extract metadata for optical duplicate marking
        info->checkBestRecord(record, totalQuality, tile, x, y); // update best record if needed
    }

    // duplicate marking for read fragments
    for (RunVector::iterator i = runFragment.begin(); i != runFragment.end(); i++) {
        offset = i->runOffset;
        BAMAlignment* record = getRead(offset);

        if ((i == runFragment.begin() || (i->libraryNameHash != ((i - 1)->libraryNameHash))) &&
            (i + 1 == runFragment.end() || (i->libraryNameHash != ((i + 1)->libraryNameHash)))) {
            continue;
        }

        if ((i == runFragment.begin() || (i->info) != ((i - 1)->info)) &&
            (i + 1 == runFragment.end() || (i->info) != ((i + 1)->info))) {
            continue;
        }

        // location and library matches
        DuplicateFragmentKey key(record, genome, i->libraryNameHash);
        foundRun = true;
        FragmentMap::iterator f = fragments.find(key);
        DuplicateMateInfo* info;
        if (f == fragments.end()) {
            fragments.put(key, DuplicateMateInfo());
            info = &fragments[key];
            //fprintf(stderr, "add %u%s/%u%s -> %d\n", key.locations[0], key.isRC[0] ? "rc" : "", key.locations[1], key.isRC[1] ? "rc" : "", mates.size());
            info->isMateMapped = (record->FLAG & SAM_MULTI_SEGMENT) != 0 && (record->FLAG & SAM_NEXT_UNMAPPED) == 0;
        }
        else {
            info = &f->value;
        }
        bool mateMapped = (record->FLAG & SAM_MULTI_SEGMENT) != 0 && (record->FLAG & SAM_NEXT_UNMAPPED) == 0;
        int totalQuality = getTotalQuality(record);
        int tile, x, y;
        getTileXY(record->read_name(), &tile, &x, &y);
        
        //
        // Prefer mapped read pairs over fragments in duplicate marking
        //
        if (mateMapped) {
            if (!info->isMateMapped) {
                info->bestReadQuality = totalQuality;
                info->setBestReadId(record->read_name());
                info->isMateMapped = true;
                info->setBestTileXY(tile, x, y);
            }
            else {
                info->checkBestRecord(record, totalQuality, tile, x, y);
            }
        }
        else {
            //
            // No best read pair found so far.
            //
            if (!info->isMateMapped) {
                info->checkBestRecord(record, totalQuality, tile, x, y);
            }
        }
    }

    if (!foundRun) {
        return; // avoid useless looping
    }

    // go back and adjust flags
    for (RunVector::iterator i = run.begin(); i != run.end(); i++) {
        if ((i == run.begin() || (i->libraryNameHash != ((i - 1)->libraryNameHash))) &&
            (i + 1 == run.end() || (i->libraryNameHash != ((i + 1)->libraryNameHash)))) {
            continue;
        }
        if ((i == run.begin() || (i->info) != ((i - 1)->info)) &&
            (i + 1 == run.end() || (i->info) != ((i + 1)->info))) {
            continue;
        }
        if ((i == run.begin() || (i->mateInfo) != ((i - 1)->mateInfo)) &&
            (i + 1 == run.end() || (i->mateInfo) != ((i + 1)->mateInfo))) {
            continue;
        }

        offset = i->runOffset;
        BAMAlignment* record = getRead(offset);

        // skip unmapped reads and reads with unmapped mates
        if (((record->FLAG & SAM_UNMAPPED) != 0) || (record->FLAG & SAM_NEXT_UNMAPPED) != 0) continue;

        DuplicateReadKey key(record, genome, i->libraryNameHash);
        MateMap::iterator m = mates.find(key);
        if (m == mates.end()) {
            continue;
        }
        DuplicateMateInfo* minfo = &m->value;
        if (!readIdsMatch(minfo->getBestReadId(), record->read_name(), record->l_read_name - 1)) {
            record->FLAG |= SAM_DUPLICATE;
        }
    }

    // clean up
    for (RunVector::iterator i = run.begin(); i != run.end(); i++) {

        // skip singletons
        if ((i == run.begin() || (i->libraryNameHash != ((i - 1)->libraryNameHash))) &&
            (i + 1 == run.end() || (i->libraryNameHash != ((i + 1)->libraryNameHash)))) {
            continue;
        }
        if ((i == run.begin() || (i->info) != ((i - 1)->info)) &&
            (i + 1 == run.end() || (i->info) != ((i + 1)->info))) {
            continue;
        }
        if ((i == run.begin() || (i->mateInfo) != ((i - 1)->mateInfo)) &&
            (i + 1 == run.end() || (i->mateInfo) != ((i + 1)->mateInfo))) {
            continue;
        }

        offset = i->runOffset;
        BAMAlignment* record = getRead(offset);

        // skip unmapped reads and reads with unmapped mates
        if (((record->FLAG & SAM_UNMAPPED) != 0) || (record->FLAG & SAM_NEXT_UNMAPPED) != 0) continue;

        DuplicateReadKey key(record, genome, i->libraryNameHash);
        MateMap::iterator m = mates.find(key);
        DuplicateMateInfo* minfo = &m->value;
        if (m != mates.end()) {
            bool isRC = (record->FLAG & SAM_REVERSE_COMPLEMENT) != 0;
            GenomeLocation loc = record->getLocation(genome);
            loc = (isRC) ? record->getUnclippedEnd(loc) : record->getUnclippedStart(loc);
            // 
            // Keep duplicate entry around till we find the mate. This allows us to match indexInFile tie breaking used in Picard MarkDup
            //
            if (loc == key.locations[1] && isRC == key.isRC[1]) {
                //fprintf(stderr, "erase %u%s/%u%s -> %d\n", key.locations[0], key.isRC[0] ? "rc" : "", key.locations[1], key.isRC[1] ? "rc" : "", mates.size());
                mates.erase(key);
            }
        }
    }

    // handle fragments
    for (RunVector::iterator i = runFragment.begin(); i != runFragment.end(); i++) {
        offset = i->runOffset;
        BAMAlignment* record = getRead(offset);

        if ((i == runFragment.begin() || (i->libraryNameHash != ((i - 1)->libraryNameHash))) &&
            (i + 1 == runFragment.end() || (i->libraryNameHash != ((i + 1)->libraryNameHash)))) {
            continue;
        }

        if ((i == runFragment.begin() || (i->info) != ((i - 1)->info)) &&
            (i + 1 == runFragment.end() || (i->info) != ((i + 1)->info))) {
            continue;
        }

        //
        // Skip unmapped reads and reads with mapped mates
        //
        if ((record->FLAG & SAM_UNMAPPED) != 0 || ((record->FLAG & SAM_MULTI_SEGMENT) != 0 && (record->FLAG & SAM_NEXT_UNMAPPED) == 0)) {
            continue;
        }

        // location and library matches
        DuplicateFragmentKey key(record, genome, i->libraryNameHash);
        FragmentMap::iterator f = fragments.find(key);

        if (f == fragments.end()) {
            continue;
        }
        DuplicateMateInfo* info = &f->value;
        if (!readIdsMatch(info->getBestReadId(), record->read_name(), record->l_read_name - 1)) {
            record->FLAG |= SAM_DUPLICATE;
        }
    }

    // clean up
    for (RunVector::iterator i = runFragment.begin(); i != runFragment.end(); i++) {
        offset = i->runOffset;
        BAMAlignment* record = getRead(offset);

        if ((i == runFragment.begin() || (i->libraryNameHash != ((i - 1)->libraryNameHash))) &&
            (i + 1 == runFragment.end() || (i->libraryNameHash != ((i + 1)->libraryNameHash)))) {
            continue;
        }

        if ((i == runFragment.begin() || (i->info) != ((i - 1)->info)) &&
            (i + 1 == runFragment.end() || (i->info) != ((i + 1)->info))) {
            continue;
        }

        //
        // Skip unmapped reads and reads with mapped mates
        //
        if ((record->FLAG & SAM_UNMAPPED) != 0 || ((record->FLAG & SAM_MULTI_SEGMENT) != 0 && (record->FLAG & SAM_NEXT_UNMAPPED) == 0)) {
            continue;
        }

        // location and library matches
        DuplicateFragmentKey key(record, genome, i->libraryNameHash);
        FragmentMap::iterator f = fragments.find(key);

        if (f != fragments.end()) {
            fragments.erase(key);
        }
    }
}


    void
BAMDupMarkFilter::onRead(BAMAlignment* lastBam, size_t lastOffset, int)
{
    //
    // Do nothing
    //
}

    int
BAMDupMarkFilter::getTotalQuality(
    BAMAlignment* bam)
{
    int result = 0;
    _uint8* p = (_uint8*) bam->qual();
    for (int i = 0; i < bam->l_seq; i++) {
        int q = *p++;
        // Picard MarkDup uses a score threshold of 15 (default)
        result += (q >= 15) ? (q != 255) * q : 0; // avoid branch?
    }
    return result;
}

    void
BAMDupMarkFilter::getTileXY(
    const char* id,
    int* o_tile,
    int* o_x,
    int* o_y)
{
    // In the 5-element format read names: tile, x and y values are 3rd, 4th and 5th elements
    // In the 7-element format (CASAVA 1.8) read names: tile, x and y values are 5th, 6th and 7th elements
    int i = 0, numColonsSeen = 0, tile = 0, x = 0, y = 0;
    char readId[120];
    strncpy(readId, id, sizeof(readId));
    readId[119] = '\0';
    int fiveElementFormatStartIndex = 0, sevenElementFormatStartIndex;
    for (i = 0; ; i++) {
        char c = readId[i];
        if (c == ':') {
            numColonsSeen++;
            if (numColonsSeen == 2) {
                fiveElementFormatStartIndex = i;
            }
            else if (numColonsSeen == 4) {
                sevenElementFormatStartIndex = i;
            }
        }
        if (c == 0 || c == ' ' || c == '/') {
            break;
        }
    }
    int fieldScanned;
    if (numColonsSeen == 4) {
        fieldScanned = sscanf(&readId[fiveElementFormatStartIndex], ":%d:%d:%d", &tile, &x, &y);
    }
    else if (numColonsSeen == 6) {
        fieldScanned = sscanf(&readId[sevenElementFormatStartIndex], ":%d:%d:%d", &tile, &x, &y);
    }
    // fprintf(stderr, "tile:%d, x:%d, y:%d\n", tile, x, y);
    *o_tile = tile;
    *o_x = x;
    *o_y = y;
}

class BAMDupMarkSupplier : public DataWriter::FilterSupplier
{
public:
    BAMDupMarkSupplier(const Genome* i_genome) :
        FilterSupplier(DataWriter::ReadFilter), genome(i_genome) {}

    virtual DataWriter::Filter* getFilter()
    { return new BAMDupMarkFilter(genome); }

    virtual void onClosing(DataWriterSupplier* supplier) {}
    virtual void onClosed(DataWriterSupplier* supplier) {}

private:
    const Genome* genome;
};

    DataWriter::FilterSupplier*
DataWriterSupplier::bamMarkDuplicates(const Genome* genome)
{
    return new BAMDupMarkSupplier(genome);
}

class BAMIndexSupplier;

class BAMIndexFilter : public BAMFilter
{
public:
    BAMIndexFilter(BAMIndexSupplier* i_supplier)
        : BAMFilter(DataWriter::ReadFilter), supplier(i_supplier) {}

protected:
    virtual void onRead(BAMAlignment* bam, size_t fileOffset, int batchIndex);

private:
    BAMIndexSupplier* supplier;
};

class BAMIndexSupplier : public DataWriter::FilterSupplier
{
public:
    BAMIndexSupplier(const char* i_indexFileName, const Genome* i_genome, GzipWriterFilterSupplier* i_gzipSupplier) :
        FilterSupplier(DataWriter::ReadFilter),
        indexFileName(i_indexFileName),
        genome(i_genome),
        gzipSupplier(i_gzipSupplier),
        lastRefId(-1),
        lastBin(0), binStart(0), lastBamEnd(0)
    {
        refs = genome ? new RefInfo[genome->getNumContigs()] : NULL;
        readCounts[0] = readCounts[1] = 0;
    }

    virtual DataWriter::Filter* getFilter()
    { return new BAMIndexFilter(this); }

    virtual void onClosing(DataWriterSupplier* supplier) {}
    virtual void onClosed(DataWriterSupplier* supplier);

private:

    friend class BAMIndexFilter;

    struct BAMChunk {
        BAMChunk() : start(0), end(0) {}
        BAMChunk(const BAMChunk& a) : start(a.start), end(a.end) {}

        _uint64 start, end;
    };
    typedef VariableSizeVector<BAMChunk> ChunkVec;
    typedef VariableSizeMap<_uint32,ChunkVec,150,MapNumericHash<_uint32>,80,-1,-2> BinMap;
    typedef VariableSizeVector<_uint64> LinearMap;
    struct RefInfo {
        BinMap bins;
        LinearMap intervals;
    };

    RefInfo* getRefInfo(int refId);

    void onRead(BAMAlignment* bam, size_t fileOffset, int batchIndex);

    void addChunk(int refId, _uint32 bin, _uint64 start, _uint64 end);

    void addInterval(int refId, int begin, int end, _uint64 fileOffset);

    const char* indexFileName;
    const Genome* genome;
    int lastRefId;
    _uint32 lastBin;
    _uint64 binStart;
    _uint64 firstBamStart;
    _uint64 lastBamEnd;
    _uint64 readCounts[2]; // mapped, unmapped
    RefInfo* refs;
    GzipWriterFilterSupplier* gzipSupplier;
};

    void
BAMIndexFilter::onRead(
    BAMAlignment* bam,
    size_t fileOffset,
    int batchIndex)
{
    supplier->onRead(bam, fileOffset, batchIndex);
}

    DataWriter::FilterSupplier*
DataWriterSupplier::bamIndex(
    const char* indexFileName,
    const Genome* genome,
    GzipWriterFilterSupplier* gzipSupplier)
{
    return new BAMIndexSupplier(indexFileName, genome, gzipSupplier);
}

    void
BAMIndexSupplier::onRead(
    BAMAlignment* bam,
    size_t fileOffset,
    int batchIndex)
{
    //fprintf(stderr, "index onRead %d:%d+%d @ %lld %d\n", bam->refID, bam->pos, bam->l_ref(), fileOffset, batchIndex);
    if (bam->refID != lastRefId) {
        if (lastRefId != -1) {
            addChunk(lastRefId, BAMAlignment::BAM_EXTRA_BIN, firstBamStart, lastBamEnd);
            addChunk(lastRefId, BAMAlignment::BAM_EXTRA_BIN, readCounts[0], readCounts[1]);
            readCounts[0] = readCounts[1] = 0;
        }
        firstBamStart = fileOffset;
    }
    readCounts[(bam->FLAG & SAM_UNMAPPED) ? 1 : 0]++;
    if (bam->refID != lastRefId || bam->bin != lastBin || lastRefId == -1) {
        addChunk(lastRefId, lastBin, binStart, fileOffset);
        lastBin = bam->bin;
        lastRefId = bam->refID;
        binStart = fileOffset;
    }
    if (! (bam->FLAG & SAM_UNMAPPED)) {
        _ASSERT(bam->pos != -1 && bam->refID != -1);
        addInterval(bam->refID, bam->pos, bam->pos + bam->l_ref() - 1, fileOffset);
    }
    lastBamEnd = fileOffset + bam->size();
}

    void
BAMIndexSupplier::onClosed(
    DataWriterSupplier* supplier)
{
    // add final chunk
    if (lastRefId != -1) {
        addChunk(lastRefId, lastBin, binStart, lastBamEnd);
        addChunk(lastRefId, BAMAlignment::BAM_EXTRA_BIN, firstBamStart, lastBamEnd);
        addChunk(lastRefId, BAMAlignment::BAM_EXTRA_BIN, readCounts[0], readCounts[1]);
    }

    // write out index file
    FILE* index = fopen(indexFileName, "wb");
    char magic[4] = {'B', 'A', 'I', 1};
    fwrite(magic, sizeof(magic), 1, index);
    _int32 n_ref = genome->getNumContigs();
    fwrite(&n_ref, sizeof(n_ref), 1, index);

    for (int i = 0; i < n_ref; i++) {
        RefInfo* info = getRefInfo(i);
        _int32 n_bin, n_intv;
        if (info == NULL) {
            n_bin = 0;
            fwrite(&n_bin, sizeof(n_bin), 1, index);
            n_intv = 0;
            fwrite(&n_intv, sizeof(n_intv), 1, index);
            continue;
        }
        n_bin = info->bins.size();
        fwrite(&n_bin, sizeof(n_bin), 1, index);
        for (BinMap::iterator j = info->bins.begin(); j != info->bins.end(); j = info->bins.next(j)) {
            _uint32 bin = j->key;
            fwrite(&bin, sizeof(bin), 1, index);
            _int32 n_chunk = (_int32) j->value.size();
            fwrite(&n_chunk, sizeof(n_chunk), 1, index);
            if (bin != BAMAlignment::BAM_EXTRA_BIN) {
                for (ChunkVec::iterator k = j->value.begin(); k != j->value.end(); k++) {
                    _uint64 chunk[2] = {gzipSupplier->toVirtualOffset(k->start), gzipSupplier->toVirtualOffset(k->end)};
                    fwrite(&chunk, sizeof(chunk), 1, index);
                }
            } else {
                _uint64 chunk[2] = {gzipSupplier->toVirtualOffset(j->value[0].start), gzipSupplier->toVirtualOffset(j->value[0].end)};
                fwrite(&chunk, sizeof(chunk), 1, index);
                chunk[0] = j->value[1].start;
                chunk[1] = j->value[1].end;
                fwrite(&chunk, sizeof(chunk), 1, index);
            }
        }
        n_intv = (_int32) info->intervals.size();
        fwrite(&n_intv, sizeof(n_intv), 1, index);
        for (LinearMap::iterator m = info->intervals.begin(); m != info->intervals.end(); m++) {
            _uint64 ioffset = gzipSupplier->toVirtualOffset(*m);
            fwrite(&ioffset, sizeof(ioffset), 1, index);
        }
    }
    fclose(index);
}

   BAMIndexSupplier::RefInfo*
BAMIndexSupplier::getRefInfo(
    int refId)
{
    return refId >= 0 && refId < genome->getNumContigs() ? &refs[refId] : NULL;
}

    void
BAMIndexSupplier::addChunk(
    int refId,
    _uint32 bin,
    _uint64 start,
    _uint64 end)
{
    RefInfo* info = getRefInfo(refId);
    if (info == NULL) {
        return;
    }
    ChunkVec* chunks = info->bins.tryFind(bin);
    if (chunks == NULL) {
        ChunkVec empty;
        info->bins.tryAdd(bin, empty, &chunks);
    }
    BAMChunk chunk;
    chunk.start = start;
    chunk.end = end;
    chunks->push_back(chunk);
}

    void
BAMIndexSupplier::addInterval(
    int refId,
    int begin,
    int end,
    _uint64 fileOffset)
{
    RefInfo* info = getRefInfo(refId);
    if (info == NULL) {
        return;
    }
    int slot = end <= 0 ? 0 : ((end - 1) / 16384);
    if (slot >= info->intervals.size()) {
        for (_int64 i = info->intervals.size(); i < slot; i++) {
            info->intervals.push_back(UINT64_MAX);
        }
        info->intervals.push_back(fileOffset);
    }
}

    bool
BgzfHeader::validate(char* buffer, size_t bytes)
{
    char* p;
    for (p = buffer; p - buffer < (_int64)bytes; ) {
        BgzfHeader* h = (BgzfHeader*) p;
        unsigned bsize = h->BSIZE() + 1;
        unsigned isize = h->ISIZE();
        if (bsize == 0 || bsize > BAM_BLOCK || isize > BAM_BLOCK ||
                bsize > max(2 * isize, isize+1000) || ! h->validate(bsize, isize)) {
            return false;
        }
        p += bsize;
    }
    return p == buffer + bytes;
}

    bool
BgzfHeader::validate(
    size_t compressed,
    size_t uncompressed)
{
    return ID1 == 0x1f && ID2 == 0x8b && CM == 8 && FLG == 4 &&
        MTIME == 0 && XFL == 0 && OS == 0 &&
        ISIZE() == uncompressed&&
        BSIZE() + 1 == compressed;
}

const int BAMReader::MAX_SEQ_LENGTH = MAX_READ_LENGTH;
const int BAMReader::MAX_RECORD_LENGTH = MAX_READ_LENGTH * 8;
