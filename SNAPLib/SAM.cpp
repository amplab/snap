/*++

Module Name:

    SAM.cpp

Abstract:

    Sequence Alignment Map (SAM) file writer and reader.

Environment:

    User mode service.

    SamWriter and SamReader (and their subclasses) aren't thread safe.

--*/

#include "stdafx.h"
#include "BigAlloc.h"
#include "Compat.h"
#include "Read.h"
#include "SAM.h"
#include "Bam.h"
#include "Tables.h"
#include "RangeSplitter.h"
#include "ParallelTask.h"
#include "Util.h"
#include "ReadSupplierQueue.h"
#include "FileFormat.h"
#include "AlignerOptions.h"
#include "directions.h"
#include "exit.h"
#include <vector>

using std::max;
using std::min;
using util::strnchr;

bool readIdsMatch(const char* id0, const char* id1)
{
    for (unsigned i = 0; ; i++) {
        char c0 = id0[i];
        char c1 = id1[i];

        if (c0 != c1) return false;
 
        // don't parse the read ID after the first space or slash, which can represent metadata (or which half of the mate pair the read is).
        if (c0 == 0 || c0 == ' ' || c0 == '/') return true;  
    }
    return true;
}

bool readIdsMatch(Read *read0, Read *read1)
{
    if (read0->getIdLength() != read1->getIdLength()) {
        return false;
    }
    for (unsigned i = 0; i < read0->getIdLength(); i++) {
        char c0 = read0->getId()[i];
        char c1 = read1->getId()[i];

        if (c0 != c1) return false;
 
        // don't parse the read ID after the first space or slash, which can represent metadata (or which half of the mate pair the read is).
        if (c0 == ' ' || c0 == '/') return true;  
    }
    return true;
}

    char *
strnchrs(char *str, char charToFind, char charToFind2, size_t maxLen) // Hokey version that looks for either of two chars
{
    for (size_t i = 0; i < maxLen; i++) {
        if (str[i] == charToFind || str[i] == charToFind2) {
            return str + i;
        }
        if (str[i] == 0) {
            return NULL;
        }
    }
    return NULL;
}

    char *
SAMReader::skipToBeyondNextFieldSeparator(char *str, const char *endOfBuffer, size_t *o_charsUntilFirstSeparator)
{
    if (NULL == str) return NULL;

    char *nextChar = str;
    while (nextChar < endOfBuffer && *nextChar != '\n' && *nextChar != '\t' && *nextChar != '\r' /* for Windows CRLF text */) {
        nextChar++;
    }

    if (NULL != o_charsUntilFirstSeparator) {
        *o_charsUntilFirstSeparator = nextChar - str;
    }

    if (nextChar >= endOfBuffer || *nextChar == '\n') {
        return NULL;
    }

    while (nextChar < endOfBuffer && ('\t' == *nextChar || '\r' == *nextChar)) {
        nextChar++;
    }

    if (nextChar >= endOfBuffer) {
        return NULL;
    }

    return nextChar;
}


    SAMReader *
SAMReader::create(
    DataSupplier* supplier,
    const char *fileName,
    int bufferCount,
    const ReaderContext& context,
    _int64 startingOffset, 
    _int64 amountOfFileToProcess)
{
    DataReader* data = supplier->getDataReader(bufferCount, maxLineLen, 0.0);
    SAMReader *reader = new SAMReader(data, context);
    reader->init(fileName, startingOffset, amountOfFileToProcess);
    return reader;
}

    void
SAMReader::readHeader(const char *fileName)
{
    // todo: allow for larger headers
    _int64 headerSize = 1024 * 1024; // 1M header max
    char* buffer = data->readHeader(&headerSize);
    if (!parseHeader(fileName, buffer, buffer + headerSize, context.genome, &headerSize, &context.headerMatchesIndex)) {
        WriteErrorMessage("SAMReader: failed to parse header on '%s'\n",fileName);
        soft_exit(1);
    }
    _ASSERT(context.header == NULL);
    char* p = new char[headerSize + 1];
    memcpy(p, buffer, headerSize);
    p[headerSize] = 0;
    context.header = p;
    context.headerBytes = context.headerLength = headerSize;
}

SAMReader::SAMReader(
    DataReader* i_data,
    const ReaderContext& i_context)
    : ReadReader(i_context), data(i_data), headerSize(-1), clipping(i_context.clipping)
{
}


//
// Implement the ReadReader form of getNextRead, which doesn't include the
// alignment results by simply throwing them away.
//
    bool
SAMReader::getNextRead(Read *readToUpdate)
{
    return getNextRead(readToUpdate, NULL, NULL, NULL, NULL, NULL, NULL);
}

    bool
SAMReader::parseHeader(
    const char *fileName, 
    char *firstLine, 
    char *endOfBuffer, 
    const Genome *genome, 
    _int64 *o_headerSize,
    bool *o_headerMatchesIndex)
{
    char *nextLineToProcess = firstLine;
    *o_headerMatchesIndex = true;
    int numSQLines = 0;
    while (NULL != nextLineToProcess && nextLineToProcess < endOfBuffer && '@' == *nextLineToProcess) {
        if (!strncmp("@SQ",nextLineToProcess,3)) {
            //
            // These lines represent sequences in the reference genome, what are
            // called "contigs" in the Genome class.  (Roughly, chromosomes or major
            // variants like some versions of the MHC genes on chr6; or more
            // particularly the things that come in different FASTA files from the
            // reference assembly).
            //
            // Verify that they actually match what's in our reference genome.
            //
            numSQLines++;
            if (nextLineToProcess + 3 >= endOfBuffer || ' ' != nextLineToProcess[3] && '\t' != nextLineToProcess[3]) {
                WriteErrorMessage("Malformed SAM file '%s' has @SQ without a following space or tab.\n",fileName);
                return false;
            }

            char *snStart = nextLineToProcess + 4;
            while (snStart < endOfBuffer && strncmp(snStart,"SN:",__min(3,endOfBuffer-snStart)) && *snStart != '\n' && *snStart != 0) {
                snStart++;
            }

            if (snStart >= endOfBuffer || *snStart == '\n' || *snStart == 0) {
                WriteErrorMessage("Malformed @SQ line doesn't have 'SN:' in file '%s'\n",fileName);
                return false;
            }

            const size_t contigNameBufferSize = 512;
            char contigName[contigNameBufferSize];
            for (unsigned i = 0; i < contigNameBufferSize && snStart+3+i < endOfBuffer; i++) {
                if (snStart[3+i] == ' ' || snStart[3+i] == '\t' || snStart[3+i] == '\n' || snStart[3+i] == 0) {
                    contigName[i] = '\0';
                } else {
                    contigName[i] = snStart[3+i];
                }
            }
            contigName[contigNameBufferSize - 1] = '\0';

            if (genome == NULL || !genome->getOffsetOfContig(contigName,NULL)) {
                *o_headerMatchesIndex = false;
            }
        } else if (!strncmp("@HD",nextLineToProcess,3) || !strncmp("@RG",nextLineToProcess,3) || !strncmp("@PG",nextLineToProcess,3) ||
            !strncmp("@CO",nextLineToProcess,3)) {
            //
            // Ignore these lines.
            //
        } else {
            WriteErrorMessage("Unrecognized header line in SAM file.\n");
            return false;
        }
		char * p = strnchr(nextLineToProcess,'\n',endOfBuffer-nextLineToProcess);
		if (p == NULL) {
            // no newline, look for null to truncate buffer
            p = (char*) memchr(nextLineToProcess, 0, endOfBuffer - nextLineToProcess);
            nextLineToProcess = p != NULL ? p + 1 : endOfBuffer;
            break;
		}
        nextLineToProcess = p + 1;
    }

    *o_headerMatchesIndex &= genome != NULL && numSQLines == genome->getNumContigs();
	*o_headerSize = nextLineToProcess - firstLine;
    return true;
}

    bool
SAMReader::parseLine(char *line, char *endOfBuffer, char *result[], size_t *linelength, size_t fieldLengths[])
{
    *linelength = 0;

    char *next = line;
    char *endOfLine = strnchr(line,'\n',endOfBuffer-line);
    if (NULL == endOfLine) {
        return false;
    }

    //
    // Skip over any leading spaces and tabs
    //
    while (next < endOfLine && (*next == ' ' || *next == '\t')) {
        next++;
    }

    for (unsigned i = 0; i < nSAMFields; i++) {
        if (NULL == next || next >= endOfLine) {
            if (i == OPT) {
                // no optional fields
                result[OPT] = NULL;
                break;
            } else {
                //
                // Too few fields.
                //
                return false;
            }
        }

        result[i] = next;
        if (i == OPT) {
            // OPT field is actually all fields until end of line
            fieldLengths[OPT] = endOfLine - next;
            break;
        }

        next = skipToBeyondNextFieldSeparator(next,endOfLine,&fieldLengths[i]);
    }

    *linelength =  endOfLine - line + 1;    // +1 skips over the \n
    return true;
}

    void
SAMReader::getReadFromLine(
    const Genome        *genome,
    char                *line, 
    char                *endOfBuffer, 
    Read                *read, 
    AlignmentResult     *alignmentResult,
    unsigned            *out_genomeLocation, 
    Direction           *direction,
    unsigned            *mapQ,
    size_t              *lineLength,
    unsigned *           flag,
    const char **        cigar,
    ReadClippingType     clipping
    )
{
    char *field[nSAMFields];
    size_t fieldLength[nSAMFields];

    if (!parseLine(line, endOfBuffer, field, lineLength, fieldLength)) {
        WriteErrorMessage( "Failed to parse SAM line:\n%.*s\n", lineLength, line);
        soft_exit(1);
    }

    //
    // We have to copy the contig name (RNAME) into its own buffer because the code in Genome expects
    // it to be a null-terminated string, while all we've got is one that's space delimited.
    //
    const size_t contigNameBufferSize = 512;
    char contigName[contigNameBufferSize];
    unsigned offsetOfContig;
    parseContigName(genome, contigName, contigNameBufferSize, &offsetOfContig, NULL, field, fieldLength);

    unsigned genomeLocation = parseLocation(offsetOfContig, field, fieldLength);

    if (NULL != out_genomeLocation) {
        *out_genomeLocation = genomeLocation;
    }

    if (fieldLength[SEQ] != fieldLength[QUAL]) {
        WriteErrorMessage("SAMReader: QUAL string unequal in length to SEQ string.\n");
        soft_exit(1);
    }

    unsigned _flag;
    const size_t flagBufferSize = 20;   // More than enough
    char flagBuffer[flagBufferSize];
    if (fieldLength[FLAG] >= flagBufferSize) {
        WriteErrorMessage("SAMReader: flag field is too long.\n");
        soft_exit(1);
    }
    memcpy(flagBuffer,field[FLAG],fieldLength[FLAG]);
    flagBuffer[fieldLength[FLAG]] = '\0';
    if (1 != sscanf(flagBuffer,"%d",&_flag)) {
        WriteErrorMessage("SAMReader: couldn't parse FLAG field.\n");
        soft_exit(1);
    }

    if (NULL != read) {
        //
        // Clip reads where the quality strings end in '#'
        //

        unsigned originalFrontClipping, originalBackClipping, originalFrontHardClipping, originalBackHardClipping;
        Read::computeClippingFromCigar(field[CIGAR], &originalFrontClipping, &originalBackClipping, &originalFrontHardClipping, &originalBackHardClipping);

        unsigned pnext = atoi(field[PNEXT]);    // Relies on atoi() returning 0 for non-numeric fields (i.e., *)

        read->init(field[QNAME],(unsigned)fieldLength[QNAME],field[SEQ],field[QUAL],(unsigned)fieldLength[SEQ], genomeLocation, atoi(field[MAPQ]), _flag, 
            originalFrontClipping, originalBackClipping, originalFrontHardClipping, originalBackHardClipping, field[RNEXT], (unsigned)fieldLength[RNEXT], pnext);
        //
        // If this read is RC in the SAM file, we need to reverse it here, since Reads are always the sense that they were as they came
        // out of the base caller.
        //

        if (_flag & SAM_REVERSE_COMPLEMENT) {
            read->becomeRC();
        }
        read->clip(clipping);

        if (field[OPT] != NULL) {
            unsigned n = (unsigned) fieldLength[OPT];
            while (n > 0 && (field[OPT][n-1] == '\n' || field[OPT][n-1] == '\r')) {
                n--;
            }
            read->setAuxiliaryData(field[OPT], n);
            for (char* p = field[OPT]; p != NULL && p < field[OPT] + fieldLength[OPT]; p = SAMReader::skipToBeyondNextFieldSeparator(p, field[OPT] + fieldLength[OPT])) {
                if (strncmp(p, "RG:Z:", 5) == 0) {
                    read->setReadGroup(READ_GROUP_FROM_AUX);
                    break;
                }
            }
        }
    }

    if (NULL != alignmentResult) {
        if (_flag & SAM_UNMAPPED) {
            *alignmentResult = NotFound;
        } else {
            if ('*' == contigName[0]) {
                WriteErrorMessage("SAMReader: mapped read didn't have RNAME filled in.\n");
                soft_exit(1);
            }
            *alignmentResult = SingleHit;   // NB: This isn't quite right, we should look at MAPQ.
        }
    }

    if (NULL != direction) {
        *direction = (_flag & SAM_REVERSE_COMPLEMENT) ? RC : FORWARD;
    }

    if (NULL != mapQ) {
        *mapQ = atoi(field[MAPQ]);
        if (*mapQ > 255) {
            WriteErrorMessage("SAMReader: MAPQ field has bogus value\n");
            soft_exit(1);
        }
    }

    if (NULL != flag) {
        *flag = _flag;
    }

    if (NULL != cigar) {
        *cigar = field[CIGAR];
    }
}

    void
SAMReader::parseContigName(
    const Genome* genome,
    char* contigName,
    size_t contigNameBufferSize,
    unsigned* o_offsetOfContig,
	int* o_indexOfContig,
    char* field[],
    size_t fieldLength[],
	unsigned rfield)
{
    if (fieldLength[rfield] >= contigNameBufferSize) {  // >= because we need a byte for the \0
        WriteErrorMessage("SAMReader: too long an RNAME.  Can't parse.\n");
        soft_exit(1);
    }
    
    memcpy(contigName,field[rfield],fieldLength[rfield]);
    contigName[fieldLength[rfield]] = '\0';

    *o_offsetOfContig = 0;
    if ('*' != contigName[0] && genome != NULL && !genome->getOffsetOfContig(contigName,o_offsetOfContig, o_indexOfContig)) {
        //WriteErrorMessage("Unable to find contig '%s' in genome.  SAM file malformed.\n",contigName);
        //soft_exit(1);
    }
}

    unsigned
SAMReader::parseLocation(
    unsigned offsetOfContig,
    char* field[],
    size_t fieldLength[],
	unsigned rfield,
	unsigned posfield)
{
    unsigned oneBasedOffsetWithinContig = 0;
    if ('*' != field[rfield][0] && '*' != field[posfield][0]) {
        //
        // We can't call sscanf directly into the mapped file, becuase it reads to the end of the
        // string even when it's satisfied all of its fields.  Since this can be gigabytes, it's not
        // really good for perf.  Instead, copy the POS field into a local buffer and null terminate it.
        //

        const unsigned posBufferSize = 20;
        char posBuffer[posBufferSize];
        if (fieldLength[posfield] >= posBufferSize) {
            WriteErrorMessage("SAMReader: POS field too long.\n");
            soft_exit(1);
        }
        memcpy(posBuffer,field[posfield],fieldLength[posfield]);
        posBuffer[fieldLength[posfield]] = '\0';
        if (0 == sscanf(posBuffer,"%d",&oneBasedOffsetWithinContig)) {
            WriteErrorMessage("SAMReader: Unable to parse position when it was expected.\n");
            soft_exit(1);
        }
        if (0 == oneBasedOffsetWithinContig) {
            WriteErrorMessage("SAMReader: Position parsed as 0 when it was expected.\n");
            soft_exit(1);
        }
        return offsetOfContig + oneBasedOffsetWithinContig - 1; // -1 is because our offset is 0 based, while SAM is 1 based.
    } else {
        return InvalidGenomeLocation;
    }
}
    
    void
SAMReader::init(
    const char *fileName,
    _int64 startingOffset,
    _int64 amountOfFileToProcess)
{
    if (! data->init(fileName)) {
        WriteErrorMessage( "Unable to read file %s\n", fileName);
        soft_exit(1);
    }

    if (0 == startingOffset) {
        readHeader(fileName);
    }

    headerSize = context.headerBytes;
    reinit(max(startingOffset, (_int64) context.headerBytes),
        amountOfFileToProcess == 0 || startingOffset >= (_int64) context.headerBytes ? amountOfFileToProcess
            : amountOfFileToProcess - (context.headerBytes - startingOffset));
}

    void
SAMReader::reinit(_int64 startingOffset, _int64 amountOfFileToProcess)
{
    _ASSERT(-1 != headerSize && startingOffset >= headerSize);  // Must call init() before reinit()
    //
    // There's no way to tell if we start at the very beginning of a read, we need to see the previous newline.
    // So, read one byte before our assigned read in case that was the terminating newline of the previous read.
    //
    if (startingOffset > headerSize) {
        startingOffset--;
        amountOfFileToProcess++;
    }
    data->reinit(startingOffset, amountOfFileToProcess);
    char* buffer;
    _int64 validBytes;
    if (!data->getData(&buffer, &validBytes)) {
        return;
    }
    if (startingOffset != headerSize) {
        char *firstNewline = strnchr(buffer,'\n',validBytes);
        if (NULL == firstNewline) {
            return;
        }

        data->advance((unsigned)(firstNewline - buffer + 1)); // +1 skips over the newline.
    }
}

    bool
SAMReader::getNextRead(
    Read *read,
    AlignmentResult *alignmentResult,
    unsigned *genomeLocation,
    Direction *direction,
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
        }
        char *newLine = strnchr(buffer, '\n', bytes);
        if (NULL == newLine) {
            //
            // There is no newline, so the line crosses the end of the buffer.
            // This should never happen since underlying reader manages overflow between chunks.
            //
            WriteErrorMessage("SAM file has too long a line, or doesn't end with a newline!  Failing.  fileOffset = %lld\n", data->getFileOffset());
            soft_exit(1);
        }

        size_t lineLength;
        read->setReadGroup(context.defaultReadGroup);
        getReadFromLine(context.genome, buffer,buffer + bytes, read, alignmentResult, genomeLocation, direction, mapQ, &lineLength, flag, cigar, clipping);
        read->setBatch(data->getBatch());
        data->advance((newLine + 1) - buffer);
    } while ((context.ignoreSecondaryAlignments && ((*flag) & SAM_SECONDARY)) ||
             (context.ignoreSupplementaryAlignments && ((*flag) & SAM_SUPPLEMENTARY)));

    return true;
}

    ReadSupplierGenerator *
SAMReader::createReadSupplierGenerator(
    const char *fileName,
    int numThreads,
    const ReaderContext& context)
{
    //
    // single-ended SAM files always can be read with the range splitter.
    //
    if (!strcmp(fileName, "-")) {
        //
        // Stdin must run from a queue, not range splitter.
        //
        ReadReader* reader;
        //
        // Because we can only have one stdin reader, we need to use a queue if we're reading from stdin
        //
        reader = SAMReader::create(DataSupplier::Stdio, "-", ReadSupplierQueue::BufferCount(numThreads), context, 0, 0);
   
        if (reader == NULL) {
            return NULL;
        }
        ReadSupplierQueue *queue = new ReadSupplierQueue(reader);
        queue->startReaders();
        return queue;
    } else {
        RangeSplitter *splitter = new RangeSplitter(QueryFileSize(fileName), numThreads, 100);
        return new RangeSplittingReadSupplierGenerator(fileName, true, numThreads, context);
    }
}
    
    PairedReadReader*
SAMReader::createPairedReader(
    const DataSupplier* supplier,
    const char *fileName,
    int bufferCount,
    _int64 startingOffset,
    _int64 amountOfFileToProcess, 
    bool quicklyDropUnpairedReads,
    const ReaderContext& context)
{
    DataSupplier *data;
    if (!strcmp("-", fileName)) {
        data = DataSupplier::Stdio;
    } else {
        data = DataSupplier::Default;
    }

    SAMReader* reader = SAMReader::create(data, fileName, bufferCount + PairedReadReader::MatchBuffers, context, 0, 0);
    if (reader == NULL) {
        return NULL;
    }
    return PairedReadReader::PairMatcher(reader, quicklyDropUnpairedReads);
}


    PairedReadSupplierGenerator *
SAMReader::createPairedReadSupplierGenerator(
    const char *fileName,
    int numThreads,
    bool quicklyDropUnpairedReads, 
    const ReaderContext& context)
{
    //
    // need to use a queue so that pairs can be matched
    //

    PairedReadReader* paired = SAMReader::createPairedReader(DataSupplier::Default, fileName,
        ReadSupplierQueue::BufferCount(numThreads), 0, 0, quicklyDropUnpairedReads, context);
    if (paired == NULL) {
        WriteErrorMessage( "Cannot create reader on %s\n", fileName);
        soft_exit(1);
    }
    ReadSupplierQueue* queue = new ReadSupplierQueue(paired);
    queue->startReaders();
    return queue;
}


const FileFormat* FileFormat::SAM[] = { new SAMFormat(false), new SAMFormat(true) };

    void
SAMFormat::getSortInfo(
    const Genome* genome,
    char* buffer,
    _int64 bytes,
	unsigned* o_location,
	unsigned* o_readBytes,
	int* o_refID,
	int* o_pos) const
{
    char* fields[SAMReader::nSAMFields];
    size_t lengths[SAMReader::nSAMFields];
    size_t lineLength;
    SAMReader::parseLine(buffer, buffer + bytes, fields, &lineLength, lengths);
    _ASSERT(lineLength < UINT32_MAX);
	if (o_readBytes != NULL) {
		*o_readBytes = (unsigned) lineLength;
	}
    if (lengths[SAMReader::POS] == 0 || fields[SAMReader::POS][0] == '*') {
		if (lengths[SAMReader::PNEXT] == 0 || fields[SAMReader::PNEXT][0] == '*') {
			if (o_location != NULL) {
				*o_location = UINT32_MAX;
			}
			if (o_refID != NULL) {
				*o_refID = -1;
			}
			if (o_pos != NULL) {
				*o_pos = 0;
			}
		} else {
			const size_t contigNameBufferSize = 512;
			char contigName[contigNameBufferSize];
			unsigned offsetOfContig;
			SAMReader::parseContigName(genome, contigName, contigNameBufferSize, &offsetOfContig, o_refID, fields, lengths, SAMReader::RNEXT);
			if (o_location != NULL) {
				*o_location = SAMReader::parseLocation(offsetOfContig, fields, lengths, SAMReader::RNEXT, SAMReader::PNEXT);
			}
		}
    } else {
        const size_t contigNameBufferSize = 512;
        char contigName[contigNameBufferSize];
        unsigned offsetOfContig;
        SAMReader::parseContigName(genome, contigName, contigNameBufferSize, &offsetOfContig, o_refID, fields, lengths);
		if (o_location != NULL) {
	        *o_location = SAMReader::parseLocation(offsetOfContig, fields, lengths);
		}
    }
}

    ReadWriterSupplier*
SAMFormat::getWriterSupplier(
    AlignerOptions* options,
    const Genome* genome,
    const Genome* transcriptome,
    const GTFReader* gtf) const
{
    DataWriterSupplier* dataSupplier;
    if (options->sortOutput) {
        size_t len = strlen(options->outputFile.fileName);
        // todo: this is going to leak, but there's no easy way to free it, and it's small...
        char* tempFileName = (char*) malloc(5 + len);
        strcpy(tempFileName, options->outputFile.fileName);
        strcpy(tempFileName + len, ".tmp");
        dataSupplier = DataWriterSupplier::sorted(this, genome, tempFileName, options->sortMemory * (1ULL << 30),
            options->numThreads, options->outputFile.fileName, NULL);
    } else {
        dataSupplier = DataWriterSupplier::create(options->outputFile.fileName);
    }
    return ReadWriterSupplier::create(this, dataSupplier, genome, transcriptome, gtf);
}

    bool
SAMFormat::writeHeader(
    const ReaderContext& context,
    char *header,
    size_t headerBufferSize,
    size_t *headerActualSize,
    bool sorted,
    int argc,
    const char **argv,
    const char *version,
    const char *rgLine)
    const
{
    char *commandLine;
	size_t commandLineSize = 0;
	for (int i = 0; i < argc; i++) {
		commandLineSize += strlen(argv[i]) + 1;	// +1 is either a space or the terminating null
	}
	commandLine = new char[commandLineSize];
	commandLine[0] = '\0';
	for (int i = 0; i < argc; i++) {
		strcat(commandLine,argv[i]);
		if (i != argc-1) {
			strcat(commandLine," ");
		}
	}

    size_t bytesConsumed = snprintf(header, headerBufferSize, "@HD\tVN:1.4\tSO:%s\n%s%s@PG\tID:SNAP\tPN:SNAP\tCL:%s\tVN:%s\n", 
		sorted ? "coordinate" : "unsorted",
        context.header == NULL ? (rgLine == NULL ? "@RG\tID:FASTQ\tSM:sample" : rgLine) : "",
        context.header == NULL ? "\n" : "",
        commandLine,version);

	delete [] commandLine;
	commandLine = NULL;
    if (bytesConsumed >= headerBufferSize) {
        WriteErrorMessage("SAMWriter: header buffer too small\n");
        return false;
    }

    if (context.header != NULL) {
		bool hasRG = false;
        for (const char* p = context.header; p < context.header + context.headerLength; ) {
            const char* newline = strnchr(p, '\n', (context.header + context.headerLength) - p);
            if (newline == NULL) {
                newline = context.header + context.headerLength;
            }
            _ASSERT(newline - p >= 3);
            // skip @HD lines, and also @SQ lines if header does not match index
			hasRG |= strncmp(p, "@RG", 3) == 0;
            if (strncmp(p, "@HD", 3) != 0 &&
                    (context.headerMatchesIndex || strncmp(p, "@SQ", 3) != 0) &&
                    strncmp(p, "@PG\tID:SNAP\t", 12) != 0) {
                if (bytesConsumed + (newline - p) + 1 >= headerBufferSize) {
                    WriteErrorMessage("SAMWriter: header buffer too small\n");
                    return false;
                }
                memcpy(header + bytesConsumed, p, (newline - p));
                * (header + bytesConsumed + (newline - p)) = '\n';
                bytesConsumed += (newline - p) + 1;
            }
            p = newline + 1;
        }
		if (! hasRG) {
			int n = snprintf(header + bytesConsumed, headerBufferSize - bytesConsumed, "%s\n",
				rgLine == NULL ? "@RG\tID:FASTQ\tSM:sample" : rgLine);
			if (n > headerBufferSize - bytesConsumed) {
				WriteErrorMessage( "SAMWriter: header buffer too small\n");
                return false;
            }
			bytesConsumed += n;
		}
    }
#ifndef SKIP_SQ_LINES
    if ((context.header == NULL || ! context.headerMatchesIndex) && context.genome != NULL) {
        // Write an @SQ line for each chromosome / contig in the genome
        const Genome::Contig *contigs = context.genome->getContigs();
        int numContigs = context.genome->getNumContigs();
        unsigned genomeLen = context.genome->getCountOfBases();
        size_t originalBytesConsumed = bytesConsumed;
        for (int i = 0; i < numContigs; i++) {
            unsigned start = contigs[i].beginningOffset;
            unsigned end = ((i + 1 < numContigs) ? contigs[i+1].beginningOffset : genomeLen) - context.genome->getChromosomePadding();
            bytesConsumed += snprintf(header + bytesConsumed, headerBufferSize - bytesConsumed, "@SQ\tSN:%s\tLN:%u\n", contigs[i].name, end - start);

            if (bytesConsumed >= headerBufferSize) {
                // todo: increase buffer size (or change to write in batch
                bytesConsumed = originalBytesConsumed;
                WriteErrorMessage("SAMWriter: header buffer too small, skipping @SQ lines\n");
                break;
            }
        }
    }
#endif // SKIP_SQ_LINES

    *headerActualSize = bytesConsumed;
    return true;
}
    
    bool
SAMFormat::createSAMLine(
    const Genome * genome,
    const Genome * transcriptome,
    const GTFReader * gtf,
    LandauVishkinWithCigar * lv,
    // output data
    char* data,
    char* quality,
    unsigned dataSize,
    const char*& contigName,
    int& contigIndex,
    int& flags,
    unsigned& positionInContig,
    int& mapQuality,
    const char*& matecontigName,
    int& mateContigIndex,
    unsigned& matePositionInContig,
    _int64& templateLength,
    unsigned& fullLength,
    const char*& clippedData,
    unsigned& clippedLength,
    unsigned& basesClippedBefore,
    unsigned& basesClippedAfter,
    // input data
    size_t& qnameLen,
    Read * read,
    AlignmentResult result, 
    unsigned genomeLocation,
    Direction direction,
    bool secondaryAlignment,
    bool isTranscriptome,
    bool useM,
    bool hasMate,
    bool firstInPair,
    Read * mate, 
    AlignmentResult mateResult,
    unsigned mateLocation,
    Direction mateDirection,
    bool mateIsTranscriptome,
    unsigned *extraBasesClippedBefore,
    unsigned *extraBasesClippedAfter)
{
    contigName = "*";
    positionInContig = 0;
    const char *cigar = "*";
    templateLength = 0;

    if (secondaryAlignment) {
        flags |= SAM_SECONDARY;
    }
    
    if (0 == qnameLen) {
         qnameLen = read->getIdLength();
    }

    //
    // If the aligner said it didn't find anything, treat it as such.  Sometimes it will emit the
    // best match that it found, even if it's not within the maximum edit distance limit (but will
    // then say NotFound).  Here, we force that to be SAM_UNMAPPED.
    //
    if (NotFound == result) {
        genomeLocation = InvalidGenomeLocation;
    }

    if (InvalidGenomeLocation == genomeLocation) {
        //
        // If it's unmapped, then always emit it in the forward direction.  This is necessary because we don't even include
        // the SAM_REVERSE_COMPLEMENT flag for unmapped reads, so there's no way to tell that we reversed it.
        //
        direction = FORWARD;
    }

    // Write the data and quality strings. If the read is reverse complemented, these need to
    // be backwards from the original read. Also, both need to be unclipped.
    clippedLength = read->getDataLength();
    fullLength = read->getUnclippedLength();
    if (fullLength > dataSize) {
        return false;
    }

    if (direction == RC) {
      for (unsigned i = 0; i < fullLength; i++) {
        data[fullLength - 1 - i] = COMPLEMENT[read->getUnclippedData()[i]];
        quality[fullLength - 1 - i] = read->getUnclippedQuality()[i];
      }
      clippedData = &data[fullLength - clippedLength - read->getFrontClippedLength()];
      basesClippedBefore = fullLength - clippedLength - read->getFrontClippedLength();
      basesClippedAfter = read->getFrontClippedLength();
    } else {
      memcpy(data, read->getUnclippedData(), read->getUnclippedLength());
      memcpy(quality, read->getUnclippedQuality(), read->getUnclippedLength());
      clippedData = read->getData();
      basesClippedBefore = read->getFrontClippedLength();
      basesClippedAfter = fullLength - clippedLength - basesClippedBefore;
    }

    int editDistance = -1;
    *extraBasesClippedAfter = 0;
    if (genomeLocation != InvalidGenomeLocation) {
        // This could be either a single hit read or a multiple hit read where we just
        // returned one location, but either way, let's print that location. We will then
        // set the quality to 60 if it was single or 0 if it was multiple. These are the
        // values the SAM FAQ suggests for aligners that don't compute confidence scores.
        if (direction == RC) {
            flags |= SAM_REVERSE_COMPLEMENT;
        }
        const Genome::Contig *contig = genome->getContigForRead(genomeLocation, read->getDataLength(), extraBasesClippedBefore);
        _ASSERT(NULL != contig && contig->length > genome->getChromosomePadding());
        if (genomeLocation + read->getDataLength() > contig->beginningOffset + contig->length - genome->getChromosomePadding()) {
            //
            // The read hangs off the end of the contig.  Soft clip it at the end.
            //
            *extraBasesClippedAfter =genomeLocation + read->getDataLength() - (contig->beginningOffset + contig->length - genome->getChromosomePadding());
        }
        genomeLocation += *extraBasesClippedBefore;

        contigName = contig->name;
        contigIndex = (int)(contig - genome->getContigs());
        positionInContig = genomeLocation - contig->beginningOffset + 1; // SAM is 1-based
        mapQuality = max(0, min(70, mapQuality));       // FIXME: manifest constant.
    } else {
        flags |= SAM_UNMAPPED;
        mapQuality = 0;
        *extraBasesClippedBefore = 0;
    }

    if (hasMate) {
        flags |= SAM_MULTI_SEGMENT;
        flags |= (firstInPair ? SAM_FIRST_SEGMENT : SAM_LAST_SEGMENT);
        if (mateLocation != InvalidGenomeLocation) {
            unsigned mateExtraBasesClippedBefore;
            const Genome::Contig *mateContig = genome->getContigForRead(mateLocation, mate->getDataLength(), &mateExtraBasesClippedBefore);
            mateLocation += mateExtraBasesClippedBefore;
            matecontigName = mateContig->name;
            mateContigIndex = (int)(mateContig - genome->getContigs());
            matePositionInContig = mateLocation - mateContig->beginningOffset + 1;

            if (mateDirection == RC) {
                flags |= SAM_NEXT_REVERSED;
            }

            if (genomeLocation == InvalidGenomeLocation) {
                //
                // The SAM spec says that for paired reads where exactly one end is unmapped that the unmapped
                // half should just have RNAME and POS copied from the mate.
                //
                contigName = matecontigName;
                contigIndex = mateContigIndex;
                matecontigName = "=";
                positionInContig = matePositionInContig;
            }

        } else {
            flags |= SAM_NEXT_UNMAPPED;
            //
            // The mate's unmapped, so point it at us.
            //
            matecontigName = "=";
            mateContigIndex = contigIndex;
            matePositionInContig = positionInContig;
        }

        if (genomeLocation != InvalidGenomeLocation && mateLocation != InvalidGenomeLocation) {
            flags |= SAM_ALL_ALIGNED;
            // Also compute the length of the whole paired-end string whose ends we saw. This is slightly
            // tricky because (a) we may have clipped some bases before/after each end and (b) we need to
            // give a signed result based on whether our read is first or second in the pair.
            _int64 myStart = genomeLocation - basesClippedBefore;
            _int64 myEnd = genomeLocation + clippedLength + basesClippedAfter;
            _int64 mateBasesClippedBefore = mate->getFrontClippedLength();
            _int64 mateBasesClippedAfter = mate->getUnclippedLength() - mate->getDataLength() - mateBasesClippedBefore;
            _int64 mateStart = mateLocation - (mateDirection == RC ? mateBasesClippedAfter : mateBasesClippedBefore);
            _int64 mateEnd = mateLocation + mate->getDataLength() + (mateDirection == FORWARD ? mateBasesClippedAfter : mateBasesClippedBefore);
			if (contigName == matecontigName) { // pointer (not value) comparison, but that's OK.
				if (myStart < mateStart) {
					templateLength = mateEnd - myStart;
				} else {
					templateLength = -(myEnd - mateStart);
				}
 			} // otherwise leave TLEN as zero.
        }

        if (contigName == matecontigName) {
            matecontigName = "=";     // SAM Spec says to do this when they're equal (and not *, which won't happen because this is a pointer, not string, compare)
        }
    }
    return true;
}

    bool
SAMFormat::writeRead(
    const Genome * genome,
    const Genome * transcriptome,
    const GTFReader *gtf,
    LandauVishkinWithCigar * lv,
    char * buffer,
    size_t bufferSpace, 
    size_t * spaceUsed,
    size_t qnameLen,
    Read * read,
    AlignmentResult result, 
    int mapQuality,
    unsigned genomeLocation,
    Direction direction,
    bool secondaryAlignment,
    bool isTranscriptome,
    unsigned tlocation,
    bool hasMate,
    bool firstInPair,
    Read * mate, 
    AlignmentResult mateResult,
    unsigned mateLocation,
    Direction mateDirection,
    bool mateIsTranscriptome,
    unsigned mateTlocation
    ) const
{
    const int MAX_READ = MAX_READ_LENGTH;
    const int cigarBufSize = MAX_READ * 2;
    char cigarBuf[cigarBufSize];

    const int cigarBufWithClippingSize = MAX_READ * 2 + 32;
    char cigarBufWithClipping[cigarBufWithClippingSize];

    int flags = 0;
    const char *contigName = "*";
    int contigIndex = -1;
    unsigned positionInContig = 0;
    const char *cigar = "*";
    const char *matecontigName = "*";
    int mateContigIndex = -1;
    unsigned matePositionInContig = 0;
    _int64 templateLength = 0;

    char data[MAX_READ];
    char quality[MAX_READ];

    const char* clippedData;
    unsigned fullLength;
    unsigned clippedLength;
    unsigned basesClippedBefore;
    unsigned extraBasesClippedBefore;   // Clipping added if we align before the beginning of a chromosome
    unsigned basesClippedAfter;
    unsigned extraBasesClippedAfter;    // Clipping added if we align off the end of a chromosome
    int editDistance = -1;

    if (! createSAMLine(genome, transcriptome, gtf, lv, data, quality, MAX_READ, contigName, contigIndex, 
        flags, positionInContig, mapQuality, matecontigName, mateContigIndex, matePositionInContig, templateLength,
        fullLength, clippedData, clippedLength, basesClippedBefore, basesClippedAfter,
        qnameLen, read, result, genomeLocation, direction, secondaryAlignment, isTranscriptome, useM,
        hasMate, firstInPair, mate, mateResult, mateLocation, mateDirection, mateIsTranscriptome, 
        &extraBasesClippedBefore, &extraBasesClippedAfter))
    {
        return false;
    }
    if (genomeLocation != InvalidGenomeLocation) {
    
        std::vector<unsigned> tokens;
        if (!isTranscriptome) {
        
            cigar = computeCigarString(genome, lv, cigarBuf, cigarBufSize, cigarBufWithClipping, cigarBufWithClippingSize, 
                                       clippedData, clippedLength, basesClippedBefore, extraBasesClippedBefore, basesClippedAfter, extraBasesClippedAfter,
                                       read->getOriginalFrontHardClipping(), read->getOriginalBackHardClipping(), genomeLocation, direction, useM, &editDistance, tokens);
                                       
        } else {
        
            cigar = computeCigarString(transcriptome, lv, cigarBuf, cigarBufSize, cigarBufWithClipping, cigarBufWithClippingSize, 
                                       clippedData, clippedLength, basesClippedBefore, extraBasesClippedBefore, basesClippedAfter, extraBasesClippedAfter,
                                       read->getOriginalFrontHardClipping(), read->getOriginalBackHardClipping(), tlocation, direction, useM, &editDistance, tokens);        
        
            //We need the pieceName for conversion             
            const Genome::Contig *transcriptomePiece = transcriptome->getContigAtLocation(tlocation);
            const char* transcriptomePieceName = transcriptomePiece->name;
            unsigned transcriptomePositionInPiece = tlocation - transcriptomePiece->beginningOffset + 1; // SAM is 1-based
 
            //Insert splice junctions
            lv->insertSpliceJunctions(gtf, tokens, transcriptomePieceName, transcriptomePositionInPiece, (char*) cigarBuf, cigarBufSize);
            cigar = cigarBuf;  
        }
    }

    // Write the SAM entry, which requires the following fields:
    //
    // 1. QNAME: Query name of the read or the read pair
    // 2. FLAG: Bitwise flag (pairing, strand, mate strand, etc.)
    // 3. RNAME: Reference sequence name
    // 4. POS: 1-Based leftmost position of clipped alignment
    // 5. MAPQ: Mapping quality (Phred-scaled)
    // 6. CIGAR: Extended CIGAR string (operations: MIDNSHP)
    // 7. MRNM: Mate reference name (‘=’ if same as RNAME)
    // 8. MPOS: 1-based leftmost mate position
    // 9. ISIZE: Inferred insert size
    // 10. SEQQuery: Sequence on the same strand as the reference
    // 11. QUAL: Query quality (ASCII-33=Phred base quality)    

    //
    // Some FASTQ files have spaces in their ID strings, which is illegal in SAM.  Just truncate them at the space.
    //
    const char *firstSpace = strnchr(read->getId(),' ',qnameLen);
    if (NULL != firstSpace) {
        qnameLen = (unsigned)(firstSpace - read->getId());
    }

    const int nmStringSize = 30;// Big enough that it won't buffer overflow regardless of the value of editDistance
    char nmString[nmStringSize];  
    snprintf(nmString, nmStringSize, "\tNM:i:%d",editDistance);

    unsigned auxLen;
    bool auxSAM;
    char* aux = read->getAuxiliaryData(&auxLen, &auxSAM);
    static bool warningPrinted = false;
    const char* readGroupSeparator = "";
    const char* readGroupString = "";
    if (aux != NULL && (! auxSAM)) {
        if (! warningPrinted) {
            WriteErrorMessage( "warning: translating optional fields from BAM->SAM not yet implemented, optional fields will not be included in output\n");
            warningPrinted = true;
        }
        if (read->getReadGroup() == READ_GROUP_FROM_AUX) {
            for (BAMAlignAux* bamAux = (BAMAlignAux*) aux; (char*) bamAux < aux + auxLen; bamAux = bamAux->next()) {
                if (bamAux->tag[0] == 'R' && bamAux->tag[1] == 'G' && bamAux->val_type == 'Z') {
                    readGroupSeparator = "\tRG:Z:";
                    readGroupString = (char*) bamAux->value();
                    break;
                }
            }
        }
        aux = NULL;
        auxLen = 0;
    }
    if (read->getReadGroup() != NULL && read->getReadGroup() != READ_GROUP_FROM_AUX) {
        readGroupSeparator = "\tRG:Z:";
        readGroupString = read->getReadGroup();
    }
    int charsInString = snprintf(buffer, bufferSpace, "%.*s\t%d\t%s\t%u\t%d\t%s\t%s\t%u\t%lld\t%.*s\t%.*s%s%.*s%s%s\tPG:Z:SNAP%s\n",
        qnameLen, read->getId(),
        flags,
        contigName,
        positionInContig,
        mapQuality,
        cigar,
        matecontigName,
        matePositionInContig,
        templateLength,
        fullLength, data,
        fullLength, quality,
        aux != NULL ? "\t" : "", auxLen, aux != NULL ? aux : "",
        readGroupSeparator, readGroupString,
        nmString);

    if (charsInString > bufferSpace) {
        //
        // Out of buffer space.
        //
        return false;
    } else if (charsInString == bufferSpace) {
      buffer[bufferSpace-1] = '\n'; // overwrite trailing null with newline
    }


    if (NULL != spaceUsed) {
        *spaceUsed = charsInString;
    }
    return true;
}

// Compute the CIGAR edit sequence string for a read against a given genome location.
// Returns this string if possible or "*" if we fail to compute it (which would likely
// be a bug due to lack of buffer space). The pointer returned may be to cigarBuf so it
// will only be valid until computeCigarString is called again.
    const char *
SAMFormat::computeCigarString(
    const Genome *              genome,
    LandauVishkinWithCigar *    lv,
    char *                      cigarBuf,
    int                         cigarBufLen,
    char *                      cigarBufWithClipping,
    int                         cigarBufWithClippingLen,
    const char *                data,
    unsigned                    dataLength,
    unsigned                    basesClippedBefore,
    unsigned                    extraBasesClippedBefore,
    unsigned                    basesClippedAfter,
    unsigned                    extraBasesClippedAfter,
    unsigned                    frontHardClipping,
    unsigned                    backHardClipping,
    unsigned                    genomeLocation,
    Direction                   direction,
	bool						useM,
    int *                       editDistance,
    std::vector<unsigned>       &tokens
)
{
    //
    // Apply the extra clipping.
    //
    genomeLocation += extraBasesClippedBefore;
    data += extraBasesClippedBefore;
    dataLength -= extraBasesClippedBefore;

    const char *reference = genome->getSubstring(genomeLocation, dataLength);
    if (NULL != reference) {
        *editDistance = lv->computeEditDistance(
                            reference,
                            dataLength - extraBasesClippedAfter + MAX_K, // Add space incase of indels.  We know there's enough, because the reference is padded.
                            data,
                            dataLength - extraBasesClippedAfter,
                            MAX_K - 1,
                            cigarBuf,
                            cigarBufLen,
						    useM,
						    tokens);
    } else {
        //
        // Fell off the end of the chromosome.
        //
        return "*";
    }

    if (*editDistance == -2) {
        WriteErrorMessage( "WARNING: computeEditDistance returned -2; cigarBuf may be too small\n");
        return "*";
    } else if (*editDistance == -1) {
        static bool warningPrinted = false;
        if (!warningPrinted) {
            WriteErrorMessage( "WARNING: computeEditDistance returned -1; this shouldn't happen\n");
            warningPrinted = true;
        }
        return "*";
    } else {
        // Add some CIGAR instructions for soft-clipping if we've ignored some bases in the read.
        char clipBefore[16] = {'\0'};
        char clipAfter[16] = {'\0'};
        char hardClipBefore[16] = {'\0'};
        char hardClipAfter[16] = {'\0'};
        if (frontHardClipping > 0) {
            snprintf(hardClipBefore, sizeof(hardClipBefore), "%uH", frontHardClipping);
            tokens.insert(tokens.begin(), 'H');
            tokens.insert(tokens.begin(), frontHardClipping);
        }
        if (basesClippedBefore + extraBasesClippedBefore > 0) {
            snprintf(clipBefore, sizeof(clipBefore), "%uS", basesClippedBefore + extraBasesClippedBefore);
            tokens.insert(tokens.begin(), 'S');
            tokens.insert(tokens.begin(), basesClippedBefore + extraBasesClippedBefore);   
        }
        if (basesClippedAfter + extraBasesClippedAfter > 0) {
            snprintf(clipAfter, sizeof(clipAfter), "%uS", basesClippedAfter + extraBasesClippedAfter);
            tokens.push_back(basesClippedAfter + extraBasesClippedAfter);
            tokens.push_back('S');
        }
        if (backHardClipping > 0) {
            snprintf(hardClipAfter, sizeof(hardClipAfter), "%uH", backHardClipping);
            tokens.push_back(backHardClipping);
            tokens.push_back('H');
        }
        snprintf(cigarBufWithClipping, cigarBufWithClippingLen, "%s%s%s%s%s", hardClipBefore, clipBefore, cigarBuf, clipAfter, hardClipAfter);

        return cigarBufWithClipping;
    }
}

    
