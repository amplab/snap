/*++

Module Name:

    FASTQ.cpp

Abstract:

    Fast FASTQ genome "query" reader.

Authors:

    Bill Bolosky, August, 2011

Environment:

    User mode service.

    This class is NOT thread safe.  It's the caller's responsibility to ensure that
    at most one thread uses an instance at any time.

Revision History:

    Adapted from Matei Zaharia's Scala implementation.

--*/

#include "stdafx.h"
#include "FASTQ.h"
#include "Compat.h"
#include "BigAlloc.h"
#include "Read.h"
#include "Util.h"
#include "exit.h"
#include "Error.h"

using std::min;
using util::strnchr;

FASTQReader::FASTQReader(
    DataReader* i_data,
    const ReaderContext& i_context)
    :
    ReadReader(i_context),
    data(i_data)
{
}

FASTQReader::~FASTQReader()
{
    delete data;
    data = NULL;
}


    FASTQReader*
FASTQReader::create(
    DataSupplier* supplier,
    const char *fileName,
    int bufferCount,
    _int64 startingOffset,
    _int64 amountOfFileToProcess,
    const ReaderContext& context)
{
    DataReader* data = supplier->getDataReader(bufferCount, maxReadSizeInBytes, 0.0, 0);
    FASTQReader* fastq = new FASTQReader(data, context);
    if (! fastq->init(fileName)) {
        WriteErrorMessage("Unable to initialize FASTQReader for file %s\n", fileName);
        soft_exit(1);
    }
    fastq->reinit(startingOffset, amountOfFileToProcess);
    return fastq;
}

    void
FASTQReader::readHeader(
    const char* fileName,
    ReaderContext& context)
{
    // no header in FQ files
    context.header = NULL;
    context.headerLength = context.headerBytes = 0;
}

    bool
FASTQReader::init(
    const char* i_fileName)
{
    fileName = i_fileName;
    return data->init(fileName);
}
    

    void
FASTQReader::reinit(
    _int64 startingOffset,
    _int64 amountOfFileToProcess)
{
    data->reinit(startingOffset, amountOfFileToProcess);
    char* buffer;
    _int64 bytes;
    if (! data->getData(&buffer, &bytes)) {
        return;
    }

    // If we're not at the start of the file, we might have the tail end of a read that someone
    // who got the previous range will process; advance past that. This is fairly tricky because
    // there can be '@' signs in the quality string (and maybe even in read names?).
    if (startingOffset != 0) {
        if (!skipPartialRecord(data)) {
            //
            // There wasn't a whole record in our range.  Skip over what we had.
            //
            data->advance(bytes);
            return;
        }
    }
}

    bool
FASTQReader::skipPartialRecord(DataReader *data)
{
    //
    // Just assume that a single FASTQ read is smaller than our buffer, so we won't exceed the buffer here.
    // Look for the pattern '{0|\n}@*\n{A|T|C|G|N}*\n+'  That is, either the beginning of the buffer or a
    // newline, followed by an '@' and some text, another newline followed by a list of bases and a newline, 
    // and then a plus.
    //
    char* buffer;
    _int64 validBytes;
    data->getData(&buffer, &validBytes);

    char *firstLineCandidate = buffer;
    if (*firstLineCandidate != '@') {
        firstLineCandidate = strnchr(buffer, '\n', validBytes) + 1;
    }

    for (;;) {
        if (firstLineCandidate - buffer >= validBytes) {
            // This happens for very small files.  
            return false;
        }

        char *secondLineCandidate = strnchr(firstLineCandidate, '\n', validBytes - (firstLineCandidate - buffer)) + 1;
        if (NULL == (secondLineCandidate-1)) {
			WriteErrorMessage("Unable to find a read in FASTQ buffer (2) at %d\n", data->getFileOffset());
            return false;
        }

        if (*firstLineCandidate != '@') {
            firstLineCandidate = secondLineCandidate;
            continue;
        }

        //
        // Scan through the second line making sure it's all bases (or 'N').  We don't have to
        // check for end-of-buffer because we know there's a null there.
        //
        char *thirdLineCandidate = secondLineCandidate;
        while (*thirdLineCandidate == 'A' || *thirdLineCandidate == 'C' || *thirdLineCandidate == 'T' || *thirdLineCandidate == 'G' ||
                *thirdLineCandidate == 'N' || *thirdLineCandidate == 'a' || *thirdLineCandidate == 'c' || *thirdLineCandidate == 't' || 
                *thirdLineCandidate == 'g' || *thirdLineCandidate == 'n') {
            thirdLineCandidate++;
        }

        if (*thirdLineCandidate == '\r') {
            //
            // CRLF text; skip the CR.
            //
            thirdLineCandidate++;
        }

        if (*thirdLineCandidate != '\n') {
            //
            // We found something that's not a base and not a newline.  It wasn't a read data (second) line.  Move up a line
            // and try again.
            //
            firstLineCandidate = secondLineCandidate;
            continue;
        }

        thirdLineCandidate++;
        if (*thirdLineCandidate != '+') {
            firstLineCandidate = secondLineCandidate;
            continue;
        }

        break;
    }

    data->advance(firstLineCandidate - buffer);
    return true;
}

//
// Try to parse a read starting at a given position pos, updating readToUpdate with it.
// Returns 0 if the parse failed or the first position past the read if it succeeds. In
// addition, if exitOnFailure is set, print a warning and exit if there is a parse error.
// (The one time we don't set this is when trying to find the first read in a chunk.)
// 
    bool
FASTQReader::getNextRead(Read *readToUpdate)
{
    //
    // Find the next newline.
    //

    char* buffer;
    _int64 validBytes;
    if (! data->getData(&buffer, &validBytes)) {
        data->nextBatch();
        if (! data->getData(&buffer, &validBytes)) {
            return false;
        }
    }
    
    _int64 bytesConsumed = getReadFromBuffer(buffer, validBytes, readToUpdate, fileName, data, context);
    if (bytesConsumed == 0) {
        return false;
    }

    data->advance(bytesConsumed);
    return true;
}

// static char LAST[100000]; static int LASTLEN = 0;

    _int64
FASTQReader::getReadFromBuffer(char *buffer, _int64 validBytes, Read *readToUpdate, const char *fileName, DataReader *data, const ReaderContext &context)
{
    //
    // Get the next four lines.
    //
    char* lines[nLinesPerFastqQuery];
    unsigned lineLengths[nLinesPerFastqQuery];
    char* scan = buffer;

    for (unsigned i = 0; i < nLinesPerFastqQuery; i++) {

        char *newLine = strnchr(scan, '\n', validBytes - (scan - buffer));
        if (NULL == newLine) {
            if (validBytes - (scan - buffer) == 1 && *scan == 0x1a && data->isEOF()) {
                // sometimes DOS files will have extra ^Z at end
                return false;
            }
            //
            // There was no next newline
            //
            if (data->isEOF()) {
                WriteErrorMessage("FASTQ file doesn't end with a newline!  Failing.  fileOffset = %lld, validBytes = %d\n",
                    data->getFileOffset(),validBytes);
                soft_exit(1);
            }
            WriteErrorMessage("FASTQ record larger than buffer size at %s:%lld\n", fileName, data->getFileOffset());
            soft_exit(1);
        }

        const size_t lineLen = newLine - scan;
        if (0 == lineLen) {
            WriteErrorMessage("Syntax error in FASTQ file: blank line.\n");
            soft_exit(1);
        }
        if (! isValidStartingCharacterForNextLine[(i + 3) % 4][*scan]) {
            WriteErrorMessage("FASTQ file %s has invalid starting character at offset %lld, line type %d, char %c\n", data->getFilename(), data->getFileOffset(), i, *scan);
            WriteErrorMessage("Line in question: '%.*s'\n", lineLen, scan);
            //WriteErrorMessage("Preceding record: '%.*s'\n", LASTLEN, LAST);
            soft_exit(1);
        }
        lines[i] = scan;
        lineLengths[i] = (unsigned) lineLen - (scan[lineLen-1] == '\r' ? 1 : 0);
        scan = newLine + (newLine[1] == '\r' ? 2 : 1);
    }

    const char *id = lines[0] + 1; // The '@' on the first line is not part of the ID
    const char* space = strnchr(id, ' ', lineLengths[0] - 1);
    readToUpdate->init(id, space != NULL ? (unsigned) (space - id) : (unsigned) lineLengths[0] - 1, lines[1], lines[3], lineLengths[1]);
    readToUpdate->clip(context.clipping);
    readToUpdate->setBatch(data->getBatch());
    readToUpdate->setReadGroup(context.defaultReadGroup);

    // memcpy(LAST, buffer, scan - buffer); LASTLEN = scan - buffer;

    return scan - buffer;

}

//
// static data & initialization
//

bool FASTQReader::isValidStartingCharacterForNextLine[FASTQReader::nLinesPerFastqQuery][256];

FASTQReader::_init FASTQReader::_initializer;


FASTQReader::_init::_init()
{
    //
    // Initialize the isValidStartingCharacterForNextLine array.
    //
    memset(isValidStartingCharacterForNextLine, 0, sizeof(isValidStartingCharacterForNextLine));

    //
    // The first line is the descriptor line and must start with an '@'
    //
    isValidStartingCharacterForNextLine[3]['@'] = true;

    //
    // The second line is the read itself and must start with a base or an
    // 'N' in either case.  A . is just a different way to encode an N.
    //
    for (const char*p = "ACTGNURYKMSWBDHVNX."; *p; p++) {
        isValidStartingCharacterForNextLine[0][*p] = true;
        isValidStartingCharacterForNextLine[0][tolower(*p)] = true;
    }

    //

    //
    // The third line is additional sequence idenfitier info and must
    // start with a '+'.
    //
    isValidStartingCharacterForNextLine[1]['+'] = true;

    //
    // Line 4 is the quality line.  It starts with a printable ascii character.
    // It would be nice to rule out the bases, N, + and @ because it might confsue the parser,
    // but some quality lines do start with those...
    //
    for (char i = '!'; i <= '~'; i++) {
        isValidStartingCharacterForNextLine[2][i] = true;
    }
}

//
// PairedInterleavedFASTQReader
//

PairedInterleavedFASTQReader::PairedInterleavedFASTQReader(
    DataReader* i_data,
    const ReaderContext& i_context) :
    data(i_data), context(i_context)
{
}

    PairedInterleavedFASTQReader* 
PairedInterleavedFASTQReader::create(DataSupplier* supplier, const char *fileName, int bufferCount, _int64 startingOffset, _int64 amountOfFileToProcess,
                                        const ReaderContext& context)
{
    DataReader* data = supplier->getDataReader(bufferCount, 2 * maxReadSizeInBytes, 0.0, 0); // 2* because we read in pairs
    PairedInterleavedFASTQReader* fastq = new PairedInterleavedFASTQReader(data, context);
    if (! fastq->init(fileName)) {
        WriteErrorMessage("Unable to initialize PairedInterleavedFASTQReader for file %s\n", fileName);
        soft_exit(1);
    }
    fastq->reinit(startingOffset, amountOfFileToProcess);
    return fastq;
}


        bool 
PairedInterleavedFASTQReader::init(const char* i_fileName)
{
    fileName = i_fileName;
    return data->init(fileName);
}

        bool 
PairedInterleavedFASTQReader::getNextReadPair(Read *read0, Read *read1)
{
   //
    // Find the next newline.
    //

    char* buffer;
    _int64 validBytes;
    if (! data->getData(&buffer, &validBytes)) {
        data->nextBatch();
        if (! data->getData(&buffer, &validBytes)) {
            return false;
        }
    }
    
    _int64 bytesConsumed = FASTQReader::getReadFromBuffer(buffer, validBytes, read0, fileName, data, context);
    if (bytesConsumed == validBytes) {
        WriteErrorMessage("Input file seems to have an odd number of reads.  Ignoring the last one.");
        return false;
    }
    bytesConsumed += FASTQReader::getReadFromBuffer(buffer + bytesConsumed, validBytes - bytesConsumed, read1, fileName, data, context);

    //
    // Validate the Read IDs.
    //
    if (read0->getIdLength() < 2 || memcmp(read0->getId() + read0->getIdLength() - 2, "/1", 2)) {
        WriteErrorMessage("PairedInterleavedFASTQReader: first read of batch doesn't have ID ending with /1: '%.*s'\n", read0->getIdLength(), read0->getId());
        soft_exit(1);
    }

    if (read1->getIdLength() < 2 || memcmp(read1->getId() + read1->getIdLength() - 2, "/2", 2)) {
        WriteErrorMessage("PairedInterleavedFASTQReader: second read of batch doesn't have ID ending with /2: '%.*s'\n", read1->getIdLength(), read1->getId());
        soft_exit(1);
    }

    data->advance(bytesConsumed);
    return true;
}

    void 
PairedInterleavedFASTQReader::reinit(_int64 startingOffset, _int64 amountOfFileToProcess)
{
    data->reinit(startingOffset, amountOfFileToProcess);
    char* buffer;
    _int64 bytes;
    if (! data->getData(&buffer, &bytes)) {
        return;
    }

    // If we're not at the start of the file, we might have the tail end of a read that someone
    // who got the previous range will process; advance past that. This is fairly tricky because
    // there can be '@' signs in the quality string (and maybe even in read names?).
    if (startingOffset != 0) {
        if (!FASTQReader::skipPartialRecord(data)) {
            return;
        }
    }

    //
    // Grab the first read from the buffer, and see if it's /1 or /2.
    //
    if (!data->getData(&buffer, &bytes)) {
        return;
    }

    Read read;
    _int64 bytesForFirstRead = FASTQReader::getReadFromBuffer(buffer, bytes, &read, fileName, data, context);
    if (read.getIdLength() < 2 || read.getId()[read.getIdLength() - 2] != '/' || (read.getId()[read.getIdLength() - 1] != '1' && read.getId()[read.getIdLength() -1] != '2') ) {
        WriteErrorMessage("PairedInterleavedFASTQReader: read ID doesn't appear to end with /1 or /2, you can't use this as a paired FASTQ file: '%.*s'\n", read.getIdLength(), read.getId());
        soft_exit(1);
    }

    if (read.getId()[read.getIdLength()-1] == '2') {
        //
        // This is the second half of a pair.  Skip it.
        //
        data->advance(bytesForFirstRead);

        //
        // Now make sure that the next read is /1.
        //
        if (!data->getData(&buffer, &bytes)) {
            //
            // This was the last read in the file.
            //
            return;
        }

        FASTQReader::getReadFromBuffer(buffer, bytes, &read, fileName, data, context);
        if (read.getIdLength() < 2 || read.getId()[read.getIdLength()-2] != '/' || read.getId()[read.getIdLength()-1] != '1') {
            WriteErrorMessage("PairedInterleavedFASTQReader: first read of pair doesn't appear to have an ID that ends in /1: '%.*s'\n", read.getIdLength(), read.getId());
            soft_exit(1);
        }
    }
}


//
// FASTQWriter
//

    FASTQWriter *
FASTQWriter::Factory(const char *filename)
{
    FILE *file = fopen(filename,"wb");
    if (NULL == file) {
        return NULL;
    }

    return new FASTQWriter(file);
}

    bool
FASTQWriter::writeRead(Read *read)
{
    size_t len = read->getIdLength() + 2 * (read->getDataLength() + 2) /* @ and + */ + 10 /* crlf + padding + null */;

    if (bufferSize - bufferOffset <= len) {
        flushBuffer();
    }

    size_t bytesUsed = snprintf(buffer + bufferOffset, bufferSize - bufferOffset, "@%.*s\n%.*s\n+\n%.*s\n", read->getIdLength(), read->getId(), read->getDataLength(), read->getData(), read->getDataLength(), read->getQuality());
 
    bufferOffset += bytesUsed;
    return true;
}


PairedFASTQReader::~PairedFASTQReader()
{
    for (int i = 0; i < 2; i++) {
        delete readers[i];
        readers[i] = NULL;
    }
}

	PairedFASTQReader *
PairedFASTQReader::create(
    DataSupplier* supplier,
    const char *fileName0,
    const char *fileName1,
    int bufferCount,
    _int64 startingOffset,
    _int64 amountOfFileToProcess,
    const ReaderContext& context)
{
    PairedFASTQReader *reader = new PairedFASTQReader;
    reader->readers[0] = FASTQReader::create(supplier, fileName0, bufferCount, startingOffset, amountOfFileToProcess, context);
    reader->readers[1] = FASTQReader::create(supplier, fileName1, bufferCount, startingOffset, amountOfFileToProcess, context);

    for (int i = 0; i < 2; i++) {
        if (NULL == reader->readers[i]) {
            delete reader;
            reader = NULL;
            return NULL;
        }
    }

    return reader;
}

    bool
PairedFASTQReader::getNextReadPair(Read *read0, Read *read1)
{
    bool worked = readers[0]->getNextRead(read0);
 
    if (readers[1]->getNextRead(read1) != worked) {
        WriteErrorMessage("PairedFASTQReader: reads of both ends responded differently.  The FASTQ files may not match properly.\n");
        soft_exit(1);
    }

    return worked;
}

    PairedReadSupplierGenerator *
PairedFASTQReader::createPairedReadSupplierGenerator(
    const char *fileName0,
    const char *fileName1,
    int numThreads,
    const ReaderContext& context,
    bool gzip)
{
    const char *fileNames[2] = {fileName0, fileName1};
    //
    // Decide whether to use the range splitter or a queue based on whether the files are the same size.
    //
    if (!strcmp("-", fileNames[0]) || !strcmp("-", fileNames[1]) || QueryFileSize(fileNames[0]) != QueryFileSize(fileNames[1]) || gzip) {
        //WriteStatusMessage("FASTQ using supplier queue\n");
        DataSupplier* dataSupplier[2];
        size_t fileSize[2];

        for (int i = 0; i < 2; i++) {
            if (!strcmp(fileNames[i], "-")) {
                fileSize[i] = 0;
                if (gzip) {
                    dataSupplier[i] = DataSupplier::GzipStdio;
                } else {
                    dataSupplier[i] = DataSupplier::Stdio;
                }
            } else {
                fileSize[i] = QueryFileSize(fileNames[i]);
                if (gzip) {
                    dataSupplier[i] = DataSupplier::GzipDefault;
                } else {
                    dataSupplier[i] = DataSupplier::Default;
                }
            }
        }
        
        
        int bufferCount = ReadSupplierQueue::BufferCount(numThreads);
        ReadReader *reader1 = FASTQReader::create(dataSupplier[0], fileName0, bufferCount, 0, fileSize[0],context);
        ReadReader *reader2 = FASTQReader::create(dataSupplier[1], fileName1, bufferCount, 0, fileSize[1],context);
        if (NULL == reader1 || NULL == reader2) {
            delete reader1;
            delete reader2;
            return NULL;
        }
        ReadSupplierQueue *queue = new ReadSupplierQueue(reader1,reader2); 
        queue->startReaders();
        return queue;
    } else {
        //WriteStatusMessage("FASTQ using range splitter\n");
        return new RangeSplittingPairedReadSupplierGenerator(fileName0, fileName1, FASTQFile, numThreads, false, context);
    }
}

    ReadSupplierGenerator *
FASTQReader::createReadSupplierGenerator(
    const char *fileName,
    int numThreads,
    const ReaderContext& context,
    bool gzip)
{
    bool isStdin = !strcmp(fileName,"-");
    if (! gzip && !isStdin) {
        //
        // Single ended uncompressed FASTQ files can be handled by a range splitter.
        //
        return new RangeSplittingReadSupplierGenerator(fileName, false, numThreads, context);
    } else {
        ReadReader* fastq;
        //
        // Because we can only have one stdin reader, we need to use a queue if we're reading from stdin
        //
        if (isStdin) {
            if (gzip) {
                fastq = FASTQReader::create(DataSupplier::GzipStdio, fileName, ReadSupplierQueue::BufferCount(numThreads), 0, 0, context);
            } else {
                fastq = FASTQReader::create(DataSupplier::Stdio, fileName, ReadSupplierQueue::BufferCount(numThreads), 0, 0, context);
            }
        } else {
            fastq = FASTQReader::create(DataSupplier::GzipDefault, fileName, ReadSupplierQueue::BufferCount(numThreads), 0, QueryFileSize(fileName), context);
        }
        if (fastq == NULL) {
            delete fastq;
            return NULL;
        }
        ReadSupplierQueue *queue = new ReadSupplierQueue(fastq);
        queue->startReaders();
        return queue;
    }
}
    
        PairedReadSupplierGenerator *
PairedInterleavedFASTQReader::createPairedReadSupplierGenerator(
    const char *fileName,
    int numThreads,
    const ReaderContext& context,
    bool gzip)
{
     bool isStdin = !strcmp(fileName,"-");
 
     if (gzip || isStdin) {
        //WriteStatusMessage("PairedInterleavedFASTQ using supplier queue\n");
        DataSupplier *dataSupplier;
        if (isStdin) {
            if (gzip) {
                dataSupplier = DataSupplier::GzipStdio;
            } else {
                dataSupplier = DataSupplier::Stdio;
            }
        } else {
            dataSupplier = DataSupplier::GzipDefault;
        }
        
        PairedReadReader *reader = PairedInterleavedFASTQReader::create(dataSupplier, fileName,
            ReadSupplierQueue::BufferCount(numThreads), 0,(stdin ? 0 : QueryFileSize(fileName)),context);
 
        if (NULL == reader ) {
            delete reader;
            return NULL;
        }
        ReadSupplierQueue *queue = new ReadSupplierQueue(reader); 
        queue->startReaders();
        return queue;
    } else {
        //WriteStatusMessage("PairedInterleavedFASTQ using range splitter\n");
        return new RangeSplittingPairedReadSupplierGenerator(fileName, NULL, InterleavedFASTQFile, numThreads, false, context);
    }
}
