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
    _int64 startingOffset,
    _int64 amountOfFileToProcess,
    const ReaderContext& context)
{
    DataReader* data = supplier->getDataReader(maxReadSizeInBytes);
    FASTQReader* fastq = new FASTQReader(data, context);
    if (! fastq->init(fileName)) {
        fprintf(stderr, "Unable to initialize FASTQReader for file %s\n", fileName);
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
        if (!skipPartialRecord()) {
            return;
        }
    }
}

    bool
FASTQReader::skipPartialRecord()
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
// This happens for very small files.  fprintf(stderr,"Unable to find a read in FASTQ buffer (1)\n");
            return false;
        }

        char *secondLineCandidate = strnchr(firstLineCandidate, '\n', validBytes - (firstLineCandidate - buffer)) + 1;
        if (NULL == (secondLineCandidate-1)) {
			fprintf(stderr,"Unable to find a read in FASTQ buffer (2) at %d\n", data->getFileOffset());
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
                *thirdLineCandidate == 'g') {
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

    //
    // Get the next four lines.
    //
    char* lines[nLinesPerFastqQuery];
    unsigned lineLengths[nLinesPerFastqQuery];
    char* scan = buffer;

    for (unsigned i = 0; i < nLinesPerFastqQuery; i++) {

        char *newLine = strnchr(scan, '\n', validBytes - (scan - buffer));
        if (NULL == newLine) {
            //
            // There was no next newline
            //
            if (data->isEOF()) {
                fprintf(stderr,"FASTQ file doesn't end with a newline!  Failing.  fileOffset = %lld, validBytes = %d\n",
                    data->getFileOffset(),validBytes);
                soft_exit(1);
            }
            fprintf(stderr, "FASTQ record larger than buffer size at %s:%lld\n", fileName, data->getFileOffset());
            soft_exit(1);
        }            //

        const size_t lineLen = newLine - scan;
        if (0 == lineLen) {
            fprintf(stderr,"Syntax error in FASTQ file: blank line.\n");
            soft_exit(1);
        }
        if (! isValidStartingCharacterForNextLine[(i + 3) % 4][*scan]) {
            fprintf(stderr, "FASTQ file has invalid starting character at offset %lld", data->getFileOffset());
            soft_exit(1);
        }
        lines[i] = scan;
        lineLengths[i] = (unsigned) lineLen - (scan[lineLen-1] == '\r' ? 1 : 0);
        scan = newLine + (newLine[1] == '\r' ? 2 : 1);
    }

    const char *id = lines[0] + 1; // The '@' on the first line is not part of the ID
    readToUpdate->init(id, (unsigned) lineLengths[0] - 1, lines[1], lines[3], lineLengths[1]);
    readToUpdate->clip(context.clipping);
    readToUpdate->setBatch(data->getBatch());
    readToUpdate->setReadGroup(context.defaultReadGroup);
    data->advance(scan - buffer);
    return true;

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
    // 'N' in either case.
    //
    isValidStartingCharacterForNextLine[0]['A'] = true;
    isValidStartingCharacterForNextLine[0]['C'] = true;
    isValidStartingCharacterForNextLine[0]['T'] = true;
    isValidStartingCharacterForNextLine[0]['G'] = true;
    isValidStartingCharacterForNextLine[0]['N'] = true;
    isValidStartingCharacterForNextLine[0]['a'] = true;
    isValidStartingCharacterForNextLine[0]['c'] = true;
    isValidStartingCharacterForNextLine[0]['t'] = true;
    isValidStartingCharacterForNextLine[0]['g'] = true;
    isValidStartingCharacterForNextLine[0]['n'] = true;

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
    size_t len = __max(read->getIdLength(),read->getDataLength());
    char *buffer = new char [len + 1];
    if (NULL == buffer) {
        fprintf(stderr,"FASTQWriter: out of memory\n");
        soft_exit(1);
    }

    memcpy(buffer,read->getId(),read->getIdLength());
    buffer[read->getIdLength()] = 0;
    bool worked = fprintf(outputFile,"@%s\n",buffer) > 0;

    memcpy(buffer,read->getData(),read->getDataLength());
    buffer[read->getDataLength()] = 0;
    worked &= fprintf(outputFile,"%s\n+\n",buffer) > 0;

    memcpy(buffer,read->getQuality(), read->getDataLength());   // The FASTQ format spec requires the same number of characters in the quality and data strings
    // Null is already in place from writing data

    worked &= fprintf(outputFile,"%s\n",buffer) > 0;

    delete [] buffer;    

    return worked;
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
    _int64 startingOffset,
    _int64 amountOfFileToProcess,
    const ReaderContext& context)
{
    PairedFASTQReader *reader = new PairedFASTQReader;
    reader->readers[0] = FASTQReader::create(supplier, fileName0,startingOffset,amountOfFileToProcess,context);
    reader->readers[1] = FASTQReader::create(supplier, fileName1,startingOffset,amountOfFileToProcess,context);

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
        fprintf(stderr,"PairedFASTQReader: reads of both ends responded differently.  The FASTQ files may not match properly.\n");
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
    //
    // Decide whether to use the range splitter or a queue based on whether the files are the same size.
    //
    if (QueryFileSize(fileName0) != QueryFileSize(fileName1) || gzip) {
        fprintf(stderr,"FASTQ using supplier queue\n");
        DataSupplier* dataSupplier = gzip ? DataSupplier::GzipDefault[false] : DataSupplier::Default[false];
        ReadReader *reader1 = FASTQReader::create(dataSupplier, fileName0,0,QueryFileSize(fileName0),context);
        ReadReader *reader2 = FASTQReader::create(dataSupplier, fileName1,0,QueryFileSize(fileName1),context);
        if (NULL == reader1 || NULL == reader2) {
            delete reader1;
            delete reader2;
            return NULL;
        }
        ReadSupplierQueue *queue = new ReadSupplierQueue(reader1,reader2); 
        queue->startReaders();
        return queue;
    } else {
        fprintf(stderr,"FASTQ using range splitter\n");
        return new RangeSplittingPairedReadSupplierGenerator(fileName0,fileName1,false,numThreads,context);
    }
}

    ReadSupplierGenerator *
FASTQReader::createReadSupplierGenerator(
    const char *fileName,
    int numThreads,
    const ReaderContext& context,
    bool gzip)
{
    if (! gzip) {
        //
        // Single ended uncompressed FASTQ files can be handled by a range splitter.
        //
        return new RangeSplittingReadSupplierGenerator(fileName, false, numThreads, context);
    } else {
        ReadReader* fastq = FASTQReader::create(DataSupplier::GzipDefault[false], fileName, 0, QueryFileSize(fileName), context);
        if (fastq == NULL) {
            delete fastq;
            return NULL;
        }
        ReadSupplierQueue *queue = new ReadSupplierQueue(fastq);
        queue->startReaders();
        return queue;
    }
}
    
