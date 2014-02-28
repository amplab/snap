/*++

Module Name:

    FASTQ.h

Abstract:

    Headers for fast FASTQ genome "query" reader.

Authors:

    Bill Bolosky, August, 2011

Environment:

    User mode service.

    This class is NOT thread safe.  It's the caller's responsibility to ensure that
    at most one thread uses an instance at any time.

Revision History:

    Adapted from Matei Zaharia's Scala implementation.

--*/

#pragma once

#include "Compat.h"
#include "Read.h"
#include "ReadSupplierQueue.h"
#include "RangeSplitter.h"
#include "DataReader.h"
#include "Error.h"

class   FASTQReader : public ReadReader {
public:

        FASTQReader(DataReader* data, const ReaderContext& i_context);

        virtual ~FASTQReader();

        static FASTQReader* create(DataSupplier* supplier,
            const char *fileName,
            int bufferCount,
            _int64 startingOffset,
            _int64 amountOfFileToProcess,
            const ReaderContext& i_context);

        static void readHeader(const char* fileName, ReaderContext& context);

        bool init(const char* i_fileName);

        static ReadSupplierGenerator *createReadSupplierGenerator(const char *fileName, int numThreads, const ReaderContext& context, bool gzip = false);
        
        virtual bool getNextRead(Read *readToUpdate);

        virtual void reinit(_int64 startingOffset, _int64 amountOfFileToProcess);
        
        virtual void holdBatch(DataBatch batch)
        { data->holdBatch(batch); }

        virtual bool releaseBatch(DataBatch batch)
        { return data->releaseBatch(batch); }
        
        static _int64 getReadFromBuffer(char *buffer, _int64 bufferSize, Read *readToUpdate, const char *fileName, DataReader *data, const ReaderContext &context);    // Returns the number of bytes consumed.

        static bool skipPartialRecord(DataReader *data);

private:

        static const int maxReadSizeInBytes = MAX_READ_LENGTH * 2 + 1000;    // Read as in sequencer read, not read-from-the-filesystem.  +1000 is for ID string, + line, newlines, etc.

        DataReader*         data;
        const char*         fileName;

        static const unsigned maxLineLen = MAX_READ_LENGTH + 500;
        static const unsigned nLinesPerFastqQuery = 4;

        static bool isValidStartingCharacterForNextLine[nLinesPerFastqQuery][256];
        static class _init
        {
        public:
            _init();
        } _initializer;
};

//
// Get read pairs from an interleaved FASTQ.  It's the same as an ordinary FASTQ reader, except that it has a different version of
// skipPartialRecord() that goes until it hits the first read in a pair.  It identifies the pairs by looking for /1 and /2 at the
// end of the read IDs.
//
class PairedInterleavedFASTQReader : public PairedReadReader {
public:
        PairedInterleavedFASTQReader(DataReader* data, const ReaderContext& i_context);

        virtual ~PairedInterleavedFASTQReader() {}

        static PairedInterleavedFASTQReader* create(DataSupplier* supplier, const char *fileName, int bufferCount, _int64 startingOffset, _int64 amountOfFileToProcess,
                                   const ReaderContext& i_context);

        static void readHeader(const char* fileName, ReaderContext& context) {
            FASTQReader::readHeader(fileName, context);
        }

        bool init(const char* i_fileName);

        static PairedReadSupplierGenerator *createPairedReadSupplierGenerator(const char *fileName, int numThreads, const ReaderContext& context, bool gzip);

        virtual bool getNextReadPair(Read *read0, Read *read1);

        virtual void reinit(_int64 startingOffset, _int64 amountOfFileToProcess);
        
        virtual void holdBatch(DataBatch batch)
        { data->holdBatch(batch); }

        virtual bool releaseBatch(DataBatch batch)
        { return data->releaseBatch(batch); }

        virtual ReaderContext* getContext()
        { return &context; }

private:

    static const int maxReadSizeInBytes = MAX_READ_LENGTH * 2 + 1000;    // Read as in sequencer read, not read-from-the-filesystem.  +1000 is for ID string, + line, newlines, etc.

        DataReader*             data;
        const char*             fileName;
        ReaderContext           context;
};

class PairedFASTQReader: public PairedReadReader {
public:
        virtual ~PairedFASTQReader();


        static PairedFASTQReader* create(DataSupplier* supplier, const char *fileName0, const char *fileName1,
                                         int bufferCount, _int64 startingOffset, _int64 amountOfFileToProcess, const ReaderContext& context);

        virtual bool getNextReadPair(Read *read0, Read *read1);

        virtual void reinit(_int64 startingOffset, _int64 amountOfFileToProcess) {
            for (int i = 0; i < 2; i++) {
                readers[i]->reinit(startingOffset,amountOfFileToProcess);
            }
        }

        virtual ReadReader *getReaderToInitializeRead(int whichHalfOfPair)
        {
            _ASSERT(0 == whichHalfOfPair || 1 == whichHalfOfPair);
            return readers[whichHalfOfPair];
        }

        static PairedReadSupplierGenerator *createPairedReadSupplierGenerator(const char *fileName0, const char *fileName1, int numThreads, const ReaderContext& context, bool gzip = false);
 
        virtual void holdBatch(DataBatch batch)
        { _ASSERT(false); /* not supported */ }

        virtual bool releaseBatch(DataBatch batch)
        { _ASSERT(false); /* not supported */ return false; }

        virtual ReaderContext* getContext()
        { return readers[0]->getContext(); }

private:

        PairedFASTQReader()
        {
            for (int i =0; i < 2; i++) {
                readers[i] = NULL;
            }
        }

        FASTQReader *readers[2];
};


class FASTQWriter { 
public:
        ~FASTQWriter() {flushBuffer(); delete [] buffer; fclose(outputFile);}
        static FASTQWriter *Factory(const char *filename);

        bool writeRead(Read *readToWrite);

private:

        void flushBuffer()
        {
            if (0 == bufferOffset) {
                return;
            }
            if (1 != fwrite(buffer, bufferOffset, 1, outputFile)) {
                WriteErrorMessage("FASTQWriter: error writing file\n");
            }

            bufferOffset = 0;
        }

        FASTQWriter(FILE *i_outputFile) : outputFile(i_outputFile) {
            bufferSize = 20 * 1024 * 1024;
            buffer = new char[bufferSize];
            bufferOffset = 0;
        }

        FILE *outputFile;

        char *buffer;
        size_t bufferSize;
        size_t bufferOffset;
};


