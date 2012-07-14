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

class   FASTQReader : public ReadReader {
public:
        virtual ~FASTQReader();


        static FASTQReader* create(const char *fileName, _int64 startingOffset, _int64 amountOfFileToProcess,
                                   ReadClippingType clipping = ClipBack);
};

class PairedFASTQReader: public PairedReadReader {
public:
        virtual ~PairedFASTQReader();


        static PairedFASTQReader* create(const char *fileName0, const char *fileName1, _int64 startingOffset, 
                                         _int64 amountOfFileToProcess, ReadClippingType clipping = ClipBack);

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


private:

        PairedFASTQReader()
        {
            for (int i =0; i < 2; i++) {
                readers[i] = NULL;
            }
        }

        FASTQReader *readers[2];
};

#ifndef _MSC_VER
class   MemMapFASTQReader : public FASTQReader {
public:
        MemMapFASTQReader(const char *fileName, _int64 startingOffset, _int64 amountOfFileToProcess, ReadClippingType i_clipping);

        virtual ~MemMapFASTQReader();

        virtual bool getNextRead(Read *readToUpdate);

        virtual void readDoneWithBuffer(unsigned *referenceCount);

        virtual void reinit(_int64 startingOffset, _int64 amountOfFileToProcess);

private:
        int fd;
        size_t fileSize;
        size_t offsetMapped;
        char *fileData;
        _uint64 pos;             // Current position within the range of the file we mapped
        _uint64 endPos;          // Where the range we were requested to parse ends; we might look one read past this
        _uint64 amountMapped;    // Where our current mmap() ends
        _uint64 lastPosMadvised; // Last position we called madvise() on to initiate a read
        ReadClippingType clipping;

        static const int maxReadSizeInBytes = 25000;
        static const int madviseSize = 4 * 1024 * 1024;

        _uint64 parseRead(_uint64 pos, Read *readToUpdate, bool exitOnFailure);

        void unmapCurrentRange();
};
#endif /* not _MSC_VER */


#ifdef _MSC_VER
class   WindowsFASTQReader : public FASTQReader {
public:
        WindowsFASTQReader(const char *fileName, _int64 startingOffset, _int64 amountOfFileToProcess, ReadClippingType i_clipping);

        virtual ~WindowsFASTQReader();

        virtual bool getNextRead(Read *readToUpdate);

        virtual void readDoneWithBuffer(unsigned *referenceCount);

        virtual void reinit(_int64 startingOffset, _int64 amountOfFileToProcess);
private:

        HANDLE hFile;

        //
        // Use several buffers so that we can run IO in parallel with parsing.
        //
        static const unsigned nBuffers = 3;
        static const unsigned bufferSize = 16 * 1024 * 1024 - 4096;

        static const int maxReadSizeInBytes = 25000;    // Read as in sequencer read, not read-from-the-filesystem.

        static const unsigned maxLineLen = 10000;
        static const unsigned nLinesPerFastqQuery = 4;

        _int64 minLineLengthSeen[nLinesPerFastqQuery];
        bool isValidStartingCharacterForNextLine[nLinesPerFastqQuery][256];

        enum BufferState {Empty, Reading, Full, UsedButReferenced};

        struct BufferInfo {
            char            *buffer;
            BufferState     state;
            DWORD           validBytes;
            DWORD           nBytesThatMayBeginARead;
            unsigned        referenceCount; // How many reads refer to this buffer?
            //
            // Some memory to hold a line that's broken over the end of this buffer and the
            // beginning of the next.
            //
            char            overflowBuffer[maxLineLen+1];
            bool            isEOF;
            unsigned        offset;     // How far has the consumer gotten?

            _int64          fileOffset;

            OVERLAPPED      lap;
        };

        BufferInfo bufferInfo[nBuffers];

        unsigned nextBufferForReader;
        unsigned nextBufferForConsumer;

        LARGE_INTEGER readOffset;
        _int64        endingOffset;
        LARGE_INTEGER fileSize;
        ReadClippingType clipping;

        bool          didInitialSkip;   // Have we skipped to the beginning of the first fastq line?  We may start in the middle of one.

        void startIo();
        void waitForBuffer(unsigned bufferNumber);

};
#endif

class FASTQWriter { 
public:
        ~FASTQWriter() {fclose(outputFile);};

        static FASTQWriter *Factory(const char *filename);

        bool writeRead(Read *readToWrite);

private:

        FASTQWriter(FILE *i_outputFile) : outputFile(i_outputFile) {}

        FILE *outputFile;
};


