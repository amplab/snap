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

class   FASTQReader : public ReadReader {
public:
        virtual ~FASTQReader();


        static FASTQReader* create(const char *fileName, _int64 startingOffset, _int64 amountOfFileToProcess,
                                   ReadClippingType clipping = ClipBack);

        static ReadSupplierGenerator *createReadSupplierGenerator(const char *fileName, int numThreads, ReadClippingType clipping = ClipBack);
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

        static PairedReadSupplierGenerator *createPairedReadSupplierGenerator(const char *fileName0, const char *fileName1, int numThreads, ReadClippingType clipping = ClipBack);
 

private:

        PairedFASTQReader()
        {
            for (int i =0; i < 2; i++) {
                readers[i] = NULL;
            }
        }

        FASTQReader *readers[2];
};

class   MemMapFASTQReader : public FASTQReader {
public:
        MemMapFASTQReader(const char *fileName, _int64 startingOffset, _int64 amountOfFileToProcess, ReadClippingType i_clipping);

        virtual ~MemMapFASTQReader();

        virtual bool getNextRead(Read *readToUpdate);

        virtual void readDoneWithBuffer(unsigned *referenceCount);

        virtual void reinit(_int64 startingOffset, _int64 amountOfFileToProcess);

private:
        FileMapper          fileMapper;

        size_t              fileSize;
        size_t              offsetMapped;
        char                *fileData;
        _uint64             pos;             // Current position within the range of the file we mapped
        _uint64             endPos;          // Where the range we were requested to parse ends; we might look one read past this
        _uint64             amountMapped;    // Where our current mmap() ends
        ReadClippingType    clipping;

        static const int maxReadSizeInBytes = 25000;
 
        _uint64 parseRead(_uint64 pos, Read *readToUpdate, bool exitOnFailure);

        void unmapCurrentRange();
};

class FASTQWriter { 
public:
        ~FASTQWriter() {fclose(outputFile);};

        static FASTQWriter *Factory(const char *filename);

        bool writeRead(Read *readToWrite);

private:

        FASTQWriter(FILE *i_outputFile) : outputFile(i_outputFile) {}

        FILE *outputFile;
};


