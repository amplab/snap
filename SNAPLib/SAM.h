/*++

Module Name:

    SAM.h

Abstract:

    Sequence Alignment Map (SAM) file writer.

Environment:

    User mode service.

    This class is NOT thread safe.  It's the caller's responsibility to ensure that
    at most one thread uses an instance at any time.

--*/

#pragma once

#include "Compat.h"
#include "LandauVishkin.h"
#include "PairedEndAligner.h"

/*
 * Output aligned reads in SAM format. See http://samtools.sourceforge.net/SAM1.pdf for details.
 */
class SAMWriter {
public:
    virtual ~SAMWriter();

    virtual bool write(Read *read, AlignmentResult result, unsigned genomeLocation, bool isRC) = 0;

    virtual bool writePair(Read *read1, Read *read2, PairedAlignmentResult *result) = 0;

    virtual bool close() = 0;

    static SAMWriter* create(const char *fileName, const Genome *genome);

    static bool generateHeader(const Genome *genome, char *header, size_t headerBufferSize, size_t *headerActualSize);
   
protected:

    //
    // Generate a SAM line into a buffer.  The line IS null terminated, but the space for the terminator
    // is not included in the returned spaceUsed value.
    //

    static bool generateSAMText(
                        Read *                      read,
                        AlignmentResult             result, 
                        unsigned                    genomeLocation, 
                        bool                        isRC, 
                        bool                        hasMate, 
                        bool                        firstInPair, 
                        Read *                      mate, 
                        AlignmentResult             mateResult, 
                        unsigned                    mateLocation,
                        bool                        isMateRC, 
                        const Genome *              genome, 
                        LandauVishkinWithCigar *    lv, 
                        char *                      buffer, 
                        size_t                      bufferSpace, 
                        size_t *                    spaceUsed = NULL
                );

private:
                    
    static const char *computeCigarString(
                        const Genome *              genome,
                        LandauVishkinWithCigar *    lv,
                        char *                      cigarBuf,
                        int                         cigarBufLen,
                        char *                      cigarBufWithClipping,
                        int                         cigarBufWithClippingLen,
                        const char *                data,
                        unsigned                    dataLength,
                        unsigned                    basesClippedBefore,
                        unsigned                    basesClippedAfter,
                        unsigned                    genomeLocation,
                        bool                        isRC
                );

};

/*
 * A SAMWriter that uses the C standard library (fwrite & co).
 */
class SimpleSAMWriter: public SAMWriter {
public:
    SimpleSAMWriter();

    virtual ~SimpleSAMWriter();

    bool open(const char* fileName, const Genome *genome);

    virtual bool write(Read *read, AlignmentResult result, unsigned genomeLocation, bool isRC);

    virtual bool writePair(Read *read0, Read *read1, PairedAlignmentResult *result);

    virtual bool close();

private:
    // Write one read's result, whether it's from a pair or not.
    void write(Read *read, AlignmentResult result, unsigned genomeLocation, bool isRC, bool hasMate, bool firstInPair,
               Read *mate, AlignmentResult mateResult, unsigned mateLocation, bool mateIsRC);

    static const int BUFFER_SIZE = 8 * 1024 * 1024;
    static const int HEADER_BUFFER_SIZE = 256 * 1024 * 1024;

    FILE *file;
    char *buffer; // For setvbuf

    const Genome *genome;

    LandauVishkinWithCigar lv;
};

#ifdef  _MSC_VER
//
// Like with SimpleSAMWriter there is one of these per thread.  Unlike SimpleSAMWriter they all
// share a single file.  Each thread maintains its own write buffer.  When the buffer is full,
// it allocates a write offset by doing an interlocked add of the writeOffset, which is shared
// among all writers of this file.
//

class WindowsSAMWriter : public SAMWriter {
public:

    WindowsSAMWriter();
    virtual ~WindowsSAMWriter();

    bool initialize(HANDLE i_hFile, const Genome *i_genome, volatile _int64 *i_nextWriteOffset);
    
    bool write(Read *read, AlignmentResult result, unsigned genomeLocation, bool isRC);

    bool writePair(Read *read0, Read *read1, PairedAlignmentResult *result);

    bool close();

private:

    bool                            startIo();
    bool                            waitForIoCompletion();

    HANDLE                          hFile;
    volatile _int64                *nextWriteOffset;

    static const size_t             BufferSize = 16 * 1024 * 1024;
    size_t                          remainingBufferSpace;
    
    unsigned                        bufferBeingCreated; // Which buffer are we generating new SAM into?
    bool                            writeOutstanding;   // Is a write pending on 1-bufferBeingCreated?
    OVERLAPPED                      lap[2];
    char                           *buffer[2];
    const Genome                   *genome;

    LandauVishkinWithCigar          lv;
};
#endif  //_MSC_VER

class ParallelSAMWriter {
public:
    virtual ~ParallelSAMWriter() {}

    static ParallelSAMWriter* create(const char *fileName, const Genome *genome, unsigned nThreads);

    virtual SAMWriter *getWriterForThread(unsigned whichThread) = 0;

    virtual bool close() = 0;
};

class SimpleParallelSAMWriter : public ParallelSAMWriter {
public:

    SimpleParallelSAMWriter() : nThreads(0), writer(NULL) {}
    virtual ~SimpleParallelSAMWriter();

    SAMWriter *getWriterForThread(unsigned whichThread) {_ASSERT(whichThread < nThreads); return writer[whichThread];}

    static SimpleParallelSAMWriter* create(const char *fileName, const Genome *genome, unsigned i_nThreads);

    virtual bool close();

private:

    unsigned             nThreads;
    SimpleSAMWriter    **writer;
};

#ifdef  _MSC_VER
class WindowsParallelSAMWriter: public ParallelSAMWriter {
public:

    WindowsParallelSAMWriter() : nThreads(0), writer(NULL), hFile(INVALID_HANDLE_VALUE), nextWriteOffset(0) {}
    virtual ~WindowsParallelSAMWriter();

    static WindowsParallelSAMWriter* create(const char *fileName, const Genome *genome, unsigned nThreads);

    SAMWriter *getWriterForThread(unsigned whichThread) {_ASSERT(whichThread < nThreads); return writer[whichThread];}

    virtual bool close();

public:

    HANDLE               hFile;
    unsigned             nThreads;
    WindowsSAMWriter   **writer;
    volatile _int64      nextWriteOffset;
};
#endif  // _MSC_VER


/*
 * Flags for the SAM file format; see http://samtools.sourceforge.net/SAM1.pdf for details.
 */
const int SAM_MULTI_SEGMENT      = 0x001; // Read had multiple segments (i.e., paired ends).
const int SAM_ALL_ALIGNED        = 0x002; // All segments of a multi-segment read were aligned.
const int SAM_UNMAPPED           = 0x004; // This segment of the read is unmapped.
const int SAM_NEXT_UNMAPPED      = 0x008; // Next segment of the read is unmapped.
const int SAM_REVERSE_COMPLEMENT = 0x010; // This segment of the read is reverse complemented.
const int SAM_NEXT_REVERSED      = 0x020; // Next segment of the read is reverse complemented.
const int SAM_FIRST_SEGMENT      = 0x040; // This is the first segment in the read.
const int SAM_LAST_SEGMENT       = 0x080; // This is the last segment in the read.
const int SAM_SECONDARY          = 0x100; // Secondary alignment for a read with multiple hits.

class SAMReader : public PairedReadReader, public ReadReader {
public:
        virtual ~SAMReader();

        virtual bool getNextRead(Read *readToUpdate);
    
        virtual bool getNextRead(Read *read, AlignmentResult *alignmentResult, unsigned *genomeLocation, bool *isRC, unsigned *mapQ,
                        unsigned *flag, const char **cigar)
        {
            return getNextRead(read,alignmentResult,genomeLocation,isRC,mapQ,flag,false,cigar);
        }

        //
        // In getNextReadPair mapQ points to an array of two unsigneds.
        //
        virtual bool getNextReadPair(Read *read1, Read *read2, PairedAlignmentResult *alignmentResult, unsigned *mapQ, const char **cigar);

        //
        // The PairedReadReader version of getNextReadPair, which throws away the alignment, mapQ and cigar values.
        //
            bool
        getNextReadPair(Read *read1, Read *read2) {
            PairedAlignmentResult pairedAlignmentResult;
            unsigned mapQ[2];

            return getNextReadPair(read1,read2,&pairedAlignmentResult,mapQ,NULL);
        }

        virtual void readDoneWithBuffer(unsigned *referenceCount) = 0;

        static SAMReader* create(const char *fileName, const Genome *genome, _int64 startingOffset, _int64 amountOfFileToProcess, 
                                 ReadClippingType clipping = ClipBack);

        virtual ReadReader *getReaderToInitializeRead(int whichHalfOfPair) {
            _ASSERT(0 == whichHalfOfPair || 1 == whichHalfOfPair);
            return this;
        }
protected:

        //
        // 0-based Field numbers for the fields within a SAM line.
        //
        static const unsigned  QNAME        = 0;
        static const unsigned  FLAG         = 1;
        static const unsigned  RNAME        = 2;
        static const unsigned  POS          = 3;
        static const unsigned  MAPQ         = 4;
        static const unsigned  CIGAR        = 5;
        static const unsigned  RNEXT        = 6;
        static const unsigned  PNEXT        = 7;
        static const unsigned  TLEN         = 8;
        static const unsigned  SEQ          = 9;
        static const unsigned  QUAL         = 10;
        static const unsigned  nSAMFields   = 11;

        // result and fieldLengths must be of size nSAMFields
        static bool parseHeader(const char *fileName, char *firstLine, char *endOfBuffer, const Genome *genome, size_t *headerSize);
        
        static bool parseLine(char *line, char *endOfBuffer, char *result[], size_t *lineLength, size_t fieldLengths[]);

        virtual bool getNextRead(Read *read, AlignmentResult *alignmentResult, 
                        unsigned *genomeLocation, bool *isRC, unsigned *mapQ, unsigned *flag, bool ignoreEndOfRange, const char **cigar) = 0;

        static void getReadFromLine(const Genome *genome, char *line, char *endOfBuffer, Read *read, AlignmentResult *alignmentResult,
                        unsigned *genomeLocation, bool *isRC, unsigned *mapQ, 
                        size_t *lineLength, unsigned *flag, unsigned **newReferenceCounts, const char **cigar, ReadClippingType clipping);
};

#ifdef  _MSC_VER

class WindowsOverlappedSAMReader : public SAMReader {
public:
        WindowsOverlappedSAMReader(ReadClippingType i_clipping);
        virtual ~WindowsOverlappedSAMReader();

        virtual void readDoneWithBuffer(unsigned *referenceCount);

        virtual void reinit(_int64 startingOffset, _int64 amountOfFileToProcess);

protected:

        virtual bool getNextRead(Read *read, AlignmentResult *alignmentResult, 
                        unsigned *genomeLocation, bool *isRC, unsigned *mapQ, unsigned *flags, bool ignoreEndOfRange, const char **cigar);

private:
        bool init(const char *fileName, const Genome *i_genome, _int64 startingOffset, _int64 amountOfFileToProcess);
        friend class SAMReader;

        static const unsigned nBuffers = 3;
        static const unsigned bufferSize = 32 * 1024 * 1024 - 4096;

        static const int maxReadSizeInBytes = 25000;    // Read as in sequencer read, not read-from-the-filesystem.
        static const unsigned maxLineLen = maxReadSizeInBytes;

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

        unsigned            nextBufferForReader;
        unsigned            nextBufferForConsumer;

        LARGE_INTEGER       readOffset;
        _int64              endingOffset;
        LARGE_INTEGER       fileSize;
        size_t              headerSize;
        ReadClippingType    clipping;

        bool                didInitialSkip;   // Have we skipped to the beginning of the first SAM line?  We may start in the middle of one.

        HANDLE              hFile;
  
        const Genome *      genome;

        void startIo();
        void waitForBuffer(unsigned bufferNumber);
};

#else   // _MSC_VER

class MemMapSAMReader : public SAMReader {
public:
        MemMapSAMReader(ReadClippingType i_clipping);
        virtual ~MemMapSAMReader();

        virtual void readDoneWithBuffer(unsigned *referenceCount);

        virtual void reinit(_int64 startingOffset, _int64 amountOfFileToProcess);

protected:
        virtual bool getNextRead(Read *read, AlignmentResult *alignmentResult, 
                        unsigned *genomeLocation, bool *isRC, unsigned *mapQ, unsigned *flags, bool ignoreEndOfRange, const char **cigar);

private:
        bool init(const char *fileName, const Genome *i_genome, _int64 startingOffset, _int64 amountOfFileToProcess);
        friend class SAMReader;

        void unmapCurrentRange();
        
        static const int maxReadSizeInBytes = 25000;
        static const int madviseSize = 4 * 1024 * 1024;

        ReadClippingType clipping;
        const Genome *genome;
        int fd;
        size_t fileSize;
        size_t headerSize;

        char *fileData;
        size_t offsetMapped;
        _uint64 pos;             // Current position within the range of the file we mapped
        _uint64 endPos;          // Where the range we were requested to parse ends; we might look one read past this
        _uint64 amountMapped;    // Where our current mmap() ends
        _uint64 lastPosMadvised; // Last position we called madvise() on to initiate a read
};


#endif  // _MSC_VER
