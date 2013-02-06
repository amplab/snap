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
#include "VariableSizeVector.h"
#include "BufferedAsync.h"
#include "Read.h"
#include "DataReader.h"

/*
 * Output aligned reads in SAM format. See http://samtools.sourceforge.net/SAM1.pdf for details.
 */
class SAMWriter {
public:
    virtual ~SAMWriter();

    virtual bool write(Read *read, AlignmentResult result, unsigned genomeLocation, bool isRC) = 0;

    virtual bool writePair(Read *read1, Read *read2, PairedAlignmentResult *result) = 0;

    virtual bool close() = 0;

    static bool generateHeader(const Genome *genome, char *header, size_t headerBufferSize, size_t *headerActualSize, bool sorted, int argc, const char **argv, const char *version, const char *rgLine);
   
    static const int HEADER_BUFFER_SIZE = 256 * 1024 * 1024;
    
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
						bool						useM,
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
                        size_t *                    spaceUsed = NULL,
                        size_t                      qnameLen = 0
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
                        bool                        isRC,
						bool						useM,
                        int                         *editDistance
                );

};

//
// There is one of these per thread.  Unlike SimpleSAMWriter they all
// share a single file.  Each thread maintains its own write buffer.  When the buffer is full,
// it allocates a write offset by doing an interlocked add of the writeOffset, which is shared
// among all writers of this file.
// Abstract base class with a subclass per platform.
//

class ThreadSAMWriter : public SAMWriter {
public:

    ThreadSAMWriter(size_t i_bufferSize, bool i_useM);
    virtual ~ThreadSAMWriter();

    bool initialize(AsyncFile* file, const Genome *i_genome, volatile _int64 *i_nextWriteOffset);
    
    bool write(Read *read, AlignmentResult result, unsigned genomeLocation, bool isRC);

    bool writePair(Read *read0, Read *read1, PairedAlignmentResult *result);

    virtual bool close();

protected:

    // hooks to allow for sorting the output buffer
    virtual void                    afterWrite(unsigned location, size_t bufferOffset, unsigned length) {}
    virtual bool                    beforeFlush(_int64 fileOffset, _int64 length) { return true; }

    bool                            startIo();

    bool                            waitForIoCompletion();

    volatile _int64                *nextWriteOffset;

    const size_t                    bufferSize;
    size_t                          remainingBufferSpace;
    
    unsigned                        bufferBeingCreated; // Which buffer are we generating new SAM into?
    AsyncFile::Writer              *writer[2];
    char                           *buffer[2];
    const Genome                   *genome;

    LandauVishkinWithCigar          lv;

	bool							useM;
};

class ParallelSAMWriter {
public:
    virtual ~ParallelSAMWriter();

    ParallelSAMWriter() : writer(NULL), file(NULL), argv(NULL), version(NULL), rgLine(NULL) {}

    static ParallelSAMWriter*       create(const char *fileName, const Genome *genome,
                                        unsigned nThreads, bool sort, size_t sortBufferMemory, bool i_useM,
                                        int argc, const char **argv, const char *version, const char *rgLine);

    static const size_t             UnsortedBufferSize = 16 * 1024 * 1024;

    virtual bool                    initialize(const char *fileName, const Genome *genome, unsigned nThreads, bool sorted);

    SAMWriter *                     getWriterForThread(int whichThread);

    virtual bool                    close();

protected:
									ParallelSAMWriter(bool i_useM, int i_argc, const char **i_argv, const char *i_version,
                                        const char *i_rgLine) : 
                                        useM(i_useM), argc(i_argc), argv(i_argv), version(i_version), rgLine(i_rgLine), writer(NULL) {}

    virtual bool                    createThreadWriters(const Genome* genome);

    AsyncFile                      *file;
    volatile _int64                 nextWriteOffset;

    ThreadSAMWriter               **writer;
    int                             nThreads;

	bool							useM;
    int                             argc;
    const char                    **argv;
    const char                     *version;
    const char                     *rgLine;
};

#pragma pack(push, 4)
struct SortEntry
{
    SortEntry() : offset(0), length(0), location(0) {}
    SortEntry(size_t i_offset, unsigned i_length, unsigned i_location)
        : offset(i_offset), length(i_length), location(i_location) {}
    size_t                      offset; // offset in file
    unsigned                    length; // number of bytes
    unsigned                    location; // location in genome
    static bool comparator(const SortEntry& e1, const SortEntry& e2)
    {
        return e1.location < e2.location;
    }
};
#pragma pack(pop)

struct SortBlock
{
    SortBlock();
    SortBlock(size_t capacity);
    SortBlock(SortBlock& other);
    void operator=(SortBlock& other);

    VariableSizeVector<SortEntry>   entries;
    size_t                          fileOffset;
    size_t                          fileBytes;
    unsigned                        index;
    BufferedAsyncReader             reader; // for reading phase
};

class SortedParallelSAMWriter;

class SortedThreadSAMWriter : public ThreadSAMWriter
{
public:

    SortedThreadSAMWriter(size_t i_bufferSize, bool useM);
    virtual ~SortedThreadSAMWriter();


    bool                            initialize(SortedParallelSAMWriter* i_parent, const Genome* i_genome);

protected:
    
    virtual void                    afterWrite(unsigned location, size_t bufferOffset, unsigned length);

    virtual bool                    beforeFlush(_int64 fileOffset, _int64 length);

private:
    SortedParallelSAMWriter*        parent;
    int                             largest; // largest location count so far
    SortBlock                       locations;
};

class SortedParallelSAMWriter : public ParallelSAMWriter
{
public:

    SortedParallelSAMWriter(size_t i_totalMemory, bool useM, int argc, const char**argv, const char *version, const char *rgLine) : 
      ParallelSAMWriter(useM, argc, argv, version, rgLine), totalMemory(i_totalMemory), tempFile(NULL), sortedFile(NULL) {}

    virtual ~SortedParallelSAMWriter() {}
    
    virtual bool                    initialize(const char *fileName, const Genome *genome, unsigned nThreads, bool sorted);

    bool                            close();

protected:
    
    virtual bool                    createThreadWriters(const Genome* genome);

private:
    
    friend class SortedThreadSAMWriter;

    void                            addLocations(SortBlock& added);

    bool                            mergeSort();

    bool                            memoryMappedSort();

    const size_t                    totalMemory;
    size_t                          headerSize;
    char*                           tempFile;
    const char*                     sortedFile;
    ExclusiveLock                   lock;
    VariableSizeVector<SortBlock>   locations;
};

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

class SAMReader : public ReadReader {
public:
        virtual ~SAMReader() {}

        SAMReader(DataReader* i_data, ReadClippingType i_clipping);

        virtual void reinit(_int64 startingOffset, _int64 amountOfFileToProcess);

        virtual bool getNextRead(Read *readToUpdate);
    
        virtual bool getNextRead(Read *read, AlignmentResult *alignmentResult, unsigned *genomeLocation, bool *isRC, unsigned *mapQ,
                        unsigned *flag, const char **cigar)
        {
            return getNextRead(read,alignmentResult,genomeLocation,isRC,mapQ,flag,false,cigar);
        }

        void releaseBefore(DataBatch batch)
        { data->releaseBefore(batch); }

        static SAMReader* create(const DataSupplier* supplier, const char *fileName, const Genome *genome, _int64 startingOffset, _int64 amountOfFileToProcess, 
                                 ReadClippingType clipping = ClipBack);
        
        static PairedReadReader* createPairedReader(const DataSupplier* supplier, const char *fileName, const Genome *genome, _int64 startingOffset, _int64 amountOfFileToProcess, 
                                 ReadClippingType clipping = ClipBack);

        static ReadSupplierGenerator *createReadSupplierGenerator(const char *fileName, int numThreads, const Genome *genome, ReadClippingType clipping = ClipBack);
        static PairedReadSupplierGenerator *createPairedReadSupplierGenerator(const char *fileName, int numThreads, const Genome *genome, ReadClippingType clipping = ClipBack);
        
        // result and fieldLengths must be of size nSAMFields
        static bool parseHeader(const char *fileName, char *firstLine, char *endOfBuffer, const Genome *genome, _int64 *o_headerSize);
        
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

        static const int maxLineLen = 2048;

        static bool parseLine(char *line, char *endOfBuffer, char *result[], size_t *lineLength, size_t fieldLengths[]);

        virtual bool getNextRead(Read *read, AlignmentResult *alignmentResult, 
                        unsigned *genomeLocation, bool *isRC, unsigned *mapQ, unsigned *flag, bool ignoreEndOfRange, const char **cigar);

        static void getReadFromLine(const Genome *genome, char *line, char *endOfBuffer, Read *read, AlignmentResult *alignmentResult,
                        unsigned *genomeLocation, bool *isRC, unsigned *mapQ, 
                        size_t *lineLength, unsigned *flag, const char **cigar, ReadClippingType clipping);

private:
        bool init(const char *fileName, const Genome *i_genome, _int64 startingOffset, _int64 amountOfFileToProcess);
        friend class SAMReader;

        DataReader*         data;
        _int64              headerSize;
        ReadClippingType    clipping;

        bool                didInitialSkip;   // Have we skipped to the beginning of the first SAM line?  We may start in the middle of one.

        const Genome *      genome;
};
