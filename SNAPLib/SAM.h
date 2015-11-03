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
#include "directions.h"
#include "Read.h"
#include "DataReader.h"
#include "FileFormat.h"

bool readIdsMatch(const char* id0, const char* id1);

bool readIdsMatch(Read *read0, Read *read1);

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
const int SAM_FAILED_QC          = 0x200; // Not passing quality controls.
const int SAM_DUPLICATE          = 0x400; // PCR or optical duplicate.
const int SAM_SUPPLEMENTARY      = 0x800; // Supplementary alignment

class SAMReader : public ReadReader {
public:
        virtual ~SAMReader() {}

        SAMReader(DataReader* i_data, const ReaderContext& i_context);

        virtual void reinit(_int64 startingOffset, _int64 amountOfFileToProcess);

        virtual bool getNextRead(Read *readToUpdate);
    
        virtual bool getNextRead(Read *read, AlignmentResult *alignmentResult, GenomeLocation *genomeLocation, Direction *direction, unsigned *mapQ,
                        unsigned *flag, const char **cigar)
        {
            return getNextRead(read, alignmentResult, genomeLocation, direction, mapQ, flag, false, cigar);
        }
        
        virtual void holdBatch(DataBatch batch)
        { data->holdBatch(batch); }

        virtual bool releaseBatch(DataBatch batch)
        { return data->releaseBatch(batch); }
        
        static SAMReader* create(DataSupplier* supplier, const char *fileName,
                int bufferCount, const ReaderContext& i_context,
                _int64 startingOffset, _int64 amountOfFileToProcess);
        
        static PairedReadReader* createPairedReader(const DataSupplier* supplier,
                const char *fileName, int bufferCount, _int64 startingOffset, _int64 amountOfFileToProcess, 
                bool quicklyDropUnpairedReads, const ReaderContext& context);

        static ReadSupplierGenerator *createReadSupplierGenerator(
            const char *fileName, int numThreads, const ReaderContext& context);

        static PairedReadSupplierGenerator *createPairedReadSupplierGenerator(
            const char *fileName, int numThreads, bool quicklyDropUnpairedReads, const ReaderContext& context);
        
        // result and fieldLengths must be of size nSAMFields
        static bool parseHeader(const char *fileName, char *firstLine, char *endOfBuffer, const Genome *genome, _int64 *o_headerSize, bool* o_headerMatchesIndex, bool *o_sawWholeHeader = NULL);
        
        static char* skipToBeyondNextFieldSeparator(char *str, const char *endOfBuffer, size_t *o_charsUntilFirstSeparator = NULL);


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
        static const unsigned  OPT          = 11;
        static const unsigned  nSAMFields   = 12;

        static const int maxLineLen = MAX_READ_LENGTH * 5;

        static bool parseLine(char *line, char *endOfBuffer, char *result[],
            size_t *lineLength, size_t fieldLengths[]);

        static size_t parseContigName(const Genome* genome, char* contigName,
            size_t contigNameBufferSize, GenomeLocation * o_locationOfContig, int* o_indexOfContig,
            char* field[], size_t fieldLength[], unsigned rfield = RNAME);  // Returns 0 on success, needed contigNameBufferSize otherwise.

        static GenomeLocation parseLocation(GenomeLocation locationOfContig, char* field[], size_t fieldLength[], unsigned rfield = RNAME, unsigned posfield = POS);

        virtual bool getNextRead(Read *read, AlignmentResult *alignmentResult, 
                        GenomeLocation *genomeLocation, Direction *direction, unsigned *mapQ, unsigned *flag, bool ignoreEndOfRange, const char **cigar);

        static void getReadFromLine(const Genome *genome, char *line, char *endOfBuffer, Read *read, AlignmentResult *alignmentResult,
                        GenomeLocation *genomeLocation, Direction *direction, unsigned *mapQ, 
                        size_t *lineLength, unsigned *flag, const char **cigar, ReadClippingType clipping);


private:
        void readHeader(const char* fileName);

		bool skipPartialHeader(_int64 *o_headerBytes);

        void init(const char *fileName, _int64 startingOffset, _int64 amountOfFileToProcess);

        DataReader*         data;
        _int64              headerSize;
        ReadClippingType    clipping;

        bool                didInitialSkip;   // Have we skipped to the beginning of the first SAM line?  We may start in the middle of one.

        friend class SAMFormat;
};

class SAMFormat : public FileFormat
{
public:
    SAMFormat(bool i_useM) : useM(i_useM) {}

    virtual void getSortInfo(const Genome* genome, char* buffer, _int64 bytes, GenomeLocation* o_location, GenomeDistance* o_readBytes, int* o_refID, int* o_pos) const;

    virtual void setupReaderContext(AlignerOptions* options, ReaderContext* readerContext) const
    { FileFormat::setupReaderContext(options, readerContext, false); }

    virtual ReadWriterSupplier* getWriterSupplier(AlignerOptions* options, const Genome* genome) const;

    virtual bool writeHeader(
        const ReaderContext& context, char *header, size_t headerBufferSize, size_t *headerActualSize,
        bool sorted, int argc, const char **argv, const char *version, const char *rgLine, bool omitSQLines) const;

    virtual bool writeRead(
        const ReaderContext& context, LandauVishkinWithCigar * lv, char * buffer, size_t bufferSpace,
        size_t * spaceUsed, size_t qnameLen, Read * read, AlignmentResult result, 
        int mapQuality, GenomeLocation genomeLocation, Direction direction, bool secondaryAlignment, int* o_addFrontClipping,
        bool hasMate = false, bool firstInPair = false, Read * mate = NULL, 
        AlignmentResult mateResult = NotFound, GenomeLocation mateLocation = 0, Direction mateDirection = FORWARD,
        bool alignedAsPair = false) const; 

    // calculate data needed to write SAM/BAM record
    // very long argument list since this was extracted from
    // original SAM record writing routine so it could be shared with BAM

    static bool
    createSAMLine(
        const Genome * genome,
        LandauVishkinWithCigar * lv,
        // output data
        char* data,
        char* quality,
        GenomeDistance dataSize,
        const char*& contigName,
        int& contigIndex,
        int& flags,
        GenomeDistance& positionInContig,
        int& mapQuality,
        const char*& mateContigName,
        int& mateContigIndex,
        GenomeDistance& matePositionInContig,
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
        GenomeLocation genomeLocation,
        Direction direction,
        bool secondaryAlignment,
        bool useM,
        bool hasMate,
        bool firstInPair,
        bool alignedAsPair,
        Read * mate, 
        AlignmentResult mateResult,
        GenomeLocation mateLocation,
        Direction mateDirection,
        GenomeDistance *extraBasesClippedBefore);

    static void computeCigar(CigarFormat cigarFormat, const Genome * genome, LandauVishkinWithCigar * lv,
        char * cigarBuf, int cigarBufLen,
        const char * data, GenomeDistance dataLength, unsigned basesClippedBefore, GenomeDistance extraBasesClippedBefore, unsigned basesClippedAfter,
        GenomeDistance *o_extraBasesClippedAfter, 
        GenomeLocation genomeLocation, bool useM, int * o_editDistance, int *o_cigarBufUsed, int * o_addFrontClipping);

private:
    static const char * computeCigarString(const Genome * genome, LandauVishkinWithCigar * lv,
        char * cigarBuf, int cigarBufLen, char * cigarBufWithClipping, int cigarBufWithClippingLen,
        const char * data, GenomeDistance dataLength, unsigned basesClippedBefore, GenomeDistance extraBasesClippedBefore, unsigned basesClippedAfter, 
        unsigned frontHardClipped, unsigned backHardClipped,
        GenomeLocation genomeLocation, Direction direction, bool useM, int * o_editDistance, int * o_addFrontClipping);

#ifdef _DEBUG
	static void validateCigarString(const Genome *genome, const char * cigarBuf, int cigarBufLen, const char *data, GenomeDistance dataLength, GenomeLocation genomeLocation, Direction direction, bool useM);
#else	// DEBUG
	inline static void validateCigarString(const Genome *genome, const char * cigarBuf, int cigarBufLen, const char *data, GenomeDistance dataLength, GenomeLocation genomeLocation, Direction direction, bool useM) {}
#endif // DEBUG


    const bool useM;
};
