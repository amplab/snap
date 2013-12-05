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
#include "GTFReader.h"

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

class SAMReader : public ReadReader {
public:
        virtual ~SAMReader() {}

        SAMReader(DataReader* i_data, const ReaderContext& i_context);

        virtual void reinit(_int64 startingOffset, _int64 amountOfFileToProcess);

        virtual bool getNextRead(Read *readToUpdate);
    
        virtual bool getNextRead(Read *read, AlignmentResult *alignmentResult, unsigned *genomeLocation, Direction *direction, unsigned *mapQ,
                        unsigned *flag, const char **cigar)
        {
            return getNextRead(read, alignmentResult, genomeLocation, direction, mapQ, flag, false, cigar);
        }

        void releaseBatch(DataBatch batch)
        { data->releaseBatch(batch); }

        static void readHeader(const char* fileName, ReaderContext& i_context);

        static SAMReader* create(DataSupplier* supplier, const char *fileName,
                const ReaderContext& i_context,
                _int64 startingOffset, _int64 amountOfFileToProcess);
        
        static PairedReadReader* createPairedReader(const DataSupplier* supplier,
                const char *fileName, _int64 startingOffset, _int64 amountOfFileToProcess, 
                bool autoRelease, const ReaderContext& context);

        static ReadSupplierGenerator *createReadSupplierGenerator(
            const char *fileName, int numThreads, const ReaderContext& context);

        static PairedReadSupplierGenerator *createPairedReadSupplierGenerator(
            const char *fileName, int numThreads, const ReaderContext& context);
        
        // result and fieldLengths must be of size nSAMFields
        static bool parseHeader(const char *fileName, char *firstLine, char *endOfBuffer, const Genome *genome, _int64 *o_headerSize, bool* o_headerMatchesIndex);
        
        static char* skipToBeyondNextRunOfSpacesAndTabs(char *str, const char *endOfBuffer, size_t *charsUntilFirstSpaceOrTab = NULL);


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

        static void parseContigName(const Genome* genome, char* contigName,
            size_t contigNameBufferSize, unsigned* o_offsetOfContig, int* o_indexOfContig,
            char* field[], size_t fieldLength[], unsigned rfield = RNAME);

        static unsigned parseLocation(unsigned offsetOfContig, char* field[], size_t fieldLength[], unsigned rfield = RNAME, unsigned posfield = POS);

        virtual bool getNextRead(Read *read, AlignmentResult *alignmentResult, 
                        unsigned *genomeLocation, Direction *direction, unsigned *mapQ, unsigned *flag, bool ignoreEndOfRange, const char **cigar);

        static void getReadFromLine(const Genome *genome, char *line, char *endOfBuffer, Read *read, AlignmentResult *alignmentResult,
                        unsigned *genomeLocation, Direction *direction, unsigned *mapQ, 
                        size_t *lineLength, unsigned *flag, const char **cigar, ReadClippingType clipping);


private:
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

    virtual bool isFormatOf(const char* filename) const;
    
    virtual void getSortInfo(const Genome* genome, char* buffer, _int64 bytes, unsigned* o_location, unsigned* o_readBytes, int* o_refID, int* o_pos) const;

    virtual ReadWriterSupplier* getWriterSupplier(AlignerOptions* options, const Genome* genome, const Genome* transcriptome, const GTFReader* gtf) const;

    virtual bool writeHeader(
        const ReaderContext& context, char *header, size_t headerBufferSize, size_t *headerActualSize,
        bool sorted, int argc, const char **argv, const char *version, const char *rgLine) const;

    virtual bool writeRead(
        const Genome * genome, const Genome * transcriptome, const GTFReader * gtf, LandauVishkinWithCigar * lv, char * buffer, size_t bufferSpace, 
        size_t * spaceUsed, size_t qnameLen, Read * read, AlignmentResult result, 
        int mapQuality, unsigned genomeLocation, Direction direction, bool isTranscriptome, unsigned tlocation,
        bool hasMate = false, bool firstInPair = false, Read * mate = NULL, 
        AlignmentResult mateResult = NotFound, unsigned mateLocation = 0, Direction mateDirection = FORWARD, bool mateIsTranscriptome = false, unsigned mateTlocation = 0) const; 

    // calculate data needed to write SAM/BAM record
    // very long argument list since this was extracted from
    // original SAM record writing routine so it could be shared with BAM

    static bool
    createSAMLine(
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
        const char*& mateContigName,
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
        unsigned *extraBasesClippedAfter);

private:
    static const char * computeCigarString(const Genome * genome, LandauVishkinWithCigar * lv,
        char * cigarBuf, int cigarBufLen, char * cigarBufWithClipping, int cigarBufWithClippingLen,
        const char * data, unsigned dataLength, unsigned basesClippedBefore, unsigned extraBasesClippedBefore, unsigned basesClippedAfter,
        unsigned extraBasesClippedAfter, unsigned frontHardCliped, unsigned backHardClipped,        
        unsigned genomeLocation, Direction direction, bool useM, int * editDistance, std::vector<unsigned> &tokens);

    const bool useM;
};
