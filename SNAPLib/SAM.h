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

        static SAMReader* create(const DataSupplier* supplier, const char *fileName,
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

        static void parsePieceName(const Genome* genome, char* pieceName,
            size_t pieceNameBufferSize, unsigned* o_offsetOfPiece, int* o_indexOfPiece,
            char* field[], size_t fieldLength[], unsigned rfield = RNAME);

        static unsigned parseLocation(unsigned offsetOfPiece, char* field[], size_t fieldLength[], unsigned rfield = RNAME, unsigned posfield = POS);

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

