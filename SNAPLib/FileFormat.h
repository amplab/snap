/*++


Module Name:

    FileFormat.h

Abstract:

    Headers for the FileFormat class for the SNAP sequencer

Authors:

    Ravi Pandya, February 2013

Environment:

    User mode service.

Revision History:

--*/

#pragma once

#include "Compat.h"
#include "Tables.h"
#include "Read.h"
#include "Genome.h"
#include "LandauVishkin.h"
#include "AlignerOptions.h"
#include "Genome.h"

//
// abstract class defining format-specific operations
// for reading and writing files of reads
//
class FileFormat
{
public:
    //
    // files
    //

    // reading
    //

    virtual void getSortInfo(const Genome* genome, char* buffer, _int64 bytes, GenomeLocation* o_location, GenomeDistance* o_readBytes, int* o_refID = NULL, int* o_pos = NULL) const = 0;

    /*

    virtual ReadReader* createReader(const DataSupplier* supplier, const char *fileName,
        const Genome *genome, _int64 startingOffset, _int64 amountOfFileToProcess, 
        ReadClippingType clipping = ClipBack) = 0;
        
    virtual PairedReadReader* createPairedReader(const DataSupplier* supplier, const char *fileName,
        const Genome *genome, _int64 startingOffset, _int64 amountOfFileToProcess, 
        ReadClippingType clipping = ClipBack) = 0;

    virtual ReadSupplierGenerator *createReadSupplierGenerator(const char *fileName, int numThreads,
        const Genome *genome, ReadClippingType clipping = ClipBack) = 0;

    virtual PairedReadSupplierGenerator *createPairedReadSupplierGenerator(const char *fileName,
        int numThreads, const Genome *genome, ReadClippingType clipping = ClipBack) = 0;
        
    // parse header and match against genome
    virtual bool parseHeader(const char *fileName, char *firstLine, char *endOfBuffer,
        const Genome *genome, _int64 *o_headerSize) = 0;
    
    */

protected:

    static const char* RGLineToAux;

    static void setupReaderContext(AlignerOptions* options, ReaderContext* readerContext, bool bam);

public:

    virtual void setupReaderContext(AlignerOptions* options, ReaderContext* readerContext) const = 0;

    //
    // writing
    //

    virtual ReadWriterSupplier* getWriterSupplier(AlignerOptions* options, const Genome* genome) const = 0;

    virtual bool writeHeader(
        const ReaderContext& context, char *header, size_t headerBufferSize, size_t *headerActualSize,
        bool sorted, int argc, const char **argv, const char *version, const char *rgLine, bool omitSQLines) const = 0;

    virtual bool writeRead(
        const ReaderContext& context, LandauVishkinWithCigar * lv, char * buffer, size_t bufferSpace,
        size_t * spaceUsed, size_t qnameLen, Read * read, AlignmentResult result,
        int mapQuality, GenomeLocation genomeLocation, Direction direction, bool secondaryAlignment, int* o_addFrontClipping,
        bool hasMate = false, bool firstInPair = false, Read * mate = NULL, 
        AlignmentResult mateResult = NotFound, GenomeLocation mateLocation = 0, Direction mateDirection = FORWARD,
        bool alignedAsPair = false) const = 0; 

    //
    // formats
    //

    static const FileFormat* SAM[2]; // 0 for =, 1 for M (useM flag)
    static const FileFormat* BAM[2];
    static const FileFormat* FASTQ;
    static const FileFormat* FASTQZ;
};
