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

    // return true if the filename is for this format based on extension etc.
    virtual bool isFormatOf(const char* filename) const = 0;
    
    // reading
    //

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

    //
    // writing
    //

    virtual bool writeHeader(
        const Genome *genome, char *header, size_t headerBufferSize, size_t *headerActualSize,
        bool sorted, int argc, const char **argv, const char *version, const char *rgLine) const = 0;

    virtual bool writeRead(
        const Genome * genome, LandauVishkinWithCigar * lv, char * buffer, size_t bufferSpace, 
        size_t * spaceUsed, size_t qnameLen, Read * read, AlignmentResult result, 
        unsigned genomeLocation, bool isRC,
        bool hasMate = false, bool firstInPair = false, Read * mate = NULL, 
        AlignmentResult mateResult = NotFound, unsigned mateLocation = 0, bool mateIsRC = false) const = 0; 

    //
    // formats
    //

    static const FileFormat* SAM[2]; // 0 for =, 1 for M (useM flag)
    static const FileFormat* BAM[2];
    static const FileFormat* FASTQ;
    static const FileFormat* FASTQZ;
};

// calculate data needed to write SAM/BAM record
// very long argument list since this was extracted from
// original SAM record writing routine so it could be shared with BAM

    bool
getSAMData(
    const Genome * genome,
    LandauVishkinWithCigar * lv,
    // output data
    char* data,
    char* quality,
    int dataSize,
    const char*& pieceName,
    int& pieceIndex,
    int& flags,
    unsigned& positionInPiece,
    int& mapQuality,
    const char*& matePieceName,
    unsigned& matePositionInPiece,
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
    bool isRC,
    bool useM,
    bool hasMate,
    bool firstInPair,
    Read * mate, 
    AlignmentResult mateResult,
    unsigned mateLocation,
    bool mateIsRC);
