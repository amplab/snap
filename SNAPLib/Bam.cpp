/*++

Module Name:

    Bam.cpp

Abstract:

    Binary Alignment Map (BAM) file writer and reader.

Environment:

    User mode service.

    BamWriter and BamReader (and their subclasses) aren't thread safe.

--*/

#include "stdafx.h"
#include "BigAlloc.h"
#include "Compat.h"
#include "Read.h"
#include "Bam.h"
#include "Tables.h"
#include "RangeSplitter.h"
#include "ParallelTask.h"

using std::max;
using std::min;

BAMReader::~BAMReader()
{
}

    bool
BAMReader::getNextReadPair(
    Read *read1,
    Read *read2,
    PairedAlignmentResult *alignmentResult,
    unsigned *mapQ,
    const char **cigar)
{
}


    void
BAMReader::readDoneWithBuffer(unsigned *referenceCount)
{
}

    BAMReader*
BAMReader::create(
    const char *fileName,
    const Genome *genome,
    _int64 startingOffset,
    _int64 amountOfFileToProcess,
    ReadClippingType clipping)
{
}

    ReadSupplierGenerator *
BAMReader::createReadSupplierGenerator(
    const char *fileName,
    int numThreads,
    const Genome *genome,
    ReadClippingType clipping)
{
}

    PairedReadSupplierGenerator *
BAMReader::createPairedReadSupplierGenerator(
    const char *fileName,
    int numThreads,
    const Genome *genome,
    ReadClippingType clipping)
{
}

    bool
BAMReader::parseHeader(
    const char *fileName,
    char *firstLine,
    char *endOfBuffer,
    const Genome *genome,
    size_t *headerSize)
{
}
        
    bool
BAMReader::parseLine(
    char *line,
    char *endOfBuffer,
    char *result[],
    size_t *lineLength,
    size_t fieldLengths[])
{
}

    bool
BAMReader::getNextRead(
    Read *read,
    AlignmentResult *alignmentResult, 
    unsigned *genomeLocation,
    bool *isRC,
    unsigned *mapQ,
    unsigned *flag,
    bool ignoreEndOfRange,
    const char **cigar)
{
}

    void
BAMReader::getReadFromLine(
    const Genome *genome,
    char *line,
    char *endOfBuffer,
    Read *read,
    AlignmentResult *alignmentResult,
    unsigned *genomeLocation,
    bool *isRC,
    unsigned *mapQ, 
    size_t *lineLength,
    unsigned *flag,
    const char **cigar,
    ReadClippingType clipping)
{
}
