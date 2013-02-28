/*++

Module Name:

    ReadWriter.cpp

Abstract:

    General file writer.

Environment:

    User mode service.

    Not thread safe.

--*/

#include "stdafx.h"
#include "BigAlloc.h"
#include "Compat.h"
#include "Read.h"
#include "SAM.h"
#include "Tables.h"
#include "RangeSplitter.h"
#include "ParallelTask.h"
#include "Util.h"
#include "ReadSupplierQueue.h"
#include "FileFormat.h"

class SimpleReadWriter : public ReadWriter
{
public:
    SimpleReadWriter(const FileFormat* i_format, DataWriter* i_writer, const Genome* i_genome)
        : format(i_format), writer(i_writer), genome(i_genome)
    {}

    ~SimpleReadWriter()
    {
        delete writer;
    }

    virtual bool writeHeader(bool sorted, int argc, const char **argv, const char *version, const char *rgLine);

    virtual bool writeRead(Read *read, AlignmentResult result, int mapQuality, unsigned genomeLocation, Direction direction);

    virtual bool writePair(Read *read0, Read *read1, PairedAlignmentResult *result);

    virtual void close();

private:
    const FileFormat* format;
    DataWriter* writer;
    const Genome* genome;
    LandauVishkinWithCigar lvc;
};

    bool
SimpleReadWriter::writeHeader(
    bool sorted,
    int argc,
    const char **argv,
    const char *version,
    const char *rgLine)
{
    char* buffer;
    size_t size;
    size_t used;

    if (! writer->getBuffer(&buffer, &size)) {
        return false;
    }

    if (! format->writeHeader(genome, buffer, size, &used, sorted, argc, argv, version, rgLine)) {
        fprintf(stderr, "Failed to write header into fresh buffer\n");
        return false;
    }

    writer->advance(used, 0);
    return true;
}

    bool
SimpleReadWriter::writeRead(
    Read *read,
    AlignmentResult result,
    int mapQuality,
    unsigned genomeLocation,
    Direction direction)
{
    char* buffer;
    size_t size;
    size_t used;
    for (int pass = 0; pass < 2; pass++) {
        if (! writer->getBuffer(&buffer, &size)) {
            return false;
        }
        if (format->writeRead(genome, &lvc, buffer, size, &used, read->getIdLength(), read, result, mapQuality, genomeLocation, direction)) {
            _ASSERT(used <= size);
            writer->advance(used, genomeLocation);
            return true;
        }
        if (pass == 1) {
            fprintf(stderr, "Failed to write into fresh buffer\n");
            exit(1);
        }
        if (! writer->nextBatch()) {
            return false;
        }
    }
    return false; // will never get here
}

    bool
SimpleReadWriter::writePair(
    Read *read0,
    Read *read1,
    PairedAlignmentResult *result)
{
    //
    // We need to write both halves of the pair into the same buffer, so that a write from
    // some other thread doesn't separate them.  So, try the writes and if either doesn't
    // work start IO and try again.
    //
    int first, second;
    Read *reads[2];
    reads[0] = read0;
    reads[1] = read1;

    if (result->location[0] <= result->location[1]) {
        first = 0;
        second = 1;
    } else {
        first = 1;
        second = 0;
    }
    size_t sizeUsed[2];

    //
    // For paired reads, we need to have the same QNAME for both of them, and it needs to be unique among all other
    // reads in the dataset.  For now, all we do is see if the read names end in /1 and /2, and if so truncate them.
    //
    size_t idLengths[2];
    idLengths[0] = read0->getIdLength();
    idLengths[1] = read1->getIdLength();
    if (idLengths[0] == idLengths[1] && idLengths[0] > 2 && read0->getId()[idLengths[0]-2] == '/' && read1->getId()[idLengths[0]-2] == '/') {
        char lastChar0, lastChar1;
        lastChar0 = read0->getId()[idLengths[0] - 1];
        lastChar1 = read1->getId()[idLengths[1] - 1];
        if ((lastChar0 == '1' || lastChar0 == '2') && (lastChar0 == '1' || lastChar1 == '2') && 
            lastChar0 != lastChar1) {
                idLengths[0] -= 2;
                idLengths[1] -= 2;
        }
    }

    for (int pass = 0; pass < 2; pass++) {
        
        char* buffer;
        size_t size;

        if (! writer->getBuffer(&buffer, &size)) {
            return false;
        }

        bool writesFit = format->writeRead(genome, &lvc, buffer, size, &sizeUsed[first],
                            idLengths[first], reads[first], result->status[first],
                            result->mapq[first], result->location[first], result->direction[first],
                            true, true, reads[second], result->status[second],
                            result->location[second], result->direction[second]);

        if (writesFit) {
            writesFit = format->writeRead(genome, &lvc, buffer + sizeUsed[first], size - sizeUsed[first], &sizeUsed[second],
                idLengths[second], reads[second], result->status[second],
                result->mapq[second], result->location[second], result->direction[second],
                true, false, reads[first], result->status[first],
                result->location[first], result->direction[first]);
            if (writesFit) {
                break;
            }
        }

        if (pass == 1) {
            fprintf(stderr,"ReadWriter: write into fresh buffer failed\n");
            return false;
        }

        if (! writer->nextBatch()) {
            return false;
        }
    }

    //
    // The strange code that determines the sort key (which uses the coordinate of the mate for unmapped reads) is because we list unmapped reads
    // with mapped mates at their mates' location so they sort together.  If both halves are unmapped, then  
    writer->advance(sizeUsed[first], result->status[first] != NotFound ? result->location[first] : ((result->status[second] != NotFound) ? result->location[second] : UINT32_MAX));
    writer->advance(sizeUsed[second], result->status[second] != NotFound ? result->location[second] : ((result->status[first] != NotFound) ? result->location[first] : UINT32_MAX));
    return true;
}

    void
SimpleReadWriter::close()
{
    writer->close();
}

class SimpleReadWriterSupplier : public ReadWriterSupplier
{
public:
    SimpleReadWriterSupplier(const FileFormat* i_format, DataWriterSupplier* i_dataSupplier, const Genome* i_genome)
        :
        format(i_format),
        dataSupplier(i_dataSupplier),
        genome(i_genome)
    {}

    ~SimpleReadWriterSupplier()
    {
        delete dataSupplier;
    }

    virtual ReadWriter* getWriter()
    {
        return new SimpleReadWriter(format, dataSupplier->getWriter(), genome);
    }

    virtual void close()
    {
        dataSupplier->close();
    }

private:
    const FileFormat* format;
    DataWriterSupplier* dataSupplier;
    const Genome* genome;
};

    ReadWriterSupplier*
ReadWriterSupplier::create(
    const FileFormat* format,
    DataWriterSupplier* dataSupplier,
    const Genome* genome)
{
    return new SimpleReadWriterSupplier(format, dataSupplier, genome);
}

