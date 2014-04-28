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
#include "exit.h"
#include "Error.h"

class SimpleReadWriter : public ReadWriter
{
public:
    SimpleReadWriter(const FileFormat* i_format, DataWriter* i_writer, const Genome* i_genome, const Genome * i_transcriptome, const GTFReader * i_gtf)
        : format(i_format), writer(i_writer), genome(i_genome), transcriptome(i_transcriptome), gtf(i_gtf)
    {}

    virtual ~SimpleReadWriter()
    {
        delete writer;
    }

    virtual bool writeHeader(const ReaderContext& context, bool sorted, int argc, const char **argv, const char *version, const char *rgLine);

    virtual bool writeRead(Read *read, AlignmentResult result, int mapQuality, unsigned genomeLocation, Direction direction, bool secondaryAlignment, bool isTranscriptome, unsigned tlocation);

    virtual bool writePair(Read *read0, Read *read1, PairedAlignmentResult *result, bool secondaryAlignment);

    virtual void close();

private:
    const FileFormat* format;
    DataWriter* writer;
    const Genome* genome;
    const Genome* transcriptome;
    const GTFReader* gtf;
    LandauVishkinWithCigar lvc;
};

    bool
SimpleReadWriter::writeHeader(
    const ReaderContext& context,
    bool sorted,
    int argc,
    const char **argv,
    const char *version,
    const char *rgLine)
{
    char* buffer;
    size_t size;
    size_t used;

	writer->inHeader(true);
    if (! writer->getBuffer(&buffer, &size)) {
        return false;
    }

    if (! format->writeHeader(context, buffer, size, &used, sorted, argc, argv, version, rgLine)) {
        WriteErrorMessage( "Failed to write header into fresh buffer\n");
        return false;
    }

    writer->advance((unsigned)used, 0);
	writer->nextBatch();
	writer->inHeader(false);
    return true;
}

    bool
SimpleReadWriter::writeRead(
    Read *read,
    AlignmentResult result,
    int mapQuality,
    unsigned genomeLocation,
    Direction direction,
    bool secondaryAlignment,
    bool isTranscriptome, 
    unsigned tlocation)
{
    char* buffer;
    size_t size;
    size_t used;
    if (result == NotFound) {
        genomeLocation = UINT32_MAX;
    }
    for (int pass = 0; pass < 2; pass++) {
        if (! writer->getBuffer(&buffer, &size)) {
            return false;
        }
        if (format->writeRead(genome, transcriptome, gtf, &lvc, buffer, size, &used, read->getIdLength(), read, result, mapQuality, genomeLocation, direction, secondaryAlignment, isTranscriptome, tlocation)) {
            _ASSERT(used <= size);

        if (used > 0xffffffff) {
            WriteErrorMessage("SimpleReadWriter:writeRead: used too big\n");
            soft_exit(1);
        }

        writer->advance((unsigned)used, genomeLocation);
            return true;
        }
        if (pass == 1) {
            WriteErrorMessage( "Failed to write into fresh buffer\n");
            soft_exit(1);
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
    PairedAlignmentResult *result,
    bool secondaryAlignment)
{
    //
    // We need to write both halves of the pair into the same buffer, so that a write from
    // some other thread doesn't separate them.  So, try the writes and if either doesn't
    // work start IO and try again.
    //
    Read *reads[2] = {read0, read1};
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
        if ((lastChar0 == '1' || lastChar0 == '2') && (lastChar1 == '1' || lastChar1 == '2') && 
            lastChar0 != lastChar1) {
                idLengths[0] -= 2;
                idLengths[1] -= 2;
        }
    }
    unsigned locations[2];
    locations[0] = result->status[0] != NotFound ? result->location[0] : UINT32_MAX;
    locations[1] = result->status[1] != NotFound ? result->location[1] : UINT32_MAX;
    int first = locations[0] > locations[1];
    int second = 1 - first;
    for (int pass = 0; pass < 2; pass++) {
        
        char* buffer;
        size_t size;

        if (! writer->getBuffer(&buffer, &size)) {
            return false;
        }

        bool writesFit = format->writeRead(genome, transcriptome, gtf, &lvc, buffer, size, &sizeUsed[first],
                            idLengths[first], reads[first], result->status[first], result->mapq[first], locations[first], result->direction[first], secondaryAlignment,
                            result->isTranscriptome[first], result->tlocation[first], true, first == 0,
                            reads[second], result->status[second], locations[second], result->direction[second], result->isTranscriptome[second], result->tlocation[second]);

        if (writesFit) {
            writesFit = format->writeRead(genome, transcriptome, gtf, &lvc, buffer + sizeUsed[first], size - sizeUsed[first], &sizeUsed[second],
                idLengths[second], reads[second], result->status[second], result->mapq[second], locations[second], result->direction[second], secondaryAlignment,
                result->isTranscriptome[second], result->tlocation[second], true, first != 0,
                reads[first], result->status[first], locations[first], result->direction[first], result->isTranscriptome[first], result->tlocation[first]);
            if (writesFit) {
                break;
            }
        }

        if (pass == 1) {
            WriteErrorMessage("ReadWriter: write into fresh buffer failed\n");
            return false;
        }

        if (! writer->nextBatch()) {
            return false;
        }
    }

    if (sizeUsed[0] > 0xffffffff || sizeUsed[1] > 0xffffffff) {
        WriteErrorMessage("SimpleReadWriter::writePair: one or the other (or both) sizeUsed too big\n");
        soft_exit(1);
    }

    //
    // The strange code that determines the sort key (which uses the coordinate of the mate for unmapped reads) is because we list unmapped reads
    // with mapped mates at their mates' location so they sort together.  If both halves are unmapped, then  
    writer->advance((unsigned)sizeUsed[first],
        locations[first] != UINT32_MAX ? locations[first] : locations[second]);

    writer->advance((unsigned)sizeUsed[second],
        locations[second] != UINT32_MAX ? locations[second] : locations[first]);
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
    SimpleReadWriterSupplier(const FileFormat* i_format, DataWriterSupplier* i_dataSupplier, const Genome* i_genome, const Genome* i_transcriptome, const GTFReader* i_gtf)
        :
        format(i_format),
        dataSupplier(i_dataSupplier),
        genome(i_genome),
        transcriptome(i_transcriptome),
        gtf(i_gtf)
    {}

    ~SimpleReadWriterSupplier()
    {
        delete dataSupplier;
    }

    virtual ReadWriter* getWriter()
    {
        return new SimpleReadWriter(format, dataSupplier->getWriter(), genome, transcriptome, gtf);
    }

    virtual void close()
    {
        dataSupplier->close();
    }

private:
    const FileFormat* format;
    DataWriterSupplier* dataSupplier;
    const Genome* genome;
    const Genome* transcriptome;
    const GTFReader* gtf;
};

    ReadWriterSupplier*
ReadWriterSupplier::create(
    const FileFormat* format,
    DataWriterSupplier* dataSupplier,
    const Genome* genome,
    const Genome* transcriptome, 
    const GTFReader* gtf)
{
    return new SimpleReadWriterSupplier(format, dataSupplier, genome, transcriptome, gtf);
}

