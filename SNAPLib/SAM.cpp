/*++

Module Name:

    SAM.cpp

Abstract:

    Sequence Alignment Map (SAM) file writer and reader.

Environment:

    User mode service.

    SamWriter and SamReader (and their subclasses) aren't thread safe.

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

using std::max;
using std::min;
using util::strnchr;


SAMWriter:: ~SAMWriter()
{
}


    bool
SAMWriter::generateHeader(const Genome *genome, char *header, size_t headerBufferSize, size_t *headerActualSize, bool sorted, int argc, const char **argv, const char *version, const char *rgLine )
{
    char *commandLine;
	size_t commandLineSize = 0;
	for (int i = 0; i < argc; i++) {
		commandLineSize += strlen(argv[i]) + 1;	// +1 is either a space or the terminating null
	}
	commandLine = new char[commandLineSize];
	commandLine[0] = '\0';
	for (int i = 0; i < argc; i++) {
		strcat(commandLine,argv[i]);
		if (i != argc-1) {
			strcat(commandLine," ");
		}
	}

    size_t bytesConsumed = snprintf(header, headerBufferSize, "@HD\tVN:1.4\tSO:%s\n%s\n@PG\tID:SNAP\tPN:SNAP\tCL:%s\tVN:%s\n", 
		sorted ? "coordinate" : "unsorted",
        rgLine == NULL ? "@RG\tID:FASTQ\tSM:sample" : rgLine,
        commandLine,version);

	delete [] commandLine;
	commandLine = NULL;
    if (bytesConsumed >= headerBufferSize) {
        fprintf(stderr,"SAMWriter: header buffer too small\n");
        return false;
    }

#ifndef SKIP_SQ_LINES
    // Write an @SQ line for each chromosome / piece in the genome
    const Genome::Piece *pieces = genome->getPieces();
    int numPieces = genome->getNumPieces();
    unsigned genomeLen = genome->getCountOfBases();
    for (int i = 0; i < numPieces; i++) {
        unsigned start = pieces[i].beginningOffset;
        unsigned end = (i + 1 < numPieces) ? pieces[i+1].beginningOffset : genomeLen;
        bytesConsumed += snprintf(header + bytesConsumed, headerBufferSize - bytesConsumed, "@SQ\tSN:%s\tLN:%u\n", pieces[i].name, end - start);

        if (bytesConsumed >= headerBufferSize) {
            fprintf(stderr,"SAMWriter: header buffer too small\n");
            return false;
        }
    }
#endif // SKIP_SQ_LINES

    *headerActualSize = bytesConsumed;
    return true;
}

// Compute the CIGAR edit sequence string for a read against a given genome location.
// Returns this string if possible or "*" if we fail to compute it (which would likely
// be a bug due to lack of buffer space). The pointer returned may be to cigarBuf so it
// will only be valid until computeCigarString is called again.
    const char *
SAMWriter::computeCigarString(
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
    int *                       editDistance
)
{
    const char *reference = genome->getSubstring(genomeLocation, dataLength);
    if (NULL != reference) {
        *editDistance = lv->computeEditDistance(
                            genome->getSubstring(genomeLocation, dataLength),
                            dataLength,
                            data,
                            dataLength,
                            MAX_K - 1,
                            cigarBuf,
                            cigarBufLen,
						    useM);
    } else {
        //
        // Fell off the end of the chromosome.
        //
        return "*";
    }

    if (*editDistance == -2) {
        fprintf(stderr, "WARNING: computeEditDistance returned -2; cigarBuf may be too small\n");
        return "*";
    } else if (*editDistance == -1) {
        static bool warningPrinted = false;
        if (!warningPrinted) {
            fprintf(stderr, "WARNING: computeEditDistance returned -1; this shouldn't happen\n");
            warningPrinted = true;
        }
        return "*";
    } else {
        // Add some CIGAR instructions for soft-clipping if we've ignored some bases in the read.
        char clipBefore[16] = {'\0'};
        char clipAfter[16] = {'\0'};
        if (basesClippedBefore > 0) {
            snprintf(clipBefore, sizeof(clipBefore), "%uS", basesClippedBefore);
        }
        if (basesClippedAfter > 0) {
            snprintf(clipAfter, sizeof(clipAfter), "%uS", basesClippedAfter);
        }
        snprintf(cigarBufWithClipping, cigarBufWithClippingLen, "%s%s%s", clipBefore, cigarBuf, clipAfter);
        return cigarBufWithClipping;
    }
}


    bool
SAMWriter::generateSAMText(
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
                bool                        mateIsRC, 
                const Genome *              genome, 
                LandauVishkinWithCigar *    lv, 
                char *                      buffer, 
                size_t                      bufferSpace, 
                size_t *                    spaceUsed,
                size_t                      qnameLen)
 {
    const int MAX_READ = 10000;
    const int cigarBufSize = MAX_READ * 2;
    char cigarBuf[cigarBufSize];

    const int cigarBufWithClippingSize = MAX_READ * 2 + 32;
    char cigarBufWithClipping[cigarBufWithClippingSize];

    int flags = 0;
    const char *pieceName = "*";
    unsigned positionInPiece = 0;
    int mapQuality = 0;
    const char *cigar = "*";
    const char *matePieceName = "*";
    unsigned matePositionInPiece = 0;
    _int64 templateLength = 0;
    
    if (0 == qnameLen) {
         qnameLen = read->getIdLength();
    }

    //
    // If the aligner said it didn't find anything, treat it as such.  Sometimes it will emit the
    // best match that it found, even if it's not within the maximum edit distance limit (but will
    // then say NotFound).  Here, we force that to be SAM_UNMAPPED.
    //
    if (NotFound == result) {
        genomeLocation = 0xffffffff;
    }

    // Write the data and quality strings. If the read is reverse complemented, these need to
    // be backwards from the original read. Also, both need to be unclipped.
    char data[MAX_READ];
    char quality[MAX_READ];
    const char *clippedData;
    unsigned clippedLength = read->getDataLength();
    unsigned fullLength = read->getUnclippedLength();
    unsigned basesClippedBefore;
    unsigned basesClippedAfter;
    if (isRC) {
      for (unsigned i = 0; i < fullLength; i++) {
        data[fullLength - 1 - i] = COMPLEMENT[read->getUnclippedData()[i]];
        quality[fullLength - 1 - i] = read->getUnclippedQuality()[i];
      }
      clippedData = &data[fullLength - clippedLength - read->getFrontClippedLength()];
      basesClippedBefore = fullLength - clippedLength - read->getFrontClippedLength();
      basesClippedAfter = read->getFrontClippedLength();
    } else {
      memcpy(data, read->getUnclippedData(), read->getUnclippedLength());
      memcpy(quality, read->getUnclippedQuality(), read->getUnclippedLength());
      clippedData = read->getData();
      basesClippedBefore = read->getFrontClippedLength();
      basesClippedAfter = fullLength - clippedLength - basesClippedBefore;
    }

    int editDistance = -1;
    if (genomeLocation != 0xFFFFFFFF) {
        // This could be either a single hit read or a multiple hit read where we just
        // returned one location, but either way, let's print that location. We will then
        // set the quality to 60 if it was single or 0 if it was multiple. These are the
        // values the SAM FAQ suggests for aligners that don't compute confidence scores.
        if (isRC) {
            flags |= SAM_REVERSE_COMPLEMENT;
        }
        const Genome::Piece *piece = genome->getPieceAtLocation(genomeLocation);
        pieceName = piece->name;
        positionInPiece = genomeLocation - piece->beginningOffset + 1; // SAM is 1-based
        cigar = computeCigarString(genome, lv, cigarBuf, cigarBufSize, cigarBufWithClipping, cigarBufWithClippingSize, 
                                   clippedData, clippedLength, basesClippedBefore, basesClippedAfter,
                                   genomeLocation, isRC, useM, &editDistance);
        mapQuality = (result == SingleHit || result == CertainHit) ? 60 : 0;
    } else {
        flags |= SAM_UNMAPPED;
    }

    if (hasMate) {
        flags |= SAM_MULTI_SEGMENT;
        flags |= (firstInPair ? SAM_FIRST_SEGMENT : SAM_LAST_SEGMENT);
        if (mateLocation != 0xFFFFFFFF) {
            const Genome::Piece *piece = genome->getPieceAtLocation(mateLocation);
            matePieceName = piece->name;
            matePositionInPiece = mateLocation - piece->beginningOffset + 1;

            if (mateIsRC) {
                flags |= SAM_NEXT_REVERSED;
            }

            if (genomeLocation == 0xFFFFFFFF) {
                //
                // The SAM spec says that for paired reads where exactly one end is unmapped that the unmapped
                // half should just have RNAME and POS copied from the mate.
                //
                pieceName = matePieceName;
                matePieceName = "=";
                positionInPiece = matePositionInPiece;
            }

        } else {
            flags |= SAM_NEXT_UNMAPPED;
            //
            // The mate's unmapped, so point it at us.
            //
            matePieceName = "=";
            matePositionInPiece = positionInPiece;
        }

        if (genomeLocation != 0xffffffff && mateLocation != 0xffffffff) {
            flags |= SAM_ALL_ALIGNED;
            // Also compute the length of the whole paired-end string whose ends we saw. This is slightly
            // tricky because (a) we may have clipped some bases before/after each end and (b) we need to
            // give a signed result based on whether our read is first or second in the pair.
            _int64 myStart = genomeLocation - basesClippedBefore;
            _int64 myEnd = genomeLocation + clippedLength + basesClippedAfter;
            _int64 mateBasesClippedBefore = mate->getFrontClippedLength();
            _int64 mateBasesClippedAfter = mate->getUnclippedLength() - mate->getDataLength() - mateBasesClippedBefore;
            _int64 mateStart = mateLocation - (mateIsRC ? mateBasesClippedAfter : mateBasesClippedBefore);
            _int64 mateEnd = mateLocation + mate->getDataLength() + (!mateIsRC ? mateBasesClippedAfter : mateBasesClippedBefore);
			if (pieceName == matePieceName) { // pointer (not value) comparison, but that's OK.
				if (myStart < mateStart) {
					templateLength = mateEnd - myStart;
				} else {
					templateLength = -(myEnd - mateStart);
				}
 			} // otherwise leave TLEN as zero.
        }

        if (pieceName == matePieceName) {
            matePieceName = "=";     // SAM Spec says to do this when they're equal (and not *, which won't happen because this is a pointer, not string, compare)
        }
    }

    // Write the SAM entry, which requires the following fields:
    //
    // 1. QNAME: Query name of the read or the read pair
    // 2. FLAG: Bitwise flag (pairing, strand, mate strand, etc.)
    // 3. RNAME: Reference sequence name
    // 4. POS: 1-Based leftmost position of clipped alignment
    // 5. MAPQ: Mapping quality (Phred-scaled)
    // 6. CIGAR: Extended CIGAR string (operations: MIDNSHP)
    // 7. MRNM: Mate reference name (‘=’ if same as RNAME)
    // 8. MPOS: 1-based leftmost mate position
    // 9. ISIZE: Inferred insert size
    // 10. SEQQuery: Sequence on the same strand as the reference
    // 11. QUAL: Query quality (ASCII-33=Phred base quality)    

    //
    // Some FASTQ files have spaces in their ID strings, which is illegal in SAM.  Just truncate them at the space.
    //
    const char *firstSpace = strnchr(read->getId(),' ',qnameLen);
    if (NULL != firstSpace) {
        qnameLen = (unsigned)(firstSpace - read->getId());
    }

    const int nmStringSize = 30;// Big enough that it won't buffer overflow regardless of the value of editDistance
    char nmString[nmStringSize];  
    if (editDistance >= 0) {
        snprintf(nmString, nmStringSize, "\tNM:i:%d",editDistance);
    } else {
        nmString[0] = '\0';
    }

    int charsInString = snprintf(buffer, bufferSpace, "%.*s\t%d\t%s\t%u\t%d\t%s\t%s\t%u\t%lld\t%.*s\t%.*s\tPG:Z:SNAP\tRG:Z:FASTQ%s\n",
        qnameLen, read->getId(),
        flags,
        pieceName,
        positionInPiece,
        mapQuality,
        cigar,
        matePieceName,
        matePositionInPiece,
        templateLength,
        fullLength, data,
        fullLength, quality,
        nmString);

    if (charsInString > bufferSpace) {
        //
        // Out of buffer space.
        //
        return false;
    } else if (charsInString == bufferSpace) {
      buffer[bufferSpace-1] = '\n'; // overwrite trailing null with newline
    }


    if (NULL != spaceUsed) {
        *spaceUsed = charsInString;
    }
    return true;
}

ThreadSAMWriter::ThreadSAMWriter(size_t i_bufferSize, bool i_useM)
    : remainingBufferSpace(i_bufferSize), bufferBeingCreated(0), bufferSize(i_bufferSize), useM(i_useM)
{
    buffer[0] = NULL;
    buffer[1] = NULL;
    writer[0] = NULL;
    writer[1] = NULL;
    genome = NULL;
    nextWriteOffset = NULL;
}

    bool
ThreadSAMWriter::initialize(AsyncFile* file, const Genome *i_genome, volatile _int64 *i_nextWriteOffset)
{
    genome = i_genome;
    nextWriteOffset = i_nextWriteOffset;
    buffer[0] = (char *)BigAlloc(bufferSize);
    buffer[1] = (char *)BigAlloc(bufferSize);
    writer[0] = file->getWriter();
    writer[1] = file->getWriter();

    if (NULL == buffer[0] || NULL == buffer[1] || NULL == writer[0] || NULL == writer[1]) {
        fprintf(stderr,"ThreadSAMWriter: failed to initialize\n");
        return false;
    }
    return true;
}
    
ThreadSAMWriter::~ThreadSAMWriter()
{
    BigDealloc(buffer[0]);
    BigDealloc(buffer[1]);
    delete writer[0];
    delete writer[1];
}

    bool
ThreadSAMWriter::close()
{
    if (remainingBufferSpace != bufferSize) {
        if (!startIo()) {
            fprintf(stderr,"WindowsSAMWriter::close(): startIo failed\n");
            return false;
        }

        if (!waitForIoCompletion()) {
            fprintf(stderr,"WindowsSAMWriter::close(): waitForIoCompletion failed\n");
            return false;
        }
    }
    bool ok = writer[0]->close();
    ok &= writer[1]->close();
    if (! ok) {
        fprintf(stderr, "WindowsSAMWriter::close closing writer failed\n");
    }

    return true;
}
    
    bool
ThreadSAMWriter::write(Read *read, AlignmentResult result, unsigned genomeLocation, bool isRC)
{
    size_t sizeUsed;
    if (!generateSAMText(read, result, genomeLocation, isRC, useM, false, true, NULL, UnknownAlignment, 0, false, genome, &lv,
            buffer[bufferBeingCreated] + bufferSize - remainingBufferSpace, remainingBufferSpace, &sizeUsed)) {

        if (!startIo()) {
            return false;
        }

        if (!generateSAMText(read, result, genomeLocation, isRC, useM, false, true, NULL, UnknownAlignment, 0, false, genome, &lv,
                buffer[bufferBeingCreated] + bufferSize - remainingBufferSpace,remainingBufferSpace, &sizeUsed)) {

            fprintf(stderr,"WindowsSAMWriter: create SAM string into fresh buffer failed\n");
            return false;
        }
    }
    size_t bufferOffset = bufferSize - remainingBufferSpace;
    remainingBufferSpace -= sizeUsed;
    afterWrite(result != NotFound ? genomeLocation : UINT32_MAX, bufferOffset, (unsigned)sizeUsed);
    return true;
}

    bool
ThreadSAMWriter::writePair(Read *read0, Read *read1, PairedAlignmentResult *result)
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

    bool writesFit = generateSAMText(reads[first], result->status[first], result->location[first], result->isRC[first], useM, true, true,
                                     reads[second], result->status[second], result->location[second], result->isRC[second],
                        genome, &lv, buffer[bufferBeingCreated] + bufferSize - remainingBufferSpace,remainingBufferSpace,&sizeUsed[first],
                        idLengths[first]);

    if (writesFit) {
        writesFit = generateSAMText(reads[second], result->status[second], result->location[second], result->isRC[second], useM, true, false,
                                    reads[first], result->status[first], result->location[first], result->isRC[first],
                        genome, &lv, buffer[bufferBeingCreated] + bufferSize - remainingBufferSpace + sizeUsed[first],remainingBufferSpace-sizeUsed[first],&sizeUsed[second],
                        idLengths[second]);
    }

    if (!writesFit) {
        if (!startIo()) {
            return false;
        }

        if (!generateSAMText(reads[first], result->status[first], result->location[first], result->isRC[first], useM, true, true,
                             reads[second], result->status[second], result->location[second], result->isRC[second],
                genome, &lv, buffer[bufferBeingCreated] + bufferSize - remainingBufferSpace,remainingBufferSpace,&sizeUsed[first],
                idLengths[first]) ||
            !generateSAMText(reads[second], result->status[second] ,result->location[second], result->isRC[second], useM, true, false,
                             reads[first], result->status[first], result->location[first], result->isRC[first],
                genome, &lv, buffer[bufferBeingCreated] + bufferSize - remainingBufferSpace + sizeUsed[first],remainingBufferSpace-sizeUsed[first],&sizeUsed[second],
                idLengths[second])) {


            fprintf(stderr,"WindowsSAMWriter: create SAM string into fresh buffer failed\n");
            return false;
        }
    }

    size_t bufferOffset[2] = {bufferSize - remainingBufferSpace, bufferSize - remainingBufferSpace + sizeUsed[first]};
    remainingBufferSpace -= (sizeUsed[0] + sizeUsed[1]);
    //
    // The strange code that determines the sort key (which uses the coordinate of the mate for unmapped reads) is because we list unmapped reads
    // with mapped mates at their mates' location so they sort together.  If both halves are unmapped, then  
    afterWrite(result->status[first] != NotFound ? result->location[first] : ((result->status[second] != NotFound) ? result->location[second] : UINT32_MAX), bufferOffset[0], (unsigned)sizeUsed[first]);
    afterWrite(result->status[second] != NotFound ? result->location[second] : ((result->status[first] != NotFound) ? result->location[first] : UINT32_MAX), bufferOffset[1], (unsigned)sizeUsed[second]);
    return true;
}

    ParallelSAMWriter *
ParallelSAMWriter::create(
    const char		*fileName,
    const Genome	*genome,
    unsigned		 nThreads,
    bool			 sort,
    size_t			 sortBufferMemory,
	bool			 useM,
    int              argc,
    const char     **argv,
    const char      *version,
    const char      *rgLine) 
{
    ParallelSAMWriter *parallelWriter = sort
        ? new SortedParallelSAMWriter(sortBufferMemory, useM, argc, argv, version, rgLine)
        : new ParallelSAMWriter(useM,argc,argv,version, rgLine);
    if (!parallelWriter->initialize(fileName, genome, nThreads, sort)) {
        fprintf(stderr, "unable to initialize parallel SAM writer\n");
        delete parallelWriter;
        return NULL;
    }
    return parallelWriter;
}

    bool

ParallelSAMWriter::initialize(const char *fileName, const Genome *genome, unsigned i_nThreads, bool sorted)
{
    file = AsyncFile::open(fileName, true);
    if (NULL == file) {
        fprintf(stderr,"Unable to create SAM file '%s'\n",fileName);
        return false;
    }

    char *headerBuffer = new char[SAMWriter::HEADER_BUFFER_SIZE];
    size_t headerActualSize;

    if (!SAMWriter::generateHeader(genome,headerBuffer,SAMWriter::HEADER_BUFFER_SIZE,&headerActualSize, sorted, argc, argv, version, rgLine)) {
        fprintf(stderr,"WindowsParallelSAMWriter: unable to generate SAM header.\n");
        delete[] headerBuffer;
        return false;
    }

    AsyncFile::Writer* hwriter = file->getWriter();

    if (NULL == hwriter) {
        fprintf(stderr,"ParallelSAMWriter: unable to create writer\n");
        delete[] headerBuffer;
        return false;
    }

    size_t bytesWritten;
    if (! hwriter->beginWrite(headerBuffer, headerActualSize, 0, &bytesWritten)) {
        fprintf(stderr,"ParallelSAMWriter: unable to write header to file\n");
        delete[] headerBuffer;
        return false;
    }

    if (! hwriter->waitForCompletion()) {
        fprintf(stderr,"ParallelSAMWriter: failed to complete\n");
        delete[] headerBuffer;
        return false;
    }
    if (! hwriter->close()) {
        fprintf(stderr, "ParallelSAMWriter: failed to close\n");
    }
    delete hwriter;
    delete[] headerBuffer;

    nextWriteOffset = headerActualSize;

    nThreads = i_nThreads;
    writer = new ThreadSAMWriter *[nThreads];

    if (!createThreadWriters(genome)) {
        fprintf(stderr,"Unable to create SAM writer.\n");
        return false;
    }

    return true;
}

    bool
ParallelSAMWriter::createThreadWriters(const Genome* genome)
{
    bool worked = true;
    for (int i = 0; i < nThreads; i++) {
        writer[i] = new ThreadSAMWriter(UnsortedBufferSize, useM);
        worked &= writer[i]->initialize(file, genome, &nextWriteOffset);
    }
    return worked;
}

ParallelSAMWriter::~ParallelSAMWriter()
{
    delete [] writer;
}

    SAMWriter*
ParallelSAMWriter::getWriterForThread(
    int whichThread)
{
    _ASSERT(whichThread < nThreads);
    return writer[whichThread];
}

    bool
ParallelSAMWriter::close()
{
    for (int i = 0; i < nThreads; i++) {
        // writers were already closed by each thread so they could flush in parallel
        delete writer[i];
        writer[i] = NULL;
    }
    if (! file->close()) {
        fprintf(stderr, "ParallelSAMWriter::close file close failed\n");
    }
    delete file;
    file = NULL;
    return true;
}

    bool
ThreadSAMWriter::startIo()
{
    //
    // It didn't fit in the buffer.  Start writing it.
    //
    _int64 writeOffset = InterlockedAdd64AndReturnNewValue(nextWriteOffset,bufferSize - remainingBufferSpace) - (bufferSize - remainingBufferSpace);
    if (!beforeFlush(writeOffset, bufferSize - remainingBufferSpace)) {
        fprintf(stderr, "ThreadSAMWriter: beforeFlush failed\n");
        return false;
    }
    if (!writer[bufferBeingCreated]->beginWrite(buffer[bufferBeingCreated], bufferSize - remainingBufferSpace, writeOffset, NULL)) {
        fprintf(stderr,"ThreadSAMWriter: WriteFile failed\n");
        return false;
    }

    //
    // If necessary, wait for the other buffer to finish writing.
    //
    if (!waitForIoCompletion()) {
        fprintf(stderr,"ThreadSAMWriter: waitForIoCompletion failed\n");
        return false;
    }
    bufferBeingCreated = 1 - bufferBeingCreated;
    remainingBufferSpace = bufferSize;

    return true;
}

    bool
ThreadSAMWriter::waitForIoCompletion()
{
    return writer[1 - bufferBeingCreated]->waitForCompletion();
}

SortBlock::SortBlock()
    : entries(0), fileOffset(0), fileBytes(0), index(0), reader()
{
}

SortBlock::SortBlock(size_t capacity)
    : entries((int)capacity), fileOffset(0), fileBytes(0), index(0), reader()
{
}

SortBlock::SortBlock(
    SortBlock& other)
{
    entries = other.entries;
    fileOffset = other.fileOffset;
    fileBytes = other.fileBytes;
    index = other.index;
    reader = other.reader;
}

    void
SortBlock::operator=(
    SortBlock& other)
{
    entries = other.entries;
    fileOffset = other.fileOffset;
    fileBytes = other.fileBytes;
    index = other.index;
    index = other.index;
    reader = other.reader;
}

SortedThreadSAMWriter::SortedThreadSAMWriter(size_t i_bufferSize, bool useM)
    : ThreadSAMWriter(i_bufferSize, useM),
        parent(NULL),
        largest(1000),
        locations(1000)
{
}

SortedThreadSAMWriter::~SortedThreadSAMWriter()
{
}

    bool
SortedThreadSAMWriter::initialize(
    SortedParallelSAMWriter* i_parent,
    const Genome* i_genome)
{
    parent = i_parent;
    return ThreadSAMWriter::initialize(parent->file, i_genome, &parent->nextWriteOffset);
}

    void
SortedThreadSAMWriter::afterWrite(
    unsigned location,
    size_t bufferOffset,
    unsigned length)
{
    SortEntry entry(bufferOffset, length, location);
    locations.entries.push_back(entry);
}

    bool
SortedThreadSAMWriter::beforeFlush(_int64 fileOffset, _int64 length)
{
    // sort buffered reads by location for later merge sort
    std::sort(locations.entries.begin(), locations.entries.end(), SortEntry::comparator);
    
    // wait for IO of other buffer to finish since we're going to sort into there
    if (! waitForIoCompletion()) {
        fprintf(stderr, "SortedThreadSAMWriter waitForIoCompletion failed\n");
        return false;
    }
    
    // copy into other buffer in sorted order & switch buffers
    unsigned target = 0;
    for (VariableSizeVector<SortEntry>::iterator i = locations.entries.begin(); i != locations.entries.end(); i++) {
        memcpy(buffer[1 - bufferBeingCreated] + target, buffer[bufferBeingCreated] + i->offset, i->length);
        i->offset = fileOffset + target;
        target += i->length;
    }
    bufferBeingCreated = 1 - bufferBeingCreated;
    
    // remember offsets for full-file sort, get a new vector of appropriate size
    locations.fileOffset = fileOffset;
    locations.fileBytes = length;
    parent->addLocations(locations);
    if (locations.entries.size() > largest) {
        largest = (locations.entries.size() * 6) / 5; // grow by 20% if it exceeds prior max
    }
    locations.entries.reserve(largest);

    return true;
}

    bool
SortedParallelSAMWriter::initialize(
    const char *fileName,
    const Genome *genome,
    unsigned i_nThreads,
    bool sorted)
{
    InitializeExclusiveLock(&lock);
    sortedFile = fileName;
    tempFile = (char*) malloc(strlen(fileName) + 5);
    strcpy(tempFile, fileName);
    strcat(tempFile, ".tmp");
    bool ok = ParallelSAMWriter::initialize(tempFile, genome, i_nThreads, sorted);
    headerSize = nextWriteOffset;
    return ok;
}
    
    bool
SortedParallelSAMWriter::createThreadWriters(const Genome* genome)
{
    size_t bufferSize = totalMemory / nThreads / 2;
    bool worked = true;
    for (int i = 0; i < nThreads; i++) {
        SortedThreadSAMWriter* w = new SortedThreadSAMWriter(bufferSize, useM);
        writer[i] = w;
        worked &= w->initialize(this, genome);
    }
    return worked;
}

class SortContext : public TaskContextBase
{
public:

    void initializeThread() {}

    void runThread();

    void finishThread(SortContext* parent) {}

    char*           source;
    size_t          bufferSize;
    unsigned        blockSize;
    RangeSplitter*  range;
    const size_t*   blockOffsets;
    size_t          entryCount;
    SortEntry*      entries;
    AsyncFile*      file;
};

void SortContext::runThread()
{
    // allocate a pair of buffers and async writers
    AsyncFile::Writer* writers[2] = {file->getWriter(), file->getWriter()};
    char* buffers[2] = {(char*) BigAlloc(bufferSize), (char*) BigAlloc(bufferSize)};
    if (buffers[0] == NULL || buffers[1] == NULL || writers[0] == NULL || writers[1] == NULL) {
        fprintf(stderr, "could not allocate write buffers\n");
        return;
    }
    int writingBuffer = 0;

    // copy blocks of source into target at desired location
    _int64 rangeStart, rangeLength;
    while (range->getNextRange(&rangeStart, &rangeLength)) {
        size_t targetOffset = blockOffsets[rangeStart];
        size_t bufferOffset = 0;
        unsigned end = min((unsigned) entryCount, (unsigned) (rangeStart + rangeLength) * blockSize );
        for (size_t read = rangeStart * blockSize; ; read++) {
            SortEntry* entry;
            if (read == end || bufferOffset + (entry = &entries[read])->length > bufferSize) {
                writers[1 - writingBuffer]->waitForCompletion();
                writers[writingBuffer]->beginWrite(buffers[writingBuffer], bufferOffset, targetOffset, NULL);
                writingBuffer = 1 - writingBuffer;
                targetOffset += bufferOffset;
                bufferOffset = 0;
                if (read == end) {
                    break;
                }
            }
            memcpy(buffers[writingBuffer] + bufferOffset, source + entry->offset, entry->length);
            bufferOffset += entry->length;
        }
        writers[1 - writingBuffer]->waitForCompletion();
    }
    delete writers[0];
    delete writers[1];
    BigDealloc(buffers[0]);
    BigDealloc(buffers[1]);
}

    bool
SortedParallelSAMWriter::close()
{
    if (! ParallelSAMWriter::close()) {
        return false;
    }
    DestroyExclusiveLock(&lock);

#ifdef _MSC_VER
    return mergeSort();
#else
#ifdef __linux__
    return mergeSort();
#else
    return memoryMappedSort(); // async io not supported on OS X
#endif
#endif
}

    bool
SortedParallelSAMWriter::mergeSort()
{
    // merge sort from temp file into sorted file
#if USE_DEVTEAM_OPTIONS
    printf("sorting...");
    _int64 start = timeInMillis();
#endif

    // first replace offset in entries with block index, and sort them all in one large array
    size_t total = 0;
    for (VariableSizeVector<SortBlock>::iterator i = locations.begin(); i != locations.end(); i++) {
        total += i->entries.size();
    }
    SortEntry* entries = (SortEntry*) BigAlloc(total * sizeof(SortEntry));
    size_t offset = 0;
    for (VariableSizeVector<SortBlock>::iterator i = locations.begin(); i != locations.end(); i++) {
        unsigned n = i->entries.size();
        size_t blockIndex = i - locations.begin();
        memcpy(entries + offset, i->entries.begin(), n * (size_t) sizeof(SortEntry));
        for (unsigned j = 0; j < n; j++) {
            entries[offset + j].offset = blockIndex;
        }
        offset += n;
    }
    std::stable_sort(entries, entries + total, SortEntry::comparator);
#if USE_DEVTEAM_OPTIONS
    printf(" %ld s\nwriting sorted reads...", (timeInMillis() - start) / 1000);
    start = timeInMillis();
#endif

    // setup - open all files, read first block, begin read for second
    AsyncFile* temp = AsyncFile::open(tempFile, false);
    const size_t ReadBufferSize = UnsortedBufferSize;
    char* buffers = (char*) BigAlloc(ReadBufferSize * 2 * locations.size());
    unsigned j = 0;
#if USE_DEVTEAM_OPTIONS
    if (timeInMillis() - start > 1000) {
        printf(" (allocated %lld Mb in %lld s)",
            ReadBufferSize * 2 * locations.size() / (2 << 20), (timeInMillis() - start) / 1000);
    }
#endif
    for (VariableSizeVector<SortBlock>::iterator i = locations.begin(); i != locations.end(); i++, j++) {
        i->reader.open(temp, i->fileOffset, i->fileBytes, ReadBufferSize, true,
            buffers + j * 2 * ReadBufferSize, buffers + (j * 2 + 1) * ReadBufferSize);
    }
    for (VariableSizeVector<SortBlock>::iterator i = locations.begin(); i != locations.end(); i++) {
        i->reader.endOpen();
    }

    // set up double-buffered output
    AsyncFile* sorted = AsyncFile::open(sortedFile, true);
	if (NULL == sorted) {
		return false;
	}
    BufferedAsyncWriter writer;
    const size_t WriteBufferSize = UnsortedBufferSize;
    if (! writer.open(sorted, WriteBufferSize)) {
        fprintf(stderr, "open sorted file for write failed\n");
        return false;
    }

    // write out header
    if (headerSize > 0) {
        void* hbuf = BigAlloc(headerSize);
        AsyncFile::Reader* hread = temp->getReader();
        bool ok = hread->beginRead(hbuf, headerSize, 0, NULL);
        ok &= hread->waitForCompletion();
        ok &= hread->close();
        delete hread;
        if (! ok ) {
            fprintf(stderr, "read header failed\n");
            return false;
        }
        if (! writer.write(hbuf, headerSize)) {
            fprintf(stderr, "write header failed\n");
            return false;
        }
        BigDealloc(hbuf);
    }

    // merge input blocks into output using pre-sorted list
    for (size_t i = 0; i < total; i++) {
        SortEntry* entry = &entries[i];
        void* buf = writer.forWrite(entry->length);
        if (buf == NULL || ! locations[entry->offset].reader.read(buf, entry->length)) {
            fprintf(stderr, "merge %s failed\n", buf ? "write" : "read");
            return false;
        }
    }

    // close everything
    BigDealloc(entries);
    bool ok = writer.close();
    ok &= sorted->close();
    delete sorted;
    _int64 readWait = 0;
    for (VariableSizeVector<SortBlock>::iterator i = locations.begin(); i != locations.end(); i++) {
        ok &= i->reader.close();
        readWait += i->reader.getWaitTimeInMillis();
    }
    ok &= temp->close();
    delete temp;
    BigDealloc(buffers);
    if (! ok) {
        printf("files did not close properly\n");
    }
    if (! DeleteSingleFile(tempFile)) {
        printf("warning: failure deleting temp file %s\n", tempFile);
    }

#if USE_DEVTEAM_OPTIONS
    printf(" %u reads in %u blocks, %lld s (%lld s read wait, %lld s write wait)\n",
        total, locations.size(), (timeInMillis() - start)/1000, readWait/1000, writer.getWaitTimeInMillis()/1000);
#endif
#ifdef PROFILE_BIGALLOC
    PrintAllocProfile();
#endif
    return true;
}

    bool
SortedParallelSAMWriter::memoryMappedSort()
{
    // sort by location and copy from temp to final file in sorted order
    // because of buffer sorting this should be a merge-sort rather than random-access
    printf("sorting...");
    _int64 start = timeInMillis();
    // todo: use in-memory merge sort instead of standard sort?
    // copy into single buffer before sorting
    size_t entryCount = 0;
    for (VariableSizeVector<SortBlock>::iterator i = locations.begin(); i != locations.end(); i++) {
        entryCount += i->entries.size();
    }
    SortEntry* entries = (SortEntry*) BigAlloc(entryCount * sizeof(SortEntry));
    size_t offset = 0;
    for (VariableSizeVector<SortBlock>::iterator i = locations.begin(); i != locations.end(); i++) {
        memcpy(entries + offset, i->entries.begin(), i->entries.size() * sizeof(SortEntry));
        offset += i->entries.size();
    }
    std::sort(entries, entries + entryCount, SortEntry::comparator);
    printf(" %ld s\n", (timeInMillis() - start) / 1000);

    printf("writing sorted reads...");
    start = timeInMillis();
    void* p;
    MemoryMappedFile* map = OpenMemoryMappedFile(tempFile, 0, QueryFileSize(tempFile), &p, false, true);
    if (map == NULL) {
        fprintf(stderr, "Could not map temporary file\n");
        // todo: just use unsorted file?
        return false;
    }

    // divide into blocks and figure out offset from base
    const unsigned blockSize = 1000;
    size_t blockCount = (entryCount + blockSize - 1) / blockSize;
    size_t* blockOffsets = new size_t[blockCount];
    offset = headerSize;
    for (unsigned block = 0; ; block++) {
        blockOffsets[block] = offset;
        if (block == blockCount - 1) {
            break;
        }
        for (unsigned i = 0; i < blockSize; i++) {
            offset += entries[block * blockSize + i].length;
        }
    }

    // setup for async permutation
    RangeSplitter range(blockCount, nThreads);
    SortContext context;
    context.totalThreads = nThreads; // todo: optimize - might be better to use 2x or more
    context.bindToProcessors = false; // disk bound, so processor affinity doesn't matter
    context.source = (char*) p;
    context.blockSize = blockSize;
    context.bufferSize = UnsortedBufferSize; // use smaller write buffers
    context.range = &range;
    context.blockOffsets = blockOffsets;
    context.entries = entries;
    context.entryCount = entryCount;
    context.file = AsyncFile::open(sortedFile, true);
    if (context.file == NULL) {
        fprintf(stderr, "could not open sorted file for write %s\n", sortedFile);
        delete[] blockOffsets;
        return false;
    }
    // write out header
    AsyncFile::Writer* hwrite = context.file->getWriter();
    hwrite->beginWrite(p, headerSize, 0, NULL);
    hwrite->waitForCompletion();
    delete hwrite;
    // run parallel tasks
    ParallelTask<SortContext> task(&context);
    task.run();

    delete [] blockOffsets;
    delete context.file;
    BigDealloc(entries);
    CloseMemoryMappedFile(map);
    if (! DeleteSingleFile(tempFile)) {
        printf("warning: failure deleting temp file %s\n", tempFile);
    }

    printf(" %lld reads, %lld bytes, %ld seconds\n",
        entryCount, offset, (timeInMillis() - start) / 1000);

    return true;
}

    void
SortedParallelSAMWriter::addLocations(
    SortBlock& added)
{
    AcquireExclusiveLock(&lock);
    locations.push_back(added);
    ReleaseExclusiveLock(&lock);
}

    char *
strnchrs(char *str, char charToFind, char charToFind2, size_t maxLen) // Hokey version that looks for either of two chars
{
    for (size_t i = 0; i < maxLen; i++) {
        if (str[i] == charToFind || str[i] == charToFind2) {
            return str + i;
        }
        if (str[i] == 0) {
            return NULL;
        }
    }
    return NULL;
}

    char *
skipToBeyondNextRunOfSpacesAndTabs(char *str, const char *endOfBuffer, size_t *charsUntilFirstSpaceOrTab = NULL)
{
    if (NULL == str) return NULL;

    char *nextChar = str;
    while (nextChar < endOfBuffer && *nextChar != ' ' && *nextChar != '\n' && *nextChar != '\t' && *nextChar != '\r' /* for Windows CRLF text */) {
        nextChar++;
    }

    if (NULL != charsUntilFirstSpaceOrTab) {
        *charsUntilFirstSpaceOrTab = nextChar - str;
    }

    if (nextChar >= endOfBuffer || *nextChar == '\n') {
        return NULL;
    }

    while (nextChar < endOfBuffer && (' ' == *nextChar || '\t' == *nextChar || '\r' == *nextChar)) {
        nextChar++;
    }

    if (nextChar >= endOfBuffer) {
        return NULL;
    }

    return nextChar;
}


    SAMReader *
SAMReader::create(
    const DataSupplier* supplier,
    const char *fileName,
    const Genome *genome, 
    _int64 startingOffset, 
    _int64 amountOfFileToProcess, 
    ReadClippingType clipping)
{
    DataReader* data = supplier->getDataReader(maxLineLen);
    SAMReader *reader = new SAMReader(data, clipping);

    if (!reader->init(fileName, genome, startingOffset, amountOfFileToProcess)) {
        //
        // Probably couldn't open the file.
        //
        delete reader;
        return NULL;
    }
    return reader;
}

SAMReader::SAMReader(
    DataReader* i_data,
    ReadClippingType i_clipping)
    : data(i_data), clipping(i_clipping)
{
}


//
// Implement the ReadReader form of getNextRead, which doesn't include the
// alignment results by simply throwing them away.
//
    bool
SAMReader::getNextRead(Read *readToUpdate)
{
    return getNextRead(readToUpdate, NULL, NULL, NULL, NULL, NULL, NULL);
}

    bool
SAMReader::parseHeader(
    const char *fileName, 
    char *firstLine, 
    char *endOfBuffer, 
    const Genome *genome, 
    _int64 *headerSize)
{
    char *nextLineToProcess = firstLine;

    while (NULL != nextLineToProcess && nextLineToProcess < endOfBuffer && '@' == *nextLineToProcess) {
        if (!strncmp("@SQ",nextLineToProcess,3)) {
            //
            // These lines represent sequences in the reference genome, what are
            // called "pieces" in the Genome class.  (Roughly, chromosomes or major
            // variants like some versions of the MHC genes on chr6; or more
            // particularly the things that come in different FASTA files from the
            // reference assembly).
            //
            // Verify that they actually match what's in our reference genome.
            //

            if (nextLineToProcess + 3 >= endOfBuffer || ' ' != nextLineToProcess[3] && '\t' != nextLineToProcess[3]) {
                fprintf(stderr,"Malformed SAM file '%s' has @SQ without a following space or tab.\n",fileName);
                return false;
            }

            char *snStart = nextLineToProcess + 4;
            while (snStart < endOfBuffer && strncmp(snStart,"SN:",__min(3,endOfBuffer-snStart)) && *snStart != '\n') {
                snStart++;
            }

            if (snStart >= endOfBuffer || *snStart == '\n') {
                fprintf(stderr,"Malformed @SQ line doesn't have 'SN:' in file '%s'\n",fileName);
                return false;
            }

            const size_t pieceNameBufferSize = 512;
            char pieceName[pieceNameBufferSize];
            for (unsigned i = 0; i < pieceNameBufferSize && snStart+3+i < endOfBuffer; i++) {
                if (snStart[3+i] == ' ' || snStart[3+i] == '\t' || snStart[3+i] == '\n') {
                    pieceName[i] = '\0';
                } else {
                    pieceName[i] = snStart[3+i];
                }
            }
            pieceName[pieceNameBufferSize - 1] = '\0';

            //if (!genome->getOffsetOfPiece(pieceName,NULL)) {
            //    fprintf(stderr,"SAM file '%s' contains sequence name '%s' that isn't in the reference genome.\n",
            //                fileName,pieceName);
            //    return false;
            //}
        } else if (!strncmp("@HD",nextLineToProcess,3) || !strncmp("@RG",nextLineToProcess,3) || !strncmp("@PG",nextLineToProcess,3) ||
            !strncmp("@CO",nextLineToProcess,3)) {
            //
            // Ignore these lines.
            //
        } else {
            fprintf(stderr,"Unrecognized header line in SAM file.\n");
            return false;
        }
        nextLineToProcess = strnchr(nextLineToProcess,'\n',endOfBuffer-nextLineToProcess) + 1;
    }

    *headerSize = nextLineToProcess - firstLine;
    return true;
}

    bool
SAMReader::parseLine(char *line, char *endOfBuffer, char *result[], size_t *linelength, size_t fieldLengths[])
{
    *linelength = 0;

    char *next = line;
    char *endOfLine = strnchr(line,'\n',endOfBuffer-line);
    if (NULL == endOfLine) {
        return false;
    }

    //
    // Skip over any leading spaces and tabs
    //
    while (next < endOfLine && (*next == ' ' || *next == '\t')) {
        next++;
    }

    for (unsigned i = 0; i < nSAMFields; i++) {
        if (NULL == next || next >= endOfLine) {
            //
            // Too few fields.
            //
            return false;
        }

        result[i] = next;

        next = skipToBeyondNextRunOfSpacesAndTabs(next,endOfLine,&fieldLengths[i]);
    }

    *linelength =  endOfLine - line + 1;    // +1 skips over the \n
    return true;
}

    void
SAMReader::getReadFromLine(
    const Genome        *genome,
    char                *line, 
    char                *endOfBuffer, 
    Read                *read, 
    AlignmentResult     *alignmentResult,
    unsigned            *genomeLocation, 
    bool                *isRC,
    unsigned            *mapQ,
    size_t              *lineLength,
    unsigned *           flag,
    const char **        cigar,
    ReadClippingType     clipping
    )
{
    char *field[nSAMFields];
    size_t fieldLength[nSAMFields];

    if (!parseLine(line, endOfBuffer, field, lineLength, fieldLength)) {
        fprintf(stderr, "Failed to parse SAM line:\n%.*s\n", lineLength, line);
        exit(1);
    }

    //
    // We have to copy the piece name (RNAME) into its own buffer because the code in Genome expects
    // it to be a null-terminated string, while all we've got is one that's space delimited.
    //
    const size_t pieceNameBufferSize = 512;
    char pieceName[pieceNameBufferSize];

    if (fieldLength[RNAME] >= pieceNameBufferSize) {  // >= because we need a byte for the \0
        fprintf(stderr,"SAMReader: too long an RNAME.  Can't parse.\n");
        exit(1);
    }
    
    memcpy(pieceName,field[RNAME],fieldLength[RNAME]);
    pieceName[fieldLength[RNAME]] = '\0';

    unsigned offsetOfPiece;
    if ('*' != pieceName[0] && !genome->getOffsetOfPiece(pieceName,&offsetOfPiece)) {
        fprintf(stderr,"Unable to find piece '%s' in genome.  SAM file malformed.\n",pieceName);
        //exit(1);
    }

    if (NULL != genomeLocation) {
        unsigned oneBasedOffsetWithinPiece = 0;
        if ('*' != pieceName[0]) {
            //
            // We can't call sscanf directly into the mapped file, becuase it reads to the end of the
            // string even when it's satisfied all of its fields.  Since this can be gigabytes, it's not
            // really good for perf.  Instead, copy the POS field into a local buffer and null terminate it.
            //

            const unsigned posBufferSize = 20;
            char posBuffer[posBufferSize];
            if (fieldLength[POS] >= posBufferSize) {
                fprintf(stderr,"SAMReader: POS field too long.\n");
                exit(1);
            }
            memcpy(posBuffer,field[POS],fieldLength[POS]);
            posBuffer[fieldLength[POS]] = '\0';
            if (0 == sscanf(posBuffer,"%d",&oneBasedOffsetWithinPiece)) {
                fprintf(stderr,"SAMReader: Unable to parse position when it was expected.\n");
                exit(1);
            }
            if (0 == oneBasedOffsetWithinPiece) {
                fprintf(stderr,"SAMReader: Position parsed as 0 when it was expected.\n");
                exit(1);
            }
            *genomeLocation = offsetOfPiece + oneBasedOffsetWithinPiece - 1; // -1 is because our offset is 0 based, while SAM is 1 based.
        } else {
            *genomeLocation = 0xffffffff;
        }
    }

    if (fieldLength[SEQ] != fieldLength[QUAL]) {
        fprintf(stderr,"SAMReader: QUAL string unequal in length to SEQ string.\n");
        exit(1);
    }

    unsigned _flag;
    const size_t flagBufferSize = 20;   // More than enough
    char flagBuffer[flagBufferSize];
    if (fieldLength[FLAG] >= flagBufferSize) {
        fprintf(stderr,"SAMReader: flag field is too long.\n");
        exit(1);
    }
    memcpy(flagBuffer,field[FLAG],fieldLength[FLAG]);
    flagBuffer[fieldLength[FLAG]] = '\0';
    if (1 != sscanf(flagBuffer,"%d",&_flag)) {
        fprintf(stderr,"SAMReader: couldn't parse FLAG field.\n");
        exit(1);
    }

    if (NULL != read) {
        //
        // Clip reads where the quality strings end in '#'
        //
        read->init(field[QNAME],(unsigned)fieldLength[QNAME],field[SEQ],field[QUAL],(unsigned)fieldLength[SEQ]);
        //
        // If this read is RC in the SAM file, we need to reverse it here, since Reads are always the sense that they were as they came
        // out of the base caller.
        //

        if (_flag & SAM_REVERSE_COMPLEMENT) {
            read->becomeRC();
        }
        read->clip(clipping);
    }

    if (NULL != alignmentResult) {
        if (_flag & SAM_UNMAPPED) {
            *alignmentResult = NotFound;
        } else {
            if ('*' == pieceName[0]) {
                fprintf(stderr,"SAMReader: mapped read didn't have RNAME filled in.\n");
                exit(1);
            }
            *alignmentResult = SingleHit;   // NB: This isn't quite right, we should look at MAPQ.
        }
    }

    if (NULL != isRC) {
        *isRC = (_flag & SAM_REVERSE_COMPLEMENT) ? true : false;
    }

    if (NULL != mapQ) {
        *mapQ = atoi(field[MAPQ]);
        if (*mapQ > 255) {
            fprintf(stderr,"SAMReader: MAPQ field has bogus value\n");
            exit(1);
        }
    }

    if (NULL != flag) {
        *flag = _flag;
    }

    if (NULL != cigar) {
        *cigar = field[CIGAR];
    }
}

    bool
SAMReader::init(
    const char *fileName,
    const Genome *i_genome,
    _int64 startingOffset,
    _int64 amountOfFileToProcess)
{
    genome = i_genome;
    
    //
    // Read and parse the header specially.
    //
    if (! data->init(fileName)) {
        return false;
    }
    headerSize = 1024 * 1024; // 1M header max
    char* buffer = data->readHeader(&headerSize);
    if (!parseHeader(fileName, buffer, buffer + headerSize, genome, &headerSize)) {
        fprintf(stderr,"SAMReader: failed to parse header on '%s'\n",fileName);
        return false;
    }

    reinit(startingOffset,amountOfFileToProcess);

    return true;
}

    void
SAMReader::reinit(_int64 startingOffset, _int64 amountOfFileToProcess)
{
    _ASSERT(0 != headerSize);  // Must call init() before reinit()
    data->reinit(
        max(headerSize, startingOffset) - 1,  // -1 is to point at the previous newline so we don't skip the first line.
        startingOffset + amountOfFileToProcess);
    char* buffer;
    _int64 validBytes;
    data->getData(&buffer, &validBytes);
    char *firstNewline = strnchr(buffer,'\n',validBytes);
    if (NULL == firstNewline) {
        return;
    }

    //
    // Parse the new first line.  If it's got SAM_MULTI_SEGMENT set and not SAM_FIRST_SEGMENT, then skip it and go with the next one.
    //
    Read read;
    AlignmentResult alignmentResult;
    unsigned genomeLocation;
    bool isRC;
    unsigned mapQ;
    size_t lineLength;
    unsigned flag;
 
    getReadFromLine(genome,firstNewline+1,buffer + validBytes,
                    &read,&alignmentResult,&genomeLocation,&isRC,&mapQ,&lineLength,
                    &flag,NULL,clipping);

    if ((flag & SAM_MULTI_SEGMENT) && !(flag & SAM_FIRST_SEGMENT)) {
        //
        // Skip this line.
        //
        data->advance((unsigned)(firstNewline + lineLength - buffer + 1)); // +1 skips over the newline.
    } else {
        data->advance((unsigned)(firstNewline - buffer + 1)); // +1 skips over the newline.
    }
}

    bool
SAMReader::getNextRead(
    Read *read,
    AlignmentResult *alignmentResult,
    unsigned *genomeLocation,
    bool *isRC,
    unsigned *mapQ, 
    unsigned *flag,
    bool ignoreEndOfRange,
    const char **cigar)
{
    char* buffer;
    _int64 bytes;
    if (! data->getData(&buffer, &bytes)) {
        return false;
    }
    char *newLine = strnchr(buffer, '\n', bytes);
    if (NULL == newLine) {
        //
        // There is no newline, so the line crosses the end of the buffer.
        // This should never happen since underlying reader manages overflow between chunks.
        //
        fprintf(stderr,"SAM file has too long a line, or doesn't end with a newline!  Failing.  fileOffset = %lld\n", data->getFileOffset());
        exit(1);
    }

    size_t lineLength;
    getReadFromLine(genome,buffer,buffer + bytes,read,alignmentResult,genomeLocation,isRC,mapQ,&lineLength,flag,cigar,clipping);
    read->setBatch(data->getBatch());
    data->advance((newLine + 1) - buffer);

    return true;
}

    ReadSupplierGenerator *
SAMReader::createReadSupplierGenerator(const char *fileName, int numThreads, const Genome *genome, ReadClippingType clipping)
{
    //
    // single-ended SAM files always can be read with the range splitter.
    //
    RangeSplitter *splitter = new RangeSplitter(QueryFileSize(fileName), numThreads, 100);
    return new RangeSplittingReadSupplierGenerator(fileName, true, clipping, numThreads, genome);
}
    
    PairedReadReader*
SAMReader::createPairedReader(
    const DataSupplier* supplier,
    const char *fileName,
    const Genome *genome,
    _int64 startingOffset,
    _int64 amountOfFileToProcess, 
    ReadClippingType clipping)
{
    SAMReader* reader = SAMReader::create(DataSupplier::Default, fileName, genome, 0, 0, clipping);
    return PairedReadReader::PairMatcher(5000, reader);
}


    PairedReadSupplierGenerator *
SAMReader::createPairedReadSupplierGenerator(const char *fileName, int numThreads, const Genome *genome, ReadClippingType clipping)
{
    //
    // need to use a queue so that pairs can be matched
    //

    ReadSupplierQueue* queue = new ReadSupplierQueue(
        SAMReader::createPairedReader(DataSupplier::Default, fileName, genome, 0, 0, clipping));
    queue->startReaders();
    return queue;
}