/*++

Module Name:

    RangeSplitter.cpp

Abstract:

    Code to split a range into pieces for multiple cores to process.  It's designed
    to handle cores that proceed at varying rates.

Authors:

    Bill Bolosky, 2011

Environment:

    User mode service.

Revision History:

    Pulled out of cSNAP.cpp to make it useful for various versions of the aligner.
    Generalized from FileSplitter to ranges, e.g. for scanning the genome /ravip/5/2012/

--*/

#include "stdafx.h"
#include "RangeSplitter.h"
#include "SAM.h"
#include "FASTQ.h"

using std::max;
using std::min;


RangeSplitter::RangeSplitter(_int64 rangeEnd_, int numThreads_, unsigned divisionSize_, _int64 rangeBegin_, unsigned minMillis_, unsigned minRangeSize_)
{
    numThreads = numThreads_;
    rangeEnd = rangeEnd_;
    rangeBegin = rangeBegin_;
    position = rangeBegin_;
    startTime = 0;                       // We'll initialize it when getNextRange() is first called
    divisionSize = divisionSize_;
    minMillis = minMillis_;
    minRangeSize = minRangeSize_;
}


bool RangeSplitter::getNextRange(_int64 *rangeStart, _int64 *rangeLength)
{
    // If there are multiple threads, start each of them off with (rangeEnd / divionSize / numThreads),
    // and then keep giving everyone 1 / (divisionSize * numThreads) of the remaining data or the amount
    // of units processed per thread in minMillis ms, whichever is bigger.
    // If there's just one thread, we give it the whole range at the beginning.

    if (startTime == 0) {
        // There's a possible "race" here if multiple threads start at the same time, but that's
        // actually OK; we just want a rough idea of when we started the processing.
        startTime = timeInMillis();
    }

    _int64 amountLeft = rangeEnd - position;
    if (amountLeft <= 0) {
        return false;
    }

    _int64 amountToTake;
    if (numThreads == 1) {
		amountToTake = rangeEnd;
    } else if (amountLeft >= (rangeEnd - rangeBegin) / divisionSize) {
        amountToTake = (rangeEnd - rangeBegin) / divisionSize / numThreads;
	    if (amountToTake == 0) {
			amountToTake = amountLeft;
	    } else {
			amountToTake = min(amountLeft, max(amountToTake, (_int64) minRangeSize));
		}
    } else {
        // Figure out units processed in minMillis ms per thread (keeping in mind we ran numThreads in total).
        _int64 unitsInMinms = (position - rangeBegin) * minMillis / max((_int64) (timeInMillis() - startTime) * numThreads, (_int64) 1);
        amountToTake = max(amountLeft / divisionSize / numThreads, unitsInMinms);
        amountToTake = max(amountToTake, (_int64) minRangeSize);  // Avoid getting really tiny amounts at the end.
    }

    _ASSERT(amountToTake > 0);
	_int64 oldPosition = position; // for debugging
    _int64 startOffset = InterlockedAdd64AndReturnNewValue(&position, amountToTake) - amountToTake;
	_ASSERT(position >= rangeBegin);
    if (startOffset >= rangeEnd) {
        // No work left to allocate.
        return false;
    }

    // Don't run past EOF if there was a race above (threads looking at amountLeft at the same time).
    amountToTake = min(amountToTake, rangeEnd - startOffset);
    _ASSERT(amountToTake > 0);

    *rangeStart = startOffset;
    *rangeLength = amountToTake;
	_ASSERT(startOffset >= rangeBegin && startOffset + amountToTake <= rangeEnd);
    return true;
}

RangeSplittingReadSupplierGenerator::RangeSplittingReadSupplierGenerator(
    const char *i_fileName,
    bool i_isSAM, 
    unsigned numThreads,
    const ReaderContext& i_context)
    : isSAM(i_isSAM), context(i_context)
{
    fileName = new char[strlen(i_fileName) + 1];
    strcpy(fileName, i_fileName);
	splitter = new RangeSplitter(QueryFileSize(fileName),numThreads, 5, 0, 200, 10*MAX_READ_LENGTH);
}

ReadSupplier *
RangeSplittingReadSupplierGenerator::generateNewReadSupplier()
{
    _int64 rangeStart, rangeLength;
    if (!splitter->getNextRange(&rangeStart, &rangeLength)) {
        return NULL;
    }

    ReadReader *underlyingReader;
    // todo: implement layered factory model
    if (isSAM) {
        underlyingReader = SAMReader::create(DataSupplier::Default[true], fileName, context, rangeStart, rangeLength);
    } else {
        underlyingReader = FASTQReader::create(DataSupplier::Default[true], fileName, rangeStart, rangeLength, context);
    }
    return new RangeSplittingReadSupplier(splitter,underlyingReader);
}

RangeSplittingReadSupplier::~RangeSplittingReadSupplier()
{
}


    Read * 
RangeSplittingReadSupplier::getNextRead()
{
    if (underlyingReader->getNextRead(&read)) {
        return &read;
    }

    _int64 rangeStart, rangeLength;
    if (!splitter->getNextRange(&rangeStart, &rangeLength)) {
        return NULL;
    }
 
    underlyingReader->reinit(rangeStart,rangeLength);
    if (!underlyingReader->getNextRead(&read)) {
        return NULL;
    }
    return &read;
}

RangeSplittingPairedReadSupplier::~RangeSplittingPairedReadSupplier()
{
}

    bool 
RangeSplittingPairedReadSupplier::getNextReadPair(Read **read1, Read **read2)
{
    *read1 = &internalRead1;
    *read2 = &internalRead2;
    if (underlyingReader->getNextReadPair(&internalRead1,&internalRead2)) {
        return true;
    }

    //
    // We need to clear out the reads, because they may contain references to the buffers in the readers.
    // These buffer reference counts get reset to 0 at reinit time, which causes problems when they're
    // still live in read.
    //

    _int64 rangeStart, rangeLength;
    if (!splitter->getNextRange(&rangeStart, &rangeLength)) {
        return false;
    }
 
    underlyingReader->reinit(rangeStart,rangeLength);
    return underlyingReader->getNextReadPair(&internalRead1, &internalRead2);
}

RangeSplittingPairedReadSupplierGenerator::RangeSplittingPairedReadSupplierGenerator(
    const char *i_fileName1, const char *i_fileName2, bool i_isSAM, unsigned numThreads, const ReaderContext& i_context) :
        isSAM(i_isSAM), context(i_context)
{
    fileName1 = new char[strlen(i_fileName1) + 1];
    strcpy(fileName1, i_fileName1); 

    if (!isSAM) {
        fileName2 = new char[strlen(i_fileName2) + 1];
        strcpy(fileName2, i_fileName2);
    } else {
        fileName2 = NULL;
    }

    splitter = new RangeSplitter(QueryFileSize(fileName1),numThreads);
}

RangeSplittingPairedReadSupplierGenerator::~RangeSplittingPairedReadSupplierGenerator()
{
    delete [] fileName1;
    delete [] fileName2;
    delete splitter;
}

    PairedReadSupplier *
RangeSplittingPairedReadSupplierGenerator::generateNewPairedReadSupplier()
{
    _int64 rangeStart, rangeLength;
    if (!splitter->getNextRange(&rangeStart, &rangeLength)) {
        return NULL;
    }

    PairedReadReader *underlyingReader;
    if (isSAM) {
        underlyingReader = SAMReader::createPairedReader(DataSupplier::Default[true], fileName1, rangeStart, rangeLength, true, context); 
    } else {
        underlyingReader = PairedFASTQReader::create(DataSupplier::Default[true], fileName1, fileName2, rangeStart, rangeLength, context);
    }

    return new RangeSplittingPairedReadSupplier(splitter,underlyingReader);
}

