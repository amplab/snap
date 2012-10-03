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

    _int64 amountToTake;
    if (numThreads == 1) {
        amountToTake = rangeEnd;
    } else if (amountLeft >= rangeEnd / divisionSize) {
        amountToTake = rangeEnd / divisionSize / numThreads;
	    if (amountToTake == 0) {
	      amountToTake = amountLeft;
	    }
    } else {
        // Figure out units processed in minMillis ms per thread (keeping in mind we ran numThreads in total).
        _int64 unitsInMinms = (position - rangeBegin) * minMillis / max((_int64) (timeInMillis() - startTime) * numThreads, (_int64) 1);
        amountToTake = max(amountLeft / divisionSize / numThreads, unitsInMinms);
        amountToTake = max(amountToTake, (_int64) minRangeSize);  // Avoid getting really tiny amounts at the end.
    }

    _ASSERT(amountToTake > 0);

    _int64 startOffset = InterlockedAdd64AndReturnNewValue(&position, amountToTake) - amountToTake;
    if (startOffset >= rangeEnd) {
        // No work left to allocate.
        return false;
    }

    // Don't run past EOF if there was a race above (threads looking at amountLeft at the same time).
    amountToTake = min(amountToTake, rangeEnd - startOffset);
    _ASSERT(amountToTake > 0);

    *rangeStart = startOffset;
    *rangeLength = amountToTake;
    return true;
}
