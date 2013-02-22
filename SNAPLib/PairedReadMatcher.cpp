/*++

Module Name:

    MultiInputReadSupplier.cp

Abstract:

    A read supplier that combines other read suppliers.  It's used when there are muliple input files to process.

Authors:

    Bill Bolosky, November, 2012

Environment:

    User mode service.

Revision History:


--*/

#pragma once
#include "stdafx.h"
#include <map>
#include "Compat.h"
#include "Read.h"
#include "DataReader.h"
#include "FixedSizeMap.h"

using std::pair;
using std::map;
using std::string;

class PairedReadMatcher: public PairedReadReader
{
public:
    PairedReadMatcher(int i_capacity, ReadReader* i_single);

    // PairedReadReader

    virtual ~PairedReadMatcher()
    { delete single; }

    virtual bool getNextReadPair(Read *read1, Read *read2);
    
    virtual void reinit(_int64 startingOffset, _int64 amountOfFileToProcess)
    { single->reinit(startingOffset, amountOfFileToProcess); }

    void releaseBefore(DataBatch batch)
    {
        if (outerReleased < batch) {
            outerReleased = batch;
            checkRelease();
        }
    }

private:

    void checkRelease()
    {
        // get min of pending, outer
        DataBatch release = pendingReleased < outerReleased ? pendingReleased : outerReleased;
        if (innerReleased < release) {
            innerReleased = release;
            single->releaseBefore(release);
        }
    }

    typedef map<string,Read> ReadMap;
    ReadMap pending;
    BatchTracker tracker;
    ReadReader* single;
    const int capacity;
    bool capacityExceeded;
    DataBatch pendingReleased; // max batch that I have released from the buffer
    DataBatch outerReleased; // max batch that has been released by others calling releaseBefore
    DataBatch innerReleased; // max batch sent to single->releaseBefore() = min(pendingReleased, outerReleased)
};

PairedReadMatcher::PairedReadMatcher(
    int i_capacity,
    ReadReader* i_single)
    : single(i_single),
    pending(),
    tracker(128),
    capacity(i_capacity),
    capacityExceeded(false),
    pendingReleased(),
    outerReleased(),
    innerReleased()
{
}

    bool
PairedReadMatcher::getNextReadPair(
    Read *read1,
    Read *read2)
{
    while (true) {
        Read one;
        if (! single->getNextRead(&one)) {
            if (pending.size() > 0) {
                fprintf(stderr, "PairedReadMatcher warning: %d unmated reads in file\n", pending.size());
            }
            return false;
        }
        // build key for pending read table, removing /1 or /2 at end
        const char* id = one.getId();
        unsigned idLength = one.getIdLength();
        if (idLength > 2 && id[idLength-2] == '/' && (id[idLength-1] == '1' || id[idLength-1] == '2')) {
            idLength -= 2;
        }
        string key(id, idLength);
        ReadMap::iterator found = pending.find(key);
        if (found == pending.end()) {
            // no match, remember it for later matching
            if (pending.size() >= capacity && ! capacityExceeded) {
                capacityExceeded = true;
                fprintf(stderr, "More than %d unmatched pending reads\n", capacity);
                // todo: deal with it; ignore for now, let the table grow...
            }
            pending[key] = one;
            tracker.addRead(one.getBatch());
            continue;
        }

        // found a match, remove it and return it
        *read1 = one;
        *read2 = found->second;
        // update reference counts for removed read batch, release if tracker requires it
        DataBatch release;
        if (tracker.removeRead(found->second.getBatch(), &release) && pendingReleased < release) {
            //printf("PairedReadMatcher::getNextReadPair pendingReleased %d\n", release.batchID);
            pendingReleased = release;
            checkRelease();
        }
        pending.erase(key);
        return true;
    }
}

// define static factory function

    PairedReadReader*
PairedReadReader::PairMatcher(
    int bufferSize,
    ReadReader* single)
{
    return new PairedReadMatcher(bufferSize, single);
}
