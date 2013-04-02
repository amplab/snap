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
#include "VariableSizeMap.h"

using std::pair;
using std::map;
using std::string;

class PairedReadMatcher: public PairedReadReader
{
public:
    PairedReadMatcher(ReadReader* i_single, bool i_autoRelease);

    // PairedReadReader

    virtual ~PairedReadMatcher();

    virtual bool getNextReadPair(Read *read1, Read *read2);
    
    virtual void reinit(_int64 startingOffset, _int64 amountOfFileToProcess)
    { single->reinit(startingOffset, amountOfFileToProcess); }

    void releaseBatch(DataBatch batch);

private:
    
    const bool autoRelease;
    ReadReader* single; // reader for single reads
    typedef map<string,Read> ReadMap;
    DataBatch batch[2]; // 0 = current, 1 = previous
    ReadMap unmatched[2]; // read id -> Read
    // used only if ! autoRelease:
    bool dependents; // true if pairs from 0->1
    ExclusiveLock lock; // exclusive access to forward/backward
    typedef VariableSizeMap<DataBatch::Key,DataBatch> BatchMap;
    BatchMap forward; // dependencies from older batch (was unmatched) -> newer batch
    BatchMap backward; // newer batch -> older batch
};

PairedReadMatcher::PairedReadMatcher(
    ReadReader* i_single,
    bool i_autoRelease)
    : single(i_single),
    forward(),
    backward(),
    dependents(false),
    autoRelease(i_autoRelease)
{
    if (! autoRelease) {
        InitializeExclusiveLock(&lock);
    }
}
    
PairedReadMatcher::~PairedReadMatcher()
{
    if (! autoRelease) {
        DestroyExclusiveLock(&lock);
    }
    delete single;
}

    bool
PairedReadMatcher::getNextReadPair(
    Read *read1,
    Read *read2)
{
    while (true) {
        Read one;
        if (! single->getNextRead(&one)) {
            int n = unmatched[0].size() + unmatched[1].size();
            if (n > 0) {
                fprintf(stderr, " warning: PairedReadMatcher%d discarding unpaired reads at eof\n", n);
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
        if (one.getBatch() != batch[0]) {
            // roll over batches
            if (unmatched[1].size() > 0) {
                fprintf(stderr, "warning: PairedReadMatcher discarding %d unpaired reads\n", unmatched[1].size());
            }
            unmatched[1] = unmatched[0];
            unmatched[0].clear();
            if (autoRelease) {
                single->releaseBatch(batch[1]);
            }
            batch[1] = batch[0];
            batch[0] = one.getBatch();
            dependents = false;
        }
        ReadMap::iterator found = unmatched[0].find(key);
        if (found != unmatched[0].end()) {
            *read2 = found->second;
            unmatched[0].erase(found);
        } else {
            // try previous batch
            found = unmatched[1].find(key);
            if (found == unmatched[1].end()) {
                // no match, remember it for later matching
	            unmatched[0][key] = one;
                continue;
            }
            // found, remember dependency
            if (autoRelease && ! dependents) {
                dependents = true;
                AcquireExclusiveLock(&lock);
                forward.put(batch[1].asKey(), batch[0]);
                backward.put(batch[0].asKey(), batch[1]);
                ReleaseExclusiveLock(&lock);
            }
            *read2 = found->second;
            read2->setBatch(batch[0]); // overwrite batch so both reads have same batch, will track deps instead
            unmatched[1].erase(found);
        }

        // found a match
        *read1 = one;
        return true;
    }
}

    void
PairedReadMatcher::releaseBatch(
    DataBatch batch)
{
    if (autoRelease) {
        return;
    }
    // only release when both forward & backward dependent batches have been released
    AcquireExclusiveLock(&lock);
    DataBatch::Key key = batch.asKey();
    DataBatch* f = forward.tryFind(key);
    bool keep = false;
    if (f != NULL) {
        DataBatch* fb = backward.tryFind(f->asKey());
        keep = fb != NULL;
        if (fb == NULL) {
            single->releaseBatch(*f);
        }
        forward.erase(key);
    }
    DataBatch* b = backward.tryFind(key);
    if (b != NULL) {
        DataBatch* bf = forward.tryFind(b->asKey());
        keep |= bf != NULL;
        if (bf == NULL) {
            single->releaseBatch(*b);
        }
        backward.erase(key);
    }
    if (! keep) {
        single->releaseBatch(batch);
    }
    ReleaseExclusiveLock(&lock);
}

// define static factory function

    PairedReadReader*
PairedReadReader::PairMatcher(
    ReadReader* single,
    bool autoRelease)
{
    return new PairedReadMatcher(single, autoRelease);
}
