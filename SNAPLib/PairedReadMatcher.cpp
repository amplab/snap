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

#include "stdafx.h"
#include <map>
#include "Compat.h"
#include "Util.h"
#include "Read.h"
#include "DataReader.h"
#include "VariableSizeMap.h"

using std::pair;

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
    typedef _uint64 StringHash;
    typedef VariableSizeMap<StringHash,Read> ReadMap;
    DataBatch batch[2]; // 0 = current, 1 = previous
    ReadMap unmatched[2]; // read id -> Read
    typedef VariableSizeMap<StringHash,ReadWithOwnMemory> OverflowMap;
    OverflowMap overflow; // read id -> Read
    _int64 overflowMatched;
    // used only if ! autoRelease:
    bool dependents; // true if pairs from 0->1
    // manage inter-batch dependencies
    // newer depends on older (i.e. only release older after newer)
    // erase forward if older released first
    // erase both if newer released first
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
    autoRelease(i_autoRelease),
    overflowMatched(0)
{
    unmatched[0] = VariableSizeMap<_uint64,Read>(10000);
    unmatched[1] = VariableSizeMap<_uint64,Read>(10000);
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
            int n2 = (int) (overflow.size() - overflowMatched);
            if (n + n2 > 0) {
                fprintf(stderr, " warning: PairedReadMatcher discarding %d+%d unpaired reads at eof\n", n, n2);
            }
            return false;
        }
        // build key for pending read table, removing /1 or /2 at end
        const char* id = one.getId();
        unsigned idLength = one.getIdLength();
        if (idLength > 2 && id[idLength-2] == '/' && (id[idLength-1] == '1' || id[idLength-1] == '2')) {
            idLength -= 2;
        }
        StringHash key = util::hash64(id, idLength);
        if (one.getBatch() != batch[0]) {
            // roll over batches
            if (unmatched[1].size() > 0) {
                //printf("warning: PairedReadMatcher overflow %d unpaired reads from %d:%d\n", unmatched[1].size(), batch[1].fileID, batch[1].batchID); //!!
                //char* buf = (char*) alloca(500);
                for (ReadMap::iterator r = unmatched[1].begin(); r != unmatched[1].end(); r = unmatched[1].next(r)) {
                    overflow.put(r->key, ReadWithOwnMemory(r->value));
                    //memcpy(buf, r->second.getId(), r->second.getIdLength());
                    //buf[r->second.getIdLength()] = 0;
                    //printf("overflow add %d:%d %s\n", batch[1].fileID, batch[1].batchID, buf);
                }
            }

            for (ReadMap::iterator i = unmatched[1].begin(); i != unmatched[1].end(); i = unmatched[1].next(i)) {
                i->value.dispose();
            }
            unmatched[1].clear();
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
            *read2 = found->value;
            //printf("current matched %d:%d->%d:%d %s\n", read2->getBatch().fileID, read2->getBatch().batchID, batch[0].fileID, batch[0].batchID, read2->getId()); //!!
            unmatched[0].erase(found->key);
        } else {
            // try previous batch
            found = unmatched[1].find(key);
            if (found == unmatched[1].end()) {
                // try overflow
                OverflowMap::iterator found2 = overflow.find(key);
                if (found2 == overflow.end()) {
                    // no match, remember it for later matching
                    unmatched[0].put(key, one);
                    //printf("unmatched add %d:%d %s\n", batch[0].fileID, batch[0].batchID, key.data()); //!!
                    continue;
                } else {
                    // copy data into read, keep in overflow table indefinitely to preserve memory
                    *read2 = * (Read*) &found2->value;
                    _ASSERT(read2->getData()[0]);
                    overflowMatched++;
                    //printf("overflow matched %d:%d %s\n", read2->getBatch().fileID, read2->getBatch().batchID, read2->getId()); //!!
                    read2->setBatch(batch[0]); // overwrite batch so both reads have same batch, will track deps instead
                }
            } else {
                // found, remember dependency
                if ((! autoRelease) && (! dependents)) {
                    dependents = true;
                    AcquireExclusiveLock(&lock);
                    //printf("add dependency %d:%d->%d:%d\n", batch[0].fileID, batch[0].batchID, batch[1].fileID, batch[1].batchID);
                    forward.put(batch[1].asKey(), batch[0]);
                    backward.put(batch[0].asKey(), batch[1]);
                    ReleaseExclusiveLock(&lock);
                }
                *read2 = found->value;
                //printf("prior matched %d:%d->%d:%d %s\n", read2->getBatch().fileID, read2->getBatch().batchID, batch[0].fileID, batch[0].batchID, read2->getId()); //!!
                read2->setBatch(batch[0]); // overwrite batch so both reads have same batch, will track deps instead
                unmatched[1].erase(found->key);
            }
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
    // case in which i'm the newer batch
    DataBatch* b = backward.tryFind(key);
    if (b != NULL) {
        DataBatch* bf = forward.tryFind(b->asKey());
        if (bf == NULL) {
            // batch I depend on already released, can release it now
            //printf("release older batch %d:%d->%d:%d\n", batch.fileID, batch.batchID, b->fileID, b->batchID);
            single->releaseBatch(*b);
        } else {
            // forget dependency so older batch can be released later
            //printf("forget newer batch dependency %d:%d->%d:%d\n", batch.fileID, batch.batchID, b->fileID, b->batchID);
            forward.erase(b->asKey());
        }
        backward.erase(key);
    }
    // case in which I'm the older batch
    DataBatch* f = forward.tryFind(key);
    if (f != NULL) {
        // someone depends on me, signal that I've been released
        //printf("keep older batch %d:%d->%d:%d\n", f->fileID, f->batchID, batch.fileID, batch.batchID);
        forward.erase(key);
    } else {
        // noone depends on me, I can be released
        //printf("release independent batch %d:%d\n", batch.fileID, batch.batchID);
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
