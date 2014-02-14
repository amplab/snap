/*++

Module Name:

    MultiInputReadSupplier.cpp

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
#include "PairedEndAligner.h"
#include "SAM.h"
#include "Error.h"

// turn on to debug matching process
//#define VALIDATE_MATCH

using std::pair;

class PairedReadMatcher: public PairedReadReader
{
public:
    PairedReadMatcher(ReadReader* i_single, bool i_autoRelease, bool i_quicklyDropUnpairedReads);

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
    DataBatch currentBatch; // for dropped reads
    bool allDroppedInCurrentBatch;
    DataBatch batch[2]; // 0 = current, 1 = previous
    bool releasedBatch[2];  // whether each batch has been released
    ReadMap unmatched[2]; // read id -> Read
    typedef VariableSizeMap<StringHash,ReadWithOwnMemory> OverflowMap;
    OverflowMap overflow; // read id -> Read
#ifdef VALIDATE_MATCH
    typedef VariableSizeMap<StringHash,char*> StringMap;
    StringMap strings;
    typedef VariableSizeMap<StringHash,int> HashSet;
    HashSet overflowUsed;
#endif
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

    bool quicklyDropUnpairedReads;
    _uint64 nReadsQuicklyDropped;

    Read localRead;
};

PairedReadMatcher::PairedReadMatcher(
    ReadReader* i_single,
    bool i_autoRelease,
    bool i_quicklyDropUnpairedReads)
    : single(i_single),
    forward(),
    backward(),
    dependents(false),
    autoRelease(i_autoRelease),
    overflowMatched(0),
    quicklyDropUnpairedReads(i_quicklyDropUnpairedReads),
    nReadsQuicklyDropped(0),
    currentBatch(0, 0), allDroppedInCurrentBatch(false)
{
    unmatched[0] = VariableSizeMap<_uint64,Read>(10000);
    unmatched[1] = VariableSizeMap<_uint64,Read>(10000);
    releasedBatch[0] = releasedBatch[1] = false;
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
    Read *outputReads[NUM_READS_PER_PAIR];
    outputReads[0] = read1;
    outputReads[1] = read2;
    int readOneToOutputRead;    // This is used to determine which of the output reads corresponds to one (the read that just came from getNextRead())
                                // That, in turn, is determined by the S/BAM flags in the read saying whether it was first-in-template.

    int skipped = 0;
    while (true) {
        if (skipped++ == 10000) {
            WriteErrorMessage( "warning: no matching read pairs in 10,000 reads, input file might be unsorted or have unexpected read id format\n");
        }

        if (! single->getNextRead(&localRead)) {
            if (allDroppedInCurrentBatch) {
                single->releaseBatch(currentBatch);
            }
            int n = unmatched[0].size() + unmatched[1].size();
            int n2 = (int) (overflow.size() - overflowMatched);
            if (n + n2 > 0) {
                WriteErrorMessage( " warning: PairedReadMatcher discarding %d+%d unpaired reads at eof\n", n, n2);
#ifdef VALIDATE_MATCH
                for (int i = 0; i < 2; i++) {
                    fprintf(stdout, "unmatched[%d]\n", i);
                    for (ReadMap::iterator j = unmatched[i].begin(); j != unmatched[i].end(); j = unmatched[i].next(j)) {
                        fprintf(stdout, "%s\n", strings[j->key]);
                    }
                }
                int printed = 0;
                fprintf(stdout, "sample of overflow\n");
                for (OverflowMap::iterator o = overflow.begin(); printed < 500 && o != overflow.end(); o = overflow.next(o)) {
                    if (NULL == overflowUsed.tryFind(o->key)) {
                        printed++;
                        fprintf(stdout, "%s\n", strings[o->key]);
                    }
                }
#endif

            }
            if (nReadsQuicklyDropped > 0) {
                WriteErrorMessage(" warning: PairedReadMatcher dropped %lld reads because they didn't have RNEXT and PNEXT filled in.\n"
                               " If your input file was generated by a single-end alignment (or this seems too big), use the -ku flag\n",
                    nReadsQuicklyDropped);
            }
            return false;
        }

        if (quicklyDropUnpairedReads) {
            if (localRead.getBatch() != currentBatch) {
                if (allDroppedInCurrentBatch) {
                    single->releaseBatch(currentBatch);
                }
                currentBatch = localRead.getBatch();
                allDroppedInCurrentBatch = true;
            }
            if (localRead.getOriginalPNEXT() == 0 || localRead.getOriginalRNEXTLength() == 1 && localRead.getOriginalRNEXT()[0] == '*') {
                nReadsQuicklyDropped++;
                skipped--;
                continue;
            }
            allDroppedInCurrentBatch = false;
        }

        if (localRead.getOriginalSAMFlags() & SAM_FIRST_SEGMENT) {
            readOneToOutputRead = 0;
        } else {
            readOneToOutputRead = 1;
        }

        // build key for pending read table, removing /1 or /2 at end
        const char* id = localRead.getId();
        unsigned idLength = localRead.getIdLength();
        // truncate at space or slash
        char* slash = (char*) memchr((void*)id, '/', idLength);
        if (slash != NULL) {
            idLength = (unsigned)(slash - id);
        }
        char* space = (char*) memchr((void*)id, ' ', idLength);
        if (space != NULL) {
            idLength = (unsigned)(space - id);
        }
        StringHash key = util::hash64(id, idLength);
#ifdef VALIDATE_MATCH
        char* s = new char[idLength+1];
        memcpy(s, id, idLength);
        s[idLength] = 0;
        char** p = strings.tryFind(key);
        if (p != NULL && strcmp(*p, s)) {
          WriteErrorMessage( "hash collision %ld of %s and %s\n", key, *p, s);
          soft_exit(1);
        }
        if (p == NULL) {
          strings.put(key, s);
        }
#endif
        if (localRead.getBatch() != batch[0]) {
            // roll over batches
            if (unmatched[1].size() > 0) {
                //fprintf(stderr, "warning: PairedReadMatcher overflow %d unpaired reads from %d:%d\n", unmatched[1].size(), batch[1].fileID, batch[1].batchID); //!!
                //char* buf = (char*) alloca(500);
                for (ReadMap::iterator r = unmatched[1].begin(); r != unmatched[1].end(); r = unmatched[1].next(r)) {
                    overflow.put(r->key, ReadWithOwnMemory(r->value));
#ifdef VALIDATE_MATCH
                    char*s2 = *strings.tryFind(r->key);
                    int len = strlen(s2);
                    _ASSERT(! strncmp(s2, r->value.getId(), len));
                    ReadWithOwnMemory* rd = overflow.tryFind(r->key);
                    _ASSERT(! strncmp(s2, rd->getId(), len));
#endif
                    //memcpy(buf, r->value.getId(), r->value.getIdLength());
                    //buf[r->value.getIdLength()] = 0;
                    //fprintf(stderr, "overflow add %d:%d %s\n", batch[1].fileID, batch[1].batchID, buf);
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
            DataBatch overflowBatch = batch[1];
            batch[1] = batch[0];
            bool releaseOverflowBatch = releasedBatch[1];
            releasedBatch[1] = releasedBatch[0];
            batch[0] = localRead.getBatch();
            releasedBatch[0] = false;
            dependents = false;
            if (releaseOverflowBatch && ! autoRelease) {
                //fprintf(stderr, "release deferred batch %d:%d\n", overflowBatch.fileID, overflowBatch.batchID);
                releaseBatch(overflowBatch);
            }
        }

        if (quicklyDropUnpairedReads && (localRead.getOriginalPNEXT() == 0 || localRead.getOriginalRNEXTLength() == 1 && localRead.getOriginalRNEXT()[0] == '*')) {
            nReadsQuicklyDropped++;
            skipped--;
            continue;
        }

        ReadMap::iterator found = unmatched[0].find(key);
        if (found != unmatched[0].end()) {
            *outputReads[1-readOneToOutputRead] = found->value;
             //fprintf(stderr, "current matched %d:%d->%d:%d %s\n", outputReads[1-readOneToOutputRead]->getBatch().fileID, outputReads[1-readOneToOutputRead]->getBatch().batchID, batch[0].fileID, batch[0].batchID, outputReads[1-readOneToOutputRead]->getId()); //!!
            unmatched[0].erase(found->key);
        } else {
            // try previous batch
            found = unmatched[1].find(key);
            if (found == unmatched[1].end()) {
                // try overflow
                OverflowMap::iterator found2 = overflow.find(key);
                if (found2 == overflow.end()) {
                    // no match, remember it for later matching
                    unmatched[0].put(key, localRead);
                    //fprintf(stderr, "unmatched add %d:%d %lx\n", batch[0].fileID, batch[0].batchID, key); //!!
                    continue;
                } else {
                    // copy data into read, keep in overflow table indefinitely to preserve memory
                    *outputReads[1-readOneToOutputRead] = * (Read*) &found2->value;
                    _ASSERT(outputReads[1-readOneToOutputRead]->getData()[0]);
                    overflowMatched++;
#ifdef VALIDATE_MATCH
                    overflowUsed.put(key, 1);
#endif
                    //fprintf(stderr, "overflow matched %d:%d %s\n", outputReads[1-readOneToOutputRead]->getBatch().fileID, outputReads[1-readOneToOutputRead]->getBatch().batchID, outputReads[1-readOneToOutputRead]->getId()); //!!
                    outputReads[1-readOneToOutputRead]->setBatch(batch[0]); // overwrite batch so both reads have same batch, will track deps instead
                }
            } else {
                // found, remember dependency
                if ((! autoRelease) && (! dependents)) {
                    dependents = true;
                    AcquireExclusiveLock(&lock);
                    //fprintf(stderr, "add dependency %d:%d->%d:%d\n", batch[0].fileID, batch[0].batchID, batch[1].fileID, batch[1].batchID);
                    forward.put(batch[1].asKey(), batch[0]);
                    backward.put(batch[0].asKey(), batch[1]);
                    ReleaseExclusiveLock(&lock);
                }

                *outputReads[1-readOneToOutputRead] = found->value;
                //fprintf(stderr, "prior matched %d:%d->%d:%d %s\n", outputReads[1-readOneToOutputRead]->getBatch().fileID, outputReads[1-readOneToOutputRead]->getBatch().batchID, batch[0].fileID, batch[0].batchID, outputReads[1-readOneToOutputRead]->getId()); //!!
                outputReads[1-readOneToOutputRead]->setBatch(batch[0]); // overwrite batch so both reads have same batch, will track deps instead
                unmatched[1].erase(found->key);
            }
        }

        // found a match
        *outputReads[readOneToOutputRead] = localRead;
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
    for (int i = 0; i < 2; i++) {
      if (batch == this->batch[i]) {
        if (! releasedBatch[i]) {
          releasedBatch[i] = true;
          //fprintf(stderr, "releaseBatch %d:%d active %d, deferred\n", batch.fileID, batch.batchID, i);
        }
        return;
      }
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
            //fprintf(stderr, "release older batch %d:%d->%d:%d\n", batch.fileID, batch.batchID, b->fileID, b->batchID);
            single->releaseBatch(*b);
        } else {
            // forget dependency so older batch can be released later
            //fprintf(stderr, "forget newer batch dependency %d:%d->%d:%d\n", batch.fileID, batch.batchID, b->fileID, b->batchID);
            forward.erase(b->asKey());
        }
        backward.erase(key);
    }
    // case in which I'm the older batch
    DataBatch* f = forward.tryFind(key);
    if (f != NULL) {
        // someone depends on me, signal that I've been released
        //fprintf(stderr, "keep older batch %d:%d->%d:%d\n", f->fileID, f->batchID, batch.fileID, batch.batchID);
        forward.erase(key);
    } else {
        // noone depends on me, I can be released
        //fprintf(stderr, "release independent batch %d:%d\n", batch.fileID, batch.batchID);
        single->releaseBatch(batch);
    }
    ReleaseExclusiveLock(&lock);
}

// define static factory function

    PairedReadReader*
PairedReadReader::PairMatcher(
    ReadReader* single,
    bool autoRelease,
    bool quicklyDropUnpairedReads)
{
    return new PairedReadMatcher(single, autoRelease, quicklyDropUnpairedReads);
}
