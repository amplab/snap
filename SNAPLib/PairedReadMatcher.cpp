/*++

Module Name:

    PairedReadMatcher.cpp

Abstract:

    Match paired-end reads coming in from different streams

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

// turn on to gather paired stats
//#define STATISTICS

using std::pair;


class PairedReadMatcher: public PairedReadReader
{
public:
    PairedReadMatcher(ReadReader* i_single, bool i_quicklyDropUnpairedReads);

    // PairedReadReader

    virtual ~PairedReadMatcher();

    virtual bool getNextReadPair(Read *read1, Read *read2);
    
    virtual void reinit(_int64 startingOffset, _int64 amountOfFileToProcess)
    { single->reinit(startingOffset, amountOfFileToProcess); }

    virtual void holdBatch(DataBatch batch)
    { single->holdBatch(batch); }

    virtual bool releaseBatch(DataBatch batch);

    virtual ReaderContext* getContext()
    { return single->getContext(); }

private:

    ReadWithOwnMemory* allocOverflowRead();
    void freeOverflowRead(ReadWithOwnMemory* read);
    
    ReadReader* single; // reader for single reads
    typedef _uint64 StringHash;
    typedef VariableSizeMap<StringHash,Read> ReadMap;
    DataBatch currentBatch; // for dropped reads
    bool allDroppedInCurrentBatch;
    DataBatch batch[2]; // 0 = current, 1 = previous
    ReadMap unmatched[2]; // read id -> Read
    typedef VariableSizeMap<PairedReadMatcher::StringHash,ReadWithOwnMemory*,150,MapNumericHash<PairedReadMatcher::StringHash>,80,0,true> OverflowMap;
    OverflowMap overflow; // read id -> Read
    typedef VariableSizeVector<ReadWithOwnMemory*> OverflowReadVector;
    OverflowReadVector blocks; // BigAlloc blocks
    static const int BlockSize = 10000; // # ReadWithOwnMemory per block
    ExclusiveLock blockLock; // protects adding to blocks list
    ReadWithOwnMemory* freeList; // head of free list, NULL if empty, use interlocked ops to update
    typedef VariableSizeMap<_uint64,OverflowReadVector*> OverflowReadReleaseMap;
    OverflowReadReleaseMap overflowRelease;
#ifdef VALIDATE_MATCH
    typedef VariableSizeMap<StringHash,char*> StringMap;
    StringMap strings;
    typedef VariableSizeMap<StringHash,int> HashSet;
    HashSet overflowUsed;
#endif
    int overflowTotal, overflowPeak;

    bool quicklyDropUnpairedReads;
    _uint64 nReadsQuicklyDropped;

    Read localRead;

#ifdef STATISTICS
    typedef struct
    {
        _int64 oldPairs; // # pairs matched from overflow
        _int64 oldBatches; // # distinct matches matched from overflow
        _int64 internalPairs; // #pairs matched within batch
        _int64 previousPairs; // #pairs matched with previous batch
        _int64 overflowPairs; // #pairs left over
        _int64 totalReads; // total reads in batch
        void clear() { memset(this, 0, sizeof(*this)); }
    } BatchStats;
    BatchStats currentStats, totalStats;
    VariableSizeMap<_int64,int> currentBatches;
#endif
};

PairedReadMatcher::PairedReadMatcher(
    ReadReader* i_single,
    bool i_quicklyDropUnpairedReads)
    : single(i_single),
    overflowTotal(0), overflowPeak(0),
    quicklyDropUnpairedReads(i_quicklyDropUnpairedReads),
    nReadsQuicklyDropped(0), freeList(NULL),
    currentBatch(0, 0), allDroppedInCurrentBatch(false)
{
    new (&unmatched[0]) VariableSizeMap<StringHash,Read>(10000);
    new (&unmatched[1]) VariableSizeMap<StringHash,Read>(10000);
    InitializeExclusiveLock(&blockLock);
#ifdef STATISTICS
    currentStats.clear();
    totalStats.clear();
#endif
}
    
PairedReadMatcher::~PairedReadMatcher()
{
    for (OverflowReadVector::iterator i = blocks.begin(); i != blocks.end(); i++) {
        BigDealloc(*i);
    }
    delete single;
	DestroyExclusiveLock(&blockLock);
}

    ReadWithOwnMemory*
PairedReadMatcher::allocOverflowRead()
{
    while (true) {
        ReadWithOwnMemory* next = freeList;
        if (next == NULL) {
            // alloc & init a new block of reads
            ReadWithOwnMemory* block = (ReadWithOwnMemory*) BigAlloc(BlockSize * sizeof(ReadWithOwnMemory));
            AcquireExclusiveLock(&blockLock);
            blocks.push_back(block);
            ReleaseExclusiveLock(&blockLock);
            for (int i = 0; i < BlockSize - 1; i++) {
                *(ReadWithOwnMemory**)&block[i] = &block[i+1];
            }
            while (true) {
                ReadWithOwnMemory* head = freeList;
                *(ReadWithOwnMemory**)&block[BlockSize-1] = head;
                if (InterlockedCompareExchangePointerAndReturnOldValue((void*volatile*)&freeList, block, head) == head) {
                    break;
                }
            }
        } else {
            ReadWithOwnMemory* nextHead = *(ReadWithOwnMemory**)next;
            if (InterlockedCompareExchangePointerAndReturnOldValue((void*volatile*)&freeList, nextHead, next) == next) {
                return next;
            }
        }
    }
}

    void
PairedReadMatcher::freeOverflowRead(
    ReadWithOwnMemory* read)
{
    while (true) {
        ReadWithOwnMemory* head = freeList;
        *(ReadWithOwnMemory**)read = head;
        if (InterlockedCompareExchangePointerAndReturnOldValue((void*volatile*)&freeList, read, head) == head) {
            return;
        }
    }
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
#ifdef USE_DEVTEAM_OPTIONS
            WriteErrorMessage("overflow total %d, peak %d\n", overflowTotal, overflowPeak);
#endif
            int n = unmatched[0].size() + unmatched[1].size() + overflow.size();
            if (n > 0) {
                WriteErrorMessage( " warning: PairedReadMatcher discarding %d unpaired reads at eof\n", n);
#ifdef USE_DEVTEAM_OPTIONS
                int printed = 0;
                char buffer[200];
                for (OverflowMap::iterator i = overflow.begin(); i != overflow.end() && printed < 10; i = overflow.next(i)) {
                    int l = min((unsigned) sizeof(buffer)-1, i->value->getIdLength());
                    memcpy(buffer, i->value->getId(), l);
                    buffer[l] = 0;
                    WriteErrorMessage("%s\n", buffer);
                    printed++;
                }
#endif
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
            single->releaseBatch(batch[0]);
            single->releaseBatch(batch[1]);
            return false;
        }

        if (quicklyDropUnpairedReads) {
            if (((localRead.getOriginalSAMFlags() & SAM_NEXT_UNMAPPED) == 0) && (localRead.getOriginalPNEXT() == 0 || localRead.getOriginalRNEXTLength() == 1 && localRead.getOriginalRNEXT()[0] == '*')) {
                nReadsQuicklyDropped++;
                skipped--;
                continue;
            }
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
#ifdef STATISTICS
            currentStats.oldBatches = currentBatches.size();
            currentStats.overflowPairs = unmatched[1].size();
            totalStats.internalPairs += currentStats.internalPairs;
            totalStats.previousPairs += currentStats.previousPairs;
            totalStats.oldBatches += currentStats.oldBatches;
            totalStats.oldPairs += currentStats.oldPairs;
            totalStats.overflowPairs += currentStats.overflowPairs;
            totalStats.totalReads += currentStats.totalReads;
            fprintf(stderr,"batch %d:%d: internal %d pairs, previous %d pairs, old %d pairs from %d batches, overflow %d pairs\n"
                "cumulative: internal %d pairs, previous %d pairs, old %d pairs from %d batches, overflow %d pairs\n",
                batch[0].fileID, batch[0].batchID, currentStats.internalPairs, currentStats.previousPairs, currentStats.oldPairs, currentStats.oldBatches, currentStats.overflowPairs,
                totalStats.internalPairs, totalStats.previousPairs, totalStats.oldPairs, totalStats.oldBatches, totalStats.overflowPairs);
            currentStats.clear();
            currentBatches.clear();
#endif
            // roll over batches
            if (unmatched[1].size() > 0) {
                // copy remaining reads into overflow map
                //fprintf(stderr,"warning: PairedReadMatcher overflow %d unpaired reads from %d:%d\n", unmatched[1].size(), batch[1].fileID, batch[1].batchID); //!!
                //char* buf = (char*) alloca(500);
                for (ReadMap::iterator r = unmatched[1].begin(); r != unmatched[1].end(); r = unmatched[1].next(r)) {
                    ReadWithOwnMemory* p = allocOverflowRead();
                    new (p) ReadWithOwnMemory(r->value);
                    _ASSERT(p->getData()[0]);
                    overflow.put(r->key, p);
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
                overflowTotal += unmatched[1].size();
                overflowPeak = max(overflow.size(), overflowPeak);
            }
            for (ReadMap::iterator i = unmatched[1].begin(); i != unmatched[1].end(); i = unmatched[1].next(i)) {
                i->value.dispose();
            }
            unmatched[1].exchange(unmatched[0]);
            unmatched[0].clear();
            single->releaseBatch(batch[1]);
            batch[1] = batch[0];
            batch[0] = localRead.getBatch();
            single->holdBatch(batch[0]);
#ifdef STATISTICS
        currentStats.totalReads++;
#endif
        }

        ReadMap::iterator found = unmatched[0].find(key);
        if (found != unmatched[0].end()) {
            *outputReads[1-readOneToOutputRead] = found->value;
             //fprintf(stderr, "current matched %d:%d->%d:%d %s\n", outputReads[1-readOneToOutputRead]->getBatch().fileID, outputReads[1-readOneToOutputRead]->getBatch().batchID, batch[0].fileID, batch[0].batchID, outputReads[1-readOneToOutputRead]->getId()); //!!
            unmatched[0].erase(found->key);
#ifdef STATISTICS
            currentStats.internalPairs++;
#endif
        } else {
            // try previous batch
            found = unmatched[1].find(key);
            if (found == unmatched[1].end()) {
                // try overflow
                OverflowMap::iterator found2 = overflow.find(key);
                if (found2 == overflow.end()) {
                    // no match, remember it for later matching
                    unmatched[0].put(key, localRead);
                    _ASSERT(localRead.getData()[0] && unmatched[0][key].getData()[0]);
                    //fprintf(stderr, "unmatched add %d:%d %lx\n", batch[0].fileID, batch[0].batchID, key); //!!
                    continue;
                } else {
                    // copy data into read, move from overflow table to release vector for current batch
                    found2->value->setBatch(batch[0]); // overwrite batch to match current
                    *outputReads[1-readOneToOutputRead] = * (Read*) found2->value;
                    _ASSERT(outputReads[1-readOneToOutputRead]->getData()[0]);
                    OverflowReadVector* v;
                    if (! overflowRelease.tryGet(batch[0].asKey(), &v)) {
                        v = new OverflowReadVector();
                        overflowRelease.put(batch[0].asKey(), v);
                        //fprintf(stderr,"overflow fetch into %d:%d\n", batch[0].fileID, batch[0].batchID);
                    }
                    v->push_back(found2->value);
                    overflow.erase(key);
                    //fprintf(stderr,"overflow matched %d:%d %s\n", read2->getBatch().fileID, read2->getBatch().batchID, read2->getId()); //!!
#ifdef VALIDATE_MATCH
                    overflowUsed.put(key, 1);
#endif
#ifdef STATISTICS
                    currentStats.oldPairs++;
                    currentBatches.put(read2->getBatch().asKey(), 1);
#endif
                }
            } else {
                // found a match in preceding batch
                *outputReads[1-readOneToOutputRead] = found->value;
                //fprintf(stderr,"prior matched %d:%d->%d:%d %s\n", found->value.getBatch().fileID, found->gvalue.etBatch().batchID, batch[0].fileID, batch[0].batchID, found->value.getId()); //!!
                unmatched[1].erase(found->key);
#ifdef STATISTICS
                currentStats.previousPairs++;
#endif
            }
        }

        // found a match
        *outputReads[readOneToOutputRead] = localRead;
        return true;
    }
}

    bool
PairedReadMatcher::releaseBatch(
    DataBatch batch)
{
    if (batch.asKey() == 0) {
        return true;
    } else if (single->releaseBatch(batch)) {
        OverflowReadVector* v = NULL;
        if (overflowRelease.tryGet(batch.asKey(), &v)) {
            // free memory for overflow reads
            //fprintf(stderr, "PairedReadMatcher release %d overflow reads for batch %d:%d\n", v->size(), batch.fileID, batch.batchID);
            for (OverflowReadVector::iterator i = v->begin(); i != v->end(); i++) {
                freeOverflowRead(*i);
            }
            delete v;
            overflowRelease.erase(batch.asKey());
        }
        return true;
    } else {
        return false;
    }
}

// define static factory function

    PairedReadReader*
PairedReadReader::PairMatcher(
    ReadReader* single,
    bool quicklyDropUnpairedReads)
{
    return new PairedReadMatcher(single, quicklyDropUnpairedReads);
}
