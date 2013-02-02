/*++

Module Name:

    GenomeIndex.h

Abstract:

    Headers for the index builder for the SNAP sequencer

Authors:

    Bill Bolosky, August, 2011

Environment:

    User mode service.

Revision History:

    Adapted from Matei Zaharia's Scala implementation.

--*/

#pragma once

#include "HashTable.h"
#include "Seed.h"
#include "Genome.h"
#include "ApproximateCounter.h"

class GenomeIndex {
public:

    //
    // run the indexer from command line arguments
    //
    static void runIndexer(int argc, const char **argv);

    //
    // Build a genome index and write it to a directory.  If you don't already have a saved index
    // the only way to get one is to build it into a directory and then load it from the directory.
    // NB: This deletes the Genome that's passed into it.
    //
    static bool BuildIndexToDirectory(const Genome *genome, int seedLen, double slack,
                                      bool computeBias, const char *directory, _uint64 overflowTableFactor,
                                      unsigned maxThreads);

    bool saveToDirectory(char *directoryName);
    static GenomeIndex *loadFromDirectory(char *directoryName);

    inline const Genome *getGenome() {return genome;}

    //
    // This looks up a seed and its reverse complement, and returns the number and list of hits for each.
    //
    void lookupSeed(Seed seed, unsigned *nHits, const unsigned **hits, unsigned *nRCHits, const unsigned **rcHits);

    //
    // Looks up a seed and its reverse complement, restricting the search to a given range of locations,
    // and returns the number and list of hits for each.
    //
    void lookupSeed(Seed seed, unsigned minLocation, unsigned maxLocation,
                    unsigned *nHits, const unsigned **hits, unsigned *nRCHits, const unsigned **rcHits);
    
    //
    // This issues a compiler prefetch for the genome data.
    //
    inline void prefetchGenomeData(unsigned genomeOffset) const {
        genome->prefetchData(genomeOffset);
    }

    inline int getSeedLength() const { return seedLen; }

    ~GenomeIndex();

    //
    // Allocate set of hash tables indexed by seeds with bias
    //
    static SNAPHashTable** allocateHashTables(unsigned* o_nTables, size_t capacity, double slack,
        int seedLen, double* biasTable = NULL);
    
private:
    
    static double GetHashTableSizeBias(unsigned whichTable, int seedSize);

    static double biasTable17[];
    static double biasTable18[];
    static double biasTable19[];
    static double biasTable20[];
    static double biasTable21[];
    static double biasTable22[];
    static double biasTable23[];

    static void ComputeBiasTable(const Genome* genome, int seedSize, double* table, unsigned maxThreads);

    struct ComputeBiasTableThreadContext {
        SingleWaiterObject              *doneObject;
        volatile int                    *runningThreadCount;
        unsigned                         genomeChunkStart;
        unsigned                         genomeChunkEnd;
        unsigned                         nHashTables;
        std::vector<ApproximateCounter> *approxCounters;
        const Genome                    *genome;
        volatile _int64                 *nBasesProcessed;
        unsigned                         seedLen;
        volatile _int64                 *validSeeds;

        ExclusiveLock                   *approximateCounterLocks;
    };

    static void ComputeBiasTableWorkerThreadMain(void *param);

    struct OverflowEntry;
    struct OverflowBackpointer;

    struct BuildHashTablesThreadContext {
        SingleWaiterObject              *doneObject;
        volatile int                    *runningThreadCount;
        unsigned                         genomeChunkStart;
        unsigned                         genomeChunkEnd;
        const Genome                    *genome;
        volatile _int64                 *nBasesProcessed;
        unsigned                         seedLen;
        volatile _int64                 *noBaseAvailable;
        volatile _int64                 *nonSeeds;
        volatile unsigned               *nextOverflowIndex;
        volatile _int64                 *bothComplementsUsed;
        GenomeIndex                     *index;
        unsigned                         nOverflowEntries;
        OverflowEntry                   *overflowEntries;
        OverflowBackpointer             *overflowBackpointers;
        unsigned                         nOverflowBackpointers;
        volatile unsigned               *nextOverflowBackpointer;
        volatile _int64                 *countOfDuplicateOverflows;

        ExclusiveLock                   *hashTableLocks;
        ExclusiveLock                   *overflowTableLock;
    };

    static void BuildHashTablesWorkerThreadMain(void *param);
    static void ApplyHashTableUpdate(BuildHashTablesThreadContext *context, unsigned whichHashTable, unsigned genomeLocation, unsigned lowBases, bool usingComplement,
                    _int64 *bothComplementsUsed, _int64 *countOfDuplicateOverflows);

    static int BackwardsUnsignedCompare(const void *, const void *);

    GenomeIndex();

    int seedLen;
    unsigned nHashTables;
    SNAPHashTable **hashTables;
    const Genome *genome;

    //
    // The overflow table is indexed by numbers > than the number of bases in the genome.
    // The hash table(s) point into the overflow table when they have a seed that's got more
    // than one instance in the genome.
    //
    unsigned overflowTableSize;
    size_t overflowTableVirtualAllocSize;
    unsigned *overflowTable;

    //
    // We have to build the overflow table in two stages.  While we're walking the genome, we first
    // assign tentative overflow table locations, and build up a list of places where each repeated
    // occurs.  Once we've read the whole thing (and so know the exact number of instances of each
    // repeated seed) we build the actual overflow table and go back and update the entries in the
    // hash table.  The next two structs hold the state while we're scanning the genome and are used
    // to build the final overflow table, and then are deleted.  They're only in the header file because
    // AddOverflowBackpointer needs them.
    //

    struct OverflowBackpointer {
        unsigned                 nextIndex;
        unsigned                 genomeOffset;
    };

    struct OverflowEntry {
        unsigned                *hashTableEntry;
        unsigned                 backpointerIndex;
        unsigned                 nInstances;
    };

    static void AddOverflowBackpointer(
                    OverflowEntry       *overflowEntries, 
                    OverflowBackpointer *overflowBackpointers,
                    unsigned             nOverflowBackpointers,
                    volatile unsigned   *nextOverflowBackpointer,
                    unsigned             genomeOffset);

    void fillInLookedUpResults(unsigned *subEntry, unsigned minLocation, unsigned maxLocation,
                               unsigned *nHits, const unsigned **hits);

};

