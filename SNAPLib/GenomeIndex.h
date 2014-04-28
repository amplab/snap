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

    const Genome *getGenome() {return genome;}

    //
    // This looks up a seed and its reverse complement, and returns the number and list of hits for each.
    // It guarantees that if the lookup succeeds that hits[-1] and rcHits[-1] are valid memory with 
    // arbirtary values.
    //
    virtual void lookupSeed(Seed seed, unsigned *nHits, const unsigned **hits, unsigned *nRCHits, const unsigned **rcHits) = 0;

    //
    // Looks up a seed and its reverse complement, restricting the search to a given range of locations,
    // and returns the number and list of hits for each.
    //
    virtual void lookupSeed(Seed seed, unsigned minLocation, unsigned maxLocation,
                    unsigned *nHits, const unsigned **hits, unsigned *nRCHits, const unsigned **rcHits) = 0;
  
    //
    // This issues a compiler prefetch for the genome data.
    //
    inline void prefetchGenomeData(unsigned genomeOffset) const {
        genome->prefetchData(genomeOffset);
    }

    inline int getSeedLength() const { return seedLen; }

    virtual ~GenomeIndex();

    //
    // run the indexer from command line arguments
    //
    static void runIndexer(int argc, const char **argv);
    static void runTranscriptomeIndexer(int argc, const char **argv);

    static GenomeIndex *loadFromDirectory(char *directoryName);

    static void printBiasTables();

protected:

    int seedLen;
    unsigned hashTableKeySize;
    unsigned nHashTables;
    const Genome *genome;

    //
    // The overflow table is indexed by numbers > than the number of bases in the genome.
    // The hash table(s) point into the overflow table when they have a seed that's got more
    // than one instance in the genome.
    //
    _uint64 overflowTableSize;
    unsigned *overflowTable;

    void *tablesBlob;   // All of the hash tables in one giant blob

    //
    // We have to build the overflow table in two stages.  While we're walking the genome, we first
    // assign tentative overflow table locations, and build up a list of places where each repeat
    // occurs.  Once we've read the whole thing (and so know the exact number of instances of each
    // repeated seed) we build the actual overflow table and go back and update the entries in the
    // hash table.
	//
	// The list of repeats works as a singly linked list, headed by the hash table entry.  The entries
	// use the index in the overflow table as links, rather than using real pointers, in order to save
	// space.  So that we can dynamically allocate overflow entries while still using indices to
	// find them, they're built in a two level table.
    //

    struct OverflowBackpointer {
        unsigned                 nextIndex;
        unsigned                 genomeOffset;
    };

	class OverflowBackpointerAnchor {
	public:
		OverflowBackpointerAnchor(unsigned maxOverflowEntries_);
		~OverflowBackpointerAnchor();

		OverflowBackpointer *getBackpointer(unsigned index);

	private:

        ExclusiveLock lock;

		static const unsigned batchSize;
		unsigned maxOverflowEntries;

		OverflowBackpointer **table;
	};


    //
    // Build a genome index and write it to a directory.  If you don't already have a saved index
    // the only way to get one is to build it into a directory and then load it from the directory.
    // NB: This deletes the Genome that's passed into it.
    //
    static bool BuildIndexToDirectory(const Genome *genome, int seedLen, double slack,
                                      bool computeBias, const char *directory,
                                      unsigned maxThreads, unsigned chromosomePaddingSize, bool forceExact, 
                                      unsigned hashTableKeySize, bool large, const char *histogramFileName);

 
    //
    // Allocate set of hash tables indexed by seeds with bias
    //
    static SNAPHashTable** allocateHashTables(unsigned* o_nTables, size_t countOfBases, double slack,
        int seedLen, unsigned hashTableKeySize, bool large, double* biasTable = NULL);
    
    static const unsigned GenomeIndexFormatMajorVersion = 4;
    static const unsigned GenomeIndexFormatMinorVersion = 0;
    
    static const unsigned largestBiasTable = 32;    // Can't be bigger than the biggest seed size, which is set in Seed.h.  Bigger than 32 means a new Seed structure.
    static const unsigned largestKeySize = 8;
    static double *hg19_biasTables[largestKeySize+1][largestBiasTable+1];
    static double *hg19_biasTables_large[largestKeySize+1][largestBiasTable+1];

    static void ComputeBiasTable(const Genome* genome, int seedSize, double* table, unsigned maxThreads, bool forceExact, unsigned hashTableKeySize, bool large);

    struct ComputeBiasTableThreadContext {
        SingleWaiterObject              *doneObject;
        volatile int                    *runningThreadCount;
        unsigned                         genomeChunkStart;
        unsigned                         genomeChunkEnd;
        unsigned                         nHashTables;
        unsigned                         hashTableKeySize;
        std::vector<ApproximateCounter> *approxCounters;
        const Genome                    *genome;
        volatile _int64                 *nBasesProcessed;
        unsigned                         seedLen;
        volatile _int64                 *validSeeds;
		bool							 large;

        ExclusiveLock                   *approximateCounterLocks;
    };

    static void ComputeBiasTableWorkerThreadMain(void *param);

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
        volatile _int64                 *seedsWithMultipleOccurrences;
        volatile _int64                 *bothComplementsUsed;
        GenomeIndex                     *index;
		OverflowBackpointerAnchor		*overflowAnchor;
        volatile unsigned               *nextOverflowBackpointer;
        volatile _int64                 *genomeLocationsInOverflowTable;
        unsigned                         hashTableKeySize;
		bool							 large;

        ExclusiveLock                   *hashTableLocks;
        ExclusiveLock                   *overflowTableLock;
    };

    struct PerHashTableBatch {
        PerHashTableBatch() : nUsed(0) {}

        static const unsigned nSeedsPerBatch = 1000;

        unsigned nUsed;

        struct Entry {
            bool usingComplement;
            _uint64 lowBases;
            unsigned genomeLocation;
        };

        Entry entries[nSeedsPerBatch];

        bool addSeed(unsigned genomeLocation, _uint64 seedLowBases, bool seedUsingComplement) {
            _ASSERT(nUsed < nSeedsPerBatch);
            entries[nUsed].lowBases = seedLowBases;
            entries[nUsed].usingComplement = seedUsingComplement;
            entries[nUsed].genomeLocation = genomeLocation;
            nUsed++;
            return nUsed >= nSeedsPerBatch;
        }

        void clear()
        {
            nUsed = 0;
        }
    };

    struct IndexBuildStats {
        IndexBuildStats() : noBaseAvailable(0), nonSeeds(0), bothComplementsUsed(0), genomeLocationsInOverflowTable(0),
                            unrecordedSkippedSeeds(0), seedsWithMultipleOccurrences(0) {}
        _int64 noBaseAvailable;
        _int64 nonSeeds;
        _int64 bothComplementsUsed;
        _int64 genomeLocationsInOverflowTable;
        _int64 seedsWithMultipleOccurrences;
        _uint64 unrecordedSkippedSeeds;
    };

    static const _int64 printPeriod;

    virtual void indexSeed(unsigned genomeLocation, Seed seed, PerHashTableBatch *batches, BuildHashTablesThreadContext *context, IndexBuildStats *stats, bool large);
    virtual void completeIndexing(PerHashTableBatch *batches, BuildHashTablesThreadContext *context, IndexBuildStats *stats, bool large);

    static void BuildHashTablesWorkerThreadMain(void *param);
    void BuildHashTablesWorkerThread(BuildHashTablesThreadContext *context);
    static void ApplyHashTableUpdate(BuildHashTablesThreadContext *context, _uint64 whichHashTable, unsigned genomeLocation, _uint64 lowBases, bool usingComplement,
                    _int64 *bothComplementsUsed, _int64 *genomeLocationsInOverflowTable, _int64 *seedsWithMultipleOccurrences, bool large);

    static int BackwardsUnsignedCompare(const void *, const void *);

    GenomeIndex();


    SNAPHashTable **hashTables;

    static unsigned AddOverflowBackpointer(
                    unsigned                     previousOverflowBackpointer,
                    OverflowBackpointerAnchor   *overflowAnchor,
                    volatile unsigned           *nextOverflowBackpointer,
                    unsigned                     genomeOffset);

    void fillInLookedUpResults(unsigned *subEntry, unsigned minLocation, unsigned maxLocation,
                               unsigned *nHits, const unsigned **hits);
};

class GenomeIndexLarge : public GenomeIndex {
public:
    //
    // This looks up a seed and its reverse complement, and returns the number and list of hits for each.
    // It guarantees that if the lookup succeeds that hits[-1] and rcHits[-1] are valid memory with 
    // arbirtary values.
    //
    void lookupSeed(Seed seed, unsigned *nHits, const unsigned **hits, unsigned *nRCHits, const unsigned **rcHits);

    //
    // Looks up a seed and its reverse complement, restricting the search to a given range of locations,
    // and returns the number and list of hits for each.
    //
    void lookupSeed(Seed seed, unsigned minLocation, unsigned maxLocation,
                    unsigned *nHits, const unsigned **hits, unsigned *nRCHits, const unsigned **rcHits);
    

    virtual ~GenomeIndexLarge();
};

class GenomeIndexSmall : public GenomeIndex {
public:
    //
    // This looks up a seed and its reverse complement, and returns the number and list of hits for each.
    // It guarantees that if the lookup succeeds that hits[-1] and rcHits[-1] are valid memory with 
    // arbirtary values.
    //
    void lookupSeed(Seed seed, unsigned *nHits, const unsigned **hits, unsigned *nRCHits, const unsigned **rcHits);

    //
    // Looks up a seed and its reverse complement, restricting the search to a given range of locations,
    // and returns the number and list of hits for each.
    //
    void lookupSeed(Seed seed, unsigned minLocation, unsigned maxLocation,
                    unsigned *nHits, const unsigned **hits, unsigned *nRCHits, const unsigned **rcHits);
    

    virtual ~GenomeIndexSmall();
};

