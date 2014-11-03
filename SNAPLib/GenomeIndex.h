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
#include "GenericFile_map.h"

class GenomeIndex {
public:
    const Genome *getGenome() {return genome;}

    //
    // This looks up a seed and its reverse complement, and returns the number and list of hits for each.
    // It guarantees that if the lookup succeeds that hits[-1] and rcHits[-1] are valid memory with 
    // arbirtary values.  The -32 version is used for indices with 32 bit genome offsets; using the version
    // that doesn't match the genome index GenomeLocation size is an error.  Check the index type with
    // doesGenomeIndexHave64BitLocations();
    //
    // The 64 bit version requires the called to supply a special location for storing a single forward or
    // reverse hit.  This is because in the hash table (but not overflow table), these may be stored in
    // 5-7 bytes in order to save space.  This means that there's no address in the hash table that can
    // be pointed to as a return value.  When only a single hit is returned, *hits == singleHit, so there's
    // no need to check on the caller's side.
    //
    void lookupSeed(Seed seed, _int64 *nHits, const GenomeLocation **hits, _int64 *nRCHits, const GenomeLocation **rcHits, GenomeLocation *singleHit, GenomeLocation *singleRCHit);
    void lookupSeed32(Seed seed, _int64 *nHits, const unsigned **hits, _int64 *nRCHits, const unsigned **rcHits);

    bool doesGenomeIndexHave64BitLocations() const {return locationSize > 4;}

    //
    // Looks up a seed and its reverse complement, restricting the search to a given range of locations,
    // and returns the number and list of hits for each.
    //
//    virtual void lookupSeed(Seed seed, unsigned minLocation, unsigned maxLocation,
//                    unsigned *nHits, const unsigned **hits, unsigned *nRCHits, const unsigned **rcHits) = 0;
  
    //
    // This issues a compiler prefetch for the genome data.
    //
    inline void prefetchGenomeData(GenomeLocation genomeLocation) const {
        genome->prefetchData(genomeLocation);
    }

    inline int getSeedLength() const { return seedLen; }

    virtual ~GenomeIndex();

    //
    // run the indexer from command line arguments
    //
    static void runIndexer(int argc, const char **argv);

    static GenomeIndex *loadFromDirectory(char *directoryName, bool map, bool prefetch);

    static void printBiasTables();

protected:

    int seedLen;
    unsigned hashTableKeySize;
    unsigned nHashTables;
    const Genome *genome;

    bool largeHashTable;
    unsigned locationSize;

    //
    // The overflow table is indexed by numbers > than the number of bases in the genome.
    // The hash table(s) point into the overflow table when they have a seed that's got more
    // than one instance in the genome.  For locationSize <= 4, the table is made of 32
    // bit entries (and pointed to by overflowTable32), otherwise it's 64 bit entries.
    //
    _uint64 overflowTableSize;
    unsigned *overflowTable32;
    _int64 *overflowTable64;
	GenericFile_map *mappedOverflowTable;

    void *tablesBlob;   // All of the hash tables in one giant blob
	GenericFile_map *mappedTables;

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
        _int64                   nextIndex;
        GenomeLocation           genomeLocation;
    };

	class OverflowBackpointerAnchor {
	public:
		OverflowBackpointerAnchor(_int64 maxOverflowEntries_);
		~OverflowBackpointerAnchor();

		OverflowBackpointer *getBackpointer(_int64 index);
		void trimTo(_int64 trimToIndex, FILE *trimFile);
		void loadFromFile(FILE *tripFile);

	private:

        ExclusiveLock lock;

		static const unsigned batchSize;
		_int64 maxOverflowEntries;

		OverflowBackpointer **table;

		static OverflowBackpointer spilledTableSlot;	// This value is used to indicate that the table slot in question has been spilled
	};


    //
    // Build a genome index and write it to a directory.  If you don't already have a saved index
    // the only way to get one is to build it into a directory and then load it from the directory.
    // NB: This deletes the Genome that's passed into it.
    //
    static bool BuildIndexToDirectory(const Genome *genome, int seedLen, double slack,
                                      bool computeBias, const char *directory,
                                      unsigned maxThreads, unsigned chromosomePaddingSize, bool forceExact, 
                                      unsigned hashTableKeySize, bool large, const char *histogramFileName,
                                      unsigned locationSize, bool smallMemory);

 
    //
    // Allocate set of hash tables indexed by seeds with bias
    //
    static SNAPHashTable** allocateHashTables(unsigned* o_nTables, GenomeDistance countOfBases, double slack,
        int seedLen, unsigned hashTableKeySize, bool large, unsigned locationSize, double* biasTable = NULL);
    
    static const unsigned GenomeIndexFormatMajorVersion = 5;
    static const unsigned GenomeIndexFormatMinorVersion = 0;
    
    static const unsigned largestBiasTable = 32;    // Can't be bigger than the biggest seed size, which is set in Seed.h.  Bigger than 32 means a new Seed structure.
    static const unsigned largestKeySize = 8;
    static double *hg19_biasTables[largestKeySize+1][largestBiasTable+1];
    static double *hg19_biasTables_large[largestKeySize+1][largestBiasTable+1];

    static void ComputeBiasTable(const Genome* genome, int seedSize, double* table, unsigned maxThreads, bool forceExact, unsigned hashTableKeySize, bool large);

    struct ComputeBiasTableThreadContext {
        SingleWaiterObject              *doneObject;
        volatile int                    *runningThreadCount;
        GenomeDistance                   genomeChunkStart;
        GenomeDistance                   genomeChunkEnd;
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
		unsigned						 nThreads;
		unsigned						 whichThread;
        SingleWaiterObject              *doneObject;
        volatile int                    *runningThreadCount;
        GenomeLocation                   genomeChunkStart;
        GenomeLocation                   genomeChunkEnd;
        const Genome                    *genome;
        volatile _int64                 *nBasesProcessed;
        unsigned                         seedLen;
        volatile _int64                 *noBaseAvailable;
        volatile _int64                 *nonSeeds;
        volatile _int64                 *seedsWithMultipleOccurrences;
        volatile _int64                 *bothComplementsUsed;
        GenomeIndex                     *index;
		OverflowBackpointerAnchor		*overflowAnchor;
        volatile _int64                 *nextOverflowBackpointer;
        volatile _int64                 *genomeLocationsInOverflowTable;
        unsigned                         hashTableKeySize;
		bool							 large;
        unsigned                         locationSize;

		//
		// The "small memory" option causes SNAP to write out the backpointer table as it's
		// built in order to save the memory it uses, because it's accessed sequentially as
		// the table is built.  In order to be able to tell when it's safe to free some of the
		// table, each thread occasionally records the largest backpointer index that it's done
		// with.  It's safe to free table chunks that are less than the least of these.  The lock
		// is used to keep more than one thread from trying to spill table at a time.
		//
		_int64							*lastBackpointerIndexUsedByThread;
		ExclusiveLock					*backpointerSpillLock;
		FILE							*backpointerSpillFile;

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
            GenomeLocation genomeLocation;
        };

        Entry entries[nSeedsPerBatch];

        bool addSeed(GenomeLocation genomeLocation, _uint64 seedLowBases, bool seedUsingComplement) {
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

    virtual void indexSeed(GenomeLocation genomeLocation, Seed seed, PerHashTableBatch *batches, BuildHashTablesThreadContext *context, IndexBuildStats *stats, bool large);
    virtual void completeIndexing(PerHashTableBatch *batches, BuildHashTablesThreadContext *context, IndexBuildStats *stats, bool large);

    static void BuildHashTablesWorkerThreadMain(void *param);
    void BuildHashTablesWorkerThread(BuildHashTablesThreadContext *context);
    static void ApplyHashTableUpdate(BuildHashTablesThreadContext *context, _uint64 whichHashTable, GenomeLocation genomeLocation, _uint64 lowBases, bool usingComplement,
                    _int64 *bothComplementsUsed, _int64 *genomeLocationsInOverflowTable, _int64 *seedsWithMultipleOccurrences, bool large);

    static int BackwardsUnsignedCompare(const void *, const void *);
    static int BackwardsInt64Compare(const void *, const void *);

    GenomeIndex();


    SNAPHashTable **hashTables;

    static _int64  AddOverflowBackpointer(
                        _int64                       previousOverflowBackpointer,
						BuildHashTablesThreadContext*context,
                        GenomeLocation               genomeLocation);

    void fillInLookedUpResults32(const unsigned *subEntry, _int64 *nHits, const unsigned **hits);
    void fillInLookedUpResults(GenomeLocation lookedUpLocation, _int64 *nHits, const GenomeLocation **hits, GenomeLocation *singleHitLocation);
};
