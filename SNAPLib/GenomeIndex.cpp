/*++

Module Name:

    GenomeIndex.cpp

Abstract:

    Index (hash table) builder for the SNAP sequencer

Authors:

    Bill Bolosky, August, 2011

Environment:

    User mode service.

Revision History:

    Adapted from Matei Zaharia's Scala implementation.

--*/

#include "stdafx.h"
#include "ApproximateCounter.h"
#include "BigAlloc.h"
#include "Compat.h"
#include "FASTA.h"
#include "FixedSizeSet.h"
#include "FixedSizeVector.h"
#include "Genome.h"
#include "GenomeIndex.h"
#include "HashTable.h"
#include "Seed.h"

using namespace std;

static const int DEFAULT_SEED_SIZE = 20;
static const double DEFAULT_SLACK = 0.3;


static void usage()
{
    fprintf(stderr,
            "Usage: snap index <input.fa> <output-dir> [<options>]\n"
            "Options:\n"
            "  -s           Seed size (default: %d)\n"
            "  -h           Hash table slack (default: %.1f)\n"
            "  -hg19        Use pre-computed table bias for hg19, which results in better speed, balance, and memory footprint but may not work for other references.\n"
            "  -Ofactor     Specify the size of the overflow space.  This will change the memory footprint, but may be needed for some genomes.\n"
            "               Larger numbers use more memory but work better with more repetitive genomes.  Smaller numbers reduce the memory\n"
            "               footprint, but may cause the index build to fail.  Making -O larger than necessary will not affect the resuting\n"
            "               index.  Factor must be between 1 and 1000, and the default is 40.\n"
            " -tMaxThreads  Specify the maximum number of threads to use. Default is the number of cores\n",
            DEFAULT_SEED_SIZE,
            DEFAULT_SLACK);
    exit(1);
}


    void
GenomeIndex::runIndexer(
    int argc,
    const char **argv)
{
    if (argc < 2) {
        usage();
    }

    const char *fastaFile = argv[0];
    const char *outputDir = argv[1];

    unsigned maxThreads = GetNumberOfProcessors();

    int seedLen = DEFAULT_SEED_SIZE;
    double slack = DEFAULT_SLACK;
    bool computeBias = true;
    _uint64 overflowTableFactor = 40;

    for (int n = 2; n < argc; n++) {
        if (strcmp(argv[n], "-s") == 0) {
            if (n + 1 < argc) {
                seedLen = atoi(argv[n+1]);
                n++;
            } else {
                usage();
            }
        } else if (strcmp(argv[n], "-h") == 0) {
            if (n + 1 < argc) {
                slack = atof(argv[n+1]);
                n++;
            } else {
                usage();
            }
        } else if (strcmp(argv[n], "-c") == 0) {
            computeBias = true;
        } else if (strcmp(argv[n], "-hg19") == 0) {
            computeBias = false;        
        } else if (argv[n][0] == '-' && argv[n][1] == 'O') {
            overflowTableFactor = atoi(argv[n]+2);
            if (overflowTableFactor < 1 || overflowTableFactor > 1000) {
                fprintf(stderr,"Overflow table factor must be between 1 and 1000 inclusive (and you need to not leave a space after '-O')\n");
                exit(1);
            }
        } else if (argv[n][0] == '-' && argv[n][1] == 't') {
            maxThreads = atoi(argv[n]+2);
            if (maxThreads < 1 || maxThreads > 100) {
                fprintf(stderr,"maxThreads must be between 1 and 100 inclusive (and you need to not leave a space after '-t')\n");
                exit(1);
            }
        } else {
            fprintf(stderr, "Invalid argument: %s\n\n", argv[n]);
            usage();
        }
    }

    if (seedLen < 16 || seedLen > 23) {
        // Right now there's some hard-coded stuff, like the seed warp table, that
        // does not work for seed lengths outside the range of 16-23.
        fprintf(stderr, "Seed length must be between 16 and 23\n");
        exit(1);
    }

    printf("Hash table slack %lf\nLoading FASTA file '%s' into memory...", slack, fastaFile);
    _int64 start = timeInMillis();
    const Genome *genome = ReadFASTAGenome(fastaFile);
    if (NULL == genome) {
        fprintf(stderr, "Unable to read FASTA file\n");
        exit(1);
    }
    printf("%llds\n", (timeInMillis() + 500 - start) / 1000);
    unsigned nBases = genome->getCountOfBases();
    if (!GenomeIndex::BuildIndexToDirectory(genome, seedLen, slack, computeBias, outputDir, overflowTableFactor, maxThreads)) {
        fprintf(stderr, "Genome index build failed\n");
        exit(1);
    }
    _int64 end = timeInMillis();
    printf("Index build and save took %llds (%lld bases/s)\n",
           (end - start) / 1000, (_int64) nBases / max((end - start) / 1000, (_int64) 1)); 
}

SNAPHashTable** GenomeIndex::allocateHashTables(
    unsigned* o_nTables,
    size_t capacity,
    double slack,
    int seedLen,
    double* biasTable)
{
    if (slack <= 0) {
        fprintf(stderr, "allocateHashTables: must have positive slack for the hash table to work.  0.3 is probably OK, 0.1 is minimal, less will wreak havoc with perf.\n");
        return false;
    }

    if (seedLen <= 0) {
        fprintf(stderr, "allocateHashTables: seedLen is too small (must be > 0, and practically should be >= 15 or so.\n");
        return false;
    }
    
    //
    // Make an array of HashTables, size depending on the seed size.  The way the index works is that we use
    // the low 32 bits (16 bases) of the seed as a hash key.  Any remaining bases are used as an index into the
    // particular hash table in question.
    //
    unsigned nHashTables = 1 << ((max(seedLen,16) - 16) * 2);

    //
    // Average size of the hash table.  We bias this later based on the actual content of the genome.
    //
    size_t hashTableSize = (size_t) ((double)capacity * (slack + 1.0) / nHashTables);
    
    SNAPHashTable **hashTables = new SNAPHashTable*[nHashTables];
    
    for (unsigned i = 0; i < nHashTables; i++) {
        //
        // Create the actual hash tables.  It turns out that the human genome is highly non-uniform in its
        // sequences of bases, so we bias the hash table sizes based on their popularity (which is emperically
        // measured.)
        //
        double bias = biasTable != NULL ? biasTable[i] : GenomeIndex::GetHashTableSizeBias(i, seedLen);
        unsigned biasedSize = (unsigned) (hashTableSize * bias);
        if (biasedSize < 100) {
            biasedSize = 100;
        }
        hashTables[i] = new SNAPHashTable(biasedSize, false);

        if (NULL == hashTables[i]) {
            fprintf(stderr,"IndexBuilder: unable to allocate HashTable %d of %d\n",i+1,nHashTables);
            exit(1);
        }
    }

    *o_nTables = nHashTables;
    return hashTables;
}

    bool
GenomeIndex::BuildIndexToDirectory(const Genome *genome, int seedLen, double slack, bool computeBias, const char *directoryName, _uint64 overflowTableFactor,
                                    unsigned maxThreads)
{
    if (mkdir(directoryName, 0777) != 0 && errno != EEXIST) {
        fprintf(stderr,"BuildIndex: failed to create directory %s\n",directoryName);
        return false;
    }

    GenomeIndex *index = new GenomeIndex();
    index->genome = NULL;   // We always delete the index when we're done, but we delete the genome first to save space during the overflow table build.

    unsigned countOfBases = genome->getCountOfBases();
    if (countOfBases > 0xfffffff0) {
        fprintf(stderr, "Genome is too big for this index.  Must be some headroom beneath 2^32 bases.\n");
        return false;
    }

    // Compute bias table sizes, unless we're using the precomputed ones hardcoded in BiasTables.cpp
    double *biasTable = NULL;
    if (computeBias) {
        unsigned nHashTables = 1 << ((max(seedLen,16) - 16) * 2);
        biasTable = new double[nHashTables];
        ComputeBiasTable(genome, seedLen, biasTable, maxThreads);
    }

    printf("Allocating memory for hash and overflow tables...");
    _int64 start = timeInMillis();
    unsigned nHashTables;
    SNAPHashTable** hashTables = index->hashTables =
        allocateHashTables(&nHashTables, countOfBases, slack, seedLen, biasTable);
    index->nHashTables = nHashTables;

    //
    // Set up the hash tables.  Each table has a key value of the lower 32 bits of the seed, and data
    // of two integers.  There is one integer each for the seed and its reverse complement (i.e., what you'd
    // get from the complementary DNA strand, A<->T and G<->C with the string order reversed).  
    // The first integer is always for the version of the seed with "lower" value, using an arbitrary
    // total order that we define in Seed.h.  Some seeds are their own reverse complements (e.g.,
    // AGCT), in which case only the first integer is used.
    //

    //
    // We want to use all 32 bits of the key space for the hash table (i.e., all possible combinations of 16 bases).
    // However, the hash table needs to have some key + data value that we promise never to use so that it'll know
    // what is a free entry.  We provide that here: -1 key and -1 in each integer.  (Well, they're unsigned, but you
    // get the idea).
    //
    unsigned unusedKeyValue = 0xffffffff;
    unsigned unusedDataValue[2];
    unusedDataValue[0] = 0xffffffff;
    unusedDataValue[1] = 0xffffffff;

    // Don't create *too* many overflow entries because they consume significant memory.
    // It would be good to grow these dynamically instead of living with a fixed-size array.
    unsigned nOverflowEntries = __min((unsigned)((_uint64)countOfBases * overflowTableFactor / (_uint64) 1000), 0xfffffffd - countOfBases);
    unsigned nOverflowBackpointers = (unsigned)__min((_uint64)countOfBases, (_uint64)nOverflowEntries * (_uint64)6);

  
    OverflowEntry *overflowEntries = new OverflowEntry[nOverflowEntries];
    if (NULL == overflowEntries) {
        fprintf(stderr,"Unable to allocate oveflow entries.\n");
        exit(1);
    }

   OverflowBackpointer *overflowBackpointers = new OverflowBackpointer[nOverflowBackpointers];
    if (NULL == overflowBackpointers) {
        fprintf(stderr,"Unable to allocate overflow backpointers\n");
        exit(1);
    }
  printf("%llds\n%d nOverflowEntries, %lld bytes, %d nOverflowBackpointers, %lld bytes\nBuilding hash tables.\n", 
        (timeInMillis() + 500 - start) / 1000,
        nOverflowEntries, (_int64)nOverflowEntries * sizeof(OverflowEntry), nOverflowBackpointers, (_int64) nOverflowBackpointers * sizeof(OverflowBackpointer));
  
    start = timeInMillis();
    volatile unsigned nextOverflowBackpointer = 0;

    volatile _int64 nonSeeds = 0;
    volatile unsigned nextOverflowIndex = 0;
    volatile _int64 countOfDuplicateOverflows = 0;     // Number of extra hits on duplicate indices.  This should come out once we implement the overflow table.
    volatile _int64 bothComplementsUsed = 0;    // Number of hash buckets where both complements are used
    volatile _int64 noBaseAvailable = 0;        // Number of places where getSubstring returned null.
    volatile _int64 nBasesProcessed = 0;
    volatile int runningThreadCount;

    SingleWaiterObject doneObject;
    CreateSingleWaiterObject(&doneObject);

    unsigned nThreads = __min(GetNumberOfProcessors(), maxThreads);
    BuildHashTablesThreadContext *threadContexts = new BuildHashTablesThreadContext[nThreads];

    ExclusiveLock *hashTableLocks = new ExclusiveLock[nHashTables];
    for (unsigned i = 0; i < nHashTables; i++) {
        InitializeExclusiveLock(&hashTableLocks[i]);
    }

    runningThreadCount = nThreads;

    unsigned nextChunkToProcess = 0;
    for (unsigned i = 0; i < nThreads; i++) {
        threadContexts[i].doneObject = &doneObject;
        threadContexts[i].genome = genome;
        threadContexts[i].genomeChunkStart = nextChunkToProcess;
        if (i == nThreads - 1) {
            nextChunkToProcess = countOfBases - seedLen - 1;
        } else {
            nextChunkToProcess += (countOfBases - seedLen) / nThreads;
        }
        threadContexts[i].genomeChunkEnd = nextChunkToProcess;
        threadContexts[i].nBasesProcessed = &nBasesProcessed;
        threadContexts[i].index = index;
        threadContexts[i].runningThreadCount = &runningThreadCount;
        threadContexts[i].seedLen = seedLen;
        threadContexts[i].noBaseAvailable = &noBaseAvailable;
        threadContexts[i].nonSeeds = &nonSeeds;
        threadContexts[i].nextOverflowIndex = &nextOverflowIndex;
        threadContexts[i].countOfDuplicateOverflows = &countOfDuplicateOverflows;
        threadContexts[i].bothComplementsUsed = &bothComplementsUsed;
        threadContexts[i].nOverflowEntries = nOverflowEntries;
        threadContexts[i].overflowEntries = overflowEntries;
        threadContexts[i].overflowBackpointers = overflowBackpointers;
        threadContexts[i].nOverflowBackpointers = nOverflowBackpointers;
        threadContexts[i].nextOverflowBackpointer = &nextOverflowBackpointer;
        threadContexts[i].countOfDuplicateOverflows = &countOfDuplicateOverflows;
        threadContexts[i].hashTableLocks = hashTableLocks;

        StartNewThread(BuildHashTablesWorkerThreadMain, &threadContexts[i]);
    }

    WaitForSingleWaiterObject(&doneObject);
    DestroySingleWaiterObject(&doneObject);

    for (unsigned j = 0; j < index->nHashTables; j++) {
        printf("HashTable[%d] has %lld used elements, loading %lld%%\n",j,(_int64)hashTables[j]->GetUsedElementCount(),
                (_int64)hashTables[j]->GetUsedElementCount() * 100 / (_int64)hashTables[j]->GetTableSize());
    }

    printf("%d(%lld%%) overflow entries, %d overflow backpointers, %d(%lld%%) duplicate overflows, %d(%lld%%) bad seeds, %d both complements used %d no string\n",
        nextOverflowIndex,
        ((_int64)nextOverflowIndex)*100 / countOfBases,
        nextOverflowBackpointer,
        countOfDuplicateOverflows,
        (_int64)countOfDuplicateOverflows * 100 / countOfBases,
        (int) nonSeeds,
        (_int64)nonSeeds *100 / countOfBases,
        bothComplementsUsed,
        noBaseAvailable);

    printf("Hash table build took %llds\nSaving genome...",(timeInMillis() + 500 - start) / 1000);
    start = timeInMillis();

    //
    // Save the genome to the directory and delete it to free up a little memory for the overflow table build.  We need all we can get.
    //
    
    const unsigned filenameBufferSize = MAX_PATH+1;
    char filenameBuffer[filenameBufferSize];
    
    snprintf(filenameBuffer,filenameBufferSize,"%s%cGenome",directoryName,PATH_SEP);
    if (!genome->saveToFile(filenameBuffer)) {
        fprintf(stderr,"GenomeIndex::saveToDirectory: Failed to save the genome itself\n");
        return false;
    }

    delete genome;
    genome = NULL;

    printf("%llds\nBuilding overflow table.\n", (timeInMillis() + 500 - start) / 1000);
    start = timeInMillis();
    fflush(stdout);

    //
    // Now build the real overflow table and simultaneously fixup the hash table entries.
    // Its format is one unsigned of the number of genome locations matching the
    // particular seed, followed by that many genome offsets.  The count of entries is
    // 3 * the entries in our overflow builder table, plus the number of surplus
    // overflows.  The 3 is one for the count, and one for each of the entries that went
    // onto the backpointer list when the overflow entry was first created.
    //
    index->overflowTableSize = nextOverflowIndex * 3 + (unsigned)countOfDuplicateOverflows;
    index->overflowTable = (unsigned *)BigAlloc(index->overflowTableSize * sizeof(*index->overflowTable),&index->overflowTableVirtualAllocSize);

    unsigned nBackpointersProcessed = 0;
    _int64 lastPrintTime = timeInMillis();

    unsigned overflowTableIndex = 0;
    for (unsigned i = 0 ; i < nextOverflowIndex; i++) {
        OverflowEntry *overflowEntry = &overflowEntries[i];
        _ASSERT(overflowEntry->nInstances >= 2);

        if (timeInMillis() - lastPrintTime > 60 * 1000) {
            printf("%d/%d duplicate seeds, %d/%d backpointers processed\n",i,nextOverflowIndex,nBackpointersProcessed,nextOverflowBackpointer-1);
            lastPrintTime = timeInMillis();
        }

        //
        // Start by patching up the hash table.
        //
        *overflowEntry->hashTableEntry = overflowTableIndex + countOfBases;

        //
        // Now fill in the overflow table. First, the count of instances of this seed in the genome.
        //
        _ASSERT(overflowTableIndex + overflowEntry->nInstances + 1 <= index->overflowTableSize);
        index->overflowTable[overflowTableIndex] = overflowEntry->nInstances;
        overflowTableIndex++;

        //
        // Followed by the actual addresses in the genome.
        //
        nBackpointersProcessed += overflowEntry->nInstances;
        for (unsigned j = 0; j < overflowEntry->nInstances; j++) {
            _ASSERT(-1 != overflowEntry->backpointerIndex);
            OverflowBackpointer *backpointer = &overflowBackpointers[overflowEntry->backpointerIndex];
            index->overflowTable[overflowTableIndex] = backpointer->genomeOffset;
            overflowTableIndex++;
            
            overflowEntry->backpointerIndex = backpointer->nextIndex;
        }

        //
        // Now sort them, because the multi thread insertion results in random order, but SNAP expects them to be in descending order.
        //
        qsort(&index->overflowTable[overflowTableIndex - overflowEntry->nInstances], overflowEntry->nInstances, sizeof(index->overflowTable[0]), BackwardsUnsignedCompare);
    }
    _ASSERT(overflowTableIndex == index->overflowTableSize);    // We used exactly what we expected to use.

    delete [] overflowEntries;
    overflowEntries = NULL;

    delete [] overflowBackpointers;
    overflowBackpointers = NULL;

    //
    // Now save out the part of the index that's independent of the genome itself.
    //
    printf("Overflow table build took %llds\nSaving genome index...", (timeInMillis() + 500 - start)/1000);
    start = timeInMillis();

    //
    // The save format is:
    //  file 'GenomeIndex' contains in order nHashTables, overflowTableSize, seedLen.
    //  File 'overflowTable' overflowTableSize bytes of the overflow table.
    //  Each hash table is saved in file base name 'GenomeIndexHash%d' where %d is the
    //  table number.
    //  And the genome itself is already saved in the same directory in its own format.
    //
    snprintf(filenameBuffer,filenameBufferSize,"%s%cGenomeIndex",directoryName,PATH_SEP);

    FILE *indexFile = fopen(filenameBuffer,"w");
    if (indexFile == NULL) {
        fprintf(stderr,"Unable to open file '%s' for write.\n",filenameBuffer);
        return false;
    }

    fprintf(indexFile,"%d %d %d",index->nHashTables, index->overflowTableSize, seedLen);

    fclose(indexFile);

    snprintf(filenameBuffer,filenameBufferSize,"%s%cOverflowTable",directoryName,PATH_SEP);
    FILE* fOverflowTable = fopen(filenameBuffer, "wb");
    if (fOverflowTable == NULL) {
        fprintf(stderr,"Unable to open overflow table file, '%s', %d\n",filenameBuffer,errno);
        return false;
    }

    const unsigned writeSize = 32 * 1024 * 1024;
    for (size_t writeOffset = 0; writeOffset < index->overflowTableSize * sizeof(*index->overflowTable); ) {
        unsigned amountToWrite = (unsigned)__min((size_t)writeSize,(size_t)index->overflowTableSize * sizeof(*index->overflowTable) - writeOffset);
#ifdef _MSC_VER
        if (amountToWrite % 4096) {
            amountToWrite = ((amountToWrite + 4095) / 4096) * 4096;
            if (amountToWrite + writeOffset > index->overflowTableVirtualAllocSize) {
                fprintf(stderr,"GenomeIndex::saveToDirectory: overflow table virtual alloc size doesn't appear to be a multiple of the page size %lld\n",
                    (_int64) index->overflowTableVirtualAllocSize);
            }
        }
#endif
        size_t amountWritten = fwrite(((char *)index->overflowTable) + writeOffset, 1, amountToWrite, fOverflowTable);
        if (amountWritten < amountToWrite) {
            fprintf(stderr,"GenomeIndex::saveToDirectory: fwrite failed, %d\n",errno);
            fclose(fOverflowTable);
            return false;
        }
        writeOffset += amountWritten;
    }
    fclose(fOverflowTable);
    fOverflowTable = NULL;

    for (unsigned i = 0; i < index->nHashTables; i++) {
        snprintf(filenameBuffer,filenameBufferSize,"%s%cGenomeIndexHash%04d",directoryName,PATH_SEP,i);
        if (!hashTables[i]->saveToFile(filenameBuffer)) {
            fprintf(stderr,"GenomeIndex::saveToDirectory: Failed to save hash table %d\n",i);
            return false;
        }
    }

    delete index;
    if (biasTable != NULL) {
        delete[] biasTable;
    }
    
    printf("%llds\n", (timeInMillis() + 500 - start) / 1000);
    
    return true;
}

//
// A comparison method for qsort that sorts unsigned ints backwards (SNAP expects them to be backwards due to 
// a historical artifact).
//
int
GenomeIndex::BackwardsUnsignedCompare(const void *first, const void *second)
{
    if (*(const unsigned *) first > *(const unsigned *)second) {
        return -1;
    } else if (*(const unsigned *) first == *(const unsigned *)second) {
        return 0;
    } else {
        return 1;
    }
}

    void
GenomeIndex::AddOverflowBackpointer(
    OverflowEntry       *overflowEntry,
    OverflowBackpointer *overflowBackpointers,
    unsigned             nOverflowBackpointers,
    volatile unsigned   *nextOverflowBackpointer,
    unsigned             genomeOffset)
{
    unsigned overflowBackpointerIndex = InterlockedIncrementAndReturnNewValue(nextOverflowBackpointer) - 1;
    if (nOverflowBackpointers <= overflowBackpointerIndex) {
        fprintf(stderr,"Ran out of overflow backpointers.  Consider using the -O switch with a larger value to increase space.\n");
        exit(1);
    }
    OverflowBackpointer *newBackpointer = &overflowBackpointers[overflowBackpointerIndex];
 
    newBackpointer->nextIndex = overflowEntry->backpointerIndex;
    newBackpointer->genomeOffset = genomeOffset;
    overflowEntry->backpointerIndex = overflowBackpointerIndex;
    overflowEntry->nInstances++;
}

    double
GenomeIndex::GetHashTableSizeBias(
    unsigned        whichTable,
    int             seedLen)
/*++

Routine Description:

    Provide the bias amount for a given hash table.

    The index uses multiple hash tables as a way to keep the hash table key size no more
    than 32 bits.  When we have seeds > 32 bits, we divide the space into 32 bit chunks and
    have one hash table for each.  The obvious thing to do is that we have n hash tables we should
    have each one be totalSize/n.  However, the human genome is highly nonuniform, and if we
    do that we have some tables overflow and others wind up underused.  So, we take the
    frequency of occurence of the various base sequences into account when we create the hash tables,
    resulting in more-or-less even loading.

    This method returns the bias to apply for a particular table/seed size combo.  Because it's empirically measured,
    it's just a giant lookup table.  If we haven't measured the particular table size, we just return
    an even weighting, which will almost certainly result in overflowing one of the tables (for 
    size 20 seeds the most popular table is 37 times bigger than the least popular).  This will run
    until you overflow, at which time the IndexBuilder will dump out the weights it's got, which will
    allow you to make a provisional entry in this table, which will get you to the final weights.

    There's one additional twist: because of duplications in the genome (and N in the reference genome)
    there will be fewer entries in the hash table than bases in the genome.  We apply this bias to the
    result as well, so the mean of all of the entries will be about 72%, not 100% like you'd otherwise
    expect.  This saves memory (or, alternately, makes the slack parameter more accurate, take your
    pick).

Arguments:

    whichTable      - which table's bias do we want
    seedLen         - what's the seed size (and hence the number of tables overall)

--*/
{
    switch (seedLen) {

        case 23: {
            _ASSERT(whichTable < 16384);
            return biasTable23[whichTable];
        }

        case 22: {
            _ASSERT(whichTable < 4096);
            return biasTable22[whichTable];
        }

        case 21: {
            _ASSERT(whichTable < 1024);
            return biasTable21[whichTable];
        }

        case 20: {
            _ASSERT(whichTable < 256);
            return biasTable20[whichTable];
        }

        case 19: {  // hg19
            _ASSERT(whichTable < 64);
            return biasTable19[whichTable];
        }
        case 18: {// hg19
            _ASSERT(whichTable < 16);
            return biasTable18[whichTable];
        }

        case 17: {
            _ASSERT(whichTable < 4);
            return biasTable17[whichTable];
        }
        default:    return 0.6; // We don't have emperical data for this type yet.
    }

    fprintf(stderr,"GetHashTableSizeBias: bogus parameters, %d, %d\n",whichTable,seedLen);
    exit(1);
    return 0;   // Just to make the compiler happy.
}

GenomeIndex::GenomeIndex() : nHashTables(0), hashTables(NULL), overflowTable(NULL), genome(NULL)
{
}

    GenomeIndex *
GenomeIndex::loadFromDirectory(char *directoryName)
{
    GenomeIndex *index = new GenomeIndex();

    const unsigned filenameBufferSize = MAX_PATH+1;
    char filenameBuffer[filenameBufferSize];
    
    snprintf(filenameBuffer,filenameBufferSize,"%s%cGenomeIndex",directoryName,PATH_SEP);

    FILE *indexFile = fopen(filenameBuffer,"r");
    if (indexFile == NULL) {
        fprintf(stderr,"Unable to open file '%s' for read.\n",filenameBuffer);
        return false;
    }

    unsigned seedLen;
    if (3 != fscanf(indexFile,"%d %d %d",&index->nHashTables, &index->overflowTableSize, &seedLen)) {
        fprintf(stderr,"GenomeIndex::LoadFromDirectory: didn't read initial values\n");
        fclose(indexFile);
        delete index;
        return NULL;
    }
    fclose(indexFile);

    if (0 == seedLen) {
        fprintf(stderr,"GenomeIndex::LoadFromDirectory: saw seed size of 0.\n");
        delete index;
        return NULL;
    }
    index->seedLen = seedLen;

    index->overflowTable = (unsigned *)BigAlloc(index->overflowTableSize * sizeof(*(index->overflowTable)),&index->overflowTableVirtualAllocSize);

    snprintf(filenameBuffer,filenameBufferSize,"%s%cOverflowTable",directoryName,PATH_SEP);
    FILE* fOverflowTable = fopen(filenameBuffer, "rb");
    if (fOverflowTable == NULL) {
        fprintf(stderr,"Unable to open overflow table file, '%s', %d\n",filenameBuffer,errno);
        delete index;
        return NULL;
    }

    const unsigned readSize = 32 * 1024 * 1024;
    for (size_t readOffset = 0; readOffset < index->overflowTableSize * sizeof(*(index->overflowTable)); ) {
        int amountToRead = (unsigned)__min((size_t)readSize,(size_t)index->overflowTableSize * sizeof(*(index->overflowTable)) - readOffset);
#ifdef _MSC_VER
        if (amountToRead % 4096) {
            amountToRead = ((amountToRead + 4095) / 4096) * 4096;
            if (amountToRead + readOffset > index->overflowTableVirtualAllocSize) {
                fprintf(stderr,"GenomeIndex::loadFromDirectory: overflow table virtual alloc size doesn't appear to be a multiple of the page size %lld\n",
                    (_int64) index->overflowTableVirtualAllocSize);
            }
        }
#endif
        int amountRead = (int)fread(((char*) index->overflowTable) + readOffset, 1, amountToRead, fOverflowTable);   
        if (amountRead < amountToRead) {
            fprintf(stderr,"GenomeIndex::loadFromDirectory: fread failed (amountToRead = %d, amountRead = %d, readOffset %lld), %d\n",amountToRead, amountRead, readOffset, errno);
            fclose(fOverflowTable);
            delete index;
            return NULL;
        }
        readOffset += amountRead;
    }
    fclose(fOverflowTable);
    fOverflowTable = NULL;

    index->hashTables = new SNAPHashTable*[index->nHashTables];

    for (unsigned i = 0; i < index->nHashTables; i++) {
        index->hashTables[i] = NULL; // We need to do this so the destructor doesn't crash if loading a hash table fails.
    }

    for (unsigned i = 0; i < index->nHashTables; i++) {
        snprintf(filenameBuffer,filenameBufferSize,"%s%cGenomeIndexHash%04d",directoryName,PATH_SEP,i);
        if (NULL == (index->hashTables[i] = new SNAPHashTable(filenameBuffer))) {
            fprintf(stderr,"GenomeIndex::loadFromDirectory: Failed to load hash table %d\n",i);
            delete index;
            return NULL;
        }
    }

    snprintf(filenameBuffer,filenameBufferSize,"%s%cGenome",directoryName,PATH_SEP);
    if (NULL == (index->genome = Genome::loadFromFile(filenameBuffer))) {
        fprintf(stderr,"GenomeIndex::loadFromDirectory: Failed to load the genome itself\n");
        delete index;
        return NULL;
    }

    return index;
}

    void
GenomeIndex::lookupSeed(Seed seed, unsigned *nHits, const unsigned **hits, unsigned *nRCHits, const unsigned **rcHits)
{
    return lookupSeed(seed, 0, 0xFFFFFFFF, nHits, hits, nRCHits, rcHits);
}

    void
GenomeIndex::lookupSeed(
    Seed              seed,
    unsigned          minLocation,
    unsigned          maxLocation,
    unsigned         *nHits,
    const unsigned  **hits,
    unsigned         *nRCHits,
    const unsigned  **rcHits)
{
    bool lookedUpComplement;

    lookedUpComplement = seed.isBiggerThanItsReverseComplement();
    if (lookedUpComplement) {
        seed = ~seed;
    }

    _ASSERT(seed.getHighBases() < nHashTables);
    unsigned lowBases = seed.getLowBases();
    unsigned *entry = hashTables[seed.getHighBases()]->Lookup(lowBases);
    if (NULL == entry) {
        *nHits = 0;
        *nRCHits = 0;
        return;
    }

    //
    // Fill in the caller's answers for the main and complement of the seed looked up.
    // Because of our hash table design, we may have had to take the complement before the
    // lookup, in which case we reverse the results so the caller gets the right thing.
    // Also, if the seed is its own reverse complement, we need to fill the same hits
    // in both return arrays.
    //
    fillInLookedUpResults((lookedUpComplement ? entry + 1 : entry), minLocation, maxLocation, nHits, hits);
    if (seed.isOwnReverseComplement()) {
      *nRCHits = *nHits;
      *rcHits = *hits;
    } else {
      fillInLookedUpResults((lookedUpComplement ? entry : entry + 1), minLocation, maxLocation, nRCHits, rcHits);
    }
}

    void
GenomeIndex::fillInLookedUpResults(
    unsigned        *subEntry,
    unsigned         minLocation,
    unsigned         maxLocation,
    unsigned        *nHits, 
    const unsigned **hits)
{
    if (*subEntry < genome->getCountOfBases()) {
        //
        // It's a singleton.
        //
        *nHits = (*subEntry >= minLocation && *subEntry <= maxLocation) ? 1 : 0;
        *hits = subEntry;
    } else if (*subEntry == 0xfffffffe) {
        //
        // It's unused, the other complement must exist.
        //
        *nHits = 0;
    } else {
        //
        // Multiple hits.  Recall that the overflow table format is first a count of
        // the number of hits for that seed, followed by the list of hits.
        //
        unsigned overflowTableOffset = *subEntry - genome->getCountOfBases();

        _ASSERT(overflowTableOffset < overflowTableSize);

        int hitCount = overflowTable[overflowTableOffset];

        _ASSERT(hitCount >= 2);
        _ASSERT(hitCount + overflowTableOffset < overflowTableSize);

        if (minLocation == 0 && maxLocation == 0xFFFFFFFF) {
            // Return all the hits
            *nHits = hitCount;
            *hits = &overflowTable[overflowTableOffset + 1];
        } else {
            // Do a binary search to find hits in the right range of locations, taking advantage
            // of the knowledge that the hit list is sorted in descending order. Specifically,
            // we'll do a binary search to find the first hit that is <= maxLocation, and then
            // a linear search forward from that to find the other ones (assuming they are few).
            unsigned *allHits = &overflowTable[overflowTableOffset + 1];
            int low = 0;
            int high = hitCount - 1;
            while (low < high - 32) {
                int mid = low + (high - low) / 2;
                if (allHits[mid] <= maxLocation) {
                    high = mid;
                } else {
                    low = mid + 1;
                }
            }
            // We stop the binary search early and do a linear scan for speed (to avoid branches).
            while (low < hitCount && allHits[low] > maxLocation) {
                low++;
            }
            int end = low;
            while (end < hitCount && allHits[end] >= minLocation) {
                end++;
            }
            *nHits = end - low;
            *hits = allHits + low;
        }
    }
}

GenomeIndex::~GenomeIndex()
{
    if (NULL != hashTables) {
        for (unsigned i = 0; i < nHashTables; i++) {
            delete hashTables[i];
            hashTables[i] = NULL;
        }
    }

    delete [] hashTables;
    hashTables = NULL;

    if (NULL != overflowTable) {
        BigDealloc(overflowTable);
        overflowTable = NULL;
    }

    delete genome;
    genome = NULL;
}

    void
GenomeIndex::ComputeBiasTable(const Genome* genome, int seedLen, double* table, unsigned maxThreads)
/**
 * Fill in table with the table size biases for a given genome and seed size.
 * We assume that table is already of the correct size for our seed size
 * (namely 4**(seedLen-16)), and just fill in the values.
 *
 * If the genome is less than 2^20 bases, we count the seeds in each table exactly;
 * otherwise, we estimate them using Flajolet-Martin approximate counters.
 */
{
    _int64 start = timeInMillis();
    printf("Computing bias table.\n");

    unsigned nHashTables = (seedLen <= 16 ? 1 : 1 << ((seedLen - 16) * 2));
    unsigned countOfBases = genome->getCountOfBases();

    static const unsigned GENOME_SIZE_FOR_EXACT_COUNT = 1 << 20;  // Needs to be a power of 2 for hash sets

    bool computeExactly = (countOfBases < GENOME_SIZE_FOR_EXACT_COUNT);
    FixedSizeVector<unsigned> numExactSeeds(nHashTables, 0);
    FixedSizeSet<_int64> exactSeedsSeen(2 * GENOME_SIZE_FOR_EXACT_COUNT);
    vector<ApproximateCounter> approxCounters(nHashTables);

    _int64 validSeeds = 0;

    if (computeExactly) {
        for (unsigned i = 0; i < (unsigned)(countOfBases - seedLen); i++) {
            if (i % 1000000 == 0) {
                printf("Bias computation: %lld / %lld\n",(_int64)i, (_int64)countOfBases);
            }
            const char *bases = genome->getSubstring(i,seedLen);
            //
            // Check it for NULL, because Genome won't return strings that cross piece (chromosome) boundaries.
            //
            if (NULL == bases) {
                continue;
            }

            //
            // We don't build seeds out of sections of the genome that contain 'N.'  If this is one, skip it.
            //
            if (!Seed::DoesTextRepresentASeed(bases, seedLen)) {
                continue;
            }

            Seed seed(bases, seedLen);
            validSeeds++;

            //
            // Figure out if we're using this base or its complement.
            //
            Seed rc = ~seed;
            bool usingComplement = seed.isBiggerThanItsReverseComplement();
            if (usingComplement) {
                seed = ~seed;       // Couldn't resist using ~ for this.
            }

            _ASSERT(seed.getHighBases() < nHashTables);
            if (computeExactly) {
                if (!exactSeedsSeen.contains(seed.getBases())) {
                    exactSeedsSeen.add(seed.getBases());
                    numExactSeeds[seed.getHighBases()]++;
                }
            } else {
                approxCounters[seed.getHighBases()].add(seed.getLowBases());
            }
        }
    } else {
        //
        // Run through the table in parallel.
        //
        unsigned nThreads = __min(GetNumberOfProcessors(), maxThreads);
        volatile int runningThreadCount = nThreads;
        volatile _int64 nBasesProcessed = 0;
        SingleWaiterObject doneObject;

        ExclusiveLock *locks;
        locks = new ExclusiveLock[nHashTables];
        for (unsigned i = 0; i < nHashTables; i++) {
            InitializeExclusiveLock(&locks[i]);
        }

        CreateSingleWaiterObject(&doneObject);

        ComputeBiasTableThreadContext *contexts = new ComputeBiasTableThreadContext[nThreads];
        unsigned nextChunkToProcess = 0;
        for (unsigned i = 0; i < nThreads; i++) {
            contexts[i].approxCounters = &approxCounters;
            contexts[i].doneObject = &doneObject;
            contexts[i].genomeChunkStart = nextChunkToProcess;
            if (i == nThreads - 1) {
                nextChunkToProcess = countOfBases - seedLen - 1;
            } else {
                nextChunkToProcess += (countOfBases - seedLen) / nThreads;
            }
            contexts[i].genomeChunkEnd = nextChunkToProcess;
            contexts[i].nHashTables = nHashTables;
            contexts[i].runningThreadCount = &runningThreadCount;
            contexts[i].genome = genome;
            contexts[i].nBasesProcessed = &nBasesProcessed;
            contexts[i].seedLen = seedLen;
            contexts[i].validSeeds = &validSeeds;
            contexts[i].approximateCounterLocks = locks;

            StartNewThread(ComputeBiasTableWorkerThreadMain, &contexts[i]);
        }
 
        WaitForSingleWaiterObject(&doneObject);
        DestroySingleWaiterObject(&doneObject);

        for (unsigned i = 0; i < nHashTables; i++) {
            DestroyExclusiveLock(&locks[i]);
        }
        delete [] locks;
    }


    double distinctSeeds = 0;
    for (unsigned i = 0; i < nHashTables; i++) {
        distinctSeeds += computeExactly ? numExactSeeds[i] : approxCounters[i].getCount();
    }

    for (unsigned i = 0; i < nHashTables; i++) {
        int count = computeExactly ? numExactSeeds[i] : approxCounters[i].getCount();
        table[i] = (count / distinctSeeds) * ((double)validSeeds / countOfBases) * nHashTables;
    }

    // printf("Bias table:\n");
    // for (unsigned i = 0; i < nHashTables; i++) {
    //     printf("%u -> %lf\n", i, table[i]);
    // }

    printf("Computed bias table in %llds\n", (timeInMillis() + 500 - start) / 1000);
}

struct PerCounterBatch {
    PerCounterBatch() : nUsed(0) {}
    static const unsigned nSeedsPerBatch = 1000;

    unsigned nUsed;
    unsigned lowBases[nSeedsPerBatch];

    bool addSeed(unsigned seedLowBases) {
        _ASSERT(nUsed < nSeedsPerBatch);
        lowBases[nUsed] = seedLowBases;
        nUsed++;
        return nUsed >= nSeedsPerBatch;
    }

    void apply(ApproximateCounter *counter) {
        for (unsigned i = 0; i < nUsed; i++) {
            counter->add(lowBases[i]);
        }
        nUsed = 0;
    }
};




    void
GenomeIndex::ComputeBiasTableWorkerThreadMain(void *param)
{
    ComputeBiasTableThreadContext *context = (ComputeBiasTableThreadContext *)param;

    unsigned countOfBases = context->genome->getCountOfBases();
    _int64 validSeeds = 0;

    //
    // Batch the insertions into the approximate counters, because otherwise we spend all of
    // our time acquiring and releasing locks.
    //
 
    PerCounterBatch *batches = new PerCounterBatch[context->nHashTables];

    _uint64 unrecordedSkippedSeeds = 0;

    for (unsigned i = context->genomeChunkStart; i < context->genomeChunkEnd; i++) {
            const char *bases = context->genome->getSubstring(i, context->seedLen);
            //
            // Check it for NULL, because Genome won't return strings that cross piece (chromosome) boundaries.
            //
            if (NULL == bases) {
                continue;
            }

            //
            // We don't build seeds out of sections of the genome that contain 'N.'  If this is one, skip it.
            //
            if (!Seed::DoesTextRepresentASeed(bases, context->seedLen)) {
                unrecordedSkippedSeeds++;
                continue;
            }

            Seed seed(bases, context->seedLen);
            validSeeds++;

            //
            // Figure out if we're using this base or its complement.
            //
            Seed rc = ~seed;
            bool usingComplement = seed.isBiggerThanItsReverseComplement();
            if (usingComplement) {
                seed = ~seed;       // Couldn't resist using ~ for this.
            }

            unsigned whichHashTable = seed.getHighBases();

            _ASSERT(whichHashTable < context->nHashTables);

            if (batches[whichHashTable].addSeed(seed.getLowBases())) {
                PerCounterBatch *batch = &batches[whichHashTable];
                AcquireExclusiveLock(&context->approximateCounterLocks[whichHashTable]);
                batch->apply(&(*context->approxCounters)[whichHashTable]);    
                ReleaseExclusiveLock(&context->approximateCounterLocks[whichHashTable]);

                _int64 basesProcessed = InterlockedAdd64AndReturnNewValue(context->nBasesProcessed, PerCounterBatch::nSeedsPerBatch + unrecordedSkippedSeeds);

                const _uint64 printBatchSize = 100000000;
                if ((_uint64)basesProcessed / printBatchSize > ((_uint64)basesProcessed - PerCounterBatch::nSeedsPerBatch - unrecordedSkippedSeeds)/printBatchSize) {
                    printf("Bias computation: %lld / %lld\n",(basesProcessed/printBatchSize)*printBatchSize, (_int64)countOfBases);
                }
                unrecordedSkippedSeeds= 0;  // We've now recorded them.
            }
    }

    for (unsigned i = 0; i < context->nHashTables; i++) {
        AcquireExclusiveLock(&context->approximateCounterLocks[i]);
        batches[i].apply(&(*context->approxCounters)[i]);
        ReleaseExclusiveLock(&context->approximateCounterLocks[i]);
    }

    delete [] batches;

    InterlockedAdd64AndReturnNewValue(context->validSeeds, validSeeds);


    if (0 == InterlockedDecrementAndReturnNewValue(context->runningThreadCount)) {
        SignalSingleWaiterObject(context->doneObject);
    }
}

struct PerHashTableBatch {
    PerHashTableBatch() : nUsed(0) {}

    static const unsigned nSeedsPerBatch = 1000;

    unsigned nUsed;

    struct Entry {
        bool usingComplement;
        unsigned lowBases;
        unsigned genomeLocation;
    };

    Entry entries[nSeedsPerBatch];

    bool addSeed(unsigned genomeLocation, unsigned seedLowBases, bool seedUsingComplement) {
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




    void 
GenomeIndex::BuildHashTablesWorkerThreadMain(void *param)
{
    BuildHashTablesThreadContext *context = (BuildHashTablesThreadContext *)param;

    unsigned countOfBases = context->genome->getCountOfBases();
    const Genome *genome = context->genome;
    unsigned seedLen = context->seedLen;
    _int64 noBaseAvailable = 0;
    _int64 nonSeeds = 0;
    _int64 bothComplementsUsed = 0;
    _int64 countOfDuplicateOverflows = 0;
 
    GenomeIndex *index = context->index;
    unsigned nHashTables = index->nHashTables;

    //
    // Batch the insertions into the hash tables, because otherwise we spend all of
    // our time acquiring and releasing locks.
    //
 
    PerHashTableBatch *batches = new PerHashTableBatch[index->nHashTables];

    _uint64 unrecordedSkippedSeeds = 0;

    for (unsigned genomeLocation = context->genomeChunkStart; genomeLocation < context->genomeChunkEnd; genomeLocation++) {

        const char *bases = genome->getSubstring(genomeLocation, seedLen);
        //
        // Check it for NULL, because Genome won't return strings that cross piece (chromosome) boundaries.
        //
        if (NULL == bases) {
            noBaseAvailable++;
            unrecordedSkippedSeeds++;
            continue;
        }

        //
        // We don't build seeds out of sections of the genome that contain 'N.'  If this is one, skip it.
        //
        if (!Seed::DoesTextRepresentASeed(bases, seedLen)) {
            nonSeeds++;
            unrecordedSkippedSeeds++;
            continue;
        }

        Seed seed(bases, seedLen);

        //
        // Figure out if we're using this base or its complement.
        //
        bool usingComplement = seed.isBiggerThanItsReverseComplement();
        if (usingComplement) {
            seed = ~seed;       // Couldn't resist using ~ for this.
        }

       unsigned whichHashTable = seed.getHighBases();
       _ASSERT(whichHashTable < index->nHashTables);
 
        if (batches[whichHashTable].addSeed(genomeLocation, seed.getLowBases(), usingComplement)) {
            AcquireExclusiveLock(&context->hashTableLocks[whichHashTable]);
            for (unsigned i = 0; i < batches[whichHashTable].nUsed; i++) {
                ApplyHashTableUpdate(context, whichHashTable, batches[whichHashTable].entries[i].genomeLocation, 
                    batches[whichHashTable].entries[i].lowBases, batches[whichHashTable].entries[i].usingComplement,
                    &bothComplementsUsed, &countOfDuplicateOverflows);
            }
            ReleaseExclusiveLock(&context->hashTableLocks[whichHashTable]);

            _int64 newNBasesProcessed = InterlockedAdd64AndReturnNewValue(context->nBasesProcessed, batches[whichHashTable].nUsed + unrecordedSkippedSeeds);
            const _int64 printPeriod = 100000000;
            if ((unsigned)(newNBasesProcessed / printPeriod) > (unsigned)((newNBasesProcessed - batches[whichHashTable].nUsed - unrecordedSkippedSeeds) / printPeriod)) {
                printf("Indexing %lld / %lld\n", (newNBasesProcessed / printPeriod) * printPeriod, countOfBases);
            }
            unrecordedSkippedSeeds = 0;
            batches[whichHashTable].clear();
        } // If we filled a batch
    } // For each genome base in our area

    //
    // Now apply the updates from the batches that were left over
    //

    for (unsigned whichHashTable = 0; whichHashTable < nHashTables; whichHashTable++) {
        AcquireExclusiveLock(&context->hashTableLocks[whichHashTable]);
        for (unsigned i = 0; i < batches[whichHashTable].nUsed; i++) {
            ApplyHashTableUpdate(context, whichHashTable, batches[whichHashTable].entries[i].genomeLocation, 
                batches[whichHashTable].entries[i].lowBases, batches[whichHashTable].entries[i].usingComplement,
                &bothComplementsUsed, &countOfDuplicateOverflows);
        }
        ReleaseExclusiveLock(&context->hashTableLocks[whichHashTable]);
    }

    InterlockedAdd64AndReturnNewValue(context->noBaseAvailable, noBaseAvailable);
    InterlockedAdd64AndReturnNewValue(context->nonSeeds, nonSeeds);
    InterlockedAdd64AndReturnNewValue(context->bothComplementsUsed, bothComplementsUsed);
    InterlockedAdd64AndReturnNewValue(context->countOfDuplicateOverflows, countOfDuplicateOverflows);

    delete [] batches;

    if (0 == InterlockedDecrementAndReturnNewValue(context->runningThreadCount)) {
        SignalSingleWaiterObject(context->doneObject);
    }

}

    void 
GenomeIndex::ApplyHashTableUpdate(BuildHashTablesThreadContext *context, unsigned whichHashTable, unsigned genomeLocation, unsigned lowBases, bool usingComplement,
                _int64 *bothComplementsUsed, _int64 *countOfDuplicateOverflows)
{
    GenomeIndex *index = context->index;
    unsigned countOfBases = context->genome->getCountOfBases();
    SNAPHashTable *hashTable = index->hashTables[whichHashTable];
    unsigned *entry = hashTable->SlowLookup(lowBases);  // use SlowLookup because we might have overflowed the table.
    if (NULL == entry) {
        //
        // We haven't yet seen either this seed or its complement.  Make a new hash table
        // entry.
        //
        unsigned newEntry[2];
        if (!usingComplement) {
            newEntry[0] = genomeLocation;
            newEntry[1] = 0xfffffffe; // Use 0xfffffffe for unused, because we gave 0xffffffff to the hash table package.
        } else{
            newEntry[0] = 0xfffffffe; // Use 0xfffffffe for unused, because we gave 0xffffffff to the hash table package.
            newEntry[1] = genomeLocation;
        }

        if (!hashTable->Insert(lowBases,newEntry)) {
            for (unsigned j = 0; j < index->nHashTables; j++) {
                printf("HashTable[%d] has %lld used elements\n",j,(_int64)index->hashTables[j]->GetUsedElementCount());
            }
            printf("IndexBuilder: exceeded size of hash table %d.\n"
                    "If you're indexing a non-human genome, make sure not to pass the -hg19 option.\n",
                    whichHashTable);
            exit(1);
        }
    } else {
        //
        // This entry already exists in the hash table.  It might just be because we've already seen the seed's complement
        // in which case we update our half of the entry.  Otherwise, it's a repeat in the genome, and we need to insert
        // it in the overflow table.
        //
        int entryIndex = usingComplement ? 1 : 0;
        if (0xfffffffe == entry[entryIndex]) {
            entry[entryIndex] = genomeLocation;
            (*bothComplementsUsed)++;
        } else if (entry[entryIndex] < countOfBases) {
            //
            // Allocate an overflow table entry.
            //
            unsigned overflowIndex = InterlockedIncrementAndReturnNewValue(context->nextOverflowIndex) - 1;
            if (overflowIndex >= context->nOverflowEntries) {
                if (0xffffffff - overflowIndex < 10) {
                    fprintf(stderr,"You've run out of index overflow address space.  Perhaps a larger seed size will help.\n");
                } else {
                    fprintf(stderr,"Index builder: Overflowed overflow table.  Consider using the -O switch with a larger value to make more space.\n");
                }
                exit(1);
            }
            //
            // And add two backpointers: one for what was already in the hash table and one
            // for the index we're processing now.
            //
            OverflowEntry *overflowEntry = &context->overflowEntries[overflowIndex];
            overflowEntry->hashTableEntry = entry + entryIndex;
            overflowEntry->nInstances = 0;
            overflowEntry->backpointerIndex = -1;

            AddOverflowBackpointer(overflowEntry, context->overflowBackpointers, context->nOverflowBackpointers, context->nextOverflowBackpointer, entry[entryIndex]);
            AddOverflowBackpointer(overflowEntry, context->overflowBackpointers, context->nOverflowBackpointers, context->nextOverflowBackpointer, genomeLocation);

            entry[entryIndex] = overflowIndex + countOfBases;
        } else {
            //
            // Stick another entry in the existing overflow bucket.
            //
            _ASSERT(entry[entryIndex] - countOfBases < context->nOverflowEntries && entry[entryIndex] >= countOfBases);
            OverflowEntry *overflowEntry = &context->overflowEntries[entry[entryIndex] - countOfBases];

            _ASSERT(overflowEntry->hashTableEntry == entry + entryIndex);

            AddOverflowBackpointer(overflowEntry, context->overflowBackpointers, context->nOverflowBackpointers, context->nextOverflowBackpointer, genomeLocation);

            (*countOfDuplicateOverflows)++;    // Should add extra stuff to the overflow table here.
        } // If the existing entry had the complement empty, needed a new overflow entry or extended an old one
    } // If new or existing entry.
}