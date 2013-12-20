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
#include "exit.h"
#include "directions.h"
#include "GTFReader.h"

using namespace std;

static const int DEFAULT_SEED_SIZE = 20;
static const double DEFAULT_SLACK = 0.3;
static const unsigned DEFAULT_PADDING = 500;
static const unsigned DEFAULT_KEY_BYTES = 4;


static void tusage()
{
    fprintf(stderr,
            "Usage: snapr transcriptome <input.gtf> <input.fa> <output-dir> [<options>]\n"
            "Options:\n"
            "  -s               Seed size (default: %d)\n"
            "  -h               Hash table slack (default: %.1f)\n"
            "  -hg19            Use pre-computed table bias for hg19, which results in better speed, balance, and memory footprint but may not work for other references.\n"
            "  -Ofactor         Specify the size of the overflow space.  This will change the memory footprint, but may be needed for some genomes.\n"
            "                   Larger numbers use more memory but work better with more repetitive genomes.  Smaller numbers reduce the memory\n"
            "                   footprint, but may cause the index build to fail.  Making -O larger than necessary will not affect the resuting\n"
            "                   index.  Factor must be between 1 and 1000, and the default is 50.\n"
            " -tMaxThreads      Specify the maximum number of threads to use. Default is the number of cores.\n"
            " -pPadding         Specify the number of Ns to put as padding between chromosomes.  This must be as large as the largest\n"
            "                   edit distance you'll ever use, and there's a performance advantage to have it be bigger than any\n"
            "                   read you'll process.  Default is %d\n"
            " -HHistogramFile   Build a histogram of seed popularity.  This is just for information, it's not used by SNAP.\n"
            " -exact            Compute hash table sizes exactly.  This will slow down index build, but may be necessary in some cases\n"
            " -keysize          The number of bytes to use for the hash table key.  Larger values increase SNAP's memory footprint, but allow larger seeds.  Default: %d\n",
            DEFAULT_SEED_SIZE,
            DEFAULT_SLACK,
            DEFAULT_PADDING,
            DEFAULT_KEY_BYTES);
    soft_exit(1);
}

static void usage()
{
    fprintf(stderr,
            "Usage: snapr index <input.fa> <output-dir> [<options>]\n"
            "Options:\n"
            "  -s               Seed size (default: %d)\n"
            "  -h               Hash table slack (default: %.1f)\n"
            "  -hg19            Use pre-computed table bias for hg19, which results in better speed, balance, and memory footprint but may not work for other references.\n"
            "  -Ofactor         Specify the size of the overflow space.  This will change the memory footprint, but may be needed for some genomes.\n"
            "                   Larger numbers use more memory but work better with more repetitive genomes.  Smaller numbers reduce the memory\n"
            "                   footprint, but may cause the index build to fail.  Making -O larger than necessary will not affect the resuting\n"
            "                   index.  Factor must be between 1 and 1000, and the default is 50.\n"
            " -tMaxThreads      Specify the maximum number of threads to use. Default is the number of cores.\n"
            " -pPadding         Specify the number of Ns to put as padding between chromosomes.  This must be as large as the largest\n"
            "                   edit distance you'll ever use, and there's a performance advantage to have it be bigger than any\n"
            "                   read you'll process.  Default is %d\n"
            " -HHistogramFile   Build a histogram of seed popularity.  This is just for information, it's not used by SNAP.\n"
            " -exact            Compute hash table sizes exactly.  This will slow down index build, but may be necessary in some cases\n"
            " -keysize          The number of bytes to use for the hash table key.  Larger values increase SNAP's memory footprint, but allow larger seeds.  Default: %d\n",
            DEFAULT_SEED_SIZE,
            DEFAULT_SLACK,
            DEFAULT_PADDING,
            DEFAULT_KEY_BYTES);
    soft_exit(1);
}


    void
GenomeIndex::runTranscriptomeIndexer(
    int argc,
    const char **argv)
{
    if (argc < 3) {
        tusage();
    }

    const char *gtfFile = argv[0];
    const char *fastaFile = argv[1];
    const char *outputDir = argv[2];

    unsigned maxThreads = GetNumberOfProcessors();

    int seedLen = DEFAULT_SEED_SIZE;
    double slack = DEFAULT_SLACK;
    bool computeBias = true;
    _uint64 overflowTableFactor = 50;
    const char *histogramFileName = NULL;
    unsigned chromosomePadding = DEFAULT_PADDING;
    bool forceExact = false;
    unsigned keySizeInBytes = DEFAULT_KEY_BYTES;

    for (int n = 3; n < argc; n++) {
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
        } else if (strcmp(argv[n], "-exact") == 0) {
            forceExact = true;
        } else if (strcmp(argv[n], "-hg19") == 0) {
            computeBias = false;
        } else if (argv[n][0] == '-' && argv[n][1] == 'H') {
            histogramFileName = argv[n] + 2;
        } else if (argv[n][0] == '-' && argv[n][1] == 'O') {
            overflowTableFactor = atoi(argv[n]+2);
            if (overflowTableFactor < 1 || overflowTableFactor > 1000) {
                fprintf(stderr,"Overflow table factor must be between 1 and 1000 inclusive (and you need to not leave a space after '-O')\n");
                soft_exit(1);
            }
        } else if (argv[n][0] == '-' && argv[n][1] == 't') {
            maxThreads = atoi(argv[n]+2);
            if (maxThreads < 1 || maxThreads > 100) {
                fprintf(stderr,"maxThreads must be between 1 and 100 inclusive (and you need to not leave a space after '-t')\n");
                soft_exit(1);
            }
        } else if (argv[n][0] == '-' && argv[n][1] == 'p') {
            chromosomePadding = atoi(argv[n]+2);
            if (0 == chromosomePadding) {
                fprintf(stderr,"Invalid chromosome padding specified, must be at least one (and in practice as large as any max edit distance you might use).\n");
                soft_exit(1);
            }
        } else if (strcmp(argv[n], "-keysize") == 0) {
            if (n + 1 < argc) {
                keySizeInBytes = atoi(argv[n+1]);
                if (keySizeInBytes < 4 || keySizeInBytes > 8) {
                    fprintf(stderr, "Key size must be between 4 and 8 inclusive\n");
                    soft_exit(1);
                }
                n++;
            } else {
                usage();
            }
        } else {
            fprintf(stderr, "Invalid argument: %s\n\n", argv[n]);
            usage();
        }
    }

    if (seedLen < 16 || seedLen > 32) {
        // Seeds are stored in 64 bits, so they can't be larger than 32 bases for now.
        fprintf(stderr, "Seed length must be between 16 and 32, inclusive\n");
        soft_exit(1);
    }

    if (seedLen < 19 && !computeBias) {
        fprintf(stderr,"For hg19, you must use seed sizes between 19 and 25 (and not specifying -hg19 won't help, it'll just take longer to fail).\n");
        soft_exit(1);
    }

    printf("Hash table slack %lf\nLoading FASTA file '%s' into memory...", slack, fastaFile);
    _int64 start = timeInMillis();
    const Genome *genome = ReadFASTAGenome(fastaFile, chromosomePadding);
    if (NULL == genome) {
        fprintf(stderr, "Unable to read FASTA file\n");
        soft_exit(1);
    }
    printf("%llds\n", (timeInMillis() + 500 - start) / 1000);
           
    //First build the transcriptome file
    GTFReader gtf;
    gtf.Load(gtfFile);
    
    //Pass in the genome to the GTF object to write out the transcriptome file
    gtf.BuildTranscriptome(genome);
    
    //Replace the input FASTA file with the newly created transcriptome file
    argv[1] = "transcriptome.fa";
    GenomeIndex::runIndexer(argc-1, argv+1);
           
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
    _uint64 overflowTableFactor = 50;
    const char *histogramFileName = NULL;
    unsigned chromosomePadding = DEFAULT_PADDING;
    bool forceExact = false;
    unsigned keySizeInBytes = DEFAULT_KEY_BYTES;

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
        } else if (strcmp(argv[n], "-exact") == 0) {
            forceExact = true;
        } else if (strcmp(argv[n], "-hg19") == 0) {
            computeBias = false;
        } else if (argv[n][0] == '-' && argv[n][1] == 'H') {
            histogramFileName = argv[n] + 2;
        } else if (argv[n][0] == '-' && argv[n][1] == 'O') {
            overflowTableFactor = atoi(argv[n]+2);
            if (overflowTableFactor < 1 || overflowTableFactor > 1000) {
                fprintf(stderr,"Overflow table factor must be between 1 and 1000 inclusive (and you need to not leave a space after '-O')\n");
                soft_exit(1);
            }
        } else if (argv[n][0] == '-' && argv[n][1] == 't') {
            maxThreads = atoi(argv[n]+2);
            if (maxThreads < 1 || maxThreads > 100) {
                fprintf(stderr,"maxThreads must be between 1 and 100 inclusive (and you need to not leave a space after '-t')\n");
                soft_exit(1);
            }
        } else if (argv[n][0] == '-' && argv[n][1] == 'p') {
            chromosomePadding = atoi(argv[n]+2);
            if (0 == chromosomePadding) {
                fprintf(stderr,"Invalid chromosome padding specified, must be at least one (and in practice as large as any max edit distance you might use).\n");
                soft_exit(1);
            }
        } else if (strcmp(argv[n], "-keysize") == 0) {
            if (n + 1 < argc) {
                keySizeInBytes = atoi(argv[n+1]);
                if (keySizeInBytes < 4 || keySizeInBytes > 8) {
                    fprintf(stderr, "Key size must be between 4 and 8 inclusive\n");
                    soft_exit(1);
                }
                n++;
            } else {
                usage();
            }
        } else {
            fprintf(stderr, "Invalid argument: %s\n\n", argv[n]);
            usage();
        }
    }

    if (seedLen < 16 || seedLen > 32) {
        // Seeds are stored in 64 bits, so they can't be larger than 32 bases for now.
        fprintf(stderr, "Seed length must be between 16 and 32, inclusive\n");
        soft_exit(1);
    }

    if (seedLen < 19 && !computeBias) {
        fprintf(stderr,"For hg19, you must use seed sizes between 19 and 25 (and not specifying -hg19 won't help, it'll just take longer to fail).\n");
        soft_exit(1);
    }

    printf("Hash table slack %lf\nLoading FASTA file '%s' into memory...", slack, fastaFile);
    _int64 start = timeInMillis();
    const Genome *genome = ReadFASTAGenome(fastaFile, chromosomePadding);
    if (NULL == genome) {
        fprintf(stderr, "Unable to read FASTA file\n");
        soft_exit(1);
    }
    printf("%llds\n", (timeInMillis() + 500 - start) / 1000);
    unsigned nBases = genome->getCountOfBases();
    if (!GenomeIndex::BuildIndexToDirectory(genome, seedLen, slack, computeBias, outputDir, overflowTableFactor, maxThreads, chromosomePadding, forceExact, keySizeInBytes, histogramFileName)) {
        fprintf(stderr, "Genome index build failed\n");
        soft_exit(1);
    }
    _int64 end = timeInMillis();
    printf("Index build and save took %llds (%lld bases/s)\n",
           (end - start) / 1000, (_int64) nBases / max((end - start) / 1000, (_int64) 1)); 
}

SNAPHashTable** GenomeIndex::allocateHashTables(
    unsigned*       o_nTables,
    size_t          capacity,
    double          slack,
    int             seedLen,
    unsigned        hashTableKeySize,
    double*         biasTable)
{
    _ASSERT(NULL != biasTable);

    BigAllocUseHugePages = false;   // Huge pages just slow down allocation and don't help much for hash table build, so don't use them.

    if (slack <= 0) {
        fprintf(stderr, "allocateHashTables: must have positive slack for the hash table to work.  0.3 is probably OK, 0.1 is minimal, less will wreak havoc with perf.\n");
        soft_exit(1);
    }

    if (seedLen <= 0) {
        fprintf(stderr, "allocateHashTables: seedLen is too small (must be > 0, and practically should be >= 15 or so.\n");
        soft_exit(1);
    }

    if (hashTableKeySize < 4 || hashTableKeySize > 8) {
        fprintf(stderr, "allocateHashTables: key size must be 4-8 inclusive\n");
        soft_exit(1);
    }

    if ((unsigned)seedLen < hashTableKeySize * 4) {
        fprintf(stderr, "allocateHashTables: key size too large for seedLen.\n");
        soft_exit(1);
    }

    if ((unsigned)seedLen > hashTableKeySize * 4 + 9) {
        fprintf(stderr, "allocateHashTables: key size too small for seeLen.\n");
        soft_exit(1);
    }
    
    //
    // Make an array of HashTables, size depending on the seed size.  The way the index works is that we use
    // the low bits of the seed as a hash key.  Any remaining bases are used as an index into the
    // particular hash table in question.  The division between "low" and "high" depends on the hash table key size.
    //
    unsigned nHashTablesToBuild = 1 << ((seedLen - hashTableKeySize * 4) * 2);

    if (nHashTablesToBuild > 256 * 1024) {
        fprintf(stderr, "allocateHashTables: key size too small for seedLen.  Try specifying -keySize and giving it a larger value.\n");
        soft_exit(1);
    }
    //
    // Average size of the hash table.  We bias this later based on the actual content of the genome.
    //
    size_t hashTableSize = (size_t) ((double)capacity * (slack + 1.0) / nHashTablesToBuild);
    
    SNAPHashTable **hashTables = new SNAPHashTable*[nHashTablesToBuild];
    
    for (unsigned i = 0; i < nHashTablesToBuild; i++) {
        //
        // Create the actual hash tables.  It turns out that the human genome is highly non-uniform in its
        // sequences of bases, so we bias the hash table sizes based on their popularity (which is emperically
        // measured), or use the estimates that we generated and passed in as "biasTable."
        //
        double bias = biasTable[i];
        unsigned biasedSize = (unsigned) (hashTableSize * bias);
        if (biasedSize < 100) {
            biasedSize = 100;
        }
        hashTables[i] = new SNAPHashTable(biasedSize, hashTableKeySize);

        if (NULL == hashTables[i]) {
            fprintf(stderr, "IndexBuilder: unable to allocate HashTable %d of %d\n", i+1, nHashTablesToBuild);
            soft_exit(1);
        }
    }

    *o_nTables = nHashTablesToBuild;
    return hashTables;
}

    bool
GenomeIndex::BuildIndexToDirectory(const Genome *genome, int seedLen, double slack, bool computeBias, const char *directoryName, _uint64 overflowTableFactor,
                                    unsigned maxThreads, unsigned chromosomePaddingSize, bool forceExact, unsigned hashTableKeySize, const char *histogramFileName)
{
    bool buildHistogram = (histogramFileName != NULL);
    FILE *histogramFile;
    if (buildHistogram) {
        histogramFile = fopen(histogramFileName, "w");
        if (NULL == histogramFile) {
            fprintf(stderr,"Unable to open histogram file '%s', skipping it.\n", histogramFileName);
            buildHistogram = false;
        }
    }


    if (mkdir(directoryName, 0777) != 0 && errno != EEXIST) {
        fprintf(stderr,"BuildIndex: failed to create directory %s\n",directoryName);
        return false;
    }

    GenomeIndex *index = new GenomeIndex();
    index->genome = NULL;   // We always delete the index when we're done, but we delete the genome first to save space during the overflow table build.

    unsigned countOfBases = genome->getCountOfBases();
    if (countOfBases > 0xfffffff0) {
        fprintf(stderr, "Genome is too big for SNAP.  Must be some headroom beneath 2^32 bases.\n");
        return false;
    }

    // Compute bias table sizes, unless we're using the precomputed ones hardcoded in BiasTables.cpp
    double *biasTable = NULL;
    if (computeBias) {
        unsigned nHashTables = 1 << ((max((unsigned)seedLen, hashTableKeySize * 4) - hashTableKeySize * 4) * 2);
        biasTable = new double[nHashTables];
        ComputeBiasTable(genome, seedLen, biasTable, maxThreads, forceExact, hashTableKeySize);
    } else {
        biasTable = hg19_biasTables[hashTableKeySize][seedLen];
        if (NULL == biasTable) {
            fprintf(stderr, "-hg19 not available for this seed length/key size pair.\n");
            return false;
        }
    }

    printf("Allocating memory for hash and overflow tables...");
    _int64 start = timeInMillis();
    unsigned nHashTables;
    SNAPHashTable** hashTables = index->hashTables =
        allocateHashTables(&nHashTables, countOfBases, slack, seedLen, hashTableKeySize, biasTable);
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
    unsigned unusedKeyValue = InvalidGenomeLocation;
    unsigned unusedDataValue[2];
    unusedDataValue[0] = InvalidGenomeLocation;
    unusedDataValue[1] = InvalidGenomeLocation;

    // Don't create *too* many overflow entries because they consume significant memory.
    // It would be good to grow these dynamically instead of living with a fixed-size array.
    unsigned nOverflowEntries = __min((unsigned)((_uint64)countOfBases * overflowTableFactor / (_uint64) 1000), 0xfffffffd - countOfBases);
    unsigned nOverflowBackpointers = (unsigned)__min((_uint64)countOfBases, (_uint64)nOverflowEntries * (_uint64)6);

  
    OverflowEntry *overflowEntries = new OverflowEntry[nOverflowEntries];
    if (NULL == overflowEntries) {
        fprintf(stderr,"Unable to allocate oveflow entries.\n");
        soft_exit(1);
    }

    OverflowBackpointer *overflowBackpointers = new OverflowBackpointer[nOverflowBackpointers];
    if (NULL == overflowBackpointers) {
        fprintf(stderr,"Unable to allocate overflow backpointers\n");
        soft_exit(1);
    }
  
    printf("%llds\n%d nOverflowEntries, %lld bytes, %u nOverflowBackpointers, %lld bytes\nBuilding hash tables.\n", 
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
        threadContexts[i].hashTableLocks = hashTableLocks;
        threadContexts[i].hashTableKeySize = hashTableKeySize;

        StartNewThread(BuildHashTablesWorkerThreadMain, &threadContexts[i]);
    }

    WaitForSingleWaiterObject(&doneObject);
    DestroySingleWaiterObject(&doneObject);

    if ((_int64)nextOverflowIndex * 3 + (_int64)countOfDuplicateOverflows + (_int64)genome->getCountOfBases() > 0xfffffff0) {
        fprintf(stderr,"Ran out of overflow table namespace. This genome cannot be indexed with this seed size.  Try a larger one.\n");
        exit(1);
    }

    size_t totalUsedHashTableElements = 0;
    for (unsigned j = 0; j < index->nHashTables; j++) {
        totalUsedHashTableElements += hashTables[j]->GetUsedElementCount();
//        printf("HashTable[%d] has %lld used elements, loading %lld%%\n",j,(_int64)hashTables[j]->GetUsedElementCount(),
//                (_int64)hashTables[j]->GetUsedElementCount() * 100 / (_int64)hashTables[j]->GetTableSize());
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

    const unsigned maxHistogramEntry = 500000;
    unsigned countOfTooBigForHistogram = 0;
    unsigned sumOfTooBigForHistogram = 0;
    unsigned largestSeed = 0;
    unsigned totalNonSingletonSeeds = 0;
    unsigned *histogram = NULL;
    if (buildHistogram) {
        histogram = new unsigned[maxHistogramEntry+1];
        for (unsigned i = 0; i <= maxHistogramEntry; i++) {
            histogram[i] = 0;
        }
    }

    unsigned overflowTableIndex = 0;
    for (unsigned i = 0 ; i < nextOverflowIndex; i++) {
        OverflowEntry *overflowEntry = &overflowEntries[i];

        _ASSERT(overflowEntry->nInstances >= 2);

        if (timeInMillis() - lastPrintTime > 60 * 1000) {
            printf("%d/%d duplicate seeds, %d/%d backpointers processed\n",i,nextOverflowIndex,nBackpointersProcessed,nextOverflowBackpointer-1);
            lastPrintTime = timeInMillis();
        }

        //
        // If we're building a histogram, update it.
        //
        if (buildHistogram) {
            totalNonSingletonSeeds += overflowEntry->nInstances;
            if (overflowEntry->nInstances > maxHistogramEntry) {
                countOfTooBigForHistogram++;
                sumOfTooBigForHistogram += overflowEntry->nInstances;
            } else {
                histogram[overflowEntry->nInstances]++;
            }
            largestSeed = __max(largestSeed, overflowEntry->nInstances);
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

    if (buildHistogram) {
        histogram[1] = (unsigned)(totalUsedHashTableElements - totalNonSingletonSeeds);
        for (unsigned i = 0; i <= maxHistogramEntry; i++) {
            if (histogram[i] != 0) {
                fprintf(histogramFile,"%d\t%d\n", i, histogram[i]);
            }
        }
        fprintf(histogramFile, "%d larger than %d with %d total genome locations, largest seed %d\n", countOfTooBigForHistogram, maxHistogramEntry, sumOfTooBigForHistogram, largestSeed);
        fclose(histogramFile);
        delete [] histogram;
    }

    //
    // Now save out the part of the index that's independent of the genome itself.
    //
    printf("Overflow table build took %llds\nSaving genome index...", (timeInMillis() + 500 - start)/1000);
    start = timeInMillis();

    //
    // The save format is:
    //  file 'GenomeIndex' contains in order major version, minor version, nHashTables, overflowTableSize, seedLen, chromosomePaddingSize.
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

    fprintf(indexFile,"%d %d %d %d %d %d %d", GenomeIndexFormatMajorVersion, GenomeIndexFormatMinorVersion, index->nHashTables, index->overflowTableSize, seedLen, chromosomePaddingSize, hashTableKeySize);

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

    snprintf(filenameBuffer,filenameBufferSize,"%s%cGenomeIndexHash",directoryName,PATH_SEP);
    FILE *tablesFile = fopen(filenameBuffer, "wb");
    if (NULL == tablesFile) {
        fprintf(stderr, "Unable to open hash table file '%s'\n", filenameBuffer);
        soft_exit(1);
    }

    for (unsigned i = 0; i < index->nHashTables; i++) {
        if (!hashTables[i]->saveToFile(tablesFile)) {
            fprintf(stderr,"GenomeIndex::saveToDirectory: Failed to save hash table %d\n",i);
            return false;
        }
    }

    fclose(tablesFile);

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
    unsigned overflowBackpointerIndex = (unsigned)InterlockedIncrementAndReturnNewValue((volatile int *)nextOverflowBackpointer) - 1;
    if (nOverflowBackpointers <= overflowBackpointerIndex) {
        fprintf(stderr,"Ran out of overflow backpointers.  Consider using the -O switch with a larger value to increase space.\n");
        soft_exit(1);
    }
    OverflowBackpointer *newBackpointer = &overflowBackpointers[overflowBackpointerIndex];
 
    newBackpointer->nextIndex = overflowEntry->backpointerIndex;
    newBackpointer->genomeOffset = genomeOffset;
    overflowEntry->backpointerIndex = overflowBackpointerIndex;
    overflowEntry->nInstances++;
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
    unsigned majorVersion, minorVersion, chromosomePadding;
    int nRead;
    if (7 != (nRead = fscanf(indexFile,"%d %d %d %d %d %d %d", &majorVersion, &minorVersion, &index->nHashTables, &index->overflowTableSize, &seedLen, &chromosomePadding, &index->hashTableKeySize))) {
        if (3 == nRead || 6 == nRead) {
            fprintf(stderr, "Indices built by versions before 0.16.28 are no longer supported.  Please rebuild your index.\n");
        } else {
            fprintf(stderr,"GenomeIndex::LoadFromDirectory: didn't read initial values\n");
        }
        fclose(indexFile);
        delete index;
        return NULL;
    }
    fclose(indexFile);

    if (majorVersion != GenomeIndexFormatMajorVersion) {
        fprintf(stderr,"This genome index appears to be from a different version of SNAP than this, and so we can't read it.  Index version %d, SNAP index format version %d\n",
            majorVersion, GenomeIndexFormatMajorVersion);
        soft_exit(1);
    }

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

    snprintf(filenameBuffer,filenameBufferSize,"%s%cGenomeIndexHash",directoryName,PATH_SEP);
    FILE *tablesFile = fopen(filenameBuffer, "rb");
    if (NULL == tablesFile) {
        fprintf(stderr,"Unable to open genome hash table file '%s'\n", filenameBuffer);
        soft_exit(1);
    }

    for (unsigned i = 0; i < index->nHashTables; i++) {
        if (NULL == (index->hashTables[i] = SNAPHashTable::loadFromFile(tablesFile))) {
            fprintf(stderr,"GenomeIndex::loadFromDirectory: Failed to load hash table %d\n",i);
            delete index;
            return NULL;
        }
    }

    fclose(tablesFile);
    tablesFile = NULL;

    snprintf(filenameBuffer,filenameBufferSize,"%s%cGenome",directoryName,PATH_SEP);
    if (NULL == (index->genome = Genome::loadFromFile(filenameBuffer, chromosomePadding))) {
        fprintf(stderr,"GenomeIndex::loadFromDirectory: Failed to load the genome itself\n");
        delete index;
        return NULL;
    }

    if ((_int64)index->genome->getCountOfBases() + (_int64)index->overflowTableSize > 0xfffffff0) {
        fprintf(stderr,"\nThis index has too many overflow entries to be valid.  Some early versions of SNAP\n"
                        "allowed building indices with too small of a seed size, and this appears to be such\n"
                        "an index.  You can no longer build indices like this, and you also can't use them\n"
                        "because they are corrupt and would produce incorrect results.  Please use an index\n"
                        "built with a larger seed size.  For hg19, the seed size must be in the range of 19-23.\n"
                        "For other reference genomes this quantity will vary.\n");
        soft_exit(1);
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

    _ASSERT(seed.getHighBases(hashTableKeySize) < nHashTables);
    _uint64 lowBases = seed.getLowBases(hashTableKeySize);
    unsigned *entry = hashTables[seed.getHighBases(hashTableKeySize)]->Lookup(lowBases);
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
    //
    // WARNING: the code in the IntersectingPairedEndAligner relies on being able to look at 
    // hits[-1].  It doesn't care about the value, but it must not be a bogus pointer.  This
    // is true with the current layout (where it will either be the hit count, the key or
    // forward pointer in the hash table entry or some intermediate hit in the case where the
    // search is constrained by minLocation/maxLocation).  If you change this, be sure to look
    // at the code and fix it.
    //
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

        if (minLocation == 0 && maxLocation == InvalidGenomeLocation) {
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
GenomeIndex::ComputeBiasTable(const Genome* genome, int seedLen, double* table, unsigned maxThreads, bool forceExact, unsigned hashTableKeySize)
/**
 * Fill in table with the table size biases for a given genome and seed size.
 * We assume that table is already of the correct size for our seed size
 * (namely 4**(seedLen-hashTableKeySize*4)), and just fill in the values.
 *
 * If the genome is less than 2^20 bases, we count the seeds in each table exactly;
 * otherwise, we estimate them using Flajolet-Martin approximate counters.
 */
{
    _int64 start = timeInMillis();
    printf("Computing bias table.\n");

    unsigned nHashTables = ((unsigned)seedLen <= (hashTableKeySize * 4) ? 1 : 1 << (((unsigned)seedLen - hashTableKeySize * 4) * 2));
    unsigned countOfBases = genome->getCountOfBases();

    static const unsigned GENOME_SIZE_FOR_EXACT_COUNT = 1 << 20;  // Needs to be a power of 2 for hash sets

    bool computeExactly = (countOfBases < GENOME_SIZE_FOR_EXACT_COUNT) || forceExact;
    if (countOfBases >= (1 << 30) && forceExact) {
        fprintf(stderr,"You can't use -exact for genomes with >= 2^30 bases.\n");
        soft_exit(1);
    }
    FixedSizeVector<int> numExactSeeds(nHashTables, 0);
    vector<ApproximateCounter> approxCounters(nHashTables);

    _int64 validSeeds = 0;

    if (computeExactly) {
        FixedSizeSet<_int64> exactSeedsSeen(2 * (forceExact ? FirstPowerOf2GreaterThanOrEqualTo(countOfBases) : GENOME_SIZE_FOR_EXACT_COUNT));
        for (unsigned i = 0; i < (unsigned)(countOfBases - seedLen); i++) {
            if (i % 10000000 == 0) {
                printf("Bias computation: %lld / %lld\n",(_int64)i, (_int64)countOfBases);
            }
            const char *bases = genome->getSubstring(i,seedLen);
            //
            // Check it for NULL, because Genome won't return strings that cross contig boundaries.
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

            _ASSERT(seed.getHighBases(hashTableKeySize) < nHashTables);
             if (!exactSeedsSeen.contains(seed.getBases())) {
                exactSeedsSeen.add(seed.getBases());
                numExactSeeds[seed.getHighBases(hashTableKeySize)]++;
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
            contexts[i].hashTableKeySize = hashTableKeySize;
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
        _uint64 count = computeExactly ? numExactSeeds[i] : approxCounters[i].getCount();
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
    _uint64 lowBases[nSeedsPerBatch];

    bool addSeed(_uint64 seedLowBases) {
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

    const _uint64 printBatchSize = 100000000;
    for (unsigned i = context->genomeChunkStart; i < context->genomeChunkEnd; i++) {

            const char *bases = context->genome->getSubstring(i, context->seedLen);
            //
            // Check it for NULL, because Genome won't return strings that cross contig boundaries.
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

            unsigned whichHashTable = seed.getHighBases(context->hashTableKeySize);

            _ASSERT(whichHashTable < context->nHashTables);

            if (batches[whichHashTable].addSeed(seed.getLowBases(context->hashTableKeySize))) {
                PerCounterBatch *batch = &batches[whichHashTable];
                AcquireExclusiveLock(&context->approximateCounterLocks[whichHashTable]);
                batch->apply(&(*context->approxCounters)[whichHashTable]);    
                ReleaseExclusiveLock(&context->approximateCounterLocks[whichHashTable]);

                _int64 basesProcessed = InterlockedAdd64AndReturnNewValue(context->nBasesProcessed, PerCounterBatch::nSeedsPerBatch + unrecordedSkippedSeeds);

                if ((_uint64)basesProcessed / printBatchSize > ((_uint64)basesProcessed - PerCounterBatch::nSeedsPerBatch - unrecordedSkippedSeeds)/printBatchSize) {
                    printf("Bias computation: %lld / %lld\n",(basesProcessed/printBatchSize)*printBatchSize, (_int64)countOfBases);
                }
                unrecordedSkippedSeeds= 0;  // We've now recorded them.
            }
    }

    for (unsigned i = 0; i < context->nHashTables; i++) {
        _int64 basesProcessed = InterlockedAdd64AndReturnNewValue(context->nBasesProcessed, batches[i].nUsed + unrecordedSkippedSeeds);

        if ((_uint64)basesProcessed / printBatchSize > ((_uint64)basesProcessed - batches[i].nUsed - unrecordedSkippedSeeds)/printBatchSize) {
            printf("Bias computation: %lld / %lld\n",(basesProcessed/printBatchSize)*printBatchSize, (_int64)countOfBases);
        }

        unrecordedSkippedSeeds = 0; // All except the first time through the loop this will be 0.

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

    const _int64 printPeriod = 100000000;
    _uint64 unrecordedSkippedSeeds = 0;

    for (unsigned genomeLocation = context->genomeChunkStart; genomeLocation < context->genomeChunkEnd; genomeLocation++) {
        const char *bases = genome->getSubstring(genomeLocation, seedLen);
        //
        // Check it for NULL, because Genome won't return strings that cross contig boundaries.
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

       unsigned whichHashTable = seed.getHighBases(context->hashTableKeySize);
       _ASSERT(whichHashTable < index->nHashTables);
 
         if (batches[whichHashTable].addSeed(genomeLocation, seed.getLowBases(context->hashTableKeySize), usingComplement)) {
            AcquireExclusiveLock(&context->hashTableLocks[whichHashTable]);
            for (unsigned i = 0; i < batches[whichHashTable].nUsed; i++) {
                ApplyHashTableUpdate(context, whichHashTable, batches[whichHashTable].entries[i].genomeLocation, 
                    batches[whichHashTable].entries[i].lowBases, batches[whichHashTable].entries[i].usingComplement,
                    &bothComplementsUsed, &countOfDuplicateOverflows);
            }
            ReleaseExclusiveLock(&context->hashTableLocks[whichHashTable]);

            _int64 newNBasesProcessed = InterlockedAdd64AndReturnNewValue(context->nBasesProcessed, batches[whichHashTable].nUsed + unrecordedSkippedSeeds);

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
        _int64 basesProcessed = InterlockedAdd64AndReturnNewValue(context->nBasesProcessed, batches[whichHashTable].nUsed + unrecordedSkippedSeeds);

        if ((_uint64)basesProcessed / printPeriod > ((_uint64)basesProcessed - batches[whichHashTable].nUsed - unrecordedSkippedSeeds)/printPeriod) {
            printf("Indexing %lld / %lld\n",(basesProcessed/printPeriod)*printPeriod, (_int64)countOfBases);
        }

        unrecordedSkippedSeeds = 0; // All except the first time through the loop this will be 0.        
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
GenomeIndex::ApplyHashTableUpdate(BuildHashTablesThreadContext *context, _uint64 whichHashTable, unsigned genomeLocation, _uint64 lowBases, bool usingComplement,
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

        if (!hashTable->Insert(lowBases, newEntry)) {
            for (unsigned j = 0; j < index->nHashTables; j++) {
                printf("HashTable[%d] has %lld used elements\n",j,(_int64)index->hashTables[j]->GetUsedElementCount());
            }
            printf("IndexBuilder: exceeded size of hash table %d.\n"
                    "If you're indexing a non-human genome, make sure not to pass the -hg19 option.  Otheriwse, use -exact or increase slack with -h.\n",
                    whichHashTable);
            soft_exit(1);
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
            unsigned overflowIndex = (unsigned)InterlockedIncrementAndReturnNewValue((volatile int *)context->nextOverflowIndex) - 1;
            if (overflowIndex >= context->nOverflowEntries) {
                if (0xffffffff - overflowIndex < 10) {
                    fprintf(stderr,"You've run out of index overflow address space.  Perhaps a larger seed size will help.\n");
                } else {
                    fprintf(stderr,"Index builder: Overflowed overflow table.  Consider using the -O switch with a larger value to make more space.\n");
                }
                soft_exit(1);
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
