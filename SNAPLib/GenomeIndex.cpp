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
#include "GenericFile.h"
#include "Genome.h"
#include "GenomeIndex.h"
#include "HashTable.h"
#include "Seed.h"
#include "exit.h"
#include "Error.h"
#include "directions.h"
#include "GTFReader.h"

using namespace std;

static const int DEFAULT_SEED_SIZE = 20;
static const double DEFAULT_SLACK = 0.3;
static const unsigned DEFAULT_PADDING = 500;
static const unsigned DEFAULT_KEY_BYTES = 4;


static void tusage()
{
    WriteErrorMessage(
           "Usage: snapr transcriptome <input.gtf> <input.fa> <output-dir> [<options>]\n"
           "Options:\n"
            "  -s               Seed size (default: %d)\n"
            "  -h               Hash table slack (default: %.1f)\n"
            "  -hg19            Use pre-computed table bias for hg19, which results in better speed, balance, and memory footprint but may not work for other references.\n"
            "  -Ofactor         This parameter is deprecated and will be ignored.\n"
            " -tMaxThreads      Specify the maximum number of threads to use. Default is the number of cores.\n"
            " -B<chars>     Specify characters to use as chromosome name terminators in the FASTA header line; these characters and anything after are\n"
            "               not part of the chromosome name.  You must specify all characters on a single -B switch.  So, for example, with -B_|,\n"
            "               the FASTA header line '>chr1|Chromosome 1' would generate a chromosome named 'chr1'.  There's a separate flag for\n"
            "               indicating that a space is a terminator.\n"
            " -bSpace       Indicates that the space character is a terminator for chromosome names (see -B above).  This may be used in addition\n"
            "               to other terminators specified by -B.  -B and -bSpace are case sensitive.\n"
            " -pPadding         Specify the number of Ns to put as padding between chromosomes.  This must be as large as the largest\n"
            "                   edit distance you'll ever use, and there's a performance advantage to have it be bigger than any\n"
            "                   read you'll process.  Default is %d\n"
            " -HHistogramFile   Build a histogram of seed popularity.  This is just for information, it's not used by SNAP.\n"
            " -exact            Compute hash table sizes exactly.  This will slow down index build, but may be necessary in some cases\n"
            " -keysize          The number of bytes to use for the hash table key.  Larger values increase SNAP's memory footprint, but allow larger seeds.  Default: %d\n"
                        " -large            Build a larger index that's a little faster, particualrly for runs with quick/inaccurate parameters.  Increases index size by\n"
                        "                   about 30%%, depending on the other index parameters and the contents of the reference genome\n"
                        ,
            DEFAULT_SEED_SIZE,
            DEFAULT_SLACK,
            DEFAULT_PADDING,
            DEFAULT_KEY_BYTES);

    soft_exit(1);
}

static void usage()
{
    WriteErrorMessage(
            "Usage: snapr index <input.fa> <output-dir> [<options>]\n"
            "Options:\n"
            "  -s               Seed size (default: %d)\n"
            "  -h               Hash table slack (default: %.1f)\n"
            "  -hg19            Use pre-computed table bias for hg19, which results in better speed, balance, and memory footprint but may not work for other references.\n"
            "  -Ofactor         This parameter is deprecated and will be ignored.\n"
            " -tMaxThreads      Specify the maximum number of threads to use. Default is the number of cores.\n"
            " -B<chars>     Specify characters to use as chromosome name terminators in the FASTA header line; these characters and anything after are\n"
            "               not part of the chromosome name.  You must specify all characters on a single -B switch.  So, for example, with -B_|,\n"
            "               the FASTA header line '>chr1|Chromosome 1' would generate a chromosome named 'chr1'.  There's a separate flag for\n"
            "               indicating that a space is a terminator.\n"
            " -bSpace       Indicates that the space character is a terminator for chromosome names (see -B above).  This may be used in addition\n"
            "               to other terminators specified by -B.  -B and -bSpace are case sensitive.\n"
            " -pPadding         Specify the number of Ns to put as padding between chromosomes.  This must be as large as the largest\n"
            "                   edit distance you'll ever use, and there's a performance advantage to have it be bigger than any\n"
            "                   read you'll process.  Default is %d\n"
            " -HHistogramFile   Build a histogram of seed popularity.  This is just for information, it's not used by SNAP.\n"
            " -exact            Compute hash table sizes exactly.  This will slow down index build, but may be necessary in some cases\n"
            " -keysize          The number of bytes to use for the hash table key.  Larger values increase SNAP's memory footprint, but allow larger seeds.  Default: %d\n"
			" -large            Build a larger index that's a little faster, particualrly for runs with quick/inaccurate parameters.  Increases index size by\n"
			"                   about 30%%, depending on the other index parameters and the contents of the reference genome\n"
			,
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
    const char *pieceNameTerminatorCharacters = NULL;
    bool spaceIsAPieceNameTerminator = false;
    const char *histogramFileName = NULL;
    unsigned chromosomePadding = DEFAULT_PADDING;
    bool forceExact = false;
    unsigned keySizeInBytes = DEFAULT_KEY_BYTES;
        bool large = false;

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
        } else if (strcmp(argv[n], "-large") == 0) {
            large = true;
        } else if (argv[n][0] == '-' && argv[n][1] == 'H') {
            histogramFileName = argv[n] + 2;
        } else if (argv[n][0] == '-' && argv[n][1] == 'O') {
            // Deprecated, ignored parameter
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
        } else if (argv[n][0] == '-' && argv[n][1] == 'B') {
            pieceNameTerminatorCharacters = argv[n] + 2;
        } else if (!strcmp(argv[n], "-bSpace")) {
            spaceIsAPieceNameTerminator = true;
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

    WriteStatusMessage("Hash table slack %lf\nLoading FASTA file '%s' into memory...", slack, fastaFile);
    BigAllocUseHugePages = false;

    _int64 start = timeInMillis();
    const Genome *genome = ReadFASTAGenome(fastaFile, pieceNameTerminatorCharacters, spaceIsAPieceNameTerminator, chromosomePadding);
    if (NULL == genome) {
        WriteErrorMessage("Unable to read FASTA file\n");        soft_exit(1);
    }
    WriteStatusMessage("%llds\n", (timeInMillis() + 500 - start) / 1000);
           
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
    const char *pieceNameTerminatorCharacters = NULL;
    bool spaceIsAPieceNameTerminator = false;
    const char *histogramFileName = NULL;
    unsigned chromosomePadding = DEFAULT_PADDING;
    bool forceExact = false;
    unsigned keySizeInBytes = DEFAULT_KEY_BYTES;
	bool large = false;

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
        } else if (strcmp(argv[n], "-large") == 0) {
            large = true;
        } else if (argv[n][0] == '-' && argv[n][1] == 'H') {
            histogramFileName = argv[n] + 2;
        } else if (argv[n][0] == '-' && argv[n][1] == 'O') {
			// Deprecated, ignored parameter
        } else if (argv[n][0] == '-' && argv[n][1] == 't') {
            maxThreads = atoi(argv[n]+2);
            if (maxThreads < 1 || maxThreads > 100) {
                WriteErrorMessage("maxThreads must be between 1 and 100 inclusive (and you need not to leave a space after '-t')\n");
                soft_exit(1);
            }
        } else if (argv[n][0] == '-' && argv[n][1] == 'p') {
            chromosomePadding = atoi(argv[n]+2);
            if (0 == chromosomePadding) {
                WriteErrorMessage("Invalid chromosome padding specified, must be at least one (and in practice as large as any max edit distance you might use).\n");
                soft_exit(1);
            }
        } else if (strcmp(argv[n], "-keysize") == 0) {
            if (n + 1 < argc) {
                keySizeInBytes = atoi(argv[n+1]);
                if (keySizeInBytes < 4 || keySizeInBytes > 8) {
                    WriteErrorMessage("Key size must be between 4 and 8 inclusive\n");
                    soft_exit(1);
                }
                n++;
            } else {
                usage();
            }
        } else if (argv[n][0] == '-' && argv[n][1] == 'B') {
            pieceNameTerminatorCharacters = argv[n] + 2;
        } else if (!strcmp(argv[n], "-bSpace")) {
            spaceIsAPieceNameTerminator = true;
        } else {
            WriteErrorMessage("Invalid argument: %s\n\n", argv[n]);
            usage();
        }
    }

    if (seedLen < 16 || seedLen > 32) {
        // Seeds are stored in 64 bits, so they can't be larger than 32 bases for now.
        WriteErrorMessage("Seed length must be between 16 and 32, inclusive\n");
        soft_exit(1);
    }

    if (seedLen < 19 && !computeBias) {
        WriteErrorMessage("For hg19, you must use seed sizes between 19 and 32, or not specify -hg19 and let SNAP compute the hash table bias.\n");
        soft_exit(1);
    }

	if (seedLen * 2 < keySizeInBytes * 8) {
		WriteErrorMessage("You must specify a smaller keysize or a larger seed size.  The seed must be big enough to fill the key\n"
			"and takes two bits per base of seed.\n");
		soft_exit(1);
	}

	if (seedLen * 2 - keySizeInBytes * 8 > 16) {
		WriteErrorMessage("You must specify a biger keysize or smaller seed len.  SNAP restricts the number of hash tables to 4^8,\n"
			"and needs 4^{excess seed len} hash tables, where excess seed len is the seed size minus the four times the key size.\n");
		soft_exit(1);
	}


    WriteStatusMessage("Hash table slack %lf\nLoading FASTA file '%s' into memory...", slack, fastaFile);

    BigAllocUseHugePages = false;

    _int64 start = timeInMillis();
    const Genome *genome = ReadFASTAGenome(fastaFile, pieceNameTerminatorCharacters, spaceIsAPieceNameTerminator, chromosomePadding);
    if (NULL == genome) {
        WriteErrorMessage("Unable to read FASTA file\n");
        soft_exit(1);
    }
    WriteStatusMessage("%llds\n", (timeInMillis() + 500 - start) / 1000);
    unsigned nBases = genome->getCountOfBases();
    if (!GenomeIndex::BuildIndexToDirectory(genome, seedLen, slack, computeBias, outputDir, maxThreads, chromosomePadding, forceExact, keySizeInBytes, 
		large, histogramFileName)) {
        WriteErrorMessage("Genome index build failed\n");
        soft_exit(1);
    }
    _int64 end = timeInMillis();
    WriteStatusMessage("Index build and save took %llds (%lld bases/s)\n",
           (end - start) / 1000, (_int64) nBases / max((end - start) / 1000, (_int64) 1)); 
}

    bool
GenomeIndex::BuildIndexToDirectory(const Genome *genome, int seedLen, double slack, bool computeBias, const char *directoryName,
                                    unsigned maxThreads, unsigned chromosomePaddingSize, bool forceExact, unsigned hashTableKeySize, 
									bool large, const char *histogramFileName)
{
    bool buildHistogram = (histogramFileName != NULL);
    FILE *histogramFile;
    if (buildHistogram) {
        histogramFile = fopen(histogramFileName, "w");
        if (NULL == histogramFile) {
            WriteErrorMessage("Unable to open histogram file '%s', skipping it.\n", histogramFileName);
            buildHistogram = false;
        }
    }


    if (mkdir(directoryName, 0777) != 0 && errno != EEXIST) {
        WriteErrorMessage("BuildIndex: failed to create directory %s\n",directoryName);
        return false;
    }

    const unsigned filenameBufferSize = MAX_PATH+1;
    char filenameBuffer[filenameBufferSize];
    
	fprintf(stderr,"Saving genome...");
	_int64 start = timeInMillis();
    snprintf(filenameBuffer,filenameBufferSize,"%s%cGenome",directoryName,PATH_SEP);
    if (!genome->saveToFile(filenameBuffer)) {
        WriteErrorMessage("GenomeIndex::saveToDirectory: Failed to save the genome itself\n");
        return false;
    }
	fprintf(stderr,"%llds\n", (timeInMillis() + 500 - start) / 1000);

	GenomeIndex *index = new GenomeIndexLarge();
    index->genome = NULL;   // We always delete the index when we're done, but we delete the genome first to save space during the overflow table build.

    unsigned countOfBases = genome->getCountOfBases();
    if (countOfBases > 0xfffffff0) {
        WriteErrorMessage("Genome is too big for SNAP.  Must be some headroom beneath 2^32 bases.\n");
        return false;
    }

    // Compute bias table sizes, unless we're using the precomputed ones hardcoded in BiasTables.cpp
    double *biasTable = NULL;
    if (computeBias) {
        unsigned nHashTables = 1 << ((max((unsigned)seedLen, hashTableKeySize * 4) - hashTableKeySize * 4) * 2);
        biasTable = new double[nHashTables];
        ComputeBiasTable(genome, seedLen, biasTable, maxThreads, forceExact, hashTableKeySize, large);
    } else {
        if (large) {
            biasTable = hg19_biasTables_large[hashTableKeySize][seedLen];
        } else {
            biasTable = hg19_biasTables[hashTableKeySize][seedLen];
        }

        if (NULL == biasTable) {
            WriteErrorMessage("-hg19 not available for this seed length/key size/small-or-large combo.\n");
            return false;
        }
    }

    WriteStatusMessage("Allocating memory for hash tables...");
    start = timeInMillis();
    unsigned nHashTables;
    SNAPHashTable** hashTables = index->hashTables =
        allocateHashTables(&nHashTables, countOfBases, slack, seedLen, hashTableKeySize, large, biasTable);
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

	OverflowBackpointerAnchor *overflowAnchor = new OverflowBackpointerAnchor(__min(0xfffffff0 - countOfBases, countOfBases));   // i.e., as much as the address space will allow.
   
    WriteStatusMessage("%llds\nBuilding hash tables.\n", (timeInMillis() + 500 - start) / 1000);
  
    start = timeInMillis();
    volatile unsigned nextOverflowBackpointer = 0;

    volatile _int64 nonSeeds = 0;
    volatile _int64 seedsWithMultipleOccurrences = 0;
    volatile _int64 genomeLocationsInOverflowTable = 0;     // Number of extra hits on duplicate indices.  This should come out once we implement the overflow table.
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
        threadContexts[i].seedsWithMultipleOccurrences = &seedsWithMultipleOccurrences;
        threadContexts[i].genomeLocationsInOverflowTable = &genomeLocationsInOverflowTable;
        threadContexts[i].bothComplementsUsed = &bothComplementsUsed;
		threadContexts[i].overflowAnchor = overflowAnchor;
        threadContexts[i].nextOverflowBackpointer = &nextOverflowBackpointer;
        threadContexts[i].hashTableLocks = hashTableLocks;
        threadContexts[i].hashTableKeySize = hashTableKeySize;
		threadContexts[i].large = large;

        StartNewThread(BuildHashTablesWorkerThreadMain, &threadContexts[i]);
    }

    WaitForSingleWaiterObject(&doneObject);
    DestroySingleWaiterObject(&doneObject);

    if (seedsWithMultipleOccurrences + genomeLocationsInOverflowTable + (_int64)genome->getCountOfBases() > 0xfffffff0) {
        WriteErrorMessage("Ran out of overflow table namespace. This genome cannot be indexed with this seed size.  Try a larger one.\n");
        exit(1);
    }

    size_t totalUsedHashTableElements = 0;
    for (unsigned j = 0; j < index->nHashTables; j++) {
        totalUsedHashTableElements += hashTables[j]->GetUsedElementCount();
//        printf("HashTable[%d] has %lld used elements, loading %lld%%\n",j,(_int64)hashTables[j]->GetUsedElementCount(),
//                (_int64)hashTables[j]->GetUsedElementCount() * 100 / (_int64)hashTables[j]->GetTableSize());
    }

    WriteStatusMessage("%lld(%lld%%) seeds occur more than once, total of %lld(%lld%%) genome locations are not unique, %d(%lld%%) bad seeds, %lld both complements used %lld no string\n",
        seedsWithMultipleOccurrences,
        ((_int64)seedsWithMultipleOccurrences)*100 / countOfBases,
        genomeLocationsInOverflowTable,
        genomeLocationsInOverflowTable * 100 / countOfBases,
        (int) nonSeeds,
        ((_int64)nonSeeds *100) / countOfBases,
        bothComplementsUsed,
        noBaseAvailable);

    WriteStatusMessage("Hash table build took %llds\n",(timeInMillis() + 500 - start) / 1000);

    //
    // We're done with the raw genome.  Delete it to save some memory.
    //
  
    delete genome;
    genome = NULL;

    WriteStatusMessage("Building overflow table.\n");
    start = timeInMillis();
    fflush(stdout);

    //
    // Now build the real overflow table and simultaneously fixup the hash table entries.
    // Its format is one unsigned of the number of genome locations matching the
    // particular seed, followed by that many genome offsets.  An overflow table
    // entry is a count followed by a reverse sorted list of genome locations.
    // For each seed with multiple occurrences in the genome, there is one count.
    // For each genome location that's not unique, there is one list entry.  So, the size
    // of the overflow table is the number of non-unique seeds plus the number of non-unique
    // genome locations.
    //
    index->overflowTableSize = seedsWithMultipleOccurrences  + genomeLocationsInOverflowTable;
    index->overflowTable = (unsigned *)BigAlloc(index->overflowTableSize * sizeof(*index->overflowTable));

 	if ((_uint64)index->overflowTableSize + (_uint64)countOfBases >= 0xfffffff0) {
		WriteErrorMessage("Not enough address space to index this genome with this seed size.  Try a larger seed size.\n");
		soft_exit(1);
	}

    unsigned nBackpointersProcessed = 0;
    _int64 lastPrintTime = timeInMillis();

    const unsigned maxHistogramEntry = 500000;
    unsigned countOfTooBigForHistogram = 0;
    unsigned sumOfTooBigForHistogram = 0;
    unsigned largestSeed = 0;
    unsigned *histogram = NULL;
    if (buildHistogram) {
        histogram = new unsigned[maxHistogramEntry+1];
        for (unsigned i = 0; i <= maxHistogramEntry; i++) {
            histogram[i] = 0;
        }
    }

	//
	// Build the overflow table by walking each of the hash tables and looking for elements to fix up.
	// Write the hash tables as we go so that we can free their memory on the fly.
	//
	snprintf(filenameBuffer,filenameBufferSize,"%s%cGenomeIndexHash",directoryName,PATH_SEP);
    FILE *tablesFile = fopen(filenameBuffer, "wb");
    if (NULL == tablesFile) {
        WriteErrorMessage("Unable to open hash table file '%s'\n", filenameBuffer);
        soft_exit(1);
    }

    size_t totalBytesWritten = 0;
    unsigned overflowTableIndex = 0;
	_uint64 duplicateSeedsProcessed = 0;

	for (unsigned whichHashTable = 0; whichHashTable < nHashTables; whichHashTable++) {
		for (unsigned whichEntry = 0; whichEntry < hashTables[whichHashTable]->GetTableSize(); whichEntry++) {
			unsigned *values = (unsigned *)hashTables[whichHashTable]->getEntryValues(whichEntry);
			for (int i = 0; i < (large ? NUM_DIRECTIONS : 1); i++) {
				if (values[i] >= countOfBases && values[i] != 0xffffffff && values[i] != 0xfffffffe) {
					//
					// This is an overflow pointer.  Fix it up.  Count the number of occurrences of this
					// seed by walking the overflow chain.
					//
					duplicateSeedsProcessed++;

					unsigned nOccurrences = 0;
					unsigned backpointerIndex = values[i] - countOfBases;
					while (backpointerIndex != 0xffffffff) {
						nOccurrences++;
						OverflowBackpointer *backpointer = overflowAnchor->getBackpointer(backpointerIndex);
						_ASSERT(overflowTableIndex + nOccurrences < index->overflowTableSize);
						index->overflowTable[overflowTableIndex + nOccurrences] = backpointer->genomeOffset;
						backpointerIndex = backpointer->nextIndex;
					}

					_ASSERT(nOccurrences > 1);

					//
					// Fill the count in as the first thing in the overflow table.
					//
                    _ASSERT(overflowTableIndex < index->overflowTableSize);
					index->overflowTable[overflowTableIndex] = nOccurrences;

					//
					// And patch the value into the hash table.
					//
					values[i] = overflowTableIndex + countOfBases;
					overflowTableIndex += 1 + nOccurrences;
                    _ASSERT(overflowTableIndex <= index->overflowTableSize);
					nBackpointersProcessed += nOccurrences;

					//
					// Sort the overflow table entries, because the paired-end aligner relies on this.  Sort them backwards, because that's
					// what it expects.  For those who are desparately curious, this is because it was originally built this way by accident
					// before there was any concept of doing binary search over a seed's hits.  When the binary search was built, it relied
					// on this.  Then, when the index build was parallelized it was easier just to preserve the old order than to change the
					// code in the aligner.  So now you know.
					//
					qsort(&index->overflowTable[overflowTableIndex -nOccurrences], nOccurrences, sizeof(index->overflowTable[0]), BackwardsUnsignedCompare);

					if (timeInMillis() - lastPrintTime > 60 * 1000) {
						WriteStatusMessage("%lld/%lld duplicate seeds, %d/%d backpointers, %d/%d hash tables processed\n", 
							duplicateSeedsProcessed, seedsWithMultipleOccurrences, nBackpointersProcessed, genomeLocationsInOverflowTable,
							whichHashTable, nHashTables);
						lastPrintTime = timeInMillis();
					}

					//
					// If we're building a histogram, update it.
					//
					if (buildHistogram) {
						if (nOccurrences > maxHistogramEntry) {
							countOfTooBigForHistogram++;
							sumOfTooBigForHistogram += nOccurrences;
						} else {
							histogram[nOccurrences]++;
						}
						largestSeed = __max(largestSeed, nOccurrences);
					}

				} // If this entry needs patching
			} // forward and RC if large table
		} // for each entry in the hash table

 		//
		// We're done with this hash table, free it to releive memory pressure.
		//
		size_t bytesWrittenThisHashTable;
        if (!hashTables[whichHashTable]->saveToFile(tablesFile, &bytesWrittenThisHashTable)) {
            WriteErrorMessage("GenomeIndex::saveToDirectory: Failed to save hash table %d\n", whichHashTable);
            return false;
        }
        totalBytesWritten += bytesWrittenThisHashTable;

		delete hashTables[whichHashTable];
		hashTables[whichHashTable] = NULL;
	} // for each hash table

    fclose(tablesFile);

    _ASSERT(overflowTableIndex == index->overflowTableSize);    // We used exactly what we expected to use.

    delete overflowAnchor;
    overflowAnchor = NULL;

    if (buildHistogram) {
        histogram[1] = (unsigned)(totalUsedHashTableElements - seedsWithMultipleOccurrences);
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
    WriteStatusMessage("Overflow table build and hash table save took %llds\nSaving overflow table...", (timeInMillis() + 500 - start)/1000);
    start = timeInMillis();


    snprintf(filenameBuffer,filenameBufferSize,"%s%cOverflowTable",directoryName,PATH_SEP);
    FILE* fOverflowTable = fopen(filenameBuffer, "wb");
    if (fOverflowTable == NULL) {
        WriteErrorMessage("Unable to open overflow table file, '%s', %d\n",filenameBuffer,errno);
        return false;
    }

    const unsigned writeSize = 32 * 1024 * 1024;
    for (size_t writeOffset = 0; writeOffset < index->overflowTableSize * sizeof(*index->overflowTable); ) {
        unsigned amountToWrite = (unsigned)__min((size_t)writeSize,(size_t)index->overflowTableSize * sizeof(*index->overflowTable) - writeOffset);
 
        size_t amountWritten = fwrite(((char *)index->overflowTable) + writeOffset, 1, amountToWrite, fOverflowTable);
        if (amountWritten < amountToWrite) {
            WriteErrorMessage("GenomeIndex::saveToDirectory: fwrite failed, %d\n",errno);
            fclose(fOverflowTable);
            return false;
        }
        writeOffset += amountWritten;
    }
    fclose(fOverflowTable);
    fOverflowTable = NULL;

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
        WriteErrorMessage("Unable to open file '%s' for write.\n",filenameBuffer);
        return false;
    }

    fprintf(indexFile,"%d %d %d %d %d %d %d %lld %d", GenomeIndexFormatMajorVersion, GenomeIndexFormatMinorVersion, index->nHashTables, 
        index->overflowTableSize, seedLen, chromosomePaddingSize, hashTableKeySize, totalBytesWritten, large ? 0 : 1); 

    fclose(indexFile);
 
    delete index;
    if (computeBias && biasTable != NULL) {
        delete[] biasTable;
    }
 
    WriteStatusMessage("%llds\n", (timeInMillis() + 500 - start) / 1000);
    
    return true;
}

SNAPHashTable** GenomeIndex::allocateHashTables(
    unsigned*       o_nTables,
    size_t          countOfBases,
    double          slack,
    int             seedLen,
    unsigned        hashTableKeySize,
	bool			large,
    double*         biasTable)
{
    _ASSERT(NULL != biasTable);

    BigAllocUseHugePages = false;   // Huge pages just slow down allocation and don't help much for hash table build, so don't use them.

    if (slack <= 0) {
        WriteErrorMessage("allocateHashTables: must have positive slack for the hash table to work.  0.3 is probably OK, 0.1 is minimal, less will wreak havoc with perf.\n");
        soft_exit(1);
    }

    if (seedLen <= 0) {
        WriteErrorMessage("allocateHashTables: seedLen is too small (must be > 0, and practically should be >= 15 or so.\n");
        soft_exit(1);
    }

    if (hashTableKeySize < 4 || hashTableKeySize > 8) {
        WriteErrorMessage("allocateHashTables: key size must be 4-8 inclusive\n");
        soft_exit(1);
    }

    if ((unsigned)seedLen < hashTableKeySize * 4) {
        WriteErrorMessage("allocateHashTables: key size too large for seedLen.\n");
        soft_exit(1);
    }

    if ((unsigned)seedLen > hashTableKeySize * 4 + 9) {
        WriteErrorMessage("allocateHashTables: key size too small for seeLen.\n");
        soft_exit(1);
    }
    
    //
    // Make an array of HashTables, size depending on the seed size.  The way the index works is that we use
    // the low bits of the seed as a hash key.  Any remaining bases are used as an index into the
    // particular hash table in question.  The division between "low" and "high" depends on the hash table key size.
    //
    unsigned nHashTablesToBuild = 1 << ((seedLen - hashTableKeySize * 4) * 2);

    if (nHashTablesToBuild > 256 * 1024) {
        WriteErrorMessage("allocateHashTables: key size too small for seedLen.  Try specifying -keySize and giving it a larger value.\n");
        soft_exit(1);
    }
    //
    // Average size of the hash table.  We bias this later based on the actual content of the genome.
    //
    size_t hashTableSize = (size_t) ((double)countOfBases * (slack + 1.0) / nHashTablesToBuild);
    
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
        hashTables[i] = new SNAPHashTable(biasedSize, hashTableKeySize, sizeof(unsigned), large ? 2 : 1, InvalidGenomeLocation);
 
        if (NULL == hashTables[i]) {
            WriteErrorMessage("IndexBuilder: unable to allocate HashTable %d of %d\n", i+1, nHashTablesToBuild);
            soft_exit(1);
        }
    }

    *o_nTables = nHashTablesToBuild;
    return hashTables;
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

    unsigned
GenomeIndex::AddOverflowBackpointer(
    unsigned                     previousOverflowBackpointer,
    OverflowBackpointerAnchor   *overflowAnchor,
    volatile unsigned           *nextOverflowBackpointer,
    unsigned                     genomeOffset)
{
    unsigned overflowBackpointerIndex = (unsigned)InterlockedIncrementAndReturnNewValue((volatile int *)nextOverflowBackpointer) - 1;
    OverflowBackpointer *newBackpointer = overflowAnchor->getBackpointer(overflowBackpointerIndex);
 
    newBackpointer->nextIndex = previousOverflowBackpointer;
    newBackpointer->genomeOffset = genomeOffset;

    return overflowBackpointerIndex;
}


GenomeIndex::GenomeIndex() : nHashTables(0), hashTables(NULL), overflowTable(NULL), genome(NULL)
{
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

    if (NULL == tablesBlob) {
        BigDealloc(tablesBlob);
        tablesBlob = NULL;
    }
}

    void
GenomeIndex::ComputeBiasTable(const Genome* genome, int seedLen, double* table, unsigned maxThreads, bool forceExact, unsigned hashTableKeySize, bool large)
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
    WriteStatusMessage("Computing bias table.\n");

    unsigned nHashTables = ((unsigned)seedLen <= (hashTableKeySize * 4) ? 1 : 1 << (((unsigned)seedLen - hashTableKeySize * 4) * 2));
    unsigned countOfBases = genome->getCountOfBases();

    static const unsigned GENOME_SIZE_FOR_EXACT_COUNT = 1 << 20;  // Needs to be a power of 2 for hash sets

    bool computeExactly = (countOfBases < GENOME_SIZE_FOR_EXACT_COUNT) || forceExact;
    if (countOfBases >= (1 << 30) && forceExact) {
        WriteErrorMessage("You can't use -exact for genomes with >= 2^30 bases.\n");
        soft_exit(1);
    }
    FixedSizeVector<int> numExactSeeds(nHashTables, 0);
    vector<ApproximateCounter> approxCounters(nHashTables);

    _int64 validSeeds = 0;

    if (computeExactly) {
        FixedSizeSet<_int64> exactSeedsSeen(2 * (forceExact ? FirstPowerOf2GreaterThanOrEqualTo(countOfBases) : GENOME_SIZE_FOR_EXACT_COUNT));
        for (unsigned i = 0; i < (unsigned)(countOfBases - seedLen); i++) {
            if (i % 10000000 == 0) {
                WriteStatusMessage("Bias computation: %lld / %lld\n",(_int64)i, (_int64)countOfBases);
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

			if (large && seed.isBiggerThanItsReverseComplement()) {
				// For large hash tables, because seeds and their reverse complements are stored
				// together, figure out which one is used for the hash table key, and use that
				// one.
				seed = ~seed;
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
			contexts[i].large = large;

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

    //printf("Bias table:\n");
    //for (unsigned i = 0; i < nHashTables; i++) {
    //    printf("%u -> %lf\n", i, table[i]);
    //}

    WriteStatusMessage("Computed bias table in %llds\n", (timeInMillis() + 500 - start) / 1000);
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
	bool large = context->large;

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

			if (large && seed.isBiggerThanItsReverseComplement()) {
				//
				// Figure out if we're using this base or its complement.
				//
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
					WriteStatusMessage("Bias computation: %lld / %lld\n",(basesProcessed/printBatchSize)*printBatchSize, (_int64)countOfBases);
				}
				unrecordedSkippedSeeds= 0;  // We've now recorded them.
			}
    }

    for (unsigned i = 0; i < context->nHashTables; i++) {
        _int64 basesProcessed = InterlockedAdd64AndReturnNewValue(context->nBasesProcessed, batches[i].nUsed + unrecordedSkippedSeeds);

        if ((_uint64)basesProcessed / printBatchSize > ((_uint64)basesProcessed - batches[i].nUsed - unrecordedSkippedSeeds)/printBatchSize) {
            WriteStatusMessage("Bias computation: %lld / %lld\n",(basesProcessed/printBatchSize)*printBatchSize, (_int64)countOfBases);
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



    void 
GenomeIndex::BuildHashTablesWorkerThreadMain(void *param)
{
    BuildHashTablesThreadContext *context = (BuildHashTablesThreadContext *)param;

    context->index->BuildHashTablesWorkerThread(context);
}
    
    void
GenomeIndex::BuildHashTablesWorkerThread(BuildHashTablesThreadContext *context)
{
    unsigned countOfBases = context->genome->getCountOfBases();
    const Genome *genome = context->genome;
    unsigned seedLen = context->seedLen;
	bool large = context->large;
 
    //
    // Batch the insertions into the hash tables, because otherwise we spend all of
    // our time acquiring and releasing locks.
    //
 
    PerHashTableBatch *batches = new PerHashTableBatch[nHashTables];
    IndexBuildStats stats;

    for (unsigned genomeLocation = context->genomeChunkStart; genomeLocation < context->genomeChunkEnd; genomeLocation++) {
        const char *bases = genome->getSubstring(genomeLocation, seedLen);
        //
        // Check it for NULL, because Genome won't return strings that cross contig boundaries.
        //
        if (NULL == bases) {
            stats.noBaseAvailable++;
            stats.unrecordedSkippedSeeds++;
            continue;
        }

        //
        // We don't build seeds out of sections of the genome that contain 'N.'  If this is one, skip it.
        //
        if (!Seed::DoesTextRepresentASeed(bases, seedLen)) {
            stats.nonSeeds++;
            stats.unrecordedSkippedSeeds++;
            continue;
        }

        Seed seed(bases, seedLen);

        indexSeed(genomeLocation, seed, batches, context, &stats, large);
    } // For each genome base in our area

    //
    // Now apply the updates from the batches that were left over
    //

    completeIndexing(batches, context, &stats, large);

    InterlockedAdd64AndReturnNewValue(context->noBaseAvailable, stats.noBaseAvailable);
    InterlockedAdd64AndReturnNewValue(context->nonSeeds, stats.nonSeeds);
    InterlockedAdd64AndReturnNewValue(context->bothComplementsUsed, stats.bothComplementsUsed);
    InterlockedAdd64AndReturnNewValue(context->genomeLocationsInOverflowTable, stats.genomeLocationsInOverflowTable);
    InterlockedAdd64AndReturnNewValue(context->seedsWithMultipleOccurrences, stats.seedsWithMultipleOccurrences);

    delete [] batches;

    if (0 == InterlockedDecrementAndReturnNewValue(context->runningThreadCount)) {
        SignalSingleWaiterObject(context->doneObject);
    }

}
    
const _int64 GenomeIndex::printPeriod = 100000000;

    void 
GenomeIndex::ApplyHashTableUpdate(BuildHashTablesThreadContext *context, _uint64 whichHashTable, unsigned genomeLocation, _uint64 lowBases, bool usingComplement,
                _int64 *bothComplementsUsed, _int64 *genomeLocationsInOverflowTable, _int64 *seedsWithMultipleOccurrences, bool large)
{
	_ASSERT(large || !usingComplement);
    GenomeIndex *index = context->index;
    unsigned countOfBases = context->genome->getCountOfBases();
    SNAPHashTable *hashTable = index->hashTables[whichHashTable];
    _ASSERT(hashTable->GetValueSizeInBytes() == 4);
    unsigned *entry = (unsigned *)hashTable->SlowLookup(lowBases);  // use SlowLookup because we might have overflowed the table.  Cast is OK because valueSize == 4
    if (NULL == entry) {
        SNAPHashTable::ValueType newEntry[2];	// We only use [0] if !large, but it doesn't hurt to declare two
		if (large) {
			//
			// We haven't yet seen either this seed or its complement.  Make a new hash table
			// entry.
			//
			if (!usingComplement) {
				newEntry[0] = genomeLocation;
				newEntry[1] = 0xfffffffe; // Use 0xfffffffe for unused, because we gave 0xffffffff to the hash table package.
			} else{
				newEntry[0] = 0xfffffffe; // Use 0xfffffffe for unused, because we gave 0xffffffff to the hash table package.
				newEntry[1] = genomeLocation;
			}
		} else {
			newEntry[0] = genomeLocation;
		}

        _ASSERT(0 != genomeLocation);

        if (!hashTable->Insert(lowBases, newEntry)) {
            for (unsigned j = 0; j < index->nHashTables; j++) {
                WriteErrorMessage("HashTable[%d] has %lld used elements\n",j,(_int64)index->hashTables[j]->GetUsedElementCount());
            }
            WriteErrorMessage("IndexBuilder: exceeded size of hash table %d.\n"
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
        if (large && 0xfffffffe == entry[entryIndex]) {
            entry[entryIndex] = genomeLocation;
            (*bothComplementsUsed)++;
        } else if (entry[entryIndex] < countOfBases) {

            (*seedsWithMultipleOccurrences)++;
            (*genomeLocationsInOverflowTable) += 2;    

            unsigned overflowIndex = AddOverflowBackpointer(0xffffffff, context->overflowAnchor, context->nextOverflowBackpointer, entry[entryIndex]);
            overflowIndex = AddOverflowBackpointer(overflowIndex, context->overflowAnchor, context->nextOverflowBackpointer, genomeLocation);

            entry[entryIndex] = overflowIndex + countOfBases;
        } else {
            //
            // Stick another entry in the existing overflow bucket.
            //

            unsigned overflowIndex = AddOverflowBackpointer(entry[entryIndex] - countOfBases, context->overflowAnchor, context->nextOverflowBackpointer, genomeLocation);
            entry[entryIndex] = overflowIndex + countOfBases;

            (*genomeLocationsInOverflowTable)++;
        } // If the existing entry had the complement empty, needed a new overflow entry or extended an old one
    } // If new or existing entry.
}

GenomeIndexLarge::~GenomeIndexLarge() {}
GenomeIndexSmall::~GenomeIndexSmall() {}

    void
GenomeIndex::indexSeed(unsigned genomeLocation, Seed seed, PerHashTableBatch *batches, BuildHashTablesThreadContext *context, IndexBuildStats *stats, bool large)
{
	bool usingComplement = large && seed.isBiggerThanItsReverseComplement();
	if (usingComplement) {
		seed = ~seed;       // Couldn't resist using ~ for this.
	}

    unsigned whichHashTable = seed.getHighBases(context->hashTableKeySize);
    _ASSERT(whichHashTable < nHashTables);
 
	if (batches[whichHashTable].addSeed(genomeLocation, seed.getLowBases(context->hashTableKeySize), usingComplement)) {
		AcquireExclusiveLock(&context->hashTableLocks[whichHashTable]);
		for (unsigned i = 0; i < batches[whichHashTable].nUsed; i++) {
			ApplyHashTableUpdate(context, whichHashTable, batches[whichHashTable].entries[i].genomeLocation, 
				batches[whichHashTable].entries[i].lowBases, batches[whichHashTable].entries[i].usingComplement,
				&stats->bothComplementsUsed, &stats->genomeLocationsInOverflowTable, &stats->seedsWithMultipleOccurrences, large);
		}
		ReleaseExclusiveLock(&context->hashTableLocks[whichHashTable]);

		_int64 newNBasesProcessed = InterlockedAdd64AndReturnNewValue(context->nBasesProcessed, batches[whichHashTable].nUsed + stats->unrecordedSkippedSeeds);

		if ((unsigned)(newNBasesProcessed / printPeriod) > (unsigned)((newNBasesProcessed - batches[whichHashTable].nUsed - stats->unrecordedSkippedSeeds) / printPeriod)) {
			WriteStatusMessage("Indexing %lld / %lld\n", (newNBasesProcessed / printPeriod) * printPeriod, context->genome->getCountOfBases());
		}
		stats->unrecordedSkippedSeeds = 0;
		batches[whichHashTable].clear();
	} // If we filled a batch
}

    void 
GenomeIndex::completeIndexing(PerHashTableBatch *batches, BuildHashTablesThreadContext *context, IndexBuildStats *stats, bool large)
{
    for (unsigned whichHashTable = 0; whichHashTable < nHashTables; whichHashTable++) {
        _int64 basesProcessed = InterlockedAdd64AndReturnNewValue(context->nBasesProcessed, batches[whichHashTable].nUsed + stats->unrecordedSkippedSeeds);

        if ((_uint64)basesProcessed / printPeriod > ((_uint64)basesProcessed - batches[whichHashTable].nUsed - stats->unrecordedSkippedSeeds)/printPeriod) {
            WriteStatusMessage("Indexing %lld / %lld\n",(basesProcessed/printPeriod)*printPeriod, (_int64)context->genome->getCountOfBases());
        }

        stats->unrecordedSkippedSeeds = 0; // All except the first time through the loop this will be 0.        
        AcquireExclusiveLock(&context->hashTableLocks[whichHashTable]);
        for (unsigned i = 0; i < batches[whichHashTable].nUsed; i++) {
            ApplyHashTableUpdate(context, whichHashTable, batches[whichHashTable].entries[i].genomeLocation, 
                batches[whichHashTable].entries[i].lowBases, batches[whichHashTable].entries[i].usingComplement,
                &stats->bothComplementsUsed, &stats->genomeLocationsInOverflowTable, 
                &stats->seedsWithMultipleOccurrences, large);
        }
        ReleaseExclusiveLock(&context->hashTableLocks[whichHashTable]);
    }
}

GenomeIndex::OverflowBackpointerAnchor::OverflowBackpointerAnchor(unsigned maxOverflowEntries_) : maxOverflowEntries(maxOverflowEntries_)
{
	_int64 roundedUpMaxOverflowEntries = ((_int64)maxOverflowEntries + batchSize - 1) / batchSize * batchSize;	// Round up to the next batch size

	table = new OverflowBackpointer *[roundedUpMaxOverflowEntries / batchSize];

	for (unsigned i = 0; i < maxOverflowEntries / batchSize; i++) {
		table[i] = NULL;
	}

    InitializeExclusiveLock(&lock);
}

GenomeIndex::OverflowBackpointerAnchor::~OverflowBackpointerAnchor()
{
	for (unsigned i = 0; i < maxOverflowEntries / batchSize; i++) {
		if (table[i] != NULL) {
			BigDealloc(table[i]);
			table[i] = NULL;
		}
	}

	delete [] table;
	table = NULL;

    DestroyExclusiveLock(&lock);
}

	GenomeIndex::OverflowBackpointer *
GenomeIndex::OverflowBackpointerAnchor::getBackpointer(unsigned index)
{
    if (index >= maxOverflowEntries) {
        WriteErrorMessage("Trying to use too many overflow entries.  You can't index this reference genome with this seed size.  Try a larger one.\n");
        soft_exit(1);
    }
	unsigned tableSlot = index / batchSize;
	if (table[tableSlot] == NULL) {
        AcquireExclusiveLock(&lock);
        if (table[tableSlot] == NULL) {
		    table[tableSlot] = (OverflowBackpointer *)BigAlloc(batchSize * sizeof(OverflowBackpointer));
		    for (unsigned i = 0; i < batchSize; i++) {
			    table[tableSlot][i].genomeOffset = 0xffffffff;
			    table[tableSlot][i].nextIndex = 0xffffffff;
		    }
        }
        ReleaseExclusiveLock(&lock);
	}
	return &table[tableSlot][index % batchSize];
}

const unsigned GenomeIndex::OverflowBackpointerAnchor::batchSize = 1024 * 1024;

    void
GenomeIndex::printBiasTables()
{
    for (int keySize = 0; keySize <= largestKeySize; keySize++) {
        for (int seedSize = 0; seedSize <= largestBiasTable; seedSize++) {
            if (NULL != hg19_biasTables_large[keySize][seedSize]) {
                printf("static double hg19_biasTable%d_%d_large[] = {\n", seedSize, keySize);
                unsigned bitsOfSeed = seedSize * 2;
                unsigned bitsOfKey = keySize * 8;   // 8 == nBits / byte
                unsigned numHashTables = 1 << (bitsOfSeed - bitsOfKey);
                for (unsigned hashTable = 0; hashTable < numHashTables; hashTable++) {
                    if (hg19_biasTables_large[keySize][seedSize][hashTable] == 0) {
                        printf("0");
                    } else if (hg19_biasTables_large[keySize][seedSize][hashTable] > 1000 || hg19_biasTables_large[keySize][seedSize][hashTable] < .01) {
                        printf("%1.2e", hg19_biasTables_large[keySize][seedSize][hashTable]);
                    } else if (hg19_biasTables_large[keySize][seedSize][hashTable] > 10) {
                        printf("%1.1f", hg19_biasTables_large[keySize][seedSize][hashTable]);
                    } else {
                        printf("%1.3f", hg19_biasTables_large[keySize][seedSize][hashTable]);
                    }

                    if (hashTable != numHashTables-1) {
                        printf(",");
                    }

                    if (hashTable % 10 == 9) {
                        printf("\n");
                    }
                } // for each hash table
                printf("\n};\n\n");
            } // if there is a bias entry for this key * seed size
        } // for each seed size
    } // for each key size

    for (int keySize = 0; keySize <= largestKeySize; keySize++) {
        for (int seedSize = 0; seedSize <= largestBiasTable; seedSize++) {
            if (NULL != hg19_biasTables[keySize][seedSize]) {
                printf("static double hg19_biasTable%d_%d[] = {\n", seedSize, keySize);
                unsigned bitsOfSeed = seedSize * 2;
                unsigned bitsOfKey = keySize * 8;   // 8 == nBits / byte
                unsigned numHashTables = 1 << (bitsOfSeed - bitsOfKey);
                for (unsigned hashTable = 0; hashTable < numHashTables; hashTable++) {
                    if (hg19_biasTables[keySize][seedSize][hashTable] == 0) {
                        printf("0");
                    } else if (hg19_biasTables[keySize][seedSize][hashTable] > 1000 || hg19_biasTables[keySize][seedSize][hashTable] < .01) {
                        printf("%1.2e", hg19_biasTables[keySize][seedSize][hashTable]);
                    } else if (hg19_biasTables[keySize][seedSize][hashTable] > 10) {
                        printf("%1.1f", hg19_biasTables[keySize][seedSize][hashTable]);
                    } else {
                        printf("%1.3f", hg19_biasTables[keySize][seedSize][hashTable]);
                    }

                    if (hashTable != numHashTables-1) {
                        printf(",");
                    }

                    if (hashTable % 10 == 9) {
                        printf("\n");
                    }
                } // for each hash table
                printf("\n};\n\n");
            } // if there is a bias entry for this key * seed size
        } // for each seed size
    } // for each key size
}

        GenomeIndex *
GenomeIndex::loadFromDirectory(char *directoryName)
{

    const unsigned filenameBufferSize = MAX_PATH+1;
    char filenameBuffer[filenameBufferSize];
    
    snprintf(filenameBuffer,filenameBufferSize,"%s%cGenomeIndex",directoryName,PATH_SEP);

    GenericFile *indexFile = GenericFile::open(filenameBuffer, GenericFile::ReadOnly);

    if (NULL == indexFile) {
        WriteErrorMessage("Unable to open file '%s' for read.\n",filenameBuffer);
        return NULL;
    }

    char indexFileBuf[1000];
    indexFile->read(indexFileBuf, sizeof(indexFileBuf));

    unsigned seedLen;
    unsigned majorVersion, minorVersion, chromosomePadding;
    int nRead;
    size_t hashTablesFileSize;
    unsigned nHashTables;
    unsigned overflowTableSize;
    unsigned hashTableKeySize;
    unsigned smallHashTable;
    if (9 != (nRead = sscanf(indexFileBuf,"%d %d %d %d %d %d %d %lld %d", &majorVersion, &minorVersion, &nHashTables, &overflowTableSize, &seedLen, &chromosomePadding, 
											&hashTableKeySize, &hashTablesFileSize, &smallHashTable))) {
        if (3 == nRead || 6 == nRead || 7 == nRead) {
            WriteErrorMessage("Indices built by versions before 1.0dev.21 are no longer supported.  Please rebuild your index.\n");
        } else {
            WriteErrorMessage("GenomeIndex::LoadFromDirectory: didn't read initial values\n");
        }
        indexFile->close();		
        delete indexFile;
        return NULL;
    }
    indexFile->close();
    delete indexFile;

    if (majorVersion != GenomeIndexFormatMajorVersion) {
        WriteErrorMessage("This genome index appears to be from a different version of SNAP than this, and so we can't read it.  Index version %d, SNAP index format version %d\n",
            majorVersion, GenomeIndexFormatMajorVersion);
        soft_exit(1);
    }

    if (0 == seedLen) {
        WriteErrorMessage("GenomeIndex::LoadFromDirectory: saw seed size of 0.\n");
        return NULL;
    }

    GenomeIndex *index;
	if (smallHashTable) {
		index = new GenomeIndexSmall();
	} else {
		index = new GenomeIndexLarge();
	}
    index->nHashTables = nHashTables;
    index->overflowTableSize = overflowTableSize;
    index->hashTableKeySize = hashTableKeySize;
    index->seedLen = seedLen;

    size_t overflowTableSizeInBytes = (size_t)index->overflowTableSize * sizeof(*(index->overflowTable));
    index->overflowTable = (unsigned *)BigAlloc(overflowTableSizeInBytes);

    snprintf(filenameBuffer,filenameBufferSize,"%s%cOverflowTable",directoryName,PATH_SEP);
    GenericFile *fOverflowTable = GenericFile::open(filenameBuffer, GenericFile::ReadOnly);

    if (NULL == fOverflowTable) {
        WriteErrorMessage("Unable to open overflow table file, '%s', %d\n",filenameBuffer,errno);
        delete index;
        return NULL;
    }

    size_t amountRead = fOverflowTable->read(index->overflowTable, overflowTableSizeInBytes);
    if (amountRead != overflowTableSizeInBytes) {
        WriteErrorMessage("Error reading overflow table, %lld != %lld bytes read.\n", amountRead, overflowTableSizeInBytes);
        soft_exit(1);
    }

    fOverflowTable->close();
    delete fOverflowTable;
    fOverflowTable = NULL;

    index->hashTables = new SNAPHashTable*[index->nHashTables];

    for (unsigned i = 0; i < index->nHashTables; i++) {
        index->hashTables[i] = NULL; // We need to do this so the destructor doesn't crash if loading a hash table fails.
    }

    snprintf(filenameBuffer,filenameBufferSize,"%s%cGenomeIndexHash",directoryName,PATH_SEP);
	GenericFile *tablesFile = GenericFile::open(filenameBuffer, GenericFile::ReadOnly);
    if (NULL == tablesFile) {
        WriteErrorMessage("Unable to open genome hash table file '%s'\n", filenameBuffer);
        soft_exit(1);
    }

    index->tablesBlob = BigAlloc(hashTablesFileSize);
    amountRead = tablesFile->read(index->tablesBlob, hashTablesFileSize);
    if (amountRead != hashTablesFileSize) {
        WriteErrorMessage("Read incorrect amount for GenomeIndexHash file, %lld != %lld\n", hashTablesFileSize, amountRead);
        soft_exit(1);
    }

    GenericFile_Blob *blobFile = GenericFile_Blob::open(index->tablesBlob, hashTablesFileSize);

    for (unsigned i = 0; i < index->nHashTables; i++) {
        if (NULL == (index->hashTables[i] = SNAPHashTable::loadFromBlob(blobFile))) {
            WriteErrorMessage("GenomeIndex::loadFromDirectory: Failed to load hash table %d\n",i);
            delete index;
            return NULL;
        }

		unsigned expectedValueCount;
		if (smallHashTable) {
			expectedValueCount = 1;
		} else {
			expectedValueCount = 2;
		}

        if (index->hashTables[i]->GetValueCount() != expectedValueCount) {
            WriteErrorMessage("Expected loaded hash table to have value count of %d, but it had %d.  Index corrupt\n", expectedValueCount, index->hashTables[i]->GetValueCount());
            delete index;
            return NULL;
        }

        if (index->hashTables[i]->GetKeySizeInBytes() != 4) {
            WriteErrorMessage("Expected loaded hash table to have key size in bytes of 4, but it had %d.  This version of SNAP is too old to work with this index.\n", index->hashTables[i]->GetKeySizeInBytes());
            delete index;
            return NULL;
        }
    }

	tablesFile->close();
	delete tablesFile;
    tablesFile = NULL;

    blobFile->close();
    delete blobFile;
    blobFile = NULL;

    snprintf(filenameBuffer,filenameBufferSize,"%s%cGenome",directoryName,PATH_SEP);
    if (NULL == (index->genome = Genome::loadFromFile(filenameBuffer, chromosomePadding))) {
        WriteErrorMessage("GenomeIndex::loadFromDirectory: Failed to load the genome itself\n");
        delete index;
        return NULL;
    }

    if ((_int64)index->genome->getCountOfBases() + (_int64)index->overflowTableSize > 0xfffffff0) {
        WriteErrorMessage("\nThis index has too many overflow entries to be valid.  Some early versions of SNAP\n"
                        "allowed building indices with too small of a seed size, and this appears to be such\n"
                        "an index.  You can no longer build indices like this, and you also can't use them\n"
                        "because they are corrupt and would produce incorrect results.  Please use an index\n"
                        "built with a larger seed size.  For hg19, the seed size must be at least 19.\n"
                        "For other reference genomes this quantity will vary.\n");
        soft_exit(1);
    }

    return index;
}

        void
GenomeIndexLarge::lookupSeed(Seed seed, unsigned *nHits, const unsigned **hits, unsigned *nRCHits, const unsigned **rcHits)
{
    return lookupSeed(seed, 0, 0xFFFFFFFF, nHits, hits, nRCHits, rcHits);
}

    void
GenomeIndexLarge::lookupSeed(
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
    _ASSERT(hashTables[seed.getHighBases(hashTableKeySize)]->GetValueSizeInBytes() == 4);
    unsigned *entry = (unsigned int *)hashTables[seed.getHighBases(hashTableKeySize)]->GetFirstValueForKey(lowBases);   // Cast OK because valueSize == 4
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
GenomeIndexSmall::lookupSeed(Seed seed, unsigned *nHits, const unsigned **hits, unsigned *nRCHits, const unsigned **rcHits)
{
    return lookupSeed(seed, 0, 0xFFFFFFFF, nHits, hits, nRCHits, rcHits);
}

    void
GenomeIndexSmall::lookupSeed(
    Seed              seed,
    unsigned          minLocation,
    unsigned          maxLocation,
    unsigned         *nHits,
    const unsigned  **hits,
    unsigned         *nRCHits,
    const unsigned  **rcHits)
{
	for (int dir = 0; dir < NUM_DIRECTIONS; dir++) {
		_ASSERT(seed.getHighBases(hashTableKeySize) < nHashTables);
		_uint64 lowBases = seed.getLowBases(hashTableKeySize);
		_ASSERT(hashTables[seed.getHighBases(hashTableKeySize)]->GetValueSizeInBytes() == 4);
		unsigned *entry = (unsigned int *)hashTables[seed.getHighBases(hashTableKeySize)]->GetFirstValueForKey(lowBases);   // Cast OK because valueSize == 4
		if (NULL == entry) {
			if (FORWARD == dir) {
				*nHits = 0;
			} else {
				*nRCHits = 0;
			}
		} else if (FORWARD == dir) {
			fillInLookedUpResults(entry, minLocation, maxLocation, nHits, hits);
		} else {
			fillInLookedUpResults(entry, minLocation, maxLocation, nRCHits, rcHits);
		}
		seed = ~seed;
    }	// For each direction
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
