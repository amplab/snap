/*++

Module Name:

    PairedAligner.cpp

Abstract:

    Functions for running the paired end aligner sub-program.


Authors:

    Hongyi Xin, June, 2016

Environment:

    User mode service.

Revision History:

    Adapted from cSNAP, which was in turn adapted from the scala prototype

--*/

//
// TODO: This is really similar to the single-end aligner overall. It would be nice
// to avoid code duplication.
//


#include "stdafx.h"
#include "options.h"
#include <time.h>
#include "Compat.h"
#include "RangeSplitter.h"
#include "GenomeIndex.h"
#include "SAM.h"
#include "TenXClusterAligner.h"
#include "Tables.h"
#include "AlignerOptions.h"
#include "AlignerContext.h"
#include "AlignerStats.h"
#include "FASTQ.h"
#include "TenXAligner.h"
#include "MultiInputReadSupplier.h"
#include "Util.h"
#include "TenXSingleAligner.h"
#include "exit.h"
#include "Error.h"

using namespace std;

using util::stringEndsWith;

static const int DEFAULT_MIN_SPACING = 50;
static const int DEFAULT_MAX_SPACING = 1000;
static const unsigned DEFAULT_MAX_BARCODE_SIZE = 60000;
static const unsigned DEFAULT_MIN_PAIRS_PER_CLUSTER = 10;
static const unsigned DEFAULT_MIN_CLUSTER_SPAN = 100000;
static const double DEFAULT_UNCLUSTERED_PENALTY = 0.000000001;
static const unsigned DEFAULT_CLUSTER_ED_COMPENSATION = 4;
static const unsigned DEFAULT_MAX_CLUSTER_SIZE = 100000;

struct TenXAlignerStats : public AlignerStats
{
    // TODO: make these constants configurable
    static const int MAX_DISTANCE = 1000;
    static const int MAX_SCORE = 15;

    _int64 sameComplement;
    _int64* distanceCounts; // histogram of distances
    // TODO: could save a bit of memory & time since this is a triangular matrix
    _int64* scoreCounts; // 2-d histogram of scores for paired ends
    static const unsigned maxMapq = 70;
    static const unsigned nTimeBuckets = 32;
    static const unsigned nHitsBuckets = 32;
    static const unsigned nLVCallsBuckets = 32;

    _int64 alignTogetherByMapqHistogram[maxMapq + 1][nTimeBuckets];
    _int64 totalTimeByMapqHistogram[maxMapq + 1][nTimeBuckets];
    _int64 nSmallHitsByTimeHistogram[nHitsBuckets][nTimeBuckets];
    _int64 nLVCallsByTimeHistogram[nLVCallsBuckets][nTimeBuckets];
    _int64 mapqByNLVCallsHistogram[maxMapq + 1][nLVCallsBuckets];
    _int64 mapqByNSmallHitsHistogram[maxMapq + 1][nHitsBuckets];

    TenXAlignerStats(AbstractStats* i_extra = NULL);

    virtual ~TenXAlignerStats();

    inline void incrementDistance(int distance) {
        distanceCounts[max(0, min(MAX_DISTANCE, distance))]++;
    }

    inline void incrementScore(int s0, int s1)
    {
        // ensure s0 <= s1, both within range
        s0 = max(0, min(MAX_SCORE, s0));
        s1 = max(0, min(MAX_SCORE, s1));
        if (s0 > s1) {
            int t = s0; s0 = s1; s1 = t;
        }
        scoreCounts[s0*(MAX_SCORE + 1) + s1]++;
    }

    inline void recordAlignTogetherMapqAndTime(unsigned mapq, _int64 timeInNanos, unsigned nSmallHits, unsigned nLVCalls) {
        int timeBucket;
        _int64 dividedTime = timeInNanos;
        for (timeBucket = 0; timeBucket < nTimeBuckets - 1; timeBucket++) {
            if (dividedTime == 0) break;
            dividedTime /= 2;
        }

        alignTogetherByMapqHistogram[mapq][timeBucket]++;
        totalTimeByMapqHistogram[mapq][timeBucket] += timeInNanos;

        int nHitsBucket;
        int dividedHits = nSmallHits;
        for (nHitsBucket = 0; nHitsBucket < nHitsBuckets; nHitsBucket++) {
            if (0 == dividedHits) break;
            dividedHits /= 2;
        }
        _ASSERT((char *)&nSmallHitsByTimeHistogram[nHitsBucket][timeBucket] < (char *)(this + 1));
        nSmallHitsByTimeHistogram[nHitsBucket][timeBucket]++;

        int nLVCallsBucket;
        int dividedLVCalls = nLVCalls;
        for (nLVCallsBucket = 0; nLVCallsBucket < nLVCallsBuckets; nLVCallsBucket++) {
            if (dividedLVCalls == 0) break;
            dividedLVCalls /= 2;
        }
        _ASSERT((char *)&nLVCallsByTimeHistogram[nLVCallsBucket][timeBucket] < (char *)(this + 1));
        nLVCallsByTimeHistogram[nLVCallsBucket][timeBucket]++;

        _ASSERT((char *)&mapqByNLVCallsHistogram[mapq][nLVCallsBucket] < (char *)(this + 1));
        mapqByNLVCallsHistogram[mapq][nLVCallsBucket]++;

        _ASSERT((char *)&mapqByNSmallHitsHistogram[mapq][nHitsBucket] < (char *)(this + 1));
        mapqByNSmallHitsHistogram[mapq][nHitsBucket]++;
    }



    virtual void add(const AbstractStats * other);

    virtual void printHistograms(FILE* output);
};

const int TenXAlignerStats::MAX_DISTANCE;
const int TenXAlignerStats::MAX_SCORE;

TenXAlignerStats::TenXAlignerStats(AbstractStats* i_extra)
    : AlignerStats(i_extra),
    sameComplement(0)
{
    int dsize = sizeof(_int64) * (MAX_DISTANCE + 1);
    distanceCounts = (_int64*)BigAlloc(dsize);
    memset(distanceCounts, 0, dsize);

    int ssize = sizeof(_int64) * (MAX_SCORE + 1)*(MAX_SCORE + 1);
    scoreCounts = (_int64*)BigAlloc(ssize);
    memset(scoreCounts, 0, ssize);

    for (unsigned mapq = 0; mapq <= maxMapq; mapq++) {
        for (unsigned timeBucket = 0; timeBucket < nTimeBuckets; timeBucket++) {
            alignTogetherByMapqHistogram[mapq][timeBucket] = 0;
            totalTimeByMapqHistogram[mapq][timeBucket] = 0;
        }
        for (unsigned smallHits = 0; smallHits < nHitsBuckets; smallHits++) {
            mapqByNSmallHitsHistogram[mapq][smallHits] = 0;
        }
        for (unsigned lvCalls = 0; lvCalls < nLVCallsBuckets; lvCalls++) {
            mapqByNLVCallsHistogram[mapq][lvCalls] = 0;
        }
    }

    for (unsigned timeBucket = 0; timeBucket < nTimeBuckets; timeBucket++) {
        for (unsigned smallHits = 0; smallHits < nHitsBuckets; smallHits++) {
            nSmallHitsByTimeHistogram[smallHits][timeBucket] = 0;
        }
        for (unsigned lvCalls = 0; lvCalls < nLVCallsBuckets; lvCalls++) {
            nLVCallsByTimeHistogram[lvCalls][timeBucket] = 0;
        }
    }

}

TenXAlignerStats::~TenXAlignerStats()
{
    BigDealloc(distanceCounts);
    BigDealloc(scoreCounts);
}

void TenXAlignerStats::add(const AbstractStats * i_other)
{
    AlignerStats::add(i_other);
    TenXAlignerStats* other = (TenXAlignerStats*)i_other;
    for (int i = 0; i < MAX_DISTANCE + 1; i++) {
        distanceCounts[i] += other->distanceCounts[i];
    }
    for (int i = 0; i < (MAX_SCORE + 1) * (MAX_SCORE + 1); i++) {
        scoreCounts[i] += other->scoreCounts[i];
    }

    for (unsigned mapq = 0; mapq <= maxMapq; mapq++) {
        for (unsigned timeBucket = 0; timeBucket < nTimeBuckets; timeBucket++) {
            alignTogetherByMapqHistogram[mapq][timeBucket] += other->alignTogetherByMapqHistogram[mapq][timeBucket];
            totalTimeByMapqHistogram[mapq][timeBucket] += other->totalTimeByMapqHistogram[mapq][timeBucket];
        }
        for (unsigned smallHits = 0; smallHits < nHitsBuckets; smallHits++) {
            mapqByNSmallHitsHistogram[mapq][smallHits] += other->mapqByNSmallHitsHistogram[mapq][smallHits];
        }
        for (unsigned lvCalls = 0; lvCalls < nLVCallsBuckets; lvCalls++) {
            mapqByNLVCallsHistogram[mapq][lvCalls] += other->mapqByNLVCallsHistogram[mapq][lvCalls];
        }
    }

    for (unsigned timeBucket = 0; timeBucket < nTimeBuckets; timeBucket++) {
        for (unsigned smallHits = 0; smallHits < nHitsBuckets; smallHits++) {
            nSmallHitsByTimeHistogram[smallHits][timeBucket] += other->nSmallHitsByTimeHistogram[smallHits][timeBucket];
        }
        for (unsigned lvCalls = 0; lvCalls < nLVCallsBuckets; lvCalls++) {
            nLVCallsByTimeHistogram[lvCalls][timeBucket] += other->nLVCallsByTimeHistogram[lvCalls][timeBucket];
        }
    }

}

void TenXAlignerStats::printHistograms(FILE* output)
{
    AlignerStats::printHistograms(output);
}

TenXAlignerOptions::TenXAlignerOptions(const char* i_commandLine)
    : AlignerOptions(i_commandLine, true),
    // indicator that secondary result overflows
    minSpacing(DEFAULT_MIN_SPACING),
    maxSpacing(DEFAULT_MAX_SPACING),

    // 10x specific parameter
    maxBarcodeSize(DEFAULT_MAX_BARCODE_SIZE),
    minPairsPerCluster(DEFAULT_MIN_PAIRS_PER_CLUSTER),
    minClusterSpan(DEFAULT_MIN_CLUSTER_SPAN),
    unclusteredPenalty(DEFAULT_UNCLUSTERED_PENALTY),
    clusterEDCompensation(DEFAULT_CLUSTER_ED_COMPENSATION),
    maxClusterNum(DEFAULT_MAX_CLUSTER_SIZE),

    // same with pairedEndAligner
    forceSpacing(false),
    noSingle(false),
    intersectingAlignerMaxHits(DEFAULT_INTERSECTING_ALIGNER_MAX_HITS),
    maxCandidatePoolSize(DEFAULT_MAX_CANDIDATE_POOL_SIZE),
    quicklyDropUnpairedReads(true)
{
}

void TenXAlignerOptions::usageMessage()
{
    AlignerOptions::usageMessage();
    WriteErrorMessage(
        "\n"
        "  -s   min and max spacing to allow between paired ends (default: %d %d).\n"
        "  -fs  force spacing to lie between min and max.\n"
        "  -H   max hits for intersecting aligner (default: %d).\n"
        "  -mcp specifies the maximum candidate pool size (An internal data structure. \n"
        "       Only increase this if you get an error message saying to do so. If you're running\n"
        "       out of memory, you may want to reduce it.  Default: %d)\n"
        "  -F b additional option to -F to require both mates to satisfy filter (default is just one)\n"
        "       If you specify -F b together with one of the other -F options, -F b MUST be second\n"
        "  -ku  Keep unpaired-looking reads in SAM/BAM input.  Ordinarily, if a read doesn't specify\n"
        "       mate information (RNEXT field is * and/or PNEXT is 0) then the code that matches reads will immdeiately\n"
        "       discard it.  Specifying this flag may cause large memory usage for some input files,\n"
        "       but may be necessary for some strangely formatted input files.  You'll also need to specify this\n"
        "       flag for SAM/BAM files that were aligned by a single-end aligner.\n"
        "  -noS no single-end aligner turned on. Just use paired-end aligner.\n"
        "\nTenXTenXTenXTenXTenXTenX    Below are 10X specific parameters    TenXTenXTenXTenXTenXTenX\n\n"
        "  -maxBar           specifies the maximum reads in a barcode. This has a direct affect on memory usage. Default: %d\n"
        "  -minClusterPairs  specifies the maximum number of pair-ended reads for a cluster to be considered legitimate. Default: %d\n"
        "  -minClusterSpan   specifies the maximum scope that SNAP detects clusters. If # of reads within minClusterSpan is\n"
        "                    smaller than minClusterReads the cluster is rejected. Default: %d\n"
        "  -UCP              specifies the unclustered probability penalty. Default: %f\n"
        "  -CED              specifies the cluster edit-distance compensation. Default: %d\n"
        "  -mCS              specifies the maximum potential cluster size. Increase it would increase the memory usage.\n"
        "                    This parameter will be taken out in the future. Default: %d\n"
        ,
        DEFAULT_MIN_SPACING,
        DEFAULT_MAX_SPACING,
        DEFAULT_INTERSECTING_ALIGNER_MAX_HITS,
        DEFAULT_MAX_CANDIDATE_POOL_SIZE,

        DEFAULT_MAX_BARCODE_SIZE,
        DEFAULT_MIN_PAIRS_PER_CLUSTER,
        DEFAULT_MIN_CLUSTER_SPAN,
        DEFAULT_UNCLUSTERED_PENALTY,
        DEFAULT_CLUSTER_ED_COMPENSATION,
        DEFAULT_MAX_CLUSTER_SIZE
        );
}

bool TenXAlignerOptions::parse(const char** argv, int argc, int& n, bool *done)
{
    *done = false;

    if (strcmp(argv[n], "-s") == 0) {
        if (n + 2 < argc) {
            minSpacing = atoi(argv[n + 1]);
            maxSpacing = atoi(argv[n + 2]);
            n += 2;
            return true;
        }
        return false;
    }
    else if (strcmp(argv[n], "-H") == 0) {
        if (n + 1 < argc) {
            intersectingAlignerMaxHits = atoi(argv[n + 1]);
            n += 1;
            return true;
        }
        return false;
    }
    else if (strcmp(argv[n], "-fs") == 0) {
        forceSpacing = true;
        return true;
    }
    else if (strcmp(argv[n], "-ku") == 0) {
        quicklyDropUnpairedReads = false;
        return true;
    }
    else if (strcmp(argv[n], "-mcp") == 0) {
        if (n + 1 < argc) {
            maxCandidatePoolSize = atoi(argv[n + 1]);
            n += 1;
            return true;
        }
        return false;
    }
    else if (strcmp(argv[n], "-F") == 0 && n + 1 < argc && strcmp(argv[n + 1], "b") == 0) {
        filterFlags |= FilterBothMatesMatch;
        n += 1;
        return true;
    }
    else if (strcmp(argv[n], "-maxBar") == 0) {
        if (n + 1 < argc) {
            maxBarcodeSize = atoi(argv[n + 1]);
            n += 1;
            return true;
        }
        return false;
    }
    else if (strcmp(argv[n], "-minClusterPairs") == 0) {
        if (n + 1 < argc) {
            minPairsPerCluster = atoi(argv[n + 1]);
            n += 1;
            return true;
        }
        return false;
    }
    else if (strcmp(argv[n], "-minClusterSpan") == 0) {
        if (n + 1 < argc) {
            minClusterSpan = atoi(argv[n + 1]);
            n += 1;
            return true;
        }
        return false;
    }
    else if (strcmp(argv[n], "-noS") == 0) {
        noSingle = true;
        return true;
    }
    else if (strcmp(argv[n], "-uCP") == 0) {
        if (n + 1 < argc) {
            unclusteredPenalty = atof(argv[n + 1]);
            n += 1;
            return true;
        }
        return false;
    }
    else if (strcmp(argv[n], "-CED") == 0) {
        if (n + 1 < argc) {
            clusterEDCompensation = atoi(argv[n + 1]);
            n += 1;
            return true;
        }
        return false;
    }
    else if (strcmp(argv[n], "-mCS") == 0) {
        if (n + 1 < argc) {
            maxClusterNum = atoi(argv[n + 1]);
            n += 1;
            return true;
        }
        return false;
    }
    return AlignerOptions::parse(argv, argc, n, done);
}

TenXAlignerContext::TenXAlignerContext(AlignerExtension* i_extension)
    : AlignerContext(0, NULL, NULL, i_extension)
{
}

bool TenXAlignerContext::initialize()
{
    AlignerContext::initialize();
    TenXAlignerOptions* options2 = (TenXAlignerOptions*)options;
    minSpacing = options2->minSpacing;
    maxSpacing = options2->maxSpacing;
    maxBarcodeSize = options2->maxBarcodeSize;
    minPairsPerCluster = options2->minPairsPerCluster;
    minClusterSpan = options2->minClusterSpan;
    forceSpacing = options2->forceSpacing;
    noSingle = options2->noSingle;
    unclusteredPenalty = options2->unclusteredPenalty;
    clusterEDCompensation = options2->clusterEDCompensation;
    maxCandidatePoolSize = options2->maxCandidatePoolSize;
    intersectingAlignerMaxHits = options2->intersectingAlignerMaxHits;
    ignoreMismatchedIDs = options2->ignoreMismatchedIDs;
    quicklyDropUnpairedReads = options2->quicklyDropUnpairedReads;
    noUkkonen = options->noUkkonen;
    noOrderedEvaluation = options->noOrderedEvaluation;

    return true;
}

AlignerStats* TenXAlignerContext::newStats()
{
    return new TenXAlignerStats();
}

void TenXAlignerContext::runTask()
{
    ParallelTask<TenXAlignerContext> task(this);
    task.run();
}


void TenXAlignerContext::runIterationThread()
{
    PreventMachineHibernationWhileThisThreadIsAlive();

    PairedReadSupplier *supplier = pairedReadSupplierGenerator->generateNewPairedReadSupplier();

    if (NULL == supplier) {
        //
        // No work for this thread to do.
        //
        return;
    }

    if (extension->runIterationThread(supplier, this)) {
        delete supplier;
        return;
    }

    if (index == NULL) { //**** not sure what's going on here
        Read *reads[NUM_READS_PER_PAIR];
        _int64 nSingleResults[2] = { 0, 0 };

        // no alignment, just input/output
        PairedAlignmentResult result;
        memset(&result, 0, sizeof(result));
        result.location[0] = result.location[1] = InvalidGenomeLocation;


        while (supplier->getNextReadPair(&reads[0], &reads[1])) {
            // Check that the two IDs form a pair; they will usually be foo/1 and foo/2 for some foo.
            if (!ignoreMismatchedIDs && !readIdsMatch(reads[0], reads[1])) {
                unsigned n[2] = { min(reads[0]->getIdLength(), 200u), min(reads[1]->getIdLength(), 200u) };
                char* p[2] = { (char*)alloca(n[0] + 1), (char*)alloca(n[1] + 1) };
                memcpy(p[0], reads[0]->getId(), n[0]); p[0][n[0]] = 0;
                memcpy(p[1], reads[1]->getId(), n[1]); p[1][n[1]] = 0;
                WriteErrorMessage("Unmatched read IDs '%s' and '%s'.  Use the -I option to ignore this.\n", p[0], p[1]);
                soft_exit(1);
            }
            stats->totalReads += 2;

            bool pass0 = options->passFilter(reads[0], result.status[0], reads[0]->getDataLength() >= minReadLength && (int)reads[0]->countOfNs() <= maxDist, false);
            bool pass1 = options->passFilter(reads[1], result.status[1], reads[1]->getDataLength() >= minReadLength && (int)reads[1]->countOfNs() <= maxDist, false);
            bool pass = (options->filterFlags & AlignerOptions::FilterBothMatesMatch)
                ? (pass0 && pass1) : (pass0 || pass1);

            if (pass) {
                stats->notFound++;
                if (NULL != readWriter) {
                    readWriter->writePairs(readerContext, reads, &result, 1, NULL, nSingleResults, true);
                }
            }
            else {
                stats->uselessReads++;
            }
        }
        delete supplier;
        return;
    }


    /*
     * Initialize some constants
     */

    int maxReadSize = MAX_READ_LENGTH;

    _int64 _maxPairedSecondaryHits_ref;
    _int64 _maxSingleSecondaryHits_ref;

    if (maxSecondaryAlignmentAdditionalEditDistance < 0) {
        _maxPairedSecondaryHits_ref = 0;
        _maxSingleSecondaryHits_ref = 0;
    }
    else {
        //
        // Since we reallocate these if they overflow, just pick a value that doesn't waste too much memory.
        //
        _maxPairedSecondaryHits_ref = 32;
        _maxSingleSecondaryHits_ref = 32;
    }

    //fprintf(stderr, "****initializing\n");

    /*
     * calculate the memory useage for reservation
     */

     // memory quota for the ClusterAlinger and it's internal SingleAligners
    size_t memoryPoolSize = TenXClusterAligner::getBigAllocatorReservation(index, maxReadSize, maxHits, index->getSeedLength(), numSeedsFromCommandLine, seedCoverage, maxDist,
        extraSearchDepth, maxCandidatePoolSize, maxSecondaryAlignmentsPerContig);

    //fprintf(stderr, "****memoryPoolSize after TenXClusterAligner reservation: %lld\n", memoryPoolSize);

    // memory quota for all the SingleAligners 
    size_t singleTenXReserve = TenXSingleAligner::getBigAllocatorReservation(index, intersectingAlignerMaxHits, maxReadSize, index->getSeedLength(),
        numSeedsFromCommandLine, seedCoverage, maxDist, extraSearchDepth, maxCandidatePoolSize,
        maxSecondaryAlignmentsPerContig);

    memoryPoolSize += singleTenXReserve * maxBarcodeSize;

    //fprintf(stderr, "****singleTenXReserve:%lld  maxBarcodeSize: %lld  memoryPoolSize: %lld\n", singleTenXReserve, maxBarcodeSize, memoryPoolSize);


    /*
     * Allocate space
     */

    BigAllocator *allocator = new BigAllocator(memoryPoolSize);

    // Allocate the single aligners pointers (single + cluster)
    TenXProgressTracker *tenXSingleTrackerArray = (TenXProgressTracker*)BigAlloc(sizeof(TenXProgressTracker) * maxBarcodeSize);

    // Allocate the clusterIsValid for clusterAligner
    bool *tenXClusterIsValid = (bool*)BigAlloc(sizeof(bool) * maxBarcodeSize);

    // Allocate the clusterCounter for clusterAligner
    unsigned *tenXClusterPairsCounter = (unsigned*)BigAlloc(sizeof(unsigned) * maxBarcodeSize);

    //fprintf(stderr, "****Before going into the loop of allocating single aligners\n");

    // Allocate and initialize erasers
    bool *clusterCounterEraser = (bool*)BigAlloc(sizeof(bool) * maxClusterNum);
    _uint8 *clusterToggleEraser = (_uint8*)BigAlloc(sizeof(_uint8) * maxClusterNum);
    for (int clusterIdx = 0; clusterIdx < maxClusterNum; clusterIdx++) {
        clusterCounterEraser[clusterIdx] = 0;
        clusterToggleEraser[clusterIdx] = false;
    }

    fprintf(stderr, "****TenX allocating\n");
    // Allocate shared cluster counter array    
    _uint8 *sharedClusterCounterAry = (_uint8*)BigAlloc(sizeof(_uint8) * maxClusterNum);

    for (int singleAlignerIdx = 0; singleAlignerIdx < maxBarcodeSize; singleAlignerIdx++) {
        // Allocate and initilazing tracker
        tenXSingleTrackerArray[singleAlignerIdx].pairNotDone = false;
        tenXSingleTrackerArray[singleAlignerIdx].singleNotDone = false;
        tenXSingleTrackerArray[singleAlignerIdx].clusterCounterAry = sharedClusterCounterAry;
        tenXSingleTrackerArray[singleAlignerIdx].clusterToggle = (bool*)BigAlloc(sizeof(bool) * maxClusterNum);
        tenXSingleTrackerArray[singleAlignerIdx].nSecondaryResults = 0;
        tenXSingleTrackerArray[singleAlignerIdx].nSingleEndSecondaryResults[0] = tenXSingleTrackerArray[singleAlignerIdx].nSingleEndSecondaryResults[1] = 0;

        // Allocate and initilazing aligner
        tenXSingleTrackerArray[singleAlignerIdx].aligner = new (allocator) TenXSingleAligner(index, maxReadSize, maxHits, maxDist, numSeedsFromCommandLine,
            seedCoverage, minSpacing, maxSpacing, intersectingAlignerMaxHits, extraSearchDepth,
            maxCandidatePoolSize, maxSecondaryAlignmentsPerContig, allocator, noUkkonen, noOrderedEvaluation, noTruncation, ignoreAlignmentAdjustmentForOm, printStatsMapQLimit, clusterEDCompensation, unclusteredPenalty,
            tenXSingleTrackerArray[singleAlignerIdx].clusterCounterAry, tenXSingleTrackerArray[singleAlignerIdx].clusterToggle);
        
        // Initialization    
        memcpy(clusterCounterEraser, tenXSingleTrackerArray[singleAlignerIdx].clusterCounterAry, sizeof(_uint8) * maxClusterNum);
        memcpy(clusterToggleEraser, tenXSingleTrackerArray[singleAlignerIdx].clusterToggle, sizeof(_uint8) * maxClusterNum);
    }

    TenXClusterAligner *aligner = new (allocator) TenXClusterAligner(
        index,
        maxReadSize,
        maxHits,
        maxDist,
        numSeedsFromCommandLine,
        seedCoverage,
        minWeightToCheck,
        forceSpacing,
        extraSearchDepth,
        noUkkonen,
        noOrderedEvaluation,
        noTruncation,
        ignoreAlignmentAdjustmentForOm,
        tenXSingleTrackerArray,
        tenXClusterIsValid,
        tenXClusterPairsCounter,
        maxBarcodeSize,
        minPairsPerCluster,
        minClusterSpan,
        unclusteredPenalty,
        clusterEDCompensation,
        minReadLength,
        maxSecondaryAlignmentsPerContig,
        printStatsMapQLimit,
        maxSecondaryAlignmentAdditionalEditDistance,
        maxSecondaryAlignments,
        allocator);

    allocator->checkCanaries();

    for (unsigned pairIdx = 0; pairIdx < maxBarcodeSize; pairIdx++) {
        // Allocate data for result arrays
        tenXSingleTrackerArray[pairIdx].results = (PairedAlignmentResult*)BigAlloc((1 + _maxPairedSecondaryHits_ref) * sizeof(PairedAlignmentResult)); // all paired results. "+1" is for the primary result
        tenXSingleTrackerArray[pairIdx].singleEndSecondaryResults = (SingleAlignmentResult*)BigAlloc(_maxSingleSecondaryHits_ref * sizeof(SingleAlignmentResult)); // all single results.
    }


    /*
     * Initialization of some data
     */

    for (unsigned pairIdx = 0; pairIdx < maxBarcodeSize; pairIdx++) {
        tenXSingleTrackerArray[pairIdx].secondaryResultBufferSize = _maxPairedSecondaryHits_ref;
        tenXSingleTrackerArray[pairIdx].singleSecondaryBufferSize = _maxSingleSecondaryHits_ref;
        tenXSingleTrackerArray[pairIdx].pairNotDone = false;
        tenXSingleTrackerArray[pairIdx].singleNotDone = false;
    }

    //fprintf(stderr, "****begin read buffering\n");

    /*
     * Buffer all the reads
     */

    ReadWriter *readWriter = this->readWriter;

#ifdef  _MSC_VER
    if (options->useTimingBarrier) {
        if (0 == InterlockedDecrementAndReturnNewValue(nThreadsAllocatingMemory)) {
            AllowEventWaitersToProceed(memoryAllocationCompleteBarrier);
        }
        else {
            WaitForEvent(memoryAllocationCompleteBarrier);
        }
    }
#endif  // _MSC_VER

    // Timing variables
    _uint64 lastReportTime = timeInMillis();
    _uint64 readsWhenLastReported = 0;
    _int64 startTime = timeInMillis();
    _int64 readFinishedTime;

    // Record time
    if (options->profile) {
        readFinishedTime = timeInMillis();
        stats->millisReading += (readFinishedTime - startTime);
    }
    // ****This is not currently useful for the moment which we only handle a single cluster
    if (AlignerOptions::useHadoopErrorMessages && stats->totalReads % 10000 == 0 && timeInMillis() - lastReportTime > 10000) {
        fprintf(stderr, "reporter:counter:SNAP,readsAligned,%lu\n", stats->totalReads - readsWhenLastReported);
        readsWhenLastReported = stats->totalReads;
        lastReportTime = timeInMillis();
    }

    _int64 nSingleResults[NUM_READS_PER_PAIR] = { 0, 0 }; // just for the sake of outputing unmapped read pairs
    unsigned totalPairsForBarcode = 0; // total legitimate read pairs of this barcode

    while (supplier->getNextReadPair(&tenXSingleTrackerArray[totalPairsForBarcode].pairedReads[0], &tenXSingleTrackerArray[totalPairsForBarcode].pairedReads[1])) {
        // A bunch of pointers delegates... For readability!
        Read **reads = tenXSingleTrackerArray[totalPairsForBarcode].pairedReads;
        bool *useful = tenXSingleTrackerArray[totalPairsForBarcode].useful;

        
        // Check that the two IDs form a pair; they will usually be foo/1 and foo/2 for some foo.
        if (!ignoreMismatchedIDs) {
            Read::checkIdMatch(reads[0], reads[1]);
        }

        stats->totalReads += 2;

        // Skip the pair if there are too many Ns and/or they're too short
        int maxDist = this->maxDist;
        useful[0] = reads[0]->getDataLength() >= minReadLength && (int)reads[0]->countOfNs() <= maxDist;
        useful[1] = reads[1]->getDataLength() >= minReadLength && (int)reads[1]->countOfNs() <= maxDist;
        if (!useful[0] && !useful[1]) {
            PairedAlignmentResult result;
            result.status[0] = NotFound;
            result.status[1] = NotFound;
            result.location[0] = InvalidGenomeLocation;
            result.location[1] = InvalidGenomeLocation;
            nSingleResults[0] = nSingleResults[1] = 0;
            result.clippingForReadAdjustment[0] = result.clippingForReadAdjustment[1] = 0;

            bool pass0 = options->passFilter(reads[0], result.status[0], true, false);
            bool pass1 = options->passFilter(reads[1], result.status[1], true, false);
            bool pass = (options->filterFlags & AlignerOptions::FilterBothMatesMatch)
                ? (pass0 && pass1) : (pass0 || pass1);

            if (pass) {
                if (NULL != readWriter) {
                    readWriter->writePairs(readerContext, reads, &result, 1, NULL, nSingleResults, true);
                }
                stats->uselessReads += 2;
            }
            else {
                stats->filtered += 2;
            }
            continue;
        }
        // at least one of the ends is good
        tenXSingleTrackerArray[totalPairsForBarcode].singleNotDone = true;
        // examine pair if both ends are good
        if (!useful[0] || !useful[1])
            tenXSingleTrackerArray[totalPairsForBarcode].pairNotDone = false;
        else
            tenXSingleTrackerArray[totalPairsForBarcode].pairNotDone = true;


        // Debugging
        if (tenXSingleTrackerArray[totalPairsForBarcode].pairNotDone) {
            char id[42];

            memcpy(id, tenXSingleTrackerArray[totalPairsForBarcode].pairedReads[0]->getId(), 41);
            id[41] = '\0';
            tenXSingleTrackerArray[totalPairsForBarcode].id[0] = string(id);
            memcpy(id, tenXSingleTrackerArray[totalPairsForBarcode].pairedReads[1]->getId(), 41);
            id[41] = '\0';
            tenXSingleTrackerArray[totalPairsForBarcode].id[1] = string(id);
        }
        // Debugging


        totalPairsForBarcode++;
        _ASSERT(totalPairsForBarcode <= maxBarcodeSize);
        // Note that useful0 and useful1 will be discarded if neither of the read of the pair is useful
    }

    //fprintf(stderr, "****begin alignment\n");

    /*
     * Align the read pairs
     */

#if     TIME_HISTOGRAM
    _int64 startTime = timeInNanos();
#endif // TIME_HISTOGRAM

    fprintf(stderr, "****TenX stage 1\n");

    // Stage 1, get seeds and find paired locations while keeping track of clusters
    bool barcodeFinished = aligner->align_first_stage(totalPairsForBarcode);
    if (barcodeFinished)
        return;

    fprintf(stderr, "****TenX stage 2A\n");
    // Stage 2, calculate ED and store paired results
    aligner->align_second_stage_clustering();

    fprintf(stderr, "****TenX stage 2B\n");
    barcodeFinished = aligner->align_second_stage_check_reallocate();
    if (!barcodeFinished) {
        for (unsigned pairIdx = 0; pairIdx < totalPairsForBarcode; pairIdx++) {
            if (tenXSingleTrackerArray[pairIdx].pairNotDone && tenXSingleTrackerArray[pairIdx].nSecondaryResults > tenXSingleTrackerArray[pairIdx].secondaryResultBufferSize) {
                BigDealloc(tenXSingleTrackerArray[pairIdx].results);
                tenXSingleTrackerArray[pairIdx].results = NULL;
                // calculate shifts here. I actually don't get it why we prefer to always make it 2-folds
                unsigned multiple = (tenXSingleTrackerArray[pairIdx].nSecondaryResults - 1) / tenXSingleTrackerArray[pairIdx].secondaryResultBufferSize + 1;
                unsigned shifter = 2;
                while (multiple > shifter) {
                    shifter *= 2;
                }
                
                tenXSingleTrackerArray[pairIdx].secondaryResultBufferSize *= shifter;
                tenXSingleTrackerArray[pairIdx].results = (PairedAlignmentResult *)BigAlloc((tenXSingleTrackerArray[pairIdx].secondaryResultBufferSize + 1) * sizeof(PairedAlignmentResult));
            }
        }
    }

    fprintf(stderr, "****TenX stage 2C\n");
    aligner->align_second_stage_generate_results();

    // Stage 3
    fprintf(stderr, "****TenX stage 3\n");
    aligner->align_third_stage();

    // Stage 4, process single mapping for unmapped reads 
    fprintf(stderr, "****TenX stage 4\n");
    if(!noSingle) {
        barcodeFinished = false;
        while (true) {
            barcodeFinished = aligner->align_forth_stage();
            if (barcodeFinished)
                break;
            for (unsigned pairIdx = 0; pairIdx < totalPairsForBarcode; pairIdx++) {
                if (tenXSingleTrackerArray[pairIdx].singleNotDone && tenXSingleTrackerArray[pairIdx].nSingleEndSecondaryResults[0] > tenXSingleTrackerArray[pairIdx].singleSecondaryBufferSize) {
                    _ASSERT(tenXSingleTrackerArray[pairIdx].nSingleEndSecondaryResults[0] > tenXSingleTrackerArray[pairIdx].singleSecondaryBufferSize);
                    BigDealloc(tenXSingleTrackerArray[pairIdx].singleEndSecondaryResults);
                    tenXSingleTrackerArray[pairIdx].singleEndSecondaryResults = NULL;
                    tenXSingleTrackerArray[pairIdx].singleSecondaryBufferSize *= 2;
                    tenXSingleTrackerArray[pairIdx].singleEndSecondaryResults = (SingleAlignmentResult *)BigAlloc(tenXSingleTrackerArray[pairIdx].singleSecondaryBufferSize * sizeof(SingleAlignmentResult));
                }
            }
        }
    }

    fprintf(stderr, "****generate output\n");

    /*
     * Output the results
     */
    
    // Record time
    _int64 alignFinishedTime;
    if (options->profile) {
        alignFinishedTime = timeInMillis();
        stats->millisAligning += (alignFinishedTime - readFinishedTime);
    }

    for (unsigned pairIdx = 0; pairIdx < totalPairsForBarcode; pairIdx++) {
        // A bunch of pointers delegates... For readability!
        Read **reads = tenXSingleTrackerArray[pairIdx].pairedReads;
        bool *useful = tenXSingleTrackerArray[pairIdx].useful;
        PairedAlignmentResult* results = tenXSingleTrackerArray[pairIdx].results;
        _int64 &nSecondaryResults = tenXSingleTrackerArray[pairIdx].nSecondaryResults;
        _int64 *nSingleEndSecondaryResults = tenXSingleTrackerArray[pairIdx].nSingleEndSecondaryResults;


#if     TIME_HISTOGRAM
        _int64 runTime = timeInNanos() - startTime;
        int timeBucket = min(30, cheezyLogBase2(runTime));
        stats->countByTimeBucket[timeBucket]++;
        stats->nanosByTimeBucket[timeBucket] += runTime;
#endif // TIME_HISTOGRAM

        if (forceSpacing && isOneLocation(results[0].status[0]) != isOneLocation(results[0].status[1])) {
            // either both align or neither do
            results[0].status[0] = results[0].status[1] = NotFound;
            results[0].location[0] = results[0].location[1] = InvalidGenomeLocation;
        }

        bool firstIsPrimary = true;
        for (int i = 0; i <= tenXSingleTrackerArray[pairIdx].secondaryResultBufferSize; i++) {  // Loop runs to <= nSecondaryResults because there's a primary result, too.
            bool pass0 = options->passFilter(reads[0], results[i].status[0], !useful[0], i != 0 || !firstIsPrimary);
            bool pass1 = options->passFilter(reads[1], results[i].status[1], !useful[1], i != 0 || !firstIsPrimary);
            bool pass = (options->filterFlags & AlignerOptions::FilterBothMatesMatch)
                ? (pass0 && pass1) : (pass0 || pass1);

            if (!pass) {
                //
                // Remove this one from the list by copying the last one here.
                //
                results[i] = results[nSecondaryResults];
                nSecondaryResults--;
                if (0 == i) {
                    firstIsPrimary = false;
                }
                i--;
            }
        }

        //
        // Now check the single secondary alignments
        //
        SingleAlignmentResult *singleResults[2] = { &tenXSingleTrackerArray[pairIdx].singleEndSecondaryResults[0], &tenXSingleTrackerArray[pairIdx].singleEndSecondaryResults[nSingleEndSecondaryResults[0]] };
        for (int whichRead = 0; whichRead < NUM_READS_PER_PAIR; whichRead++) {
            for (int whichAlignment = 0; whichAlignment < nSingleEndSecondaryResults[whichRead]; whichAlignment++) {
                if (!options->passFilter(reads[whichRead], singleResults[whichRead][whichAlignment].status, false, true)) {
                    singleResults[whichRead][whichAlignment] = singleResults[whichRead][nSingleEndSecondaryResults[whichRead] - 1];
                    nSingleEndSecondaryResults[whichRead]--;
                    whichAlignment--;
                }
            }
        }

        if (NULL != readWriter) {
            readWriter->writePairs(readerContext, reads, results, nSecondaryResults + 1, singleResults, nSingleEndSecondaryResults, firstIsPrimary);
        }

        // ****Not sure about all these stats. It's a legacy from the normal pairEndMapper. But now it's cluster based so it doesn't seem right no more. Whowever still wants meaningful stats from this, you need to fix this.
        if (options->profile) {
            startTime = timeInMillis();
            stats->millisWriting += (startTime - alignFinishedTime);
        }

        stats->extraAlignments += nSecondaryResults + (firstIsPrimary ? 0 : 1); // If first isn't primary, it's secondary.

        if (firstIsPrimary) {
            updateStats((TenXAlignerStats*)stats, reads[0], reads[1], results, useful[0], useful[1]);
        }
        else {
            stats->filtered += 2;
        }
    }   // while we have a read pair

    stats->lvCalls = aligner->getLocationsScored();

    //fprintf(stderr, "****begin cleanup\n");

    /*
     * Deallocate and clean up
     */

    allocator->checkCanaries();

    BigDealloc(tenXClusterIsValid);
    BigDealloc(tenXClusterPairsCounter);

    for (unsigned pairIdx = 0; pairIdx < maxBarcodeSize; pairIdx++) {
        BigDealloc(tenXSingleTrackerArray[pairIdx].results);
        tenXSingleTrackerArray[pairIdx].results = NULL;

        BigDealloc(tenXSingleTrackerArray[pairIdx].singleEndSecondaryResults);
        tenXSingleTrackerArray[pairIdx].singleEndSecondaryResults = NULL;

        BigDealloc(tenXSingleTrackerArray[pairIdx].clusterToggle);
        tenXSingleTrackerArray[pairIdx].clusterToggle = NULL;
    }

    for (unsigned singleAlignerIdx = 0; singleAlignerIdx < maxBarcodeSize; singleAlignerIdx++) {
        tenXSingleTrackerArray[singleAlignerIdx].aligner->~TenXSingleAligner();
    }
    
    aligner->~TenXClusterAligner();
    
    BigDealloc(tenXSingleTrackerArray);
    BigDealloc(sharedClusterCounterAry);
    BigDealloc(clusterCounterEraser);
    BigDealloc(clusterToggleEraser);
    
    delete supplier;

    delete allocator;
}


void TenXAlignerContext::updateStats(TenXAlignerStats* stats, Read* read0, Read* read1, PairedAlignmentResult* result, bool useful0, bool useful1)
{
    bool useful[2] = { useful0, useful1 };

    // Update stats
    for (int r = 0; r < 2; r++) {
        if (useful[r]) {
            if (isOneLocation(result->status[r])) {
                stats->singleHits++;
            }
            else if (result->status[r] == MultipleHits) {
                stats->multiHits++;
            }
            else {
                _ASSERT(result->status[r] == NotFound);
                stats->notFound++;
            }
            // Add in MAPQ stats
            if (result->status[r] != NotFound) {
                int mapq = result->mapq[r];
                _ASSERT(mapq >= 0 && mapq <= AlignerStats::maxMapq);
                stats->mapqHistogram[mapq]++;
            }
        }
        else {
            stats->uselessReads++;
        }

    }

    if (result->direction[0] == result->direction[1]) {
        stats->sameComplement++;
    }

    if (isOneLocation(result->status[0]) && isOneLocation(result->status[1])) {
        stats->incrementDistance(abs((int)(result->location[0] - result->location[1])));
        stats->incrementScore(result->score[0], result->score[1]);
    }

    if (result->fromAlignTogether) {
        stats->recordAlignTogetherMapqAndTime(__max(result->mapq[0], result->mapq[1]), result->nanosInAlignTogether, result->nSmallHits, result->nLVCalls);
    }

    if (result->alignedAsPair) {
        stats->alignedAsPairs += 2; // They are a pair, after all.  Hence, +2.
    }
}

void
TenXAlignerContext::typeSpecificBeginIteration()
{
    if (1 == options->nInputs) {
        //
        // We've only got one input, so just connect it directly to the consumer.
        //
        pairedReadSupplierGenerator = options->inputs[0].createPairedReadSupplierGenerator(options->numThreads, quicklyDropUnpairedReads, readerContext);
    }
    else {
        //
        // We've got multiple inputs, so use a MultiInputReadSupplier to combine the individual inputs.
        //
        PairedReadSupplierGenerator **generators = new PairedReadSupplierGenerator *[options->nInputs];
        printf("PairedReadSupplierGenerator\n");
        // use separate context for each supplier, initialized from common
        for (int i = 0; i < options->nInputs; i++) {
            ReaderContext context(readerContext);
            generators[i] = options->inputs[i].createPairedReadSupplierGenerator(options->numThreads, quicklyDropUnpairedReads, context);
        }
        pairedReadSupplierGenerator = new MultiInputPairedReadSupplierGenerator(options->nInputs, generators);
        printf("MultiInputPairedReadSupplierGenerator\n");
    }
    ReaderContext* context = pairedReadSupplierGenerator->getContext();
    readerContext.header = context->header;
    readerContext.headerBytes = context->headerBytes;
    readerContext.headerLength = context->headerLength;
    readerContext.headerMatchesIndex = context->headerMatchesIndex;
}
void
TenXAlignerContext::typeSpecificNextIteration()
{
    if (readerContext.header != NULL) {
        delete[] readerContext.header;
        readerContext.header = NULL;
        readerContext.headerLength = readerContext.headerBytes = 0;
        readerContext.headerMatchesIndex = false;
    }
    delete pairedReadSupplierGenerator;
    pairedReadSupplierGenerator = NULL;
}
