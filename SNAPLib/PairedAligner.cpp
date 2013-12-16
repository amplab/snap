/*++

Module Name:

    PairedAligner.cpp

Abstract:

    Functions for running the paired end aligner sub-program.


Authors:

    Matei Zaharia, February, 2012

Environment:
`
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
#include "Range.h"
#include "SAM.h"
#include "ChimericPairedEndAligner.h"
#include "Tables.h"
#include "WGsim.h"
#include "AlignerOptions.h"
#include "AlignerContext.h"
#include "AlignerStats.h"
#include "FASTQ.h"
#include "PairedAligner.h"
#include "MultiInputReadSupplier.h"
#include "Util.h"
#include "IntersectingPairedEndAligner.h"
#include "exit.h"
#include "AlignmentFilter.h"

using namespace std;

using util::stringEndsWith;

static const int DEFAULT_MIN_SPACING = 50;
static const int DEFAULT_MAX_SPACING = 1000;

struct PairedAlignerStats : public AlignerStats
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

    _int64 alignTogetherByMapqHistogram[maxMapq+1][nTimeBuckets];
    _int64 totalTimeByMapqHistogram[maxMapq+1][nTimeBuckets];
    _int64 nSmallHitsByTimeHistogram[nHitsBuckets][nTimeBuckets];
    _int64 nLVCallsByTimeHistogram[nLVCallsBuckets][nTimeBuckets];
    _int64 mapqByNLVCallsHistogram[maxMapq+1][nLVCallsBuckets];
    _int64 mapqByNSmallHitsHistogram[maxMapq+1][nHitsBuckets];

    PairedAlignerStats(AbstractStats* i_extra = NULL);

    virtual ~PairedAlignerStats();

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
        scoreCounts[s0*(MAX_SCORE+1)+s1]++;
    }

    inline void recordAlignTogetherMapqAndTime(unsigned mapq, _int64 timeInNanos, unsigned nSmallHits, unsigned nLVCalls) {
        int timeBucket;
        _int64 dividedTime = timeInNanos;
        for (timeBucket = 0; timeBucket < nTimeBuckets-1; timeBucket++) {
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

const int PairedAlignerStats::MAX_DISTANCE;
const int PairedAlignerStats::MAX_SCORE;

PairedAlignerStats::PairedAlignerStats(AbstractStats* i_extra)
    : AlignerStats(i_extra),
    sameComplement(0)
{
    int dsize = sizeof(_int64) * (MAX_DISTANCE+1);
    distanceCounts = (_int64*)BigAlloc(dsize);
    memset(distanceCounts, 0, dsize);

    int ssize = sizeof(_int64) * (MAX_SCORE+1)*(MAX_SCORE+1);
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

PairedAlignerStats::~PairedAlignerStats()
{
    BigDealloc(distanceCounts);
    BigDealloc(scoreCounts);
}

void PairedAlignerStats::add(const AbstractStats * i_other)
{
    AlignerStats::add(i_other);
    PairedAlignerStats* other = (PairedAlignerStats*) i_other;
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

void PairedAlignerStats::printHistograms(FILE* output)
{
    AlignerStats::printHistograms(output);
}

PairedAlignerOptions::PairedAlignerOptions(const char* i_commandLine)
    : AlignerOptions(i_commandLine, true),
    minSpacing(DEFAULT_MIN_SPACING),
    maxSpacing(DEFAULT_MAX_SPACING),
    forceSpacing(false),
    intersectingAlignerMaxHits(DEFAULT_INTERSECTING_ALIGNER_MAX_HITS),
    maxCandidatePoolSize(DEFAULT_MAX_CANDIDATE_POOL_SIZE)
{
}

void PairedAlignerOptions::usageMessage()
{
    AlignerOptions::usageMessage();
    printf(
        "  -s   min and max spacing to allow between paired ends (default: %d %d)\n"
        "  -fs  force spacing to lie between min and max\n"
        "  -H   max hits for intersecting aligner (default: %d)\n"
        "  -mcp specifies the maximum candidate pool size (An internal data structure. \n"
        "       Only increase this if you get an error message saying to do so. If you're running\n"
        "       out of memory, you may want to reduce it.  Default: %d)\n",
        DEFAULT_MIN_SPACING,
        DEFAULT_MAX_SPACING,
        DEFAULT_INTERSECTING_ALIGNER_MAX_HITS,
        DEFAULT_MAX_CANDIDATE_POOL_SIZE);
}

bool PairedAlignerOptions::parse(const char** argv, int argc, int& n, bool *done)
{
    *done = false;

    if (strcmp(argv[n], "-s") == 0) {
        if (n + 2 < argc) {
            minSpacing = atoi(argv[n+1]);
            maxSpacing = atoi(argv[n+2]);
            n += 2;
            return true;
        } 
        return false;
    } else if (strcmp(argv[n], "-H") == 0) {
        if (n + 1 < argc) {
            intersectingAlignerMaxHits = atoi(argv[n+1]);
            n += 1;
            return true;
        } 
        return false;
    } else if (strcmp(argv[n], "-fs") == 0) {
        forceSpacing = true;
        return true;
    } else if (strcmp(argv[n], "-mcp") == 0) {
        if (n + 1 < argc) {
            maxCandidatePoolSize = atoi(argv[n+1]);
            n += 1;
            return true;
        } 
        return false;
    }
    return AlignerOptions::parse(argv, argc, n, done);
}

PairedAlignerContext::PairedAlignerContext(AlignerExtension* i_extension)
    : AlignerContext( 0,  NULL, NULL, i_extension)
{
}

AlignerOptions* PairedAlignerContext::parseOptions(int i_argc, const char **i_argv, const char *i_version, unsigned *argsConsumed)
{
    argc = i_argc;
    argv = i_argv;
    version = i_version;

    PairedAlignerOptions* options = new PairedAlignerOptions(
        "snapr paired <genome-dir> <transcriptome-dir> <annotation> <input file(s)> <read2.fq> [<options>]\n"
        "   where <input file(s)> is a list of files to process.  FASTQ\n"
        "   files must come in pairs, since each read end is in a separate file.");
    options->extra = extension->extraOptions();
    if (argc < 4) {
        fprintf(stderr, "Too few parameters\n");
        options->usage();
    }

    options->indexDir = argv[0];
    options->transcriptomeDir = argv[1];
    options->annotation = argv[2];
    //
    // Figure out how many inputs there are.  All options begin with a '-', so count the
    // args until we hit an option.  FASTQ files come in pairs, and each pair only counts
    // as one input.
    //
    int nInputs = 0;
    bool foundFirstHalfOfFASTQ = false;
    for (int i = 3; i < argc; i++) {
        if (argv[i][0] == '-' || argv[i][0] == ',' && argv[i][1] == '\0') {
                break;
        }
        
        if (stringEndsWith(argv[i],".sam") || stringEndsWith(argv[i],".bam")) {
            if (foundFirstHalfOfFASTQ) {
                fprintf(stderr,"For the paired aligner, FASTQ files must come in pairs.  I found SAM/BAM file '%s' after first half FASTQ file '%s'.\n",
                    argv[i],argv[i-1]);
                soft_exit(1);
            }
            nInputs++;
        } else {
            if (foundFirstHalfOfFASTQ) {
                nInputs++;
            }
            foundFirstHalfOfFASTQ = !foundFirstHalfOfFASTQ;
        }
    }
    if (foundFirstHalfOfFASTQ) {
        fprintf(stderr,"For the paired aligner, FASTQ files must come in pairs.  The last one is unmatched.\n");
        soft_exit(1);
    }
    if (0 == nInputs) {
        fprintf(stderr,"Didn't see any input files\n");
        options->usage();
    }
    //
    // Now build the input array.
    //
    options->nInputs = nInputs;
    options->inputs = new SNAPInput[nInputs];
    int i;
    int whichInput = 0;
    for (i = 3; i < argc; i++) {
        if (argv[i][0] == '-') {
                break;
        }

        if (stringEndsWith(argv[i],".sam") || stringEndsWith(argv[i],".bam")) {
            _ASSERT(!foundFirstHalfOfFASTQ);
            options->inputs[whichInput].fileType = stringEndsWith(argv[i],".sam") ? SAMFile : BAMFile;
            options->inputs[whichInput].fileName = argv[i];
            whichInput++;
        } else {
            if (foundFirstHalfOfFASTQ) {
                options->inputs[whichInput].fileType =
                    stringEndsWith(argv[i],".gzip") || stringEndsWith(argv[i], ".gz") // todo: more suffixes?
                        ? GZipFASTQFile : FASTQFile;
                options->inputs[whichInput].secondFileName = argv[i];
                whichInput++;
            } else {
                options->inputs[whichInput].fileName = argv[i];
            }
            foundFirstHalfOfFASTQ = !foundFirstHalfOfFASTQ;
        }
    }

    for (/* i initialized by previous loop*/; i < argc; i++) {
        bool done;
        int oldI = i;
        if (!options->parse(argv, argc, i, &done)) {
            fprintf(stderr, "Didn't understand options starting at %s\n", argv[oldI]);
            options->usage();
        }

        if (done) {
            i++;    // For the ',' arg
            break;
        }
    }

    if (options->maxDist.end + options->extraSearchDepth >= MAX_K) {
        fprintf(stderr,"You specified too large of a maximum edit distance combined with extra search depth.  The must add up to less than %d.\n", MAX_K);
        fprintf(stderr,"Either reduce their sum, or change MAX_K in LandauVishkin.h and recompile.\n");
        soft_exit(1);
    }
        
    *argsConsumed = i;
    return options;
}

void PairedAlignerContext::initialize()
{
    AlignerContext::initialize();
    PairedAlignerOptions* options2 = (PairedAlignerOptions*) options;
    minSpacing = options2->minSpacing;
    maxSpacing = options2->maxSpacing;
    forceSpacing = options2->forceSpacing;
    maxCandidatePoolSize = options2->maxCandidatePoolSize;
    intersectingAlignerMaxHits = options2->intersectingAlignerMaxHits;
    ignoreMismatchedIDs = options2->ignoreMismatchedIDs;
}

AlignerStats* PairedAlignerContext::newStats()
{
    return new PairedAlignerStats();
}

void PairedAlignerContext::runTask()
{
    ParallelTask<PairedAlignerContext> task(this);
    task.run();
}

void PairedAlignerContext::runIterationThread()
{
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
	if (index == NULL) {
        // no alignment, just input/output
        Read *read0;
        Read *read1;
        PairedAlignmentResult result;
        memset(&result, 0, sizeof(result));
        result.location[0] = result.location[1] = InvalidGenomeLocation;
        
        while (supplier->getNextReadPair(&read0,&read1)) {
            // Check that the two IDs form a pair; they will usually be foo/1 and foo/2 for some foo.
            if (!ignoreMismatchedIDs && !readIdsMatch(read0, read1)) {
                unsigned n[2] = {min(read0->getIdLength(), 200u), min(read1->getIdLength(), 200u)};
                char* p[2] = {(char*) alloca(n[0] + 1), (char*) alloca(n[1] + 1)};
                memcpy(p[0], read0->getId(), n[0]); p[0][n[0]] = 0;
                memcpy(p[1], read1->getId(), n[1]); p[1][n[1]] = 0;
                fprintf(stderr, "Unmatched read IDs '%s' and '%s'.  Use the -I option to ignore this.\n", p[0], p[1]);
                soft_exit(1);
            }
            stats->totalReads += 2;
            writePair(read0, read1, &result);
        }
        delete supplier;
        return;
    }

    int maxReadSize = MAX_READ_LENGTH;
    size_t g_memoryPoolSize = IntersectingPairedEndAligner::getBigAllocatorReservation(index, intersectingAlignerMaxHits, maxReadSize, index->getSeedLength(), 
                                                                numSeedsFromCommandLine, seedCoverage, maxDist, extraSearchDepth, maxCandidatePoolSize);

    g_memoryPoolSize += ChimericPairedEndAligner::getBigAllocatorReservation(index, maxReadSize, maxHits, index->getSeedLength(), numSeedsFromCommandLine, seedCoverage, maxDist,
                                                    extraSearchDepth, maxCandidatePoolSize);

    BigAllocator *g_allocator = new BigAllocator(g_memoryPoolSize);
    
    IntersectingPairedEndAligner *g_intersectingAligner = new (g_allocator) IntersectingPairedEndAligner(index, maxReadSize, maxHits, maxDist, numSeedsFromCommandLine, 
                                                                seedCoverage, minSpacing, maxSpacing, intersectingAlignerMaxHits, extraSearchDepth, 
                                                                maxCandidatePoolSize, g_allocator);


    ChimericPairedEndAligner *g_aligner = new (g_allocator) ChimericPairedEndAligner(
        index,
        maxReadSize,
        maxHits,
        maxDist,
        numSeedsFromCommandLine,
        seedCoverage,
        forceSpacing,
        extraSearchDepth,
        g_intersectingAligner,
        g_allocator);

    g_allocator->checkCanaries();

    BigAllocator *c_allocator = NULL;
    IntersectingPairedEndAligner *c_intersectingAligner = NULL;
    ChimericPairedEndAligner *c_aligner = NULL;

    if (contamination != NULL) {

      //Contamination database for paired end reads
      size_t c_memoryPoolSize = IntersectingPairedEndAligner::getBigAllocatorReservation(contamination, intersectingAlignerMaxHits, maxReadSize, contamination->getSeedLength(),
                                                                  numSeedsFromCommandLine, seedCoverage, maxDist, extraSearchDepth, maxCandidatePoolSize);

      c_memoryPoolSize += ChimericPairedEndAligner::getBigAllocatorReservation(contamination, maxReadSize, maxHits, contamination->getSeedLength(), numSeedsFromCommandLine, seedCoverage, maxDist, extraSearchDepth, maxCandidatePoolSize);

      c_allocator = new BigAllocator(c_memoryPoolSize);

      
      c_intersectingAligner = new (c_allocator) IntersectingPairedEndAligner(contamination, maxReadSize, maxHits, maxDist, numSeedsFromCommandLine,
                                                                  seedCoverage, minSpacing, maxSpacing, intersectingAlignerMaxHits, extraSearchDepth,
                                                                  maxCandidatePoolSize, c_allocator);

      c_aligner = new (c_allocator) ChimericPairedEndAligner(
          contamination,
          maxReadSize,
          maxHits, 
          maxDist,
          numSeedsFromCommandLine,
          seedCoverage,
          forceSpacing,
          extraSearchDepth,
          c_intersectingAligner,
          c_allocator);

      c_allocator->checkCanaries();
     
    }

    unsigned singleAlignerMaxHits = 300;
    unsigned p_numSeedsFromCommandLine = 0;
    float p_seedCoverage = maxReadSize / (index->getSeedLength()*2);
    BigAllocator *t_allocator = new BigAllocator(BaseAligner::getBigAllocatorReservation(true, singleAlignerMaxHits, maxReadSize, transcriptome->getSeedLength(), numSeedsFromCommandLine, seedCoverage));

    BaseAligner *t_aligner = new (t_allocator) BaseAligner(
            transcriptome,
            singleAlignerMaxHits,
            maxDist,
            maxReadSize,
            numSeedsFromCommandLine,            
            seedCoverage,
            extraSearchDepth,
            NULL,               // LV (no need to cache in the single aligner)
            NULL,               // reverse LV
            stats,
            t_allocator);

    t_allocator->checkCanaries();
    t_aligner->setExplorePopularSeeds(options->explorePopularSeeds);
    t_aligner->setStopOnFirstHit(options->stopOnFirstHit);

    BigAllocator *p_allocator = new BigAllocator(BaseAligner::getBigAllocatorReservation(true, singleAlignerMaxHits, maxReadSize, index->getSeedLength(), p_numSeedsFromCommandLine, p_seedCoverage));

    BaseAligner *p_aligner = new (p_allocator) BaseAligner(
            index,
            singleAlignerMaxHits,
            maxDist,
            maxReadSize,
            p_numSeedsFromCommandLine,
            p_seedCoverage, //seedCoverage = 0 to ensure seeds across entire read for partial alignments
            extraSearchDepth,
            NULL,               // LV (no need to cache in the single aligner)
            NULL,               // reverse LV
            stats,
            p_allocator);

    p_allocator->checkCanaries();
    p_aligner->setExplorePopularSeeds(options->explorePopularSeeds);
    p_aligner->setStopOnFirstHit(options->stopOnFirstHit);

    ReadWriter *readWriter = this->readWriter;

#ifdef  _MSC_VER
    if (options->useTimingBarrier) {
        if (0 == InterlockedDecrementAndReturnNewValue(nThreadsAllocatingMemory)) {
            AllowEventWaitersToProceed(memoryAllocationCompleteBarrier);
        } else {
            WaitForEvent(memoryAllocationCompleteBarrier);
        }
    }
#endif  // _MSC_VER

    // Align the reads.
    Read *read0;
    Read *read1;
    while (supplier->getNextReadPair(&read0,&read1)) {
        // Check that the two IDs form a pair; they will usually be foo/1 and foo/2 for some foo.
        if (!ignoreMismatchedIDs) {
			Read::checkIdMatch(read0, read1);
		}

        stats->totalReads += 2;

        // Skip the pair if there are too many Ns or 2s.
        int maxDist = this->maxDist;
        bool useful0 = read0->getDataLength() >= 50 && (int)read0->countOfNs() <= maxDist;
        bool useful1 = read1->getDataLength() >= 50 && (int)read1->countOfNs() <= maxDist;

        //Quality filtering
        bool quality0 = read0->qualityFilter(options->minPercentAbovePhred, options->minPhred, options->phredOffset);
        bool quality1 = read1->qualityFilter(options->minPercentAbovePhred, options->minPhred, options->phredOffset);
        
        if ((!useful0 && !useful1) || (!quality0 || !quality0)) {
            PairedAlignmentResult result;
            result.isTranscriptome[0] = false;
            result.isTranscriptome[1] = false;
            result.status[0] = NotFound;
            result.status[1] = NotFound;
            result.location[0] = InvalidGenomeLocation;
            result.location[1] = InvalidGenomeLocation;
            writePair(read0, read1, &result);
            continue;
        } else {
            // Here one the reads might still be hopeless, but maybe we can align the other.
            stats->usefulReads += (useful0 && useful1) ? 2 : 1;
        }

        PairedAlignmentResult result, contaminantResult;
        result.isTranscriptome[0] = false;
        result.isTranscriptome[1] = false;

        //Make users setting
        AlignmentFilter filter(read0, read1, index->getGenome(), transcriptome->getGenome(), gtf, minSpacing, maxSpacing, options->confDiff, options->maxDist.start, index->getSeedLength(), p_aligner);

        const unsigned maxHitsToGet = 1000;
        unsigned loc0, loc1;
        Direction rc0, rc1, temp;
        int score0, score1;
        int mapq0, mapq1;
        AlignmentResult status0, status1;
        
        int       transcriptome_multiHitsFound0;
        unsigned  transcriptome_multiHitLocations0[maxHitsToGet];
        bool      transcriptome_multiHitRCs0[maxHitsToGet];
        int       transcriptome_multiHitScores0[maxHitsToGet]; 
        int       transcriptome_multiHitsFound1;
        unsigned  transcriptome_multiHitLocations1[maxHitsToGet];
        bool      transcriptome_multiHitRCs1[maxHitsToGet];
        int       transcriptome_multiHitScores1[maxHitsToGet]; 
        
        t_aligner->setReadId(0);      
        status0 = t_aligner->AlignRead(read0, &loc0, &rc0, &score0, &mapq0, 0, 0, 0, maxHitsToGet, &transcriptome_multiHitsFound0, (unsigned*)&transcriptome_multiHitLocations0, (bool*)&transcriptome_multiHitRCs0, (int*)&transcriptome_multiHitScores0);            
       
        t_aligner->setReadId(1);      
        status0 = t_aligner->AlignRead(read1, &loc1, &rc1, &score1, &mapq1, 0, 0, 0, maxHitsToGet, &transcriptome_multiHitsFound1, (unsigned*)&transcriptome_multiHitLocations1, (bool*)&transcriptome_multiHitRCs1, (int*)&transcriptome_multiHitScores1);            
 
        //Add reads to filter
        for (int i = 0; i < transcriptome_multiHitsFound0; ++i) {
            filter.AddAlignment(transcriptome_multiHitLocations0[i], transcriptome_multiHitRCs0[i], transcriptome_multiHitScores0[i], 0, true, false);
        }

        for (int i = 0; i < transcriptome_multiHitsFound1; ++i) {
            filter.AddAlignment(transcriptome_multiHitLocations1[i], transcriptome_multiHitRCs1[i], transcriptome_multiHitScores1[i], 0, true, true);
        }    
    
        g_aligner->align(read0, read1, &result);
        filter.AddAlignment(result.location[0], result.direction[0], result.score[0], result.mapq[0], false, false);
        filter.AddAlignment(result.location[1], result.direction[1], result.score[1], result.mapq[1], false, true);
     
        //Perform the primary filtering of all aligned reads
        unsigned status = filter.Filter(&result);

        //If the read is still unaligned
        if ((result.status[0] == NotFound) && (result.status[1] == NotFound)) {
        
          //If the contamination database is present
          if (c_aligner != NULL) {

            c_aligner->align(read0, read1, &contaminantResult);
            if ((contaminantResult.status[0] != NotFound) && (contaminantResult.status[1] != NotFound)) {

              c_filter->AddAlignment(contaminantResult.location[0], string(read0->getId(), read0->getIdLength()), string(read0->getData(), read0->getDataLength()));
              c_filter->AddAlignment(contaminantResult.location[1], string(read1->getId(), read1->getIdLength()), string(read1->getData(), read1->getDataLength()));

            }
          }
        }

        if (forceSpacing && isOneLocation(result.status[0]) != isOneLocation(result.status[1])) {
            // either both align or neither do
            result.status[0] = result.status[1] = NotFound;
            result.location[0] = result.location[1] = InvalidGenomeLocation;
        }
#if 0       // cheese
        if (result.score[0] + result.score[1] >= 5) {
            double divisor = __max(1,((result.score[0] + result.score[1]) *2.0) / 5.0);
            if (result.mapq[0] < 50) {
                result.mapq[0] /= 2;
            }
            if (result.mapq[1] < 50) {
                result.mapq[1] /= 2;
            }
        }
#endif // 0

        writePair(read0, read1, &result);

        updateStats((PairedAlignerStats*) stats, read0, read1, &result);
    }

    stats->lvCalls = g_aligner->getLocationsScored();

    g_allocator->checkCanaries();

    g_aligner->~ChimericPairedEndAligner();
    delete supplier;

    g_intersectingAligner->~IntersectingPairedEndAligner();
    delete g_allocator;

    t_aligner->~BaseAligner();
    delete t_allocator;

    p_aligner->~BaseAligner();
    delete p_allocator;

    if (c_allocator != NULL) {
        c_allocator->checkCanaries();
        c_aligner->~ChimericPairedEndAligner();
        c_intersectingAligner->~IntersectingPairedEndAligner();
        delete c_allocator;
    }

}

void PairedAlignerContext::writePair(Read* read0, Read* read1, PairedAlignmentResult* result)
{
    if (readWriter != NULL && (options->passFilter(read0, result->status[0]) || options->passFilter(read1, result->status[1]))) {
        readWriter->writePair(read0, read1, result);
    }
}

void PairedAlignerContext::updateStats(PairedAlignerStats* stats, Read* read0, Read* read1, PairedAlignmentResult* result)
{
    // Update stats
    for (int r = 0; r < 2; r++) {
        bool wasError = false;
        if (computeError && result->status[r] != NotFound) {
            wasError = wgsimReadMisaligned((r == 0 ? read0 : read1), result->location[r], index, options->misalignThreshold);
        }
        if (isOneLocation(result->status[r])) {
            stats->singleHits++;
            stats->errors += wasError ? 1 : 0;
        } else if (result->status[r] == MultipleHits) {
            stats->multiHits++;
        } else {
            _ASSERT(result->status[r] == NotFound);
            stats->notFound++;
        }
        // Add in MAPQ stats
        if (result->status[r] != NotFound) {
            int mapq = result->mapq[r];
            _ASSERT(mapq >= 0 && mapq <= AlignerStats::maxMapq);
            stats->mapqHistogram[mapq]++;
            stats->mapqErrors[mapq] += wasError ? 1 : 0;
        }
    }

    if (result->direction[0] == result->direction[1]) {
        stats->sameComplement++;
    }

    if (isOneLocation(result->status[0]) && isOneLocation(result->status[1])) {
        stats->incrementDistance(abs((int) (result->location[0] - result->location[1])));
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
PairedAlignerContext::typeSpecificBeginIteration()
{
    if (1 == options->nInputs) {
        //
        // We've only got one input, so just connect it directly to the consumer.
        //
        options->inputs[0].readHeader(readerContext);
        pairedReadSupplierGenerator = options->inputs[0].createPairedReadSupplierGenerator(options->numThreads, readerContext);
    } else {
        //
        // We've got multiple inputs, so use a MultiInputReadSupplier to combine the individual inputs.
        //
        PairedReadSupplierGenerator **generators = new PairedReadSupplierGenerator *[options->nInputs];
        // use separate context for each supplier, initialized from common
        for (int i = 0; i < options->nInputs; i++) {
            ReaderContext context(readerContext);
            options->inputs[i].readHeader(context);
            generators[i] = options->inputs[i].createPairedReadSupplierGenerator(options->numThreads, context);
        }
        pairedReadSupplierGenerator = new MultiInputPairedReadSupplierGenerator(options->nInputs,generators);
    }
}
    void 
PairedAlignerContext::typeSpecificNextIteration()
{
    delete pairedReadSupplierGenerator;
    pairedReadSupplierGenerator = NULL;
}
