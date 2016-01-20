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
#include "SAM.h"
#include "ChimericPairedEndAligner.h"
#include "Tables.h"
#include "AlignerOptions.h"
#include "AlignerContext.h"
#include "AlignerStats.h"
#include "FASTQ.h"
#include "PairedAligner.h"
#include "MultiInputReadSupplier.h"
#include "Util.h"
#include "IntersectingPairedEndAligner.h"
#include "exit.h"
#include "Error.h"

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
    maxCandidatePoolSize(DEFAULT_MAX_CANDIDATE_POOL_SIZE),
    quicklyDropUnpairedReads(true)
{
}

void PairedAlignerOptions::usageMessage()
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
        ,
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
    } else if (strcmp(argv[n], "-ku") == 0) {
        quicklyDropUnpairedReads = false;
        return true;
    } else if (strcmp(argv[n], "-mcp") == 0) {
        if (n + 1 < argc) {
            maxCandidatePoolSize = atoi(argv[n+1]);
            n += 1;
            return true;
        } 
        return false;
    } else if (strcmp(argv[n], "-F") == 0 && n + 1 < argc && strcmp(argv[n + 1],"b") == 0) {
        filterFlags |= FilterBothMatesMatch;
        n += 1;
        return true;
    }
    return AlignerOptions::parse(argv, argc, n, done);
}

PairedAlignerContext::PairedAlignerContext(AlignerExtension* i_extension)
    : AlignerContext( 0,  NULL, NULL, i_extension)
{
}

bool PairedAlignerContext::initialize()
{
    AlignerContext::initialize();
    PairedAlignerOptions* options2 = (PairedAlignerOptions*) options;
    minSpacing = options2->minSpacing;
    maxSpacing = options2->maxSpacing;
    forceSpacing = options2->forceSpacing;
    maxCandidatePoolSize = options2->maxCandidatePoolSize;
    intersectingAlignerMaxHits = options2->intersectingAlignerMaxHits;
    ignoreMismatchedIDs = options2->ignoreMismatchedIDs;
    quicklyDropUnpairedReads = options2->quicklyDropUnpairedReads;
    noUkkonen = options->noUkkonen;
    noOrderedEvaluation = options->noOrderedEvaluation;

	return true;
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

    Read *reads[NUM_READS_PER_PAIR];
    int nSingleResults[2] = { 0, 0 };

 	if (index == NULL) {
        // no alignment, just input/output
        PairedAlignmentResult result;
        memset(&result, 0, sizeof(result));
        result.location[0] = result.location[1] = InvalidGenomeLocation;

         
        while (supplier->getNextReadPair(&reads[0],&reads[1])) {
            // Check that the two IDs form a pair; they will usually be foo/1 and foo/2 for some foo.
            if (!ignoreMismatchedIDs && !readIdsMatch(reads[0], reads[1])) {
                unsigned n[2] = {min(reads[0]->getIdLength(), 200u), min(reads[1]->getIdLength(), 200u)};
                char* p[2] = {(char*) alloca(n[0] + 1), (char*) alloca(n[1] + 1)};
                memcpy(p[0], reads[0]->getId(), n[0]); p[0][n[0]] = 0;
                memcpy(p[1], reads[1]->getId(), n[1]); p[1][n[1]] = 0;
                WriteErrorMessage( "Unmatched read IDs '%s' and '%s'.  Use the -I option to ignore this.\n", p[0], p[1]);
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
            } else {
                stats->uselessReads++;
            }
        }
        delete supplier;
        return;
    }

    int maxReadSize = MAX_READ_LENGTH;
    size_t memoryPoolSize = IntersectingPairedEndAligner::getBigAllocatorReservation(index, intersectingAlignerMaxHits, maxReadSize, index->getSeedLength(), 
                                                                numSeedsFromCommandLine, seedCoverage, maxDist, extraSearchDepth, maxCandidatePoolSize,
                                                                maxSecondaryAlignmentsPerContig);

    memoryPoolSize += ChimericPairedEndAligner::getBigAllocatorReservation(index, maxReadSize, maxHits, index->getSeedLength(), numSeedsFromCommandLine, seedCoverage, maxDist,
        extraSearchDepth, maxCandidatePoolSize, maxSecondaryAlignmentsPerContig);

    unsigned maxPairedSecondaryHits;
    unsigned maxSingleSecondaryHits;

    if (maxSecondaryAlignmentAdditionalEditDistance < 0) {
        maxPairedSecondaryHits = 0;
        maxSingleSecondaryHits = 0;
    } else {
        maxPairedSecondaryHits = IntersectingPairedEndAligner::getMaxSecondaryResults(numSeedsFromCommandLine, seedCoverage, maxReadSize, maxHits, index->getSeedLength(), minSpacing, maxSpacing);
        maxSingleSecondaryHits = ChimericPairedEndAligner::getMaxSingleEndSecondaryResults(numSeedsFromCommandLine, seedCoverage, maxReadSize, maxHits, index->getSeedLength());
    }

    memoryPoolSize += (1 + maxPairedSecondaryHits) * sizeof(PairedAlignmentResult) + maxSingleSecondaryHits * sizeof(SingleAlignmentResult);

    BigAllocator *allocator = new BigAllocator(memoryPoolSize);
    
    IntersectingPairedEndAligner *intersectingAligner = new (allocator) IntersectingPairedEndAligner(index, maxReadSize, maxHits, maxDist, numSeedsFromCommandLine, 
                                                                seedCoverage, minSpacing, maxSpacing, intersectingAlignerMaxHits, extraSearchDepth, 
                                                                maxCandidatePoolSize, maxSecondaryAlignmentsPerContig ,allocator, noUkkonen, noOrderedEvaluation, noTruncation);


    ChimericPairedEndAligner *aligner = new (allocator) ChimericPairedEndAligner(
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
        intersectingAligner,
		minReadLength,
        maxSecondaryAlignmentsPerContig,
        allocator);

    allocator->checkCanaries();

    PairedAlignmentResult *results = (PairedAlignmentResult *)allocator->allocate((1 + maxPairedSecondaryHits) * sizeof(*results)); // 1 + is for the primary result
    SingleAlignmentResult *singleSecondaryResults = (SingleAlignmentResult *)allocator->allocate(maxSingleSecondaryHits * sizeof(*singleSecondaryResults));

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
    _uint64 lastReportTime = timeInMillis();
    _uint64 readsWhenLastReported = 0;

    while (supplier->getNextReadPair(&reads[0], &reads[1])) {
        // Check that the two IDs form a pair; they will usually be foo/1 and foo/2 for some foo.
        if (!ignoreMismatchedIDs) {
            Read::checkIdMatch(reads[0], reads[1]);
        }

        stats->totalReads += 2;

        if (AlignerOptions::useHadoopErrorMessages && stats->totalReads % 10000 == 0 && timeInMillis() - lastReportTime > 10000) {
            fprintf(stderr, "reporter:counter:SNAP,readsAligned,%lu\n", stats->totalReads - readsWhenLastReported);
            readsWhenLastReported = stats->totalReads;
            lastReportTime = timeInMillis();
        }

        // Skip the pair if there are too many Ns and/or they're too short
        int maxDist = this->maxDist;
        bool useful0 = reads[0]->getDataLength() >= minReadLength && (int)reads[0]->countOfNs() <= maxDist;
        bool useful1 = reads[1]->getDataLength() >= minReadLength && (int)reads[1]->countOfNs() <= maxDist;
        if (!useful0 && !useful1) {
            PairedAlignmentResult result;
            result.status[0] = NotFound;
            result.status[1] = NotFound;
            result.location[0] = InvalidGenomeLocation;
            result.location[1] = InvalidGenomeLocation;
            nSingleResults[0] = nSingleResults[1] = 0;

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



#if     TIME_HISTOGRAM
        _int64 startTime = timeInNanos();
#endif // TIME_HISTOGRAM

        int nSecondaryResults;
        int nSingleSecondaryResults[2];

        aligner->align(reads[0], reads[1], results, maxSecondaryAlignmentAdditionalEditDistance, maxPairedSecondaryHits, &nSecondaryResults, results + 1,
            maxSingleSecondaryHits, maxSecondaryAlignments, &nSingleSecondaryResults[0], &nSingleSecondaryResults[1], singleSecondaryResults);

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
        for (int i = 0; i <= nSecondaryResults; i++) {  // Loop runs to <= nSecondaryResults because there's a primary result, too.
            bool pass0 = options->passFilter(reads[0], results[i].status[0], !useful0, i != 0 || !firstIsPrimary);
            bool pass1 = options->passFilter(reads[1], results[i].status[1], !useful1, i != 0 || !firstIsPrimary);
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
        SingleAlignmentResult *singleResults[2] = { singleSecondaryResults, singleSecondaryResults + nSingleSecondaryResults[0] };
        for (int whichRead = 0; whichRead < NUM_READS_PER_PAIR; whichRead++) {
            for (int whichAlignment = 0; whichAlignment < nSingleSecondaryResults[whichRead]; whichAlignment++) {
                if (!options->passFilter(reads[whichRead], singleResults[whichRead][whichAlignment].status, false, true)) {
                    singleResults[whichRead][whichAlignment] = singleResults[whichRead][nSingleSecondaryResults[whichRead] - 1];
                    nSingleSecondaryResults[whichRead]--;
                    whichAlignment--;
                }
            }
        }

        if (NULL != readWriter) {
            readWriter->writePairs(readerContext, reads, results, nSecondaryResults + 1, singleResults, nSingleSecondaryResults, firstIsPrimary);
        }

        stats->extraAlignments += nSecondaryResults + (firstIsPrimary ? 0 : 1); // If first isn't primary, it's secondary.
        if (firstIsPrimary) {
            updateStats((PairedAlignerStats*)stats, reads[0], reads[1], &results[0], useful0, useful1);
        } else {
            stats->filtered += 2;
        }
    }   // while we have a read pair

    stats->lvCalls = aligner->getLocationsScored();

    allocator->checkCanaries();

    aligner->~ChimericPairedEndAligner();
    delete supplier;

    intersectingAligner->~IntersectingPairedEndAligner();
    delete allocator;
}


void PairedAlignerContext::updateStats(PairedAlignerStats* stats, Read* read0, Read* read1, PairedAlignmentResult* result, bool useful0, bool useful1)
{
	bool useful[2] = { useful0, useful1 };

    // Update stats
    for (int r = 0; r < 2; r++) {
		if (useful[r]) {
			if (isOneLocation(result->status[r])) {
				stats->singleHits++;
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
            }
        } else {
            stats->uselessReads++;
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
        pairedReadSupplierGenerator = options->inputs[0].createPairedReadSupplierGenerator(options->numThreads, quicklyDropUnpairedReads, readerContext);
    } else {
        //
        // We've got multiple inputs, so use a MultiInputReadSupplier to combine the individual inputs.
        //
        PairedReadSupplierGenerator **generators = new PairedReadSupplierGenerator *[options->nInputs];
        // use separate context for each supplier, initialized from common
        for (int i = 0; i < options->nInputs; i++) {
            ReaderContext context(readerContext);
            generators[i] = options->inputs[i].createPairedReadSupplierGenerator(options->numThreads, quicklyDropUnpairedReads, context);
        }
        pairedReadSupplierGenerator = new MultiInputPairedReadSupplierGenerator(options->nInputs,generators);
    }
    ReaderContext* context = pairedReadSupplierGenerator->getContext();
    readerContext.header = context->header;
    readerContext.headerBytes = context->headerBytes;
    readerContext.headerLength = context->headerLength;
    readerContext.headerMatchesIndex = context->headerMatchesIndex;
}
    void 
PairedAlignerContext::typeSpecificNextIteration()
{
    if (readerContext.header != NULL) {
        delete [] readerContext.header;
        readerContext.header = NULL;
        readerContext.headerLength = readerContext.headerBytes = 0;
        readerContext.headerMatchesIndex = false;
    }
    delete pairedReadSupplierGenerator;
    pairedReadSupplierGenerator = NULL;
}
