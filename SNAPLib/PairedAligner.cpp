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
#include "SmarterPairedEndAligner.h"
#include "Tables.h"
#include "WGsim.h"
#include "GoodRandom.h"
#include "AlignerOptions.h"
#include "AlignerContext.h"
#include "AlignerStats.h"
#include "FASTQ.h"
#include "PairedAligner.h"
#include "ThirdPairedEndAligner.h"
#include "MultiInputReadSupplier.h"
#include "Util.h"
#include "BloomPairedEndAligner.h"
#include "IntersectingPairedEndAligner.h"
#include "exit.h"

using namespace std;

using util::stringEndsWith;

static const int DEFAULT_MIN_SPACING = 100;
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

    virtual void add(AbstractStats * other);

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
}

PairedAlignerStats::~PairedAlignerStats()
{
    BigDealloc(distanceCounts);
    BigDealloc(scoreCounts);
}

void PairedAlignerStats::add(AbstractStats * i_other)
{
    AlignerStats::add(i_other);
    PairedAlignerStats* other = (PairedAlignerStats*) i_other;
    for (int i = 0; i < MAX_DISTANCE + 1; i++) {
        distanceCounts[i] += other->distanceCounts[i];
    }
    for (int i = 0; i < (MAX_SCORE + 1) * (MAX_SCORE + 1); i++) {
        scoreCounts[i] += other->scoreCounts[i];
    }
}

void PairedAlignerStats::printHistograms(FILE* output)
{
    AlignerStats::printHistograms(output);
    // print all non-zeros
    fprintf(output, "\ndistance\tpairs\n");
    for (int i = 0; i <= MAX_DISTANCE; i++) {
        if (distanceCounts[i] != 0) {
            fprintf(output, "%d\t%lld\n", i, distanceCounts[i]);
        }
    }

    fprintf(output, "\nscores\n");
    int max1 = MAX_SCORE;
    bool found = false;
    while (max1 > 0) {
        for (int i = 0; i <= max1 && !found; i++) {
            found = scoreCounts[i * (MAX_SCORE + 1) + max1] != 0;
        }
        if (found) {
            break;
        }
        max1--;
    }
    for (int s1 = 0; s1 <= max1; s1++) {
        for (int s0 = 0; s0 <= s1; s0++) {
            fprintf(output, "%lld%s", scoreCounts[s0 * (MAX_SCORE + 1) + s1], s0 < s1 ? "\t" : "\n");
        }
    }
}

static bool readIdsMatch(Read *read0, Read *read1)
{
    if (read0->getIdLength() != read1->getIdLength()) {
        return false;
    }
    for (unsigned i = 0; i < read0->getIdLength(); i++) {
        char c0 = read0->getId()[i];
        char c1 = read1->getId()[i];

        if (c0 != c1) return false;
 
        // don't parse the read ID after the first space or slash, which can represent metadata (or which half of the mate pair the read is).
        if (c0 == ' ' || c0 == '/') return true;  
    }
    return true;
}

PairedAlignerOptions::PairedAlignerOptions(const char* i_commandLine)
    : AlignerOptions(i_commandLine, true),
    minSpacing(DEFAULT_MIN_SPACING),
    maxSpacing(DEFAULT_MAX_SPACING)
{
}

void PairedAlignerOptions::usageMessage()
{
    AlignerOptions::usageMessage();
    printf(
        "  -s   min and max spacing to allow between paired ends (default: %d %d)\n",
        DEFAULT_MIN_SPACING,
        DEFAULT_MAX_SPACING);
}

bool PairedAlignerOptions::parse(const char** argv, int argc, int& n)
{
    if (strcmp(argv[n], "-s") == 0) {
        if (n + 2 < argc) {
            minSpacing = atoi(argv[n+1]);
            maxSpacing = atoi(argv[n+2]);
            n += 2;
            return true;
        } 
        return false;
    }
    return AlignerOptions::parse(argv, argc, n);
}

PairedAlignerContext::PairedAlignerContext(AlignerExtension* i_extension)
    : AlignerContext( 0,  NULL, NULL, i_extension)
{
}

AlignerOptions* PairedAlignerContext::parseOptions(int i_argc, const char **i_argv, const char *i_version)
{
    argc = i_argc;
    argv = i_argv;
    version = i_version;

    PairedAlignerOptions* options = new PairedAlignerOptions(
        "snap paired <index-dir> <input file(s)> <read2.fq> [<options>]\n"
        "   where <input file(s)> is a list of files to process.  FASTQ\n"
        "   files must come in pairs, since each read end is in a separate file.");
    options->extra = extension->extraOptions();
    if (argc < 2) {
        options->usage();
    }

    options->indexDir = argv[0];
    //
    // Figure out how many inputs there are.  All options begin with a '-', so count the
    // args until we hit an option.  FASTQ files come in pairs, and each pair only counts
    // as one input.
    //
    int nInputs = 0;
    bool foundFirstHalfOfFASTQ = false;
    for (int i = 1; i < argc; i++) {
        if (argv[i][0] == '-') {
                break;
        }
        
        if (stringEndsWith(argv[i],".sam")) {
            if (foundFirstHalfOfFASTQ) {
                fprintf(stderr,"For the paired aligner, FASTQ files must come in pairs.  I found SAM file '%s' after first half FASTQ file '%s'.\n",
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
        options->usage();
    }
    //
    // Now build the input array.
    //
    options->nInputs = nInputs;
    options->inputs = new SNAPInput[nInputs];
    int i;
    int whichInput = 0;
    for (i = 1; i < argc; i++) {
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
    _ASSERT(whichInput == nInputs);

    for (/* i initialized by previous loop*/; i < argc; i++) {
        if (!options->parse(argv, argc, i)) {
            options->usage();
        }
    }

    return options;
}

void PairedAlignerContext::initialize()
{
    AlignerContext::initialize();
    PairedAlignerOptions* options2 = (PairedAlignerOptions*) options;
    minSpacing = options2->minSpacing;
    maxSpacing = options2->maxSpacing;
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

    int maxReadSize = 10000;
#if     1
    SmarterPairedEndAligner *aligner = new SmarterPairedEndAligner(
            index,
            maxReadSize,
            confDiff,
            maxHits,
            maxDist,
            numSeeds,
            minSpacing,
            maxSpacing,
            adaptiveConfDiff);
#endif  // 0
    BigAllocator *allocator = new BigAllocator(100 * 1024 * 1024);
//    IntersectingPairedEndAligner *aligner = new IntersectingPairedEndAligner(index, maxReadSize, maxHits, maxDist, numSeeds, minSpacing, maxSpacing, allocator);
            
    /*BigAllocator *allocator = new BigAllocator(ThirdPairedEndAligner::getBigAllocatorReservation(index, maxHits, maxReadSize, index->getSeedLength(), numSeeds));
    ThirdPairedEndAligner *aligner = new(allocator) ThirdPairedEndAligner(index, maxReadSize, maxHits, maxDist, numSeeds, minSpacing, maxSpacing, allocator);
    allocator->assertAllMemoryUsed();*/

    ReadWriter *readWriter = this->readWriter;

    // Align the reads.
    Read *read0;
    Read *read1;
    _int64 readNum = 0;
    while (supplier->getNextReadPair(&read0,&read1)) {
        if (1 != selectivity && GoodFastRandom(selectivity-1) != 0) {
            //
            // Skip this read.
            //
            continue;
        }

#ifdef PROFILE
        readNum++;
        bool record = (readNum % 1000 == 0);
        _int64 start;
        if (record) {
            start = timeInNanos();
        }
#endif
            
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

        // Skip the pair if there are too many Ns or 2s.
        int maxDist = this->maxDist;
        bool useful0 = read0->getDataLength() >= 50 && (int)read0->countOfNs() <= maxDist;
        bool useful1 = read1->getDataLength() >= 50 && (int)read1->countOfNs() <= maxDist;
        if (!useful0 && !useful1) {
            PairedAlignmentResult result;
            result.status[0] = NotFound;
            result.status[1] = NotFound;
            result.location[0] = 0xFFFFFFFF;
            result.location[1] = 0xFFFFFFFF;
            writePair(read0, read1, &result);
            continue;
        } else {
            // Here one the reads might still be hopeless, but maybe we can align the other.
            stats->usefulReads += (useful0 && useful1) ? 2 : 1;
        }

        PairedAlignmentResult result;
        aligner->align(read0, read1, &result);

        writePair(read0, read1, &result);

        updateStats((PairedAlignerStats*) stats, read0, read1, &result);

#ifdef PROFILE
        if (record) {
            _int64 end = timeInNanos();
            printf("%d %lld %.*s %s %s %lld\n",
                context->threadNum,
                readNum,
                read0.getIdLength(), read0.getId(),
                AlignmentResultToString(result.status[0]),
                AlignmentResultToString(result.status[1]),
                end - start);
        }
#endif

    }
    //printf("Time in s: %lld: thread ran out of work.  Last range was %8lld bytes in %4lldms, starting at %10lld.  Total %4d ranges and %10lld bytes.\n",timeInMillis() / 1000, rangeLength, timeInMillis() - rangeStartTime, rangeStart, totalRanges, totalBytes);

//    aligner->~ThirdPairedEndAligner();
//    delete allocator;

    delete aligner;
    delete supplier;
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
}

    void 
PairedAlignerContext::typeSpecificBeginIteration()
{
    if (1 == options->nInputs) {
        //
        // We've only got one input, so just connect it directly to the consumer.
        //
        pairedReadSupplierGenerator = options->inputs[0].createPairedReadSupplierGenerator(options->numThreads, index->getGenome(), options->clipping);
    } else {
        //
        // We've got multiple inputs, so use a MultiInputReadSupplier to combine the individual inputs.
        //
        PairedReadSupplierGenerator **generators = new PairedReadSupplierGenerator *[options->nInputs];
        for (int i = 0; i < options->nInputs; i++) {
            generators[i] = options->inputs[i].createPairedReadSupplierGenerator(options->numThreads, index->getGenome(), options->clipping);
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
