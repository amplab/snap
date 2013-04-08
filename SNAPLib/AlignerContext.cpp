/*++


Module Name:

    AlignerContext.cpp

Abstract:

    Common parameters for running single & paired alignment.

Authors:

    Ravi Pandya, May, 2012

Environment:
`
    User mode service.

Revision History:

    Integrated from SingleAligner.cpp & PairedAligner.cpp

--*/

#include "stdafx.h"
#include "Compat.h"
#include "options.h"
#include "AlignerOptions.h"
#include "AlignerContext.h"
#include "AlignerStats.h"
#include "BaseAligner.h"
#include "FileFormat.h"
#include "exit.h"

using std::max;
using std::min;

//
// Save the index & index directory globally so that we don't need to reload them on multiple runs.
//
GenomeIndex *g_index = NULL;
char *g_indexDirectory = NULL;

AlignerContext::AlignerContext(int i_argc, const char **i_argv, const char *i_version, AlignerExtension* i_extension)
    :
    index(NULL),
    writerSupplier(NULL),
    options(NULL),
    stats(NULL),
    extension(i_extension != NULL ? i_extension : new AlignerExtension()),
    readWriter(NULL),
    argc(i_argc),
    argv(i_argv),
    version(i_version),
    perfFile(NULL)
{
}

AlignerContext::~AlignerContext()
{
    delete extension;
    if (NULL != perfFile) {
        fclose(perfFile);
    }
}

void AlignerContext::runAlignment(int argc, const char **argv, const char *version, unsigned *argsConsumed)
{
    options = parseOptions(argc, argv, version, argsConsumed);
    initialize();
    extension->initialize();
    
    if (! extension->skipAlignment()) {
        printStatsHeader();
        do {

            beginIteration();

            runTask();
            
            finishIteration();
            
            printStats();

        } while (nextIteration());
    }

    extension->finishAlignment();
}

    void
AlignerContext::initializeThread()
{
    stats = newStats(); // separate copy per thread
    stats->extra = extension->extraStats();
    readWriter = writerSupplier != NULL ? writerSupplier->getWriter() : NULL;
    extension = extension->copy();
}

    void
AlignerContext::runThread()
{
    extension->beginThread();
    runIterationThread();
    if (readWriter != NULL) {
        readWriter->close();
        delete readWriter;
    }
    extension->finishThread();
}
    
    void
AlignerContext::finishThread(AlignerContext* common)
{
    common->stats->add(stats);
    delete stats;
    stats = NULL;
    delete extension;
    extension = NULL;
}

    void
AlignerContext::initialize()
{
    if (g_indexDirectory == NULL || strcmp(g_indexDirectory, options->indexDir) != 0) {
        delete g_index;
        g_index = NULL;
        delete g_indexDirectory;
        g_indexDirectory = new char [strlen(options->indexDir) + 1];
        strcpy(g_indexDirectory, options->indexDir);

        printf("Loading index from directory... ");
        fflush(stdout);
        _int64 loadStart = timeInMillis();
        index = GenomeIndex::loadFromDirectory((char*) options->indexDir);
        if (index == NULL) {
            fprintf(stderr, "Index load failed, aborting.\n");
            soft_exit(1);
        }
        g_index = index;

        _int64 loadTime = timeInMillis() - loadStart;
        printf("%llds.  %u bases, seed size %d\n",
            loadTime / 1000, index->getGenome()->getCountOfBases(), index->getSeedLength());
    } else {
        index = g_index;
    }

    if (options->outputFileTemplate != NULL && (options->maxHits.size() > 1 || options->maxDist.size() > 1 || options->numSeeds.size() > 1
                || options->confDiff.size() > 1 || options->adaptiveConfDiff.size() > 1)) {
        fprintf(stderr, "WARNING: You gave ranges for some parameters, so SAM files will be overwritten!\n");
    }

    confDiff_ = options->confDiff.start;
    maxHits_ = options->maxHits.start;
    maxDist_ = options->maxDist.start;
    numSeeds_ = options->numSeeds.start;
    adaptiveConfDiff_ = options->adaptiveConfDiff.start;

    if (options->perfFileName != NULL) {
        perfFile = fopen(options->perfFileName,"a");
        if (NULL == perfFile) {
            fprintf(stderr,"Unable to open perf file '%s'\n", options->perfFileName);
            soft_exit(1);
        }
    }

    DataSupplier::ThreadCount = options->numThreads;
}

    void
AlignerContext::printStatsHeader()
{
    printf("ConfDif\tMaxHits\tMaxDist\tMaxSeed\tConfAd\t%%Used\t%%Unique\t%%Multi\t%%!Found\t%%Error\t%%Pairs\tReads/s\n");
}

    void
AlignerContext::beginIteration()
{
    writerSupplier = NULL;
    if (NULL != options->outputFileTemplate) {
        const FileFormat* format = 
            FileFormat::SAM[0]->isFormatOf(options->outputFileTemplate) ? FileFormat::SAM[options->useM] :
            FileFormat::BAM[0]->isFormatOf(options->outputFileTemplate) ? FileFormat::BAM[options->useM] :
            NULL;
        if (format != NULL) {
            writerSupplier = format->getWriterSupplier(options, index->getGenome());
            ReadWriter* headerWriter = writerSupplier->getWriter();
            headerWriter->writeHeader(options->sortOutput, argc, argv, version, options->rgLineContents);
            headerWriter->close();
            delete headerWriter;
        } else {
            fprintf(stderr, "warning: no output, unable to determine format of output file %s\n", options->outputFileTemplate);
        }
    }

    alignStart = timeInMillis();

    clipping = options->clipping;
    totalThreads = options->numThreads;
    computeError = options->computeError;
    bindToProcessors = options->bindToProcessors;
    maxDist = maxDist_;
    maxHits = maxHits_;
    numSeeds = numSeeds_;
    confDiff = confDiff_;
    adaptiveConfDiff = adaptiveConfDiff_;

    if (stats != NULL) {
        delete stats;
    }
    stats = newStats();
    stats->extra = extension->extraStats();
    extension->beginIteration();
    typeSpecificBeginIteration();
}

    void
AlignerContext::finishIteration()
{
    extension->finishIteration();

    if (NULL != writerSupplier) {
        writerSupplier->close();
        delete writerSupplier;
        writerSupplier = NULL;
    }

    alignTime = /*timeInMillis() - alignStart -- use the time from ParallelTask.h, that may exclude memory allocation time*/ time;
}

    bool
AlignerContext::nextIteration()
{
     typeSpecificNextIteration();
     if ((adaptiveConfDiff_ += options->adaptiveConfDiff.step) > options->adaptiveConfDiff.end) {
        adaptiveConfDiff_ = options->adaptiveConfDiff.start;
        if ((numSeeds_ += options->numSeeds.step) > options->numSeeds.end) {
            numSeeds_ = options->numSeeds.start;
            if ((maxDist_ += options->maxDist.step) > options->maxDist.end) {
                maxDist_ = options->maxDist.start;
                if ((maxHits_ += options->maxHits.step) > options->maxHits.end) {
                    maxHits_ = options->maxHits.start;
                    if ((confDiff_ += options->confDiff.step) > options->confDiff.end) {
                        return false;
                    }
                }
            }
        }
    }

    return true;
}

    void
AlignerContext::printStats()
{
    double usefulReads = max((double) stats->usefulReads, 1.0);
    char errorRate[16];
    if (options->computeError) {
        _int64 numSingle = max(stats->singleHits, (_int64) 1);
        snprintf(errorRate, sizeof(errorRate), "%0.3f%%", (100.0 * stats->errors) / numSingle);
    } else {
        snprintf(errorRate, sizeof(errorRate), "-");
    }
    printf("%d\t%d\t%d\t%d\t%d\t%0.2f%%\t%0.2f%%\t%0.2f%%\t%0.2f%%\t%s\t%0.2f%%\t%.0f (at: %lld)\n",
            confDiff_, maxHits_, maxDist_, numSeeds_, adaptiveConfDiff_,
            100.0 * usefulReads / max(stats->totalReads, (_int64) 1),
            100.0 * stats->singleHits / usefulReads,
            100.0 * stats->multiHits / usefulReads,
            100.0 * stats->notFound / usefulReads,
            errorRate,
            100.0 * stats->alignedAsPairs / usefulReads,
            (1000.0 * usefulReads) / max(alignTime, (_int64) 1), alignTime);
    if (NULL != perfFile) {
        fprintf(perfFile, "%d\t%d\t%d\t%d\t%d\t%0.2f%%\t%0.2f%%\t%0.2f%%\t%0.2f%%\t%s\t%0.2f%%\t%.0f\n",
                confDiff_, maxHits_, maxDist_, numSeeds_, adaptiveConfDiff_,
                100.0 * usefulReads / max(stats->totalReads, (_int64) 1),
                100.0 * stats->singleHits / usefulReads,
                100.0 * stats->multiHits / usefulReads,
                100.0 * stats->notFound / usefulReads,
                errorRate,
                100.0 * stats->alignedAsPairs / usefulReads,
                (1000.0 * usefulReads) / max(alignTime, (_int64) 1));

        for (AlignerStats::ThreadPerfEntry *threadEntry = stats->threadEntry; NULL != threadEntry; threadEntry = threadEntry->next) {
            fprintf(perfFile, "%d\t%d\t%lld\t%lld\t0x%llx\t%llx\t%llx\t%llx\t%llx\n", threadEntry->threadNumber, threadEntry->threadId, threadEntry->nReads, threadEntry->lvCalls, threadEntry->candidateEntries,
                threadEntry->hashAnchor[0], threadEntry->hashAnchor[1], threadEntry->alignerObject, threadEntry->stackPointer);
        }
        fprintf(perfFile,"\n");
    }
    // Running counts to compute a ROC curve (with error rate and %aligned above a given MAPQ)
    double totalAligned = 0;
    double totalErrors = 0;
    for (int i = AlignerStats::maxMapq; i >= 0; i--) {
        totalAligned += stats->mapqHistogram[i];
        totalErrors += stats->mapqErrors[i];
        double truePositives = (totalAligned - totalErrors) / max(stats->totalReads, (_int64) 1);
        double falsePositives = totalErrors / totalAligned;
        if (i <= 10 || i % 2 == 0 || i == 69) {
            printf("%d\t%d\t%d\t%.3f\t%.2E\n", i, stats->mapqHistogram[i], stats->mapqErrors[i], truePositives, falsePositives);
        }
    }

    stats->printHistograms(stdout);

#ifdef  TIME_STRING_DISTANCE
    printf("%llds, %lld calls in BSD noneClose, not -1\n",  stats->nanosTimeInBSD[0][1]/1000000000, stats->BSDCounts[0][1]);
    printf("%llds, %lld calls in BSD noneClose, -1\n",      stats->nanosTimeInBSD[0][0]/1000000000, stats->BSDCounts[0][0]);
    printf("%llds, %lld calls in BSD close, not -1\n",      stats->nanosTimeInBSD[1][1]/1000000000, stats->BSDCounts[1][1]);
    printf("%llds, %lld calls in BSD close, -1\n",          stats->nanosTimeInBSD[1][0]/1000000000, stats->BSDCounts[1][0]);
    printf("%llds, %lld calls in Hamming\n",                stats->hammingNanos/1000000000,         stats->hammingCount);
#endif  // TIME_STRING_DISTANCE

    extension->printStats();
}

