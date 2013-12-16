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
GenomeIndex *t_index = NULL;
char *t_indexDirectory = NULL;
GenomeIndex *c_index = NULL;
char *c_indexDirectory = NULL;

AlignerContext::AlignerContext(int i_argc, const char **i_argv, const char *i_version, AlignerExtension* i_extension)
    :
    index(NULL),
    transcriptome(NULL),
    contamination(NULL),
    c_filter(NULL),
    gtf(NULL),
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

    if (NULL != index) {
        delete index;
    }
    
    if (NULL != transcriptome) {
        delete transcriptome;
    }

    if (NULL != contamination) {
        delete contamination;
    }
    
    if (NULL != gtf) {
        delete gtf;
    }

    if (NULL != c_filter) {
        delete c_filter;
    }

    delete extension;
    if (NULL != perfFile) {
        fclose(perfFile);
    }
}

void AlignerContext::runAlignment(int argc, const char **argv, const char *version, unsigned *argsConsumed)
{
    options = parseOptions(argc, argv, version, argsConsumed);
#ifdef _MSC_VER
    useTimingBarrier = options->useTimingBarrier;
#endif
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
    PrintBigAllocProfile();
    PrintWaitProfile();
    
    //Process all the intervals and read counts
    gtf->AnalyzeReadIntervals();
    gtf->WriteReadCounts();

    if (c_filter != NULL) {
      //c_filter->Print();
      c_filter->Write();
    }
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

        if (strcmp(options->indexDir, "-") != 0) {
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
            printf("no alignment, input/output only\n");
        }
    } else {
        index = g_index;
    }
    
    if (t_indexDirectory == NULL || strcmp(t_indexDirectory, options->transcriptomeDir) != 0) {
    
        delete t_index;
        t_index = NULL;
        delete t_indexDirectory;
        t_indexDirectory = new char [strlen(options->transcriptomeDir) + 1];
        strcpy(t_indexDirectory, options->transcriptomeDir);

        if (strcmp(options->transcriptomeDir, "-") != 0) {
            printf("Loading transcriptome index from directory... ");
            fflush(stdout);
            _int64 loadStart = timeInMillis();
            transcriptome = GenomeIndex::loadFromDirectory((char*) options->transcriptomeDir);
            if (transcriptome == NULL) {
                fprintf(stderr, "Transcriptome index load failed, aborting.\n");
                soft_exit(1);
            }
            t_index = transcriptome;
            
            _int64 loadTime = timeInMillis() - loadStart;
            printf("%llds.  %u bases, seed size %d\n",
                loadTime / 1000, transcriptome->getGenome()->getCountOfBases(), transcriptome->getSeedLength());
        } else {
            printf("no alignment, input/output only\n");
        }
    } else {
        transcriptome = t_index;
    }
 
    //Load the contamination database
    if (options->contaminationDir != NULL) {
      if (c_indexDirectory == NULL || strcmp(c_indexDirectory, options->contaminationDir) != 0) {

        delete c_index;
        c_index = NULL;
        delete c_indexDirectory;
        c_indexDirectory = new char [strlen(options->contaminationDir) + 1];
        strcpy(c_indexDirectory, options->contaminationDir);

        if (strcmp(options->contaminationDir, "-") != 0) {
            printf("Loading contamination index from directory... ");
            fflush(stdout);
            _int64 loadStart = timeInMillis();
            contamination = GenomeIndex::loadFromDirectory((char*) options->contaminationDir);
            if (contamination == NULL) {
                fprintf(stderr, "Contamination index load failed, aborting.\n");
                soft_exit(1);
            }
            c_index = contamination;

            _int64 loadTime = timeInMillis() - loadStart;
            printf("%llds.  %u bases, seed size %d\n",
                loadTime / 1000, contamination->getGenome()->getCountOfBases(), contamination->getSeedLength());

            if (c_filter != NULL) {
              delete c_filter;
            }
            //Create the contamination filter
            c_filter = new ContaminationFilter(contamination->getGenome(), options->outputFileTemplate);

          } else {
            printf("no alignment, input/output only\n");
          }
      } else {
          contamination = c_index;
      }
    }

    //Create GTFReader and load the annotation
    gtf = new GTFReader(options->outputFileTemplate);
    gtf->Load(options->annotation);
    
    if (options->outputFileTemplate != NULL && (options->maxHits.size() > 1 || options->maxDist.size() > 1)) {
        fprintf(stderr, "WARNING: You gave ranges for some parameters, so SAM files will be overwritten!\n");
    }

    maxHits_ = options->maxHits.start;
    maxDist_ = options->maxDist.start;
    extraSearchDepth = options->extraSearchDepth;

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
    printf("MaxHits\tMaxDist\t%%Used\t%%Unique\t%%Multi\t%%!Found\t%%Error\t%%Pairs\tlvCalls\tNumReads\tReads/s\n");
}

    void
AlignerContext::beginIteration()
{
    writerSupplier = NULL;
    alignStart = timeInMillis();
    clipping = options->clipping;
    totalThreads = options->numThreads;
    computeError = options->computeError;
    bindToProcessors = options->bindToProcessors;
    maxDist = maxDist_;
    maxHits = maxHits_;
    numSeedsFromCommandLine = options->numSeedsFromCommandLine;
    seedCoverage = options->seedCoverage;
    if (stats != NULL) {
        delete stats;
    }
    stats = newStats();
    stats->extra = extension->extraStats();
    extension->beginIteration();
    
    readerContext.clipping = options->clipping;
    readerContext.defaultReadGroup = options->defaultReadGroup;
    readerContext.genome = index != NULL ? index->getGenome() : NULL;
    readerContext.transcriptome = transcriptome != NULL ? transcriptome->getGenome() : NULL;
    readerContext.gtf = gtf != NULL ? gtf : NULL;
    readerContext.paired = false;
    readerContext.ignoreSecondaryAlignments = options->ignoreSecondaryAlignments;
    DataSupplier::ExpansionFactor = options->expansionFactor;
	readerContext.header = NULL;
	readerContext.headerLength = 0;
	readerContext.headerBytes = 0;

    typeSpecificBeginIteration();

    if (NULL != options->outputFileTemplate) {
        const FileFormat* format = 
            FileFormat::SAM[0]->isFormatOf(options->outputFileTemplate) ? FileFormat::SAM[options->useM] :
            FileFormat::BAM[0]->isFormatOf(options->outputFileTemplate) ? FileFormat::BAM[options->useM] :
            NULL;
        if (format != NULL) {
            writerSupplier = format->getWriterSupplier(options, readerContext.genome, readerContext.transcriptome, readerContext.gtf);
            ReadWriter* headerWriter = writerSupplier->getWriter();
            headerWriter->writeHeader(readerContext, options->sortOutput, argc, argv, version, options->rgLineContents);
            headerWriter->close();
            delete headerWriter;
        } else {
            fprintf(stderr, "warning: no output, unable to determine format of output file %s\n", options->outputFileTemplate);
        }
    }
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
    if ((maxDist_ += options->maxDist.step) > options->maxDist.end) {
        maxDist_ = options->maxDist.start;
        if ((maxHits_ += options->maxHits.step) > options->maxHits.end) {
            return false;
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
    printf("%d\t%d\t%0.2f%%\t%0.2f%%\t%0.2f%%\t%0.2f%%\t%s\t%0.2f%%\t%lld\t%lld\t%.0f (at: %lld)\n",
            maxHits_, maxDist_, 
            100.0 * usefulReads / max(stats->totalReads, (_int64) 1),
            100.0 * stats->singleHits / usefulReads,
            100.0 * stats->multiHits / usefulReads,
            100.0 * stats->notFound / usefulReads,
            errorRate,
            100.0 * stats->alignedAsPairs / usefulReads,
            stats->lvCalls,
            stats->totalReads,
            (1000.0 * usefulReads) / max(alignTime, (_int64) 1), 
            alignTime);
    if (NULL != perfFile) {
        fprintf(perfFile, "%d\t%d\t%0.2f%%\t%0.2f%%\t%0.2f%%\t%0.2f%%\t%s\t%0.2f%%\t%lld\t%lld\tt%.0f\n",
                maxHits_, maxDist_, 
                100.0 * usefulReads / max(stats->totalReads, (_int64) 1),
                100.0 * stats->singleHits / usefulReads,
                100.0 * stats->multiHits / usefulReads,
                100.0 * stats->notFound / usefulReads,
                stats->lvCalls,
                errorRate,
                100.0 * stats->alignedAsPairs / usefulReads,
                stats->totalReads,
                (1000.0 * usefulReads) / max(alignTime, (_int64) 1));

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
//            printf("%d\t%d\t%d\t%.3f\t%.2E\n", i, stats->mapqHistogram[i], stats->mapqErrors[i], truePositives, falsePositives);
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
