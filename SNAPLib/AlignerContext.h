/*++

Module Name:

    AlignerContext.h

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

#pragma once
#include "stdafx.h"
#include "Genome.h"
#include "Range.h"
#include "RangeSplitter.h"
#include "AlignerOptions.h"
#include "AlignerStats.h"
#include "ParallelTask.h"

class AlignerExtension;

/*++
    Common context state shared across threads during alignment process
--*/
class AlignerContext : public TaskContextBase
{
public:

    AlignerContext(int i_argc, const char **i_argv, const char *i_version, AlignerExtension* i_extension = NULL);

    ~AlignerContext();

    // running alignment

    void runAlignment(int argc, const char **argv, const char *version);
    
    // ParallelTask template

    void initializeThread();

    void runThread();

    void finishThread(AlignerContext* common);

protected:
    
    void printStatsHeader();
    
    void printStats();
    
    void beginIteration();
    
    void finishIteration();
    
    // advance to next iteration in range, return false when past end
    bool nextIteration();
    
    // overrideable by concrete single/paired alignment subclasses
    
    // parse options from the command line
    virtual AlignerOptions* parseOptions(int argc, const char **argv, const char *version) = 0;
    
    // initialize from options
    virtual void initialize();

    // new stats object
    virtual AlignerStats* newStats() = 0;
    
    // instantiate and run a parallel task
    virtual void runTask() = 0;

    // run single thread within single iteration
    virtual void runIterationThread() = 0;

    // common state across all threads
    GenomeIndex        *index;
    ParallelSAMWriter  *parallelSamWriter;
    _int64              alignStart;
    _int64              alignTime;
    RangeSplitter       fileSplitterState;
    AlignerOptions     *options;
    AlignerStats       *stats;
    AlignerExtension   *extension;
    friend class AlignerContext2;
    unsigned            maxDist;
    int                 numSeeds;
    int                 maxHits;
    int                 confDiff;
    int                 adaptiveConfDiff;
    bool                computeError;
    const char         *inputFilename;
    bool                inputFileIsFASTQ;   // Else SAM
    RangeSplitter      *fileSplitter;
    unsigned            selectivity;
    bool                detailedStats;
    ReadClippingType    clipping;
    int                 argc;
    const char        **argv;
    const char         *version;

    // iteration variables
    int                 confDiff_;
    int                 maxHits_;
    int                 maxDist_;
    int                 numSeeds_;
    int                 adaptiveConfDiff_;

    // Per-thread context state used during alignment process
    SAMWriter          *samWriter;
};

// abstract class for extending base context

class AlignerExtension
{
public:

    virtual ~AlignerExtension() {}

    virtual AbstractOptions* extraOptions() { return NULL; }

    virtual AbstractStats* extraStats() { return NULL; }

    virtual bool skipAlignment() { return false; }

    virtual void initialize() {}

    virtual void beginIteration() {}

    virtual AlignerExtension* copy() { return new AlignerExtension(); }

    virtual void beginThread() {}

    virtual void finishThread() {}

    virtual void finishIteration() {}

    virtual void printStats() {}

    virtual void finishAlignment() {}
};
