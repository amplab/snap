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
#include "RangeSplitter.h"
#include "AlignerOptions.h"
#include "AlignerStats.h"
#include "ParallelTask.h"
#include "GenomeIndex.h"

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

    void runAlignment(int argc, const char **argv, const char *version, unsigned *nArgsConsumed);
    
    // ParallelTask template

    void initializeThread();

    void runThread();

    void finishThread(AlignerContext* common);

    void printStats();
    
    void beginIteration();
    
    void finishIteration();
    
    // advance to next iteration in range, return false when past end
    bool nextIteration();
    
    // overrideable by concrete single/paired alignment subclasses
    
    // parse options from the command line
    AlignerOptions* parseOptions(int argc, const char **argv, const char *version, unsigned *argsConsumed, bool paired);
    
    // initialize from options
    virtual bool initialize();

    // new stats object
    virtual AlignerStats* newStats() = 0;
    
    // instantiate and run a parallel task
    virtual void runTask() = 0;

    // run single thread within single iteration
    virtual void runIterationThread() = 0;

    virtual void typeSpecificBeginIteration() = 0;
    virtual void typeSpecificNextIteration() = 0;

    virtual bool isPaired() = 0;

    friend class AlignerContext2;
 
    // common state across all threads
    GenomeIndex                         *index;
    ReadWriterSupplier                  *writerSupplier;
    ReaderContext                        readerContext;
    _int64                               alignStart;
    _int64                               alignTime;
    AlignerOptions                      *options;
    AlignerStats                        *stats;
    AlignerExtension                    *extension;
    unsigned                             maxDist;
    unsigned                             numSeedsFromCommandLine;
    double                               seedCoverage;
    unsigned                             minWeightToCheck;
    int                                  maxHits;
    bool                                 detailedStats;
    ReadClippingType                     clipping;
    unsigned                             extraSearchDepth;
    int                                  argc;
    const char                         **argv;
    const char                          *version;
    FILE                                *perfFile;
    bool                                 noUkkonen;
    bool                                 noOrderedEvaluation;
	bool								 noTruncation;
    int                                  maxSecondaryAlignmentAdditionalEditDistance;
	int									 maxSecondaryAlignments;
    int                                  maxSecondaryAlignmentsPerContig;
	unsigned							 minReadLength;


    // iteration variables
    int                 maxHits_;
    int                 maxDist_;

    // Per-thread context state used during alignment process
    ReadWriter         *readWriter;
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

	virtual bool runIterationThread(PairedReadSupplier* supplier, AlignerContext* threadContext) { return false; }
	
	virtual bool runIterationThread(ReadSupplier* supplier, AlignerContext* threadContext) { return false; }

    virtual void finishThread() {}

    virtual void finishIteration() {}

    virtual void printStats() {}

    virtual void finishAlignment() {}
};
