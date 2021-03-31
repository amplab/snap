/*++

Module Name:

    PairedAligner.h

Abstract:

    Header for the paired end aligner scaffolding

Authors:

    Matei Zaharia, February, 2012

Environment:

    User mode service.

Revision History:

    Adapted from cSNAP, which was in turn adapted from the scala prototype

--*/

#pragma once
#include "stdafx.h"
#include "AlignerContext.h"
#include "ReadSupplierQueue.h"

struct PairedAlignerStats;

static const int DEFAULT_BATCH_SIZE_IS_ESTIMATION = 256 * 1024;
// Parameters from BWA-MEM
static const int OUTLIER_BOUND = 2;
static const int MAPPING_BOUND = 3;
static const int MAX_STDDEV = 4;

class PairedAlignerContext : public AlignerContext
{
public:

    PairedAlignerContext(AlignerExtension* i_extension = NULL);
    
protected:

    // AlignerContext
    
    virtual bool initialize();

    virtual AlignerStats* newStats();
    
    virtual void runTask();
    
    virtual void runIterationThread();

    // for subclasses

    virtual void updateStats(PairedAlignerStats* stats, Read* read0, Read* read1, PairedAlignmentResult* result, bool useful0, bool useful1);

    bool isPaired() {return true;}

    void computeSpacingDist(GenomeDistance* pairedEndSpacing, int* minSpacing, int* maxSpacing, double* avg, double* stddev);

protected:

    virtual void typeSpecificBeginIteration();
    virtual void typeSpecificNextIteration();

    PairedReadSupplierGenerator *pairedReadSupplierGenerator;
 
    int                 minSpacing;
    int                 maxSpacing;
    bool                forceSpacing;
    unsigned            intersectingAlignerMaxHits;
    unsigned            maxCandidatePoolSize;
    const char         *fastqFile1;
    bool                ignoreMismatchedIDs;
    bool                quicklyDropUnpairedReads;
    bool                inferSpacing;
    int                 maxSeedsSingleEnd;
    int                 minScoreRealignment;
    int                 minScoreGapRealignmentALT;
    int                 minAGScoreImprovement;
    bool                enableHammingScoringBaseAligner;

    GenomeDistance pairedEndSpacing[DEFAULT_BATCH_SIZE_IS_ESTIMATION];      // Spacing between reads aligned as FR pairs
    static int compareBySpacing(const void *first_, const void *second_);

	friend class AlignerContext2;
};

struct PairedAlignerOptions : public AlignerOptions
{
    PairedAlignerOptions(const char* i_commandLine);

    virtual void usageMessage();

    virtual bool parse(const char** argv, int argc, int& n, bool *done);

    virtual bool isPaired() { return true; }

    int         minSpacing;
    int         maxSpacing;
    bool        forceSpacing;
    unsigned    intersectingAlignerMaxHits;
    unsigned    maxCandidatePoolSize;
    bool        quicklyDropUnpairedReads;
    bool        inferSpacing;
    int         maxSeedsSingleEnd;
    int         minScoreRealignment;
    int         minScoreGapRealignmentALT;
    int         minAGScoreImprovement;
    bool        enableHammingScoringBaseAligner;
};
