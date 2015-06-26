/*++

Module Name:

    SingleAligner.cpp

Abstract:

    Functions for running the single end aligner sub-program.

Authors:

    Matei Zaharia, February, 2012

Environment:
`
    User mode service.

Revision History:

    Adapted from cSNAP, which was in turn adapted from the scala prototype

--*/

#pragma once
#include "stdafx.h"
#include "AlignerContext.h"
#include "ReadSupplierQueue.h"

struct PairedAlignerStats;

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
};
