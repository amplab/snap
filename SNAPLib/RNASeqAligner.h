/*++

Module Name:

    SingleAligner.cpp

Abstract:

    Functions for running the single end aligner sub-program.

Authors:

    Andrew Magis, March 2013
    
    Adapted from PairedAligner.h by Matei Zaharia

Environment:
`
    User mode service.

Revision History:

    Adapted from cSNAP, which was in turn adapted from the scala prototype

--*/

#pragma once
#include "stdafx.h"
#include "AlignerContext.h"

struct RNASeqAlignerStats;

class RNASeqAlignerContext : public AlignerContext
{
public:

    RNASeqAlignerContext(AlignerExtension* i_extension = NULL);
    
protected:

    // AlignerContext
    
    virtual AlignerOptions* parseOptions(int argc, const char **argv, const char *version);

    virtual void initialize();

    virtual AlignerStats* newStats();
    
    virtual void runTask();
    
    virtual void runIterationThread();

    // for subclasses

    virtual void writePair(Read* read0, Read* read1, PairedAlignmentResult* result);

    virtual void updateStats(RNASeqAlignerStats* stats, Read* read0, Read* read1, PairedAlignmentResult* result);

protected:

    int                 minSpacing;
    int                 maxSpacing;
    const char         *fastqFile1;
    bool                ignoreMismatchedIDs;
};

struct RNASeqAlignerOptions : public AlignerOptions
{
    RNASeqAlignerOptions(const char* i_commandLine);

    virtual void usageMessage();

    virtual bool parse(const char** argv, int argc, int& n);

    int minSpacing;
    int maxSpacing;
    const char* fastqFile1;
};
