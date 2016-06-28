/*++

Module Name:

TenX.h

Abstract:

Headers for the 10x aligner

Authors:

Hongyi Xin and Bill Bolosky, May, 2016

Environment:

User mode service.

Revision History:

Cloned and modified the paired end aligner

--*/

#pragma once
#include "stdafx.h"
#include "AlignerContext.h"
#include "ReadSupplierQueue.h"

struct TenXAlignerStats;

class TenXAlignerContext : public AlignerContext
{
public:

    TenXAlignerContext(AlignerExtension* i_extension = NULL);

protected:

    // AlignerContext

    virtual bool initialize();

    virtual AlignerStats* newStats();

    virtual void runTask();

    virtual void runIterationThread();

    // for subclasses

    virtual void updateStats(TenXAlignerStats* stats, Read* read0, Read* read1, PairedAlignmentResult* result, bool useful0, bool useful1);

    bool isPaired() { return true; }

protected:

    virtual void typeSpecificBeginIteration();
    virtual void typeSpecificNextIteration();

    PairedReadSupplierGenerator *pairedReadSupplierGenerator;

    int                 minSpacing;
    int                 maxSpacing;
    int                 maxClusterSize;
    int                 minReadsPerCluster;
    bool                forceSpacing;
    unsigned            intersectingAlignerMaxHits;
    unsigned            maxCandidatePoolSize;
    const char         *fastqFile1;
    bool                ignoreMismatchedIDs;
    bool                quicklyDropUnpairedReads;

    friend class AlignerContext2;
};

struct TenXAlignerOptions : public AlignerOptions
{
    TenXAlignerOptions(const char* i_commandLine);

    virtual void usageMessage();

    virtual bool parse(const char** argv, int argc, int& n, bool *done);

    virtual bool isPaired() { return true; }

    int         minSpacing;
    int         maxSpacing;
    int         maxClusterSize;
    int         minReadsPerCluster;
    bool        forceSpacing;
    unsigned    intersectingAlignerMaxHits;
    unsigned    maxCandidatePoolSize;
    bool        quicklyDropUnpairedReads;
};
