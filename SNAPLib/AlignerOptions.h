/*++

Module Name:

    AlignerOptions.h

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
#include "options.h"
#include "Range.h"
#include "Genome.h"
#include "SAM.h"
#include "RangeSplitter.h"

struct AbstractOptions
{
    virtual void usageMessage() = 0;

    virtual bool parse(char** argv, int argc, int& n) = 0;
};

struct AlignerOptions : public AbstractOptions
{
    AlignerOptions(const char* i_commandLine, bool forPairedEnd = false);

    const char         *commandLine;
    const char         *indexDir;
    int                 numThreads;
    Range               maxDist;
    Range               numSeeds;
    Range               maxHits;
    Range               confDiff;
    Range               adaptiveConfDiff;
    bool                computeError;
    bool                bindToProcessors;
    bool                ignoreMismatchedIDs;
    unsigned            selectivity;
    char               *samFileTemplate;
    bool                doAlignerPrefetch;
    char               *inputFilename;
    bool                inputFileIsFASTQ;   // Else SAM
    ReadClippingType    clipping;
    AbstractOptions    *extra; // extra options

    void usage();

    virtual void usageMessage();

    virtual bool parse(char** argv, int argc, int& n);
};
