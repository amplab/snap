/*++

Module Name:

    SingleAligner.cpp

Abstract:

    Functions for running the single end aligner sub-program.

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
#include "AlignerStats.h"
#include "ReadSupplierQueue.h"

class SingleAlignerContext : public AlignerContext
{
public:

    SingleAlignerContext(AlignerExtension* i_extension = NULL);
    

protected:

    // AlignerContext overrides

    virtual AlignerOptions* parseOptions(int i_argc, const char **i_argv, const char *i_version, unsigned *argsUsed);

    virtual AlignerStats* newStats();
    
    virtual void runTask();
    
    virtual void runIterationThread();

    virtual void typeSpecificBeginIteration();
    virtual void typeSpecificNextIteration();

    // for subclasses

    virtual void writeRead(Read* read, AlignmentResult result, unsigned location, Direction direction, bool isTranscriptome, unsigned tlocation, int score, int mapq);

    virtual void updateStats(AlignerStats* stats, Read* read, AlignmentResult result, unsigned location, int score, int mapq, bool wasError);

    //RangeSplittingReadSupplierGenerator   *readSupplierGenerator;

    ReadSupplierGenerator *readSupplierGenerator;

	friend class AlignerContext2;
};
