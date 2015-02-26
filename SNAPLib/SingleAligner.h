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
#include "AlignmentResult.h"

class SingleAlignerContext : public AlignerContext
{
public:

    SingleAlignerContext(AlignerExtension* i_extension = NULL);
    

protected:

    // AlignerContext overrides

    virtual AlignerStats* newStats();
    
    virtual void runTask();
    
    virtual void runIterationThread();

    virtual void typeSpecificBeginIteration();
    virtual void typeSpecificNextIteration();

    // for subclasses

    virtual void updateStats(AlignerStats* stats, Read* read, AlignmentResult result, int score, int mapq);

    //RangeSplittingReadSupplierGenerator   *readSupplierGenerator;

    ReadSupplierGenerator *readSupplierGenerator;

	friend class AlignerContext2;

    bool isPaired() {return false;}
};
