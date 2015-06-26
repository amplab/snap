/*++

Module Name:

    PairedEndAligner.h

Abstract:

    Superclass for implementations of paired end aligners

Authors:

    Matei Zaharia, December, 2011

Environment:

    User mode service.

Revision History:

--*/

#pragma once

#include "AlignmentResult.h"
#include "directions.h"
#include "LandauVishkin.h"
#include "Read.h"



/**
 * Abstract interface for paired-end aligners.
 */
class PairedEndAligner
{
public:
    virtual ~PairedEndAligner() {}
    
    virtual void align(
        Read                  *read0,
        Read                  *read1,
        PairedAlignmentResult *result,
        int                    maxEditDistanceForSecondaryResults,
        int                    secondaryResultBufferSize,
        int                   *nSecondaryResults,
        PairedAlignmentResult *secondaryResults,             // The caller passes in a buffer of secondaryResultBufferSize and it's filled in by align()
        int                    singleSecondaryBufferSize,
        int                    maxSecondaryAlignmentsToReturn,
        int                   *nSingleEndSecondaryResultsForFirstRead,
        int                   *nSingleEndSecondaryResultsForSecondRead,
        SingleAlignmentResult *singleEndSecondaryResults     // Single-end secondary alignments for when the paired-end alignment didn't work properly
        ) = 0;

    virtual void setLandauVishkin(
        LandauVishkin<1>        *landauVishkin,
        LandauVishkin<-1>       *reverseLandauVishkin)
    {
    }

    virtual _int64 getLocationsScored() const  = 0;
};
