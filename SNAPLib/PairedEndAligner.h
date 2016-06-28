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

#ifndef __PAIREDENDALIGNER_H__

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
    
    virtual bool align(
        Read                  *read0,
        Read                  *read1,
        PairedAlignmentResult *result,
        int                    maxEditDistanceForSecondaryResults,
        _int64                 secondaryResultBufferSize,
        _int64                *nSecondaryResults,
        PairedAlignmentResult *secondaryResults,             // The caller passes in a buffer of secondaryResultBufferSize and it's filled in by align()
        _int64                 singleSecondaryBufferSize,
        _int64                 maxSecondaryAlignmentsToReturn,
        _int64                *nSingleEndSecondaryResultsForFirstRead,
        _int64                *nSingleEndSecondaryResultsForSecondRead,
        SingleAlignmentResult *singleEndSecondaryResults     // Single-end secondary alignments for when the paired-end alignment didn't work properly
        ) = 0; // return value is false iff there wasn't enough space in the secondary buffers

    virtual void setLandauVishkin(
        LandauVishkin<1>        *landauVishkin,
        LandauVishkin<-1>       *reverseLandauVishkin)
    {
    }

    virtual _int64 getLocationsScored() const  = 0;
};

#endif // __PAIREDENDALIGNER_H__
