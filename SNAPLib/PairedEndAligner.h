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
#include "AffineGap.h"
#include "AffineGapVectorized.h"
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
        PairedAlignmentResult* result,
        PairedAlignmentResult* firstALTResult,
        int                    maxEditDistanceForSecondaryResults,
        _int64                 secondaryResultBufferSize,
        _int64                *nSecondaryResults,
        PairedAlignmentResult *secondaryResults,             // The caller passes in a buffer of secondaryResultBufferSize and it's filled in by align()
        _int64                 singleSecondaryBufferSize,
        _int64                 maxSecondaryAlignmentsToReturn,
        _int64                *nSingleEndSecondaryResultsForFirstRead,
        _int64                *nSingleEndSecondaryResultsForSecondRead,
        SingleAlignmentResult *singleEndSecondaryResults,     // Single-end secondary alignments for when the paired-end alignment didn't work properly
        _int64                 maxPairedCandidatesForAffineGapBufferSize,
        _int64                *nPairedCandidatesForAffineGap,
        PairedAlignmentResult *pairedCandidatesForAffineGap,
        _int64                 maxSingleCandidatesForAffineGapBufferSize,
        _int64                *nSingleCandidatesForAffineGapFirstRead,
        _int64                *nSingleCandidatesForAffineGapSecondRead,
        SingleAlignmentResult *singleCandidatesForAffineGap,
        int                    maxK
    ) = 0; // return value is false iff there wasn't enough space in the secondary buffers

    virtual void setLandauVishkin(
        LandauVishkin<1>        *landauVishkin,
        LandauVishkin<-1>       *reverseLandauVishkin)
    {
    }

    virtual void setAffineGap(
        AffineGapVectorized<1>        *affineGap,
        AffineGapVectorized<-1>       *reverseAffineGap)
        // AffineGap<1>        *affineGap,
        // AffineGap<-1>       *reverseAffineGap)
    {
    }

    virtual _int64 getLocationsScored() const  = 0;
};
