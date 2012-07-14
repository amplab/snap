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

#include "Aligner.h"

struct PairedAlignmentResult {
    AlignmentResult status[2];  // SingleHit or CertainHit if aligned, MultipleHit if matches DB
                                // but not confidently aligned, or NotFound.

    unsigned location[2];       // Genome location of each read.
    
    bool isRC[2];               // Did we match the reverse complement? In general the two reads should have
                                // opposite orientations because they're part of the same original fragment,
                                // but it seems possible for a piece of the genome to get cut cleanly and flip
                                // in a translocation event, which would cause both ends of a fragment aligning
                                // there to be in the same orientation w.r.t. the reference genome.

    int score[2];               // score of each end if matched
};

/**
 * Abstract interface for paired-end aligners.
 */
class PairedEndAligner
{
public:
    virtual ~PairedEndAligner();
    
    virtual void align(
        Read                  *read0,
        Read                  *read1,
        PairedAlignmentResult *result) = 0;
};
