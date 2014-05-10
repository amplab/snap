/*++

Module Name:

    AlignmentResult.h

Abstract:

    Header for SNAP genome alignment results

Authors:

    Bill Bolosky, May, 2014

Environment:

    User mode service.

    This class is NOT thread safe.  It's the caller's responsibility to ensure that
    at most one thread uses an instance at any time.

Revision History:

   Pulled out of other places in SNAP

--*/

#pragma once
#include "Genome.h"
#include "directions.h"

enum AlignmentResult {NotFound, SingleHit, MultipleHits, UnknownAlignment}; // BB: Changed Unknown to UnknownAlignment because of a conflict w/Windows headers

bool isAValidAlignmentResult(AlignmentResult result);


inline const char *AlignmentResultToString(AlignmentResult result) {
    switch (result) {
        case NotFound: return "NotFound";
        case SingleHit: return "SingleHit";
        case MultipleHits: return "MultipleHits";
        case UnknownAlignment: return "Unknown";
        default: return "Unknown alignment result type";
    }
}

struct SingleAlignmentResult {
    AlignmentResult status;  

    GenomeLocation  location;       // Aligned genome location.
    Direction       direction;      // Did we match the reverse complement? 
    int             score;          // score of each end if matched

    int             mapq;           // mapping quality, encoded like a Phred score (but as an integer, not ASCII Phred + 33).
};

// Does an AlignmentResult represent a single location?
inline bool isOneLocation(AlignmentResult result) {
    return result == SingleHit;
}
