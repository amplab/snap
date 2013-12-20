/*++

Module Name:

    Aligner.h

Abstract:

    Header for SNAP genome aligner

Authors:

    Bill Bolosky, August, 2011

Environment:

    User mode service.

    This class is NOT thread safe.  It's the caller's responsibility to ensure that
    at most one thread uses an instance at any time.

Revision History:

    Adapted from Matei Zaharia's Scala implementation.

--*/

#pragma once
#include "options.h"
#include "Genome.h"
#include "GenomeIndex.h"
#include "Histogram.h"
#include "Read.h"
#include "Seed.h"
#include "directions.h"

inline const char *AlignmentResultToString(AlignmentResult result) {
    switch (result) {
        case NotFound: return "NotFound";
        case SingleHit: return "SingleHit";
        case MultipleHits: return "MultipleHits";
        case UnknownAlignment: return "Unknown";
        default: return "Unknown alignment result type";
    }
}

// Does an AlignmentResult represent a single location?
inline bool isOneLocation(AlignmentResult result) {
    return result == SingleHit;
}

extern bool doAlignerPrefetch;

class Aligner {
    public:
   
    virtual ~Aligner() {}

        virtual AlignmentResult
    AlignRead(
        Read        *read,
        unsigned    *genomeLocation,
        Direction   *hitDirection,
        int         *finalScore = NULL,
        int         *mapq = NULL) = 0;

    virtual _int64 getNHashTableLookups() const = 0;
    virtual _int64 getLocationsScored() const  = 0;
    virtual _int64 getNHitsIgnoredBecauseOfTooHighPopularity() const = 0;
    virtual _int64 getNReadsIgnoredBecauseOfTooManyNs() const = 0;
    virtual _int64 getNIndelsMerged() const = 0;

    virtual void addIgnoredReads(_int64 newlyIgnoredReads) = 0;

    virtual const char *getRCTranslationTable() const = 0;

    virtual int getMaxK() const = 0;

    virtual const char *getName() const = 0;

};
