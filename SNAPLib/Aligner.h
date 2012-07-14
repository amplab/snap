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

enum AlignmentResult {NotFound, CertainHit, SingleHit, MultipleHits, UnknownAlignment}; // BB: Changed Unknown to UnknownAlignment because of a conflict w/Windows headers

inline const char *AlignmentResultToString(AlignmentResult result) {
    switch (result) {
        case NotFound: return "NotFound";
        case CertainHit: return "SingleHit";    // Just for now
        case SingleHit: return "SingleHit";
        case MultipleHits: return "MultipleHits";
        case UnknownAlignment: return "Unknown";
        default: return "Unknown alignment result type";
    }
}

// Does an AlignmentResult represent a single location?
inline bool isOneLocation(AlignmentResult result) {
    return result == SingleHit || result == CertainHit;
}

extern bool doAlignerPrefetch;

class Aligner {
    public:
   
    virtual ~Aligner() {}

        virtual AlignmentResult
    AlignRead(
        Read        *read,
        unsigned    *genomeLocation,
        bool        *hitIsRC,
        int         *finalScore = NULL) = 0;

        virtual void
    ComputeHitDistribution(
        Read        *read,
        unsigned     correctGenomeLocation,
        bool         correctHitIsRC,
        unsigned    *hitCountBySeed,
        unsigned    *rcHitCountBySeed,
        unsigned    &nSeedsApplied,
        unsigned    &nRCSeedsApplied,
        unsigned    *hitsCountsContainingCorrectLocation) = 0;

    virtual _int64 getNHashTableLookups() const = 0;
    virtual _int64 getLocationsScored() const  = 0;
    virtual _int64 getNHitsIgnoredBecauseOfTooHighPopularity() const = 0;
    virtual _int64 getNReadsIgnoredBecauseOfTooManyNs() const = 0;
    virtual _int64 getNIndelsMerged() const = 0;

    virtual void addIgnoredReads(_int64 newlyIgnoredReads) = 0;

#if     MAINTAIN_HISTOGRAMS
    virtual const Histogram *getLVHistogram() const  = 0;
    virtual const Histogram *getLookupHistogram() const = 0;
    virtual const Histogram *getLVHistogramForMulti() const = 0;
    virtual const Histogram *getLVHistogramWhenBestFound() const = 0;
#endif  // MAINTAIN_HISTOGRAMS


    virtual const char *getRCTranslationTable() const = 0;

    virtual int getMaxK() const = 0;
    virtual int getConfDiff() const = 0;

    virtual const char *getName() const = 0;

};
