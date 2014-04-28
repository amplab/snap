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

struct SingleAlignmentResult {
    AlignmentResult status;  

    unsigned location;          // Aligned genome location.
    Direction direction;        // Did we match the reverse complement? 
    int score;                  // score of each end if matched

    int mapq;                   // mapping quality, encoded like a Phred score (but as an integer, not ASCII Phred + 33).
    bool isTranscriptome;       //is this alignment a transcriptome alignment or genome alignment
    unsigned tlocation;         //pos of the transcriptome alignment
};

// Does an AlignmentResult represent a single location?
inline bool isOneLocation(AlignmentResult result) {
    return result == SingleHit;
}

extern bool doAlignerPrefetch;

class Aligner {
    public:
   
    virtual ~Aligner() {}

        virtual void
    AlignRead(
        Read                    *read,
        SingleAlignmentResult   *primaryResult,
        int                      maxEditDistanceForSecondaryResults,
        int                      secondaryResultBufferSize,
        int                     *nSecondaryResults,
        SingleAlignmentResult   *secondaryResults             // The caller passes in a buffer of secondaryResultBufferSize and it's filled in by AlignRead()
    ) = 0;      // Retun value is true if there was enough room in the secondary alignment buffer for everything that was found.

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
