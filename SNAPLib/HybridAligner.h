/*++

Module Name:

    HybridAligner.h

Abstract:

    Header for SNAP genome aligner

Authors:

    Bill Bolosky, October, 2011

Environment:

    User mode service.

    This class is NOT thread safe.  It's the caller's responsibility to ensure that
    at most one thread uses an instance at any time.

Revision History:

    Adapted from Matei Zaharia's Scala implementation.

--*/

#pragma once

#include "options.h"
#include "Read.h"
#include "Seed.h"
#include "Genome.h"
#include "GenomeIndex.h"
#include "Histogram.h"
#include "Aligner.h"
#include "IntersectingAligner.h"
#include "BaseAligner.h"

class HybridAligner: public Aligner {
public:

    HybridAligner(
        GenomeIndex *genomeIndex, 
        unsigned     confDiff, 
        unsigned     maxHitsToConsider, 
        unsigned     maxK,
        unsigned     maxReadSize,
        unsigned     maxSeedsToUse,
        unsigned     lvCutoff);

        AlignmentResult
    AlignRead(
        Read        *read,
        unsigned    *genomeLocation,
        bool        *hitIsRC,
        int         *score = NULL);

    _int64 getNReadsIgnoredBecauseOfTooManyNs() const {return baseAligner->getNReadsIgnoredBecauseOfTooManyNs();}
    int getMaxK() const {return baseAligner->getMaxK();}
    const char *getRCTranslationTable() const {return baseAligner->getRCTranslationTable();}
    int getConfDiff() const {return baseAligner->getConfDiff();}
    _int64 getNHashTableLookups() {return baseAligner->getNHashTableLookups() + intersectingAligner->getNHashTableLookups();}
    _int64 getLocationsScored() {return baseAligner->getLocationsScored() + intersectingAligner->getLocationsScored();}
    _int64 getNHitsIgnoredBecauseOfTooHighPopularity() const {return baseAligner->getNHitsIgnoredBecauseOfTooHighPopularity();}
    _int64 getNIndelsMerged() const {return baseAligner->getNIndelsMerged() + intersectingAligner->getNIndelsMerged();}
    _int64 getNHashTableLookups() const {return baseAligner->getNHashTableLookups() + intersectingAligner->getNHashTableLookups();}
    _int64 getLocationsScored() const {return baseAligner->getLocationsScored() + intersectingAligner->getLocationsScored();}
    void addIgnoredReads(_int64 newlyIgnoredReads) {}

        void
    ComputeHitDistribution(
        Read        *read,
        unsigned     correctGenomeLocation,
        bool         correctHitIsRC,
        unsigned    *hitCountBySeed,
        unsigned    *rcHitCountBySeed,
        unsigned    &nSeedsApplied,
        unsigned    &nRCSeedsApplied,
        unsigned    *hitsCountsContainingCorrectLocation)
    {
        return baseAligner->ComputeHitDistribution(read,correctGenomeLocation,correctHitIsRC,hitCountBySeed,rcHitCountBySeed,nSeedsApplied,nRCSeedsApplied,
                            hitsCountsContainingCorrectLocation);
    }

    const char *getName() const {return "Hybrid Aligner";}

private:

    BaseAligner             *baseAligner;
    IntersectingAligner     *intersectingAligner;
};



