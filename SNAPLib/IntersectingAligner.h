/*++

Module Name:

    IntersectingAligner.h

Abstract:

    Header for SNAP genome aligner

Authors:

    Bill Bolosky, September, 2011

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
#include "LandauVishkin.h"

class IntersectingAligner: public Aligner {
public:

    IntersectingAligner(
        GenomeIndex *i_genomeIndex, 
        unsigned     i_confDiff, 
        unsigned     i_maxHitsToConsider, 
        unsigned     i_maxK,
        unsigned     i_maxReadSize,
        unsigned     i_maxSeedsToUse,
        unsigned     i_lvCutoff);

        AlignmentResult
    AlignRead(
        Read        *read,
        unsigned    *genomeLocation,
        bool        *hitIsRC,
        int         *score = NULL);

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
        fprintf(stderr,"ComputeHitDistribution doesn't work for IntersectingAligner.\n");
        exit(1);
    }
    

    _int64 getNHashTableLookups() const {return nHashTableLookups;}
    _int64 getLocationsScored() const {return nLocationsScored;}
    _int64 getNHitsIgnoredBecauseOfTooHighPopularity() const {return nHitsIgnoredBecauseOfTooHighPopularity;}
    _int64 getNReadsIgnoredBecauseOfTooManyNs() const {return nReadsIgnoredBecauseOfTooManyNs;}
    _int64 getNIndelsMerged() const {return nIndelsMerged;}

#if     MAINTAIN_HISTOGRAMS
    const Histogram *getLVHistogram() const {return lvHistogram;}
    const Histogram *getLookupHistogram() const {return lookupHistogram;}
    const Histogram *getLVHistogramForMulti() const {return lvHistogramForMulti;}
    const Histogram *getLVHistogramWhenBestFound() const {return lvCountWhenBestFound;}
#endif  // MAINTAIN_HISTOGRAMS


    const char *getRCTranslationTable() const {return rcTranslationTable;}

    ~IntersectingAligner();

    inline int getMaxK() const {return maxK;}
    inline int getConfDiff() const {return confDiff;}

    const char *getName() const {return "Intersecting Aligner";}

    struct LookupResult {
        unsigned        nHits;
        bool            isRC;
        const unsigned *hits;
        unsigned        offset; // The offset of the seed within the read.  This is corrected for during intersection.
    };

    void addIgnoredReads(_int64 newlyIgnoredReads) {};

private:
    LandauVishkin landauVishkin;

    char rcTranslationTable[256];


        AlignmentResult
    AlignReadInternal(
        unsigned    seedOffset,
        Read        *read,
        unsigned    *genomeLocation,
        bool        *hitIsRC,
        int         *score);

    static const unsigned maxIntersections = 4;

#if     MAINTAIN_HISTOGRAMS
    Histogram   *lvHistogram;
    Histogram   *lookupHistogram;
    Histogram   *lvHistogramForMulti;
    Histogram   *lvCountWhenBestFound;
#endif  // MAINTAIN_HISTOGRAMS

    _int64 nHashTableLookups;
    _int64 nLocationsScored;
    _int64 nHitsIgnoredBecauseOfTooHighPopularity;
    _int64 nReadsIgnoredBecauseOfTooManyNs;
    _int64 nIndelsMerged;

    //
    // A bitvector indexed by offset in the read indicating whether this seed is used.
    // This is here to avoid doing a memory allocation in the aligner.
    //
    BYTE *seedUsed;
    BYTE *seedUsedAsAllocated;  // Use this for deleting seedUsed.

    inline bool IsSeedUsed(unsigned indexInRead) const {
        return (seedUsed[indexInRead / 8] & (1 << (indexInRead % 8))) != 0;
    }

    inline void SetSeedUsed(unsigned indexInRead) {
        seedUsed[indexInRead / 8] |= (1 << (indexInRead % 8));
    }

    LookupResult *lookupResults;
    int nLookupResults;

    LookupResult *rcLookupResults;
    int nRCLookupResults;

    //
    // A big array of ulongs that we can use for intersection results.
    //
    unsigned *intersectionResultSpace;
    unsigned intersectionResultsSpaceSize;
    unsigned usedEntriesInIntersectionResultSpace;

    struct IntersectionResult {
        int             nMisses;        // How many results isn't it in?
        unsigned        nHits;          // What's the size of the intersection?
        bool            isRC;
        unsigned        *hits;
    };

    IntersectionResult *intersectionResults;
    unsigned nIntersectionResults;

    IntersectionResult *rcIntersectionResults;
    unsigned nRCIntersectionResults;

    //
    // This is really a local to computeIntersection, but we allocate it here to avoid having to
    // do a memory allocation/deallocation there.
    //
    int *intersectionOffsets;
    int *alreadyIntersectedOffsets;

    void
    computeIntersection(
        IntersectionResult      *alreadyIntersected,
        int                      nAlreadyIntersected,
        int                      totalSetsRepresentedByAlreadyIntersected,
        int                      maxMisses,
        LookupResult            *lookups,
        int                      nLookups,
        IntersectionResult      *results,
        bool                    *truncatedResult);  // Set iff we got more than maxHitsToConsider in any result set
        
    static const unsigned UnusedScoreValue = 0xffff;

    const Genome *genome;
    GenomeIndex *genomeIndex;
    unsigned seedLen;
    int confDiff;
    unsigned maxHitsToConsider;
    int maxK;
    int maxReadSize;
    int maxSeedsToUse; // Max number of seeds to look up in the hash table
    unsigned lvCutoff;

    char *rcReadData;

    unsigned nTable[256];

    inline unsigned getWrappedNextSeedToTest(unsigned wrapCount) {
            if (0 == wrapCount) {
                return 0;
            }
            switch (seedLen) {
                case 23: {
                    switch (wrapCount) {
                        case 1: return 12;
                        case 2: return 6;
                        case 3: return 17;
                        case 4: return 3;
                        case 5: return 9;
                        case 6: return 20;
                        case 7: return 14;
                        case 8: return 1;
                        case 9: return 4;
                        case 10: return 7;
                        case 11: return 10;
                        case 12: return 15;
                        case 13: return 18;
                        case 14: return 21;
                        case 15: return 4;
                        case 16: return 2;
                        case 17: return 5;
                        case 18: return 11;
                        case 19: return 16;
                        case 20: return 19;
                        case 21: return 22;
                        case 22: return 8;
                    }
                }
                case 22: {
                    switch (wrapCount) {
                        case 1: return 11;
                        case 2: return 6;
                        case 3: return 16;
                        case 4: return 3;
                        case 5: return 9;
                        case 6: return 14;
                        case 7: return 19;
                        case 8: return 2;
                        case 9: return 7;
                        case 10: return 12;
                        case 11: return 17;
                        case 12: return 20;
                        case 13: return 4;
                        case 14: return 1;
                        case 15: return 10;
                        case 16: return 13;
                        case 17: return 15;
                        case 18: return 18;
                        case 19: return 21;
                        case 20: return 5;
                        case 21: return 8;
                        default: _ASSERT(!"NOTREACHED");
                    }
                }

                case 21: {
                    switch (wrapCount) {
                        case 1: return 11;
                        case 2: return 6;
                        case 3: return 16;
                        case 4: return 3;
                        case 5: return 9;
                        case 6: return 13;
                        case 7: return 17;
                        case 8: return 18;
                        case 9: return 2;
                        case 10: return 5;
                        case 11: return 8;
                        case 12: return 15;
                        case 13: return 20;
                        case 14: return 1;
                        case 15: return 4;
                        case 16: return 7;
                        case 17: return 10;
                        case 18: return 12;
                        case 19: return 14;
                        case 20: return 19;
                        default: _ASSERT(!"NOTREACHED");
                    }
                }
                case 20: {
                    switch (wrapCount) {
                        case 1: return 10;
                        case 2: return 5;
                        case 3: return 15;
                        case 4: return 2;
                        case 5: return 7;
                        case 6: return 12;
                        case 7: return 17;
                        case 8: return 3;
                        case 9: return 9;
                        case 10: return 11;
                        case 11: return 13;
                        case 12: return 19;
                        case 13: return 1;
                        case 14: return 4;
                        case 15: return 6;
                        case 16: return 8;
                        case 17: return 14;
                        case 18: return 18;
                        case 19: return 16;
                        default: _ASSERT(!"NOTREACHED");
                    }
                }

                case 19: {
                    switch (wrapCount) {
                        case 1: return 10;
                        case 2: return 4;
                        case 3: return 14;
                        case 4: return 2;
                        case 5: return 6;
                        case 6: return 8;
                        case 7: return 12;
                        case 8: return 16;
                        case 9: return 18;
                        case 10: return 1;
                        case 11: return 3;
                        case 12: return 5;
                        case 13: return 7;
                        case 14: return 9;
                        case 15: return 11;
                        case 16: return 13;
                        case 17: return 15;
                        case 18: return 17;
                        default: _ASSERT(!"NOTREACHED");
                    }
                }

                case 18: {
                    switch (wrapCount) {
                        case 1: return 9;
                        case 2: return 4;
                        case 3: return 13;
                        case 4: return 2;
                        case 5: return 6;
                        case 6: return 11;
                        case 7: return 15;
                        case 8: return 1;
                        case 9: return 3;
                        case 10: return 5;
                        case 11: return 7;
                        case 12: return 8;
                        case 13: return 10;
                        case 14: return 12;
                        case 15: return 14;
                        case 16: return 16;
                        case 17: return 17;
                        default: _ASSERT(!"NOTREACHED");
                    }
                }

                case 17: {
                    switch (wrapCount) {
                        case 1: return 8;
                        case 2: return 4;
                        case 3: return 12;
                        case 4: return 2;
                        case 5: return 6;
                        case 6: return 10;
                        case 7: return 14;
                        case 8: return 1;
                        case 9: return 3;
                        case 10: return 5;
                        case 11: return 7;
                        case 12: return 9;
                        case 13: return 11;
                        case 14: return 13;
                        case 15: return 15;
                        case 16: return 16;
                        default: _ASSERT(!"NOTREACHED");
                    }
                }

                case 16: {
                    switch (wrapCount) {
                        case 1: return 8;
                        case 2: return 4;
                        case 3: return 12;
                        case 4: return 2;
                        case 5: return 6;
                        case 6: return 10;
                        case 7: return 14;
                        case 8: return 1;
                        case 9: return 3;
                        case 10: return 5;
                        case 11: return 7;
                        case 12: return 9;
                        case 13: return 11;
                        case 14: return 13;
                        case 15: return 15;
                        default: _ASSERT(!"NOTREACHED");
                    }
                } // inner switch
                default: fprintf(stderr,"IntersectingAligner: Not set up to run with this seed size\n"); exit(1);
            } // outer switch
        } 
};
