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
#include "Read.h"
#include "Seed.h"
#include "Genome.h"
#include "GenomeIndex.h"
#include "Histogram.h"
#include "Aligner.h"
#include "LandauVishkin.h"

#define USE_CAREFUL_BASE_STATUS     0   // Turning this off speeds up the aligner code path at the cost of a few possible extra lookups

class OldAligner: public Aligner {
public:

    OldAligner(
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
        int         *finalScore = NULL);

    //
    // Statistics gathering.
    //
        void
    ComputeHitDistribution(
        Read        *read,
        unsigned     correctGenomeLocation,
        bool         correctHitIsRC,
        unsigned    *hitCountBySeed,
        unsigned    *rcHitCountBySeed,
        unsigned    &nSeedsApplied,
        unsigned    &nRCSeedsApplied,
        unsigned    *hitsCountsContainingCorrectLocation);

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

    ~OldAligner();

    inline int getMaxK() const {return maxK;}
    inline int getConfDiff() const {return confDiff;}

    const char *getName() const {return "Old Aligner";}

    void addIgnoredReads(_int64 newlyIgnoredReads) {}

private:

    LandauVishkin landauVishkin;

    char rcTranslationTable[256];

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

#if     USE_CAREFUL_BASE_STATUS
    struct BaseStatus {
        unsigned        nSeedsIncludingThisBase;
        unsigned        nRCSeedsIncludingThisBase;
    };
    //
    // A table indexed by offset in the Read that describes what we've done with the (potential) seed that
    // starts at this location.
    //
    BaseStatus *bases;  // This is here just to avoid dynamically allocating this on every call to AlignRead
#endif  // USE_CAREFUL_BASE_STATUS

    struct Candidate {
        unsigned        genomeLocation;
        unsigned        weight;
        unsigned        minPossibleScore;
        //
        // In cases of insertions and deletions we might see
        // several nearby hits, because the insertion or deletion
        // messes up re-aligning the seed match with the read.  If that
        // happens, then we just add the new seed to the existing one
        // by updating "minRange" and "maxRange" to include the new
        // hit, and than scoring in that range.
        unsigned        minRange;
        unsigned        maxRange;
        unsigned        score;
        bool            scored;
        bool            isRC;

        inline bool operator>(OldAligner::Candidate &peer) {
            if (genomeLocation == peer.genomeLocation) {
                return isRC > peer.isRC;
            }
            return genomeLocation > peer.genomeLocation;
        }

        inline bool operator>=(OldAligner::Candidate &peer) {
            if (genomeLocation == peer.genomeLocation) {
                return isRC >= peer.isRC;
            }
            return genomeLocation > peer.genomeLocation;
        }

        inline bool operator<(OldAligner::Candidate &peer) {
            if (genomeLocation == peer.genomeLocation) {
                return isRC < peer.isRC;
            }
            return genomeLocation < peer.genomeLocation;
        }

        inline bool operator<=(OldAligner::Candidate &peer) {
            if (genomeLocation == peer.genomeLocation) {
                return isRC <= peer.isRC;
            }
            return genomeLocation < peer.genomeLocation;
        }

        inline bool operator==(OldAligner::Candidate &peer) {
            return genomeLocation == peer.genomeLocation && isRC == peer.isRC;
        }

        inline bool operator!=(OldAligner::Candidate &peer) {
            return genomeLocation != peer.genomeLocation || isRC != peer.isRC;
        }
    };

    static const unsigned UnusedScoreValue = 0xffff;

        bool
    score(
        bool             forceResult,
        Read            *read,
        Read            *rcRead,
        AlignmentResult *result,
        int             *finalScore,
        unsigned        *singleHitGenomeLocation,
        bool            *hitIsRC,
        unsigned         nSeedsApplied,
        unsigned         nRCSeedsApplied,
        unsigned         mostSeedsContainingAnyParticularBase,
        unsigned         mostRCSeedsContainingAnyParticularBase,
        unsigned        &lowestPossibleScoreOfAnyUnseenLocation,
        unsigned        &lowestPossibleRCScoreOfAnyUnseenLocation,
        Candidate       *candidates,
        unsigned        &bestScore,
        unsigned        &bestScoreGenomeOffset,
        unsigned        &secondBestScore,
        unsigned        &lvScores,
        unsigned        &lvScoresAfterBestFound,
        bool             anyHitsSkippedBecauseOfTooHighPopularityThisRun);
    
    inline void clearCandidates() {nCandidates = 0;}
    Candidate *findCandidate(unsigned genomeLocation, bool isRC);
    Candidate *allocateNewCandidate(unsigned genomeLocation, bool isRC);
    void removeCandidateAtIndex(int index);

    const Genome *genome;
    GenomeIndex *genomeIndex;
    unsigned seedLen;
    unsigned confDiff;
    unsigned maxHitsToConsider;
    unsigned maxK;
    unsigned maxReadSize;
    unsigned maxSeedsToUse; // Max number of seeds to look up in the hash table
    unsigned lvCutoff;

    int nCandidates;
    Candidate *candidates;

    char *rcReadData;

    unsigned nTable[256];

#if     DBG
    void AssertCandidatesAreSorted();
#else   // DBG
    inline void AssertCandidatesAreSorted() {}
#endif  // DBG

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
                default: fprintf(stderr,"OldAligner: Not set up to run with this seed size\n"); exit(1);
            } // outer switch
        } 
};
