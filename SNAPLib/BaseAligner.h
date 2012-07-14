/*++

Module Name:

    BaseAligner.h

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

#include "Aligner.h"
#include "LandauVishkin.h"
#include "BoundedStringDistance.h"
#include "BigAlloc.h"

class BaseAligner: public Aligner {
public:

    BaseAligner(
        GenomeIndex    *i_genomeIndex, 
        unsigned        i_confDiff, 
        unsigned        i_maxHitsToConsider, 
        unsigned        i_maxK,
        unsigned        i_maxReadSize,
        unsigned        i_maxSeedsToUse,
        unsigned        i_lvCutoff,
        unsigned        i_adaptiveConfDiffThreshold,
        LandauVishkin  *i_landauVishkin = NULL);

    virtual ~BaseAligner();

        virtual AlignmentResult
    AlignRead(
        Read        *read,
        unsigned    *genomeLocation,
        bool        *hitIsRC,
        int         *finalScore = NULL);
        
    //
    // A richer version of AlignRead that allows for searching near a given location.
    // If searchRadius is not 0, constrain the search to distance maxSpacing around
    // searchLocation in the orientation given by searchRC.
    //
        AlignmentResult
    AlignRead(
        Read        *read,
        unsigned    *genomeLocation,
        bool        *hitIsRC,
        int         *finalScore,
        unsigned     searchRadius,       // If non-zero, constrain search around searchLocation in direction searchRC.
        unsigned     searchLocation,
        bool         searchRC);
        
    //
    // A richer version of AlignRead that allows for searching near a given location, as well as returning
    // multiple hits if the best hits are within distance confDiff of each other.
    //
        AlignmentResult
    AlignRead(
        Read        *read,
        unsigned    *genomeLocation,
        bool        *hitIsRC,
        int         *finalScore,
        unsigned     searchRadius,       // If non-zero, constrain search around searchLocation in direction searchRC.
        unsigned     searchLocation,
        bool         searchRC,
        int          maxHitsToGet,       // If maxHitsToGet > 1, output up to this many hits within confDiff of the best
        int         *multiHitsFound,     // inside multiHitLocations / RCs instead of returning MultipleHits right away.
        unsigned    *multiHitLocations,
        bool        *multiHitRCs,
        int         *multiHitScores);

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
    void addIgnoredReads(_int64 newlyIgnoredReads) {nReadsIgnoredBecauseOfTooManyNs += newlyIgnoredReads;}

#if     MAINTAIN_HISTOGRAMS
    const Histogram *getLVHistogram() const {return lvHistogram;}
    const Histogram *getLookupHistogram() const {return lookupHistogram;}
    const Histogram *getLVHistogramForMulti() const {return lvHistogramForMulti;}
    const Histogram *getLVHistogramWhenBestFound() const {return lvCountWhenBestFound;}
#endif  // MAINTAIN_HISTOGRAMS


    const char *getRCTranslationTable() const {return rcTranslationTable;}

    inline int getMaxK() const {return maxK;}
    inline int getConfDiff() const {return confDiff;}

    inline void setMaxK(int maxK_) {maxK = maxK_;}

    inline void setReadId(int readId_) {readId = readId_;}

    const char *getName() const {return "Base Aligner";}

    inline bool checkedAllSeeds() {return popularSeedsSkipped == 0;}

    void *operator new(size_t size) {return BigAlloc(size);}
    void operator delete(void *ptr) {BigDealloc(ptr);}

    inline bool getExplorePopularSeeds() {return explorePopularSeeds;}
    inline void setExplorePopularSeeds(bool newValue) {explorePopularSeeds = newValue;}

    inline bool getStopOnFirstHit() {return stopOnFirstHit;}
    inline void setStopOnFirstHit(bool newValue) {stopOnFirstHit = newValue;}

private:
    LandauVishkin *landauVishkin;
    BoundedStringDistance<> *bsd;
    bool ownLandauVishkin;

    // Maximum distance to merge candidates that differe in indels over.
    // This can't be bigger than 32, else some bitvectors overflow.
    // TODO(matei): this seems to work better when we make it lower; why?
    static const unsigned maxMergeDist = 15; 

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

    struct Candidate {
        Candidate() {init();}
        void init();

        _int64          scoredInEpoch;
        unsigned        score;
    };

    struct HashTableElement {
        HashTableElement();
        void init();

        //
        // Doubly linked list for the weight buckets.
        //
        HashTableElement    *weightNext;
        HashTableElement    *weightPrev;

        //
        // Singly linked list for the hash table buckets.
        //
        HashTableElement    *next;

        _uint64      candidatesUsed;

        unsigned             baseGenomeLocation;
        unsigned             weight;
        unsigned             lowestPossibleScore;
        unsigned             bestScore;
        bool                 isRC;
        bool                 allExtantCandidatesScored;

        Candidate            candidates[maxMergeDist * 2];
    };

    //
    // Clearing out all of the pointers in the hash tables is expensive relative to running
    // an alignment, because usually the table is much bigger than the number of entries in it.
    // So, we avoid that expense by simply not clearing out the table at all.  Instead, along with
    // the pointers we keep an epoch number.  There's a corresponding epoch number in the
    // BaseAligner object, and if the two differ then the hash table bucket is empty.  We increment
    // the epoch number in the BaseAligner at the beginning of each alignment, thus effectively
    // clearing the hash table from the last run.
    //
    struct HashTableAnchor {
        HashTableElement *element;
        _int64            epoch;
    };

    _int64 hashTableEpoch;

    unsigned nUsedHashTableElements;
    unsigned hashTableElementPoolSize;
    HashTableElement *hashTableElementPool;

    const HashTableElement emptyHashTableElement;

    unsigned candidateHashTablesSize;
    HashTableAnchor *candidateHashTable[2]; // 0 is normal, 1 reverse complement
    
    HashTableElement *weightLists;
    unsigned highestUsedWeightList;

    static inline unsigned hash(unsigned key) {
#if     1
        key = key * 131;    // Believe it or not, we spend a long time computing the hash, so we're better off with more table entries and a dopey function.
#else   // 1
        //
        // Hash the key.  Use the hash finalizer from the 64 bit MurmurHash3, http://code.google.com/p/smhasher/wiki/MurmurHash3,
        // which is public domain code.
        //
    
        key ^= key >> 16; 
        key *= 0x85ebca6b; 
        key ^= key >> 13; 
        key *= 0xc2b2ae35; 
        key ^= key >> 16;
#endif  // 1
        return key;
    }


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
        unsigned        &bestScoreGenomeLocation,
        unsigned        &secondBestScore,
        unsigned        &secondBestScoreGenomeLocation,
        bool            &secondBestScoreIsRC,
        unsigned        &scoreLimit,
        unsigned        &lvScores,
        unsigned        &lvScoresAfterBestFound,
        unsigned         maxHitsToGet);
    
    void clearCandidates();

    void findCandidate(unsigned genomeLocation, bool isRC, Candidate **candidate, HashTableElement **hashTableElement);
    void allocateNewCandidate(unsigned genomeLoation, bool isRC, unsigned lowestPossibleScore, Candidate **candidate, HashTableElement **hashTableElement);
    void incrementWeight(HashTableElement *element);

    void fillHitsFound(unsigned maxHitsToGet, int *multiHitsFound, 
                       unsigned *multiHitLocations, bool *multiHitRCs, int *multiHitScores);

    const Genome *genome;
    GenomeIndex *genomeIndex;
    unsigned seedLen;
    unsigned confDiff;
    unsigned maxHitsToConsider;
    unsigned maxK;
    unsigned maxReadSize;
    unsigned maxSeedsToUse; // Max number of seeds to look up in the hash table
    unsigned lvCutoff;
    unsigned adaptiveConfDiffThreshold; // Increase confDiff by 1 if this many seeds are repetitive.

    int nCandidates;
    Candidate *candidates;

    char *rcReadData;

    unsigned nTable[256];

    int readId;
    
    // Store the best hits at a given edit distance, as well as their number
    static const int MAX_MULTI_HITS_TO_GET = 512;
    unsigned hitCount[MAX_K];
    unsigned hitLocations[MAX_K][MAX_MULTI_HITS_TO_GET];
    bool hitRCs[MAX_K][MAX_MULTI_HITS_TO_GET];

    // How many overly popular (> maxHits) seeds we skipped this run
    unsigned popularSeedsSkipped;

    bool explorePopularSeeds; // Whether we should explore the first maxHits hits even for overly
                              // popular seeds (useful for filtering reads that come from a database
                              // with many very similar sequences).

    bool stopOnFirstHit;      // Whether to stop the first time a location matches with less than
                              // maxK edit distance (useful when using SNAP for filtering only).

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
                default: fprintf(stderr,"BaseAligner: Not set up to run with this seed size\n"); exit(1);
            } // outer switch
        } 
};
