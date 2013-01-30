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
#include "ProbabilityDistance.h"
#include "SimilarityMap.h"
#include "AlignerStats.h"
#include "directions.h"

class BaseAligner: public Aligner {
public:

    BaseAligner(
        GenomeIndex    *i_genomeIndex, 
        unsigned        i_confDiff, 
        unsigned        i_maxHitsToConsider, 
        unsigned        i_maxK,
        unsigned        i_maxReadSize,
        unsigned        i_maxSeedsToUse,
        unsigned        i_adaptiveConfDiffThreshold,
        LandauVishkin<>*i_landauVishkin = NULL,
        SimilarityMap  *i_similarityMap = NULL,
        AlignerStats   *i_stats = NULL,
        BigAllocator    *allocator = NULL);

    virtual ~BaseAligner();

        virtual AlignmentResult
    AlignRead(
        Read        *inputRead,
        unsigned    *genomeLocation,
        Direction   *hitDirection,
        int         *finalScore = NULL,
        int         *mapq = NULL);
        
    //
    // A richer version of AlignRead that allows for searching near a given location.
    // If searchRadius is not 0, constrain the search to distance maxSpacing around
    // searchLocation in the orientation given by searchRC.
    //
        AlignmentResult
    AlignRead(
        Read        *inputRead,
        unsigned    *genomeLocation,
        Direction   *hitDirection,
        int         *finalScore,
        int         *mapq,
        unsigned     searchRadius,       // If non-zero, constrain search around searchLocation in direction searchRC.
        unsigned     searchLocation,
        Direction    searchDirection);
        
    //
    // A richer version of AlignRead that allows for searching near a given location, as well as returning
    // multiple hits if the best hits are within distance confDiff of each other.
    //
        AlignmentResult
    AlignRead(
        Read        *inputRead,
        unsigned    *genomeLocation,
        Direction   *hitDirection,
        int         *finalScore,
        int         *mapq,
        unsigned     searchRadius,       // If non-zero, constrain search around searchLocation in direction searchRC.
        unsigned     searchLocation,
        Direction    searchDirection,
        int          maxHitsToGet,       // If maxHitsToGet > 1, output up to this many hits within confDiff of the best
        int         *multiHitsFound,     // inside multiHitLocations / RCs instead of returning MultipleHits right away.
        unsigned    *multiHitLocations,
        Direction   *multiHitDirections,
        int         *multiHitScores);

    //
    // Statistics gathering.
    //
        void
    ComputeHitDistribution(
        Read        *read,
        unsigned     correctGenomeLocation,
        Direction    correctHitDirection,
        unsigned    *hitCountBySeed[NUM_DIRECTIONS],
        unsigned    *nSeedsApplied[NUM_DIRECTIONS],
        unsigned    *hitsContainingCorrectLocation);

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

    void *operator new(size_t size, BigAllocator *allocator) {_ASSERT(size == sizeof(BaseAligner)); return allocator->allocate(size);}
    void operator delete(void * ptr, BigAllocator *allocator) {/* do nothing; the owner of the BigAlloctor is responsible for cleanup*/}

    inline bool getExplorePopularSeeds() {return explorePopularSeeds;}
    inline void setExplorePopularSeeds(bool newValue) {explorePopularSeeds = newValue;}

    inline bool getStopOnFirstHit() {return stopOnFirstHit;}
    inline void setStopOnFirstHit(bool newValue) {stopOnFirstHit = newValue;}

    static size_t getBigAllocatorReservation(bool ownLandauVishkin, unsigned maxHitsToConsider, unsigned maxReadSize, unsigned seedLen, unsigned maxSeedsToUse);

private:

    bool hadBigAllocator;

#if     defined(USE_BOUNDED_STRING_DISTANCE)
    BoundedStringDistance<> *boundedStringDist;
#endif  // bsd or LV

    LandauVishkin<> *landauVishkin;
    LandauVishkin<-1> *reverseLandauVishkin;
    bool ownLandauVishkin;

    ProbabilityDistance *probDistance;

    // Maximum distance to merge candidates that differe in indels over.
    // This can't be bigger than 32, else some bitvectors overflow.
    // TODO(matei): this seems to work better when we make it lower; why?
    static const unsigned maxMergeDist = 31; 

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

        unsigned        score;
        int             seedOffset;
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

        _uint64             candidatesUsed;    // Really candidates we still need to score
        _uint64             candidatesScored;

        unsigned             baseGenomeLocation;
        unsigned             weight;
        unsigned             lowestPossibleScore;
        unsigned             bestScore;
        Direction            direction;
        bool                 allExtantCandidatesScored;
        double               matchProbabilityForBestScore;

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
    HashTableAnchor *candidateHashTable[NUM_DIRECTIONS];
    
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

    // MAPQ parameters, currently not set to match Mason.  Using #define because VC won't allow "static const double".
#define SNP_PROB  0.001
#define GAP_OPEN_PROB  0.001
#define GAP_EXTEND_PROB  0.5

    //
    // Storage that's used during a call to AlignRead, but that's also needed by the
    // score function.  Since BaseAligner is single threaded, it's easier just to make
    // them member variables than to pass them around.
    //
    unsigned lowestPossibleScoreOfAnyUnseenLocation[NUM_DIRECTIONS];
    unsigned mostSeedsContainingAnyParticularBase[NUM_DIRECTIONS];
    unsigned nSeedsApplied[NUM_DIRECTIONS];
    unsigned bestScore;
    unsigned bestScoreGenomeLocation;
    unsigned secondBestScore;
    unsigned secondBestScoreGenomeLocation;
    int      secondBestScoreDirection;
    unsigned scoreLimit;
    unsigned lvScores;
    unsigned lvScoresAfterBestFound;
    double probabilityOfAllCandidates;
    double probabilityOfBestCandidate;
    int firstPassSeedsNotSkipped[NUM_DIRECTIONS];
    unsigned smallestSkippedSeed[NUM_DIRECTIONS];
    unsigned biggestClusterScored;
    unsigned highestWeightListChecked;
    bool usedHammingThisAlignment;

    double totalProbabilityByDepth[AlignerStats::maxMaxHits];
    void updateProbabilityMass();

        bool
    score(
        bool             forceResult,
        Read            *read[NUM_DIRECTIONS],
        AlignmentResult *result,
        int             *finalScore,
        unsigned        *singleHitGenomeLocation,
        Direction       *hitDirection,
        Candidate       *candidates,
        unsigned         maxHitsToGet,
        int             *mapq);
    
    void clearCandidates();

    bool findElement(unsigned genomeLocation, Direction direction, HashTableElement **hashTableElement);
    void findCandidate(unsigned genomeLocation, Direction direction, Candidate **candidate, HashTableElement **hashTableElement);
    void allocateNewCandidate(unsigned genomeLoation, Direction direction, unsigned lowestPossibleScore, int seedOffset, Candidate **candidate, HashTableElement **hashTableElement);
    void incrementWeight(HashTableElement *element);

    void fillHitsFound(unsigned maxHitsToGet, int *multiHitsFound, 
                       unsigned *multiHitLocations, Direction *multiHitDirections, int *multiHitScores);

    const Genome *genome;
    GenomeIndex *genomeIndex;
    SimilarityMap *similarityMap;
    unsigned seedLen;
    unsigned confDiff;
    unsigned maxHitsToConsider;
    unsigned maxK;
    unsigned maxReadSize;
    unsigned maxSeedsToUse; // Max number of seeds to look up in the hash table
    unsigned adaptiveConfDiffThreshold; // Increase confDiff by 1 if this many seeds are repetitive.

    int nCandidates;
    Candidate *candidates;

    char *rcReadData;
    char *rcReadQuality;       // This is allocated along with rcReadData in a single BigAlloc call, so don't free it.
    char *reversedRead[NUM_DIRECTIONS];
    char *reversedGenomeData;  // This will also be allocated along with reversedRead

    unsigned nTable[256];

    int readId;
    
    // Store the best hits at a given edit distance, as well as their number
    static const int MAX_MULTI_HITS_TO_GET = 512;
    unsigned  hitCount[MAX_K];
    unsigned  hitLocations[MAX_K][MAX_MULTI_HITS_TO_GET];
    Direction hitDirections[MAX_K][MAX_MULTI_HITS_TO_GET];

    // How many overly popular (> maxHits) seeds we skipped this run
    unsigned popularSeedsSkipped;

    bool explorePopularSeeds; // Whether we should explore the first maxHits hits even for overly
                              // popular seeds (useful for filtering reads that come from a database
                              // with many very similar sequences).

    bool stopOnFirstHit;      // Whether to stop the first time a location matches with less than
                              // maxK edit distance (useful when using SNAP for filtering only).

    AlignerStats *stats;

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
