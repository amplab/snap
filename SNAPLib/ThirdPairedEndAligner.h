/*++

Module Name:

    ThirdPairedEndAligner.h

Abstract:

    A paired-end aligner that works much like the base aligner.

Authors:

    Bill Bolosky, February, 2013

Environment:

    User mode service.

Revision History:

--*/

#pragma once

#include "PairedEndAligner.h"
#include "BaseAligner.h"
#include "BigAlloc.h"
#include "directions.h"
#include "LandauVishkin.h"

#if     0


class ThirdPairedEndAligner : public PairedEndAligner
{
public:
    ThirdPairedEndAligner(
        GenomeIndex  *index_,
        unsigned      maxReadSize_,
        unsigned      maxHits_,
        unsigned      maxK_,
        unsigned      maxSeeds_,
        unsigned      minSpacing_,                 // Minimum distance to allow between the two ends.
        unsigned      maxSpacing_,                 // Maximum distance to allow between the two ends.
        BigAllocator  *allocator); 
    
    virtual ~ThirdPairedEndAligner();
    
    virtual void align(
        Read                  *read0,
        Read                  *read1,
        PairedAlignmentResult *result);

    void *operator new(size_t size, BigAllocator *allocator) {_ASSERT(size == sizeof(ThirdPairedEndAligner)); return allocator->allocate(size);}
    void operator delete(void *ptr, BigAllocator *allocator) {/*Do nothing, the owner of the allocator is responsible for freeing memory*/}

    static size_t getBigAllocatorReservation(GenomeIndex * index, unsigned maxHitsToConsider, unsigned maxReadSize, unsigned seedLen, unsigned maxSeedsToUse);

private:

    ThirdPairedEndAligner() {}  // This is for the counting allocator, it doesn't build a useful object

    void clearCandidates();

    static const int NUM_READS_PER_PAIR = 2;    // This is just to make it clear what the array subscripts are, it doesn't ever make sense to change

    void allocateDynamicMemory(BigAllocator *allocator, unsigned maxReadSize, unsigned maxHitsToConsider, unsigned maxSeedsToUse);

    GenomeIndex *   index;
    const Genome *  genome;
    unsigned        genomeSize;
    unsigned        maxReadSize;
    unsigned        maxHits;
    unsigned        maxK;
    unsigned        maxSeeds;
    unsigned        minSpacing;
    unsigned        maxSpacing;
    unsigned        seedLen;
    unsigned        distanceToSearchBeyondBestScore;

    char *rcReadData[NUM_READS_PER_PAIR];                   // the reverse complement of the data for each read
    char *rcReadQuality[NUM_READS_PER_PAIR];                // the reversed quality strings for each read
    unsigned readLen[NUM_READS_PER_PAIR];

    Read *reads[NUM_READS_PER_PAIR][NUM_DIRECTIONS];        // These are the reads that are provided in the align call, together with their reverse complements, which are computed.

    char *reversedRead[NUM_READS_PER_PAIR][NUM_DIRECTIONS]; // The reversed data for each read for forward and RC.  This is used in the backwards LV

    LandauVishkin<> *landauVishkin;
    LandauVishkin<-1> *reverseLandauVishkin;

    char rcTranslationTable[256];
    unsigned nTable[256];

    BaseAligner *baseAligner;

    BYTE *seedUsed;

    inline bool IsSeedUsed(unsigned indexInRead) const {
        return (seedUsed[indexInRead / 8] & (1 << (indexInRead % 8))) != 0;
    }

    inline void SetSeedUsed(unsigned indexInRead) {
        seedUsed[indexInRead / 8] |= (1 << (indexInRead % 8));
    }

    void alignWithBaseAligner(Read *read0, Read *read1, PairedAlignmentResult *result, int maxMapq);

    //
    // We keep track of the hash table hits in two ways.  First, each hit is recorded in a Candidate structure, which keeps track of the hit for each read, and is used
    // to coalsece nearby hits (so that we don't treat close alignments with indels as different and mess up our mapq).  Candidates track 64 offsets each.
    // The second way is with SuperGroups, which are used to find potential mate pairs.  A SuperGroup is just a bit vector, each bit of which indicates whether a given
    // Candidate (structure) exists.
    // 
    // Both candidates and SuperGroups are stored in hash tables.  There's one hash table for each direction of each read.
    //
    // Like in the BaseAligner, the hash tables aren't cleared at the beginning of each alignment, because it would be too slow.  Instead, each hash table entry contains
    // an epoch number, and any entries that have old epochs are virtually empty.  That way, clearing the hash tables just means incrementing the epoch number.
    //
    unsigned hashTableEpochNumber;

    unsigned candidateGroupHashTableSize;
    unsigned superGroupHashTableSize;

    static const unsigned GroupSpan = 64;                     // Can't be > 64, since we use a _uint64 for a bit vector
    static const unsigned SuperGroupSpan = 64 * GroupSpan;    // Again, can't be > 64 Groups, since we use _uint64 bit vector

    struct SuperGroup {
        inline void init() {usedGroups = 0;}

        unsigned baseOffset;
        _uint64 usedGroups;
        SuperGroup *next;
    };

    struct SuperGroupHashTableAnchor {
        unsigned             hashTableEpoch;
        SuperGroup          *supergroups;
    };

    SuperGroupHashTableAnchor *superGroupHashTable[NUM_READS_PER_PAIR][NUM_DIRECTIONS];

    bool isThereAMateCandidate(unsigned genomeLocation, unsigned whichRead, Direction dir);
    bool isThereACandidateInRange(unsigned minLoc, unsigned maxLox, unsigned whichRead, Direction dir);

    struct Candidate {
        unsigned seedOffset;
        unsigned score;
    };

    struct CandidateGroup {
        void init(unsigned baseGenomeLocation_, unsigned whichRead_, Direction direction_) 
        {
            baseGenomeLocation = baseGenomeLocation_;
            weight = 0;
            bestScore= -1;
            whichRead = whichRead_;
            direction = direction_;
            allExtantCandidatesScored = false;   // Since there soon will be one
            matchProbabilityForBestScore = 0;
            //
            // Don't bother with the hash table pointer, it just gets immediately overwritten anyway.
            //
            weightListNext = weightListPrev = NULL;
            usedCandidates = 0;
            scoredCandidates = 0;
            hasAKnownMate = false;
        }

        void removeFromWeightList()
        {
            if (NULL == weightListNext) {
                _ASSERT(weightListPrev == NULL);
                return;
            }
            weightListNext->weightListPrev = weightListPrev;
            weightListPrev->weightListNext = weightListNext;

            weightListNext = weightListPrev = NULL;
        }

        unsigned    baseGenomeLocation;
        unsigned    weight;
        unsigned    bestPossibleScore;
        unsigned    bestScore;
        unsigned    whichRead;
        Direction   direction;
        bool        allExtantCandidatesScored;
        bool        hasAKnownMate;       // Set if we know that there's seed hit on the other read in the other direction that's in the right range
        double      matchProbabilityForBestScore;
        unsigned    bestScoreGenomeLocation;
        //
        // Singly linked list for hash table (don't need doubly linked, because we never remove)
        //
        CandidateGroup *next;

        //
        // Doubly linked list for weight bucket.
        //
        CandidateGroup *weightListNext;
        CandidateGroup *weightListPrev;

        //
        // A singly linked list of candidates used during scoring.
        //
        CandidateGroup *scoringNext;

        _uint64 usedCandidates;         // bit vector of used candidates
        _uint64 scoredCandidates;       // bit vector of scored candidate

        Candidate candidates[GroupSpan];
    };

    struct CandidateGroupHashTableAnchor {
        unsigned        hashTableEpoch;
        CandidateGroup  *candidateGroups;
    };

    CandidateGroupHashTableAnchor *candidateGroupHashTable[NUM_READS_PER_PAIR][NUM_DIRECTIONS];
    CandidateGroup *candidateGroupPool;
    SuperGroup *superGroupPool;

    CandidateGroup *weightLists;
    unsigned maxWeightListUsed;
    
    unsigned candidateGroupPoolSize;
    unsigned nUsedCandidateGroups;
    unsigned superGroupPoolSize;
    unsigned nUsedSuperGroups;

    //
    // The probability of the best pair we've found so far.
    //
    double probabilityOfBestPair;

    //
    // "Local probability" means the probability that each end is correct given that the pair itself is correct.
    // Consider the example where there's exactly one decent match for one read, but the other one has several
    // that are all within the correct range for the first one.  Then the local probability for the second read
    // is lower than the first.  The overall probability of an alignment then is 
    // pairProbability * localProbability/ allPairProbability.
    //
    double localBestPairProbability[NUM_READS_PER_PAIR];
    double probabilityOfAllPairs;
    unsigned bestPairScore;
    unsigned scoreLimit;
    unsigned bestResultGenomeLocation[NUM_READS_PER_PAIR];
    Direction bestResultDirection[NUM_READS_PER_PAIR];
    unsigned bestResultScore[NUM_READS_PER_PAIR];
    unsigned popularSeedsSkipped[NUM_READS_PER_PAIR];

    void findCandidateAndCreateIfNotExtant(unsigned whichRead, Direction dir, unsigned genomeLocation, CandidateGroup **group, Candidate **candidate, unsigned bestPossibleScoreIfNew);
    SuperGroup *findSuperGroupAndCreateIfNotExtant(unsigned whichRead, Direction dir, unsigned genomeLocation);

    static inline unsigned hash(unsigned value, unsigned tableSize) {return (value * 131) % tableSize;} // Very quick-n-dirty hash function.  Speed of hashing is more important than distribution here.

    bool score(
        bool                     forceResult,
        Read                    *read[NUM_READS_PER_PAIR][NUM_DIRECTIONS],
        PairedAlignmentResult   *result);

    unsigned bestPassSeedsNotSkipped[NUM_READS_PER_PAIR][NUM_DIRECTIONS];
    void addCandidateGroupsInRangeToScoringList(CandidateGroup **listHead, unsigned baseGenomeLocation, bool lookBelow, unsigned whichRead, Direction dir, unsigned *bestPossibleScoreOfAnyMate);
    _uint64 getSupergroupBitMaskForRange(unsigned supergroupBaseOffset, unsigned minLoc, unsigned maxLoc);
    CandidateGroup *lookupCandidateGroup(unsigned genomeLocation, unsigned whichRead, Direction dir);
    void scoreGroup(CandidateGroup *group, unsigned whichRead, Direction direction);
};

#endif  // 0
