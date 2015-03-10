/*++

Module Name:

    IntersectingPairedEndAligner.h

Abstract:

    A paired-end aligner based on set intersections to narrow down possible candidate locations.

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
#include "FixedSizeMap.h"

const unsigned DEFAULT_INTERSECTING_ALIGNER_MAX_HITS = 2000;
const unsigned DEFAULT_MAX_CANDIDATE_POOL_SIZE = 1000000;

class IntersectingPairedEndAligner : public PairedEndAligner
{
public:
    IntersectingPairedEndAligner(
        GenomeIndex  *index_,
        unsigned      maxReadSize_,
        unsigned      maxHits_,
        unsigned      maxK_,
        unsigned      maxSeedsFromCommandLine_,
        double        seedCoverage_,
        unsigned      minSpacing_,                 // Minimum distance to allow between the two ends.
        unsigned      maxSpacing_,                 // Maximum distance to allow between the two ends.
        unsigned      maxBigHits_,
        unsigned      extraSearchDepth_,
        unsigned      maxCandidatePoolSize,
        int           maxSecondaryAlignmentsPerContig_,
        BigAllocator  *allocator,
        bool          noUkkonen_,
        bool          noOrderedEvaluation_,
		bool		  noTruncation_);

     static unsigned getMaxSecondaryResults(unsigned numSeedsFromCommandLine, double seedCoverage, unsigned maxReadSize, unsigned maxHits, unsigned seedLength, unsigned minSpacing, unsigned maxSpacing)
     {
        unsigned maxSeedsToUse;
        if (0 != numSeedsFromCommandLine) {
            maxSeedsToUse = numSeedsFromCommandLine;
        } else {
            maxSeedsToUse = (unsigned)(maxReadSize * seedCoverage / seedLength);
        }

		//
		// The number of hits we can conceivably get is for each seed a result for every hit, times every possible pair for that hit.  The possible pairs
		// run from min to max distance on either side, but if they're within max merge distance then they'll be merged.
		//
        return NUM_DIRECTIONS * maxHits * maxSeedsToUse * (maxSpacing - minSpacing + 1 + maxMergeDistance - 1) / maxMergeDistance * 2;
     }

     void setLandauVishkin(
        LandauVishkin<1> *landauVishkin_,
        LandauVishkin<-1> *reverseLandauVishkin_) 
    {
        landauVishkin = landauVishkin_;
        reverseLandauVishkin = reverseLandauVishkin_;
    }
    
    virtual ~IntersectingPairedEndAligner();
    
    virtual void align(
        Read                  *read0,
        Read                  *read1,
        PairedAlignmentResult *result,
        int                    maxEditDistanceForSecondaryResults,
        int                    secondaryResultBufferSize,
        int                   *nSecondaryResults,
        PairedAlignmentResult *secondaryResults,             // The caller passes in a buffer of secondaryResultBufferSize and it's filled in by align()
        int                    singleSecondaryBufferSize,
        int                    maxSecondaryResultsToReturn,
        int                   *nSingleEndSecondaryResultsForFirstRead,
        int                   *nSingleEndSecondaryResultsForSecondRead,
        SingleAlignmentResult *singleEndSecondaryResults     // Single-end secondary alignments for when the paired-end alignment didn't work properly
        );

    static size_t getBigAllocatorReservation(GenomeIndex * index, unsigned maxBigHitsToConsider, unsigned maxReadSize, unsigned seedLen, unsigned maxSeedsFromCommandLine, 
                                             double seedCoverage, unsigned maxEditDistanceToConsider, unsigned maxExtraSearchDepth, unsigned maxCandidatePoolSize,
                                             int maxSecondaryAlignmentsPerContig);

    void *operator new(size_t size, BigAllocator *allocator) {_ASSERT(size == sizeof(IntersectingPairedEndAligner)); return allocator->allocate(size);}
    void operator delete(void *ptr, BigAllocator *allocator) {/* do nothing.  Memory gets cleaned up when the allocator is deleted.*/}

    void *operator new(size_t size) {return BigAlloc(size);}
    void operator delete(void *ptr) {BigDealloc(ptr);}

    virtual _int64 getLocationsScored() const {
         return nLocationsScored;
     }


private:

    IntersectingPairedEndAligner() {}  // This is for the counting allocator, it doesn't build a useful object

    static const int NUM_SET_PAIRS = 2;         // A "set pair" is read0 FORWARD + read1 RC, or read0 RC + read1 FORWARD.  Again, it doesn't make sense to change this.

    void allocateDynamicMemory(BigAllocator *allocator, unsigned maxReadSize, unsigned maxBigHitsToConsider, unsigned maxSeedsToUse, 
                               unsigned maxEditDistanceToConsider, unsigned maxExtraSearchDepth, unsigned maxCandidatePoolSize,
                               int maxSecondaryAlignmentsPerContig);

    GenomeIndex *   index;
    const Genome *  genome;
    GenomeDistance  genomeSize;
    unsigned        maxReadSize;
    unsigned        maxHits;
    unsigned        maxBigHits;
    unsigned        extraSearchDepth;
    unsigned        maxK;
    unsigned        numSeedsFromCommandLine;
    double          seedCoverage;
    static const unsigned MAX_MAX_SEEDS = 30;
    unsigned        minSpacing;
    unsigned        maxSpacing;
    unsigned        seedLen;
    bool            doesGenomeIndexHave64BitLocations;
    _int64          nLocationsScored;
    bool            noUkkonen;
    bool            noOrderedEvaluation;
	bool			noTruncation;

	static const unsigned        maxMergeDistance;
	
	//
	// It's a template, because we 
    // have different sizes of genome locations depending on the hash table format.  So, GL must be unsigned or GenomeLocation
    //    
    template<class GL> struct HashTableLookup {
        unsigned        seedOffset;
        _int64          nHits;
        const GL  *     hits;
        unsigned        whichDisjointHitSet;

        //
        // We keep the hash table lookups that haven't been exhaused in a circular list.
        //
        HashTableLookup<GL> *nextLookupWithRemainingMembers;
        HashTableLookup<GL> *prevLookupWithRemainingMembers;

        //
        // State for handling the binary search of a location in this lookup.
        // This would ordinarily be stack local state in the binary search
        // routine, but because a) we want to interleave the steps of the binary
        // search in order to allow cache prefetches to have time to execute;
        // and b) we don't want to do dynamic memory allocation (really at all),
        // it gets stuck here.
        //
        int limit[2];   // The upper and lower limits of the current binary search in hits
        GL maxGenomeLocationToFindThisSeed;
        
        //
        // A linked list of lookups that haven't yet completed this binary search.  This is a linked
        // list with no header element, so testing for emptiness needs to happen at removal time.
        // It's done that way to avoid a comparison for list head that would result in a hard-to-predict
        // branch.
        //
        HashTableLookup<GL> *nextLookupForCurrentBinarySearch;
        HashTableLookup<GL> *prevLookupForCurrentBinarySearch;

        _int64           currentHitForIntersection;

        //
        // A place for the hash table to write in singletons.  We need this because when the hash table is
        // built with > 4 byte genome locations, it usually doesn't store 8 bytes, so we need to
        // provide the lookup function a place to write the result.  Since we need one per
        // lookup, it goes here.
        //
        GL singletonGenomeLocation[2];  // The [2] is because we need to look one before sometimes, and that allows space
    };
    
    //
    // A set of seed hits, represented by the lookups that came out of the big hash table.  It can be over 32 or
    // 64 bit indices, but its external interface is always 64 bits (it extends on the way out if necessary).
    //
    class HashTableHitSet {
    public:
        HashTableHitSet() {}
        void firstInit(unsigned maxSeeds_, unsigned maxMergeDistance_, BigAllocator *allocator, bool doesGenomeIndexHave64BitLocations_);

        //
        // Reset to empty state.
        //
        void init();

        //
        // Record a hash table lookup.  All recording must be done before any
        // calls to getNextHitLessThanOrEqualTo.  A disjoint hit set is a set of hits
		// that don't share any bases in the read.  This is interesting because the edit
		// distance of a read must be at least the number of seeds that didn't hit for
		// any disjoint hit set (because there must be a difference in the read within a
		// seed for it not to hit, and since the reads are disjoint there can't be a case
		// where the same difference caused two seeds to miss).
        //
        void recordLookup(unsigned seedOffset, _int64 nHits, const unsigned *hits, bool beginsDisjointHitSet);
        void recordLookup(unsigned seedOffset, _int64 nHits, const GenomeLocation *hits, bool beginsDisjointHitSet);

        //
        // This efficiently works through the set looking for the next hit at or below this address.
        // A HashTableHitSet only allows a single iteration through its address space per call to
        // init().
        //
        bool    getNextHitLessThanOrEqualTo(GenomeLocation maxGenomeLocationToFind, GenomeLocation *actualGenomeLocationFound, unsigned *seedOffsetFound);

        //
        // Walk down just one step, don't binary search.
        //
        bool getNextLowerHit(GenomeLocation *genomeLocation, unsigned *seedOffsetFound);


        //
        // Find the highest genome address.
        //
        bool    getFirstHit(GenomeLocation *genomeLocation, unsigned *seedOffsetFound);

		unsigned computeBestPossibleScoreForCurrentHit();

        //
        // This is bit of storage that the 64 bit lookup needs in order to extend singleton hits into 64 bits, since they may be
        // stored in the index in fewer.
        //
        GenomeLocation *getNextSingletonLocation()
        {
            return &lookups64[nLookupsUsed].singletonGenomeLocation[1];
        }


    private:
        struct DisjointHitSet {
            unsigned countOfExhaustedHits;
            unsigned missCount;
        };

        int                                 currentDisjointHitSet;
        DisjointHitSet  *                   disjointHitSets;
        HashTableLookup<unsigned> *         lookups32;
        HashTableLookup<GenomeLocation> *   lookups64;
        HashTableLookup<unsigned>           lookupListHead32[1];
        HashTableLookup<GenomeLocation>     lookupListHead64[1];
        unsigned                            maxSeeds;
        unsigned                            nLookupsUsed;
        GenomeLocation                      mostRecentLocationReturned;
		unsigned		                    maxMergeDistance;
        bool                                doesGenomeIndexHave64BitLocations;
    };

    HashTableHitSet *                       hashTableHitSets[NUM_READS_PER_PAIR][NUM_DIRECTIONS];

    int                                     countOfHashTableLookups[NUM_READS_PER_PAIR];
    _int64                                  totalHashTableHits[NUM_READS_PER_PAIR][NUM_DIRECTIONS];
    _int64                                  largestHashTableHit[NUM_READS_PER_PAIR][NUM_DIRECTIONS];
    unsigned                                readWithMoreHits;
    unsigned                                readWithFewerHits;

    //
    // A location that's been scored (or waiting to be scored).  This is needed in order to do merging
    // of close-together hits and to track potential mate pairs.
    //
    struct HitLocation {
        GenomeLocation  genomeLocation;
        int             genomeLocationOffset;   // This is needed because we might get an offset back from scoring (because it's really scoring a range).
        unsigned        seedOffset;
        bool            isScored;           // Mate pairs are sometimes not scored when they're inserted, because they
        unsigned        score;
        unsigned        maxK;               // The maxK that this was scored with (we may need to rescore if we need a higher maxK and score is -1)
        double          matchProbability;
		unsigned		bestPossibleScore;

        //
        // We have to be careful in the case where lots of offsets in a row match well against the read (think
        // about repetitive short sequences, i.e., ATTATTATTATT...).  We want to merge the close ones together,
        // but if the repetitive sequence extends longer than maxMerge, we don't want to just slide the window
        // over the whole range and declare it all to be one.  There is really no good definition for the right
        // thing to do here, so instead all we do is that when we declare two candidates to be matched we
        // pick one of them to be the match primary and then coalesce all matches that are within maxMatchDistance
        // of the match primary.  No one can match with any of the locations in the set that's beyond maxMatchDistance
        // from the set primary.  This means that in the case of repetitve sequences that we'll declare locations
        // right next to one another not to be matches.  There's really no way around this while avoiding
        // matching things that are possibly much more than maxMatchDistance apart.
        //
        GenomeLocation  genomeLocationOfNearestMatchedCandidate;
    };


    char *rcReadData[NUM_READS_PER_PAIR];                   // the reverse complement of the data for each read
    char *rcReadQuality[NUM_READS_PER_PAIR];                // the reversed quality strings for each read
    unsigned readLen[NUM_READS_PER_PAIR];

    Read *reads[NUM_READS_PER_PAIR][NUM_DIRECTIONS];        // These are the reads that are provided in the align call, together with their reverse complements, which are computed.
    Read rcReads[NUM_READS_PER_PAIR][NUM_DIRECTIONS];

    char *reversedRead[NUM_READS_PER_PAIR][NUM_DIRECTIONS]; // The reversed data for each read for forward and RC.  This is used in the backwards LV

    LandauVishkin<> *landauVishkin;
    LandauVishkin<-1> *reverseLandauVishkin;

    char rcTranslationTable[256];
    unsigned nTable[256];

    BYTE *seedUsed;

    inline bool IsSeedUsed(_int64 indexInRead) const {
        return (seedUsed[indexInRead / 8] & (1 << (indexInRead % 8))) != 0;
    }

    inline void SetSeedUsed(_int64 indexInRead) {
        seedUsed[indexInRead / 8] |= (1 << (indexInRead % 8));
    }

    //
    // "Local probability" means the probability that each end is correct given that the pair itself is correct.
    // Consider the example where there's exactly one decent match for one read, but the other one has several
    // that are all within the correct range for the first one.  Then the local probability for the second read
    // is lower than the first.  The overall probability of an alignment then is 
    // pairProbability * localProbability/ allPairProbability.
    //
    double localBestPairProbability[NUM_READS_PER_PAIR];

    void scoreLocation(
            unsigned             whichRead,
            Direction            direction,
            GenomeLocation       genomeLocation,
            unsigned             seedOffset,
            unsigned             scoreLimit,
            unsigned            *score,
            double              *matchProbability,
            int                 *genomeLocationOffset   // The computed offset for genomeLocation (which is needed because we scan several different possible starting locations)
    );

    //
    // These are used to keep track of places where we should merge together candidate locations for MAPQ purposes, because they're sufficiently
    // close in the genome.
    //
    struct MergeAnchor {
        double          matchProbability;
        GenomeLocation  locationForReadWithMoreHits;
        GenomeLocation  locationForReadWithFewerHits;
        int             pairScore;

        void init(GenomeLocation locationForReadWithMoreHits_, GenomeLocation locationForReadWithFewerHits_, double matchProbability_, int pairScore_) {
            locationForReadWithMoreHits = locationForReadWithMoreHits_;
            locationForReadWithFewerHits = locationForReadWithFewerHits_;
            matchProbability = matchProbability_;
            pairScore = pairScore_;
        }

        //
        // Returns whether this candidate is a match for this merge anchor.
        //
        bool doesRangeMatch(GenomeLocation newMoreHitLocation, GenomeLocation newFewerHitLocation) {
            GenomeDistance deltaMore = DistanceBetweenGenomeLocations(locationForReadWithMoreHits, newMoreHitLocation);
            GenomeDistance deltaFewer = DistanceBetweenGenomeLocations(locationForReadWithFewerHits, newFewerHitLocation);

            return deltaMore < 50 && deltaFewer < 50;
        }


        //
        // Returns true and sets oldMatchProbability if this should be eliminated due to a match.
        //
        bool checkMerge(GenomeLocation newMoreHitLocation, GenomeLocation newFewerHitLocation, double newMatchProbability, int newPairScore, 
                        double *oldMatchProbability); 
    };

    //
    // We keep track of pairs of locations to score using two structs, one for each end.  The ends for the read with fewer hits points into
    // a list of structs for the end with more hits, so that we don't need one stuct for each pair, just one for each end, and also so that 
    // we don't need to score the mates more than once if they could be paired with more than one location from the end with fewer hits.
    //

    struct ScoringMateCandidate {
        //
        // These are kept in arrays in decreasing genome order, one for each set pair, so you can find the next largest location by just looking one
        // index lower, and vice versa.
        //
        double                  matchProbability;
        GenomeLocation          readWithMoreHitsGenomeLocation;
        unsigned                bestPossibleScore;
        unsigned                score;
        unsigned                scoreLimit;             // The scoreLimit with which score was computed
        unsigned                seedOffset;
        int                     genomeOffset;

        void init(GenomeLocation readWithMoreHitsGenomeLocation_, unsigned bestPossibleScore_, unsigned seedOffset_) {
            readWithMoreHitsGenomeLocation = readWithMoreHitsGenomeLocation_;
            bestPossibleScore = bestPossibleScore_;
            seedOffset = seedOffset_;
            score = -2;
            scoreLimit = -1;
            matchProbability = 0;
            genomeOffset = 0;
        }
    };

    struct ScoringCandidate {
        ScoringCandidate *      scoreListNext;              // This is a singly-linked list
        MergeAnchor *           mergeAnchor;
        unsigned                scoringMateCandidateIndex;  // Index into the array of scoring mate candidates where we should look 
        GenomeLocation          readWithFewerHitsGenomeLocation;
        unsigned                whichSetPair;
        unsigned                seedOffset;

        unsigned                bestPossibleScore;

        void init(GenomeLocation readWithFewerHitsGenomeLocation_, unsigned whichSetPair_, unsigned scoringMateCandidateIndex_, unsigned seedOffset_,
                  unsigned bestPossibleScore_, ScoringCandidate *scoreListNext_)
        {
            readWithFewerHitsGenomeLocation = readWithFewerHitsGenomeLocation_;
            whichSetPair = whichSetPair_;
            _ASSERT(whichSetPair < NUM_SET_PAIRS);  // You wouldn't think this would be necessary, but...
            scoringMateCandidateIndex = scoringMateCandidateIndex_;
            seedOffset = seedOffset_;
            bestPossibleScore = bestPossibleScore_;
            scoreListNext = scoreListNext_;
            mergeAnchor = NULL;
         }
    };

    //
    // A pool of scoring candidates.  For each alignment call, we free them all by resetting lowestFreeScoringCandidatePoolEntry to 0,
    // and then fill in the content when they're initialized.  This means that for alignments with few candidates we'll be using the same
    // entries over and over, so they're likely to be in the cache.  We have maxK * maxSeeds * 2 of these in the pool, so we can't possibly run
    // out.  We rely on their being allocated in descending genome order within a set pair.
    //
    ScoringCandidate *scoringCandidatePool;
    unsigned scoringCandidatePoolSize;
    unsigned lowestFreeScoringCandidatePoolEntry;

    //
    // maxK + 1 lists of Scoring Candidates.  The lists correspond to bestPossibleScore for the candidate and its best mate.
    //

    ScoringCandidate    **scoringCandidates;

    //
    // The scoring mates.  The each set scoringCandidatePoolSize / 2.
    //
    ScoringMateCandidate * scoringMateCandidates[NUM_SET_PAIRS];
    unsigned lowestFreeScoringMateCandidate[NUM_SET_PAIRS];

    //
    // Merge anchors.  Again, we allocate an upper bound number of them, which is the same as the number of scoring candidates.
    //
    MergeAnchor *mergeAnchorPool;
    unsigned firstFreeMergeAnchor;
    unsigned mergeAnchorPoolSize;


    struct HitsPerContigCounts {
        _int64  epoch;              // Rather than zeroing this whole array every time, we just bump the epoch number; results with an old epoch are considered zero
        int     hits;
    };

    HitsPerContigCounts *hitsPerContigCounts;   // How many alignments are we reporting for each contig.  Used to implement -mpc, otheriwse unallocated.
    int maxSecondaryAlignmentsPerContig;
    _int64 contigCountEpoch;
};
