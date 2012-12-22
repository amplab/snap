#pragma once

#include "options.h"
#include "Compat.h"
#include "Genome.h"
#include "FixedSizeMap.h"
#include "FixedSizeVector.h"
#include "Tables.h"

//
// A cluster of similar regions in a genome.
//
class RegionCluster {
public:
    struct Member {
        unsigned location;
        char editDistanceToConsensus;

        Member() {}

        Member(unsigned location_): location(location_) {}
    };

    RegionCluster(): consensusString(NULL) {}

    ~RegionCluster();

    void initialize(unsigned id_, unsigned numMembers_, unsigned stringLength_);

    // Call this method after filling in the members array to compute consensus string, mini-seeds,
    // and various statistics required for alignment.
    void computeMemberInfo(const Genome *genome);

    unsigned id;
    FixedSizeVector<Member> members;
    unsigned numMembers;
    unsigned stringLength;
    char *consensusString;
    FixedSizeVector<short> consensusMiniSeeds;

    // How far in edit distance is the farthest member from the consensus string? This can be used
    // to prune the whole cluster based on an edit distance calculation to the consensus string.
    unsigned maxEditDistanceToConsensus;

    int *firstMemberIndexByDistance;  // Used to search members list within a range of distances.

    static const int MINI_SEED_LENGTH = 5;  // Number of bases to use to compute mini-seeds
    static const int MINI_SEED_BITS = 2 * MINI_SEED_LENGTH; // Number of bits to store mini-seeds 
                                                            // (could be less if we hash them)

    // Hash data[0..MINI_SEED_LENGTH] as a 16-bit integer
    static inline short hashMiniSeed(const char *data) {
        short bases = 0;
        for (int i = 0; i < MINI_SEED_LENGTH; i++) {
            bases |= (BASE_VALUE_NO_N[data[i]] << (2 * i));
        }
        return bases;
    }

    // Sort members by edit distance to consensus
    struct MemberComparator {
        inline bool operator() (const Member& m1, const Member& m2) {
            if (m1.editDistanceToConsensus != m2.editDistanceToConsensus) {
                return m1.editDistanceToConsensus < m2.editDistanceToConsensus;
            } else {
                return m1.location < m2.location;
            }
        }
    };

    void *operator new[](size_t size) {return BigAlloc(size);}
    void operator delete[](void *ptr) {BigDealloc(ptr);}
};
