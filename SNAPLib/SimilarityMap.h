#pragma once

#include "options.h"
#include "BigAlloc.h"
#include "Compat.h"
#include "RegionCluster.h"
#include "Genome.h"

//
// Holds a precomputed index of similar regions within a genome. The regions are grouped into
// "clusters" of substrings of a given length, so that each location in the genome is part of
// a unique cluster, or part of no cluster if it is not similar to any other substring. Each
// cluster's ID is the location of the first substring inside it. The special cluster ID
// 0xFFFFFFFF means "no cluster".
//
class SimilarityMap {
public:
    static const unsigned NO_CLUSTER = 0xFFFFFFFF;

    struct ClusterInfo {
        unsigned clusterId;
        RegionCluster *cluster;

        void *operator new[](size_t size) {return BigAlloc(size);}
        void operator delete[](void *ptr) {BigDealloc(ptr);}
    };

    SimilarityMap(unsigned genomeLength_, unsigned numClusters_, int mergeDistance_);

    ~SimilarityMap();

    void addCluster(unsigned clusterId, RegionCluster *cluster);

    //
    // Return the cluster info for a particular genome location.
    //
    inline ClusterInfo getClusterInfo(unsigned location) { return clusterInfo[location]; }

    inline unsigned getNumClusterMembers(unsigned location) { return numClusterMembers[location]; }

    inline unsigned getNumClusters() { return numClusters; }

    inline unsigned getMergeDistance() { return mergeDistance; }

    //
    // Load a similarity map from a file created by SimFinder.
    //
    static SimilarityMap *load(const char *filename, const Genome *genome, bool computeMemberInfo = true);

private:
    unsigned genomeLength;
    unsigned numClusters;
    int mergeDistance;

    ClusterInfo *clusterInfo;         // indexed by genome location
    unsigned *numClusterMembers;      // indexed by genome location
};
