#include "stdafx.h"
#include "SimilarityMap.h"
#include "LandauVishkin.h"


SimilarityMap::SimilarityMap(unsigned genomeLength_, unsigned numClusters_, int mergeDistance_)
    : genomeLength(genomeLength_), numClusters(numClusters_), mergeDistance(mergeDistance_)
{
    // MATEI: Skipping these for now because they take a lot of memory and we don't use them yet
    //clusterInfo = new ClusterInfo[genomeLength];
    //for (unsigned i = 0; i < genomeLength; i++) {
    //    clusterInfo[i].clusterId = NO_CLUSTER;
    //    clusterInfo[i].cluster = NULL;
    //}

    numClusterMembers = new unsigned[genomeLength];
    for (unsigned i = 0; i < genomeLength; i++) {
        numClusterMembers[i] = 1;
    }
}


SimilarityMap::~SimilarityMap()
{
    //for (unsigned i = 0; i < genomeLength; i++) {
    //    if (clusterInfo[i].clusterId == i) {
    //        delete clusterInfo[i].cluster;  // Only deletes each cluster once (on its first element)
    //    }
    //}
    //delete[] clusterInfo;
    
    delete[] numClusterMembers;
}


void SimilarityMap::addCluster(unsigned newClusterId, RegionCluster *cluster)
{
    for (int i = 0; i < cluster->members.size(); i++) {
        clusterInfo[cluster->members[i].location].clusterId = newClusterId;
        clusterInfo[cluster->members[i].location].cluster = cluster;
    }
}


SimilarityMap *SimilarityMap::load(const char *filename, const Genome *genome, bool computeMemberInfo)
{
    FILE *in = fopen(filename, "r");
    if (in == NULL) {
        fprintf(stderr, "Failed to open '%s' for reading\n", filename);
        return NULL;
    }

    unsigned genomeLen;
    unsigned numClusters;
    unsigned stringLen;
    int mergeDistance;
    if (fscanf(in, "%u %u %u %d", &genomeLen, &numClusters, &stringLen, &mergeDistance) != 4) {
        fprintf(stderr, "Could not read similarity map header from '%s'\n", filename);
        return NULL;
    }

    SimilarityMap *simMap = new SimilarityMap(genomeLen, numClusters, mergeDistance);
    //RegionCluster *clusters = new RegionCluster[numClusters];

    for (unsigned clusterNum = 0; clusterNum < numClusters; clusterNum++) {
        unsigned numMembers;
        if (fscanf(in, "%u", &numMembers) != 1) {
            fprintf(stderr, "Could not parse similarity map (EOF before seeing all clusters)\n");
            return NULL;
        }
        unsigned firstMember;
        if (fscanf(in, "%u", &firstMember) != 1) {
            fprintf(stderr, "Could not parse similarity map (EOF in cluster member list)\n");
            return NULL;
        }

        simMap->numClusterMembers[firstMember] = numMembers;
        for (unsigned i = 1; i < numMembers; i++) {
            unsigned member;
            if (fscanf(in, "%u", &member) != 1) {
                fprintf(stderr, "Could not parse similarity map (EOF in cluster member list)\n");
                return NULL;
            }
            simMap->numClusterMembers[member] = numMembers;
        }

        // MATEI: Skip initialization of RegionCluster objects for now
        //RegionCluster *cluster = &(clusters[clusterNum]);
        //cluster->initialize(firstMember, numMembers, stringLen);
        //cluster->members.push_back(RegionCluster::Member(firstMember));
        //for (unsigned i = 1; i < numMembers; i++) {
        //    unsigned member;
        //    if (fscanf(in, "%u", &member) != 1) {
        //        fprintf(stderr, "Could not parse similarity map (EOF in cluster member list)\n");
        //        return NULL;
        //    }
        //    cluster->members.push_back(RegionCluster::Member(member));
        //}
        //if (computeMemberInfo) {
        //    cluster->computeMemberInfo(genome);
        //}
        //simMap->addCluster(firstMember, cluster);
    }

    return simMap;
}
