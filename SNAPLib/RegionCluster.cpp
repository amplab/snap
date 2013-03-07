#include "stdafx.h"
#include "RegionCluster.h"
#include "LandauVishkin.h"

using std::sort;
using std::string;
using std::vector;


void RegionCluster::initialize(unsigned id_, unsigned numMembers_, unsigned stringLength_)
{
    id = id_;
    numMembers = numMembers_;
    stringLength = stringLength_;
    members.reserve(numMembers);
    consensusString = NULL;
    firstMemberIndexByDistance = NULL;
}


RegionCluster::~RegionCluster()
{
    if (consensusString != NULL) {
        delete[] consensusString;
    }
    if (firstMemberIndexByDistance != NULL) {
        delete[] firstMemberIndexByDistance;
    }
}


void RegionCluster::computeMemberInfo(const Genome *genome)
{
  static LandauVishkin<> lv; // TODO: Not thread-safe!

    // Compute the consensus string
    consensusString = new char[stringLength];

    FixedSizeVector<const char*> memberStrings(numMembers);
    for (int i = 0; i < numMembers; i++) {
        const char *str = genome->getSubstring(members[i].location, stringLength);
        if (str == NULL) {
            fprintf(stderr, "Genome.getSubstring returned NULL on cluster member %u\n",
                    members[i].location);
            soft_exit(1);
        }
        memberStrings.push_back(str);
    }

    int counts[128];  // Count the number of occurrences of each string in each position

    for (int pos = 0; pos < stringLength; pos++) {
        counts['A'] = counts['C'] = counts['G'] = counts['T'] = 0;
        for (int i = 0; i < numMembers; i++) {
            counts[memberStrings[i][pos]]++;
        }
        // Find the most popular character at position i in the members
        char best = 'A';
        int bestCount = counts['A'];
        if (counts['C'] > bestCount) {
            best = 'C'; bestCount = counts['C'];
        }
        if (counts['G'] > bestCount) {
            best = 'G'; bestCount = counts['G'];
        }
        if (counts['T'] > bestCount) {
            best = 'T'; bestCount = counts['T'];
        }
        consensusString[pos] = best;
    }

    // Find the maximum edit distance from the consensus string to any member
    maxEditDistanceToConsensus = 0;
    for (int i = 0; i < numMembers; i++) {
        int distance = lv.computeEditDistance(
            consensusString, stringLength, memberStrings[i], stringLength, stringLength);
        members[i].editDistanceToConsensus = distance;
        if (distance > maxEditDistanceToConsensus) {
            maxEditDistanceToConsensus = distance;
        }
    }

    // Sort the members by edit distance to consensus
    sort(members.begin(), members.end(), MemberComparator());

    // Fix up the member strings array because we've re-sorted our members
    for (int i = 0; i < numMembers; i++) {
        memberStrings[i] = genome->getSubstring(members[i].location, stringLength);
    }

    // Figure out first index at each edit distance
    firstMemberIndexByDistance = new int[maxEditDistanceToConsensus + 2];
    memset(firstMemberIndexByDistance, 0, (maxEditDistanceToConsensus + 2) * sizeof(int));
    int prevDistance = 0;
    for (int i = 0; i < numMembers; i++) {
        if (members[i].editDistanceToConsensus > prevDistance) {
            for (int d = prevDistance + 1; d <= members[i].editDistanceToConsensus; d++) {
                firstMemberIndexByDistance[d] = i;
            }
            prevDistance = members[i].editDistanceToConsensus;
        }
    }
    firstMemberIndexByDistance[maxEditDistanceToConsensus + 1] = members.size();   // Sentinel

    // Compute the consensus string's mini-seeds
    int miniSeedsPerMember = stringLength / MINI_SEED_LENGTH;
    consensusMiniSeeds.reserve(miniSeedsPerMember);
    for (int i = 0; i < miniSeedsPerMember; i++) {
        consensusMiniSeeds.push_back(hashMiniSeed(consensusString + (MINI_SEED_LENGTH * i)));
    }
}
