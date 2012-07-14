/*++

Module Name:

    SimplePairedEndAligner.h

Abstract:

    First cut at a paired-end aligner.

Authors:

    Matei Zaharia, December, 2011

Environment:

    User mode service.

Revision History:

--*/

#pragma once

#include "PairedEndAligner.h"
#include "BaseAligner.h"

class SimplePairedEndAligner : public PairedEndAligner
{
public:
    SimplePairedEndAligner(
        GenomeIndex  *_index,
        unsigned      _maxReadSize,
        unsigned      _confDiff,
        unsigned      _maxHits,
        unsigned      _maxK,
        unsigned      _maxSeeds,
        unsigned      _maxSpacing,                 // Maximum distance to allow between the two ends.
        unsigned      _adaptiveConfDiffThreshold,  // Increase confDiff if this many seeds in the read have multiple hits.
        unsigned      _multiHitConfDiff,           // Difference between best and second-best pairs if multiple pairs align.
        unsigned      _multiHitConfDiffThreshold); // Increase multiHitConfDiff if there are more than this many good pairs.
    
    virtual ~SimplePairedEndAligner();
    
    virtual void align(
        Read                  *read0,
        Read                  *read1,
        PairedAlignmentResult *result);

private:
    static const int MAX_HITS_PER_READ = 1024;  // If both reads are multi-hit, how many locations to get for each?
    
    GenomeIndex  *index;
    
    unsigned      maxReadSize;
    
    unsigned      confDiff;
    unsigned      maxHits;
    unsigned      maxK;
    unsigned      maxSeeds;
    unsigned      maxSpacing;
    unsigned      adaptiveConfDiffThreshold;
    unsigned      multiHitConfDiff;
    unsigned      multiHitConfDiffThreshold;
    
    BaseAligner  *singleAligner;

    // Multi-hit state shared across methods (kind of ugly)
    int numMultiHits[2];
    unsigned multiHitLocs[2][MAX_HITS_PER_READ];
    bool multiHitRCs[2][MAX_HITS_PER_READ];
    int multiHitScores[2][MAX_HITS_PER_READ];

    void handleMultiHits1(Read *read0, Read *read1, PairedAlignmentResult *result);
    
    void handleMultiHits2(Read *read0, Read *read1, PairedAlignmentResult *result);
};
