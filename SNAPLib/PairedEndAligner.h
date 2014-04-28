/*++

Module Name:

    PairedEndAligner.h

Abstract:

    Superclass for implementations of paired end aligners

Authors:

    Matei Zaharia, December, 2011

Environment:

    User mode service.

Revision History:

--*/

#pragma once

#include "Aligner.h"
#include "directions.h"
#include "LandauVishkin.h"

const int NUM_READS_PER_PAIR = 2;    // This is just to make it clear what the array subscripts are, it doesn't ever make sense to change

struct PairedAlignmentResult {
    AlignmentResult status[NUM_READS_PER_PAIR]; // SingleHit or CertainHit if aligned, MultipleHit if matches DB
                                                // but not confidently aligned, or NotFound.

    unsigned location[NUM_READS_PER_PAIR];      // Genome location of each read.
    
    Direction direction[NUM_READS_PER_PAIR];    // Did we match the reverse complement? In general the two reads should have
                                                // opposite orientations because they're part of the same original fragment,
                                                // but it seems possible for a piece of the genome to get cut cleanly and flip
                                                // in a translocation event, which would cause both ends of a fragment aligning
                                                // there to be in the same orientation w.r.t. the reference genome.

    int score[NUM_READS_PER_PAIR];              // score of each end if matched

    int mapq[NUM_READS_PER_PAIR];               // mapping quality of each end, encoded like a Phred score (but as an integer, not ASCII Phred + 33).

    bool isTranscriptome[NUM_READS_PER_PAIR];   //is this alignment a transcriptome alignment or genome alignment
    unsigned tlocation[NUM_READS_PER_PAIR];     //pos of the trancriptome alignment
    //char flag[NUM_READS_PER_PAIR];              //alignment flag

    bool fromAlignTogether;                     // Was this alignment created by aligning both reads together, rather than from some combination of single-end aligners?
    bool alignedAsPair;                         // Were the reads aligned as a pair, or separately?
    _int64 nanosInAlignTogether;
    unsigned nLVCalls;
    unsigned nSmallHits;
};

/**
 * Abstract interface for paired-end aligners.
 */
class PairedEndAligner
{
public:
    virtual ~PairedEndAligner() {}
    
    virtual void align(
        Read                  *read0,
        Read                  *read1,
        PairedAlignmentResult *result,
        int                    maxEditDistanceForSecondaryResults,
        int                    secondaryResultBufferSize,
        int                   *nSecondaryResults,
        PairedAlignmentResult *secondaryResults,             // The caller passes in a buffer of secondaryResultBufferSize and it's filled in by align()
        int                    singleSecondaryBufferSize,
        int                   *nSingleEndSecondaryResultsForFirstRead,
        int                   *nSingleEndSecondaryResultsForSecondRead,
        SingleAlignmentResult *singleEndSecondaryResults     // Single-end secondary alignments for when the paired-end alignment didn't work properly
        ) = 0;

    virtual void setLandauVishkin(
        LandauVishkin<1>        *landauVishkin,
        LandauVishkin<-1>       *reverseLandauVishkin)
    {
    }

    virtual _int64 getLocationsScored() const  = 0;
};
