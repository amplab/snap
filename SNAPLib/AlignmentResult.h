/*++

Module Name:

    AlignmentResult.h

Abstract:

    Header for SNAP genome alignment results

Authors:

    Bill Bolosky, May, 2014

Environment:

    User mode service.

    This class is NOT thread safe.  It's the caller's responsibility to ensure that
    at most one thread uses an instance at any time.

Revision History:

   Pulled out of other places in SNAP

--*/

#pragma once
#include "Genome.h"
#include "directions.h"

class Read;

enum AlignmentResult {NotFound, SingleHit, MultipleHits, UnknownAlignment}; // BB: Changed Unknown to UnknownAlignment because of a conflict w/Windows headers

bool isAValidAlignmentResult(AlignmentResult result);


inline const char *AlignmentResultToString(AlignmentResult result) {
    switch (result) {
        case NotFound: return "NotFound";
        case SingleHit: return "SingleHit";
        case MultipleHits: return "MultipleHits";
        case UnknownAlignment: return "Unknown";
        default: return "Unknown alignment result type";
    }
}

struct SingleAlignmentResult {
	AlignmentResult status;

    GenomeLocation  location;	// Aligned genome location.
    Direction       direction;	// Did we match the reverse complement? 
    int             score;		// score of each end if matched

    int             mapq;		// mapping quality, encoded like a Phred score (but as an integer, not ASCII Phred + 33).

    static int compareByContigAndScore(const void *first, const void *second);      // qsort()-style compare routine
    static int compareByScore(const void *first, const void *second);               // qsort()-style compare routine
};

// Does an AlignmentResult represent a single location?
inline bool isOneLocation(AlignmentResult result) {
    return result == SingleHit;
}

const int NUM_READS_PER_PAIR = 2;    // This is just to make it clear what the array subscripts are, it doesn't ever make sense to change

struct PairedAlignmentResult {
	AlignmentResult status[NUM_READS_PER_PAIR]; // SingleHit or CertainHit if aligned, MultipleHit if matches DB
	// but not confidently aligned, or NotFound.

	GenomeLocation location[NUM_READS_PER_PAIR];// Genome location of each read.

	Direction direction[NUM_READS_PER_PAIR];    // Did we match the reverse complement? In general the two reads should have
	// opposite orientations because they're part of the same original fragment,
	// but it seems possible for a piece of the genome to get cut cleanly and flip
	// in a translocation event, which would cause both ends of a fragment aligning
	// there to be in the same orientation w.r.t. the reference genome.

	int score[NUM_READS_PER_PAIR];              // score of each end if matched

	int mapq[NUM_READS_PER_PAIR];               // mapping quality of each end, encoded like a Phred score (but as an integer, not ASCII Phred + 33).

	bool fromAlignTogether;                     // Was this alignment created by aligning both reads together, rather than from some combination of single-end aligners?
	bool alignedAsPair;                         // Were the reads aligned as a pair, or separately?
	_int64 nanosInAlignTogether;
	unsigned nLVCalls;
	unsigned nSmallHits;

    static int compareByContigAndScore(const void *first, const void *second);      // qsort()-style compare routine
    static int compareByScore(const void *first, const void *second);               // qsort()-style compare routine
};