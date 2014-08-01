/*++

Module Name:

AlignmentResult.cpp

Abstract:

SNAP genome alignment results

Authors:

Bill Bolosky & Xuya Wang, July, 2014

Environment:

User mode service.

This class is NOT thread safe.  It's the caller's responsibility to ensure that
at most one thread uses an instance at any time.

Revision History:

--*/

#include "stdafx.h"
#include "AlignmentResult.h"
#include "Read.h"

GenomeLocation correctLocationForSoftClipping(AlignmentResult status, GenomeLocation location, Direction direction, Read *read, const Genome *genome)
{
	if (NotFound == status || UnknownAlignment == status) {
		return location;
	}

	GenomeLocation newLocation;
	if (FORWARD == direction && read->getOriginalFrontClipping() != 0) {
		newLocation = location - read->getOriginalFrontClipping();
	} else if (RC == direction && read->getOriginalBackClipping() != 0) {
		newLocation = location - read->getOriginalBackClipping();
	} else {
		return location;
	}

	if (GenomeLocationAsInt64(newLocation) < 0) {
		return location;
	}

	const Genome::Contig *contig = genome->getContigAtLocation(location);
	if (contig->beginningLocation > newLocation) {
		//
		// This would push the alignment into a different contig.  Skip the fixup.
		//
		return location;
	}

	return newLocation;
}

	void
SingleAlignmentResult::correctAlignmentForSoftClipping(Read *read, const Genome *genome)
{
	location = correctLocationForSoftClipping(status, location, direction, read, genome);
}

	void
PairedAlignmentResult::correctAlignmentForSoftClipping(Read *read[], const Genome *genome)
{
	for (Direction direction = FORWARD; direction < NUM_DIRECTIONS; direction++) {
		location[direction] = correctLocationForSoftClipping(status[direction], location[direction], direction, read[direction], genome);
	}
}
