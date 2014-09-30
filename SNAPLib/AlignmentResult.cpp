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
	if (FORWARD == direction && read->getFrontClippedLength() != 0) {
		newLocation = location - read->getFrontClippedLength();
	} else if (RC == direction && read->getBackClippedLength() != 0) {
		newLocation = location - read->getBackClippedLength();
	} else {
		return location;
	}

	if (GenomeLocationAsInt64(newLocation) < 0) {
		return location;
	}

	const Genome::Contig *contig = genome->getContigAtLocation(location);
	if (contig == NULL || contig->beginningLocation > newLocation) {
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
	for (unsigned whichRead = 0; whichRead < NUM_READS_PER_PAIR; whichRead++) {
		location[whichRead] = correctLocationForSoftClipping(status[whichRead], location[whichRead], direction[whichRead], read[whichRead], genome);
	}
}
