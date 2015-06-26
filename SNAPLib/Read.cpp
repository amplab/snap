/*++


Module Name:

    Read.cpp

Abstract:

   Read class for the SNAP sequencer

Authors:

    Bill Bolosky, May, 2012

Environment:

    User mode service.

Revision History:

--*/

#include "stdafx.h"
#include "Read.h"
#include "SAM.h"
#include "Error.h"

	bool
isAValidAlignmentResult(
	AlignmentResult result)
{
	return result == NotFound || result == SingleHit || result == MultipleHits || result == UnknownAlignment;
}

	void
Read::checkIdMatch(
	Read* read0,
	Read* read1)
{
    if (!readIdsMatch(read0, read1)) {
        unsigned n[2] = {min(read0->getIdLength(), 200u), min(read1->getIdLength(), 200u)};
        char* p[2] = {(char*) alloca(n[0] + 1), (char*) alloca(n[1] + 1)};
        memcpy(p[0], read0->getId(), n[0]); p[0][n[0]] = 0;
        memcpy(p[1], read1->getId(), n[1]); p[1][n[1]] = 0;
        WriteErrorMessage("Unmatched read IDs '%s' and '%s'.  Use the -I option to ignore this.\n", p[0], p[1]);
        soft_exit(1);
    }
}

    
const unsigned Read::localBufferLength = 3 * MAX_READ_LENGTH;
const unsigned DEFAULT_MIN_READ_LENGTH = 50;
