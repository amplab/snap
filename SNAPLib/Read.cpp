/*++


Module Name:

    Read.h

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


bool isAValidAlignmentResult(AlignmentResult result) {return result == NotFound || result == SingleHit || result == MultipleHits || result == UnknownAlignment;}
