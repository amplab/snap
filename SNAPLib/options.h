/*++

Module Name:

    options.h

Abstract:

    Compile time options for cSNAP

Authors:

    Bill Bolosky, September, 2011

Environment:

    User mode service.

    This class is NOT thread safe.  It's the caller's responsibility to ensure that
    at most one thread uses an instance at any time.

Revision History:

    Adapted from Matei Zaharia's Scala implementation.

--*/

#pragma once

#define MAINTAIN_HISTOGRAMS     0

#ifndef USE_DEVTEAM_OPTIONS
#define USE_DEVTEAM_OPTIONS     1 // Options only for the development team, not for the release version
#endif
