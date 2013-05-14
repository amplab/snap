/*++


Module Name:

	GoodRandom.h

Abstract:

    Cryptographically secure random number generator that's fast because it batches getting its random
    source from CryptGenRand().  It's careful not to bias results when the request range isn't a power
    of two.

Authors:

    Bill Bolosky, April, 2011

Environment:

    Not thread safe.

--*/
#pragma once

#include "Compat.h"

unsigned MinBytesToStore(_uint64 maxValue);

//
// This returns a cryptographically secure random number from {0..maxValue} (NOT {0..maxValue-1}!)
// It's computationally efficient, too.
//
_uint64 GoodFastRandom(_uint64 maxValue);
