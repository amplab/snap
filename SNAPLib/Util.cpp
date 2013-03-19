/*++

Module Name:

    Util.cpp

Abstract:

    Generic support routines that don't seem to belong elsewhere.

Authors:

    Bill Bolosky, March, 2013

Environment:

    User mode service.

Revision History:

    Factored from other places

--*/

#include "stdafx.h"
#include "Util.h"


int FirstPowerOf2GreaterThanOrEqualTo(int value)
{
    int highestBitSet;
    for (highestBitSet = 0; highestBitSet <= 30; highestBitSet++) { // Only go to 31, since this is signed
        if (!(value & ~((1 << highestBitSet) - 1))) {
            highestBitSet -= 1;
            break;
        }
    }
    if (1 << highestBitSet == value) return value;
    return 1 << (highestBitSet + 1);
}
