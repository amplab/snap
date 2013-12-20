/*++

Module Name:

    directions.h

Abstract:

    Definitions for basic read directions (forward & reverse compliment)

Authors:

    Bill Bolosky, January, 2013

Environment:

    User mode service.

Revision History:

    
--*/

#pragma once

const int NUM_DIRECTIONS = 2;    // Forward and reverse compliment

typedef int Direction;

const int FORWARD = 0;
const int RC = 1;

inline Direction OppositeDirection(Direction direction) {
    _ASSERT(FORWARD == direction || RC == direction);
    return 1-direction;
}