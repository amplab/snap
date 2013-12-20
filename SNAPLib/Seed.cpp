/*++

Module Name:

    Seed.cpp

Abstract:

    Code to handle seeds in the SNAP sequencer.

Authors:

    Bill Bolosky, August, 2011

Environment:
`
    User mode service.

Revision History:

    Adapted from Matei Zaharia's Scala implementation.

--*/

#include "stdafx.h"
#include "Seed.h"

    bool
Seed::DoesTextRepresentASeed(const char *textBases, unsigned seedLen)
{
    for (unsigned i = 0; i < seedLen; i++) {
        switch (textBases[i]) {
            case 'A':
            case 'G':
            case 'C':
            case 'T':
                    break;
            default: return false;
        }
    }
    return true;
}
    
    Seed
Seed::fromBases(
    _int64 bases,
    int seedLength)
{
    _int64 rc = 0;
    _int64 b = bases;
    for (int i = 0; i < seedLength; i++) {
        rc = (rc << 2) | ((b & 3) ^ 3);
        b = b >> 2;
    }
    return Seed(bases, rc);
}