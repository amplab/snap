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
