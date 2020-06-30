/*++

Module Name:

    Main.cpp

Abstract:

   Just call into the command processor from SNAPLib

Authors:

    Matei Zaharia & Bill Bolosky, February, 2012

Environment:

    User mode service.

Revision History:

    Adapted from cSNAP, which was in turn adapted from the scala prototype

--*/

#include "stdafx.h"
#include "CommandProcessor.h"

//The SNAP_VERSION string moved to SNAPLib\CommandProcessor.cpp

_int64 BJBProgramStartTimeInMillis;
int main(int argc, const char **argv)
{
    BJBProgramStartTimeInMillis = timeInMillis();
	ProcessTopLevelCommands(argc, argv);
    fprintf(stderr, "Total run time %lld\n", (timeInMillis() + 500 - BJBProgramStartTimeInMillis) / 1000);
}
