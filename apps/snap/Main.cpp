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

int main(int argc, const char **argv)
{
	ProcessTopLevelCommands(argc, argv);
}
