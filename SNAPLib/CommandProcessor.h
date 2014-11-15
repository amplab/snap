/*++

Module Name:

CommandProcessor.h

Abstract:

Header for running the top-level commands of SNAP

Authors:

Bill Bolosky, November, 2014

Environment:
`
User mode service.

Revision History:

--*/

#pragma once
#include "Compat.h"

extern void ProcessTopLevelCommands(int argc, const char **argv);

extern NamedPipe *CommandPipe;
extern const char *CommandExecutedString;	// Sent back along the command pipe to indicate that the whole thing is done and SNAPCommand should exit
