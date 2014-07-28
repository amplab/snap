/*++

Module Name:

    exit.cpp

Abstract:

    SNAP soft exit function

Authors:

    Bill Bolosky, February, 2013

Environment:

    User mode service.

Revision History:


--*/

#include "stdafx.h"
#include "exit.h"
#include "Error.h"

//
// This exists solely as a place to set a breakpoint when debugging SNAP.  It gets called both from the
// soft_exit function and also from a few places in the code where we want to exit without printing the
// warning message (like after printing the usage string), which was causing confusion.
//
void soft_exit_no_print(int n)
{
    exit(n);
}

void soft_exit_function(int n, const char *fileName, int lineNum)
{
    WriteErrorMessage("SNAP exited with exit code %d from line %d of file %s\n", n, lineNum, fileName);
    soft_exit_no_print(n);
}
