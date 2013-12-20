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

void soft_exit_function(int n, const char *fileName, int lineNum)
{
    fprintf(stderr,"SNAP exited with exit code %d from line %d of file %s\n", n, lineNum, fileName);
    exit(n);
}
