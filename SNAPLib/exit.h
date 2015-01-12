/*++

Module Name:

    exit.h

Abstract:

    Header for SNAP soft exit function

Authors:

    Bill Bolosky, February, 2013

Environment:

    User mode service.

Revision History:


--*/

#pragma once

#define soft_exit(n) soft_exit_function(n, __FILE__, __LINE__)

void soft_exit_no_print(int n);

void soft_exit_function(int n, const char *fileName, int lineNum);
