/*++

Module Name:

    Error.h

Abstract:

    Header for SNAP error-message writer

Authors:

    Bill Bolosky, Feburary, 2014

Environment:

    User mode service.

Revision History:


--*/

#pragma once
#include "Compat.h"

    void
WriteErrorMessage(const char *message, ...);

    void
WriteStatusMessage(const char *message, ...);

    void
WriteProgressCounter(const char *counterName, _int64 increment);