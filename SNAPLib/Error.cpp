/*++

Module Name:

    Error.cpp

Abstract:

    SNAP error-message writer

Authors:

    Bill Bolosky, Feburary, 2014

Environment:

    User mode service.

Revision History:


--*/

#include "stdafx.h"
#include "Compat.h"
#include "Error.h"
#include "AlignerOptions.h"

void
WriteMessageToFile(FILE *file, const char *message, va_list args)
{
    if (AlignerOptions::useHadoopErrorMessages) {
        const size_t bufferSize = 10240;
        char buffer[bufferSize];
        vsnprintf(buffer, bufferSize - 1, message, args);
        fprintf(stderr,"reporter:status:%s", buffer);           // Always use stderr in Hadoop mode, regardless of whether this is an error
        fprintf(stderr, "%s", buffer);                          // And also print without the prefix, so it shows up in both logs
    } else {
        if (AlignerOptions::outputToStdout && stdout == file) {
            vfprintf(stderr, message, args);
        } else {
            vfprintf(file, message, args);
        }
    }
}

    void
WriteErrorMessage(const char *message, ...)
{
    va_list args;
    va_start(args, message);
    WriteMessageToFile(stderr, message, args);
}

    void
WriteStatusMessage(const char *message, ...)
{
    va_list args;
    va_start(args, message);
    WriteMessageToFile(stdout, message, args);
}

    void
WriteProgressCounter(const char *counterName, _int64 increment)
{
    if (!AlignerOptions::useHadoopErrorMessages) {
        //
        // No counters unless in Hadoop mode.
        //
        return;
    }

    fprintf(stderr,"reporter:counter:SNAP,%s,%lld\n", counterName, increment);
}