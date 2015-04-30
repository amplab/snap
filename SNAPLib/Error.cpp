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
#include "CommandProcessor.h"

	void
WriteMessageToFile(FILE *file, const char *message)
{
    if (AlignerOptions::useHadoopErrorMessages) {
        fprintf(stderr,"reporter:status:%s", message);           // Always use stderr in Hadoop mode, regardless of whether this is an error
        fprintf(stderr, "%s", message);                          // And also print without the prefix, so it shows up in both logs
		fflush(stderr);
    } else {
        if (AlignerOptions::outputToStdout && stdout == file) {
	  fprintf(stderr, "%s", message);
			fflush(stderr);
        } else {
	  fprintf(file, "%s", message);
			fflush(file);
        }
    }
}


    void
WriteErrorMessage(const char *message, ...)
{
    va_list args;
    va_start(args, message);
    const size_t bufferSize = 10240;
    char buffer[bufferSize];
    vsnprintf(buffer, bufferSize - 1, message, args);
    WriteMessageToFile(stderr, buffer);
	if (NULL != CommandPipe) {
	  WriteToNamedPipe(CommandPipe, buffer);
	}
}

    void
WriteStatusMessage(const char *message, ...)
{
    va_list args;
    va_start(args, message);
    const size_t bufferSize = 10240;
    char buffer[bufferSize];
    vsnprintf(buffer, bufferSize - 1, message, args);
    WriteMessageToFile(stdout, buffer);
	if (NULL != CommandPipe) {
	  WriteToNamedPipe(CommandPipe, buffer);
	}
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
	fflush(stderr);
}
