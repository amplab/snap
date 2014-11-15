/*++

Module Name:

SNAPCommand.cpp

Abstract:

Send a command to SNAP running daemon mode

Authors:

Matei Zaharia & Bill Bolosky, February, 2012

Environment:

User mode service.

Revision History:

Adapted from cSNAP, which was in turn adapted from the scala prototype

--*/

#include "stdafx.h"
#include "Compat.h"
#include "exit.h"
#include "CommandProcessor.h"

void
usage()
{
	fprintf(stderr, "usage: SNAPCommand {-p PipeName} <command to send to SNAP>\n");
	fprintf(stderr, "Send command 'exit' to SNAP to have the server process exit.\n");
	soft_exit_no_print(1);
}


int main(int argc, const char **argv) 
{
	if (argc < 2) {
		usage();
	}

	const char *pipeName;
	int startingArg;

	if (strcmp(argv[1], "-p") == 0) {
		if (argc < 4) usage();

		pipeName = argv[2];
		startingArg = 3;
	} else {
		pipeName = DEFAULT_NAMED_PIPE_NAME;
		startingArg = 1;
	}

	NamedPipe *serverPipe = OpenNamedPipe(pipeName, false);

	if (NULL == serverPipe) {
		fprintf(stderr, "Unable to open pipe to server\n");
		soft_exit(1);
	}

	//
	// Send down the args.
	//
	char argcBuffer[100];
	sprintf(argcBuffer, "%d", argc - startingArg + 1);	// +1 is for the command name, argv[0]
	if (!WriteToNamedPipe(serverPipe, argcBuffer)) {
		fprintf(stderr, "Unable to send arg count to server\n");
		soft_exit(1);
	}

	if (!WriteToNamedPipe(serverPipe, argv[0])) {
		fprintf(stderr, "Error sending arg '%s' to server\n", argv[0]);
		soft_exit(1);
	}

	for (int i = startingArg; i < argc; i++) {
		if (!WriteToNamedPipe(serverPipe, argv[i])) {
			fprintf(stderr, "Error sending arg '%s' to server\n", argv[i]);
			soft_exit(1);
		}
	}

	//
	// Now process the results from SNAP, printing them out and waiting for the terminator.
	//
	const size_t outputBufferSize = 100000;
	char outputBuffer[outputBufferSize];

	while (ReadFromNamedPipe(serverPipe, outputBuffer, outputBufferSize)) {
		if (strcmp(outputBuffer, CommandExecutedString) == 0) {
			soft_exit_no_print(0);
		}
		printf("%s", outputBuffer);
		fflush(stdout);
	}

	fprintf(stderr, "Error reading from server pipe\n");
	soft_exit(1);

	return 1;
}
