/*++

Module Name:

CommandProcessor.cpp

Abstract:

Code for running the top-level commands of SNAP

Authors:

Bill Bolosky, November, 2014

Environment:
`
User mode service.

Revision History:

Pulled from the main program and expanded to handle daemon mode

--*/

#include "stdafx.h"
#include "options.h"
#include "FASTA.h"
#include "GenomeIndex.h"
#include "SingleAligner.h"
#include "PairedAligner.h"
#include "exit.h"
#include "SeedSequencer.h"
#include "AlignerOptions.h"
#include "CommandProcessor.h"
#include "Error.h"
#include "Compat.h"

const char *SNAP_VERSION = "1.0beta.23";

static void usage()
{
	WriteErrorMessage(
		"Usage: snap-aligner <command> [<options>]\n"
		"Commands:\n"
		"   index    build a genome index\n"
		"   single   align single-end reads\n"
		"   paired   align paired-end reads\n"
		"   daemon   run in daemon mode--accept commands remotely\n"
		"Type a command without arguments to see its help.\n");
}

void ProcessNonDaemonCommands(int argc, const char **argv) {
	if (strcmp(argv[1], "index") == 0) {
		if (CommandPipe == NULL) {
			GenomeIndex::runIndexer(argc - 2, argv + 2);
		} else {
			//
			// The error cases in index build don't really free memory properly, so we just don't allows it in daemon mode.
			//
			WriteErrorMessage("The index command is not available in daemon mode.  Please run 'snap-aligner index' directly.\n");
		}
	} else if (strcmp(argv[1], "single") == 0 || strcmp(argv[1], "paired") == 0) {
		for (int i = 1; i < argc; /* i is increased below */) {
			unsigned nArgsConsumed;
			if (strcmp(argv[i], "single") == 0) {
				SingleAlignerContext single;
				single.runAlignment(argc - i, argv + i, SNAP_VERSION, &nArgsConsumed);
			} else if (strcmp(argv[i], "paired") == 0) {
				PairedAlignerContext paired;
				paired.runAlignment(argc - i, argv + i, SNAP_VERSION, &nArgsConsumed);
			} else {
				fprintf(stderr, "Invalid command: %s\n\n", argv[i]);
				usage();
				return;
			}
			_ASSERT(nArgsConsumed > 0);
			i += nArgsConsumed;
		}
	} else {
		WriteErrorMessage("Invalid command: %s\n\n", argv[1]);
		usage();
	}
}

static void daemonUsage()
{
	fprintf(stderr, "Usage: snap-aligner daemon [Named pipe name]\n");
	soft_exit_no_print(1);    // Don't use soft_exit, it's confusing people to get an "error" message after the usage
}

void RunDaemonMode(int argc, const char **argv)
{
	if (argc < 2 || argc > 3) {
		daemonUsage();
	}

	printf("SNAP in daemon mode, waiting for commands to execute\n");

	const char *pipeName = argc == 3 ? argv[2] : DEFAULT_NAMED_PIPE_NAME;
	CommandPipe = OpenNamedPipe(pipeName, true);

	if (NULL == CommandPipe) {
		WriteErrorMessage("Unable to open named pipe for command IO.\n");
		soft_exit(1);
	}

	const size_t commandBufferSize = 10000;	// Yes, this is fixed size, no it's not a buffer overflow.  The named pipe reader just quits if it's too long.
	char commandBuffer[commandBufferSize];

	//
	// Format of commands is argc (in ascii) followed by argc arguments, each in one line.
	//
	for (;;) {
		if (!ReadFromNamedPipe(CommandPipe, commandBuffer, commandBufferSize)) {
			CloseNamedPipe(CommandPipe);
			CommandPipe = NULL;
			WriteStatusMessage("Named pipe closed.  Exiting\n");
			soft_exit_no_print(0);
		}

		int argc = atoi(commandBuffer);
		if (0 == argc) {
			WriteErrorMessage("Expected argument count on named pipe, got '%s'; ignoring.\n", commandBuffer);
		} else {
			char **argv = new char*[argc];
			for (int i = 0; i < argc; i++) {
				argv[i] = new char[commandBufferSize];
				if (!ReadFromNamedPipe(CommandPipe, argv[i], commandBufferSize)) {
					CloseNamedPipe(CommandPipe);
					CommandPipe = NULL;
					WriteStatusMessage("Error reading argument #%d from named pipe.\n", i);
					soft_exit(1);
				}
			} // for each arg

			if (argc > 1 && strcmp(argv[1], "exit") == 0) {
				WriteStatusMessage("SNAP server exiting by request\n");
				WriteToNamedPipe(CommandPipe, CommandExecutedString);
				soft_exit_no_print(1);
			}

			printf("Executing command: ");
			for (int i = 1; i < argc; i++) {
				printf("%s ", argv[i]);
			}
			printf("\n");

			ProcessNonDaemonCommands(argc, (const char **) argv);

			printf("\n");

			for (int i = 0; i < argc; i++) {
				delete[] argv[i];
				argv[i] = NULL;
			}
			delete[] argv;
			argv = NULL;
		}
		WriteToNamedPipe(CommandPipe, CommandExecutedString);
	}
}

void ProcessTopLevelCommands(int argc, const char **argv)
{
	fprintf(stderr, "Welcome to SNAP version %s.\n\n", SNAP_VERSION);       // Can't use WriteStatusMessage, because we haven't parsed args yet to determine if -hdp is specified.  Just stick with stderr.

	InitializeSeedSequencers();

	if (argc < 2) {
		usage();
		soft_exit_no_print(1);
	}

	if (strcmp(argv[1], "daemon") == 0) {
		RunDaemonMode(argc, argv);
	} else {
		ProcessNonDaemonCommands(argc, argv);
	}
}

NamedPipe *CommandPipe = NULL;
const char *CommandExecutedString = "***SNAP Command completed execution***";
