//
// FixSplitNamesorted.cpp : Splitting namesorted SAM files doens't result in exact matches because of
// things like supplementary alignments.  This program takes a set of mostly-identical splits
// and produces overflow files for each input that contain the reads that don't match.
//

#include <iostream>
#include <windows.h>

void usage()
{
	fprintf(stderr, "usage: FixSplitNamesorted input1 input2 ... inputN\n");
	fprintf(stderr, "'input' is the base name of a namesorted, split SAM file.  So, for instance, \\machine\d$\dir\hg001.snap.namesorted\n");
	fprintf(stderr, "This program will then read the inputs with .nnnnnnnn appended to the names and write into outputs\n");
	fprintf(stderr, "<file>.extra any reads that aren't in all of the inputs\n");
	exit(1);
}

int main(int argc, const char **argv)
{
	if (argc < 3) usage();

	int nInputs = argc - 1;

	HANDLE* inputHandles = new HANDLE[nInputs];
	HANDLE* outputHandles = new HANDLE[nInputs];

	bool* inputsCompleted = new bool[nInputs];

	for (int i = 0; i < nInputs; i++) {
		size_t outputFilenameLen = strlen(argv[i + 1]) + 7;// +7 is ".extra" and \0
		char* outputFilename = new char[outputFilenameLen];
		sprintf_s(outputFilename, outputFilenameLen, "%s.extra", argv[i - 1]);

		outputHandles[i] = CreateFile(outputFilename, GENERIC_READ, FILE_SHARE_READ, NULL, CREATE_ALWAYS, 0, NULL);

		if (INVALID_HANDLE_VALUE == outputHandles[i]) {
			fprintf(stderr, "Unable to open '%s', %d\n", outputFilename, GetLastError());
			exit(1);
		}

		delete[] outputFilename;

		inputsCompleted[i] = false;
	} // for each input 

	for (int fileIndex = 0; true; fileIndex++) { // Exit is from mid-loop when there are no more input files.
		for (int i = 0; i < nInputs; i++) {

		}
	}
}
