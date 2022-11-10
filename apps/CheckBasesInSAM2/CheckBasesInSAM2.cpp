#include <stdio.h>
#include <Windows.h>


void usage()
{
	fprintf(stderr, "usage: CheckBasesInSAM2 inputfile.sam\n");
	exit(1);
} // usage

HANDLE hInputFile;
__int64 fileSize;
volatile LONG64 TotalLinesProcessed = 0;
const LONG64 LinesPerDot = 10000000;
volatile LONG nThreadsRunning = 0;
HANDLE hDone;


struct ThreadContext {
	__int64			startingOffset;
	__int64			amountToProcess;
};

const int ChunkSize = 128 * 1024 * 1024;	// Map and process 128MB at a time

bool ValidCharacters[256];

void WorkerThread(PVOID pvContext) 
{
	ThreadContext* context = (ThreadContext*)pvContext;

	__int64 amountRemaining = context->amountToProcess;

	while (context->amountToProcess > 0) {
		__int64 linesProcessed = 0;

		int amountThisChunk = (int)__min(ChunkSize, amountRemaining);	// cast OK because of min
		int amountToMap = (int)__min(amountThisChunk + 64 * 1024, fileSize - context->startingOffset - (context->amountToProcess - amountRemaining));	// To allow for a line that runs off the end of the mapped chunk.

		HANDLE hMapping = CreateFileMapping(hInputFile, NULL, PAGE_READONLY, 0, 0, NULL);
		if (NULL == hMapping) {
			fprintf(stderr, "CreateFileMapping failed, %d\n", GetLastError());
			exit(1);
		}

		LARGE_INTEGER liOffset;
		liOffset.QuadPart = context->startingOffset + context->amountToProcess - amountRemaining;

		const char *mappedFile = (const char *)MapViewOfFile(hMapping, FILE_MAP_READ, liOffset.HighPart, liOffset.LowPart, 
			(liOffset.QuadPart + amountToMap) >= fileSize ? 0 : amountToMap);
		if (NULL == mappedFile) {
			fprintf(stderr, "MapViewOfFile failed, %d\n", GetLastError());
			exit(1);
		}

		const char* currentChar = mappedFile;
		if (liOffset.QuadPart != 0) {
			//
			// Go to the first newline so we start at the beginning of a line.
			//
			while (currentChar < mappedFile + amountToMap && *currentChar != '\n') {
				currentChar++;
			}
			currentChar++;
		} // If we're not at the beginning of the file

		while (currentChar < mappedFile + amountThisChunk) {	// Run to the size of the chunk, not to the amount mapped
			const char* lineStart = currentChar;
			bool printLine = false;
			char badChar = ' ';

			if (*currentChar != '@') {	// We just skip header lines by falling through 
				int nTabsSeen = 0;
				while (currentChar < mappedFile + amountToMap && nTabsSeen < 9) {
					if (*currentChar == '\t') {
						nTabsSeen++;
					}
					currentChar++;
				}

				while (currentChar < mappedFile + amountToMap) {
					if (*currentChar == '\t') {
						break;
					}

					if (!ValidCharacters[*currentChar]) {
						printLine = true;
						badChar = *currentChar;
						break;
					}

					currentChar++;
				} // while we're in SEQ
			} // If it's not a header line

			//
			// Now go to the end of the line
			while (currentChar < mappedFile + amountToMap && *currentChar != '\n') {
				currentChar++;
			}
			currentChar++;

			if (printLine) {
				//
				// Probably should lock here, but with luck this is so rare that it doesn't matter
				//
				fprintf(stderr, "\nFound line with bad seq char '%c': ", badChar);
				for (const char* printChar = lineStart; printChar < currentChar; printChar++) {
					fprintf(stderr, "%c", *printChar);
				}
				fprintf(stderr, "\n");
			}

			linesProcessed++;
		} // While we have lines in this chunk to handle

		if (!UnmapViewOfFile(mappedFile)) {
			fprintf(stderr, "UnmapViewOfFile failed, %d\n", GetLastError());
			exit(1);
		}

		CloseHandle(hMapping);

		LONG64 newLinesProcessed = InterlockedAdd64(&TotalLinesProcessed, linesProcessed);
		if (newLinesProcessed / LinesPerDot != (newLinesProcessed - linesProcessed) / LinesPerDot) {
			printf(".");
		}

		amountRemaining -= amountThisChunk;
	} // while we have chunks to do

	if (0 == InterlockedDecrement(&nThreadsRunning)) {
		SetEvent(hDone);
	}
} // WorkerThread

int main(int argc, const char** argv)
{
	if (argc != 2) {
		usage();
	}

	for (int i = 0; i < 256; i++) {
		ValidCharacters[i] = false;
	}

	ValidCharacters['A'] = true;
	ValidCharacters['a'] = true;
	ValidCharacters['T'] = true;
	ValidCharacters['t'] = true;
	ValidCharacters['C'] = true;
	ValidCharacters['c'] = true;
	ValidCharacters['G'] = true;
	ValidCharacters['g'] = true;
	ValidCharacters['N'] = true;
	ValidCharacters['n'] = true;
	ValidCharacters['U'] = true;
	ValidCharacters['u'] = true;
	ValidCharacters['R'] = true;
	ValidCharacters['r'] = true;
	ValidCharacters['Y'] = true;
	ValidCharacters['y'] = true;
	ValidCharacters['K'] = true;
	ValidCharacters['k'] = true;
	ValidCharacters['M'] = true;
	ValidCharacters['m'] = true;
	ValidCharacters['S'] = true;
	ValidCharacters['s'] = true;
	ValidCharacters['W'] = true;
	ValidCharacters['w'] = true;
	ValidCharacters['B'] = true;
	ValidCharacters['b'] = true;
	ValidCharacters['D'] = true;
	ValidCharacters['d'] = true;
	ValidCharacters['H'] = true;
	ValidCharacters['h'] = true;
	ValidCharacters['V'] = true;
	ValidCharacters['v'] = true;

	const char* inputFilename = argv[1];
	hInputFile = CreateFileA(inputFilename, GENERIC_READ, FILE_SHARE_READ, NULL, OPEN_EXISTING, FILE_FLAG_SEQUENTIAL_SCAN, NULL);

	if (INVALID_HANDLE_VALUE == hInputFile) {
		fprintf(stderr, "Open of '%s' failed, %d\n", inputFilename, GetLastError());
		return 1;
	}

	BY_HANDLE_FILE_INFORMATION fileInfo[1];
	if (!GetFileInformationByHandle(hInputFile, fileInfo)) {
		fprintf(stderr, "GetFileInformationBYHandle failed, %d\n", GetLastError());
		return 1;
	}

	LARGE_INTEGER liFileSize;
	liFileSize.HighPart = fileInfo->nFileSizeHigh;
	liFileSize.LowPart = fileInfo->nFileSizeLow;
	fileSize = liFileSize.QuadPart;

	hDone = CreateEvent(NULL, TRUE, FALSE, NULL);
	if (NULL == hDone) {
		fprintf(stderr, "Unable to CreateEvent, %d\n", GetLastError());
	}

	printf("Progress (1 dot/10M lines): ");

	//
	// One thread for debugging.
	//
	ThreadContext context[1];
	context->amountToProcess = fileSize;
	context->startingOffset = 0;

	WorkerThread(context);

	CloseHandle(hInputFile);


} // main