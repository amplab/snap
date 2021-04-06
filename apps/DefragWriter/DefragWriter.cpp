// DefragWriter.cpp : This file contains the 'main' function. Program execution begins and ends there.
//

#include <iostream>
#include <Windows.h>

int main(int argc, char **argv)
{
	if (argc != 2) {
		fprintf(stderr, "usage: DefragWriter outputFilename\n");
		fprintf(stderr, "copies stdin to the output file writing in 1 GB chunks\n");
		return 1;
	}

	HANDLE hInputFile = GetStdHandle(STD_INPUT_HANDLE);
	if (hInputFile == INVALID_HANDLE_VALUE) {
		fprintf(stderr, "GetStdHandle(STD_INPUT_HANDLE) failed, %d\n", GetLastError());
		return GetLastError();
	}

	HANDLE hOutputFile = CreateFileA(argv[1], GENERIC_WRITE, FILE_SHARE_READ, NULL, CREATE_ALWAYS, 0, NULL);
	if (hOutputFile == INVALID_HANDLE_VALUE) {
		fprintf(stderr, "Unable to open %s, %d\n", argv[1], GetLastError());
		return GetLastError();
	}

	const DWORD bufferSize = 1024 * 1024 * 1024;
	char* buffer = (char*)malloc(bufferSize);

	DWORD bufferOffset = 0;

	DWORD bytesRead;

	while (ReadFile(hInputFile, buffer + bufferOffset, bufferSize - bufferOffset, &bytesRead, NULL)) {
		bufferOffset += bytesRead;
		if (bufferOffset == bufferSize) {
			DWORD writeOffset = 0;

			while (writeOffset < bufferSize) {
				DWORD bytesWritten;
				if (!WriteFile(hOutputFile, buffer + writeOffset, bufferSize - writeOffset, &bytesWritten, NULL)) {
					fprintf(stderr, "Error writing output file: %d\n", GetLastError());
					return GetLastError();
				}

				writeOffset += bytesWritten;
			} // while we have something to write

			bufferOffset = 0;
		} // if we have a full buffer
	} // while we read stuff

	if (bufferOffset != 0) {
		//
		// Write the rump of the buffer
		//
		DWORD writeOffset = 0;

		while (writeOffset < bufferOffset) {
			DWORD bytesWritten;
			if (!WriteFile(hOutputFile, buffer + writeOffset, bufferOffset - writeOffset, &bytesWritten, NULL)) {
				fprintf(stderr, "Error writing output file: %d\n", GetLastError());
				return GetLastError();
			}

			writeOffset += bytesWritten;
		} // while we have something to write
	} // if we have leftover stuff in the buffer

	CloseHandle(hOutputFile);
}

