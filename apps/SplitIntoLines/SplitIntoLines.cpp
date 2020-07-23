// SplitIntoLines.cpp : This file contains the 'main' function. Program execution begins and ends there.
//

#include <iostream>
#include <Windows.h>

void usage()
{
    fprintf(stderr, "usage: SplitIntoLines -nLines inputFile {outputFileBase}\n");
    exit(1);
}



int main(int argc, char **argv)
{
    if (argc != 3 && argc != 4) 
    {
        usage();
    }

    if (argv[1][0] != '-')
    {
        usage();
    }

    int nLinesPerChunk = atoi(argv[1] + 1);
    if (nLinesPerChunk < 1) {
        fprintf(stderr, "Must be at least one line per output file\n");
        exit(1);
    }

    HANDLE hFile = CreateFile(argv[2], GENERIC_READ, FILE_SHARE_READ, NULL, OPEN_EXISTING, FILE_FLAG_SEQUENTIAL_SCAN, NULL);
    if (INVALID_HANDLE_VALUE == hFile) 
    {
        fprintf(stderr, "Unable to open '%s', %d\n", argv[2], GetLastError());
        exit(1);
    }

    HANDLE hMapping = CreateFileMapping(hFile, NULL, PAGE_READONLY, 0, 0, NULL);
    if (NULL == hMapping) 
    {
        fprintf(stderr, "CreateFileMapping failed, %d\n", GetLastError());
        exit(1);
    }

    PVOID fileData = MapViewOfFile(hMapping, FILE_MAP_READ, 0, 0, 0);
    if (NULL == fileData) 
    {
        fprintf(stderr, "MapViewOfFile failed, %d\n", GetLastError());
        exit(1);
    }

    LARGE_INTEGER liFileSize;
    if (!GetFileSizeEx(hFile, &liFileSize)) {
        fprintf(stderr, "GetFileSizeEx failed, %d\n", GetLastError());
        exit(1);
    }

    _int64 fileSize = liFileSize.QuadPart;
    char* endOfFile = (char*)fileData + fileSize;

    char* linePointer = (char *)fileData;

    int outputFileNumber = 0;

    const char* outputFileBaseName = "x";

    if (argc > 3) {
        outputFileBaseName = argv[3];
    }

    int nOutputFiles = 0;

    while (linePointer < endOfFile) {

        char* beginningOfThisChunk = linePointer;

        for (int i = 0; i < nLinesPerChunk; i++)
        {
            while (linePointer < endOfFile && *(linePointer++) != '\n') {}  // There's certainly a faster way to do this.
        }

        size_t outputFilenameBufferSize = strlen(outputFileBaseName) + 10; // 1 for . 8 for digits and 1 for null termination
        char* outputFilename = new char[outputFilenameBufferSize];
        sprintf_s(outputFilename, outputFilenameBufferSize, "%s.%08d", outputFileBaseName, nOutputFiles);

        HANDLE hOutputFile = CreateFile(outputFilename, GENERIC_WRITE, 0, NULL, CREATE_ALWAYS, FILE_FLAG_SEQUENTIAL_SCAN, NULL);
        if (INVALID_HANDLE_VALUE == hOutputFile) {
            fprintf(stderr, "Unable to open '%s', %d\n", outputFilename, GetLastError());
            exit(1);
        }

        nOutputFiles++;

        char* writePointer = beginningOfThisChunk;

        while (writePointer < linePointer) 
        {
            DWORD amountToWrite = (int)(__min(1024 * 1024 * 1024, linePointer - writePointer));
            DWORD amountWritten;
            if (!WriteFile(hOutputFile, writePointer, amountToWrite, &amountWritten, NULL)) {
                fprintf(stderr, "WriteFile failed on '%s', %d\n", outputFilename, GetLastError());
                exit(1);
            }

            if (amountWritten < 1) {
                fprintf(stderr, "WriteFile didn't write any data\n");
                exit(1);
            }

            writePointer += amountWritten;
        }

        CloseHandle(hOutputFile);
        delete[] outputFilename;

    } // while we have data.

    UnmapViewOfFile(fileData);
    CloseHandle(hMapping);
    CloseHandle(hFile);
}

