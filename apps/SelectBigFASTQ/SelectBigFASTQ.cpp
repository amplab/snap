// SelectBigFASTQ.cpp : This file contains the 'main' function. Program execution begins and ends there.
//

#define _CRT_RAND_S
#define _CRT_SECURE_NO_WARNINGS
#include <iostream>
#include <Windows.h>
#include <stdlib.h>

void usage()
{
    fprintf(stderr, "usage: SelectBigFASTQ nReadsToWrite inputFile1 {inputFile2} outputFile1 {outputFile2}\n");
    fprintf(stderr, "       nReadsToWrite is per output file.\n");
    exit(1);
}

int main(int argc, char **argv)
{
    if (argc != 4 && argc != 6) usage();

    int nReadsToWrite = atoi(argv[1]);
    if (nReadsToWrite < 1) usage();

    int nInputFiles = (argc == 4) ? 1 : 2;

    const char** inputFileNames = new const char* [nInputFiles];
    inputFileNames[0] = argv[2];
    
    const char** outputFileNames = new const char* [nInputFiles];

    if (nInputFiles == 1) {
        outputFileNames[0] = argv[3];
    } else {
        inputFileNames[1] = argv[3];
        outputFileNames[0] = argv[4];
        outputFileNames[1] = argv[5];
    }

    //
    // Estimate the number of reads by reading the first 10 and counting the bytes.
    //

    size_t nBytesInFirst10Reads = 0;
    FILE* inputForEstimate = fopen(inputFileNames[0], "r");
    const int lineBufferSize = 10240;
    char lineBuffer[lineBufferSize];

    for (int i = 0; i < 40; i++) {
        fgets(lineBuffer, lineBufferSize, inputForEstimate);
        nBytesInFirst10Reads += strlen(lineBuffer) + 2; //+1 is because of Windows CRLF.  fgets strips them.  If we're reading standard (\n only) text, then this will cause us to overallocate the lines array, which is fine.
    }

    fclose(inputForEstimate);

    HANDLE inputFileHandles[2];
    HANDLE fileMappingObjects[2];
    const char* fileMappings[2];

    for (int whichFile = 0; whichFile < nInputFiles; whichFile++) {
        inputFileHandles[whichFile] = CreateFileA(inputFileNames[whichFile], GENERIC_READ, FILE_SHARE_READ, NULL, OPEN_EXISTING, FILE_FLAG_SEQUENTIAL_SCAN, NULL);

        if (INVALID_HANDLE_VALUE == inputFileHandles[whichFile]) {
            fprintf(stderr, "Unable to open input file %s, %d\n", inputFileNames[whichFile], GetLastError());
            exit(1);
        }

        fileMappingObjects[whichFile] = CreateFileMapping(inputFileHandles[whichFile], NULL, PAGE_READONLY, 0, 0, NULL);
        if (NULL == fileMappingObjects[whichFile]) {
            fprintf(stderr, "Failed to create a file mapping object, %d\n", GetLastError());
            exit(1);
        }

        fileMappings[whichFile] = (const char *)MapViewOfFileEx(fileMappingObjects[whichFile], FILE_MAP_READ, 0, 0, 0, NULL);
        if (NULL == fileMappings[whichFile]) {
            fprintf(stderr, "Unable to MapViewOfFileEx, %d\n", GetLastError());
            exit(1);
        }
    }
    
    LARGE_INTEGER *fileSizes = new LARGE_INTEGER[nInputFiles];
    for (int whichFile = 0; whichFile < nInputFiles; whichFile++) {
        if (!GetFileSizeEx(inputFileHandles[0], &fileSizes[whichFile])) {
            fprintf(stderr, "Unable to GetFileSizeEx, %d\n", GetLastError());
            exit(1);
        }
    }

    _int64 approximateNReads = fileSizes[0].QuadPart * 10 / nBytesInFirst10Reads * 11 / 10; // The *11/10 makes this a 10% overestimate.

    const char*** readStarts = new const char** [nInputFiles];
    for (int whichFile = 0; whichFile < nInputFiles; whichFile++) {
        readStarts[whichFile] = new const char* [approximateNReads];
    }

    size_t* fileOffsets = new size_t[nInputFiles];
    _int64 nReads = 0;

    for (int whichFile = 0; whichFile < nInputFiles; whichFile++) {
        fileOffsets[whichFile] = 0;
    }

    printf("Scanning FASTQ file(s), one dot/million reads/file: "); fflush(stdout);

    while (fileOffsets[0] < fileSizes[0].QuadPart) {
        if (nReads >= approximateNReads) {
            fprintf(stderr, "Estimated too few reads and we're about to overflow the buffer.  FIXME!!\n");
            exit(1);
        }

        for (int whichFile = 0; whichFile < nInputFiles; whichFile++) {
            readStarts[whichFile][nReads] = fileMappings[whichFile] + fileOffsets[whichFile];

            for (int nLineBreaks = 0; nLineBreaks < 4; nLineBreaks++) { // The 4 lines of a FASTQ read
                //
                // Scan for end-of-line
                //
                while (fileMappings[whichFile][fileOffsets[whichFile]] != '\n') fileOffsets[whichFile]++;
                fileOffsets[whichFile]++;
                if (fileOffsets[whichFile] < fileSizes[whichFile].QuadPart && fileMappings[whichFile][fileOffsets[whichFile]] == '\r') fileOffsets[whichFile]++;
            }
        }

        nReads++;
        if (nReads % 1000000 == 0) {
            printf(".");
            fflush(stdout);
        }
    } // while we still have space

    printf("\n%lld reads\n", nReads);

    FILE** outputFiles = new FILE * [nInputFiles];
    for (int whichFile = 0; whichFile < nInputFiles; whichFile++) {
        outputFiles[whichFile] = fopen(outputFileNames[whichFile], "w");
        if (NULL == outputFiles[whichFile]) {
            fprintf(stderr, "Unable to open output file %s\n", outputFileNames[whichFile]);
            exit(1);
        }
        fileOffsets[whichFile] = 0;
    }

    if (nReads >= MAXUINT) {
        fprintf(stderr, "More than MAXUINT reads.  Rewrite the code to generate a 64 bit random number.\n");
        exit(1);
    }

    printf("One dot per percentile written: ");
    int lastPercentile = 0;
    for (int readToWrite = 0; readToWrite < __min(nReadsToWrite, nReads); readToWrite++) {
        unsigned randomNumber;
        int retVal;
        if (0 != (retVal = rand_s(&randomNumber))) {
            fprintf(stderr, "rand_s failed, %d\n", retVal);
            exit(1);
        }

        unsigned readIndex = randomNumber % (nReads - readToWrite);

        for (int whichFile = 0; whichFile < nInputFiles; whichFile++) {
            const char* readPointer = readStarts[whichFile][readIndex];
            for (int whichLine = 0; whichLine < 4; whichLine++) {
                while (*readPointer != '\n') {
                    putc(*readPointer, outputFiles[whichFile]);
                    readPointer++;
                }
                readPointer++;  // Skip the \n

                if (readPointer < fileMappings[whichFile] + fileSizes[whichFile].QuadPart && *readPointer == '\r') readPointer++;

                putc('\n', outputFiles[whichFile]);
            }

            readStarts[whichFile][readIndex] = readStarts[whichFile][nReads - readToWrite - 1]; // Overwrite the pointer with the last one, which we'll ignore after we loop.
        }

        if (lastPercentile != readToWrite * 100 / nReadsToWrite) {
            lastPercentile = readToWrite * 100 / nReadsToWrite;
            printf(".");
            fflush(stdout);
        }
    }

    printf("\n");
    fflush(stdout);

    for (int whichFile = 0; whichFile < nInputFiles; whichFile++) {
        fclose(outputFiles[whichFile]);
    }

}
