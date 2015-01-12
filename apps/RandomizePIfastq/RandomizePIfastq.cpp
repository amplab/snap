/*++

Module Name:

    RandomizePIfastq.cpp

Abstract:

   Cheezy little app that takes a paired, interleaved FASTQ and more-or-less randomizes it.

Authors:

    Bill Bolosky, September, 2013

Revision History:

    
--*/

#include "stdafx.h"
#include "GoodRandom.h"

class IntermediateFile {
public:
    IntermediateFile(const char *fileName, size_t i_bufferSize);
    ~IntermediateFile();

 
    void addLines(char *lines, size_t size);
    void close();

private:
       void flush();


    HANDLE       hFile;
    size_t       bufferSize;
    size_t       bufferUsed;
    char        *buffer;
};

IntermediateFile::IntermediateFile(const char *fileName, size_t i_bufferSize) : bufferSize(i_bufferSize)
{
    hFile = CreateFile(fileName, GENERIC_READ|GENERIC_WRITE, FILE_SHARE_READ, NULL, CREATE_ALWAYS, FILE_ATTRIBUTE_NORMAL, NULL);
    if (INVALID_HANDLE_VALUE == hFile) {
        fprintf(stderr,"Unable to open intermediate file '%s', %d\n", fileName, GetLastError());
        return;
    }
    buffer = new char[bufferSize];
    bufferUsed = 0;
}

IntermediateFile::~IntermediateFile()
{
    close();
}

    void
IntermediateFile::addLines(char *lines, size_t size)
{
    if (size + bufferUsed > bufferSize) {
        flush();
    }

    memcpy(buffer + bufferUsed, lines, size);

    bufferUsed += size;
}

    void
IntermediateFile::close()
{
    flush();
    CloseHandle(hFile);
    hFile = NULL;
}

    void
IntermediateFile::flush()
{
    DWORD bytesWritten;
    if (!WriteFile(hFile, buffer, bufferUsed, &bytesWritten, NULL)) {
        fprintf(stderr,"WriteFile failed, %d\n", GetLastError());
        return;
    }
    bufferUsed = 0;
}

    void
usage()
{
    fprintf(stderr,"usage: RandomizePIfastq inputFile outputFile\n");
}

const int inputBufferSize = 100 * 1024 * 1024;
char *inputBuffer[2];

    void
main(int argc, const char **argv)
{
    for (int i = 0; i < 2; i++) {
        inputBuffer[i] = new char[inputBufferSize];
    }
    const int nFiles = 1000;

    if (3 != argc) {
        usage();
        return;
    }

    HANDLE hInputFile = CreateFile(argv[1],GENERIC_READ,FILE_SHARE_READ,NULL,OPEN_EXISTING, FILE_ATTRIBUTE_NORMAL,NULL);
    if (INVALID_HANDLE_VALUE == hInputFile) {
        fprintf(stderr,"Unable to open '%s' for input, error %d\n", argv[1], GetLastError());
        return;
    }

    IntermediateFile **files = new IntermediateFile *[nFiles];
    for (int i = 0; i < nFiles; i++) {
        char fileName[100];
        sprintf(fileName,"piFQ.piece.%04d", i);
        files[i] = new IntermediateFile(fileName, 10 * 1024 * 1024);
    }


    const int nLinesPerItem = 8;

    int whichInputBuffer = 0;
    DWORD usedBytes = 0;
    DWORD validBytes = 0;
    bool done = false;
    for (;;) {
        DWORD leftoverBytes = validBytes - usedBytes;   // i.e., how much was in the old buffer minus what we used.
        memcpy(inputBuffer[whichInputBuffer], inputBuffer[1 - whichInputBuffer] + usedBytes, leftoverBytes);
        DWORD bytesRead;
        if (!ReadFile(hInputFile, inputBuffer[whichInputBuffer] + leftoverBytes, inputBufferSize - leftoverBytes, &bytesRead, NULL)) {
            fprintf(stderr,"Read error on input file %d\n", GetLastError());
            return;
        }

        if (0 == bytesRead) {
            if (0 != leftoverBytes) {
                fprintf(stderr,"Input file doesn't seem to have a multiple of %d lines\n", nLinesPerItem);
            }
            break;
        }

        validBytes = leftoverBytes + bytesRead; // i.e., what we copied from the last one, plus what's new
        usedBytes = 0;

        for (;;) {
 
            char *nextItem = inputBuffer[whichInputBuffer] + usedBytes;

            int nNewLines = 0;
            while (nNewLines < nLinesPerItem && usedBytes < validBytes) {
                if (inputBuffer[whichInputBuffer][usedBytes] == '\n') {
                    nNewLines++;
                }
                usedBytes++;    // NB: deliberately eating the \n
            }

            if (nNewLines < nLinesPerItem) {
                usedBytes = nextItem - inputBuffer[whichInputBuffer];
                break;
            }

            _uint64 whichFile = GoodFastRandom(nFiles - 1);
            files[whichFile]->addLines(nextItem, usedBytes - (nextItem - inputBuffer[whichInputBuffer]));
        }

        whichInputBuffer = 1 - whichInputBuffer;
    }

    for (int i = 0; i < nFiles; i++) {
        delete files[i];    // This also flushes them
    }
}