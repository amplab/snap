
#include <iostream>
#include <windows.h>



void usage()
{
    fprintf(stderr, "usage: BillSplit2 inputFilename outputFilenameBase nLinesPerFragment\n");
    fprintf(stderr, "output files will be named <outputFilenameBase>.n where n is a number with leading zeroes.\n");

    exit(1);
}


void WriteToOutputFile(const char* outputFilenameBase, int* pNextOutputFileNumber, HANDLE* phOutputFile, const char* begin, const char* end)
{
    if (begin == end) {
        return;
    }

    if (*phOutputFile == INVALID_HANDLE_VALUE) {
        size_t outputFilenameBufferSize = strlen(outputFilenameBase) + 1 + 8 + 1;  // .nnnnnnnn\0
        char* outputFilename = new char[outputFilenameBufferSize];

        snprintf(outputFilename, outputFilenameBufferSize, "%s.%08d", outputFilenameBase, *pNextOutputFileNumber);

        *phOutputFile = CreateFileA(outputFilename, GENERIC_WRITE, FILE_SHARE_READ, NULL, CREATE_ALWAYS, FILE_FLAG_SEQUENTIAL_SCAN, NULL);
        if (INVALID_HANDLE_VALUE == *phOutputFile || NULL == *phOutputFile) {
            fprintf(stderr, "Unable to open output file '%s', %d\n", outputFilename, GetLastError());
            exit(1);
        }

        delete[] outputFilename;

        (*pNextOutputFileNumber)++;
    }

    DWORD amountLeftToWrite = (DWORD)(end - begin);
    size_t amountWritten = 0;
    DWORD nBytesWritten;

    while (amountLeftToWrite > 0) {
        if (!WriteFile(*phOutputFile, begin + amountWritten, amountLeftToWrite, &nBytesWritten, NULL)) {
            fprintf(stderr, "WriteFile failed, %d\n", GetLastError());
            exit(1);
        }

        if (nBytesWritten < 1) {
            fprintf(stderr, "WriteFile succeeded but wrote no data.");
            exit(1);
        }

        amountLeftToWrite -= nBytesWritten;
        nBytesWritten += nBytesWritten;
    }
} // WriteOutputFile

int main(int argc, const char** argv)
{
    if (4 != argc) usage();

    const char* inputFilename = argv[1];
    const char* outputFilenameBase = argv[2];
    _int64 nLinesPerFragment = _atoi64(argv[3]);

    if (nLinesPerFragment < 1) usage();

    HANDLE hInputFile = CreateFileA(inputFilename, GENERIC_READ, FILE_SHARE_READ, NULL, OPEN_EXISTING, FILE_FLAG_SEQUENTIAL_SCAN, NULL);

    if (hInputFile == NULL || hInputFile == INVALID_HANDLE_VALUE) {
        fprintf(stderr, "Unable to open input file %s, %d.\n", inputFilename, GetLastError());
        exit(1);
    }


    const char* writeStartLine = NULL;
    const char* currentLine = NULL;
    const char* endAddress = NULL;
    HANDLE hOutputFile = INVALID_HANDLE_VALUE;
    _int64 nLinesInCurrentOutputFile = 0;
    int nextOutputFileNumber = 0;

    _int64 nBytesRead = 0;

    DWORD maxReadSize = 512 * 1024 * 1024;
    char* buffer = new char[maxReadSize];

    LARGE_INTEGER liFileSize;
    if (!GetFileSizeEx(hInputFile, &liFileSize)) {
        fprintf(stderr, "Unable to get input file size, %d\n", GetLastError());
        exit(1);
    }

    _int64 fileSize = liFileSize.QuadPart;

    while (nBytesRead < fileSize) {
        if (endAddress == NULL) {
            DWORD bytesReadThisRead;
            if (!ReadFile(hInputFile, buffer, maxReadSize, &bytesReadThisRead, NULL)) {
                fprintf(stderr, "ReadFile failed, %d\n", GetLastError());
                exit(1);
            }

            if (bytesReadThisRead == 0) {
                fprintf(stderr, "ReadFile completed successfully with no bytes read\n");
                exit(1);
            }

            endAddress = buffer + bytesReadThisRead;

            currentLine = writeStartLine = buffer;
        }

        while (currentLine < endAddress && *currentLine != '\n') {
            currentLine++;
            nBytesRead++;
        }

        if (currentLine < endAddress) {
            nLinesInCurrentOutputFile++;

            currentLine++;  // Skip the \n
            nBytesRead++;

            if (nLinesInCurrentOutputFile >= nLinesPerFragment) {
                WriteToOutputFile(outputFilenameBase, &nextOutputFileNumber, &hOutputFile, writeStartLine, currentLine);
                CloseHandle(hOutputFile);
                hOutputFile = INVALID_HANDLE_VALUE;

                nLinesInCurrentOutputFile = 0;
                writeStartLine = currentLine;
            }
        } else {
            //
            // Done with this input chunk.  Write it.
            //
            WriteToOutputFile(outputFilenameBase, &nextOutputFileNumber, &hOutputFile, writeStartLine, endAddress);

            currentLine = endAddress = writeStartLine = NULL;
        } 
    } // still have file to go

    if (nLinesInCurrentOutputFile != 0) {
        WriteToOutputFile(outputFilenameBase, &nextOutputFileNumber, &hOutputFile, writeStartLine, endAddress);
    }

    if (hOutputFile != INVALID_HANDLE_VALUE) {
        CloseHandle(hOutputFile);
    }
    return 0;
} // main