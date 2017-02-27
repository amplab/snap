// BillWriter.cpp : Defines the entry point for the console application.
//

#include "stdafx.h"

void soft_exit(int code)    // Just a place to set a breakpoint
{
    exit(code);
}

class BillWriteFile
{
public:

    BillWriteFile()
    {
        buffer = new char[bufferSize];
        bufferUsed = 0;
        hFile = INVALID_HANDLE_VALUE;
    }

    bool open(char *filename)
    {
        if (INVALID_HANDLE_VALUE != hFile) {
            fprintf(stderr, "Can only open an BillWriteFile that's not already open.\n");
            soft_exit(1);
        }
        hFile = CreateFile(filename, GENERIC_READ | GENERIC_WRITE, FILE_SHARE_READ, NULL, CREATE_ALWAYS, FILE_ATTRIBUTE_NORMAL | FILE_FLAG_SEQUENTIAL_SCAN, NULL);
        if (INVALID_HANDLE_VALUE == hFile) {
            return false;
        }

        return true;
    }

    bool preallocate(LARGE_INTEGER preallocationSize)
    {
        LARGE_INTEGER newFilePointer;
        LARGE_INTEGER zero;
        zero.QuadPart = 0;

        if (!SetFilePointerEx(hFile, preallocationSize, &newFilePointer, FILE_BEGIN) || !SetEndOfFile(hFile) || !SetFilePointerEx(hFile, zero, &newFilePointer, FILE_BEGIN))
        {
            return false;
        }

        return true;
    }

    bool write(char *data, DWORD length)
    {
        if (length + bufferUsed > bufferSize && bufferUsed > 0) {
            if (!flushBuffer()) {
                return false;
            }
        }

        if (length >= bufferSize) {
            //
            // Don't buffer this giant write, just write it directly.
            //
            DWORD nBytesWritten;
            if (!WriteFile(hFile, data, length, &nBytesWritten, NULL)) {
                return false;
            }

            if (nBytesWritten != length) {
                fprintf(stderr, "Didn't write the expected number of bytes (2), %d != %d\n", length, nBytesWritten);
                soft_exit(1);
            }
            return true;
        }

        memcpy(buffer + bufferUsed, data, length);
        bufferUsed += length;

        return true;
    }

    bool close()
    {
        if (INVALID_HANDLE_VALUE == hFile) {
            fprintf(stderr, "Can't close unopened file\n");
            soft_exit(1);
        }

        if (!flushBuffer()) {
            return false;
        }

        CloseHandle(hFile);
        hFile = INVALID_HANDLE_VALUE;

        return true;
    }

    ~BillWriteFile()
    {
        if (INVALID_HANDLE_VALUE != hFile) {
            fprintf(stderr, "Must close BillWriteFile before destroying the object.\n");
            soft_exit(1);
        }

        delete[] buffer;
    }

private:

    bool flushBuffer()
    {
        if (bufferUsed == 0) {
            return true;
        }

        DWORD nBytesWritten;

        if (!WriteFile(hFile, buffer, bufferUsed, &nBytesWritten, NULL)) {
            return false;
        }

        if (nBytesWritten != bufferUsed) {
            fprintf(stderr, "Didn't write the expected number of bytes, %d != %d\n", bufferUsed, nBytesWritten);
            soft_exit(1);
        }

        bufferUsed = 0;

        return true;
    }

    HANDLE hFile;
    char *buffer;
    DWORD bufferUsed;


    static const size_t bufferSize = 64 * 1024 * 1024;
};

void usage()
{
    fprintf(stderr, "Usage: BillWriter outputFile\n");
    fprintf(stderr, "Reads stdin and writes it in big chunks to outputFile.\n");
    soft_exit(1);
}

int main(int argc, char* argv[])
{
    if (2 != argc) usage();

    HANDLE hInput = GetStdHandle(STD_INPUT_HANDLE);

    if (NULL == hInput) {
        fprintf(stderr, "Couldn't get stdin handle, %d\n", GetLastError());
        soft_exit(1);
    }

    BillWriteFile outputFile;

    if (!outputFile.open(argv[1])) {
        fprintf(stderr, "Unable to open output file '%s', %d\n", argv[1], GetLastError());
        soft_exit(1);
    }

    const int inputBufferSize = 1024 * 1024;
    char *inputBuffer = new char[inputBufferSize];

    for (;;) {
        DWORD nBytesRead;
        if (!ReadFile(hInput, inputBuffer, inputBufferSize, &nBytesRead, NULL)) {
            if (GetLastError() == ERROR_BROKEN_PIPE) {
                nBytesRead = 0; // Just EOF
            } else {
                fprintf(stderr, "Unable to read input, %d\n", GetLastError());
                soft_exit(1);
            }
        }

        if (0 == nBytesRead) {
            outputFile.close();
            soft_exit(0);
        }

        if (!outputFile.write(inputBuffer, nBytesRead)) {
            fprintf(stderr, "Error writing output file, %d\n", GetLastError());
            soft_exit(1);
        }
    }

    //NOTREACHED

    return 0;
}

