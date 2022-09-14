//
// FixNovoHg1.cpp : For reasons I don't understand Novoalign's hg1 output has
// some extra newlines in it.  This gets rid of them
//

#include <iostream>
#include <Windows.h>

int main(int argc, const char **argv)
{
    _int64 startTime = GetTickCount64();

    if (3 != argc) {
        fprintf(stderr, "usage: FixNovoHg1 input.sam output.sam\n");
        exit(1);
    }

    HANDLE hInputFile = CreateFile(argv[1], GENERIC_READ, FILE_SHARE_READ, NULL, OPEN_EXISTING, FILE_FLAG_SEQUENTIAL_SCAN, NULL);
    if (INVALID_HANDLE_VALUE == hInputFile) {
        fprintf(stderr, "Unable to to open '%s' for read, %d\n", argv[1], GetLastError());
        exit(1);
    }

    HANDLE hOutputFile = CreateFile(argv[2], GENERIC_WRITE, FILE_SHARE_READ, NULL, CREATE_NEW, FILE_FLAG_SEQUENTIAL_SCAN, NULL);
    if (INVALID_HANDLE_VALUE == hOutputFile) {
        fprintf(stderr, "Unable to create '%s' for write.  The file must not already exist.  Error %d\n", argv[2], GetLastError());
        exit(1);
    }

    const size_t bufferSize = 1024 * 1024 * 1024;
    char* buffer = new char[bufferSize];
    bool previousByteWasCR = false;

    _int64 nExcised = 0;
    _int64 totalBytesRead = 0;

    DWORD bytesInBuffer;
    for (;;) {
        if (!ReadFile(hInputFile, buffer, bufferSize, &bytesInBuffer, NULL)) {
            fprintf(stderr, "ReadFile failed, %d\n", GetLastError());
            exit(1);
        }

        if (bytesInBuffer == 0) {
            //
            // EOF.
            //
            break;
        }

        totalBytesRead += bytesInBuffer;

        DWORD startOffset = 0;
        for (DWORD currentOffset = 0; currentOffset < bytesInBuffer; currentOffset++) {
            if (buffer[currentOffset] == '\r') {
                previousByteWasCR = true;
            } else if (buffer[currentOffset] == '\t' && previousByteWasCR) {
                //
                // Write the buffer from the last write position (startOffset) to just before
                // the \r (i.e., two bytes less than where we are).
                //
                if (currentOffset > 1) {
                    DWORD bytesToWrite = currentOffset - 2 - startOffset;
                    DWORD bytesWritten;
                    if (!WriteFile(hOutputFile, buffer + startOffset, currentOffset - 2 - startOffset, &bytesWritten, NULL)) {
                        fprintf(stderr, "WriteFile failed, %d\n", GetLastError());
                        exit(1);
                    }

                    if (bytesWritten != bytesToWrite) {
                        fprintf(stderr, "WriteFile didn't write the whole buffer, and I didn't implement support for partial writes.\n");
                        exit(1);
                    }
                }
                startOffset = currentOffset - 1;
                previousByteWasCR = false;

                nExcised++;
            } else {
                previousByteWasCR = false;
            }
        }

        //
        // Write the last chunk.  If the last byte was CR, don't write it.
        //
        DWORD bytesToWrite = bytesInBuffer - startOffset;
        if (bytesToWrite > 0 && previousByteWasCR) {
            bytesToWrite--;
        }

        if (bytesToWrite > 0) {
            DWORD bytesWritten;
            if (!WriteFile(hOutputFile, buffer + startOffset, bytesToWrite, &bytesWritten, NULL)) {
                fprintf(stderr, "WriteFile failed, %d\n", GetLastError());
                exit(1);
            }

            if (bytesWritten != bytesToWrite) {
                fprintf(stderr, "WriteFile only wrote %d bytes when we requesed %d\n", bytesToWrite, bytesWritten);
                exit(1);
            }
        }

    } // for ever

    CloseHandle(hInputFile);
    CloseHandle(hOutputFile);

    delete[] buffer;

    printf("Processed %lld bytes in %llds and excised %lld spurious carriage returns\n", totalBytesRead, (GetTickCount64() - startTime) / 1000, nExcised);
    
}

