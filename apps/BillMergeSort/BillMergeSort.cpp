// BIllMergeSort.cpp : Defines the entry point for the console application.
//

#include "stdafx.h"

void soft_exit(int code)    // Just a place to set a breakpoint
{
    exit(code);
}


void usage()
{
    fprintf(stderr, "Usage: BillMergeSort inputFile outputFile\n");
    soft_exit(1);
}


struct LineData {
    char *lineStart;
    int lineLength;
};

class BillsBufferedInputFile {
public:
    BillsBufferedInputFile(HANDLE hFile_) {
        hFile = hFile_;
        atEOF = false;
        validBytesInBuffer = usedBytesInBuffer = 0;
    }

    inline char getNextChar() {
        if (atEOF) {
            fprintf(stderr, "Tried to getNextChar at EOF.\n");
            soft_exit(1);
        }

        if (usedBytesInBuffer >= validBytesInBuffer) {
            fprintf(stderr, "Called getNextChar() with an empty buffer.  You must call atEOF before each call to getNextChar()\n");
            soft_exit(1);
        }

        char value = inputBuffer[usedBytesInBuffer];
        usedBytesInBuffer++;
        return value;
    }

    inline bool isAtEOF() {
        if (atEOF) {
            return true;
        }
        if (usedBytesInBuffer >= validBytesInBuffer) {
            advanceBuffer();
        }
        return atEOF;
    }

private:
    static const int inputBufferSize = 16 * 1024 * 1024;
    char inputBuffer[inputBufferSize];

    bool atEOF;
    int validBytesInBuffer;
    int usedBytesInBuffer;

    void advanceBuffer() {
        DWORD nBytesRead;
        if (!ReadFile(hFile, inputBuffer, inputBufferSize, &nBytesRead, NULL)) {
            fprintf(stderr, "Error reading intermediate file, %d\n", GetLastError());
            soft_exit(1);
        }

        validBytesInBuffer = nBytesRead;
        if (0 == nBytesRead) {
            atEOF = true;
            CloseHandle(hFile);
        }

        usedBytesInBuffer = 0;
    }

    HANDLE hFile;
};


int compareLines(const void * line1Pointer, const void * line2Pointer)
{
    //
    // The parameters are pointers to LineData *s, and we want to compare the strings that they point to.
    //
    LineData *line1 = (LineData *)line1Pointer;
    LineData *line2 = (LineData *)line2Pointer;

    return strncmp(line1->lineStart, line2->lineStart, __min(line1->lineLength, line2->lineLength));    // NB: the lines contain the newline character, so substrings aren't a problem with strncmp
}

int longestLineSeen = 0;

bool ReadLineFromIntermediateFile(LineData *lineData, BillsBufferedInputFile *inputFile)
{
    lineData->lineLength = 0;

    while (!inputFile->isAtEOF() && lineData->lineLength < longestLineSeen) {
        lineData->lineStart[lineData->lineLength] = inputFile->getNextChar();
        if (lineData->lineStart[lineData->lineLength] == '\n') {
            lineData->lineStart[lineData->lineLength + 1] = '\0';   // Remember, there's an extra character for the null terminator in the buffer.
            lineData->lineLength++;
            return true;
        }
        lineData->lineLength++;
    }

    if (inputFile->isAtEOF()) {
        if (lineData->lineLength == 0) {
            return false;
        }

        lineData->lineStart[lineData->lineLength] = '\0';
        lineData->lineLength++;
        return true;
    }

    fprintf(stderr, "Found a line that got longer since it was written.  Bizarre, or code bug.  lineData->lineLength = %d\n", lineData->lineLength);
    soft_exit(1);
    return false;
}

SYSTEM_INFO systemInfo;
struct WorkerThreadState {
    LineData *lineBuffers;
    size_t nLineBuffers;
    size_t lineBuffersMerged;   // This is used in the merge phase, so it's not really thread state, but it's convenient to put it here.
    HANDLE hDoneEvent;
    volatile long *nThreadsStillRunning;
};

DWORD WINAPI SortWorkerThreadMain(void *pvWorkerThreadState)
{
    WorkerThreadState *state = (WorkerThreadState *)pvWorkerThreadState;

    qsort(state->lineBuffers, state->nLineBuffers, sizeof(*state->lineBuffers), compareLines);

    if (0 == InterlockedDecrement(state->nThreadsStillRunning)) {
        if (!SetEvent(state->hDoneEvent)) {
            fprintf(stderr, "Unable to SetEvent, %d\n", GetLastError());
            soft_exit(1);
        }
        // state is now invalid
    }

    return 0;
}

void parallelInMemoryMergeSort(LineData *lineBuffers, size_t nLineBuffers)
{
    if (nLineBuffers < 2) {
        return;
    }

    //
    // Fire off threads to sort chunks of the lines.  We want at most one thread per processor, and also at least 1 GB of data/thread to avoid much false sharing near the boundaries.
    //

    volatile long nThreadsRunning = 0;
    long nThreads = 0;
    size_t usedLineBuffers = 0;

    WorkerThreadState *threadState = new WorkerThreadState[systemInfo.dwNumberOfProcessors];
    HANDLE hEvent = CreateEvent(NULL, TRUE, FALSE, NULL);

    if (INVALID_HANDLE_VALUE == hEvent) {
        fprintf(stderr, "Unable to create event, %d\n", GetLastError());
        soft_exit(1);
    }

    while (usedLineBuffers < nLineBuffers) {
        if ((DWORD)nThreads >= systemInfo.dwNumberOfProcessors) {
            fprintf(stderr, "Code bug: tried to start too many threads.\n");
            soft_exit(1);
        }

        size_t endForNextThread = __min(usedLineBuffers + nLineBuffers / systemInfo.dwNumberOfProcessors + 1, nLineBuffers);

        size_t minSize = 1024 * 1024 * 1024;

        if ((size_t)(lineBuffers[endForNextThread].lineStart - lineBuffers[usedLineBuffers].lineStart) < minSize && endForNextThread < nLineBuffers - 1) {
            //
            // Just guesstimate.  We don't need to be really accurate.
            //
            size_t dataInRange = lineBuffers[endForNextThread].lineStart - lineBuffers[usedLineBuffers].lineStart;
            endForNextThread = __min(nLineBuffers, 
                endForNextThread + ((minSize - dataInRange) / (dataInRange / (endForNextThread - usedLineBuffers))));
        }

        threadState[nThreads].lineBuffers = lineBuffers + usedLineBuffers;
        threadState[nThreads].nLineBuffers = endForNextThread - usedLineBuffers;
        threadState[nThreads].hDoneEvent = hEvent;
        threadState[nThreads].nThreadsStillRunning = &nThreadsRunning;

        usedLineBuffers = endForNextThread;

        nThreads++;
    }

    nThreadsRunning = nThreads;

    for (long i = 0; i < nThreads; i++) {
        DWORD threadId;
        HANDLE hThread = CreateThread(NULL, 0, SortWorkerThreadMain, threadState + i, 0, &threadId);
        if (INVALID_HANDLE_VALUE == hThread) {
            fprintf(stderr, "Unable to create thread, %d\n", GetLastError());
            soft_exit(1);
        }
    }

    if (WAIT_OBJECT_0 != WaitForSingleObject(hEvent, INFINITE)) {
        fprintf(stderr, "WaitForSingleObject failed, %d\n", GetLastError());
        soft_exit(1);
    }

    //
    // Now merge.
    //
    LineData *mergedLineData = new LineData[nLineBuffers];

    int nThreadsWithDataRemaining = 0;
    for (int i = 0; i < nThreads; i++) {
        if (threadState[i].nLineBuffers > 0) {
            threadState[i].lineBuffersMerged = 0;
            nThreadsWithDataRemaining++;
        }
    }

    size_t nMergedLines = 0;

    while (nThreadsWithDataRemaining > 0) {
        int threadWithNextLine = -1;

        for (int i = 0; i < nThreads; i++) {
            if (threadState[i].lineBuffersMerged < threadState[i].nLineBuffers &&
                (-1 == threadWithNextLine || compareLines(threadState[i].lineBuffers + threadState[i].lineBuffersMerged, threadState[threadWithNextLine].lineBuffers + threadState[threadWithNextLine].lineBuffersMerged) < 0)) {
                threadWithNextLine = i;
            }
        }

        if (-1 == threadWithNextLine) {
            fprintf(stderr, "Code bug: didn't find merge point.\n");
            soft_exit(1);
        }

        WorkerThreadState *thisState = &threadState[threadWithNextLine];

        mergedLineData[nMergedLines] = thisState->lineBuffers[thisState->lineBuffersMerged];
        nMergedLines++;

        thisState->lineBuffersMerged++;
        if (thisState->lineBuffersMerged >= thisState->nLineBuffers) {
            nThreadsWithDataRemaining--;
        }
    }

    if (nMergedLines != nLineBuffers) {
        fprintf(stderr, "Code bug: didn't merge the right number of lines.\n");
        soft_exit(1);
    }
    
    memcpy(lineBuffers, mergedLineData, sizeof(*lineBuffers) * nLineBuffers);
    delete[] mergedLineData;
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

int main(int argc, char* argv[])
{
    _int64 startTickCount = GetTickCount64();

    if (argc != 3) usage();

    GetSystemInfo(&systemInfo); // This is the rare Win32 API that literally can't fail.

    MEMORYSTATUSEX memoryStatus;
    memoryStatus.dwLength = sizeof(memoryStatus);

    if (!GlobalMemoryStatusEx(&memoryStatus)) {
        fprintf(stderr, "GlobalMemoryStatusEx failed, %d\n", GetLastError());
        soft_exit(1);
    }

    _int64 memoryToUse = (memoryStatus.ullAvailPhys * 7) / 10;

    if (memoryToUse < 1024 * 1024 * 1024) {
        fprintf(stderr, "GlobalMemoryStatusEx reports less than a gigabyte of free memory, quitting. (%lld)\n", memoryToUse);
        soft_exit(1);
    }

    HANDLE hInputFile = CreateFile(argv[1], GENERIC_READ, FILE_SHARE_READ, NULL, OPEN_EXISTING, 0, NULL);

    if (INVALID_HANDLE_VALUE == hInputFile) {
        fprintf(stderr, "Unable to open file '%s' for read, %d\n", argv[1], GetLastError());
        soft_exit(1);
    }

    FILE_STANDARD_INFO fileInfo;

    if (!GetFileInformationByHandleEx(hInputFile, FileStandardInfo, &fileInfo, sizeof(fileInfo))) {
        fprintf(stderr, "Unable to get file information, %d\n", GetLastError());
        soft_exit(1);
    }

    if (0 == fileInfo.EndOfFile.QuadPart) {
        HANDLE hOutputFile = CreateFile(argv[2], GENERIC_WRITE, FILE_SHARE_READ, NULL, CREATE_ALWAYS, FILE_ATTRIBUTE_NORMAL, NULL);
        if (INVALID_HANDLE_VALUE == hOutputFile) {
            fprintf(stderr, "Error opening output file for empty input, %d\n", GetLastError());
            soft_exit(1);
        }
        CloseHandle(hOutputFile);
        soft_exit(0);
    }

    HANDLE hMapping = CreateFileMapping(hInputFile, NULL, PAGE_READONLY, 0, 0, NULL);

    if (NULL == hMapping) {
        fprintf(stderr, "Unable to CreateFileMapping, %d\n", GetLastError());
        soft_exit(1);
    }

    //
    // Figure out how many chunks we'll need.
    //
    _int64 nChunks = (fileInfo.EndOfFile.QuadPart + memoryToUse - 1) / memoryToUse;

    if (0 == nChunks) {
        fprintf(stderr, "Computed 0 chunks.  Since we special-cased the empty input file, this is a code bug.  FileSize = %lld, memoryToUse = %lld\n", fileInfo.EndOfFile.QuadPart, memoryToUse);
        soft_exit(1);
    }
    
    _int64 chunkSize = (fileInfo.EndOfFile.QuadPart + nChunks - 1) / nChunks;

    _int64 lineBufferSize = memoryToUse / 80;  // Just a stab at it.
    LineData *lines = new LineData[lineBufferSize];

    LARGE_INTEGER amountOfFileProcessed;
    amountOfFileProcessed.QuadPart = 0;

    _int64 totalLines = 0;

    printf("Processing file in %lld chunks of size about %lldGB\n", nChunks, (fileInfo.EndOfFile.QuadPart / nChunks + 512 * 1024 * 1024) / (1024 * 1024 * 1024));

    size_t intermediateFileNameBufferSize = strlen(argv[2]) + 20;
    char *intermediateFileName = new char[intermediateFileNameBufferSize];

    _int64 tickCount;

    for (_int64 whichChunk = 0; whichChunk < nChunks; whichChunk++) {
        int nLinesInChunk = 0;

        //
        // Map a chunk plus a something much bigger than a line, or to EOF, whichever is smaller.
        //

        //
        // We need to map in PAGE_SIZE offsets, so round down the mapping offset.
        //
        LARGE_INTEGER mappingOffset;
        mappingOffset.QuadPart = amountOfFileProcessed.QuadPart / systemInfo.dwAllocationGranularity * systemInfo.dwAllocationGranularity;

        size_t mappedSize = __min(chunkSize + 10 * 1024 * 1024, fileInfo.EndOfFile.QuadPart - mappingOffset.QuadPart);


        char *mappedBase = (char *)MapViewOfFile(hMapping, FILE_MAP_READ, mappingOffset.HighPart, mappingOffset.LowPart, mappedSize);
        if (NULL == mappedBase) {
            fprintf(stderr, "Unable to map chunk, %d\n", GetLastError());
            soft_exit(1);
        }
        char *chunkBase = mappedBase + (amountOfFileProcessed.QuadPart - mappingOffset.QuadPart);

        char *mappingEnd = mappedBase + mappedSize;

        char *nextLineToProcess = chunkBase;

        tickCount = GetTickCount64();
        printf("Finding line breaks for chunk %lld...", whichChunk); fflush(stdout);

        while (nextLineToProcess < mappingEnd && nextLineToProcess - chunkBase < chunkSize) {
            if (nLinesInChunk == lineBufferSize) {
                //
                // Expand the line start buffer by 30%.
                //

                _int64 newLineBufferSize = (lineBufferSize * 13) / 10;
                LineData *newLineBuffer = new LineData[newLineBufferSize];

                memcpy(newLineBuffer, lines, sizeof(*newLineBuffer)* nLinesInChunk);
                delete[] lines;
                lines = newLineBuffer;

                lineBufferSize = newLineBufferSize;
            }

            lines[nLinesInChunk].lineStart = nextLineToProcess;

            while (nextLineToProcess < mappingEnd && *nextLineToProcess != '\n') {
                nextLineToProcess++;
            }

            if (nextLineToProcess < mappingEnd) nextLineToProcess++;  // Skip over the newline

            if (nextLineToProcess - lines[nLinesInChunk].lineStart > MAXINT) {
                fprintf(stderr, "Input file contains too big of a line (>= 2GB)\n");
                soft_exit(1);
            }

            lines[nLinesInChunk].lineLength = (int)(nextLineToProcess - lines[nLinesInChunk].lineStart);
            longestLineSeen = __max(longestLineSeen, lines[nLinesInChunk].lineLength);
            nLinesInChunk++;
            totalLines++;
        }

        printf("%llds\nSorting chunk %lld (%ld lines)...", (GetTickCount64() - tickCount + 500) / 1000, whichChunk, nLinesInChunk); fflush(stdout);

        tickCount = GetTickCount64();

        //qsort(lines, nLinesInChunk, sizeof(*lines), compareLines);
        parallelInMemoryMergeSort(lines, nLinesInChunk);

        printf("%llds\nWriting out chunk %lld...", (GetTickCount64() - tickCount + 500) / 1000, whichChunk); fflush(stdout);
        tickCount = GetTickCount64();

        BillWriteFile intermediateFile;

        sprintf_s(intermediateFileName, intermediateFileNameBufferSize, "%s.%lld", argv[2], whichChunk);

        if (!intermediateFile.open(intermediateFileName)) {
            fprintf(stderr, "Unable to open intermediate file '%s', %d\n", intermediateFileName, GetLastError());
            soft_exit(1);
        }

        //
        // Pre-allocate the intermediate file to avoid fragmentation.
        //
        LARGE_INTEGER intermediateFileSize;
        intermediateFileSize.QuadPart = nextLineToProcess - chunkBase;
        if (!intermediateFile.preallocate(intermediateFileSize)) {
            fprintf(stderr, "Unable to preallocate intermediate file, %d.\n", GetLastError());
            soft_exit(1);
        }

        for (int whichLine = 0; whichLine < nLinesInChunk; whichLine++) {
            if (!intermediateFile.write(lines[whichLine].lineStart, lines[whichLine].lineLength)) {
                fprintf(stderr, "Error writing intermediate file '%s', %d\n", intermediateFileName, GetLastError());
                soft_exit(1);
            }
        }

        amountOfFileProcessed.QuadPart += nextLineToProcess - chunkBase;

        if (!intermediateFile.close()) {
            fprintf(stderr, "Unable to close intermediate file, %d\n", GetLastError());
            soft_exit(1);
        }

        if (!UnmapViewOfFile(mappedBase)) {
            fprintf(stderr, "Unable to UnmapViewOfFile, %d\n", GetLastError());
            soft_exit(1);
        }

        printf("%llds\n", (GetTickCount64() - tickCount + 500) / 1000);
    } // for each chunk.

    if (amountOfFileProcessed.QuadPart != fileInfo.EndOfFile.QuadPart) {
        fprintf(stderr, "Didn't reach EOF.  Probably a program bug.  %lld != %lld\n", amountOfFileProcessed.QuadPart, fileInfo.EndOfFile.QuadPart);
        soft_exit(1);
    }
    
    CloseHandle(hMapping);
    CloseHandle(hInputFile);

    if (nChunks < 2) {
        sprintf_s(intermediateFileName, intermediateFileNameBufferSize, "%s.0", argv[2]);
        if (!MoveFileEx(intermediateFileName, argv[2], MOVEFILE_REPLACE_EXISTING)) {
            fprintf(stderr, "Unable to move single intermediate file to output file, %d\n", GetLastError());
            soft_exit(1);
        }
        soft_exit(0);
    }

    //
    // Now do the merge.
    //
    printf("Merging into output file..."); fflush(stdout);
    tickCount = GetTickCount64();

    BillsBufferedInputFile **intermediateFiles = new BillsBufferedInputFile *[nChunks];
    for (_int64 whichChunk = 0; whichChunk < nChunks; whichChunk++) {
        sprintf_s(intermediateFileName, intermediateFileNameBufferSize, "%s.%lld", argv[2], whichChunk);
        HANDLE hIntermediateFile = CreateFile(intermediateFileName, GENERIC_READ, FILE_SHARE_READ, NULL, OPEN_EXISTING, FILE_FLAG_SEQUENTIAL_SCAN, NULL);
        if (INVALID_HANDLE_VALUE == hIntermediateFile) {
            fprintf(stderr, "Unable to reopen intermediate file '%s', %d\n", intermediateFileName, GetLastError());
            soft_exit(1);
        }

        intermediateFiles[whichChunk] = new BillsBufferedInputFile(hIntermediateFile);
    }

    LineData *intermediateFileNextLines = new LineData[nChunks];
    bool *intermediateFileFinished = new bool[nChunks];
    _int64 nIntermediateFilesRemaining = nChunks;
    for (_int64 whichChunk = 0; whichChunk < nChunks; whichChunk++) {
        intermediateFileNextLines[whichChunk].lineStart = new char[longestLineSeen + 1];
        intermediateFileFinished[whichChunk] = false;

        if (!ReadLineFromIntermediateFile(&intermediateFileNextLines[whichChunk], intermediateFiles[whichChunk])) {
            intermediateFileFinished[whichChunk] = true;
            nIntermediateFilesRemaining--;
        }
    }

    BillWriteFile outputFile;
    if (!outputFile.open(argv[2])) {
        fprintf(stderr, "Unable to open output file '%s', %d\n", argv[2], GetLastError());
        soft_exit(1);
    }

    if (!outputFile.preallocate(fileInfo.EndOfFile)) {
        fprintf(stderr, "Unable to preallocate output file, %d.\n", GetLastError());
        soft_exit(1);
    }

    while (nIntermediateFilesRemaining > 0) {
        int nextLineToUse = -1;

        for (int whichChunk = 0; whichChunk < nChunks; whichChunk++) {
            if (!intermediateFileFinished[whichChunk]) {
                if (-1 == nextLineToUse) {
                    nextLineToUse = whichChunk;
                }
                else if (compareLines(&intermediateFileNextLines[whichChunk], &intermediateFileNextLines[nextLineToUse]) < 0) {
                    nextLineToUse = whichChunk;
                }
            }
        }

        if (-1 == nextLineToUse) {
            fprintf(stderr, "Code bug: didn't find a line even though we should have one.\n");
            soft_exit(1);
        }

        if (!outputFile.write(intermediateFileNextLines[nextLineToUse].lineStart, intermediateFileNextLines[nextLineToUse].lineLength)) {
            fprintf(stderr, "Error writing output file, %d\n", GetLastError());
            soft_exit(1);
        }

        if (!ReadLineFromIntermediateFile(&intermediateFileNextLines[nextLineToUse], intermediateFiles[nextLineToUse])) {
            intermediateFileFinished[nextLineToUse] = true;
            nIntermediateFilesRemaining--;
        }
    }

    if (!outputFile.close()) {
        fprintf(stderr, "Error closing output file, %d\n", GetLastError());
        soft_exit(1);
    }

    for (_int64 whichChunk = 0; whichChunk < nChunks; whichChunk++) {
        sprintf_s(intermediateFileName, intermediateFileNameBufferSize, "%s.%lld", argv[2], whichChunk);
        if (!DeleteFile(intermediateFileName)) {
            fprintf(stderr, "Unable to delete intermediate file '%s', %d\n", intermediateFileName, GetLastError());
            //
            // Just keep trying the rest anyway.
            //
        }
    }

    printf("%llds\nProcessed a total of %lld lines and %lld bytes in %llds, %lld MB/s\n", (GetTickCount64() - tickCount + 500) / 1000, totalLines, fileInfo.EndOfFile.QuadPart, (GetTickCount64() - startTickCount + 500) / 1000,
        fileInfo.EndOfFile.QuadPart / (GetTickCount64() - tickCount) / 1024);
}

