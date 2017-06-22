// GenerateConsolodatedExtractedReads.cpp : Defines the entry point for the console application.
//

#include "stdafx.h"
#include "Compat.h"

//
// String these together to hold the output from a run of samtools, or a set of them.
//

struct ReadsBuffer {
    char *buffer;
    size_t size;
    ReadsBuffer *next;
};

const size_t filenameBufferSize = 500;
struct SamtoolsRun {
    size_t size;
    _int64 fileOffset;  // Only gets filled in once it's written.
    char   *filename;
    ReadsBuffer *readsBuffers;

    SamtoolsRun *next;  // For chaining these together
};

SamtoolsRun *runsReadyToWrite = NULL;


CRITICAL_SECTION criticalSection[1];
size_t dataSizeOnWriteQueue = 0;

const size_t MinWriteSize = 128 * 1024 * 1024;
const size_t MaxWriteQueueSize = 512 * 1024 * 1024;

HANDLE runsReadyToWriteQueueHasEnoughDataEvent;
HANDLE runsReadyToWriteNotTooFullEvent;
HANDLE writerDoneEvent;

volatile long nextSamtoolsRunToProcess = 0;  // Which is to say index into inputLines
long nInputLines;
char **inputLines;

volatile long nThreadsRunning = 16;  // Default parallelism

char *samtoolsPathname = "c:\\bolosky\\bin\\samtools.exe";

HANDLE hOutputFile;
FILE *indexFile;

void usage()
{
    fprintf(stderr, "usage: GenerateConsolodatedExtractedReads InputScriptFile OutputFilenameBase {nThreads {samtoolsPathname}}\n");
    fprintf(stderr, "Take a script of samtools runs and executes it, combining the results into a combined file with an index.\n");
    fprintf(stderr, "This is mostly to avoid creating tens of millions of individual files and making the filesystem unhappy.\n");
    fprintf(stderr, "Specifying 0 threads means to use one per core on the machine (which is the default anyway).\n");

    exit(1);
}

int __cdecl
CompareSamtoolsRuns(const void *ppRun1, const void *ppRun2)
{
    const SamtoolsRun *run1 = *(SamtoolsRun **)ppRun1;
    const SamtoolsRun *run2 = *(SamtoolsRun **)ppRun2;

    return strcmp(run1->filename, run2->filename);
}

DWORD __stdcall
 WriterThreadMain(void *param)
{
     SamtoolsRun *runsWritten = NULL;
     int nRunsWritten = 0;
     size_t nBytesWritten = 0;
     size_t offsetInOutputFile = 0;

     _int64 startTime = GetTickCount64();

     for (;;) {
        if (0 != nThreadsRunning) {
            if (WAIT_OBJECT_0 != WaitForSingleObject(runsReadyToWriteQueueHasEnoughDataEvent, INFINITE)) {
                fprintf(stderr, "WriterThread: WaitForSingleObject failed, %d\n", GetLastError());
                exit(1);
            }
        }

        //
        // Enter the critical section and grab what work's available, if any.
        //
        EnterCriticalSection(criticalSection);
        if (nThreadsRunning == 0 && runsReadyToWrite == NULL) {
            LeaveCriticalSection(criticalSection);
            break;
        }

        SamtoolsRun *runsToWrite = NULL;
        size_t sizeOfRunsToWrite;
        runsToWrite = runsReadyToWrite;
        runsReadyToWrite = NULL;
        sizeOfRunsToWrite = dataSizeOnWriteQueue;
        dataSizeOnWriteQueue = 0;
        ResetEvent(runsReadyToWriteQueueHasEnoughDataEvent);
        SetEvent(runsReadyToWriteNotTooFullEvent);

        LeaveCriticalSection(criticalSection);

        if (NULL != runsToWrite) {
            char *consolodatedOutputBuffer = new char[sizeOfRunsToWrite];
            size_t offsetInConsolodatedOutputBuffer = 0;

            size_t dataSeen = 0;
            for (SamtoolsRun *thisRun = runsToWrite; NULL != thisRun; thisRun = thisRun->next) {
                nRunsWritten++;
                thisRun->fileOffset = offsetInOutputFile;
                offsetInOutputFile += thisRun->size;
                dataSeen += thisRun->size;

                for (ReadsBuffer *thisBuffer = thisRun->readsBuffers; NULL != thisBuffer; thisBuffer = thisBuffer->next) {
                    memcpy(consolodatedOutputBuffer + offsetInConsolodatedOutputBuffer, thisBuffer->buffer, thisBuffer->size);
                    offsetInConsolodatedOutputBuffer += thisBuffer->size;
                }
            }

            if (dataSeen != sizeOfRunsToWrite) {
                fprintf(stderr, "Didn't see expected amount of data in write buffers, %lld != %lld\n", dataSeen, sizeOfRunsToWrite);
                exit(1);
            }

            size_t sizeLeftToWrite = sizeOfRunsToWrite;
			size_t amountOfConsolodatedOutputBufferUsed = 0;

            while (0 != sizeLeftToWrite) {
                DWORD bytesWritten;
                DWORD bytesToWrite;
                if (sizeLeftToWrite <= 1024 * 1024 * 1024) {
                    bytesToWrite = (DWORD)sizeLeftToWrite;
                } else {
                    bytesToWrite = 1024 * 1024 * 1024;
                }

                if (!WriteFile(hOutputFile, consolodatedOutputBuffer + amountOfConsolodatedOutputBufferUsed, bytesToWrite, &bytesWritten, NULL)) {
                    fprintf(stderr, "Write to output file failed, %d\n", GetLastError());
                    exit(1);
                }

                if (0 == bytesWritten) {
                    fprintf(stderr, "Empty write on output file!\n");
                    exit(1);
                }

                sizeLeftToWrite -= bytesWritten;
				amountOfConsolodatedOutputBufferUsed += bytesWritten;
            }

            nBytesWritten += sizeOfRunsToWrite;

            delete[] consolodatedOutputBuffer;
            consolodatedOutputBuffer = NULL;

            //
            // Now run through the runs and release their memory.
            //
            for (SamtoolsRun *thisRun = runsToWrite; TRUE /* exit via break below */; thisRun = thisRun->next) {
                while (thisRun->readsBuffers != NULL) {
                    ReadsBuffer *bufferToDelete = thisRun->readsBuffers;
                    thisRun->readsBuffers = bufferToDelete->next;
                    delete[] bufferToDelete->buffer;
                    delete bufferToDelete;
                }

                //
                // Remember the runs we just wrote (if this is the last run on the list)
                //
                if (thisRun->next == NULL) {
                    thisRun->next = runsWritten;
                    runsWritten = runsToWrite;
                    break;
                }
            }
        }  // If we grabbed any work.
     } // Loop while there are still live collectors or work to do.

     CloseHandle(hOutputFile);
     hOutputFile = NULL;

     //
     // Write the index.  It's just a list of what the filename would have been followed by the offset and size in the file of the data, sorted by filename order, in ascii.
     // We start by sorting the input.
     //

     SamtoolsRun **allRuns = new SamtoolsRun *[nRunsWritten];
     for (int i = 0; i < nRunsWritten; i++) {
         allRuns[i] = runsWritten;
         runsWritten = runsWritten->next;
     }

     _ASSERT(NULL == runsWritten);

     qsort(allRuns, nRunsWritten, sizeof(*allRuns), CompareSamtoolsRuns);

     for (int i = 0; i < nRunsWritten; i++) {
         fprintf(indexFile, "%s\t%lld\t%lld\n", allRuns[i]->filename, allRuns[i]->fileOffset, allRuns[i]->size);
     }

     printf("Wrote data for %d samtools runs, %lld bytes in %llds\n", nRunsWritten, nBytesWritten, (GetTickCount64() - startTime + 500) / 1000);

     SetEvent(writerDoneEvent);
     return 0;
}

DWORD __stdcall
ProcessorThreadMain(void *param)
{
    SECURITY_ATTRIBUTES sAttr;
    sAttr.bInheritHandle = TRUE;
    sAttr.lpSecurityDescriptor = NULL;
    sAttr.nLength = sizeof(sAttr);

    int nRunsThisThread = 0;

    for (int lineToProcess = InterlockedIncrement(&nextSamtoolsRunToProcess) - 1; lineToProcess < nInputLines; lineToProcess = InterlockedIncrement(&nextSamtoolsRunToProcess) - 1) {
        if (lineToProcess != 0 && lineToProcess % 1000 == 0) {
            if (lineToProcess % 100000 == 0) {
                printf("\n%d:", lineToProcess);
            }
            printf(".");
            fflush(stdout);
        }
        if (WAIT_OBJECT_0 != WaitForSingleObject(runsReadyToWriteNotTooFullEvent, INFINITE)) {
            fprintf(stderr, "Waiting for write queue not to be overfull failed, %d\n", GetLastError());
            exit(1);
        }

        char *line = inputLines[lineToProcess];

        //
        // Line should be of the form "samtools ... > outputFilename".
        //
        char *redirectChar = strchr(line, '>');

        if (redirectChar == NULL) {
            fprintf(stderr, "Unable to parse samtools line '%s', line number %d\n", line, lineToProcess + 1); // +1 converts to 1 based, which most text editors use
            continue;
        }

        *redirectChar = '\0';
        char *outputFilename = redirectChar + 1;
        while (*outputFilename == ' ') {
            outputFilename++;
        }

        //
        // There's a common problem caused by samtools being case sensitive to the contig name, which generally causes problems with
        // sex chromosomes, because we'll generate chrx and it will expect chrX.  So, if we get an error the first time through we'll
        // capitalize the last letter of the contig name and retry once.
        //
        SamtoolsRun *run;
        BOOL worked = false;
        for (int samtoolsRetryCount = 0; samtoolsRetryCount < 2; samtoolsRetryCount++) {
            HANDLE hStdOutRead, hStdOutWrite;
            if (!CreatePipe(&hStdOutRead, &hStdOutWrite, &sAttr, 0)) {
                fprintf(stderr, "CreatePipe failed, %d\n", GetLastError());
                exit(1);
            }

            if (!SetHandleInformation(hStdOutRead, HANDLE_FLAG_INHERIT, 0)) {   // Make sure the child doesn't inhert the read handle, only the write one.
                fprintf(stderr, "SetHandleInformation failed, %d\n", GetLastError());
                exit(1);
            }

            HANDLE hStdErrRead, hStdErrWrite;
            if (!CreatePipe(&hStdErrRead, &hStdErrWrite, &sAttr, 0)) {
                fprintf(stderr, "CreatePipe failed, %d\n", GetLastError());
                exit(1);
            }

            if (!SetHandleInformation(hStdErrRead, HANDLE_FLAG_INHERIT, 0)) {   // Make sure the child doesn't inhert the read handle, only the write one.
                fprintf(stderr, "SetHandleInformation failed, %d\n", GetLastError());
                exit(1);
            }

            PROCESS_INFORMATION piProcInfo;
            STARTUPINFO siStartInfo;

            ZeroMemory(&piProcInfo, sizeof(PROCESS_INFORMATION));

            ZeroMemory(&siStartInfo, sizeof(STARTUPINFO));
            siStartInfo.cb = sizeof(STARTUPINFO);
            siStartInfo.hStdError = hStdErrWrite;
            siStartInfo.hStdOutput = hStdOutWrite;
            siStartInfo.hStdInput = NULL;
            siStartInfo.dwFlags |= STARTF_USESTDHANDLES;

            BOOL createProcessWorked;
            for (int retryCount = 0; retryCount <= 20; retryCount++) {
                createProcessWorked = CreateProcess(samtoolsPathname, line, NULL, NULL, TRUE, 0, NULL, NULL, &siStartInfo, &piProcInfo);
                if (createProcessWorked) {
                    break;
                }

                if (ERROR_PATH_NOT_FOUND == GetLastError()) {
                    fprintf(stderr, "Path not found attempting to start samtools at '%s', aborting.\n", samtoolsPathname);
                    exit(1);
                }

                fprintf(stderr, "CreateProcess failed, %d, pausing and retrying, retry count %d\n", GetLastError(), retryCount);

                Sleep(4 << retryCount); // when retryCount is 20, this is a little more than an hour.
            }
            if (!createProcessWorked) {
                fprintf(stderr, "Too many retries on CreateProcess(%s...), giving up.\n", line);
                exit(1);
            }
            CloseHandle(piProcInfo.hThread);    // Don't need this one.
            piProcInfo.hThread = NULL;

            //
            // Have to close the far ends of these pipes, else we don't get EOF when the child process exits, because we could still conceivably write to our
            // handle to the pipe.
            //
            CloseHandle(hStdOutWrite);
            CloseHandle(hStdErrWrite);

            size_t readsBufferSize = 16 * 1024; // We grow this exponentially as we get more output.

            nRunsThisThread++;

            run = new SamtoolsRun;
            run->filename = outputFilename;
            run->size = 0;
            run->readsBuffers = new ReadsBuffer;
            run->readsBuffers->buffer = new char[readsBufferSize];
            run->readsBuffers->size = 0;
            run->readsBuffers->next = NULL;

            ReadsBuffer *readsBuffer = run->readsBuffers;

            //
            // Now read the output from samtools.
            //
            DWORD bytesRead;
            for (;;) {
                if (!ReadFile(hStdOutRead, readsBuffer->buffer + readsBuffer->size, readsBufferSize - readsBuffer->size, &bytesRead, NULL)) {
                    if (GetLastError() == ERROR_BROKEN_PIPE) {
                        break;
                    } else {
                        fprintf(stderr, "Read from stdout pipe to samtools failed, %d\n", GetLastError());
                        exit(1);
                    }
                }

                if (0 == bytesRead) {
                    break;
                }

                readsBuffer->size += bytesRead;
                run->size += bytesRead;
                if (readsBuffer->size == readsBufferSize) {
                    //
                    // Allocate a new buffer.
                    //
                    readsBufferSize *= 2;
                    readsBuffer->next = new ReadsBuffer;
                    readsBuffer = readsBuffer->next;

                    readsBuffer->buffer = new char[readsBufferSize];
                    readsBuffer->size = 0;
                    readsBuffer->next = NULL;
                }
            }

            //
            // See if there's anything on stderr.
            //
            const size_t stdErrorBufferSize = 10 * 1024;
            char stdErrorBuffer[stdErrorBufferSize + 1];
            bool anythingWrittenToStdError = false;
            while (ReadFile(hStdErrRead, stdErrorBuffer, stdErrorBufferSize, &bytesRead, NULL)) {
                if (0 == bytesRead) {
                    break;
                }
                stdErrorBuffer[bytesRead] = '\0';
                if (!anythingWrittenToStdError) {
                    if (0 != samtoolsRetryCount) {
                        fprintf(stderr, "stderr output for line %d '%s': %s", lineToProcess + 1, line, stdErrorBuffer);
                    }
                    anythingWrittenToStdError = true;
                } else if (0 != samtoolsRetryCount) {
                    fprintf(stderr, "%s", stdErrorBuffer);
                }
            }

            //
            // Clean up the samtools process.
            //
            WaitForSingleObject(piProcInfo.hProcess, INFINITE); // Wait for the process to really exit.

            DWORD exitCode;
            if (!GetExitCodeProcess(piProcInfo.hProcess, &exitCode)) {
                fprintf(stderr, "GetProcessExitCode failed, %d\n", GetLastError());
                CloseHandle(hStdOutRead);
                CloseHandle(hStdErrRead);
                CloseHandle(piProcInfo.hProcess);
            } else {

                CloseHandle(hStdOutRead);
                CloseHandle(hStdErrRead);
                CloseHandle(piProcInfo.hProcess);

                if (exitCode != 0) {
                    fprintf(stderr, "\nsamtools run failed with exit code %d, line number %d, command line: %s\n", exitCode, lineToProcess + 1, line);
                } else if (0 == samtoolsRetryCount && anythingWrittenToStdError) {
                    while (NULL != run->readsBuffers) {
                        ReadsBuffer *toDelete = run->readsBuffers;
                        run->readsBuffers = toDelete->next;
                        delete[] toDelete->buffer;
                        delete toDelete;
                    }
                    delete run;

                    //
                    // Now fix up line to have a cap right before the : in the fourth item in the line or change chrun to chrUn.  Get there by scanning for the third space.
                    //
                    char *linePtr = line;
                    for (int spaceCount = 0; spaceCount < 3; spaceCount++) {
                        linePtr = strchr(linePtr + 1, ' ');
                        if (NULL == linePtr) {
                            fprintf(stderr, "Couldn't find enough spaces in input line to capitalize for retry: %s\n", line);
                            break;
                        }
                    }

                    //
                    // See if this starts with chrun_.  If so, capitalize the U and the two characters after the underscore.
                    //
                    if (!strncmp(" chrun_", linePtr, 7) && strlen(linePtr) > 8) {
                        linePtr[4] = 'U';

                        for (int i = 7; i <= 8; i++) {
                            linePtr[i] = toupper(linePtr[i]);
                        }
                    } else {
                        //
                        // Now find the colon.
                        //
                        linePtr = strchr(linePtr, ':');
                        if (NULL == linePtr) {
                            fprintf(stderr, "Unable to find colon in input line: %s\n", line);
                            samtoolsRetryCount++;   // Don't retry
                            break;
                        }

                        //
                        // Back up to the character before the colon, and capitalize it if possible.
                        // 
                        linePtr--;

                        if (*linePtr >= 'a' && *linePtr <= 'z') {
                            *linePtr += 'A' - 'a';
                        } else {
                            fprintf(stderr, "Not retrying because character before the colon wasn't a lower case letter.  Line: %s\n", line);
                            samtoolsRetryCount++;   // Don't retry
                            break;
                        }
                    }
                } else {
                    //
                    // It worked, no need to retry.
                    //
                    worked = true;
                    break;
                }
            }
        }

        //
        // And add our buffer to the queue for the writer.
        //
        if (worked) {
            EnterCriticalSection(criticalSection);
            run->next = runsReadyToWrite;
            runsReadyToWrite = run;
            dataSizeOnWriteQueue += run->size;

            if (dataSizeOnWriteQueue >= MinWriteSize) {
                SetEvent(runsReadyToWriteQueueHasEnoughDataEvent);

                if (dataSizeOnWriteQueue >= MaxWriteQueueSize) {
                    ResetEvent(runsReadyToWriteNotTooFullEvent);
                }
            }
            run = NULL;
            LeaveCriticalSection(criticalSection);
        }
    }

    if (InterlockedDecrement(&nThreadsRunning) == 0) {
        //
        // If we're the last processor, wake up the writer even if there's not enough data left.
        //
        SetEvent(runsReadyToWriteQueueHasEnoughDataEvent);
    }

    return 0;
}

int main(int argc, char* argv[])
{

    if (argc < 3 || argc > 5) usage();

    nThreadsRunning = GetNumberOfProcessors();

    if (argc > 3) {
        int arg = atoi(argv[3]);
        if (arg > 0) {
            nThreadsRunning = arg;
        }
    }

    if (argc > 4) {
        samtoolsPathname = argv[4];
    }


    //
    // Read in the input file.
    //
    struct InputFileLine {
        char *data;
        InputFileLine *next;
    } *inputFileLines = NULL;

    const size_t inputFileLineBufferSize = 1000;
    char inputFileLineBuffer[inputFileLineBufferSize + 1];

    FILE *inputFile;
    if (!strcmp(argv[1], "-")) {
        inputFile = stdin;
    } else {
        inputFile = fopen(argv[1], "r");
        if (NULL == inputFile) {
            fprintf(stderr, "Unable to open input file '%s'\n", argv[1]);
            return 1;
        }
    }

    size_t inputFileSize = 0;

    char *nextLine;
    while (NULL != (nextLine = fgets(inputFileLineBuffer, inputFileLineBufferSize, inputFile))) {
        InputFileLine *nextLine = new InputFileLine();
        inputFileLineBuffer[inputFileLineBufferSize] = '\0';    // The buffer is one bigger than "size"
        nextLine->data = new char[strlen(inputFileLineBuffer) + 1];
        strcpy(nextLine->data, inputFileLineBuffer);
        nextLine->next = inputFileLines;    // Yes, this revereses the order, but we run them in parallel out of order anyway.
        inputFileLines = nextLine;
        inputFileSize += strlen(inputFileLineBuffer);
    }
    fclose(inputFile);


    char *inputFileData = new char[inputFileSize + 1];
    inputFileData[inputFileSize] = '\0';
    size_t bytesProcessed = 0;
    while (NULL != inputFileLines) {
        size_t lineSize = strlen(inputFileLines->data);
        if (bytesProcessed + lineSize > inputFileSize) {
            fprintf(stderr, "Code bug!  Would copy over the end of the buffer.\n");
            return 1;
        }
        memcpy(inputFileData + bytesProcessed, inputFileLines->data, lineSize);
        bytesProcessed += lineSize;

        InputFileLine *dying = inputFileLines;
        inputFileLines = dying->next;
        delete[] dying->data;
        delete dying;
    }

    if (bytesProcessed != inputFileSize) {
        fprintf(stderr, "Got a different number of bytes while processing the input file, %lld != %lld.  Code bug.\n", bytesProcessed, inputFileSize);
        return 1;
    }
    inputFileData[inputFileSize] = '\0';


    hOutputFile = CreateFile(argv[2], GENERIC_READ | GENERIC_WRITE, FILE_SHARE_READ, NULL, CREATE_ALWAYS, FILE_ATTRIBUTE_NORMAL, NULL);
    if (INVALID_HANDLE_VALUE == hOutputFile) {
        fprintf(stderr, "Unable to open output file '%s', %d\n", argv[2], GetLastError());
        exit(1);
    }

    char *indexFilename = new char[strlen(argv[2]) + 7];
    sprintf(indexFilename, "%s.index", argv[2]);
    indexFile = fopen(indexFilename, "w");
    if (NULL == indexFile) {
        fprintf(stderr, "Unable to open index file '%s'\n", indexFilename);
        exit(1);
    }

    InitializeCriticalSection(criticalSection);
    runsReadyToWriteQueueHasEnoughDataEvent = CreateEvent(NULL, TRUE, FALSE, NULL);   // Manual reset event that starts clear (nothing is on the queue)
    runsReadyToWriteNotTooFullEvent = CreateEvent(NULL, TRUE, TRUE, NULL);      // Manual reset event that starts set (the queue isn't overfull)
    writerDoneEvent = CreateEvent(NULL, TRUE, FALSE, NULL);      // Manual reset event that starts clear (the writer isn't done before we start)

    if (NULL == runsReadyToWriteQueueHasEnoughDataEvent || NULL == runsReadyToWriteNotTooFullEvent || NULL == writerDoneEvent) {
        fprintf(stderr, "CreateEvent failed\n");
        exit(1);
    }


    //
    // Now find the line breaks in the input file.  Make two passes: one to count them and then a second one to record them.
    // Writing this makes me realize how much easier this would be in C#.  :-(
    //

    nInputLines = 0;
    for (char *ptr = inputFileData; ptr < inputFileData + inputFileSize; ptr++) {
        if (*ptr == '\r') {
            //
            // Just blow away the annoying cr in crlf text.
            //
            *ptr = '\0';
        }
        if (*ptr == '\n') {
            nInputLines++;
        }
    }

    inputLines = new char *[nInputLines];
    inputLines[0] = inputFileData;

    int whichInputLine = 1;
    for (char *ptr = inputFileData; ptr < inputFileData + inputFileSize; ptr++) {
        if (*ptr == '\n') {
            if (whichInputLine >= nInputLines) {
                if (whichInputLine > nInputLines) {
                    fprintf(stderr, "Strange..saw more input lines the second time around.\n");
                    exit(1);
                }
                //
                // Just ignore the last one, since we're recording the next line at the previous one's \n
                //
            } else {
                inputLines[whichInputLine] = ptr + 1;
            }
            *ptr = '\0';    // Make each line its own string.
            whichInputLine++;
        }
    }

    if (1 == nInputLines && !strcmp(inputLines[0], "This script intentionally left blank.")) {
        nInputLines = 0;
    } else if (0 == nInputLines) {
        fprintf(stderr, "Empty script file (which doesn't happen on purpose).\n");
        exit(1);
    }

    printf("Extracting %d samples\n", nInputLines);
    fflush(stdout);

    //
    // Now launch the worker threads.  First the writer, then the producers.
    //
    DWORD writerThreadId;
    if (!CreateThread(NULL, 0, WriterThreadMain, NULL, 0, &writerThreadId)) {
        fprintf(stderr, "Unable to create writer thread, %d\n", GetLastError());
        exit(1);
    }

    long nThreadsToStart = nThreadsRunning;
    for (int i = 0; i < nThreadsToStart; i++) {
        DWORD processorThreadId;
        if (!CreateThread(NULL, 0, ProcessorThreadMain, NULL, 0, &processorThreadId)) {
            fprintf(stderr, "Unable to create processor thread, %d\n", GetLastError());
            exit(1);
        }
    }

    if (WAIT_OBJECT_0 != WaitForSingleObject(writerDoneEvent, INFINITE)) {
        fprintf(stderr, "Wait for overall done object failed, %d\n", GetLastError());
        exit(1);
    }

    fprintf(indexFile, "**done**\n");
    fclose(indexFile);

	return 0;
}

