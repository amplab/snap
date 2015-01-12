/*++

Module Name:

    wc.cpp

Abstract:

   Version of the standard wc (word count) program that uses 64 bit counters.

Authors:

    Bill Bolosky, September, 2013

Revision History:

    
--*/

#include "stdafx.h"
#include "Compat.h"
#include "exit.h"

void usage()
{
    fprintf(stderr,"usage: wc [-lwc] [files]\n");
    soft_exit(1);
}

struct InputFile {
    InputFile() : lines(0), words(0), chars(0), next(NULL) {}

    char *fileName;
    _uint64  lines, words, chars;
    InputFile *next;
};

SingleWaiterObject allThreadsDone;
volatile _int64 nRunningThreads;

void WorkerThreadMain(void *context)
{
    InputFile *inputFile = (InputFile *)context;

    FILE *file;

    if (!strcmp(inputFile->fileName, "-")) {
        file = stdin;
    } else {
        file = fopen(inputFile->fileName, "rb");
        if (NULL == file) {
            fprintf(stderr,"wc: unable to open input file '%s'\n", inputFile->fileName);
            soft_exit(1);
        }
    }

    const size_t bufferSize = 8 * 1024 * 1024;    // A decent disk IO size
    unsigned char *buffer = new unsigned char[bufferSize];

    _uint64 lines = 0, words = 0, chars = 0;

    //
    // Rather than using conditional branches, use lookup tables.
    //
    int isSeparator[256];
    int isLineBreak[256];
    int isWordPart[256];

    for (int x = 0; x < 256; x++) {
        isSeparator[x] = isLineBreak[x] = isWordPart[x] = 0;
    }

    //
    // Characters that come between words.  Personally, I'd include punctuation, but the old wc doesn't.
    //
    isSeparator[' '] = 1;
    isSeparator['\t'] = 1;
    isSeparator['\n'] = 1;
    isSeparator['\r'] = 1;

    //
    // Characters that end a line.  We don't include "\r" because then CRLF text would look like it has twice
    // as many lines as it does.
    //
    isLineBreak['\n'] = 1;

    //
    // Characters that are parts of words.
    //
    for (int x = 'a'; x <= 'z'; x++) {
        isWordPart[x] = 1;
    }

    for (int x = 'A'; x <= 'Z'; x++) {
        isWordPart[x] = 1;
    }

    for (int x = '0'; x <= '9'; x++) {
        isWordPart[x] = 1;
    }

    //
    // Now process the input file.
    //
    _uint64 notInAWord = 1;

    size_t validBytes;
    while (0 != (validBytes = fread(buffer, 1, bufferSize, file))) {
        chars += validBytes;

        for (size_t i = 0; i < validBytes; i++) {
            unsigned char nextChar = buffer[i];

            lines += isLineBreak[nextChar];

            //
            // Use branch-free logic to compute the word count.
            //
            words += notInAWord * isWordPart[nextChar];    // if (notInAWord && isWordPart) words++

            notInAWord = 1 - ((1 - isSeparator[nextChar]) * (1 - notInAWord)); // notInAWord = isSeparator || notInAWord 
            notInAWord = notInAWord * (1 - isWordPart[nextChar]); // notInAWord = notInAWord && !isWordPart
        } // for each char of the buffer
    }

    if (!feof(file)) {
        fprintf(stderr,"Error reading file '%s'\n", inputFile->fileName);
        soft_exit(1);
    }

    inputFile->chars = chars;
    inputFile->words = words;
    inputFile->lines = lines;

    fclose(file);
    delete [] buffer;
 
    if (0 == InterlockedAdd64AndReturnNewValue(&nRunningThreads, -1)) {
        SignalSingleWaiterObject(&allThreadsDone);
    }
}

void printOutputLine(
    _uint64     chars,
    _uint64     words,
    _uint64     lines,
    const char *fileName,
    bool        printChars,
    bool        printWords,
    bool        printLines)
{
    printf("\t");   // Not sure why, but unix wc does this, so I'll keep it.
    if (printLines) {
        printf("%llu\t", lines);
    }
    if (printWords) {
        printf("%llu\t", words);
    }
    if (printChars) {
        printf("%llu\t", chars);
    }
    printf("%s\n", fileName);
}


int main(int argc, char* argv[])
{
    bool cmdLinePrintChars = false, cmdLinePrintWords = false, cmdLinePrintLines = false, seenStdin = false;

    InputFile *inputFiles = NULL;
    InputFile *lastInputFile = NULL;

    _uint64 nInputFiles = 0;

    for (int i = 1; i < argc; i++) {
        if (argv[i][0] == '-' && argv[i][1] != '\0') {
            //
            // Option.
            //
            for (size_t j = 1; j < strlen(argv[i]); j++) {
                switch (argv[i][j]) {
                    case 'l': cmdLinePrintLines = true; break;
                    case 'w': cmdLinePrintWords = true; break;
                    case 'c': cmdLinePrintChars = true; break;
                    default: usage();
                } // switch
            } // for each char in options
        } else {
            nInputFiles++;

            if (!strcmp(argv[i], "-")) {
                if (seenStdin) {
                    fprintf(stderr,"Can't specify '-' (read input from stdin) more than once\n");
                    soft_exit(1);
                }
                seenStdin = true;
            }
            InputFile *inputFile = new InputFile;
            inputFile->fileName = argv[i];
            if (lastInputFile == NULL) {
                inputFiles = lastInputFile = inputFile;
            } else {
                lastInputFile->next = inputFile;
                lastInputFile = inputFile;
            }
        }
    } // for all args

    if (0 == nInputFiles) {
        //
        // Read from stdin
        //
        InputFile *inputFile = new InputFile;
        inputFile->fileName = "-";
        inputFiles = lastInputFile = inputFile;
        nInputFiles = 1;
    }

    //
    // Fire up a thread for each file.
    //
    CreateSingleWaiterObject(&allThreadsDone);
    nRunningThreads = nInputFiles;

    InputFile *inputFile = inputFiles;
    while (NULL != inputFile) {
        StartNewThread(WorkerThreadMain, inputFile);
        inputFile = inputFile->next;
    }
    WaitForSingleWaiterObject(&allThreadsDone);

	_uint64 chars = 0, words = 0, lines = 0;

    bool printChars, printWords, printLines;

    if (cmdLinePrintChars || cmdLinePrintWords || cmdLinePrintLines) {
        printChars = cmdLinePrintChars;
        printWords = cmdLinePrintWords;
        printLines = cmdLinePrintLines;
    } else {
        printChars = printWords = printLines = true;
    }

    inputFile = inputFiles;
    while (NULL != inputFile) {
        printOutputLine(inputFile->chars, inputFile->words, inputFile->lines, inputFile->fileName, printChars, printWords, printLines);

        chars += inputFile->chars;
        words += inputFile->words;
        lines += inputFile->lines;

        inputFile = inputFile->next;
    }

    if (nInputFiles > 1) {
        printOutputLine(chars, words, lines, "Totals", printChars, printWords, printLines);
    }
}

