/*++

Module Name:

    FASTQ.cpp

Abstract:

    Fast FASTQ genome "query" reader.

Authors:

    Bill Bolosky, August, 2011

Environment:

    User mode service.

    This class is NOT thread safe.  It's the caller's responsibility to ensure that
    at most one thread uses an instance at any time.

Revision History:

    Adapted from Matei Zaharia's Scala implementation.

--*/

#include "stdafx.h"
#include "FASTQ.h"
#include "Compat.h"
#include "BigAlloc.h"
#include "Read.h"

using std::min;


FASTQReader::~FASTQReader()
{
}


    FASTQReader*
FASTQReader::create(const char *fileName, _int64 startingOffset, _int64 amountOfFileToProcess, ReadClippingType clipping)
{
#ifdef _MSC_VER
  return new WindowsFASTQReader(fileName, startingOffset, amountOfFileToProcess, clipping);
#else
  return new MemMapFASTQReader(fileName, startingOffset, amountOfFileToProcess, clipping);
#endif
}


#ifndef _MSC_VER

MemMapFASTQReader::MemMapFASTQReader(const char* fileName, _int64 startingOffset, _int64 amountOfFileToProcess, ReadClippingType i_clipping)
{
    fd = open(fileName, O_RDONLY);
    if (fd == -1) {
	fprintf(stderr, "Failed to open %s\n", fileName);
	exit(1);
    }

    struct stat sb;
    int r = fstat(fd, &sb);
    if (r == -1) {
	fprintf(stderr, "Failed to stat %s\n", fileName);
	exit(1);
    }
    fileSize = sb.st_size;

    fileData = NULL;

    clipping = i_clipping;

    reinit(startingOffset, amountOfFileToProcess);
}
    

void MemMapFASTQReader::reinit(_int64 startingOffset, _int64 amountOfFileToProcess) {
    unmapCurrentRange();

    _int64 misalignment = (startingOffset % getpagesize());
    _int64 alignedOffset = startingOffset - misalignment;

    size_t amountToMap = min((_uint64) amountOfFileToProcess + misalignment + maxReadSizeInBytes,
                             (_uint64) fileSize - alignedOffset);
    //printf("Going to map %llu bytes starting at %lld (amountOfFile=%lld)\n", amountToMap, alignedOffset, amountOfFileToProcess);

    fileData = (char *) mmap(NULL, amountToMap, PROT_READ, MAP_SHARED, fd, alignedOffset);
    if (fileData == MAP_FAILED) {
	fprintf(stderr, "mmap failed on FASTQ file\n");
	exit(1);
    }

    int r = madvise(fileData, min((size_t) madviseSize, amountToMap), MADV_WILLNEED);
    _ASSERT(r == 0);
    lastPosMadvised = 0;

    pos = misalignment;
    endPos = pos + amountOfFileToProcess;
    offsetMapped = alignedOffset;
    amountMapped = amountToMap;

    // If we're not at the start of the file, we might have the tail end of a read that someone
    // who got the previous range will process; advance past that. This is fairly tricky because
    // there can be '@' signs in the quality string (and maybe even in read names?). We'll do it
    // by trying to load a read starting at each offset until we get one that works.
    if (startingOffset != 0) {
        Read r;
        while (pos < endPos && parseRead(pos, &r, false) == 0) {
            pos++;
        }
        //printf("Skipped %lld (%lld) bytes to get to a read\n", pos-misalignment, pos);
    }
}


void MemMapFASTQReader::unmapCurrentRange()
{
    if (fileData != NULL) {
        munmap(fileData, amountMapped);
        fileData = NULL;
    }
}


MemMapFASTQReader::~MemMapFASTQReader()
{
    unmapCurrentRange();
    close(fd);
}

    bool
MemMapFASTQReader::getNextRead(Read *readToUpdate)
{
    _uint64 newPos = parseRead(pos, readToUpdate, true);
    if (newPos == 0) {
        return false;
    } else {
        pos = newPos;
        return true;
    }
}

//
// Try to parse a read starting at a given position pos, updating readToUpdate with it.
// Returns 0 if the parse failed or the first position past the read if it succeeds. In
// addition, if exitOnFailure is set, print a warning and exit if there is a parse error.
// (The one time we don't set this is when trying to find the first read in a chunk.)
// 
    _uint64
MemMapFASTQReader::parseRead(_uint64 pos, Read *readToUpdate, bool exitOnFailure)
{
    if (pos >= endPos) {
        return 0;
    }

    if (fileData[pos] != '@') {
        if (exitOnFailure) {
            fprintf(stderr, "Syntax error in FASTQ file at offset %llu. "
                "Expected line to start with '@' but it didn't\n", (_uint64) (pos + offsetMapped));
            exit(1);
        } else {
            return 0;
        }
    }

    char *id = fileData + pos + 1; // The '@' is not part of the ID
    while (pos < amountMapped && fileData[pos] != '\n') {
        pos++;
    }
    if (pos == amountMapped) {
        if (exitOnFailure) {
            fprintf(stderr, "Syntax error in FASTQ file at offset %llu. "
                "Reached end of range before a complete read.\n", (_uint64) (pos + offsetMapped));
            exit(1);
        } else {
            return 0;
        }
    }
    unsigned idLen = (unsigned) (fileData + pos - id);
    pos++;

    char *data = fileData + pos;
    while (pos < amountMapped && fileData[pos] != '\n') {
        pos++;
    }
    if (pos == amountMapped) {
        if (exitOnFailure) {
            fprintf(stderr, "Syntax error in FASTQ file at offset %llu. "
                "Reached end of range before a complete read.\n", (_uint64) (pos + offsetMapped));
            exit(1);
        } else {
            return 0;
        }
    }
    unsigned dataLen = (unsigned) (fileData + pos - data);
    pos++;

    if (pos >= amountMapped || fileData[pos] != '+') {
        if (exitOnFailure) {
            fprintf(stderr, "Syntax error in FASTQ file at offset %llu. "
                    "Expected line starting with '+'.\n", (_uint64) (pos + offsetMapped));
            exit(1);
        } else {
            return 0;
        }
    }
    pos++;

    while (pos < amountMapped && fileData[pos] != '\n') {
        pos++;
    }
    if (pos == amountMapped) {
        if (exitOnFailure) {
            fprintf(stderr, "Syntax error in FASTQ file at offset %llu. "
                    "Reached end of range before a complete read.\n", (_uint64) (pos + offsetMapped));
            exit(1);
        } else {
            return 0;
        }
    }
    pos++;

    char *quality = fileData + pos;
    if (pos + dataLen > amountMapped) {
        if (exitOnFailure) {
            fprintf(stderr, "Syntax error in FASTQ file at offset %llu. "
                    "Quality string is not long enough.", (_uint64) (pos + offsetMapped));
            exit(1);
        } else {
            return 0;
        }
    } else if (pos + dataLen == amountMapped) {
      pos = amountMapped;
    } else {
      if (fileData[pos + dataLen] != '\n') {
          if (exitOnFailure) {
              fprintf(stderr, "Syntax error in FASTQ file at offset %llu. "
                      "Quality string is of the wrong length.", (_uint64) (pos + offsetMapped));
              exit(1);
          } else {
              return 0;
          }
      }
      pos += dataLen + 1;
    }

    // Call madvise() to (a) start reading more bytes if we're past half our current
    // range and (b) tell the OS we won't need any stuff we've read in the past
    if (pos > lastPosMadvised + madviseSize / 2) {
        _uint64 offset = lastPosMadvised + madviseSize;
        _uint64 len = (offset > amountMapped ? 0 : min(amountMapped - offset, (_uint64) madviseSize));
        if (len > 0) {
            // Start reading new range
            int r = madvise(fileData + offset, len, MADV_WILLNEED);
            _ASSERT(r == 0);
        }
        if (lastPosMadvised > 0) {
          // Unload the range we had before our current one
          int r = madvise(fileData + lastPosMadvised - madviseSize, madviseSize, MADV_DONTNEED);
          _ASSERT(r == 0);
        }
        lastPosMadvised = offset;
    }

    readToUpdate->init(id, idLen, data, quality, dataLen);
    readToUpdate->clip(clipping);
    return pos;
}


    void
MemMapFASTQReader::readDoneWithBuffer(unsigned *referenceCount)
{
    // Do nothing, since file is memory-mapped.
}

#endif /* not _MSC_VER */


#ifdef _MSC_VER

WindowsFASTQReader::WindowsFASTQReader(const char *fileName, _int64 startingOffset, 
                                        _int64 amountOfFileToProcess, ReadClippingType i_clipping)
{
    didInitialSkip = false;
    clipping = i_clipping;
    //
    // Initilize the buffer info struct.
    //
    for (unsigned i = 0 ; i < nBuffers; i++) {
        bufferInfo[i].buffer = (char *)BigAlloc(bufferSize + 1);    // +1 gives us a place to put a terminating null
        if (NULL == bufferInfo[i].buffer) {
            fprintf(stderr,"FASTQ Reader: unable to allocate IO buffer\n");
            exit(1);
        }

        bufferInfo[i].buffer[bufferSize] = 0;       // The terminating null.
        
        bufferInfo[i].lap.hEvent = CreateEvent(NULL,TRUE,FALSE,NULL);
        if (NULL == bufferInfo[i].lap.hEvent) {
            fprintf(stderr,"Unable to create event for FASTQ reader\n");
            exit(1);
        }

        bufferInfo[i].state = Empty;
        bufferInfo[i].isEOF = false;
        bufferInfo[i].offset = 0;
        bufferInfo[i].referenceCount = 0;
    }

    nextBufferForReader = 0;
    nextBufferForConsumer = 0;

    //
    // Open the file.
    //
    hFile = CreateFile(fileName,GENERIC_READ,FILE_SHARE_READ,NULL,OPEN_EXISTING,FILE_FLAG_OVERLAPPED/*|FILE_FLAG_SEQUENTIAL_SCAN*/,NULL);
    if (INVALID_HANDLE_VALUE == hFile) {
        fprintf(stderr,"Unable to open FASTQ file '%s', %d\n",fileName,GetLastError());
        exit(1);
    }

    //
    // Get the file size so we can compute our ranges.
    //
    if (!GetFileSizeEx(hFile,&fileSize)) {
        fprintf(stderr,"FASTQ reader: unable to get file size, %d\n",GetLastError());
        exit(1);
    }

    readOffset.QuadPart = startingOffset;
    if (0 == amountOfFileToProcess) {
        if (0 != startingOffset) {
            fprintf(stderr,"FASTQ: if amount to process is 0, so must be starting offset (which means process the whole file)\n");
            exit(1);
        }
        endingOffset = fileSize.QuadPart;
    } else {
        endingOffset = startingOffset + amountOfFileToProcess;
    }

    //
    // Initialize the line length hints.    BB: This is causing problems for short reads, so we'll just skip the optimization for now.
    //
    minLineLengthSeen[0] = 1;
    minLineLengthSeen[1] = 1;
    minLineLengthSeen[2] = 1;               // This is the line that just contains a "+"
    minLineLengthSeen[3] = 1;

    //
    // And the isValidStartingCharacterForNextLine array.
    //
    for (unsigned i = 0; i < nLinesPerFastqQuery; i++) {
        for (unsigned j = 0; j < 256; j++) {
            isValidStartingCharacterForNextLine[i][j] = false;
        }
    }

    //
    // The first line is the descriptor line and must start with an '@'
    //
    isValidStartingCharacterForNextLine[3]['@'] = true;

    //
    // The second line is the read itself and must start with a base or an
    // 'N'.
    //
    isValidStartingCharacterForNextLine[0]['A'] = true;
    isValidStartingCharacterForNextLine[0]['C'] = true;
    isValidStartingCharacterForNextLine[0]['T'] = true;
    isValidStartingCharacterForNextLine[0]['G'] = true;
    isValidStartingCharacterForNextLine[0]['N'] = true;

    //
    // The third line is additional sequence idenfitier info and must
    // start with a '+'.
    //
    isValidStartingCharacterForNextLine[1]['+'] = true;

    //
    // Line 4 is the quality line.  It starts with a printable ascii character.
    // We're ruling out the bases, N, + and @ because it might confsue the parser.
    //
    for (char i = 'a'; i <= 'z'; i++) {
        isValidStartingCharacterForNextLine[2][i] = true;
    }
    for (char i = 'A'; i <= 'Z'; i++) {
        isValidStartingCharacterForNextLine[2][i] = true;
    }
#if     0   // Some reads use these letters as quality, we'll just have to deal with them.
    isValidStartingCharacterForNextLine[2]['A'] = false;
    isValidStartingCharacterForNextLine[2]['T'] = false;
    isValidStartingCharacterForNextLine[2]['C'] = false;
    isValidStartingCharacterForNextLine[2]['G'] = false;
    isValidStartingCharacterForNextLine[2]['N'] = false;
#endif  // 0

    for (char i = '0'; i <= '9'; i++) {
        isValidStartingCharacterForNextLine[2][i] = true;
    }
    isValidStartingCharacterForNextLine[2]['-'] = true;
    isValidStartingCharacterForNextLine[2]['!'] = true;
    isValidStartingCharacterForNextLine[2]['#'] = true;
    isValidStartingCharacterForNextLine[2]['$'] = true;
    isValidStartingCharacterForNextLine[2]['%'] = true;
    isValidStartingCharacterForNextLine[2]['^'] = true;
    isValidStartingCharacterForNextLine[2]['&'] = true;
    isValidStartingCharacterForNextLine[2]['*'] = true;
    isValidStartingCharacterForNextLine[2]['('] = true;
    isValidStartingCharacterForNextLine[2][')'] = true;
    isValidStartingCharacterForNextLine[2]['_'] = true;
    isValidStartingCharacterForNextLine[2]['='] = true;
    isValidStartingCharacterForNextLine[2]['{'] = true;
    isValidStartingCharacterForNextLine[2]['}'] = true;
    isValidStartingCharacterForNextLine[2]['|'] = true;
    isValidStartingCharacterForNextLine[2]['\\'] = true;
    isValidStartingCharacterForNextLine[2]['|'] = true;
    isValidStartingCharacterForNextLine[2]['['] = true;
    isValidStartingCharacterForNextLine[2][']'] = true;
    isValidStartingCharacterForNextLine[2][';'] = true;
    isValidStartingCharacterForNextLine[2][':'] = true;
    isValidStartingCharacterForNextLine[2]['\"'] = true;
    isValidStartingCharacterForNextLine[2]['\''] = true;
    isValidStartingCharacterForNextLine[2]['?'] = true;
    isValidStartingCharacterForNextLine[2]['/'] = true;
    isValidStartingCharacterForNextLine[2]['>'] = true;
    isValidStartingCharacterForNextLine[2]['.'] = true;
    isValidStartingCharacterForNextLine[2]['<'] = true;
    isValidStartingCharacterForNextLine[2][','] = true;
    isValidStartingCharacterForNextLine[2]['~'] = true;
    isValidStartingCharacterForNextLine[2]['`'] = true;
    isValidStartingCharacterForNextLine[2]['@'] = true; // You've gotta worry about this one.
    isValidStartingCharacterForNextLine[2]['+'] = true; // This, too

    //
    // Kick off the IO.
    //
    startIo();
}
    void 
WindowsFASTQReader::reinit(_int64 startingOffset, _int64 amountOfFileToProcess)
{
    didInitialSkip = false;
    //
    // Reinitilize the buffer info struct.
    //
    for (unsigned i = 0 ; i < nBuffers; i++) {
        bufferInfo[i].state = Empty;
        bufferInfo[i].isEOF = false;
        bufferInfo[i].offset = 0;
        bufferInfo[i].referenceCount = 0;
    }

    nextBufferForReader = 0;
    nextBufferForConsumer = 0;

    readOffset.QuadPart = startingOffset;
    endingOffset = startingOffset + amountOfFileToProcess;

    startIo();
}


WindowsFASTQReader::~WindowsFASTQReader()
{
    for (unsigned i = 0; i < nBuffers; i++) {
        BigDealloc(bufferInfo[i].buffer);
        bufferInfo[i].buffer = NULL;
        CloseHandle(bufferInfo[i].lap.hEvent);
    }
    CloseHandle(hFile);
}


    void
WindowsFASTQReader::startIo()
{
    //
    // Launch reads on whatever buffers are ready.
    //
    while (bufferInfo[nextBufferForReader].state == Empty) {
        BufferInfo *info = &bufferInfo[nextBufferForReader];

        if (readOffset.QuadPart >= fileSize.QuadPart || readOffset.QuadPart >= endingOffset + maxReadSizeInBytes) {
            info->validBytes = 0;
            info->nBytesThatMayBeginARead = 0;
            info->isEOF = true;
            info->state = Full;
            SetEvent(info->lap.hEvent);
            return;
        }

        unsigned amountToRead;
        if (fileSize.QuadPart - readOffset.QuadPart > bufferSize && endingOffset + maxReadSizeInBytes - readOffset.QuadPart > bufferSize) {
            amountToRead = bufferSize;

            if (readOffset.QuadPart + amountToRead > endingOffset) {
                info->nBytesThatMayBeginARead = (unsigned)(endingOffset - readOffset.QuadPart);
            } else {
                info->nBytesThatMayBeginARead = amountToRead;
            }
        } else {
            amountToRead = (unsigned)__min(fileSize.QuadPart - readOffset.QuadPart,endingOffset+maxReadSizeInBytes - readOffset.QuadPart);
            if (endingOffset <= readOffset.QuadPart) {
                //
                // We're only reading this for overflow buffer.
                //
                info->nBytesThatMayBeginARead = 0;
            } else {
                info->nBytesThatMayBeginARead = __min(amountToRead,(unsigned)(endingOffset - readOffset.QuadPart));    // Don't begin a read past endingOffset
            }
            info->isEOF = true;
        }

        _ASSERT(amountToRead > info->nBytesThatMayBeginARead || !info->isEOF || fileSize.QuadPart == readOffset.QuadPart + amountToRead);
        ResetEvent(info->lap.hEvent);
        info->lap.Offset = readOffset.LowPart;
        info->lap.OffsetHigh = readOffset.HighPart;
        info->fileOffset = readOffset.QuadPart;
         
        if (!ReadFile(
                hFile,
                info->buffer,
                amountToRead,
                &info->validBytes,
                &info->lap)) {

            if (GetLastError() != ERROR_IO_PENDING) {
                fprintf(stderr,"FASTQReader::startIo(): readFile failed, %d\n",GetLastError());
                exit(1);
            }
        }

        readOffset.QuadPart += amountToRead;
        info->state = Reading;
        info->offset = 0;

        nextBufferForReader = (nextBufferForReader + 1) % nBuffers;
    }
}

    bool
WindowsFASTQReader::getNextRead(Read *readToUpdate)
{
    BufferInfo *info = &bufferInfo[nextBufferForConsumer];
    waitForBuffer(nextBufferForConsumer);
    if (info->isEOF && info->offset >= info->validBytes) {
        //
        // EOF.
        //
        return false;
    }

    if (info->offset >= info->nBytesThatMayBeginARead && info->validBytes != info->nBytesThatMayBeginARead) {
        //
        // We have more valid bytes, but it's not ours to handle.  This only happens at the end of our work.
        //
        return false;
    }

    //
    // Get the next four lines.
    //
    char *lines[nLinesPerFastqQuery];
    unsigned line1DataLen;

    unsigned *referenceCounts[2];
    referenceCounts[0] = &info->referenceCount;
    referenceCounts[1] = NULL;

    for (unsigned i = 0; i < nLinesPerFastqQuery; i++) {
        if (info->state != Full) {
            waitForBuffer(nextBufferForConsumer);
        }

        if (!didInitialSkip) {
            //
            // Just assume that a single FASTQ read is smaller than our buffer, so we won't exceed the buffer here.
            // Look for the pattern '{0|\n}@*\n{A|T|C|G|N}*\n+'  That is, either the beginning of the buffer or a
            // newline, followed by an '@' and some text, another newline followed by a list of bases and a newline, 
            // and then a plus.
            //
            _ASSERT(0 == info->offset);

            char *firstLineCandidate = info->buffer;
            if (*firstLineCandidate != '@') {
                firstLineCandidate = strchr(info->buffer,'\n') + 1;
            }

            for (;;) {
                if (firstLineCandidate - info->buffer >= info->nBytesThatMayBeginARead) {
// This happens for very small files.  fprintf(stderr,"Unable to find a read in FASTQ buffer (1)\n");
                    return false;
                }

                char *secondLineCandidate = strchr(firstLineCandidate,'\n') + 1;
                if (NULL == (secondLineCandidate-1)) {
					fprintf(stderr,"Unable to find a read in FASTQ buffer (2) at %d\n", (firstLineCandidate - info->buffer) + info->fileOffset);
                    return false;
                }

                if (*firstLineCandidate != '@') {
                    firstLineCandidate = secondLineCandidate;
                    continue;
                }

                //
                // Scan through the second line making sure it's all bases (or 'N').  We don't have to
                // check for end-of-buffer because we know there's a null there.
                //
                char *thirdLineCandidate = secondLineCandidate;
                while (*thirdLineCandidate == 'A' || *thirdLineCandidate == 'C' || *thirdLineCandidate == 'T' || *thirdLineCandidate == 'G' ||
                       *thirdLineCandidate == 'N') {
                    thirdLineCandidate++;
                }

                if (*thirdLineCandidate == '\r') {
                    //
                    // CRLF text; skip the CR.
                    //
                    thirdLineCandidate++;
                }

                if (*thirdLineCandidate != '\n') {
                    //
                    // We found something that's not a base and not a newline.  It wasn't a read data (second) line.  Move up a line
                    // and try again.
                    //
                    firstLineCandidate = secondLineCandidate;
                    continue;
                }

                thirdLineCandidate++;
                if (*thirdLineCandidate != '+') {
                    firstLineCandidate = secondLineCandidate;
                    continue;
                }

                break;
            }

            info->offset = (unsigned)(firstLineCandidate - info->buffer);   // cast OK because buffer size < 2^32.
            didInitialSkip = true;
        }

        //
        // Find the next newline.  Recall that the buffer is always NULL-terminated.
        // We start by skipping a certain distance into the line to avoid spending
        // a long time in strchr.  We then check that the next line we find
        // starts with the right thing.  If we didn't, then fall back on
        // strchr and update minLineLength.
        //

        char *newLine;
        if (info->offset + minLineLengthSeen[i]+1 >= info->validBytes) {
            //
            // We'd be past end-of-buffer, take the slow path.
            //
            newLine = strchr(info->buffer + info->offset,'\n');
        } else {
            newLine = strchr(info->buffer + info->offset + minLineLengthSeen[i],'\n');
            if (NULL != newLine && !isValidStartingCharacterForNextLine[i][*(newLine +1)] && newLine + 1 - info->buffer != info->validBytes) {
                //
                // We probably overshot.  Fall back on a full scan and update
                // minLineLength.
                //
                newLine = strchr(info->buffer + info->offset,'\n');
                if (NULL != newLine) {
                    if (!isValidStartingCharacterForNextLine[i][*(newLine+1)]) {
                        fprintf(stderr,"FASTQ Parsing error: %c(%d) isn't a valid starting character for line %d\n",*(newLine+1),(unsigned)*(newLine+1),(i+1)%nLinesPerFastqQuery);
                        exit(1);
                    }
                    minLineLengthSeen[i] = __min(newLine - (info->buffer + info->offset), minLineLengthSeen[i]);
                }
            }
        }

        if (NULL == newLine) {
            //
            // There was no next newline; read into the next buffer.
            //
            if (info->isEOF) {
                fprintf(stderr,"FASTQ file doesn't end with a newline!  Failing.  fileOffset = %lld, offset = %d, validBytes = %d, nBytesThatMayBeginARead %d\n",
                    info->fileOffset,info->offset,info->validBytes,info->nBytesThatMayBeginARead);
                exit(1);
            }
            if (bufferInfo[(nextBufferForConsumer + 1) %nBuffers].state != Full) {
                waitForBuffer((nextBufferForConsumer + 1) % nBuffers);
            }

            unsigned amountFromOldBuffer = info->validBytes - info->offset;

            lines[i] = info->overflowBuffer;
            if (amountFromOldBuffer > maxLineLen) {
                fprintf(stderr,"Error parsing FASTQ file.  Either it has very long text lines or it is corrupt.\n");
                exit(1);
            }

            memcpy(lines[i],info->buffer + info->offset, info->validBytes - info->offset);
            info->state = UsedButReferenced;        // The consumer is no longer using this buffer, but it's still referecned by Read(s)
            info->referenceCount++;
            info->offset = info->validBytes;

            nextBufferForConsumer = (nextBufferForConsumer + 1) % nBuffers;
            info = &bufferInfo[nextBufferForConsumer];
            referenceCounts[1] = &info->referenceCount;

            newLine = strchr(info->buffer,'\n');
            _ASSERT(NULL != newLine);
            *newLine = 0;
            if (newLine - info->buffer + 1 + amountFromOldBuffer > maxLineLen) {
                fprintf(stderr,"Error parsing FASTQ file(2).  Either it has very long text lines or it is corrupt.\n");
                exit(1);
            }
            memcpy(lines[i] + amountFromOldBuffer, info->buffer, newLine - info->buffer + 1);

            if (1 == i) {
                line1DataLen = (unsigned)(amountFromOldBuffer + newLine - info->buffer);
            }
        } else {
            lines[i] = info->buffer + info->offset;

            if (1 == i) {
                line1DataLen = (unsigned)(newLine - (info->buffer + info->offset));
            }
        }
        if (newLine > info->buffer + info->offset && *(newLine-1) == '\r') {
            //
            // Kill CR in CRLF text.
            //
            *(newLine-1) = 0;
            if (1 == i) {
                line1DataLen--;
            }
        } else {
            *newLine = 0;
        }

        info->offset = (unsigned)((newLine - info->buffer) + 1);
    }

    info->referenceCount++;

    if ('@' != lines[0][0]) {
        fprintf(stderr,"Syntax error in FASTQ file.  Expected line to start with '@', but it didn't: '%s'\n",lines[0]);
        exit(1);
    }

    if ('+' != lines[2][0]) {
        fprintf(stderr,"Syntax error in FASTQ file.  Expected line to start with '+', but it didn't: '%s'\n",lines[2]);
        exit(1);
    }

    const char *id = lines[0] + 1; // The '@' on the first line is not part of the ID
    readToUpdate->init(id, (unsigned) strlen(id), lines[1], lines[3], line1DataLen, referenceCounts);
    readToUpdate->clip(clipping);

    return true;
}

    void
WindowsFASTQReader::waitForBuffer(unsigned bufferNumber)
{
    BufferInfo *info = &bufferInfo[bufferNumber];

    if (info->state == Full) {
        return;
    }

    if (info->state == UsedButReferenced) {
        fprintf(stderr,"Overlapped buffer manager: waiting for buffer that's in UsedButReferenced.  Almost certainly a bug.\n");
        exit(1);
    }

    if (info->state != Reading) {
        startIo();
    }

    if (!GetOverlappedResult(hFile,&info->lap,&info->validBytes,TRUE)) {
        fprintf(stderr,"Error reading FASTQ file, %d\n",GetLastError());
        exit(1);
    }

    info->state = Full;
    info->buffer[info->validBytes] = 0;
    ResetEvent(info->lap.hEvent);
}


    void
WindowsFASTQReader::readDoneWithBuffer(unsigned *referenceCount)
{
    if (0 != *referenceCount) {
        return;
    }

    BufferInfo *info = NULL;
    for (unsigned i = 0; i < nBuffers; i++) {
        if (&bufferInfo[i].referenceCount == referenceCount) {
            info = &bufferInfo[i];
            break;
        }
    }

    _ASSERT(NULL != info);

    if (info->state == UsedButReferenced) {
        info->state = Empty;

        startIo();
    }
}

#endif /* _MSC_VER */

    FASTQWriter *
FASTQWriter::Factory(const char *filename)
{
    FILE *file = fopen(filename,"wb");
    if (NULL == file) {
        return NULL;
    }

    return new FASTQWriter(file);
}

    bool
FASTQWriter::writeRead(Read *read)
{
    size_t len = __max(read->getIdLength(),read->getDataLength());
    char *buffer = new char [len + 1];
    if (NULL == buffer) {
        fprintf(stderr,"FASTQWriter: out of memory\n");
        exit(1);
    }

    memcpy(buffer,read->getId(),read->getIdLength());
    buffer[read->getIdLength()] = 0;
    bool worked = fprintf(outputFile,"@%s\n",buffer) > 0;

    memcpy(buffer,read->getData(),read->getDataLength());
    buffer[read->getDataLength()] = 0;
    worked &= fprintf(outputFile,"%s\n+\n",buffer) > 0;

    memcpy(buffer,read->getQuality(), read->getDataLength());   // The FASTQ format spec requires the same number of characters in the quality and data strings
    // Null is already in place from writing data

    worked &= fprintf(outputFile,"%s\n",buffer) > 0;

    delete [] buffer;    

    return worked;
}


PairedFASTQReader::~PairedFASTQReader()
{
    for (int i = 0; i < 2; i++) {
        delete readers[i];
        readers[i] = NULL;
    }
}

	PairedFASTQReader *
PairedFASTQReader::create(const char *fileName0, const char *fileName1, _int64 startingOffset, _int64 amountOfFileToProcess, ReadClippingType clipping)
{
    PairedFASTQReader *reader = new PairedFASTQReader;
    reader->readers[0] = FASTQReader::create(fileName0,startingOffset,amountOfFileToProcess,clipping);
    reader->readers[1] = FASTQReader::create(fileName1,startingOffset,amountOfFileToProcess,clipping);

    for (int i = 0; i < 2; i++) {
        if (NULL == reader->readers[i]) {
            delete reader;
            reader = NULL;
            return NULL;
        }
    }

    return reader;
}

    bool
PairedFASTQReader::getNextReadPair(Read *read0, Read *read1)
{
    if (!readers[0]->getNextRead(read0)) {
        return false;
    }

    if (!readers[1]->getNextRead(read1)) {
        fprintf(stderr,"PairedFASTQReader: failed to read mate.  The FASTQ files may not match properly.\n");
        exit(1);
    }

    return true;
}
