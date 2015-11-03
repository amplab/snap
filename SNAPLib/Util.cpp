/*++

Module Name:

    Util.cpp

Abstract:

    Generic support routines that don't seem to belong elsewhere.

Authors:

    Bill Bolosky, March, 2013

Environment:

    User mode service.

Revision History:

    Factored from other places

--*/

#include "stdafx.h"
#include "Util.h"
#include "Error.h"
#include "GenericFile.h"


_int64 FirstPowerOf2GreaterThanOrEqualTo(_int64 value)
{
    int highestBitSet;
    for (highestBitSet = 0; highestBitSet <= 62; highestBitSet++) { // Only go to 63, since this is signed
        if (!(value & ~((((_int64)1) << highestBitSet) - 1))) {
            highestBitSet -= 1;
            break;
        }
    }
    if (((_int64)1) << highestBitSet == value) return value;
    return ((_int64)1) << (highestBitSet + 1);
}

int cheezyLogBase2(_int64 value)
{
    int retVal = 0;
    value /= 2; // Since 2^0 = 1; we'll also define cheezyLogBase2(x) = 0 where x<= 0.
    while (value > 0) {
        retVal++;
        value /= 2;
    }
    return retVal;
}

    void
util::memrevcpy(
    void* dst,
    const void* src,
    size_t bytes)
{
    size_t dwords = bytes >> 3;
    _uint64* p = (_uint64*) dst;
    const _uint64* q = (const _uint64*) ((const char*)src + bytes - 8);
    for (size_t i = 0; i < dwords; i++) {
        *p++ = ByteSwapUI64(*q--);
    }
    int left = (int) (bytes & 7);
    if (left > 0) {
        char* p2 = (char*) p;
        const char* q2 = (left - 1) + (const char*) src;
        for (int i = 0; i < left; i++) {
            *p2++ = *q2--;
        }
    }
}

NWaiter::NWaiter(size_t n)
{
	_signalsRequired = n;
	_signalsReceived = 0;
	InitializeExclusiveLock(&_lock);
	CreateEventObject(&_waiter);
}

NWaiter::~NWaiter()
{
	DestroyExclusiveLock(&_lock);
	DestroyEventObject(&_waiter);
}


void NWaiter::wait()
{
	while (true) {
		bool done;

		AcquireExclusiveLock(&_lock);
		done = (_signalsReceived >= _signalsRequired);
		ReleaseExclusiveLock(&_lock);

		if (done)
			return;
		else {
			WaitForEvent(&_waiter);
		}
	}
}

void NWaiter::signal()
{
	AcquireExclusiveLock(&_lock);
	_signalsReceived += 1;
	ReleaseExclusiveLock(&_lock);

	AllowEventWaitersToProceed(&_waiter);
}

char *FormatUIntWithCommas(_uint64 val, char *outputBuffer, size_t outputBufferSize)
{
	//
	// First, figure out the number of digits.
	//
	unsigned nDigits = 0;
	_uint64 x = val;
	while (x > 0) {
		nDigits++;
		x = x / 10;
	}

	if (0 == nDigits) {
		//
		// Special case for the value 0 (which, I suppose if the world was rational, would be represented by the empty string.  :-))
		//
		_ASSERT(0 == val);
		nDigits = 1;
	}

	int nCommas = (nDigits - 1) / 3;

	if (outputBufferSize < nDigits + nCommas + 1) {
		WriteErrorMessage("Internal error: too small buffer for FormatUIntWithCommas, value %lld, outputBufferSize %lld\n", val, outputBufferSize);
		if (outputBufferSize > 0) {
			*outputBuffer = 0;
		} else {
			soft_exit(1);
		}
		return outputBuffer;
	}

	//
	// Now build up the string backwards.
	//
	size_t offset = nDigits + nCommas;
	outputBuffer[offset] = '\0';

	if (0 == val) {
		outputBuffer[0] = '0';
		return outputBuffer;
	}

	x = val;
	while (x > 0) {
		char tempBuf[5];
		if (x > 999) {
			sprintf(tempBuf, ",%03lld", x % 1000);
			_ASSERT(strlen(tempBuf) == 4);
		} else {
			sprintf(tempBuf, "%d", x);
		}

		_ASSERT(offset >= strlen(tempBuf));
		offset -= strlen(tempBuf);
		memcpy(outputBuffer + offset, tempBuf, strlen(tempBuf));
		x /= 1000;
	}

	return outputBuffer;
}

//
// Version of fgets that dynamically (re-)allocates the buffer to be big enough to fit the whole line
//

class FgetsObject
{
public:
    virtual char *fgets(char *s, int size) = 0;
};

char *genericReallocatingFgets(char **buffer, int *io_bufferSize, FgetsObject *getsObject)
{
    if (*io_bufferSize == 0) {
        //
        // Just pick a decent size to start out with.
        //
        *io_bufferSize = 128;
        *buffer = new char[*io_bufferSize];
    }

    int offset = 0;
    for (;;) { // loop with middle exit

        if (NULL == getsObject->fgets((*buffer) + offset, *io_bufferSize - offset)) {    // Recall that if this clips, then the next time it just keeps reading from the clip.
            //
            // EOF or error.
            //
            return NULL;
        }

        int len = (int)strlen(*buffer);
        if (len != *io_bufferSize - 1 || (*buffer)[len - 1] == '\n') {
            //
            // It fit.
            //
            return *buffer;
        }

        //
        // Double the buffer and retry.
        //

        if (((_int64)*io_bufferSize) * 2 > 0x7fffffff) {
            WriteErrorMessage("Trying to fgets() a string bigger than 2^31 bytes.  Perhaps you've supplied an incorrect input file.");
            soft_exit(1);
        }
        int newBufferSize = *io_bufferSize * 2;
        char *newBuffer = new char[newBufferSize];
        offset = len;
        memcpy(newBuffer, *buffer, offset);
        delete[] *buffer;
        *buffer = newBuffer;
        *io_bufferSize = newBufferSize;
    }

    //NOTREACHED
    return NULL;    // Just to keep the compiler happy
}

class StdioFGetsObject : public FgetsObject
{
public:
    StdioFGetsObject(FILE *_file) : file(_file) {}
    virtual char *fgets(char *s, int size) {
        return ::fgets(s, size, file);
    }
private:
    FILE *file;
};

char *reallocatingFgets(char **buffer, int *io_bufferSize, FILE *stream)
{
    StdioFGetsObject fgetsObject(stream);
    return genericReallocatingFgets(buffer, io_bufferSize, &fgetsObject);
}

class GenericFileFGetsObject : public FgetsObject
{
public:
    GenericFileFGetsObject(GenericFile *_file) : file(_file) {}
    virtual char *fgets(char *s, int size) {
        return file->gets(s, size);
    }
private:
    GenericFile *file;
};

char *reallocatingFgetsGenericFile(char **buffer, int *io_bufferSize, GenericFile *file)
{
    GenericFileFGetsObject fgetsObject(file);
    return genericReallocatingFgets(buffer, io_bufferSize, &fgetsObject);
}
