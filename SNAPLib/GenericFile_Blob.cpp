/*++

Module Name:

    GenericFile_Blob.cpp

Abstract:

    Generic IO class for SNAP that can read from an in-memory blob.

Authors:

    Bill Bolosky, March, 2014

Environment:

    User mode service.

Revision History:


--*/

#include "stdafx.h"
#include "GenericFile_Blob.h"

GenericFile_Blob::GenericFile_Blob(void *i_blob, size_t i_blobSize) : blob((char *)i_blob), readPointer((char *)i_blob), blobSize(i_blobSize), blobEnd((char *)i_blob + i_blobSize)
{
}

GenericFile_Blob *
GenericFile_Blob::open(void *i_blob, size_t i_blobSize)
{
	GenericFile_Blob *file = new GenericFile_Blob(i_blob, i_blobSize);
 
    return file;
}

size_t
GenericFile_Blob::read(void *ptr, size_t count)
{
    size_t bytesReturned;
    void *base = mapAndAdvance(count, &bytesReturned);

    memcpy(ptr, base, bytesReturned);
    return bytesReturned;
}

int
GenericFile_Blob::advance(long long offset)
{
    if (offset < 0) {
        if (readPointer - blob > -1 * offset) {
            return 1;
        }
        readPointer += offset;
        return 0;
    }

    long long amountRemaining = blobEnd - readPointer;
    if (amountRemaining < offset) {
        return EOF;
    }
    readPointer += offset;

    return 0;
}


	void 
GenericFile_Blob::close()
{
    blob = readPointer = blobEnd = NULL;
    blobSize = 0;
}

GenericFile_Blob::~GenericFile_Blob()
{
    close();
}

	int 
GenericFile_Blob::getchar()
{
    if (readPointer >= blobEnd) {
        return EOF;
    }

    unsigned char c = *(unsigned char *)readPointer;
    readPointer++;

    return (int)c;
}

char *GenericFile_Blob::gets(char *buf, size_t count)
{
	return _gets_impl(buf, count);
}

    //
    // This gives a pointer into the blob rather than copying it.
    // It's the caller's responsibility to assure that the blob doesn't
    // go away while this pointer's still in use.
    //
    void *
GenericFile_Blob::mapAndAdvance(size_t count, size_t *bytesReturned)
{
    size_t amountRemaining = blobEnd - readPointer;

    if (count > amountRemaining) {
        *bytesReturned = amountRemaining;
    } else {
        *bytesReturned = count;
    }

    void *retVal = readPointer;
    readPointer += *bytesReturned;

    return retVal;
}
    size_t
GenericFile_Blob::getAmountUsed()
{
    return readPointer - blob;
}
