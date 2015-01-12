/*++

Module Name:

    GenericFile_Blob.h

Abstract:

    Generic IO class for SNAP that can read from an in-memory blob.

Authors:

    Bill Bolosky, March, 2014

Environment:

    User mode service.

Revision History:


--*/

#pragma once

#include "GenericFile.h"

class GenericFile_Blob: public GenericFile
{
public:
    //
    // This object does not take ownership of the blob, so it's the caller's
    // responsibility to assure that it continues to exist as long as the
    // GenericFile_Blob object does, and to free the memory as necessary.
    // In addition, if anyone calls mapAndAdvance, the caller must
    // assure that the blob continues to exist until all uses of that
    // pointer are finished as well.
    //
	static GenericFile_Blob *open(void *i_blob, size_t i_blobSize);

	virtual size_t read(void *ptr, size_t count);
	virtual int getchar();
	virtual char *gets(char *buf, size_t count);
	virtual int advance(long long offset);
	virtual void close();
	virtual ~GenericFile_Blob();
    virtual size_t getAmountUsed();



    //
    // This gives a pointer into the blob rather than copying it.
    // It's the caller's responsibility to assure that the blob doesn't
    // go away while this pointer's still in use.
    //
    void *mapAndAdvance(size_t count, size_t *bytesReturned);

protected:
	GenericFile_Blob(void *i_blob, size_t i_blobSize);

private:
	char    *blob;
    char    *blobEnd;
    char    *readPointer;
    size_t   blobSize;
};
