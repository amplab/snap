/*++

Module Name:

    GenericFile_stdio.h

Abstract:

    Generic IO class for SNAP that can read from stdio.

Authors:

    Jeremy Elson, February 2014

Environment:

    User mode service.

Revision History:


--*/

#pragma once

#include "GenericFile.h"

class GenericFile_stdio : public GenericFile
{
public:
	static GenericFile_stdio *open(const char *filename, Mode mode);
	static GenericFile_stdio *open(const char *filename); // no Mode means ReadOnly
	virtual size_t read(void *ptr, size_t count);
	virtual int getchar();
	virtual char *gets(char *buf, size_t count);
	virtual int advance(long long offset);
	virtual ~GenericFile_stdio();
	virtual void close();

private:
	GenericFile_stdio();

	FILE *_file;
};
