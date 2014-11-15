/*++

Module Name:

GenericFile_map.h

Abstract:

Generic IO class for SNAP that can map an input file.

Authors:

Bill Bolosky, September, 2014

Environment:

User mode service.

Revision History:


--*/

#pragma once

#include "GenericFile_Blob.h"
#include "Compat.h"

class GenericFile_map : public GenericFile_Blob
{
public:
	static GenericFile_map *open(const char *filename);
	virtual ~GenericFile_map();
	virtual _int64 prefetch();	// Ignore the return value, it's just to trick the compiler into not optimizing it away.
	virtual void close();

private:
	GenericFile_map(MemoryMappedFile *i_mappedFile, void *i_contents, size_t i_fileSize);

	MemoryMappedFile *mappedFile;
	const char *contents;
	size_t fileSize;
};
