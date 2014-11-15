/*++

Module Name:

GenericFile_map.cpp

Abstract:

Generic IO class for SNAP that can map an input file.

Authors:

Bill Bolosky, September, 2014

Environment:

User mode service.

Revision History:


--*/

#include "stdafx.h"
#include "GenericFile_map.h"
#include "Error.h"
#include "exit.h"

GenericFile_map *GenericFile_map::open(const char *filename)
{
	size_t fileSize = QueryFileSize(filename);
	void *contents;
	MemoryMappedFile *mappedFile = OpenMemoryMappedFile(filename, 0, fileSize, &contents);

	return new GenericFile_map(mappedFile, contents, fileSize);
}

GenericFile_map::GenericFile_map(MemoryMappedFile *i_mappedFile, void *i_contents, size_t i_fileSize) : mappedFile(i_mappedFile), contents((const char *)i_contents), fileSize(i_fileSize), GenericFile_Blob(i_contents, i_fileSize)
{
}

	void
GenericFile_map::close()
{
	if (NULL != mappedFile) {
		CloseMemoryMappedFile(mappedFile);
		mappedFile = NULL;

		GenericFile_Blob::close();
	}
}

GenericFile_map::~GenericFile_map()
{
	close();
}

	_int64
GenericFile_map::prefetch()
{
    AdviseMemoryMappedFilePrefetch(mappedFile);
	int pageSize = getpagesize();

	_int64 total = 0;

	for (size_t offset = 0; offset < fileSize / sizeof (_int64); offset += 4 ) {
		total += ((_int64 *)contents)[offset];
	}

	return total;		// We're returning this just to keep the compiler from optimizing away the whole thing.
}
