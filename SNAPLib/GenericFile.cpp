/*++

Module Name:

    GenericFile.cpp

Abstract:

    Generic IO class for SNAP that can read from either the filesystem or HDFS.

Authors:

    Jeremy Elson, February 2014

Environment:

    User mode service.

Revision History:


--*/

#include "stdafx.h"
#include <string.h>

#include "Compat.h"
#include "GenericFile.h"
#include "GenericFile_stdio.h"

#ifdef SNAP_HDFS
# include "GenericFile_HDFS.h"
#endif

const char *GenericFile::HDFS_PREFIX = "hdfs:/";

GenericFile::GenericFile()
{
	_filename = NULL;
}

GenericFile::~GenericFile()
{
	free(_filename);
}

GenericFile *GenericFile::open(const char *filename, Mode mode)
{
	if (NULL == filename) {
		return NULL;
	}

	GenericFile *retval = NULL;

	if (0 == strncmp(filename, HDFS_PREFIX, strlen(HDFS_PREFIX))) {
#ifdef SNAP_HDFS
		retval = GenericFile_HDFS::open(filename, mode);
#else
		fprintf(stderr, "SNAP not compiled with HDFS support. Set HADOOP_HOME and recompile.\n");
		retval = NULL;
#endif
	} else {
        retval = GenericFile_stdio::open(filename, mode);
	}

	if (NULL != retval) {
		retval->_filename = strdup(filename);
		retval->_mode = mode;
	}

	return retval;
}

// gets -- read until a newline. Based on the K&R implementation.
char *GenericFile::_gets_impl(char *buf, size_t count)
{
	int c;
	char *next;

	if (count == 0) {
		return NULL;
	}

	next = buf;
	while (--count > 0 && (c = getchar()) != EOF) {
		// put the input char into the current pointer position, then increment it.
		// if a newline is encountered, break
		if ((*next++ = c) == '\n')
			break;
	}

	*next = '\0';
	return (c == EOF && next == buf) ? NULL : buf;
}

	_int64
GenericFile::prefetch()
{
	const size_t ioSize = 128 * 1024 * 1024;

	char *buffer = new char[ioSize];

	for (;;) {
		if (0 == read(buffer, ioSize)) {
			delete[] buffer;
			return 0;
		}
	}
}