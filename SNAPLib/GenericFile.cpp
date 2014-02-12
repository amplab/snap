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

#include "GenericFile.h"
#include "GenericFile_HDFS.h"
#include "GenericFile_stdio.h"

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
		retval = GenericFile_HDFS::open(filename, mode);
	} else {
		retval = GenericFile_stdio::open(filename, mode);
	}

	if (NULL != retval) {
		retval->_filename = _strdup(filename);
		retval->_mode = mode;
	}

	return retval;
}



