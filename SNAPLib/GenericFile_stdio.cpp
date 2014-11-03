/*++

Module Name:

    GenericFile_stdio.cpp

Abstract:

    Generic IO class for SNAP that can read from stdio.

Authors:

    Jeremy Elson, February 2014

Environment:

    User mode service.

Revision History:


--*/

#include "stdafx.h"
#include "Compat.h"
#include "GenericFile_stdio.h"
#include "Error.h"

GenericFile_stdio::GenericFile_stdio()
{
	_file = NULL;
}

GenericFile_stdio::~GenericFile_stdio()
{
}

GenericFile_stdio *GenericFile_stdio::open(const char *filename, Mode mode)
{
	GenericFile_stdio *retval = new GenericFile_stdio();
	retval->_mode = mode;

	const char *fMode = NULL;

	switch (mode) {
	case ReadOnly:
		fMode = "rb";
		break;
	case WriteOnly:
		fMode = "wb";
		break;
	}

	retval->_file = fopen(filename, fMode);

	if (retval->_file == NULL) {
		delete retval;
		return NULL;
	}

	return retval;
}

GenericFile_stdio *GenericFile_stdio::open(const char *filename)
{
        return open(filename, ReadOnly);
}

size_t GenericFile_stdio::read(void *ptr, size_t count)
{

	return fread(ptr, 1, count, _file);

}

int GenericFile_stdio::getchar()
{
	return fgetc(_file);
}

char *GenericFile_stdio::gets(char *buf, size_t count)
{
	return fgets(buf, (int) count, _file);
}

int GenericFile_stdio::advance(long long offset)
{
	return _fseek64bit(_file, offset, SEEK_CUR);
}
 
void GenericFile_stdio::close()
{
	fclose(_file);
}
