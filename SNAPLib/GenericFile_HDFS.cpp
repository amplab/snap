/*++

Module Name:

    GenericFile_HDFS.cpp

Abstract:

    Generic IO class for SNAP that can read from HDFS.

Authors:

    Jeremy Elson, February 2014

Environment:

    User mode service.

Revision History:


--*/

#include "stdafx.h"

#ifdef SNAP_HDFS

#include "Compat.h"
#include "GenericFile.h"
#include "GenericFile_HDFS.h"

GenericFile_HDFS::GenericFile_HDFS()
{
	_fs = NULL;
	_file = NULL;
}

GenericFile_HDFS::~GenericFile_HDFS()
{
}


GenericFile_HDFS *GenericFile_HDFS::open(const char *filename, Mode mode)
{
	GenericFile_HDFS *retval = new GenericFile_HDFS();

#ifdef _MSC_VER
	staticLibInit();
#endif

	retval->_fs = hdfsConnect("default", 0);

	if (NULL == retval->_fs) {
		fprintf(stderr, "can't open HDFS");
		goto fail;
	}

	switch (mode) {
	case ReadOnly:
		retval->_file = hdfsOpenFile(retval->_fs, filename, O_RDONLY, 0, 0, 0);
		break;
	case WriteOnly:
		retval->_file = hdfsOpenFile(retval->_fs, filename, O_WRONLY | O_CREAT, 0, 0, 0);
		break;
	default:
		fprintf(stderr, "GenericFile_HDFS::open::unknown file mode");
		break;
	}

	if (0 == retval->_file) {
		fprintf(stderr, "couldn't open hdfs file");
		goto fail;
	}

	return retval;

fail:
	if (NULL != retval->_fs) {
		hdfsDisconnect(retval->_fs);
	}

	delete retval;
	return NULL;
	
}

size_t GenericFile_HDFS::read(void *ptr, size_t count)
{
	size_t totalRead = 0;

	while (true) {
		// HDFS takes read arguments as signed 32-bit ints, so do smaller
		// reads if we got a larger argument
		size_t desiredRead = count - totalRead;
		tSize readSize;

		if (desiredRead > INT32_MAX) {
			readSize = INT32_MAX;
		} else {
			readSize = (tSize) desiredRead;
		}
			
		tSize retval = hdfsRead(_fs, _file, ((char *) ptr) + totalRead, readSize);

		if (retval < 0) {
			perror("hdfsRead");
			return totalRead;
		} else {
			totalRead += retval;

			if (retval == 0 || totalRead >= count) {
				return totalRead;
			}
		}
	}
}

// ridiculously slow implementaiton of getChar.
int GenericFile_HDFS::getchar()
{
	char buf[1];

	if (1 == read(buf, sizeof(buf))) {
		return buf[0];
	} else {
		return EOF;
	}
}

// gets -- read until a newline. Based on the K&R implementation, but
// a zillion times slower because we're going all the way out to the JVM
// for each character. We can buffer locally if perf hurts too much.
char *GenericFile_HDFS::gets(char *buf, size_t count)
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

int GenericFile_HDFS::advance(long offset)
{
	if (offset == 0)
		return 0;

	tOffset currOffset = hdfsTell(_fs, _file);
	return hdfsSeek(_fs, _file, currOffset + offset);
}

void GenericFile_HDFS::close()
{
	if (_mode == GenericFile::WriteOnly) {
		if (hdfsFlush(_fs, _file)) {
			fprintf(stderr, "Failed to flush %s!\n", _filename);
		}
	}

	hdfsCloseFile(_fs, _file);
}


#endif // SNAP_HDFS
