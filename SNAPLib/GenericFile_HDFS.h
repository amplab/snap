/*++

Module Name:

    GenericFile_HDFS.h

Abstract:

    Generic IO class for SNAP that can read from HDFS.

Authors:

    Jeremy Elson, February 2014

Environment:

    User mode service.

Revision History:


--*/

#pragma once

#include "GenericFile.h"

#include "pdclibhdfs/inc/hdfs.h"

class GenericFile_HDFS : public GenericFile
{
public:
	static GenericFile_HDFS *open(const char *filename, Mode mode);

	virtual size_t read(void *ptr, size_t count);
	virtual char *gets(char *buf, size_t count);
	virtual int advance(long offset);
	virtual void close();
	virtual ~GenericFile_HDFS();

	int getchar();

private:
	GenericFile_HDFS();

	hdfsFS _fs;
	hdfsFile _file;
};
