/*++

Module Name:

    GenericFile.h

Abstract:

    Generic IO class for SNAP that can read from either stdlib or HDFS.

Authors:

    Jeremy Elson, February 2014

Environment:

    User mode service.

Revision History:


--*/

#pragma once
#include "Compat.h"

class GenericFile
{
public:
	enum Mode
	{
		ReadOnly,
		WriteOnly,
	};

	static const char *HDFS_PREFIX;

	// Factory that returns either:
    //   * a GenericFile_HDFS object if the filename starts with "hdfs://"
    //   * a GenericFile_stdio object otherwise
    static GenericFile *open(const char *fileName, Mode mode);

	// Read 'count' bytes into the memory pointed at by 'ptr'.
    // Returns the actual number of bytes read, or -1 on error.
	virtual size_t read(void *ptr, size_t count) = 0;

    // Gets a single character from the stream, like stdio getc.
    // On success, returns the char promoted to an int.	
    // On failure, returns EOF.
    virtual int getchar() = 0;

	// Gets a string from the file and store it as a c string in 'str' until (num-1)
	// characters have been read or either a newline or eod-of-file is reached,
	// whichever happens first.
	virtual char *gets(char *buf, size_t count) = 0;

	// Advance forward or back by byteOffset bytes in the file.
	virtual int advance(long long byteOffset) = 0;

    // Close the file.
	virtual void close() = 0;

	// Return the name of the file.
	char *getFilename() { return _filename; }

	virtual ~GenericFile();

	virtual _int64 prefetch();	// Ignore the return value, it's just to trick the compiler into not optimizing it away.

protected:
	char *_gets_impl(char *buf, size_t count);
	GenericFile();
	Mode _mode;
	char *_filename;
};
