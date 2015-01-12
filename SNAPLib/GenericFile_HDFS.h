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
#include "Util.h"

#include "hdfs.h"

#include <vector>

class GenericFile_HDFS : public GenericFile
{
public:
	static GenericFile_HDFS *open(const char *filename, Mode mode);

	virtual size_t read(void *ptr, size_t count);
	virtual int getchar();
	virtual char *gets(char *buf, size_t count);
	virtual int advance(long long offset);
    int seek(long long offset);
	virtual void close();
	virtual ~GenericFile_HDFS();

private:
	// private constructor -- must use factory
	GenericFile_HDFS();

	// private methods and data	
	static int _initFlag;
	static int _staticInit();
	static ExclusiveLock _staticLock;
	size_t _readMultiThreaded(void *ptr, size_t count);
	size_t _readSingleThreaded(void *ptr, size_t count);

	// this is static because of an apparent bug in the HDFS library
	// that prevents clients from holding more than one handle. If you
	// open two connections to the same filesystem, then close one,
	// the other is also closed. Ugh. As a work-around we'll just keep
	// one global.
	// See: https://issues.apache.org/jira/browse/HDFS-925
	static hdfsFS _fs; 
	
	hdfsFile _file;

	class _HdfsWorkItem {
	public:
		void *ptr;
		tOffset startOffset;
		size_t count;

		_HdfsWorkItem(void *ptrArg, tOffset startOffsetArg, size_t countArg)
		{
			this->ptr = ptrArg;
			this->startOffset = startOffsetArg;
			this->count = countArg;
		}
	};

	class _HdfsWorkQueue {
	public:
		_HdfsWorkQueue(GenericFile_HDFS *gFile, void *ptr, tOffset startOffset, size_t count);
		size_t size();
		_HdfsWorkItem *getOne();
		GenericFile_HDFS *getFile() { return _gFile; }
		void createNWaiter(size_t n);
		void wait();
		void signalThreadDone();
		void signalError() { _error = true; }
		bool isErrorThrown() { return _error; }
		~_HdfsWorkQueue();

	private:
		GenericFile_HDFS *_gFile;
		ExclusiveLock _workQueueLock;
		NWaiter *_nWaiter;
		std::vector<_HdfsWorkItem *> _workItemList;
		bool _error;
	};

	static void _hdfsReaderThread(void *workQueueObject);
	static const size_t _MAX_READ_THREADS = 32;
	static const size_t _READ_CHUNK_SIZE = 16*1024*1024; // 16MB
};

