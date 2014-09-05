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
#include "Error.h"
#include "GenericFile.h"
#include "GenericFile_HDFS.h"
#include "Util.h"

ExclusiveLock GenericFile_HDFS::_staticLock;
hdfsFS GenericFile_HDFS::_fs;

// this goes last so staticInit runs after others are initialized;
// it has side-effects that initialize the other members
int GenericFile_HDFS::_initFlag = GenericFile_HDFS::_staticInit(); 

int GenericFile_HDFS::_staticInit()
{
	InitializeExclusiveLock(&_staticLock);
	SetExclusiveLockWholeProgramScope(&_staticLock);

	_fs = NULL;
	srand(unsigned(time(NULL)));

	return 1;
}

GenericFile_HDFS::GenericFile_HDFS()
{
	_file = NULL;
}

GenericFile_HDFS::~GenericFile_HDFS()
{
}

GenericFile_HDFS *GenericFile_HDFS::open(const char *filename, Mode mode)
{
	GenericFile_HDFS *retval = new GenericFile_HDFS();
	bool noFs = false;

	AcquireExclusiveLock(&_staticLock);

	if (NULL == _fs) {
		// note, because we're now holding this globally due to HDFS' ridiculous
		// bug (see:  https://issues.apache.org/jira/browse/HDFS-925), this never gets closed.
		// sorry. feel free to refcount it and hdfsDisconnect() it in ::close(). --JE
		_fs = hdfsConnect("default", 0);

		if (NULL == _fs) {
			fprintf(stderr, "can't open HDFS");
			noFs = true;
		}
	}

	ReleaseExclusiveLock(&_staticLock);

	if (noFs) {
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
	delete retval;
	return NULL;
	
}

// for debugging; ignore
static char *startingBytes(char *ptr, size_t totalRead)
{
	static char buf[100];

	for (size_t i = 0; i < min(totalRead, 10); i++)
		sprintf(buf+i*3, "%02x:", (unsigned char) ptr[i]);

	return buf;
}

size_t GenericFile_HDFS::_readSingleThreaded(void *ptr, size_t count)
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

		// WriteErrorMessage("reading 0x%x from 0x%llx --> %x\n", readSize, hdfsTell(_fs, _file), (char *) ptr+totalRead);

		tSize retval = hdfsRead(_fs, _file, ((char *) ptr) + totalRead, readSize);

		if (retval < 0) {
			perror("hdfsRead");
			return totalRead;
		} else {
			totalRead += retval;

			if (retval == 0 || totalRead >= count) {
				// WriteErrorMessage("...read at 0x%x, %ld bytes, starts with %s\n", ptr, totalRead, startingBytes((char *)ptr, totalRead));
				return totalRead;
			}
		}
	}
}


/////////////////
/// Implementation of multi-threaded work queue for HDFS reads
////////////////

GenericFile_HDFS::_HdfsWorkQueue::_HdfsWorkQueue(GenericFile_HDFS *gFile, void *void_ptr, tOffset startOffset, size_t count)
{
	InitializeExclusiveLock(&_workQueueLock);
	_error = false;

	this->_gFile = gFile;

	// cast to char ptr necessary to do pointer math
	char *ptr = (char *) void_ptr;

	while (count > 0) {
		long readSize = count >= _READ_CHUNK_SIZE ? _READ_CHUNK_SIZE : (long) count;

		_workItemList.push_back(new _HdfsWorkItem(ptr, startOffset, readSize));

		ptr += readSize;
		startOffset += readSize;
		count -= readSize;
	};

	std::random_shuffle(_workItemList.begin(), _workItemList.end());
}

size_t GenericFile_HDFS::_HdfsWorkQueue::size()
{
	AcquireExclusiveLock(&_workQueueLock);
	size_t retval = _workItemList.size();
	ReleaseExclusiveLock(&_workQueueLock);
	return retval;
}

GenericFile_HDFS::_HdfsWorkItem *GenericFile_HDFS::_HdfsWorkQueue::getOne()
{
	_HdfsWorkItem *retval = NULL;

	// If an error has been thrown, don't bother with the rest of the work
	if (isErrorThrown()) {
		return NULL;
	}

	AcquireExclusiveLock(&_workQueueLock);

	if (_workItemList.empty()) {
		retval = NULL;
	} else {		
		retval = _workItemList.back();
		_workItemList.pop_back();
	}

	ReleaseExclusiveLock(&_workQueueLock);
	return retval;
}

void GenericFile_HDFS::_HdfsWorkQueue::createNWaiter(size_t n)
{
	_nWaiter = new NWaiter(n);
}

void GenericFile_HDFS::_HdfsWorkQueue::signalThreadDone()
{
	_nWaiter->signal();
}

void GenericFile_HDFS::_HdfsWorkQueue::wait()
{
	_nWaiter->wait();
}

GenericFile_HDFS::_HdfsWorkQueue::~_HdfsWorkQueue()
{
	if (!_workItemList.empty()) {
		WriteErrorMessage("HDFS work queue not empty when destroyed! This is a bug. A grave one.");
		exit(1);
	}

	DestroyExclusiveLock(&_workQueueLock);
	delete _nWaiter;
}


size_t GenericFile_HDFS::_readMultiThreaded(void *ptr, size_t count)
{
	_HdfsWorkQueue *workQueue = new _HdfsWorkQueue(this, ptr, hdfsTell(_fs, _file), count);
	size_t numThreads = workQueue->size();

	if (numThreads > _MAX_READ_THREADS) {
		numThreads = _MAX_READ_THREADS;
	}

	// WriteErrorMessage("Reading %lld bytes with %d threads (%lld work queue items)\n", count, numThreads, workQueue->size());

	// Tell the work queue how many threads we're creating
	workQueue->createNWaiter(numThreads);

	for (size_t i = 0; i < numThreads; i++) {
		StartNewThread(_hdfsReaderThread, workQueue);
	}

	// Wait for all threads to finish
	workQueue->wait();

	size_t retval;

	if (workQueue->isErrorThrown()) {
		retval = 0;
	} else {
		retval = count;
	}

	delete workQueue;
	return retval;
}

void GenericFile_HDFS::_hdfsReaderThread(void *workQueueVoid)
{
	_HdfsWorkQueue *workQueue = static_cast<_HdfsWorkQueue*>(workQueueVoid);

	if (NULL == workQueue) {
		WriteErrorMessage("HDFS reader thread didn't get a valid work queue!\n");
		exit(-1); // can't even signal the thread is done
	}

	// we need a separate file structure so that we can seek independently
	GenericFile_HDFS *localFile = GenericFile_HDFS::open(workQueue->getFile()->getFilename(), ReadOnly);

	if (NULL == localFile) {
		WriteErrorMessage("HDFS reader thread could not open file!\n");
		workQueue->signalError();
		goto done;
	} 

	while (true) {
		_HdfsWorkItem *nextItem = workQueue->getOne();
		
		if (NULL == nextItem)
			goto done;
		
		localFile->seek(nextItem->startOffset);

		// WriteErrorMessage("readerThread: reading 0x%x from %x --> %x\n", nextItem->count, nextItem->startOffset, nextItem->ptr);

		size_t retval = localFile->_readSingleThreaded(nextItem->ptr, nextItem->count);
		
		if (retval != nextItem->count) {
			WriteErrorMessage("HDFS read error: starting at offset %lld, tried to read %lld, only got %lld\n",
				nextItem->startOffset, nextItem->count, retval);
			workQueue->signalError();
		}
	}

done:
	if (NULL != localFile) {
		localFile->close();
		delete localFile;
	}

	workQueue->signalThreadDone();
}


size_t GenericFile_HDFS::read(void *ptr, size_t count)
{
	if (count <= _READ_CHUNK_SIZE || _MAX_READ_THREADS == 0) {
		return _readSingleThreaded(ptr, count);
	} else {
		size_t retval = _readMultiThreaded(ptr, count);

		// move the file pointer forward on the "main" file,
		// because the multi threaded reader will open its own file pointers
		advance(count);

		return retval;
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

// very slow because we're going all the way out to the JVM
// for each character. We can buffer locally if perf hurts too much.
char *GenericFile_HDFS::gets(char *buf, size_t count)
{
	return _gets_impl(buf, count);
}

int GenericFile_HDFS::advance(long long offset)
{
	if (offset == 0)
		return 0;

	tOffset currOffset = hdfsTell(_fs, _file);
	return hdfsSeek(_fs, _file, currOffset + offset);
}

int GenericFile_HDFS::seek(long long offset)
{
	return hdfsSeek(_fs, _file, offset);
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
