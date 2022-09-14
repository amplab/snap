// SelectSAMByAlignmentTime.cpp : This file contains the 'main' function. Program execution begins and ends there.
//

#include <iostream>
#include <stdio.h>
#include <windows.h>


int minTimeToExclude;
int maxTimeToExclude = INT32_MAX;
CRITICAL_SECTION CriticalSection[1];
volatile _int64 NextFileOffsetToProcess = 0;
_int64 fileSize;
HANDLE blockReadyEvent;	// Set by a worker when it adds an output block, cleared by the main thread (writer) when it doesn't have the next block ready
HANDLE readerThrottle;

const char* inputFileName;

_int64 totalLinesProcessed = 0;
_int64 totalLinesSkipped = 0;

const size_t threadChunkSize = 100 * 1024 * 1024;	// How much a thread processes at at time.
const size_t readAheadAmount = (size_t)10 * 1024 * 1024 * 1024;

//
// This is just copied from SNAPLib
//
int cheezyLogBase2(_int64 value)
{
	int retVal = 0;
	value /= 2; // Since 2^0 = 1; we'll also define cheezyLogBase2(x) = 0 where x<= 0.
	while (value > 0) {
		retVal++;
		value /= 2;
	}
	return retVal;
}

class TimingHistogram {
public:
	TimingHistogram() {}

	void addSample(_int64 time) {
		int logTime = cheezyLogBase2(time);
		int whichBucket = logTime >= nBuckets ? nBuckets - 1 : logTime;

		bucket[whichBucket].count++;
		bucket[whichBucket].totalTimeInUs += time;
	}

	void printHistogram(FILE* outfile) {
		fprintf(outfile, "time (us)\tcount\ttotal time (us)\tcumulative count\tcumulative time (us)\tpdf count\tpdf time\tcdf count\tcdf time\n");

		_int64 totalCount = 0, totalTimeInUs = 0;
		int maxUsedBucket = 0;
		for (int i = 0; i < nBuckets; i++) {
			totalCount += bucket[i].count;
			totalTimeInUs += bucket[i].totalTimeInUs;

			if (bucket[i].count > 0) {
				maxUsedBucket = i;
			}
		}

		_int64 cumulativeCount = 0, cumulativeTimeInUs = 0;
		for (int i = 0; i <= maxUsedBucket; i++) {
			cumulativeCount += bucket[i].count;
			cumulativeTimeInUs += bucket[i].totalTimeInUs;

			fprintf(outfile, "%lld\t%lld\t%lld\t%lld\t%lld\t%f\t%f\t%f\t%f\t\n",
				(_int64)1 << i, bucket[i].count, bucket[i].totalTimeInUs, cumulativeCount, cumulativeTimeInUs,
				(double)bucket[i].count / totalCount, (double)bucket[i].totalTimeInUs / totalTimeInUs,
				(double)cumulativeCount / totalCount, (double)cumulativeTimeInUs / totalTimeInUs);
		}

		fprintf(outfile, "\ntotal count %lld, total time in us %lld\n", totalCount, totalTimeInUs);
	} // printHistogram

	void mergeInto(TimingHistogram* peer) {
		for (int i = 0; i < nBuckets; i++) {
			peer->bucket[i].count += bucket[i].count;
			peer->bucket[i].totalTimeInUs += bucket[i].totalTimeInUs;
		}
	} // mergeInto


private:
	static const int nBuckets = 26;

	struct Bucket {
		Bucket() : count(0), totalTimeInUs(0) {}

		_int64		count;
		_int64		totalTimeInUs;
	};

	Bucket bucket[nBuckets];
}; // TimingHistogram

class LinearHistogram {
public:
	LinearHistogram(int minBucket_, int maxBucket_) {
		minBucket = minBucket_;
		maxBucket = maxBucket_;

		if (minBucket > maxBucket) {
			fprintf(stderr, "LinearHistogram: minBucket can't be bigger than maxBucket\n");
			exit(1);
		}

		bucket = new Bucket[maxBucket - minBucket + 1];
	} // ctor

	~LinearHistogram() {
		delete[] bucket;
	}

	void addValue(int whichBucket, _int64 timeInUs) {
		if (whichBucket < minBucket || whichBucket > maxBucket) {
			fprintf(stderr, "LinearHistogram: value outside range.\n");
			exit(1);
		}

		bucket[whichBucket - minBucket].count++;
		bucket[whichBucket - minBucket].totalTimeInUs += timeInUs;
	} // addValue

	void mergeInto(LinearHistogram* peer) {
		if (peer->minBucket != minBucket || peer->maxBucket != maxBucket) {
			fprintf(stderr, "trying to merge LinearHistograms of different shapes\n");
			exit(1);
		}

		for (int i = 0; i <= maxBucket - minBucket; i++) {
			peer->bucket[i].count += bucket[i].count;
			peer->bucket[i].totalTimeInUs += bucket[i].totalTimeInUs;
		}
	} // mergeInto

	void printHistogram(FILE* outfile) {
		fprintf(outfile, "value\tcount\ttotal time (us)\tcumulative count\tcumulative time (us)\tpdf count\tpdf time\tcdf count\tcdf time\n");

		_int64 totalCount = 0, totalTimeInUs = 0;
		int maxUsedBucket = 0;
		for (int i = 0; i <= maxBucket - minBucket; i++) {
			totalCount += bucket[i].count;
			totalTimeInUs += bucket[i - minBucket].totalTimeInUs;

			if (bucket[i].count > 0) {
				maxUsedBucket = i;
			}
		}

		_int64 cumulativeCount = 0, cumulativeTimeInUs = 0;
		for (int i = 0; i <= maxUsedBucket; i++) {
			cumulativeCount += bucket[i].count;
			cumulativeTimeInUs += bucket[i].totalTimeInUs;

			fprintf(outfile, "%d\t%lld\t%lld\t%lld\t%lld\t%f\t%f\t%f\t%f\t\n",
				i + minBucket, bucket[i].count, bucket[i].totalTimeInUs, cumulativeCount, cumulativeTimeInUs,
				(double)bucket[i].count / totalCount, (double)bucket[i].totalTimeInUs / totalTimeInUs,
				(double)cumulativeCount / totalCount, (double)cumulativeTimeInUs / totalTimeInUs);
		}

		fprintf(outfile, "\ntotal count %lld, total time in us %lld\n", totalCount, totalTimeInUs);
	} // printHistogram

private:
	int minBucket, maxBucket;

	struct Bucket {
		Bucket() : count(0), totalTimeInUs(0) {}

		_int64		count;
		_int64		totalTimeInUs;
	};

	Bucket* bucket;

}; //

TimingHistogram g_timingHistogram[1];

class OutputBlock {
public:
	char* data;
	size_t dataSize;
	_int64 fileOffset;	// In the input file.  This will be a multiple of threadChunkSize

	OutputBlock* next;
	OutputBlock* prev;

	OutputBlock()
	{
		data = NULL;
		dataSize = 0;
		fileOffset = -1;
		next = prev = this;
	}

	~OutputBlock()
	{
		delete[] data;
	}

	void AddToList(OutputBlock *listHeader) {
		//
		// The list is in fileOffset order.
		//

		_ASSERT(fileOffset > listHeader->fileOffset); // Else we could infinite loop.

		OutputBlock *blockToAddAfter = listHeader->prev;

		while (blockToAddAfter->fileOffset > fileOffset) {
			blockToAddAfter = blockToAddAfter->prev;
		}

		prev = blockToAddAfter;
		next = blockToAddAfter->next;
		prev->next = this;
		next->prev = this;
	}

	void remove() {
		prev->next = next;
		next->prev = prev;

		prev = next = NULL;
	}
};

OutputBlock BlocksReadyToWrite[1];


unsigned GetNumberOfProcessors()
{
	SYSTEM_INFO systemInfo[1];
	GetSystemInfo(systemInfo);

	return systemInfo->dwNumberOfProcessors;
}

DWORD WorkerThread(PVOID param) {
	_int64 linesProcessed = 0;
	_int64 linesSkipped = 0;
	
	int readBufferSize = threadChunkSize + 100000;
	char* readBuffer = new char[readBufferSize];

	HANDLE hInputFile = CreateFileA(inputFileName, GENERIC_READ, FILE_SHARE_READ, NULL, OPEN_EXISTING, FILE_FLAG_SEQUENTIAL_SCAN, NULL);
	if (INVALID_HANDLE_VALUE == hInputFile) {
		fprintf(stderr, "Unable to open %s, %d\n", inputFileName, GetLastError());
		exit(1);
	}

	TimingHistogram timingHistogram[1];

	for (;;) {
		_int64 readOffset = InterlockedAdd64(&NextFileOffsetToProcess, threadChunkSize) - threadChunkSize;
		if (readOffset > fileSize) {
			InterlockedAdd64(&totalLinesProcessed, linesProcessed);
			InterlockedAdd64(&totalLinesSkipped, linesSkipped);

			EnterCriticalSection(CriticalSection);
			timingHistogram->mergeInto(g_timingHistogram);
			LeaveCriticalSection(CriticalSection);
			
			return 0;
		}

		if (WAIT_OBJECT_0 != WaitForSingleObject(readerThrottle, INFINITE)) {
			fprintf(stderr, "WaitForSingleObject on throttle failed, %d\n", GetLastError());
			exit(1);
		}

		LARGE_INTEGER liFileOffset;
		liFileOffset.QuadPart = readOffset;

		if (!SetFilePointerEx(hInputFile, liFileOffset, NULL, SEEK_SET)) {
			fprintf(stderr, "Unable to SetFilePointer, %d\n", GetLastError());
			exit(1);
		}

		DWORD bytesRead;
		if (!ReadFile(hInputFile, readBuffer, readBufferSize, &bytesRead, NULL)) {
			fprintf(stderr, "ReadFile failed, %d\n", GetLastError());
			exit(1);
		}

		if (bytesRead != readBufferSize && bytesRead != (DWORD)(fileSize - readOffset)) {
			fprintf(stderr, "Read unexpected number of bytes %d != %d at offset %lld\n", bytesRead, readBufferSize, readOffset);
			exit(1);
		}

		_int64 stopOffset = readOffset + threadChunkSize;	// Run to the first \n after this offset or EOF

		OutputBlock* block = new OutputBlock;
		block->fileOffset = readOffset;
		block->data = new char[threadChunkSize + 10000];	// Won't necessarily use it all, but this is a bound (and really, a buffer overflow, but this should never have hostile input).

		const char* line = readBuffer;
		if (readOffset != 0) {
			while (line < readBuffer + bytesRead && *line != '\n') {
				line++;
			}

			line++;	// Consume the newline

			if (line >= readBuffer + bytesRead) {
				//
				// The last chunk didn't include a whole line.  We're done.
				//
				delete block;

				InterlockedAdd64(&totalLinesProcessed, linesProcessed);
				InterlockedAdd64(&totalLinesSkipped, linesSkipped);

				EnterCriticalSection(CriticalSection);
				timingHistogram->mergeInto(g_timingHistogram);
				LeaveCriticalSection(CriticalSection);

				return 0;
			}
		}

		while (line <= readBuffer + threadChunkSize && line < readBuffer + bytesRead) {
			const char* endOfLine = line;
			while (endOfLine < readBuffer + bytesRead && *endOfLine != '\n') {	// Don't check for end of chunk, we always consume the first newline in the next chunk.
				endOfLine++;
			}
			endOfLine++;	// So it points at the next line

			if (*line == '@') {
				//
				// A header line.  Copy it identically.
				//
				memcpy(block->data + block->dataSize, line, endOfLine - line);
				block->dataSize += endOfLine - line;
			} else {
				//
				// Look for the AT tag, "\tAT:i:"
				//
				linesProcessed++; // don't count headers

				bool foundAT = false;
				for (const char* current = line; current < endOfLine - 7; current++) {	// -7 leaves space for \tAT:i: at least one digit and then something after (tab or newline)
					if (current[0] == '\t' && current[1] == 'A' && current[2] == 'T' && current[3] == ':' && current[4] == 'i' && current[5] == ':') {
						int alignmentTime = atoi(current + 6);

						if (alignmentTime < minTimeToExclude || alignmentTime > maxTimeToExclude) {
							memcpy(block->data + block->dataSize, line, endOfLine - line);
							block->dataSize += endOfLine - line;
						} else {
							linesSkipped++;
						}

						timingHistogram->addSample(alignmentTime);

						foundAT = true;
						break;	// out of the loop looking for AT
					} // if we're at AT:i:
				} // walking the line

				if (!foundAT) {
					const int bytesToPrint = 500;
					char buffer[bytesToPrint+1];
					memcpy(buffer, line, bytesToPrint);
					buffer[bytesToPrint] = '\0';
					fprintf(stderr, "Found line without an AT tag, first bit: %s\n", buffer);
				}


			} // not a header line

			line = endOfLine;
		} // while we have data in this chunk

		if (block->dataSize > 0) {
			EnterCriticalSection(CriticalSection);
			block->AddToList(BlocksReadyToWrite);
			SetEvent(blockReadyEvent);
			if (readOffset > BlocksReadyToWrite->next->fileOffset + readAheadAmount) {
				ResetEvent(readerThrottle);
			}
			LeaveCriticalSection(CriticalSection);
		} else {
			delete block;
		}
	} // for ever

	return 0;	// NOTREACHED but needed for the compiler to be happy
}

void usage()
{
	fprintf(stderr, "usage: SelectSAMByAlignmentTime inputFile outputFile minTimeToExclude {maxTimeToExclude}\n");
	exit(1);
} // usage

int main(int argc, char **argv)
{
	if (argc != 4 && argc != 5) {
		usage();
	}

	minTimeToExclude = atoi(argv[3]);
	if (argc >= 5) {
		maxTimeToExclude = atoi(argv[4]);
	}

	if (minTimeToExclude < 0 || minTimeToExclude >= maxTimeToExclude) {
		fprintf(stderr, "min time to exclude can't be negative and must be less than max time to exclude\n");
		exit(1);
	}

	inputFileName = argv[1];

	HANDLE hInputFile = CreateFileA(inputFileName, GENERIC_READ, FILE_SHARE_READ, NULL, OPEN_EXISTING, FILE_FLAG_SEQUENTIAL_SCAN, NULL);
	if (INVALID_HANDLE_VALUE == hInputFile) {
		fprintf(stderr, "Unable to open %s, %d\n", inputFileName, GetLastError());
		exit(1);
	}

	LARGE_INTEGER liFileSize;
	if (!GetFileSizeEx(hInputFile, &liFileSize)) {
		fprintf(stderr, "Unable to GetFileSizeEx(), %d\n", GetLastError());
		exit(1);
	}

	CloseHandle(hInputFile);

	fileSize = liFileSize.QuadPart;

	InitializeCriticalSection(CriticalSection);

	blockReadyEvent = CreateEvent(NULL, TRUE, FALSE, NULL);
	if (NULL == blockReadyEvent) {
		fprintf(stderr, "Unable to CreateEvent, %d\n", GetLastError());
		exit(1);
	}

	readerThrottle = CreateEvent(NULL, TRUE, TRUE, NULL);
	if (NULL == blockReadyEvent) {
		fprintf(stderr, "Unable to CreateEvent, %d\n", GetLastError());
		exit(1);
	}


	HANDLE hOutputFile = CreateFileA(argv[2], GENERIC_WRITE, FILE_SHARE_READ, NULL, CREATE_ALWAYS, FILE_ATTRIBUTE_NORMAL, NULL);
	if (INVALID_HANDLE_VALUE == hOutputFile) {
		fprintf(stderr, "Unable to open output file %s, %d\n", argv[2], GetLastError());
		exit(1);
	}

	unsigned nProcessors = GetNumberOfProcessors();

	for (int i = 0; i < nProcessors; i++) {
		if (!CreateThread(NULL, 0, WorkerThread, NULL, 0, NULL)) {
			fprintf(stderr, "CreateThread failed, %d\n", GetLastError());
			exit(1);
		}
	}

	_int64 nextBlockOffset = 0;

	_int64 startTime = GetTickCount64();

	printf("Processing file (1 dot/GB): ");
	_int64 lastDot = 0;

	while (nextBlockOffset < fileSize) {
		DWORD waitResult = WaitForSingleObject(blockReadyEvent, INFINITE);
		if (waitResult != WAIT_OBJECT_0) {
			fprintf(stderr, "WaitForSingleObject failed, %d\n", GetLastError());
			exit(1);
		}

		for (;;) {
			OutputBlock* blockToWrite = NULL;

			EnterCriticalSection(CriticalSection);
			if (BlocksReadyToWrite->next->fileOffset == nextBlockOffset) {
				blockToWrite = BlocksReadyToWrite->next;
				blockToWrite->remove();
			}

			if (blockToWrite == NULL) {
				ResetEvent(blockReadyEvent);
			} else if (blockToWrite->fileOffset + readAheadAmount > NextFileOffsetToProcess) {
				//
				// We're not too far behind.  Let the readers go.
				//
				SetEvent(readerThrottle);
			}
			LeaveCriticalSection(CriticalSection);

			if (blockToWrite == NULL) {
				break;	// Go back to sleep waiting for work
			}

			DWORD nBytesWritten;
			if (!WriteFile(hOutputFile, blockToWrite->data, blockToWrite->dataSize, &nBytesWritten, NULL)) {
				fprintf(stderr, "WriteFile failed, %d\n", GetLastError());
				exit(1);
			}

			nextBlockOffset += threadChunkSize;
			delete blockToWrite;

			if (nextBlockOffset >= lastDot + 1024 * 1024 * 1024) {
				printf(".");
				lastDot += 1024 * 1024 * 1024;
			}
		} // for ever (write run without sleeping)
	} // while we have work

	printf("\nSkipped %lld of %lld reads in %llds\n", totalLinesSkipped, totalLinesProcessed, (GetTickCount64() - startTime) / 1000);

	g_timingHistogram->printHistogram(stdout);

	CloseHandle(hOutputFile);

} // main

