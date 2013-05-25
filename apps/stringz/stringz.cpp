// stringz.cpp : Defines the entry point for the console application.
//

#include "stdafx.h"
#include "GoodRandom.h"
#include <math.h>


size_t bigStringSize = 1000000;
unsigned smallStringSize = 11;
char *bigString = NULL;
unsigned alphabetSize = 4;
double differenceProb = .001;
unsigned unitDifferenceTable[256];
char *alphabetTable = NULL;
double *probByEditDistance  = NULL;

void
GenerateBigString()
{
	bigString = new char[bigStringSize];
	alphabetTable = new char[alphabetSize];
	if (alphabetSize == 4) {
		//
		// Special case 4 to be the genetic bases.
		//
		alphabetTable[0] = 'A';
		alphabetTable[1] = 'C';
		alphabetTable[2] = 'G';
		alphabetTable[3] = 'T';
	} else {
		for (unsigned i = 0; i < alphabetSize; i++) {
			alphabetTable[i] = 'A' + i;
		}
	}

	for (size_t i = 0; i < bigStringSize; i++) {
		bigString[i] = alphabetTable[GoodFastRandom(alphabetSize-1)];
	}
}

void Hamming(const char *str1, const char *str2, size_t len, double *matchProb, unsigned *editDistance)
{
	*editDistance = 0;

	for (size_t i = 0; i < len; i++) {
		char xorValue = str1[i] ^ str2[i];				// This is 0 if they're equal, and not if they're not.
		*editDistance += unitDifferenceTable[xorValue];	// Ditto for edit distance
	}

    *matchProb = probByEditDistance[*editDistance];
}

void
usage()
{
	fprintf(stderr,"usage: stringz nStrings [smallStringSize [bigStringSize [differenceProb [alphabetSize]]]]\n");
	exit(1);
}

const unsigned max_MAPQ = 70;

struct PerMAPQ {
	size_t		count;
	size_t		errors;

	PerMAPQ() : count(0), errors(0) {}
};

struct MAPQHistogram {
	PerMAPQ		perMAPQ[max_MAPQ+1];

	void operator+=(MAPQHistogram &peer)
	{
		for (unsigned i = 0; i <= max_MAPQ; i++) {
			perMAPQ[i].count += peer.perMAPQ[i].count;
			perMAPQ[i].errors += peer.perMAPQ[i].errors;
		}
	}
};

void runTest(size_t count, MAPQHistogram *result)
{
	MAPQHistogram histogram;

	char *smallString = new char[smallStringSize] + 1;
	smallString[smallStringSize] = 0;	// Not that anything looks at this.
	unsigned oneHundredPer = (unsigned)(100/differenceProb);	// The odds of a difference are 100 per this many.  It lets us use GoodFastRandom to determine whether to introduce a difference.

	for (size_t i = 0; i < count; i++) {
		//
		// Select a random offset in the big string.
		//
		size_t bigStringOffset = (size_t)GoodFastRandom(bigStringSize - smallStringSize);

		//
		// Mutate the string
		//
		for (unsigned j = 0; j < smallStringSize; j++) {
			if (GoodFastRandom(oneHundredPer-1) < 100) {
				smallString[j] = alphabetTable[(bigString[j+bigStringOffset] + 1)%alphabetSize];
			} else {
				smallString[j] = bigString[j+bigStringOffset];
			}
		}

		size_t bestOffset = -1;
		unsigned bestEditDistance = smallStringSize + 1;
		double bestProbability = 0;
		double totalProbability = 0;
		size_t bestLocation = bigStringSize;

		for (size_t j = 0; j < bigStringSize - smallStringSize; j++) {
			double matchProb;
			unsigned editDistance;
			Hamming(smallString, bigString + j, smallStringSize, &matchProb, &editDistance);
			if (editDistance < bestEditDistance) {
				bestEditDistance = editDistance;
				bestProbability = matchProb;
				bestLocation = j;
			}
			totalProbability += matchProb;
		}

		//
		// Compute MAPQ and update the histogram
		//

		unsigned MAPQ;
		if (bestProbability == totalProbability) {
			MAPQ = max_MAPQ;
		} else {
			MAPQ = __min(max_MAPQ-1, (int)(-10 * log10(1 - bestProbability/totalProbability)));
		}

		histogram.perMAPQ[MAPQ].count++;
		if (bestLocation != bigStringOffset) {
			histogram.perMAPQ[MAPQ].errors++;
		}
	}

	*result = histogram;
}

struct ThreadContext {
	MAPQHistogram histogram;
	size_t		  count;
};

volatile int nThreadsRunning;
SingleWaiterObject threadsDone;

void WorkerThreadMain(void *param)
{
	ThreadContext *context = (ThreadContext *)param;
	runTest(context->count, &context->histogram);

	if (0 == InterlockedDecrementAndReturnNewValue(&nThreadsRunning)) {
		SignalSingleWaiterObject(&threadsDone);
	}
}


void main(int argc, char* argv[])
{
	if (argc < 2 || argc > 6) usage();

	size_t nStringz;
	if (1 != sscanf(argv[1],"%lld", &nStringz) || nStringz <= 0) {
		usage();
	}

	if (argc > 2) {
		if (1 != sscanf(argv[2], "%d", &smallStringSize) || smallStringSize <= 0) usage();
		if (argc > 3) {
			if (1 != sscanf(argv[3], "%d", &bigStringSize) || bigStringSize <= 0) usage();
			if (argc > 4) {
				if (1 != sscanf(argv[4], "%lf", &differenceProb) || differenceProb <= 0.0) usage();
				if (argc > 5) {
					if (1 != sscanf(argv[5], "%d", &alphabetSize) || alphabetSize <=0 || alphabetSize >= 26) usage();
				}
			}
		}
	}

	if (smallStringSize > bigStringSize) usage();

	for (int i = 1; i < 256; i++) {
		unitDifferenceTable[i] = 1;
	}
	unitDifferenceTable[0] = 0;

    probByEditDistance = new double[smallStringSize+1];
    for (unsigned i = 0; i <= smallStringSize; i++) {
        probByEditDistance[i] = 1.0;
        for (unsigned j = 0; j < i; j++) {
            probByEditDistance[i] *= differenceProb/(alphabetSize - 1);
        }
        for (unsigned j = i+1; j <= smallStringSize; j++) {
            probByEditDistance[i] *= (1 - differenceProb)/*/(alphabetSize - 1)*/;
        }
    }

#ifdef _MSC_VER
    if (!SetPriorityClass(GetCurrentProcess(), IDLE_PRIORITY_CLASS)) {
        fprintf(stderr,"SetPriorityClass failed, %d\n", GetLastError());
    }
#endif // _MSC_VER


	GenerateBigString();

	unsigned nThreads = GetNumberOfProcessors();
	ThreadContext *threadContexts = new ThreadContext[nThreads];

	_int64 start = timeInMillis();

	nThreadsRunning = nThreads;
	CreateSingleWaiterObject(&threadsDone);

	for (unsigned i = 0; i < nThreads - 1; i++) {
		threadContexts[i].count = nStringz / nThreads;
		StartNewThread(WorkerThreadMain, &threadContexts[i]);
	}
	threadContexts[nThreads-1].count = nStringz - (nStringz/nThreads) * (nThreads - 1);
	StartNewThread(WorkerThreadMain, &threadContexts[nThreads - 1]);

	WaitForSingleWaiterObject(&threadsDone);

	MAPQHistogram result;

	for (unsigned i = 0; i < nThreads; i++) {
		result += threadContexts[i].histogram;
	}

	for (unsigned i = 0; i <= max_MAPQ; i++) {
		printf("%d\t%lld\t%lld\n", i, result.perMAPQ[i].count, result.perMAPQ[i].errors);
	}

	_int64 stop = timeInMillis();
	printf("\nProcessed %lld stringz in %llds, %lld stringz/s\n", nStringz, (stop - start + 500) / 1000, (nStringz * 1000) / (stop - start));


}

