/*++

Module Name:

KMerAligner.cpp

Abstract:

Functions for running the single end aligner sub-program.

Authors:

Tracy Ballinger & Bill Bolosky, July, 2014

Environment:
`
User mode service.

Revision History:

Adapted from cSNAP, which was in turn adapted from the scala prototype

--*/

#include "stdafx.h"
#include "options.h"
#include "BaseAligner.h"
#include "Compat.h"
#include "RangeSplitter.h"
#include "GenomeIndex.h"
#include "Range.h"
#include "SAM.h"
#include "Tables.h"
#include "AlignerContext.h"
#include "AlignerOptions.h"
#include "FASTQ.h"
#include "Util.h"
#include "SingleAligner.h"
#include "MultiInputReadSupplier.h"
#include "KMerAligner.h"

using namespace std;
using util::stringEndsWith;



KMerAlignerContext::KMerAlignerContext(AlignerExtension* i_extension)
: extension(this), AlignerContext(0, NULL, NULL, &extension)
{
	KMerAlignerExtension *extension = new KMerAlignerExtension(this);
}

AlignerStats*
KMerAlignerContext::newStats()
{
	return new AlignerStats();
}

void
KMerAlignerContext::runTask()
{
	ParallelTask<KMerAlignerContext> task(this);
	task.run();
}


struct SeedAndInputFile {
	Seed		seed;
	unsigned	whichInputFile;
};

void GetHashTableAndKeyFromSeed(unsigned seedLength, Seed seed, unsigned hashTableBits, _uint64 *o_key, unsigned *o_whichHashTable)
{
	*o_whichHashTable = (unsigned)(seed.getBases() >> (seedLength * 2 - hashTableBits));
	*o_key = seed.getBases() & ((((_uint64)1) << (seedLength * 2 - hashTableBits))-1);
}

class PendingInserts {
public:
	PendingInserts(SNAPHashTable *i_hashTable, ExclusiveLock *i_hashTableLock, unsigned i_nInputFiles, unsigned nHashTables, unsigned i_seedLength) : 
		hashTable(i_hashTable), hashTableLock(i_hashTableLock), usedCount(0), nInputFiles(i_nInputFiles), seedLength(i_seedLength)
	{
		hashTableBits = cheezyLogBase2(nHashTables);

		inserts = new SeedAndInputFile[batchSize];
		insertArray = new SNAPHashTable::ValueType[nInputFiles];
		for (unsigned i = 0; i < nInputFiles; i++) {
			insertArray[i] = 0;
		}
	}

	~PendingInserts() {
		delete[] inserts;
	}

	void recordSeed(Seed seed, unsigned whichInputFile)
	{
		_ASSERT(usedCount < batchSize);
		inserts[usedCount].seed = seed;
		inserts[usedCount].whichInputFile = whichInputFile;
		usedCount++;
		if (usedCount == batchSize) {
			flush();
		}
	}

	void flush()
	{
		AcquireExclusiveLock(hashTableLock);
		for (unsigned i = 0; i < usedCount; i++) {
			unsigned short *hashTableValues;
			_uint64 key;
			unsigned whichHashTable;
			GetHashTableAndKeyFromSeed(seedLength, inserts[i].seed, hashTableBits, &key, &whichHashTable);
			if (NULL == (hashTableValues = (unsigned short *)hashTable->GetFirstValueForKey(key))) {
				insertArray[inserts[i].whichInputFile] = 1;
				if (!hashTable->Insert(key, insertArray)) {
					WriteErrorMessage("Hash table full\n");
					soft_exit(1);
				}
				insertArray[inserts[i].whichInputFile] = 0;
			} else {
				if (hashTableValues[inserts[i].whichInputFile] < 0xfffe) {
					hashTableValues[inserts[i].whichInputFile]++;
				}
			}
		}
		ReleaseExclusiveLock(hashTableLock);

		usedCount = 0;
	}

private:
	static const unsigned batchSize;

	unsigned						usedCount;
	unsigned						nInputFiles;
	SeedAndInputFile	*			inserts;
	SNAPHashTable		*			hashTable;
	ExclusiveLock		*			hashTableLock;
	SNAPHashTable::ValueType	*	insertArray;
	unsigned						hashTableBits;
	unsigned						seedLength;
};

const unsigned PendingInserts::batchSize = 1000;


void
KMerAlignerContext::runIterationThread()
{
	PreventMachineHibernationWhileThisThreadIsAlive();

	unsigned seedLength = index->getSeedLength();
	unsigned hashTableBits = cheezyLogBase2(nHashTables);

	PendingInserts **pendingInserts = new PendingInserts *[nHashTables];
	for (unsigned i = 0; i < nHashTables; i++) {
		pendingInserts[i] = new PendingInserts(hashTables[i], &hashTableLocks[i], options->nInputs, nHashTables, seedLength);
	}

	for (unsigned whichInputFile = 0; whichInputFile < options->nInputs; whichInputFile++) {
		ReadSupplier * readSupplier = readSupplierGenerators[whichInputFile]->generateNewReadSupplier();
		if (NULL == readSupplier) {
			continue;
		}

		Read *read;
		while (NULL != (read = readSupplier->getNextRead())) {
			stats->totalReads++;

			if (read->getDataLength() < seedLength) {
				continue;
			}

			stats->usefulReads++;

			for (unsigned nextSeedToTest = 0; nextSeedToTest < read->getDataLength() - seedLength; nextSeedToTest++) {
#if 0
				const unsigned prefetchDepth = 0;
				if (doAlignerPrefetch && nextSeedToTest + prefetchDepth < read->getDataLength() - seedLength && Seed::DoesTextRepresentASeed(read->getData() + nextSeedToTest + prefetchDepth, seedLength)) {
					Seed seed(read->getData() + nextSeedToTest + prefetchDepth, seedLength);
					index->prefetchLookup(seed);
				}

#endif // 0
				if (!Seed::DoesTextRepresentASeed(read->getData() + nextSeedToTest, seedLength)) {
					continue;
				}

				Seed seed(read->getData() + nextSeedToTest, seedLength);

				_uint64 key;
				unsigned whichHashTable;
				GetHashTableAndKeyFromSeed(seedLength, seed, hashTableBits, &key, &whichHashTable);

				if (index->doesGenomeIndexHave64BitLocations()) {
					_int64 nHits[NUM_DIRECTIONS];
					const GenomeLocation *hits[NUM_DIRECTIONS];
					GenomeLocation hitBuffer;

					index->lookupSeed(seed, &nHits[FORWARD], &hits[FORWARD], &nHits[RC], &hits[RC], &hitBuffer, &hitBuffer);

					for (Direction direction = FORWARD; direction < NUM_DIRECTIONS; direction++) {
						if (0 == hits[direction]) {
							pendingInserts[whichHashTable]->recordSeed(seed, whichInputFile);
						}
						seed = ~seed;
					}
				}
				else {
					_int64 nHits[NUM_DIRECTIONS];
					const unsigned *hits[NUM_DIRECTIONS];

					index->lookupSeed32(seed, &nHits[FORWARD], &hits[FORWARD], &nHits[RC], &hits[RC]);

					if (0 == nHits[FORWARD] + nHits[RC]) {
						if (seed.isBiggerThanItsReverseComplement()) {
							seed = ~seed;
						}

						GetHashTableAndKeyFromSeed(seedLength, seed, hashTableBits, &key, &whichHashTable);
						pendingInserts[whichHashTable]->recordSeed(seed, whichInputFile);
					}
				} // if 34/64 index
			} // for each seed in the read
		} // for each read in the input file

		delete readSupplier;
	} // for each input file

	for (unsigned i = 0; i < nHashTables; i++) {
		pendingInserts[i]->flush();
	}
}


void
KMerAlignerContext::updateStats(
AlignerStats* stats,
Read* read,
AlignmentResult result,
int score,
int mapq)
{
}

const unsigned KMerAlignerContext::hashTableSize = 45000000000;

void
KMerAlignerContext::typeSpecificBeginIteration()
{
	unsigned seedSize = index->getSeedLength();

	nHashTables = 1 << (seedSize * 2) % 8;
	keyBytes = (seedSize * 2) / 8;

	if (nHashTables < 64) {
		nHashTables *= 256;
		keyBytes--;
	}

	hashTables = new SNAPHashTable*[nHashTables];
	hashTableLocks = new ExclusiveLock[nHashTables];

	for (unsigned i = 0; i < nHashTables; i++) {
		hashTables[i] = new SNAPHashTable(hashTableSize / nHashTables * (i == 0 ? 20 : 1), keyBytes, 2, options->nInputs, 0xffff);
		InitializeExclusiveLock(&hashTableLocks[i]);
	}

	readSupplierGenerators = new ReadSupplierGenerator *[options->nInputs];
	for (unsigned i = 0; i < options->nInputs; i++) {
		ReaderContext context(readerContext);
		readSupplierGenerators[i] = options->inputs[i].createReadSupplierGenerator(options->numThreads, context);
	}

	outputFile = fopen("f:\\temp\\kmers.txt", "w");
}
void
KMerAlignerContext::typeSpecificNextIteration()
{
}

void KMerAlignerExtension::finishAlignment()
{
	context->finish();
}

void KMerAlignerContext::finish()
{
	unsigned hashTableBits = cheezyLogBase2(nHashTables);
	unsigned seedLen = index->getSeedLength();
	_uint64 singletons = 0;
	_uint64 hits = 0;

	for (_uint64 whichHashTable = 0; whichHashTable < nHashTables; whichHashTable++) {
		for (_uint64 whichEntry = 0; whichEntry < hashTableSize / nHashTables; whichEntry++) {
			void *entry = hashTables[whichHashTable]->getEntry(whichEntry);
			if (!hashTables[whichHashTable]->doesEntryHaveInvalidValue(entry)) {
				_uint64 totalCount = 0;
				for (unsigned whichFile = 0; whichFile < options->nInputs; whichFile++) {
					totalCount += hashTables[whichHashTable]->getValueFromEntry(entry, whichFile);
				}

				if (totalCount > 1) {
					hits++;

					_uint64 seedBinary = hashTables[whichHashTable]->getKey(entry);
					seedBinary |= (whichHashTable << (_uint64)(seedLen * 2 - hashTableBits));
					Seed seed;
					for (int whichBase = 0; whichBase < seedLen; whichBase++) {
						seed.setBase(seedLen - whichBase, seedLen, (seedBinary >> (whichBase * 2)) & 3);
					}

					char seedString[100];
					seed.toString(seedString, seedLen);
					seedString[seedLen] = '\0';
					fprintf(outputFile, "%08lld\t%s", totalCount, seedString);
					for (unsigned whichFile = 0; whichFile < options->nInputs; whichFile++) {
						fprintf(outputFile, "\t%d", hashTables[whichHashTable]->getValueFromEntry(entry, whichFile));
					}
					fprintf(outputFile, "\n");
				}
				else if (1 == totalCount) {
					singletons++;
				}
			}
		}
	}


	fprintf(stderr, "%lld hits, %lld singletons\n", hits, singletons);
	fclose(outputFile);
}