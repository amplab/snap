/*++

Module Name:

KMerAligner.cpp

Abstract:

Function for running the kmer counter functionality.

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
#include "CountingHash.h"

using namespace std;
using util::stringEndsWith;


static const unsigned DEFAULT_HASHTABLE_SIZE = 1; // 40;  //The hashtable size in Gb.
static const unsigned DEFAULT_BLOOM_SIZE = 1; // 40; // The size of the bloomfilter in Gb.  


KMerAlignerContext::KMerAlignerContext(AlignerExtension* i_extension)
: extension(this), AlignerContext(0, NULL, NULL, &extension)
{
	KMerAlignerExtension *extension = new KMerAlignerExtension(this);
}

void KMerAlignerContext::initialize()
{
	AlignerContext::initialize();
	KmerOptions* options2 = (KmerOptions*)options;
	usingBloom = options2->usingBloom; 
	outputDir = options2->outputDir;
	prettyPrint = options2->prettyPrint;
	const unsigned filenameBufferSize = MAX_PATH + 1;
	char filenameBuffer[filenameBufferSize];
	if (prettyPrint){
		snprintf(filenameBuffer, filenameBufferSize, "%s%ckmer_counts.txt", outputDir, PATH_SEP);
	}
	else {
		snprintf(filenameBuffer, filenameBufferSize, "%s%ckmer_counts.bkmer", outputDir, PATH_SEP);
	}
	outputFile = fopen(filenameBuffer, "w");
	seedLength = options2->seedLength;
	computeBias = options2->computeBias;
	addCounts = options2->addCounts; 
	joinCounts = options2->joinCounts; 
	hashTableCapacity = (_uint64)(1000000000 * unsigned(options2->hashTableGb));
	bloomFilterSize = (_uint64) (1000000000 * unsigned(options2->bloomFilterGb));
	if (!options2->usingIndex){
		index = NULL;
	}
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

unsigned KMerAlignerContext::nInputs; 
unsigned KMerAlignerContext::keyBytes; 

struct SeedAndInputFile {
	Seed		seed;
	unsigned	whichInputFile;
};

void GetHashTableAndKeyFromSeed(unsigned seedLength, Seed seed, unsigned hashTableBits, _uint64 *o_key, unsigned *o_whichHashTable)
{
	*o_whichHashTable = (unsigned)(seed.getBases() >> (seedLength * 2 - hashTableBits));
	*o_key = seed.getBases() & ((((_uint64)1) << (seedLength * 2 - hashTableBits))-1);
}

void KMerAlignerContext::printKmerHashTable(SNAPHashTable * hashTable, unsigned hashTableBits, unsigned whichHashTable, FILE * fileToPrintTo)
{
	_uint64 singletons = 0;
	_uint64 hits = 0; 

	for (_uint64 whichEntry = 0; whichEntry < hashTable->tableSize; whichEntry++) {
//	for (_uint64 whichEntry = 0; whichEntry < 100; whichEntry++) {
		void *entry = hashTable->getEntry(whichEntry);
		if (!hashTable->doesEntryHaveInvalidValue(entry)) {
			_uint64 totalCount = 0;
			for (unsigned whichFile = 0; whichFile < nInputs; whichFile++) {
				totalCount += hashTable->getValueFromEntry(entry, whichFile);
			}

			if (totalCount > 0) {
				hits++;

				_uint64 seedBinary = hashTable->getKey(entry);
				_uint64 hashBinary = ((_uint64) whichHashTable << (_uint64)((seedLength * 2) - hashTableBits)); //
				seedBinary |= ((_uint64) whichHashTable << (_uint64)((seedLength * 2) - hashTableBits));
				Seed seed;
				for (int whichBase = 0; whichBase < seedLength; whichBase++) {
					seed.setBase(seedLength - (whichBase +1), seedLength, (seedBinary >> (whichBase * 2)) & 3);
				}

				char seedString[100];
				seed.toString(seedString, seedLength);
				seedString[seedLength] = '\0';
				fprintf(fileToPrintTo, "%08lld\t%s", totalCount, seedString);
				for (unsigned whichFile = 0; whichFile < nInputs; whichFile++) {
					fprintf(fileToPrintTo, "\t%d", hashTable->getValueFromEntry(entry, whichFile));
				}
				fprintf(fileToPrintTo, "\n");
				if (1 == totalCount){
					singletons++;
				}
			}
			else if (1 == totalCount) {
				singletons++;
			}
		}
	}
	fprintf(stderr, "%lld hits, %lld singletons\n", hits, singletons);
}

class PendingInserts {
public:
	PendingInserts(KMerAlignerContext * i_kmerAligner, unsigned i_whichHashTable, unsigned i_seedLength): 
		kmerAligner(i_kmerAligner), usedCount(0), whichHashTable(i_whichHashTable), hashDumpCounter(0)
/*	PendingInserts(SNAPHashTable *i_hashTable, ExclusiveLock *i_hashTableLock, unsigned i_nInputFiles, unsigned i_nHashTables, unsigned i_seedLength, bool usingBloom, const char * i_outputDir, unsigned *i_hashDumpCounter) :
		hashTable(i_hashTable),
		hashTableLock(i_hashTableLock),
		usedCount(0),
		nInputFiles(i_nInputFiles),  
		seedLength(i_seedLength),  
		outputDir(i_outputDir),
		hashDumpCounter(i_hashDumpCounter)*/
	{
		unsigned nHashTables = i_kmerAligner->nHashTables;
		hashTable = i_kmerAligner->hashTables[whichHashTable]; 
		hashTableLock = &(i_kmerAligner->hashTableLocks[whichHashTable]); 
		nInputFiles = i_kmerAligner->nInputs;
		seedLength = i_kmerAligner->seedLength; 
		usingBloom = i_kmerAligner->usingBloom; 
		outputDir = i_kmerAligner->outputDir; 

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

	bool recordSeed(Seed seed, unsigned whichInputFile)
	{
		_ASSERT(usedCount < batchSize);
		inserts[usedCount].seed = seed;
		inserts[usedCount].whichInputFile = whichInputFile;
		usedCount++;
		if (usedCount == batchSize) {
			flush();
		}
		return true;
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
				insertArray[inserts[i].whichInputFile] = usingBloom ? 2:1 ;  //If we are using a bloomfilter, the initial count will be 2.
				// if it passes the bloom filter, add it to the hash. 
				if (!hashTable->Insert(key, insertArray)) {
					WriteErrorMessage("Hash table %d is full\n", whichHashTable);
					//dump the hash table if it is full
					hashTable->sortHash(SNAPHashTable::keyCompare);
					const unsigned filenameBufferSize = MAX_PATH + 1;
					char filenameBuffer[filenameBufferSize];
					size_t bytesWrittenThisHashTable;
					FILE * dumpToFile; 
					if (kmerAligner->prettyPrint){
						snprintf(filenameBuffer, filenameBufferSize, "%s%ckmerHashTable%d_%d.txt", outputDir, PATH_SEP, whichHashTable, *hashDumpCounter);
						dumpToFile = fopen(filenameBuffer, "w");
						kmerAligner->printKmerHashTable(hashTable, hashTableBits, whichHashTable, dumpToFile);
					}
					else{
						//make the tableSize the usedElementCount so we don't write all of the empty slots. 
						hashTable->tableSize = hashTable->usedElementCount; 
						snprintf(filenameBuffer, filenameBufferSize, "%s%ckmerHashTable%d_%d", outputDir, PATH_SEP, whichHashTable, *hashDumpCounter);
						dumpToFile = fopen(filenameBuffer, "w");
						hashTable->saveToFile(dumpToFile, &bytesWrittenThisHashTable);
					}
					fclose(dumpToFile); 
					*hashDumpCounter++; 
					//save parameters for creating a new hashtable. 
					unsigned int i_tableSize = hashTable->tableSize; 
					unsigned int i_keySizeInBytes = hashTable->GetKeySizeInBytes(); 
					unsigned int i_valueSizeInBytes = hashTable->valueSizeInBytes;
					unsigned int i_valueCount = hashTable->valueCount;
					unsigned int i_invalidValue = hashTable->invalidValueValue; 
					delete hashTable;
					//create a new hashtable. 
					hashTable = new SNAPHashTable(i_tableSize, i_keySizeInBytes, i_valueSizeInBytes, i_valueCount, i_invalidValue);
					if (!hashTable->Insert(key, insertArray)){
						WriteErrorMessage("New hashtable shouldn't be full");
						soft_exit(1);
					}
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
	KMerAlignerContext		*		kmerAligner; 
	static const unsigned batchSize;
	unsigned						usedCount;
	unsigned						whichHashTable;
	unsigned						nInputFiles;
	unsigned						nHashTables; 
	SeedAndInputFile	*			inserts;
	SNAPHashTable		*			hashTable;
	ExclusiveLock		*			hashTableLock;
	SNAPHashTable::ValueType	*	insertArray;
	unsigned						hashTableBits;
	unsigned						seedLength; 
	unsigned		*				hashDumpCounter;
	bool			 				usingBloom; 
	const char *					outputDir; 
	

};

const unsigned PendingInserts::batchSize = 1000;
bool computeBias = false; 

KmerOptions::KmerOptions(const char* i_commandLine)
: AlignerOptions(i_commandLine, false), 
	usingBloom(true),
	usingIndex(true),
	outputDir(".\\"),
	seedLength(DEFAULT_SEED_SIZE),
	computeBias(false),
	hashTableGb(DEFAULT_HASHTABLE_SIZE),
	bloomFilterGb(DEFAULT_BLOOM_SIZE),
	prettyPrint(false),
	joinCounts(false), 
	addCounts(false) 
{
}

void KmerOptions::usageMessage()
{
	WriteErrorMessage(
		"Usage: \nsnap kmer <indexDir><inputFile(s)>[<options>] where <input file(s) is a list of files to process.\n"
		"Options:\n"
		"  -D   Directory name to dump output files (full tables and final kmer counts).  Default is the current working directory.\n"
		"  -noG   Don't filter kmers against the given index in indexDir.  In this case, use \"-\" as indexDir.\n"
		"  -s   Seed size (If a genome index is used, the seed size will be the size used to build the given index.\n"
		"       Otherwise, the defaults is %d\n"
		"  -h   Hash table slack (default: %.1f)\n"
		"  -hg19   Use pre-computed table bias for hg19, which results in better speed, balance, and memory footprint but may not work for other references.\n"
		"  -cb  Compute bias for tables before counting.  If neither cp or hg19 options are used, the hash tables will have equal sizes.\n"
		"  -hs  Hash tables size in Gb. (default: %d)\n"
		"  -bfs Bloom filter size in Gb. (default: %d)\n"
		"  -t   Number of threads (default is one per core)\n"
		"  -b   Bind each thread to its processor (off by default)\n"
		"  -sn  Sort output file by the total count of the kmer. Files are normally sorted by kmer sequence. (not implemented yet)\n"
		"  -sm  Memory to use for sorting in Gb\n"
		"  -nobf   Don't filter out singletons using a bloom filter.\n",  
		"  -merge  Join the given kmer count files by extending the array counts. (input files should be in .snapkmer format).\n"
		"  -add    Join the given .snapkmer files by adding the counts.  (Input files should be in .snapkmer format).\n"
		"  -p   Print the output in text (tsv) format.\n",
		DEFAULT_SEED_SIZE, 
		DEFAULT_SLACK, 
		DEFAULT_HASHTABLE_SIZE, 
		DEFAULT_BLOOM_SIZE);
}

bool KmerOptions::parse(const char** argv, int argc, int& n, bool *done)
{
	*done = false;

	if (strcmp(argv[n], "-D") == 0) {
		if (n + 1 < argc) {
			outputDir = argv[n + 1];
			n += 1;
			return true;
		}
		return false;
	}
	else if (strcmp(argv[n], "-noG") == 0) {
		if (n + 1 < argc) {
			usingIndex = false;
			return true;
		}
		return false;
	}
	else if (strcmp(argv[n], "-s") == 0) {
		if (n + 1 < argc) {
			seedLength = atoi(argv[n + 1]);
			if (usingIndex){
				WriteErrorMessage(
					"Although you've specified the kmer length, the seedlength used for the given genome will be used.\n"
					"If they are the same, then rest easy.  If they are different, you can either count kmers without\n"
					"filtering against the given genome index, (use \"-\" in place of the index directory name), or \n"
					"remake the index with the desired seedlength.\n");
			}
			n += 1; 
			return true;
		}
		return false; 
	}
	else if (strcmp(argv[n], "-hg19") == 0) {
		computeBias = false;
		return true; 
	}
	else if (strcmp(argv[n], "-cb") == 0) {
		computeBias = true;
		return true; 
	}
	else if (strcmp(argv[n], "-hs") == 0) {
		if (n + 1 < argc){
			hashTableGb = atoi(argv[n + 1]);
			n += 1; 
			return true;
		}
		return false; 
	}
	else if (strcmp(argv[n], "-bfs") == 0){
		if (n + 1 < argc){
			bloomFilterGb = atoi(argv[n + 1]);
			n += 1;
			return true;
		}
		return false;
	}
	else if (strcmp(argv[n], "-sn") == 0){
		WriteErrorMessage("Sorting by kmer count is not yet implemented.\n");
		return true; 
	}
	else if (strcmp(argv[n], "-merge") == 0){
		if (addCounts){
			WriteErrorMessage("-merge and -add are mutually exclusive.  Please pick one.\n");
			soft_exit(1); 
		}
		joinCounts = true; 
		return true; 
	}
	else if (strcmp(argv[n], "-add") == 0){
		if (joinCounts){
			WriteErrorMessage("-merge and -add are mutually exclusive.  Please pick one.\n");
			soft_exit(1);
		}
		addCounts = true; 
		return true; 
	}
	else if (strcmp(argv[n], "-p") == 0){
		prettyPrint = true;
		return true; 
	}
	return AlignerOptions::parse(argv, argc, n, done);
}


void
KMerAlignerContext::runIterationThread()
{
	PreventMachineHibernationWhileThisThreadIsAlive();

	//unsigned seedLength = index->getSeedLength();
	unsigned hashTableBits = cheezyLogBase2(nHashTables);

	PendingInserts **pendingInserts = new PendingInserts *[nHashTables];
	for (unsigned i = 0; i < nHashTables; i++) {
		//pendingInserts[i] = new PendingInserts(hashTables[i], &hashTableLocks[i], options->nInputs, nHashTables, seedLength, bloomFilters[i], outputDir, &hashDumpCounter);
		pendingInserts[i] = new PendingInserts(this, i, seedLength);
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
			unsigned datalength = read->getDataLength(); 
			for (unsigned nextSeedToTest = 0; nextSeedToTest < read->getDataLength() - seedLength; nextSeedToTest++) {
#if 0
				const unsigned prefetchDepth = 0;
				if (doAlignerPrefetch && nextSeedToTest + prefetchDepth < read->getDataLength() - seedLength && Seed::DoesTextRepresentASeed(read->getData() + nextSeedToTest + prefetchDepth, seedLength)) {
					Seed seed(read->getData() + nextSeedToTest + prefetchDepth, seedLength);
					index->prefetchLookup(seed);
				}

#endif // 0
				const char* seedtext = read->getData() + nextSeedToTest; 
				if (!Seed::DoesTextRepresentASeed(read->getData() + nextSeedToTest, seedLength)) {
					continue;
				}

				Seed seed(read->getData() + nextSeedToTest, seedLength);
				if (seed.isBiggerThanItsReverseComplement()) {
					seed = ~seed;
				}

				_uint64 key;
				unsigned whichHashTable;
				GetHashTableAndKeyFromSeed(seedLength, seed, hashTableBits, &key, &whichHashTable);
				
				/* before filtering against our SNAP index of the reference genome, 
				use a Bloom filter to make sure that the seed occurs more than once.
				Singleton kmers are probably just noise in the sequencing run. */
				if (usingBloom){
					_uint64 bloomcount = bloomFilters[whichHashTable]->increment(key);
					if (bloomcount == 0) {  //If this is the first time we've seen this key, then skip it.
						continue;
					}
				}
				/*Filter against some reference genome */
				if (NULL != index){
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
							//						if (seed.isBiggerThanItsReverseComplement()) {
							//						seed = ~seed;
							//				}
							pendingInserts[whichHashTable]->recordSeed(seed, whichInputFile);
						}
					} // if 34/64 index
				}  // if filtering against the index. 
				else{
					pendingInserts[whichHashTable]->recordSeed(seed, whichInputFile); 
				}
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

int KMerAlignerContext::kmerKeyCompare(const void *first, const void *second)
{
	const unsigned* key1 = (const unsigned*)first; 
	const unsigned * key2 = (const unsigned*)second; 
	return *key1 - *key2;  
}

int KMerAlignerContext::kmerCountCompare(const void *first, const void *second)
{
/*	const unsigned *ifirst = (const unsigned*)first;
	const unsigned *isecond = (const unsigned*)second; 
	const unsigned firstvals = *ifirst << keyBytes; 
	int totalCount1 = 0; 
	for (int i = 0; i < nInputs; i++){
		
	}

	for (unsigned whichFile = 0; whichFile < nInputs; whichFile++) {
		totalCount += hashTable->getValueFromEntry(entry, whichFile);
	}
*/
	return 1; 
}

void
KMerAlignerContext::typeSpecificBeginIteration()
{
	nInputs = options->nInputs;

	if (NULL != index){
		seedLength = index->getSeedLength();
	} //else seedLength is already defined (default =20, or given in command line)
	unsigned seedSize = seedLength; // index->getSeedLength();
	nHashTables = 1 << (seedSize * 2) % 8;
	keyBytes = (seedSize * 2) / 8;

	if (nHashTables < 64) {
		nHashTables *= 256;
		keyBytes--;
	}
	hashDumpCounter = 0; 
	//Get bias tables for kmer frequencies (just use pre-computed hg19 bias for now). 
	double *biasTable = NULL;
	int hashTableKeySize = keyBytes;
	if (!computeBias) {
		biasTable = GenomeIndex::hg19_biasTables[hashTableKeySize][seedLength];
	}
	if (NULL == biasTable) {
		WriteErrorMessage("-hg19 not available for this seed length/key size/small-or-large combo.  Computing bias tables the hard way.\n");
		computeBias = true;
	}
	if (computeBias) {
		biasTable = new double[nHashTables];
		WriteErrorMessage("computing bias needs to be implemented.\n");
		soft_exit(1);
	}

	//make a bloom filter
	if (usingBloom){
		bloomFilters = new CountingHash*[nHashTables];
		WriteStatusMessage("Allocating bloomFilters...");
		for (unsigned i = 0; i < nHashTables; i++) {
			double bias = biasTable[i];
			_uint64 capacity = (unsigned)((bloomFilterSize / nHashTables) * bias);
			if (capacity < 100){
				capacity = 100; 
			}
			_uint64 keyCount = (0.5 - capacity * log(2) * log(2)) / log(0.0001);
			bloomFilters[i] = new CountingHash(keyCount, 1);  //max count as 1 makes this a bloom filter 
		}
		WriteStatusMessage("bloomFilters allocated!\n");
	}
	else {
		for (unsigned i = 0; i < nHashTables; i++) {
			bloomFilters[i] = NULL;
		}
	}

	hashTables = new SNAPHashTable*[nHashTables];
	hashTableLocks = new ExclusiveLock[nHashTables];
	WriteStatusMessage("Allocating HashTables...");
	hashTableSize = (unsigned)(hashTableCapacity / (keyBytes + 2 * nInputs)); 
	for (unsigned i = 0; i < nHashTables; i++) {
		double bias = biasTable[i];
		unsigned biasedSize = (unsigned)((hashTableSize / nHashTables) * bias);
		if (biasedSize < 100){
			biasedSize = 100; 
		}
		//biasedSize = 1000000; 
		hashTables[i] = new SNAPHashTable(biasedSize, keyBytes, 2, options->nInputs, 0xffff);
		InitializeExclusiveLock(&hashTableLocks[i]);
	}
	WriteStatusMessage("HashTables allocated!\n");
	readSupplierGenerators = new ReadSupplierGenerator *[options->nInputs];
	for (unsigned i = 0; i < options->nInputs; i++) {
		ReaderContext context(readerContext);
		readSupplierGenerators[i] = options->inputs[i].createReadSupplierGenerator(options->numThreads, context);
	}

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
	_uint64 singletons = 0;
	_uint64 hits = 0;
	size_t totalBytesWritten = 0; 
	
	for (_uint64 whichHashTable = 0; whichHashTable < nHashTables; whichHashTable++) {
		hashTables[whichHashTable]->sortHash(SNAPHashTable::keyCompare); 
		if (prettyPrint){
			printKmerHashTable(hashTables[whichHashTable], hashTableBits, whichHashTable, outputFile);
		}
		else{ 
			size_t bytesWrittenThisHashTable;
			hashTables[whichHashTable]->tableSize = hashTables[whichHashTable]->usedElementCount; 
			if (!hashTables[whichHashTable]->saveToFile(outputFile, &bytesWrittenThisHashTable)) {
				WriteErrorMessage("GenomeIndex::saveToDirectory: Failed to save hash table %d\n", whichHashTable);
				exit(1);
			}
			totalBytesWritten += bytesWrittenThisHashTable;
		}
	}
	fclose(outputFile);
	// print bloom filter stats
	if (usingBloom){
		for (_uint64 whichHashTable = 0; whichHashTable < nHashTables; whichHashTable++) {
			double pError = bloomFilters[whichHashTable]->calculatePerror();
			fprintf(stderr, "For bloomfilter %i, bitsSet: %lld, entries: %lld, error rate: %e\n", whichHashTable, bloomFilters[whichHashTable]->getBitsSet(), bloomFilters[whichHashTable]->getKeyCount(), pError);
		}
	}
}


