/*++

Module Name:

	KMerAligner.h

Abstract:

Functions for running the kmer counter sub-program.

Authors:

Tracy Ballinger & Bill Bolosky, July, 2014

Environment:

User mode service.

Revision History:


--*/

#pragma once
#include "stdafx.h"
#include "AlignerContext.h"
#include "AlignerStats.h"
#include "ReadSupplierQueue.h"
#include "AlignmentResult.h"
#include "SingleAligner.h"
#include "CountingHash.h"

class KMerAlignerContext;

class KMerAlignerExtension : public AlignerExtension
{
public:
	KMerAlignerExtension(KMerAlignerContext *i_context) : context(i_context) {}
	virtual void finishAlignment();

private:
	KMerAlignerContext *context;

};


class KMerAlignerContext : public AlignerContext
{
public:

	KMerAlignerContext(AlignerExtension* i_extension = NULL);

	void finish(); 
	void sortKmerHash(SNAPHashTable hashtable, bool orderByCount);
	void printKmerHashTable(SNAPHashTable * hashTable, unsigned hashTableBits, unsigned whichHashTable, FILE * fileToPrintTo); 

protected:

	// AlignerContext overrides

	virtual void initialize();

	virtual AlignerStats* newStats();

	virtual void runTask();

	virtual void runIterationThread();

	virtual void typeSpecificBeginIteration();
	virtual void typeSpecificNextIteration();

	// for subclasses

	virtual void writeRead(Read* read, const SingleAlignmentResult &result, bool secondaryAlignment) {
		WriteErrorMessage("Shouldn't be here\n");
		soft_exit(1);
	}

	virtual void updateStats(AlignerStats* stats, Read* read, AlignmentResult result, int score, int mapq);

	//RangeSplittingReadSupplierGenerator   *readSupplierGenerator;

	ReadSupplierGenerator **readSupplierGenerators;

	friend class AlignerContext2;

	friend class PendingInserts; 

	bool isPaired() { return false; }

	static int kmerKeyCompare(const void *, const void *); 
	static int kmerCountCompare(const void *, const void*);


private:

	unsigned				nHashTables;
	static unsigned			nInputs; 
	static unsigned			keyBytes;
	SNAPHashTable **		hashTables;
	ExclusiveLock *			hashTableLocks;
	KMerAlignerExtension	extension;
	const char *			outputDir; 
	FILE					*outputFile;
	bool					usingBloom; //flag for whether we're using a bloom filter or not. 
	CountingHash**			bloomFilters;  //one bloomFilter per hashtable
	_uint64					bloomFilterSize;
	bool					usingIndex; //flag for whether we're filtering against an index or not. 
	_uint64					hashTableCapacity; //size of hash table in bits
	unsigned				hashTableSize; //number of entries hash table can have
	unsigned				hashDumpCounter; 
	bool					computeBias; 
	bool					prettyPrint; 
	int						seedLength;
	bool					joinCounts; 
	bool					addCounts; 
};

struct KmerOptions : public AlignerOptions
{
	KmerOptions(const char* i_commandLine);

	virtual void usageMessage();

	virtual bool parse(const char** argv, int argc, int& n, bool *done);

	virtual bool isPaired() { return false; }

	bool			usingBloom;
	bool			usingIndex; 
	const char *	outputDir;
	int				seedLength; 
	bool			computeBias; 
	int				hashTableGb; 
	int				bloomFilterGb; 
	bool			prettyPrint;
	bool			joinCounts;
	bool			addCounts;
};
