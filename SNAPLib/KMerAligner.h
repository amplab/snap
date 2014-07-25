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

protected:

	// AlignerContext overrides

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

	bool isPaired() { return false; }


private:

	unsigned				nHashTables;
	unsigned				keyBytes;
	SNAPHashTable **		hashTables;
	ExclusiveLock *			hashTableLocks;
	KMerAlignerExtension	extension;
	FILE					*outputFile;

	static const unsigned hashTableSize;
};
