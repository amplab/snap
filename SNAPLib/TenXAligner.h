/*++

Module Name:

TenX.h

Abstract:

Headers for the 10x aligner

Authors:

Hongyi Xin and Bill Bolosky, May, 2016

Environment:

User mode service.

Revision History:

Cloned and modified the paired end aligner

--*/

#pragma once
#include "stdafx.h"
#include "AlignerContext.h"
#include "ReadSupplierQueue.h"

struct TenXAlignerStats;

class TenXAlignerContext : public AlignerContext
{
public:

	TenXAlignerContext(AlignerExtension* i_extension = NULL);

protected:

	// AlignerContext

	virtual bool initialize();

	virtual AlignerStats* newStats();

	virtual void runTask();

	virtual void runIterationThread();

	// for subclasses

	virtual void updateStats(TenXAlignerStats* stats, Read* read0, Read* read1, PairedAlignmentResult* result, bool useful0, bool useful1);

	MappingMode_t getMappingMode() { return m_tenx; }

protected:

	virtual void typeSpecificBeginIteration();
	virtual void typeSpecificNextIteration();

	PairedReadSupplierGenerator *pairedReadSupplierGenerator;

	int					minSpacing;
	int					maxSpacing;
	_uint64				maxBarcodeSize;
	unsigned			minPairsPerCluster;
	_uint64				minClusterSpan;
	bool				forceSpacing;
	bool				noSingle;
	unsigned			intersectingAlignerMaxHits;
	unsigned			maxCandidatePoolSize;
	const char			*fastqFile1;
	bool				ignoreMismatchedIDs;
	bool				quicklyDropUnpairedReads;
	double				unclusteredPenalty;
	unsigned			clusterEDCompensation;
	unsigned			maxClusterNum;

	friend class AlignerContext2;
};

struct TenXAlignerOptions : public AlignerOptions
{
	TenXAlignerOptions(const char* i_commandLine);

	virtual void usageMessage();

	virtual bool parse(const char** argv, int argc, int& n, bool *done);

	virtual bool isPaired() { return true; }

	int					minSpacing;
	int					maxSpacing;
	_uint64				maxBarcodeSize; // make it 64 bit so that the carry won't overflow...
	unsigned			minPairsPerCluster;
	_uint64				minClusterSpan;
	bool				forceSpacing;
	bool				noSingle;
	unsigned			intersectingAlignerMaxHits;
	unsigned			maxCandidatePoolSize;
	bool				quicklyDropUnpairedReads;
	double				unclusteredPenalty;
	unsigned			clusterEDCompensation;
	unsigned			maxClusterNum;
};
