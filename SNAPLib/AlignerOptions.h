/*++

Module Name:

    AlignerOptions.h

Abstract:

    Common parameters for running single & paired alignment.

Authors:

    Ravi Pandya, May, 2012

Environment:
`
    User mode service.

Revision History:

    Integrated from SingleAligner.cpp & PairedAligner.cpp

--*/

#pragma once

#include "stdafx.h"
#include "options.h"
#include "Range.h"
#include "Genome.h"
#include "RangeSplitter.h"

#define MAPQ_LIMIT_FOR_SINGLE_HIT 10

struct AbstractOptions
{
    virtual void usageMessage() = 0;

    virtual bool parse(const char** argv, int argc, int& n, bool *done) = 0;
};

enum FileType {UnknownFileType, SAMFile, FASTQFile, BAMFile, GZipFASTQFile, CRAMFile};  // As more as needed

struct SNAPInput {
    SNAPInput() : fileName(NULL), secondFileName(NULL), fileType(UnknownFileType) {}
    const char          *fileName;
    const char          *secondFileName;
    FileType             fileType;

    void readHeader(ReaderContext& context);
    PairedReadSupplierGenerator *createPairedReadSupplierGenerator(int numThreads, const ReaderContext& context);
    ReadSupplierGenerator *createReadSupplierGenerator(int numThreads, const ReaderContext& context);
};

struct AlignerOptions : public AbstractOptions
{
    AlignerOptions(const char* i_commandLine, bool forPairedEnd = false);

    const char         *commandLine;
    const char         *indexDir;
    const char         *transcriptomeDir;
    const char         *contaminationDir;
    const char         *annotation;
    const char         *similarityMapFile;
    int                 numThreads;
    Range               maxDist;
    unsigned            numSeedsFromCommandLine;
    double              seedCoverage;       // Exclusive with numSeeds; this is readSize/seedSize
    bool                seedCountSpecified; // Has either -n or -sc been specified?  This bool is used to make sure they're not both specified on the command line
    Range               maxHits;
    bool                computeError;
    bool                bindToProcessors;
    bool                ignoreMismatchedIDs;
    const char         *outputFileTemplate;
    int                 nInputs;
    SNAPInput          *inputs;
    ReadClippingType    clipping;
    bool                sortOutput;
    bool                noIndex;
    bool                noDuplicateMarking;
    bool                noQualityCalibration;
    unsigned            sortMemory; // total output sorting buffer size in Gb
    unsigned            filterFlags;
    bool                explorePopularSeeds;
    bool                stopOnFirstHit;
	bool				useM;	// Should we generate CIGAR strings using = and X, or using the old-style M?
    unsigned            gapPenalty; // if non-zero use gap penalty aligner
    unsigned            misalignThreshold; // For error reporting: min distance from real location to mark a read as misaligned
    AbstractOptions    *extra; // extra options
    const char         *rgLineContents;
    const char         *perfFileName;
    bool                useTimingBarrier;
    unsigned            extraSearchDepth;
    const char         *defaultReadGroup; // if not specified in input
    bool                ignoreSecondaryAlignments;
    bool                preserveClipping;
    float               expansionFactor;

    //Added this back in for now
    unsigned            confDiff;

    //Quality filtering options
    float               minPercentAbovePhred;
    unsigned            minPhred;
    unsigned            phredOffset;

    void usage();

    virtual void usageMessage();

    virtual bool parse(const char** argv, int argc, int& n, bool *done);

    enum FilterFlags
    {
        FilterUnaligned =           0x0001,
        FilterSingleHit =           0x0002,
        FilterMultipleHits =        0x0004,
    };

    bool passFilter(Read* read, AlignmentResult result);
    
    virtual bool isPaired() { return false; }
};
