/*++

Module Name:

    AlignerOptions.h

Abstract:

    Common parameters for running single & paired alignment.

Authors:

    Ravi Pandya, May, 2012

Environment:

    User mode service.

Revision History:

    Integrated from SingleAligner.cpp & PairedAligner.cpp

--*/

#pragma once

#include "stdafx.h"
#include "options.h"
#include "Genome.h"
#include "Read.h"

#define INSTRUMENTATION_FOR_PAPER 0 // Turn this on to generate raw data about hit sets and their intersections for the paper

#if INSTRUMENTATION_FOR_PAPER

#define MAX_HIT_SIZE_LOG_2  15  // This is for the instrumentation
extern _int64 g_alignmentTimeByHitCountsOfEachSeed[MAX_HIT_SIZE_LOG_2+1][MAX_HIT_SIZE_LOG_2+1];  // In the paired-end aligner, if you have seeds A and B with hit set sizes |A| and |B| then the total time in ns gets added into g_alignmentTimeByHitCountsOfEachSeed[log2(|A|)][log2(|B|)]
extern _int64 g_alignmentCountByHitCountsOfEachSeed[MAX_HIT_SIZE_LOG_2 + 1][MAX_HIT_SIZE_LOG_2 + 1];  // Same as above, but just add one per time.
extern _int64 g_scoreCountByHitCountsOfEachSeed[MAX_HIT_SIZE_LOG_2 + 1][MAX_HIT_SIZE_LOG_2 + 1];
extern _int64 g_setIntersectionSizeByHitCountsOfEachSeed[MAX_HIT_SIZE_LOG_2 + 1][MAX_HIT_SIZE_LOG_2 + 1];
extern _int64 g_100xtotalRatioOfSetIntersectionSizeToSmallerSeedHitCountByCountsOfEachSeed[MAX_HIT_SIZE_LOG_2 + 1][MAX_HIT_SIZE_LOG_2 + 1];
extern _int64 g_totalSizeOfSmallerHitSet;
extern _int64 g_totalSizeOfSetIntersection;
extern _int64 g_alignmentsWithMoreThanOneCandidateWhereTheBestCandidateIsScoredFirst;
extern _int64 g_alignmentsWithMoreThanOneCandidate;

#endif // INSTRUMENTATION_FOR_PAPER

#define MAPQ_LIMIT_FOR_SINGLE_HIT 10

struct AbstractOptions
{
    virtual void usageMessage() = 0;

    virtual bool parse(const char** argv, int argc, int& n, bool *done) = 0;
};

enum FileType {UnknownFileType, SAMFile, FASTQFile, BAMFile, InterleavedFASTQFile, CRAMFile};  // Add more as needed

struct SNAPFile {
	SNAPFile() : fileName(NULL), secondFileName(NULL), fileType(UnknownFileType), isStdio(false), omitSQLines(false) {}
    const char          *fileName;
    const char          *secondFileName;
    FileType             fileType;
    bool                 isCompressed;
    bool                 isStdio;           // Only applies to the first file for two-file inputs
	bool				 omitSQLines;		// Special (formerly) undocumented option for Charles Chiu's group.  Mostly a bad idea.

    PairedReadSupplierGenerator *createPairedReadSupplierGenerator(int numThreads, bool quicklyDropUnpairedReads, const ReaderContext& context);
    ReadSupplierGenerator *createReadSupplierGenerator(int numThreads, const ReaderContext& context);
    static bool generateFromCommandLine(const char **args, int nArgs, int *argsConsumed, SNAPFile *snapFile, bool paired, bool isInput);
};

struct AlignerOptions : public AbstractOptions
{
    AlignerOptions(const char* i_commandLine, bool forPairedEnd = false);

    const char         *commandLine;
    const char         *indexDir;
    const char         *similarityMapFile;
    int                 numThreads;
    unsigned            maxDist;
    float               maxDistFraction;
    unsigned            maxDistForIndels;
    unsigned            numSeedsFromCommandLine;
    double              seedCoverage;       // Exclusive with numSeeds; this is readSize/seedSize
    bool                seedCountSpecified; // Has either -n or -sc been specified?  This bool is used to make sure they're not both specified on the command line
    unsigned            maxHits;
    int                 minWeightToCheck;
    bool                bindToProcessors;
    bool                ignoreMismatchedIDs;
    SNAPFile            outputFile;
    int                 nInputs;
    SNAPFile           *inputs;
    ReadClippingType    clipping;
    bool                sortOutput;
    bool                noIndex;
    bool                noDuplicateMarking;
    bool                noQualityCalibration;   // This doesn't appear to be used.  
    unsigned            sortMemory; // total output sorting buffer size in Gb
    unsigned            filterFlags;
    bool                explorePopularSeeds;
    bool                stopOnFirstHit;
	bool				useM;	// Should we generate CIGAR strings using = and X, or using the old-style M?
    unsigned            gapPenalty; // if non-zero use gap penalty aligner
    AbstractOptions    *extra; // extra options
    const char         *rgLineContents;
    const char         *perfFileName;
    bool                useTimingBarrier;
    unsigned            extraSearchDepth;
    const char         *defaultReadGroup; // if not specified in input
    bool                ignoreSecondaryAlignments; // on input, default true
    int                 maxSecondaryAlignmentAdditionalEditDistance;
	int					maxSecondaryAlignments;
    int                 maxSecondaryAlignmentsPerContig;
    int                 flattenMAPQAtOrBelow;
    bool                preserveClipping;
    float               expansionFactor;
    bool                noUkkonen;
    bool                noOrderedEvaluation;
	bool				noTruncation;
    bool                useAffineGap;
    bool                useSoftClipping;
    unsigned            matchReward;
    unsigned            subPenalty;
    unsigned            gapOpenPenalty;
    unsigned            gapExtendPenalty;
	unsigned			minReadLength;
	bool				mapIndex;
	bool				prefetchIndex;
    size_t              writeBufferSize;
    bool                dropIndexBeforeSort;
    bool                killIfTooSlow;
    const char *        sortIntermediateDirectory;
    bool                profile;
    bool                profileAffineGap;
    bool                ignoreAlignmentAdjustmentsForOm;
    bool                emitInternalScore;
    char                internalScoreTag[3];
	bool				altAwareness;
    int                 maxScoreGapToPreferNonALTAlignment;
    bool                emitALTAlignments;
    
    static bool         useHadoopErrorMessages; // This is static because it's global (and I didn't want to push the options object to every place in the code)
    static bool         outputToStdout;         // Likewise

    void usage();

    virtual void usageMessage();

    virtual bool parse(const char** argv, int argc, int& n, bool *done);

    enum FilterFlags
    {
        FilterUnaligned =           0x0001,
        FilterSingleHit =           0x0002,
        FilterMultipleHits =        0x0004,
        FilterBothMatesMatch =      0x0008,
		FilterTooShort =            0x0010
    };

    bool passFilter(Read* read, AlignmentResult result, bool tooShort, bool secondaryAlignment);
    
    virtual bool isPaired() { return false; }
};
