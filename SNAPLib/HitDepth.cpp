/*++

Module Name:

    HitDepth.cpp

Abstract:

    Special function to find the seed with the minimum hit depth to align each
    locus in a set of contigs.

Authors:

    Bill Bolosky, March, 2022

Environment:

    User mode service.

Revision History:

--*/

#include "stdafx.h"
#include "Compat.h"
#include "Histogram.h"
#include "exit.h"
#include "AlignerOptions.h"
#include "GenomeIndex.h"

#if HIT_DEPTH_COUNTING
void CountHitDepthUsage()
{
    fprintf(stderr, "This is a special function to compute data for a paper.\n");
    fprintf(stderr, "It looks at every locus in a set of contigs and finds the seed\n");
    fprintf(stderr, "with the fewest hits that contains the correct alignment across\n");
    fprintf(stderr, "a range of seed sizes.  It is intended to show some concept of\n");
    fprintf(stderr, "'difficulty' of aligning different portions of the genome.\n\n");
    fprintf(stderr, "usage:\n");
    fprintf(stderr, "  depth index-filename-base minSeedSize maxSeedSize seedSizeForBaseAlignment outputFile {contigFile}\n");
    fprintf(stderr, "        index-filename-base is the base pathname for a set of SNAP indices of various seed sizes.\n");
    fprintf(stderr, "                            the individual index names are the base concatenated with the seed size.\n");
    fprintf(stderr, "        minSeedSize and maxSeedSize are the min and max of the range of seed sizes to test\n");
    fprintf(stderr, "        seedSizeForBaseAlignment is the seed size used to align the reads to determine the 'correct' alignment\n");
    fprintf(stderr, "        ouputFile is the output filename\n");
    fprintf(stderr, "        contigFile is a file with a list of contigs to compute (one per line).  If it is not specified then\n");
    fprintf(stderr, "                   the standard list of primary contigs from hg38 are used (i.e., chr1-chr22, chrX, chrY, and chrM)\n");

    soft_exit(1);
} // CountHitDepthUsage

char* contigsToUse[] = {
    "chr1",
    "chr2",
    "chr3",
    "chr4",
    "chr5",
    "chr6",
    "chr7",
    "chr8",
    "chr9",
    "chr10",
    "chr11",
    "chr12",
    "chr13",
    "chr14",
    "chr15",
    "chr16",
    "chr17",
    "chr18",
    "chr19",
    "chr20",
    "chr21",
    "chr22",
    "chrX",
    "chrY",
    "chrM"
};

GenomeIndex *g_HitDepthGenomeIndex = NULL;

void CountHitDepthLoadIndex(const char* indexFileNameBase, int seedSize) 
{
    printf("Loading index for seed size %d...", seedSize);

    _int64 startLoadTime = timeInMillis();

    size_t indexFileNameBufferSize = strlen(indexFileNameBase) + 20;// 20 is more than the largest int + 1 for \0, so plenty of space
    char* indexFileName = new char[indexFileNameBufferSize];

    snprintf(indexFileName, indexFileNameBufferSize, "%s%d", indexFileNameBase, seedSize);

    g_HitDepthGenomeIndex = GenomeIndex::loadFromDirectory(indexFileName, true, true);

    if (g_HitDepthGenomeIndex == NULL) {
        fprintf(stderr, "\nIndex load failed\n");
        soft_exit(1);
    }

    printf("%llds\n", (timeInMillis() - startLoadTime + 500) / 1000);
}

void CountHitDepth(int argc, const char** argv)
{
    if (argc != 5 && argc != 6) {
        CountHitDepthUsage();
    }

    const char* indexFileNameBase = argv[0];

    int minSeedSize = atoi(argv[1]);
    int maxSeedSize = atoi(argv[2]);

    if (minSeedSize <= 0 || maxSeedSize < minSeedSize) {
        fprintf(stderr, "Min seed size must be strictly positive and no greater than max seed size\n");
        soft_exit(1);
    }

    int seedSizeForBaseAlignment = atoi(argv[3]);
    if (seedSizeForBaseAlignment < minSeedSize || seedSizeForBaseAlignment > maxSeedSize) {
        fprintf(stderr, "Seed size for base alignment must be in the range of min seed size - max seed size\n");
        soft_exit(1);
    }

    if (argc == 6) {
        fprintf(stderr, "Contig file not yet implemented\n");
        soft_exit(1);
    }

    FILE* outputFile = fopen(argv[4], "w");
    if (NULL == outputFile) {
        fprintf(stderr, "Unable to open output file\n");
        soft_exit(1);
    }

    CountHitDepthLoadIndex(indexFileNameBase, seedSizeForBaseAlignment);

} // CountHitDepth
#endif // HIT_DEPTH_COUNTING
