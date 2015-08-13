// CountReadsCovering.cpp : Defines the entry point for the console application.
//

#include "stdafx.h"
#include "exit.h"
#include "Genome.h"
#include "pool.h"
#include "avl2.h"
#include "BigAlloc.h"
#include "DataReader.h"
#include "Read.h"
#include "bam.h"
#include "sam.h"
#include "Util.h"


void usage()
{
    fprintf(stderr, "usage: CountReadsCovering index knownGenesFile inputFile outputFile\n");
    fprintf(stderr, "Counts the reads that cover each gene and exon in the knownGenesFile.\n");
    soft_exit(1);
}

struct Isoform;

struct Exon {
    Exon() : isoform(NULL), mappedCount(0), start(0), end(0) {}

    Isoform *           isoform;
    volatile __int64    mappedCount;
    GenomeLocation      start;
    GenomeLocation      end;            // First base NOT covered by this Exon (i.e,. start + length)
};

struct Isoform {
    Isoform() : exons(NULL), mappedCount(0), knownGeneLine(NULL), next(0) {}

    unsigned            nExons;
    Exon *              exons;
    volatile _int64     mappedCount;
    char *              knownGeneLine;
    Isoform *           next;   // List of all of them for generating the final output.
};

class Region { // A chunk of genome space with the same set of Isoforms mapping to it.
public:

    Region(GenomeLocation _start, GenomeLocation _end, unsigned _nExons, unsigned _nExonSpace) : start(_start), end(_end), nExons(_nExons), nExonSpace(_nExonSpace), next(NULL)
    {
        exons = new Exon *[nExonSpace];
    }

    Region(const Region *peer) : start(peer->start), end(peer->end), nExons(peer->nExons), nExonSpace(peer->nExonSpace), next(NULL)
    {
        exons = new Exon *[nExonSpace];
        for (unsigned i = 0; i < nExons; i++) {
            exons[i] = peer->exons[i];
        }
    }
    GenomeLocation      start;
    GenomeLocation      end;
    Region *            next;           // All regions in GenomeLocation order
    Exon **             exons;          // First base NOT covered by this Exon (i.e,. start + length)
    unsigned            nExons;
    unsigned            nExonSpace;     // How big is exons; this is >= nExons, and allows us to have allocation hysteresis

    GenomeLocation getKey() {
        return start;
    }

    void addExon(Exon *exon) {
        if (nExonSpace == nExons) {
            //
            // Expand the space.
            //
            Exon **newExons = new Exon*[nExonSpace * 2];
            for (unsigned i = 0; i < nExons; i++) {
                newExons[i] = exons[i];
            }
            delete[] exons;
            exons = newExons;
            nExonSpace *= 2;
        }

        exons[nExons] = exon;
        nExons++;
    }
};

AVLTree<Region, GenomeLocation> *regions;
volatile int nRunningThreads;
ReadSupplierGenerator *readSupplierGenerator;
SingleWaiterObject allThreadsDone;
volatile _int64 TotalReads = 0;
volatile _int64 MappedReads = 0;
volatile _int64 IsoformHits = 0;
volatile _int64 ExonHits = 0;

//
// Definitions of the knownGene fields
//
const int kg_name       = 0;
const int kg_chrom      = 1;
const int kg_strand     = 2;
const int kg_txStart    = 3;
const int kg_txEnd      = 4;
const int kg_cdsStart   = 5;
const int kg_cdsEnd     = 6;
const int kg_exonCount  = 7;
const int kg_exonStarts = 8;    // This is a comma separated list
const int kg_exonEnds   = 9;    // This is a comma separated list
const int kg_proteinID  = 10;
const int kg_alignID    = 11;

int ComparePointers(const void *one, const void *two)
{
    if (one < two) return -1;
    if (one == two) return 0;
    return 1;
}

void ReaderThread(void *param)
{
    ReadSupplier *readSupplier = readSupplierGenerator->generateNewReadSupplier();

    _int64 totalReads = 0;
    _int64 mappedReads = 0;
    _int64 isoformHits = 0;
    _int64 exonHits = 0;

    int BJBNs = 0;

    Read *read;
    while (NULL != (read = readSupplier->getNextRead())) {
        const int maxExons = 1000;
        Exon *exons[maxExons];
        Isoform *isoforms[maxExons];
 
        int nExons = 0;
        totalReads++;
        if (read->getOriginalSAMFlags() & SAM_UNMAPPED || read->getOriginalAlignedLocation() == InvalidGenomeLocation) {
            continue;
        }
        mappedReads++;

        //
        // Run through the read and handle each chunk separated by an N in the cigar string.
        //
        GenomeLocation chunkStart = read->getOriginalAlignedLocation();
        GenomeLocation chunkEnd = chunkStart;
        for (int i = 0; i <= read->getOriginalNBAMCigar(); i++) {
            _uint32 cigarOp = read->getOriginalBAMCigar()[i] & 0xf;
            _uint32 cigarCount = (read->getOriginalBAMCigar()[i] >> 4);
            if (i == read->getOriginalNBAMCigar() || cigarOp == BAM_CIGAR_N) {
                if (cigarOp == BAM_CIGAR_N) {
                    BJBNs++;
                }
                for (Region *region = regions->findFirstLessThanOrEqualTo(chunkStart); NULL != region && region->start < chunkEnd; region = region->next) {
                    for (unsigned i = 0; i < region->nExons; i++) {
                        if (nExons < maxExons) {
                            exons[nExons] = region->exons[i];
                            isoforms[nExons] = region->exons[i]->isoform;
                        }
                        nExons++;
                    }
                }

                if (nExons > maxExons) {
                    WriteErrorMessage("Not enough space for exons, %d > %d\n", nExons, maxExons);
                    continue;
                }

                if (i == read->getOriginalNBAMCigar()) {
                    break;
                }

                chunkStart = chunkEnd + cigarCount;
                chunkEnd = chunkStart;
            } // If this is the end of a chunk

            if (cigarOp != BAM_CIGAR_S && cigarOp != BAM_CIGAR_H && cigarOp != BAM_CIGAR_N) {
                chunkEnd += cigarCount;
            }
        } // for each cigar op.

        qsort(exons, nExons, sizeof(Exon *), ComparePointers);
        qsort(isoforms, nExons, sizeof(Isoform *), ComparePointers);

        for (int i = 0; i < nExons; i++) {
            if (i == 0 || exons[i] != exons[i - 1]) {
                exonHits++;
                InterlockedAdd64AndReturnNewValue(&exons[i]->mappedCount, 1);
            }

            if (i == 0 || isoforms[i] != isoforms[i - 1]) {
                isoformHits++;
                InterlockedAdd64AndReturnNewValue(&isoforms[i]->mappedCount, 1);
            }
        }
    }

    InterlockedAdd64AndReturnNewValue(&TotalReads, totalReads);
    InterlockedAdd64AndReturnNewValue(&MappedReads, mappedReads);
    InterlockedAdd64AndReturnNewValue(&IsoformHits, isoformHits);
    InterlockedAdd64AndReturnNewValue(&ExonHits, exonHits);

    if (0 == InterlockedDecrementAndReturnNewValue(&nRunningThreads)) {
        SignalSingleWaiterObject(&allThreadsDone);
    }
}

void checkRegions()
{
    GenomeLocation previousEnd = 0;
    for (Region *region = regions->findMin(); NULL != region; region = regions->findFirstGreaterThan(region)) {
        _ASSERT(region->start >= previousEnd);
        for (unsigned i = 0; i < region->nExons; i++) {
            _ASSERT(region->exons[i]->start <= region->start && region->exons[i]->end >= region->end);  // the region is a subset of the exon
        }
        previousEnd = region->end;
    }
}

const char *translateAlternateChromName(const char *inputName)
{
    if (!strcmp(inputName, "GL000191.1")) return "chr1_gl000191_random";
    if (!strcmp(inputName, "GL000192.1")) return "chr1_gl000192_random";
    if (!strcmp(inputName, "GL000193.1")) return "chr4_gl000193_random";
    if (!strcmp(inputName, "GL000194.1")) return "chr4_gl000194_random";
    if (!strcmp(inputName, "GL000195.1")) return "chr7_gl000195_random";
    if (!strcmp(inputName, "GL000205.1")) return "chr17_gl000205_random";
    if (!strcmp(inputName, "GL000209.1")) return "chr19_gl000209_random";
    if (!strcmp(inputName, "GL000211.1")) return "chrUn_gl000211";
    if (!strcmp(inputName, "GL000212.1")) return "chrUn_gl000212";
    if (!strcmp(inputName, "GL000213.1")) return "chrUn_gl000213";
    if (!strcmp(inputName, "GL000214.1")) return "chrUn_gl000214";
    if (!strcmp(inputName, "GL000218.1")) return "chrUn_gl000218";
    if (!strcmp(inputName, "GL000219.1")) return "chrUn_gl000219";
    if (!strcmp(inputName, "GL000220.1")) return "chrUn_gl000220";
    if (!strcmp(inputName, "GL000229.1")) return "chrUn_gl000229";
    if (!strcmp(inputName, "GL000237.1")) return "chrUn_gl000237";
    if (!strcmp(inputName, "GL000241.1")) return "chrUn_gl000241";

    return NULL;
}

int main(int argc, char* argv[])
{
    _int64 start = timeInMillis();

    if (5 != argc) {
        usage();
    }

    const size_t inputBufferLength = 128 * 1024;
    char inputBuffer[inputBufferLength];
    char headerBuffer[inputBufferLength];

    FILE *knownGenesFile = fopen(argv[2], "r");
    if (NULL == knownGenesFile) {
        fprintf(stderr, "Unable to open '%s'\n", argv[2]);
        soft_exit(1);
    }

    static const char *genomeSuffix = "Genome";
    size_t filenameLen = strlen(argv[1]) + 1 + strlen(genomeSuffix) + 1;
    char *fileName = new char[strlen(argv[1]) + 1 + strlen(genomeSuffix) + 1];
    snprintf(fileName, filenameLen, "%s%c%s", argv[1], PATH_SEP, genomeSuffix);
    const Genome *genome = Genome::loadFromFile(fileName, 0);
    if (NULL == genome) {
        fprintf(stderr, "Unable to load genome from file '%s'\n", fileName);
        return -1;
    }
    delete[] fileName;
    fileName = NULL;

    //
    // Open the input file.
    //
    int nThreads = GetNumberOfProcessors();
nThreads = 1; // BJB
    DataSupplier::ThreadCount = nThreads;

    ReaderContext readerContext;
    readerContext.clipping = NoClipping;
    readerContext.defaultReadGroup = "";
    readerContext.genome = genome;
    readerContext.ignoreSecondaryAlignments = true;
    readerContext.ignoreSupplementaryAlignments = true;
    readerContext.header = NULL;
    readerContext.headerLength = 0;
    readerContext.headerBytes = 0;

    if (NULL != strrchr(argv[3], '.') && !_stricmp(strrchr(argv[3], '.'), ".bam")) {
        readSupplierGenerator = BAMReader::createReadSupplierGenerator(argv[3], nThreads, readerContext);
    }
    else {
        readSupplierGenerator = SAMReader::createReadSupplierGenerator(argv[3], nThreads, readerContext);
    }

    regions = new AVLTree<Region, GenomeLocation>();

    //
    // Run through the knownGenes file and build the regions tree.  
    //
    fgets(headerBuffer, inputBufferLength, knownGenesFile);  // Save the header line

    int lineNumber = 2; // Because we just skipped 1

    Isoform *isoforms = NULL;
    Isoform *lastIsoform = NULL;    // This strange thing is a way to do a singly linked list that's in order.

    while (fgets(inputBuffer, inputBufferLength, knownGenesFile)) {
        //
        // Split the line into fields as separated by tabs (which, in C# would be inputLine.Split('\t').  Sigh.)
        //
        int whichField = 0;
        const int nFields = 12;
        char *fields[nFields];


        Isoform *isoform = new Isoform;
        if (NULL == isoforms) {
            isoforms = isoform;
        } else {
            lastIsoform->next = isoform;
        }
        isoform->next = NULL;
        lastIsoform = isoform;


        isoform->knownGeneLine = new char[strlen(inputBuffer) + 1];
        strcpy(isoform->knownGeneLine, inputBuffer);

        for (char *ptr = inputBuffer; NULL != ptr; ptr = strchr(ptr, '\t')) {
            if (0 != whichField) {
                //
                // Turn the tab into a null terminator for the previous field.
                //
                *ptr = '\0';
                ptr++;
            }

            if (whichField >= nFields) {
                fprintf(stderr, "Too few fields on line %d\n", lineNumber);
                soft_exit(1);
            }

            fields[whichField] = ptr;
            whichField++;
        }

        if (whichField != nFields) {
            fprintf(stderr, "Too few fields on line %d\n", lineNumber);
            soft_exit(1);
        }


        isoform->nExons = atoi(fields[kg_exonCount]);

        if (isoform->nExons <= 0) {
            fprintf(stderr, "Isoform %s has zero (or negative) exon count.  Ignoring.  Line %d.\n", fields[kg_name], lineNumber);
            delete isoform;
            continue;
        }

        isoform->exons = new Exon[isoform->nExons];

        GenomeLocation contigBase;
        if (!genome->getLocationOfContig(fields[kg_chrom], &contigBase)) {
            //
            // Try stripping off the leading "chr" if it has one.
            //
            if (strlen(fields[kg_chrom]) < 3 || fields[kg_chrom][0] != 'c' || fields[kg_chrom][1] != 'h' || fields[kg_chrom][2] != 'r' || !genome->getLocationOfContig(fields[kg_chrom] + 3, &contigBase)) {
                //
                // Change chrM into MT
                //
                if (strlen(fields[kg_chrom]) != 4 || fields[kg_chrom][0] != 'c' || fields[kg_chrom][1] != 'h' || fields[kg_chrom][2] != 'r' || fields[kg_chrom][3] != 'M' || !genome->getLocationOfContig("MT", &contigBase)) {
                    //
                    // Try the alternate naming.
                    //
                    const char *alternateName = translateAlternateChromName(fields[kg_chrom]);
                    if (NULL == alternateName || !genome->getLocationOfContig(alternateName, &contigBase)) {
                        //
                        // And ignore long things that didn't otherwise work; they're probably minor haplotypes that aren't in this reference (especially if it's grch37-lite, which specifically has them stripped).
                        //
                        if (strlen(fields[kg_chrom]) < 5) {
                            fprintf(stderr, "Unable to find contig '%s', line %d\n", fields[kg_chrom], lineNumber);
                            soft_exit(1);
                        }
                    }
                }
            }
        }

        char *startPtr = fields[kg_exonStarts];
        char *endPtr = fields[kg_exonEnds];

        for (unsigned i = 0; i < isoform->nExons; i++) {
            if ((char *)1 == startPtr || (char *)1 == endPtr) {
                fprintf(stderr, "Too few exons on line %d\n", lineNumber);
                soft_exit(1);
            }

            isoform->exons[i].start = contigBase + atoi(startPtr) - 1;  // -1 because coordinates are 1-based
            isoform->exons[i].end = contigBase + atoi(endPtr) - 1;
            isoform->exons[i].isoform = isoform;
            if (isoform->exons[i].end <= isoform->exons[i].start) {
                fprintf(stderr, "Zero or negative length exon on line %d\n", lineNumber);
                soft_exit(1);
            }

            startPtr = strchr(startPtr + 1, ',') + 1;
            endPtr = strchr(endPtr + 1, ',') + 1;
        }

        //
        // Now insert the isoform's exons into the region map.
        //
        for (unsigned i = 0; i < isoform->nExons; i++) {
            Exon *exon = &isoform->exons[i];
            GenomeLocation coveredUpTo = exon->start;

            while (coveredUpTo < exon->end) {
                Region *region = regions->findFirstLessThanOrEqualTo(coveredUpTo);

                if (NULL == region || region->end <= coveredUpTo) {
                    //
                    // There's no region covering the first part of what we need.  Create a new one.  It will extend
                    // to the end of the exon or the start of the next region, whichever's first.
                    //
                    Region *nextRegion = regions->findFirstGreaterThan(coveredUpTo);

                    GenomeLocation end;
                    if (NULL != nextRegion && nextRegion->start < exon->end) {
                        end = nextRegion->start;
                    } else {
                        end = exon->end;
                    }
                    Region *newRegion = new Region(coveredUpTo, end, 1, 2);
                    newRegion->exons[0] = exon;

                    regions->insert(newRegion);
                    coveredUpTo = newRegion->end;
                } else {
                    //
                    // The region before us covers our start.
                    //
                    if (region->start < coveredUpTo) {
                        //
                        // It has part hanging off the beginning.  This should only happen if we're the first part of the exon.  Chop it into two, but don't add us to
                        // the second half; instead that'll happen the next time around the loop.
                        //
                        _ASSERT(coveredUpTo == exon->start);
                        Region *newRegion = new Region(region);

                        region->end = coveredUpTo;
                        newRegion->start = coveredUpTo;
                        regions->insert(newRegion);
                    } else {
                        //
                        // It starts exactly where we do.
                        //
                        _ASSERT(region->start == coveredUpTo);
                        if (region->end > exon->end) {
                            //
                            // But ends beyond where we do.  Split it at our and and add ourself to the first half.
                            //
                            Region *newRegion = new Region(region);
                            newRegion->start = exon->end;
                            region->end = exon->end;
                            regions->insert(newRegion);
                            region->addExon(exon);
                        } else {
                            //
                            // Just add ourself to it.  It ends before us or exactly where we end.
                            //
                            region->addExon(exon);
                        }
                    }
                    coveredUpTo = region->end;
                } // region doesn't cover our start


            } // while we have space in this exon
        } // for each exon

        lineNumber++;
    } // for each isoform (i.e., non-header line in known genes file)

#ifdef _DEBUG
    checkRegions();
#endif // _DEBUG

    //
    // Now run through the tree and link the regions together.
    //
    Region *prevRegion = NULL;
    for (Region *region = regions->findMin(); NULL != region; region = regions->findFirstGreaterThan(region)) {
        if (NULL != prevRegion) {
            prevRegion->next = region;
        }
        prevRegion = region;
    }

    //
    // Start up the reader threads.
    //
    CreateSingleWaiterObject(&allThreadsDone);

    nRunningThreads = nThreads;

    for (int i = 0; i < nThreads; i++) {
        StartNewThread(ReaderThread, NULL);
    }

    WaitForSingleWaiterObject(&allThreadsDone);

    if (0 == MappedReads) {
        WriteErrorMessage("No reads mapped.  Something is wrong, not writing output.\n");
        soft_exit(1);
    }

    FILE *outputFile = fopen(argv[4], "w");
    if (NULL == outputFile) {
        fprintf(stderr, "Unable to open output file '%s'\n", argv[4]);
        soft_exit(1);
    }

    char *newLine = strchr(headerBuffer, '\n');
    if (NULL != newLine) {
        *newLine = '\0';
    }
    fprintf(outputFile, "%s\tExonHits(%lld)\tIsoformHits(%lld)\tFracIsoformHits\tv1.0\n", headerBuffer, ExonHits, IsoformHits);

    for (Isoform *isoform = isoforms; NULL != isoform; isoform = isoform->next) {
        if (NULL != (newLine = strchr(isoform->knownGeneLine, '\n'))) {
            *newLine = '\0';
        }

        fprintf(outputFile, "%s\t", isoform->knownGeneLine);
        for (unsigned i = 0; i < isoform->nExons; i++) {
            fprintf(outputFile, "%lld,", isoform->exons[i].mappedCount);
        }
        fprintf(outputFile, "\t%lld\t%0.9f\n", isoform->mappedCount, (double)isoform->mappedCount / (double)IsoformHits);
    }

    const int numberMaxWidth = 20;
    char readCount[numberMaxWidth];
    char mappedReadCount[numberMaxWidth];
    char isoformHitCount[numberMaxWidth];
    char exonHitCount[numberMaxWidth];
    char time[numberMaxWidth];

    extern char *FormatUIntWithCommas(_uint64 val, char *outputBuffer, size_t outputBufferSize);

    printf("%s reads (%s mapped), %s isoform hits, %s exon hits in %ss\n",
        FormatUIntWithCommas(TotalReads, readCount, numberMaxWidth),
        FormatUIntWithCommas(MappedReads, mappedReadCount, numberMaxWidth),
        FormatUIntWithCommas(IsoformHits, isoformHitCount, numberMaxWidth),
        FormatUIntWithCommas(ExonHits, exonHitCount, numberMaxWidth),
        FormatUIntWithCommas((timeInMillis() - start + 500) / 1000, time, numberMaxWidth));

	return 0;
}

