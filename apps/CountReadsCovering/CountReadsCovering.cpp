// CountReadsCovering.cpp : Defines the entry point for the console application.
//

#include "stdafx.h"
#include "exit.h"
#include "Genome.h"


void usage()
{
    fprintf(stderr, "usage: CountReadsCovering index listOfLocationsToCheck inputFile\n");
    soft_exit(1);
}

int main(int argc, char* argv[])
{
    if (4 != argc) {
        usage();
    }

    FILE *listOfLocationsFile = fopen(argv[2], "r");
    if (NULL == listOfLocationsFile) {
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
    // First, count the lines in the listOfLocations so we know how much memory to allocate.
    //
    int nLocationsToCheck = 0;
    const int inputBufferSize = 1024;
    char inputBuffer[inputBufferSize];

    while (NULL != fgets(inputBuffer, inputBufferSize, listOfLocationsFile)) {
        nLocationsToCheck++;
    }

    GenomeLocation *locationsToCheck = new GenomeLocation[nLocationsToCheck];

    rewind(listOfLocationsFile);
    for (int i = 0; i < nLocationsToCheck; i++) {
        if (NULL == fgets(inputBuffer, inputBufferSize, listOfLocationsFile)) {
            fprintf(stderr, "Error reading listOfLocations file second time\n");
            soft_exit(1);
        }

        char *tab = strchr(inputBuffer, '\t');
        if (NULL == tab) {
            fprintf(stderr, "Incorrect format for listOfLocations like '%s'; should be contig offset\n", inputBuffer);
            soft_exit(1);
        }

        *tab = 0;   // Make two strings.
        GenomeLocation contigBase;
        if (!genome->getLocationOfContig(inputBuffer, &contigBase)) {
            fprintf(stderr, "Unable to find contig '%s', listOfLocations line %n\n", inputBuffer, i + 1);
            soft_exit(1);
        }

        int offset = atoi(tab + 1);
        if (0 == offset) {
            fprintf(stderr, "Failed to parse offset in listOfLocations line %n\n", i);
            soft_exit(1);
        }

        locationsToCheck[i] = contigBase + offset - 1;  // -1 because we're 0 based and S/BAM is 1-based.
    }

    now sort locationsToCheck.


	return 0;
}

