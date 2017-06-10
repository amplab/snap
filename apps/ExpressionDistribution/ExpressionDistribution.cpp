// ExpressionDistribution.cpp : Defines the entry point for the console application.
//

#include "stdafx.h"
#include "DataReader.h"

struct ListOfFilenames {
    char *filename;
	char *mappedBaseCountFilename;
    ListOfFilenames *next;
};

struct ListOfChromosomes;

struct Disease {
    Disease() : diseaseName(NULL), allCountFiles(NULL), listOfChromosomes(NULL), nTotalFiles(0), nCorruptFiles(0), nMissingFiles(0) {}
    //~Disease();
    char *diseaseName;
    ListOfFilenames *allCountFiles;
    ListOfChromosomes *listOfChromosomes;

    int nTotalFiles;
    int nCorruptFiles;
    int nMissingFiles;
};

struct StdDevState {
    StdDevState() : sum(0), sumOfSquares(0), n(0) {}
    double sum;
    double sumOfSquares;
    _int64 n;

    double getStdDev() {
        if (0 == n) return 0;

        return sqrt(n * sumOfSquares - sum * sum) / (double)n;
    }

    double getMean() {
        if (0 == n) return 0;
        return sum / (double)n;
    }

    void recordObservation(double observation)
    {
        sum += observation;
        sumOfSquares += observation * observation;
        n++;
    }
};

struct Chromosome {
    Chromosome(int _size)
    {
        size = _size;
        stdDevState = new StdDevState*[size+1]; // it's 1-based, hence +1
        memset(stdDevState, 0, sizeof(*stdDevState) * (size + 1));
    }

    ~Chromosome()
    {
        for (int i = 0; i <= size; i++) {
            if (stdDevState[i]) {
                delete stdDevState[i];
            }
        }
        delete[] stdDevState;
    }

    char *name;
    int size;
    StdDevState **stdDevState;
};

struct ListOfChromosomes {
    ListOfChromosomes() : chromosome(NULL), next(NULL) {};
    ~ListOfChromosomes()
    {
        delete chromosome;
        delete next;    // Tail recursive.
    }
    Chromosome *chromosome;
    ListOfChromosomes * next;
};

class CompressedLineReader {
public:
    CompressedLineReader(char *filename);
    ~CompressedLineReader();
    char *getNextLine(int *o_lineLength);
    bool getNextLine(char *lineBuffer, size_t lineBufferSize);

private:
    int amountToAdvance;
    bool readyForNextBatch;
    DataReader *reader;
};

CompressedLineReader::CompressedLineReader(char *filename)
{
    DataSupplier *gzipSupplier = DataSupplier::Gzip(DataSupplier::Default);
    reader = gzipSupplier->getDataReader(10, 100 * 1024, 5, 10 * 1024 * 1024);

    reader->init(filename);
    _int64 headerIOSize = 1024;
    reader->readHeader(&headerIOSize);
    reader->reinit(0, 0);

    amountToAdvance = 0;
    readyForNextBatch = false;
}

CompressedLineReader::~CompressedLineReader()
{
    delete reader;
}

char *
CompressedLineReader::getNextLine(int *o_lineLength)
{
    reader->advance(amountToAdvance);
    if (readyForNextBatch) {
        reader->nextBatch();
    }

    char *dataBuffer;
    _int64 validBytes;
    _int64 bytesToStart;
    if (!reader->getData(&dataBuffer, &validBytes, &bytesToStart)) {
        return NULL;    // EOF
    }

    int bytesToEOL;
    int EOLBytes;   // The count of crlf bytes

    for (bytesToEOL = 0; bytesToEOL < validBytes && dataBuffer[bytesToEOL] != '\r' && dataBuffer[bytesToEOL] != '\n'; bytesToEOL++) {
    }

    if (bytesToEOL < validBytes && dataBuffer[bytesToEOL] == '\r') {
        EOLBytes = 1;
    }
    else {
        EOLBytes = 0;
    }

    if (dataBuffer[bytesToEOL + EOLBytes] == '\n') {
        EOLBytes++;
    }

    amountToAdvance = bytesToEOL + EOLBytes;
    readyForNextBatch = amountToAdvance >= bytesToStart;

    *o_lineLength = bytesToEOL;
    return dataBuffer; 
}

bool
CompressedLineReader::getNextLine(char *lineBuffer, size_t lineBufferSize)
{
    int lineSize;
    char *data = getNextLine(&lineSize);

    if (NULL == data) {
        return false;
    }

    if (lineSize >= lineBufferSize) {
        return false;
    }

    memcpy(lineBuffer, data, lineSize);
    lineBuffer[lineSize] = '\0';

    return true;
}

void ProcessAllcountFile(char *filename, char *mappedBaseCountFilename, Disease *disease)
{
    static _int64 nAllocated = 0;
    static _int64 nHit = 0;

    int lineNumber = 0;

	_int64 mappedBaseCount;

	FILE *mappedBaseCountFile = fopen(mappedBaseCountFilename, "r");
	if (NULL == mappedBaseCountFile) {
		fprintf(stderr, "Unable to open mapped base count file %s\n", mappedBaseCountFilename);
		disease->nCorruptFiles++;
		return;	// Don't goto done, because we haven't yet created reader
	}

	if (1 != fscanf(mappedBaseCountFile, "%lld\t%*s", &mappedBaseCount)) {
		fprintf(stderr, "Unable to parse mapped base count file %s\n", mappedBaseCountFilename);
		disease->nCorruptFiles++;
		fclose(mappedBaseCountFile);
		return;
	}

	fclose(mappedBaseCountFile);

	if (mappedBaseCount < 1000000) {
		fprintf(stderr, "Mapped base count from file %s does not pass sanity check: %lld\n", mappedBaseCountFilename, mappedBaseCount);
		disease->nCorruptFiles++;
		return;
	}

    HANDLE hFile = CreateFile(filename, GENERIC_READ, FILE_SHARE_READ, NULL, OPEN_EXISTING, FILE_ATTRIBUTE_NORMAL, NULL);
    if (INVALID_HANDLE_VALUE == hFile) {
        disease->nMissingFiles++;
        fprintf(stderr, "File %s is missing or can't be opened,  %d\n", filename, GetLastError());
        return;
    }
    CloseHandle(hFile);

    CompressedLineReader *reader = new CompressedLineReader(filename);

    const int lineBufferSize = 1000;
    char lineBuffer[lineBufferSize];

    if (!reader->getNextLine(lineBuffer, lineBufferSize)) {
        fprintf(stderr, "Missing header in file %s\n", filename);
        disease->nCorruptFiles++;
        goto done;
    }
    lineNumber++;

    //
    // Second line is the total and mapped read count
    //
    if (!reader->getNextLine(lineBuffer, lineBufferSize)) {
        fprintf(stderr, "Missing mapped read count in file %s\n", filename);
        disease->nCorruptFiles++;
        goto done;
    }
    lineNumber++;

    _int64 nMappedReads, nLowQualityReads, nMitochondrialReads, nTotalReads;

    if (4 != sscanf(lineBuffer, "%lld mapped high quality nuclear reads, %lld low quality reads, %lld mitochondrial reads, %lld total reads", &nMappedReads, &nLowQualityReads, &nMitochondrialReads, &nTotalReads)) {
        fprintf(stderr, "Unable to parse count line '%s' in %s\n", lineBuffer, filename);
        disease->nCorruptFiles++;
        goto done;
    }

    if (nMappedReads == 0 || nMappedReads + nLowQualityReads + nMitochondrialReads > nTotalReads || (nMappedReads + nLowQualityReads + nMitochondrialReads) * 4 < nTotalReads) {
        fprintf(stderr, "Read count '%s' doesn't pass sanity check for file %s\n", lineBuffer, filename);
        disease->nCorruptFiles++;
        goto done;
    }

    if (!reader->getNextLine(lineBuffer, lineBufferSize) || !reader->getNextLine(lineBuffer, lineBufferSize)) { // Blank line & NumContigs line
        fprintf(stderr, "Error reading NumContigs line in %s\n", filename);
        disease->nCorruptFiles++;
        goto done;
    }
    lineNumber += 2;

    int numContigs;
    if (1 != sscanf(lineBuffer, "NumContigs: %d\n", &numContigs) || numContigs < 1 || numContigs > 1000) {
        fprintf(stderr, "Unable to parse num contigs line '%s' from %s\n", lineBuffer, filename);
        disease->nCorruptFiles++;
        goto done;
    }
    
    if (!reader->getNextLine(lineBuffer, lineBufferSize)) {
        fprintf(stderr, "Unable to read ContigName header line in %s\n", filename);
        disease->nCorruptFiles++;
        goto done;
    }
    lineNumber++;

    lineBuffer[lineBufferSize - 1] = '\0';

    for (int i = 0; i < numContigs; i++) {
        if (!reader->getNextLine(lineBuffer, lineBufferSize)) {
            fprintf(stderr, "File truncated reading contigs, %s\n", filename);
            disease->nCorruptFiles++;
            goto done;
        }
        lineNumber++;

        char chromosomeName[lineBufferSize];
        int chromosomeSize;

        if (2 != sscanf(lineBuffer, "%[^\t]\t%d", chromosomeName, &chromosomeSize) || chromosomeSize <= 0) {
            fprintf(stderr, "Unable to parse contig line '%s' in %s\n", lineBuffer, filename);
            disease->nCorruptFiles++;
            goto done;
        }

        if (chromosomeName[0] != 'c' || chromosomeName[1] != 'h' || chromosomeName[2] != 'r') {
            //
            // Prepend "chr"
            //
            sprintf(lineBuffer, "chr%s", chromosomeName);
            strcpy(chromosomeName, lineBuffer);
        }

        //
        // See if we already have it.
        //
        ListOfChromosomes *chromosomeEntry;
        for (chromosomeEntry = disease->listOfChromosomes; NULL != chromosomeEntry; chromosomeEntry = chromosomeEntry->next) {
            if (!strcmp(chromosomeEntry->chromosome->name, chromosomeName)) {
                if (chromosomeEntry->chromosome->size < chromosomeSize) {
                    fprintf(stderr, "Need to expand chromosome %s, %d > %d, %s\n", chromosomeName, chromosomeSize, chromosomeEntry->chromosome->size, filename);
                    disease->nCorruptFiles++;   // Not really corrupt, more like I have to write this code.  Whatever.
                    goto done;
                }
                break;
            }
        }

        if (NULL == chromosomeEntry) {
            chromosomeEntry = new ListOfChromosomes;
            chromosomeEntry->chromosome = new Chromosome(chromosomeSize);
            chromosomeEntry->chromosome->name = new char[strlen(chromosomeName) + 1];
            strcpy(chromosomeEntry->chromosome->name, chromosomeName);

            chromosomeEntry->next = disease->listOfChromosomes;
            disease->listOfChromosomes = chromosomeEntry;
        }
    }

    //
    // Finally, process the actual meat of the file.  The format is that a chromosome starts with >chromosomeName, and then each line is of one of three formats:
    //  offset\tmappedCount
    //  mappedCount
    //  xRepeatCount
    //
    // And the last line of the file contains "**done**".
    //
    // The first one means that that offset in the chromosome has that many reads mapped to it.  The second means that the next address after the one most recently referenced has
    // that mapped count, and the final one (which has a literal "x) means that the next RepeatCount offsets have the same mappedCount.
    //

    Chromosome *chromosome = NULL;
    int nextOffset = 0;
    int currentMappedCount = 0;
    while (reader->getNextLine(lineBuffer, lineBufferSize)) {
        lineNumber++;

        if (lineBuffer[0] == '>') {
            char chromosomeName[lineBufferSize];
            if (lineBuffer[1] != 'c' || lineBuffer[2] != 'h' || lineBuffer[3] != 'r') {
                strcpy(chromosomeName, "chr");
                strcat(chromosomeName, lineBuffer + 1);
            }
            else {
                strcpy(chromosomeName, lineBuffer + 1);
            }

            chromosome = NULL;
            for (ListOfChromosomes *chromosomeEntry = disease->listOfChromosomes; NULL != chromosomeEntry; chromosomeEntry = chromosomeEntry->next) {
                if (!strcmp(chromosomeEntry->chromosome->name, chromosomeName)) {
                    chromosome = chromosomeEntry->chromosome;
                    break;
                }
            }

            if (NULL == chromosome) {
                fprintf(stderr, "Unable to find chromosome %s in file %s\n", chromosomeName, filename);
                disease->nCorruptFiles++;
                goto done;
            }
            //printf("%s\n", chromosome->name); fflush(stdout);
            nextOffset = 0;
            currentMappedCount = 0;
            continue;
        }

        if (NULL == chromosome) {
            fprintf(stderr, "Body of file %s doesn't start with chromosome line, instead it's %s\n", filename, lineBuffer);
            disease->nCorruptFiles++;
            goto done;
        }

        if (!strcmp(lineBuffer, "**done**")) {
            //
            // Finished successfully.
            //
            goto done;
        }

       char *tab = strchr(lineBuffer, '\t');
       int repeatCount;
       if (NULL != tab) {
           //
           // It's format 1: address\tmappedReadCount
           //
           int offset;
           if (2 != sscanf(lineBuffer, "%x\t%x", &offset, &currentMappedCount)) {
               fprintf(stderr, "Error parsing data line %s in file %s\n", lineBuffer, filename);
               disease->nCorruptFiles++;
               goto done;
           }

           if (offset <= nextOffset) {
               fprintf(stderr, "Saw unexpectedly low offset %d in file %s chromosome %s (expected more than %d)\n", offset, filename, chromosome->name, nextOffset);
               disease->nCorruptFiles++;
               goto done;
           }
           repeatCount = 1;
           nextOffset = offset;
       }
       else if (lineBuffer[0] == 'x') {
           if (1 != sscanf(lineBuffer, "x%x", &repeatCount) || repeatCount < 2) {
               fprintf(stderr, "Unable to parse line %s (or too small repeatCount), file %s\n", lineBuffer, filename);
               disease->nCorruptFiles++;
               goto done;
           }
           repeatCount -= 1;    // Because the first one is already done.
       }
       else {
           if (1 != sscanf(lineBuffer, "%x", &currentMappedCount) || currentMappedCount < 1) {
               fprintf(stderr, "Couldn't parse line %s (or too small mapped count) in file %s\n", lineBuffer, filename);
               disease->nCorruptFiles++;
               goto done;
           }
           repeatCount = 1;
       }

       if (0 == nextOffset) {
           fprintf(stderr, "Chromosome didn't start with address, chromosome %s file %s\n", chromosome->name, filename);
           disease->nCorruptFiles++;
           goto done;
       }

       for (int i = 0; i < repeatCount; i++) {
           if (nextOffset > chromosome->size) {
               fprintf(stderr, "Mapped area ran off end of chromosome %s in file %s, offset %d\n", chromosome->name, filename, nextOffset);
               disease->nCorruptFiles++;
               goto done;
           }

           if (chromosome->stdDevState[nextOffset] == NULL) {
               chromosome->stdDevState[nextOffset] = new StdDevState();
               if ((++nAllocated) % 200000000 == 0) {
                   //printf("%lldM allocated\n", nAllocated / 1000000);
               }
           }
           else {
               if ((++nHit) % 200000000 == 0) {
                   //printf("%lldM hit\n", nHit / 1000000);
               }
           }
           chromosome->stdDevState[nextOffset]->recordObservation((double)currentMappedCount / (double)mappedBaseCount);
           nextOffset++;
       } // for each repeat count
    } // While we have an input line.

    //
    // If we get here, we fell off the end of the file without seeing **done**, so it's been truncated.
    //
    lineBuffer[lineBufferSize - 1] = '\0';
    fprintf(stderr, "Truncated file %s, last line (%d): %s\n", filename, lineNumber, lineBuffer);
    disease->nCorruptFiles++;


done:
    delete reader;
}


void ProcessDisease(Disease *disease)
{
    printf("Disease: %s\n", disease->diseaseName);
    disease->listOfChromosomes = NULL;
    int nFiles = 0;
    _int64 start = timeInMillis();

    for (ListOfFilenames *allcountFile = disease->allCountFiles; NULL != allcountFile; allcountFile = allcountFile->next) {
        nFiles++;
        //printf("%d: %s\t(%llds)\n", nFiles, allcountFile->filename, (timeInMillis() - start) / 1000); fflush(stdout);
        ProcessAllcountFile(allcountFile->filename, allcountFile->mappedBaseCountFilename, disease);
    }
}


void usage()
{
    fprintf(stderr, "usage: ExpressionDistribution casesFilename expressionFilesDirectory projectColumn tumorRNAAllcountColumn tumorRNAMappedBaseCountColumn <diseaseName>\n");
    exit(1);
}

int main(int argc, char* argv[])
{

    FILE *casesFile = NULL;
    
    int retryCount = 0;
    
    if (argc != 6 && argc != 7) {
        usage();
    }

    int projectColumn;
    int tumorRNAAllcountColumn;
	int tumorRNAMappedBaseCountColumn;

    projectColumn = atoi(argv[3]);
    tumorRNAAllcountColumn = atoi(argv[4]);
	tumorRNAMappedBaseCountColumn = atoi(argv[5]);

    if (projectColumn < 0|| tumorRNAAllcountColumn < 0 || tumorRNAMappedBaseCountColumn < 0) {
        fprintf(stderr, "Invalid project, allcount or mapped base count columns.\n");
        return -1;
    }

    const size_t maxColumn = 100;

    int highestInterestingColumnNumber = __max(projectColumn, __max(tumorRNAAllcountColumn, tumorRNAMappedBaseCountColumn));

    if (highestInterestingColumnNumber >= maxColumn) {
        fprintf(stderr, "One or both of the column numbers are bigger than permitted.\n");
        return -1;
    }

    while (NULL == (casesFile = fopen(argv[1], "r"))) {
        if (retryCount++ == 100) {
            fprintf(stderr, "Gave up opening cases file ('%s') after 100 tries\n", argv[1]);
            return -1;
        }
        Sleep(1000);
    }

    const int maxDiseases = 100;
    int nDiseases = 0;
    Disease *diseases = new Disease[maxDiseases];

    const size_t bufferSize = 10000;
    char inputBuffer[bufferSize];

    fgets(inputBuffer, bufferSize, casesFile);    // Skip the header line

    bool seenDone = false;

    while (fgets(inputBuffer, bufferSize, casesFile)) {

        //
        // fgets leaves the newline around.  Kill it if it's there.
        //
        char *newline = strchr(inputBuffer, '\n');
        if (NULL != newline) {
            *newline = '\0';
        }

        if (seenDone) {
            fprintf(stderr, "cases file continues beyond **done**.\n");
            return -1;
        }

        if (!strcmp(inputBuffer, "**done**")) {
            seenDone = true;
            continue;
        }

        // format is tab separated columns.  The project ID is of the form TCGA-DISEASE
        char *fields[maxColumn];

        char *nextField = inputBuffer;
        for (int i = 0; i <= highestInterestingColumnNumber; i++) {
            if (nextField == NULL) {
                fprintf(stderr, "Truncated cases file, line has too few fields.\n");
                return -1;
            }

            fields[i] = nextField;

            char *tab = strchr(nextField, '\t');
            if (NULL == tab) {
                nextField = NULL;
            } else {
                *tab = '\0';
                nextField = tab + 1;
            }
        }

        //
        // Null terminate the highest column string (if it isn't already).
        //
        char *charToCheck;
        for (charToCheck = fields[highestInterestingColumnNumber]; *charToCheck != '\t' && *charToCheck != '\0'; charToCheck++) {
            // This loop body intentionally left blank
        }
        *charToCheck = '\0';

        char *diseaseName;
        if (strlen(fields[projectColumn]) < 6 || fields[projectColumn][0] != 'T' || fields[projectColumn][1] != 'C' || fields[projectColumn][2] != 'G' || fields[projectColumn][3] != 'A' || fields[projectColumn][4] != '-') {
            fprintf(stderr, "Invalid project name (%s), doesn't start with TCGA-", fields[projectColumn]);
            return -1;
        }

        diseaseName = fields[projectColumn] + 5;
        for (char *charToLower = diseaseName; *charToLower != '\0'; charToLower++) {
            *charToLower = tolower(*charToLower);
        }

        Disease *disease = NULL;
        for (int i = 0; i < nDiseases; i++) {
            if (!strcmp(diseases[i].diseaseName, diseaseName)) {
                disease = &diseases[i];
                break;
            }
        }

        if (NULL == disease) {
            if (nDiseases == maxDiseases) {
                fprintf(stderr, "Ran out of diseases!\n");
                return -1;
            }
            disease = &diseases[nDiseases];
            nDiseases++;

            disease->diseaseName = new char[strlen(diseaseName) + 1];
            strcpy(disease->diseaseName, diseaseName);
        }


        char *tumorAllcountsFileField = fields[tumorRNAAllcountColumn];
		char *tumorRNAMappedBaseCountFileField = fields[tumorRNAMappedBaseCountColumn];

        if (strcmp(tumorAllcountsFileField, "") && strcmp(tumorRNAMappedBaseCountFileField, "")) {
            //
            // It's not null, add it to the disease.
            //
            ListOfFilenames *entry = new ListOfFilenames;
            entry->filename = new char[strlen(tumorAllcountsFileField) + 1];
            strcpy(entry->filename, tumorAllcountsFileField);
			entry->mappedBaseCountFilename = new char[strlen(tumorRNAMappedBaseCountFileField) + 1];
			strcpy(entry->mappedBaseCountFilename, tumorRNAMappedBaseCountFileField);

            entry->next = disease->allCountFiles;
            disease->allCountFiles = entry;
            disease->nTotalFiles++;
        }
    } // while there's another line in the cases file

    if (!seenDone) {
        fprintf(stderr, "Truncated cases file: does not end with **done**\n");
        return -1;
    }

    fclose(casesFile);
    casesFile = NULL;

    _int64 start = timeInMillis();
    for (int i = 0; i < nDiseases; i++) {
        Disease *disease = &diseases[i];
        if (argc == 7 && strcmp(argv[6], disease->diseaseName)) {
            continue;
        }
        ProcessDisease(disease);
        printf("%s has %d total, %d missing and %d corrupt files\n", disease->diseaseName, disease->nTotalFiles, disease->nMissingFiles, disease->nCorruptFiles);
        int nGoodFiles = disease->nTotalFiles - disease->nMissingFiles - disease->nCorruptFiles;

        char fileNameBuffer[1000];
        sprintf(fileNameBuffer, "%s\\expression_%s", argv[2], disease->diseaseName);
        FILE *outputFile = fopen(fileNameBuffer, "w");
        if (NULL == outputFile) {
            fprintf(stderr, "Unable to open output file '%s'\n", argv[2]);
        } else {
            for (ListOfChromosomes *entry = disease->listOfChromosomes; NULL != entry; entry = entry->next) {
                Chromosome *chromosome = entry->chromosome;
                fprintf(outputFile, "%s\n", chromosome->name);
                for (int i = 0; i <= chromosome->size; i++) {
                    if (chromosome->stdDevState[i] != 0 && chromosome->stdDevState[i]->n * 10 >= nGoodFiles) {  // Only bother with ones where there's at least 10% with at least one read mapped
                        _int64 oldN = chromosome->stdDevState[i]->n;
                        chromosome->stdDevState[i]->n = nGoodFiles; // Account for the ones with zero reads mapped here
                        fprintf(outputFile, "%d\t%lld\t%e\t%e\n", i, oldN, chromosome->stdDevState[i]->getMean(), chromosome->stdDevState[i]->getStdDev());
                    }
                }
            }
            fclose(outputFile);
        }

        delete disease->listOfChromosomes;  // This deletes the whole thing in one fell swoop.
        disease->listOfChromosomes = NULL;

    }
    fprintf(stderr, "Took %lldm\n", (timeInMillis() - start + 30000)/60000);

	return 0;
}

