// ExpressionDistribution.cpp : Defines the entry point for the console application.
//

#include "stdafx.h"
#include "DataReader.h"

struct ListOfFilenames {
    char *filename;
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

void ProcessAllcountFile(char *filename, Disease *disease)
{
    static _int64 nAllocated = 0;
    static _int64 nHit = 0;

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

    //
    // Second line is the total and mapped read count
    //
    if (!reader->getNextLine(lineBuffer, lineBufferSize)) {
        fprintf(stderr, "Missing mapped read count in file %s\n", filename);
        disease->nCorruptFiles++;
        goto done;
    }

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

    for (int i = 0; i < numContigs; i++) {
        if (!reader->getNextLine(lineBuffer, lineBufferSize)) {
            fprintf(stderr, "File truncated reading contigs, %s\n", filename);
            disease->nCorruptFiles++;
            goto done;
        }

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
           chromosome->stdDevState[nextOffset]->recordObservation((double)currentMappedCount / (double)nMappedReads);
           nextOffset++;
       }


    } // While we have an input line.

    //
    // If we get here, we fell off the end of the file without seeing **done**, so it's been truncated.
    //
    fprintf(stderr, "Truncated file %s\n", filename);
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
        ProcessAllcountFile(allcountFile->filename, disease);
    }
}



int main(int argc, char* argv[])
{

    FILE *experimentsFile = NULL;
    
    int retryCount = 0;
    while (NULL == (experimentsFile = fopen("\\\\bolosky\\f$\\temp\\expression\\experiments.txt", "r"))) {
        if (retryCount++ == 100) {
            fprintf(stderr, "Gave up opening experiments file after 100 tries\n");
            return 1;
        }
        Sleep(1000);
    }

    const int maxDiseases = 100;
    int nDiseases = 0;
    Disease *diseases = new Disease[maxDiseases];

    const size_t bufferSize = 10000;
    char inputBuffer[bufferSize];

    fgets(inputBuffer, bufferSize, experimentsFile);    // Skip the header line

    while (fgets(inputBuffer, bufferSize, experimentsFile)) {
        // format is tab separated columns.  Column 0 is disease_abbr, column 16 tumorAllcountFile

        const int diseaseNameBufferSize = 100;
        char diseaseNameBuffer[diseaseNameBufferSize];
        memcpy(diseaseNameBuffer, inputBuffer, __min(diseaseNameBufferSize, strlen(inputBuffer)));
        diseaseNameBuffer[diseaseNameBufferSize - 1] = '\0';
        char *tab = strchr(diseaseNameBuffer, '\t');
        if (NULL == tab) {
            fprintf(stderr, "Unparsable experiment line %s\n", inputBuffer);
            return 1;
        }

        *tab = '\0';

        if (!strcmp(diseaseNameBuffer, "ov")) {
            char *reference = tab + 1;
            tab = strchr(reference, '\t');
            if (NULL == tab) {
                fprintf(stderr, "Unparsable ov line %s", inputBuffer);
                return 1;
            }
            *tab = '\0';
            if (!strcmp(reference, "ncbi36_bccagsc_variant")) {
                strcpy(diseaseNameBuffer, "ov_hg18");
            }
        }

        Disease *disease = NULL;
        for (int i = 0; i < nDiseases; i++) {
            if (!strcmp(diseases[i].diseaseName, diseaseNameBuffer)) {
                disease = &diseases[i];
                break;
            }
        }

        if (NULL == disease) {
            if (nDiseases == maxDiseases) {
                fprintf(stderr, "Ran out of diseases!\n");
                return 1;
            }
            disease = &diseases[nDiseases];
            nDiseases++;

            disease->diseaseName = new char[strlen(diseaseNameBuffer) + 1];
            strcpy(disease->diseaseName, diseaseNameBuffer);
        }


        char *tumorAllcountsFileField = inputBuffer;
        for (int i = 1; i <= 16; i++) {
            tumorAllcountsFileField = strchr(tumorAllcountsFileField, '\t');
            if (NULL == tumorAllcountsFileField) {
                fprintf(stderr, "Couldn't find enough tabs in experiment line %s\n", inputBuffer);
                return 1;
            }
            tumorAllcountsFileField++; // Skip over the tab
        }

        char *nextTab = strchr(tumorAllcountsFileField, '\t');
        if (NULL == nextTab) {
            fprintf(stderr, "Can't find trailing tab in experiments like %s\n", inputBuffer);
            return 1;
        }
        *nextTab = '\0';

        if (strcmp(tumorAllcountsFileField, "")) {
            //
            // It's not null, add it to the disease.
            //
            ListOfFilenames *entry = new ListOfFilenames;
            entry->filename = new char[strlen(tumorAllcountsFileField) + 1];
            strcpy(entry->filename, tumorAllcountsFileField);
            entry->next = disease->allCountFiles;
            disease->allCountFiles = entry;
            disease->nTotalFiles++;
        }
    } // while there's another line in the experiments file
    fclose(experimentsFile);
    experimentsFile = NULL;

    _int64 start = timeInMillis();
    for (int i = 0; i < nDiseases; i++) {
        Disease *disease = &diseases[i];
        if (argc == 2 && strcmp(argv[1], disease->diseaseName)) {
            continue;
        }
        ProcessDisease(disease);
        printf("%s has %d total, %d missing and %d corrupt files\n", disease->diseaseName, disease->nTotalFiles, disease->nMissingFiles, disease->nCorruptFiles);
        int nGoodFiles = disease->nTotalFiles - disease->nMissingFiles - disease->nCorruptFiles;

        char fileNameBuffer[1000];
        sprintf(fileNameBuffer, "\\\\msr-genomics-1\\d$\\temp\\expression_%s", disease->diseaseName);
        FILE *outputFile = fopen(fileNameBuffer, "w");
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

        delete disease->listOfChromosomes;  // This deletes the whole thing in one fell swoop.
        disease->listOfChromosomes = NULL;

    }
    fprintf(stderr, "Took %lldms\n", timeInMillis() - start);



	return 0;
}

