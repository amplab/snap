// RemoveSlashesFromReadNames.cpp : Remove the /1 and /2 from paired-end reads in a sam file
//

#include <iostream>
#include <windows.h>


void usage()
{
    fprintf(stderr, "usage: RemoveSlashesFromReadNames input.sam output.sam\n");
    exit(1);
}
int main(int argc, const char **argv)
{
    if (argc != 3) {
        usage();
    }

    FILE* inputFile;
    fopen_s(&inputFile, argv[1], "r");
    if (NULL == inputFile) {
        fprintf(stderr, "Unable to open %s for input\n", argv[1]);
        exit(1);
    }

    FILE* outputFile;
    fopen_s(&outputFile, argv[2], "w");
    if (NULL == outputFile) {
        fprintf(stderr, "Unable to open %s for output\n", argv[2]);
        exit(1);
    }

    int bufferSize = 1024 * 1024; // should be plenty
    char* buffer = new char[bufferSize];

    while (NULL != fgets(buffer, bufferSize - 1, inputFile)) {
        if (strchr(buffer, '\n') == NULL) {
            fprintf(stderr, "Found a line bigger than buffer (!)\n");
            exit(1);
        }

        if (buffer[0] == '@') {
            //
            // Header line.
            //
            fprintf(outputFile, "%s", buffer);
        } else {
            char *firstTab = strchr(buffer, '\t');
            if (firstTab == NULL) {
                fprintf(stderr, "Malformed SAM line (no tabs): %s\n", buffer);
                exit(1);
            }

            if (firstTab - buffer > 2 && (firstTab[-1] == '1' || firstTab[-1] == '2') && firstTab[-2] == '/') {
                //
                // Write out the whole line except for the /1 or /2
                //
                fwrite(buffer, 1, firstTab - buffer - 2, outputFile);
                fwrite(firstTab, 1, strlen(buffer) - (firstTab - buffer), outputFile);
            } else {
                //
                // Nothing to cut.  Write the whole line.
                //
                fwrite(buffer, 1, strlen(buffer), outputFile);
            }

        }
    } // while we have input

    fclose(inputFile);
    fclose(outputFile);

} // main

