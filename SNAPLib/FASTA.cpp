/*++

Module Name:

    FASTA.cpp

Abstract:

    FASTA reader

Authors:

    Bill Bolosky, August, 2011

Environment:

    User mode service.

Revision History:

    Adapted from Matei Zaharia's Scala implementation.

--*/

#include "stdafx.h"
#include "Compat.h"
#include "FASTA.h"


namespace {

_uint64 getFileSize(const char *fileName)
{
#ifdef _MSC_VER
    // On Windows, use GetFileSizeEx
    HANDLE hFile = CreateFile(fileName,GENERIC_READ,FILE_SHARE_READ,NULL,OPEN_EXISTING,NULL,NULL);

    if (INVALID_HANDLE_VALUE == hFile) {
        fprintf(stderr,"Unable to open FASTA file '%s', %d\n",fileName,GetLastError());
        return NULL;
    }

    LARGE_INTEGER fileSize;
    if (!GetFileSizeEx(hFile,&fileSize)) {
        fprintf(stderr,"Unable to get file size of FASTA file '%s', %d\n",fileName,GetLastError());
        CloseHandle(hFile);
        return NULL;
    }

    CloseHandle(hFile);
    hFile = NULL;
    if (0 == fileSize.QuadPart) {
        fprintf(stderr,"FASTA file '%s' appears to be empty\n",fileName);
        return NULL;
    }
    return fileSize.QuadPart;
#else
    // On Unix, use stat
    int fd = open(fileName, O_RDONLY);
    _ASSERT(fd != -1);
    struct stat sb;
    int r = fstat(fd, &sb);
    _ASSERT(r != -1);
    _uint64 fileSize = sb.st_size;
    close(fd);
    return fileSize;
#endif
}

}


    const Genome *
ReadFASTAGenome(const char *fileName)
{
    //
    // We need to know a bound on the size of the genome before we create the Genome object.
    // A bound is the number of bytes in the FASTA file, because we store at most one base per
    // byte.  Get the file size to use for this bound.
    //
    _uint64 fileSize = getFileSize(fileName);

    if (fileSize >> 32 != 0) {
        fprintf(stderr,"This tool only works with genomes with 2^32 bases or fewer.\n");
        return NULL;
    }

    Genome *genome = new Genome((unsigned) fileSize, (unsigned)fileSize);

    FILE *fastaFile = fopen(fileName, "r");
    if (fastaFile == NULL) {
        fprintf(stderr,"Unable to open FASTA file '%s' (even though we already got its size)\n",fileName);
        delete genome;
        return NULL;
    }

    const size_t lineBufferSize = 4096;
    char lineBuffer[lineBufferSize];

    while (NULL != fgets(lineBuffer,lineBufferSize,fastaFile)) {
        if (lineBuffer[0] == '>') {
            lineBuffer[strlen(lineBuffer) - 1] = '\0';   // Remove the trailing newline from fgets
            genome->startPiece(lineBuffer+1);
        } else {
            //
            // Convert it to upper case and truncate the newline before adding it to the genome.
            //

            char *newline = strchr(lineBuffer,'\n');
            if (NULL != newline) {
                *newline = 0;
            }

            //
            // But convert any 'N' to 'n'.  This is so we don't match the N from the genome with N
            // in reads (where we just do a straight text comparison.
            //
            size_t lineLen = strlen(lineBuffer);

			for (unsigned i = 0; i < lineLen; i++) {
              lineBuffer[i] = toupper(lineBuffer[i]);
            }

			for (unsigned i = 0; i < lineLen; i++) {
                if ('N' == lineBuffer[i]) {
                    lineBuffer[i] = 'n';
                }
            }
            genome->addData(lineBuffer);
        }
    }

    fclose(fastaFile);
    return genome;
}

//
// TODO: Reduce code duplication with the mutator.
// 
bool AppendFASTAGenome(const Genome *genome, FILE *fasta, const char *prefix="")
{
    int nPieces = genome->getNumPieces();
    const Genome::Piece *pieces = genome->getPieces();
    for (int i = 0; i < nPieces; ++i) {
        const Genome::Piece &piece = pieces[i];
        unsigned start = piece.beginningOffset;
        unsigned end = i + 1 < nPieces ? pieces[i + 1].beginningOffset : genome->getCountOfBases();
        unsigned size = end - start;
        const char *bases = genome->getSubstring(start, size);
        
        fprintf(fasta, ">%s%s\n", prefix, piece.name);
        fwrite(bases, 1, size, fasta);
        fputc('\n', fasta);
    }
    return !ferror(fasta);
}

    bool
AppendFASTADiploidGenome(const DiploidGenome *diploidGenome, FILE *fasta)
{
	return AppendFASTAGenome(diploidGenome->getGenome(true), fasta, diploidFASTASexPrefix(false)) &&
		AppendFASTAGenome(diploidGenome->getGenome(false), fasta, diploidFASTASexPrefix(true));			// getGenome takes isFemale, diploidFASTAPrefix isMale.
}

// 
// TODO: Factor out this boilerplate, to take a function object.
// 
bool WriteFASTADiploidGenome(const DiploidGenome *diploidGenome, const char *fileName)
{
    FILE *file = fopen(fileName, "wb");
    if (!file) {
        fprintf(stderr, "Can't write FASTA file '%s'\n", fileName);
        return false;
    }
    bool ok = AppendFASTADiploidGenome(diploidGenome, file);
    fclose(file);
    return ok;
}
