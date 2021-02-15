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
#include "Error.h"
#include "exit.h"
#include "Util.h"

using namespace std;

    bool
IsContigALT(
	const char		*contigName,
	GenomeDistance	 contigSize,
	const char* const*opt_in_alt_names,
	int				 opt_in_alt_names_count,
	const char* const*opt_out_alt_names,
	int				 opt_out_alt_names_count,
	GenomeDistance	 maxSizeForAutomaticALT,
    bool             autoALT)
{

	for (int i = 0; i < opt_out_alt_names_count; i++) {
		if (!_stricmp(opt_out_alt_names[i], contigName)) {
			return false;
		}
	} // opt out

	if (contigSize <= maxSizeForAutomaticALT) {
		return true;
	}

	for (int i = 0; i < opt_in_alt_names_count; i++) {
		if (!_stricmp(opt_in_alt_names[i], contigName)) {
			return true;
		} // match
	} // opt in

    if (autoALT && ((strlen(contigName) > 4 && !_stricmp(contigName + strlen(contigName) - 4, "_alt")) || 
                    (strlen(contigName) > 3 && (contigName[0] == 'H' || contigName[0] == 'h') && (contigName[1] == 'L' || contigName[1] == 'l') && (contigName[2] == 'A' || contigName[2] == 'a') && contigName[3] == '-')))  {
        return true;
    }

    return false;
} // MarkALTContigIfAppropriate


//
// There are several ways of specifying ALT contigs.  There is an opt-in list of ALTs, an opt-out list of regular chromosomes (these must be mutually
// exclusive), and a size cutoff below which is contig is an ALT.  The opt-in and opt-out lists supersede the size cutoff.
//

struct ContigLine {
    char* bases;
    ContigLine* next;

    ContigLine() {
        bases = NULL;
        next = NULL;
    }
};

struct RawContigData {
    ContigLine* lines;
    ContigLine* lastLine;
    GenomeDistance totalSize;
    char* name;
    RawContigData* next;
    int contigNumber;   // Where this contig is in the original FASTA file

    RawContigData() {
        lines = NULL;
        lastLine = NULL;
        totalSize = 0;
        name = NULL;
        next = NULL;
        contigNumber = -1;
    }

    void addLine(const char* bases) {
        ContigLine*line = new ContigLine();
        size_t size = strlen(bases) + 1;
        line->bases = new char[size];
        strncpy(line->bases, bases, size);
        if (lines == NULL) {
            lines = lastLine = line;
        } else {
            lastLine->next = line;
            lastLine = line;
        }
    } // addLine

    ~RawContigData() {
        delete[] name;

        while (lines != NULL) {
            ContigLine* toDelete = lines;
            lines = lines->next;
            delete[] toDelete->bases;
            delete toDelete;
        }
    }
}; // RawContigData

   void
AddRawContigToList(
    RawContigData*   currentContig,
    RawContigData**  altContigs,
    RawContigData**  regularContigs,
    const char* const* opt_in_alt_names,
    int				 opt_in_alt_names_count,
    const char* const* opt_out_alt_names,
    int				 opt_out_alt_names_count,
    GenomeDistance	 maxSizeForAutomaticALT,
    bool             autoALT)
{
    if (IsContigALT(currentContig->name, currentContig->totalSize, opt_in_alt_names, opt_in_alt_names_count, opt_out_alt_names,
        opt_out_alt_names_count, maxSizeForAutomaticALT, autoALT)) {
        currentContig->next = *altContigs;
        *altContigs = currentContig;
    }  else {
        currentContig->next = *regularContigs;
        *regularContigs = currentContig;
    }
}

   void
AddContigToGenome(
    RawContigData   *contig,
    Genome          *genome,
    char            *paddingBuffer)
{
    genome->addData(paddingBuffer);
    genome->startContig(contig->name);

    for (ContigLine* line = contig->lines; NULL != line; line = line->next) {
        genome->addData(line->bases);
    }
} // AddContigToGenome

   void
ReverseContigList(RawContigData** head) 
{
    if (*head == NULL) {
        return;
    }

    RawContigData* cur = *head;
    RawContigData* prev = NULL;

    for (;;) {
        RawContigData* next = cur->next;
        cur->next = prev;

        if (next == NULL) {
            *head = cur;
            return;
        }

        prev = cur;
        cur = next;
    } // for ever
} // ReverseContigList

    const Genome *
ReadFASTAGenome(
	const char		*fileName,
	const char		*pieceNameTerminatorCharacters,
	bool			 spaceIsAPieceNameTerminator,
	unsigned		 chromosomePaddingSize,
	const char* const*opt_in_alt_names,
	int				 opt_in_alt_names_count,
	const char* const*opt_out_alt_names,
	int				 opt_out_alt_names_count,
	GenomeDistance	 maxSizeForAutomaticALT,
    bool             autoALT)
{
    //
    // We need to know a bound on the size of the genome before we create the Genome object.
    // A bound is the number of bytes in the FASTA file, because we store at most one base per
    // byte. We count the bytes as we read them rather than getting the file size so that
    // we can deal with inputs that are redirected to pipes.
    //
    _int64 fileSize = 0;
    bool isValidGenomeCharacter[256];

    for (int i = 0; i < 256; i++) {
        isValidGenomeCharacter[i] = false;
    }

    isValidGenomeCharacter['A'] = isValidGenomeCharacter['T'] = isValidGenomeCharacter['C'] = isValidGenomeCharacter['G'] = isValidGenomeCharacter['N'] = true;
    isValidGenomeCharacter['a'] = isValidGenomeCharacter['t'] = isValidGenomeCharacter['c'] = isValidGenomeCharacter['g'] = isValidGenomeCharacter['n'] = true;

    FILE *fastaFile = fopen(fileName, "r");
    if (fastaFile == NULL) {
        WriteErrorMessage("Unable to open FASTA file '%s' (even though we already got its size)\n",fileName);
        return NULL;
    }

    int lineBufferSize = 0;
    char *lineBuffer;
 
    //
    // Read in the raw data and split it up by chromosome.  We're going to rearrange it to put the ALT chromosomes at the end so that it's
    // quick to test whether a GenomeLocation is ALT or not.
    //
    unsigned nContigs = 0;

    RawContigData* regularContigs = NULL;
    RawContigData* altContigs = NULL;
    RawContigData* currentContig = NULL;

    bool warningIssued = false;
    bool inAContig = false;

    int nextContigNumber = 0;

    while (NULL != reallocatingFgets(&lineBuffer, &lineBufferSize, fastaFile)) {
        fileSize += strlen(lineBuffer);
        if (lineBuffer[0] == '>') {
            nContigs++;

            if (inAContig) {
                AddRawContigToList(currentContig, &altContigs, &regularContigs, opt_in_alt_names, opt_in_alt_names_count, opt_out_alt_names, opt_out_alt_names_count, maxSizeForAutomaticALT, autoALT);
			}

            currentContig = new RawContigData();

            inAContig = true;

            //
            // Now supply the contig name.
            //
            if (NULL != pieceNameTerminatorCharacters) {
                for (int i = 0; i < strlen(pieceNameTerminatorCharacters); i++) {
                    char *terminator = strchr(lineBuffer+1, pieceNameTerminatorCharacters[i]);
                    if (NULL != terminator) {
                        *terminator = '\0';
                    }
                }
            }

            if (spaceIsAPieceNameTerminator) {
                char *terminator = strchr(lineBuffer, ' ');
                if (NULL != terminator) {
                    *terminator = '\0';
                }
                terminator = strchr(lineBuffer, '\t');
                if (NULL != terminator) {
                    *terminator = '\0';
                }
            }

            char *terminator = strchr(lineBuffer, '\n');
            if (NULL != terminator) {
                *terminator = '\0';
            }

            terminator = strchr(lineBuffer, '\r');
            if (NULL != terminator) {
                *terminator = '\0';
            }

            size_t nameLength = strlen(lineBuffer + 1) + 1;

            currentContig->name = new char[nameLength];
            strncpy(currentContig->name, lineBuffer + 1, nameLength);
            currentContig->contigNumber = nextContigNumber;
            nextContigNumber++;

        } else {
            if (!inAContig) {
                WriteErrorMessage("\nFASTA file doesn't begin with a contig name (i.e., the first line doesn't start with '>').\n");
                soft_exit(1);
            }

            //
            // Convert it to upper case and truncate the newline before adding it to the genome.
            //

            char *newline = strchr(lineBuffer, '\n');
            if (NULL != newline) {
                *newline = 0;
            }

            size_t lineLen = strlen(lineBuffer);

			for (unsigned i = 0; i < lineLen; i++) {
              lineBuffer[i] = toupper(lineBuffer[i]);
            }

			currentContig->totalSize += lineLen;

			for (unsigned i = 0; i < lineLen; i++) {
                if (!isValidGenomeCharacter[(unsigned char)lineBuffer[i]]) {
                    if (!warningIssued) {
                        WriteErrorMessage("\nFASTA file contained a character that's not a valid base (or N): '%c', full line '%s'; \nconverting to 'N'.  This may happen again, but there will be no more warnings.\n", lineBuffer[i], lineBuffer);
                        warningIssued = true;
                    }
                    lineBuffer[i] = 'N';
                }
            }

            currentContig->addLine(lineBuffer);
        }
    }

	if (!inAContig) {
		WriteErrorMessage("The FASTA file was empty.");
		return NULL;
	}

    AddRawContigToList(currentContig, &altContigs, &regularContigs, opt_in_alt_names, opt_in_alt_names_count, opt_out_alt_names, opt_out_alt_names_count, maxSizeForAutomaticALT, autoALT);

    //
    // AddRawContigToList reversed them.  Reverse them again so that they're in the same order as the FASTA, except that all ALTs follow all non-ALTs.
    // That's necessary to have the test for ALT be a simple comparison.  We fix that at sort time so they come out in the original order.
    //
    ReverseContigList(&altContigs);
    ReverseContigList(&regularContigs);

    Genome* genome = new Genome(fileSize + ((_int64)nContigs + 1) * (size_t)chromosomePaddingSize, fileSize + ((_int64)nContigs + 1) * (size_t)chromosomePaddingSize, chromosomePaddingSize, nContigs + 1, altContigs != NULL, true);

    char* paddingBuffer = new char[(GenomeDistance)chromosomePaddingSize + 1];
    for (unsigned i = 0; i < chromosomePaddingSize; i++) {
        paddingBuffer[i] = 'n';
    }
    paddingBuffer[chromosomePaddingSize] = '\0';

    while (regularContigs != NULL) {
        AddContigToGenome(regularContigs, genome, paddingBuffer);
        RawContigData* toDelete = regularContigs;
        regularContigs = regularContigs->next;
        delete toDelete;
    }

    while (altContigs != NULL) {
        AddContigToGenome(altContigs, genome, paddingBuffer);
        genome->markContigALT(altContigs->name);
        RawContigData* toDelete = altContigs;
        altContigs = altContigs->next;
        delete toDelete;
    }

    //
    // And finally add padding at the end of the genome.
    //
    genome->addData(paddingBuffer);
    genome->fillInContigLengths();
    genome->sortContigsByName();

    fclose(fastaFile);
    delete [] paddingBuffer;
    delete [] lineBuffer;
    return genome;
}

//
// TODO: Reduce code duplication with the mutator.
//
bool AppendFASTAGenome(const Genome *genome, FILE *fasta, const char *prefix="")
{
    int nContigs = genome->getNumContigs();
    const Genome::Contig *contigs = genome->getContigs();
    for (int i = 0; i < nContigs; ++i) {
        const Genome::Contig &contig = contigs[i];
        GenomeLocation start = contig.beginningLocation;
        GenomeLocation end = i + 1 < nContigs ? contigs[i + 1].beginningLocation : genome->getCountOfBases();
        GenomeDistance size = end - start;
        const char *bases = genome->getSubstring(start, size);

        fprintf(fasta, ">%s%s\n", prefix, contig.name);
        fwrite(bases, 1, size, fasta);
        fputc('\n', fasta);
    }
    return !ferror(fasta);
}
