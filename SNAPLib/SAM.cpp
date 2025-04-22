/*++

Module Name:

    SAM.cpp

Abstract:

    Sequence Alignment Map (SAM) file writer and reader.

Environment:

    User mode service.

    SamWriter and SamReader (and their subclasses) aren't thread safe.

--*/

#include "stdafx.h"
#include "BigAlloc.h"
#include "Compat.h"
#include "Read.h"
#include "SAM.h"
#include "Bam.h"
#include "Tables.h"
#include "RangeSplitter.h"
#include "ParallelTask.h"
#include "Util.h"
#include "ReadSupplierQueue.h"
#include "FileFormat.h"
#include "AlignerOptions.h"
#include "directions.h"
#include "exit.h"

using std::max;
using std::min;
using util::strnchr;

bool readIdsMatch(const char* id0, const char* id1, size_t len)
{
    const char* id0Base = id0;
    const char* id0End = id0 + len;
    while (true) {
        _uint64 x;
        x = *((_uint64*)id0) ^ *((_uint64*)id1);
        if (x) {
            unsigned long zeroes;
            CountTrailingZeroes(x, zeroes);
            zeroes >>= 3;
            return ((size_t)(id0 - id0Base) + (size_t)zeroes) >= len;
        }
        id0 += 8;
        if (id0 >= id0End) {
            return true;
        }
        id1 += 8;
    }
}

bool readIdsMatch(const char* id0, const char* id1, _int64* innerLoopCount)
{
    for (unsigned i = 0; ; i++) {
        char c0 = id0[i];
        char c1 = id1[i];

        (*innerLoopCount)++;

        if (c0 != c1) return false;
 
        // don't parse the read ID after the first space or slash, which can represent metadata (or which half of the mate pair the read is).
        if (c0 == 0 || c0 == ' ' || c0 == '/') return true;  
    }
    return true;
}

bool readIdsMatch(Read *read0, Read *read1)
{
    if (read0->getIdLength() != read1->getIdLength()) {
        return false;
    }

    for (unsigned i = 0; i < read0->getIdLength(); i++) {
        char c0 = read0->getId()[i];
        char c1 = read1->getId()[i];

        if (c0 != c1) return false;
 
        // don't parse the read ID after the first space or slash, which can represent metadata (or which half of the mate pair the read is).
        if (c0 == ' ' || c0 == '/') return true;  
    }

    return true;
}

    char *
strnchrs(char *str, char charToFind, char charToFind2, size_t maxLen) // Hokey version that looks for either of two chars
{
    for (size_t i = 0; i < maxLen; i++) {
        if (str[i] == charToFind || str[i] == charToFind2) {
            return str + i;
        }
        if (str[i] == 0) {
            return NULL;
        }
    }
    return NULL;
}

    char *
SAMReader::skipToBeyondNextFieldSeparator(char *str, const char *endOfBuffer, size_t *o_charsUntilFirstSeparator)
{
    if (NULL == str) return NULL;

    char *nextChar = str;
    while (nextChar < endOfBuffer && *nextChar != '\n' && *nextChar != '\t' && *nextChar != '\r' /* for Windows CRLF text */) {
        nextChar++;
    }

    if (NULL != o_charsUntilFirstSeparator) {
        *o_charsUntilFirstSeparator = nextChar - str;
    }

    if (nextChar >= endOfBuffer || *nextChar == '\n') {
        return NULL;
    }

    while (nextChar < endOfBuffer && ('\t' == *nextChar || '\r' == *nextChar)) {
        nextChar++;
    }

    if (nextChar >= endOfBuffer) {
        return NULL;
    }

    return nextChar;
}


    SAMReader *
SAMReader::create(
    DataSupplier* supplier,
    const char *fileName,
    int bufferCount,
    const ReaderContext& context,
    _int64 startingOffset, 
    _int64 amountOfFileToProcess)
{
    DataReader* data = supplier->getDataReader(bufferCount, maxLineLen, 0.0, 0);
    SAMReader *reader = new SAMReader(data, context);
    reader->init(fileName, startingOffset, amountOfFileToProcess);
    return reader;
}

    void
SAMReader::readHeader(const char *fileName)
{
    // todo: allow for larger headers
    _int64 headerSize = 512 * 1024; // 1M header initially (it's doubled before we use it)
	_int64 oldHeaderSize = 0;
 
	char* buffer;
	bool sawWholeHeader;
	do {
		headerSize *= 2;
		buffer = data->readHeader(&headerSize);
		if (oldHeaderSize >= headerSize) {
			//
			// No new data, we hit EOF
			//
			return;
		}
		oldHeaderSize = headerSize;

		if (!parseHeader(fileName, buffer, buffer + headerSize, context.genome, &headerSize, &context.headerMatchesIndex, &sawWholeHeader, NULL, NULL, &numRGLines, &rgLines, &rgLineOffsets)) {
			WriteErrorMessage("SAMReader: failed to parse header on '%s'\n", fileName);
			soft_exit(1);
		}
	} while (!sawWholeHeader);
    _ASSERT(context.header == NULL);
    char* p = new char[headerSize + 1];
    memcpy(p, buffer, headerSize);
    p[headerSize] = 0;
    context.header = p;
    context.headerBytes = context.headerLength = headerSize;
    context.numRGLines = numRGLines;
    context.rgLines = rgLines;
    context.rgLineOffsets = rgLineOffsets;
}

SAMReader::SAMReader(
    DataReader* i_data,
    const ReaderContext& i_context)
    : ReadReader(i_context), data(i_data), headerSize(-1), clipping(i_context.clipping), numRGLines(0), rgLines(NULL), rgLineOffsets(NULL)
{
}


//
// Implement the ReadReader form of getNextRead, which doesn't include the
// alignment results by simply throwing them away.
//
    bool
SAMReader::getNextRead(Read *readToUpdate)
{
    return getNextRead(readToUpdate, NULL, NULL, NULL, NULL, NULL, NULL);
}

    bool
SAMReader::parseHeader(
    const char *fileName, 
    char *firstLine, 
    char *endOfBuffer, 
    const Genome *genome, 
    _int64 *o_headerSize,
    bool *o_headerMatchesIndex,
	bool *o_sawWholeHeader,
    int *o_n_ref,
    GenomeLocation **o_ref_locations,
    int *o_n_rg,
    char **o_rgLines,
    size_t **o_rgLineOffsets)
{
    _ASSERT((NULL == o_n_ref) == (NULL == o_ref_locations));    // both or neither are NULL, not one of each

    char *nextLineToProcess = firstLine;
    *o_headerMatchesIndex = true;
    int numSQLines = 0;
    int numRGLines = 0;
    int n_ref_slots = 4096;
    size_t n_rg_slots = 128;
    size_t rg_slot_size = 1024;
    size_t rg_total_size = n_rg_slots * rg_slot_size;
    GenomeLocation *ref_locations = NULL;
    size_t lineBufferSize = 1024;
    char* lineBuffer = new char[lineBufferSize];
    char* rgLines = NULL;
    size_t* rgLineOffsets = NULL;

    if (NULL != o_ref_locations) {
        ref_locations = (GenomeLocation *)BigAlloc(sizeof(GenomeLocation)* n_ref_slots);
    }

    if (NULL != o_rgLines) {
        rgLines = new char[rg_total_size];
        rgLineOffsets = new size_t[n_rg_slots];
        rgLineOffsets[0] = 0;
    }

    while (NULL != nextLineToProcess && nextLineToProcess < endOfBuffer && '@' == *nextLineToProcess) {
		//
		// Make sure we have the complete line.
		//
		bool foundCompleteLine = false;
        char* endOfLine;

		for (char *c = nextLineToProcess; c < endOfBuffer; c++) {
			if (*c == '\n') {
				foundCompleteLine = true;
                if ((size_t)(c - nextLineToProcess) + 1 >= lineBufferSize) {
                    delete [] lineBuffer;
                    lineBufferSize = max(lineBufferSize * 2, (size_t)(c - nextLineToProcess) + 2);
                    lineBuffer = new char[lineBufferSize];
                }
                memcpy(lineBuffer, nextLineToProcess, c - nextLineToProcess);
                lineBuffer[c - nextLineToProcess] = '\0';
                endOfLine = &lineBuffer[c - nextLineToProcess + 1];
                break;
			}
		}
		if (!foundCompleteLine) {
			*o_sawWholeHeader = false;
            delete[] lineBuffer;
			return true;	// Parsed OK, but incomplete
		}

        if (!strncmp("@SQ",lineBuffer,3)) {
            //
            // These lines represent sequences in the reference genome, what are
            // called "contigs" in the Genome class.  (Roughly, chromosomes or major
            // variants like some versions of the MHC genes on chr6; or more
            // particularly the things that come in different FASTA files from the
            // reference assembly).
            //
            // Verify that they actually match what's in our reference genome.
            //
            if (lineBuffer + 3 >= endOfLine || (' ' != lineBuffer[3] && '\t' != lineBuffer[3])) {
                WriteErrorMessage("Malformed SAM file '%s' has @SQ without a following space or tab.\n", fileName);
                delete[] lineBuffer;
                return false;
            }

            char *snStart = lineBuffer + 4;
            while (snStart < endOfLine && strncmp(snStart,"SN:",__min(3, endOfLine -snStart)) && *snStart != '\n' && *snStart != 0) {
                snStart++;
            }

            if (snStart >= endOfLine || *snStart == '\n' || *snStart == 0) {
                WriteErrorMessage("Malformed @SQ line doesn't have 'SN:' in file '%s'\n",fileName);
                delete[] lineBuffer;
                return false;
            }

            size_t contigNameBufferSize = 512;
            char *contigName = new char[contigNameBufferSize];
            for (unsigned i = 0; snStart+3+i < endOfLine; i++) {
                if (i >= contigNameBufferSize) {
                    //
                    // Get more buffer.
                    //
                    size_t newSize = contigNameBufferSize * 2;
                    char *newBuffer = new char[newSize];
                    memcpy(newBuffer, contigName, contigNameBufferSize);
                    delete[] contigName;
                    contigName = newBuffer;
                    contigNameBufferSize = newSize;
                }
                if (snStart[3+i] == ' ' || snStart[3+i] == '\t' || snStart[3+i] == '\n' || snStart[3+i] == 0) {
                    contigName[i] = '\0';
                } else {
                    contigName[i] = snStart[3+i];
                }
            }
            contigName[contigNameBufferSize - 1] = '\0';

            GenomeLocation contigBase = InvalidGenomeLocation;
            if (NULL != genome) {
                if (!genome->getLocationOfContig(contigName, &contigBase)) {
                    //
                    // Try stripping off the leading "chr" if it has one.
                    //
                    if (strlen(contigName) < 3 || contigName[0] != 'c' || contigName[1] != 'h' || contigName[2] != 'r' || !genome->getLocationOfContig(contigName + 3, &contigBase)) {
                        //
                        // Change chrM into MT
                        //
                        if (strlen(contigName) != 4 || contigName[0] != 'c' || contigName[1] != 'h' || contigName[2] != 'r' || contigName[3] != 'M' || !genome->getLocationOfContig("MT", &contigBase)) {
                            //
                            // Try prefixing chr to short names.
                            //
                            const size_t maxShortNameSize = 2;
                            char prefixedName[maxShortNameSize + 4];
                            if (strlen(contigName) <= maxShortNameSize) {
                                sprintf(prefixedName, "chr%s", contigName);
                                genome->getLocationOfContig(prefixedName, &contigBase);
                            }
                        }
                    }
                }
            }

            //
            // If we have not seen the contig before we cannot copy over header.
            // Also if we are running with ALT contigs, don't copy over header even if the contig name matches, since the contig's ALT status may have changed from the previous run
            //
            if (InvalidGenomeLocation == contigBase || (genome->isGenomeLocationALT(contigBase))) {
                *o_headerMatchesIndex = false;
            }

            if (NULL != o_ref_locations) {
                if (numSQLines >= n_ref_slots) {
                    GenomeLocation *new_ref_locations = (GenomeLocation *)BigAlloc(sizeof(GenomeLocation)* n_ref_slots * 2);
                    memcpy(new_ref_locations, ref_locations, sizeof(GenomeLocation)* n_ref_slots);
                    BigDealloc(ref_locations);
                    ref_locations = new_ref_locations;
                    n_ref_slots *= 2;
                }
                ref_locations[numSQLines] = contigBase;
            }

            numSQLines++;
            delete[] contigName;
        } else if (!strncmp("@HD", nextLineToProcess, 3) || !strncmp("@RG", nextLineToProcess, 3) || !strncmp("@PG", nextLineToProcess, 3) ||
            !strncmp("@CO",nextLineToProcess,3)) {

            if (!strncmp("@RG", nextLineToProcess, 3)) {
                if (nextLineToProcess + 3 >= endOfBuffer || (' ' != nextLineToProcess[3] && '\t' != nextLineToProcess[3])) {
                    WriteErrorMessage("Malformed SAM file '%s' has @RG without a following space or tab.\n", fileName);
                    return false;
                }

                char* rgStart = nextLineToProcess + 4;
                char* rgSlot = new char[rg_slot_size];
                int bytesUsed = 0;
                bool foundID = false;
                for (int i = 0; rgStart + i < endOfBuffer; i++) {
                    if (bytesUsed >= rg_slot_size) {
                        //
                        // Get more buffer.
                        //
                        size_t newSize = (size_t)rg_slot_size * 2;
                        char* newBuffer = new char[newSize];
                        memcpy(newBuffer, rgSlot, rg_slot_size);
                        delete[] rgSlot;
                        rgSlot = newBuffer;
                        rg_slot_size = newSize;
                    }
                    if (rgStart[i] == 'I') {
                        if (rgStart + i + 2 >= endOfBuffer) {
                            delete[] rgSlot;
                            return false;
                        }
                        if (rgStart[i + 1] == 'D' && rgStart[i + 2] == ':') {
                            for (int j = i + 3; rgStart + j < endOfBuffer; j++) {
                                if (bytesUsed >= rg_slot_size) {
                                    //
                                    // Get more buffer.
                                    //
                                    size_t newSize = (size_t)rg_slot_size * 2;
                                    char* newBuffer = new char[newSize];
                                    memcpy(newBuffer, rgSlot, rg_slot_size);
                                    delete[] rgSlot;
                                    rgSlot = newBuffer;
                                    rg_slot_size = newSize;
                                }
                                if (rgStart[j] == '\t' || rgStart[j] == 0 || rgStart[j] == '\n' || rgStart[j] == ' ') {
                                    foundID = true;
                                    rgSlot[bytesUsed++] = '\t';
                                    break;
                                }
                                else {
                                    rgSlot[bytesUsed++] = rgStart[j];
                                }
                            }
                        }
                    }
                    if (foundID) break;
                }

                if (foundID) {
                    bool foundLB = false;
                    for (int i = 0; rgStart + i < endOfBuffer; i++) {
                        if (bytesUsed >= rg_slot_size) {
                            //
                            // Get more buffer.
                            //
                            size_t newSize = (size_t)rg_slot_size * 2;
                            char* newBuffer = new char[newSize];
                            memcpy(newBuffer, rgSlot, rg_slot_size);
                            delete[] rgSlot;
                            rgSlot = newBuffer;
                            rg_slot_size = newSize;
                        }

                        if (rgStart[i] == 'L') {
                            if (rgStart + i + 2 >= endOfBuffer) {
                                delete[] rgSlot;
                                return false;
                            }

                            if (rgStart[i + 1] == 'B' && rgStart[i + 2] == ':') {
                                for (int j = i + 3; rgStart + j < endOfBuffer; j++) {
                                    if (bytesUsed >= rg_slot_size) {
                                        //
                                        // Get more buffer.
                                        //
                                        size_t newSize = (size_t)rg_slot_size * 2;
                                        char* newBuffer = new char[newSize];
                                        memcpy(newBuffer, rgSlot, rg_slot_size);
                                        delete[] rgSlot;
                                        rgSlot = newBuffer;
                                        rg_slot_size = newSize;
                                    }

                                    if (rgStart[j] == '\t' || rgStart[j] == 0 || rgStart[j] == '\n' || rgStart[j] == ' ') {
                                        rgSlot[bytesUsed++] = '\t';
                                        foundLB = true;
                                        break;
                                    } else {
                                        rgSlot[bytesUsed++] = rgStart[j];
                                    }
                                } // for j
                            } //if B:
                        } // L

                        if (foundLB) {
                            break;
                        }
                    } // for i
                } // if foundID

                rgSlot[bytesUsed++] = '\0';

                numRGLines++;

                if (NULL != o_rgLines) {
                    if ((_int64)numRGLines + 1 >= (_int64)n_rg_slots) {
                        char* newRGLines = new char[n_rg_slots * rg_slot_size * 2];
                        memcpy(newRGLines, rgLines, sizeof(char) * n_rg_slots * rg_slot_size);
                        size_t* newRGLineOffsets = new size_t[n_rg_slots * 2];
                        memcpy(newRGLineOffsets, rgLineOffsets, sizeof(size_t) * n_rg_slots);
                        delete[] rgLines;
                        delete[] rgLineOffsets;
                        rgLines = newRGLines;
                        rgLineOffsets = newRGLineOffsets;
                        n_rg_slots *= 2;
                    }
                    rgLineOffsets[numRGLines] = rgLineOffsets[numRGLines - 1] + strlen(rgSlot);
                    strcpy(rgLines + rgLineOffsets[numRGLines - 1], rgSlot);
                }

                delete[] rgSlot;
            } // RG
        } else {
            WriteErrorMessage("Unrecognized header line in SAM file.\n");
            delete[] lineBuffer;
            return false;
        }

		char * p = strnchr(nextLineToProcess,'\n',endOfBuffer-nextLineToProcess);
		if (p == NULL) {
            // no newline, look for null to truncate buffer
            p = (char*) memchr(nextLineToProcess, 0, endOfBuffer - nextLineToProcess);
            nextLineToProcess = p != NULL ? p + 1 : endOfBuffer;
            break;
		}
        nextLineToProcess = p + 1;
    }

    *o_headerMatchesIndex &= genome != NULL && numSQLines == genome->getNumContigs();
	*o_headerSize = nextLineToProcess - firstLine;
	if (NULL != o_sawWholeHeader) {
		*o_sawWholeHeader = nextLineToProcess < endOfBuffer;
	}

    if (NULL != o_ref_locations) {
        *o_n_ref = numSQLines;
        *o_ref_locations = ref_locations;
    }

    delete[] lineBuffer;

    if (NULL != o_rgLines) {
        *o_n_rg = numRGLines;
        *o_rgLines = rgLines;
        *o_rgLineOffsets = rgLineOffsets;
    }

    return true;
}

    bool
SAMReader::parseLine(char *line, char *endOfBuffer, char *result[], size_t *linelength, size_t fieldLengths[])
{
    *linelength = 0;

    char* prev = line;
    char *next = line;
    char *endOfLine = strnchr(line,'\n',endOfBuffer-line);
    if (NULL == endOfLine) {
        return false;
    }

    //
    // Skip over any leading spaces and tabs
    //
    while (next < endOfLine && (*next == ' ' || *next == '\t')) {
        next++;
    }

    for (unsigned i = 0; i < nSAMFields; i++) {
        if (NULL == next || next >= endOfLine) {
            if (i == OPT) {
                // no optional fields
                result[OPT] = NULL;
                break;
            } else {
                //
                // Too few fields.
                //
                return false;
            }
        }

        result[i] = next;
        if (i == OPT) {
            // OPT field is actually all fields until end of line
            fieldLengths[OPT] = endOfLine - next;
            break;
        }

        next = skipToBeyondNextFieldSeparator(next,endOfLine,&fieldLengths[i]);
    }

    *linelength =  endOfLine - line + 1;    // +1 skips over the \n
    return true;
}

    void
SAMReader::getReadFromLine(
    const Genome        *genome,
    char                *line, 
    char                *endOfBuffer, 
    Read                *read, 
    AlignmentResult     *alignmentResult,
    GenomeLocation      *out_genomeLocation, 
    Direction           *direction,
    unsigned            *mapQ,
    size_t              *lineLength,
    unsigned *           flag,
    const char **        cigar,
    ReadClippingType     clipping,
    char                *rgLines,
    int                  numRGLines,
    size_t              *rgLineOffsets
    )
{
    char *field[nSAMFields];
    size_t fieldLength[nSAMFields];

    if (!parseLine(line, endOfBuffer, field, lineLength, fieldLength)) {
        WriteErrorMessage( "Failed to parse SAM line:\n%.*s\n", lineLength, line);
        soft_exit(1);
    }

    //
    // We have to copy the contig name (RNAME) into its own buffer because the code in Genome expects
    // it to be a null-terminated string, while all we've got is one that's space delimited.
    //
    const size_t contigNameBufferSize = 512;
    char contigNameBuffer[contigNameBufferSize];
    char *contigName = contigNameBuffer;
    size_t neededSize;
    GenomeLocation locationOfContig;
    if (0 != (neededSize = parseContigName(genome, contigName, contigNameBufferSize, &locationOfContig, NULL, field, fieldLength))) {
        contigName = new char[neededSize];
        if (0 != parseContigName(genome, contigName, neededSize, &locationOfContig, NULL, field, fieldLength)) {
            WriteErrorMessage("SAMReader::getReadFromLine: reallocated contigName was still too small\n");
            soft_exit(1);
        }
    }

    GenomeLocation genomeLocation;
    
    if (InvalidGenomeLocation != locationOfContig) {
        genomeLocation = parseLocation(locationOfContig, field, fieldLength);
    } else {
        genomeLocation = InvalidGenomeLocation;
    }

    if (NULL != out_genomeLocation) {
        *out_genomeLocation = genomeLocation;
    }

    if (fieldLength[SEQ] != fieldLength[QUAL]) {
        WriteErrorMessage("SAMReader: QUAL string unequal in length to SEQ string.\n");
        soft_exit(1);
    }

    unsigned _flag;
    const size_t flagBufferSize = 20;   // More than enough
    char flagBuffer[flagBufferSize];
    if (fieldLength[FLAG] >= flagBufferSize) {
        WriteErrorMessage("SAMReader: flag field is too long.\n");
        soft_exit(1);
    }
    memcpy(flagBuffer,field[FLAG],fieldLength[FLAG]);
    flagBuffer[fieldLength[FLAG]] = '\0';
    if (1 != sscanf(flagBuffer,"%d",&_flag)) {
        WriteErrorMessage("SAMReader: couldn't parse FLAG field.\n");
        soft_exit(1);
    }

    if (NULL != read) {
        //
        // Clip reads where the quality strings end in '#'
        //

        unsigned originalFrontClipping, originalBackClipping, originalFrontHardClipping, originalBackHardClipping;
        Read::computeClippingFromCigar(field[CIGAR], &originalFrontClipping, &originalBackClipping, &originalFrontHardClipping, &originalBackHardClipping);

        unsigned pnext = atoi(field[PNEXT]);    // Relies on atoi() returning 0 for non-numeric fields (i.e., *)

        read->init(field[QNAME],(unsigned)fieldLength[QNAME],field[SEQ],field[QUAL],(unsigned)fieldLength[SEQ], genomeLocation, atoi(field[MAPQ]), _flag, 
            originalFrontClipping, originalBackClipping, originalFrontHardClipping, originalBackHardClipping, field[RNEXT], (unsigned)fieldLength[RNEXT], pnext);

        //
        // If this read is RC in the SAM file, we need to reverse it here, since Reads are always the sense that they were as they came
        // out of the base caller.
        //

        if (_flag & SAM_REVERSE_COMPLEMENT) {
            read->becomeRC();
        }
        read->clip(clipping);

        if (field[OPT] != NULL) {
            unsigned n = (unsigned) fieldLength[OPT];
            while (n > 0 && (field[OPT][n-1] == '\n' || field[OPT][n-1] == '\r')) {
                n--;
            }

            read->setAuxiliaryData(field[OPT], n);
            char* rgFromAux = NULL;
            size_t rgFromAuxLen = 0;
            for (char* p = field[OPT]; p != NULL && p < field[OPT] + fieldLength[OPT]; p = SAMReader::skipToBeyondNextFieldSeparator(p, field[OPT] + fieldLength[OPT])) {
                if (strncmp(p, "RG:Z:", 5) == 0) {
                    rgFromAux = p + 5;
                    p = SAMReader::skipToBeyondNextFieldSeparator(p, field[OPT] + fieldLength[OPT]);
                    rgFromAuxLen = (p != NULL) ? p - 1 - rgFromAux : field[OPT] + fieldLength[OPT] - rgFromAux;
                    read->setReadGroup(READ_GROUP_FROM_AUX);
                    break;
                }
            }

            // LB
            if (rgLines != NULL) {
                // get library for read group
                for (int i = 0; i < numRGLines; i++) {
                    char* rgStart = rgLines + rgLineOffsets[i];
                    char* rgEnd = strchr(rgStart, '\t');
                    size_t rgLen = rgEnd - rgStart;
                    if (rgFromAuxLen != rgLen) continue;
                    if (!strncmp(rgFromAux, rgStart, rgLen)) {
                        char* lbStart = rgEnd + 1;
                        char* lbEnd = strchr(lbStart, '\t');
                        if (lbEnd != NULL) {
                            read->setLibrary(lbStart);
                            read->setLibraryLength((int)(lbEnd - lbStart));
                        }
                    }
                }
            } // rgLines != NULL
        } // field[OPT] != NULL
    } // NULL != read

    if (NULL != alignmentResult) {
        if (_flag & SAM_UNMAPPED) {
            *alignmentResult = NotFound;
        } else {
            if ('*' == contigName[0]) {
                WriteErrorMessage("SAMReader: mapped read didn't have RNAME filled in.\n");
                soft_exit(1);
            }
            *alignmentResult = SingleHit;   // NB: This isn't quite right, we should look at MAPQ.
        }
    }

    if (NULL != direction) {
        *direction = (_flag & SAM_REVERSE_COMPLEMENT) ? RC : FORWARD;
    }

    if (NULL != mapQ) {
        *mapQ = atoi(field[MAPQ]);
        if (*mapQ > 255) {
            WriteErrorMessage("SAMReader: MAPQ field has bogus value\n");
            soft_exit(1);
        }
    }

    if (NULL != flag) {
        *flag = _flag;
    }

    if (NULL != cigar) {
        *cigar = field[CIGAR];
    }

    if (contigName != contigNameBuffer) {
        delete[] contigName;
    }
}

    size_t
SAMReader::parseContigName(
    const Genome* genome,
    char* contigName,
    size_t contigNameBufferSize,
    GenomeLocation* o_locationOfContig,
    InternalContigNum* o_indexOfContig,
    char* field[],
    size_t fieldLength[],
	unsigned rfield)
{
    if (fieldLength[rfield] >= contigNameBufferSize) {  // >= because we need a byte for the \0
        return fieldLength[rfield] + 1; // +1 for trailing null
    }
    
    memcpy(contigName,field[rfield],fieldLength[rfield]);
    contigName[fieldLength[rfield]] = '\0';

    *o_locationOfContig = InvalidGenomeLocation;
    if ('*' != contigName[0] && genome != NULL && !genome->getLocationOfContig(contigName, o_locationOfContig, o_indexOfContig)) {
        //WriteErrorMessage("Unable to find contig '%s' in genome.  SAM file malformed.\n",contigName);
        //soft_exit(1);
    }
    return 0;
}

    GenomeLocation
SAMReader::parseLocation(
    GenomeLocation locationOfContig,
    char* field[],
    size_t fieldLength[],
	unsigned rfield,
	unsigned posfield,
    int *o_pos)
{
    unsigned oneBasedOffsetWithinContig = 0;
    if ('*' != field[rfield][0] && '*' != field[posfield][0]) {
        //
        // We can't call sscanf directly into the mapped file, becuase it reads to the end of the
        // string even when it's satisfied all of its fields.  Since this can be gigabytes, it's not
        // really good for perf.  Instead, copy the POS field into a local buffer and null terminate it.
        //

        const unsigned posBufferSize = 20;
        char posBuffer[posBufferSize];
        if (fieldLength[posfield] >= posBufferSize) {
            WriteErrorMessage("SAMReader: POS field too long.\n");
            soft_exit(1);
        }

        memcpy(posBuffer,field[posfield],fieldLength[posfield]);
        posBuffer[fieldLength[posfield]] = '\0';
        if (0 == sscanf(posBuffer,"%d",&oneBasedOffsetWithinContig)) {
            WriteErrorMessage("SAMReader: Unable to parse position when it was expected.\n");
            soft_exit(1);
        }

        if (0 == oneBasedOffsetWithinContig) {
            WriteErrorMessage("SAMReader: Position parsed as 0 when it was expected.\n");
            soft_exit(1);
        }

        if (NULL != o_pos) {
            *o_pos = oneBasedOffsetWithinContig;
        }

        return locationOfContig + oneBasedOffsetWithinContig - 1; // -1 is because our offset is 0 based, while SAM is 1 based.
    } else {
        if (NULL != o_pos) {
            *o_pos = 0;
        }

        return InvalidGenomeLocation;
    }
}
    
    void
SAMReader::init(
    const char *fileName,
    _int64 startingOffset,
    _int64 amountOfFileToProcess)
{
    if (! data->init(fileName)) {
        WriteErrorMessage( "Unable to read file %s\n", fileName);
        soft_exit(1);
    }

    if (0 == startingOffset) {
        readHeader(fileName);
    }

    headerSize = context.headerBytes;
    reinit(max(startingOffset, (_int64) context.headerBytes),
        amountOfFileToProcess == 0 || startingOffset >= (_int64) context.headerBytes ? amountOfFileToProcess
            : amountOfFileToProcess - (context.headerBytes - startingOffset));
}

    void
SAMReader::reinit(_int64 startingOffset, _int64 amountOfFileToProcess)
{
    _ASSERT(-1 != headerSize && startingOffset >= headerSize);  // Must call init() before reinit()
    //
    // There's no way to tell if we start at the very beginning of a read, we need to see the previous newline.
    // So, read one byte before our assigned read in case that was the terminating newline of the previous read.
    //
    if (startingOffset > headerSize) {
        startingOffset--;
        amountOfFileToProcess++;
    }
    data->reinit(startingOffset, amountOfFileToProcess);
    char* buffer;
    _int64 validBytes;
    if (!data->getData(&buffer, &validBytes)) {
        return;
    }
    if (startingOffset != headerSize) {
        char *firstNewline = strnchr(buffer,'\n',validBytes);
        if (NULL == firstNewline) {
            return;
        }

        data->advance((unsigned)(firstNewline - buffer + 1)); // +1 skips over the newline.
    }
}

    bool
SAMReader::getNextRead(
    Read *read,
    AlignmentResult *alignmentResult,
    GenomeLocation *genomeLocation,
    Direction *direction,
    unsigned *mapQ, 
    unsigned *flag,
    bool ignoreEndOfRange,
    const char **cigar)
{
    unsigned local_flag;
    if (NULL == flag) {
        flag = &local_flag;
    }

    do {
        char* buffer;
        _int64 bytes;
        if (! data->getData(&buffer, &bytes)) {
            data->nextBatch();
            if (! data->getData(&buffer, &bytes)) {
                return false;
            }
        }

        char *newLine = strnchr(buffer, '\n', bytes);
        if (NULL == newLine) {
            //
            // There is no newline, so the line crosses the end of the buffer.
            // This should never happen since underlying reader manages overflow between chunks.
            //
            WriteErrorMessage("SAM file has too long a line, or doesn't end with a newline!  Failing.  fileOffset = %lld\n", data->getFileOffset());
            soft_exit(1);
        }

        size_t lineLength;
        read->setReadGroup(context.defaultReadGroup);
        getReadFromLine(context.genome, buffer,buffer + bytes, read, alignmentResult, genomeLocation, direction, mapQ, &lineLength, flag, cigar, clipping,
            context.rgLines, context.numRGLines, context.rgLineOffsets);
        read->setBatch(data->getBatch());
        data->advance((newLine + 1) - buffer);

    } while ((context.ignoreSecondaryAlignments && ((*flag) & SAM_SECONDARY)) ||
             (context.ignoreSupplementaryAlignments && ((*flag) & SAM_SUPPLEMENTARY)));

    return true;
}

    ReadSupplierGenerator *
SAMReader::createReadSupplierGenerator(
    const char *fileName,
    int numThreads,
    const ReaderContext& context)
{
    //
    // single-ended SAM files always can be read with the range splitter, unless reading from stdin, which needs a queue
    //
    if (!strcmp(fileName, "-")) {
        //
        // Stdin must run from a queue, not range splitter.
        //
        ReadReader* reader;
        //
        // Because we can only have one stdin reader, we need to use a queue if we're reading from stdin
        //
        reader = SAMReader::create(DataSupplier::Stdio, "-", ReadSupplierQueue::BufferCount(numThreads), context, 0, 0);
   
        if (reader == NULL) {
            return NULL;
        }

        ReadSupplierQueue *queue = new ReadSupplierQueue(reader);
        queue->startReaders();
        return queue;
    } else {
        RangeSplitter *splitter = new RangeSplitter(QueryFileSize(fileName), numThreads, 100);
        return new RangeSplittingReadSupplierGenerator(fileName, true, numThreads, context);
    }
}
    
    PairedReadReader*
SAMReader::createPairedReader(
    const DataSupplier* supplier,
    const char *fileName,
    int bufferCount,
    _int64 startingOffset,
    _int64 amountOfFileToProcess, 
    bool quicklyDropUnpairedReads,
    const ReaderContext& context)
{
    DataSupplier *data;
    if (!strcmp("-", fileName)) {
        data = DataSupplier::Stdio;
    } else {
        data = DataSupplier::Default;
    }

    SAMReader* reader = SAMReader::create(data, fileName, bufferCount + PairedReadReader::MatchBuffers, context, 0, 0);
    if (reader == NULL) {
        return NULL;
    }
    return PairedReadReader::PairMatcher(reader, quicklyDropUnpairedReads);
}


    PairedReadSupplierGenerator *
SAMReader::createPairedReadSupplierGenerator(
    const char *fileName,
    int numThreads,
    bool quicklyDropUnpairedReads, 
    const ReaderContext& context)
{
    //
    // need to use a queue so that pairs can be matched
    //

    PairedReadReader* paired = SAMReader::createPairedReader(DataSupplier::Default, fileName,
        ReadSupplierQueue::BufferCount(numThreads), 0, 0, quicklyDropUnpairedReads, context);
    if (paired == NULL) {
        WriteErrorMessage( "Cannot create reader on %s\n", fileName);
        soft_exit(1);
    }
    ReadSupplierQueue* queue = new ReadSupplierQueue(paired);
    queue->startReaders();
    return queue;
}


const FileFormat* FileFormat::SAM[] = { new SAMFormat(false), new SAMFormat(true) };

    void
SAMFormat::getSortInfo(
    const Genome* genome,
    char* buffer,
    _int64 bytes,
	GenomeLocation* o_location,
	GenomeDistance* o_readBytes,
	OriginalContigNum* o_refID,
	int* o_pos) const
{
    char* fields[SAMReader::nSAMFields];
    size_t lengths[SAMReader::nSAMFields];
    size_t lineLength;
    SAMReader::parseLine(buffer, buffer + bytes, fields, &lineLength, lengths);
    _ASSERT(lineLength < UINT32_MAX);
	if (o_readBytes != NULL) {
		*o_readBytes = (unsigned) lineLength;
	}

    GenomeLocation locationBuffer;
    if (o_location == NULL) {
        o_location = &locationBuffer;   // We need to call parseLocation to figure out pos, so if we're not returning the value just store it on the stack and discard
    }

    //
    // Fill these in to indicate whether we're sorting by the aligned location or by next (which we do for unaligned reads so they land by their mate pair).
    //
    int posField;
    int rField;

    if (lengths[SAMReader::RNAME] == 0 || fields[SAMReader::RNAME][0] == '*') {
        if (lengths[SAMReader::RNEXT] == 0 || fields[SAMReader::RNEXT][0] == '*' || fields[SAMReader::RNEXT][0] == '=') { // The check for = means it's unaligned because we only get here if RNAME is '*'
            if (o_location != NULL) {
                *o_location = InvalidGenomeLocation;
            }

            if (o_refID != NULL) {
                *o_refID = OriginalContigNum(INT32_MAX);    // So that it sorts to the end
            }

            if (o_pos != NULL) {
                *o_pos = 0;
            }

            return;
        } else {
            posField = SAMReader::PNEXT;
            rField = SAMReader::RNEXT;
        }
    } else {
        posField = SAMReader::POS;
        rField = SAMReader::RNAME;
    }
  
	const size_t contigNameBufferSize = 512;        // We do a static buffer with reallocation so that in the usual case there is no dynamic memory allocation.  If you have enormous contig names, you'll just run a little slower
    char contigNameBuffer[contigNameBufferSize];
    char *contigName = contigNameBuffer;
	GenomeLocation locationOfContig;
    size_t neededSize;
    InternalContigNum internalContigNum;
    if (0 != (neededSize = SAMReader::parseContigName(genome, contigName, contigNameBufferSize, &locationOfContig, &internalContigNum, fields, lengths, rField))) {
        //
        // Need a bigger buffer.
        //
        contigName = new char[neededSize];
        if (0 != SAMReader::parseContigName(genome, contigName, neededSize, &locationOfContig, &internalContigNum, fields, lengths, rField)) {
            WriteErrorMessage("SAMFormat::getSortInfo: reallocated buffer size is still too small.\n"); // This really shouldn't happen
            soft_exit(1);
        }
    }

    *o_location = SAMReader::parseLocation(locationOfContig, fields, lengths, rField, posField, o_pos);

    if (o_refID != NULL) {
        *o_refID = genome->getContigByInternalNumber(internalContigNum)->originalContigNumber;
    }

    if (contigName != contigNameBuffer) {
        delete[] contigName;
    }
}

// which @RG line fields to put in aux data of every read
const char* FileFormat::RGLineToAux = "IDLBPLPUSM";

    void
FileFormat::setupReaderContext(
    AlignerOptions* options,
    ReaderContext* readerContext,
    bool bam)
{
    if (options->rgLineContents == NULL || *options->rgLineContents == '\0') {
        readerContext->defaultReadGroupAux = "";
        readerContext->defaultReadGroupAuxLen = 0;
        return;
    }
    char* buffer = new char[strlen(options->rgLineContents) * 3]; // can't expend > 2x
    const char* from = options->rgLineContents;
    char* to = buffer;

    // skip @RG
    _ASSERT(strncmp(from, "@RG", 3) == 0);
    while (*from && *from != '\t') {
        from++;
    }

    while (*from) {
        if (!(from[0] == '\t' && from[1] && from[1] != '\t' && from[2] && from[2] != '\t' && from[3] == ':')) {
            WriteErrorMessage("Invalid @RG line: %s\n", options->rgLineContents);
            soft_exit(1);
        }

        bool keep = false;
        bool isID = false;
        for (const char* a = RGLineToAux; *a; a += 2) {
            if (from[1] == a[0] && from[2] == a[1]) {
                keep = true;
                isID = from[1] == 'I' && from[2] == 'D';
                break;
            }
        }

        if (keep) {
            if (bam) {
                BAMAlignAux* aux = (BAMAlignAux*)to;
                aux->tag[0] = isID ? 'R' : from[1];
                aux->tag[1] = isID ? 'G' : from[2];
                aux->val_type = 'Z';
                from += 4; // skip \tXX:
                to = (char*)aux->value();
                while (*from && *from != '\t') {
                    *to++ = *from++;
                }
                *to++ = 0;
            } else {
                // turn \tXX: into \tXX:Z:, change ID to RG
                *to++ = *from++;
                if (isID) {
                    *to++ = 'R';
                    *to++ = 'G';
                    from += 2;
                } else {
                    *to++ = *from++;
                    *to++ = *from++;
                }
                *to++ = *from++;
                *to++ = 'Z';
                *to++ = ':';
                // copy string attribute
                while (*from && *from != '\t') {
                    *to++ = *from++;
                }
            }
        } else {
            from += 4;
            while (*from && *from != '\t') {
                from++;
            }
        }
    }
    readerContext->defaultReadGroupAux = buffer;
    readerContext->defaultReadGroupAuxLen = (int) (to - buffer);
}

    ReadWriterSupplier*
SAMFormat::getWriterSupplier(
    AlignerOptions* options,
    const Genome* genome) const
{
    DataWriterSupplier* dataSupplier;
    if (options->sortOutput) {
        DataWriter::FilterSupplier* filters = NULL;
        if (!options->noDuplicateMarking) {
            filters = DataWriterSupplier::samMarkDuplicates(genome);
        }
        dataSupplier = DataWriterSupplier::sorted(this, genome, DataWriterSupplier::generateSortIntermediateFilePathName(options), options->sortMemory * (1ULL << 30),
            options->numThreads, options->outputFile.fileName, filters, options->writeBufferSize, options->emitInternalScore, options->internalScoreTag);
    } else {
        dataSupplier = DataWriterSupplier::create(options->outputFile.fileName, options->writeBufferSize, options->emitInternalScore, options->internalScoreTag);
    }

    return ReadWriterSupplier::create(this, dataSupplier, genome, options->killIfTooSlow, options->emitInternalScore, options->internalScoreTag, options->ignoreAlignmentAdjustmentsForOm,
        options->matchReward, options->subPenalty, options->gapOpenPenalty, options->gapExtendPenalty, options->attachAlignmentTimes);
}

    bool
SAMFormat::writeHeader(
    const ReaderContext& context,
    char *header,
    size_t headerBufferSize,
    size_t *headerActualSize,
    bool sorted,
    int argc,
    const char **argv,
    const char *version,
    const char *rgLine,
	bool omitSQLines)	// Hacky option for Charles
    const
{
    char *commandLine;
	size_t commandLineSize = 0;
	for (int i = 0; i < argc; i++) {
		commandLineSize += strlen(argv[i]) + 1;	// +1 is either a space or the terminating null
	}

	commandLine = new char[commandLineSize];
	commandLine[0] = '\0';
	for (int i = 0; i < argc; i++) {
		strcat(commandLine,argv[i]);
		if (i != argc-1) {
			strcat(commandLine," ");
		}
	}

    size_t bytesConsumed = snprintf(header, headerBufferSize, "@HD\tVN:1.6\t%s\n%s%s@PG\tID:SNAP\tPN:SNAP\tCL:%s\tVN:%s\n", 
		sorted ? "SO:coordinate" : "GO:query",
        context.header == NULL ? (rgLine == NULL ? "@RG\tID:FASTQ\tSM:sample" : rgLine) : "",
        context.header == NULL ? "\n" : "",
        commandLine,version);

	delete [] commandLine;
	commandLine = NULL;
    if (bytesConsumed >= headerBufferSize) {
        //WriteErrorMessage("SAMWriter: header buffer too small\n");
        return false;
    }

    if (context.header != NULL) {
		bool hasRG = false;
        for (const char* p = context.header; p < context.header + context.headerLength; ) {
            const char* newline = strnchr(p, '\n', (context.header + context.headerLength) - p);
            if (newline == NULL) {
                newline = context.header + context.headerLength;
            }

            _ASSERT(newline - p >= 3);
            // skip @HD lines, and also @SQ lines if header does not match index
			hasRG |= strncmp(p, "@RG", 3) == 0;

            if (strncmp(p, "@HD", 3) != 0 &&
                    (context.headerMatchesIndex || strncmp(p, "@SQ", 3) != 0) &&
                    strncmp(p, "@PG\tID:SNAP\t", 12) != 0) {

                if (bytesConsumed + (newline - p) + 1 >= headerBufferSize) {
                    //WriteErrorMessage("SAMWriter: header buffer too small\n");
                    return false;
                }

                memcpy(header + bytesConsumed, p, (newline - p));
                * (header + bytesConsumed + (newline - p)) = '\n';
                bytesConsumed += (newline - p) + 1;
            }
            p = newline + 1;
        }

		if (! hasRG) {
			int n = snprintf(header + bytesConsumed, headerBufferSize - bytesConsumed, "%s\n",
				rgLine == NULL ? "@RG\tID:FASTQ\tSM:sample" : rgLine);
			if (n > headerBufferSize - bytesConsumed) {
				//WriteErrorMessage( "SAMWriter: header buffer too small\n");
                return false;
            }
			bytesConsumed += n;
		}
    }

    if ((context.header == NULL || ! context.headerMatchesIndex) && context.genome != NULL && !omitSQLines) {
        // Write an @SQ line for each chromosome / contig in the genome
        int numContigs = context.genome->getNumContigs();
        GenomeDistance genomeLen = context.genome->getCountOfBases();
        size_t originalBytesConsumed = bytesConsumed;

        for (int i = 0; i < numContigs; i++) {
            const Genome::Contig* contig = context.genome->getContigByOriginalContigNumber(OriginalContigNum(i));
            bytesConsumed += snprintf(header + bytesConsumed, headerBufferSize - bytesConsumed, "@SQ\tSN:%s\tLN:%llu%s\n", contig->name, contig->length - context.genome->getChromosomePadding(), contig->isALT ? "\tAH:*":"");

            if (bytesConsumed >= headerBufferSize) {
                // todo: increase buffer size (or change to write in batch
                bytesConsumed = originalBytesConsumed;
                //WriteErrorMessage("SAMWriter: header buffer too small, skipping @SQ lines\n");
                return false;
            }
        }
    }

    *headerActualSize = bytesConsumed;
    return true;
}

    void
SAMFormat::fillMateInfo(
    const Genome * genome,
    int& flags,
    Read * read,
    GenomeLocation genomeLocation,
    Direction direction,
    const char*& contigName,
    OriginalContigNum *contigIndex,
    GenomeDistance& positionInContig,
    _int64& templateLength,
    unsigned basesClippedBefore,
    bool firstInPair,
    bool alignedAsPair,
    Read * mate,
    GenomeLocation mateLocation,
    Direction mateDirection,
    const char*& matecontigName,
    OriginalContigNum *mateContigIndex,
    GenomeDistance& matePositionInContig,
    unsigned mateBasesClippedBefore,
    int myRefSpanFromCigar,
    int mateRefSpanFromCigar)
{
    flags |= SAM_MULTI_SEGMENT;
    flags |= (firstInPair ? SAM_FIRST_SEGMENT : SAM_LAST_SEGMENT);

    GenomeDistance extraBasesClippedBefore, mateExtraBasesClippedBefore;

    if (mateLocation != InvalidGenomeLocation) {
        const Genome::Contig *mateContig = genome->getContigForRead(mateLocation, mate->getDataLength(), &mateExtraBasesClippedBefore);
        mateLocation += mateExtraBasesClippedBefore;
        matecontigName = mateContig->name;
        *mateContigIndex = mateContig->originalContigNumber;
        matePositionInContig = mateLocation - mateContig->beginningLocation + 1;

        if (mateDirection == RC) {
            flags |= SAM_NEXT_REVERSED;
        }

        if (genomeLocation == InvalidGenomeLocation) {
            //
            // The SAM spec says that for paired reads where exactly one end is unmapped that the unmapped
            // half should just have RNAME and POS copied from the mate.
            //
            contigName = matecontigName;
            *contigIndex = *mateContigIndex;
            matecontigName = "=";
            positionInContig = matePositionInContig;
        }

    } else {
        flags |= SAM_NEXT_UNMAPPED;
        //
        // The mate's unmapped, so point it at us.
        //
        matecontigName = "=";
        *mateContigIndex = *contigIndex;
        matePositionInContig = positionInContig;
    }

    if (genomeLocation != InvalidGenomeLocation && mateLocation != InvalidGenomeLocation) {
        if (alignedAsPair) {
            flags |= SAM_ALL_ALIGNED;
        }

        // Add in any clipping from running off the end of reference
        const Genome::Contig* contig = genome->getContigForRead(genomeLocation, read->getDataLength(), &extraBasesClippedBefore);
        genomeLocation += extraBasesClippedBefore;

        // Also compute the length of the whole paired-end string whose ends we saw. This is slightly
        // tricky because (a) we may have clipped some bases before/after each end and (b) we need to
        // give a signed result based on whether our read is first or second in the pair.
        GenomeLocation myStart = genomeLocation - basesClippedBefore - extraBasesClippedBefore;
        GenomeLocation myEnd = genomeLocation + myRefSpanFromCigar;
        GenomeLocation mateStart = mateLocation - mateBasesClippedBefore - mateExtraBasesClippedBefore;
        GenomeLocation mateEnd = mateLocation + mateRefSpanFromCigar;

        *contigIndex = contig->originalContigNumber;

        if (myStart < mateStart) {
            if (direction == FORWARD) {
                if (mateDirection == RC) {
                    templateLength = mateEnd - myStart; // FR
                } else {
                    templateLength = mateStart - myStart; // FF
                }
            } else {
                if (mateDirection == FORWARD) {
                    templateLength = mateStart - myEnd; // RF
                } else {
                    templateLength = mateEnd - myEnd; // RR
                }
            }
        } else {
            if (direction == RC) {
                if (mateDirection == FORWARD) {
                    templateLength = -(myEnd - mateStart);
                } else {
                    templateLength = -(myEnd - mateEnd);
                }
            } else {
                if (mateDirection == FORWARD) {
                    templateLength = -(myStart - mateStart);
                } else {
                    templateLength = -(myStart - mateEnd);
                }
            }
        }
    }

    if (contigName == matecontigName) {
        matecontigName = "=";     // SAM Spec says to do this when they're equal (and not *, which won't happen because this is a pointer, not string, compare)
    }
}

    bool
SAMFormat::createSAMLine(
    const Genome * genome,
    // output data
    char* data,
    char* quality,
    GenomeDistance dataSize,
    const char*& contigName,
    OriginalContigNum *contigIndex,
    int& flags,
    GenomeDistance& positionInContig,
    int& mapQuality,
    const char*& matecontigName,
    OriginalContigNum* mateContigIndex,
    GenomeDistance& matePositionInContig,
    _int64& templateLength,
    unsigned& fullLength,
    const char*& clippedData,
    const char*& clippedQuality,
    unsigned& clippedLength,
    unsigned& basesClippedBefore,
    unsigned& basesClippedAfter,
    unsigned& mateBasesClippedBefore,
    unsigned& mateBasesClippedAfter,
    const char*& FASTQComment,
    unsigned &FASTQCommentLength,
    // input data
    size_t& qnameLen,
    Read * read,
    AlignmentResult result,
    GenomeLocation genomeLocation,
    Direction direction,
    bool secondaryAlignment,
    bool supplementaryAlignment,
    bool useM,
    bool hasMate,
    bool firstInPair,
    bool alignedAsPair,
    Read * mate,
    AlignmentResult mateResult,
    GenomeLocation mateLocation,
    Direction mateDirection,
    GenomeDistance *extraBasesClippedBefore,
    int bpClippedBefore,
    int bpClippedAfter,
    int mateBpClippedBefore,
    int mateBpClippedAfter)
{
    contigName = "*";
    *contigIndex = OriginalContigNum(-1);
    positionInContig = 0;
    const char *cigar = "*";
    templateLength = 0;

    if (secondaryAlignment) {
        flags |= SAM_SECONDARY;
    }

    if (supplementaryAlignment) {
        flags |= SAM_SUPPLEMENTARY;
    }

    if (0 == qnameLen) {
        qnameLen = read->getIdLength();
    }

    if (read->getFASTQComment() == NULL || read->getFASTQCommentLength() == 0) {
        FASTQComment = "";
        FASTQCommentLength = 0;
    } else {
        FASTQComment = read->getFASTQComment();
        FASTQCommentLength = read->getFASTQCommentLength();
    }

    //
    // If the aligner said it didn't find anything, treat it as such.  Sometimes it will emit the
    // best match that it found, even if it's not within the maximum edit distance limit (but will
    // then say NotFound).  Here, we force that to be SAM_UNMAPPED.
    //
    if (NotFound == result) {
        genomeLocation = InvalidGenomeLocation;
    }

    if (InvalidGenomeLocation == genomeLocation) {
        //
        // If it's unmapped, then always emit it in the forward direction.  This is necessary because we don't even include
        // the SAM_REVERSE_COMPLEMENT flag for unmapped reads, so there's no way to tell that we reversed it.
        //
        direction = FORWARD;
    }

    // Write the data and quality strings. If the read is reverse complemented, these need to
    // be backwards from the original read. Also, both need to be unclipped.
    clippedLength = read->getDataLength();
    fullLength = read->getUnclippedLength();
    if (fullLength > dataSize) {
        return false;
    }

    if (direction == RC) {
        for (unsigned i = 0; i < fullLength; i++) {
            data[fullLength - 1 - i] = COMPLEMENT[read->getUnclippedData()[i]];
            quality[fullLength - 1 - i] = read->getUnclippedQuality()[i];
        }

        clippedData = &data[fullLength - clippedLength - read->getFrontClippedLength()];
        clippedQuality = &quality[fullLength - clippedLength - read->getFrontClippedLength()];
        basesClippedBefore = fullLength - clippedLength - read->getFrontClippedLength();
        basesClippedAfter = read->getFrontClippedLength();
    } else {
        memcpy(data, read->getUnclippedData(), read->getUnclippedLength());
        memcpy(quality, read->getUnclippedQuality(), read->getUnclippedLength());
        clippedData = read->getData();
        clippedQuality = read->getQuality();
        basesClippedBefore = read->getFrontClippedLength();
        basesClippedAfter = fullLength - clippedLength - basesClippedBefore;
    }

    // Add soft-clipping from seed extension
    basesClippedBefore += bpClippedBefore;
    basesClippedAfter += bpClippedAfter;
    clippedData += bpClippedBefore;
    clippedQuality += bpClippedBefore;
    clippedLength -= (bpClippedBefore + bpClippedAfter);

    int editDistance = -1;
    if (genomeLocation != InvalidGenomeLocation) {
        if (direction == RC) {
            flags |= SAM_REVERSE_COMPLEMENT;
        }

        const Genome::Contig *contig = genome->getContigForRead(genomeLocation, read->getDataLength(), extraBasesClippedBefore);
        _ASSERT(NULL != contig && contig->length > genome->getChromosomePadding());
        genomeLocation += *extraBasesClippedBefore;

        contigName = contig->name;
        *contigIndex = contig->originalContigNumber;
        positionInContig = genomeLocation - contig->beginningLocation + 1; // SAM is 1-based
        mapQuality = max(0, min(70, mapQuality));       // FIXME: manifest constant.
    } else {
        flags |= SAM_UNMAPPED;
        flags &= ~SAM_REVERSE_COMPLEMENT;
        mapQuality = 0;
        *extraBasesClippedBefore = 0;
        *contigIndex = OriginalContigNum(-1);
        positionInContig = 0;
    }

    return true;
}

    bool 
SAMFormat::writePairs(
    const ReaderContext& context,
    LandauVishkinWithCigar * lv,
    AffineGapVectorizedWithCigar * ag,
    bool useAffineGap,
    char * buffer,
    size_t bufferSpace,
    size_t * spaceUsed,
    size_t * qnameLen,
    Read ** reads,
    GenomeLocation * locations,
    PairedAlignmentResult * result,
    bool isSecondary,
    bool emitInternalScore,
    char *internalScoreTag,
    bool attachAlignmentTime,
    int * writeOrder,
    int* cumulativePositiveAddFrontClipping,
    bool * secondReadLocationChanged,
    bool * outOfSpace) const
{

    const int MAX_READ = MAX_READ_LENGTH;
    const int cigarBufSize = MAX_READ * 2;
    char cigarBuf[NUM_READS_PER_PAIR][cigarBufSize];

    const int cigarBufWithClippingSize = MAX_READ * 2 + 32;
    char cigarBufWithClipping[NUM_READS_PER_PAIR][cigarBufWithClippingSize];

    int flags[NUM_READS_PER_PAIR] = {0, 0};
    const char *contigName[NUM_READS_PER_PAIR] = {"*", "*"};
    OriginalContigNum contigIndex[NUM_READS_PER_PAIR] = { OriginalContigNum(-1), OriginalContigNum(-1)};
    GenomeDistance positionInContig[NUM_READS_PER_PAIR] = {0, 0};
    const char *mateContigName[NUM_READS_PER_PAIR] = {"*", "*"};
    OriginalContigNum mateContigIndex[NUM_READS_PER_PAIR] = { OriginalContigNum(-1), OriginalContigNum(-1)};
    GenomeDistance matePositionInContig[NUM_READS_PER_PAIR] = {0, 0};
    const char *cigar[NUM_READS_PER_PAIR] = {"*", "*"};
    _int64 templateLength[NUM_READS_PER_PAIR] = {0, 0};
    int refSpanFromCigar[NUM_READS_PER_PAIR] = {0, 0};

    char data[NUM_READS_PER_PAIR][MAX_READ];
    char quality[NUM_READS_PER_PAIR][MAX_READ];

    const char* clippedData[NUM_READS_PER_PAIR];
    const char* clippedQuality[NUM_READS_PER_PAIR];
    unsigned fullLength[NUM_READS_PER_PAIR];
    unsigned clippedLength[NUM_READS_PER_PAIR];
    unsigned basesClippedBefore[NUM_READS_PER_PAIR];
    unsigned basesClippedAfter[NUM_READS_PER_PAIR];
    GenomeDistance extraBasesClippedBefore[NUM_READS_PER_PAIR];   // Clipping added if we align before the beginning of a chromosome
    int editDistance[NUM_READS_PER_PAIR] = {-1, -1};
    const char* FASTQComment[NUM_READS_PER_PAIR];
    unsigned FASTQCommentLength[NUM_READS_PER_PAIR];

    // Create SAM entry and compute CIGAR
    for (int firstOrSecond = 0; firstOrSecond < NUM_READS_PER_PAIR; firstOrSecond++) {
        int whichRead = writeOrder[firstOrSecond];
        Read* read = reads[whichRead];
        bool firstInPair = writeOrder[firstOrSecond] == 0;

        int addFrontClipping;
        do {
            addFrontClipping = 0;

            if (!createSAMLine(context.genome, data[whichRead], quality[whichRead], MAX_READ, contigName[whichRead], &contigIndex[whichRead],
                flags[whichRead], positionInContig[whichRead], result->mapq[whichRead], contigName[1 - whichRead], &contigIndex[1 - whichRead],
                positionInContig[1 - whichRead], templateLength[whichRead],
                fullLength[whichRead], clippedData[whichRead], clippedQuality[whichRead], clippedLength[whichRead], basesClippedBefore[whichRead], basesClippedAfter[whichRead],
                basesClippedBefore[1 - whichRead], basesClippedAfter[1 - whichRead], FASTQComment[whichRead], FASTQCommentLength[whichRead], qnameLen[whichRead], reads[whichRead], 
                result->status[whichRead], locations[whichRead], result->direction[whichRead], isSecondary, result->supplementary[whichRead], useM,
                true, firstInPair, result->alignedAsPair, reads[1 - whichRead], result->status[1 - whichRead], locations[1 - whichRead], result->direction[1 - whichRead], 
                &extraBasesClippedBefore[whichRead], result->basesClippedBefore[whichRead], result->basesClippedAfter[whichRead], 
                result->basesClippedBefore[1 - whichRead], result->basesClippedAfter[1 - whichRead]))
            {
                return false;
            }

            if (locations[whichRead] != InvalidGenomeLocation) {
                if (useAffineGap && (result->usedAffineGapScoring[whichRead] || result->score[whichRead] > 0)) {
                    cigar[whichRead] = computeCigarString(context.genome, ag, cigarBuf[whichRead], cigarBufSize, cigarBufWithClipping[whichRead], cigarBufWithClippingSize,
                        clippedData[whichRead], clippedQuality[whichRead], clippedLength[whichRead], result->score[whichRead], basesClippedBefore[whichRead], extraBasesClippedBefore[whichRead], basesClippedAfter[whichRead],
                        read->getOriginalFrontHardClipping(), read->getOriginalBackHardClipping(), locations[whichRead], result->direction[whichRead], useM,
                        &editDistance[whichRead], &addFrontClipping, &refSpanFromCigar[whichRead]);

                    if (addFrontClipping != 0) {
                        *secondReadLocationChanged = firstOrSecond == 1;
                        const Genome::Contig *originalContig = context.genome->getContigAtLocation(locations[whichRead]);
                        const Genome::Contig *newContig = context.genome->getContigAtLocation(locations[whichRead] + addFrontClipping);
                        if (newContig != originalContig || NULL == newContig || locations[whichRead] + addFrontClipping > originalContig->beginningLocation + originalContig->length - context.genome->getChromosomePadding()) {
                            //
                            // Altering this would push us over a contig boundary.  Just give up on the read.
                            //
                            result->status[whichRead] = NotFound;
                            result->location[whichRead] = InvalidGenomeLocation;
                            locations[whichRead] = InvalidGenomeLocation;
                            cigar[whichRead] = "*";
                            editDistance[whichRead] = -1;
                            result->direction[whichRead] = FORWARD;
                        } else {
                            if (addFrontClipping < 0) { // Insertion (soft-clip)
                                cumulativePositiveAddFrontClipping[firstOrSecond] += addFrontClipping;
                                if (result->direction[whichRead] == FORWARD) {
                                    reads[whichRead]->setAdditionalFrontClipping(-cumulativePositiveAddFrontClipping[firstOrSecond]);
                                } else {
                                    reads[whichRead]->setAdditionalBackClipping(-cumulativePositiveAddFrontClipping[firstOrSecond]);
                                }
                            } else { // Deletion
                                locations[whichRead] += addFrontClipping;
                            }
                        }
                    }
                } else {
                    cigar[whichRead] = computeCigarString(context.genome, lv, cigarBuf[whichRead], cigarBufSize, cigarBufWithClipping[whichRead], cigarBufWithClippingSize,
                        clippedData[whichRead], clippedLength[whichRead], basesClippedBefore[whichRead], extraBasesClippedBefore[whichRead], basesClippedAfter[whichRead], 
                        read->getOriginalFrontHardClipping(), read->getOriginalBackHardClipping(), locations[whichRead], result->direction[whichRead], useM,
                        &editDistance[whichRead], &addFrontClipping, &refSpanFromCigar[whichRead]);

                    if (addFrontClipping != 0) {
                        *secondReadLocationChanged = firstOrSecond == 1;
                        const Genome::Contig *originalContig = context.genome->getContigAtLocation(locations[whichRead]);
                        const Genome::Contig *newContig = context.genome->getContigAtLocation(locations[whichRead] + addFrontClipping);
                        if (newContig != originalContig || NULL == newContig || locations[whichRead] + addFrontClipping > originalContig->beginningLocation + originalContig->length - context.genome->getChromosomePadding()) {
                            //
                            // Altering this would push us over a contig boundary.  Just give up on the read.
                            //
                            result->status[whichRead] = NotFound;
                            result->location[whichRead] = InvalidGenomeLocation;
                            locations[whichRead] = InvalidGenomeLocation;
                            cigar[whichRead] = "*";
                            editDistance[whichRead] = -1;
                            result->direction[whichRead] = FORWARD;
                        } else {
                            if (addFrontClipping > 0) {
                                cumulativePositiveAddFrontClipping[firstOrSecond] += addFrontClipping;
                                reads[whichRead]->setAdditionalFrontClipping(cumulativePositiveAddFrontClipping[firstOrSecond]);
                            }
                            locations[whichRead] += addFrontClipping;
                        }
                    }
                }
            }
        } while (addFrontClipping != 0);
	}

    // Fill mate information
    for (int firstOrSecond = 0; firstOrSecond < NUM_READS_PER_PAIR; firstOrSecond++) {
        int whichRead = writeOrder[firstOrSecond];
        bool firstInPair = writeOrder[firstOrSecond] == 0;
        fillMateInfo(context.genome, flags[whichRead], reads[whichRead], locations[whichRead], result->direction[whichRead], 
            contigName[whichRead], &contigIndex[whichRead], positionInContig[whichRead], templateLength[whichRead], basesClippedBefore[whichRead],
            firstInPair, result->alignedAsPair, reads[1 - whichRead], locations[1 - whichRead], result->direction[1 - whichRead],
            mateContigName[whichRead], &mateContigIndex[whichRead], matePositionInContig[whichRead], basesClippedBefore[1 - whichRead],
            refSpanFromCigar[whichRead], refSpanFromCigar[1 - whichRead]);
    }

    // Write the SAM entry, which requires the following fields:
    //
    // 1. QNAME: Query name of the read or the read pair
    // 2. FLAG: Bitwise flag (pairing, strand, mate strand, etc.)
    // 3. RNAME: Reference sequence name
    // 4. POS: 1-Based leftmost position of clipped alignment
    // 5. MAPQ: Mapping quality (Phred-scaled)
    // 6. CIGAR: Extended CIGAR string (operations: MIDNSHP)
    // 7. MRNM: Mate reference name (‘=’ if same as RNAME)
    // 8. MPOS: 1-based leftmost mate position
    // 9. ISIZE: Inferred insert size
    // 10. SEQQuery: Sequence on the same strand as the reference
    // 11. QUAL: Query quality (ASCII-33=Phred base quality)    

    for (int firstOrSecond = 0; firstOrSecond < NUM_READS_PER_PAIR; firstOrSecond++) {
        int whichRead = writeOrder[firstOrSecond];
        Read* read = reads[whichRead];
        //
        // Some FASTQ files have spaces in their ID strings, which is illegal in SAM.  Just truncate them at the space.
        //
        const char *firstSpace = strnchr(read->getId(),' ',qnameLen[whichRead]);
        if (NULL != firstSpace) {
            qnameLen[whichRead] = (unsigned)(firstSpace - read->getId());
        }

        const int nmStringSize = 30;// Big enough that it won't buffer overflow regardless of the value of editDistance
        char nmString[nmStringSize];  
        snprintf(nmString, nmStringSize, "\tNM:i:%d",editDistance[whichRead]);

        unsigned auxLen;
        bool auxSAM;
        char* aux = read->getAuxiliaryData(&auxLen, &auxSAM);
        static bool warningPrinted = false;
        const char* readGroupSeparator = "";
        const char* readGroupString = "";
        if (aux != NULL && (! auxSAM)) {
            if (! warningPrinted) {
                WriteErrorMessage( "warning: translating optional fields from BAM->SAM not yet implemented, optional fields will not be included in output\n");
                warningPrinted = true;
            }
            if (read->getReadGroup() == READ_GROUP_FROM_AUX) {
                for (BAMAlignAux* bamAux = (BAMAlignAux*) aux; (char*) bamAux < aux + auxLen; bamAux = bamAux->next()) {
                    if (bamAux->tag[0] == 'R' && bamAux->tag[1] == 'G' && bamAux->val_type == 'Z') {
                        readGroupSeparator = "\tRG:Z:";
                        readGroupString = (char*) bamAux->value();
                        break;
                    }
                }
            }
            aux = NULL;
            auxLen = 0;
        }

        const char* rglineAux = "";
        int rglineAuxLen = 0;
        if (read->getReadGroup() != NULL && read->getReadGroup() != READ_GROUP_FROM_AUX) {
            if (*readGroupString == 0 || strcmp(readGroupString, context.defaultReadGroup) == 0) {
                readGroupSeparator = "";
                readGroupString = "";
                rglineAux = context.defaultReadGroupAux;
                rglineAuxLen = context.defaultReadGroupAuxLen;
            } else {
                readGroupSeparator = "\tRG:Z:";
                readGroupString = read->getReadGroup();
            }
        }

        const int internalScoreBufferSize = 100;    // Should be plenty for \tXX:i:%d
        char internalScoreBuffer[internalScoreBufferSize];
        if (emitInternalScore) {
            int charsInInternalScore = snprintf(internalScoreBuffer, internalScoreBufferSize - 1, "\t%s:i:%d", internalScoreTag, (flags[whichRead] & SAM_UNMAPPED) ? -1 : result->scorePriorToClipping[whichRead]);
            if (charsInInternalScore >= internalScoreBufferSize) {
                WriteErrorMessage("SAMFormat::writeRead overran internal buffer for internal score tag, which is kind of surprising.  %d\n", charsInInternalScore);
            }
        } else {
            internalScoreBuffer[0] = '\0';
        }

        const int alignmentTimeBufferSize = 100;    // Should be plenty for \tAT:i:%d
        char alignmentTimeBuffer[alignmentTimeBufferSize];
        if (attachAlignmentTime) {
            int alignmentTimeInMicroseconds;
            if (result->alignmentTimeInNanoseconds / 1000 > MAXINT32) {
                alignmentTimeInMicroseconds = MAXINT32;
            } else {
                alignmentTimeInMicroseconds = (int)(result->alignmentTimeInNanoseconds / 1000);
            }
            int charsInAlignmentTime = snprintf(alignmentTimeBuffer, alignmentTimeBufferSize - 1, "\tAT:i:%d", alignmentTimeInMicroseconds);
            if (charsInAlignmentTime >= alignmentTimeBufferSize) {
                WriteErrorMessage("SAMFormat::writeRead overran internal buffer for alignment time tag, which is kind of surprising.  %d\n", charsInAlignmentTime);
            }
        } else {
            alignmentTimeBuffer[0] = '\0';
        }


        // QS
        int mqs = 0;
        _uint8* p = (_uint8*)quality[1 - whichRead];
        for (unsigned i = 0; i < fullLength[1 - whichRead]; i++) {
            int q = *p++;
            q -= '!';
            // Picard MarkDup uses a score threshold of 15 (default)
            mqs += (q >= 15) ? (q != 255) * q : 0; // avoid branch?
        }
        const int mqsStringSize = 30;// Big enough that it won't buffer overflow regardless of the value of mateQualityScore
        char mqsString[nmStringSize];
        snprintf(mqsString, mqsStringSize, "\tQS:i:%d", mqs);

        // LB
        const int libraryStringSize = 512;
        char libraryString[libraryStringSize];
        const char* library = read->getLibrary();
        if (library != NULL) {
            size_t libraryLength = read->getLibraryLength();
            if (libraryLength >= libraryStringSize) {
                WriteErrorMessage("LB field too long\n");
                soft_exit(1);
            }
            snprintf(libraryString, libraryStringSize, "\tLB:Z:%.*s", (int)libraryLength, library);
        } else {
            libraryString[0] = '\0';
        }

        int charsInString = snprintf(buffer, bufferSpace, "%.*s\t%d\t%s\t%llu\t%d\t%s\t%s\t%llu\t%d\t%.*s\t%.*s%s%.*s%s%s\tPG:Z:SNAP%s%.*s%s%s%s%s%s%.*s\n",
            (unsigned)qnameLen[whichRead], read->getId(),
            flags[whichRead],
            contigName[whichRead],
            positionInContig[whichRead],
            result->mapq[whichRead],
            cigar[whichRead],
            mateContigName[whichRead],
            matePositionInContig[whichRead],
            (_int32) templateLength[whichRead],
            fullLength[whichRead], data[whichRead],
            fullLength[whichRead], quality[whichRead],
            aux != NULL ? "\t" : "", auxLen, aux != NULL ? aux : "",
            readGroupSeparator, readGroupString,
            nmString, rglineAuxLen, rglineAux,
            internalScoreBuffer,
            alignmentTimeBuffer,
            mqsString,
            libraryString,
            (FASTQCommentLength[whichRead] == 0) ? "" : "\t",
            FASTQCommentLength[whichRead],
            FASTQComment[whichRead]);

        if (charsInString > bufferSpace) {
            //
            // Out of buffer space.
            //
            *outOfSpace = true;
            return false;
        } else if (charsInString == bufferSpace) {
            buffer[bufferSpace-1] = '\n'; // overwrite trailing null with newline
        }

        if (NULL != spaceUsed) {
            spaceUsed[firstOrSecond] = charsInString;
        }

        buffer += charsInString;
        bufferSpace -= charsInString;
    }
    return true;
}

    bool
SAMFormat::writeRead(
    const ReaderContext& context,
    LandauVishkinWithCigar * lv,
    char * buffer,
    size_t bufferSpace, 
    size_t * spaceUsed,
    size_t qnameLen,
    Read * read,
    AlignmentResult result, 
    int mapQuality,
    GenomeLocation genomeLocation,
    Direction direction,
    bool secondaryAlignment,
    bool supplementaryAlignment,
    int * o_addFrontClipping,
    int internalScore,
    bool emitInternalScore,
    char *internalScoreTag,
    bool attachAlignmentTime,
    _int64 alignmentTimeInNanoseconds,
    int bpClippedBefore,
    int bpClippedAfter,
    bool hasMate,
    bool firstInPair,
    Read * mate, 
    AlignmentResult mateResult,
    GenomeLocation mateLocation,
    Direction mateDirection,
    bool alignedAsPair,
    int mateBpClippedBefore,
    int mateBpClippedAfter
    ) const
{
    const int MAX_READ = MAX_READ_LENGTH;
    const int cigarBufSize = MAX_READ * 2;
    char cigarBuf[cigarBufSize];

    const int cigarBufWithClippingSize = MAX_READ * 2 + 32;
    char cigarBufWithClipping[cigarBufWithClippingSize];

    int flags = 0;
    const char *contigName = "*";
    OriginalContigNum contigIndex = OriginalContigNum(-1);
    GenomeDistance positionInContig = 0;
    const char *cigar = "*";
    const char *matecontigName = "*";
    OriginalContigNum mateContigIndex = OriginalContigNum(-1);
    GenomeDistance matePositionInContig = 0;
    _int64 templateLength = 0;
    int refSpanFromCigar = 0;

    char data[MAX_READ];
    char quality[MAX_READ];

    const char* clippedData;
    const char* clippedQuality;
    unsigned fullLength;
    unsigned clippedLength;
    unsigned basesClippedBefore, mateBasesClippedBefore;
    GenomeDistance extraBasesClippedBefore;   // Clipping added if we align before the beginning of a chromosome
    unsigned basesClippedAfter, mateBasesClippedAfter;
    int editDistance = -1;

    const char* FASTQComment;
    unsigned FASTQCommentLength;

    *o_addFrontClipping = 0;

	if (!createSAMLine(context.genome, data, quality, MAX_READ, contigName, &contigIndex,
        flags, positionInContig, mapQuality, matecontigName, &mateContigIndex, matePositionInContig, templateLength,
        fullLength, clippedData, clippedQuality, clippedLength, basesClippedBefore, basesClippedAfter, mateBasesClippedBefore, mateBasesClippedAfter,
        FASTQComment, FASTQCommentLength, qnameLen, read, result, genomeLocation, direction, secondaryAlignment, supplementaryAlignment, useM,
        hasMate, firstInPair, alignedAsPair, mate, mateResult, mateLocation, mateDirection, 
        &extraBasesClippedBefore, bpClippedBefore, bpClippedAfter, mateBpClippedBefore, mateBpClippedAfter))
    {
        return false;
    }

	if (genomeLocation != InvalidGenomeLocation) {
		cigar = computeCigarString(context.genome, lv, cigarBuf, cigarBufSize, cigarBufWithClipping, cigarBufWithClippingSize,
			clippedData, clippedLength, basesClippedBefore, extraBasesClippedBefore, basesClippedAfter, 
			read->getOriginalFrontHardClipping(), read->getOriginalBackHardClipping(), genomeLocation, direction, useM,
			&editDistance, o_addFrontClipping, &refSpanFromCigar);
		if (*o_addFrontClipping != 0) {
			return false;
		}
	}


    // Write the SAM entry, which requires the following fields:
    //
    // 1. QNAME: Query name of the read or the read pair
    // 2. FLAG: Bitwise flag (pairing, strand, mate strand, etc.)
    // 3. RNAME: Reference sequence name
    // 4. POS: 1-Based leftmost position of clipped alignment
    // 5. MAPQ: Mapping quality (Phred-scaled)
    // 6. CIGAR: Extended CIGAR string (operations: MIDNSHP)
    // 7. MRNM: Mate reference name (‘=’ if same as RNAME)
    // 8. MPOS: 1-based leftmost mate position
    // 9. ISIZE: Inferred insert size
    // 10. SEQQuery: Sequence on the same strand as the reference
    // 11. QUAL: Query quality (ASCII-33=Phred base quality)    

    //
    // Some FASTQ files have spaces in their ID strings, which is illegal in SAM.  Just truncate them at the space.
    //
    const char *firstSpace = strnchr(read->getId(),' ',qnameLen);
    if (NULL != firstSpace) {
        qnameLen = (unsigned)(firstSpace - read->getId());
    }

    const int nmStringSize = 30;// Big enough that it won't buffer overflow regardless of the value of editDistance
    char nmString[nmStringSize];  
    snprintf(nmString, nmStringSize, "\tNM:i:%d",editDistance);

    unsigned auxLen;
    bool auxSAM;
    char* aux = read->getAuxiliaryData(&auxLen, &auxSAM);
    static bool warningPrinted = false;
    const char* readGroupSeparator = "";
    const char* readGroupString = "";
    if (aux != NULL && (! auxSAM)) {
        if (! warningPrinted) {
            WriteErrorMessage( "warning: translating optional fields from BAM->SAM not yet implemented, optional fields will not be included in output\n");
            warningPrinted = true;
        }
        if (read->getReadGroup() == READ_GROUP_FROM_AUX) {
            for (BAMAlignAux* bamAux = (BAMAlignAux*) aux; (char*) bamAux < aux + auxLen; bamAux = bamAux->next()) {
                if (bamAux->tag[0] == 'R' && bamAux->tag[1] == 'G' && bamAux->val_type == 'Z') {
                    readGroupSeparator = "\tRG:Z:";
                    readGroupString = (char*) bamAux->value();
                    break;
                }
            }
        }
        aux = NULL;
        auxLen = 0;
    }
    const char* rglineAux = "";
    int rglineAuxLen = 0;
    if (read->getReadGroup() != NULL && read->getReadGroup() != READ_GROUP_FROM_AUX) {
        if (*readGroupString == 0 || strcmp(readGroupString, context.defaultReadGroup) == 0) {
            readGroupSeparator = "";
            readGroupString = "";
            rglineAux = context.defaultReadGroupAux;
            rglineAuxLen = context.defaultReadGroupAuxLen;
        } else {
            readGroupSeparator = "\tRG:Z:";
            readGroupString = read->getReadGroup();
        }
    }

    const int internalScoreBufferSize = 100;    // Should be plenty for \tXX:i:%d
    char internalScoreBuffer[internalScoreBufferSize];
    if (emitInternalScore) {
        int charsInInternalScore = snprintf(internalScoreBuffer, internalScoreBufferSize - 1, "\t%s:i:%d", internalScoreTag, (flags & SAM_UNMAPPED) ? -1 : internalScore);
        if (charsInInternalScore >= internalScoreBufferSize) {
            WriteErrorMessage("SAMFormat::writeRead overran internal buffer for internal score tag, which is kind of surprising.  %d\n", charsInInternalScore);
        }
    } else {
        internalScoreBuffer[0] = '\0';
    }

    const int alignmentTimeBufferSize = 100;    // Should be plenty for AT:i:123456789
    char alignmentTimeBuffer[alignmentTimeBufferSize];
    if (attachAlignmentTime) {
        _int64 alignmentTimeInMicroseconds = (alignmentTimeInNanoseconds + 500) / 1000;
        if (alignmentTimeInMicroseconds >= MAXINT32) {  // MAXINT is about 2 billion, so this would be ~35 minutes for one read (pair)
            alignmentTimeInMicroseconds = 0;
        };
        int charsInAligmentTime = snprintf(alignmentTimeBuffer, alignmentTimeBufferSize - 1, "\tAT:i:%d", (int)alignmentTimeInMicroseconds);
        if (charsInAligmentTime >= alignmentTimeBufferSize) {
            WriteErrorMessage("SAMFormat::writeRead overran internal buffer for alignment time tag, which is kind of surprising.  %d\n", charsInAligmentTime);
        }
    } else {
        alignmentTimeBuffer[0] = '\0';
    }

    int charsInString = snprintf(buffer, bufferSpace, "%.*s\t%d\t%s\t%llu\t%d\t%s\t%s\t%llu\t%lld\t%.*s\t%.*s%s%.*s%s%s\tPG:Z:SNAP%s%.*s%s%s%s%.*s\n",
        (unsigned)qnameLen, read->getId(),
        flags,
        contigName,
        positionInContig,
        mapQuality,
        cigar,
        matecontigName,
        matePositionInContig,
        templateLength,
        fullLength, data,
        fullLength, quality,
        aux != NULL ? "\t" : "", auxLen, aux != NULL ? aux : "",
        readGroupSeparator, readGroupString,
        nmString, rglineAuxLen, rglineAux,
        internalScoreBuffer,
        alignmentTimeBuffer,
        (FASTQCommentLength == 0) ? "" : "\t",
        FASTQCommentLength,
        FASTQComment
       );

    if (charsInString > bufferSpace) {
        //
        // Out of buffer space.
        //
        return false;
    } else if (charsInString == bufferSpace) {
      buffer[bufferSpace-1] = '\n'; // overwrite trailing null with newline
    }


    if (NULL != spaceUsed) {
        *spaceUsed = charsInString;
    }
    return true;
}

    bool
SAMFormat::writeRead(
    const ReaderContext& context,
    AffineGapVectorizedWithCigar * ag,
    char * buffer,
    size_t bufferSpace,
    size_t * spaceUsed,
    size_t qnameLen,
    Read * read,
    AlignmentResult result,
    int mapQuality,
    GenomeLocation genomeLocation,
    Direction direction,
    bool secondaryAlignment,
    bool supplementaryAlignment,
    int * o_addFrontClipping,
    int score,
    int internalScore,
    bool emitInternalScore,
    char *internalScoreTag,
    bool attachAlignmentTime,
    _int64 alignmentTimeInNanoseconds,
    int bpClippedBefore,
    int bpClippedAfter,
    bool hasMate,
    bool firstInPair,
    Read * mate,
    AlignmentResult mateResult,
    GenomeLocation mateLocation,
    Direction mateDirection,
    bool alignedAsPair,
    int mateBpClippedBefore,
    int mateBpClippedAfter
) const
{
    const int MAX_READ = MAX_READ_LENGTH;
    const int cigarBufSize = MAX_READ * 2;
    char cigarBuf[cigarBufSize];

    const int cigarBufWithClippingSize = MAX_READ * 2 + 32;
    char cigarBufWithClipping[cigarBufWithClippingSize];

    int flags = 0;
    const char *contigName = "*";
    OriginalContigNum contigIndex = -1;
    GenomeDistance positionInContig = 0;
    const char *cigar = "*";
    const char *matecontigName = "*";
    OriginalContigNum mateContigIndex = -1;
    GenomeDistance matePositionInContig = 0;
    _int64 templateLength = 0;
    int refSpanFromCigar = 0;

    char data[MAX_READ];
    char quality[MAX_READ];

    const char* clippedData;
    const char* clippedQuality;
    unsigned fullLength;
    unsigned clippedLength;
    unsigned basesClippedBefore, mateBasesClippedBefore;
    GenomeDistance extraBasesClippedBefore;   // Clipping added if we align before the beginning of a chromosome
    unsigned basesClippedAfter, mateBasesClippedAfter;
    int editDistance = -1;

    const char* FASTQComment;
    unsigned FASTQCommentLength;

    *o_addFrontClipping = 0;

    if (!createSAMLine(context.genome, data, quality, MAX_READ, contigName, &contigIndex,
        flags, positionInContig, mapQuality, matecontigName, &mateContigIndex, matePositionInContig, templateLength,
        fullLength, clippedData, clippedQuality, clippedLength, basesClippedBefore, basesClippedAfter, mateBasesClippedBefore, mateBasesClippedAfter,
        FASTQComment, FASTQCommentLength, qnameLen, read, result, genomeLocation, direction, secondaryAlignment, supplementaryAlignment, useM,
        hasMate, firstInPair, alignedAsPair, mate, mateResult, mateLocation, mateDirection,
        &extraBasesClippedBefore, bpClippedBefore, bpClippedAfter, mateBpClippedBefore, mateBpClippedAfter))
    {
        return false;
    }

    if (extraBasesClippedBefore != 0) {
        *o_addFrontClipping = (int)extraBasesClippedBefore;
        return false;
    }

    if (genomeLocation != InvalidGenomeLocation) {
        cigar = computeCigarString(context.genome, ag, cigarBuf, cigarBufSize, cigarBufWithClipping, cigarBufWithClippingSize,
            clippedData, clippedQuality, clippedLength, score, basesClippedBefore, extraBasesClippedBefore, basesClippedAfter,
            read->getOriginalFrontHardClipping(), read->getOriginalBackHardClipping(), genomeLocation, direction, useM,
            &editDistance, o_addFrontClipping, &refSpanFromCigar);
        // Uncomment for debug
        if (editDistance == -1) {
            const char* read_data = read->getUnclippedData();
            const char* readId = read->getId();
            for (unsigned i = 0; i < read->getIdLength(); ++i) {
                printf("%c", readId[i]);
            }
            printf(",");
            for (unsigned i = 0; i < read->getUnclippedLength(); ++i) {
                printf("%c", read_data[i]);
            }
            printf("\n");
        }
        if (*o_addFrontClipping != 0) {
            return false;
        }
    }


    // Write the SAM entry, which requires the following fields:
    //
    // 1. QNAME: Query name of the read or the read pair
    // 2. FLAG: Bitwise flag (pairing, strand, mate strand, etc.)
    // 3. RNAME: Reference sequence name
    // 4. POS: 1-Based leftmost position of clipped alignment
    // 5. MAPQ: Mapping quality (Phred-scaled)
    // 6. CIGAR: Extended CIGAR string (operations: MIDNSHP)
    // 7. MRNM: Mate reference name (‘=’ if same as RNAME)
    // 8. MPOS: 1-based leftmost mate position
    // 9. ISIZE: Inferred insert size
    // 10. SEQQuery: Sequence on the same strand as the reference
    // 11. QUAL: Query quality (ASCII-33=Phred base quality)    

    //
    // Some FASTQ files have spaces in their ID strings, which is illegal in SAM.  Just truncate them at the space.
    //
    const char *firstSpace = strnchr(read->getId(), ' ', qnameLen);
    if (NULL != firstSpace) {
        qnameLen = (unsigned)(firstSpace - read->getId());
    }

    const int nmStringSize = 30;// Big enough that it won't buffer overflow regardless of the value of editDistance
    char nmString[nmStringSize];
    snprintf(nmString, nmStringSize, "\tNM:i:%d", editDistance);

    unsigned auxLen;
    bool auxSAM;
    char* aux = read->getAuxiliaryData(&auxLen, &auxSAM);
    static bool warningPrinted = false;
    const char* readGroupSeparator = "";
    const char* readGroupString = "";
    if (aux != NULL && (!auxSAM)) {
        if (!warningPrinted) {
            WriteErrorMessage("warning: translating optional fields from BAM->SAM not yet implemented, optional fields will not be included in output\n");
            warningPrinted = true;
        }
        if (read->getReadGroup() == READ_GROUP_FROM_AUX) {
            for (BAMAlignAux* bamAux = (BAMAlignAux*)aux; (char*)bamAux < aux + auxLen; bamAux = bamAux->next()) {
                if (bamAux->tag[0] == 'R' && bamAux->tag[1] == 'G' && bamAux->val_type == 'Z') {
                    readGroupSeparator = "\tRG:Z:";
                    readGroupString = (char*)bamAux->value();
                    break;
                }
            }
        }
        aux = NULL;
        auxLen = 0;
    }
    const char* rglineAux = "";
    int rglineAuxLen = 0;
    if (read->getReadGroup() != NULL && read->getReadGroup() != READ_GROUP_FROM_AUX) {
        if (*readGroupString == 0 || strcmp(readGroupString, context.defaultReadGroup) == 0) {
            readGroupSeparator = "";
            readGroupString = "";
            rglineAux = context.defaultReadGroupAux;
            rglineAuxLen = context.defaultReadGroupAuxLen;
        }
        else {
            readGroupSeparator = "\tRG:Z:";
            readGroupString = read->getReadGroup();
        }
    }
    const int internalScoreBufferSize = 100;    // Should be plenty for \tXX:i:%d
    char internalScoreBuffer[internalScoreBufferSize];
    if (emitInternalScore) {
        int charsInInternalScore = snprintf(internalScoreBuffer, internalScoreBufferSize - 1, "\t%s:i:%d", internalScoreTag, (flags & SAM_UNMAPPED) ? -1 : internalScore);
        if (charsInInternalScore >= internalScoreBufferSize) {
            WriteErrorMessage("SAMFormat::writeRead overran internal buffer for internal score tag, which is kind of surprising.  %d\n", charsInInternalScore);
        }
    } else {
        internalScoreBuffer[0] = '\0';
    }

    const int alignmentTimeBufferSize = 100;    // Should be plenty for AT:i:123456789
    char alignmentTimeBuffer[alignmentTimeBufferSize];
    if (attachAlignmentTime) {
        _int64 alignmentTimeInMicroseconds = (alignmentTimeInNanoseconds + 500) / 1000;
        if (alignmentTimeInMicroseconds >= MAXINT32) {  // MAXINT is about 2 billion, so this would be ~35 minutes for one read (pair)
            alignmentTimeInMicroseconds = 0;
        };
        int charsInAligmentTime = snprintf(alignmentTimeBuffer, alignmentTimeBufferSize - 1, "\tAT:i:%d", (int)alignmentTimeInMicroseconds);
        if (charsInAligmentTime >= alignmentTimeBufferSize) {
            WriteErrorMessage("SAMFormat::writeRead overran internal buffer for alignment time tag, which is kind of surprising.  %d\n", charsInAligmentTime);
        }
    } else {
        alignmentTimeBuffer[0] = '\0';
    }

    int charsInString = snprintf(buffer, bufferSpace, "%.*s\t%d\t%s\t%llu\t%d\t%s\t%s\t%llu\t%lld\t%.*s\t%.*s%s%.*s%s%s\tPG:Z:SNAP%s%.*s%s%s%s%.*s\n",
        (unsigned)qnameLen, read->getId(),
        flags,
        contigName,
        positionInContig,
        mapQuality,
        cigar,
        matecontigName,
        matePositionInContig,
        templateLength,
        fullLength, data,
        fullLength, quality,
        aux != NULL ? "\t" : "", auxLen, aux != NULL ? aux : "",
        readGroupSeparator, readGroupString,
        nmString, rglineAuxLen, rglineAux,
        internalScoreBuffer,
        alignmentTimeBuffer,
        (FASTQCommentLength == 0) ? "" : "\t",
        FASTQCommentLength,
        FASTQComment);

    if (charsInString > bufferSpace) {
        //
        // Out of buffer space.
        //
        return false;
    }
    else if (charsInString == bufferSpace) {
        buffer[bufferSpace - 1] = '\n'; // overwrite trailing null with newline
    }


    if (NULL != spaceUsed) {
        *spaceUsed = charsInString;
    }
    return true;
}

//
// Common cigar string computation between SAM and BAM formats.
//
    void
SAMFormat::computeCigar(
    CigarFormat cigarFormat, 
    const Genome * genome, 
    LandauVishkinWithCigar * lv,
    char * cigarBuf, 
    int cigarBufLen, 
    const char * data, 
    GenomeDistance dataLength, 
    unsigned basesClippedBefore, 
    GenomeDistance extraBasesClippedBefore, 
    unsigned basesClippedAfter,
    GenomeDistance *o_extraBasesClippedAfter, 
    GenomeLocation genomeLocation, 
    bool useM, 
    int * o_editDistance, 
    int *o_cigarBufUsed, 
    int * o_addFrontClipping)
{
    if (dataLength > INT32_MAX - MAX_K) {
        dataLength = INT32_MAX - MAX_K;
    }

    int netIndel;
    *o_extraBasesClippedAfter = 0;

    //
    // Apply the extra clipping.
    //
    genomeLocation += extraBasesClippedBefore;
    data += extraBasesClippedBefore;
    dataLength -= extraBasesClippedBefore;

    const Genome::Contig *contig = genome->getContigAtLocation(genomeLocation);

    if (genomeLocation + dataLength > contig->beginningLocation + contig->length - genome->getChromosomePadding()) {
        //
        // The read hangs off the end of the contig.  Soft clip it at the end.  This is a tentative amount that assumes no net indels in the
        // mapping, we'll refine it later if needed.
        //
        *o_extraBasesClippedAfter = genomeLocation + dataLength - (contig->beginningLocation + contig->length - genome->getChromosomePadding());
    } else {
        *o_extraBasesClippedAfter = 0;
    }

    const char *reference = genome->getSubstring(genomeLocation, dataLength);
    if (NULL == reference) {
        //
        // Fell off the end of the contig.
        //
        *o_editDistance = 0;
        *o_addFrontClipping = 0;
        *o_cigarBufUsed = 0;
        *cigarBuf = '*';
        return;
    }

    *o_editDistance = lv->computeEditDistanceNormalized(
        reference,
        (int)(dataLength - *o_extraBasesClippedAfter + MAX_K), // Add space incase of indels.  We know there's enough, because the reference is padded.
        data,
        (int)(dataLength - *o_extraBasesClippedAfter),
        MAX_K - 1,
        cigarBuf,
        cigarBufLen,
        useM,
        cigarFormat,
        o_cigarBufUsed,
        o_addFrontClipping,
        &netIndel);

    if (*o_addFrontClipping != 0) {
        //
        // On this path, there really isn't a returned cigar string, it's sort of like an exception.  We're going up a level and
        // trying a different alignment.
        //
        return;
    }

    //
    // Normally, we'd be done.  However, if the amount that we would clip at the end of the read because of hanging off of the end
    // of the contig changed, then we need to recompute.  In some cases this is an iterative processess as we add or remove bits
    // of read.  
    //
    GenomeDistance newExtraBasesClippedAfter = __max(0, genomeLocation + dataLength + netIndel - (contig->beginningLocation + contig->length - genome->getChromosomePadding()));
    for (GenomeDistance pass = 0; pass < dataLength; pass++) {
        if (newExtraBasesClippedAfter == *o_extraBasesClippedAfter) {
            return;
        }

        *o_extraBasesClippedAfter = newExtraBasesClippedAfter;

        *o_editDistance = lv->computeEditDistanceNormalized(
            reference,
            (int)(dataLength - *o_extraBasesClippedAfter + MAX_K), // Add space incase of indels.  We know there's enough, because the reference is padded.
            data,
            (int)(dataLength - *o_extraBasesClippedAfter),
            MAX_K - 1,
            cigarBuf,
            cigarBufLen,
            useM,
            cigarFormat,
            o_cigarBufUsed,
            o_addFrontClipping,
            &netIndel);

        newExtraBasesClippedAfter = __max(0, genomeLocation + dataLength + netIndel - (contig->beginningLocation + contig->length - genome->getChromosomePadding()));
     }

     _ASSERT(!"cigar computation didn't converge");
    *o_extraBasesClippedAfter = newExtraBasesClippedAfter;
}

//
// Common cigar string computation between SAM and BAM formats.
//
    void
SAMFormat::computeCigar(
    CigarFormat cigarFormat,
    const Genome * genome,
    AffineGapVectorizedWithCigar * ag,
    char * cigarBuf,
    int cigarBufLen,
    const char * data,
    const char * quality,
    GenomeDistance dataLength,
    int score,
    unsigned basesClippedBefore,
    GenomeDistance extraBasesClippedBefore,
    unsigned basesClippedAfter,
    GenomeDistance *o_extraBasesClippedAfter,
    GenomeLocation genomeLocation,
    bool useM,
    int * o_editDistance,
    int *o_cigarBufUsed,
    int * o_addFrontClipping,
    int *o_backClippingMissedByLV)
{
    if (dataLength > INT32_MAX - MAX_K) {
        dataLength = INT32_MAX - MAX_K;
    }

    int netIndel = 0; // FIXME: Check how to use it with affine gap
    *o_extraBasesClippedAfter = 0;

    //
    // Apply the extra clipping.
    //
    genomeLocation += extraBasesClippedBefore;
    data += extraBasesClippedBefore;
    dataLength -= extraBasesClippedBefore;

    const Genome::Contig *contig = genome->getContigAtLocation(genomeLocation);

    if (genomeLocation + dataLength > contig->beginningLocation + contig->length - genome->getChromosomePadding()) {
        //
        // The read hangs off the end of the contig.  Soft clip it at the end.  This is a tentative amount that assumes no net indels in the
        // mapping, we'll refine it later if needed.
        //
        *o_extraBasesClippedAfter = genomeLocation + dataLength - (contig->beginningLocation + contig->length - genome->getChromosomePadding());
    } else {
        *o_extraBasesClippedAfter = 0;
    }

    const char *reference = genome->getSubstring(genomeLocation, dataLength);
    if (NULL == reference) {
        //
        // Fell off the end of the contig.
        //
        *o_editDistance = 0;
        *o_addFrontClipping = 0;
        *o_cigarBufUsed = 0;
        *cigarBuf = '*';
        return;
    }

    *o_editDistance = ag->computeGlobalScoreNormalized(
        reference,
        (int)(dataLength - *o_extraBasesClippedAfter + MAX_K), // Add space incase of indels.  We know there's enough, because the reference is padded.
        data,
        quality,
        (int)(dataLength - *o_extraBasesClippedAfter),
        score,
        cigarBuf,
        cigarBufLen,
        useM,
        cigarFormat,
        o_cigarBufUsed,
        o_addFrontClipping,
        &netIndel,
        o_backClippingMissedByLV);

    if (*o_addFrontClipping != 0) {
        //
        // On this path, there really isn't a returned cigar string, it's sort of like an exception.  We're going up a level and
        // trying a different alignment.
        //
        return;
    }

    //
    // Normally, we'd be done.  However, if the amount that we would clip at the end of the read because of hanging off of the end
    // of the contig changed, then we need to recompute.  In some cases this is an iterative processess as we add or remove bits
    // of read.  
    //
    GenomeDistance newExtraBasesClippedAfter = __max(0, genomeLocation + dataLength + netIndel - (contig->beginningLocation + contig->length - genome->getChromosomePadding()));
    for (GenomeDistance pass = 0; pass < dataLength; pass++) {
        if (newExtraBasesClippedAfter <= *o_extraBasesClippedAfter) {
            return;
        }

        *o_extraBasesClippedAfter = newExtraBasesClippedAfter;

        *o_editDistance = ag->computeGlobalScoreNormalized(
            reference,
            (int)(dataLength - *o_extraBasesClippedAfter + MAX_K), // Add space incase of indels.  We know there's enough, because the reference is padded.
            data,
            quality,
            (int)(dataLength - *o_extraBasesClippedAfter),
            score,
            cigarBuf,
            cigarBufLen,
            useM,
            cigarFormat,
            o_cigarBufUsed,
            o_addFrontClipping,
            &netIndel,
            o_backClippingMissedByLV);

        newExtraBasesClippedAfter = __max(0, genomeLocation + dataLength + netIndel - (contig->beginningLocation + contig->length - genome->getChromosomePadding()));
    }

	WriteErrorMessage("cigar computation didn't converge: data:%.*s\n", dataLength, data);
    _ASSERT(!"cigar computation didn't converge");
    *o_extraBasesClippedAfter = newExtraBasesClippedAfter;
}

// Compute the CIGAR edit sequence string for a read against a given genome location.
// Returns this string if possible or "*" if we fail to compute it (which would likely
// be a bug due to lack of buffer space). The pointer returned may be to cigarBuf so it
// will only be valid until computeCigarString is called again.
    const char *
SAMFormat::computeCigarString(
    const Genome *              genome,
    LandauVishkinWithCigar *    lv,
    char *                      cigarBuf,
    int                         cigarBufLen,
    char *                      cigarBufWithClipping,
    int                         cigarBufWithClippingLen,
    const char *                data,
    GenomeDistance              dataLength,
    unsigned                    basesClippedBefore,
    GenomeDistance              extraBasesClippedBefore,
    unsigned                    basesClippedAfter,
    unsigned                    frontHardClipping,
    unsigned                    backHardClipping,
    GenomeLocation              genomeLocation,
    Direction                   direction,
	bool						useM,
    int *                       o_editDistance,
    int *                       o_addFrontClipping,
    int *                       o_refSpan
)
{
    GenomeDistance extraBasesClippedAfter;
    int cigarBufUsed;

    computeCigar(COMPACT_CIGAR_STRING, genome, lv, cigarBuf, cigarBufLen, data, dataLength, basesClippedBefore,
        extraBasesClippedBefore, basesClippedAfter, &extraBasesClippedAfter, genomeLocation, useM,
        o_editDistance, &cigarBufUsed, o_addFrontClipping);

    if (*o_addFrontClipping != 0) {
        return NULL;
    }

    if (*o_editDistance == -2) {
        WriteErrorMessage( "WARNING: computeEditDistance returned -2; cigarBuf may be too small\n");
        return "*";
    } else if (*o_editDistance == -1) {
        static bool warningPrinted = false;
        if (!warningPrinted) {
            WriteErrorMessage( "WARNING: computeEditDistance returned -1; this shouldn't happen. Read %.*s\n", dataLength, data);
            warningPrinted = true;
        }
        return "*";
    } else {
        // Add some CIGAR instructions for soft-clipping if we've ignored some bases in the read.
        char clipBefore[16] = {'\0'};
        char clipAfter[16] = {'\0'};
        char hardClipBefore[16] = {'\0'};
        char hardClipAfter[16] = {'\0'};

        if (frontHardClipping > 0) {
            snprintf(hardClipBefore, sizeof(hardClipBefore), "%uH", frontHardClipping);
        }

        if (basesClippedBefore + extraBasesClippedBefore > 0) {
            snprintf(clipBefore, sizeof(clipBefore), "%lluS", basesClippedBefore + extraBasesClippedBefore);
        }

        if (basesClippedAfter + extraBasesClippedAfter > 0) {
            snprintf(clipAfter, sizeof(clipAfter), "%lluS", basesClippedAfter + extraBasesClippedAfter);
        }

        if (backHardClipping > 0) {
            snprintf(hardClipAfter, sizeof(hardClipAfter), "%uH", backHardClipping);
        }

        snprintf(cigarBufWithClipping, cigarBufWithClippingLen, "%s%s%s%s%s", hardClipBefore, clipBefore, cigarBuf, clipAfter, hardClipAfter);

		validateCigarString(genome, cigarBufWithClipping, cigarBufWithClippingLen, 
			data - basesClippedBefore, dataLength + ((size_t)basesClippedBefore + (size_t)basesClippedAfter), genomeLocation + extraBasesClippedBefore, direction, useM);

        *o_refSpan = 0;
        getRefSpanFromCigar(cigarBufWithClipping, cigarBufWithClippingLen, o_refSpan);

        return cigarBufWithClipping;
    }
}

// Compute the CIGAR edit sequence string for a read against a given genome location.
// Returns this string if possible or "*" if we fail to compute it (which would likely
// be a bug due to lack of buffer space). The pointer returned may be to cigarBuf so it
// will only be valid until computeCigarString is called again.
    const char *
SAMFormat::computeCigarString(
    const Genome *              genome,
    AffineGapVectorizedWithCigar *ag,
    char *                      cigarBuf,
    int                         cigarBufLen,
    char *                      cigarBufWithClipping,
    int                         cigarBufWithClippingLen,
    const char *                data,
    const char *                quality,
    GenomeDistance              dataLength,
    int                         score,
    unsigned                    basesClippedBefore,
    GenomeDistance              extraBasesClippedBefore,
    unsigned                    basesClippedAfter,
    unsigned                    frontHardClipping,
    unsigned                    backHardClipping,
    GenomeLocation              genomeLocation,
    Direction                   direction,
    bool						useM,
    int *                       o_editDistance,
    int *                       o_addFrontClipping,
    int *                       o_refSpan
)
{
    GenomeDistance extraBasesClippedAfter;
    int cigarBufUsed;
    int backClippingMissedByLV = 0;

    computeCigar(COMPACT_CIGAR_STRING, genome, ag, cigarBuf, cigarBufLen, data, quality, dataLength, score, basesClippedBefore,
        extraBasesClippedBefore, basesClippedAfter, &extraBasesClippedAfter, genomeLocation, useM,
        o_editDistance, &cigarBufUsed, o_addFrontClipping, &backClippingMissedByLV);

    if (*o_addFrontClipping != 0) {
        return NULL;
    }

    if (*o_editDistance == -2) {
        WriteErrorMessage("WARNING: computeGlobalScore returned -2; cigarBuf may be too small\n");
        return "*";
    } else if (*o_editDistance == -1) {
        static bool warningPrinted = false;
        if (!warningPrinted) {
            WriteErrorMessage("WARNING: computeGlobalScore returned -1; this shouldn't happen\n");
            warningPrinted = true;
        }
        return "*";
    } else {
        // There may be a better way to do this. For now, whenever we see tail insertions, soft-clip them
        basesClippedAfter += backClippingMissedByLV;
        dataLength -= backClippingMissedByLV;
        
        // if (dataLength > 101) {
        //     WriteErrorMessage("computeCigarString: Data length negative for read:%.*s, isRC:%d\n", 101, data - basesClippedBefore, direction);
        //     soft_exit(1);
        // }

        // Add some CIGAR instructions for soft-clipping if we've ignored some bases in the read.
        char clipBefore[16] = { '\0' };
        char clipAfter[16] = { '\0' };
        char hardClipBefore[16] = { '\0' };
        char hardClipAfter[16] = { '\0' };

        if (frontHardClipping > 0) {
            snprintf(hardClipBefore, sizeof(hardClipBefore), "%uH", frontHardClipping);
        }

        if (basesClippedBefore + extraBasesClippedBefore > 0) {
            snprintf(clipBefore, sizeof(clipBefore), "%lluS", basesClippedBefore + extraBasesClippedBefore);
        }

        if (basesClippedAfter + extraBasesClippedAfter > 0) {
            snprintf(clipAfter, sizeof(clipAfter), "%lluS", basesClippedAfter + extraBasesClippedAfter);
        }

        if (backHardClipping > 0) {
            snprintf(hardClipAfter, sizeof(hardClipAfter), "%uH", backHardClipping);
        }

        snprintf(cigarBufWithClipping, cigarBufWithClippingLen, "%s%s%s%s%s", hardClipBefore, clipBefore, cigarBuf, clipAfter, hardClipAfter);

        validateCigarString(genome, cigarBufWithClipping, cigarBufWithClippingLen,
            data - basesClippedBefore, dataLength + (basesClippedBefore + basesClippedAfter), genomeLocation + extraBasesClippedBefore, direction, useM);

        *o_refSpan = 0;
        getRefSpanFromCigar(cigarBufWithClipping, cigarBufWithClippingLen, o_refSpan);

        return cigarBufWithClipping;
    }
}

    void 
SAMFormat::getRefSpanFromCigar(const char * cigarBuf, int cigarBufLen, int* refSpan) 
{
    if (cigarBufLen == 0) return;
    const char *nextChunkOfCigar = cigarBuf;
    unsigned len;
    char op;
    int fieldsScanned = sscanf(nextChunkOfCigar, "%d%c", &len, &op);
    if (op != 'S' && op != 'H') {
        *refSpan += len;
    }
    //
    // Now scan over the current op.
    //
    while ('0' <= *nextChunkOfCigar && '9' >= *nextChunkOfCigar) {
        nextChunkOfCigar++;
    }

    nextChunkOfCigar++;
    while ('\0' != *nextChunkOfCigar) {
        fieldsScanned = sscanf(nextChunkOfCigar, "%d%c", &len, &op);
        if (op != 'I') {
            *refSpan += len;
        }
        //
        // Now scan over the current op.
        //
        while ('0' <= *nextChunkOfCigar && '9' >= *nextChunkOfCigar) {
            nextChunkOfCigar++;
        }
        nextChunkOfCigar++;
    }
}

    void
SAMFormat::printRead(const Genome* genome, Read* read, GenomeLocation location)
{
    if (read) {
        const char* read_data = read->getUnclippedData();
        const char* readId = read->getId();

        for (unsigned i = 0; i < read->getIdLength(); ++i) {
            printf("%c", readId[i]);
        }

        printf(",");

        for (unsigned i = 0; i < read->getUnclippedLength(); ++i) {
            printf("%c", read_data[i]);
        }

        printf("\n");
        printf("Aligned location %s:%llu\n",
            genome->getContigAtLocation(location)->name,
            location - genome->getContigAtLocation(location)->beginningLocation);
    }
}

// #ifdef _DEBUG
	void 
SAMFormat::validateCigarString(
	const Genome *genome, const char * cigarBuf, int cigarBufLen, const char *data, GenomeDistance dataLength, GenomeLocation genomeLocation, Direction direction, bool useM, Read* read)
{
	const char *nextChunkOfCigar = cigarBuf;
	GenomeDistance offsetInData = 0;
	const char *reference = genome->getSubstring(genomeLocation, dataLength);

	if (NULL == reference) {
		WriteErrorMessage("validateCigarString: couldn't look up genome data for location %lld, data length %lld, read ID %.*s cigar string %.*s\n", GenomeLocationAsInt64(genomeLocation), dataLength,
                            read->getIdLength(), read->getId(), cigarBufLen, cigarBuf);
        printRead(genome, read, genomeLocation);
		soft_exit(1);
	}

	GenomeDistance offsetInReference = 0;
	bool sawNonH = false;	// This is to make sure that the clipping types (H & S) occur only at the beginning or end of the cigar string.
	bool sawTailS = false;	// Did we see a S
	bool sawLeadingS = false;	// Have we seen the soft clip at the front of the cigar string?
	bool sawTrailingH = false;
	char previousOp = '\0';	// Make sure that we don't have two consecutive ops of the same type that should be merged
	bool sawXorM = false;
	bool lastItemWasIndel = false;

	//
	// First check to see that it's null terminated
	//
	bool nullTerminated = false;
	for (size_t offset = 0; offset < cigarBufLen; offset++) {
		if ('\0' == cigarBuf[offset]) {
			nullTerminated = true;
			break;
		}
	}

	if (!nullTerminated) {
		WriteErrorMessage("validateCigarString: non-null-terminated or overflow cigar string: '%.*s'\n", cigarBufLen, cigarBuf);
        printRead(genome, read, genomeLocation);
        soft_exit(1);
	}

	const Genome::Contig *contig = genome->getContigAtLocation(genomeLocation);
	if (NULL == contig) {
		WriteErrorMessage("validateCigarString: read alignment location isn't in a chromosome, genomeLocation %lld\n", GenomeLocationAsInt64(genomeLocation));
        printRead(genome, read, genomeLocation);
        soft_exit(1);
	}

	if (genomeLocation >= contig->beginningLocation + contig->length - genome->getChromosomePadding()) {
		WriteErrorMessage("validateCigarString: alignment location is in genome padding: %lld, contig name %s, base %lld, len %lld, padding size %d\n",
			GenomeLocationAsInt64(genomeLocation), contig->name, GenomeLocationAsInt64(contig->beginningLocation), contig->length, genome->getChromosomePadding());
        printRead(genome, read, genomeLocation);
        soft_exit(1);
	}

	while ('\0' != *nextChunkOfCigar) {
		unsigned len;
		char op;
		int fieldsScanned = sscanf(nextChunkOfCigar, "%d%c", &len, &op);
		if (2 != fieldsScanned) {
			WriteErrorMessage("validateCigarString: didn't scan two fields here '%s' in overall cigar string '%s'\n", nextChunkOfCigar, cigarBuf);
            printRead(genome, read, genomeLocation);
            soft_exit(1);
		}

		if (0 == len) {
			WriteErrorMessage("validateCigarString: got zero length field here '%s' in overall cigar string '%s'\n", nextChunkOfCigar, cigarBuf);
            printRead(genome, read, genomeLocation);
            soft_exit(1);
		}

		if (op != 'H' && sawTailS) {
			WriteErrorMessage("validateCigarString: saw incorrect op type after what should have been the terminal soft or hard clipping here '%s', in overall cigar string '%s'\n",
				nextChunkOfCigar, cigarBuf);
            printRead(genome, read, genomeLocation);
            soft_exit(1);
		}

		if (sawTrailingH) {
			WriteErrorMessage("validateCigarString: saw op after what should have been the terminal hard clip here '%s' in overall cigar '%s'\n", nextChunkOfCigar, cigarBuf);
            printRead(genome, read, genomeLocation);
            soft_exit(1);
		}

		if (op == previousOp) {
			WriteErrorMessage("validateCigarString: saw consecutive ops of the same type '%c' here '%s' in overall cigar '%s'\n", op, nextChunkOfCigar, cigarBuf);
            printRead(genome, read, genomeLocation);
            soft_exit(1);
		}

		switch (op) {
			case 'M': 
			{
				if (!useM) {
					WriteErrorMessage("validateCigarString: generated an M when we were supposed to use X and = here '%s' in overall cigar string '%s'\n", nextChunkOfCigar, cigarBuf);
                    printRead(genome, read, genomeLocation);
                    soft_exit(1);
				}
				offsetInData += len;
				sawNonH = true;
				sawXorM = true;
				lastItemWasIndel = false;
				break;
			}

			case 'X': 
			case '=': 
			{
				if (useM) {
					WriteErrorMessage("validateCigarString: generated an %c when were supposed to use M here '%s' in overall cigar string '%s'\n", op, nextChunkOfCigar, cigarBuf);
                    printRead(genome, read, genomeLocation);
                    soft_exit(1);
				}

				if (len + offsetInData > dataLength) {
					WriteErrorMessage("validateCigarString: cigar string overflowed read length, here '%s', overall cigar '%s'\n", nextChunkOfCigar, cigarBuf);
                    printRead(genome, read, genomeLocation);
                    soft_exit(1);
				}

				for (unsigned offset = 0; offset < len; offset++) {
					if ((data[offset + offsetInData] == reference[offset + offsetInReference]) == ('X' == op)) {
						WriteErrorMessage("validateCigarString: saw a (non-)matching base in an %c range, offset %d, offsetInData %lld, offsetInReference %lld, data '%.*s', reference '%.*s', here '%s', overall cigar '%s'\n",
							op, offset, offsetInData, offsetInReference, dataLength, data, dataLength, reference, nextChunkOfCigar, cigarBuf);
                        printRead(genome, read, genomeLocation);
                        soft_exit(1);
					}
				}

				offsetInData += len;
				offsetInReference += len;
				sawNonH = true;
				sawXorM = true;
				lastItemWasIndel = false;
				break;
			}

			case 'I': 
			{
				//
				// Insertion uses up bases in the read but not in the reference.
				//
				if (len + offsetInData > dataLength) {
					WriteErrorMessage("validateCigarString: insertion pushes cigar string overlength, here '%s' in overall cigar '%s'\n", nextChunkOfCigar, cigarBuf);
                    printRead(genome, read, genomeLocation);
                    soft_exit(1);
				}

				if (!sawXorM) {
					WriteErrorMessage("validateCigarString: cigar string started with I (after clipping) here '%s' in overall cigar '%s'\n", nextChunkOfCigar, cigarBuf);
                    printRead(genome, read, genomeLocation);
                    soft_exit(1);
				}

                if (previousOp == 'D') {
                    WriteErrorMessage("validateCigarString: cigar string had D immediately followed by I here '%'s in overall cigar '%s'\n", nextChunkOfCigar, cigarBuf);
                    printRead(genome, read, genomeLocation);
                    soft_exit(1);
                }

				offsetInData += len;
				sawNonH = true;
				lastItemWasIndel = true;
				break;
			}

			case 'D':
			{
				if (!sawXorM) {
					WriteErrorMessage("validateCigarString: cigar string started with D (after clipping) here '%s' in overall cigar '%s'\n", nextChunkOfCigar, cigarBuf);
                    printRead(genome, read, genomeLocation);
                    soft_exit(1);
				}
						
                if (previousOp == 'I') {
                    WriteErrorMessage("validateCigarString: cigar string had I immediately followed by D here '%'s in overall cigar '%s'\n", nextChunkOfCigar, cigarBuf);
                    printRead(genome, read, genomeLocation);
                    soft_exit(1);
                }

                //
				// D uses up bases in the reference but not the read.
				//
				offsetInReference += len;
				sawNonH = true;
				lastItemWasIndel = true;
				break;
			}

			case 'N':
			case 'P':
			{
				WriteErrorMessage("validateCigarString: saw valid op type '%c' that SNAP shouldn't generate, here '%s' in overall cigar string '%s'\n", op, nextChunkOfCigar, cigarBuf);
                printRead(genome, read, genomeLocation);
                soft_exit(1);
			}

			case 'H':
			{
				//
				// Hard clip bases do not occur in the read string at all.  All we can validate is that this is the first or last thing in the cigar string.
				//
				if (nextChunkOfCigar == cigarBuf) {
					//
					// First thing, this is OK.
					//
					break;
				}
				sawTrailingH = true;
				break;
			}

			case 'S':
			{
				if (sawNonH) {
					sawTailS = true;
				}
                sawNonH = true;
 				offsetInData += len;
				break;
			}


			default: {
				WriteErrorMessage("validateCigarString: got unrecognized cigar op '%c', here '%s' in overall string '%s'\n", op, nextChunkOfCigar, cigarBuf);
                printRead(genome, read, genomeLocation);
                soft_exit(1);
			}
		}

		previousOp = op;
		//
		// Now scan over the current op.
		//
		while ('0' <= *nextChunkOfCigar && '9' >= *nextChunkOfCigar) {
			nextChunkOfCigar++;
		}
		if (*nextChunkOfCigar != op) {
			WriteErrorMessage("validateCigarString: bug in validation code; expected op '%c', got '%c' at '%s' in '%s'\n", op, *nextChunkOfCigar, nextChunkOfCigar, cigarBuf);
            printRead(genome, read, genomeLocation);
            soft_exit(1);
		}
		nextChunkOfCigar++;
	}

	if (offsetInData != dataLength) {
		WriteErrorMessage("validateCigarString: Didn't consume entire read data, got %lld of %lld, cigar '%s'\n", offsetInData, dataLength, cigarBuf);
        printRead(genome, read, genomeLocation);
        soft_exit(1);
	}

	if (lastItemWasIndel) {
		WriteErrorMessage("validateCigarString: cigar string ended with indel '%s'\n", cigarBuf);
        printRead(genome, read, genomeLocation);
        soft_exit(1);
	}

    //
    // Make sure none of the non-soft-clipped part of the read is mapped onto padding.
    //
    if (genomeLocation + offsetInReference > contig->beginningLocation + contig->length - genome->getChromosomePadding()) {
        WriteErrorMessage("validateCigarString: alignment runs into contig padding: %lld, contig name %s, base %lld, len %lld, padding size %d, offsetInReference %lld\n",
            GenomeLocationAsInt64(genomeLocation), contig->name, GenomeLocationAsInt64(contig->beginningLocation), contig->length, genome->getChromosomePadding(), offsetInReference);
        printRead(genome, read, genomeLocation);
        soft_exit(1);
    }
}

// #endif // _DEBUG

struct OffsetInfo {
    size_t offset;
    bool isDuplicate;

    OffsetInfo() : offset(0), isDuplicate(false) {}
    OffsetInfo(size_t offset_, bool isDuplicate_) {
        offset = offset_;
        isDuplicate = isDuplicate_;
    }
    ~OffsetInfo() {}
};

class SAMFilter : public DataWriter::Filter
{
public:
    SAMFilter(DataWriter::FilterType i_type) : Filter(i_type), offsets(1000), header(false) {}

    virtual ~SAMFilter() {
        offsets.clear();
    }

    virtual void inHeader(bool flag)
    {
        header = flag;
    }

    virtual void onAdvance(DataWriter* writer, size_t batchOffset, char* data, GenomeDistance bytes, GenomeLocation location);

    virtual size_t onNextBatch(DataWriter* writer, size_t offset, size_t bytes, bool lastBatch = false, bool* needMoreBuffer = NULL, size_t* fromBufferUsed = NULL) = 0;

protected:

    void getAlignmentInfo(const Genome* genome, char* buffer, _int64 bytes, SAMAlignment* sam);

    void updateSAMLine(char* toBuffer, size_t* toUsed, char* fromBuffer, size_t bytes);

    static int getTotalQuality(char* qual, int l_seq);

    static void getTileXY(const char* id, int* o_tile, int* o_x, int* o_y);

    VariableSizeVector<OffsetInfo> offsets;
    DataWriter* currentWriter;
    char* currentBuffer;
    size_t currentBufferBytes; // # of valid bytes
    size_t currentOffset; // logical file offset of beginning of current buffer

private:
    bool header;
};

struct SAMDupMarkEntry
{
    SAMDupMarkEntry() : libraryNameHash(0), runOffset(0), flag(0), qual(0), mateQual(0), mateInfo(0), mateLocation(InvalidGenomeLocation), info(0), tlen(0), qName(NULL), qNameLength(0) {}

    SAMDupMarkEntry(size_t libraryNameHash_, size_t runOffset_, int flag_, _int32 qual_, _int32 mateQual_, _uint64 mateInfo_, GenomeLocation mateLocation_, _uint64 info_, _int64 tlen_, char* qName_, size_t qNameLength_)
        : libraryNameHash(libraryNameHash_),
          runOffset(runOffset_),
          flag(flag_),
          qual(qual_),
          mateQual(mateQual_),
          mateInfo(mateInfo_),
          mateLocation(mateLocation_),
          info(info_),
          tlen(tlen_),
          qName(qName_),
          qNameLength(qNameLength_)
    {
    }

    bool operator<(const SAMDupMarkEntry& b) const {
        return (libraryNameHash < b.libraryNameHash) ||
            ((libraryNameHash == b.libraryNameHash) && (info < b.info)) ||
            ((libraryNameHash == b.libraryNameHash && info == b.info && mateInfo < b.mateInfo));
    }

    static size_t hash(char* str) // We could probably use a better hash function, but I doubt it matters
    {
        size_t value = 0x123456789abcdef0;
        char* libraryName = str;

        for (int i = 0; i < strlen(str); i++)
        {
            value = (value * 131) ^ libraryName[i];
        }

        return value;
    } // hash

    size_t libraryNameHash; // hash value of library name string LB:Z:xxxx
    size_t runOffset; // read offset in run
    int flag; // SAM flag
    _int32 qual; // read quality score
    _int32 mateQual; // mate quality score (only considering bases with phred score >= 15)
    _uint64 mateInfo; // mate information has: mateLocation (calculated from 32-bit signed TLEN) and matedirection
    GenomeLocation mateLocation; // actual mateLocation calculated from RNEXT and PNEXT
    _uint64 info; // read information has: location and direction
    _int64 tlen;
    char* qName;
    size_t qNameLength;
};

    void
SAMFilter::getAlignmentInfo(
    const Genome* genome, 
    char* buffer, 
    _int64 bytes,
    SAMAlignment* sam)
{
    char* fields[SAMReader::nSAMFields];
    size_t lengths[SAMReader::nSAMFields];
    size_t lineLength;
    if (!SAMReader::parseLine(buffer, buffer + bytes, fields, &lineLength, lengths)) {
        WriteErrorMessage("Failed to parse SAM line:\n%.*s\n", lineLength, buffer);
        soft_exit(1);
    }
    _ASSERT(lineLength < UINT32_MAX);

    //
    // We have to copy the contig name (RNAME) into its own buffer because the code in Genome expects
    // it to be a null-terminated string, while all we've got is one that's space delimited.
    //
    const size_t contigNameBufferSize = 512;
    char contigNameBuffer[contigNameBufferSize];
    char *contigName = contigNameBuffer;
    size_t neededSize;
    GenomeLocation locationOfContig;
    if (0 != (neededSize = SAMReader::parseContigName(genome, contigName, contigNameBufferSize, &locationOfContig, NULL, fields, lengths))) {
        contigName = new char[neededSize];
        if (0 != SAMReader::parseContigName(genome, contigName, neededSize, &locationOfContig, NULL, fields, lengths)) {
            WriteErrorMessage("SAMFilter::getAlignmentInfo: reallocated contigName was still too small\n");
            soft_exit(1);
        }
    }

    const size_t mateContigNameBufferSize = 512;
    char mateContigNameBuffer[mateContigNameBufferSize];
    char* mateContigName = mateContigNameBuffer;
    GenomeLocation mateLocationOfContig;
    if (0 != (neededSize = SAMReader::parseContigName(genome, mateContigName, mateContigNameBufferSize, &mateLocationOfContig, NULL, fields, lengths, 6))) {
        mateContigName = new char[neededSize];
        if (0 != SAMReader::parseContigName(genome, mateContigName, neededSize, &mateLocationOfContig, NULL, fields, lengths, 6)) {
            WriteErrorMessage("SAMFilter::getAlignmentInfo: reallocated mateContigName was still too small\n");
            soft_exit(1);
        }
    }

    GenomeLocation genomeLocation, mateGenomeLocation;
    
    if (InvalidGenomeLocation != locationOfContig) {
        genomeLocation = SAMReader::parseLocation(locationOfContig, fields, lengths);
    } else {
        genomeLocation = InvalidGenomeLocation;
    }

    if (InvalidGenomeLocation != mateLocationOfContig) {
        mateGenomeLocation = SAMReader::parseLocation(mateLocationOfContig, fields, lengths, 6, 7);
    }
    else {
        mateGenomeLocation = InvalidGenomeLocation;
    }

    // FLAG
    int _flag;
    const size_t flagBufferSize = 20;   // More than enough
    char flagBuffer[flagBufferSize];
    if (lengths[SAMReader::FLAG] >= flagBufferSize) {
        WriteErrorMessage("SAMFilter::getAlignmentInfo flag field is too long. %.*s\n", lengths[SAMReader::FLAG], fields[SAMReader::FLAG]);
        soft_exit(1);
    }
    memcpy(flagBuffer, fields[SAMReader::FLAG], lengths[SAMReader::FLAG]);
    flagBuffer[lengths[SAMReader::FLAG]] = '\0';
    if (1 != sscanf(flagBuffer, "%d", &_flag)) {
        WriteErrorMessage("SAMFilter::getAlignmentInfo couldn't parse FLAG field.\n");
        soft_exit(1);
    }

    // TLEN
    _int64 tlen = atoll(fields[SAMReader::TLEN]);

    int mateQual = 0;
    size_t libraryNameHash = 0;
    if (fields[SAMReader::OPT] != NULL) {
        unsigned n = (unsigned)lengths[SAMReader::OPT];
        while (n > 0 && (fields[SAMReader::OPT][n - 1] == '\n' || fields[SAMReader::OPT][n - 1] == '\r')) {
            n--;
        }
        for (char* p = fields[SAMReader::OPT]; p != NULL && p < fields[SAMReader::OPT] + lengths[SAMReader::OPT]; p = SAMReader::skipToBeyondNextFieldSeparator(p, fields[SAMReader::OPT] + lengths[SAMReader::OPT])) {
            if (strncmp(p, "QS:i:", 5) == 0) {
                const size_t mateQualBufferSize = 20; // should be more than enough
                char mateQualBuffer[mateQualBufferSize];
                char* mateQualValueStart = p + 5;
                p = SAMReader::skipToBeyondNextFieldSeparator(p, fields[SAMReader::OPT] + lengths[SAMReader::OPT]);
                size_t mateQualValueLen = (p != NULL) ? (p - 1 - mateQualValueStart) : (fields[SAMReader::OPT] + lengths[SAMReader::OPT] - mateQualValueStart); // fields separator at p - 1
                if (mateQualValueLen >= mateQualBufferSize) {
                    WriteErrorMessage("SAMFilter::getAlignmentInfo QS field is too long.\n");
                    soft_exit(1);
                }
                memcpy(mateQualBuffer, mateQualValueStart, mateQualValueLen);
                mateQualBuffer[mateQualValueLen] = '\0';
                if (1 != sscanf(mateQualBuffer, "%d", &mateQual)) {
                    WriteErrorMessage("SAMFilter::getAlignmentInfo couldn't parse QS field.\n");
                    soft_exit(1);
                }
            }
            if (p != NULL && p < fields[SAMReader::OPT] + lengths[SAMReader::OPT]) {
                if (strncmp(p, "LB:Z:", 5) == 0) {
                    const size_t libraryNameBufferSize = 512; // should be more than enough
                    char libraryNameBuffer[libraryNameBufferSize];
                    char* libraryNameStart = p + 5;
                    p = SAMReader::skipToBeyondNextFieldSeparator(p, fields[SAMReader::OPT] + lengths[SAMReader::OPT]);
                    size_t libraryNameLen = (p != NULL) ? (p - 1 - libraryNameStart) : (fields[SAMReader::OPT] + lengths[SAMReader::OPT] - libraryNameStart); // fields separator at p - 1
                    if (libraryNameLen >= libraryNameBufferSize) {
                        WriteErrorMessage("SAMFilter::getAlignmentInfo LB field is too long.\n");
                        soft_exit(1);
                    }
                    memcpy(libraryNameBuffer, libraryNameStart, libraryNameLen);
                    libraryNameBuffer[libraryNameLen] = '\0';
                    libraryNameHash = SAMDupMarkEntry::hash(libraryNameBuffer);
                }
            }
        }
    }

    if (sam != NULL) {
        sam->flag = _flag;
        sam->tlen = tlen;
        sam->location = genomeLocation;
        sam->mateLocation = mateGenomeLocation;
        sam->unclippedStartLocation = sam->getUnclippedStart(genomeLocation, fields[SAMReader::CIGAR], (int)lengths[SAMReader::CIGAR]);
        sam->unclippedEndLocation = sam->getUnclippedEnd(genomeLocation, fields[SAMReader::CIGAR], (int)lengths[SAMReader::CIGAR]);
        sam->setReadName(fields[SAMReader::QNAME], lengths[SAMReader::QNAME]);
        sam->qual = getTotalQuality(fields[SAMReader::QUAL], (int)lengths[SAMReader::QUAL]);
        sam->mateQual = mateQual;
        sam->libraryNameHash = libraryNameHash;
        sam->flagOffset = fields[SAMReader::FLAG];
    }

    if (contigName != contigNameBuffer) {
        delete[] contigName;
    }

    if (mateContigName != mateContigNameBuffer) {
        delete[] mateContigName;
    }
}

    void 
SAMFilter::updateSAMLine(
    char* toBuffer,
    size_t* toUsed,
    char* fromBuffer,
    size_t bytes)
{
    char* fields[SAMReader::nSAMFields];
    size_t lengths[SAMReader::nSAMFields];
    size_t lineLength;

    if (!SAMReader::parseLine(fromBuffer, fromBuffer + bytes, fields, &lineLength, lengths)) {
        WriteErrorMessage("Failed to parse SAM line:\n%.*s\n", lineLength, fromBuffer);
        soft_exit(1);
    }
    _ASSERT(lineLength < UINT32_MAX);
    
    int flag;
    const size_t flagBufferSize = 20;   // More than enough
    char flagBuffer[flagBufferSize];
    if (lengths[SAMReader::FLAG] >= flagBufferSize) {
        WriteErrorMessage("SAMFilter::updateSAMLine flag field is too long.\n");
        soft_exit(1);
    }
    memcpy(flagBuffer, fields[SAMReader::FLAG], lengths[SAMReader::FLAG]);
    flagBuffer[lengths[SAMReader::FLAG]] = '\0';
    if (1 != sscanf(flagBuffer, "%d", &flag)) {
        WriteErrorMessage("SAMFilter::updateSAMLine couldn't parse FLAG field.\n");
        soft_exit(1);
    }

    flag |= SAM_DUPLICATE;

    char* prev = fields[SAMReader::OPT];
    //
    // We remove either the QS: or LB: fields at the end of the SAM record to make space for the new flag. Otherwise we will run out of space in buffer.
    // XXX: This is just awful and wrong.  We must fix this.
    //
    size_t lengthOPTExcludingLastField = 0;
    bool seenPG = false;
    for (char* p = fields[SAMReader::OPT]; p != NULL && p < fields[SAMReader::OPT] + lengths[SAMReader::OPT]; p = SAMReader::skipToBeyondNextFieldSeparator(p, fields[SAMReader::OPT] + lengths[SAMReader::OPT])) {
        lengthOPTExcludingLastField += (p - prev);
        if (!strncmp(p, "PG:Z:SNAP", 9)) {
            seenPG = true;
        }
        prev = p;
    }

    //
    // Trim a trailing \t if it's there.
    //
    if (lengthOPTExcludingLastField > 0 && fields[SAMReader::OPT][lengthOPTExcludingLastField - 1] == '\t') {
        lengthOPTExcludingLastField--;
    }

    lengths[SAMReader::OPT] = lengthOPTExcludingLastField;



    // FIXME: needs to be changed if more mandatory SAM fields are added
    int charsInString;
    if (seenPG) {
        charsInString = snprintf(toBuffer + *toUsed, lineLength, "%.*s\t%d\t%.*s\t%.*s\t%.*s\t%.*s\t%.*s\t%.*s\t%.*s\t%.*s\t%.*s\t%.*s\n",
            (int)lengths[SAMReader::QNAME], fields[SAMReader::QNAME],
            flag,
            (int)lengths[SAMReader::RNAME], fields[SAMReader::RNAME],
            (int)lengths[SAMReader::POS], fields[SAMReader::POS],
            (int)lengths[SAMReader::MAPQ], fields[SAMReader::MAPQ],
            (int)lengths[SAMReader::CIGAR], fields[SAMReader::CIGAR],
            (int)lengths[SAMReader::RNEXT], fields[SAMReader::RNEXT],
            (int)lengths[SAMReader::PNEXT], fields[SAMReader::PNEXT],
            (int)lengths[SAMReader::TLEN], fields[SAMReader::TLEN],
            (int)lengths[SAMReader::SEQ], fields[SAMReader::SEQ],
            (int)lengths[SAMReader::QUAL], fields[SAMReader::QUAL],
            (int)lengths[SAMReader::OPT], fields[SAMReader::OPT]);
    } else {
        charsInString = snprintf(toBuffer + *toUsed, lineLength, "%.*s\t%d\t%.*s\t%.*s\t%.*s\t%.*s\t%.*s\t%.*s\t%.*s\t%.*s\t%.*s\tPG:Z:SNAP\t%.*s\n",
            (int)lengths[SAMReader::QNAME], fields[SAMReader::QNAME],
            flag,
            (int)lengths[SAMReader::RNAME], fields[SAMReader::RNAME],
            (int)lengths[SAMReader::POS], fields[SAMReader::POS],
            (int)lengths[SAMReader::MAPQ], fields[SAMReader::MAPQ],
            (int)lengths[SAMReader::CIGAR], fields[SAMReader::CIGAR],
            (int)lengths[SAMReader::RNEXT], fields[SAMReader::RNEXT],
            (int)lengths[SAMReader::PNEXT], fields[SAMReader::PNEXT],
            (int)lengths[SAMReader::TLEN], fields[SAMReader::TLEN],
            (int)lengths[SAMReader::SEQ], fields[SAMReader::SEQ],
            (int)lengths[SAMReader::QUAL], fields[SAMReader::QUAL],
            (int)lengths[SAMReader::OPT], fields[SAMReader::OPT]);
    }

    *toUsed += charsInString;

    if (charsInString > bytes || charsInString > lineLength) {
        WriteErrorMessage("SAMFilter::updateSAMLine: charsInString %d is larger than bytes %lld, linelength %lld\n", charsInString, bytes, lineLength);
        soft_exit(1);
    }
}


    void
SAMFilter::onAdvance(
    DataWriter* writer,
    size_t batchOffset,
    char* data,
    GenomeDistance bytes,
    GenomeLocation location)
{
    if (! header) {
        OffsetInfo oinfo(batchOffset, false);
        offsets.push_back(oinfo);
    }
}

    int
SAMFilter::getTotalQuality(
    char* qual,
    int l_seq)
{
    int result = 0;
    _uint8* p = (_uint8*) qual;
    for (int i = 0; i < l_seq; i++) {
        int q = *p++;
        q -= '!';
        // Picard MarkDup uses a score threshold of 15 (default)
        result += (q >= 15) ? (q != 255) * q : 0; // avoid branch?
    }
    return result;
}

    void
SAMFilter::getTileXY(
    const char* id,
    int* o_tile,
    int* o_x,
    int* o_y)
{
    // In the 5-element format read names: tile, x and y values are 3rd, 4th and 5th elements
    // In the 7-element format (CASAVA 1.8) read names: tile, x and y values are 5th, 6th and 7th elements
    int i = 0, numColonsSeen = 0, tile = 0, x = 0, y = 0;
    char readId[120];
    strncpy(readId, id, sizeof(readId));
    readId[119] = '\0';
    int fiveElementFormatStartIndex = 0, sevenElementFormatStartIndex;
    for (i = 0; ; i++) {
        char c = readId[i];
        if (c == ':') {
            numColonsSeen++;
            if (numColonsSeen == 2) {
                fiveElementFormatStartIndex = i;
            }
            else if (numColonsSeen == 4) {
                sevenElementFormatStartIndex = i;
            }
        }
        if (c == 0 || c == ' ' || c == '/') {
            break;
        }
    }
    int fieldScanned;
    if (numColonsSeen == 4) {
        fieldScanned = sscanf(&readId[fiveElementFormatStartIndex], ":%d:%d:%d", &tile, &x, &y);
    }
    else if (numColonsSeen == 6) {
        fieldScanned = sscanf(&readId[sevenElementFormatStartIndex], ":%d:%d:%d", &tile, &x, &y);
    }
    // fprintf(stderr, "tile:%d, x:%d, y:%d\n", tile, x, y);
    *o_tile = tile;
    *o_x = x;
    *o_y = y;
}

//
// One of these per possible duplicate read (matches on location, next location, RC, next RC).
//
struct DuplicateReadKey
{
    DuplicateReadKey()
    { memset(this, 0, sizeof(DuplicateReadKey)); }

    DuplicateReadKey(SAMDupMarkEntry* sam, const Genome *genome, size_t _libraryHash)
    {
        if (sam == NULL) {
            locations[0] = locations[1] = InvalidGenomeLocation;
            isRC[0] = isRC[1] = false;
            libraryHash = 0;
        } else {
            isRC[0] = (sam->flag & SAM_REVERSE_COMPLEMENT) != 0;
            isRC[1] = (sam->flag & SAM_NEXT_REVERSED) != 0;
            locations[0] = sam->info >> 1;
            locations[1] = locations[0] + sam->tlen;
            libraryHash = _libraryHash;
            if (((((_uint64) GenomeLocationAsInt64(locations[0])) << 1) | (isRC[0] ? 1 : 0)) > ((((_uint64) GenomeLocationAsInt64(locations[1])) << 1) | (isRC[1] ? 1 : 0))) {
                const GenomeLocation t = locations[1];
                locations[1] = locations[0];
                locations[0] = t;
                const bool f = isRC[1];
                isRC[1] = isRC[0];
                isRC[0] = f;
            }
        }
    }

    bool operator==(const DuplicateReadKey& b) const
    {
        return libraryHash == b.libraryHash && locations[0] == b.locations[0] && locations[1] == b.locations[1] &&
            isRC[0] == b.isRC[0] && isRC[1] == b.isRC[1];
    }

    bool operator!=(const DuplicateReadKey& b) const
    {
        return ! ((*this) == b);
    }

    bool operator<(const DuplicateReadKey& b) const
    {
        return libraryHash < b.libraryHash ||
            (libraryHash == b.libraryHash && locations[0] < b.locations[0]) ||
            (locations[0] == b.locations[0] &&
                (locations[1] < b.locations[1] ||
                    (locations[1] == b.locations[1] &&
                        isRC[0] * 2 + isRC[1] <  b.isRC[0] *2 + b.isRC[1])));
    }


    // required for use as a key in VariableSizeMap template
    DuplicateReadKey(int x)
    { locations[0] = locations[1] = x; isRC[0] = isRC[1] = false; }
    bool operator==(int x) const
    { return locations[0] == (_uint32) x && locations[1] == (_uint32) x; }
    bool operator!=(int x) const
    { return locations[0] != (_uint32) x || locations[1] != (_uint32) x; }
    operator _uint64()
    { return ((_uint64) (GenomeLocationAsInt64(locations[1]) ^ (isRC[1] ? 1 : 0))) << 32 | (_uint64) (GenomeLocationAsInt64(locations[0]) ^ (isRC[0] ? 1 : 0)); }

    GenomeLocation locations[NUM_READS_PER_PAIR];
    bool isRC[NUM_READS_PER_PAIR];
    size_t libraryHash; // hash value of library name string LB:Z:xxxx
};

//
// Fragments duplicate marking does not require mate information
//
struct DuplicateFragmentKey
{
    DuplicateFragmentKey()
    {
        memset(this, 0, sizeof(DuplicateFragmentKey));
    }

    DuplicateFragmentKey(SAMDupMarkEntry* sam, const Genome* genome, size_t _libraryHash)
    {
        if (sam == NULL) {
            location = InvalidGenomeLocation;
            isRC = false;
            libraryHash = 0;
        }
        else {
            isRC = (sam->flag & SAM_REVERSE_COMPLEMENT) != 0;
            location = sam->info >> 1;
            libraryHash = _libraryHash;
        }
    }

    bool operator==(const DuplicateFragmentKey& b) const
    {
        return libraryHash == b.libraryHash && location == b.location && isRC == b.isRC;
    }

    bool operator!=(const DuplicateFragmentKey& b) const
    {
        return !((*this) == b);
    }

    bool operator<(const DuplicateFragmentKey& b) const
    {
        return libraryHash < b.libraryHash ||
            (libraryHash == b.libraryHash && location < b.location) ||
            (location == b.location && isRC < b.isRC);
    }

    // required for use as a key in VariableSizeMap template
    DuplicateFragmentKey(int x)
    {
        location = x; isRC = false;
    }
    bool operator==(int x) const
    {
        return location == (_uint32)x;
    }
    bool operator!=(int x) const
    {
        return location != (_uint32)x;
    }
    operator _uint64()
    {
        return (_uint64)(GenomeLocationAsInt64(location) ^ (isRC ? 1 : 0));
    }

    GenomeLocation location;
    bool isRC;
    size_t libraryHash;
};

struct DuplicateMateInfo
{
    DuplicateMateInfo()
    {
        memset(this, 0, sizeof(DuplicateMateInfo));
    }

    ~DuplicateMateInfo()
    {
    }

    bool isMateMapped;
    int bestReadQuality; // total quality of first/both best reads
    char bestReadId[120];
    
    //
    // The below are useful for marking optical duplicates. Helps break ties in Illumina reads
    //
    int tile;
    int x;
    int y;

    void setBestReadId(const char* id) { strncpy(bestReadId, id, sizeof(bestReadId)); }
    const char* getBestReadId() { return bestReadId; }
    void setBestTileXY(int tile_, int x_, int y_) { tile = tile_; x = x_; y = y_; }
    void getBestTileXY(int* tile_, int* x_, int* y_) { *tile_ = tile; *x_ = x; *y_ = y; }

    void checkBestRecord(SAMDupMarkEntry* sam, int totalQuality_, int tile_, int x_, int y_, bool isDuplicate) {

        if (isDuplicate) {
            return;
        }
        if (totalQuality_ > bestReadQuality) {
            bestReadQuality = totalQuality_;
            setBestReadId(sam->qName);
            setBestTileXY(tile_, x_, y_);
        }
        else if (totalQuality_ == bestReadQuality) {
            if (tile_ < tile) {
                bestReadQuality = totalQuality_;
                setBestReadId(sam->qName);
                setBestTileXY(tile_, x_, y_);
            }
            else if (tile_ == tile) {
                if (x_ < x) {
                    bestReadQuality = totalQuality_;
                    setBestReadId(sam->qName);
                    setBestTileXY(tile_, x_, y_);
                }
                else if (x_ == x) {
                    if (y_ < y) {
                        bestReadQuality = totalQuality_;
                        setBestReadId(sam->qName);
                        setBestTileXY(tile_, x_, y_);
                    }
                }
            }
        }
    }
};

class SAMDupMarkFilter : public SAMFilter
{
public:
    SAMDupMarkFilter(const Genome* i_genome) :
        SAMFilter(DataWriter::DupMarkFilter),
        genome(i_genome), runOffset(0), runLocation(InvalidGenomeLocation), prevRunLocation(InvalidGenomeLocation), runCount(0), mates(), fragments(), buffer(NULL), bufferSize(0), bufferUsed(0)
    {
    }

    ~SAMDupMarkFilter()
    {
        if (buffer != NULL) {
            BigDealloc(buffer);
            buffer = NULL;
        }

        run.clear();
        runFragment.clear();
    }

    virtual size_t onNextBatch(DataWriter* writer, size_t offset, size_t bytes, bool lastBatch = false, bool* needMoreBuffer = NULL, size_t* fromBufferUsed = NULL);

    void dupMarkBatch(size_t runStartIndex);

    void updateBuffer(size_t lastByte);

private:

    const Genome* genome;
    size_t runOffset; // offset in file of first read in run
    GenomeLocation runLocation; // location in genome
    GenomeLocation prevRunLocation; // location in genome
    int runCount; // number of aligned reads

    typedef VariableSizeMap<DuplicateReadKey,DuplicateMateInfo,150,MapNumericHash<DuplicateReadKey>,70,0,-2> MateMap;
    typedef VariableSizeMap<DuplicateFragmentKey, DuplicateMateInfo, 150, MapNumericHash<DuplicateFragmentKey>, 70, 0, -2> FragmentMap;
    typedef VariableSizeVector<SAMDupMarkEntry> RunVector;

    RunVector run; // used for paired-end duplicate marking
    RunVector runFragment; // used for single-end duplicate marking 
    MateMap mates;
    FragmentMap fragments;
    char* buffer; // store results after duplicate marking here
    size_t bufferSizeInit; // initial buffer size. Shrink back buffer to this size after MarkDup
    size_t bufferSize;
    size_t bufferUsed;

};

    void
SAMDupMarkFilter::updateBuffer(size_t lastByte)
{
    // fprintf(stderr, "Buffer used: %lld, Buffer size: %lld\n", bufferUsed, bufferSize);
    //
    // Update flags for duplicate records
    //
    size_t nRecords = offsets.size();
    for (size_t i = 0; i < nRecords; i++) {
        if (offsets[i].offset >= lastByte) {
            break;
        }
        size_t recordSize = (i != nRecords - 1) ? offsets[i + 1].offset - offsets[i].offset : currentBufferBytes - offsets[i].offset;
        if (!offsets[i].isDuplicate) {
            memcpy(buffer + bufferUsed, currentBuffer + offsets[i].offset, recordSize);
            bufferUsed += recordSize;
        } else {
            updateSAMLine(buffer, &bufferUsed, currentBuffer + offsets[i].offset, recordSize);
        }
    }
}

    size_t
SAMDupMarkFilter::onNextBatch(
    DataWriter* writer,
    size_t offset,
    size_t bytes,
    bool lastBatch,
    bool* needMoreBuffer,
    size_t* fromBytesUsed)
{
    // 
    // Nothing to write
    //
    if (bytes == 0 || *needMoreBuffer) {
        return 0;
    }

    size_t currentBufferSize;
    bool ok = writer->getBatch(-1, &currentBuffer, &currentBufferSize, NULL, NULL, &currentBufferBytes, &currentOffset);
    if (!ok) {
        WriteErrorMessage("Error writing to output file\n");
        soft_exit(1);
    }
    currentWriter = writer;

    if (buffer == NULL) {
        buffer = (char*)BigAlloc(currentBufferSize);
        if (buffer == NULL) {
            WriteErrorMessage("Unable to allocate %lld bytes for SAM MarkDup buffer\n", currentBufferSize);
            soft_exit(1);
        }
        bufferSize = currentBufferSize;
        bufferSizeInit = bufferSize;
    }
    else if (currentBufferSize > bufferSize) {
        BigDealloc(buffer);
        buffer = (char*)BigAlloc(currentBufferSize);
        if (buffer == NULL) {
            WriteErrorMessage("Unable to allocate %lld bytes for SAM MarkDup buffer\n", currentBufferSize);
            soft_exit(1);
        }
        bufferSize = currentBufferSize;
    }
    bufferUsed = 0;

    size_t next_i = 0, unfinishedRunStart = 0;
    bool foundNextBatchStart = false;
    size_t nextBatchStart = 0;
    runCount = 0;
    runOffset = 0;
    runLocation = InvalidGenomeLocation;
    run.clear();
    runFragment.clear();

    SAMAlignment firstSam, lastSam;
    size_t nRecords = offsets.size();

    for (size_t i = 0; i < nRecords; i = next_i) {
        size_t recordSize = (i != nRecords - 1) ? offsets[i + 1].offset - offsets[i].offset : currentBufferBytes - offsets[i].offset;
        getAlignmentInfo(genome, currentBuffer + offsets[i].offset, recordSize, &lastSam);
        GenomeLocation logicalLocation = lastSam.unclippedStartLocation;

        //
        // Initialize run
        //
        if (runLocation == InvalidGenomeLocation) {
            runCount = 1;
            runLocation = logicalLocation;
            runOffset = currentOffset + offsets[i].offset;
            unfinishedRunStart = i;
            next_i = i + 1;

            if ((lastSam.flag & SAM_SECONDARY) == 0 && (lastSam.flag & SAM_SUPPLEMENTARY) == 0) {
                bool isRC = (lastSam.flag & SAM_REVERSE_COMPLEMENT) != 0;
                bool isMateRC = (lastSam.flag & SAM_NEXT_REVERSED) != 0;
                GenomeLocation myLoc = isRC ? lastSam.unclippedEndLocation : lastSam.unclippedStartLocation;
                GenomeLocation mateLoc = myLoc + lastSam.tlen;

                SAMDupMarkEntry entry(lastSam.libraryNameHash, i, lastSam.flag, lastSam.qual, lastSam.mateQual, 
                    (((_uint64)GenomeLocationAsInt64(mateLoc)) << 1) | (isMateRC ? 1 : 0), lastSam.mateLocation,
                    (((_uint64)GenomeLocationAsInt64(myLoc)) << 1) | (isRC ? 1 : 0), lastSam.tlen, lastSam.qName, lastSam.qNameLength);

                SAMDupMarkEntry entryFragment(lastSam.libraryNameHash, i, lastSam.flag, lastSam.qual, 0, 0, 0,
                    (((_uint64)GenomeLocationAsInt64(myLoc)) << 1) | (isRC ? 1 : 0), 0, lastSam.qName, lastSam.qNameLength);

                // don't need to look at mate information in single-ended datasets
                if ((lastSam.flag & SAM_MULTI_SEGMENT) != 0) {
                    run.push_back(entry);
                }

                runFragment.push_back(entryFragment);
            }
        } else {
            // 
            // Track the read from which we need to start the next run. Next run starts at nextBatchStart
            // 
            if (!foundNextBatchStart && logicalLocation > runLocation + MAX_READ_LENGTH + MAX_K) {
                nextBatchStart = i;
                foundNextBatchStart = true;
                next_i = i + 1;
            }
            //
            // Add read to run
            // 
            if (logicalLocation <= runLocation + (2 * (MAX_READ_LENGTH + MAX_K))) {
                runCount++;
                next_i = i + 1;

                if ((lastSam.flag & SAM_SECONDARY) == 0 && (lastSam.flag & SAM_SUPPLEMENTARY) == 0) {
                    bool isRC = (lastSam.flag & SAM_REVERSE_COMPLEMENT) != 0;
                    bool isMateRC = (lastSam.flag & SAM_NEXT_REVERSED) != 0;
                    GenomeLocation myLoc = isRC ? lastSam.unclippedEndLocation : lastSam.unclippedStartLocation;
                    GenomeLocation mateLoc = myLoc + lastSam.tlen;

                    SAMDupMarkEntry entry(lastSam.libraryNameHash, i, lastSam.flag, lastSam.qual, lastSam.mateQual,
                        (((_uint64)GenomeLocationAsInt64(mateLoc)) << 1) | (isMateRC ? 1 : 0), lastSam.mateLocation,
                        (((_uint64)GenomeLocationAsInt64(myLoc)) << 1) | (isRC ? 1 : 0), lastSam.tlen, lastSam.qName, lastSam.qNameLength);

                    SAMDupMarkEntry entryFragment(lastSam.libraryNameHash, i, lastSam.flag, lastSam.qual, 0, 0, 0,
                        (((_uint64)GenomeLocationAsInt64(myLoc)) << 1) | (isRC ? 1 : 0), 0, lastSam.qName, lastSam.qNameLength);

                    // don't need to look at mate information in single-ended datasets
                    if ((lastSam.flag & SAM_MULTI_SEGMENT) != 0) {
                        run.push_back(entry);
                    }

                    runFragment.push_back(entryFragment);
                }

            } else {
                // 
                // We are done with the run. Begin marking duplicates
                // 
                dupMarkBatch(unfinishedRunStart);

                //
                // Reset run parameters in preparation for next run
                // 
                runCount = 0;
                runOffset = 0;
                runLocation = InvalidGenomeLocation;
                run.clear();
                runFragment.clear();
                next_i = nextBatchStart;
                foundNextBatchStart = false;
            }
        }
    } // end offsets

    size_t bytesRead = min<long long>(bytes, currentBufferBytes);

    //
    // Run could potentially span multiple buffers
    //
    if (runCount > 1) { // runs of size <= 1 cannot have duplicates
        if (lastBatch) {
            dupMarkBatch(unfinishedRunStart);

            updateBuffer(bytesRead);

            // Copy over updated buffer
            memcpy(currentBuffer + offsets[0].offset, buffer, bufferUsed);

            offsets.clear();

        } else if (offsets[unfinishedRunStart].offset == 0) { // we did not yet mark duplicates for this run
            //
            // If we have a different run from what we have seen before, 
            // simply mark all duplicates in the run we currently have
            //
            if (runLocation != prevRunLocation) {
                dupMarkBatch(unfinishedRunStart);

                updateBuffer(bytesRead);

                // Copy over updated buffer
                memcpy(currentBuffer + offsets[0].offset, buffer, bufferUsed);

                offsets.clear();
            } else {
                *needMoreBuffer = true;
                return 0;
            }
        } else {
            bytesRead = offsets[unfinishedRunStart].offset;
            updateBuffer(bytesRead);

            // Copy over updated buffer
            memcpy(currentBuffer + offsets[0].offset, buffer, bufferUsed);

            //
            // Copy over read offsets for those reads that will be duplicate marked in the next batch
            //
            VariableSizeVector<OffsetInfo> nextBatchOffsets;
            for (_int64 i = unfinishedRunStart; i < offsets.size(); i++) {
                OffsetInfo oinfo(offsets[i].offset - offsets[unfinishedRunStart].offset, offsets[i].isDuplicate);
                nextBatchOffsets.push_back(oinfo);
            }

            offsets.clear();

            for (_int64 i = 0; i < nextBatchOffsets.size(); i++) {
                offsets.push_back(nextBatchOffsets[i]);
            }
            nextBatchOffsets.clear();

        }
    } // runcount > 1
    else {
        if (runCount > 0) {
            //
            // Commit any duplicate marking updates
            //
            updateBuffer(bytesRead);
            // Copy over updated buffer
            memcpy(currentBuffer + offsets[0].offset, buffer, bufferUsed);
        }
        offsets.clear();
    }

    prevRunLocation = runLocation;
    currentWriter = NULL;
    currentBuffer = NULL;
    currentBufferBytes = 0;
    currentOffset = 0;
    *fromBytesUsed = bytesRead;

    //
    // Shrink buffers to reduce memory consumption
    //
    if (bufferSize > bufferSizeInit) {
        size_t newBufferSize = bufferSizeInit;
        char* newBuffer = (char*)BigAlloc(newBufferSize);
        if (newBuffer == NULL) {
            WriteErrorMessage("Unable to allocate %lld bytes for MarkDup buffer\n", newBufferSize);
            soft_exit(1);
        }
        BigDealloc(buffer);
        buffer = newBuffer;
        bufferSize = newBufferSize;
    }

    return (runCount > 0) ? bufferUsed : bytesRead; // return bytes written in current batch.
}

    void
SAMDupMarkFilter::dupMarkBatch(size_t runStartIndex) {

    size_t nRecords = offsets.size();

    if (run.size() == 0 && runFragment.size() == 0) {
        return;
    }
    // ensure that duplicates reads are adjacent
    std::stable_sort(run.begin(), run.end());
    std::stable_sort(runFragment.begin(), runFragment.end());

    bool foundRun = false;
    size_t offsetIndex = 0;
    for (RunVector::iterator i = run.begin(); i != run.end(); i++) {
      
        // skip unmapped reads and reads with unmapped mates
        if (((i->flag & SAM_UNMAPPED) != 0) || (i->flag & SAM_NEXT_UNMAPPED) != 0) continue;

        bool prevRecordMismatch = (i == run.begin()) || ((i->libraryNameHash != ((i - 1)->libraryNameHash)) ||
            (i->info != (i - 1)->info) ||
            (i->mateInfo != (i - 1)->mateInfo));

        bool nextRecordMismatch = (i + 1 == run.end()) || ((i->libraryNameHash != ((i + 1)->libraryNameHash)) ||
            (i->info != (i + 1)->info) ||
            (i->mateInfo != (i + 1)->mateInfo));

        // adjacent entries with different library names, different location/orientation, different mate location/orientation
        if (prevRecordMismatch && nextRecordMismatch) {
            continue;
        }

        size_t offsetIndex = i->runOffset;
        foundRun = true;
        DuplicateReadKey key(i, genome, i->libraryNameHash);
        MateMap::iterator f = mates.find(key);
        DuplicateMateInfo* info;
        if (f == mates.end()) {
            mates.put(key, DuplicateMateInfo());
            info = &mates[key];
            //fprintf(stderr, "add %u%s/%u%s -> %d\n", key.locations[0], key.isRC[0] ? "rc" : "", key.locations[1], key.isRC[1] ? "rc" : "", mates.size());
            info->isMateMapped = true;
        } else {
            info = &f->value;
        }
        int totalQuality = i->qual;
        if ((i->flag & SAM_MULTI_SEGMENT) != 0) {
            totalQuality += i->mateQual;
        }
        int tile, x, y;
        getTileXY(i->qName, &tile, &x, &y); // parse read name and extract metadata for optical duplicate marking
        info->checkBestRecord(i, totalQuality, tile, x, y, offsets[offsetIndex].isDuplicate); // update best record if needed
    }

    // duplicate marking for read fragments
    for (RunVector::iterator i = runFragment.begin(); i != runFragment.end(); i++) {

        //
        // Skip unmapped reads
        //
        if ((i->flag & SAM_UNMAPPED) != 0) {
            continue;
        }

        if ((i == runFragment.begin() || (i->libraryNameHash != ((i - 1)->libraryNameHash))) &&
            (i + 1 == runFragment.end() || (i->libraryNameHash != ((i + 1)->libraryNameHash)))) {
            continue;
        }

        if ((i == runFragment.begin() || (i->info) != ((i - 1)->info)) &&
            (i + 1 == runFragment.end() || (i->info) != ((i + 1)->info))) {
            continue;
        }

        size_t offsetIndex = i->runOffset;
        
        // location and library matches
        foundRun = true;
        DuplicateFragmentKey key(i, genome, i->libraryNameHash);
        FragmentMap::iterator f = fragments.find(key);
        DuplicateMateInfo* info;
        if (f == fragments.end()) {
            fragments.put(key, DuplicateMateInfo());
            info = &fragments[key];
            //fprintf(stderr, "add %u%s/%u%s -> %d\n", key.locations[0], key.isRC[0] ? "rc" : "", key.locations[1], key.isRC[1] ? "rc" : "", mates.size());
            info->isMateMapped = (i->flag & SAM_MULTI_SEGMENT) != 0 && (i->flag & SAM_NEXT_UNMAPPED) == 0;
        } else {
            info = &f->value;
        }

        bool mateMapped = (i->flag & SAM_MULTI_SEGMENT) != 0 && (i->flag & SAM_NEXT_UNMAPPED) == 0;
        int totalQuality = i->qual;
        int tile, x, y;
        getTileXY(i->qName, &tile, &x, &y);
        
        //
        // Prefer mapped read pairs over fragments in duplicate marking
        //
        if (mateMapped) {
            if (!info->isMateMapped) {
                info->bestReadQuality = totalQuality;
                info->setBestReadId(i->qName);
                info->isMateMapped = true;
                info->setBestTileXY(tile, x, y);
            } else {
                info->checkBestRecord(i, totalQuality, tile, x, y, offsets[offsetIndex].isDuplicate);
            }
        } else {
            //
            // No best read pair found so far.
            //
            if (!info->isMateMapped) {
                info->checkBestRecord(i, totalQuality, tile, x, y, offsets[offsetIndex].isDuplicate);
            }
        }
    }

    if (!foundRun) {
        return; // avoid useless looping
    }

    // go back and adjust flags
    for (RunVector::iterator i = run.begin(); i != run.end(); i++) {

        // skip unmapped reads and reads with unmapped mates
        if (((i->flag & SAM_UNMAPPED) != 0) || (i->flag & SAM_NEXT_UNMAPPED) != 0) continue;

        bool prevRecordMismatch = (i == run.begin()) || ((i->libraryNameHash != ((i - 1)->libraryNameHash)) ||
            (i->info != (i - 1)->info) ||
            (i->mateInfo != (i - 1)->mateInfo));

        bool nextRecordMismatch = (i + 1 == run.end()) || ((i->libraryNameHash != ((i + 1)->libraryNameHash)) ||
            (i->info != (i + 1)->info) ||
            (i->mateInfo != (i + 1)->mateInfo));


        // adjacent entries with different library names, different location/orientation, different mate location/orientation
        if (prevRecordMismatch && nextRecordMismatch) {
            continue;
        }

        size_t offsetIndex = i->runOffset;

        DuplicateReadKey key(i, genome, i->libraryNameHash);
        MateMap::iterator m = mates.find(key);
        if (m == mates.end()) {
            continue;
        }

        DuplicateMateInfo* minfo = &m->value;
        if (!readIdsMatch(minfo->getBestReadId(), i->qName, i->qNameLength)) {
            i->flag |= SAM_DUPLICATE;
            offsets[offsetIndex].isDuplicate |= true;
        }
    }

    // clean up
    for (RunVector::iterator i = run.begin(); i != run.end(); i++) {

        bool prevRecordMismatch = (i == run.begin()) || ((i->libraryNameHash != ((i - 1)->libraryNameHash)) ||
            (i->info != (i - 1)->info) ||
            (i->mateInfo != (i - 1)->mateInfo));

        bool nextRecordMismatch = (i + 1 == run.end()) || ((i->libraryNameHash != ((i + 1)->libraryNameHash)) ||
            (i->info != (i + 1)->info) ||
            (i->mateInfo != (i + 1)->mateInfo));

        // adjacent entries with different library names, different location/orientation, different mate location/orientation
        if (prevRecordMismatch && nextRecordMismatch) {
            continue;
        }

        size_t offsetIndex = i->runOffset;

        // skip unmapped reads and reads with unmapped mates
        if (((i->flag & SAM_UNMAPPED) != 0) || (i->flag & SAM_NEXT_UNMAPPED) != 0) continue;

        DuplicateReadKey key(i, genome, i->libraryNameHash);
        MateMap::iterator m = mates.find(key);
        DuplicateMateInfo* minfo = &m->value;
        if (m != mates.end()) {
            bool isRC = (i->flag & SAM_REVERSE_COMPLEMENT) != 0;
            GenomeLocation loc = i->info >> 1;

            GenomeLocation nextLoc = i->mateLocation;
            GenomeDistance spacing = nextLoc > loc ? nextLoc - loc : loc - nextLoc;
            // 
            // Keep duplicate entry around till we find the mate. This allows us to match indexInFile tie breaking used in Picard MarkDup.
            // TODO: We have not ensured this when the read and its mate are mapped more than INT_MAX apart, since TLEN is 32-bit signed.
            //       The SAM spec says TLEN = 0 for reads mapped to different chromosomes. In these cases we must store the offset of the
            //       mate in an auxiliary tag.
            //
            if (spacing > INT_MAX || (loc == key.locations[1] && isRC == key.isRC[1]) ||
                (loc == key.locations[0] && isRC == key.isRC[0] && (loc > nextLoc))) {
                //fprintf(stderr, "erase %u%s/%u%s -> %d\n", key.locations[0], key.isRC[0] ? "rc" : "", key.locations[1], key.isRC[1] ? "rc" : "", mates.size());
                mates.erase(key);
            }
        }
    }

    // handle fragments
    for (RunVector::iterator i = runFragment.begin(); i != runFragment.end(); i++) {

        if ((i == runFragment.begin() || (i->libraryNameHash != ((i - 1)->libraryNameHash))) &&
            (i + 1 == runFragment.end() || (i->libraryNameHash != ((i + 1)->libraryNameHash)))) {
            continue;
        }

        if ((i == runFragment.begin() || (i->info) != ((i - 1)->info)) &&
            (i + 1 == runFragment.end() || (i->info) != ((i + 1)->info))) {
            continue;
        }

        size_t offsetIndex = i->runOffset;

        //
        // Skip unmapped reads and reads with mapped mates
        //
        if ((i->flag & SAM_UNMAPPED) != 0 || ((i->flag & SAM_MULTI_SEGMENT) != 0 && (i->flag & SAM_NEXT_UNMAPPED) == 0)) {
            continue;
        }

        // location and library matches
        DuplicateFragmentKey key(i, genome, i->libraryNameHash);
        FragmentMap::iterator f = fragments.find(key);

        if (f == fragments.end()) {
            continue;
        }
        DuplicateMateInfo* info = &f->value;
        if (!readIdsMatch(info->getBestReadId(), i->qName, i->qNameLength)) {
            i->flag |= SAM_DUPLICATE;
            offsets[offsetIndex].isDuplicate |= true;
        }
    }

    // clean up
    fragments.clear();

    run.clear();
    runFragment.clear();
}

class SAMDupMarkSupplier : public DataWriter::FilterSupplier
{
public:
    SAMDupMarkSupplier(const Genome* i_genome) :
        FilterSupplier(DataWriter::ReadFilter), genome(i_genome) {}

    virtual DataWriter::Filter* getFilter()
    {
        return new SAMDupMarkFilter(genome);
    }

    virtual void onClosing(DataWriterSupplier* supplier) {}
    virtual void onClosed(DataWriterSupplier* supplier) {}

private:
    const Genome* genome;
};

    DataWriter::FilterSupplier*
DataWriterSupplier::samMarkDuplicates(const Genome* genome)
{
    return new SAMDupMarkSupplier(genome);
}
