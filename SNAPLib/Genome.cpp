/*++

Module Name:

    Geonome.cpp

Abstract:

    Genome class for the SNAP sequencer

Authors:

    Bill Bolosky, August, 2011

Environment:

    User mode service.

Revision History:

    Adapted from Matei Zaharia's Scala implementation.

--*/

#include "stdafx.h"
#include "Genome.h"
#include "GenericFile.h"
#include "GenericFile_map.h"
#include "Compat.h"
#include "BigAlloc.h"
#include "exit.h"
#include "Error.h"
#include "Util.h"

Genome::Genome(GenomeDistance i_maxBases, GenomeDistance nBasesStored, unsigned i_chromosomePadding, unsigned i_maxContigs)
: maxBases(i_maxBases), minLocation(0), maxLocation(i_maxBases), chromosomePadding(i_chromosomePadding), maxContigs(i_maxContigs), mappedFile(NULL)
{
    bases = ((char *) BigAlloc(nBasesStored + 2 * N_PADDING)) + N_PADDING;
    if (NULL == bases) {
        WriteErrorMessage("Genome: unable to allocate memory for %llu bases\n", GenomeLocationAsInt64(maxBases));
        soft_exit(1);
    }

    // Add N's for the N_PADDING bases before and after the genome itself
    memset(bases - N_PADDING, 'n', N_PADDING);
    memset(bases + nBasesStored, 'n', N_PADDING);

    nBases = 0;

    nContigs = 0;
    contigs = (Contig *)BigAlloc(sizeof(Contig) * maxContigs);
    contigNumberByOriginalOrder = (InternalContigNum *)BigAlloc(sizeof(InternalContigNum) * maxContigs);
    contigsByName = NULL;
}

    void
Genome::addData(const char *data, GenomeDistance len)
{
    if (nBases + len > GenomeLocationAsInt64(maxBases)) {
        WriteErrorMessage("Tried to write beyond allocated genome size (or tried to write into a genome that was loaded from a file).\n"
                          "Size = %lld\n", GenomeLocationAsInt64(maxBases));
        soft_exit(1);
    }

    memcpy(bases + nBases,data,len);
    nBases += (unsigned)len;
}

    void
Genome::addData(const char *data)
{
    addData(data, strlen(data));
}

    void
Genome::startContig(const char *contigName, int originalContigNumber)
{
    if (nContigs == maxContigs) {
        //
        // Reallocate (maybe we're sequencing a tree that's got lots of chromosomes).
        //
        int newMaxContigs = maxContigs * 2;
        Contig *newContigs = new Contig[newMaxContigs];
        if (NULL == newContigs) {
            WriteErrorMessage("Genome: unable to reallocate contig array to size %d\n", newMaxContigs);
            soft_exit(1);
        }
        for (int i = 0; i < nContigs; i++) {
            newContigs[i] = contigs[i];
        }

        delete [] contigs;
        contigs = newContigs;
        maxContigs = newMaxContigs;
    }

    contigs[nContigs].beginningLocation = nBases;
    size_t len = strlen(contigName) + 1;
    contigs[nContigs].name = new char[len];
    contigs[nContigs].nameLength = (unsigned)len-1;

    strncpy(contigs[nContigs].name,contigName,len);
    contigs[nContigs].name[len-1] = '\0';
	contigs[nContigs].isALT = false;
    contigs[nContigs].originalContigNumber = originalContigNumber;

    nContigs++;
}


	void 
Genome::markContigALT(
	const char			*contigName)
{
	for (int i = 0; i < nContigs; i++) {
		if (!strcmp(contigs[i].name, contigName)) {
			contigs[i].isALT = true;
			return;
		} // if it matches
	} // for each contig
} // markContigALT

    int
Genome::getContigNumByName(
    const char* contigName) const
{
    for (int i = 0; i < nContigs; i++) {
        if (!strcmp(contigs[i].name, contigName)) {
            return i;
        } // if it matches
    }
    return -1;
}

    GenomeLocation
Genome::getBeginningLocation(
    const char* contigName) const
{
    int contigNum = getContigNumByName(contigName);
    return contigs[contigNum].beginningLocation;
}

    void
Genome::markContigLiftover(
    char *contigName,
    char** altLiftoverContigNames,
    unsigned* altLiftoverContigFlags,
    char** altLiftoverProjContigNames,
    unsigned* altLiftoverProjContigOffsets,
    char** altLiftoverProjCigar,
    int nAltLiftover)
{
    int contigNum = getContigNumByName(contigName);
    for (int i = 0; i < nAltLiftover; i++) {
        if (NULL != altLiftoverContigNames[i]) {
            if (!strcmp(altLiftoverContigNames[i], contigName)) {
                int projContigNum = getContigNumByName(altLiftoverProjContigNames[i]);
                contigs[contigNum].projBeginningLocation = contigs[projContigNum].beginningLocation + altLiftoverProjContigOffsets[i] - 1;
                contigs[contigNum].isProjRC = (altLiftoverContigFlags[i] & 16) != 0; // bit 4 for rc
                contigs[contigNum].projCigar = altLiftoverProjCigar[i];
                return;
            }
        }
    }
}

Genome::~Genome()
{
    for (int i = 0; i < nContigs; i++) {
        delete [] contigs[i].name;
        contigs[i].name = NULL;
        delete [] contigs[i].projCigar;
        contigs[i].projCigar = NULL;
        delete [] contigs[i].projCigarOps;
        contigs[i].projCigarOps = NULL;
    }

    BigDealloc(contigs);
    BigDealloc(contigNumberByOriginalOrder);

    if (contigsByName) {
        delete [] contigsByName;
    }
    contigs = NULL;

    if (NULL != mappedFile) {
        mappedFile->close();
        delete mappedFile;
    }
    else {
        BigDealloc(bases - N_PADDING);
    }
}

// Flags for the options field in the header
#define	GENOME_FLAG_ALT_CONTIGS_MARKED		            0x1 // This should always be true now, we dropped support for indices that don't mark ALTs in 1.0.4.  (That doesn't mean that they have to have alts marked, it's just the format.)

// Flags for the per-contig options
#define GENOME_FLAG_CONTIG_IS_ALT			0x1

#define GENOME_FLAG_ALT_PROJ_CONTIG_IS_RC   0x1
    bool
Genome::saveToFile(const char *fileName) const
{
    //
    // Save file format is the number of bases, the number of contigs and flags followed by
    //  the contigs themselves, rounded up to 4K, followed by the bases.
    //

    FILE *saveFile = fopen(fileName,"wb");
    if (saveFile == NULL) {
        WriteErrorMessage("Genome::saveToFile: unable to open file '%s'\n",fileName);
        return false;
    } 

    fprintf(saveFile,"%lld %d %d\n",nBases, nContigs, GENOME_FLAG_ALT_CONTIGS_MARKED);	
    char *curChar = NULL;

    for (int i = 0; i < nContigs; i++) {
        for (int n = 0; n < strlen(contigs[i].name); n++)
        {
             curChar = contigs[i].name + n;
             if (*curChar == ' '){ *curChar = '_'; }
        }

        fprintf(saveFile, "%lld %x %d %lld %x %d %d %s %s\n", GenomeLocationAsInt64(contigs[i].beginningLocation), (contigs[i].isALT ? GENOME_FLAG_CONTIG_IS_ALT : 0), OriginalContigNumToInt(contigs[i].originalContigNumber), GenomeLocationAsInt64(contigs[i].projBeginningLocation),
            contigs[i].isProjRC ? GENOME_FLAG_ALT_PROJ_CONTIG_IS_RC : 0, (int) strlen(contigs[i].name), contigs[i].projCigar != NULL ? (int) strlen(contigs[i].projCigar) : 1, contigs[i].name, (contigs[i].projCigar != NULL) ? contigs[i].projCigar : "*");
    }

	//
	// Write it out in (big) chunks.  For whatever reason, fwrite with really big sizes seems not to
	// work as well as one would like.
	//
	const size_t max_chunk_size = 1 * 1024 * 1024 * 1024;	// 1 GB (or GiB for the obsessively precise)

	size_t bases_to_write = nBases;
	size_t bases_written = 0;
	while (bases_to_write > 0) {
		size_t bases_this_write = __min(bases_to_write, max_chunk_size);
		if (bases_this_write != fwrite(bases + bases_written, 1, bases_this_write, saveFile)) {
			WriteErrorMessage("Genome::saveToFile: fwrite failed\n");
			fclose(saveFile);
			return false;
		}
		bases_to_write -= bases_this_write;
		bases_written += bases_this_write;
	}

	_ASSERT(bases_written == nBases);

    fclose(saveFile);
    return true;
}

    void
Genome::setUpContigNumbersByOriginalOrder()
{
    //
    // This assumes it's already allocated.
    //

    memset(contigNumberByOriginalOrder, -1, sizeof(InternalContigNum) * nContigs);

    for (int i = 0; i < nContigs; i++) {
        OriginalContigNum originalContigNumber = contigs[i].originalContigNumber;

        if (OriginalContigNumToInt(originalContigNumber) < 0 || OriginalContigNumToInt(originalContigNumber) >= nContigs || InternalContigNumToInt(contigNumberByOriginalOrder[OriginalContigNumToInt(originalContigNumber)]) != -1) {
            WriteErrorMessage("Something's wrong with the original contig numbers in the genome.  If rebuilding the index doesn't help, please submit a bug report\n");
            soft_exit(1);
        }

        contigNumberByOriginalOrder[OriginalContigNumToInt(originalContigNumber)] = InternalContigNum(i);
    } // for each contig
} // setUpContigNumbersByOriginalOrder

    const Genome *
Genome::loadFromFile(const char *fileName, unsigned chromosomePadding, GenomeLocation minLocation, GenomeDistance length, bool map)
{    
    GenericFile *loadFile;
    GenomeDistance nBases;
    unsigned nContigs;

    if (!openFileAndGetSizes(fileName, &loadFile, &nBases, &nContigs, map)) {
        //
        // It already printed an error.  Just fail.
        //
        return NULL;
    }

    GenomeLocation maxLocation(nBases); 

    if (0 == length) {
        length = maxLocation - minLocation;
    } else {
        //
        // Don't let length go beyond nBases.
        //
        length = __min(length, maxLocation - minLocation);
        maxLocation = minLocation + length;
    }

    Genome *genome = new Genome(nBases, length, chromosomePadding, nContigs);
   
    genome->nBases = nBases;
    genome->nContigs = genome->maxContigs = nContigs;
    genome->contigs = (Contig *)BigAlloc(sizeof(Contig) * nContigs);
    genome->contigNumberByOriginalOrder = (InternalContigNum *)BigAlloc(sizeof(InternalContigNum) * nContigs);

    genome->minLocation = minLocation;
    if (GenomeLocationAsInt64(minLocation) >= nBases) {
        WriteErrorMessage("Genome::loadFromFile: specified minOffset %u >= nBases %u\n", GenomeLocationAsInt64(minLocation), nBases);
        soft_exit(-1);
    }

    genome->maxLocation = maxLocation;

    int contigNameBufferSize = 0;
    char *contigNameBuffer = NULL;
    size_t n;
    char *curName;
    for (unsigned i = 0; i < nContigs; i++) {
        if (NULL == reallocatingFgetsGenericFile(&contigNameBuffer, &contigNameBufferSize, loadFile)) {	 
            WriteErrorMessage("Unable to read contig description\n");
            delete[] contigNameBuffer;
            delete genome;
            return NULL;
        }

        int spacesFound = 0;
        n = 0;
        for (;;) {
            if (contigNameBuffer[n] == ' ') {
                spacesFound++;
                if (spacesFound >= 7) {
                    contigNameBuffer[n] = '\0';
                    break;
                } // If we've found all three
            } // if it's a space

            n++;

            if (n >= contigNameBufferSize) {
                WriteErrorMessage("Corrupt contig line in genome file found.  Try rebuilding the index.  If that doesn't work, submit a bug report.  Line: '%s'\n", contigNameBuffer);
                soft_exit(1);
            }   // if we got to the end without finding enough spaces
        } // for ever

        _int64 contigStart;
        _int64 projContigStart;
        int contigNameSize;
        int cigarLength;
        int originalContigNumber = -1;
        int flags;
        int altProjRC;

        if (7 != sscanf(contigNameBuffer, "%lld %x %d %lld %x %d %d", &contigStart, &flags, &originalContigNumber, &projContigStart, &altProjRC, &contigNameSize, &cigarLength)) {
            WriteErrorMessage("Unable to parse contig start flags or original contig number in genome file '%s', '%s%'\n", fileName, contigNameBuffer);
            soft_exit(1);
        }

        genome->contigs[i].beginningLocation = GenomeLocation(contigStart);
        genome->contigs[i].projBeginningLocation = GenomeLocation(projContigStart);
	    contigNameBuffer[n] = ' '; 
	    n++; // increment n so we start copying at the position after the space

		genome->contigs[i].isALT = (flags & GENOME_FLAG_CONTIG_IS_ALT) != 0;
        genome->contigs[i].internalContigNumber = InternalContigNum(i);
        genome->contigs[i].originalContigNumber = OriginalContigNum(originalContigNumber);

        genome->contigs[i].name = new char[contigNameSize + 1];
        genome->contigs[i].nameLength = (unsigned)contigNameSize;
	    curName = genome->contigs[i].name;
	    for (int pos = 0; pos < contigNameSize; pos++) {
	      curName[pos] = contigNameBuffer[pos + n];
	    }
        curName[contigNameSize] = '\0';

        n += contigNameSize;
        contigNameBuffer[n] = ' ';
        n++;

        if (cigarLength > 1) {
            genome->contigs[i].projCigar = new char[cigarLength + 1];
            for (int pos = 0; pos < cigarLength; pos++) {
                genome->contigs[i].projCigar[pos] = contigNameBuffer[pos + n];
            }
            genome->contigs[i].projCigar[cigarLength] = '\0';
            genome->contigs[i].isProjRC = (altProjRC & GENOME_FLAG_ALT_PROJ_CONTIG_IS_RC) != 0;
            genome->contigs[i].projCigarOps = new ProjCigarOpFormat[cigarLength];
            char* nextChunkOfCigar = genome->contigs[i].projCigar;
            ProjCigarOpFormat* cigarOps = genome->contigs[i].projCigarOps;
            int scannedOps = 0;
            while ('\0' != *nextChunkOfCigar) {
                int fieldsScanned = sscanf(nextChunkOfCigar, "%d%c", &cigarOps[scannedOps].count, &cigarOps[scannedOps].action);
                scannedOps++;
                while ('0' <= *nextChunkOfCigar && '9' >= *nextChunkOfCigar) {
                    nextChunkOfCigar++;
                }
                nextChunkOfCigar++;
                _ASSERT(scannedOps <= cigarLength);
            }
            genome->contigs[i].countProjCigarOps = scannedOps;
        }
    } // for each contig



    if (0 != loadFile->advance(GenomeLocationAsInt64(minLocation))) {
        WriteErrorMessage("Genome::loadFromFile: _fseek64bit failed\n");
        soft_exit(1);
    }

    size_t readSize;
    if (map) {
        GenericFile_map *mappedFile = (GenericFile_map *)loadFile;
        BigDealloc(genome->bases - N_PADDING);
        genome->bases = (char *)mappedFile->mapAndAdvance(length, &readSize);
        genome->mappedFile = mappedFile;
        mappedFile->prefetch();
    } else {
        readSize = loadFile->read(genome->bases, length);

        loadFile->close();
        delete loadFile;
        loadFile = NULL;
    }

	if (length != readSize) {
		WriteErrorMessage("Genome::loadFromFile: fread of bases failed; wanted %u, got %d\n", length, readSize);
		delete loadFile;
		delete genome;
		return NULL;
	}
	
	genome->fillInContigLengths();
    genome->sortContigsByName();
    genome->setUpContigNumbersByOriginalOrder();
    delete[] contigNameBuffer;
    return genome;
}

    bool
contigComparator(
    const Genome::Contig& a,
    const Genome::Contig& b)
{
    return strcmp(a.name, b.name) < 0;
}

    void
Genome::sortContigsByName()
{
    if (contigsByName) {
        delete [] contigsByName;
    }
    contigsByName = new Contig[nContigs];
    memcpy(contigsByName, contigs, nContigs * sizeof(Contig));
    std::sort(contigsByName, contigsByName + nContigs, contigComparator);

    GenomeLocation highestNonALTContig = 0;
    GenomeLocation lowestALTContig = LLONG_MAX;

    for (int i = 0; i < nContigs; i++) {
        if (contigs[i].isALT) {
            lowestALTContig = __min(lowestALTContig, contigs[i].beginningLocation);
        } else {
            highestNonALTContig = __max(highestNonALTContig, contigs[i].beginningLocation + contigs[i].length);
        }
    }

    if (highestNonALTContig < lowestALTContig && lowestALTContig != LLONG_MAX) {
        WriteErrorMessage("ALT and non-ALT regions overlap.  This is a code bug, it shouldn't happen with any FASTA file.");
        soft_exit(1);
    }

    genomeLocationOfFirstALTContig = lowestALTContig;
}

    bool
Genome::openFileAndGetSizes(const char *filename, GenericFile **file, GenomeDistance *nBases, unsigned *nContigs, bool map)
{
	if (map) {
		*file = GenericFile_map::open(filename);
	} else {
		*file = GenericFile::open(filename, GenericFile::ReadOnly);
	}

    if (*file == NULL) {
        WriteErrorMessage("Genome::openFileAndGetSizes: unable to open file '%s'\n",filename);
        return false;
    } 

    char linebuf[2000];
    char *retval = (*file)->gets(linebuf, sizeof(linebuf));

    unsigned flags = 0;
    if (NULL == retval || 3 != sscanf(linebuf,"%lld %d %d\n", nBases, nContigs, &flags)) {
        (*file)->close();
        delete *file;
        *file = NULL;
        WriteErrorMessage("Genome::openFileAndGetSizes: unable to read header\n");
        return false;
    }

    if (!(flags & GENOME_FLAG_ALT_CONTIGS_MARKED)) {
        WriteErrorMessage("Genome doesn't have the ALT contigs marked flag set.  Either it is corrupt or this is a code bug.  Rebuild the index and if this recurs please submit a bug report.\n");
        soft_exit(1);
    }

    return true;
}


    bool 
Genome::getSizeFromFile(const char *fileName, GenomeDistance *nBases, unsigned *nContigs)
{
    GenericFile *file;
    GenomeDistance localNBases;
    unsigned localnContigs;
    
    if (!openFileAndGetSizes(fileName,&file, nBases ? nBases : &localNBases, nContigs ? nContigs : &localnContigs, false)) {
        return false;
    }

    file->close();
    delete file;
    return true;
}


    bool
Genome::getLocationOfContig(const char *contigName, GenomeLocation *location, InternalContigNum * index) const
{
    if (contigsByName) {
        int low = 0;
        int high = nContigs - 1;
        while (low <= high) {
            int mid = (low + high) / 2;
            int c = strcmp(contigsByName[mid].name, contigName);
            if (c == 0) {
                if (location != NULL) {
                    *location = contigsByName[mid].beginningLocation;
                }
                if (index != NULL) {
                    *index = InternalContigNumToInt(contigsByName[mid].internalContigNumber);
                }
                return true;
            } else if (c < 0) {
                low = mid + 1;
            } else {
                high = mid - 1;
            }
        }
        return false;
    }

    for (int i = 0; i < nContigs; i++) {
        if (!strcmp(contigName,contigs[i].name)) {
            if (NULL != location) {
                *location = contigs[i].beginningLocation;
            }
			if (index != NULL) {
				*index = InternalContigNumToInt(i);
			}
            return true;
        }
    }
    return false;
}

Genome::Contig ContigForInvalidGenomeLocation;

    const Genome::Contig *
Genome::getContigAtLocation(GenomeLocation location) const
{
    if (location == InvalidGenomeLocation) {
        return &ContigForInvalidGenomeLocation;
    }
    _ASSERT(location < nBases);
    int low = 0;
    int high = nContigs - 1;
    while (low <= high) {
        int mid = (low + high) / 2;
        if (contigs[mid].beginningLocation <= location &&
                (mid == nContigs-1 || contigs[mid+1].beginningLocation > location)) {
            return &contigs[mid];
        } else if (contigs[mid].beginningLocation <= location) {
            low = mid + 1;
        } else {
            high = mid - 1;
        }
    }
    return NULL; // Should not be reached
}

    int
Genome::getContigNumAtLocation(GenomeLocation location) const
{
    const Contig *contig = getContigAtLocation(location);
    return (int)(contig - contigs);
}

    const Genome::Contig *
Genome::getNextContigAfterLocation(GenomeLocation location) const
{
    _ASSERT(location < nBases);
    if (nContigs > 0 && location < contigs[0].beginningLocation) {
            return &contigs[0];
    }

    int low = 0;
    int high = nContigs - 1;
    while (low <= high) {
        int mid = (low + high) / 2;
        if (contigs[mid].beginningLocation <= location &&
                (mid == nContigs-1 || contigs[mid+1].beginningLocation > location)) {
            if (mid >= nContigs - 1) {
                //
                // This location landed in the last contig, so return NULL for the next one.
                //
                return NULL;
            } else {
                return &contigs[mid+1];
            }
        } else if (contigs[mid].beginningLocation <= location) {
            low = mid + 1;
        } else {
            high = mid - 1;
        }
    }
    return NULL; // Should not be reached
}

int Genome::getProjContigNumAtLocation(GenomeLocation location) const {
    int contig = getContigNumAtLocation(location);
    GenomeLocation projLocation = contigs[contig].projBeginningLocation;
    return getContigNumAtLocation(projLocation);
}

GenomeDistance Genome::getContigProjSpan(int contig) const {
    GenomeDistance offset = 0;
    ProjCigarOpFormat* cigar = contigs[contig].projCigarOps;
    int nCigarOps = contigs[contig].countProjCigarOps;
    if (nCigarOps > 0) {
        char action = cigar[0].action;
        int count = cigar[0].count;
        if (action != 'S' && action != 'H') {
            offset += count;
        }
    }
    for (int i = 1; i < nCigarOps; ++i) {
        char action = cigar[i].action;
        int count = cigar[i].count;
        if (action == 'M' || action == 'I') {
            offset += count;
        }
    }
    return offset;
}

GenomeLocation Genome::getProjLocation(GenomeLocation location, int alnSpan) const {
    int contig = getContigNumAtLocation(location);
    bool isProjRC = contigs[contig].isProjRC;
    GenomeDistance offset = !isProjRC ? location - contigs[contig].beginningLocation : getContigProjSpan(contig) - (location - contigs[contig].beginningLocation + alnSpan);
    GenomeDistance projOffset = 0;
    ProjCigarOpFormat* cigar = contigs[contig].projCigarOps;
    int cigarOpCount = contigs[contig].countProjCigarOps;
    int start = 0;
    if (cigarOpCount > 0) {
        char action = cigar[0].action;
        if (action == 'S' || action == 'H') {
            start = 1;
        }
    }
    // From bwa-postalt.js: https://github.com/lh3/bwa/blob/master/bwakit/bwa-postalt.js
    GenomeDistance x = 0, y = 0; // x - REF, y - ALT
    for (int i = start; i < cigarOpCount; ++i) {
        char action = cigar[i].action;
        int count = cigar[i].count;
        if (action == 'M') {
            if (y <= offset && offset < y + count) {
                projOffset = x + (offset - y);
                break;
            }
            x += count;
            y += count;
        }
        else if (action == 'D') {
            x += count;
        }
        else if (action == 'I') {
            if (y <= offset && offset < y + count) {
                projOffset = x;
                break;
        }
            y += count;
        }
        else if (action == 'S' || action == 'H') {
            if (y <= offset && offset < y + count) {
                projOffset = -1;
                break;
            }
            y += count;
        }
    }
    if (projOffset == -1) {
        return InvalidGenomeLocation;
    }
    return contigs[contig].projBeginningLocation + projOffset;
}

bool Genome::isProjContigRC(GenomeLocation location) const {
    int contig = getContigNumAtLocation(location);
    return contigs[contig].isProjRC;
}

GenomeDistance DistanceBetweenGenomeLocations(GenomeLocation locationA, GenomeLocation locationB) 
{
    if (locationA > locationB) return locationA - locationB;
    return locationB - locationA;
}

void Genome::fillInContigLengths()
{
    if (nContigs == 0) return;

    for (int i = 0; i < nContigs - 1; i++) {
        contigs[i].length =  contigs[i+1].beginningLocation - contigs[i].beginningLocation;
    }

    contigs[nContigs-1].length = nBases - GenomeLocationAsInt64(contigs[nContigs-1].beginningLocation);
}

const Genome::Contig *Genome::getContigForRead(GenomeLocation location, unsigned readLength, GenomeDistance *extraBasesClippedBefore) const 
{
    const Contig *contig = getContigAtLocation(location);

    //
    // Sometimes, a read aligns before the beginning of a chromosome (imagine prepending a few bases to the read).
    // In that case, we want to handle it by soft-clipping the bases off of the beginning of the read.  We detect it
    // here by looking to see if the aligned location plus the read length crosses a contig boundary.  It also might
    // happen that it is aligned before the first contig, in which case contig will be NULL.
    //
     if (NULL == contig || location + readLength > contig->beginningLocation + contig->length) {
        //
        // We should never align over the end of a chromosome, only before the beginning.  So move this into the next
        // chromosome.
        //
        contig = getNextContigAfterLocation(location);
        _ASSERT(NULL != contig);
        _ASSERT(contig->beginningLocation > location && contig->beginningLocation < location + readLength);
        *extraBasesClippedBefore = contig->beginningLocation - location;
    } else {
        *extraBasesClippedBefore = 0;
    }

    return contig;
}

int 
Genome::getNumALTContigs() const 
{
	int nALTContigs = 0;

	for (int i = 0; i < nContigs; i++) {
		if (contigs[i].isALT) {
			nALTContigs++;
		}
	}

	return nALTContigs;
}

    char* 
Genome::genomeLocationInStringForm(GenomeLocation location, char* buffer, size_t bufferSize) const
{
    const Contig* contig = getContigAtLocation(location);

    const char* contigName;
    GenomeDistance offsetInContig;

    if (contig == NULL) {
        contigName = "UnknownContig";
        offsetInContig = GenomeLocationAsInt64(location);
    } else {
        contigName = contig->name;
        offsetInContig = location - contig->beginningLocation;
    }

    snprintf(buffer, bufferSize, "%s:", contigName);
    size_t usedSize = strlen(buffer);

    const size_t uintBufferSize = 200;
    char uintBuffer[uintBufferSize];

    snprintf(buffer + usedSize, bufferSize - usedSize, "%s", FormatUIntWithCommas(offsetInContig, uintBuffer, uintBufferSize));

    return buffer;
} // genomeLocationInStringForm


GenomeLocation InvalidGenomeLocation;   // Gets set on genome build/load
