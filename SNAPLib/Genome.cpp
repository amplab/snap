/*++

Module Name:

    geonome.cpp

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
: maxBases(i_maxBases), minLocation(0), maxLocation(i_maxBases), chromosomePadding(i_chromosomePadding), maxContigs(i_maxContigs),
  mappedFile(NULL)
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
    contigs = new Contig[maxContigs];
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
Genome::startContig(const char *contigName)
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

    nContigs++;
}


Genome::~Genome()
{
    BigDealloc(bases - N_PADDING);
    for (int i = 0; i < nContigs; i++) {
        delete [] contigs[i].name;
        contigs[i].name = NULL;
    }

    delete [] contigs;
    if (contigsByName) {
        delete [] contigsByName;
    }
    contigs = NULL;

	if (NULL != mappedFile) {
		mappedFile->close();
		delete mappedFile;
	}
}


    bool
Genome::saveToFile(const char *fileName) const
{
    //
    // Save file format is (in binary) the number of bases, the number of contigs, followed by
    //  the contigs themselves, rounded up to 4K, followed by the bases.
    //

    FILE *saveFile = fopen(fileName,"wb");
    if (saveFile == NULL) {
        WriteErrorMessage("Genome::saveToFile: unable to open file '%s'\n",fileName);
        return false;
    } 

    fprintf(saveFile,"%lld %d\n",nBases, nContigs);
    char *curChar = NULL;

    for (int i = 0; i < nContigs; i++) {
        for (int n = 0; n < strlen(contigs[i].name); n++){
         curChar = contigs[i].name + n;
         if (*curChar == ' '){ *curChar = '_'; }
        }
        fprintf(saveFile,"%lld %s\n",contigs[i].beginningLocation, contigs[i].name);
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

    Genome *genome = new Genome(nBases, length, chromosomePadding);
   
    genome->nBases = nBases;
    genome->nContigs = genome->maxContigs = nContigs;
    genome->contigs = new Contig[nContigs];
    genome->minLocation = minLocation;
    if (GenomeLocationAsInt64(minLocation) >= nBases) {
        WriteErrorMessage("Genome::loadFromFile: specified minOffset %u >= nBases %u\n", GenomeLocationAsInt64(minLocation), nBases);
        soft_exit(-1);
    }

    genome->maxLocation = maxLocation;

    int contigNameBufferSize = 0;
    char *contigNameBuffer = NULL;
    unsigned n;
    size_t contigSize;
    char *curName;
    for (unsigned i = 0; i < nContigs; i++) {
        if (NULL == reallocatingFgetsGenericFile(&contigNameBuffer, &contigNameBufferSize, loadFile)) {	 
            WriteErrorMessage("Unable to read contig description\n");
            delete[] contigNameBuffer;
            delete genome;
            return NULL;
        }

        for (n = 0; n < (unsigned)contigNameBufferSize; n++) {
	        if (contigNameBuffer[n] == ' ') {
	            contigNameBuffer[n] = '\0'; 
	            break;
	        }
	    }

        _int64 contigStart;
        if (1 != sscanf(contigNameBuffer, "%lld", &contigStart)) {
            WriteErrorMessage("Unable to parse contig start in genome file '%s', '%s%'\n", fileName, contigNameBuffer);
            soft_exit(1);
        }
        genome->contigs[i].beginningLocation = GenomeLocation(contigStart);
	    contigNameBuffer[n] = ' '; 
	    n++; // increment n so we start copying at the position after the space
	    contigSize = strlen(contigNameBuffer + n) - 1; //don't include the final \n
        genome->contigs[i].name = new char[contigSize + 1];
        genome->contigs[i].nameLength = (unsigned)contigSize;
	    curName = genome->contigs[i].name;
	    for (unsigned pos = 0; pos < contigSize; pos++) {
	      curName[pos] = contigNameBuffer[pos + n];
	    }
        curName[contigSize] = '\0';
    } // for each contig

    if (0 != loadFile->advance(GenomeLocationAsInt64(minLocation))) {
        WriteErrorMessage("Genome::loadFromFile: _fseek64bit failed\n");
        soft_exit(1);
    }

    size_t readSize;
	if (map) {
		GenericFile_map *mappedFile = (GenericFile_map *)loadFile;
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

    if (NULL == retval || 2 != sscanf(linebuf,"%lld %d\n", nBases, nContigs)) {
        (*file)->close();
        delete *file;
        *file = NULL;
        WriteErrorMessage("Genome::openFileAndGetSizes: unable to read header\n");
        return false;
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
Genome::getLocationOfContig(const char *contigName, GenomeLocation *location, int * index) const
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
                    *index = mid;
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
				*index = i;
			}
            return true;
        }
    }
    return false;
}


    const Genome::Contig *
Genome::getContigAtLocation(GenomeLocation location) const
{
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

GenomeLocation InvalidGenomeLocation;   // Gets set on genome build/load