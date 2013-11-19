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
#include "Compat.h"
#include "BigAlloc.h"
#include "exit.h"

Genome::Genome(unsigned i_maxBases, unsigned nBasesStored, unsigned i_chromosomePadding)
    : maxBases(i_maxBases), minOffset(0), maxOffset(i_maxBases), chromosomePadding(i_chromosomePadding)
{
    bases = ((char *) BigAlloc(nBasesStored + 2 * N_PADDING)) + N_PADDING;
    if (NULL == bases) {
        fprintf(stderr,"Genome: unable to allocate memory for %llu bases\n",(_int64)maxBases);
        soft_exit(1);
    }

    // Add N's for the N_PADDING bases before and after the genome itself
    memset(bases - N_PADDING, 'n', N_PADDING);
    memset(bases + nBasesStored, 'n', N_PADDING);

    nBases = 0;

    maxContigs = 32; // A power of two that's bigger than the usual number of chromosomes, so we don't have to
                    // reallocate in practice.

    nContigs = 0;
    contigs = new Contig[maxContigs];
    contigsByName = NULL;
}

    void
Genome::addData(const char *data, size_t len)
{
    if ((size_t)nBases + len > maxBases) {
        fprintf(stderr,"Tried to write beyond allocated genome size (or tried to write into a genome that was loaded from a file).\n");
        fprintf(stderr,"Size = %lld\n",(_int64)maxBases);
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
            fprintf(stderr,"Genome: unable to reallocate contig array to size %d\n",newMaxContigs);
            soft_exit(1);
        }
        for (int i = 0; i < nContigs; i++) {
            newContigs[i] = contigs[i];
        }

        delete [] contigs;
        contigs = newContigs;
        maxContigs = newMaxContigs;
    }

    contigs[nContigs].beginningOffset = nBases;
    size_t len = strlen(contigName) + 1;
    contigs[nContigs].name = new char[len];
    if (NULL == contigs[nContigs].name) {
        fprintf(stderr,"Unable to allocate space for contig name\n");
        soft_exit(1);
    }

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
        fprintf(stderr,"Genome::saveToFile: unable to open file '%s'\n",fileName);
        return false;
    } 

    fprintf(saveFile,"%d %d\n",nBases,nContigs);
    char *curChar = NULL;

    for (int i = 0; i < nContigs; i++) {
        for (int n = 0; n < strlen(contigs[i].name); n++){
         curChar = contigs[i].name + n;
         if (*curChar == ' '){ *curChar = '_'; }
        }
        fprintf(saveFile,"%d %s\n",contigs[i].beginningOffset,contigs[i].name);
    }

    if (nBases != fwrite(bases,1,nBases,saveFile)) {
        fprintf(stderr,"Genome::saveToFile: fwrite failed\n");
        fclose(saveFile);
        return false;
    }

    fclose(saveFile);
    return true;
}

    const Genome *
Genome::loadFromFile(const char *fileName, unsigned chromosomePadding, unsigned i_minOffset, unsigned length)
{    
    FILE *loadFile;
    unsigned nBases,nContigs;

    if (!openFileAndGetSizes(fileName,&loadFile,&nBases,&nContigs)) {
        //
        // It already printed an error.  Just fail.
        //
        return NULL;
    }

    if (0 == length) {
        length = nBases - i_minOffset;
    } else {
        //
        // Don't let length go beyond nBases.
        //
        length = __min(length,nBases - i_minOffset);
    }

    Genome *genome = new Genome(nBases,length, chromosomePadding);
   
    genome->nBases = nBases;
    genome->nContigs = genome->maxContigs = nContigs;
    genome->contigs = new Contig[nContigs];
    genome->minOffset = i_minOffset;
    if (i_minOffset >= nBases) {
        fprintf(stderr,"Genome::loadFromFile: specified minOffset %u >= nBases %u\n",i_minOffset,nBases);
    }

 

    genome->maxOffset = i_minOffset + length;

    static const unsigned contigNameBufferSize = 512;
    char contigNameBuffer[contigNameBufferSize];
    unsigned n;
    size_t contigSize;
    char *curName;
    for (unsigned i = 0; i < nContigs; i++) {
        if (NULL == fgets(contigNameBuffer, contigNameBufferSize, loadFile)){
	  
	  fprintf(stderr,"Unable to read contig description\n");
            delete genome;
            return NULL;
        }

	for (n = 0; n < contigNameBufferSize; n++){
	  if (contigNameBuffer[n] == ' ') {
	    contigNameBuffer[n] = '\0'; 
	    break;
	  }
	}

    genome->contigs[i].beginningOffset = atoi(contigNameBuffer);
	contigNameBuffer[n] = ' '; 
	n++; // increment n so we start copying at the position after the space
	contigSize = strlen(contigNameBuffer + n) - 1; //don't include the final \n
        genome->contigs[i].name = new char[contigSize + 1];
	curName = genome->contigs[i].name;
	for (unsigned pos = 0; pos < contigSize; pos++) {
	  curName[pos] = contigNameBuffer[pos + n];
	}
        curName[contigSize] = '\0';
    }

    //
    // Skip over the miserable \n that gets left in the file.
    //
    /*  char newline;
    if (1 != fread(&newline,1,1,loadFile)) {
        fprintf(stderr,"Genome::loadFromFile: Unable to read expected newline\n");
        delete genome;
        return NULL;
    }

    if (newline != 10) {
        fprintf(stderr,"Genome::loadFromFile: Expected newline to be 0x0a, got 0x%02x\n",newline);
        delete genome;
        return NULL;
    }
    */

    if (0 != _fseek64bit(loadFile,i_minOffset,SEEK_CUR)) {
        fprintf(stderr,"Genome::loadFromFile: _fseek64bit failed\n");
        soft_exit(1);
    }

    if (length != fread(genome->bases,1,length,loadFile)) {
        fprintf(stderr,"Genome::loadFromFile: fread of bases failed\n");
        fclose(loadFile);
        delete genome;
        return NULL;
    }

    fclose(loadFile);
    genome->fillInContigLengths();
    genome->sortContigsByName();
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
Genome::openFileAndGetSizes(const char *filename, FILE **file, unsigned *nBases, unsigned *nContigs)
{
    *file = fopen(filename,"rb");
    if (*file == NULL) {
        fprintf(stderr,"Genome::openFileAndGetSizes: unable to open file '%s'\n",filename);
        return false;
    } 

    if (2 != fscanf(*file,"%d %d\n",nBases,nContigs)) {
        fclose(*file);
        *file = NULL;
        fprintf(stderr,"Genome::openFileAndGetSizes: unable to read header\n");
        return false;
    }
    return true;
}


    bool 
Genome::getSizeFromFile(const char *fileName, unsigned *nBases, unsigned *nContigs)
{
    FILE *file;
    unsigned localNBases, localnContigs;
    
    if (!openFileAndGetSizes(fileName,&file,nBases ? nBases : &localNBases, nContigs ? nContigs : &localnContigs)) {
        return false;
    }

    fclose(file);
    return true;
}


    bool
Genome::getOffsetOfContig(const char *contigName, unsigned *offset, int * index) const
{
    if (contigsByName) {
        int low = 0;
        int high = nContigs - 1;
        while (low <= high) {
            int mid = (low + high) / 2;
            int c = strcmp(contigsByName[mid].name, contigName);
            if (c == 0) {
                if (offset != NULL) {
                    *offset = contigsByName[mid].beginningOffset;
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
            if (NULL != offset) {
                *offset = contigs[i].beginningOffset;
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
Genome::getContigAtLocation(unsigned location) const
{
    _ASSERT(location < nBases);
    int low = 0;
    int high = nContigs - 1;
    while (low <= high) {
        int mid = (low + high) / 2;
        if (contigs[mid].beginningOffset <= location &&
                (mid == nContigs-1 || contigs[mid+1].beginningOffset > location)) {
            return &contigs[mid];
        } else if (contigs[mid].beginningOffset <= location) {
            low = mid + 1;
        } else {
            high = mid - 1;
        }
    }
    return NULL; // Should not be reached
}

    const Genome::Contig *
Genome::getNextContigAfterLocation(unsigned location) const
{
    _ASSERT(location < nBases);
    if (nContigs > 0 && location < contigs[0].beginningOffset) {
            return &contigs[0];
    }

    int low = 0;
    int high = nContigs - 1;
    while (low <= high) {
        int mid = (low + high) / 2;
        if (contigs[mid].beginningOffset <= location &&
                (mid == nContigs-1 || contigs[mid+1].beginningOffset > location)) {
            if (mid >= nContigs - 1) {
                //
                // This location landed in the last contig, so return NULL for the next one.
                //
                return NULL;
            } else {
                return &contigs[mid+1];
            }
        } else if (contigs[mid].beginningOffset <= location) {
            low = mid + 1;
        } else {
            high = mid - 1;
        }
    }
    return NULL; // Should not be reached
}

//
// Makes a copy of a Genome, but with only one of the sex chromosomes.
//
// The fate of the mitochondrion is that of the X chromosome.
//
    Genome *
Genome::copy(bool copyX, bool copyY, bool copyM) const
{
    Genome *newCopy = new Genome(getCountOfBases(),getCountOfBases(), chromosomePadding);

    if (NULL == newCopy) {
        fprintf(stderr,"Genome::copy: failed to allocate space for copy.\n");
        return NULL;
    }

    const Genome::Contig *currentContig = NULL;
    const Genome::Contig *nextContig = getContigAtLocation(0);

    unsigned offsetInReference = 0;
    while (offsetInReference < getCountOfBases()) {
        if (NULL != nextContig && offsetInReference >= nextContig->beginningOffset) {
            //
            // Start of a new contig.  See if we want to skip it.
            //
            currentContig = nextContig;
            nextContig = getNextContigAfterLocation(offsetInReference + 1);
            if ((!copyX && !strcmp(currentContig->name,"chrX")) ||
                (!copyY && !strcmp(currentContig->name,"chrY")) ||
                (!copyM && !strcmp(currentContig->name,"chrM"))) {
                //
                // Yes, skip over this contig.
                //
                nextContig = getNextContigAfterLocation(offsetInReference + 1);
                if (NULL == nextContig) {
                    //
                    // The chromosome that we're skipping was the last one, so we're done.
                    //
                    break;
                } else {
                    offsetInReference = nextContig->beginningOffset;
                    continue;
                }
            } // If skipping this chromosome

            newCopy->startContig(currentContig->name);
        } // If new contig beginning

        const size_t maxCopySize = 10000;
        char dataBuffer[maxCopySize + 1];

        unsigned amountToCopy = maxCopySize;
        if (nextContig && nextContig->beginningOffset < offsetInReference + amountToCopy) {
            amountToCopy = nextContig->beginningOffset - offsetInReference;
        }

        if (getCountOfBases() < offsetInReference + amountToCopy) {
            amountToCopy = getCountOfBases() - offsetInReference;
        }

        memcpy(dataBuffer,getSubstring(offsetInReference,amountToCopy), amountToCopy);
        dataBuffer[amountToCopy] = '\0';

        newCopy->addData(dataBuffer);

        offsetInReference += amountToCopy;
    }

    newCopy->fillInContigLengths();
    newCopy->sortContigsByName();
    return newCopy;
}

unsigned DistanceBetweenGenomeLocations(unsigned locationA, unsigned locationB) 
{
    unsigned largerGenomeOffset = __max(locationA, locationB);
    unsigned smallerGenomeOffset = __min(locationA, locationB);

    return largerGenomeOffset - smallerGenomeOffset;
}

void Genome::fillInContigLengths()
{
    if (nContigs == 0) return;

    for (int i = 0; i < nContigs - 1; i++) {
        contigs[i].length =  contigs[i+1].beginningOffset - contigs[i].beginningOffset;
    }

    contigs[nContigs-1].length = nBases - contigs[nContigs-1].beginningOffset;
}

const Genome::Contig *Genome::getContigForRead(unsigned location, unsigned readLength, unsigned *extraBasesClippedBefore) const 
{
    const Contig *contig = getContigAtLocation(location);

    //
    // Sometimes, a read aligns before the beginning of a chromosome (imagine prepending a few bases to the read).
    // In that case, we want to handle it by soft-clipping the bases off of the beginning of the read.  We detect it
    // here by looking to see if the aligned location plus the read length crosses a contig boundary.  It also might
    // happen that it is aligned before the first contig, in which case contig will be NULL.
    //
     if (NULL == contig || location + readLength > contig->beginningOffset + contig->length) {
        //
        // We should never align over the end of a chromosome, only before the beginning.  So move this into the next
        // chromosome.
        //
        contig = getNextContigAfterLocation(location);
        _ASSERT(NULL != contig);
        _ASSERT(contig->beginningOffset > location && contig->beginningOffset < location + readLength);
        *extraBasesClippedBefore = contig->beginningOffset - location;
    } else {
        *extraBasesClippedBefore = 0;
    }

    return contig;
}