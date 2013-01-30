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

Genome::Genome(unsigned i_maxBases, unsigned nBasesStored) : maxBases(i_maxBases), minOffset(0), maxOffset(i_maxBases)
{
    bases = ((char *) BigAlloc(nBasesStored + 2 * N_PADDING)) + N_PADDING;
    if (NULL == bases) {
        fprintf(stderr,"Genome: unable to allocate memory for %llu bases\n",(_int64)maxBases);
        exit(1);
    }

    // Add N's for the N_PADDING bases before and after the genome itself
    memset(bases - N_PADDING, 'N', N_PADDING);
    memset(bases + nBasesStored, 'N', N_PADDING);

    nBases = 0;

    maxPieces = 32; // A power of two that's bigger than the usual number of chromosomes, so we don't have to
                    // reallocate in practice.

    nPieces = 0;
    pieces = new Piece[maxPieces];
}

    void
Genome::addData(const char *data, size_t len)
{
    if ((size_t)nBases + len > maxBases) {
        fprintf(stderr,"Tried to write beyond allocated genome size (or tried to write into a genome that was loaded from a file).\n");
        fprintf(stderr,"Size = %lld\n",(_int64)maxBases);
        exit(1);
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
Genome::startPiece(const char *pieceName)
{
    if (nPieces == maxPieces) {
        //
        // Reallocate (maybe we're sequencing a tree that's got lots of chromosomes).
        //
        int newMaxPieces = maxPieces * 2;
        Piece *newPieces = new Piece[newMaxPieces];
        if (NULL == newPieces) {
            fprintf(stderr,"Genome: unable to reallocate piece array to size %d\n",newMaxPieces);
            exit(1);
        }
        for (int i = 0; i < nPieces; i++) {
            newPieces[i] = pieces[i];
        }

        delete [] pieces;
        pieces = newPieces;
        maxPieces = newMaxPieces;
    }

    pieces[nPieces].beginningOffset = nBases;
    size_t len = strlen(pieceName) + 1;
    pieces[nPieces].name = new char[len];
    if (NULL == pieces[nPieces].name) {
        fprintf(stderr,"Unable to allocate space for piece name\n");
        exit(1);
    }

    strncpy(pieces[nPieces].name,pieceName,len);
    pieces[nPieces].name[len-1] = '\0';

    nPieces++;
}


Genome::~Genome()
{
    BigDealloc(bases - N_PADDING);
    for (int i = 0; i < nPieces; i++) {
        delete [] pieces[i].name;
        pieces[i].name = NULL;
    }

    delete [] pieces;
    pieces = NULL;
}


    bool
Genome::saveToFile(const char *fileName) const
{
    //
    // Save file format is (in binary) the number of bases, the number of pieces, followed by
    //  the pieces themselves, rounded up to 4K, followed by the bases.
    //

    FILE *saveFile = fopen(fileName,"wb");
    if (saveFile == NULL) {
        fprintf(stderr,"Genome::saveToFile: unable to open file '%s'\n",fileName);
        return false;
    } 

    fprintf(saveFile,"%d %d\n",nBases,nPieces);
    char *curChar = NULL;

    for (int i = 0; i < nPieces; i++) {
        for (int n = 0; n < strlen(pieces[i].name); n++){
         curChar = pieces[i].name + n;
         if (*curChar == ' '){ *curChar = '_'; }
        }
        fprintf(saveFile,"%d %s\n",pieces[i].beginningOffset,pieces[i].name);
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
Genome::loadFromFile(const char *fileName, unsigned i_minOffset, unsigned length)
{    
    FILE *loadFile;
    unsigned nBases,nPieces;

    if (!openFileAndGetSizes(fileName,&loadFile,&nBases,&nPieces)) {
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

    Genome *genome = new Genome(nBases,length);
   
    genome->nBases = nBases;
    genome->nPieces = genome->maxPieces = nPieces;
    genome->pieces = new Piece[nPieces];
    genome->minOffset = i_minOffset;
    if (i_minOffset >= nBases) {
        fprintf(stderr,"Genome::loadFromFile: specified minOffset %u >= nBases %u\n",i_minOffset,nBases);
    }

 

    genome->maxOffset = i_minOffset + length;

    static const unsigned pieceNameBufferSize = 512;
    char pieceNameBuffer[pieceNameBufferSize];
    unsigned n;
    size_t pieceSize;
    char *curName;
    for (unsigned i = 0; i < nPieces; i++) {
        if (NULL == fgets(pieceNameBuffer, pieceNameBufferSize, loadFile)){
	  
	  fprintf(stderr,"Unable to read piece description\n");
            delete genome;
            return NULL;
        }

	for (n = 0; n < pieceNameBufferSize; n++){
	  if (pieceNameBuffer[n] == ' ') {
	    pieceNameBuffer[n] = '\0'; 
	    break;
	  }
	}

    genome->pieces[i].beginningOffset = atoi(pieceNameBuffer);
	pieceNameBuffer[n] = ' '; 
	n++; // increment n so we start copying at the position after the space
	pieceSize = strlen(pieceNameBuffer + n) - 1; //don't include the final \n
        genome->pieces[i].name = new char[pieceSize + 1];
	curName = genome->pieces[i].name;
	for (unsigned pos = 0; pos < pieceSize; pos++) {
	  curName[pos] = pieceNameBuffer[pos + n];
	}
        curName[pieceSize] = '\0';
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
        exit(1);
    }

    if (length != fread(genome->bases,1,length,loadFile)) {
        fprintf(stderr,"Genome::loadFromFile: fread of bases failed\n");
        fclose(loadFile);
        delete genome;
        return NULL;
    }

    fclose(loadFile);
    return genome;
}

    bool
Genome::openFileAndGetSizes(const char *filename, FILE **file, unsigned *nBases, unsigned *nPieces)
{
    *file = fopen(filename,"rb");
    if (*file == NULL) {
        fprintf(stderr,"Genome::openFileAndGetSizes: unable to open file '%s'\n",filename);
        return false;
    } 

    if (2 != fscanf(*file,"%d %d\n",nBases,nPieces)) {
        fclose(*file);
        *file = NULL;
        fprintf(stderr,"Genome::openFileAndGetSizes: unable to read header\n");
        return false;
    }
    return true;
}


    bool 
Genome::getSizeFromFile(const char *fileName, unsigned *nBases, unsigned *nPieces)
{
    FILE *file;
    unsigned localNBases, localNPieces;
    
    if (!openFileAndGetSizes(fileName,&file,nBases ? nBases : &localNBases, nPieces ? nPieces : &localNPieces)) {
        return false;
    }

    fclose(file);
    return true;
}


    bool
Genome::getOffsetOfPiece(const char *pieceName, unsigned *offset) const
{
    for (int i = 0; i < nPieces; i++) {
        if (!strcmp(pieceName,pieces[i].name)) {
            if (NULL != offset) {
                *offset = pieces[i].beginningOffset;
            }
            return true;
        }
    }
    return false;
}


    const Genome::Piece *
Genome::getPieceAtLocation(unsigned location) const
{
    _ASSERT(location < nBases);
    int low = 0;
    int high = nPieces - 1;
    while (low <= high) {
        int mid = (low + high) / 2;
        if (pieces[mid].beginningOffset <= location &&
                (mid == nPieces-1 || pieces[mid+1].beginningOffset > location)) {
            return &pieces[mid];
        } else if (pieces[mid].beginningOffset <= location) {
            low = mid + 1;
        } else {
            high = mid - 1;
        }
    }
    return NULL; // Should not be reached
}

    const Genome::Piece *
Genome::getNextPieceAfterLocation(unsigned location) const
{
    _ASSERT(location < nBases);
    int low = 0;
    int high = nPieces - 1;
    while (low <= high) {
        int mid = (low + high) / 2;
        if (pieces[mid].beginningOffset <= location &&
                (mid == nPieces-1 || pieces[mid+1].beginningOffset > location)) {
            if (mid >= nPieces - 1) {
                //
                // This location landed in the last piece, so return NULL for the next one.
                //
                return NULL;
            } else {
                return &pieces[mid+1];
            }
        } else if (pieces[mid].beginningOffset <= location) {
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
    Genome *newCopy = new Genome(getCountOfBases(),getCountOfBases());

    if (NULL == newCopy) {
        fprintf(stderr,"Genome::copy: failed to allocate space for copy.\n");
        return NULL;
    }

    const Genome::Piece *currentPiece = NULL;
    const Genome::Piece *nextPiece = getPieceAtLocation(0);

    unsigned offsetInReference = 0;
    while (offsetInReference < getCountOfBases()) {
        if (NULL != nextPiece && offsetInReference >= nextPiece->beginningOffset) {
            //
            // Start of a new piece.  See if we want to skip it.
            //
            currentPiece = nextPiece;
            nextPiece = getNextPieceAfterLocation(offsetInReference + 1);
            if ((!copyX && !strcmp(currentPiece->name,"chrX")) ||
                (!copyY && !strcmp(currentPiece->name,"chrY")) ||
                (!copyM && !strcmp(currentPiece->name,"chrM"))) {
                //
                // Yes, skip over this piece.
                //
                nextPiece = getNextPieceAfterLocation(offsetInReference + 1);
                if (NULL == nextPiece) {
                    //
                    // The chromosome that we're skipping was the last one, so we're done.
                    //
                    break;
                } else {
                    offsetInReference = nextPiece->beginningOffset;
                    continue;
                }
            } // If skipping this chromosome

            newCopy->startPiece(currentPiece->name);
        } // If new piece beginning

        const size_t maxCopySize = 10000;
        char dataBuffer[maxCopySize + 1];

        unsigned amountToCopy = maxCopySize;
        if (nextPiece && nextPiece->beginningOffset < offsetInReference + amountToCopy) {
            amountToCopy = nextPiece->beginningOffset - offsetInReference;
        }

        if (getCountOfBases() < offsetInReference + amountToCopy) {
            amountToCopy = getCountOfBases() - offsetInReference;
        }

        memcpy(dataBuffer,getSubstring(offsetInReference,amountToCopy), amountToCopy);
        dataBuffer[amountToCopy] = '\0';

        newCopy->addData(dataBuffer);

        offsetInReference += amountToCopy;
    }

    return newCopy;
}


    DiploidGenome *
DiploidGenome::Factory(const Genome *referenceGenome, bool isMale)
{
    //
    // This is just the special case where we use the same reference genome for both parents.
    //
    return Factory(referenceGenome, referenceGenome, isMale);
}

    DiploidGenome *
DiploidGenome::Factory(const Genome *motherReference, const Genome *fatherReference, bool isMale)
{
    //
    // We have to create parent genomes that have the right set of chromosomes.  Make ParentA as the mother (i.e., include chrX and not
    // chrY).
    //
    Genome *mom = motherReference->copyGenomeOneSex(false, true);   // Always get an X chromosome from mom
    if (NULL == mom) {
        fprintf(stderr,"DiploidGenome::Factory: Unable to allocate mom's genome.\n");
        return NULL;
    }

    Genome *dad = fatherReference->copyGenomeOneSex(isMale, false);  // Get dad's Y if this is a male, otherwise dad's X
    if (NULL == dad) {
        fprintf(stderr,"DiploidGenome::Factory: Unable to allocate dad's genome.\n");
        delete mom;
        
        return NULL;
    }

    DiploidGenome *diploidGenome = new DiploidGenome(mom, dad, true);
    if (NULL == diploidGenome) {
        fprintf(stderr,"DiploidGenome::Factory: Unable to allocate container.\n");  // Sometimes you write code that you just know will never execute.
        delete mom;
        delete dad;
        return NULL;
    }

    return diploidGenome;
}



    DiploidGenome *
DiploidGenome::Factory(const Genome *motherReference, const Genome *fatherReference)
{
    return new DiploidGenome(motherReference, fatherReference, false);
}

    DiploidGenome *
DiploidGenome::CopyingFactory(const Genome *motherReference, const Genome *fatherReference)
{
    Genome *motherCopy = motherReference->copy();
    if (NULL == motherCopy) {
        fprintf(stderr,"DiploidGenome: unable to copy mother's genome\n");
        return NULL;
    }

    Genome *fatherCopy = fatherReference->copy();
    if (NULL == fatherCopy) {
        fprintf(stderr,"DiploidGenome: unable to copy father's genome\n");
        delete motherCopy;
        return NULL;
    }

    DiploidGenome *diploidGenome = new DiploidGenome(motherCopy, fatherCopy, true);
    if (NULL == diploidGenome) {
        fprintf(stderr,"DiploidGenome: unable to allocate container.\n");
        delete motherCopy;
        delete fatherCopy;
        return NULL;
    }

    return diploidGenome;
}

    const char *
DiploidGenome::FilenameBase = "DiploidParent";


    DiploidGenome *
DiploidGenome::loadFromDirectory(const char *directoryName)
{
    size_t fileNameLen = strlen(directoryName) + 1 + strlen(FilenameBase) + 2;  // +1 is PATH_SEP, +2 is A/B and \0
    char *filename = new char[fileNameLen];

    if (NULL == filename) {
        fprintf(stderr,"DiploidGenome::loadFromDirectory: failed to allocate filename buffer.\n");
        return NULL;
    }

    sprintf(filename,"%s%c%s%c",directoryName,PATH_SEP,FilenameBase,'A');
    const Genome *motherGenome = Genome::loadFromFile(filename);

    if (NULL == motherGenome) {
        fprintf(stderr,"DiploidGenome::loadFromDirectory: failed to load mother's genome.\n");
        delete [] filename;
        return NULL;
    }

    sprintf(filename,"%s%c%s%c",directoryName,PATH_SEP,FilenameBase,'B');
    const Genome *fatherGenome = Genome::loadFromFile(filename);

    if (NULL == fatherGenome) {
        fprintf(stderr,"DiploidGenome::loadFromDirectory: failed to load father's genome.\n");
        delete motherGenome;
        delete [] filename;
        return NULL;
    }

    delete [] filename;

    DiploidGenome *diploidGenome = new DiploidGenome(motherGenome, fatherGenome, true);
    if (NULL == diploidGenome) {
        fprintf(stderr,"DiploidGenome::loadFromDirectory: failed to allocate DiploidGenome container.\n");
        delete motherGenome;
        delete fatherGenome;
        return NULL;
    }

    return diploidGenome;
}

    bool
DiploidGenome::saveToDirectory(const char *directoryName) const
{

    size_t fileNameLen = strlen(directoryName) + 1 + strlen(FilenameBase) + 2;  // +1 is PATH_SEP, +2 is A/B and \0
    char *filename = new char[fileNameLen];

    if (NULL == filename) {
        fprintf(stderr,"DiploidGenome::saveToDirectory: failed to allocate filename buffer.\n");
        return false;
    }

    sprintf(filename,"%s%c%s%c",directoryName,PATH_SEP,FilenameBase,'A');
    if (!motherGenome->saveToFile(filename)) {
        fprintf(stderr,"DiploidGenome::saveToDirectory: failed to save mother's genome.\n");
        delete [] filename;
        return false;
    }

    sprintf(filename,"%s%c%s%c",directoryName,PATH_SEP,FilenameBase,'B');
    if (!fatherGenome->saveToFile(filename)) {
        fprintf(stderr,"DiploidGenome::saveToDirectory: failed to save father's genome.\n");
        delete motherGenome;
        delete [] filename;
        return false;
    }

    delete [] filename;
    return true;
}

    const Genome *
DiploidGenome::getGenome(bool fromMother) const
{
    if (fromMother) return motherGenome;
    return fatherGenome;
}

DiploidGenome::~DiploidGenome()
{
    if (ownsGenomes) {
        delete motherGenome;
        delete fatherGenome;
    }
}
