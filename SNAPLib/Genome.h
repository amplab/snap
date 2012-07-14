/*++

Module Name:

    geonome.h

Abstract:

    Genome class for the SNAP sequencer

Authors:

    Bill Bolosky, August, 2011

Environment:

    User mode service.

Revision History:

    Adapted from Matei Zaharia's Scala implementation.

--*/

#pragma once
#include <xmmintrin.h>
#include "Compat.h"

class Genome {
public:
        //
        // Methods for building a genome.
        //

        //
        // Create a new genome.  It's got a limit on the number of bases, but there's no need to
        // store that many, it's just an upper bound.  It will, of course, use memory proportional
        // to the bound.
        //
        Genome(
            unsigned             i_maxBases,
            unsigned             nBasesStored);

        void startPiece(
            const char          *pieceName);

        void addData(
            const char          *data);

        void addData(const char *data, size_t len);


        ~Genome();

        //
        // Methods for loading a genome from a file, and saving one to a file.  When you save and
        // then load a genome the space used is reduced from the max that was reserved when it was
        // first created to the amount actually used (rounded up to a page size).  However, saved
        // and loaded genomes can't be added to, they're read only.
        //
        // minOffset and length are used to read in only a part of a whole genome.
        //
        static const Genome *loadFromFile(const char *fileName, unsigned i_minOffset = 0, unsigned length = 0);  
                                                                  // This loads from a genome save
                                                                  // file, not a FASTA file.  Use
                                                                  // FASTA.h for FASTA loads.

        static bool getSizeFromFile(const char *fileName, unsigned *nBases, unsigned *nPieces);

        bool saveToFile(const char *fileName) const;

        //
        // Methods to read the genome.
        //
        inline const char *getSubstring(size_t offset, size_t lengthNeeded) const {
            // TODO: Don't return substrings that span piece (chromosome) boundaries
            if (offset + lengthNeeded > nBases) {
                return NULL;
            }

            _ASSERT(offset >= minOffset && offset + lengthNeeded <= maxOffset); // If the caller asks for a genome slice, it's only legal to look within it.

            return bases + (offset-minOffset);
        }

        inline unsigned getCountOfBases() const {return nBases;}

        bool getOffsetOfPiece(const char *pieceName, unsigned *offset) const;

        inline void prefetchData(unsigned genomeOffset) const {
            _mm_prefetch(bases + genomeOffset,_MM_HINT_T2);
            _mm_prefetch(bases + genomeOffset + 64,_MM_HINT_T2);
        }

        struct Piece {
            unsigned     beginningOffset;
            char        *name;
        };

        inline const Piece *getPieces() const { return pieces; }

        inline int getNumPieces() const { return nPieces; }

        const Piece *getPieceAtLocation(unsigned location) const;
        const Piece *getNextPieceAfterLocation(unsigned location) const;

        Genome *copy() const {return copy(true,true,true);}
        Genome *copyGenomeOneSex(bool useY, bool useM) const {return copy(!useY,useY,useM);}
private:

        //
        // The actual genome.
        char        *bases;
        unsigned     nBases;
        unsigned     maxBases;

        unsigned     minOffset;
        unsigned     maxOffset;

        //
        // A genome is made up of a bunch of pieces, typically chromosomes.  Pieces have names,
        // which are stored here.
        //
        int          nPieces;
        int          maxPieces;

        Piece       *pieces;

        Genome *copy(bool copyX, bool copyY, bool copyM) const;

        static bool openFileAndGetSizes(const char *filename, FILE **file, unsigned *nBases, unsigned *nPieces);

};

class DiploidGenome {
public:
        static DiploidGenome *Factory(const Genome *referenceGenome, bool isMale);
        static DiploidGenome *Factory(const Genome *motherReference, const Genome *fatherReference,
                                      bool isMale);

        //
        // This version of Factory just takes references to the parent Genomes, and it's
        // the caller's responsibility to make sure they continue to exist.
        //        
        static DiploidGenome *Factory(const Genome *motherReference, const Genome *fatherReference);

        //
        // The copying factory makes copies of the underlying genomes.
        //
        static DiploidGenome *CopyingFactory(const Genome *motherReference,
                                             const Genome *fatherReference);

        static DiploidGenome *loadFromDirectory(const char *directoryName);


        bool saveToDirectory(const char *directoryName) const;
        const Genome *getGenome(bool fromMother) const;
        const unsigned getCountOfBases(bool fromMother) const {
            return fromMother ? motherGenome->getCountOfBases() : fatherGenome->getCountOfBases();
        }
        const _int64 getCountOfBasesBothHalves() const {
            return (_int64)motherGenome->getCountOfBases() +
                (_int64)fatherGenome->getCountOfBases();
        }

        ~DiploidGenome();
private:

        DiploidGenome(const Genome *i_motherGenome, const Genome *i_fatherGenome,
                      bool i_ownsGenomes) : 
            motherGenome(i_motherGenome), fatherGenome(i_fatherGenome),
            ownsGenomes(i_ownsGenomes) {} 

        const Genome *motherGenome;
        const Genome *fatherGenome;

        bool ownsGenomes;

        static Genome *copyGenomeOneSex(const Genome *referenceGenome, bool useY);

        static const char *FilenameBase;
};
