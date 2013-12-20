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

const unsigned InvalidGenomeLocation = 0xffffffff;

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
            unsigned             nBasesStored,
            unsigned             i_chromosomePadding);

        void startContig(
            const char          *contigName);

        void addData(
            const char          *data);

        void addData(const char *data, size_t len);

        const unsigned getChromosomePadding() const {return chromosomePadding;}

        ~Genome();

        //
        // Methods for loading a genome from a file, and saving one to a file.  When you save and
        // then load a genome the space used is reduced from the max that was reserved when it was
        // first created to the amount actually used (rounded up to a page size).  However, saved
        // and loaded genomes can't be added to, they're read only.
        //
        // minOffset and length are used to read in only a part of a whole genome.
        //
        static const Genome *loadFromFile(const char *fileName, unsigned chromosomePadding, unsigned i_minOffset = 0, unsigned length = 0);  
                                                                  // This loads from a genome save
                                                                  // file, not a FASTA file.  Use
                                                                  // FASTA.h for FASTA loads.

        static bool getSizeFromFile(const char *fileName, unsigned *nBases, unsigned *nContigs);

        bool saveToFile(const char *fileName) const;

        //
        // Methods to read the genome.
        //
        inline const char *getSubstring(size_t offset, size_t lengthNeeded) const {
            if (offset > nBases || offset + lengthNeeded > nBases + N_PADDING) {
                // The first part of the test is for the unsigned version of a negative offset.
                return NULL;
            }

            if (lengthNeeded <= chromosomePadding) {
                return bases + (offset - minOffset);
            }

            _ASSERT(offset >= minOffset && offset + lengthNeeded <= maxOffset + N_PADDING); // If the caller asks for a genome slice, it's only legal to look within it.

            if (lengthNeeded == 0) {
                return bases + (offset - minOffset);
            }

            //
            // See if the substring crosses a contig (chromosome) boundary.  If so, disallow it.
            //

            if (nContigs > 100) {
                //
                // Start by special casing the last contig (it makes the rest of the code easier).
                //
                if (contigs[nContigs - 1].beginningOffset <= offset) {
                    //
                    // Because it starts in the last contig, it's OK because we already checked overflow
                // of the whole genome.
                //
                return bases + (offset-minOffset);
            }
    
                int min = 0;
                int max = nContigs - 2;
                while (min <= max) {
                    int i = (min + max) / 2;
                    if (contigs[i].beginningOffset <= offset) {
                        if (contigs[i+1].beginningOffset > offset) {
                            if (contigs[i+1].beginningOffset <= offset + lengthNeeded - 1) {
                                return NULL;    // This crosses a contig boundary.
                            } else {
                                return bases + (offset-minOffset);
                            }
                        } else {
                            min = i+1;
                    }
                } else {
                        max = i-1;
                    }
                }
    
                _ASSERT(false && "NOTREACHED");
                return NULL;
            } else {
                //
                // Use linear rather than binary search for small numbers of contigs, because binary search
                // confuses the branch predictor, and so is slower even though it uses many fewer instructions.
                //
                for (int i = 0 ; i < nContigs; i++) {
                    if (offset + lengthNeeded - 1 >= contigs[i].beginningOffset) {
                        if (offset < contigs[i].beginningOffset) {
                            return NULL;        // crosses a contig boundary.
                        } else {
                            return bases + (offset-minOffset);
                        }
                    }
                }
                _ASSERT(false && "NOTREACHED");
                return NULL;
            }
        } 

        inline unsigned getCountOfBases() const {return nBases;}

        bool getOffsetOfContig(const char *contigName, unsigned *offset, int* index = NULL) const;

        inline void prefetchData(unsigned genomeOffset) const {
            _mm_prefetch(bases + genomeOffset,_MM_HINT_T2);
            _mm_prefetch(bases + genomeOffset + 64,_MM_HINT_T2);
        }

        struct Contig {
            unsigned     beginningOffset;
            unsigned     length;
            char        *name;
        };

        inline const Contig *getContigs() const { return contigs; }

        inline int getNumContigs() const { return nContigs; }

        const Contig *getContigAtLocation(unsigned location) const;
        const Contig *getContigForRead(unsigned location, unsigned readLength, unsigned *extraBasesClippedBefore) const;
        const Contig *getNextContigAfterLocation(unsigned location) const;

        Genome *copy() const {return copy(true,true,true);}
        Genome *copyGenomeOneSex(bool useY, bool useM) const {return copy(!useY,useY,useM);}

        //
        // These are only public so creators of new genomes (i.e., FASTA) can use them.
        // 
        void    fillInContigLengths();
        void    sortContigsByName();

private:

        static const int N_PADDING = 100; // Padding to add on either end of the genome to allow substring reads past it

        //
        // The actual genome.
        char        *bases;       // Will point to offset N_PADDING in an array of nBases + 2 * N_PADDING
        unsigned     nBases;
        unsigned     maxBases;

        unsigned     minOffset;
        unsigned     maxOffset;

        //
        // A genome is made up of a bunch of contigs, typically chromosomes.  Contigs have names,
        // which are stored here.
        //
        int          nContigs;
        int          maxContigs;

        Contig      *contigs;    // This is always in order (it's not possible to express it otherwise in FASTA).

        Contig      *contigsByName;
        Genome *copy(bool copyX, bool copyY, bool copyM) const;

        static bool openFileAndGetSizes(const char *filename, FILE **file, unsigned *nBases, unsigned *nContigs);


        const unsigned chromosomePadding;
};

unsigned DistanceBetweenGenomeLocations(unsigned locationA, unsigned locationB);
