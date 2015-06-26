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
#include "Compat.h"
#include "GenericFile.h"
#include "GenericFile_map.h"

//
// We have two different classes to represent a place in a genome and a distance between places in a genome.
// In reality, they're both just 64 bit ints, but the classes are set up to encourage the user to keep
// in mind the difference.  So, a genome location might be something
// like "chromosome 12, base 12345" which would be represented in (0-based) genome coordinates as some 
// 64 bit int that's the base of cheomosome 12 + 12344 (one less because we're 0-based and the nomenclature
// uses 1-based).
// In the non-debug build, GenomeLocation is just defined as an _int64, so that no matter how dumb the compiler
// can't possibly screw it up.  However, in the debug build you get the happy type checking.
//

typedef _int64 GenomeDistance;

#ifdef _DEBUG

class GenomeLocation {
public:
    GenomeLocation(_int64 i_location) : location(i_location) {}
    GenomeLocation(const GenomeLocation &peer) : location(peer.location) {}
    GenomeLocation() {location = -1;}

    inline GenomeLocation operator=(const GenomeLocation &peer) {
        location = peer.location;
        return *this;
    }

    inline GenomeLocation operator=(const _int64 value) {
       location = value;
       return *this;
    }

    inline GenomeLocation operator++() {
        location++;
        return *this;
    }

    inline GenomeLocation operator--() {
        location--;
        return *this;
    }

    // The postfix versions
    inline GenomeLocation operator++(int foo) {
        location++;
        return *this - 1;
    }

    inline GenomeLocation operator--(int foo) {
        location--;
        return *this + 1;
    }

    inline bool operator==(const GenomeLocation &peer) const {
        return location == peer.location;
    }

    inline bool operator>=(const GenomeLocation &peer) const {
        return location >= peer.location;
    }

    inline bool operator>(const GenomeLocation &peer) const {
        return location > peer.location;
    }

    inline bool operator<=(const GenomeLocation &peer) const {
        return location <= peer.location;
    }

    inline bool operator<(const GenomeLocation &peer) const {
        return location < peer.location;
    }

    inline bool operator!=(const GenomeLocation &peer) const {
        return location != peer.location;
    }

    inline GenomeLocation operator+(const GenomeDistance distance) const {
        GenomeLocation retVal(location + distance);
        return retVal;
    }

    inline GenomeDistance operator-(const GenomeLocation &otherLoc) const {
        return location - otherLoc.location;
    }

    inline GenomeLocation operator-(const GenomeDistance distance) const {
        return location - distance;
    }

    inline GenomeLocation operator+=(const GenomeDistance distance)  {
        location += distance;
        return *this;
    }

    inline GenomeLocation operator-=(const GenomeDistance distance) {
        location -= distance;
        return *this;
    }

    _int64         location;
};

inline _int64 GenomeLocationAsInt64(GenomeLocation genomeLocation) {
    return genomeLocation.location;
}

inline unsigned GenomeLocationAsInt32(GenomeLocation genomeLocation) {
    _ASSERT(genomeLocation.location <= 0xffffffff && genomeLocation.location >= 0);
    return (unsigned)genomeLocation.location;
}
#else   // _DEBUG
typedef _int64 GenomeLocation;

inline _int64 GenomeLocationAsInt64(GenomeLocation genomeLocation)
{
    return genomeLocation;
}

inline unsigned GenomeLocationAsInt32(GenomeLocation genomeLocation) {
    _ASSERT(genomeLocation <= 0xffffffff && genomeLocation>= 0);    // One might wonder about the value of an _ASSERT in code that's only non-_DEBUG.  Think of it as an uppity comment.  :-)
    return (unsigned)genomeLocation;
}

#endif // _DEBUG

typedef _int64 GenomeDistance;

extern GenomeLocation InvalidGenomeLocation;

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
            GenomeDistance          i_maxBases,
            GenomeDistance          nBasesStored,
            unsigned                i_chromosomePadding,
            unsigned                maxContigs = 32);

        void startContig(
            const char          *contigName);

        void addData(
            const char          *data);

        void addData(const char *data, GenomeDistance len);

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
        static const Genome *loadFromFile(const char *fileName, unsigned chromosomePadding, GenomeLocation i_minLocation = 0, GenomeDistance length = 0, bool map = false);
                                                                  // This loads from a genome save
                                                                  // file, not a FASTA file.  Use
                                                                  // FASTA.h for FASTA loads.

        static bool getSizeFromFile(const char *fileName, GenomeDistance *nBases, unsigned *nContigs);

        bool saveToFile(const char *fileName) const;

        //
        // Methods to read the genome.
        //
		inline const char *getSubstring(GenomeLocation location, GenomeDistance lengthNeeded) const {
			if (location > nBases || location + lengthNeeded > nBases + N_PADDING) {
				// The first part of the test is for the unsigned version of a negative offset.
				return NULL;
			}

			// If we're in the padding, then the base will be an n, and we can't short circuit.  Recall that we use lower case n in the reference so it won't match with N in the read.
			if (lengthNeeded <= chromosomePadding && bases[GenomeLocationAsInt64(location)] != 'n') {
				return bases + (location - minLocation);
			}

			_ASSERT(location >= minLocation && location + lengthNeeded <= maxLocation + N_PADDING); // If the caller asks for a genome slice, it's only legal to look within it.

			if (lengthNeeded == 0) {
				return bases + (location - minLocation);
			}

			const Contig *contig = getContigAtLocation(location);
			if (NULL == contig) {
				return NULL;
			}

			_ASSERT(contig->beginningLocation <= location && contig->beginningLocation + contig->length >= location);
			if (contig->beginningLocation + contig->length <= location + lengthNeeded) {
				return NULL;
			}

			return bases + (location - minLocation);
		}

        inline GenomeDistance getCountOfBases() const {return nBases;}

        bool getLocationOfContig(const char *contigName, GenomeLocation *location, int* index = NULL) const;

        inline void prefetchData(GenomeLocation genomeLocation) const {
            _mm_prefetch(bases + GenomeLocationAsInt64(genomeLocation), _MM_HINT_T2);
            _mm_prefetch(bases + GenomeLocationAsInt64(genomeLocation) + 64, _MM_HINT_T2);
        }

        struct Contig {
            Contig() : beginningLocation(InvalidGenomeLocation), length(0), nameLength(0), name(NULL) {}
            GenomeLocation     beginningLocation;
            GenomeDistance     length;
            unsigned           nameLength;
            char              *name;
        };

        inline const Contig *getContigs() const { return contigs; }

        inline int getNumContigs() const { return nContigs; }

        const Contig *getContigAtLocation(GenomeLocation location) const;
        const Contig *getContigForRead(GenomeLocation location, unsigned readLength, GenomeDistance *extraBasesClippedBefore) const;
        const Contig *getNextContigAfterLocation(GenomeLocation location) const;
        int getContigNumAtLocation(GenomeLocation location) const;    // Returns the contig number, which runs from 0 .. getNumContigs() - 1.

// unused        Genome *copy() const {return copy(true,true,true);}
// unused        Genome *copyGenomeOneSex(bool useY, bool useM) const {return copy(!useY,useY,useM);}

        //
        // These are only public so creators of new genomes (i.e., FASTA) can use them.
        //
        void    fillInContigLengths();
        void    sortContigsByName();

private:

        static const int N_PADDING = 100; // Padding to add on either end of the genome to allow substring reads past it

        //
        // The actual genome.
        char                *bases;       // Will point to offset N_PADDING in an array of nBases + 2 * N_PADDING
        GenomeDistance       nBases;
        GenomeLocation       maxBases;

        GenomeLocation       minLocation;
        GenomeLocation       maxLocation;

        //
        // A genome is made up of a bunch of contigs, typically chromosomes.  Contigs have names,
        // which are stored here.
        //
        int          nContigs;
        int          maxContigs;

        Contig      *contigs;    // This is always in order (it's not possible to express it otherwise in FASTA).

        Contig      *contigsByName;
        Genome *copy(bool copyX, bool copyY, bool copyM) const;

        static bool openFileAndGetSizes(const char *filename, GenericFile **file, GenomeDistance *nBases, unsigned *nContigs, bool map);

        const unsigned chromosomePadding;

		GenericFile_map *mappedFile;
};

GenomeDistance DistanceBetweenGenomeLocations(GenomeLocation locationA, GenomeLocation locationB);
inline bool genomeLocationIsWithin(GenomeLocation locationA, GenomeLocation locationB, GenomeDistance distance)
{
    return DistanceBetweenGenomeLocations(locationA, locationB) <= distance;
}
