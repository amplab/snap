/*++

Module Name:

    Geonome.h

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

//
// A wrapper for contig nums in the order of the orginal FASTA.
// It's done this way to make its type different from internal contig nums.
//
class OriginalContigNum {
public:
    OriginalContigNum(int contigNumFromFASTA_) {
        contigNumFromFASTA = contigNumFromFASTA_;
    }

    OriginalContigNum() {
        contigNumFromFASTA = -1;
    }

    //
    // We use unsigned comparisons because we want -1, which is unaligned,
    // to be bigger than any ordinary contig.  We have to use -1 because
    // it gets written into BAM files and the BAM spec requires it.
    //
    bool operator==(const OriginalContigNum& peer) const {
        return (unsigned)contigNumFromFASTA == (unsigned)peer.contigNumFromFASTA;
    }

    bool operator!=(const OriginalContigNum& peer) const {
        return (unsigned)contigNumFromFASTA != (unsigned)peer.contigNumFromFASTA;
    }

    bool operator<(const OriginalContigNum& peer) const {
        return (unsigned)contigNumFromFASTA < (unsigned)peer.contigNumFromFASTA;
    }

    bool operator<=(const OriginalContigNum& peer) const {
        return (unsigned)contigNumFromFASTA <= (unsigned)peer.contigNumFromFASTA;
    }

    bool operator>(const OriginalContigNum& peer) const {
        return (unsigned)contigNumFromFASTA > (unsigned)peer.contigNumFromFASTA;
    }

    bool operator>=(const OriginalContigNum& peer) const {
        return (unsigned)contigNumFromFASTA >= (unsigned)peer.contigNumFromFASTA;
    }

private:

    int     contigNumFromFASTA;

    friend int OriginalContigNumToInt(const OriginalContigNum originalContigNum);

};

inline int OriginalContigNumToInt(const OriginalContigNum originalContigNum) {
    return originalContigNum.contigNumFromFASTA;
}

//
// A wrapper for contig nums in the order of our index (i.e., with ALTs later than others).
// It's done this way to make its type different from original contig nums.
//
class InternalContigNum {
public:
    InternalContigNum(int contigNumFromIndex_)
    {
        contigNumFromIndex = contigNumFromIndex_;
    }

    InternalContigNum() {
        contigNumFromIndex = -1;
    }

    bool operator==(const InternalContigNum& peer) const {
        return contigNumFromIndex == peer.contigNumFromIndex;
    }

    bool operator!=(const InternalContigNum& peer) const {
        return contigNumFromIndex != peer.contigNumFromIndex;
    }

    bool operator<(const InternalContigNum& peer) const {
        return contigNumFromIndex < peer.contigNumFromIndex;
    }

    bool operator<=(const InternalContigNum& peer) const {
        return contigNumFromIndex <= peer.contigNumFromIndex;
    }

    bool operator>(const InternalContigNum& peer) const {
        return contigNumFromIndex > peer.contigNumFromIndex;
    }

    bool operator>=(const InternalContigNum& peer) const {
        return contigNumFromIndex >= peer.contigNumFromIndex;
    }

private:

    int     contigNumFromIndex;
    friend int InternalContigNumToInt(const InternalContigNum internalContigNum);

};

inline int InternalContigNumToInt(const InternalContigNum internalContigNum) {
    return internalContigNum.contigNumFromIndex;
}

#else   // _DEBUG
typedef _int64 GenomeLocation;

inline _int64 GenomeLocationAsInt64(GenomeLocation genomeLocation)
{
    return genomeLocation;
}

inline unsigned GenomeLocationAsInt32(GenomeLocation genomeLocation) {
    _ASSERT(genomeLocation <= 0xffffffff && genomeLocation >= 0);    // One might wonder about the value of an _ASSERT in code that's only non-_DEBUG.  Think of it as an uppity comment.  :-)
    return (unsigned)genomeLocation;
}

typedef int InternalContigNum;
typedef unsigned OriginalContigNum; // Unsigned so that -1 is big rather than small.  See comment in the _DEBUG version for details.
inline int OriginalContigNumToInt(const OriginalContigNum originalContigNum) { return originalContigNum; }
inline int InternalContigNumToInt(const InternalContigNum internalContigNum) { return internalContigNum; }

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
            unsigned                maxContigs);

		void startContig(
			const char          *contigName,
            int                  originalContigNumber);

		void markContigALT(
			const char			*contigName);

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
        static const Genome *loadFromFile(
								const char *fileName, 
								unsigned chromosomePadding, 
								GenomeLocation i_minLocation = 0, 
								GenomeDistance length = 0, 
								bool map = false);
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

        bool getLocationOfContig(const char *contigName, GenomeLocation *location, InternalContigNum* index = NULL) const;

        inline void prefetchData(GenomeLocation genomeLocation) const {
            _mm_prefetch(bases + GenomeLocationAsInt64(genomeLocation), _MM_HINT_T2);
            _mm_prefetch(bases + GenomeLocationAsInt64(genomeLocation) + 64, _MM_HINT_T2);
        }

        struct Contig {
            Contig() : beginningLocation(InvalidGenomeLocation), length(0), nameLength(0), name(NULL), isALT(false), originalContigNumber(OriginalContigNum(-1)), internalContigNumber(InternalContigNum(-1)) {}
            GenomeLocation     beginningLocation;
            GenomeDistance     length;                  // This includes the padding
            unsigned           nameLength;
            InternalContigNum  internalContigNumber;    // The contig number in the order in the index; this corresponds to 
            OriginalContigNum  originalContigNumber;    // Where was this contig in the FASTA file?  We reorder them because of ALTs, but then undo that at sort time.
            char              *name;
			bool			   isALT;
        };

        //inline const Contig *getContigs() const { return contigs; }

        inline const Contig* getContigByInternalNumber(InternalContigNum internalContigNum) const {
            return &contigs[InternalContigNumToInt(internalContigNum)];
        }

        inline const Contig* getContigByOriginalContigNumber(OriginalContigNum originalContigNum) const {
            return getContigByInternalNumber(contigNumberByOriginalOrder[OriginalContigNumToInt(originalContigNum)]);
        }

        // This maps the original contig number to the index of the reordered contig, i.e., the index in contigs[].  It's the inverse mapping of contigs[i].originalContigNumber
        inline const InternalContigNum* getContigNumbersByOriginalOrder() const { return contigNumberByOriginalOrder; }    

        inline int getNumContigs() const { return nContigs; }

		int getNumALTContigs() const;

        const Contig *getContigAtLocation(GenomeLocation location) const;
        const Contig *getContigForRead(GenomeLocation location, unsigned readLength, GenomeDistance *extraBasesClippedBefore) const;
        const Contig *getNextContigAfterLocation(GenomeLocation location) const;
        int getContigNumAtLocation(GenomeLocation location) const;    // Returns the contig number, which runs from 0 .. getNumContigs() - 1.

        //
        // These are only public so creators of new genomes (i.e., FASTA) can use them.
        //
        void    fillInContigLengths();
        void    sortContigsByName();

        inline bool isGenomeLocationALT(GenomeLocation location) const {
            return location >= genomeLocationOfFirstALTContig;
        }

        char* genomeLocationInStringForm(GenomeLocation location, char *buffer, size_t bufferSize) const;

        void setUpContigNumbersByOriginalOrder();

private:

        static const int N_PADDING = 1000; // Padding to add on either end of the genome to allow substring reads past it

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


        Contig              *contigsByName;
        InternalContigNum   *contigNumberByOriginalOrder;
 
        static bool openFileAndGetSizes(const char *filename, GenericFile **file, GenomeDistance *nBases, unsigned *nContigs, bool map);

        const unsigned chromosomePadding;

		GenericFile_map *mappedFile;

        GenomeLocation  genomeLocationOfFirstALTContig;
};

//
// A wrapper for GenomeLocation that does comparison by the original order of contigs from the FASTA that generated the index, not
// by GenomeLocation space.  The trick here is that the index *almost* preserves the original order.  The only thing that it does
// is to move all of the ALT contigs to the end, so that the aligner can test whether a GenomeLocation is ALT by just comparing it
// with the cutoff.  So, to compare two genome locations in the original order if they're both non-ALT or both ALT then you can
// still just compare the raw location.  If they're one of each then you need to go to the genome and look up the contig structure,
// and from there you can just look at the original contig number to make the comparison.
//
class GenomeLocationOrderedByOriginalContigs
{
    GenomeLocation  location;
    const Genome    *genome;

public:
    GenomeLocationOrderedByOriginalContigs(GenomeLocation i_location, const Genome* i_genome) : location(i_location), genome(i_genome) {}

    inline bool operator>=(GenomeLocationOrderedByOriginalContigs& peer)
    {
        if (location == InvalidGenomeLocation) {
            return true;
        }

        if (peer.location == InvalidGenomeLocation) {
            return false;
        }

        if (genome->isGenomeLocationALT(location) == genome->isGenomeLocationALT(peer.location))
        {
            return location >= peer.location;
        }

        return genome->getContigAtLocation(location)->originalContigNumber >= genome->getContigAtLocation(peer.location)->originalContigNumber;
    }

    inline bool operator>(GenomeLocationOrderedByOriginalContigs& peer)
    {
        if (location == peer.location) {
            return false;
        }

        if (location == InvalidGenomeLocation) {
            return true;
        }

        if (peer.location == InvalidGenomeLocation) {
            return false;
        }

        if (genome->isGenomeLocationALT(location) == genome->isGenomeLocationALT(peer.location))
        {
            return location > peer.location;
        }

        return genome->getContigAtLocation(location)->originalContigNumber > genome->getContigAtLocation(peer.location)->originalContigNumber;
    }

    inline bool operator<(GenomeLocationOrderedByOriginalContigs& peer)
    {
        if (location == peer.location) {
            return false;
        }

        if (location == InvalidGenomeLocation) {
            return false;
        }

        if (peer.location == InvalidGenomeLocation) {
            return true;
        }

        if (genome->isGenomeLocationALT(location) == genome->isGenomeLocationALT(peer.location))
        {
            return location < peer.location;
        }

        return genome->getContigAtLocation(location)->originalContigNumber < genome->getContigAtLocation(peer.location)->originalContigNumber;
    }

    inline bool operator<=(GenomeLocationOrderedByOriginalContigs& peer)
    {
        if (location == peer.location) {
            return true;
        }

        if (location == InvalidGenomeLocation) {
            return false;
        }

        if (peer.location == InvalidGenomeLocation) {
            return true;
        }

        if (genome->isGenomeLocationALT(location) == genome->isGenomeLocationALT(peer.location))
        {
            return location <= peer.location;
        }

        return genome->getContigAtLocation(location)->originalContigNumber <= genome->getContigAtLocation(peer.location)->originalContigNumber;
    }

    inline bool operator==(GenomeLocationOrderedByOriginalContigs& peer)
    {
        return location == peer.location;
    }

    inline bool operator!=(GenomeLocationOrderedByOriginalContigs& peer)
    {
        return location != peer.location;
    }
};

class ContigAndPos {
public:
    ContigAndPos(OriginalContigNum originalContigNum_, int pos_) {
        originalContigNum = originalContigNum_;
        pos = pos_;
    }

    ContigAndPos() {
        originalContigNum = -1;
        pos = 0;
    }

    inline bool operator==(const ContigAndPos& peer) const {
        return originalContigNum == peer.originalContigNum && pos == peer.pos;
    }

    inline bool operator!=(const ContigAndPos& peer) const {
        return originalContigNum != peer.originalContigNum || pos != peer.pos;
    }

    inline bool operator>=(const ContigAndPos& peer) const {
        if (originalContigNum == peer.originalContigNum) {
            return pos >= peer.pos;
        }

        return originalContigNum > peer.originalContigNum;
    }

    inline bool operator>(const ContigAndPos& peer) const {
        if (originalContigNum == peer.originalContigNum) {
            return pos > peer.pos;
        }

        return originalContigNum > peer.originalContigNum;
    }

    inline bool operator<=(const ContigAndPos& peer) const {
        if (originalContigNum == peer.originalContigNum) {
            return pos <= peer.pos;
        }

        return originalContigNum < peer.originalContigNum;
    }

    inline bool operator<(const ContigAndPos& peer) const {
        if (originalContigNum == peer.originalContigNum) {
            return pos < peer.pos;
        }

        return originalContigNum < peer.originalContigNum;
    }

    inline OriginalContigNum getContig() { return originalContigNum; }
    inline int getPos() { return pos; }

private:
    OriginalContigNum originalContigNum;
    int pos;
};

GenomeDistance DistanceBetweenGenomeLocations(GenomeLocation locationA, GenomeLocation locationB);
inline bool genomeLocationIsWithin(GenomeLocation locationA, GenomeLocation locationB, GenomeDistance distance)
{
    return DistanceBetweenGenomeLocations(locationA, locationB) <= distance;
}
