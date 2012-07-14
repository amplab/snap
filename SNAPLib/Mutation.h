#pragma once

#include "stdafx.h"
#include "Genome.h"

// 
// Header for code involving mutation.
// Based on VCF terminology: see <http://www.1000genomes.org/node/101>.
// 

namespace sim {

typedef std::string Chrom;      // A piece of the genome, e.g., the mitochondrion.

typedef unsigned Pos;           // Position in a chrom.  WARNING: may be 0-based or 1-based.  See:
                                // <http://www.biostars.org/post/show/6373/what-are-the-advantagesdisadvantages-of-one-based-vs-zero-based-genome-coordinate-systems/>

typedef std::string Allele;

class Locus {
public:
    Locus(const Chrom &i_chrom, Pos i_pos) : chrom(i_chrom), pos(i_pos) {}
    const Chrom &getChrom() const {
        return chrom;
    }
    const Pos getPos() const {
        return pos;
    }
private:
    Chrom chrom;
    Pos pos;
};

const int SEXES = 2;            // Convenient approximation.
enum Sex {FEMALE, MALE};        // Sex is parity of number of X or Y chroms.

#define FORSEX for (int sex = 0; sex < SEXES; ++sex)

typedef char Base;

//
// A genotype at a particular locus.
//
class Gtype {
public:
    Gtype(const Allele &allele0, const Allele &allele1) {
        addAllele(allele0);
        addAllele(allele1);
    }
    explicit Gtype(const std::vector<const char *> alleles) {
        FORSEX
            addAllele(alleles[sex]);
    }
    explicit Gtype(const std::string &gtype) {
        sz_t slashPos = gtype.find('/');
        _ASSERT(slashPos != std::string::npos);
        addAllele(gtype, 0, slashPos);
        addAllele(gtype, slashPos + 1, gtype.size());
    }
    const Allele *getAllele(Sex sex) const {
        return &alleles[sex];
    }
    bool homozygous() const {
        return alleles[0] == alleles[1];
    }
    const Gtype swapped() const {
        Gtype ret = *this;
        std::swap(ret.alleles[0], ret.alleles[1]);
        return ret;
    }
    const Gtype &merge(const Gtype &that) {
        FORSEX
            alleles[sex] += that.alleles[sex];
        return *this;
    }
    static Gtype homozygousFactory(const Allele &allele) {
        return Gtype(allele, allele);
    }
private:
    typedef std::string::size_type sz_t;
    void addAllele(const Allele &allele) {
        alleles.push_back(allele == "-" ? "" : allele);
    }
    void addAllele(const std::string &gtype, sz_t first, sz_t last) {
        addAllele(gtype.substr(first, last));
    }
    friend bool operator<(const Gtype &, const Gtype &);
    std::vector<Allele> alleles;
};

inline bool operator<(const Gtype &gtype1, const Gtype &gtype2)
{
    return gtype1.alleles < gtype2.alleles;
}

inline unsigned getOffset(const Locus &locus, const Genome &genome)
{
    unsigned offset;
    if (!genome.getOffsetOfPiece((locus.getChrom()).c_str(), &offset)) {
        fprintf(stderr, "Error: no piece '%s' in genome.\n", locus.getChrom().c_str());
        exit(1);
    }
    return offset + locus.getPos();
}

//
// Replacement of a reference allele with an alternative something.
// 
class AbstractAlleleChange {
public:
    AbstractAlleleChange(const Allele &i_ref) : ref(i_ref) {}
    unsigned numReplaced() const {
        return (unsigned)ref.size();
    }
    const Allele &getRef() const {
        return ref;
    }
    bool isInsertion() const {
        return ref.empty();
    }
protected:
    Allele ref;
};

//
// Replacement of a reference allele with an alternative allele.
// 
class AlleleChange : public AbstractAlleleChange {
public:
    AlleleChange(const Allele &i_ref, const Allele &i_alt) :
        AbstractAlleleChange(i_ref), alt(i_alt) {}
    const Allele &getAlt() const {
        return alt;
    }
    bool isDeletion() const {
        return alt.empty();
    }
    const AlleleChange &merge(const AlleleChange &that) {
        ref += that.ref;
        alt += that.alt;
        return *this;
    }
    static AlleleChange trivial(const Locus &locus, const Genome &genome) {
        unsigned offset = getOffset(locus, genome);
        Allele refBase(genome.getSubstring(offset, 1), 1);
        return AlleleChange(refBase, refBase);
    }
private:
    Allele alt;
};

//
// Replacement of a reference allele with an alternative genotype (i.e., allele pair).
//
// TODO: reduce code duplication with 'AlleleChange', by templatizing 'AbstractAlleleChange'.
// 
class DiploidAlleleChange : public AbstractAlleleChange {
public:
    DiploidAlleleChange(const Allele &i_ref, const Gtype &i_alt) :
        AbstractAlleleChange(i_ref), alt(i_alt) {}
    const Gtype &getAlt() const {
        return alt;
    }
    bool isDeletion() const {
        FORSEX
            if (alt.getAllele(Sex(sex))->empty())
                return true;
        return false;
    }
    const DiploidAlleleChange &merge(const DiploidAlleleChange &that) {
        ref += that.ref;
        alt.merge(that.alt);
        return *this;
    }
    static DiploidAlleleChange trivial(const Locus &locus, const Genome &genome) {
        unsigned offset = getOffset(locus, genome);
        Allele refBase(genome.getSubstring(offset, 1), 1);
        return DiploidAlleleChange(refBase, Gtype::homozygousFactory(refBase));
    }
private:
    Gtype alt;
};

}
