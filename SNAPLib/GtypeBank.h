#pragma once

#include <functional>
#include "stdafx.h"
#include "Mutation.h"
#include "Mutator.h"
#include "Chrom.h"

namespace sim {

// 
// Probability distribution over genotypes.
//
// The total probability mass might be
// slightly different from 1.0 due to rounding in the input probabilities, or
// significantly less than 1.0 for a male X chromosome, where we ignore heterozygous genotypes.
//
// For probabilities, we use dbSNP's frequency precision: thousandths.
//
class GtypeDist {
public:
    GtypeDist() : probabilityMass(0.0), homozygosity(true) {}
    const Gtype *sample() const;
    void update(const Gtype &gtype, double probability) {
        dist[gtype] = gtype.homozygous() ? probability : dist[gtype.swapped()] = probability / 2.0;
        probabilityMass += probability;
        if (!gtype.homozygous())
            homozygosity = false;
    }
    bool empty() const {
        return probabilityMass < 0.001;
    }
    bool homozygous() const {
        return homozygosity;
    }
    GtypeDist homozygousVersion() const {
        GtypeDist homo;
        for (Dist::const_iterator p = dist.begin(); p != dist.end(); ++p) {
            const Gtype &gtype = p->first;
            if (gtype.homozygous())
                homo.update(gtype, p->second);
        }
        return homo;
    }
    //
    // How likely is a sampled allele to satisfy a given predicate?
    //
    template <class Pred> double alleleProbability(Pred pred) const {
        _ASSERT(!empty());
        double prob = 0;
        for (Dist::const_iterator p = dist.begin(); p != dist.end(); ++p)
            FORSEX
                if (pred(*p->first.getAllele(Sex(sex))))
                    prob += p->second;
        return prob / (2.0 * probabilityMass);
    }
    //
    // Does a given predicate hold for all alleles in this genotype distribution?
    //
    template <class Pred> bool forAlleles(Pred pred) const {
        for (Dist::const_iterator p = dist.begin(); p != dist.end(); ++p)
            FORSEX
                if (!pred(*p->first.getAllele(Sex(sex))))
                    return false;
        return true;
    }
private:
    typedef std::map<Gtype, double> Dist;
    Dist dist;
    double probabilityMass;
    bool homozygosity;
};

//
// A SNP in the generalized sense of dbSNP; a probability distribution over
// arbitrary replacements of a fixed segment of a reference chromosome.
// 
struct SNP {
    //
    // How likely is a copy of the reference allele to mutate?
    // 
    double variation() const {
        return gtypeDist.alleleProbability(std::bind2nd(std::not_equal_to<Allele>(), ref));
    }
    bool isSnp();
    std::string ref;
    GtypeDist gtypeDist;
};

//
// Genotype statistics.
// 
class GtypeStats {
public:
    GtypeStats() : totalSNPs(0), variation(0.0) {}
    void add(const SNP &snp, int ploidy) {
        _ASSERT(0 <= ploidy && ploidy <= 2);
        ++totalSNPs;
        variation += snp.variation() * ploidy;
    }
    void show(FILE *out=stdout) const {
        fprintf(out,
                "  total SNPs (in the general sense of dbSNP): %u\n"
                "  expected mutations: %g\n", totalSNPs, variation);
    }
private:
    Pos totalSNPs;
    double variation;
};

//
// A collection of probability distributions over genotypes.
// 
// (The abbreviation "gtype" comes from dbSNP's XML Genotype Exchange format.)
// 
class GtypeBank {
public:
    explicit GtypeBank(Sex sex) : chromSexer(sex) {}
    //
    // Add a SNP unless it would overlap with an existing SNP.
    // 
    bool addSNP(const std::string &chrom, Pos pos, const SNP &snp) {
        _ASSERT(!snp.gtypeDist.empty());
        SNPs *onChrom = &snps[chrom];
        int ploidy = chromSexer.sex(chrom).size();
        switch (ploidy) {
        case 2:
            return addSNP(pos, snp, onChrom, ploidy);
        case 1:
            const GtypeDist &gtypeDist = snp.gtypeDist;
            if (gtypeDist.homozygous()) {
                return addSNP(pos, snp, onChrom, ploidy);
            } else {
                SNP homo = {snp.ref, gtypeDist.homozygousVersion()};
                return !homo.gtypeDist.empty() && addSNP(pos, homo, onChrom, ploidy);
            }
        }
    }
    bool addSureSNP(const Chrom &chrom, Pos pos, const Allele &ref, const Gtype &gtype) {
        SNP snp = {ref};
        snp.gtypeDist.update(gtype, 1.0);
        return addSNP(chrom, pos, snp);
    }
    const ChromSexer *getChromSexer() {
        return &chromSexer;
    }

    //
    // By sampling all SNPs, create a deterministic generator of fake genomes of fixed sex.
    // 
    // Dump the diff from 'baseGenome' in 'vcfFileName'.
    // 
    void stageMutations(const Genome &baseGenome,
                        const std::string &vcfFileName, DiploidMutator *) const;
    DiploidMutator makeDiploidMutator(const Genome &baseGenome,
                                      const std::string &vcfFileName) const {
        DiploidMutator diploidMutator;
        stageMutations(baseGenome, vcfFileName, &diploidMutator);
        return diploidMutator;
    }

    void showStats(FILE *out=stdout) const {
        stats.show(out);
    }
private:
    typedef std::map<Pos, SNP> SNPs;
    bool addSNP(Pos pos, const SNP &snp, SNPs *snps, int ploidy);
    std::map<std::string, SNPs> snps;
    ChromSexer chromSexer;
    GtypeStats stats;
};

}
