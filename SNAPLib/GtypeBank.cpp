#define DEBUG

#include "GtypeBank.h"
#include "Util.h"
#include "GoodRandom.h"
#include "Chrom.h"
#include "VCF.h"

using std::string;
using std::vector;
using std::map;
using util::die;

namespace sim {

namespace {

bool isSingleBase(const Allele &allele) {
    return allele.size() == 1;
}

}

    bool
SNP::isSnp()
{
    return isSingleBase(ref) && gtypeDist.forAlleles(isSingleBase);
}

    const Gtype *
sim::GtypeDist::sample() const
{
    _ASSERT(!empty());
    double draw = GoodFastRandom(probabilityMass * 1000 - 1) / 1000.0;
    double lower = 0.0;
    for (Dist::const_iterator p = dist.begin(); p != dist.end(); ++p) {
        double upper = lower + p->second;
        if (lower <= draw && draw < upper)
            return &p->first;
        lower = upper;
    }
    die("Bad genotype probability distribution.");
}

namespace {

typedef map<Pos, SNP> SNPs;

void addDiploidSNPs(const string &chrom, const SNPs &snps, const Genome &genome,
                    DiploidMutator *diploidMutator, VCFWriter *vcfWriter)
{
    DiploidVCFBuffer vcfBuffer(genome, vcfWriter);
    for (SNPs::const_iterator p = snps.begin(); p != snps.end(); ++p) {
        const SNP &snp = p->second;
        const Gtype &gtype = *snp.gtypeDist.sample();
        FORSEX
            diploidMutator->addMutation(chrom, p->first, snp.ref.size(),
                                        *gtype.getAllele(Sex(sex)), Sex(sex));
        vcfBuffer.addMutation(Locus(chrom, p->first),
                              DiploidAlleleChange(snp.ref, gtype));
    }
}

//
// TODO: Refactor for more encapsulation and less code duplication.
// 
void addHaploidSNPs(const string &chrom, const SNPs &snps, const Genome &genome, Sex parent,
                    DiploidMutator *diploidMutator, VCFWriter *vcfWriter)
{
    HaploidVCFBuffer vcfBuffer(genome, vcfWriter);
    for (SNPs::const_iterator p = snps.begin(); p != snps.end(); ++p) {
        const SNP &snp = p->second;
        const Gtype &gtype = *snp.gtypeDist.sample();
        const Allele &allele = *gtype.getAllele(parent);
        diploidMutator->addMutation(chrom, p->first, snp.ref.size(), allele, parent);
        vcfBuffer.addMutation(Locus(chrom, p->first), AlleleChange(snp.ref, allele));
    }
}

} // namespace

    void
GtypeBank::stageMutations(const Genome &baseGenome, const string &vcfFileName,
                          DiploidMutator *diploidMutator) const
{
    VCFWriter vcfWriter(fopen(vcfFileName.c_str(), "wb"));
    for (map<string, SNPs>::const_iterator p = snps.begin(); p != snps.end(); ++p) {
        const string &chrom = p->first;
        const vector<Sex> &sexes = chromSexer.sex(chrom);
        switch (sexes.size()) {
        case 2:
            addDiploidSNPs(chrom, p->second, baseGenome, diploidMutator, &vcfWriter);
            break;
        case 1:
            addHaploidSNPs(chrom, p->second, baseGenome, sexes[0], diploidMutator, &vcfWriter);
            break;
        default:
            _ASSERT(sexes.empty());
        }
    }
}

namespace {

//
// If 'snp' were mapped to 'prev', would 'next' fall after the reference allele?
//
// Assume the reference allele is nonempty.
// 
bool after(Pos prev, const SNP&snp, Pos next)
{
    return prev + snp.ref.size() <= next;
}

typedef std::map<Pos, SNP> SNPs;

//
// Return whether 'snp' wouldn't overlap with an existing SNP in 'snps' if added at 'pos'.
// 
bool wouldFit(Pos pos, const SNP &snp, const SNPs &snps)
{
   SNPs::const_iterator next = snps.upper_bound(pos);
   SNPs::const_reverse_iterator prev(next); // <http://stackoverflow.com/a/9503696>
   return (next == snps.end() || after(pos, snp, next->first)) &&
       (prev == snps.rend() || after(prev->first, prev->second, pos));
}

} // namespace

    bool
GtypeBank::addSNP(Pos pos, const SNP &snp, SNPs *snps, int ploidy)
{
    bool add = wouldFit(pos, snp, *snps);
    if (add) {
        (*snps)[pos] = snp;
        stats.add(snp, ploidy);
    }
    return add;
}

} // namespace sim
