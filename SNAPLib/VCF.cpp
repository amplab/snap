#include "stdafx.h"
#include "Util.h"
#include "VCF.h"

using std::string;
using util::joinWithSep;

    void
sim::VCFWriter::addLine(const Locus &locus, const Allele &ref, const Gtype &gtype)
{
    std::vector<Allele> distinctAlleles(1, ref); // Inefficient.
    std::vector<char> alleleIndices(SEXES, '0');
    FORSEX {
        const Allele &allele = *gtype.getAllele(Sex(sex));
        util::addIfAbsent(&distinctAlleles, allele);
        alleleIndices[sex] += (char)util::findIndex(BEGEND(distinctAlleles), allele);
    }
    if (distinctAlleles.size() > 1) {
        startLine(locus, ref);
        string alleleString = joinWithSep(++distinctAlleles.begin(), distinctAlleles.end(), ',');
        endLine(alleleString, joinWithSep(BEGEND(alleleIndices), '|'));
    }
}
