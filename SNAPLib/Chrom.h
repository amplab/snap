#pragma once

#include "stdafx.h"
#include "Genome.h"
#include "Mutation.h"

//
// Chromosome utilities.
// 

namespace sim {

//
// Ad hoc helper: return whether 'pieceName' denotes 'chrom'.
// 
inline bool matchesChromName(const std::string &pieceName, const std::string &chrom)
{
    std::string prefix = "chr";
    _ASSERT(!chrom.find(prefix));
    bool mitochondrial = chrom == "chrMT";
    if (pieceName == chrom || mitochondrial && pieceName == "chrM")
        return true;            // Matches the format Jesse first worked with.
    std::string longForm =
        mitochondrial ? "mitochondrion" : "chromosome " + chrom.substr(prefix.size());
    return pieceName.find(longForm) != std::string::npos; // Matches GRCh37.p5.
}

//
// Return the first piece corresponding to 'chrom' in 'genome'; NULL for failure.
// 
inline const Genome::Piece *getChrom(const Genome &genome, const std::string &chrom)
{
    for (int i = 0; i < genome.getNumPieces(); ++i) {
        const Genome::Piece *piece = genome.getPieces() + i;
        if (matchesChromName(piece->name, chrom))
            return piece;
    }
    return NULL;
}

class ChromSexer {
public:
    explicit ChromSexer(Sex i_child) : cache(4), child(i_child) {
        FORSEX {
            cache[sex].push_back(Sex(sex));
            cache[2].push_back(Sex(sex));
        }
    }
    //
    // Does a child of given 'sex' inherit 'chrom' from mom, dad, both, or neither?
    //
    const std::vector<Sex> &sex(const std::string &chrom) const {
        return cache[getIndex(chrom)];
    }
private:
    std::vector<std::vector<Sex> > cache;
    Sex child;
    int getIndex(const std::string &chrom) const {
        if (matchesChromName(chrom, "chrX") && child || matchesChromName(chrom, "chrMT"))
            return FEMALE;
        if (matchesChromName(chrom, "chrY"))
            return child ? MALE : 3;
        return 2;
    }
};

} // namespace sim
