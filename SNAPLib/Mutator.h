#pragma once

#include <map>
#include "stdafx.h"
#include "Genome.h"
#include "Mutation.h"

namespace sim {

class Mutator;           // TODO: Hide.

class DiploidMutator
{
public:
    DiploidMutator();
    ~DiploidMutator();

    //
    // Create a synthetic genome from a reference genome.
    //
    const DiploidGenome *mutate(const Genome &genome, Sex sex) const;

    //
    // Register a mutation in 0-based chromosome coordinates.
    // WARNING: 'pos' must be an offset from the start of the chromosome, not the genome.
    //
    void addMutation(const std::string &chrom, Pos pos, unsigned numReplaced,
                     const std::string &replacement, Sex sex);

    //
    // Chromosome coordinate conversion: mutated to unmutated.
    // 
    Pos getTranslation(Sex sex, const std::string &chrom, Pos pos) const;

private:
    Mutator *mutators; // Mom and dad, in that order.
};

}
