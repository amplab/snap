#include "stdafx.h"
#include "Mutator.h"
#include "Util.h"
#include "GoodRandom.h"
#include "OffsetMap.h"

using std::vector;
using std::string;
using std::map;
using util::getOrElse;
using sim::Allele;
using sim::Base;
using sim::Pos;

namespace {

struct Mutation {
    unsigned numReplaced;   // How many bases this genome replaces.
    Allele replacement;   // May be empty, for a pure deletion.
};

//
// Helper for mutating a single chromosome.
// 
class ChromMutator
{
public:
    ChromMutator() : madeOffsetMap(false) {}
    //
    // Append to `output' a mutant of reference sequence ['refBegin', 'refEnd').
    //
    void mutate(const Base *refBegin, const Base *refEnd, Genome *output) const;
    void addMutation(Pos pos, const Mutation &mutation) {
        Mutations::const_reverse_iterator prev(mutations.upper_bound(pos)); // <http://stackoverflow.com/a/9503696>
        if (prev != mutations.rend() && pos < prev->first + prev->second.numReplaced) {
            printf("Ignoring overlapping mutation at 0-based chromosomal %u\n", pos);
        } else {
            mutations[pos] = mutation;
            madeOffsetMap = false;
        }
    }
    Pos getTranslation(Pos pos) const {
        if (!madeOffsetMap)
            makeOffsetMap();
        return offsetMap.getTranslation(pos);
    }
private:
    typedef map<Pos, Mutation> Mutations;
    Mutations mutations;
    mutable sim::OffsetMap offsetMap;
    mutable bool madeOffsetMap;
    void makeOffsetMap() const;
};

    void
ChromMutator::mutate(const Base *refBegin, const Base *refEnd, Genome *output) const
{
    const Base *b = refBegin;       // Next base to look at.
    for (Mutations::const_iterator m = mutations.begin(); m != mutations.end(); ++m) {
        const Mutation &mutation = m->second;
        const Base *mBegin = refBegin + m->first;
        _ASSERT(mBegin >= b);
        // printf("DEBUG: copying %lu unmutated bases starting with %c\n", mBegin - b, *b);
        output->addData(b, mBegin - b);
        // printf("DEBUG: m at position %u\n", m->first);
        output->addData(mutation.replacement.c_str());
        b = mBegin + mutation.numReplaced;
    }
    _ASSERT(b <= refEnd);       // TODO: Improve error handling.
    output->addData(b, refEnd - b);
}
    
    void
ChromMutator::makeOffsetMap() const
{
    for (Mutations::const_iterator m = mutations.begin(); m != mutations.end(); ++m) {
        const Mutation &mutation = m->second;
        offsetMap.addMutation(m->first, mutation.numReplaced, (unsigned)mutation.replacement.size());
    }
    madeOffsetMap = true;
}

}

namespace sim {

class Mutator
{
public:
    Mutator() : basesGained(0) {}
    void addMutation(const string& chrom, Pos pos, unsigned numReplaced,
                     const Allele &replacement) {
        Mutation mutation = {numReplaced, replacement};
        chromMutators[chrom].addMutation(pos, mutation);
        unsigned numReplacements = (unsigned)replacement.size();
        if (numReplacements > numReplaced) // Overcounting is okay.
            basesGained += numReplacements - numReplaced;
    }
    const Genome *mutate(const Genome &genome) const;
    Pos getTranslation(const string &chrom, Pos pos) const {
        return getOrElse(chromMutators, chrom).getTranslation(pos);
    }
private:
    unsigned basesGained;
    map<string, ChromMutator> chromMutators;
};

    const Genome *
Mutator::mutate(const Genome &genome) const
{
  Genome *newGenome = new Genome(genome.getCountOfBases() + basesGained,genome.getCountOfBases() + basesGained);
    int nPieces = genome.getNumPieces();
    const Genome::Piece *pieces = genome.getPieces();
    for (int i = 0; i < nPieces; ++i) {
        const Genome::Piece &piece = pieces[i];
        unsigned start = piece.beginningOffset;
        unsigned end = i + 1 < nPieces ? pieces[i + 1].beginningOffset : genome.getCountOfBases();
        unsigned size = end - start;
        const Base *bases = genome.getSubstring(start, size);
        const char *name = piece.name;
        // printf("DEBUG: mutating %s of size %u, from %u to %u\n", name, size, start, end);
        newGenome->startPiece(name);
        getOrElse(chromMutators, string(name)).mutate(bases, bases + size, newGenome);
    }
    return newGenome;
}

DiploidMutator::DiploidMutator() {
    mutators = new Mutator[SEXES];
}

DiploidMutator::~DiploidMutator() {
    delete[] mutators;
}

    const DiploidGenome *
DiploidMutator::mutate(const Genome &genome, Sex sex) const
{
    const DiploidGenome *diploidGenome = DiploidGenome::Factory(&genome, sex == MALE);
    const Genome *mutated[SEXES];
    FORSEX
        mutated[sex] = mutators[sex].mutate(*diploidGenome->getGenome(sex == FEMALE));
    return DiploidGenome::Factory(mutated[0], mutated[1]);
}

    void
DiploidMutator::addMutation(const string &chrom, Pos pos, unsigned numReplaced,
                            const string &replacement, Sex sex)
{
    mutators[sex].addMutation(chrom, pos, numReplaced, replacement);
}

    Pos
DiploidMutator::getTranslation(Sex sex, const std::string &chrom, Pos pos) const
{
    return mutators[sex].getTranslation(chrom, pos);
}

}
