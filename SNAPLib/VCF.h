#pragma once

#include <queue>
#include "Mutation.h"
#include "Compat.h"
#include "Genome.h"

namespace sim {

// 
// Write a VCF file.
// 
// Variant Call Format version 4.0: <http://www.1000genomes.org/node/101>.
//
// Polymorphisms fed to this class must be in order and nonoverlapping, with nonempty alleles.
// 
class VCFWriter {
public:
    explicit VCFWriter(FILE *i_output) : output(i_output) {
        fprintf(output,
                "##fileformat=VCFv4.0\n"
                "##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">\n"
                "##source=JessePrototype\n"
                "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tNA00001\n");
    }
    ~VCFWriter() {
        fclose(output);
    }
    //
    // Write a nontrivial haploid or diploid mutation.
    // Internally, this function converts from 0-based to 1-based chromosome coordinates.
    // 
    template<class AlleleChange>
    void addLine(const Locus &locus, const AlleleChange &alleleChange) {
        addLine(locus, alleleChange.getRef(), alleleChange.getAlt());
    }
    void addLine(const Locus &locus, const Allele &ref, const Allele &alt) {
        if (alt != ref) {
            startLine(locus, ref);
            endLine(alt, "1");
        }
    }
    void addLine(const Locus &locus, const Allele &ref, const Gtype &gtype);

private:
    FILE *output;
    void startLine(const Locus &locus, const Allele &ref) {
        fprintf(output, "%s\t%u\t.\t%s\t", locus.getChrom().c_str(), locus.getPos() + 1,
                ref.c_str());
    }
    void endLine(const std::string &alleles, const std::string &gtype) {
        static const char *boilerplate = "200\tPASS\t.\tGT"; // Huge Phred score.
        fprintf(output, "%s\t%s\t%s\n", alleles.c_str(), boilerplate, gtype.c_str());
        _ASSERT(!ferror(output));
    }
};

//
// Hide the fact that VCF alleles can't be empty.
// 
template<class AlleleChange> class VCFBuffer {
public:
    VCFBuffer(const Genome &i_genome, VCFWriter *i_vcfWriter) :
        genome(i_genome), vcfWriter(i_vcfWriter) {}
    ~VCFBuffer() {
        flush();
    }
    //
    // Add a mutation at 'locus'.
    //
    // 'locus' must be valid (0-based) in the genome used to construct this 'VCFBuffer'.
    // Each allele must be an ACGT string, possibly empty.
    // Mutations must be added in order, without overlapping.
    // 
    void addMutation(const Locus &locus, const AlleleChange &alleleChange) {
        _ASSERT(!alleleChange.isInsertion() || !alleleChange.isDeletion());
        Line line(locus, alleleChange);
        if (!alleleChange.isInsertion() && !alleleChange.isDeletion()) {
            pushLine(line);
        } else {
            if (alleleChange.isInsertion()) {  // In wgsim's stdout, an insertion follows the base.
                pushLine(trivialLine(locus, genome).merge(line));
            } else if (alleleChange.isDeletion()) {
                if (!locus.getPos())
                    fprintf(stderr,
                            "WARNING: Jesse made me ignore deletion at chromosome start.\n");
                else if (!buffer.empty() && buffer.back().justBefore(locus)) {
                    buffer.back().merge(line);
                } else {
                    Locus shiftedLocus(locus.getChrom(), locus.getPos() - 1);
                    pushLine(trivialLine(shiftedLocus, genome).merge(line));
                }
            }
        }
    }
    //
    // Write any unwritten VCF lines.
    // 
    void flush() {
        while (!buffer.empty()) {
            writeLine(buffer.front());
            buffer.pop();
        }
    }
private:
    struct Line {
        Line(const Locus &i_locus, const AlleleChange &i_alleleChange) :
            locus(i_locus), alleleChange(i_alleleChange) {}
        Locus locus;
        AlleleChange alleleChange;
        bool justBefore(const Locus &nextLocus) const {
            return locus.getChrom() == nextLocus.getChrom() &&
                locus.getPos() + alleleChange.numReplaced() == nextLocus.getPos();
        }
        const Line &merge(const Line &that) {
            alleleChange.merge(that.alleleChange);
            return *this;
        }
    };
    static Line trivialLine(const Locus &locus, const Genome &genome) {
        return Line(locus, AlleleChange::trivial(locus, genome));
    }
    void writeLine(const Line &line) {
        vcfWriter->addLine(line.locus, line.alleleChange);
    }
    void pushLine(const Line &line) {
        flush();
        buffer.push(line);
    }
    std::queue<Line> buffer;
    const Genome &genome;
    VCFWriter *vcfWriter;
};

typedef VCFBuffer<AlleleChange> HaploidVCFBuffer;
typedef VCFBuffer<DiploidAlleleChange> DiploidVCFBuffer;

}
