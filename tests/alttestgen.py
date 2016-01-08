# alttestgen.py
#
# generate data files for alt-contig test
#
# test is alttest.py
#

import sys
import os
import shutil
import subprocess
import random

import pandas as pd

BASES = "ACTG"
RCBASES = {"A":"T", "T":"A", "C":"G", "G":"C"}

def random_bases(n):
    result = ""
    for i in range(n):
        result = result + random.choice(BASES)
    return result

def random_mutate(seq, p = 0.02):
    for i in range(len(seq)):
        if random.random() <= p:
            b = "ACTG".find(seq[i:i+1])
            seq = seq[:i] + random.choice(BASES[:b] + BASES[b+1:]) + seq[i + 1:]
    return seq

def rc(seq):
    result = ""
    for c in seq:
        result = RCBASES[c] + result
    return result

class Read:
    def __init__(self, id, chr, pos, seq, qual=None):
        self.id = id
        self.chr = chr
        self.pos = pos
        self.seq = seq
        self.qual = qual

    def __str__(self):
        return "Read({}, {}, {}, {})".format(self.id, self.chr, self.pos, self.seq)

class Contig:
    def __init__(self, name, accession, seq, isAlt=False, parent=None, parentLoc = 0, isAltRC=False):
        self.name = name
        self.accession = accession
        self.seq = seq
        self.isAlt = isAlt
        self.parent = parent
        self.parentLoc = parentLoc
        self.isAltRC = isAltRC

    def __str__(self):
        return "Contig({}, {}, {}, {}, {}, {}, {})".format(
            self.name, self.accession, self.seq, 'alt' if self.isAlt else 'ref',
            self.parent, self.parentLoc, 'rc' if self.isAltRC else '')

class Genome:
    def __init__(self, contigs={}):
        self.contigs = contigs
    
    def add(self, contig):
        self.contigs[contig.name] = contig

    def add_alt(self, name, accession, parent, start, stop, isRC=False, pmut=0.01):
        pc = self.contigs[parent]
        altseq = random_mutate(pc.seq[start:stop], pmut)
        if (isRC):
            altseq = rc(altseq)
        self.add(Contig(name, accession, altseq, True, parent, start, isRC))

    def make_read(self, chr, pos, isReverse=False, len=100, pmut=.01, id=None):
        if id == None:
            id = "r{:05d}_{}_{}_{}".format(random.randint(0,99999), chr, pos+1, ('r' if isReverse else 'f'))
        seq = random_mutate(self.contigs[chr].seq[pos:pos + len], pmut)
        return Read(id, chr, pos, seq)

    def make_pair(self, chr1, pos1, chr2, pos2, len=100, pmut=.01):
        id = "r{:05d}_{}_{}_{}_{}".format(random.randint(0,99999), chr1, pos1 + 1, chr2, pos2 + 1)
        r1 = self.make_read(chr1, pos1, False, len, pmut, id + "/1")
        r2 = self.make_read(chr2, pos2, True, len, pmut, id + "/2")
        return [r1, r2]

    def to_sam_pair(self, read1, read2):
        rc1 = 1 if self.contigs[read1.chr].isAltRC else 0
        rc2 = 0 if self.contigs[read2.chr].isAltRC else 1
        r1 = "{}\t{}\t{}\t{}\t{}\t{}M\t{}\t{}\t{}\t{}\t{}\n".format(
            read1.id, 67+16*rc1+32*rc2, read1.chr, read1.pos + 1, 60, len(read1.seq), read2.chr,
            read2.pos + 1, abs(read1.pos - read2.pos + len(read2.seq)), read1.seq, (['ABCD','DCBA'][rc1]*int(len(read1.seq)/4+1))[:len(read1.seq)])
        return r1 + "{}\t{}\t{}\t{}\t{}\t{}M\t{}\t{}\t{}\t{}\t{}\n".format(
            read2.id, 131+16*rc2+32*rc1, read2.chr, read2.pos + 1, 60, len(read2.seq), read1.chr,
            read1.pos + 1, abs(read1.pos - read2.pos + len(read2.seq)), read2.seq, (['ABCD','DCBA'][rc2]*int(len(read2.seq)/4+1))[:len(read2.seq)])

    def write_fasta(self, filename):
        with open(filename, 'w') as file:
            for write_alts in [False, True]:
                cnames = self.contigs.keys()
                cnames.sort()
                for cname in cnames:
                    contig = self.contigs[cname]
                    if contig.isAlt == write_alts:
                        file.write(">{}|gb|{}\n".format(contig.name, contig.accession))
                        LINE_LEN=100
                        for i in range(0, len(contig.seq), LINE_LEN):
                            file.write("{}\n".format(contig.seq[i:i+LINE_LEN]))

    def write_alts(self, filename):
        with open(filename, 'w') as file:
            file.write("#alt_scaf_acc\tparent_acc\tori\talt_scaf_start\talt_scaf_stop\tparent_start\tparent_stop\talt_start_tail\talt_stop_tail\n")
            for contig in self.contigs.values():
                if contig.isAlt:
                    file.write("{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\n".format(
                        contig.accession, self.contigs[contig.parent].accession, '-' if contig.isAltRC else '+',
                        1, len(contig.seq), 1 + contig.parentLoc, contig.parentLoc + len(contig.seq), 0, 0))

g = Genome()
seq = random_bases(7000)
seq = seq + random_mutate(seq[5000:6000]) + random_bases(1000) + random_mutate(seq[5000:6000]) + random_bases(1000)
g.add(Contig("chr1", "C01", seq))
g.add_alt("chr1a", "C01A", "chr1", 1000, 2000)
g.add_alt("chr1b", "C01B", "chr1", 3000, 4000, True)
g.add_alt("chr1c", "C01C", "chr1", 5000, 6000)
g.add_alt("chr1d", "C01D", "chr1", 7000, 8000, True)
g.add_alt("chr1e", "C01E", "chr1", 9000, 10000)
g.add_alt("chr1f", "C01F", "chr1", 9000, 10000)
g.write_fasta("test.fa")
g.write_alts("test_alts.txt")

with open("test.sam", "w") as file:
    for i in [0, 100, 600, 800, 900]:
        for j in range(5):
            start = j * 2000
            chralt = ['chr1a', 'chr1b', 'chr1c', 'chr1d', 'chr1e'][j]
            ialt = 900 - i if g.contigs[chralt].isAltRC else i
            [r1, r2] = g.make_pair('chr1', i + start, chralt, ialt)
            file.write(g.to_sam_pair(r1,r2))
            [r1, r2] = g.make_pair('chr1', i + start, 'chr1', i + start + 1000)
            file.write(g.to_sam_pair(r1,r2))
            [r1, r2] = g.make_pair(chralt, ialt, 'chr1', i + start + 2000)
            file.write(g.to_sam_pair(r1,r2))
            [r1, r2] = g.make_pair('chr1', i + start + 1000, 'chr1', i + start + 2000)
            file.write(g.to_sam_pair(r1,r2))
