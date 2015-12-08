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

    def to_sam_pair(self, other):
        r1 = "{}\t{}\t{}\t{}\t{}\t{}M\t{}\t{}\t{}\t{}\t{}\n".format(
            self.id, 99, self.chr, self.pos, 60, len(self.seq), other.chr,
            other.pos, abs(self.pos - other.pos + len(other.seq)), self.seq, 'A'*len(self.seq))
        return r1 + "{}\t{}\t{}\t{}\t{}\t{}M\t{}\t{}\t{}\t{}\t{}\n".format(
            other.id, 147, other.chr, other.pos, 60, len(other.seq), self.chr,
            self.pos, abs(self.pos - other.pos + len(other.seq)), other.seq, 'A'*len(other.seq))

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

    def add_alt(self, name, accession, parent, start, stop, isRC=False, pmut=0.1):
        pc = self.contigs[parent]
        altseq = random_mutate(pc.seq[start:stop], pmut)
        if (isRC):
            altseq = rc(altseq)
        self.add(Contig(name, accession, altseq, True, parent, start, isRC))

    def get_seq(self, chr, start, end):
        return self.contigs[chr].seq[start:end]

    def make_read(self, chr, pos, isRC=False, len=100, pmut=.02, id=None):
        if id == None:
            id = "r{:05d}_{}_{}_{}".format(random.randint(0,99999), chr, pos, ('r' if isRC else 'f'))
        return Read(id, chr, pos, random_mutate(self.get_seq(chr, pos, pos + len), pmut))

    def make_pair(self, chr1, pos1, chr2, pos2, len=100, pmut=.02):
        id = "r{:05d}_{}_{}_{}_{}".format(random.randint(0,99999), chr1, pos1, chr2, pos2)
        r1 = self.make_read(chr1, pos1, False, len, pmut, id + "/1")
        r2 = self.make_read(chr2, pos2, True, len, pmut, id + "/2")
        return [r1, r2]

    def write_fasta(self, filename):
        with open(filename, 'w') as file:
            for contig in self.contigs.values():
                file.write(">{}|gb|{}\n".format(contig.name, contig.accession))
                for i in range(0, len(contig.seq), 80):
                    file.write("{}\n".format(contig.seq[i:i+80]))

    def write_alts(self, filename):
        with open(filename, 'w') as file:
            file.write("#alt_scaf_acc\tparent_acc\tori\talt_scaf_start\talt_scaf_stop\tparent_start\tparent_stop\talt_start_tail\talt_stop_tail\n")
            for contig in self.contigs.values():
                if contig.isAlt:
                    file.write("{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\n".format(
                        contig.accession, self.contigs[contig.parent].accession, '-' if contig.isAltRC else '+',
                        1, len(contig.seq), 1 + contig.parentLoc, contig.parentLoc + len(contig.seq), 0, 0))

g = Genome()
g.add(Contig("chr1", "C01", random_bases(2000)))
g.add_alt("chr1a", "C01A", "chr1", 500, 1500)
g.write_fasta("test.fa")
g.write_alts("test_alts.txt")

with open("test.sam", "w") as file:
    for i in range(0, 101, 10):
        [r1, r2] = g.make_pair('chr1', 500 + i, 'chr1a', i)
        file.write(r1.to_sam_pair(r2))
