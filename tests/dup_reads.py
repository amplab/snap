# dup_reads.py
#
# create duplicate reads
#

import sys
import random

def readread(f):
    result = [f.readline(),f.readline(),f.readline(),f.readline()]
    if (result[0] and (result[0][0] != "@" or result[2][0] != "+" or len(result[1]) != len(result[3]))):
        sys.stderr.write("invalid fasta file near %s" % (result[0]))
        exit(1)
    return result

def writeread(f, r):
    for i in range(4):
        f.write(r[i])

if (len(sys.argv) < 4 or len(sys.argv) > 5):
    print "usage: %s <# of duplicate reads> <max duplication> read1.fq [read2.fq]" % sys.argv[p]
    exit(1)

dupcount = int(sys.argv[1])
maxdup = int(sys.argv[2])

in1 = open(sys.argv[3], "r")
out1 = open("dup_" + sys.argv[3], "w")
paired = len(sys.argv) >= 5
if paired:
    in2 = open(sys.argv[4], "r")
    out2 = open("dup_" + sys.argv[4], "w")

for i in range(0, dupcount):
    r1 = readread(in1)
    if paired:
        r2 = readread(in2)
    ndup = random.randint(2,maxdup)
    for j in range(0, ndup):
        writeread(out1, ["@dup%d_%s" % (j, r1[0][1:]), r1[1], r1[2], r1[3]])
        if paired:
            writeread(out2, ["@dup%d_%s" % (j, r2[0][1:]), r2[1], r2[2], r2[3]])
