# datatest.py
#
# Run data i/o tests on SNAP
#
# There are 3 possibilities for input:
# FQ |SAM | BAM
#
# There are 2 references datatest.fa and datatest2.fa (with an extra refseq)
#
# There are 2 possibilities for output:
# SAM | BAM
#
# There are four possible reference files for output
#   FQ | SAM input file (BAM is like SAM)
#   datatest | datatest2 reference file
#
# Input & reference files are stored in datatest
# Temp files are put in datatest/temp
#

import sys
import os
import shutil
import subprocess

if len(sys.argv) != 4:
    print "usage: %s data_dir snap-aligner bin_dir" % sys.argv[0]
    exit(1)

data = sys.argv[1]
snap = sys.argv[2]
bin = sys.argv[3]
quick = False

run = "" # declare global

def _f(name):
    return os.path.normpath(data + "/" + name.replace("#", run))

def _ff(names):
    return [_f(x) for x in names]
    
def runit(args, tag, strict=False, stdout=None, stdin=None):
    print "> %s" % ' '.join(args)
    fout = _f(("temp/stdout-%s" % tag)) if stdout == None else stdout
    ferr = _f("temp/stderr-%s" % tag)
    retcode = subprocess.call(args, stdout=open(fout, "w"), stderr=open(ferr, "w"))
    if retcode != 0:
        print "Run %s exited with %d" % (tag.replace("#", run), retcode)
        print open(ferr, "r").read(),
        if strict:
            exit(1)
        return False
    else:
        return True

# setup data & temp directories with all needed input files

temp = os.path.normpath(data + "/temp")
if not quick:
    if os.path.exists(temp):
        shutil.rmtree(temp)
    os.mkdir(temp)
    # create indexex
    runit([snap, "index", _f("datatest.fa"), _f("temp/datatest.idx"), "-c"], "snap-index", strict=True)
    runit([snap, "index", _f("datatest2.fa"), _f("temp/datatest2.idx"), "-c", "-O500"], "snap-index", strict=True)

runs = 0
succeeded = 0
for input_format in ["fq", "bam", "sam"]:
    for index in ["datatest", "datatest2"]:
        for output_format in ["sam", "bam"]:
            runs += 1
            temps = [] # temporary files, deleted on success
            # build & run snap command line
            outfile = _f("temp/%s-%s.%s" % (input_format, index, output_format))
            args = [snap, "single", _f("temp/%s.idx" % index), _f("datatest." + input_format), "-t", "1", "-rg", "group1", "-o", outfile]
            run = "%s-%s-%s" % (input_format, index, output_format)
            ok = runit(args, "snap-#")
            temps.append(outfile)
            if not ok: continue
            # validate output
            ok = runit(["java", "-jar", _f("ValidateSamFile.jar"), "input=" + outfile, "output=" + _f("temp/validate-#")], "validate-#")
            temps.append(_f("temp/validate-#"))
            if not ok: continue
            # translate output to sam if needed
            check_output = outfile
            if output_format == "bam":
                check_output = outfile + ".sam"
                ok = runit([bin + "samtools", "view", "-h", outfile, "-o", check_output], "bam2sam-#")
                temps.append(check_output)
                if not ok: continue
            # compare with reference file
            # don't yet translate attrs from bam<->sam
            if input_format == "fq" or input_format == output_format:
                # remove @PG line from output file
                checkfile = _f("temp/nopg-#.sam")
                temps.append(checkfile)
                ok = runit([bin + "grep", "-v", "@PG", check_output],"grep-#", stdout=checkfile)
                if not ok: continue
                reffile = _f("correct-%s-%s.sam" % ((input_format if input_format != "bam" else "sam"), index))
                ok = runit([bin + "diff", checkfile, reffile], "diff-#")
                if not ok: continue
            # delete temp files
            for f in temps:
                os.remove(f)
            succeeded += 1

print "completed %d runs, %d failures" % (runs, runs - succeeded)
