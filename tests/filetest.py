# filetest.py
#
# Run file i/o tests on SNAP
#
# There are 4x2 possibilities for input:
# FQ | FQZ | SAM | BAM
# Single | paired
#
# There is a compressed reference in data/xx.fa.gz
#
# There are 3 input files stored:
#   xx_in1.fq.gz xx_in2.fq.gz xx_in12.bam
# The others are created from these using gunzip and samtools
# FQ/FQZ single : in1
# FQ/FQZ paired : in1 & in2
# BAM/SAM single/paired : in12
#
# There are 2x2 possibilities for output:
# SAM | BAM
# Unsorted | sorted
#
# These are validated against 4 output files using samtools & diff
#   xx_ref1.bam xx_sorted_ref1.bam xx_ref12.bam xx_sorted_ref12.bam
#
# Input & output files are stored in data
# Temp files are put in data/temp
#

import sys
import os
import shutil
import subprocess

if (len(sys.argv) < 5 or len(sys.argv) > 6):
    print "usage: %s data_dir file_base snap-aligner bin_dir [-quick]" % sys.argv[0]
    exit(1)

data = sys.argv[1]
template = sys.argv[2]
snap = sys.argv[3]
bin = sys.argv[4]
quick = len(sys.argv) == 6

run = "" # declare global

def _f(name):
    return os.path.normpath(data + "/" + name.replace("^", template).replace("#", run))

def _ff(names):
    return [_f(x) for x in names]
    
def runit(args, tag, strict=False, stdout=None, stdin=None):
    print "$ %s%s%s" % (' '.join(args), "" if stdin == None else " < "+stdin, "" if stdout == None else " > "+stdout)
    fout = _f(("temp/stdout-%s" % tag)) if stdout == None else stdout
    ferr = _f("temp/stderr-%s" % tag)
    retcode = subprocess.call(args, stdout=open(fout, "w"), stderr=open(ferr, "w"), stdin=(None if stdin==None else open(stdin, "r")))
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

	runit([bin + "gzip", "-d", "-c", _f("^.fa.gz")], "gunzip-ref", strict=True, stdout=_f("temp/^.fa"))
	runit([bin + "gzip", "-d", "-c", _f("^_in1.fq.gz")], "gunzip1", strict=True, stdout=_f("temp/^_in1.fq"))
	runit([bin + "gzip", "-d", "-c", _f("^_in2.fq.gz")], "gunzip2", strict=True, stdout=_f("temp/^_in2.fq"))
	runit([bin + "samtools", "view", "-h", _f("^_in12.bam"), "-o", _f("temp/^_in12.sam")], "samtools1", strict=True)
	runit([bin + "samtools", "view", _f("^_ref1.bam"), "-o", _f("temp/^_ref1u.sam")], "samtools2", strict=True)
	runit([bin + "sort", _f("temp/^_ref1u.sam")], "namesort-ref1", stdout=_f("temp/^_ref1.sam"), strict=True)
	os.remove(_f("temp/^_ref1u.sam"))
	runit([bin + "samtools", "view", _f("^_sorted_ref1.bam"), "-o", _f("temp/^_sorted_ref1.sam")], "samtools3", strict=True)
	runit([bin + "samtools", "view", _f("^_ref12.bam"), "-o", _f("temp/^_ref12u.sam")], "samtools4", strict=True)
	runit([bin + "sort", _f("temp/^_ref12u.sam")], "namesort-ref12", stdout=_f("temp/^_ref12.sam"), strict=True)
	os.remove(_f("temp/^_ref12u.sam"))
	runit([bin + "samtools", "view", _f("^_sorted_ref12.bam"), "-o", _f("temp/^_sorted_ref12.sam")], "samtools5", strict=True)
	# create index
	runit([snap, "index", _f("temp/^.fa"), _f("temp/^.idx")], "snap-index", strict=True)

inputs = {
    "fq": [["temp/^_in1.fq"], ["temp/^_in1.fq", "temp/^_in2.fq"]],
    "fq.gz": [["^_in1.fq.gz"], ["^_in1.fq.gz", "^_in2.fq.gz"]],
    "sam": [["temp/^_in12.sam"], ["temp/^_in12.sam"]],
    "bam": [["^_in12.bam"], ["^_in12.bam"]]
}

formats = {
    "fq": "-fastq",
    "fq.gz": "-compressedFastq",
    "sam": "-sam",
    "bam": "-bam"
}

runs = 0
failed = 0
for input_format in ["fq", "fq.gz", "bam", "sam"]:
    for paired in [0, 1]:
        for input_format_2 in ["", "fq", "fq.gz", "bam", "sam"]:
            test_output_format = ["sam", "bam"] if input_format == "fq" and input_format_2 == "" else ["bam"]
            test_sorted = [0, 1] if (input_format == "fq" or input_format == "bam") and input_format_2 == "" else [0]
            for output_format in test_output_format:
                for sorted in test_sorted:
                    test_pipe = [0, 1] if sorted == 0 and (input_format == "bam" or input_format == "sam" or paired == 0) and input_format_2 == "" else [0]
                    for pipe in test_pipe:
                        test_threads = ["1", "4"] if input_format == "fq" and output_format == "bam" and sorted == 1 and pipe == 0 else ["4"]
                        for threads in test_threads:
                            runs += 1
                            temps = [] # temporary files, deleted on success
                            # build & run snap command line
                            args = [snap, ["single", "paired"][paired], _f("temp/^.idx")]
                            args = args + ["-t", threads, "-b", "-="]
                            if pipe == 0:
                                args = args + _ff(inputs[input_format][paired])
                                if input_format_2 != "":
                                    args = args + _ff(inputs[input_format_2][paired])
                            else:
                                args = args + [formats[input_format], "-"]
                            if sorted:
                                args.append("-so")
                            run = "%s-%s-%s-%s-%s-%s-%s" % (input_format, input_format_2, ["single", "paired"][paired], output_format, ["unsorted", "sorted"][sorted], threads, ["file", "pipe"][pipe])
                            outfile = _f("temp/output-#." + output_format)
                            if pipe == 0:
                                args = args + ["-o", outfile]
                            else:
                                args = args + ["-o", formats[output_format], "-"]
                            ok = runit(args, "snap-#") if pipe ==0 else runit(args, "snap-#", stdin=_f(inputs[input_format][paired][0]), stdout=outfile)
                            temps.append(outfile)
                            if not ok:
                                failed += 1
                                continue
                            # validate output
                            ok = runit(["java", "-jar", _f("ValidateSamFile.jar"), "input=" + outfile, "output=" + _f("temp/validate-#")], "validate-#")
                            temps.append(_f("temp/validate-#"))
                            if not ok:
                                failed += 1
                                continue
                            if False:
                                # todo: need more sophisticated diff
                                # remove header, convert to SAM if needed
                                samfile = _f("temp/output-nh-#.sam")
                                ok = runit([bin + "samtools", "view"] + {"bam":[], "sam":["-S"]}[output_format] + [outfile, "-o", samfile], "samtools-view-#")
                                temps.append(samfile)
                                if not ok:
                                    failed += 1
                                    continue
                                outfile = samfile
                                # sort by name if not sorted by coordinate
                                if sorted == 0:
                                    sortfile = _f("temp/output-namesort-#.sam")
                                    ok = runit([bin + "sort", outfile], "namesort-#", stdout=sortfile)
                                    temps.append(sortfile)
                                    if not ok:
                                        failed += 1
                                        continue
                                    outfile = sortfile
                                # compare with reference file
                                use12 = "12" if paired or input_format == "sam" or input_format == "bam" else "1"
                                reffile = _f("temp/^%s_ref%s.sam" % (["", "_sorted"][sorted], use12))
                                ok = runit([bin + "diff", outfile, reffile], "diff-#")
                                if not ok:
                                    failed += 1
                                    continue
                            # delete temp files
                            for f in temps:
                                os.remove(f)
print "completed %d runs, %d failures" % (runs, failed)
