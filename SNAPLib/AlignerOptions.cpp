/*++

Module Name:

    AlignerOptions.cpp

Abstract:

    Common parameters for running single & paired alignment.
    Common parameters for running single & paired alignment.

Authors:

    Ravi Pandya, May, 2012

Environment:

    User mode service.

Revision History:

    Integrated from SingleAligner.cpp & PairedAligner.cpp

--*/

#include "stdafx.h"
#include "options.h"
#include "AlignerOptions.h"
#include "FASTQ.h"
#include "SAM.h"
#include "Bam.h"
#include "exit.h"
#include "Error.h"
#include "BaseAligner.h"
#include "CommandProcessor.h"

AlignerOptions::AlignerOptions(
    const char* i_commandLine,
    bool forPairedEnd)
    :
    commandLine(i_commandLine),
    indexDir(NULL),
    similarityMapFile(NULL),
    numThreads(GetNumberOfProcessors()),
    bindToProcessors(true),
    ignoreMismatchedIDs(false),
    clipping(ClipBack),
    sortOutput(false),
    noIndex(false),
    noDuplicateMarking(false),
    noQualityCalibration(false),
    sortMemory(0),
    filterFlags(0),
    explorePopularSeeds(false),
    stopOnFirstHit(false),
	useM(true),
    gapPenalty(0),
	extra(NULL),
    rgLineContents("@RG\tID:FASTQ\tPL:Illumina\tPU:pu\tLB:lb\tSM:sm"),
    perfFileName(NULL),
    useTimingBarrier(false),
    extraSearchDepth(1),
    defaultReadGroup("FASTQ"),
    seedCountSpecified(false),
    minWeightToCheck(1),
    numSeedsFromCommandLine(0),
    ignoreSecondaryAlignments(true),
    maxSecondaryAlignmentAdditionalEditDistance(-1),
	maxSecondaryAlignments(0x7fffffff),
    maxSecondaryAlignmentsPerContig(-1),    // -1 means don't limit
    preserveClipping(false),
    expansionFactor(1.0),
    noUkkonen(false),
    noOrderedEvaluation(false),
	noTruncation(false),
    useAffineGap(true),
    matchReward(1),
    subPenalty(4),
    gapOpenPenalty(6),
    gapExtendPenalty(1),
    minAGScore(30),
	minReadLength(DEFAULT_MIN_READ_LENGTH),
    maxDistFraction(0.0),
	mapIndex(false),
	prefetchIndex(false),
    writeBufferSize(16 * 1024 * 1024),
    dropIndexBeforeSort(false),
    killIfTooSlow(false),
    sortIntermediateDirectory(NULL),
    profile(false),
    ignoreAlignmentAdjustmentsForOm(true),
    emitInternalScore(false),
	altAwareness(true),
    emitALTAlignments(false),
    maxScoreGapToPreferNonALTAlignment(2)
{
    if (forPairedEnd) {
        maxDist                 = 32;
        seedCoverage            = 0;
        numSeedsFromCommandLine = 8;
        maxHits                 = 300;
     } else {
        maxDist                 = 14;
        numSeedsFromCommandLine = 25;
        maxHits                 = 300;
		seedCoverage			= 0;
    }

    initializeLVProbabilitiesToPhredPlus33();
}

    void
AlignerOptions::usage()
{
    usageMessage();
}

    void
        AlignerOptions::usageMessage()
    {
        WriteErrorMessage(
            "Usage: \n%s\n"
            "Options:\n"
            "  -o   filename  output alignments to filename in SAM or BAM format, depending on the file extension or\n"
            "       explicit type specifier (see below).  Use a dash with an explicit type specifier to write to\n"
            "       stdout, so for example -o -sam - would write SAM output to stdout\n"
            "  -d   maximum edit distance allowed per read or pair (default: %d)\n"
            "  -n   number of seeds to use per read\n"
            "  -sc  Seed coverage (i.e., readSize/seedSize).  Floating point.  Exclusive with -n.  (default uses -n)\n"
            "  -h   maximum hits to consider per seed (default: %d)\n"
            "  -ms  minimum seed matches per location (default: %d)\n"
            "  -t   number of threads (default is one per core)\n"
            "  -b-  Don't bind each thread to its processor (note the double dash)\n"
            "  -P   disables cache prefetching in the genome; may be helpful for machines\n"
            "       with small caches or lots of cores/cache\n"
            "  -so  sort output file by alignment location\n"
            "  -sm  memory to use for sorting in Gb\n"
            " -sid  Specifies the sort intermediate directory.  When SNAP is sorting, it aligns the reads in the order in which they come in, and writes\n"
            "       the aligned reads in batches to a temporary file.  When the aligning is done, it does a merge sort from the temporary file into the\n"
            "       final output file.  By default, the intermediate file is in the same directory as the output file, but for performance or space\n"
            "       reasons, you might want to put it elsewhere.  If so, use this option.\n"            
            "  -x   explore some hits of overly popular seeds (useful for filtering)\n"
            "  -S   suppress additional processing (sorted BAM output only)\n"
            "       i=index, d=duplicate marking\n"
            "  -f   stop on first match within edit distance limit (filtering mode)\n"
            "  -F   filter output (a=aligned only, s=single hit only (MAPQ >= %d), u=unaligned only, l=long enough to align (see -mrl))\n"
            "  -E   an alternate (and fully general) way to specify filter options.  Emit only these types s = single hit (MAPQ >= %d), m = multiple hit (MAPQ < %d),\n"
            "       x = not long enough to align, u = unaligned, b = filter must apply to both ends of a paired-end read.  Combine the letters after\n"
            "       -E, so for example -E smu will emit all reads that aren't too short/have too many Ns (because it leaves off l).  -E smx is the same\n"
            "       as -F a, -E ux is the same as -F u, and so forth.\n"
            "       When filtering in paired-end mode (either with -F or -E) unless you specify the b flag a read will be emitted if it's mate pair passes the filter\n"
            "       Even if the read itself does not.  If you specify b mode, then a read will be emitted only if it and its partner both pass the filter.\n"
#if     USE_DEVTEAM_OPTIONS
            "  -I   ignore IDs that don't match in the paired-end aligner\n"
#ifdef  _MSC_VER    // Only need this on Windows, since memory allocation is fast on Linux
            "  -B   Insert barrier after per-thread memory allocation to improve timing accuracy\n"
#endif  // _MSC_VER
#endif  // USE_DEVTEAM_OPTIONS
            "  -Cxx must be followed by two + or - symbols saying whether to clip low-quality\n"
            "       bases from front and back of read respectively; default: back only (-C-+)\n"
            "  -=   use the new style CIGAR strings with = and X rather than M.  The opposite of -M\n"
            "  -pf  specify the name of a file to contain the run speed\n"
            "  -hp  Indicates to use huge pages (this may speed up alignment and slow down index load).\n"
            "  -D   Specifies the extra search depth (the edit distance beyond the best hit that SNAP uses to compute MAPQ).  Default %d\n"
            "  -rg  Specify the default read group if it is not specified in the input file\n"
            "  -R   Specify the entire read group line for the SAM/BAM output.  This must include an ID tag.  If it doesn't start with\n"
            "       '@RG' SNAP will add that.  Specify tabs by \\t.  Two backslashes will generate a single backslash.\n"
            "       backslash followed by anything else is illegal.  So, '-R @RG\\tID:foo\\tDS:my data' would generate reads\n"
            "       with defualt tag foo, and an @RG line that also included the DS:my data field.\n"
            "  -sa  Include reads from SAM or BAM files with the secondary (0x100) or supplementary (0x800) flag set; default is to drop them.\n"
            "  -om  Output multiple alignments.  Takes as a parameter the maximum extra edit distance relative to the best alignment\n"
            "       to allow for secondary alignments\n"
            " -omax Limit the number of alignments per read generated by -om.  This means that if -om would generate more\n"
            "       than -omax secondary alignments, SNAP will write out only the best -omax of them, where 'best' means\n"
            "       'with the lowest edit distance'.  Ties are broken arbitrarily.\n"
            "  -mpc Limit the number of alignments generated by -om to this many per contig (chromosome/FASTA entry);\n"
            "       'mpc' means 'max per contig; default unlimited.  This filter is applied prior to -omax.  The primary alignment\n"
            "       is counted.\n"
            "  -pc  Preserve the soft clipping for reads coming from SAM or BAM files\n"
            "  -xf  Increase expansion factor for BAM and GZ files (default %.1f)\n"
            "  -hdp Use Hadoop-style prefixes (reporter:status:...) on error messages, and emit hadoop-style progress messages\n"
            "  -mrl Specify the minimum read length to align, reads shorter than this (after clipping) stay unaligned.  This should be\n"
            "       a good bit bigger than the seed length or you might get some questionable alignments.  Default %d\n"
            "  -map Use file mapping to load the index rather than reading it.  This might speed up index loading in cases\n"
            "       where SNAP is run repatedly on the same index, and the index is larger than half of the memory size\n"
            "       of the machine.  On some operating systems, loading an index with -map is much slower than without if the\n"
            "       index is not in memory.  You might consider adding -pre to prefetch the index into system cache when loading\n"
            "       with -map when you don't expect the index to be in cache.\n"
            "  -pre Prefetch the index into system cache.  This is only meaningful with -map, and only helps if the index is not\n"
            "       already in memory and your operating system is slow at reading mapped files (i.e., some versions of Linux,\n"
            "       but not Windows).\n"
            "  -lp  Run SNAP at low scheduling priority (Only implemented on Windows)\n"
#ifdef LONG_READS
            "  -dp  Edit distance as a percentage of read length (single only, overrides -d)\n"
#endif
            "  -nu  No Ukkonen: don't reduce edit distance search based on prior candidates. This option is purely for\n"
            "       evaluating the performance effect of using Ukkonen's algorithm rather than Smith-Waterman, and specifying\n"
            "       it will slow down execution without improving the alignments.\n"
            "  -no  No Ordering: don't order the evalutation of reads so as to select more likely candidates first.  This option\n"
            "       is purely for evaluating the performance effect of the read evaluation order, and specifying it will slow\n"
            "       down execution without improving alignments.\n"
            "  -nt  Don't truncate searches based on missed seed hits.  This option is purely for evaluating the performance effect\n"
            "       of candidate truncation, and specifying it will slow down execution without improving alignments.\n"
            " -wbs  Write buffer size in megabytes.  Don't specify this unless you've gotten an error message saying to make it bigger.  Default 16.\n"
            "  -di  Drop the index after aligning and before sorting.  This frees up memory for the sort at the expense of not having the index loaded for your next run.\n"
            " -kts  Kill if too slow.  Monitor our progress and kill ourself if we're not moving fast enough.  This is intended for use on machines\n"
            "       with limited memory, where some alignment tasks will push SNAP into paging, and take disproportinaltely long.  This allows the script\n"
            "       to move on to the next alignment.  Only works when generating output, and not during the sort phase.  If you're running out of memory\n"
            "       sorting, try using -di.\n"

            " -pro  Profile alignment to give you an idea of how much time is spent aligning and how much waiting for IO\n"
            "  -ae  Apply the end-of-contig soft clipping before the -om processing rather than after it.  A read that's soft clipped because of hanging off one end or the other\n"
            "       of a contig does not have a penalty in its NM tag, but it does in SNAP's internal scoring.  This flag says to use the NM value for -om processing\n"
            "       rather than SNAP's internal score.\n"
            "  -is  Write SNAP's internal score for an alignment into the output.  The value following -is specifies the tag to use, and must be a two letter\n"
            "       value starting with X, Y or Z.  So, -is ZQ will cause SNAP to write ZQ:i:3 on a read with internal score 3.  Generally, the internal scores\n"
            "       are the same as the NM values, except that they contain penalties for soft clipping reads that hang over the end of contigs (but not for\n"
            "       soft clipping that's due to # quality scores or that was present in the input SAM/BAM file and retained due to -pc)\n"
            "  -G-  disable affine gap scoring (default: true)\n"
            "       Scoring parameters (works only when -G- is not used)\n"
            "           cost for match -gm (default: %u)\n"
            "           cost for substitution -gs (default: %u)\n"
            "           cost for opening a gap -go (default: %u)\n"
            "           cost for extending a gap -ge (default: %u)\n"
            "           minimum score threshold for read -gt (default: %u)\n"
            "  -A-  Disable ALT awareness.  The default is to try to map reads to the primary assembly and only to choose ALT alignments when they're much better,\n"
            "       and to compute MAPQ for non-ALT alignments using only non-ALT hits. This flag disables that behavior and results in ALT-oblivious behavior.\n"
            "  -ea  Emit ALT alignments.  When the aligner is ALT aware (i.e., -A- isn't specified) if it finds an ALT alignment that would have been\n"
            "       the primary alignment if -A- had been specified but isn't without -A-, SNAP will emit the read with the supplementary alignment\n"
            "       flag set and MAPQ computed across all potential mappings, both primary and ALT\n"
            "  -asg Maximum score gap to prefer a non-ALT alignment.  If the best non-ALT alignment is more than this much worse than the best ALT alignment\n"
            "       emit the ALT alignment as the primary result rather than as a supplementary result. (default: %u)\n"
		,
            commandLine,
            maxDist,
            maxHits,
			minWeightToCheck,
            MAPQ_LIMIT_FOR_SINGLE_HIT, MAPQ_LIMIT_FOR_SINGLE_HIT, MAPQ_LIMIT_FOR_SINGLE_HIT,
			extraSearchDepth,
			expansionFactor,
			DEFAULT_MIN_READ_LENGTH,
            matchReward,
            subPenalty,
            gapOpenPenalty,
            gapExtendPenalty,
            minAGScore,
            maxScoreGapToPreferNonALTAlignment
        );

    if (extra != NULL) {
        extra->usageMessage();
    }

    WriteErrorMessage("\n\n"
                "You may process more than one alignment without restarting SNAP, and if possible without reloading\n"
                "the index.  In order to do this, list on the command line all of the parameters for the first\n"
                "alignment, followed by a comma (separated by a space from the other parameters) followed by the\n"
                "parameters for the next alignment (including single or paired).  You may have as many of these\n"
                "as you please.  If two consecutive alignments use the same index, it will not be reloaded.\n"
                "So, for example, you could do 'snap-aligner single hg19-20 foo.fq -o foo.sam , paired hg19-20 end1.fq end2.fq -o paired.sam'\n"
                "and it would not reload the index between the single and paired alignments.\n",
                "SNAP doesn't parse the options for later runs until the earlier ones have completed, so if you make\n"
                "an error in one, it may take a while for you to notice.  So, be careful (or check back shortly after\n"
                "you think each run will have completed).\n\n");

    WriteErrorMessage("When specifying an input or output file, you can simply list the filename, in which case\n"
                      "SNAP will infer the type of the file from the file extension (.sam or .bam for example),\n"
                      "or you can explicitly specify the file type by preceding the filename with one of the\n"
                      " following type specifiers (which are case sensitive):\n"
                      "    -fastq\n"
                      "    -compressedFastq\n"
                      "    -sam\n"
                      "    -bam\n"
                      "    -pairedFastq\n"
                      "    -pairedInterleavedFastq\n"
                      "    -pairedCompressedInterleavedFastq\n"
                      "\n"
                      "So, for example, you could specify -bam input.file to make SNAP treat input.file as a BAM file,\n"
                      "even though it would ordinarily assume a FASTQ file for input or a SAM file for output when it\n"
                      "doesn't recoginize the file extension.\n"
                      "In order to use a file name that begins with a '-' and not have SNAP treat it as a switch, you must\n"
                      "explicitly specify the type.  But really, that's just confusing and you shouldn't do it.\n"
                      "Input and output may also be from/to stdin/stdout. To do that, use a - for the input or output file\n"
                      "name and give an explicit type specifier.  So, for example, \n"
                      "snap-aligner single myIndex -fastq - -o -sam -\n"
                      "would read FASTQ from stdin and write SAM to stdout.\n"
    );
}

    bool
        AlignerOptions::parse(
            const char** argv,
            int argc,
            int& n,
            bool* done)
    {
        *done = false;

        if (strcmp(argv[n], "-d") == 0) {
            if (n + 1 < argc) {
                maxDist = atoi(argv[n + 1]);
                n++;
                return true;
            }
        }
        else if (strcmp(argv[n], "-n") == 0) {
            if (n + 1 < argc) {
                if (seedCountSpecified) {
                    WriteErrorMessage("-sc and -n are mutually exclusive.  Please use only one.\n");
                    return false;
                }
                seedCountSpecified = true;
                numSeedsFromCommandLine = atoi(argv[n + 1]);
                n++;
                return true;
            }
        }
        else if (strcmp(argv[n], "-sc") == 0) {
            if (n + 1 < argc) {
                if (seedCountSpecified) {
                    WriteErrorMessage("-sc and -n are mutually exclusive.  Please use only one.\n");
                    return false;
                }
                seedCountSpecified = true;
                seedCoverage = atof(argv[n + 1]);
                numSeedsFromCommandLine = 0;
                n++;
                return true;
            }
        }
        else if (strcmp(argv[n], "-ms") == 0) {
            if (n + 1 < argc) {
                minWeightToCheck = (unsigned)atoi(argv[n + 1]);
                if (minWeightToCheck > 1000) {
                    fprintf(stderr, "-ms must be between 1 and 1000\n");
                    return false;
                }
                n++;
                return true;
            }
        }
        else if (strcmp(argv[n], "-h") == 0) {
            if (n + 1 < argc) {
                maxHits = atoi(argv[n + 1]);
                n++;
                return true;
            }
        }
        else if (strcmp(argv[n], "-c") == 0) { // conf diff is deprecated, but we just ignore it rather than throwing an error.
            if (n + 1 < argc) {
                n++;
                return true;
            }
        }
        else if (strcmp(argv[n], "-a") == 0) { // adaptive conf diff is deprecated, but we just ignore it rather than throwing an error.
            if (n + 1 < argc) {
                n++;
                return true;
            }
        }
        else if (strcmp(argv[n], "-t") == 0) {
            if (n + 1 < argc) {
                numThreads = atoi(argv[n + 1]);
                n++;
                return true;
            }
        }
        else if (strcmp(argv[n], "-o") == 0) {
            int argsConsumed;
            if (!SNAPFile::generateFromCommandLine(argv + n + 1, argc - n - 1, &argsConsumed, &outputFile, false, false)) {
                WriteErrorMessage("Must have a file specifier after -o\n");
                return false;
            }
            if (outputFile.isStdio) {
                AlignerOptions::outputToStdout = true;
            }
            n += argsConsumed;
            return true;
        } else if (strcmp(argv[n], "-P") == 0) {
            doAlignerPrefetch = false;
            return true;
        } else if (strcmp(argv[n], "-kts") == 0) {
            killIfTooSlow = true;
            return true;
        } else if (strcmp(argv[n], "-b") == 0) {
            bindToProcessors = true;
            return true;
        } else if (strcmp(argv[n], "--b") == 0 || strcmp(argv[n], "-b-") == 0) {
            bindToProcessors = false;
            return true;
        } else if (strcmp(argv[n], "-so") == 0) {
            sortOutput = true;
            return true;
        } else if (strcmp(argv[n], "-map") == 0) {
            mapIndex = true;
            return true;
        } else if (strcmp(argv[n], "-pre") == 0) {
            prefetchIndex = true;
            return true;
        } else if (strcmp(argv[n], "-ae") == 0) {
            ignoreAlignmentAdjustmentsForOm = false;
            return true;
        } else if (strcmp(argv[n], "-S") == 0) {
            if (n + 1 < argc) {
                n++;
                for (const char* p = argv[n]; *p; p++) {
                    switch (*p) {
                    case 'i':
                        noIndex = true;
                        break;
                    case 'd':
                        noDuplicateMarking = true;
                        break;
#if 0   // This is unused, so I took it out of the arg parsing to avoid unusual behavior if you were to do something like -S i d q, expecting q to be an input file.
                    case 'q':
                        noQualityCalibration = true;
                        break;

#endif // 0
                    default:
                        return false;
                    }
                }
                return true;
            }
        } else if (strcmp(argv[n], "-sm") == 0) {
            if (n + 1 < argc && argv[n + 1][0] >= '0' && argv[n + 1][0] <= '9') {
                sortMemory = atoi(argv[n + 1]);
                n++;
                return true;
            }
        } else if (strcmp(argv[n], "-is") == 0) {
            if (n + 1 >= argc || strlen(argv[n + 1]) != 2 || argv[n + 1][0] < 'X' || argv[n + 1][0] > 'Z' || argv[n + 1][1] < 'A' || argv[n + 1][1] > 'Z') {
                WriteErrorMessage("-is switch must be followed by two letter tag that consists of X, Y, or Z and a capital letter.\n");
                return false;
            }
            emitInternalScore = true;
            strcpy(internalScoreTag, argv[n + 1]);
            n++;
            return true;
        }
        else if (strcmp(argv[n], "-F") == 0) {
            if (n + 1 < argc) {
                n++;
                if (strcmp(argv[n], "a") == 0) {
                    if (0 != filterFlags) {
                        WriteErrorMessage("Specified -F %s after a previous -F or -E option.  Choose one (or put -F b after -F %s)\n", argv[n], argv[n]);
                        return false;
                    }
                    filterFlags = FilterSingleHit | FilterMultipleHits | FilterTooShort;
                }
                else if (strcmp(argv[n], "s") == 0) {
                    if (0 != filterFlags) {
                        WriteErrorMessage("Specified -F %s after a previous -F or -E option.  Choose one (or put -F b after -F %s)\n", argv[n], argv[n]);
                        return false;
                    }
                    filterFlags = FilterSingleHit | FilterTooShort;
                }
                else if (strcmp(argv[n], "u") == 0) {
                    if (0 != filterFlags) {
                        WriteErrorMessage("Specified -F %s after a previous -F or -E  option.  Choose one (or put -F b after -F %s)\n", argv[n], argv[n]);
                        return false;
                    }
                    filterFlags = FilterUnaligned | FilterTooShort;
                }
                else if (strcmp(argv[n], "l") == 0) {
                    if (0 != filterFlags) {
                        WriteErrorMessage("Specified -F %s after a previous -F or -E  option.  Choose one (or put -F b after -F %s)\n", argv[n], argv[n]);
                        return false;
                    }
                    filterFlags = FilterSingleHit | FilterMultipleHits | FilterUnaligned;
                }
                else if (strcmp(argv[n], "b") == 0) {
                    // ignore paired-end option(s)
                }
                else {
                    WriteErrorMessage("Unknown option type after -F: %s\n", argv[n]);
                    return false;
                }
                return true;
            }
        }
        else if (strcmp(argv[n], "-E") == 0) {
            if (n + 1 < argc) {
                if (0 != filterFlags) {
                    WriteErrorMessage("You can have only one -F and/or -E switch (excepting -F b)\n");
                    return false;
                }
                n++;
                for (int whichChar = 0; whichChar < strlen(argv[n]); whichChar++) {
                    switch (argv[n][whichChar]) {
                    case 's':   filterFlags |= FilterSingleHit; break;
                    case 'm':   filterFlags |= FilterMultipleHits; break;
                    case 'x':   filterFlags |= FilterTooShort; break;
                    case 'u':   filterFlags |= FilterUnaligned; break;
                    case 'b':   filterFlags |= FilterBothMatesMatch; break;
                    default:    WriteErrorMessage("Unrecognized filter type after -E '%c'; must be one of smxub\n", argv[n][whichChar]); return false;
                    }
                }
                return true;
            }
        }
        else if (strcmp(argv[n], "-x") == 0) {
            explorePopularSeeds = true;
            return true;
        }
        else if (strcmp(argv[n], "-f") == 0) {
            stopOnFirstHit = true;
            return true;
#if     USE_DEVTEAM_OPTIONS
        }
        else if (strcmp(argv[n], "-I") == 0) {
            ignoreMismatchedIDs = true;
            return true;
#ifdef  _MSC_VER
        }
        else if (strcmp(argv[n], "-B") == 0) {
            useTimingBarrier = true;
            return true;
#endif  // _MSC_VER
#endif  // USE_DEVTEAM_OPTIONS
        }
        else if (strcmp(argv[n], "-M") == 0) {
            useM = true;
            return true;
        }
        else if (strcmp(argv[n], "-=") == 0) {
            useM = false;
            return true;
        }
        else if (strcmp(argv[n], "-sa") == 0) {
            ignoreSecondaryAlignments = false;
            return true;
        }
        else if (strcmp(argv[n], "-om") == 0) {
            if (n + 1 >= argc) {
                WriteErrorMessage("-om requires an additional value\n");
                return false;
            }
            //
            // Check that the parameter is actually numeric.  This is to avoid having someone do "-om -anotherSwitch" and
            // having the additional switch silently consumed here.
            //
            if (argv[n + 1][0] < '0' || argv[n + 1][0] > '9') {
                WriteErrorMessage("-om requires a numerical parameter.\n");
                return false;
            }
            maxSecondaryAlignmentAdditionalEditDistance = atoi(argv[n + 1]);

            n++;

            return true;
        }
        else if (strcmp(argv[n], "-omax") == 0) {
            if (n + 1 >= argc) {
                WriteErrorMessage("-omax requires an additional value\n");
                return false;
            }

            maxSecondaryAlignments = atoi(argv[n + 1]);

            if (maxSecondaryAlignments <= 0) {
                WriteErrorMessage("-omax must be strictly positive\n");
            }

            n++;

            return true;
        }
        else if (strcmp(argv[n], "-sid") == 0) {
            if (n + 1 >= argc) {
                WriteErrorMessage("-sid requires an additional value\n");
                return false;
            }

            if (argv[n + 1][0] == '-') {
                WriteErrorMessage("The directory name after -sid must not start with a dash (it's just too confusing when compared with a command line switch)\n");
                return false;
            }

            sortIntermediateDirectory = argv[n + 1];
            n++;
            return true;
        }
        else if (strcmp(argv[n], "-mpc") == 0) {
            if (n + 1 >= argc) {
                WriteErrorMessage("-mpc requires an additional value\n");
                return false;
            }

            maxSecondaryAlignmentsPerContig = atoi(argv[n + 1]);

            if (maxSecondaryAlignmentsPerContig <= 0) {
                WriteErrorMessage("-mpc must be strictly positive\n");
                return false;
            }

            n++;

            return true;
        }
        else if (strcmp(argv[n], "-wbs") == 0) {
            if (n + 1 >= argc) {
                WriteErrorMessage("-wbs requires an additional value\n");
                return false;
            }
            //
            // Check that the parameter is actually numeric.  This is to avoid having someone do "-wbs -anotherSwitch" and
            // having the additional switch silently consumed here.
            //
            if (argv[n + 1][0] < '0' || argv[n + 1][0] > '9') {
                WriteErrorMessage("-wbs requires a numerical parameter.\n");
                return false;
            }
            writeBufferSize = atoi(argv[n + 1]) * 1024 * 1024;

            if (writeBufferSize <= 0) {
                WriteErrorMessage("-wbs must be bigger than zero");
                return false;
            }

            n++;

            return true;
        }
        else if (strcmp(argv[n], "-xf") == 0) {
            if (n + 1 < argc) {
                n++;
                expansionFactor = (float)atof(argv[n]);
                return expansionFactor > 0;
            }
        }
        else if (strcmp(argv[n], "-pc") == 0) {
            preserveClipping = true;
            return true;
        }
        else if (strcmp(argv[n], "-pro") == 0) {
            profile = true;
            return true;
        } else if (strcmp(argv[n], "-G") == 0) {
            useAffineGap = true;
            return true;
        } else if (strcmp(argv[n], "-G-") == 0) {
            useAffineGap = false;
            return true;
        } else if (strcmp(argv[n], "-gm") == 0) {
            if (n + 1 >= argc) {
                WriteErrorMessage("-gm requires an additional value\n");
                return false;
            }
            n++;
            matchReward = atoi(argv[n]);
            if (matchReward <= 0) {
                WriteErrorMessage("-gm must be greater than zero");
                return false;
            }
            return true;
        } else if (strcmp(argv[n], "-gs") == 0) {
            if (n + 1 >= argc) {
                WriteErrorMessage("-gs requires an additional value\n");
                return false;
            }
            n++;
            subPenalty = atoi(argv[n]);
            if (subPenalty <= 0) {
                WriteErrorMessage("-gs must be greater than zero");
                return false;
            }
            return true;
        } else if (strcmp(argv[n], "-go") == 0) {
            if (n + 1 >= argc) {
                WriteErrorMessage("-go requires an additional value\n");
                return false;
            }
            n++;
            gapOpenPenalty = atoi(argv[n]);
            if (gapOpenPenalty <= 0) {
                WriteErrorMessage("-go must be greater than zero");
                return false;
            }
            return true;
        } else if (strcmp(argv[n], "-ge") == 0) {
            if (n + 1 >= argc) {
                WriteErrorMessage("-ge requires an additional value\n");
                return false;
            }
            n++;
            gapExtendPenalty = atoi(argv[n]);
            if (gapExtendPenalty <= 0) {
                WriteErrorMessage("-ge must be greater than zero");
                return false;
            }
            return true;
        } else if (strcmp(argv[n], "-gt") == 0) {
            if (n + 1 >= argc) {
                WriteErrorMessage("-gt requires an additional value\n");
                return false;
            }
            n++;
            minAGScore = atoi(argv[n]);
            if (minAGScore <= 0) {
                WriteErrorMessage("-gt must be greater than zero");
                return false;
            }
            return true;
        } else if (strcmp(argv[n], "-mrl") == 0) {
            if (n + 1 < argc) {
                n++;
                minReadLength = atoi(argv[n]);
                return minReadLength > 0;
            }
        } else if (strcmp(argv[n], "-dp") == 0) {
            if (n + 1 < argc) {
                n++;
                maxDistFraction = (float)(0.01 * atof(argv[n]));
                return (!isPaired()) && maxDistFraction > 0.0 && maxDistFraction < 1.0;
            }
        } else if (strcmp(argv[n], "-R") == 0) {
            if (n + 1 < argc) {
                //
                // Check the line for sanity.  It must consist either of @RG\t<fields> or just <fields> (in which
                // case we add the @RG part).  It must contain a field called ID. Fields are separated by tabs.
                // We don't require that the fields be things that are listed in the SAM spec, however, because
                // new ones might be added.
                //
                bool needsRG = !(argv[n + 1][0] == '@' && argv[n + 1][1] == 'R' && argv[n + 1][2] == 'G' && argv[n + 1][3] == '\\' && argv[n + 1][4] == 't');
                const unsigned buflen = (unsigned)(strlen(argv[n + 1]) + 1 + (needsRG ? 4 : 0));
                char* buffer = new char[buflen];
                char* copyToPtr = buffer;
                const char* copyFromPtr = argv[n + 1];
                if (needsRG) {
                    memcpy(copyToPtr, "@RG\t", 4);
                    copyToPtr += 4;
                }

                //
                // First copy the line, converting \t into tabs.
                //

                bool pendingBackslash = false;

                while (*copyFromPtr != '\0') {
                    if (pendingBackslash) {
                        if (*copyFromPtr == 't' || *copyFromPtr == '\\') {
                            _ASSERT((unsigned)(copyToPtr - buffer) < buflen);
                            *copyToPtr = (*copyFromPtr == 't') ? '\t' : '\\';
                            copyToPtr++;
                            copyFromPtr++;
                            pendingBackslash = false;
                        }
                        else {
                            WriteErrorMessage("Unrecognized escape character in -R parameter.  A backslash must be followed by a t or another backslash.\n");
                            return false;
                        }
                    }
                    else {
                        //
                        // Emit the character literally unless it's a backslash.
                        //
                        pendingBackslash = *copyFromPtr == '\\';

                        if (!pendingBackslash) {
                            _ASSERT((unsigned)(copyToPtr - buffer) < buflen);
                            *copyToPtr = *copyFromPtr;
                            copyToPtr++;
                        }

                        copyFromPtr++;
                    }
                } // while
                _ASSERT((unsigned)(copyToPtr - buffer) < buflen);
                *copyToPtr = '\0';  // Null terminate the string

                //
                // Now run through the line looking for <tab>ID:..., and use that to set the default read group.
                //
                int bytesAlong = 0;
                defaultReadGroup = NULL;
                for (int i = 0; NULL == defaultReadGroup && i < strlen(buffer); i++) {
                    switch (bytesAlong) {
                    case 0:
                        if (buffer[i] == '\t') {
                            bytesAlong = 1;
                        }
                        break;

                    case 1:
                        if (buffer[i] == 'I') {
                            bytesAlong = 2;
                        }
                        else {
                            bytesAlong = 0;
                        }
                        break;

                    case 2:
                        if (buffer[i] == 'D') {
                            bytesAlong = 3;
                        }
                        else {
                            bytesAlong = 0;
                        }
                        break;

                    case 3:
                        if (buffer[i] == ':') {
                            if (NULL != defaultReadGroup) {
                                WriteErrorMessage("read group string specified with -R contained more than one ID field.\n");
                                return false;
                            }
                            //
                            // The ID tag starts at i+1.
                            //
                            int idTagSize = 0;
                            for (idTagSize = 0; buffer[i + 1 + idTagSize] != '\t' && buffer[i + 1 + idTagSize] != '\0'; idTagSize++) {
                                // This loop body intentionally left blank.
                            }

                            if (0 == idTagSize) {
                                WriteErrorMessage("The ID tag on the read group line specified by -R must not be empty\n");
                                return false;
                            }
                            char* newReadGroup = new char[(size_t)idTagSize + 1]; // +1 for null.
                            memcpy(newReadGroup, buffer + i + 1, idTagSize);
                            newReadGroup[idTagSize] = '\0';
                            defaultReadGroup = newReadGroup; // +1 for null.

                        }
                        else {
                            bytesAlong = 0;
                        }
                        break;

                    default:
                        WriteErrorMessage("Invalid bytesAlong = %d", bytesAlong);
                        soft_exit(1);
                    } // switch
                } // for

                if (NULL == defaultReadGroup) {
                    WriteErrorMessage("The string specified after -R must include an ID field.\n");
                    return false;
                }

                rgLineContents = buffer;    // This leaks, but so what?
                n++;
                return true;

            }
            else {
                WriteErrorMessage("-R requires a value");
                return false;
            }
        } else if (strcmp(argv[n], "-pf") == 0) {
            if (n + 1 < argc) {
                perfFileName = argv[n + 1];
                n++;
                return true;
            }
            else {
                WriteErrorMessage("Must specify the name of the perf file after -pf\n");
            }
        } else if (strcmp(argv[n], "-rg") == 0) {
            if (n + 1 < argc) {
                char* newReadGroup = new char[strlen(argv[n + 1]) + 1];
                strcpy(newReadGroup, argv[n + 1]);
                defaultReadGroup = newReadGroup;
                n++;
                static const char* format = "@RG\tID:%s\tPL:Illumina\tPU:pu\tLB:lb\tSM:sm";
                char* s = new char[1 + strlen(defaultReadGroup) + strlen(format)];
                sprintf(s, format, defaultReadGroup);
                rgLineContents = s;
                return true;
            }
            else {
                WriteErrorMessage("Must specify the default read group after -rg\n");
            }
        } else if (strcmp(argv[n], "--hp") == 0) {
            BigAllocUseHugePages = false;
            return true;
        } else if (strcmp(argv[n], "-hp") == 0) {
            BigAllocUseHugePages = true;
            return true;
        } else if (strcmp(argv[n], "-hdp") == 0) {
            AlignerOptions::useHadoopErrorMessages = true;
            return true;
        } else if (strcmp(argv[n], "-lp") == 0) {
            SetToLowSchedulingPriority();
            return true;
        } else if (strcmp(argv[n], "-nu") == 0) {
            noUkkonen = true;
            return true;
        } else if (strcmp(argv[n], "-no") == 0) {
            noOrderedEvaluation = true;
            return true;
        } else if (strcmp(argv[n], "-nt") == 0) {
            noTruncation = true;
            return true;
        } else if (strcmp(argv[n], "-di") == 0) {
            dropIndexBeforeSort = true;
            return true;
        } else if (strcmp(argv[n], "-D") == 0) {
            if (n + 1 < argc) {
                extraSearchDepth = atoi(argv[n + 1]);
                n++;
                return true;
            } else {
                WriteErrorMessage("Must specify the desired extra search depth after -D\n");
            }
        } else if (strlen(argv[n]) >= 2 && '-' == argv[n][0] && 'C' == argv[n][1]) {
            if (strlen(argv[n]) != 4 || '-' != argv[n][2] && '+' != argv[n][2] ||
                '-' != argv[n][3] && '+' != argv[n][3]) {

                WriteErrorMessage("Invalid -C argument.\n\n");
                return false;
            }

            if ('-' == argv[n][2]) {
                if ('-' == argv[n][3]) {
                    clipping = NoClipping;
                } else {
                    clipping = ClipBack;
                }
            } else {
                if ('-' == argv[n][3]) {
                    clipping = ClipFront;
                } else {
                    clipping = ClipFrontAndBack;
                }
            }
            return true;
        }
        else if (strcmp(argv[n], "-A-") == 0) {
            altAwareness = false;
            return true;
        } else if (strcmp(argv[n], "-ea") == 0) {
            emitALTAlignments = true;
            return true;
        } else if (strcmp(argv[n], "-asg") == 0) {
            if (n + 1 < argc) {
                maxScoreGapToPreferNonALTAlignment = atoi(argv[n + 1]);
                n++;
                return true;
            } else {
                WriteErrorMessage("Must specify the ALT score gap after -ag\n");
                return false;
            }
        } else if (strcmp(argv[n], ",") == 0) {
            //
            // End of args for this run.
            //
            *done = true;
            return true;
        } else if (extra != NULL) {
            return extra->parse(argv, argc, n, done);
        }
    return false;
}

    bool
AlignerOptions::passFilter(
    Read* read,
    AlignmentResult result,
	bool tooShort,
    bool secondaryAlignment)
{
    if (filterFlags == 0) {
        return true;
    }

    if (tooShort) {
        return (filterFlags & FilterTooShort) != 0;
    }

    if (result == MultipleHits && secondaryAlignment && (filterFlags & FilterSingleHit)) { // Don't filter out secondary alignments for low MAPQ
        return true;
    }

    switch (result) {
    case NotFound:
    case UnknownAlignment:
        return (filterFlags & FilterUnaligned) != 0;
    case SingleHit:
        return (filterFlags & FilterSingleHit) != 0;
    case MultipleHits:
        return (filterFlags & FilterMultipleHits) != 0;
    default:
        return false; // shouldn't happen!
    }
}

    PairedReadSupplierGenerator *
SNAPFile::createPairedReadSupplierGenerator(int numThreads, bool quicklyDropUnpairedReads, const ReaderContext& context)
{
    _ASSERT(fileType == SAMFile || fileType == BAMFile || fileType == InterleavedFASTQFile || secondFileName != NULL); // Caller's responsibility to check this

    switch (fileType) {
    case SAMFile:
        return SAMReader::createPairedReadSupplierGenerator(fileName, numThreads, quicklyDropUnpairedReads, context);
        
    case BAMFile:
        return BAMReader::createPairedReadSupplierGenerator(fileName,numThreads, quicklyDropUnpairedReads, context);

    case FASTQFile:
        return PairedFASTQReader::createPairedReadSupplierGenerator(fileName, secondFileName, numThreads, context, isCompressed);

    case InterleavedFASTQFile:
        return PairedInterleavedFASTQReader::createPairedReadSupplierGenerator(fileName, numThreads, context, isCompressed);
        
    default:
        _ASSERT(false);
        WriteErrorMessage("SNAPFile::createPairedReadSupplierGenerator: invalid file type (%d)\n", fileType);
        soft_exit(1);
        return NULL;
    }
}

    ReadSupplierGenerator *
SNAPFile::createReadSupplierGenerator(int numThreads, const ReaderContext& context)
{
    _ASSERT(secondFileName == NULL);
    switch (fileType) {
    case SAMFile:
        return SAMReader::createReadSupplierGenerator(fileName, numThreads, context);
        
    case BAMFile:
        return BAMReader::createReadSupplierGenerator(fileName,numThreads, context);

    case FASTQFile:
        return FASTQReader::createReadSupplierGenerator(fileName, numThreads, context, isCompressed);

    default:
        _ASSERT(false);
        WriteErrorMessage("SNAPFile::createReadSupplierGenerator: invalid file type (%d)\n", fileType);
        soft_exit(1);
        return NULL;
    }
}

        bool 
SNAPFile::generateFromCommandLine(const char **args, int nArgs, int *argsConsumed, SNAPFile *snapFile, bool paired, bool isInput)
{
    snapFile->fileName = NULL;
    snapFile->secondFileName = NULL;
    snapFile->isCompressed = false;
    *argsConsumed = 0;
    snapFile->isStdio = false;

    if (0 == nArgs) {
        return false;
    }

    //
    // Check to see if this is an explicit file type.
    //
    if ('-' == args[0][0] && '\0' != args[0][1]) { // starts with - but isn't just a - (which means to use stdio without a type specifier)
        if (1 == nArgs) {
            return false;
        }

        if (!strcmp(args[1], "-")) {
            snapFile->isStdio = true;
        }

        if (!strcmp(args[0], "-fastq") || !strcmp(args[0], "-compressedFastq")) {
            if (!isInput) {
                WriteErrorMessage("%s is not a valid output file type.\n", args[0]);
				return false;
            }

            if (paired && nArgs < 3) {
                WriteErrorMessage("Expected a pair of fastQ files, but instead just got one\n");
				return false;
            }

 
            snapFile->isCompressed = !strcmp(args[0], "-compressedFastq");

            if (paired) {
				if (nArgs < 3) {
					WriteErrorMessage("paired FASTQ requires two consecutive input files, and the last item on your command line is the first half of a FASTQ pair.\n");
					return false;
				}
                snapFile->fileType = FASTQFile;
                snapFile->secondFileName = args[2];
                if (!strcmp("-", args[2])) {
                    if (snapFile->isStdio) {
                        WriteErrorMessage("Can't have both halves of paired FASTQ files be stdin ('-').  Did you mean to use the interleaved FASTQ type?\n");
						return false;
                    }
                    snapFile->isStdio = true;
                }
                *argsConsumed = 3;
            } else {
               snapFile->fileType = FASTQFile;
               *argsConsumed = 2;
            }
        } else if (!strcmp(args[0], "-sam")) {
            snapFile->fileType = SAMFile;
            *argsConsumed = 2;
		} else if (!strcmp(args[0], "-samNoSQ") && !isInput) {	// No header is only valid for output file types
			snapFile->fileType = SAMFile;
			snapFile->omitSQLines = true;
			*argsConsumed = 2;
		} else if (!strcmp(args[0], "-bam")) {
            snapFile->fileType = BAMFile;
            snapFile->isCompressed = true;
            *argsConsumed = 2;
        } else if (!strcmp(args[0], "-pairedInterleavedFastq") || !strcmp(args[0], "-pairedCompressedInterleavedFastq")) {
            if (!paired) {
                WriteErrorMessage("Specified %s for a single-end alignment.  To treat it as single-end, just use ordinary fastq (or compressed fastq, as appropriate)\n", args[0]);
				return false;
            }

            snapFile->fileType = InterleavedFASTQFile;
            snapFile->isCompressed = !strcmp(args[0], "-pairedCompressedInterleavedFastq");
            *argsConsumed = 2;
        } else {
            //
            // starts with '-' but isn't a file-type specifier
            //
            return false;
        }
        snapFile->fileName = args[1];
        return true;
    }

    //
    // Just a filename.  Infer the type.
    //

    *argsConsumed = 1;
    snapFile->fileName = args[0];
    snapFile->isStdio = '-' == args[0][0] && '\0' == args[0][1];

    if (util::stringEndsWith(args[0], ".sam")) {
        snapFile->fileType = SAMFile;
        snapFile->isCompressed = false;
    } else if (util::stringEndsWith(args[0], ".bam")) {
        snapFile->fileType = BAMFile;
        snapFile->isCompressed = true;
    } else if (!isInput) {
        //
        // No default output file type.
        //
        WriteErrorMessage("You specified an output file with name '%s', which doesn't end in .sam or .bam, and doesn't have an explicit type\n"
                          "specifier.  There is no default output file type.  Consider doing something like '-o -bam %s'\n", args[0], args[0]);
		return false;
    } else if (util::stringEndsWith(args[0], ".fq") || util::stringEndsWith(args[0], ".fastq") ||
        util::stringEndsWith(args[0], ".fq.gz") || util::stringEndsWith(args[0], ".fastq.gz") ||
        util::stringEndsWith(args[0], ".fq.gzip") || util::stringEndsWith(args[0], ".fastq.gzip")) {

        // 
        // It's a fastq input file (either by default or because it's got a .fq or .fastq extension, we don't
        // need to check).  See if it's also compressed.
        //
        snapFile->fileType= FASTQFile;
        if (util::stringEndsWith(args[0], ".gz") || util::stringEndsWith(args[0], ".gzip")) {
            snapFile->isCompressed = true;
        } else {
            snapFile->isCompressed = false;
        }

        snapFile->isStdio = !strcmp(args[0], "-");

        if (paired) {
			if (nArgs < 2) {
				WriteErrorMessage("paired FASTQ requires two input files, and the last item on your command line is the first half of a FASTQ pair.\n");
				return false;
			}
			snapFile->secondFileName = args[1];
            if (!strcmp(args[1], "-")) {
                if (snapFile->isStdio) {
                    WriteErrorMessage("Can't have both halves of paired FASTQ files be stdin ('-').  Did you mean to use the interleaved FASTQ type?\n");
					return false;
                }
				if (CommandPipe != NULL) {
					WriteErrorMessage("You may not write to stdout in daemon mode\n");
					return false;
				}
                snapFile->isStdio = true;
            }

            *argsConsumed = 2;
        }
    } else {
        if (snapFile->isStdio) {
            WriteErrorMessage("Stdio IO always requires an explicit file type.  So, for example, do 'snap-aligner single index-directory -fastq -' to read FASTQ from stdin\n");
        } else {
            WriteErrorMessage("Unknown file type for file name '%s', please specify file type with -fastq, -sam, -bam, etc.\n", snapFile->fileName);
        }
		return false;
    }

    return true;
}

    bool         
AlignerOptions::useHadoopErrorMessages= false;

    bool
AlignerOptions::outputToStdout = false;

