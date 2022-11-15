@rem total number of reads in simulator output
set /a total=460767070
@rem n reads in this sample.  Divide before multiply to avoid int32 overflow (at the loss of a tiny bit of precision)
set /a n=%total%/30*%1
cd /d d:\temp
\gdc\bin\SampleFastqWithKnownSize %n% %total% d:\sequence\InfoGainSim\hg3_neat_shortreads_30x_fixed_read1.fq d:\sequence\InfoGainSim\hg3_neat_shortreads_30x_fixed_read2.fq d:\temp\hg3_neat_shortreads_sampled_%1x_1.fq d:\temp\hg3_neat_shortreads_sampled_%1x_2.fq
\gdc\bin\snap.exe paired d:\sequence\indices\Homo_sapiens_assembly38-24-large-liftover d:\temp\hg3_neat_shortreads_sampled_%1x_1.fq d:\temp\hg3_neat_shortreads_sampled_%1x_2.fq -so -sm 20 -o d:\temp\hg3_neat_shortreads_sampled_%1x.bam

wsl /mnt/d/gdc/bin/run_haplotype_caller_sampled.sh %1

wsl /mnt/d/gdc/bin/run_happy_sampled.sh %1
