@rem total number of reads in simulator output
set /a total=229111160*2
@rem n reads in this sample.  Divide before multiply to avoid int32 overflow (at the loss of a tiny bit of precision)
set /a n=%total%/30*%1
cd /d d:\temp
\gdc\bin\SampleFastqWithKnownSize %n% %total% d:\sequence\InfoGainSim\simulated_hg3_200bp_paired_30x_1.fastq d:\sequence\InfoGainSim\simulated_hg3_200bp_paired_30x_2.fastq d:\temp\simulated_hg3_sampled_%1x_1.fq d:\temp\simulated_hg3_sampled_%1x_2.fq
\gdc\bin\snap.exe paired d:\sequence\indices\Homo_sapiens_assembly38-24-large-liftover d:\temp\simulated_hg3_sampled_%1x_1.fq d:\temp\simulated_hg3_sampled_%1x_2.fq -so -sm 20 -o d:\temp\simulated_hg3_sampled_%1x.bam

wsl /mnt/d/gdc/bin/run_haplotype_caller_sampled %1

wsl /mnt/d/gdc/bin/run_happy_mason_sampled %1
