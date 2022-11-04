rem first arg is long read coverage, second is short

rd /s /q d:\temp\InfogainTemp
md d:\temp\InfogainTemp
cd /d d:\temp\InfogainTemp
d:\gdc\bin\SimulateLongReadsAsFASTA.exe d:\sequence\InfoGainSim\hg003.fasta d:\temp\InfogainTemp\longreads_%1x 10000 %1 -n

wsl /mnt/d/gdc/bin/minimap2 -t 16 -a /mnt/d/sequence/indices/Homo_sapiens_assembly38-long.mmi /mnt/d/temp/InfoGainTemp/longreads_%1x.fastq -o /mnt/d/temp/InfogainTemp/longreads_%1x.sam

d:\gdc\bin\FASTAFromSAM.exe longreads_%1x.sam longread_

for %%r in (longread_*fasta) do call d:\temp\run_infogain_alignment %%r %2


cd /d d:\temp
rd /s /q d:\temp\InfogainTemp