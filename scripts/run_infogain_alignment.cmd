rd /s /q snapTempIndex
d:\gdc\bin\snap.exe index %1 snapTempIndex -qq

@rem mason takes input in pairs.  So total bases (i.e., read size of 10K) times short read coverage divided by read *pair* size
set /a nPairs = 10000 * %2 / 400
wsl /mnt/d/gdc/bin/mason_simulator -ir /mnt/d/temp/InfogainTemp/%1 -n %npairs% -o /mnt/d/temp/InfogainTemp/shortreads_1.fastq -or /mnt/d/temp/InfogainTemp/shortreads_2.fastq --num-threads 16 --illumina-read-length 200

d:\gdc\bin\snap.exe paired snapTempIndex shortreads_1.fastq shortreads_2.fastq -o %1.sam -qq











