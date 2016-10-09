del samtoolsForVcf.cmd
copy \\gcr\scratch\b99\bolosky\samtools.exe .
\\gcr\scratch\b99\bolosky\GenerateScriptFromVariants %1 %2 samtoolsForVcf.cmd %3
\\gcr\scratch\b99\bolosky\GenerateConsolodatedExtractedReads samtoolsForVcf.cmd %4 0 .\samtools.exe
del samtoolsForVcf.cmd
del samtools.exe

