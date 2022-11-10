date
mkdir /mnt/d/temp/HaplotypeCallerTemp/simulated_hg3_$1x
rm /mnt/d/temp/HaplotypeCallerTemp/simulated_hg3_$1x/*
cd ~/gatk
cat /mnt/d/gdc/chromosomes.txt | parallel -j 12 ./gatk HaplotypeCaller -R /mnt/d/sequence/genomes/Homo_sapiens_assembly38.fasta -I /mnt/d/temp/simulated_hg3_sampled_$1x.bam -O /mnt/d/temp/HaplotypeCallerTemp/simulated_hg3_$1x/simulated_hg3_$1x.{}.g.vcf.gz -L {} -ERC GVCF -pcr-indel-model NONE --dont-use-soft-clipped-bases true --dbsnp /mnt/d/sequence/genomes/Homo_sapiens_assembly38.dbsnp138.vcf --java-options '-DGATK_STACKTRACE_ON_USER_EXCEPTION=true'
bcftools concat -o /mnt/d/temp/simulated_hg3_$1x.vcf /mnt/d/temp/HaplotypeCallerTemp/simulated_hg3_$1x/simulated_hg3_$1x.*.g.vcf.gz 
rm -rf /mnt/d/temp/HaplotypeCallerTemp/simulated_hg3_$1x
date
