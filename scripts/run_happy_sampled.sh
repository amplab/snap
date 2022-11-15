date
bcftools annotate -x FORMAT/AD -O z -o /mnt/d/temp/hg3_neat_shortreads_sampled_$1x.NoAd.vcf.gz /mnt/d/temp/hg3_neat_shortreads_sampled_$1x.vcf
bcftools view -c 1 -O z -o /mnt/d/temp/hg3_neat_shortreads_sampled_$1x.NoAd.NoHomRef.vcf.gz /mnt/d/temp/hg3_neat_shortreads_sampled_$1x.NoAd.vcf.gz
cd ~/hap.py
mkdir /mnt/d/temp/output.hg3_neat_shortreads_sampled_$1x
python hap.py-install/bin/hap.py /mnt/d/sequence/ground-truth-renamed/hg003.vcf.gz /mnt/d/temp/hg3_neat_shortreads_sampled_$1x.NoAd.NoHomRef.vcf.gz -r /mnt/d/sequence/genomes/Homo_sapiens_assembly38.fasta -o /mnt/d/temp/output.hg3_neat_shortreads_sampled_$1x/simulated_hg3_$1x.concordance --engine=vcfeval --roc QUAL
cd /mnt/d/temp/output.hg3_neat_shortreads_sampled_$1x
tar cvf /mnt/d/temp/shg3_neat_shortreads_sampled_$1x.concordance.tar .
cd ~
rm -rf /mnt/d/temp/hg3_neat_shortreads_sampled_$1x.output /mnt/d/temp/hg3_neat_shortreads_sampled_$1x.NoAd.vcf.gz /mnt/d/temp/hg3_neat_shortreads_sampled_$1x.NoAd.NoHomRef.vcf.gz /mnt/d/temp/output.hg3_neat_shortreads_sampled_$1x
