# mkdir ./gatk_tmp
# gatk CreateSequenceDictionary -R /home/wayne/Project/SC/Sanger/0911/HA-1_res/ref/HA-1-2-C07_110.fa
# samtools faidx /home/wayne/Project/SC/Sanger/0911/HA-1_res/ref/HA-1-2-C07_110.fa
bwa mem -t 24 /home/wayne/Project/SC/Sanger/0911/HA-1_res/ref/HA-1-2-C07_110.fa 0911/HA-1_res/fq/HA-1-2-C07_110.fq -R '@RG\tID:2T\tPL:illumina\tPU:2T\tLB:2T\tSM:2T\tCN:VectorBuilder'|samtools sort -@ 24 -o BL.sort.bam -
gatk MarkDuplicates -I BL.sort.bam -O ./HA-1-2-C07_110_marked_duplicates.bam -M metrics.txt
gatk  HaplotypeCaller -R ~/Project/SC/Sanger/0911/HA-1_res/ref/HA-1-2-C07_110.fa -I ./HA-1-2-C07_110_marked_duplicates.bam -O ./gatk_tmp/BL.raw.vcf --tmp-dir .//gatk_tmp  --pcr-indel-model AGGRESSIVE --min-base-quality-score 20 --native-pair-hmm-threads 24 --sample-name 2T
# gatk GenotypeGVCFs -R ~/Project/SC/Sanger/0911/HA-1_res/ref/HA-1-2-C07_110.fa --variant  .//gatk_tmp/BL.raw.gvcf.gz -O  .//gatk_tmp/BL.variants.vcf.gz --tmp-dir .//gatk_tmp/
# gatk VariantFiltration -R ~/Project/SC/Sanger/0911/HA-1_res/ref/HA-1-2-C07_110.fa -V .//gatk_tmp/BL.variants.vcf.gz -O .//gatk_tmp/filtered.vcf
