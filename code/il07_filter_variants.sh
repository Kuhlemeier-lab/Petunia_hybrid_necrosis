#!/bin/bash

#SBATCH --mail-user=xxx

# Mail on NONE, BEGIN, END, FAIL, REQUEUE, ALL
#SBATCH --mail-type=ALL

#SBATCH --job-name="7filt"

# Runtime and memory
#SBATCH --time=02:10:00  # 4h
#SBATCH --mem-per-cpu=8G  # 6Gb
#SBATCH --cpus-per-task=4

# Partition
#SBATCH --partition=epyc2
#SBATCH --account=xxx
#SBATCH --chdir=/xxx/necrosis
#SBATCH --output=code/il07_filter_variants_%A.out
#SBATCH --error=code/il07_filter_variants_%A.err
#################################
chdir=/xxx/necrosis
scdir=/xxx/necrosis

module load vital-it
module add UHTS/Analysis/GenomeAnalysisTK/4.0.4.0
module add UHTS/Analysis/picard-tools/2.18.11

module list 2>&1

echo -e "#### Filter variants for quality and select
 those that match the phenotype observed.

# Apply hard filters to SNPs and INDELs separately,
#  extract variants following the phenotype,
#  summarise variant counts in file summary.txt.
#  `date`
#  current dir: `pwd`
###\t####\n
"

genome=${scdir}/data/genomes/Peax403.fasta

# Filter variants
echo -e "## Hard filter SNP variants 
"

# note that the online docs for gatk 4.0.4.0 report different
# argument names (filterExpression) but in fact they are wrong.
GenomeAnalysisTK VariantFiltration \
   -R ${genome} \
   -V ${scdir}/data/raw/il_variants/raw_snp.vcf \
   --filter-expression "QD<2.0" \
   --filter-name "QD" \
   --filter-expression "FS>60.0" \
   --filter-name "FS" \
   --filter-expression "SOR>3.0" \
   --filter-name "SOR" \
   --filter-expression "MQ<40.0" \
   --filter-name "MQ" \
   --filter-expression "MQRankSum<-12.5" \
   --filter-name "MQRankSum" \
   --filter-expression "ReadPosRankSum<-8.0" \
   --filter-name "ReadPosRankSum" \
   -O ${scdir}/data/raw/il_variants/filtered_snp.vcf

echo -e "Hard filter done.
Select variants...
"

GenomeAnalysisTK SelectVariants \
   -R ${genome} \
   -V ${scdir}/data/raw/il_variants/filtered_snp.vcf \
   --exclude-filtered true \
   -O ${scdir}/data/raw/il_variants/selected_snp.vcf

echo -e "Select done.
## Hard filter INDEL variants
"

GenomeAnalysisTK VariantFiltration \
   -R ${genome} \
   -V ${scdir}/data/raw/il_variants/raw_indel.vcf \
   --filter-expression "QD<2.0" \
   --filter-name "QD" \
   --filter-expression "FS>200.0" \
   --filter-name "FS" \
   --filter-expression "SOR>10.0" \
   --filter-name "SOR" \
   --filter-expression "ReadPosRankSum<-20.0" \
   --filter-name "ReadPosRankSum" \
   -O ${scdir}/data/raw/il_variants/filtered_indel.vcf

echo -e "Hard filter done.
Select variants...
"

GenomeAnalysisTK SelectVariants \
   -R ${genome} \
   -V ${scdir}/data/raw/il_variants/filtered_indel.vcf \
   --exclude-filtered true \
   -O ${scdir}/data/raw/il_variants/selected_indel.vcf

# merge SNP and INDEL variants in a single file
echo -e "Selection done.
Merging SNPs and INDELs...
"

picard-tools MergeVcfs \
    I=${scdir}/data/raw/il_variants/selected_snp.vcf \
    I=${scdir}/data/raw/il_variants/selected_indel.vcf \
    O=${scdir}/data/raw/il_variants/clean.vcf

echo -e "Hard filter, select, merge SNP and INDEL done
"
