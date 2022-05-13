#!/bin/bash

#SBATCH --mail-user=xxx

# Mail on NONE, BEGIN, END, FAIL, REQUEUE, ALL
#SBATCH --mail-type=ALL

#SBATCH --job-name="6gvcf"

# Runtime and memory
#SBATCH --time=24:00:00  # 10h 
#SBATCH --mem-per-cpu=8G  # 12Gb 
#SBATCH --cpus-per-task=6

# Partition
#SBATCH --partition=epyc2
#SBATCH --account=ips_ck
#SBATCH --chdir=/xxx/necrosis
#SBATCH --output=code/il06_combineGvcf_%A.out
#SBATCH --error=code/il06_combineGvcf_%A.err
#################################
chdir=/xxx/necrosis
scdir=/xxx/necrosis

echo -e "#### Combine and then genotype gvcf files
## `date`
## current dir: `pwd`
####\t####\n"

module load vital-it
module add UHTS/Analysis/GenomeAnalysisTK/4.0.4.0
module add UHTS/Analysis/picard-tools/2.18.11
module add R/3.4.2

module list 2>&1

echo -e "\n## Combine GVCFs\n"

genome=${scdir}/data/genomes/Peax402INV.fasta

GenomeAnalysisTK CombineGVCFs \
     -R ${genome} \
     --variant ${scdir}/data/raw/il_variants/kmh1.g.vcf \
     --variant ${scdir}/data/raw/il_variants/kmh2.g.vcf \
     --variant ${scdir}/data/raw/il_variants/kmh3.g.vcf \
     --variant ${scdir}/data/raw/il_variants/kmh4.g.vcf \
     --variant ${scdir}/data/raw/il_variants/kmh5.g.vcf \
     --variant ${scdir}/data/raw/il_variants/kmh6.g.vcf \
     --variant ${scdir}/data/raw/il_variants/kmh7.g.vcf \
     --variant ${scdir}/data/raw/il_variants/kmh8.g.vcf \
     --variant ${scdir}/data/raw/il_variants/kmh9.g.vcf \
     --variant ${scdir}/data/raw/il_variants/kmh10.g.vcf \
     --variant ${scdir}/data/raw/il_variants/kmh11.g.vcf \
     --variant ${scdir}/data/raw/il_variants/kmh12.g.vcf \
     --variant ${scdir}/data/raw/il_variants/kmh13.g.vcf \
     --variant ${scdir}/data/raw/il_variants/kmh14.g.vcf \
     -O ${scdir}/data/raw/il_variants/cohort.g.vcf

echo -e "\n## Combine done.\nGenotype GVCF\n"

GenomeAnalysisTK GenotypeGVCFs \
     -R ${genome} \
     -V ${scdir}/data/raw/il_variants/cohort.g.vcf \
     -O ${scdir}/data/raw/il_variants/raw.vcf \
     --seconds-between-progress-updates 600

# remove cohort.g.vcf
rm ${scdir}/data/raw/il_variants/cohort.g.vcf ${scdir}/data/raw/il_variants/cohort.g.vcf.idx

echo -e "\nCombine and genotype done.\n"

echo -e "Separate SNPs from INDELs.\n"

picard-tools SplitVcfs \
    I=${scdir}/data/raw/il_variants/raw.vcf \
    SNP_OUTPUT=${scdir}/data/raw/il_variants/raw_snp.vcf \
    INDEL_OUTPUT=${scdir}/data/raw/il_variants/raw_indel.vcf \
    STRICT=false

echo -e "\nExtracting a table for quality parameters plotting.\n"
GenomeAnalysisTK VariantsToTable \
     -R ${genome} \
     -V ${scdir}/data/raw/il_variants/raw_snp.vcf \
     -O ${scdir}/data/raw/il_variants/raw_snp.vcf.table \
     -F CHROM \
     -F POS \
     -F QD \
     -F FS \
     -F SOR \
     -F MQ \
     -F MQRankSum \
     -F ReadPosRankSum \
     --show-filtered

GenomeAnalysisTK  VariantsToTable \
    -R ${genome} \
    -V ${scdir}/data/raw/il_variants/raw_indel.vcf \
    -O ${scdir}/data/raw/il_variants/raw_indel.vcf.table \
    -F CHROM \
    -F POS \
    -F QD \
    -F FS \
    -F SOR \
    -F ReadPosRankSum \
    --show-filtered

echo -e "\nraw_snp.vcf.table and raw_indel.vcf.table were saved\n"

# call plot_distribution.R
echo -e "Make distribution plots of SNPs and INDELs.\n"

if [ ! -d ${chdir}/figures/exploratory/il_variants ]; then
	mkdir -p ${chdir}/figures/exploratory/il_variants
fi
# with GATK 4.1.3 I get a trailing tab on the header line of the table,
# so I remove it beforehand to avoid problems in R.
sed -i "s/\t$//" ${scdir}/data/raw/il_variants/raw_indel.vcf.table
sed -i "s/\t$//" ${scdir}/data/raw/il_variants/raw_snp.vcf.table
# requires optparse library
Rscript ${chdir}/code/plot_vcfq_distribution.R -i ${scdir}/data/raw/il_variants/raw_indel.vcf.table --vcftype INDEL -o ${chdir}/figures/exploratory/il_variants
Rscript ${chdir}/code/plot_vcfq_distribution.R -i ${scdir}/data/raw/il_variants/raw_snp.vcf.table --vcftype SNP -o ${chdir}/figures/exploratory/il_variants
echo -e "Plots done. They are saved in figures/exploratory/il_variants .\n"

rm ${scdir}/data/raw/il_variants/raw_snp.vcf.table
rm ${scdir}/data/raw/il_variants/raw_indel.vcf.table
