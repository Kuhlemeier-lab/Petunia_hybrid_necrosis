#!/bin/bash

#SBATCH --mail-user=xxx
#SBATCH --mail-type=ALL
#SBATCH --job-name="7plotsnp"

# Runtime and memory
#SBATCH --time=24:00:00  # 2'
#SBATCH --mem-per-cpu=6G  # 0.8Gb
#SBATCH --cpus-per-task=1

# Partition
#SBATCH --partition=epyc2

#SBATCH --account=xxx
#SBATCH --chdir=/xxx/necrosis
#SBATCH --output=code/bsa07_SNPselect_qualityplot_%A.out
#SBATCH --error=code/bsa07_SNPselect_qualityplot_%A.err

#################################
chdir=/xxx/necrosis
scdir=/xxx/necrosis

module add vital-it
module add UHTS/Analysis/GenomeAnalysisTK/4.0.4.0
module load R/3.4.2

echo -e "## `date`\n## Select SNP variants.\n"

GenomeAnalysisTK SelectVariants \
    -R ${scdir}/data/genomes/Peax403.fasta \
    -V ${scdir}/data/raw/bsa_variants/BSA_raw.vcf \
    -O ${scdir}/data/raw/bsa_variants/BSA_SNP_raw.vcf \
    --select-type-to-include SNP

echo -e "\nExtracting a table for quality parameters plotting.\n"
GenomeAnalysisTK VariantsToTable \
     -R ${scdir}/data/genomes/Peax403.fasta \
     -V ${scdir}/data/raw/bsa_variants/BSA_SNP_raw.vcf \
     -O ${scdir}/data/raw/bsa_variants/BSA_SNP_raw.vcf.table \
     -F CHROM \
     -F POS \
     -F QD \
     -F FS \
     -F SOR \
     -F MQ \
     -F MQRankSum \
     -F ReadPosRankSum \
     --show-filtered

echo -e "## `date`\n## Plot quality of raw SNP.\n"

${chdir}/code/plot_vcfq_distribution.R -i ${scdir}/data/raw/bsa_variants/BSA_SNP_raw.vcf.table --vcftype SNP -o ${chdir}/figures/exploratory/bsa_variants
