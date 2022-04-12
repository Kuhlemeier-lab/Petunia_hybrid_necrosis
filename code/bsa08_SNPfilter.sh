#!/bin/bash

#SBATCH --mail-user=xxx
#SBATCH --mail-type=ALL
#SBATCH --job-name="9filter"

# Runtime and memory
#SBATCH --time=12:00:00  # 4'
#SBATCH --mem-per-cpu=6G  # 2.7Gb
#SBATCH --cpus-per-task=4

# Partition
#SBATCH --partition=epyc2

#SBATCH --account=xxx
#SBATCH --chdir=/xxx/necrosis
#SBATCH --output=code/bsa08_SNPfilter_%A.out
#SBATCH --error=code/bsa08_SNPfilter_%A.err

#################################
chdir=/xxx/necrosis
scdir=/xxx/necrosis

module add vital-it
module add UHTS/Analysis/GenomeAnalysisTK/4.0.4.0
module add UHTS/Analysis/vcftools/0.1.15

echo -e "## `date`\n## Filter variants."

GenomeAnalysisTK VariantFiltration \
    -R ${scdir}/data/genomes/Peax402INV.fasta \
    -V ${scdir}/data/raw/bsa_variants/BSA_SNP_raw.vcf \
    -window 10 \
    -cluster 3 \
    --filter-expression "QD<2.0" \
    --filter-name "QD" \
    --filter-expression "FS>30.0" \
    --filter-name "FS" \
    -O ${scdir}/data/raw/bsa_variants/BSA_SNP_raw_gatkfiltered.vcf

GenomeAnalysisTK SelectVariants \
    -R ${scdir}/data/genomes/Peax402INV.fasta \
    -V ${scdir}/data/raw/bsa_variants/BSA_SNP_raw_gatkfiltered.vcf \
    --restrict-alleles-to BIALLELIC \
    --select-type-to-include SNP \
    --exclude-filtered \
    -O ${scdir}/data/raw/bsa_variants/BSA_SNP_biallelic_gatkselected.vcf

echo -e "## `date`\n## Done."

# filter for depth at least 100 and thin for 100bp
vcftools --vcf ${scdir}/data/raw/bsa_variants/BSA_SNP_biallelic_gatkselected.vcf \
  --minDP 100 \
  --max-missing 1 \
  --remove-filtered-geno-all \
  --recode \
  --out ${scdir}/data/raw/bsa_variants/BSA_SNP_biallelic_gatkselected_minDP100

vcftools --vcf ${scdir}/data/raw/bsa_variants/BSA_SNP_biallelic_gatkselected_minDP100.recode.vcf \
  --thin 100 \
  --recode \
  --out ${scdir}/data/raw/bsa_variants/BSA_SNP_biallelic_gatkselected_minDP100_thin100

