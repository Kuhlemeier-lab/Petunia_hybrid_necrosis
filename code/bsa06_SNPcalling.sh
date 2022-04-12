#!/bin/bash

#SBATCH --mail-user=xxx
#SBATCH --mail-type=ALL
#SBATCH --job-name="callvar"

# Runtime and memory
#SBATCH --time=48:00:00  # 22h
#SBATCH --mem-per-cpu=8G  # 4.5Gb
#SBATCH --cpus-per-task=4

# Partition
#SBATCH --partition=epyc2

#SBATCH --account=xxx
#SBATCH --chdir=/xxx/necrosis
#SBATCH --output=code/bsa06_SNPcalling_%A.out
#SBATCH --error=code/bsa06_SNPcalling_%A.err

#################################
chdir=/xxx/necrosis
scdir=/xxx/necrosis

if [ ! -d ${scdir}/data/raw/bsa_variants ]; then
  mkdir -p ${scdir}/data/raw/bsa_variants
fi

module add vital-it
module add UHTS/Analysis/GenomeAnalysisTK/4.0.4.0
module add UHTS/Analysis/vcftools/0.1.15

echo -e "## `date`\n## Call variants.\n"

GenomeAnalysisTK HaplotypeCaller \
    -R ${scdir}/data/genomes/Peax402INV.fasta \
    -I ${scdir}/data/raw/bsa_aligned_reads/SRR13907523.sp.bam \
    -I ${scdir}/data/raw/bsa_aligned_reads/SRR13907524.sp.bam \
    --dont-use-soft-clipped-bases \
    --pcr-indel-model NONE \
    -stand-call-conf 30 \
    --native-pair-hmm-threads 4 \
    -O ${scdir}/data/raw/bsa_variants/BSA_raw.vcf 
