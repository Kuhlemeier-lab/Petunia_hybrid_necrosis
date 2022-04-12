#!/bin/bash

#SBATCH --mail-user=xxx
#SBATCH --mail-type=ALL
#SBATCH --job-name="markdup"

# Runtime and memory
#SBATCH --time=24:00:00  #  
#SBATCH --mem-per-cpu=12G  # .

#SBATCH --array=23,24

# Partition
#SBATCH --partition=epyc2

#SBATCH --account=xxx
#SBATCH --chdir=/xxx/necrosis
#SBATCH --output=code/bsa05_markdup_%A_%a.out
#SBATCH --error=code/bsa05_markdup_%A_%a.err

#################################
chdir=/xxx/necrosis
scdir=/xxx/necrosis

module add vital-it
module add UHTS/Analysis/picard-tools/2.18.11
module add UHTS/Analysis/GenomeAnalysisTK/4.0.4.0

echo -e "## `date`\n## Mark duplicated reads and split N cigar reads..."

cd ${scdir}/data/raw/bsa_aligned_reads

picard-tools MarkDuplicates \
    I=SRR139075${SLURM_ARRAY_TASK_ID}_Aligned.sortedByCoord_readgroup.out.bam \
    O=SRR139075${SLURM_ARRAY_TASK_ID}.md.bam \
    CREATE_INDEX=true \
    VALIDATION_STRINGENCY=STRICT \
    M=SRR139075${SLURM_ARRAY_TASK_ID}.output.metrics \
    TMP_DIR=${scdir}/data/raw/bsa_aligned_reads

GenomeAnalysisTK SplitNCigarReads \
    -R ${scdir}/data/genomes/Peax402INV.fasta \
    -I SRR139075${SLURM_ARRAY_TASK_ID}.md.bam \
    -O SRR139075${SLURM_ARRAY_TASK_ID}.sp.bam \
    --TMP_DIR ${scdir}/data/raw/bsa_aligned_reads

echo -e " Done!!.\n `date`"
