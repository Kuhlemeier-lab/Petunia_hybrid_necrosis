#!/bin/bash

#SBATCH --mail-user=xxx

# Mail on NONE, BEGIN, END, FAIL, REQUEUE, ALL
#SBATCH --mail-type=ALL

#SBATCH --job-name="3align"

# Runtime and memory
#SBATCH --time=08:00:00  #  6h
#SBATCH --mem-per-cpu=8G  #  20Gb total
#SBATCH --cpus-per-task=4

# Array
#SBATCH --array=23,24

# Partition
#SBATCH --partition=epyc2

#SBATCH --account=xxx
#SBATCH --chdir=/xxx/necrosis
#SBATCH --output=code/bsa03_STARmapping_%A_%a.out
#SBATCH --error=code/bsa03_STARmapping_%A_%a.err

#################################
chdir=/xxx/necrosis
scdir=/xxx/necrosis

# Create a data/raw/aligned_reads folder if not already existing

if [ ! -d ${scdir}/data/raw/bsa_aligned_reads ]; then
	mkdir -p ${scdir}/data/raw/bsa_aligned_reads
fi

module load vital-it/7
module load UHTS/Aligner/STAR/2.6.0c

echo -e "## `date`\n## Reads mapping...\n "

STAR --genomeDir ${scdir}/data/genomes \
  --runThreadN 4 \
  --readFilesIn ${scdir}/data/raw/bsa_rawreads/SRR139075${SLURM_ARRAY_TASK_ID}_1_trim.fq.gz ${scdir}/data/raw/bsa_rawreads/SRR139075${SLURM_ARRAY_TASK_ID}_2_trim.fq.gz \
  --sjdbOverhang 149 \
  --outFilterType BySJout \
  --outFilterMultimapNmax 20 \
  --outSAMtype BAM SortedByCoordinate \
  --twopassMode Basic --outFileNamePrefix ${scdir}/data/raw/bsa_aligned_reads/SRR139075${SLURM_ARRAY_TASK_ID}_ \
  --readFilesCommand zcat --genomeLoad NoSharedMemory

echo -e "Reads mapping done.\n `date`"

