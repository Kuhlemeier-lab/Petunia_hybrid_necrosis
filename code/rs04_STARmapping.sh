#!/bin/bash

#SBATCH --mail-user=xxx

# Mail on NONE, BEGIN, END, FAIL, REQUEUE, ALL
#SBATCH --mail-type=ALL

#SBATCH --job-name="4align"

# Runtime and memory
#SBATCH --time=02:00:00  # 10'
#SBATCH --mem-per-cpu=4G  # 26Gb
#SBATCH --cpus-per-task=16

# Array
#SBATCH --array=1-9

# Partition
#SBATCH --partition=epyc2

#SBATCH --account=xxx
#SBATCH --chdir=/xxx/necrosis
#SBATCH --output=code/rs04_STARmapping_%A_%a.out
#SBATCH --error=code/rs04_STARmapping_%A_%a.err

#################################
chdir=/xxx/necrosis
scdir=/xxx/necrosis

# Create a data/raw/aligned_reads folder if not already existing

if [ ! -d ${scdir}/data/raw/rs_aligned_reads ]; then
	mkdir -p ${scdir}/data/raw/rs_aligned_reads
fi

module load vital-it/7
module load UHTS/Aligner/STAR/2.6.0c

echo -e "## `date`\n## Reads mapping...\n "

STAR --genomeDir ${scdir}/data/genomes/Pax_125bp \
    --runThreadN 16 \
    --readFilesIn ${scdir}/data/raw/rs_rawreads/kmh${SLURM_ARRAY_TASK_ID}_trim.fastq.gz \
    --outFilterType BySJout \
    --outFilterMultimapNmax 20 \
    --outSAMtype BAM SortedByCoordinate \
    --twopassMode Basic \
    --outFileNamePrefix ${scdir}/data/raw/rs_aligned_reads/kmh${SLURM_ARRAY_TASK_ID}_trim_STAR_ --readFilesCommand zcat --genomeLoad NoSharedMemory

echo -e "Reads mapping done.\n `date`"
