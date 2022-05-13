#!/bin/bash

#SBATCH --mail-user=xxx

# Mail on NONE, BEGIN, END, FAIL, REQUEUE, ALL
#SBATCH --mail-type=ALL

#SBATCH --job-name="2fq_trim"

# Runtime and memory
#SBATCH --time=12:00:00  # 30'
#SBATCH --mem-per-cpu=8G  # 4Gb
#SBATCH --cpus-per-task=4

# Array
#SBATCH --array=1-14

# Partition
#SBATCH --partition=epyc2

#SBATCH --account=xxx
#SBATCH --chdir=/xxx/necrosis
#SBATCH --output=code/il02_fastqc_trim_%A_%a.out
#SBATCH --error=code/il02_fastqc_trim_%A_%a.err

#################################
chdir=/xxx/necrosis
scdir=/xxx/necrosis

echo -e "#### Quality and trimming of raw reads
## `date`
## current dir: `pwd`
####\t####\n"

module load vital-it
module load UHTS/Quality_control/fastqc/0.11.7 
module load UHTS/Analysis/trimmomatic/0.36 

module list 2>&1

# assess quality of each raw read file with fastqc
# R1 is the forward reads, R2 the reverse
echo -e "\n\n## Assessing fastq R1 quality..."
cd ${scdir}/data/raw/il_rawreads

fastqc -t 4 kmh${SLURM_ARRAY_TASK_ID}_R1.fastq.gz

echo -e "\n## `date`
## Assessing fastq R2 quality..."
fastqc -t 4 kmh${SLURM_ARRAY_TASK_ID}_R2.fastq.gz

echo -e "Assessing quality done.\n"

# trimming
echo -e "## `date`\n## Trimming "
# since I noted that file kmh1 has some quality issue in the last part of the reads (after 145 bp), I decided to crop it at 145.
if [ "${SLURM_ARRAY_TASK_ID}" = "1" ] # if file is kmh1 then crop the reads, otherwise don't.
then
    trimmomatic PE -threads 16 \
        -phred33 \
        kmh${SLURM_ARRAY_TASK_ID}_R1.fastq.gz kmh${SLURM_ARRAY_TASK_ID}_R2.fastq.gz kmh"${SLURM_ARRAY_TASK_ID}"_R1.fastq.gz.trimmed kmh"${SLURM_ARRAY_TASK_ID}"_R1.fastq.gz.unpaired kmh"${SLURM_ARRAY_TASK_ID}"_R2.fastq.gz.trimmed kmh"${SLURM_ARRAY_TASK_ID}"_R2.fastq.gz.unpaired \
        LEADING:3 \
        TRAILING:3 \
        SLIDINGWINDOW:4:15 \
        MINLEN:100 \
        CROP:145
else
    trimmomatic PE -threads 16 \
        -phred33 \
        kmh${SLURM_ARRAY_TASK_ID}_R1.fastq.gz kmh${SLURM_ARRAY_TASK_ID}_R2.fastq.gz kmh"${SLURM_ARRAY_TASK_ID}"_R1.fastq.gz.trimmed kmh"${SLURM_ARRAY_TASK_ID}"_R1.fastq.gz.unpaired kmh"${SLURM_ARRAY_TASK_ID}"_R2.fastq.gz.trimmed kmh"${SLURM_ARRAY_TASK_ID}"_R2.fastq.gz.unpaired \
        LEADING:3 \
        TRAILING:3 \
        SLIDINGWINDOW:4:15 \
        MINLEN:100
fi

echo -e "Trimming done.\n `date`"

fastqc -t 4 kmh${SLURM_ARRAY_TASK_ID}_R1.fastq.gz.trimmed
fastqc -t 4 kmh${SLURM_ARRAY_TASK_ID}_R2.fastq.gz.trimmed
