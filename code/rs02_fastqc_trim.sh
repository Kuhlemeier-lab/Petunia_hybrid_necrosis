#!/bin/bash

#SBATCH --mail-user=xxx

# Mail on NONE, BEGIN, END, FAIL, REQUEUE, ALL
#SBATCH --mail-type=ALL

#SBATCH --job-name="02fastq"

# Runtime and memory
#SBATCH --time=02:00:00  # 24'
#SBATCH --mem-per-cpu=4G  # 2Gb
#SBATCH --cpus-per-task=4

# Array
#SBATCH --array=1-9

# Partition
#SBATCH --partition=epyc2

#SBATCH --account=xxx
#SBATCH --chdir=/xxx/necrosis
#SBATCH --output=code/rs02_fastqc_trim_%A_%a.out
#SBATCH --error=code/rs02_fastqc_trim_%A_%a.err

#################################
chdir=/xxx/necrosis
scdir=/xxx/necrosis

echo -e "#### Quality of raw reads
## `date`
## current dir: `pwd`
####\t####\n"

module load vital-it
module load UHTS/Quality_control/fastqc/0.11.7
module load UHTS/Analysis/trimmomatic/0.36

# assess quality of each raw read file with fastqc
echo -e "\n\n## Assessing fastq quality..."
cd ${scdir}/data/raw/rs_rawreads

fastqc -t 4 kmh${SLURM_ARRAY_TASK_ID}.fastq.gz

echo -e "Assessing quality done.\n Trim reads."

trimmomatic SE \
  -threads 4 \
  -phred33 \
  kmh${SLURM_ARRAY_TASK_ID}.fastq.gz \
  kmh${SLURM_ARRAY_TASK_ID}_trim.fastq.gz \
  LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36 \
  ILLUMINACLIP:/software/UHTS/Analysis/trimmomatic/0.36/adapters/TruSeq3-SE.fa:2:30:3

echo -e "Trimming done.\n `date`"

# assess quality of each raw read file with fastqc
echo -e "\n\n## Assessing fastq quality..."

fastqc -t 4 kmh${SLURM_ARRAY_TASK_ID}_trim.fastq.gz

echo -e "Assessing quality of trimmed reads done.\n"
