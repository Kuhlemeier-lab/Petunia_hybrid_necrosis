#!/bin/bash

#SBATCH --mail-user=xxx

# Mail on NONE, BEGIN, END, FAIL, REQUEUE, ALL
#SBATCH --mail-type=ALL

#SBATCH --job-name="01fastq"

# Runtime and memory
#SBATCH --time=08:00:00  #
#SBATCH --mem-per-cpu=4Gb  # 
#SBATCH --cpus-per-task=4

#SBATCH --array=23,24

# Partition
#SBATCH --partition=epyc2

#SBATCH --account=xxx
#SBATCH --chdir=/xxx/necrosis
#SBATCH --output=code/bsa01_quality_trim_%A_%a.out
#SBATCH --error=code/bsa01_quality_trim_%A_%a.err

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
cd ${scdir}/data/raw/bsa_rawreads

fastqc -t 4 SRR*${SLURM_ARRAY_TASK_ID}_*.fq.gz
echo -e "Assessing quality done.\n Trim reads."

trimmomatic PE \
  -threads 4 \
  -phred33 \
  SRR139075${SLURM_ARRAY_TASK_ID}_1.fq.gz \
  SRR139075${SLURM_ARRAY_TASK_ID}_2.fq.gz \
  SRR139075${SLURM_ARRAY_TASK_ID}_1_trim.fq.gz \
  SRR139075${SLURM_ARRAY_TASK_ID}_1_unpaired.fq.gz \
  SRR139075${SLURM_ARRAY_TASK_ID}_2_trim.fq.gz \
  SRR139075${SLURM_ARRAY_TASK_ID}_2_unpaired.fq.gz \
  LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36
echo -e "Trimming done.\n `date`"

echo -e "#### Quality of trimmed reads
## `date`
## current dir: `pwd`
####\t####\n"

# assess quality of each raw read file with fastqc
echo -e "\n\n## Assessing fastq quality..."

for i in {1..2}
do
  fastqc -t 4 SRR139075${SLURM_ARRAY_TASK_ID}_${i}_trim.fq.gz
done

