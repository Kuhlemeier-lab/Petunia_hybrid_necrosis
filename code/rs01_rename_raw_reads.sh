#!/bin/bash

#SBATCH --mail-user=xxx

# Mail on NONE, BEGIN, END, FAIL, REQUEUE, ALL
#SBATCH --mail-type=ALL

#SBATCH --job-name="1renameraw"

# Runtime and memory
#SBATCH --time=02:00:00  # 10', time it actually took
#SBATCH --mem-per-cpu=4G  # 520Kb, memory it actually used
#SBATCH --cpus-per-task=4

# Partition
#SBATCH --partition=epyc2

#SBATCH --account=xxx
#SBATCH --chdir=/xxx/necrosis
#SBATCH --output=code/rs01_rename_raw_reads_%A.out
#SBATCH --error=code/rs01_rename_raw_reads_%A.err

#################################
chdir=/xxx/necrosis
scdir=/xxx/necrosis

echo -e "#### Rename raw reads
## `date`
## current dir: `pwd`
####\t####\n"

# Move to data/raw/raw_reads
cd ${scdir}/data/raw/rs_rawreads

# rename raw reads files to cut the useless part of the file name and allow for job array

mv SRR13861746_1.fastq.gz kmh1.fastq.gz
mv SRR13861745_1.fastq.gz kmh2.fastq.gz
mv SRR13861744_1.fastq.gz kmh3.fastq.gz
mv SRR13861743_1.fastq.gz kmh4.fastq.gz
mv SRR13861742_1.fastq.gz kmh5.fastq.gz
mv SRR13861741_1.fastq.gz kmh6.fastq.gz
mv SRR13861740_1.fastq.gz kmh7.fastq.gz
mv SRR13861739_1.fastq.gz kmh8.fastq.gz
mv SRR13861738_1.fastq.gz kmh9.fastq.gz

echo -e "Done."
