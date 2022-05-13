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
#SBATCH --output=code/il01_rename_raw_reads_%A.out
#SBATCH --error=code/il01_rename_raw_reads_%A.err

#################################
chdir=/xxx/necrosis
scdir=/xxx/necrosis

echo -e "#### Rename raw reads
## `date`
## current dir: `pwd`
####\t####\n"

# Move to data/raw/raw_reads
cd ${scdir}/data/raw/il_rawreads

# rename raw reads files to cut the useless part of the file name and allow for job array

mv SRR13809747_1.fastq.gz kmh1_R1.fastq.gz
mv SRR13809746_1.fastq.gz kmh2_R1.fastq.gz
mv SRR13809741_1.fastq.gz kmh3_R1.fastq.gz
mv SRR13809740_1.fastq.gz kmh4_R1.fastq.gz
mv SRR13809739_1.fastq.gz kmh5_R1.fastq.gz
mv SRR13809738_1.fastq.gz kmh6_R1.fastq.gz
mv SRR13809737_1.fastq.gz kmh7_R1.fastq.gz
mv SRR13809736_1.fastq.gz kmh8_R1.fastq.gz
mv SRR13809735_1.fastq.gz kmh9_R1.fastq.gz
mv SRR13809734_1.fastq.gz kmh10_R1.fastq.gz
mv SRR13809745_1.fastq.gz kmh11_R1.fastq.gz
mv SRR13809744_1.fastq.gz kmh12_R1.fastq.gz
mv SRR13809743_1.fastq.gz kmh13_R1.fastq.gz
mv SRR13809742_1.fastq.gz kmh14_R1.fastq.gz

mv SRR13809747_2.fastq.gz kmh1_R2.fastq.gz
mv SRR13809746_2.fastq.gz kmh2_R2.fastq.gz
mv SRR13809741_2.fastq.gz kmh3_R2.fastq.gz
mv SRR13809740_2.fastq.gz kmh4_R2.fastq.gz
mv SRR13809739_2.fastq.gz kmh5_R2.fastq.gz
mv SRR13809738_2.fastq.gz kmh6_R2.fastq.gz
mv SRR13809737_2.fastq.gz kmh7_R2.fastq.gz
mv SRR13809736_2.fastq.gz kmh8_R2.fastq.gz
mv SRR13809735_2.fastq.gz kmh9_R2.fastq.gz
mv SRR13809734_2.fastq.gz kmh10_R2.fastq.gz
mv SRR13809745_2.fastq.gz kmh11_R2.fastq.gz
mv SRR13809744_2.fastq.gz kmh12_R2.fastq.gz
mv SRR13809743_2.fastq.gz kmh13_R2.fastq.gz
mv SRR13809742_2.fastq.gz kmh14_R2.fastq.gz

echo -e "Done."
