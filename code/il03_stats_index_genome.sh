#!/bin/bash

#SBATCH --mail-user=xxx

# Mail on NONE, BEGIN, END, FAIL, REQUEUE, ALL
#SBATCH --mail-type=ALL

#SBATCH --job-name="3stats"

# Runtime and memory
#SBATCH --time=02:00:00  #1.5h
#SBATCH --mem-per-cpu=6G  # 2Gb
#SBATCH --cpus-per-task=4

# Partition
#SBATCH --partition=epyc2
#SBATCH --account=ips_ck
#SBATCH --chdir=/xxx/necrosis
#SBATCH --output=code/il03_stats_index_genome_%A.out
#SBATCH --error=code/il03_stats_index_genome_%A.err
#################################
chdir=/xxx/necrosis
scdir=/xxx/necrosis


echo -e "#### Stats on trimmed reads and genome index
## `date`
## current dir: `pwd`
####\t####\n"

module load vital-it
module load SequenceAnalysis/ea-utils/1.1.2		#fastq-stats
module add UHTS/Aligner/bwa/0.7.17 
module add UHTS/Analysis/samtools/1.8
module add UHTS/Analysis/picard-tools/2.18.11
module list 2>&1

## Raw reads quality statistics
#echo -e "## Calculation of quality stats of raw reads\n"
#cd ${scdir}/data/raw/il_rawreads
#echo -e "Stats on raw reads\n`date`\n" > stats_raw_reads.txt
#for file in *.fastq.gz  # for each raw file in the folder, calculates stats and store them in stats_raw_reads.txt
#do
#	echo -e "\n################################# \n $file" >> stats_raw_reads.txt
#	fastq-stats $file >> stats_raw_reads.txt
#done
#echo -e "Quality stats done and stored in data/raw/il_rawreads/stats_raw_reads.txt.\n"
## Clean reads quality statistics
#echo -e "## Calculation of quality stats of trimmed reads\n"
#cd ${scdir}/data/raw/il_rawreads
#echo -e "Stats on trimmed reads\n`date`\n" > stats_trimmed_reads.txt
#for file in *.fastq.gz.trimmed  # for each trimmed file in the folder, calculates stats and store them in stats_trimmed_reads.txt
#do
#	echo -e "\n################################# \n $file" >> stats_trimmed_reads.txt
#	fastq-stats $file >> stats_trimmed_reads.txt
#done
#echo -e "Quality stats done and stored in data/raw/raw_reads/stats_trimmed_reads.txt.\n"
#
## Index genome
genome=${scdir}/data/genomes/Peax402INV.fasta
#
#bwa index ${genome}
#samtools faidx ${genome}
picard-tools CreateSequenceDictionary R=${genome} O=${genome}.dict
echo -e "Indexing done.\n"
