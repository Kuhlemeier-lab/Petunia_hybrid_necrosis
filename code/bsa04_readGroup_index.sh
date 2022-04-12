#!/bin/bash

#SBATCH --mail-user=xxx
#SBATCH --mail-type=ALL
#SBATCH --job-name="groupIndex"

# Runtime and memory
#SBATCH --time=24:00:00  # 20h 
#SBATCH --mem-per-cpu=12G  # 170 Mb

# Partition
#SBATCH --partition=epyc2

#SBATCH --account=xxx
#SBATCH --chdir=/xxx/necrosis
#SBATCH --output=code/bsa04_readGroup_index_%A.out
#SBATCH --error=code/bsa04_readGroup_index_%A.err

#################################
chdir=/xxx/necrosis
scdir=/xxx/necrosis

module add vital-it
module add UHTS/Analysis/picard-tools/2.18.11
module add UHTS/Analysis/samtools/1.8

echo -e "## `date`\n## Preparation of index and dictionary..."

picard-tools CreateSequenceDictionary R=${scdir}/data/genomes/Peax402INV.fasta O=${scdir}/data/genomes/Peax402INV.dict

samtools faidx ${scdir}/data/genomes/Peax402INV.fasta

echo -e "Indexing and dictionary done.\n"


echo -e "## `date`\n## Add read groups IDs "

cd ${scdir}/data/raw/bsa_aligned_reads

picard-tools AddOrReplaceReadGroups \
    RGLB=1 \
    RGPL=Illumina \
    RGPU=1 \
    RGSM=23 \
    I=SRR13907523_Aligned.sortedByCoord.out.bam \
    O=SRR13907523_Aligned.sortedByCoord_readgroup.out.bam

picard-tools AddOrReplaceReadGroups \
    RGLB=1 \
    RGPL=Illumina \
    RGPU=2 \
    RGSM=24 \
    I=SRR13907524_Aligned.sortedByCoord.out.bam \
    O=SRR13907524_Aligned.sortedByCoord_readgroup.out.bam

echo -e "Read groups added.\n `date`"

