#!/bin/bash

#SBATCH --mail-user=xxx
#SBATCH --mail-type=ALL
#SBATCH --job-name="5counts"

# Runtime and memory
#SBATCH --time=02:00:00  # 10'
#SBATCH --mem-per-cpu=4G  # 500Mb
#SBATCH --cpus-per-task=16
# Partition
#SBATCH --partition=epyc2

#SBATCH --account=xxx
#SBATCH --chdir=/xxx/necrosis
#SBATCH --output=code/rs05_featureCounts_%A.out
#SBATCH --error=code/rs05_featureCounts_%A.err

#################################
chdir=/xxx/necrosis
scdir=/xxx/necrosis

if [ ! -d ${scdir}/data/raw/rs_counts ]; then
	mkdir -p ${scdir}/data/raw/rs_counts
fi

module load vital-it
module add UHTS/Analysis/subread/1.6.0;

# move to folder containing the aligned reads
cd ${scdir}/data/raw/rs_aligned_reads

echo -e "## `date`\n## Counting starts "

featureCounts -T 16 \
    -t gene \
    -s 2 \
    -g ID \
    -a ${scdir}/data/genomes/peaxi162AQ_Peax402INV.cds.gff \
    -o ${scdir}/data/raw/rs_counts/rs_counts.txt kmh1_trim_STAR_Aligned.sortedByCoord.out.bam kmh2_trim_STAR_Aligned.sortedByCoord.out.bam kmh3_trim_STAR_Aligned.sortedByCoord.out.bam kmh4_trim_STAR_Aligned.sortedByCoord.out.bam kmh5_trim_STAR_Aligned.sortedByCoord.out.bam kmh6_trim_STAR_Aligned.sortedByCoord.out.bam kmh7_trim_STAR_Aligned.sortedByCoord.out.bam kmh8_trim_STAR_Aligned.sortedByCoord.out.bam kmh9_trim_STAR_Aligned.sortedByCoord.out.bam

echo -e "Counting done.\n `date`"
