#!/bin/bash

#SBATCH --mail-user=xxx

# Mail on NONE, BEGIN, END, FAIL, REQUEUE, ALL
#SBATCH --mail-type=ALL

#SBATCH --job-name="4align"

# Runtime and memory
#SBATCH --time=12:00:00  # 6.5h 
#SBATCH --mem-per-cpu=24G  # 18Gb 
#SBATCH --cpus-per-task=4
# Array
#SBATCH --array=1-14
# Partition
#SBATCH --partition=epyc2
#SBATCH --account=xxx
#SBATCH --chdir=/xxx/necrosis
#SBATCH --output=code/il04_align_%A_%a.out
#SBATCH --error=code/il04_align_%A_%a.err
#################################
chdir=/xxx/necrosis
scdir=/xxx/necrosis

echo -e "#### Align reads, mark duplicates
## `date`
## current dir: `pwd`
####\t####\n"

module load vital-it
module add UHTS/Aligner/bwa/0.7.17
module add UHTS/Analysis/samtools/1.8
#module add Development/java/1.8.0_121
module add UHTS/Analysis/picard-tools/2.9.0 

module list 2>&1

if [ ! -d ${scdir}/data/raw/il_alignedreads ]; then
  mkdir ${scdir}/data/raw/il_alignedreads
fi

genome=${scdir}/data/genomes/Peax404.fasta

# align reads with BWA
bwa mem -t 4 \
  -M ${genome} \
  ${scdir}/data/raw/il_rawreads/kmh${SLURM_ARRAY_TASK_ID}_R1.fastq.gz.trimmed \
  ${scdir}/data/raw/il_rawreads/kmh${SLURM_ARRAY_TASK_ID}_R2.fastq.gz.trimmed > ${scdir}/data/raw/il_alignedreads/kmh${SLURM_ARRAY_TASK_ID}.sam

# compress the aligned reads into an unsorted bam file
samtools view -b -@ 4 -t ${genome}.fai ${scdir}/data/raw/il_alignedreads/kmh${SLURM_ARRAY_TASK_ID}.sam > ${scdir}/data/raw/il_alignedreads/kmh${SLURM_ARRAY_TASK_ID}_unsorted.bam

# sort the reads in the bam file
samtools sort -@ 4 -o ${scdir}/data/raw/il_alignedreads/kmh${SLURM_ARRAY_TASK_ID}.bam ${scdir}/data/raw/il_alignedreads/kmh${SLURM_ARRAY_TASK_ID}_unsorted.bam

# remove the .sam and unsorted.bam files because no longer needed.
rm ${scdir}/data/raw/il_alignedreads/kmh${SLURM_ARRAY_TASK_ID}_unsorted.bam
rm ${scdir}/data/raw/il_alignedreads/kmh${SLURM_ARRAY_TASK_ID}.sam
echo -e "\n## Done aligning.\n"

# Add sample name in bam files
echo -e "\n## Edit name bam\n"
picard-tools AddOrReplaceReadGroups \
  I=${scdir}/data/raw/il_alignedreads/kmh${SLURM_ARRAY_TASK_ID}.bam \
  O=${scdir}/data/raw/il_alignedreads/kmh${SLURM_ARRAY_TASK_ID}_ed.bam \
  RGID=${SLURM_ARRAY_TASK_ID} \
  RGLB=lib${SLURM_ARRAY_TASK_ID} \
  RGPL=illumina \
  RGPU=unit1 \
  RGSM=kmh${SLURM_ARRAY_TASK_ID}
echo -e "Edit sample names in bam files done.\n"

echo -e "\n## Mark duplicates \n"
picard-tools MarkDuplicates \
  I=${scdir}/data/raw/il_alignedreads/kmh${SLURM_ARRAY_TASK_ID}_ed.bam \
  O=${scdir}/data/raw/il_alignedreads/kmh${SLURM_ARRAY_TASK_ID}_ed_markdup.bam \
  M=${scdir}/data/raw/il_alignedreads/kmh${SLURM_ARRAY_TASK_ID}_ed_markdup_metrics.txt \
  MAX_FILE_HANDLES= 1000 \
  TAGGING_POLICY=All \
  CREATE_INDEX=true

echo -e "\n## Collect WGS metrics.\n"
picard-tools CollectWgsMetrics \
  I=${scdir}/data/raw/il_alignedreads/kmh${SLURM_ARRAY_TASK_ID}_ed_markdup.bam \
  O=${scdir}/data/raw/il_alignedreads/kmh${SLURM_ARRAY_TASK_ID}_ed_markdup_wgsmetrics.txt \
  R=${genome}

echo -e "Count mapped reads\n"
samtools view -F 0x4 ${scdir}/data/raw/il_alignedreads/kmh${SLURM_ARRAY_TASK_ID}_ed_markdup.bam | cut -f 1 | sort | uniq | wc -l

# remove old bam file
rm ${scdir}/data/raw/il_alignedreads/kmh${SLURM_ARRAY_TASK_ID}_ed.bam ${scdir}/data/raw/il_alignedreads/kmh${SLURM_ARRAY_TASK_ID}_ed.bam.bai
