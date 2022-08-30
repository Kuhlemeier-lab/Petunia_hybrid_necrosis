#!/bin/bash

#SBATCH --mail-user=xxx

# Mail on NONE, BEGIN, END, FAIL, REQUEUE, ALL
#SBATCH --mail-type=ALL

#SBATCH --job-name="5call"

# Runtime and memory
#SBATCH --time=72:00:00  # 36h
#SBATCH --mem-per-cpu=8G  # 6Gb
#SBATCH --cpus-per-task=6

# Array
#SBATCH --array=1-14
# Partition
#SBATCH --partition=epyc2
#SBATCH --account=xxx
#SBATCH --chdir=/xxx/necrosis
#SBATCH --output=code/il05_call_variants_%A_%a.out
#SBATCH --error=code/il05_call_variants_%A_%a.err
#################################
chdir=/xxx/necrosis
scdir=/xxx/necrosis

echo -e "#### Call variants 
## `date`
## current dir: `pwd`
####\t####\n"

module load vital-it
module add UHTS/Analysis/GenomeAnalysisTK/4.0.4.0 
module list 2>&1

genome=${scdir}/data/genomes/Peax403.fasta

if [ ! -d ${scdir}/data/raw/il_variants ]; then
  mkdir ${scdir}/data/raw/il_variants
fi

echo -e "\n## Call variants\n"
GenomeAnalysisTK HaplotypeCaller \
  -R ${genome} \
  -I ${scdir}/data/raw/il_alignedreads/kmh${SLURM_ARRAY_TASK_ID}_ed_markdup.bam \
  -ERC GVCF \
  --do-not-run-physical-phasing \
  -O ${scdir}/data/raw/il_variants/kmh${SLURM_ARRAY_TASK_ID}.g.vcf \
  --native-pair-hmm-threads 6
echo -e "Done call variants.\n"
