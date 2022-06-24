#!/bin/bash

#SBATCH --mail-user=xxx

# Mail on NONE, BEGIN, END, FAIL, REQUEUE, ALL
#SBATCH --mail-type=ALL

#SBATCH --job-name="8annvcf"

# Runtime and memory
#SBATCH --time=02:00:00 # 1'
#SBATCH --mem-per-cpu=4G # 1Mb  #3Gb
#SBATCH --cpus-per-task=4

# Partition
#SBATCH --partition=epyc2

#SBATCH --account=xxx
#SBATCH --chdir=/xxx/necrosis
#SBATCH --output=code/il08b_annotate_vcf_%A.out
#SBATCH --error=code/il08b_annotate_vcf_%A.err
#################################
chdir=/xxx/necrosis
scdir=/xxx/necrosis


echo -e "#### Annotate the final variants file.
##  `date`
##  current dir: `pwd`
####\t####\n
"
module add vital-it
module add UHTS/Analysis/snpEff/4.3t
module list 2>&1

database="Peax402INV"
se_config=${chdir}/code/snpEff-vitalit.config

snpEff -t \
    -v ${database} \
    -c ${se_config} \
    ${scdir}/data/raw/il_variants/clean.vcf > ${scdir}/data/raw/il_variants/clean.annotated.vcf

snpEff -t \
    -v ${database} \
    -c ${se_config} \
    ${scdir}/data/raw/il_variants/clean_phenotype1_necroticIsAlternate.vcf > ${scdir}/data/raw/il_variants/clean_phenotype1_necroticIsAlternate.annotated.vcf

snpEff -t \
    -v ${database} \
    -c ${se_config} \
    ${scdir}/data/raw/il_variants/clean_phenotype2_necroticIsReference.vcf > ${scdir}/data/raw/il_variants/clean_phenotype2_necroticIsReference.annotated.vcf
