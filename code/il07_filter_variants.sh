#!/bin/bash

#SBATCH --mail-user=xxx

# Mail on NONE, BEGIN, END, FAIL, REQUEUE, ALL
#SBATCH --mail-type=ALL

#SBATCH --job-name="7filt"

# Runtime and memory
#SBATCH --time=02:10:00  # 4h
#SBATCH --mem-per-cpu=8G  # 6Gb
#SBATCH --cpus-per-task=4

# Partition
#SBATCH --partition=epyc2
##SBATCH --qos=job_epyc2_debug
#SBATCH --account=ips_ck
#SBATCH --chdir=/xxx/necrosis
#SBATCH --output=code/il07_filter_variants_%A.out
#SBATCH --error=code/il07_filter_variants_%A.err
#################################
chdir=/xxx/necrosis
scdir=/xxx/necrosis

module load vital-it
module add UHTS/Analysis/GenomeAnalysisTK/4.0.4.0
module add UHTS/Analysis/picard-tools/2.18.11

module list 2>&1

echo -e "#### Filter variants for quality and select
 those that match the phenotype observed.

# Apply hard filters to SNPs and INDELs separately,
#  extract variants following the phenotype,
#  summarise variant counts in file summary.txt.
#  `date`
#  current dir: `pwd`
###\t####\n
"

genome=${scdir}/data/genomes/Peax402INV.fasta

# Filter variants
echo -e "## Hard filter SNP variants 
"

# note that the online docs for gatk 4.0.4.0 report different
# argument names (filterExpression) but in fact they are wrong.
GenomeAnalysisTK VariantFiltration \
   -R ${genome} \
   -V ${scdir}/data/raw/il_variants/raw_snp.vcf \
   --filter-expression "QD<2.0" \
   --filter-name "QD" \
   --filter-expression "FS>60.0" \
   --filter-name "FS" \
   --filter-expression "SOR>3.0" \
   --filter-name "SOR" \
   --filter-expression "MQ<40.0" \
   --filter-name "MQ" \
   --filter-expression "MQRankSum<-12.5" \
   --filter-name "MQRankSum" \
   --filter-expression "ReadPosRankSum<-8.0" \
   --filter-name "ReadPosRankSum" \
   -O ${scdir}/data/raw/il_variants/filtered_snp.vcf

echo -e "Hard filter done.
Select variants...
"

GenomeAnalysisTK SelectVariants \
   -R ${genome} \
   -V ${scdir}/data/raw/il_variants/filtered_snp.vcf \
   --exclude-filtered true \
   -O ${scdir}/data/raw/il_variants/selected_snp.vcf

echo -e "Select done.
## Hard filter INDEL variants
"

GenomeAnalysisTK VariantFiltration \
   -R ${genome} \
   -V ${scdir}/data/raw/il_variants/raw_indel.vcf \
   --filter-expression "QD<2.0" \
   --filter-name "QD" \
   --filter-expression "FS>200.0" \
   --filter-name "FS" \
   --filter-expression "SOR>10.0" \
   --filter-name "SOR" \
   --filter-expression "ReadPosRankSum<-20.0" \
   --filter-name "ReadPosRankSum" \
   -O ${scdir}/data/raw/il_variants/filtered_indel.vcf

echo -e "Hard filter done.
Select variants...
"

GenomeAnalysisTK SelectVariants \
   -R ${genome} \
   -V ${scdir}/data/raw/il_variants/filtered_indel.vcf \
   --exclude-filtered true \
   -O ${scdir}/data/raw/il_variants/selected_indel.vcf

# merge SNP and INDEL variants in a single file
echo -e "Selection done.
Merging SNPs and INDELs...
"

picard-tools MergeVcfs \
    I=${scdir}/data/raw/il_variants/selected_snp.vcf \
    I=${scdir}/data/raw/il_variants/selected_indel.vcf \
    O=${scdir}/data/raw/il_variants/clean.vcf

echo -e "Hard filter, select, merge SNP and INDEL done
"

# Extract variants associated to phenotype

echo -e "
# Selection of variants following the phenotype.
"

# the sample name (in the quotes after getGenotype) corresponds to names given in script il04.
# Here we select as if first sample, kmh1 which is necrotic is homozygous for the alternate allele,
GenomeAnalysisTK SelectVariants \
   -R ${genome} \
   --restrict-alleles-to BIALLELIC \
   -V ${scdir}/data/raw/il_variants/clean.vcf \
   --selectExpressions 'vc.getGenotype("kmh1").isHomVar()' \
   --selectExpressions 'vc.getGenotype("kmh2").isHomVar()' \
   --selectExpressions 'vc.getGenotype("kmh3").isHomRef()' \
   --selectExpressions 'vc.getGenotype("kmh4").isHomRef()' \
   --selectExpressions 'vc.getGenotype("kmh5").isHomVar()' \
   --selectExpressions 'vc.getGenotype("kmh6").isHomRef()' \
   --selectExpressions 'vc.getGenotype("kmh7").isHomRef()' \
   --selectExpressions 'vc.getGenotype("kmh8").isHomVar()' \
   --selectExpressions 'vc.getGenotype("kmh9").isHomRef()' \
   --selectExpressions 'vc.getGenotype("kmh10").isHomRef()' \
   --selectExpressions 'vc.getGenotype("kmh11").isHomRef()' \
   --selectExpressions 'vc.getGenotype("kmh12").isHomRef()' \
   --selectExpressions 'vc.getGenotype("kmh13").isHomRef()' \
   --selectExpressions 'vc.getGenotype("kmh14").isHomRef()' \
   -O ${scdir}/data/raw/il_variants/clean_phenotype1_necroticIsAlternate.vcf

# Now we select as if kmh1 (necrotic) is homozygous reference.
GenomeAnalysisTK SelectVariants \
    -R ${genome} \
    --restrict-alleles-to BIALLELIC \
    -V ${scdir}/data/raw/il_variants/clean.vcf \
    --selectExpressions 'vc.getGenotype("kmh1").isHomRef()' \
    --selectExpressions 'vc.getGenotype("kmh2").isHomRef()' \
    --selectExpressions 'vc.getGenotype("kmh3").isHomVar()' \
    --selectExpressions 'vc.getGenotype("kmh4").isHomVar()' \
    --selectExpressions 'vc.getGenotype("kmh5").isHomRef()' \
    --selectExpressions 'vc.getGenotype("kmh6").isHomVar()' \
    --selectExpressions 'vc.getGenotype("kmh7").isHomVar()' \
    --selectExpressions 'vc.getGenotype("kmh8").isHomRef()' \
    --selectExpressions 'vc.getGenotype("kmh9").isHomVar()' \
    --selectExpressions 'vc.getGenotype("kmh10").isHomVar()' \
    --selectExpressions 'vc.getGenotype("kmh11").isHomVar()' \
    --selectExpressions 'vc.getGenotype("kmh12").isHomVar()' \
    --selectExpressions 'vc.getGenotype("kmh13").isHomVar()' \
    --selectExpressions 'vc.getGenotype("kmh14").isHomVar()' \
    -O ${scdir}/data/raw/il_variants/clean_phenotype2_necroticIsReference.vcf

echo -e "Selection done."

# Count variants in different steps and save them in summary.txt
echo -e "
# Summary of variant counts..."

echo -e "####    Summary of variants
## `date`
" > ${scdir}/data/raw/il_variants/summary.txt

# make and empty temporary file to store the summary values before formatting
echo "" > ${scdir}/data/raw/il_variants/summary.tmp

echo -e "FILE\tCOUNT" >> ${scdir}/data/raw/il_variants/summary.tmp
count=$(grep -vc '^#' ${scdir}/data/raw/il_variants/raw.vcf)
echo -e "raw.vcf,${count}" >> ${scdir}/data/raw/il_variants/summary.tmp
count=$(grep -vc '^#' ${scdir}/data/raw/il_variants/raw_snp.vcf)
echo -e "raw_snp.vcf,${count}" >> ${scdir}/data/raw/il_variants/summary.tmp
count=$(grep -vc '^#' ${scdir}/data/raw/il_variants/raw_indel.vcf)
echo -e "raw_indel.vcf,${count}" >> ${scdir}/data/raw/il_variants/summary.tmp
count=$(grep -vc '^#' ${scdir}/data/raw/il_variants/selected_snp.vcf)
echo -e "selected_snp.vcf,${count}" >> ${scdir}/data/raw/il_variants/summary.tmp
count=$(grep -vc '^#' ${scdir}/data/raw/il_variants/selected_indel.vcf)
echo -e "selected_indel.vcf,${count}" >> ${scdir}/data/raw/il_variants/summary.tmp
count=$(grep -vc '^#' ${scdir}/data/raw/il_variants/clean.vcf)
echo -e "clean.vcf,${count}" >> ${scdir}/data/raw/il_variants/summary.tmp
count=$(grep -Pc '1/1.*1/1.*1/1.*1/1.*1/1.*1/1.*1/1.*1/1.*1/1.*1/1.*1/1.*1/1.*1/1.*1/1' ${scdir}/data/raw/il_variants/clean.vcf)
echo -e "all samples homozygous alternate (1/1),${count}" >> ${scdir}/data/raw/il_variants/summary.tmp
count=$(grep -cvP "^#" ${scdir}/data/raw/il_variants/clean_phenotype1_necroticIsAlternate.vcf)
echo -e "clean_phenotype1_necroticIsAlternate.vcf,${count}" >> ${scdir}/data/raw/il_variants/summary.tmp
# for the next step, we checked the order of the samples in the vcf file,
# which is kmh: 1,10,11,12,13,14,2,3,4,5,6,7,8,9
count=$(grep -cP '1/1.*0/0.*0/0.*0/0.*0/0.*0/0.*1/1.*0/0.*0/0.*1/1.*0/0.*0/0.*1/1.*0/0' ${scdir}/data/raw/il_variants/clean.vcf)
echo -e "grep count on phenotype1,${count}" >> ${scdir}/data/raw/il_variants/summary.tmp
count=$(grep -cPv "^#" ${scdir}/data/raw/il_variants/clean_phenotype2_necroticIsReference.vcf)
echo -e "clean_phenotype2_necroticIsReference.vcf,${count}" >> ${scdir}/data/raw/il_variants/summary.tmp
count=$(grep -cP '0/0.*1/1.*1/1.*1/1.*1/1.*1/1.*0/0.*1/1.*1/1.*0/0.*1/1.*1/1.*0/0.*1/1' ${scdir}/data/raw/il_variants/clean.vcf)
echo -e "grep count on phenotype2,${count}" >> ${scdir}/data/raw/il_variants/summary.tmp
cat ${scdir}/data/raw/il_variants/summary.tmp | column -s "," -t >> ${scdir}/data/raw/il_variants/summary.txt  # format first part of the summary table
echo -e "\n### Phenotype 1 Necrotis is alternate ###" >> ${scdir}/data/raw/il_variants/summary.txt
echo -e "VARIANTS\tSCF" >> ${scdir}/data/raw/il_variants/summary.txt
grep -oP '^(Chr|Scaffold_)\d+' ${scdir}/data/raw/il_variants/clean_phenotype1_necroticIsAlternate.vcf | sort | uniq -c >> ${scdir}/data/raw/il_variants/summary.txt
echo -e "\n### Phenotype 2 Necrotic is reference ###" >> ${scdir}/data/raw/il_variants/summary.txt
echo -e "VARIANTS\tSCF" >> ${scdir}/data/raw/il_variants/summary.txt
grep -oP '^(Chr|Scaffold_)\d+' ${scdir}/data/raw/il_variants/clean_phenotype2_necroticIsReference.vcf | sort | uniq -c >> ${scdir}/data/raw/il_variants/summary.txt

# remove summary.tmp
rm ${scdir}/data/raw/il_variants/summary.tmp

echo -e "Summary done. Results are stored in ${scdir}/data/raw/il_variants/summary.txt.

End at `date`"
