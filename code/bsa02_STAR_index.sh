#!/bin/bash

#SBATCH --mail-user=xxx

# Mail on NONE, BEGIN, END, FAIL, REQUEUE, ALL
#SBATCH --mail-type=ALL

#SBATCH --job-name="2index"

# Runtime and memory
#SBATCH --time=08:00:00  #
#SBATCH --mem-per-cpu=32Gb  # 
#SBATCH --cpus-per-task=4

# Partition
#SBATCH --partition=epyc2

#SBATCH --account=xxx
#SBATCH --chdir=/xxx/necrosis
#SBATCH --output=code/bsa02_STAR_index_%A.out
#SBATCH --error=code/bsa02_STAR_index_%A.err

#################################
chdir=/xxx/necrosis
scdir=/xxx/necrosis

module load vital-it
module load UHTS/Aligner/STAR/2.6.0c

# Build a genome index for STAR aligner

## link the original genome file from a location in the cluster to the current project folder
ln -s /xxx/genomes/Peax402INV.fasta ${scdir}/data/genomes/Peax402INV.fasta
## link the annotation file
ln -s /xxx/genomes/peaxi162AQ_Peax402INV.cds.gff ${scdir}/data/genomes/peaxi162AQ_Peax402INV.cds.gff

echo -e "## `date`\n## Genome indexing\n"

STAR --runMode genomeGenerate \
  --runThreadN 4 \
  --genomeDir ${scdir}/data/genomes \
  --genomeFastaFiles ${scdir}/data/genomes/Peax402INV.fasta \
  --sjdbGTFfile ${scdir}/data/genomes/peaxi162AQ_Peax402INV.cds.gff \
  --sjdbOverhang 149 \
  --sjdbGTFtagExonParentTranscript Parent \
  --genomeChrBinNbits 16 --limitGenomeGenerateRAM 24000000000

echo -e "Indexing done.\n `date`"
