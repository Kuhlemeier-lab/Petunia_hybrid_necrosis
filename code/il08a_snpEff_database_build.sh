#!/bin/bash

#SBATCH --mail-user=xxx
# Mail on NONE, BEGIN, END, FAIL, REQUEUE, ALL
#SBATCH --mail-type=ALL

#SBATCH --job-name="snpeffDB"

# Runtime and memory
#SBATCH --time=02:00:00  # 1'
#SBATCH --mem-per-cpu=4G  # 2Gb

# Partition
#SBATCH --partition=epyc2

#SBATCH --account=xxx
#SBATCH --chdir=/xxx/necrosis
#SBATCH --output=code/il08a_snpEff_database_build_%A.out
#SBATCH --error=code/il08a_snpEff_database_build_%A.err
#################################
chdir=/xxx/necrosis
scdir=/xxx/necrosis

module add vital-it
module add UHTS/Analysis/snpEff/4.3t

genome=${scdir}/data/genomes/Peax402INV.fasta
annotation=${scdir}/data/genomes/peaxi162AQ_Peax402INV.cds.gff

# specify a local config file
se_config=${chdir}/code/snpEff-vitalit.config

if [ ! -d ${scdir}/data/genomes/snpEff/Peax402INV ]; then
	mkdir -p ${scdir}/data/genomes/snpEff/Peax402INV
fi

cd ${scdir}/data/genomes/snpEff/Peax402INV

cp ${genome} sequences.fa
cp ${annotation} genes.gff

head ${se_config}

snpEff build \
    -gff3 -v Peax402INV \
    -c ${se_config}
