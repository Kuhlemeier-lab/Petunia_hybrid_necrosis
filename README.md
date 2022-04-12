# Petunia axillaris X P. exserta hybrid necrosis

Code used in the study of Petunia axillaris X P. exserta hybrid necrosis

## BSR-Seq

To identify the regions of the genome asssociated with the presence and abscence of the necrosis phenotype, we applied a BSR-Seq approach. This is similar to a BSA (bulk segregant analysis) but is performed on RNA-seq data. The approach is very similar to that used by Soyk et al. 2017 [doi.org/10.1038/ng.3733](https://doi.org/10.1038/ng.3733).

We produced two F2 populations, coming from the crosses P. axillaris N X P. exserta and P. exserta X P. axillaris N.
The plants from these populations were then phenotyped for the presence and absence of necrotic symptoms and were grouped accordingly.
Note that in the article we only present the results from the cross P. axillaris N X P. exserta. The raw reads for the other cross are nonetheless available on SRA.

We performed RNA extraction (see article materials and methods for details), and pooled the necrotic individuals together in a sample, and the healthy individual in another sample.

### Sequencing

Was performed by Novogene, with the Eukaryotic RNA-seq (polyA enriched library prep).
- paired end reads
- 150 bp read length
- 250-300 bp insert cDNA library
- Sanger Illumina 1.9 encoding

### Raw reads

Are available on NCBI BioProject [PRJNA708139](https://www.ncbi.nlm.nih.gov/bioproject/PRJNA708139), under SRA accessions:
- SRR13907521: healthy individuals, P. exserta X P. axillaris N
- SRR13907522: necrotic individuals, P. exserta X P. axillaris N
- SRR13907523: healthy individuals, P. axillaris N X P. exserta
- SRR13907524: necrotic individuals, P. axillaris N X P. exserta

Forward reads are numbered 1, reverse are numbered 2.

Read quality asessment performed with fastqc.
Reads trimmed with trimmomatic parameters `LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36`.
Script [bsa01_quality_trim.sh](code/bsa01_quality_trim.sh).

Read numbers before and after cleaning listed in [read_alignment_stats.csv](data/read_alignment_stats.csv).

### Read alignment

Genome index with STAR: [bsa02_STAR_index.sh](code/bsa02_STAR_index.sh).

Alignment with [bsa03_STARmapping.sh](code/bsa03_STARmapping.sh).

```
--sjdbOverhang 149
--outFilterType BySJout
--outFilterMultimapNmax 20
--twopassMode Basic
```

Aligned reads number in [read_alignment_stats.csv](data/read_alignment_stats.csv).

We then add the read group info in each bam file with [bsa04_readGroup_index.sh](code/bsa04_readGroup_index.sh).

And mark duplicated reads with [bsa05_markdup.sh](code/bsa05_markdup.sh).

### Variant calling

We called variants using GATK v4.0.4.0, using the haplotypeCaller and filtering out variants based on quality values because Petunia does not have a high quality variants database to use BQSR.

Genome dictionary created in a previous step with [bsa04_readGroup_index.sh](code/bsa04_readGroup_index.sh), same script to add read group info in the bam files.

We called variants with [bsa06_SNPcalling.sh](code/bsa06_SNPcalling.sh), with parameters

```
--dont-use-soft-clipped-bases
--pcr-indel-model NONE
--stand-call-conf 30
--native-pair-hmm-threads 4
```

We then plot the distribution of the quality values of the variants to know if the GATK suggested filters are fine. This is done with [bsa07_SNPselect_qualityplot.sh](code/bsa07_SNPselect_qualityplot.sh) which uses [plot_vcfq_distribution.R](code/plot_vcfq_distribution.R) to make the [plots](data/snp_quality.pdf).

We then apply the filters for quality with [bsa08_SNPfilter.sh](code/bsa08_SNPfilter.sh). Filters:

```
-window 10 
-cluster 3 
--filter-expression "QD<2.0" 
--filter-name "QD" 
--filter-expression "FS>30.0" 
--filter-name "FS" 
```

Then we select variants passing filters and only the biallelic SNPs. We then keep only variants where the read depth is at least 100 and then we thin the dataset to 100bp.

Number of variants in each dataset is listed in [variant_numbers.csv](data/variant_numbers.csv).

The final set of variants is available in [BSA_SNP_biallelic_gatkselected_minDP100_thin100.recode.vcf](data/BSA_SNP_biallelic_gatkselected_minDP100_thin100.recode.vcf).

### Bulk segregant analysis

Is performed in R, with script [bsa09_analysis.R](code/bsa09_analysis.R).




## IL shallow sequencing

Low coverage sequencing performed to define the boundaries of the introgression region in the IL5 line (IL2 in Hermann et al 2013, https://doi.org/10.1016/j.cub.2013.03.069).


### Sequencing

Library preparation and sequencing were performed by the Next Generation Sequencing platform of the University of Bern, with settings:
- whole genome sequencing library
- Paired end
- **length?**
- Illumina HiSeq 3000


### Raw reads

Have been uploaded to NCBI SRA under BioProject [PRJNA705072](https://www.ncbi.nlm.nih.gov/bioproject/PRJNA705072).

- SRR13809747: sample KMH1
- SRR13809746: sample KMH2
- SRR13809741: sample KMH3
- SRR13809740: sample KMH4
- SRR13809739: sample KMH5
- SRR13809738: sample KMH6
- SRR13809737: sample KMH7
- SRR13809736: sample KMH8
- SRR13809735: sample KMH9
- SRR13809734: sample KMH10
- SRR13809745: sample KMH11
- SRR13809744: sample KMH12
- SRR13809743: sample KMH13
- SRR13809742: sample KMH14

### Alignment

### Variant calling

## RNAseq

To identify genes differentially expressed in the necrotic and healthy plants from the introgression lines we performed an RNAseq experiment.
The samples used are plants constituting the progeny of a single selfed plant heterozygous for the introgression IL2 (see article materials and methods for details).
In the progeny of this heterozygous plants we have some plants homozygous exserta, homozygous axillaris and heterozygous in the introgression.
Three plants per genotype were selected, and one leaf per plant was collected. The tissue was collected from leaves of the homozygous axillaris IL when they just started displaying necrotic symptoms, and equivalent tissue was collected from the other genotypes.

### Sequencing

Was performed at the [Lausanne Genomics Technologies Facility](https://wp.unil.ch/gtf/).
- Single end reads
- 125 bp long
- Illumina HiSeq
- TruSeq stranded RNA library prep

### Raw reads

Have been uploaded to NCBI SRA under BioProject [PRJNA705649](https://www.ncbi.nlm.nih.gov/bioproject/PRJNA705649).

- SRR13861738: heterozygous introgression, replicate 3
- SRR13861739: heterozygous introgression, replicate 2
- SRR13861740: heterozygous introgression, replicate 1
- SRR13861741: homozygous exserta introgression, replicate 3
- SRR13861742: homozygous exserta introgression, replicate 2
- SRR13861743: homozygous exserta introgression, replicate 1
- SRR13861744: homozygous axillaris introgression, replicate 3
- SRR13861745: homozygous axillaris introgression, replicate 2
- SRR13861746: homozygous axillaris introgression, replicate 1

### Alignment


### Read count

### Differential expression analysis

### Allele-specific expression analysis

### GO term analysis


## Software versions

fastqc/0.11.7

trimmomatic/0.36

STAR/2.6.0c

picard-tools/2.18.11

samtools/1.8

vcftools/0.1.15

GenomeAnalysisTK/4.0.4.0

R on the computing cluster: R/3.4.2


