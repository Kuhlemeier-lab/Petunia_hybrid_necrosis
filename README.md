# Petunia axillaris X P. exserta hybrid necrosis

Code used in the study of Petunia axillaris X P. exserta hybrid necrosis, associated with **publication**.

Authors: Chaobin Li, Mathieu Hanemian, Marta Binaghi, Cris Kuhlemeier

Author of this page: Marta Binaghi

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

Read quality assessment performed with fastqc.
Reads trimmed with trimmomatic parameters `LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36`.
Script [bsa01_quality_trim.sh](code/bsa01_quality_trim.sh).

Read numbers before and after cleaning listed in [bsa_read_alignment_stats.csv](data/bsa_read_alignment_stats.csv).

### Read alignment

Genome index with STAR: [bsa02_STAR_index.sh](code/bsa02_STAR_index.sh).

Alignment with [bsa03_STARmapping.sh](code/bsa03_STARmapping.sh).

```
--sjdbOverhang 149
--outFilterType BySJout
--outFilterMultimapNmax 20
--twopassMode Basic
```

Aligned reads number in [bsa_read_alignment_stats.csv](data/bsa_read_alignment_stats.csv).

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

We then plot the distribution of the quality values of the variants to know if the GATK suggested filters are fine. This is done with [bsa07_SNPselect_qualityplot.sh](code/bsa07_SNPselect_qualityplot.sh) which uses [plot_vcfq_distribution.R](code/plot_vcfq_distribution.R) to make the [plots](data/bsa_snp_quality.pdf).

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

Number of variants in each dataset is listed in [bsa_variant_numbers.csv](data/bsa_variant_numbers.csv).

The final set of variants is available in [BSA_SNP_biallelic_gatkselected_minDP100_thin100.recode.vcf](data/BSA_SNP_biallelic_gatkselected_minDP100_thin100.recode.vcf).

### Bulk segregant analysis

Is performed in R, with script [bsa09_analysis.R](code/bsa09_analysis.R).

The plots shown in the manuscript are obtained with script [bsa10_plots_manuscript.R](code/bsa10_plots_manuscript.R).


---

## IL shallow sequencing

Low coverage sequencing performed to define the boundaries of the introgression region in the IL5 line (IL2 in Hermann et al 2013, https://doi.org/10.1016/j.cub.2013.03.069).


### Sequencing

Library preparation and sequencing were performed by the Next Generation Sequencing platform of the University of Bern, with settings:
- whole genome sequencing library
- Paired end
- 150 bp long reads
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

Forward reads are numbered 1, reverse are numbered 2.

Reads renamed with [il01_rename_raw_reads.sh](code/il01_rename_raw_reads.sh).

Read quality assessment performed with fastqc. Reads trimmed with trimmomatic parameters `LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:100`.
Script [il02_fastqc_trim.sh](code/il02_fastqc_trim.sh). Summary of raw and trimmed reads obtained with fastq-stats, [il03_stats_index_genome.sh](code/il03_stats_index_genome.sh).

Read numbers before and after cleaning listed in [il_read_alignment_stats.csv](data/il_read_alignment_stats.csv).


### Alignment

Genome index done in script [il03_stats_index_genome.sh](code/il03_stats_index_genome.sh), alignment performed with BWA MEM, with default parameters in script [il04_align.sh](code/il04_align.sh).

The proportion of reads aligned is shown in the table [il_read_alignment_stats.csv](data/il_read_alignment_stats.csv).

### Variant calling

Is performed with GATK 4.0.4.0, following GATK best practices for organisms without a database of high quality variants. Variants are called in ERC mode.
Variant calling in script [il05_call_variants.sh](code/il05_call_variants.sh), then the single samples are combined in a vcf file with [il06_combineGvcfs.sh](il06_combineGvcfs.sh). The same script is also extracting the quality values to plot them and verify that the hard filters of GATK are suitable for the dataset. The plotting is done with an R script [plot_vcfq_distribution.R](code/plot_vcfq_distribution.R). The plots are shown in [il_snp_quality.pdf](data/il_snp_quality.pdf) and [il_indel_quality.pdf](data/il_indel_quality.pdf).

I then apply the hard filters to the variants in script [il07_filter_variants.sh](code/il07_filter_variants.sh). In this script we also perform some additional steps. In particular we filter the variants in order to select only those that are associated with the necrotic and healthy phenotype. To do so, we simply require the variant to be homozygous reference in the necrotic samples and homozygous alternate in the healthy samples. We also perform the opposite selection, just to check.
The resulting dataset of variants which follow the phenotype of the samples (necrotic is homozygous reference) is in [il_clean_phenotype2_necroticIsReference.vcf](data/il_clean_phenotype2_necroticIsReference.vcf). **to add**

We then also annotate the variants using SnpEff. **to do**



---

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

- bwa/0.7.17
- fastqc/0.11.7
- fastq-stats ea-utils/1.1.2
- GenomeAnalysisTK/4.0.4.0
- picard-tools/2.18.11 or 2.9.0, see scripts to know which version was used
- R on the computing cluster: R/3.4.2
- R on the local machine: R/3.3.3
- samtools/1.8
- STAR/2.6.0c
- trimmomatic/0.36
- vcftools/0.1.15


## R libraries used

For the versions, see at the bottom of the R scripts in the [code](code/) folder.

- vcfUtils custom package, available here [vcfUtils.tar.gz](code/vcfUtils.tar.gz)
- tidyr Wickham H, Girlich M (2022). tidyr: Tidy Messy Data. https://tidyr.tidyverse.org, https://github.com/tidyverse/tidyr
- dplyr Wickham H, François R, Henry L, Müller K (2022). dplyr: A Grammar of Data Manipulation. https://dplyr.tidyverse.org, https://github.com/tidyverse/dplyr
- ggplot2 Wickham H (2016). ggplot2: Elegant Graphics for Data Analysis. Springer-Verlag New York. ISBN 978-3-319-24277-4, https://ggplot2.tidyverse.org
- optparse https://cran.r-project.org/package=optparse


