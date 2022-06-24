# Petunia axillaris X P. exserta hybrid necrosis

Code used in the study of Petunia axillaris X P. exserta hybrid necrosis, associated with **publication**.

Authors: Chaobin Li, Mathieu Hanemian, Marta Binaghi, Cris Kuhlemeier

Author of this page: Marta Binaghi

Genome sequence and annotation used in all the analyses: **to add**

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
The resulting datasets of variants which follow the phenotype of the samples (necrotic is homozygous reference) is in [il_clean_phenotype2_necroticIsReference.vcf](data/il_clean_phenotype2_necroticIsReference.vcf). The opposite selection is available in [il_clean_phenotype1_necroticIsAlternate.vcf](data/il_clean_phenotype1_necroticIsAlternate.vcf).
**to add**

### Variant annotation MAYBE WILL NOT BE INCLUDED

We annotate the variants using [SnpEff](https://pcingola.github.io/SnpEff/).

The genome of Petunia axillaris in the version that we used is not available as pre-built in snpEff so we used script [il08a_snpEff_database_build.sh](code/il08a_snpEff_database_build.sh) to build a custom database starting from the fasta and gff file. Notice that if you want to do this, you need to modify the config file accordingly. See instructions in the snpEff documentation.

We then annotate the vcf files with script [il08b_annotate_vcf.sh](code/il08b_annotate_vcf.sh).

---

## RNAseq

To identify genes differentially expressed in the necrotic and healthy plants from the introgression lines we performed an RNAseq experiment.
The samples used are plants constituting the progeny of a single selfed plant heterozygous for the introgression IL2 (see article materials and methods for details).
In the progeny of this heterozygous plants we have some plants homozygous exserta, homozygous axillaris and heterozygous in the introgression.
Three plants per genotype were selected, and one leaf per plant was collected. The tissue was collected from leaves of the homozygous axillaris IL when they just started displaying necrotic symptoms, and equivalent tissue was collected from the other genotypes.

RNA was extracted with **to add**

The RNAseq data were used to perform a differential expression (DE) analysis and the results were then used to perform a GO term analysis.

### Sequencing

Was performed at the [Lausanne Genomics Technologies Facility](https://wp.unil.ch/gtf/).
- Single end reads
- 125 bp long
- Illumina HiSeq
- TruSeq stranded RNA library prep

### Raw reads

Have been uploaded to NCBI SRA under BioProject [PRJNA705649](https://www.ncbi.nlm.nih.gov/bioproject/PRJNA705649).

- SRR13861738: heterozygous introgression, replicate 3 (kmh9)
- SRR13861739: heterozygous introgression, replicate 2 (kmh8)
- SRR13861740: heterozygous introgression, replicate 1 (kmh7)
- SRR13861741: homozygous exserta introgression, replicate 3 (kmh6)
- SRR13861742: homozygous exserta introgression, replicate 2 (kmh5)
- SRR13861743: homozygous exserta introgression, replicate 1 (kmh4)
- SRR13861744: homozygous axillaris introgression, replicate 3 (kmh3)
- SRR13861745: homozygous axillaris introgression, replicate 2 (kmh2)
- SRR13861746: homozygous axillaris introgression, replicate 1 (kmh1)

Raw reads are renamed with script [rs01_rename_raw_reads.sh](code/rs01_rename_raw_reads.sh).

Quality control and trimming in script [rs02_fastqc_trim.sh](code/rs02_fastqc_trim.sh). Parameters for trimming:

```
trimmomatic SE \
  -threads 4 \
  -phred33 \
  kmh${SLURM_ARRAY_TASK_ID}.fastq.gz \
  kmh${SLURM_ARRAY_TASK_ID}_trim.fastq.gz \
  LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36 \
  ILLUMINACLIP:TruSeq3-SE.fa:2:30:3
```


The read numbers before and after trimming are reported in [rs_read_alignment_stats.csv](data/rs_read_alignment_stats.csv).

### Alignment

Note that I build another STAR index for the genome because in the genomeGenerate step they require the read length and here it is 125 bp. So I use script [rs03_STAR_genome_index.sh](code/rs03_STAR_genome_index.sh) to redo the genome index. Parameters:

```
STAR --runMode genomeGenerate \
  --runThreadN 4 \
  --genomeDir ${scdir}/data/genomes/Pax_125bp \
  --genomeFastaFiles ${scdir}/data/genomes/Peax402INV.fasta \
  --sjdbGTFfile ${scdir}/data/genomes/peaxi162AQ_Peax402INV.cds.gff \
  --sjdbOverhang 124 \
  --sjdbGTFtagExonParentTranscript Parent \
  --genomeChrBinNbits 16 --limitGenomeGenerateRAM 24000000000
```

The alignment is then performed with [rs04_STARmapping.sh](code/rs04_STARmapping.sh). Aligned read numbers reported in [rs_read_alignment_stats.csv](data/rs_read_alignment_stats.csv). Alignment performed with parameters:

```
STAR --genomeDir data/genomes/Pax_125bp \
    --runThreadN 16 \
    --readFilesIn data/raw/rs_rawreads/kmh${SLURM_ARRAY_TASK_ID}_trim.fastq.gz \
    --outFilterType BySJout \
    --outFilterMultimapNmax 20 \
    --outSAMtype BAM SortedByCoordinate \
    --twopassMode Basic \
    --outFileNamePrefix data/raw/rs_aligned_reads/kmh${SLURM_ARRAY_TASK_ID}_trim_STAR_ \
    --readFilesCommand zcat \
    --genomeLoad NoSharedMemory
```

### Read count

Is performed with subread featureCounts function, in script [rs05_featureCounts.sh](code/rs05_featureCounts.sh). Parameters are:

```
featureCounts -T 16 \
    -t gene \
    -s 2 \
    -g ID \
    -a ${scdir}/data/genomes/peaxi162AQ_Peax402INV.cds.gff \
    -o ${scdir}/data/raw/rs_counts/rs_counts.txt kmh1_trim_STAR_Aligned.sortedByCoord.out.bam kmh2_trim_STAR_Aligned.sortedByCoord.out.bam kmh3_trim_STAR_Aligned.sortedByCoord.out.bam kmh4_trim_STAR_Aligned.sortedByCoord.out.bam kmh5_trim_STAR_Aligned.sortedByCoord.out.bam kmh6_trim_STAR_Aligned.sortedByCoord.out.bam kmh7_trim_STAR_Aligned.sortedByCoord.out.bam kmh8_trim_STAR_Aligned.sortedByCoord.out.bam kmh9_trim_STAR_Aligned.sortedByCoord.out.bam
```

### Differential expression analysis

The DE analysis is done with DESeq2 in the R markdown script [rs06_DE_analysis.Rmd](code/rs06_DE_analysis.Rmd). The output is in [rs06_DE_analysis.html](code/rs06_DE_analysis.html).

The results including raw gene counts, normalised gene counts and DE of the comparison axillaris VS exserta are available in table [rs_DE_results_axVSex.csv](data/rs_DE_results_axVSex.csv).



## Software versions

- bwa/0.7.17
- fastqc/0.11.7
- fastq-stats ea-utils/1.1.2
- GenomeAnalysisTK/4.0.4.0
- picard-tools/2.18.11 or 2.9.0, see scripts to know which version was used
- R on the computing cluster: R/3.4.2
- R on the local machine: R/3.3.3
- samtools/1.8
- SnpEff/4.3t
- STAR/2.6.0c
- subread/1.6.0
- trimmomatic/0.36
- vcftools/0.1.15


## R libraries used

For the versions, see at the bottom of the R scripts in the [code](code/) folder.

- DESeq2 Love MI, Huber W, Anders S (2014). “Moderated estimation of fold change and dispersion for RNA-seq data with DESeq2.” Genome Biology, 15, 550. doi: 10.1186/s13059-014-0550-8. 
- dplyr Wickham H, François R, Henry L, Müller K (2022). dplyr: A Grammar of Data Manipulation. https://dplyr.tidyverse.org, https://github.com/tidyverse/dplyr
- ggplot2 Wickham H (2016). ggplot2: Elegant Graphics for Data Analysis. Springer-Verlag New York. ISBN 978-3-319-24277-4, https://ggplot2.tidyverse.org
- optparse https://cran.r-project.org/package=optparse
- pheatmap https://cran.r-project.org/package=pheatmap
- PoiClaClu https://cran.r-project.org/package=PoiClaClu
- RColorBrewer https://cran.r-project.org/package=RColorBrewer
- tidyr Wickham H, Girlich M (2022). tidyr: Tidy Messy Data. https://tidyr.tidyverse.org, https://github.com/tidyverse/tidyr
- vcfUtils custom package, available here [vcfUtils.tar.gz](code/vcfUtils.tar.gz)


