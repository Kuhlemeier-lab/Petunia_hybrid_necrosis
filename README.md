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
- 250-300 bp inster cDNA library
- Sanger Illumina 1.9 encoding

### Raw reads

Are available on NCBI BioProject [PRJNA708139](https://www.ncbi.nlm.nih.gov/bioproject/PRJNA708139), under SRA accessions:
- SRR13907521: healthy individuals, P. exserta X P. axillaris N
- SRR13907522: necrotic individuals, P. exserta X P. axillaris N
- SRR13907523: healthy individuals, P. axillaris N X P. exserta
- SRR13907524: necrotic individuals, P. axillaris N X P. exserta

Forward reads are numbered 1, reverse are numbered 2.

### Read alignment

### Variant calling

### Bulk segregant analysis



## IL shallow sequencing

### Sequencing


### Raw reads

Have been uploaded to NCBI SRA under BioProject [](https://www.ncbi.nlm.nih.gov/bioproject/)


## RNAseq

To identify genes differentially expressed in the necrotic and healthy plants we performed an RNAseq experiment.
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

Have been uploaded to NCBI SRA under BioProject [](https://www.ncbi.nlm.nih.gov/bioproject/).

- SRR13861738:
- SRR13861739:
- SRR13861740:
- SRR13861741:
- SRR13861742:
- SRR13861743:
- SRR13861744:
- SRR13861745:
- SRR13861746:

### Alignment


### Read count

### Differential expression analysis

### Allele-specific expression analysis

### GO term analysis


