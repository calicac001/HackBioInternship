---
title: HackBio Internship - Stage 2 
subtitle: RNAseq - Plants
author: Chloe Nichole Calica
date: September 23, 2025
---

# Project Background
Berkowitz et al. (2021) investigated how different tissues of the Arabidopsis thaliana leaf—epidermis, mesophyll, and vasculature—respond to multiple abiotic stresses. In our case, we will be determining the genes that respond to UV stress in the vasculature.

Specifically, you are to perform differential expression analysis in the vasculature of UV treated tissues with Water (control).

|Replicate | Control | UV-C Treated |
|----------|---------|--------------|
|1 |SRR12808527 |SRR12808497|
|2 |SRR12808528 |SRR12808498|
|3 |SRR12808529 |SRR12808499|

## Goals
* Complete Bash and R code for each step
* Understanding and Interpretation of results
* List of top 100 differentially expressed genes
* List of top 5 enriched pathways

# Methods
This section outlines the steps performed in order to obtain, process, and analyze the associated RNASeq data.

## Dataset
* Using the provided `SRR` accession codes, datasets were found in [SRA Explorer](https://sra-explorer.info/) for easy download.
* Following script was executed to download the datasets in the terminal:

> download.sh
```
#!/usr/bin/env bash

# Create dir for datasets
mkdir -p datasets
cd datasets/

# Control
curl -L ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR128/027/SRR12808527/SRR12808527.fastq.gz -o control_1.fastq.gz
curl -L ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR128/028/SRR12808528/SRR12808528.fastq.gz -o control_2.fastq.gz
curl -L ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR128/029/SRR12808529/SRR12808529.fastq.gz -o control_3.fastq.gz

# UV-C Treatment
curl -L ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR128/097/SRR12808497/SRR12808497.fastq.gz -o uvc_1.fastq.gz
curl -L ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR128/098/SRR12808498/SRR12808498.fastq.gz -o uvc_2.fastq.gz
curl -L ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR128/099/SRR12808499/SRR12808499.fastq.gz -o uvc_3.fastq.gz
```

* Note that output files were renamed to reflect sample conditions (control vs uvc) and replicate number for easier reference later on.

## Preprocessing
Before any downstream analyses, the datasets first undergo some preprocessing steps to check and adjust their quality.

### Quality Control (QC)
QC is performed to check the raw reads for problems such as low-quality bases, leftover adapter sequences, or weird patterns.

FastQC was used on the raw reads to create the quality control reports while MultiQC was used to aggregate the reports into one output.

> qc.sh
```
#! /bin/bash

# create directories for output
mkdir -p qc_reports

for filename in datasets/*.fastq.gz; do
	fastqc -o qc_reports/ $filename
done

multiqc qc_reports/
```

#### Initial QC Results
Initial QC of the raw reads indicated apparent failures in per-base sequence content, sequence duplication levels, and overrepresented sequences (Fig. 1). These patterns are common in RNA-seq data due to highly abundant transcripts as well as hexamer priming bias.The detailed QC report can be found [here](result/multiqc_report.html). 

[MultiQC Heatmap Summary of Initial QC](result/multiqc_summary.png)
**Figure 1. Summary Heatmap of the Initial QC Reports for the _A. thaliana_ Samples**

### Trimming
Trimming of adapters and low-quality bases was performed using FastP to remove technical artifacts while preserving biologically relevant sequences.”

> trim.sh
```
#!/bin/bash

#trim raw reads
mkdir -p trimmed

for filename in datasets/*.fastq.gz; do
    base=$(basename "$filename" .fastq.gz)
    fastp -i "$filename" -o "trimmed/${base}_trim.fastq.gz" -h "trimmed/${base}_report.html"
done
```

#### Re-QC
After trimming, the quality of the trimmed reads is checked by doing another QC step. Using a modified `qc.sh`, another round of FastQC was performed on the trimmed reads.

> re_qc.sh
```
#! /bin/bash

# create directories for output
mkdir -p qc_reports_2

for filename in trimmed/*.fastq.gz; do
	fastqc -o qc_reports_2/ $filename
done

multiqc qc_reports_2/
```

#### Re-QC Results
Detailed results of the MultiQC after trimming are shown [here](result/multiqc_report_trimmed.html). Most warnings and failures were not removed by the default FastP (Fig. 2), specifically in "Overrepresented Sequences," "Per-base Sequence Content," and "Sequence Duplication Levels."

Examination of the overrepresented sequences by BLAST revealed that they primarily corresponded to rRNA and chloroplast sequences, reflecting highly abundant biological transcripts rather than technical artifacts. The per-base sequence content bias, observed mainly at the first few bases of reads, is consistent with random hexamer priming commonly used in RNA-seq library preparation. Similarly, the high duplication levels largely result from these highly expressed transcripts. These features are typical in RNA-seq data and do not reflect poor library quality. Consequently, no further trimming of polyA tails or the first several bases was performed.

![MultiQC Summary of Trimmed](result/multiqc_summary_trimmed.png)
**Figure 2.  Summary Heatmap of the  QC Reports for the Trimmed Reads**

## Transcriptome Mapping
Despite the FastQC “fail” flags shown in the previous section, the trimmed reads are suitable for downstream transcriptome mapping with STAR, as the aligner tolerates soft-clipping of polyA tails and minor sequence bias at read starts, ensuring accurate alignment of mRNA reads.

### Reference Genome
Before mapping, the reference genome for _A. thaliana_ was obtained from Ensembl.

> get_genome.sh
```
#!/bin/bash

# Genome FASTA
wget -O a_thaliana_genome.fa.gz https://ftp.ensemblgenomes.ebi.ac.uk/pub/plants/release-62/fasta/arabidopsis_thaliana/dna/Arabidopsis_thaliana.TAIR10.dna.toplevel.fa.gz

# Unzip
gunzip a_thaliana_genome.fa.gz
```

### Mapping with STAR
RNA-seq reads were mapped to the Arabidopsis thaliana reference genome (TAIR10) using STAR, a splice-aware aligner optimized for transcriptome data. The raw reads in compressed FASTQ format were processed directly by STAR, which decompressed them on-the-fly using gunzip. A genome index was first generated from the reference to facilitate accurate alignment The output was sorted and stored as coordinate-sorted BAM files for downstream analyses, including transcript quantification and differential expression.

> mapping.sh
```
#! /bin/bash

# create genome index directory
STAR --runMode genomeGenerate --genomeDir genomeIndex --genomeFastaFiles a_thaliana_genome.fa

# create directories for output
mkdir -p mapped

for infile in trimmed/*.fastq.gz ; do
	outfile=$(basename "$infile" .fastq.gz)
	STAR --genomeDir genomeIndex \
        --readFilesIn $infile \
        --readFilesCommand gunzip -c \
        --outFileNamePrefix mapped/$outfile \
        --outSAMtype BAM SortedByCoordinate \
        --outSAMattributes All
done

```

## Counting Transcripts
Following mapping, transcript abundance was quantified by counting reads mapped to annotated genomic features. Using the coordinate-sorted BAM files generated by STAR, reads overlapping exons were assigned to genes based on the reference GTF annotation. This approach allows for accurate estimation of gene-level expression by accounting for reads spanning multiple exons and properly handling strand-specific information when applicable. The resulting count matrix ([`counts.txt`](result/counts.txt)) serve as the basis for differential expression analysis.

> count_reads.sh
```
#!/bin/bash

# Annotation GTF
wget -O a_thaliana_annotation.gff3.gz https://ftp.ensemblgenomes.ebi.ac.uk/pub/plants/release-62/gff3/arabidopsis_thaliana/Arabidopsis_thaliana.TAIR10.62.gff3.gz

# Unzip
gunzip a_thaliana_annotation.gff3.gz

# Run featureCounts
featureCounts -O -t gene -g ID -a a_thaliana_annotation.gff3 -o counts.txt mapped/*.bam
```

## Differential Expression Analysis (DEA)
DEA was performed to identify genes showing significant changes in expression between experimental conditions. Gene-level counts derived from STAR-aligned reads were used as input, allowing for robust statistical comparison across samples while accounting for variability in sequencing depth and biological replicates.

### Metadata
Although the datasets were renamed for easier identification, a metadata file was still constructed to be used for DEA. Properly structuring the metadata ensures that downstream analyses can accurately model the experimental design and interpret differential expression results.

**Table 1. Metadata Used for DEA**
|sample|condition|
|------|---------|
|control_1|control|
|control_2|control|
|control_3|control|
|uvc_1|uvc_treatment|
|uvc_2|uvc_treatment|
|uvc_3|uvc_treatment|

### Performing DESeq2
Following read alignment and transcript counting, the gene-level count data were imported into R for downstream analysis. Differential expression was carried out using DESeq2, which models RNA-seq count data while accounting for differences in sequencing depth and variability across samples. Experimental conditions from the metadata file were incorporated into the model, allowing for accurate identification of genes with significant expression changes between sample groups.

### Functional Enrichment

# Results
## Differentially Expressed Genes
## Enriched Pathways

# Discussion

# Professional Profile
* Github: https://github.com/calicac001/HackBioInternship/NGS_Fundamentals/Stage3_RNA_Sequencing
* LinkedIn: www.linkedin.com/in/chloe-nichole-calica-3abb18262
