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
curl -L ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR128/099/SRR12808499/SRR12808499.fastq.gz -o uvcc_3.fastq.gz
```

* Note that output files were renamed to reflect sample conditions (control vs uvc) and replicate number for easier reference later on.

## Preprocessing
Before any downstream analyses, the datasets first undergo some preprocessing steps to ensure that their quality is compatible with the next steps.

### Quality Control (QC)
QC is performed to check the raw reads for problems such as low-quality bases, leftover adapter sequences, or weird patterns.

The program FastQC was used to create QC reports for the raw reads while MultiQC was used to aggregate the reports into one output.

> qc.sh
```
#! /bin/bash

# create directories for output
mkdir -p qc

for filename in  datasets/*.fq; do
	fastqc -o qc/ $filename
done

multiqc datasets/
```

#### Initial QC Results


### Trimming

#### Re-QC
#### Re-QC Results

## Transcriptome Mapping
### Reference Genome
### Mapping with STAR

## Counting Transcripts
### Genome Annotation File
### FeatureCounts


## Differential Expression Analysis (DEA)
### Metadata
### Performing DESeq2
### Functional Enrichment

# Results
## Count Data Preview
## Differentially Expressed Genes
## Enriched Pathways

# Discussion

# Professional Profile
* Github: https://github.com/calicac001/HackBioInternship/NGS_Fundamentals/Stage3_RNA_Sequencing
* LinkedIn: www.linkedin.com/in/chloe-nichole-calica-3abb18262
