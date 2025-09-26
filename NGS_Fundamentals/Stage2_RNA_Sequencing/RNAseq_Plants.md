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

> a_thaliana_dea.R
```
# Set Working Directory
setwd("C:/Users/chloe/Desktop/Projects/HackBio/NGS_Fundamentals/Stage2_RNA_Sequencing")

# Install and Load libraries
if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")

if (!requireNamespace("DESeq2", quietly = TRUE)) BiocManager::install("DESeq2")
if (!requireNamespace("pheatmap", quietly = TRUE)) install.packages("pheatmap")

library(DESeq2)
library(pheatmap)

# Files
at_count <- read.delim('data/counts.txt', header = T)
at_meta <- read.delim('data/metadata.tsv', header = T)
at_meta$condition <- as.factor(at_meta$condition)

# Keep only the count data
raw_counts <- at_count[, at_meta$sample]
rownames(raw_counts) <- at_count$Geneid

# Create DeSeq dataset
dds <- DESeqDataSetFromMatrix(countData = raw_counts,
                              colData = at_meta,
                              design = ~condition)

# Perform DEA
dds <- DESeq(dds)
dea_results <- results(dds)

# Remove those with no pval then look at pval distribution
plot(density(x = na.omit(dea_results$pvalue)), main="")

# Plot DE genes
plot(x = dea_results$log2FoldChange, 
     y = -log10(dea_results$padj),
     pch = 19, 
     col = 'grey',
     ylim = c(0,20),
     ylab = 'Adjusted P-Value',
     xlab = 'Log2 FC')

abline(v = c(-2, 2), h = -log10(0.05), lwd = 0.5, lty = 2)

# Upregulated Genes
upregulated <- subset(dea_results, padj < 0.05 & log2FoldChange > 2)
points(upregulated$log2FoldChange,
       y = -log10(upregulated$padj), 
       pch = 19,
       col = 'salmon')

# Downregulated
downregulated <- subset(dea_results, padj < 0.05 & log2FoldChange < -2)
points(downregulated$log2FoldChange,
       y = -log10(downregulated$padj), 
       pch = 19,
       col = 'lightblue')

# Merge genes in one dataset
degs <- rbind(raw_counts[rownames(upregulated),], 
              raw_counts[rownames(downregulated),])
pheatmap(degs, 
         cluster_rows = F,
         cluster_cols = F,
         show_rownames = F,
         scale = 'row',
         show_colnames = T)

rownames(upregulated)
rownames(downregulated)

#exporting the files
write.csv(upregulated, 'result/upregulated.csv')
write.csv(downregulated, 'result/downregulated.csv')
write.csv(raw_counts, 'result/raw_counts.csv')
```

### Functional Enrichment
Differentially expressed genes from the RNA-seq analysis were analyzed using [ShinyGO v0.85](https://bioinformatics.sdstate.edu/go/) to identify enriched KEGG pathways. Significantly enriched pathways were selected based on an FDR-adjusted p-value cutoff of 0.05. Key statistics such as the number of genes involved, fold enrichment, and pathway coverage were recorded.

# Results
To understand the transcriptional changes captured by the RNA-seq dataset, we first identified differentially expressed genes (DEGs) between the experimental conditions. These DEGs were then further analyzed to uncover biological pathways and processes that are significantly affected, providing insight into the molecular mechanisms underlying the observed changes.

## Differentially Expressed Genes
Analysis of the RNA-seq dataset identified a total of 2,659 differentially expressed genes, of which 1,858 were upregulated and 801 were downregulated under the experimental conditions. The most significant DEGs were predominantly upregulated and included genes involved in transcription regulation, stress response, and secondary metabolism. Notably, the top 15 DEGs (Table 2) feature WRKY33 and WRKY40, key transcription factors in plant stress signaling, as well as genes encoding FAD-binding Berberine family proteins, cytochrome P450s, uridine diphosphate glycosyltransferases, and glutathione S-transferases, reflecting potential roles in defense and metabolic processes. Several highly upregulated genes remain uncharacterized, suggesting areas for future functional investigation.
_Note: A table of the top 100 DEGs can be found [here](result/top100deg.md) as a markdown table while all annotated DEGs can be found in this [csv file](result/degs.csv)._

**Table 2. Top 15 Differentially Expressed Genes for the _A. thaliana_ dataset**
| Ensembl Gene ID | Up/Down Regulated | Adjusted P-value | Symbol    | Description                                                                 |
|-----------------|-------------------|------------------|----------|-----------------------------------------------------------------------------|
| AT4G20830       | Up                | 1.9528E-220      | AT4G20830| FAD-binding Berberine family protein                                        |
| AT2G38470       | Up                | 8.6855E-192      | WRKY33   | WRKY DNA-binding protein 33                                                 |
| AT4G20835       | Up                | 2.0983E-191      | AT4G20835| Uncharacterized protein                                                     |
| AT3G26830       | Up                | 1.6487E-186      | PAD3     | Cytochrome P450 superfamily protein                                         |
| AT1G05680       | Up                | 2.7223E-163      | UGT74E2  | Uridine diphosphate glycosyltransferase 74E2                                |
| AT2G41105       | Up                | 4.2769E-157      | AT2G41105| Uncharacterized protein                                                     |
| AT2G18193       | Up                | 1.0317E-155      | AT2G18193| P-loop containing nucleoside triphosphate hydrolases superfamily protein    |
| AT4G02380       | Up                | 9.0478E-147      | SAG21    | Senescence-associated 21                                                    |
| AT1G19180       | Up                | 4.1741E-144      | JAZ1     | Jasmonate-zim-domain protein 1                                              |
| AT2G47000       | Up                | 1.8222E-143      | ABCB4    | ATP binding cassette subfamily B4                                           |
| AT2G41100       | Up                | 4.4672E-133      | TCH3     | Calcium-binding EF hand family protein                                      |
| AT1G72520       | Up                | 1.7717E-129      | LOX4     | PLAT/LH2 domain-containing lipoxygenase family protein                      |
| AT1G02930       | Up                | 3.2366E-124      | GSTF6    | Glutathione S-transferase 6                                                 |
| AT3G63380       | Up                | 7.5781E-121      | AT3G63380| ATPase E1-E2 type family protein / haloacid dehalogenase-like hydrolase fam |
| AT1G80840       | Up                | 3.9099E-119      | WRKY40   | WRKY DNA-binding protein 40                                                 |

## Enriched Pathways
KEGG pathway enrichment analysis revealed several pathways significantly overrepresented among the differentially expressed genes (Table X). Notably, genes involved in photosynthesis-antenna proteins (FDR = 4.9 × 10⁻², fold enrichment = 3.6), phenylalanine, tyrosine, and tryptophan biosynthesis (FDR = 2.0 × 10⁻³, fold enrichment = 3.1), and plant-pathogen interaction (FDR = 3.1 × 10⁻¹⁰, fold enrichment = 2.7) were enriched. Additional enriched pathways included glutathione metabolism and MAPK signaling in plants, indicating potential stress response and defense mechanisms captured by the dataset.

**Table 3. Enriched Pathways for the _A. thaliana_ dataset**
| Enrichment FDR | nGenes | Pathway Genes | Fold Enrichment | Pathways                                    |
|----------------|--------|---------------|----------------|---------------------------------------------|
| 4.9E-02        | 7      | 22            | 3.6            | Photosynthesis-antenna proteins              |
| 2.0E-03        | 16     | 56            | 3.1            | Phenylalanine tyrosine and tryptophan biosynthesis |
| 3.1E-10        | 56     | 216           | 2.7            | Plant-pathogen interaction                   |
| 1.5E-02        | 21     | 102           | 2.2            | Glutathione metabolism                       |
| 3.3E-03        | 29     | 144           | 2.1            | MAPK signaling pathway-plant                 |


# Discussion
The differential expression analysis revealed a predominance of upregulated genes in response to UV treatment in the vasculature, suggesting an active transcriptional response to stress. Many of the top DEGs, including **WRKY33, WRKY40, PAD3, and LOX4**, are known regulators of stress signaling and defense responses, consistent with prior observations that vascular tissues mount distinct transcriptional programs under abiotic stresses (Berkowitz et al., 2021). The KEGG pathway enrichment analysis further supports this, highlighting overrepresented pathways such as **phenylalanine, tyrosine, and tryptophan biosynthesis**, **plant-pathogen interactions**, and **MAPK signaling**, which are linked to secondary metabolism, reactive oxygen species signaling, and defense activation. These findings suggest that UV stress in the vasculature triggers both protective metabolic responses and signaling cascades that may coordinate stress adaptation, aligning with the tissue-specific responses previously described in Arabidopsis leaves. The presence of uncharacterized yet highly upregulated genes indicates that additional, potentially novel, components of the UV stress response remain to be explored in the vasculature.

While our analysis identifies a set of differentially expressed genes and enriched pathways in response to UV stress in the vasculature, these findings should be interpreted with some caution. Validation through independent experiments, such as qRT-PCR or additional RNA-seq datasets, would strengthen the reliability of the conclusions. Gene expression may also be influenced by other biological or environmental factors, including circadian rhythms or developmental stage, which could act as confounders. Furthermore, while the thresholds used for defining DEGs (log2 fold change > 2, adjusted p-value < 0.05) prioritize highly responsive genes, they may exclude genes with smaller but biologically relevant changes. Finally, the biological significance of enriched pathways could be explored more comprehensively through cross-referencing with literature and other functional databases, helping to confirm that the observed transcriptional responses are specifically linked to UV stress.
