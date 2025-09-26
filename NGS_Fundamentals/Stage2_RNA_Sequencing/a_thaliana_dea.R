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




