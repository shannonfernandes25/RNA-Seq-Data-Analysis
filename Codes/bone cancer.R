setwd("D:/project/bone_cancer")

library(DESeq2)

counts <- as.matrix(read.csv("counts.csv", row.names=1))
coldata <- read.csv("metadata.csv", row.names=1)

dds <- DESeqDataSetFromMatrix(countData = counts,
                              colData = coldata,
                              design = ~treatment)

dds <- DESeq(dds)

results <- results(dds)
head(results)


plotMA(results, main="DESeq2 MA Plot")

library(ggplot2)
results$log2FoldChange <- as.numeric(results$log2FoldChange)
ggplot(results, aes(x=log2FoldChange, y=-log10(padj))) +
  geom_point(alpha=0.4) + xlim(c(-5, 5)) + ylim(c(0, 10))


library(pheatmap)
topGenes <- head(order(results$padj), 20)
normalizedCounts <- counts(dds, normalized=TRUE)
pheatmap(normalizedCounts[topGenes,])

significant_genes <- results[which(results$padj<0.05),]
significant_genes


write.csv(significant_genes,"significant_genes.csv")

# PCA Analysis 
rld <- rlog(dds) 
pca <- prcomp(t(assay(rld))) 
pca_data <- data.frame(pca$x, treatment = colData(dds)$treatment) 

# PCA Plot using ggplot2 
pca_plot <- ggplot(pca_data, aes(x = PC1, y = PC2, color = treatment)) + 
  geom_point(alpha=0.7, size=3) + ggtitle("PCA Plot") + 
  xlab("Principal Component 1") + ylab("Principal Component 2") + 
  theme_minimal() 

print(pca_plot)
