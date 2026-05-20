#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 3) {
  stop("Usage: Rscript deseq2_analysis.R <counts_matrix> <sample_sheet> <out_dir>", call. = FALSE)
}

library(DESeq2)
library(ggplot2)
library(pheatmap)
library(ggrepel)

counts_path  <- args[1]
metadata_path <- args[2]
out_dir       <- args[3]

dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

counts_raw <- read.table(counts_path, header = TRUE, skip = 1, row.names = 1, check.names = FALSE)
count_matrix <- counts_raw[, 6:ncol(counts_raw)]
colnames(count_matrix) <- gsub(".sortedByCoord.out.bam$", "", colnames(count_matrix))

metadata <- read.csv(metadata_path, header = TRUE, row.names = 1)
count_matrix <- count_matrix[, rownames(metadata)]

dds <- DESeqDataSetFromMatrix(countData = count_matrix,
                              colData = metadata,
                              design = ~ condition)

keep <- rowSums(counts(dds)) >= 10
dds <- dds[keep,]

dds <- DESeq(dds)
res <- results(dds, contrast = c("condition", "Tumor", "Normal"))
res_lfc <- lfcShrink(dds, coef = 2, type = "ashr")

write.csv(as.data.frame(res_lfc), file = file.path(out_dir, "deseq2_deg_results.csv"))

vsd <- vst(dds, blind = FALSE)
pca_data <- plotPCA(vsd, intgroup = "condition", returnData = TRUE)
percentVar <- round(100 * attr(pca_data, "percentVar"))

ggplot(pca_data, aes(PC1, PC2, color = condition)) +
  geom_point(size = 4, alpha = 0.8) +
  xlab(paste0("PC1: ", percentVar[1], "% variance")) +
  ylab(paste0("PC2: ", percentVar[2], "% variance")) +
  theme_minimal(base_size = 14) +
  scale_color_manual(values = c("Normal" = "#2c3e50", "Tumor" = "#e74c3c")) +
  labs(title = "Principal Component Analysis (PCA)", subtitle = "Tumor vs Normal Expression Profiles")
ggsave(file.path(out_dir, "pca_plot.png"), width = 7, height = 6, dpi = 300)

df_res <- as.data.frame(res_lfc)
df_res$significance <- "Not Significant"
df_res$significance[df_res$log2FoldChange > 2 & df_res$padj < 0.01] <- "Upregulated"
df_res$significance[df_res$log2FoldChange < -2 & df_res$padj < 0.01] <- "Downregulated"

df_res$gene <- rownames(df_res)
top_genes <- head(df_res[df_res$significance != "Not Significant", ][order(df_res$padj[df_res$significance != "Not Significant"]), ], 10)$gene

ggplot(df_res, aes(x = log2FoldChange, y = -log10(padj), color = significance)) +
  geom_point(alpha = 0.6, size = 1.5) +
  scale_color_manual(values = c("Upregulated" = "#e74c3c", "Downregulated" = "#3498db", "Not Significant" = "#bdc3c7")) +
  geom_text_repel(data = subset(df_res, gene %in% top_genes), aes(label = gene), size = 4, box.padding = 0.4, color = "black") +
  theme_classic(base_size = 14) +
  geom_vline(xintercept = c(-2, 2), linetype = "dashed", color = "gray40") +
  geom_hline(yintercept = -log10(0.01), linetype = "dashed", color = "gray40") +
  labs(title = "Differential Gene Expression Volcano Plot", x = "log2 Fold Change", y = "-log10 Adjusted P-value")
ggsave(file.path(out_dir, "volcano_plot.png"), width = 8, height = 7, dpi = 300)

cat("DESeq2 execution complete. Statistics and plots generated successfully.\n")
