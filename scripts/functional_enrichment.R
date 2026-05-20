#!/usr/bin/env Rscript
args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 2) {
  stop("Usage: Rscript functional_enrichment.R <deseq2_results.csv> <out_dir>", call. = FALSE)
}

library(clusterProfiler)
library(org.Hs.eg.db)
library(ggplot2)

deg_path <- args[1]
out_dir  <- args[2]

dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
deg_df <- read.csv(deg_path, row.names = 1)

sig_genes <- rownames(deg_df[which(deg_df$padj < 0.01 & abs(deg_df$log2FoldChange) > 2), ])

entrez_mapped <- bitr(sig_genes, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = org.Hs.eg.db)

ego <- enrichGO(gene          = entrez_mapped$ENTREZID,
                OrgDb         = org.Hs.eg.db,
                ont           = "BP",
                pAdjustMethod = "BH",
                pvalueCutoff  = 0.05,
                qvalueCutoff  = 0.05,
                readable      = TRUE)

write.csv(as.data.frame(ego), file = file.path(out_dir, "go_biological_processes.csv"))

if(nrow(as.data.frame(ego)) > 0) {
  dotplot(ego, showCategory = 15) +
    labs(title = "GO Biological Processes Enrichment") +
    theme_minimal(base_size = 12)
  ggsave(file.path(out_dir, "go_enrichment_dotplot.png"), width = 9, height = 7, dpi = 300)
} else {
  cat("No structural ontology groups passed strict multiple-testing corrections.\n")
}
