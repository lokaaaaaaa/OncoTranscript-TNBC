# Benchmarks

Measured on the TCGA-BRCA discovery run (240 samples, 50M PE reads / sample,
GRCh38 + GENCODE v44). Hardware: 1× compute node, 32 vCPU, 128 GB RAM, NVMe scratch.

## Per-stage runtime & memory

| Stage              | Wall time | Peak RSS | Notes |
|--------------------|-----------|----------|-------|
| FastQC             | 12 min    | 1.1 GB   | embarrassingly parallel |
| Trim Galore        | 9 min     | 0.8 GB   | quality + adapter trim |
| STAR 2-pass align  | 1 h 47 m  | 14.2 GB  | dominant stage; see notes |
| featureCounts      | 18 min    | 2.4 GB   | gene-level, stranded |
| DESeq2             | 6 min     | 3.1 GB   | `~ plate + subtype` |
| clusterProfiler GSEA| 4 min    | 2.0 GB   | Hallmark + KEGG + GO:BP |
| RF classifier      | 2 min     | 0.9 GB   | 5-fold CV, 500 trees |
| **Total**          | **~3 h**  | 14.2 GB peak |  |

Smoke test (`tests/test_data/`, 4 samples × 200k reads): **~4 min** end-to-end
on a 2021 M1 laptop under Docker.

## Scientific metrics

| Metric                              | Value       |
|-------------------------------------|-------------|
| DEGs (padj<0.05, \|log2FC\|>1)      | **1,247**   |
| GSEA enriched Hallmark pathways     | 18 (FDR<0.05) |
| Top NES — IFN_GAMMA_RESPONSE        | +2.41 (padj 1.4e-6) |
| Top NES — ESTROGEN_RESPONSE_EARLY   | -2.62 (padj 9.0e-7) |
| RF 5-fold CV AUROC (TCGA)           | **0.969 ± 0.021** |
| RF AUROC on GSE58135                | **0.87**    |
| RF accuracy / F1 on GSE58135        | 0.76 / 0.88 |
| IFN-γ signature survival (Cox)      | HR 0.46 (0.24–0.87), log-rank p=0.018 |

Raw metric files: [`results/run_2024-11-12/ml/metrics.json`](results/run_2024-11-12/ml/metrics.json),
[`results/run_2024-11-12/gsea/hallmark_gsea.tsv`](results/run_2024-11-12/gsea/hallmark_gsea.tsv).
