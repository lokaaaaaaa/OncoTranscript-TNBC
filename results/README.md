# Example results

Outputs from the run logged at `run_2024-11-12/` on the full TCGA-BRCA
discovery cohort (n=240). Committed as a reference so the figures and tables
in the README are verifiable without re-running the pipeline.

```
results/
└── run_2024-11-12/
    ├── qc/                    FastQC + STAR alignment summary
    │   └── multiqc_summary.tsv
    ├── counts/                featureCounts (excerpt)
    │   └── featurecounts_excerpt.tsv
    ├── deseq2/                DESeq2 results (top 2000 by padj)
    │   └── tnbc_vs_hrpos.tsv
    ├── gsea/                  clusterProfiler GSEA tables
    │   └── hallmark_gsea.tsv
    ├── ml/                    classifier metrics + feature importances
    │   └── metrics.json
    ├── figures/               PNGs used in the README
    └── pipeline_info/         Nextflow trace / timeline
        └── execution_trace.txt
```

A full run regenerates all of these — see `BENCHMARKS.md` for expected wall
times.
