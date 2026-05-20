# Future Work

Planned extensions to OncoTranscript-TNBC, roughly ordered by priority.

## Near-term

- **Single-cell deconvolution** (CIBERSORTx, quanTIseq) to attribute the IFN-γ signal to specific tumor-infiltrating immune populations — the bulk RNA-seq signal is almost certainly driven by TILs but cannot be confirmed without deconvolution.
- **Survival association** of the top-5 classifier features against TCGA-BRCA clinical follow-up using multivariable Cox PH, adjusting for stage, age, and PAM50 subtype.
- **SHAP interpretability** on the Random Forest classifier to replace the current feature-importance scores with model-level SHAP values that account for feature interactions.

## Medium-term

- **Multi-cohort meta-analysis** spanning TCGA-BRCA + METABRIC + GSE96058 with ComBat-seq batch correction.
- **Graph-regularised classifier** that injects STRING-DB PPI priors as a regularisation term, replacing the current independent-feature Random Forest.
- **Nextflow Tower / Seqera** integration for cloud execution (AWS Batch, Google Cloud Life Sciences).

## Long-term

- **nf-core compatibility** layer to enable submission to the nf-core pipeline registry.
- **PAM50 subtyping** via genefu to replace IHC-based subtype labels with molecular subtypes.
- **Multi-modal integration** — combine RNA-seq signatures with DNA methylation and copy-number profiles from TCGA.
