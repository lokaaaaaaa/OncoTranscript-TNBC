# Implementation Notes

A running log of the practical bits that don't fit in the README. Kept honest
on purpose — the workflow is not magic and a few stages were genuinely annoying
to get right.

## STAR was the bottleneck

STAR 2-pass on GRCh38 + GENCODE v44 was by far the most expensive stage.

- Peak RSS hovered around **14 GB per sample**; anything below ~16 GB available
  RAM would OOM-kill the second pass.
- I originally tried `--genomeLoad LoadAndKeep` to share the index across
  parallel jobs, but on the shared cluster the index pages got evicted between
  tasks and runtimes actually got worse. Reverted to `NoSharedMemory` and let
  Nextflow handle parallelism via `maxForks`.
- Wall time per sample (50M PE reads, 8 threads): ~22 min alignment + ~4 min
  sorting/indexing.

## DESeq2 — batch surprised me

The TCGA discovery cohort spans multiple sequencing plates. A first naive
`~ subtype` design produced ~3,400 "DE" genes, of which a suspicious chunk were
pseudogenes and lincRNAs that loaded almost entirely on plate. Adding plate as
a covariate (`~ plate + subtype`) collapsed that to **1,247** genes and the
GSEA result became much cleaner — the IFN-γ and E2F signals strengthened while
the spurious lincRNA signal disappeared. Lesson: always check `colData` before
fitting.

`apeglm` shrinkage matters for the classifier too — ranking by unshrunken LFC
put a few high-variance, low-count genes at the top that didn't replicate in
GSE58135.

## GSE58135 quirks

- The GEO series is a mix of paired tumor/normal and tumor-only samples. The
  sample sheet ships only the **TNBC and HR+ primary tumors** (n=84) and drops
  the paired normals to keep the validation comparable to the discovery setup.
- Subtype labels in the series matrix are free-text ("triple negative",
  "ER+/PR+/HER2-", etc.). They're normalised in
  `data/metadata/sample_sheet.csv` (pre-normalised; see the cohort notes in the README) — worth re-checking if you swap the
  cohort.

## Classifier — validation-set drop is expected

5-fold CV on TCGA gives AUROC ≈ 0.94. On GSE58135 it drops to ≈ 0.91. That
~0.03 gap is the cohort-shift cost (different library prep, RIN distribution,
batch). I deliberately did **not** retune hyperparameters on GSE58135 — that
would defeat the purpose of an independent validation cohort.

## Things I'd do differently

- Switch STAR → `salmon` for quantification. STAR alignment is overkill when
  the downstream is gene-level DE; salmon would cut wall time by ~5x.
- Use `pyDESeq2` instead of shelling out to R from Nextflow. The R↔Python
  serialisation through TSV is the ugliest part of the pipeline.
- Track DAGs in MLflow rather than Nextflow's static `dag.html`.

## Limitations (be honest with reviewers)

- **Bulk RNA-seq only.** The IFN-γ signal almost certainly comes from
  tumor-infiltrating lymphocytes, but bulk data can't say so directly.
  Single-cell or CIBERSORTx deconvolution is the right next step (see
  `FUTURE_WORK.md`).
- **TCGA-BRCA is a US-centric cohort.** Generalisation to other ancestries is
  untested.
- **Subtype labels are immunohistochemistry-based**, not PAM50. A handful of
  "TNBC" samples are likely luminal-AR or HER2-low under a molecular
  classification.
- **Survival analysis is exploratory.** n=120 with limited follow-up is
  under-powered for definitive HR estimates; treat the KM curve as a
  hypothesis, not a claim.
