# Smoke-test data

4-sample subset used by CI to validate samplesheet parsing, Nextflow configuration, and Python dependency imports in under 2 minutes — no STAR index or real FASTQ files required.

## Files

| File | Description |
|------|-------------|
| `sample_sheet.csv` | 4 samples: 2 TNBC + 2 HRpos, one from TCGA, one from GSE58135 |

## Running locally

```bash
bash tests/pipeline_tests/test_pipeline.sh
```

## Full pipeline run (requires references)

```bash
nextflow run main.nf \
  -profile docker \
  --samplesheet data/metadata/sample_sheet.csv \
  --genome_fasta references/GRCh38.primary_assembly.fa \
  --gtf references/gencode.v44.annotation.gtf \
  --star_index references/star_idx/ \
  --outdir results/
```
