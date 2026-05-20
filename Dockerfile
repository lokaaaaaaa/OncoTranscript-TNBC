FROM condaforge/mambaforge:23.3.1-1 as builder

LABEL maintainer="dlokaveenasri@gmail.com" \
      description="OncoTranscript-TNBC: STAR + DESeq2 + RF pipeline for TCGA-BRCA TNBC vs HR+" \
      version="1.0.0"

COPY environment.yml /tmp/environment.yml

RUN mamba env create -f /tmp/environment.yml && mamba clean -afy

ENV PATH /opt/conda/envs/rnaseq-pipeline-env/bin:$PATH
RUN echo "source activate rnaseq-pipeline-env" > ~/.bashrc

WORKDIR /workspace
