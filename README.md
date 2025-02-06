# Transcript-level RNA-seq data analysis using human T2T reference

This repository contains scripts to reproduce the analysis in the following workflow paper:

**Detecting differential transcript expression and usage in RNA-seq experiments with Salmon and edgeR using the human T2T reference**
Xueyi Dong, Junli Nie, Gordon K. Smyth, Yunshun Chen

## Example RNA-seq data availability

The example data used in this protocol is the Illumina RNA-seq data from [Dong *et al*](https://doi.org/10.1038/s41592-023-02026-3), which is available from Gene Expression Omnibus (GEO) under accession number [GSE172421](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE172421).

## Overview of the repository

* [`data`](data): This folder is where any data required to run the workflow should be stored. RNA-seq reads should be saved into `data/reads`, while reference genome and gene annotation files should be downloaded into `data/reference`. Here we provided the sample information and experimental design spreadsheet ([`data/targets.txt`](data/targets.txt)) for the example data.

* [`setup`](setup): This folder contains scripts for the experimental setup, including downloading and preparing data files and installing required software packages.

* [`workflow`](workflow): This folder contains scripts for all the steps in this analysis workflow.

## Instructions for running the workflow

All the scripts in this repository should be run from the project root directory.

To make sure the workflow can be reproduced, the users should follow the following order:

1. Clone this repository and navigate to the directory of the repository.

2. Install required software tools :
  -   Follow the instructions in the "Equipment setup" section in the paper to install SRA Toolkit, GffRead, Salmon and R.
  -   Run [`install_R_packages.R`](setup/install_R_packages.R) in R to install required R packages.

3. Download and prepare data:
  -   Download [human T2T-CHM13v2.0 reference genome sequence](https://s3-us-west-2.amazonaws.com/human-pangenomics/T2T/CHM13/assemblies/analysis_set/chm13v2.0_maskedY.fa.gz) and [annotation file](https://s3-us-west-2.amazonaws.com/human-pangenomics/T2T/CHM13/assemblies/annotation/chm13.draft_v2.0.gene_annotation.gff3) into the directory `data/reference`.
  -   Run [`download_sra.slurm`](setup/download_sra.slurm), [`get_fastq.slurm`](setup/get_fastq.slurm) and [`merge_tech_batch.sh`](setup/merge_tech_batch.sh) in order to download and prepare the example RNA-seq reads data.

4. Run the scripts of each step of the workflow under the [workflow](workflow) folder in order.

## Notes

* Due to the stochastic nature of Salmon’s quasi-mapping algorithm and Gibbs resampling, the result of each Salmon quantification run can be slightly different. This difference will impact downstream analysis such that different transcripts may be filtered out and slightly different number of differential expression or usage transcripts may be detected.
