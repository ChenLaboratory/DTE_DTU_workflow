# Transcript-level RNA-seq data analysis using Salmon, edgeR v4 and the human T2T reference

This repository contains scripts to reproduce the analysis in the following workflow paper:

**Differential transcript expression and differential transcript usage using Salmon, edgeR v4 and the human T2T reference**

Xueyi Dong, Lizhong Chen, Junli Nie, Gordon K. Smyth, Yunshun Chen

## Example RNA-seq data availability

The example data used in this protocol is the Illumina RNA-seq data from [Dong *et al*](https://doi.org/10.1038/s41592-023-02026-3), which is available from Gene Expression Omnibus (GEO) under accession number [GSE172421](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE172421).

## Overview of the repository

* [`data`](data): This folder is where any data required to run the workflow should be stored. RNA-seq reads should be saved into `data/reads`, while reference genome and gene annotation files should be downloaded into `data/reference`. Here we provided the sample information and experimental design spreadsheet ([`data/targets.txt`](data/targets.txt)) for the example data.

* [`setup`](setup): This folder contains scripts for the experimental setup, including downloading and preparing data files and installing required software packages.

* [`workflow`](workflow): This folder contains scripts for all the steps in this analysis workflow.

## System requirements

### Operating system and hardware

This protocol is designed to be run on Linux operating system.  We recommand at least 20GB of RAM.

This protocol has been tested on Red Hat Enterprise Linux 9.3 (Plow) operating system.

### Software

The protocol is dependent on the following software:

* [SRA Toolkit software](https://github.com/ncbi/sra-tools) (version 3.1.0 or later)

* [GffRead software](https://github.com/gpertea/gffread) (version 0.12.7 or later)

* [Salmon software](https://combine-lab.github.io/salmon/) (version 1.10.0 or later) 

* [R](https://www.r-project.org) (version 4.5.2 or later)

* R packages:

    *	edgeR version 4.8.0

    * limma version 3.66.0

    * rtracklayer version 1.70.0

    * RColorBrewer version 1.1-3

    * ggplot2 version 4.0.1

    * Gviz version 1.54.0

    * pheatmap version 1.0.13

    * readr version 2.1.6

    * jsonlite version 2.0.0

The protocol has been tested with SRA Toolkit version 3.1.0, GffRead version 0.12.7, Salmon version 1.10.0, R version 4.5.2. The R session information from our test run can be found in the "Session information" section of the [expected output](results/expected_workflow_output.html).

## Instructions for running the workflow

### Running the workflow on example data

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

### Using the workflow on your data

When you want to use our workflow on your own data, we recommend the following:

* Choose a suitable version of the reference genome and annotation for your data.

* Prepare a target file to save your experimental design and sample information. The format can be found in this file: [data/targets.txt](data/targets.txt).

* Adjust the design matrix (stage 4, step 24) and contrasts (stage 5, step 29) according to your experimental design. We recommend the following article as a guide on how to set up your design and contrasts properly: **Law et al., A guide to creating design matrices for gene expression experiments, F1000Research, 2020**, [DOI: 10.12688/f1000research.27893.1](https://doi.org/10.12688/f1000research.27893.1).

## Expected output

The output from the test run can be found in [this file](results/expected_workflow_output.html). This file also includes the elapsed time of each R-based stage and the R session information in the test run.

## Notes

* Due to the stochastic nature of Salmon’s quasi-mapping algorithm and Gibbs resampling, the result of each Salmon quantification run can be slightly different. This difference will impact downstream analysis such that different transcripts may be filtered out and a slightly different number of differential expression or usage transcripts may be detected.
