# Transcript-level RNA-seq data analysis using human T2T reference

This repository contains scripts to reproduce the analysis in the following workflow paper:

**Detecting differential transcript expression and usage in RNA-seq experiments with Salmon and edgeR using the human T2T reference**
Xueyi Dong, Junli Nie, Gordon K. Smyth, Yunshun Chen

## Example RNA-seq data availability

The example data used in this protocol is the Illumina RNA-seq data from [Dong *et al*](https://doi.org/10.1038/s41592-023-02026-3](https://doi.org/10.1038/s41592-023-02026-3), which is available from Gene Expression Omnibus (GEO) under accession number [GSE172421](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE172421).

## Overview of the repository

* [`data`](data): This folder is where any data required to run the workflow should be stored. RNA-seq reads should be saved into `data/reads`, while reference genome and gene annotation files should be downloaded into `data/reference`. Here we provided the sample information and experimental design spreadsheet ([`data/targets.txt`](data/targets.txt)) for the example data.

* [`setup`](setup): This folder contains scripts for the experimental setup, including downloading and preparing data files and installing required software packages.

* [`workflow`](workflow): This folder contains scripts for all the steps in this analysis workflow.

## Instructions for running the workflow

All the scripts in this repository should be run from the project root directory.

To make sure the workflow can be reproduced, the users should follow the following order:

1. Clone this repository

2. 
