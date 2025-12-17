#!/bin/bash
#SBATCH --time=4:00:00
#SBATCH --cpus-per-task=1
#SBATCH --mem=2G
#SBATCH --job-name=sra_download

# NOTE: This script is designed to be used on Slurm. sra-toolkit is required to be pre-installed as a Slurm module.
# NOTE: if you want to run it directory on a UNIX environment without using Slurm, make sure sra-toolkit is installed and added to your PATH and skip line 12.

# USAGE: sbatch get_fastq.slurm <PROJECT_DIR> <SAMPLE> # replace <PROJECT_DIR> with your project directory and replace <SAMPLE> with the SRA accession of one sample

module load sra-toolkit

cd $1

# Extract fastq
fasterq-dump -O data/reads data/reads/$2/`ls ./$2`
gzip data/reads/$2\_1.fastq
gzip data/reads/$2\_2.fastq
