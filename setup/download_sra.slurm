#!/bin/bash
#SBATCH --time=48:00:00
#SBATCH --cpus-per-task=1
#SBATCH --mem=20G
#SBATCH --job-name=sra_download

# NOTE: This script is designed to be used on Slurm. sra-toolkit is required to be pre-installed as a Slurm module.
# NOTE: if you want to run it directory on a UNIX environment without using Slurm, make sure sra-toolkit is installed and added to your PATH and skip line 13.

# USAGE: sbatch download_sra.slurm <PROJECT_DIR> # replace <PROJECT_DIR> with your project directory


module load sra-toolkit

cd $1

# Download 
prefetch  --option-file setup/SraAccList.csv --output-directory data/reads
