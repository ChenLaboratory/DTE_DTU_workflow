#!/bin/bash
#SBATCH --time=48:00:00
#SBATCH --cpus-per-task=2
#SBATCH --mem=20G
#SBATCH --job-name=sra_download
#SBATCH --mail-user=dong.x@wehi.edu.au
#SBATCH --mail-type=BEGIN,END,FAIL

module load sra-toolkit

cd /vast/projects/lupus/xueyi/NatProt_LUAD

## Download 
prefetch  --option-file setup/SraAccList.csv --output-directory data/reads
