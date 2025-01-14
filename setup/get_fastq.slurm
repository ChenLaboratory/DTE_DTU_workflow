#!/bin/bash
#SBATCH --time=4:00:00
#SBATCH --cpus-per-task=1
#SBATCH --mem=2G
#SBATCH --job-name=sra_download
#SBATCH --mail-user=dong.x@wehi.edu.au
#SBATCH --mail-type=BEGIN,END,FAIL

module load sra-toolkit

cd /vast/projects/lupus/xueyi/NatProt_LUAD

## Extract fastq
	fasterq-dump -O /vast/projects/lupus/xueyi/NatProt_LUAD/data/reads /vast/projects/lupus/xueyi/NatProt_LUAD/data/reads/$1/`ls ./$1`
	gzip /vast/projects/lupus/xueyi/NatProt_LUAD/data/reads/$1\_1.fastq
	gzip /vast/projects/lupus/xueyi/NatProt_LUAD/data/reads/$1\_2.fastq
