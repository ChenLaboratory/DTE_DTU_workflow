#!/bin/bash
#SBATCH --time=24:00:00
#SBATCH --cpus-per-task=8
#SBATCH --mem=64G
#SBATCH --mail-user=dong.x@wehi.edu.au
#SBATCH --mail-type=END,FAIL

module load apptainer

# pull the docker image
apptainer pull docker://xueyidong/dte_dtu_workflow:latest

# download reference genome and gene annotation
echo "Setup 1: download reference genome and gene annotation"
setup/download_annotation.sh

# download RNA-seq data from SRA
echo "Setup 2: download RNA-seq data from SRA"
apptainer exec dte_dtu_workflow_latest.sif bash setup/download_and_prepare_sra.sh

# prepare the data: merge technical replicates
echo "Setup 3: prepare the data: merge technical replicates"
setup/merge_tech_batch.sh

# workflow steps
echo "Stage 0: preparation of transcriptome index "
apptainer exec dte_dtu_workflow_latest.sif bash workflow/0_prepare_index.sh

echo "Stage 1: transcript quantification"
apptainer exec dte_dtu_workflow_latest.sif bash workflow/1_salmon_quantify.sh

echo "Stage 2 to 6 are run in R. Here we knit the Rmd report."
apptainer exec dte_dtu_workflow_latest.sif Rscript -e "rmarkdown::render('workflow/workflow.Rmd')"