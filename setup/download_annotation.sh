# downloading and preparing annotation resources
mkdir -p data/reference

wget -P data/reference -nv \
  https://s3-us-west-2.amazonaws.com/human-pangenomics/T2T/CHM13/assemblies/analysis_set/chm13v2.0_maskedY.fa.gz
gunzip data/reference/chm13v2.0_maskedY.fa.gz

wget  -P data/reference -nv \
  https://s3-us-west-2.amazonaws.com/human-pangenomics/T2T/CHM13/assemblies/annotation/chm13.draft_v2.0.gene_annotation.gff3

# won't need to run the next line because the file is already there if you have cloned the git repo
# wget -P data https://raw.githubusercontent.com/ChenLaboratory/DTE_DTU_workflow/refs/heads/main/data/targets.txt -nv