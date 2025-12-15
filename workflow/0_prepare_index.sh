# run the next line if gffread and salmon haven't been added to path; replace <$1> with the correct directory
# export PATH="$1/software/bin/:$PATH"

# 1.	Prepare transcriptome sequence. Extract the transcript sequences of all transcripts in the gene annotation GFF file using GffRead
gffread -w data/reference/transcripts.fa \
  -g data/reference/chm13v2.0_maskedY.fa \
  data/reference/chm13.draft_v2.0.gene_annotation.gff3

# 2.	Combine transcriptome and genome sequence in a new file as gentrome
cat data/reference/transcripts.fa \
  data/reference/chm13v2.0_maskedY.fa \
  > data/reference/gentrome.fa

# 3.	Extract the name of the genome targets to create a list of the decoys
grep "^>" data/reference/chm13v2.0_maskedY.fa \
  | cut -d " " -f 1\
  > data/reference/decoys.txt
sed -i.bak -e 's/>//g' data/reference/decoys.txt

# 4.	Build decoy-aware transcriptome index
salmon index -t data/reference/gentrome.fa \
  -d data/reference/decoys.txt \
  -p 8 \
  -i data/reference/salmon_index