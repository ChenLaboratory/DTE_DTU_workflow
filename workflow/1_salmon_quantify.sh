# run the next line if salmon hasn't been added to path; replace <$1> with the correct directory
# export PATH="$1/software/bin/:$PATH"

for SAMPLE in SRR1428605{7..9} SRR14286066 SRR14286069 SRR31761441_merged; do
  mkdir -p results/salmon_output/$SAMPLE
  salmon quant -i data/reference/salmon_index\
    -l A \
    -1 data/reads/$SAMPLE\_1.fastq.gz \
    -2 data/reads/$SAMPLE\_2.fastq.gz \
    -o results/salmon_output/$SAMPLE \
    -p  8 \
    --numGibbsSamples 50
done