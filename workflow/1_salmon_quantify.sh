# run the next line if salmon hasn't been added to path; replace <$1> with the correct directory
# export PATH="$1/software/bin/:$PATH"

for SAMPLE in `ls data/reads/*1.fastq.gz \
  | xargs -n 1 basename \
  | sed 's/_1.fastq.gz$//'`;do
  mkdir -p results/salmon_output/$SAMPLE
  salmon quant -i data/reference/salmon_index\
    -l A \
    -1 data/reads/$SAMPLE\_R1.fastq.gz \
    -2 data/reads/$SAMPLE\_R2.fastq.gz \
    -o results/salmon_output/$SAMPLE \
    -p  8 \
    --numGibbsSamples 50
done