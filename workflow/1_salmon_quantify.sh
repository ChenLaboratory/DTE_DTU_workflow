export PATH="/vast/projects/lupus/xueyi/NatProt_LUAD/software/bin/:$PATH"

mkdir -p results/salmon_output
for SAMPLE in `ls data/reads/*1.fastq.gz | xargs -n 1 basename| sed 's/_1.fastq.gz$//'`
do
  echo $SAMPLE
  date '+%A %W %Y %X'
  mkdir -p results/salmon_output/$SAMPLE
  salmon quant -i data/reference/salmon_index\
  -l A \
  -1 data/reads/$SAMPLE\_R1.fastq.gz \
  -2 data/reads/$SAMPLE\_R2.fastq.gz \
  --validateMappings -o results/salmon_output/$SAMPLE -p 8 --numGibbsSamples 35
done

date '+%A %W %Y %X'