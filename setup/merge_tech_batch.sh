zcat  data/reads/SRR31761441_1.fastq.gz  \
  data/reads/SRR31761442_1.fastq.gz \
  | gzip > data/reads/SRR31761441_merged_1.fastq.gz
zcat  data/reads/SRR31761441_2.fastq.gz  \
  data/reads/SRR31761442_2.fastq.gz \
  | gzip >  data/reads/SRR31761441_merged_2.fastq.gz
rm  data/reads/SRR31761441_1.fastq.gz \
  data/reads/SRR31761441_2.fastq.gz \
  data/reads/SRR31761442_1.fastq.gz \
  data/reads/SRR31761442_2.fastq.gz