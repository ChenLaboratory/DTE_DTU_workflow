for SAMPLE in SRR1428605{7..9} SRR14286066 SRR14286069 SRR3176144{1..2} 
  do 
 	prefetch  $SAMPLE --output-directory data/reads 
 	fasterq-dump -O data/reads data/reads/$SAMPLE/ 
 	 pigz -p 8 data/reads/$SAMPLE\_1.fastq 
 	 pigz -p 8 data/reads/$SAMPLE\_2.fastq 
  done 