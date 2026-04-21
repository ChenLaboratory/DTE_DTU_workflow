library(rtracklayer)
gff <- import("data/reference/chm13.draft_v2.0.gene_annotation.gff3")

transcripts <- gff[gff$type == "transcript", ]

transcript_annotation <- data.frame(
  gene_id = transcripts$gene_id,
  transcript_id = transcripts$ID,
  ensembl = substr(transcripts$source_transcript, 1, 15),
  gene = substr(transcripts$source_gene, 1, 15),
  symbol = transcripts$gene_name,
  gene_biotype = transcripts$gene_biotype,
  transcript_biotype = transcripts$transcript_biotype,
  chromosome = as.character(seqnames(transcripts))
)

symbols <- unlist(
  lapply(split(transcript_annotation, transcript_annotation$gene_id),
         function(df) {
           x <- df$symbol
           if (length(unique(x)) > 1) {
             unique(x[!grepl("^MSTRG", x)])[1]
           } else {
             unique(x)
           }
         })
)
transcript_annotation$symbol <- symbols[match(transcript_annotation$gene_id, 
                                              names(symbols))]

write.csv(transcript_annotation, "data/reference/transcript_annotation.csv", 
          row.names = FALSE)