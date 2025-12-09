# print current time for timing purpose
print(Sys.time())

library(rtracklayer)
gff <- import("data/reference/chm13.draft_v2.0.gene_annotation.gff3")

transcripts <- gff[gff$type == "transcript", ]

transcript_annotation <- data.frame(
  transcript_id = transcripts$ID,
  gene_id = transcripts$gene_id,
  ensembl = substr(transcripts$source_transcript, 1, 15),
  gene = substr(transcripts$source_gene, 1, 15),
  symbol = transcripts$gene_name,
  gene_biotype = transcripts$gene_biotype,
  transcript_biotype = transcripts$transcript_biotype,
  chromosome = as.character(seqnames(transcripts))
)

transcript_annotation$symbol <- ave(
  transcript_annotation$symbol, 
  transcript_annotation$gene_id, 
  FUN = function(x) {
    if (length(unique(x)) > 1) {
      unique(x[!grepl("^MSTRG", x)])
    } else {
      x
    }
  }
)

write.csv(transcript_annotation, "data/reference/transcript_annotation.csv", 
          row.names = FALSE)
print(Sys.time())