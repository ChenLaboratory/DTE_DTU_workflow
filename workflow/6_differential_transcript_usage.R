# model fitting using step 4 script
# source("workflow/4_fit_model.R")

# test for specific comparisons
# make contrast 
contrast <- makeContrasts(HCC827 - H1975, levels = design)
# test for differential transcript usage for the specified contrast
ds <- diffSplice(fit, contrast = contrast, 
                 geneid = "gene_id", exonid = "transcript_id")
# generate a list of all differentially spliced genes (simes test)
ts_simes <- topSplice(ds, 
                      test = "simes", number = Inf)
# # add gene symbol to the list for annotation
# ts_simes$symbol <- transcript_annotation$symbol[
#   match(ts_simes$gene_id, transcript_annotation$gene_id)]
# print the number of significant differentially spliced genes (simes test)
table(ts_simes$FDR < 0.05)
# print top differentially spliced genes (simes test)
topSplice(ds, test = "simes")
# generate a list of all differentially spliced transcripts
ts_transcript <- topSplice(ds, 
                              test = "t", number = Inf)
# print the number of significant differentially spliced transcripts
table(ts_transcript$FDR < 0.05)
# print top differentially spliced transcripts
topSplice(ds, test = "t")

#-------------- Results visualization

# Select example gene ZNF880 (CHM13_G0029458)

# make heatmaps to visualize the transcript usage
# extract counts of specific genes from samples in the comparison
m1 <- which(y$genes$gene_id == "CHM13_G0029458")
m2 <- colSums(t(design) * as.numeric(contrast))!= 0
dat <- as.data.frame(y$counts[m1, m2])
library(pheatmap)
annotation_col = data.frame(
  group = y$samples$group
)
ann_colors = list(
  group = c(H1975 = "red", HCC827 = "blue")
)
rownames(annotation_col) <- colnames(y)
pdf("results/figure/DTU_heatmap_scaled_count.pdf", height = 3, width = 5)
pheatmap(dat, scale = "column", cluster_cols = FALSE, 
         annotation_col = annotation_col, annotation_colors = ann_colors,
         main="ZNF880")
dev.off()


# make bar plots to visualize the transcript usage
dat <- dat * y$genes$Overdispersion[m1]
# transform extracted counts data into long format to prepare for visualization
dat$transcript <- rownames(dat)
dat <- reshape(dat, varying = list(names(dat)[-ncol(dat)]),
               v.names = "count", timevar = "sample",
               times = colnames(dat)[-ncol(dat)], idvar = "transcript",
               direction = "long")
# add sample information
dat$group <- y$samples$group[match(dat$sample, rownames(y$samples))]
# add transcript information
transcript_colors <- brewer.pal(
  sum(transcript_annotation$gene_id == "CHM13_G0029458"),
                                "Set2")
names(transcript_colors) <- transcript_annotation$transcript_id[
  transcript_annotation$gene_id == "CHM13_G0029458"]
# make stacked proportional bar plot to visualize the transcript usage
library(ggplot2)
pdf("results/figure/DTU_barplot.pdf", height = 3, width = 5)
ggplot(dat, aes(x = sample, y = count, fill = transcript)) +
  geom_bar(stat = "identity", position = position_fill()) +
  facet_grid(cols = vars(group), scales = "free") +
  labs(y = "transcript usage proportion", title = "ZNF880 transcript usage")+
  theme_bw() +
  theme(axis.text.x = element_text(angle = 30, hjust = 1)) +
  scale_fill_manual(values = transcript_colors)
dev.off()

# Visualize annotation transcripts
gff <- rtracklayer::import("data/reference/chm13.draft_v2.0.gene_annotation.gff3")
library(Gviz)
annotation_DTU <- gff[gff$gene_id == "CHM13_G0029458" & gff$type == "exon"]
annotation_DTU$transcript <- annotation_DTU$transcript_id
annotation_DTU <- annotation_DTU[order(annotation_DTU$transcript_id)]
axisTrack <- GenomeAxisTrack()
grtrack <- GeneRegionTrack(
  annotation_DTU,
  transcriptAnnotation = "transcript",
  fill = transcript_colors[as.character(annotation_DTU$transcript_id)]
)
# pdf("results/figure/DTU_transcript_model.pdf", height = 3, width = 8)
plotTracks(list(axisTrack, grtrack))
# dev.off()