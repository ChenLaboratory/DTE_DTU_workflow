# model fitting using step 4 script
# source("workflow/4_fit_model.R")

# test for specific comparisons
# make contrast 
contrast <- makeContrasts(HCC827 - H1975, levels = design)

# test for differential transcript usage for the specified contrast
ds <- diffSplice(fit, contrast = contrast, 
                 geneid = "gene_id", exonid = "transcript_id")

# generate a list of all differentially spliced genes (simes test)
ts_simes <- topSplice(ds, test = "simes", number = Inf)
# print the number of significant differentially spliced genes (simes test)
table(ts_simes$FDR < 0.05)
# print top differentially spliced genes (simes test)
topSplice(ds, test = "simes")

# generate a list of all differentially spliced transcripts
ts_transcript <- topSplice(ds, test = "t", number = Inf)
# print the number of significant differentially spliced transcripts
table(ts_transcript$FDR < 0.05)
# print top differentially spliced transcripts
topSplice(ds, test = "t")

# calculate TPM and expression proportion of transcripts
TPM <- tpm(y, effective.tx.length = fit$genes$EffectiveLength,
           rta.overdispersion = fit$genes$Overdispersion)
TPMProp <- tpmProp(TPM, geneid = y$genes$gene_id)
# Select example gene ZNF880 (CHM13_G0029458)
g <- which(y$genes$gene_id == "CHM13_G0029458")
TPM[g, ]
TPMProp[g, ]

# Results visualization
# make heatmaps to visualize the transcript usage
library(pheatmap)
annotation_col = data.frame(group = y$samples$group)
ann_colors = list(group = c(H1975 = "red", HCC827 = "blue"))
rownames(annotation_col) <- colnames(y)
# Figure 4b
# pdf("results/figure/DTU_heatmap_scaled_TPM.pdf", height = 4, width = 5)
pheatmap(log(TPM[g, ] + 1), scale = "column", cluster_cols = FALSE, 
         annotation_col = annotation_col, annotation_colors = ann_colors,
         main="ZNF880")
# dev.off()


# make bar plots to visualize the transcript usage
dat <- as.data.frame(TPMProp[g, ])
# transform expression proportion data into long format to prepare for visualization
dat$transcript <- rownames(dat)
dat <- reshape(dat, varying = list(names(dat)[-ncol(dat)]),
               v.names = "proportion", timevar = "sample",
               times = colnames(dat)[-ncol(dat)], idvar = "transcript",
               direction = "long")
# add sample information
dat$group <- y$samples$group[match(dat$sample, rownames(y$samples))]
# add transcript information
library(RColorBrewer)
transcript_colors <- brewer.pal(
  sum(transcript_annotation$gene_id == "CHM13_G0029458"), "Set2")
names(transcript_colors) <- transcript_annotation$transcript_id[
  transcript_annotation$gene_id == "CHM13_G0029458"]
# make stacked proportional bar plot to visualize the transcript usage
library(ggplot2)
# Figure 4c
# pdf("results/figure/DTU_barplot.pdf", height = 4, width = 5)
ggplot(dat, aes(x = sample, y = proportion, fill = transcript)) +
  geom_bar(stat = "identity") +
  facet_grid(cols = vars(group), scales = "free") +
  labs(y = "transcript usage proportion", title = "ZNF880 transcript usage")+
  theme_bw() +
  theme(axis.text.x = element_text(angle = 30, hjust = 1)) +
  scale_fill_manual(values = transcript_colors)
# dev.off()

# Visualize annotation transcripts
# Skip the next line if code from stage 2 has been run
gff <- rtracklayer::import("data/reference/chm13.draft_v2.0.gene_annotation.gff3")

annotation_DTU <- gff[gff$gene_id == "CHM13_G0029458" & gff$type == "exon"]
annotation_DTU$transcript <- annotation_DTU$transcript_id
annotation_DTU <- annotation_DTU[order(annotation_DTU$transcript_id)]

library(Gviz)
axisTrack <- GenomeAxisTrack()
grtrack <- GeneRegionTrack(
  annotation_DTU,
  transcriptAnnotation = "transcript",
  fill = transcript_colors[as.character(annotation_DTU$transcript_id)]
)
# Figure 4d
# pdf("results/figure/DTU_transcript_model.pdf", height = 3, width = 8)
plotTracks(list(axisTrack, grtrack))
# dev.off()