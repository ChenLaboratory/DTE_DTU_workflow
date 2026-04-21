# model fitting using step 4 script
# source("workflow/4_fit_model.R")

# test for specific comparisons
contrast <- makeContrasts(HCC827 - H1975, levels = design)

res <- glmQLFTest(fit, contrast = contrast)

# test for DE
is_de <- decideTests(res)
summary(is_de)

# top DE results
topTags(res)

# save complete DTE results table into TSV file
tt <- topTags(res, n = Inf)
write.table(tt, file = "results/DTE_results.tsv", sep = "\t")

# Figure 4a
# pdf("results/figure/MD_small.pdf", height = 5, width = 4)
plotMD(res, main = "HCC827 vs H1975", cex = 0.3)
# dev.off()

# volcano plot
library(ggplot2)
# library(ggrepel)
fdr_cutoff <- 0.05
tt$table$de_status <- ifelse(tt$table$FDR > fdr_cutoff, "NotSig", 
                         ifelse(tt$table$logFC > 0, "Up", "Down"))

# Create a column for highlighting top 2 DE transcripts
tt$table$top2 <- ""
tt$table$top5[1:2] <- tt$table$transcript_id[1:2]
# pdf("results/figure/volcano_small.pdf", height = 5, width = 4)
ggplot(tt$table, aes(x = logFC, y = -log10(FDR),  label = top5)) +
  geom_point(aes(colour = de_status), size = 0.3) +
  geom_text(hjust = -0.05) +
  scale_colour_manual(values = c(
    "Up" = "red",
    "Down" = "blue",
    "NotSig" = "black"
  )) +
  geom_hline(yintercept = -log10(fdr_cutoff), linetype = "dashed")+
  theme_classic()
# dev.off()
 

