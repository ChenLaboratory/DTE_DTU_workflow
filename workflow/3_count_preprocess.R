library(edgeR)
library(limma)
library(stringr)
library(RColorBrewer)

salmon_dir <- "results/salmon_output"
samples <- list.files(salmon_dir)

counts <- catchSalmon(file.path(salmon_dir, samples))

y <- DGEList(counts = counts$counts / counts$annotation$Overdispersion, 
             genes = counts$annotation)

colnames(y) <- sub("results/salmon_output/", "", colnames(y))

y$samples$group <- factor(rep(c("H1975", "HCC827"), each = 3))

transcript_annotation <- read.csv("data/reference/transcript_annotation.csv")

m <- match(rownames(y), transcript_annotation$transcript_id)
y$genes <- cbind(y$genes, transcript_annotation[m, ])

# filt <- rowSums(y$counts) > 10
filt <- filterByExpr(y, min.count = 3, min.total.count = 10)
table(filt)
y <- y[filt, , keep.lib.sizes=FALSE]

# keep1 <- y$genes$chromosome != "chrY"
# keep2 <- y$genes$gene_biotype %in% c("protein_coding", "lncRNA", "miRNA", "snRNA", "snoRNA")
# keep <- keep1 & keep2
# table(keep)
# y <- y[keep, , keep.lib.sizes = FALSE]

keep <- !grepl("pseudogene", y$genes$gene_biotype)
table(keep)
y <- y[keep, , keep.lib.sizes = FALSE]

y <- normLibSizes(y)
y$samples

# col <- brewer.pal(6, "Set2")
col <- c("red", "blue")
# pdf("results/figure/MDS.pdf", height = 6, width = 6)
plotMDS(y, col = col[as.numeric(y$samples$group)], pch = 16)
legend("topleft", levels(y$samples$group), text.col = col)
# dev.off()

# pdf("results/figure/MDS_wide.pdf", height = 4, width = 6)
# par(mar = c(5, 4, 4, 8))
# plotMDS(y, col = col[as.numeric(y$samples$group)], pch = 16)
# legend("topright", inset = c(-0.4, 0), levels(y$samples$group), text.col = col,
#        , xpd = TRUE)
# dev.off()