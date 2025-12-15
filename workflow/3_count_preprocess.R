library(edgeR)
library(limma)

salmon_dir <- "results/salmon_output"
samples <- list.files(salmon_dir)
samples

targets <- read.delim("data/targets.txt", row.names = 1)
targets

counts <- catchSalmon(file.path(salmon_dir, samples))
y <- DGEList(counts = counts$counts / counts$annotation$Overdispersion, 
             genes = counts$annotation, group = factor(targets$CellLine))

colnames(y) <- rownames(targets)

transcript_annotation <- read.csv("data/reference/transcript_annotation.csv")

m <- match(rownames(y), transcript_annotation$transcript_id)
y$genes <- cbind(y$genes, transcript_annotation[m, ])

filt <- filterByExpr(y, min.count = 3, min.total.count = 10)
table(filt)
y <- y[filt, , keep.lib.sizes=FALSE]

keep <- !grepl("pseudogene", y$genes$gene_biotype)
table(keep)
y <- y[keep, , keep.lib.sizes = FALSE]

y <- normLibSizes(y)
y$samples


col <- c("red", "blue")
# Figure 3a
# pdf("results/figure/MDS.pdf", height = 6, width = 6)
plotMDS(y, col = col[as.numeric(y$samples$group)], pch = 16)
legend("topleft", levels(y$samples$group), text.col = col)
# dev.off()