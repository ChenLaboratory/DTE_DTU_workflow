# Data preprocess using step 3 script
# source("workflow/3_count_preprocess.R")

# create design
design <- model.matrix(~ 0 + y$samples$group)
colnames(design) <- sub("y$samples$group", "", colnames(design), fixed = TRUE)

# model fitting
y <- estimateDisp(y, design = design)

fit <- glmQLFit(y, design = design)
# pdf("results/figure/BCV_QL.pdf", height = 8, width = 5)
pdf("results/figure/BCV.pdf", height = 4, width = 5)
# par(mfrow = c(2, 1))
par(mgp = c(2, 1, 0)) 
plotBCV(y)
dev.off()
pdf("results/figure/QL.pdf", height = 4, width = 5)
par(mgp = c(2, 1, 0)) 
plotQLDisp(fit)
dev.off()
