# Data preprocess using step 3 script
# source("workflow/3_count_preprocess.R")

# create design
design <- model.matrix(~ 0 + y$samples$group)
colnames(design) <- sub("y$samples$group", "", colnames(design), fixed = TRUE)

# model fitting
y <- estimateDisp(y, design = design)

# Figure 3b
# pdf("results/figure/BCV.pdf", height = 4, width = 5)
# par(mgp = c(2, 1, 0)) 
plotBCV(y)
# dev.off()

fit <- glmQLFit(y, design = design)

# Figure 3c
# pdf("results/figure/QL.pdf", height = 4, width = 5)
# par(mgp = c(2, 1, 0)) 
plotQLDisp(fit)
# dev.off()