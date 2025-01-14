# Data preprocess using step 3 script
# source("workflow/3_count_preprocess.R")

# create design
design <- model.matrix(~ 0 + y$samples$group)
colnames(design) <- sub("y$samples$group", "", colnames(design), fixed = TRUE)

# model fitting
y <- estimateDisp(y, design = design)

fit <- glmQLFit(y, design = design, legacy = FALSE)
# pdf("results/figure/BCV_QL.pdf", height = 6, width = 5)
# par(mfrow = c(2, 1))
plotBCV(y)
plotQLDisp(fit)
# dev.off()
